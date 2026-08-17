Module RayTracing
    use constants
    use Quadrature
    use IO
    use Geometry
    use BasicMathFunctions
    use Sources
    Implicit None
    
    ! NOTES
    ! For AddAngularFluence I will need to use a different set of raytracing nodes, because
    ! doing the quadrature then mass-eliminating does not get you the unc. at the spatial points
    ! unless the unc. is linear in the element (which we explicitly assume it isn't)
    ! This error propagates because you need to do EI and PI
    ! Idea for AddAngularFluence is to trace to points that are infinitesimally closer to the
    ! interior of the element. Essentially taking:
    ! r = r^{e}_{k} - \epsilon * (r^{e}_{CM} - r^{e}_{k}) / |r^{e}_{CM} - r^{e}_{k}|
    ! This will prevent any ambiguity for, e.g., beams traveling perfectly along a surface
    ! (which is still possible but dramatically less likely, if this is not done then there
    ! are meshes + beams where this will ALWAYS happen multiple times)
    
Contains

Subroutine RayTraceBeam (beam)
    Implicit None
    
    !  ================================================================
    !    This subroutine determines the path lengths of a beam in the 
    !    different materials in the mesh. The optical path length can 
    !    only be determined with cross sections (specifically, the 
    !    attenuation coefficient), so that is done elsewhere.
    !  ================================================================
    
    Type (ExternalBeam),             Intent (InOut) :: beam
    
    Integer                                         :: e, eb, f, i, iRay, im2sd, k, kq, m, q, sdof
    Real (KREAL),        Parameter                  :: TOL = 1.0e-8_KREAL ! EPSILON(ONE)
    Integer                                         :: Nm
    Integer                                         :: NE
    Integer                                         :: NENK
    Integer                                         :: NK
    Integer                                         :: NR
    Integer,             Allocatable                :: qtmp     (:)
    Integer,             Allocatable                :: einb     (:)
    Real (KREAL)                                    :: rs       (2)
    Real (KREAL)                                    :: vs       (2)
    Real (KREAL)                                    :: cqs      (2)
    Real (KREAL)                                    :: wqs      (2)
    Real (KREAL)                                    :: R0       (3)
    Real (KREAL)                                    :: Rp       (3)
    Real (KREAL)                                    :: r        (3)
    Real (KREAL)                                    :: wq       (4)
    Real (KREAL),        Allocatable                :: pathlen  (:)
    Real (KREAL)                                    :: re       (3, 4)
    Real (KREAL)                                    :: cq       (3, 4)
    Real (KREAL)                                    :: v        (3, 4)
    
    Nm   = beam%mesh%Nm
    NE   = beam%mesh%NE
    NENK = beam%mesh%NENK
    
    ALLOCATE(pathlen(Nm))
    
    !  =================================================================
    !    Prepare the rays object with necessary addressing information  
    !  =================================================================
    
    NR = COUNT(beam%quadmask) ! Number of quadrature points for which a ray is needed (= number of rays)
    
    ALLOCATE(beam%rays(NR))
    
    ! Allocate pathlen
    do iRay = 1, NR
        ALLOCATE(beam%rays(iRay)%pathlen(Nm), source = ZERO)
    end do
    
    ALLOCATE(einb, source=PACK([(e,e=1,NE)], mask=beam%elemmask))
    
    !  =====================
    !    Begin ray tracing  
    !  =====================
    
    if (beam%mesh%slab) then
        ! Slab case is easy
        
        do iRay = 1, NR
            beam%rays(iRay)%q = iRay ! Since all slab points are in the beam
        end do
        
        call PrepLinearQuadrature (cqs, wqs)
        
        do iRay = 1, NR
            ! Visit a ray (mesh sdof)
            
            do e = 1, (iRay + 1) / 2 - 1 ! Integer division gives the element in which q resides
                m = beam%mesh%el(e)%mat
                beam%rays(iRay)%pathlen(m) = beam%rays(iRay)%pathlen(m) + beam%mesh%el(e)%vol ! "Volume" is actually length of the element
            end do
            
            ! Now add the coordinate of the quad point relative to the element starting point
            e = (iRay + 1) / 2
            m = beam%mesh%el((iRay + 1) / 2)%mat
            
            rs(1) = beam%mesh%rg(1,e)
            rs(2) = beam%mesh%rg(1,e + 1)
            vs = MapToLinearQuad (cqs, rs)
            
            if (odd(iRay)) then
                beam%rays(iRay)%pathlen(m) = beam%rays(iRay)%pathlen(m) + (vs(1) - beam%mesh%rg(1,e))
            else
                beam%rays(iRay)%pathlen(m) = beam%rays(iRay)%pathlen(m) + (vs(2) - beam%mesh%rg(1,e))
            end if
            
        end do
        
    else if (streq(beam%fldgeo%angdist, 'spherical')) then
        ! Spherical case: Ray origin is constant
        R0 = beam%fldgeo%origin
        
        call PrepTetrahedralQuadrature (cq, wq)
        
        iRay = 0
        do eb = 1, SIZE(einb)
            ! ----- HERE IS WHERE THE TETRAHEDRAL CASE IS ENFORCED -----
            e = einb(eb)
            
            re(1:3, 1:4) = beam%mesh%rg(1:3, beam%mesh%el(e)%node(1:4))
            v = MapToTetrahedralQuad (cq, re)
            
            do kq = 1, 4
                q = beam%mesh%offset(e) + kq
                if (.not. beam%quadmask(q)) cycle
                
                ! If you made it to this point, you are a valid ray, so iterate the ray index
                iRay = iRay + 1
                
                ! Assign the q index of this ray
                beam%rays(iRay)%q = q
                
                ! Express the node relative to the origin
                r = v(1:3,kq) - R0
                
                ! Ray trace
                call TraceVector (Nm, R0, r, beam%mesh, TOL, pathlen)
                
                beam%rays(iRay)%pathlen = pathlen
                
            end do
            ! -----------------------------------------------------------
        end do
        
    else
        ! Planar case: Ray origin can be interpreted one of two ways:
        !              1. It is arbitrarily far from the mesh, along the beam axis
        !              2. It is given by back-projecting from the quadrature point to the plane of the beam
        ! But the first one doesn't really provide any useful identities that can be used below so we'll just
        ! do the 2nd. We take the plane of the source to be given by the beam origin and beam axis.
        R0 = beam%fldgeo%origin
        
        call PrepTetrahedralQuadrature (cq, wq)
        
        iRay = 0
        do eb = 1, SIZE(einb)
            ! ----- HERE IS WHERE THE TETRAHEDRAL CASE IS ENFORCED -----
            e = einb(eb)
            
            re(1:3, 1:4) = beam%mesh%rg(1:3, beam%mesh%el(e)%node(1:4))
            v = MapToTetrahedralQuad (cq, re)
            
            do kq = 1, 4
                q = beam%mesh%offset(e) + kq
                if (.not. beam%quadmask(q)) cycle
                
                ! If you made it to this point, you are a valid ray, so iterate the ray index
                iRay = iRay + 1
                
                ! Assign the q index of this ray
                beam%rays(iRay)%q = q
                
                ! Construct the origin by projecting onto the basis e1 and e2 then taking the beam axis component to zero
                ! We take Rp in the beam coordinate system as R0 + (a, b, 0),
                ! where a = (v - R0) * e1, b = (v - R0) * e2.
                ! Use r as a temporary for v - R0
                r  = v(1:3,kq) - R0
                Rp = R0                                              &
                   + DOT_PRODUCT(r, beam%fldgeo%e1) * beam%fldgeo%e1 &
                   + DOT_PRODUCT(r, beam%fldgeo%e2) * beam%fldgeo%e2
                
                ! Express the node relative to the origin
                r = v(1:3,kq) - Rp
                
                ! Ray trace
                call TraceVector (Nm, Rp, r, beam%mesh, TOL, pathlen)
                
                beam%rays(iRay)%pathlen = pathlen
                
            end do
            ! -----------------------------------------------------------
        end do
        
    end if
    
End Subroutine

Subroutine TraceVector (Nm, R0, r, mesh, TOL, pathlen)
    Implicit None
    Integer,           Intent (In)  :: Nm           ! Number of materials in the mesh
    Real (KREAL),      Intent (In)  :: R0      (3)  ! Ray origin
    Real (KREAL),      Intent (In)  :: r       (3)  ! Point to be traced
    Type (MeshClass),  Intent (In)  :: mesh         ! Mesh
    Real (KREAL),      Intent (In)  :: TOL          ! Numerical tolerance
    
    Real (KREAL),      Intent (Out) :: pathlen (Nm)
    
    Integer                         :: e, f, fg, freg, i, il, k0, k1, k2, mat, reg
    Logical                         :: success
    Integer                         :: nkeep
    Integer                         :: sl
    Real (KREAL)                    :: length
    Real (KREAL)                    :: det
    Real (KREAL)                    :: u
    Real (KREAL)                    :: v
    Real (KREAL)                    :: t
    Real (KREAL)                    :: vert0 (3)
    Real (KREAL)                    :: vert1 (3)
    Real (KREAL)                    :: vert2 (3)
    Real (KREAL)                    :: edge1 (3)
    Real (KREAL)                    :: edge2 (3)
    Real (KREAL)                    :: pvec  (3)
    Real (KREAL)                    :: tvec  (3)
    Real (KREAL)                    :: qvec  (3)
    Real (KREAL)                    :: n     (3)
    Real (KREAL)                    :: nref  (3)
    Real (KREAL),      Allocatable  :: llist (:)
    Character (LEN=:), Allocatable  :: errmess
    
    pathlen = ZERO
    length  = NORM2(r)
    success = .FALSE.
    
    do reg = 1, SIZE(mesh%region)
        
        if (Nm > 1) then
            ALLOCATE(llist(0))
            mat = mesh%region(reg)%mat
        else
            mat = 1 ! I think this doesn't get used though
        end if
        
        ! Find the interface face which intersects the ray connecting r and R0
        do freg = 1, SIZE(mesh%region(reg)%fg)
            fg = mesh%region(reg)%fg(freg)
            
            k0 = mesh%face(fg)%node(1)
            k1 = mesh%face(fg)%node(2)
            k2 = mesh%face(fg)%node(3)
            
            vert0 = mesh%rg(:,k0)
            vert1 = mesh%rg(:,k1)
            vert2 = mesh%rg(:,k2)
            
            edge1 = vert1 - vert0
            edge2 = vert2 - vert0
            
            ! INSTEAD OF CHECKING BACKFACING STUFF, Can't I just limit the length to norm2(r - R0)??? Although an issue is that nodes on the outgoing boundary would get 0 pathlength... Certainly treatable...
            e   = mesh%face(fg)%sharedby(1)
            f   = FINDLOC(mesh%el(e)%face, fg, DIM=1) ! Could be slow, store this in mesh%face(fg)??? ! THIS IS VERY SLOW I'D BET
            n   = mesh%el(e)%n(:,f)
            
            ! Check if face is backfacing AND Nm = 1. If so, cycle
            det = DOT_PRODUCT(n, r)
            if (Nm == 1) then
                if (det > -TOL) cycle ! Change this back to det > TOL? This is weeding out parallel faces
            else
                if (ABS(det) < TOL) cycle ! Weed out parallel faces.. They're technically possible even in the quadrature scheme but not likely, and in any case we'll choose to ignore them
            end if
            
            nref = CROSS_PRODUCT(edge1, edge2)
            
            if (DOT_PRODUCT(n, nref) < ZERO) then
                ! Swap the definitions of edges 1 and edges 2, because the way they were defined gives the wrong sign to
                ! following values (Could I avoid this somehow? Choice of edge1 and edge2 definitely matters)
                nref  = edge1
                edge1 = edge2
                edge2 = nref
            end if
            ! shouldve documented what i was doing below...
            pvec = CROSS_PRODUCT(r, edge2)
            
            det = DOT_PRODUCT(edge1, pvec)
            
            if (ABS(det) < TOL) cycle
            
            tvec = R0 - vert0
            
            u = DOT_PRODUCT(tvec, pvec)
            
            qvec = CROSS_PRODUCT(tvec, edge1)
            
            v = DOT_PRODUCT(r, qvec)
            
            if (det > ZERO) then
                if (u < -TOL .or. u > det + TOL) cycle
                if (v < -TOL .or. u + v > det + TOL) cycle
            else
                if (u > TOL .or. u < det - TOL) cycle
                if (v > TOL .or. u + v < det - TOL) cycle
            end if
            
            t = DOT_PRODUCT(edge2, qvec) / det
            
            if (Nm == 1) then
                success = .TRUE.
                exit
            else
                if (t > ONE) cycle ! Do not keep the outgoing boundary, which will otherwise be recorded with t > 1
                llist = [llist, t] ! t is just intersection in units of length
            end if
            
        end do
        
        if (Nm == 1) exit ! If Nm == 1 and you reach this point you'll definitely hit the error condition anyway
            ! If there is just one intersection, then check that it is nonzero
            ! If it is close to zero then just leave
            ! NOTE that if tracing QUADRATURE rays then there will be no points precisely on a face,
            ! so we can do this
        if (SIZE(llist) == 1) then
            if (ABS(llist(1)) < 1.0e-12_KREAL) then
                DEALLOCATE(llist)
                cycle ! Region contributes nothing
            end if
        else if (SIZE(llist) == 0) then
            DEALLOCATE(llist)
            cycle ! Region contributes nothing
        end if
        
        call quicksort (llist)
        
        ! Edge case - remove any values that are identical. This means a ray happened to pass dead center through multiple faces
        nkeep = 1
        
        do i = 2, SIZE(llist)
            if (ABS(llist(i) - llist(nkeep)) > 1.0e-12_KREAL) then ! 1.0e-12 is somewhat arbitrary. It's supposed to be smaller than ~Angstrom but larger than precision error
                nkeep = nkeep + 1
                llist(nkeep) = llist(i)
            end if
        end do
        
        llist = llist(1:nkeep)
        
        sl = SIZE(llist) / 2 ! Integer division to take floor
        
        do il = 1, sl
            pathlen(mat) = pathlen(mat) + (llist(2 * il) - llist(2 * il - 1))
        end do
        
        if (ODD(SIZE(llist))) then
            ! If llist is odd then this ray is within this region.
            ! Add the last path
            pathlen(mat) = pathlen(mat) + (ONE - llist(SIZE(llist))) ! In units of length
        end if
        
        DEALLOCATE(llist)
        
    end do
    
    if (Nm == 1) then
        if (success) then
            pathlen(1) = (ONE - t) * length ! We take 1 - t to find the fracton of length that occurred within the mesh
        else
            errmess = 'RayTracing.f90: TraceVector: While ray tracing no intersections were found. ' //&
                      'Check your mesh file for errors.'
            call stophere (errmess)
        end if
    else
        pathlen = pathlen * length
    end if
    
End Subroutine

!Subroutine TraceVectorHeterogeneous (R0, r, mesh, TOL, ells) ! RENAME
!    Implicit None
!    Real (KREAL),     Intent (In)  :: R0  (3) ! Beam origin
!    Real (KREAL),     Intent (In)  :: r (3) ! Beam direction
!    Type (MeshClass), Intent (In)  :: mesh   ! Mesh
!    Real (KREAL),     Intent (In)  :: TOL
!
!    Real (KREAL),     Intent (Out) :: ells (Nm)
!
!    Integer                        :: f, i, mat
!    Logical                        :: success
!    Integer                        :: N
!    Integer                        :: sizeofmat
!    Real (KREAL)                   :: ell
!    Real (KREAL)                   :: det
!    Real (KREAL)                   :: u
!    Real (KREAL)                   :: v
!    Real (KREAL)                   :: tval
!    Real (KREAL)                   :: t0
!    Real (KREAL)                   :: tN
!    Real (KREAL)                   :: vert0  (3)
!    Real (KREAL)                   :: vert1  (3)
!    Real (KREAL)                   :: vert2  (3)
!    Real (KREAL)                   :: edge1  (3)
!    Real (KREAL)                   :: edge2  (3)
!    Real (KREAL)                   :: pvec   (3)
!    Real (KREAL)                   :: tvec   (3)
!    Real (KREAL)                   :: qvec   (3)
!    Real (KREAL)                   :: Normal (3)
!    Real (KREAL),     Allocatable  :: t      (:)
!    Real (KREAL),     Allocatable  :: B_t    (:)
!    Character (100)                :: errmess
!    
!    ! Algorithm:
!    ! We only trace through faces that represent an interface between materials (surface faces count,
!    ! as they are essentially an interface between a material and vacuum).
!    ! We just find the intersections (if more than one exist) with these interface faces, and then
!    ! we sort the distances and use that sorted array to determine the travel distances in the various
!    ! materials. This is much faster, of course, than tracing through every face in the mesh.
!    
!    ! TO FIND INTERFACES:
!    ! Visit each element and check its neighbors. 
!    ! If the element is on the bdy, it's an interface (with vacuum).
!    ! If a neighbor has a different element, that face is an interface
!    
!    ! Yes that makes sense. Just outline all interface elements, and then ray trace through those,
!    ! then organize the travel distances, as this subroutine does. Perfect, just need to translate Cfkf and matfc.
!    
!    t0 = 0
!    tN = 1
!    
!    ells = 0
!    ell = norm2(r)
!    success = .TRUE. ! Not a mistake. Later on taking success = success .AND. (something).
!    do mat = 1, Nm ! Restructure? Does this loop permit multiple intersections? I think so, see the B_t stuff.
!        ALLOCATE(t(0))
!        
!        sizeofmat = findloc(matfc(:,mat), 0, dim = 1) - 1 ! Why the -1??
!        
!        do f = 1, sizeofmat
!            vert0 = rglobal(:,Cfkf(matfc(f,mat),1))
!            vert1 = rglobal(:,Cfkf(matfc(f,mat),3)) ! SWAPPED ! WHAT DID THIS MEAN??
!            vert2 = rglobal(:,Cfkf(matfc(f,mat),2)) ! SWAPPED ! WHAT DID THIS MEAN??
!
!            edge1 = vert1 - vert0
!            edge2 = vert2 - vert0
!
!            Normal = cross_product(edge1,edge2) ! Do I have normals already??? I think so
!            det    = DOT_PRODUCT(Normal, r)
!            ! NOTES: No culling since material interfaces have normals with direction up to +/- 1,
!            ! also removing parallel faces.
!            if (abs(det) < TOL) cycle
!
!            pvec = cross_product(r, edge2)
!
!            det = dot_product(edge1, pvec)
!
!            if (det < TOL) cycle ! If stuff doesn't work, investigate this step
!
!            tvec = R0 - vert0
!
!            u = dot_product(tvec, pvec)
!
!            if (u < -TOL .or. u > det + TOL) cycle
!
!            qvec = cross_product(tvec, edge1)
!
!            v = dot_product(r, qvec)
!
!            if (v < -TOL .or. u + v > det + TOL) cycle
!
!            tval = dot_product(edge2, qvec)/det
!
!            allocate(B_t(size(t)+1))
!            B_t(1:size(t)) = t
!            B_t(size(t)+1) = tval
!            call move_alloc(B_t, t)
!        end do
!
!        success = success .AND. size(t) > 0
!
!        call qsort(t)
!
!        ! I used stupid conventions again. When not tired, redo conventions.
!        N = size(t) + 1
!
!        allocate(B_t(N+2))
!        B_t(1) = t0
!        B_t(2:N+1) = t
!        B_t(N+2) = tN
!        call move_alloc(B_t, t)
!
!        ells(mat) = ell
!        do i = 1, floor(0.5*(N+1))
!            ells(mat) = ells(mat) - (t(2*i-1 + 1)-t(2*i-2 + 1))*ell ! +1 to shift indexing with t(0) stuff
!        end do
!
!        deallocate(t)
!    end do
!
!    if (.not. success) then
!        errmess = 'Error while ray tracing mesh - no intersections found. Check your mesh files'
!        call stophere (errmess)
!    end if
!    
!End Subroutine

End Module