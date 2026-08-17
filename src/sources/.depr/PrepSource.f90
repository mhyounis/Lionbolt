Submodule (Sources) submodPrepSource
Contains
Module Subroutine PrepSource (self)
    
    ! DEPRECATED BUT CONTAINS A LOT OF CODE THAT WAS TRICKY TO WRITE. MAY BE USED IN THE FUTURE.
    ! Nominally gave a storage-saving representation of the external beam source, that could be
    ! transported to give the uncollided angular fluence.
    ! This however is prone to ray artifacts and is generally not superior to raytracing.
    
    use constants
    use IO
    use Quadrature
    use ShapeFunctions
    use FEAInnerProducts
    use Jacobians
    use Sources
    Implicit None
    Class (ExternalBeam), Intent (InOut) :: self
    
    Integer                              :: b, bloc, e, eb, f, fg, i, k, kb, kf, kg, kq
    Real (KREAL),         Parameter      :: TOL = 1.0e-12_KREAL
    Logical                              :: proceed
    Logical                              :: spherical = .TRUE.
    Logical                              :: anyfaces
    Integer                              :: d
    Integer                              :: idx
    Integer                              :: NE
    Integer                              :: NEb
    Integer                              :: NFe
    Integer                              :: NENK
    Integer                              :: NI
    Integer                              :: NKe
    Integer                              :: NKf
    Integer                              :: NQ
    Integer,              Allocatable    :: einb    (:)
    Integer,              Allocatable    :: tmprows (:)
    Integer                              :: f2k     (4, 3) ! Maps from an element's face f and face's local node index kf to the element's local node index. This has been constructed solely for tetrahedra
    Real (KREAL)                         :: varsigma
    Real (KREAL)                         :: TWOA
    Real (KREAL)                         :: tmp1
    Real (KREAL)                         :: tmp2
    Real (KREAL)                         :: tmp3
    Real (KREAL)                         :: tmp4
    Real (KREAL)                         :: k0      (3)
    Real (KREAL)                         :: R0      (3)
    Real (KREAL)                         :: n       (3)
    Real (KREAL)                         :: kv      (3)
    Real (KREAL)                         :: wq      (3)
    Real (KREAL)                         :: Iek     (8)
    Real (KREAL)                         :: Ifk     (8)
    Real (KREAL),         Allocatable    :: tmpvals (:)
    Real (KREAL)                         :: cq      (2, 3)
    Real (KREAL)                         :: r       (3, 3)
    Real (KREAL)                         :: v       (3, 3)
    Real (KREAL)                         :: IM      (4, 4)
    
    ! A note about construction of the boundary angular fluence ---
    ! In the planar beam case (including slab, which is always planar), BoundaryFluence
    ! is truly a boundary fluence and contains no angular distribution.
    ! In the spherical beam case, it is an angular fluence.
    ! Neither case is normalized over the discrete ordinates set, so that is done separately.
    ! 
    ! In the planar beam case we store the angular distribution normalization (which is independent of position)
    ! in tmp4.
    ! Then, tmp3 contains the angular distribution evaluated at a given discrete ordinate, including the normalization factor
    ! tmp2 is a shape function
    ! tmp1 is the boundary fluence
    ! The result of the integrand is thus tmp1 * tmp2 * tmp3
    ! 
    ! In the spherical beam case, tmp4 is neglected
    ! tmp3 contains ONLY the angular distribution normalization, which is position-dependent
    ! tmp2 is a shape function
    ! tmp1 is the boundary ANGULAR fluence without normalization
    ! The result of the integrand is thus tmp1 * tmp2 * tmp3 again.
    
    !  ===============
    !    Assignments  
    !  ===============
    
    d    = SIZE(self%mesh%rg, dim=1)
    NE   = self%mesh%NE
    NEb  = COUNT(self%elemmask)
    NENK = self%mesh%NENK
    NI   = self%angular%NI
    
    f2k(1,:) = [1, 2, 3]
    f2k(2,:) = [1, 2, 4]
    f2k(3,:) = [1, 3, 4]
    f2k(4,:) = [2, 3, 4]
    
    !  ===================
    !    Pre-Allocations  
    !  ===================
    
    if (.not. self%mesh%slab) then
        call PrepTriangularQuadrature (cq, wq)
    end if
    call InitializeLegendreCoefficients (self%angular%L)
    
    ! If this has already been called, it is expected that the user
    ! changed something in mesh or angular and thus wants to rebuild the prep
    if (ALLOCATED(self%rows)) DEALLOCATE(self%rows)
    if (ALLOCATED(self%vals)) DEALLOCATE(self%vals)
    
    ! Going to proceed with the plan that I pre-allocate rows and vals to NENK, then truncate
    ! to the correct Nb as I visit angle-by-angle and element-by-element and construct the face contributions.
    ! I am doing this instead of NI * NENK for a balance between storage and performance.
    ALLOCATE(self%rows(0))
    ALLOCATE(self%vals(0))
    ALLOCATE(tmprows(NENK))
    ALLOCATE(tmpvals(NENK))
    
    ALLOCATE(einb, source=PACK([(e,e=1,NE)], mask=self%elemmask))
    
    !  ======================================
    !    Integrate and assign source values  
    !  ======================================
    
    if (streq(self%fldgeo%angdist, 'planar')) then
        spherical = .FALSE.
        k0   = self%fldgeo%axis
        tmp4 = self%AngularNorm ([ZERO, ZERO, ZERO]) ! Position doesn't matter
    else if (streq(self%fldgeo%angdist, 'spherical')) then
        R0 = self%fldgeo%origin
    end if
    
    b = 0
    
    ! FOR NOW --- Keeping slab and general case separate. There is definitely a way to merge them,
    ! but ultimately the point of typically keeping them merged is for diagnosis of issues with 
    ! general case. That said, the old code (which was merged) was validated for slab but not 
    ! general. So that's not always the best way to diagnose the general case.
    
    if (.not. self%mesh%slab) then
        
        do i = 1, NI
            
            kv(1:3) = self%angular%k(1:3,i)
            
            ! Make tmp3 by multiplying tmp4 by the actual angular distribution
            if (.not. spherical) then
                k0 = self%fldgeo%axis
                
                tmp3 = tmp4 * EvalBeamAngDist (kv, k0)
            end if
            
            tmprows = 0
            tmpvals = ZERO
            bloc    = 0
            do eb = 1, NEb
                e = einb(eb)
                
                NKe = self%mesh%el(e)%NK
                NFe = self%mesh%el(e)%NF
                idx = self%mesh%el(e)%idx
                
                ! Initialize the integrals for this element
                Iek = ZERO
                
                anyfaces = .FALSE.
                
                do f = 1, NFe
                    
                    !  -----------------------------------------
                    !    Check if this face is in the boundary
                    !  -----------------------------------------
                    
                    fg = self%mesh%el(e)%face(f)
                    
                    proceed = self%mesh%face(fg)%bdy
                    
                    if (.not. proceed) cycle
                    
                    !  ------------------------------------------------------
                    !    Check if this face is facing the discrete ordinate  
                    !  ------------------------------------------------------
                    
                    n(1:d) = self%mesh%el(e)%n(1:d,f)
                    
                    varsigma = DOT_PRODUCT(kv(1:d), n(1:d))
                    
                    proceed = varsigma < TOL
                    
                    if (.not. proceed) cycle
                    
                    ! Indicate that at least one face in this element is nontrivial
                    anyfaces = .TRUE.
                    
                    !  -------------------------
                    !    Evaluate the integral  
                    !  -------------------------
                    
                    TWOA = TWO * self%mesh%face(fg)%area ! This is technically the surface Jacobian vector's magnitude. Specific to triangle (face of tetrahedra)
                    NKf  = self%mesh%face(fg)%NK
                    NQ   = NKf
                    
                    r(1:3, 1:NKf) = self%mesh%rg(1:3, self%mesh%face(fg)%node(1:NKf))
                    v(1:3, 1:NQ)  = MapToTriangularQuad (cq, r) ! NOTE - for slab, the position you send to BoundaryFluence doesn't actually matter, at all.
                    
                    Ifk = ZERO
                    do kq = 1, NQ
                        tmp1 = BoundaryFluence (self%fldgeo, v(1:3, kq), kv)
                        
                        if (ISNAN(tmp1)) cycle ! Fluence returns NaN when the point is not in the field
                        
                        tmp1 = tmp1
                        
                        do kf = 1, NKf ! I could loop this outside of kq if I stored phi values. Consider doing that
                            k = f2k(f, kf)
                            
                            tmp2 = TriangularSF (kf, cq(1:2, kq))
                            
                            Ifk(k) = Ifk(k) + wq(kq) * tmp1 * tmp2
                        end do
                    end do
                    
                    Ifk(1:NKe) = Ifk(1:NKe) * TWOA
                    
                    Iek(1:NKe) = Iek(1:NKe) + varsigma * Ifk(1:NKe)
                    
                end do
                
                if (.not. anyfaces) cycle ! If no faces are on the boundary and facing the discrete ordinate, just cycle 
                
                !  ----------------------------------------------------------
                !    Mass-eliminate the integral to finally form the values  
                !  ----------------------------------------------------------
                
                IM = self%mesh%IMr(idx)%m / self%mesh%el(e)%detJv
                ! Don't forget the minus sign
                tmpvals(bloc + 1 : bloc + NKe) = tmpvals(bloc + 1 : bloc + NKe) - MATMUL(IM(1:NKe, 1:NKe), Iek(1:NKe))
                
                ! Assign the sadof to b
                do k = 1, NKe
                    tmprows(bloc + k) = self%mesh%offset(e) + k + (i - 1) * NENK
                end do
                
                ! Iterate the local (per-angle) row index
                bloc = bloc + NKe
                
            end do
            
            ! Iterate the total row index and update the row and vals arrays
            b = b + bloc
            
            self%rows = [self%rows, tmprows(1:bloc)]
            self%vals = [self%vals, tmpvals(1:bloc)]
            
        end do
        
    else
        
        k0(1) = COS(self%fldgeo%axis(1) * PI / 180.0_KREAL)
        
        do i = 1, NI
            
            kv(1) = self%angular%k(1,i)
            
            tmp3 = tmp4 * EvalBeamAngDist (kv(1), k0(1))
            ! tmp3 = ONE / self%angular%w(i)
            ! tmp3 = ONE
            
            tmprows = 0
            tmpvals = ZERO
            bloc    = 0
            do eb = 1, NEb
                e = einb(eb)
                
                NKe = self%mesh%el(e)%NK
                NFe = self%mesh%el(e)%NF
                idx = self%mesh%el(e)%idx
                
                ! ! Just use Iek as a dummy variable for the fluences
                ! Iek = ZERO
                ! do k = 1, NKe
                !     kg = self%mesh%el(e)%node(k)
                !     
                !     r(1,1) = self%mesh%rg(1, kg)
                !     
                !     Iek(k) = BoundaryFluence (self%fldgeo, r(1:3,1), kv, override=.TRUE.)
                ! end do
                
                anyfaces = .FALSE.
                
                do f = 1, NFe
                    
                    !  -----------------------------------------
                    !    Check if this face is in the boundary
                    !  -----------------------------------------
                    
                    fg = self%mesh%el(e)%face(f)
                    
                    proceed = self%mesh%face(fg)%bdy
                    
                    if (.not. proceed) cycle
                    
                    !  ------------------------------------------------------
                    !    Check if this face is facing the discrete ordinate  
                    !  ------------------------------------------------------
                    
                    n(1) = self%mesh%el(e)%n(1,f)
                    
                    varsigma = kv(1) * n(1)
                    
                    proceed = varsigma < TOL
                    
                    if (.not. proceed) cycle
                    
                    ! Indicate that at least one face in this element is nontrivial
                    anyfaces = .TRUE.
                    
                    tmpvals(bloc + 1 : bloc + NKe) = tmpvals(bloc + 1 : bloc + NKe) &
                                                   - varsigma * tmp3 * self%mesh%I(e)%IFM(f)%m(1:2,f) ! MATMUL(self%mesh%I(e)%IFM(f)%m, Iek(1:NKe))
                    
                end do
                
                if (.not. anyfaces) cycle
                
                ! Assign the sadof to b
                do k = 1, NKe
                    tmprows(bloc + k) = self%mesh%offset(e) + k + (i - 1) * NENK
                end do
                
                ! Iterate the local (per-angle) row index
                bloc = bloc + NKe
                
            end do
            
            ! Iterate the total row index and update the row and vals arrays
            b = b + bloc ! DON'T NEED THIS RIGHT?
            
            self%rows = [self%rows, tmprows(1:bloc)]
            self%vals = [self%vals, tmpvals(1:bloc)]
            
        end do
        
    end if
    
    ! MHY LATER - STILL FINDING SOME IDENTICALLY ZERO ENTRIES IN THE SOURCE. Not an accuracy problem but its not efficient
    ! Investigate proceed = varsigma < TOL, try flipping TOL to - TOL.
    
    ! do b = 1, SIZE(self%vals)
    !     print *, self%rows(b), self%vals(b)
    ! end do
    
    call DestroyLegendreCoefficients ()
    
End Subroutine
End Submodule