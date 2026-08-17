Submodule (Sources) submodFCSSource
Contains
Module Subroutine FCSSource (self, OpK, s)
    use constants
    use Parallelism
    use IO
    use profiler
    use Quadrature
    use Polynomials
    use ShapeFunctions
    use LinearOperatorClass
    use ScatteringOperator
    use Sources
    Implicit None
    Class (ExternalBeam),    Intent (InOut) :: self
    Type (ScatteringOp),     Intent (In)    :: OpK
    
    Type (SpaceAngleVector), Intent (InOut) :: s
    
    Integer                                 :: e, eb, ell, i, idx, im2sd, k, kq, m, q, sdof
    Integer                                 :: NE
    Integer                                 :: NKe
    Integer                                 :: istart
    Integer,                 Allocatable    :: einb    (:)
    Real (KREAL)                            :: angfac
    Real (KREAL)                            :: kk0
    Real (KREAL),            Allocatable    :: phi0q   (:)
    Real (KREAL),            Allocatable    :: phi0ek  (:)
    Real (KREAL)                            :: cqs     (2)
    Real (KREAL)                            :: wqs     (2)
    Real (KREAL)                            :: kv      (3)
    Real (KREAL)                            :: k0      (3)
    Real (KREAL)                            :: R0      (3)
    Real (KREAL)                            :: r       (3)
    Real (KREAL)                            :: wq      (4)
    Real (KREAL)                            :: tmp8    (8)
    Real (KREAL)                            :: re      (3,4)
    Real (KREAL)                            :: cq      (3,4)
    Real (KREAL)                            :: v       (3,4)
    Real (KREAL),            Allocatable    :: u       (:,:)
    
    if (.not. ALLOCATED(s%v)) then
        s = SpaceAngleVector (self%mesh, self%angular)
    end if
    if (.not. OpK%PHYSREADY) then
        call stophere ('FCSSource.f90: FCSSource: K operator provided has not had its physics built. You must build this.')
    end if
    phi0q = self%FluenceQuadrature (OpK%XS%t) ! This is important - it uses the attenuation assigned to this group. So you cannot do this FCS Source stuff for off-diagonal entries. (You don't ever need to right?)
                                              ! Perhaps delta down will have something to say about this?
    call InitializeLegendreCoefficients (self%angular%L)
    
    NE = self%mesh%NE
    
    if (self%mesh%slab) then
        
        ALLOCATE(phi0ek(self%mesh%NENK), source=ZERO)
        ALLOCATE(u(2,2))
        
        call PrepLinearQuadrature (cqs, wqs)
        
        do kq = 1, 2
            do k = 1, 2
                u(k,kq) = LineSF (k, cqs(kq))
            end do
        end do
        
        ! Get the fluence for elements and nodes by integrating with shape functions
        ! using the quadrature node info.
        !$OMP PARALLEL DO DEFAULT (SHARED) &
        !$OMP PRIVATE (e, NKe, k, sdof, kq, q, idx, istart)
        do e = 1, NE
            NKe = self%mesh%el(e)%NK
            
            do k = 1, NKe
                sdof = self%mesh%offset(e) + k
                
                phi0ek(sdof) = ZERO
                
                do kq = 1, NKe
                    q = self%mesh%offset(e) + kq
                    if (.not. self%quadmask(q)) cycle
                    
                    phi0ek(sdof) = phi0ek(sdof) + wqs(kq) * u(k,kq) * phi0q(q)
                end do
                
            end do
            
            ! Now, you matmul with I1inv, i.e., matmul with IMr and divide by detJv
            ! However, we don't use 1/detJv because we ought to have multiplied all Iek
            ! by detJv. Cancels out if we do neither.
            idx    = self%mesh%el(e)%idx
            istart = self%mesh%offset(e)
            phi0ek(istart + 1 : istart + NKe) = &
                MATMUL(self%mesh%IMr(idx)%m, phi0ek(istart + 1 : istart + NKe))
            
        end do
        !$OMP END PARALLEL DO
        
        ! Now append the angular distribution
        k0(1) = COS(self%fldgeo%axis(1) * PI / 180.0_KREAL)
        !$OMP PARALLEL DO DEFAULT (SHARED) &
        !$OMP PRIVATE (i, kv, ell, angfac, m, im2sd, sdof)
        do i = 1, self%angular%NI
            kv(1) = self%angular%k(1,i)
            
            if (self%mesh%Nm == 1) then
                angfac = ZERO
                do ell = 0, self%angular%L
                    angfac = angfac &
                           + HALF * (2 * ell + 1) * SimpleHorner (x=kv(1), c=coeff(ell,0:ell)) &
                           * SimpleHorner (x=k0(1), c=coeff(ell,0:ell)) * OpK%XS%s(ell, 1)
                end do
                s%v(1 + (i - 1) * self%mesh%NENK : i * self%mesh%NENK) = &
                    s%v(1 + (i - 1) * self%mesh%NENK : i * self%mesh%NENK) + phi0ek * angfac
            else
                do m = 1, self%mesh%Nm
                    angfac = ZERO
                    do ell = 0, self%angular%L
                        angfac = angfac &
                               + HALF * (2 * ell + 1) * SimpleHorner (x=kv(1), c=coeff(ell,0:ell)) &
                               * SimpleHorner (x=k0(1), c=coeff(ell,0:ell)) * OpK%XS%s(ell, m)
                    end do
                    
                    do im2sd = 1, SIZE(self%mesh%mat2sd(m)%v)
                        sdof = self%mesh%mat2sd(m)%v(im2sd)
                        s%v(sdof + (i - 1) * self%mesh%NENK) = &
                            s%v(sdof + (i - 1) * self%mesh%NENK) + phi0ek(sdof) * angfac
                    end do
                end do
            end if
        end do
        !$OMP END PARALLEL DO
        
    else if (streq(self%fldgeo%angdist, 'planar')) then
        
        ALLOCATE(einb, source=PACK([(e,e=1,NE)], mask=self%elemmask))
        ALLOCATE(u(4,4))
        
        ALLOCATE(phi0ek(self%mesh%NENK), source=ZERO)
        
        call PrepTetrahedralQuadrature (cq, wq)
        
        ! ----- HERE IS WHERE THE TETRAHEDRAL CASE IS ENFORCED -----
        do kq = 1, 4
            do k = 1, 4
                u(k,kq) = TetrahedralSF (k, cq(:,kq))
            end do
        end do
        
        ! Get the fluence for elements and nodes by integrating with shape functions
        ! using the quadrature node info.
        !$OMP PARALLEL DO DEFAULT (SHARED) &
        !$OMP PRIVATE (eb, e, NKe, k, sdof, kq, q, idx, istart)
        do eb = 1, SIZE(einb)
            e   = einb(eb)
            NKe = self%mesh%el(e)%NK
            
            do k = 1, NKe
                sdof = self%mesh%offset(e) + k
                
                do kq = 1, 4
                    q = self%mesh%offset(e) + kq
                    if (.not. self%quadmask(q)) cycle
                    
                    ! note, phi0ek was already initialized as zero
                    phi0ek(sdof) = phi0ek(sdof) + wq(kq) * u(k,kq) * phi0q(q)
                end do
                
            end do
            
            idx    = self%mesh%el(e)%idx
            istart = self%mesh%offset(e)
            phi0ek(istart + 1 : istart + NKe) = &
                MATMUL(self%mesh%IMr(idx)%m, phi0ek(istart + 1 : istart + NKe))
            
        end do
        !$OMP END PARALLEL DO
        
        ! Now append the angular distribution
        k0 = self%fldgeo%axis
        ! MHY LATER - just a code cleanup thing, I should rethink the defaults
        !$OMP PARALLEL DO DEFAULT (SHARED) &
        !$OMP PRIVATE (i, kv, kk0, ell, angfac, m, im2sd, sdof)
        do i = 1, self%angular%NI
            kv  = self%angular%k(1:3,i)
            kk0 = DOT_PRODUCT(kv, k0)
            
            if (self%mesh%Nm == 1) then
                angfac = ZERO
                do ell = 0, self%angular%L
                    angfac = angfac &
                           + (2 * ell + 1) * SimpleHorner (x=kk0, c=coeff(ell,0:ell)) * OpK%XS%s(ell, 1) / FOURPI
                    if (ISNAN(angfac)) then
                        print *, i, ell
                    end if
                end do
                s%v(1 + (i - 1) * self%mesh%NENK : i * self%mesh%NENK) = &
                    s%v(1 + (i - 1) * self%mesh%NENK : i * self%mesh%NENK) + phi0ek * angfac
            else
                do m = 1, self%mesh%Nm
                    angfac = ZERO
                    do ell = 0, self%angular%L
                        angfac = angfac &
                           + (2 * ell + 1) * SimpleHorner (x=kk0, c=coeff(ell,0:ell)) * OpK%XS%s(ell, m) / FOURPI
                    end do
                    
                    do im2sd = 1, SIZE(self%mesh%mat2sd(m)%v)
                        sdof = self%mesh%mat2sd(m)%v(im2sd)
                        s%v(sdof + (i - 1) * self%mesh%NENK) = &
                            s%v(sdof + (i - 1) * self%mesh%NENK) + phi0ek(sdof) * angfac
                    end do
                end do
            end if
        end do
        !$OMP END PARALLEL DO
        
    else if (streq(self%fldgeo%angdist, 'spherical')) then
        ! Spherical case is a little trickier because the mesh nodes determine k0.
        ! You do indeed need to account for this in the spatial integration
        ! call profile ('FCSSource')
        ALLOCATE(einb, source=PACK([(e,e=1,NE)], mask=self%elemmask))
        ALLOCATE(u(4,4))
        
        R0 = self%fldgeo%origin
        
        call PrepTetrahedralQuadrature (cq, wq)
        
        ! ----- HERE IS WHERE THE TETRAHEDRAL CASE IS ENFORCED -----
        do kq = 1, 4
            do k = 1, 4
                u(k,kq) = TetrahedralSF (k, cq(:,kq))
            end do
        end do
        
        !$OMP PARALLEL DO DEFAULT (SHARED) &
        !$OMP PRIVATE (i, eb, e, NKe, m, re, v, k, sdof, tmp8, kq, q, k0, kv, kk0) &
        !$OMP PRIVATE (angfac, ell, idx, istart)
        do i = 1, self%angular%NI
            
            do eb = 1, SIZE(einb)
                e   = einb(eb)
                NKe = self%mesh%el(e)%NK
                m   = self%mesh%el(e)%mat
                
                re(1:3, 1:4) = self%mesh%rg(1:3, self%mesh%el(e)%node(1:4)) ! Here too
                
                v = MapToTetrahedralQuad (cq, re)
                
                do k = 1, NKe
                    sdof = self%mesh%offset(e) + k
                    
                    tmp8(k) = ZERO
                    
                    do kq = 1, 4 ! ----- HERE IS WHERE THE TETRAHEDRAL CASE IS ENFORCED -----
                        q = self%mesh%offset(e) + kq ! Here too
                        if (.not. self%quadmask(q)) cycle
                        
                        k0  = (v(:,kq) - R0) / NORM2(v(:,kq) - R0)
                        kv  = self%angular%k(1:3,i)
                        kk0 = DOT_PRODUCT(k0, kv)
                        angfac = ZERO
                        do ell = 0, self%angular%L
                            ! MHY LATER - If I refactor the delta down stuff I may have to revisit this
                            angfac = angfac &
                                   + (2 * ell + 1) * SimpleHorner (x=kk0, c=coeff(ell,0:ell)) * OpK%XS%s(ell, m) / FOURPI
                        end do
                        
                        tmp8(k) = tmp8(k) + wq(kq) * u(k,kq) * phi0q(q) * angfac ! angfac here includes scattering XS, unlike in AddAngularFluence
                    end do
                    
                end do
                
                idx = self%mesh%el(e)%idx
                tmp8(1:NKe) = MATMUL(self%mesh%IMr(idx)%m, tmp8(1:NKe))
                
                ! Now contribute to s for this element and angle
                istart = self%mesh%offset(e) + (i - 1) * self%mesh%NENK 
                s%v(istart + 1 : istart + NKe) = s%v(istart + 1 : istart + NKe) + tmp8(1:NKe)
            end do
            
        end do
        !$OMP END PARALLEL DO
        ! call profile ('FCSSource')
    end if
    
    call DestroyLegendreCoefficients ()
    
End Subroutine
End Submodule