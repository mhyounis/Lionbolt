Submodule (Sources) submodAddAngularFluence
Contains
Module Subroutine AddAngularFluence (self, attn, unc)
    ! Stores the spatial profile of the uncollided fluence.
    ! This can be re-used if energy iteration is being performed,
    ! preventing the neeed to continually reconstruct these values.
    use constants
    use Parallelism
    use IO
    use profiler
    use Quadrature
    use ShapeFunctions
    use Sources
    use LinearOperatorClass
    Implicit None
    Class (ExternalBeam),    Intent (InOut) :: self
    Real (KREAL),            Intent (In)    :: attn (:)
    
    Type (SpaceAngleVector), Intent (InOut) :: unc
    
    Integer                                 :: e, eb, i, idx, k, kq, q, sdof
    Integer                                 :: NE
    Integer                                 :: NKe
    Integer                                 :: istart
    Integer,                 Allocatable    :: einb   (:)
    Real (KREAL)                            :: angfac
    Real (KREAL),            Allocatable    :: phi0q  (:)
    Real (KREAL),            Allocatable    :: phi0ek (:)
    Real (KREAL)                            :: cqs    (2)
    Real (KREAL)                            :: wqs    (2)
    Real (KREAL)                            :: kv     (3)
    Real (KREAL)                            :: k0     (3)
    Real (KREAL)                            :: R0     (3)
    Real (KREAL)                            :: r      (3)
    Real (KREAL)                            :: wq     (4)
    Real (KREAL)                            :: tmp8   (8)
    Real (KREAL)                            :: re     (3,4)
    Real (KREAL)                            :: cq     (3,4)
    Real (KREAL)                            :: v      (3,4)
    Real (KREAL),            Allocatable    :: u      (:,:)
    
    if (.not. ALLOCATED(unc%v)) then
        unc = SpaceAngleVector (self%mesh, self%angular)
    end if
    phi0q = self%FluenceQuadrature (attn)
    
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
        !$OMP PRIVATE (i, kv, angfac)
        do i = 1, self%angular%NI
            kv(1) = self%angular%k(1,i)
            
            angfac = SlabAngularDistribution (kv(1), k0(1))
            unc%v(1 + (i - 1) * self%mesh%NENK : i * self%mesh%NENK) = &
                unc%v(1 + (i - 1) * self%mesh%NENK : i * self%mesh%NENK) + phi0ek * angfac
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
                
                phi0ek(sdof) = ZERO
                
                do kq = 1, NKe
                    q = self%mesh%offset(e) + kq
                    if (.not. self%quadmask(q)) cycle
                    
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
        !$OMP PARALLEL DO DEFAULT (SHARED) &
        !$OMP PRIVATE (i, kv, angfac)
        do i = 1, self%angular%NI
            kv = self%angular%k(1:3,i)
            
            angfac = GenAngularDistribution (kv, k0)
            unc%v(1 + (i - 1) * self%mesh%NENK : i * self%mesh%NENK) = &
                unc%v(1 + (i - 1) * self%mesh%NENK : i * self%mesh%NENK) + phi0ek * angfac
        end do
        !$OMP END PARALLEL DO
        
    else if (streq(self%fldgeo%angdist, 'spherical')) then
        ! Spherical case is a little trickier because the mesh nodes determine k0.
        ! You do indeed need to account for this in the spatial integration
        ! call profile ('AddAngularFluence')
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
        !$OMP PRIVATE (i, eb, e, NKe, re, v, k, sdof, tmp8, kq, q, k0, kv, angfac) &
        !$OMP PRIVATE (idx, istart)
        do i = 1, self%angular%NI
            
            do eb = 1, SIZE(einb)
                e   = einb(eb)
                NKe = self%mesh%el(e)%NK
                
                re(1:3, 1:4) = self%mesh%rg(1:3, self%mesh%el(e)%node(1:4)) ! Here too
                
                v = MapToTetrahedralQuad (cq, re)
                
                do k = 1, NKe
                    sdof = self%mesh%offset(e) + k
                    
                    tmp8(k) = ZERO
                    
                    do kq = 1, NKe ! ----- HERE IS WHERE THE TETRAHEDRAL CASE IS ENFORCED -----
                        q = self%mesh%offset(e) + kq ! Here too
                        if (.not. self%quadmask(q)) cycle
                        
                        k0     = (v(:,kq) - R0) / NORM2(v(:,kq) - R0)
                        kv     = self%angular%k(1:3,i)
                        angfac = GenAngularDistribution (kv, k0)
                        
                        tmp8(k) = tmp8(k) + wq(kq) * u(k,kq) * phi0q(q) * angfac
                    end do
                    
                end do
                
                idx = self%mesh%el(e)%idx
                tmp8(1:NKe) = MATMUL(self%mesh%IMr(idx)%m, tmp8(1:NKe))
                
                ! Now contribute to unc for this element and angle
                istart = self%mesh%offset(e) + (i - 1) * self%mesh%NENK 
                unc%v(istart + 1 : istart + NKe) = unc%v(istart + 1 : istart + NKe) + tmp8(1:NKe)
            end do
            
        end do
        !$OMP END PARALLEL DO
        ! call profile ('AddAngularFluence')
        
    end if
    
    call DestroyLegendreCoefficients ()
    
End Subroutine
End Submodule