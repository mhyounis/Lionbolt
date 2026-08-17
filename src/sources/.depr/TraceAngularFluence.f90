Submodule (Sources) submodTraceAngularFluence
Contains
Module Subroutine TraceAngularFluence (self, attn, x)
    use constants
    use Quadrature
    use ShapeFunctions
    use Geometry
    use AngularSpace
    use LinearOperatorClass
    use Sources
    Implicit None
    Class (ExternalBeam),    Intent (InOut) :: self
    Real (KREAL),            Intent (In)    :: attn (:)
    
    Type (SpaceAngleVector), Intent (InOut) :: x
    
    Integer                                 :: e, eb, i, iRay, k, q, kq
    Logical                                 :: spherical
    Integer                                 :: NENK
    Integer                                 :: NI
    Integer                                 :: NR
    Integer                                 :: uNI
    Integer                                 :: ios
    Integer                                 :: NKe
    Integer                                 :: NQ
    Integer,                 Allocatable    :: einb   (:)
    Real (KREAL)                            :: tmp1
    Real (KREAL)                            :: tmp2
    Real (KREAL)                            :: tmp3
    Real (KREAL)                            :: tmp4
    Real (KREAL)                            :: OPL
    Real (KREAL)                            :: Beer
    Real (KREAL)                            :: k0  (3)
    Real (KREAL)                            :: R0  (3)
    Real (KREAL)                            :: kv  (3)
    Real (KREAL)                            :: wq  (4)
    Real (KREAL)                            :: Iek (8)
    Real (KREAL)                            :: cq  (3,4)
    Real (KREAL)                            :: r   (3,4)
    Real (KREAL)                            :: v   (3,4)
    
    ! A note about construction of the boundary angular fluence ---
    ! I hate to repeat code in different routines, especially such unintuitive code,
    ! but we use the same construction as in PrepSource.f90/PrepSource. See that if you're
    ! interested.
    
    if (.not. ALLOCATED(x%v)) then
        call stophere ('TraceAngularFluence.f90: TraceAngularFluence: The SpaceAngleVector ' // & 
                       'object sent to this routine must be pre-allocated.')
    else
        ! Initialize as zero
        x%v = ZERO
    end if
    
    call InitializeLegendreCoefficients (self%angular%L)
    
    ! Ray trace if needed
    if (.not. ALLOCATED(self%rays)) call self%RayTrace ()
    
    ! Assignments
    NENK = self%mesh%NENK
    NI   = self%angular%NI
    NR   = SIZE(self%rays)
    
    call PrepTetrahedralQuadrature (cq, wq)
    
    ALLOCATE(einb, source=PACK([(e,e=1,self%mesh%NE)], mask=self%elemmask))
    
    ! In the planar case, we will separately form the volume-integrated fluence and the angular distribution
    ! and just append it later
    if (streq(self%fldgeo%angdist, 'planar')) then
        spherical = .FALSE.
        k0 = self%fldgeo%axis
        tmp4 = self%AngularNorm ([ZERO, ZERO, ZERO]) ! Position doesn't matter
        
        uNI = 1 ! So that you just form one angular distribution
    else if (streq(self%fldgeo%angdist, 'spherical')) then
        R0 = self%fldgeo%origin
        
        uNI = NI
    end if
    
    ! Now, you visit each element in the beam and construct the integral
    do i = 1, uNI
        
        if (spherical) kv = self%angular%k(1:3,i)
        
        ! Angular offset
        ios = (i - 1) * NENK
        
        iRay = 0
        do eb = 1, SIZE(einb)
            e = einb(eb)
            
            NKe = self%mesh%el(e)%NK
            NQ  = NKe ! Temporary, works for tetrahedral elements only
            
            r(1:3, 1:4) = self%mesh%rg(1:3, self%mesh%el(e)%node(1:4))
            v = MapToTetrahedralQuad (cq, r)
            
            Iek = ZERO
            do kq = 1, 4
                if (.not. self%quadmask(self%mesh%offset(e) + kq)) cycle
                
                iRay = iRay + 1
                
                tmp1 = BoundaryFluence (self%fldgeo, v(1:3, kq), kv, override=.TRUE.)
                
                if (spherical) tmp3 = self%AngularNorm (v(1:3, kq))
                
                ! Beer's law factor --- exponential of negative optical path length
                ! Optical path length in my structure is given by summing pathlen over the materials, 
                ! weighted by the attenuation coefficient
                Beer = EXP(- DOT_PRODUCT(attn, self%rays(iRay)%pathlen))
                
                tmp1 = tmp1 * tmp3 * Beer
                
                do k = 1, NKe
                    tmp2 = TetrahedralSF (k, cq(1:3, kq))
                    
                    Iek(k) = Iek(k) + wq(kq) * tmp1 * tmp2
                end do
                
            end do
            
            ! Now, you matmul with I1inv, i.e., matmul with IMr and divide by detJv
            ! However, we don't use 1/detJv because we ought to have multiplied all Iek
            ! by detJv. Cancels out if we do neither.
            x%v(self%mesh%offset(e) + 1 + ios : self%mesh%offset(e) + NKe + ios) = MATMUL(self%mesh%IMr(4)%m, Iek(1:NKe))
        end do
    end do
    
    ! Construct the angular distribution for the planar case
    if (.not. spherical) then
        k0 = self%fldgeo%axis
        
        do i = 1, NI
            ! First make angular fluence isotropic (only 1 angle was populated)
            x%v(1 + (i - 1) * NENK : NENK + (i - 1) * NENK ) = x%v(1 : NENK)
            
            kv(1:3) = self%angular%k(1:3,i)
            x%v(1 + (i - 1) * NENK : NENK + (i - 1) * NENK ) = tmp4 * EvalBeamAngDist (kv, k0)
        end do
    end if
    
    ! Finally scale by number of particles in beam
    x%v = x%v * self%n
    
    call DestroyLegendreCoefficients ()
    
End Subroutine
End Submodule