Program TransportMMS
    use constants
    use HDF5Tree
    use NittanyAPI
    use LionboltAPI
    Implicit None
    
    !                               ####################
    !                                   TransportMMS    
    !                                    Written by     
    !                                 Muhsin H. Younis  
    !                               ####################
    
    !  ==================================================================================
    !    This driver file performs the Method of Manufactured Solutions (MMS)
    !    to validate the inversion of the transport operator in Lionbolt.
    !    Basically, given any manually defined source s that satisfies vacuum boundary 
    !    conditions, it will compute T^{-1}s, where T is the transport operator, 
    !    -\mathbf{\hat{k}}cdot\nabla + \sigma for some constant attenuation coefficient 
    !    \sigma.
    !    
    !    It will output a single angular fluence array to a particle named 'virtual.'
    !    
    !    The below should demonstrate the simplicity and user-friendliness of the
    !    Lionbolt API. Besides the declarations and the definition of the source, all
    !    it takes is ~10 lines of code to carry out this rather complex operation from
    !    scratch.
    !  ==================================================================================
    
    Type (MeshClass),    Target :: mesh
    Type (AngularClass), Target :: angular
    Type (XSType),       Target :: XS
    Type (TransportOp)          :: OpT
    
    Type (SpaceAngleVector)     :: s
    Type (SpaceAngleVector)     :: strue
    
    Integer                     :: e, i, k, sdof
    Real (KREAL)                :: x
    Real (KREAL)                :: y
    Real (KREAL)                :: z
    Real (KREAL)                :: r  (3)
    Real (KREAL)                :: kv (3)
    Real (KREAL), Parameter     :: A     = TEN
    Real (KREAL), Parameter     :: sigma = 1.0e-10_KREAL
    
    mesh    = MeshClass ('simple_cube.msh')
    call mesh%Scale([TWO, TWO, TWO])
    call mesh%PostProcess ()
    angular = AngularClass (15, 'GMRES', .FALSE.)
    
    ALLOCATE(XS%t(1))
    XS%t = sigma
    
    call OpT%Build (mesh, angular, XS)
    
    s     = SpaceAngleVector (mesh, angular)
    strue = SpaceAngleVector (mesh, angular)
    
    s%v     = ZERO
    strue%v = ZERO
    do i = 1, angular%NI
        kv = angular%k(:,i)
        do e = 1, mesh%NE
            do k = 1, 4
                sdof = mesh%offset(e) + k
                r    = mesh%rg(1:3, mesh%connectivity(sdof))
                
                x = r(1)
                y = r(2)
                z = r(3)
                
                ! Spatial Gaussian solution, uniform in angle
                s%v(sdof + (i - 1) * mesh%NENK) = &
                    (- TWO * A * DOT_PRODUCT(kv, r) + XS%t(1)) * EXP(- A * DOT_PRODUCT(r, r))
                
                strue%v(sdof + (i - 1) * mesh%NENK) = EXP(- A * DOT_PRODUCT(r, r))
                
                ! ! Spatial Guassian with angular-dependent amplitude.
                ! s%v(sdof + (i - 1) * mesh%NENK) = &
                !     (- TWO * A * DOT_PRODUCT(kv, r)**2 + ONE + XS%t(1) * DOT_PRODUCT(kv, r)) &
                !      * EXP(- A * DOT_PRODUCT(r, r))
                ! 
                ! strue%v(sdof + (i - 1) * mesh%NENK) = DOT_PRODUCT(kv, r) * EXP(- A * DOT_PRODUCT(r, r))
            end do
        end do
    end do
    
    call OpT%MatInv (s)
    
    ! MHY LATER - need to make this more user friendly, and all accessible without core library
    call OpenHDF5             ('TransportMMS.h5')
    call WriteProblemTypeH5   ('general')
    call WriteMeshH5          (mesh)
    call OpenParticleH5       ('virtual')
    call AppendOrdinatesH5    (angular%w, angular%k)
    call AppendAngularFluence (s%v, mesh%NENK, .FALSE., 1) ! Super ugly and un-friendly
    call CloseParticleH5      ('virtual')
    call OpenParticleH5       ('exact')
    call AppendOrdinatesH5    (angular%w, angular%k)
    call AppendAngularFluence (strue%v, mesh%NENK, .FALSE., 1)
    call CloseParticleH5      ('exact')
    call CloseHDF5            ()
    
    print *, 'TransportMMS.f90 - Successfully terminated.'
    
End Program