Program HetScatterMMS
    use constants
    use HDF5Tree
    use NittanyAPI
    use LionboltAPI
    Implicit None
    
    Type (MeshClass),    Target :: mesh
    Type (AngularClass), Target :: angular
    Type (XSType),       Target :: XS
    Type (ScatteringOp)         :: OpK
    
    Type (SpaceAngleVector)     :: s
    Type (SpaceAngleVector)     :: strue
    
    Integer                     :: e, ell, i, k, sdof
    Integer                     :: L
    Real (KREAL)                :: r  (3)
    Real (KREAL)                :: kv (3)
    Real (KREAL)                :: tmp
    Real (KREAL)                :: S2S1 = ONE
    
    mesh = MeshClass ('/home/myounis/work/meshes/mat1-mat2-mat1_pencil-beam.msh')
    call mesh%Translate([- HALF, - HALF, - HALF])
    call mesh%Scale([TWO, TWO, TWO])
    call mesh%PostProcess ()
    
    L = 15
    angular = AngularClass (L=L, solver='GMRES', slab=.FALSE.)
    
    ALLOCATE(XS%s(0:L,2))
    
    do ell = 0, L
        XS%s(ell,1) = ScatteringMoment (ell)
        XS%s(ell,2) = S2S1 * XS%s(ell,1)
    end do
    
    call OpK%Build (mesh, angular, XS)
    
    s     = SpaceAngleVector (mesh, angular)
    strue = SpaceAngleVector (mesh, angular)
    
    s%v = ZERO
    do i = 1, angular%NI
        kv = angular%k(:,i)
        do e = 1, mesh%NE
            do k = 1, 4
                sdof = mesh%offset(e) + k
                r    = mesh%rg(1:3, mesh%connectivity(sdof))
                
                tmp  = NORM2(r + kv)
                
                s%v(sdof + (i - 1) * mesh%NENK) = EXP(DOT_PRODUCT(r,kv))
                if (tmp == ZERO) then
                    strue%v(sdof + (i - 1) * mesh%NENK) = FOURPI
                else
                    strue%v(sdof + (i - 1) * mesh%NENK) = TWOPI * (EXP(tmp) - EXP(-tmp)) / tmp
                end if
                
                if (mesh%el(e)%mat == 2) then
                    strue%v(sdof + (i - 1) * mesh%NENK) = strue%v(sdof + (i - 1) * mesh%NENK) * S2S1
                end if
                
            end do
        end do
    end do
    
    call OpK%MatVec (s)
    
    ! CONSIDER LOOKING AT FLUENCES TOO
    
    ! MHY LATER - need to make this more user friendly, and all accessible without core library
    call OpenHDF5             ('HetScatterMMS.h5')
    call WriteProblemTypeH5   ('general')
    call WriteMeshH5          (mesh)
    call OpenParticleH5       ('virtual')
    call AppendOrdinatesH5    (angular%w, angular%k)
    call AppendAngularFluence (s%v, mesh%NENK, .FALSE., 1) ! Super ugly and un-friendly
    call CloseParticleH5      ('virtual')
    call OpenParticleH5       ('exact')
    call AppendOrdinatesH5    (angular%w, angular%k)
    call AppendAngularFluence (strue%v, mesh%NENK, .FALSE., 1) ! Super ugly and un-friendly
    call CloseParticleH5      ('exact')
    call CloseHDF5            ()
    
    print *, 'HetScatterMMS.f90 - Successfully terminated.'
    
    Contains
    
    Function ScatteringMoment (n) Result (f)
        Implicit None
        Integer, Intent (In) :: n
        
        Real (KREAL)         :: f
        
        Integer              :: k, m
        Real (KREAL)         :: tmp
        Real (KREAL)         :: S1
        Real (KREAL)         :: S2
        Real (KQUAD)         :: fQ
        Real (KQUAD)         :: S1Q
        Real (KQUAD)         :: S2Q
        
        if (n >= 14) then
            f = ZERO
            return
        end if
        
        S1  = ZERO
        S2  = ZERO
        tmp = ONE
        do m = 0, n
            S1 = S1 + tmp * (-1)**m
            S2 = S2 + tmp
            
            if (m == n) exit
            
            tmp = tmp * (n + m + 1) * (n - m) / (2 * (m + 1))
        end do
        
        S1Q = REAL(S1, KIND=KQUAD)
        S2Q = REAL(S2, KIND=KQUAD)
        fQ  = S1Q * EXP(REAL(ONE, KIND=KQUAD)) - S2Q * (-1)**n * EXP(- REAL(ONE, KIND=KQUAD))
        
        f = TWOPI * REAL(fQ, KIND=KREAL)
        
    End Function
    
End Program