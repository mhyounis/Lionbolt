Module FEAInnerProducts
    use constants
    use types
    use BasicMathFunctions
    use Quadrature
    use ShapeFunctions
    Implicit None

! Uses N = 4 Gauss-Legendre quadrature in each dimension
Integer,                   Private :: Nq
Real (KREAL), Allocatable, Private :: xq (:)
Real (KREAL), Allocatable, Private :: wq (:)

Contains

Subroutine InitializeFEAQuadrature (N)
    Implicit None
    Integer, Intent (In) :: N
    
    Nq = N
    
    call PrepGaussLegendreQuadrature (Nq, xq, wq)
    
End Subroutine

Subroutine DestroyFEAQuadrature ()
    Implicit None
    
    DEALLOCATE(xq)
    DEALLOCATE(wq)
    
End Subroutine

Subroutine FEAMassRef (idx, I)
    Implicit None
    !  =====================================================
    !    Computes the mass integral in a reference element
    !  =====================================================
    Integer,                   Intent (In)  :: idx ! Element type index (according to GMSH)
    
    Real (KREAL), Allocatable, Intent (Out) :: I (:,:)
    
    Integer                                 :: ii, j, k, kk, kp
    Integer                                 :: NK
    Real (KREAL)                            :: x
    Real (KREAL)                            :: y
    Real (KREAL)                            :: z
    Real (KREAL)                            :: uk
    Real (KREAL)                            :: ukp
    Real (KREAL)                            :: r (3)
    
    ! COULD USE SYMMETRY OF RESULTING INTEGRAL. THESE ONLY NEED TO BE DONE ONCE THOUGH SO PROBABLY NOT WORTH IT
    
    select case (idx)
    case (1)
        NK = 2
        ALLOCATE(I(NK, NK))
        
        I = ZERO
        
        do ii = 1, Nq
            z = xq(ii)
            do kp = 1, NK
                ukp = LineSF(kp, z)
                do k = 1, NK
                    uk = LineSF(k, z)
                    
                    I(k,kp) = I(k,kp) + wq(ii) * uk * ukp
                end do
            end do
        end do
        
    case (4)
        NK = 4
        ALLOCATE(I(NK, NK))
        
        I = ZERO
        
        do ii = 1, Nq
            x = HALF * xq(ii) + HALF
            do j = 1, Nq
                y = HALF * (ONE - x) * xq(j) + HALF * (ONE - x)
                do kk = 1, Nq
                    z = HALF * (ONE - x - y) * xq(kk) + HALF * (ONE - x - y)
                    
                    r = [x, y, z]
                    do kp = 1, NK
                        ukp = TetrahedralSF(kp, r)
                        do k = 1, NK
                            uk = TetrahedralSF(k, r)
                            
                            I(k,kp) = I(k,kp) + wq(ii) * wq(j) * wq(kk) * EIGTH * (1 - x) * (1 - x - y) * uk * ukp
                        end do
                    end do
                    
                end do
            end do
        end do
        
    case (5)
        NK = 8
        ALLOCATE(I(NK, NK))
        
        I = ZERO

        do ii = 1, Nq
            x = xq(ii)
            do j = 1, Nq
                y = xq(j)
                do kk = 1, Nq
                    z = xq(kk)
                    
                    r = [x, y, z]
                    do kp = 1, NK
                        ukp = HexahedralSF(kp, r)
                        do k = 1, NK
                            uk = HexahedralSF(k, r)
                            
                            I(k,kp) = I(k,kp) + wq(ii) * wq(j) * wq(kk) * uk * ukp
                        end do
                    end do
                    
                end do
            end do
        end do
        
    case (6)
        NK = 6
        print *, 'wip12124'
        stop
    case (7)
        NK = 5
        print *, 'wip12547'
        stop
    end select
    
End Subroutine

Subroutine FEAConvectionRef (idx, I)
    Implicit None
    !  ============================================================
    !    Computes the convection integrals in a reference element
    !  ============================================================
    Integer,                   Intent (In)  :: idx ! Element type index (according to GMSH)
    
    Real (KREAL), Allocatable, Intent (Out) :: I (:,:,:)
    
    Integer                                 :: ii, j, k, kk, kp
    Integer                                 :: NK
    Real (KREAL)                            :: x
    Real (KREAL)                            :: y
    Real (KREAL)                            :: z
    Real (KREAL)                            :: ukp
    Real (KREAL)                            :: Duk  (3)
    Real (KREAL)                            :: r    (3)
    
    select case (idx)
    case (1)
        NK = 2
        ALLOCATE(I(NK, NK, 1))
        
        I = ZERO
        
        do ii = 1, Nq
            z = xq(ii)
            do kp = 1, NK
                ukp = LineSF(kp, z)
                do k = 1, NK
                    Duk(1) = GradLineSF(k, z)
                    
                    I(k,kp,1) = I(k,kp,1) + wq(ii) * ukp * Duk(1)
                end do
            end do
        end do
    
    case (4)
        NK = 4
        ALLOCATE(I(NK, NK, 3))
        
        I = ZERO

        do ii = 1, Nq
            x = HALF * xq(ii) + HALF
            do j = 1, Nq
                y = HALF * (ONE - x) * xq(j) + HALF * (ONE - x)
                do kk = 1, Nq
                    z = HALF * (ONE - x - y) * xq(kk) + HALF * (ONE - x - y)
                    
                    r = [x, y, z]
                    do kp = 1, NK
                        ukp = TetrahedralSF(kp, r)
                        do k = 1, NK
                            Duk = GradTetrahedralSF(k, r)
                            
                            I(k,kp,:) = I(k,kp,:) + wq(ii) * wq(j) * wq(kk) * EIGTH * (ONE - x) * (ONE - x - y) * ukp * Duk
                        end do
                    end do
                end do
            end do
        end do
        
    case (5)
        NK = 8
        ALLOCATE(I(NK, NK, 3))
        
        I = ZERO
        
        do ii = 1, Nq
            x = xq(ii)
            do j = 1, Nq
                y = xq(j)
                do kk = 1, Nq
                    z = xq(kk)
                    
                    r = [x, y, z]
                    do kp = 1, NK
                        ukp = HexahedralSF(kp, r)
                        do k = 1, NK
                            Duk = GradHexahedralSF(k, r)
                            
                            I(k,kp,:) = I(k,kp,:) + wq(ii) * wq(j) * wq(kk) * ukp * Duk
                        end do
                    end do
                    
                end do
            end do
        end do
        
    case (6)
        NK = 6
        print *, 'wip12124111'
        stop
    case (7)
        NK = 5
        print *, 'wip12547111'
        stop
    end select
    
End Subroutine

Subroutine FEAFaceMassRef (NK, I)
    Implicit None
    !  ===============================================================
    !    Computes the mass integral on a face in a reference element
    !  ===============================================================
    Integer,                   Intent (In)  :: NK ! Number of nodes on a face
    
    Real (KREAL), Allocatable, Intent (Out) :: I (:,:)
    
    Integer                                 :: ii, j, k, kp
    Real (KREAL)                            :: x
    Real (KREAL)                            :: y
    Real (KREAL)                            :: uk
    Real (KREAL)                            :: ukp
    Real (KREAL)                            :: r (2)
    
    select case (NK)
    case (1)
        ALLOCATE(I(NK, NK))
        
        I(1, 1) = ONE
        
    case (3)
        ALLOCATE(I(NK, NK))
        
        I = ZERO
        
        do ii = 1, Nq
            x = HALF * xq(ii) + HALF
            do j = 1, Nq
                y = HALF * (ONE - x) * xq(j) + HALF * (ONE - x)
                
                r = [x, y]
                do kp = 1, NK
                    ukp = TriangularSF(kp, r)
                    do k = 1, NK
                        uk = TriangularSF(k, r)
                        
                        I(k,kp) = I(k,kp) + wq(ii) * wq(j) * FOURTH * (ONE - x) * uk * ukp
                    end do
                end do
            end do
        end do
        
    case (4)
        ALLOCATE(I(NK, NK))
        
        I = ZERO
        
        do ii = 1, Nq
            x = xq(ii)
            do j = 1, Nq
                y = xq(j)
                
                r = [x, y]
                do kp = 1, NK
                    ukp = RectangularSF(kp, r)
                    do k = 1, NK
                        uk = RectangularSF(k, r)
                        
                        I(k,kp) = I(k,kp) + wq(ii) * wq(j) * uk * ukp
                    end do
                end do
            end do
        end do
        
    end select
    
End Subroutine

Subroutine MEFEAConvection (IM, ICr, JvinvT, detJv, IC)
    Implicit None
    ! Determine mass-eliminated convective inner product. Array must be pre-allocated
    Real (KREAL), Intent (In)  :: IM     (:,:)
    Real (KREAL), Intent (In)  :: ICr    (:,:,:)
    Real (KREAL), Intent (In)  :: JvinvT (3,3)
    Real (KREAL), Intent (In)  :: detJv
    
    Real (KREAL), Intent (Out) :: IC     (:,:,:)
    
    Integer                    :: dir, k, kp
    Integer                    :: d
    Integer                    :: NK
    
    d  = SIZE(ICr, DIM=3)
    NK = SIZE(ICr, DIM=1)
    
    do kp = 1, NK
        do k = 1, NK
            IC(k,kp,1:d) = MATMUL(JvinvT(1:d,1:d), ICr(k,kp,1:d))
        end do
    end do
    
    do dir = 1, d
        IC(1:NK,1:NK,dir) = detJv * MATMUL(IM, IC(1:NK,1:NK,dir))
    end do
    
End Subroutine

Subroutine MEFEAFaceMass (IM, f2k, IFMr, normJs, IFM)
    Implicit None
    ! Determine mass-eliminated face mass inner product. Array must be pre-allocated
    Real (KREAL), Intent (In)  :: IM     (:,:)
    Integer,      Intent (In)  :: f2k    (:)
    Real (KREAL), Intent (In)  :: IFMr   (:,:)
    Real (KREAL), Intent (In)  :: normJs
    
    Real (KREAL), Intent (Out) :: IFM    (:,:)
    
    Integer                    :: k, k1, k2, kp
    Integer                    :: NK
    Integer                    :: NKf
    
    NK  = SIZE(IM, DIM=1)
    NKf = SIZE(IFMr, DIM=1)
    
    IFM = ZERO
    
    do kp = 1, NKf
        k2 = f2k(kp)
        do k = 1, NKf
            k1 = f2k(k)
            
            IFM(k1,k2) = IFMr(k,kp)
        end do
    end do
    
    IFM = MATMUL(IM, IFM) * normJs
    ! IFM = IFM * normJs
    
End Subroutine

Subroutine MEFEAInterfaceMass (IM, f2k, enodes, epnodes, IFMr, normJs, IIM)
    Implicit None
    ! Determine mass-eliminated interface mass inner product. Array must be pre-allocated
    Real (KREAL), Intent (In)  :: IM      (:,:)
    Integer,      Intent (In)  :: f2k     (:)
    Integer,      Intent (In)  :: enodes  (:)
    Integer,      Intent (In)  :: epnodes (:)
    Real (KREAL), Intent (In)  :: IFMr    (:,:)
    Real (KREAL), Intent (In)  :: normJs
    
    Real (KREAL), Intent (Out) :: IIM     (:,:)
    
    Integer                    :: k, kp, k1, k2, kg, kgp
    Integer                    :: NK
    Integer                    :: NKf
    
    NK  = SIZE(IM, DIM=1)
    NKf = SIZE(IFMr, DIM=1)
    
    IIM = ZERO
    
    do k = 1, NKf
        k1 = f2k(k)
        kg = enodes(k1)
        do kp = 1, NKf
            kgp = enodes(f2k(kp))
            k2 = FINDLOC(epnodes, kgp, DIM=1)
            
            IIM(k1,k2) = IFMr(k,kp)
        end do
    end do
    
    IIM = MATMUL(IM, IIM) * normJs
    
End Subroutine

End Module FEAInnerProducts