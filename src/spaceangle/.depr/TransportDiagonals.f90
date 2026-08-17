Subroutine TransportDiagonals (NE, NI, NK, el2mat, sl, khat, attn, n, IP, T)
    ! CURRENTLY DEPRECATED. MAY FIND USE SOON / EVENTUALLY.
    use constants
    use types
    use NittanyAPI
    use Geometry
    use IO
    use BasicMathFunctions
    Implicit None
    Integer, Intent (In) :: NE
    Integer, Intent (In) :: NI
    Integer, Intent (In) :: NK (:)
    Integer, Intent (In) :: el2mat (:)
    Integer, Intent (In) :: sl (:,:)
    Real (KREAL), Intent (In) :: khat (:,:)
    Real (KREAL),         Intent (In)    :: attn  (:)   ! (mat)
    Type (CRealM), Intent (In) :: n (:)
    Type (InnerProductType), Intent (In) :: IP (:)
    
    Type (CRealM),        Intent (InOut) :: T     (:,:) ! (es, i) ! MUST BE PRE-ALLOCATED
    
    Integer                              :: dir, e, es, f, i, k
    Integer                              :: info
    Integer                              :: mat
    Integer                              :: NF
    Integer                              :: d
    Integer                              :: nwork
    Integer :: NKe
    Real (KREAL)                         :: kdotn
    Real (KREAL)                         :: kvec    (3)
    Real (KREAL),         Allocatable    :: work    (:)
    Type (CIntV)                         :: ipiv    (8)   ! ipiv for BLAS routine
    Real (KREAL)                         :: ident (2,2)
    
    ident = ZERO
    ident(1,1) = ONE
    ident(2,2) = ONE
    
    d = SIZE(khat, DIM=1)
    
    ALLOCATE(ipiv(2)%v(2)) ! 1D linear elements
    ALLOCATE(ipiv(4)%v(4)) ! Tetrahedra
    ALLOCATE(ipiv(8)%v(8)) ! Hexahedra
    ALLOCATE(ipiv(6)%v(6)) ! Prism
    ALLOCATE(ipiv(5)%v(5)) ! Pyramid
    
    ! CAN PARALLELIZE OVER es AS WELL IF I DON'T ASSIGN KHAT...
    do i = 1, NI
        kvec(1:d) = khat(1:d,i)
        do es = 1, NE
            e = sl(es,i)
            
            NKe = NK(e)
            mat = el2mat(e)
            
            !  --------------
            !    Initialize
            !  --------------
            
            T(es,i)%m = ZERO
            
            !  ----------
            !    M term  
            !  ----------
            
            do k = 1, NKe
                T(es,i)%m(k,k) = attn(mat)
            end do
            
            !  ----------
            !    G term  
            !  ----------
            
            do dir = 1, d
                T(es,i)%m = T(es,i)%m - kvec(dir) * IP(e)%IC(1:NKe,1:NKe,dir)
            end do
            
            !  -------------
            !    F_up term
            !  -------------
            ! ALL AROUND CAN BE MADE MORE EFFICIENT...
            NF = SIZE(IP(e)%IFM) ! Inefficient???
            do f = 1, NF
                kdotn = DOT_PRODUCT(kvec(1:d), n(e)%m(:,f))
                if (kdotn < ZERO) cycle
                T(es,i)%m = T(es,i)%m + kdotn * IP(e)%IFM(f)%m
            end do
            
            call dgetrf (NKe, NKe, T(es,i)%m, NKe, ipiv(NKe)%v, info)
            
            nwork = -1
            ALLOCATE(work(1))
            call dgetri (NKe, T(es,i)%m, NKe, ipiv(NKe)%v, work, nwork, info)
            nwork = INT(work(1))
            DEALLOCATE(work)
            ALLOCATE(work(nwork))
            
            call dgetri (NKe, T(es,i)%m, NKe, ipiv(NKe)%v, work, nwork, info)
            
            DEALLOCATE(work)
            
        end do
    end do
    
    ! Deprecated slab code. faster of course but slab is meant for development only anyway and having the logic for slab and general be merged is better for validating general
    !else
    !    do i = 1, NI
    !        mu = prec%k(1,i)
    !        
    !        varsigma1 = MERGE(mu, ZERO, mu > ZERO)
    !        varsigma2 = MERGE(mu, ZERO, mu < ZERO)
    !        
    !        do es = 1, NE
    !            e = prec%sl(es,i)
    !        
    !            invdz = ONE / prec%vol(e)
    !            
    !            NK = 2
    !            
    !            mat = prec%el2mat(e)
    !            
    !            A = (-FOUR * varsigma2 + THREE * mu      ) * invdz + attn(mat)
    !            B = ( THREE * mu       - TWO * varsigma1 ) * invdz
    !            C = ( TWO * varsigma2  - THREE * mu      ) * invdz
    !            D = (-THREE * mu       + FOUR * varsigma1) * invdz + attn(mat)
    !            
    !            select case (solver)
    !            case (1)
    !                invdet = ONE / (A * D - B * C)
    !                
    !                T(es,i)%m(1,1) =  invdet * D
    !                T(es,i)%m(1,2) = -invdet * B
    !                T(es,i)%m(2,1) = -invdet * C
    !                T(es,i)%m(2,2) =  invdet * A
    !            case (2)
    !                T(es,i)%m(1,1) = A
    !                T(es,i)%m(1,2) = B
    !                T(es,i)%m(2,1) = C
    !                T(es,i)%m(2,2) = D
    !            end select
    !        end do
    !    end do
    !end if
    
End Subroutine