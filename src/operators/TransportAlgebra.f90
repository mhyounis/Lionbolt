Module TransportAlgebra
    use constants
    use Parallelism
    use types
    use Geometry
    use SpaceAngleInterface
    Implicit None
    
    Private :: LocalLU
    
    Interface ApplyTransport
        Module Procedure SNMGXS
        ! Module Procedure PNMGXS ! Should probably be a way to unify these, with indexing or with more objects/types/classes. obviously treat that later.
    End Interface
    
Contains
    
    Subroutine SNMGXS ()
        Implicit None
        
        print *, 'TransportOp%MatVec: WIP'
        stop
        
        ! if i recall correctly from wiscobolt this should not require sweep ordering
        ! I should be able to visit an element and do diagonal + neighboring faces easily
        
    End Subroutine
    
    Subroutine Sweep (d, NE, NI, NK, el2idx, el2mat, sl, attn, khat, n, IP, IAF, s)
        Implicit None
        Integer,                 Intent (In)    :: d              ! Dimension of physical space
        Integer,                 Intent (In)    :: NE             ! Number of elements
        Integer,                 Intent (In)    :: NI             ! Number of discrete ordinates
        Integer,                 Intent (In)    :: NK     (NE)    ! Number of nodes in an element (e)
        Integer,                 Intent (In)    :: el2idx (NE)    ! Gives FEM index from element
        Integer,                 Intent (In)    :: el2mat (NE)    ! Gives material index from element
        Integer,                 Intent (In)    :: sl     (NE,NI) ! Sweep list (es,i)
        Real (KREAL),            Intent (In)    :: attn   (:)     ! Attenuation coefficients (mat)
        Real (KREAL),            Intent (In)    :: khat   (d,NI)  ! Discrete ordinates set
        Type (CRealM),           Intent (In)    :: n      (NE)    ! Normal vectors (e)%m(dir, f)
        Type (InnerProductType), Intent (In)    :: IP     (7)     ! Mesh inner products
        Type (UpstreamSys),      Intent (In)    :: IAF            ! Describes off-diagonal elements of the transport matrix (outer index, i.e., mesh element to element)
        
        Real (KREAL),            Intent (InOut) :: s      (:)   ! Enters as source, leaves as solution
        
        Integer                                 :: e, es, f, fup, i, idx, j, k, le, mat, row
        Integer                                 :: ips
        Integer                                 :: ipe
        Integer                                 :: ipsf
        Integer                                 :: ipef
        Integer                                 :: NKe
        Integer                                 :: iafuf0
        Integer                                 :: nup
        Real (KREAL)                            :: kdotn
        Real (KREAL)                            :: lk    (d)
        Real (KREAL)                            :: tmp8  (8)
        Real (KREAL)                            :: Tdiag (8, 8)
        
        !  If I do GPU parallelization, I can further break up the elements into sweep blocks using:
        !  do si = 1, Ns(i)
        !  sumofNs = SUM(Ns(1:i-1))
        !  do es = ss(si + sumofNs) + 1, ss(si + sumofNs + 1)
        
        !$OMP PARALLEL DO DEFAULT (SHARED) &
        !$OMP PRIVATE (i, idx, le, es, e, fup, f, kdotn, nup, lk, NKe, mat, ips, ipe, ipsf, ipef, tmp8, Tdiag, k, row, j)
        do i = 1, NI
            lk = khat(:,i)
            
            ! Visit its dependent elements one-by-one
            do es = 1, NE
                e   = sl(es, i)
                NKe = NK(e)
                idx = el2idx(e)
                le  = IP(idx)%elmap(e)
                mat = el2mat(e)
                
                ! Now define the range of k values in the collapsed-indexed arrays to which this element corresponds
                ips = IAF%SAoffset(es,i) + 1
                ipe = IAF%SAoffset(es + 1,i)
                
                ! Visit each contributing face
                iafuf0 = IAF%UFoffset(es,i) ! Just storing this value to avoid having to grab it so many times
                nup    = IAF%UFoffset(es + 1,i) - iafuf0
                do fup = 1, nup
                    f = IAF%upfaces(iafuf0 + fup)
                    
                    ipsf = IAF%ipsf(iafuf0 + fup)
                    ipef = IAF%ipef(iafuf0 + fup)
                    
                    kdotn = DOT_PRODUCT(lk, n(e)%m(:,f)) ! n should really also follow the idx thing, since the only reason it's stored like this is the NFe dimension may change
                    
                    tmp8(1:NKe) = - kdotn * MATMUL(IP(idx)%IIM(:,:,f,le), s(ipsf:ipef))
                    s(ips:ipe)  = s(ips:ipe) + tmp8(1:NKe)
                end do
                
                ! Apply inverse of the diagonal block
                tmp8(1:NKe) = s(ips:ipe)
                call ApplyInvertedDiagonal (NKe, attn(mat), lk, n(e)%m, IP(idx)%IC(:,:,:,le), IP(idx)%IFM(:,:,:,le), Tdiag, tmp8)
                s(ips:ipe) = tmp8(1:NKe)
                
            end do
            
        end do
        !$OMP END PARALLEL DO
        
    End Subroutine
    
    Subroutine ApplyInvertedDiagonal (NKe, attn, khat, n, IC, IFM, Tdiag, tmp8)
        Implicit None
        Integer,      Intent (In)    :: NKe
        Real (KREAL), Intent (In)    :: attn
        Real (KREAL), Intent (In)    :: khat  (:)
        Real (KREAL), Intent (In)    :: n     (:,:)
        Real (KREAL), Intent (In)    :: IC    (:,:,:)
        Real (KREAL), Intent (In)    :: IFM   (:,:,:)
        Real (KREAL), Intent (InOut) :: Tdiag (8, 8)
        Real (KREAL), Intent (InOut) :: tmp8  (8)
        
        Integer                      :: dir, f, k
        Real (KREAL)                 :: kdotn
        ! Integer                       :: ipiv (8) ! TEMPORARY
        ! Integer                       :: info     ! TEMPORARY
        
        ! ALL AROUND COULD BE MADE MORE EFFICIENT...
        
        Tdiag(1:NKe, 1:NKe) = ZERO
        
        !  ----------
        !    M term  
        !  ----------
        
        do k = 1, NKe
            Tdiag(k,k) = attn
        end do
        
        !  ----------
        !    G term  
        !  ----------
        
        do dir = 1, SIZE(khat)
            Tdiag(1:NKe,1:NKe) = Tdiag(1:NKe,1:NKe) - khat(dir) * IC(:,:,dir)
        end do
        
        !  -------------
        !    F_up term  
        !  -------------
        
        do f = 1, SIZE(IFM, DIM=3)
            kdotn = DOT_PRODUCT(khat, n(:,f))
            if (kdotn < ZERO) cycle
            Tdiag(1:NKe,1:NKe) = Tdiag(1:NKe,1:NKe) + kdotn * IFM(:,:,f)
        end do
        
        call LocalLU (NKe, Tdiag, tmp8)
        
        ! ! MHY LATER - TESTING
        ! if (attn  / abs(Tdiag(1,1) - attn) > 1.0e-6_KREAL) then
        !     call LocalLU (NKe, Tdiag, tmp8)
        ! else
        !     ! For unstable systems use a more stable method
        !     ! If this doesn't get hit for air in 3D meshes then I really shouldn't use it.
        !     ! Double check criterion too...
        !     ! What I'll do is try the air calculation with ONLY this turned on.
        !     ! If this doesn't make a difference I can profile the two techniques here separately.
        !     ! If this is slower then just comment it out and go with LocalLU.
        !     call dgetrf (NKe, NKe, Tdiag(1:NKe,1:NKe), NKe, ipiv, info)
        !     call dgetrs ('N', NKe, 1, Tdiag(1:NKe,1:NKe), NKe, ipiv, tmp8(1:NKe), NKe, info)
        ! end if
        ! ! ! !
        
    End Subroutine
    
    Subroutine LocalLU (n, A, b)
        ! Just a small, quick LU inversion routine. Could use this elsewhere too. We'll see
        ! Re-wrote from LU_decomposition, originally in wiscobolt/math.f08.
        ! Kept here for speed of access via ApplyInvertedDiagonal.
        ! MHY LATER - probably use this in FiniteElementAnalysis too. Oh yea that's a great idea
        ! Also consider demanding that the temporary array x (8) is stored somewhere permanently.
        ! ALSO NOTE - If cross sections ever go to zero, then this routine will fail,
        !             because the matrix A will be poorly conditioned.
        !             As well as some other cases
        Implicit None
        Integer,      Intent (In)    :: n
        Real (KREAL), Intent (InOut) :: A (8, 8)
        Real (KREAL), Intent (InOut) :: b (8)
        
        Integer                      :: i, j, k
        Integer                      :: pivot
        Real (KREAL)                 :: maxval
        Real (KREAL)                 :: t
        Real (KREAL)                 :: factor
        Real (KREAL)                 :: x (8)
        
        ! Could I/should I do some of this matrix-less? 8 x 8 is not that expensive
        
        do k = 1, n - 1
            
            maxval = ABS(A(k,k))
            pivot  = k
            
            do i = k + 1, n
                if (ABS(A(i,k)) > maxval) then
                    maxval = ABS(A(i,k))
                    pivot  = i
                end if
            end do
            
            ! Swap rows if necessary
            if (pivot /= k) then
                A([k,pivot],1:8) = A([pivot,k],1:8)
                
                t        = b(k)
                b(k)     = b(pivot)
                b(pivot) = t
            end if
            
            ! Eliminate below pivot
            do i = k + 1, n
                
                factor = A(i,k) / A(k,k)
                A(i,k) = factor
                
                do j = k + 1, n
                    A(i,j) = A(i,j) - factor * A(k,j)
                end do
                
                b(i) = b(i) - factor * b(k)
                
            end do
            
        end do
        
        ! Back substitution
        x(n) = b(n) / A(n,n)
        do i = n - 1, 1, -1
            x(i) = b(i)
            do j = i + 1, n
                x(i) = x(i) - A(i,j) * x(j)
            end do
            x(i) = x(i) / A(i,i)
        end do
        
        ! Copy solution back to RHS
        b = x
        
    End Subroutine
    
End Module TransportAlgebra