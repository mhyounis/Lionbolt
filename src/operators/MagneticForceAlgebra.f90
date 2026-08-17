Module MagneticForceAlgebra
    use constants
    use types
    use Geometry
    use Parallelism
    Implicit None
    
    Interface ApplyMagneticForce
        Module Procedure SNMGXS
        ! Module Procedure PNMGXS ! Should probably be a way to unify these, with indexing or with more objects/types/classes. obviously treat that later.
    End Interface
    
    Contains
    
    Subroutine SNMGXS (NENK, NI, invp, A, B, k, xw, x)
        Implicit None
        Integer,              Intent (In)    :: NENK
        Integer,              Intent (In)    :: NI
        Real (KREAL),         Intent (In)    :: invp
        Real (KREAL),         Intent (In)    :: A   (:,:)
        Real (KREAL),         Intent (In)    :: B   (:,:)
        Real (KREAL),         Intent (In)    :: k   (:,:)
        Real (KREAL),         Intent (InOut) :: xw  (:,:)
        
        Real (KREAL), Target, Intent (InOut) :: x   (:)
        
        Integer                              :: i, ip
        Integer                              :: d
        Real (KREAL), Pointer                :: xwp (:,:)
        
        d = SIZE(k, dim=1)
        
        ! Map to pointer
        xwp(1:NENK, 1:NI) => x
        
        ! Form A * x
        call dgemm ('N','N', NENK, NI, NI, ONE, xwp, NENK, A(1,1), NI, ZERO, xw, NENK)
        
        xwp = xw
        
        ! Note to self - dgemm will REPLACE xw if third to last argument is ZERO.
        
        ! Form B * x
        
    End Subroutine
    
End Module