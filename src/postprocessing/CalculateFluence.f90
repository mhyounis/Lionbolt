Submodule (PostProcInterface) submodCalculateFluence
Contains
Module Subroutine CalculateFluence (s, w, fl)
    use constants
    Implicit None
    Real (KREAL), Intent (In)  :: s  (:)
    Real (KREAL), Intent (In)  :: w  (:)
    
    Real (KREAL), Intent (Out) :: fl (:)
    
    Integer                    :: i
    Integer                    :: NE
    Integer                    :: NI
    Integer                    :: NENK
    Integer                    :: start
    
    NE   = SIZE(fl)
    NI   = SIZE(w)
    NENK = SIZE(s) / NI
    
    fl = ZERO
    do i = 1, NI
        start = (i - 1) * NENK
        
        fl(1:NENK) = fl(1:NENK)                         &
                   + w(i) * s(1 + start : NENK + start)
    end do
    
End Subroutine
End Submodule