Submodule (SpaceAngleInterface) submodEncodeSweep
Contains
Module Subroutine EncodeSweep (encode, NE, NK, NENK, NI, offset, SAoffset, sl, t, s)
    use constants
    Implicit None
    Logical,      Intent (In)    :: encode                ! Whether to encode (TRUE) or decode (FALSE) the sweep
    Integer,      Intent (In)    :: NE
    Integer,      Intent (In)    :: NK       (NE)
    Integer,      Intent (In)    :: NENK
    Integer,      Intent (In)    :: NI
    Integer,      Intent (In)    :: offset   (NE + 1)     ! Generic offset
    Integer,      Intent (In)    :: SAoffset (NE + 1, NI) ! Space-angle offset
    Integer,      Intent (In)    :: sl       (NE, NI)     ! Sweep list
    
    Real (KREAL), Intent (InOut) :: t        (:)          ! Pre-allocated temporary array
    Real (KREAL), Intent (InOut) :: s        (:)
    
    Integer                      :: e, es, i
    Integer                      :: NKe
    Integer                      :: os1
    Integer                      :: os2
    Integer                      :: ios
    
    ! MHY LATER - really want to do something to optimize this. Not sure. It is already decently fast I think
    ! Think pointers?
    
    if (encode) then
        !$OMP PARALLEL DO DEFAULT (SHARED) PRIVATE (i, ios, es, os1, os2, e, NKe) SCHEDULE (STATIC)
        do i = 1, NI
            ios = (i - 1) * NENK
            
            do es = 1, NE
                os2 = SAoffset(es,i)
                
                e   = sl(es,i)
                NKe = NK(e)
                os1 = offset(e) + ios
                
                t(os2 + 1 : os2 + NKe) = s(os1 + 1 : os1 + NKe)
            end do
            
        end do
        !$OMP END PARALLEL DO
    else
        !$OMP PARALLEL DO DEFAULT (SHARED) PRIVATE (i, ios, es, os1, os2, e, NKe) SCHEDULE (STATIC)
        do i = 1, NI
            ios = (i - 1) * NENK
            
            do es = 1, NE
                os2 = offset(es) + ios
                
                e   = sl(es,i)
                NKe = NK(es)
                os1 = SAoffset(e,i) ! MHY LATER - seems wrong but works for case of all-tetrahedra mesh. Could be worth looking into flipping some stuff to get agreement in general case, if there's an issue
                
                t(os1 + 1 : os1 + NKe) = s(os2 + 1 : os2 + NKe)
            end do
            
        end do
        !$OMP END PARALLEL DO
    end if
    
    s = t
    
End Subroutine
End Submodule