Module ScatteringAlgebra
    use constants
    use types
    Implicit None
    
    Interface ApplyScattering
        Module Procedure SNMGXS
        Module Procedure PNMGXS
    End Interface
    
    Contains
    
    Subroutine SNMGXS (Nm, NI, NENK, YYS, mat2sd, xw, x)
        Implicit None
        Integer,              Intent (In)    :: Nm
        Integer,              Intent (In)    :: NI
        Integer,              Intent (In)    :: NENK
        Real (KREAL),         Intent (In)    :: YYS    (:,:,:) ! Angular scattering kernel                                      ! (i', i, mat)
        Type (CIntV),         Intent (In)    :: mat2sd (:)     ! Gives the set of spatial d.o.f.s corresponding to material mat ! (mat)
        Real (KREAL),         Intent (InOut) :: xw     (:,:)   ! Workspace array, rank two for re-ordering
        
        Real (KREAL), Target, Intent (InOut) :: x (:)          ! In as angular fluence, out as source
        
        Integer                              :: im2sd, sdof, mat
        Real (KREAL), Pointer                :: xwp  (:,:)
        Real (KREAL)                         :: xtmp (NI)
        Real (KREAL)                         :: ytmp (NI)
        
        ! Map to pointer
        xwp(1:NENK, 1:NI) => x
        
        if (Nm == 1) then
            ! Execute matmul using workspace array
            !call openblas_set_num_threads(16)
            call dgemm ('N','N', NENK, NI, NI, ONE, xwp, NENK, YYS(1,1,1), NI, ZERO, xw, NENK)
            !call openblas_set_num_threads(1)
        else
            ! MHY LATER - Huge opportunity for optimization.
            do mat = 1, Nm
                do im2sd = 1, SIZE(mat2sd(mat)%v)
                    sdof = mat2sd(mat)%v(im2sd)
                    
                    xtmp = xwp(sdof, 1:NI)
                    call dgemv ('T', NI, NI, ONE, YYS(:,:,mat), NI, xtmp, 1, ZERO, ytmp, 1) ! Transpose here because of how YYS is defined, i.e., YYS(i', i, m)
                    
                    xw(sdof, 1:NI) = ytmp ! This construct could be used to avoid xw entirely... Generalize to Nm = 1 too?
                end do
            end do
            
        end if
        
        ! Map back to x via pointer
        xwp = xw
        
    End Subroutine
    
    Subroutine PNMGXS ()
        Implicit None
        
    End Subroutine
    
    Subroutine ApplyDeltaDown (Nm, NI, NENK, XSdd, mat2sd, x)
        Implicit None
        Integer,      Intent (In)    :: Nm
        Integer,      Intent (In)    :: NI
        Integer,      Intent (In)    :: NENK
        Real (KREAL), Intent (In)    :: XSdd   (:)     ! Material-dependent delta down cross section                          ! (mat)
        Type (CIntV), Intent (In)    :: mat2sd (:)     ! Gives the set of true elements x nodes corresponding to material mat ! (mat)
        
        Real (KREAL), Intent (InOut) :: x      (:)     ! In as angular fluence, out as source
        
        Integer                      :: i, ip, im2sd, mat
        
        ! In this routine we use YYS(i,i',mat) = XSdd(mat) * delta(i,i')
        
        if (Nm == 1) then
            call dscal (NENK * NI, XSdd(1), x, 1) ! Fast vector scaling
        else
            do mat = 1, Nm
                do im2sd = 1, SIZE(mat2sd(mat)%v)
                    ip = mat2sd(mat)%v(im2sd)
                    
                    ! MHY LATER - EXTREMELY SLOW, ANOTHER OPTIMIZATION OPPORTUNITY
                    ! idea here is to visit the sdof in every angle and make a point-wise multiplication.
                    do i = 1, NI
                        x(ip + (i - 1) * NENK) = x(ip + (i - 1) * NENK) * XSdd(mat)
                    end do
                end do
            end do
        end if
        
    End Subroutine
    
End Module ScatteringAlgebra