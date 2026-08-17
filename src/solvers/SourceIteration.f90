Submodule (SolversInterface) submodSourceIteration
Contains
Subroutine SourceIteration (OpT, OpK, s)
    use constants
    use profiler
    use LinearOperatorClass
    Implicit None
    Class (LinearOperator),  Intent (InOut) :: OpT   ! Must send these in to avoid circular dependency with BoltzmannOperator. (While I don't have to declare the Op as BoltzmannOp, I would need to use BoltzmannOp's T and K attributes, which are not part of the abstract LinearOperator class)
    Class (LinearOperator),  Intent (InOut) :: OpK   ! Unfortunate. But this whole routine could be described as a more general leading/lagging operator thing
    
    Type (SpaceAngleVector), Intent (InOut) :: s ! This must be sent in as T^{-1}s ! Should consider not making this generally the case, and just having this routine take a Tprecond logical. Make sure to keep this routine a general lagging routine
    
    Integer                                 :: iter
    Real (KREAL),            Parameter      :: rTOL = 1.0e-8_KREAL  ! Relative tolerance
    Real (KREAL),            Parameter      :: aTOL = 1.0e-12_KREAL ! Absolute tolerance
    Integer                                 :: N
    Real (KREAL)                            :: norms ! Norm of source initially
    Real (KREAL)                            :: resnorm
    Type (SpaceAngleVector)                 :: w     ! Temporary vector to store the contribution of the Mth iterate
    
    N = SIZE(s%v)
    
    ALLOCATE(w%v(N))
    
    write (iuout,'(7X,A)')  'Solver : SI'
    write (iuout,'(10X,A)') 'Leading Operator : ' // OpT%name
    write (iuout,'(10X,A)') 'Lagging Operator : ' // OpK%name
    
    write (iuout,'(10X,A)') '-------------------------------------'
    write (iuout,'(10X,A)') ' ITERATE           RESIDUAL NORM     '
    write (iuout,'(10X,A)') '---------     -----------------------'
    
    ! Initialize outside of the loop, because of how residual is checked.
    iter = 1
    
    w%v = s%v
    norms = NORM2(s%v)
    
    ! Now calculate the contribution of this iterate
    
    call OpK%MatVec (w)
    call OpT%MatInv (w)
    
    ! Add it to running tally
    s%v = s%v + w%v
    
    write (iuout,'(13X,A)') '  1         ----------------------' ! No residual check for first iterate
    
    do
        iter = iter + 1
        
        ! Residual is checked first because the residual of the Mth iterate is actually the contribution of
        ! the (M-1)th iterate, so that can be checked immediately.
        
        resnorm = NORM2(w%v)
        
        write (iuout,'(13X,I3,7X,ES24.16)') iter, resnorm
        
        if (resnorm <= aTOL + norms * rTOL) then
            write (iuout,*)                 ''
            write (iuout,'(10X,A)')         'CONVERGED.'
            write (iuout,'(10X,A,I3)')      'Total number of iterates : ', iter
            write (iuout,'(10X,A,ES24.16)') 'Final Residual Norm      : ', resnorm
            write (iuout,*)                 ''
            return
        end if
        
        ! Now calculate the contribution of this iterate
        if (gtiming) call Timestamp ('Scattering')
        call OpK%MatVec (w)
        if (gtiming) call Timestamp ('Scattering')
        
        if (gtiming) call Timestamp ('Sweeping')
        call OpT%MatInv (w)
        if (gtiming) call Timestamp ('Sweeping')
        
        ! Add it to running tally
        s%v = s%v + w%v
        
    end do
    
End Subroutine
End Submodule