Module SPARSKITWrappers
    !  =======================================================
    !    This module contains wrappers to SPARSKIT routines.
    !    
    !    Given a linear operator object, a source, and a
    !    solution vector, it relies on SPARSKIT
    !    implementations of a given solver to solve the
    !    desired linear system.
    !  =======================================================
    use constants
    use profiler
    use IO
    use LinearOperatorClass
    ! Other file dependencies:
    ! iters.f
    Implicit None
    
Contains

Subroutine GMRESWrapper (m, s, Op, x, Precond, max_iter)
    Implicit None
    Integer,                           Intent (In)    :: m         ! Size of Krylov subspace (restart parameter)
    Type (SpaceAngleVector),           Intent (In)    :: s         ! Source
    Class (LinearOperator),            Intent (InOut) :: Op        ! Linear operator (according to implemented operator wrappers in LinearOperator class. See LinearOperators.f90)
    
    Type (SpaceAngleVector),           Intent (InOut) :: x         ! Guess and solution
    
    Class (LinearOperator),  Optional, Intent (InOut) :: Precond   ! Preconditioning operator ! WIP. Implemented, I think, but not tested nor validated.
    Integer,                 Optional, Intent (In)    :: max_iter  ! Maximum number of GMRES iterates
    
    Integer                                           :: iter
    Logical,                 Save                     :: firstrun = .TRUE.
    Logical                                           :: doprecond
    Integer                                           :: n
    Integer                                           :: wsize
    Integer                                           :: ipar (16)
    Real (KREAL)                                      :: fpar (16)
    Real (KREAL),            Allocatable              :: w    (:)
    Type (SpaceAngleVector)                           :: tx        ! Temporary w array for matvec ! MHY LATER - see if you can forego this. It gets allocated once every time a MatInv is called
    Character (512)                                   :: errmess
    
    if (firstrun) then
        write (iuout,'(7X,A)') 'Solver : GMRES'
        write (iuout,'(7X,A)') '[Ref. : Y. Saad & M. H. Schultz, SIAM J. Sci. Comput., 7, 856 (1986)]'
        firstrun = .FALSE.
    else
        write (iuout,'(7X,A)') 'Solver : GMRES'
    end if
    write (iuout,'(10X,A)') 'Operator : ' // Op%name
    
    doprecond = PRESENT(Precond)
    
    if (gmemprofiling) call RSSLogger ('--- IN GMRESWrapper --- Before allocating workspace')
    
    n     = SIZE(x%v)
    wsize = (n + 3) * (m + 2) + (m + 1) * m / 2
    ALLOCATE(w(wsize))
    ALLOCATE(tx%v(n))
    
    if (gmemprofiling) call RSSLogger ('--- IN GMRESWrapper --- After allocating workspace')
    
    if (wsize < 0) then
        errmess = 'GMRES with restart has been called, but the workspace array has a size exceeding ' // &
                  'the maximum single-precision integer value. This is a developmental issue and should ' // &
                  'be resolved soon. For now, make your mesh more coarse or lower your angular discretization.'
        call stophere (errmess)
    end if
    
    ipar(1) = 0
    if (doprecond) then
        ipar(2) = 1
    else
        ipar(2) = 0
    end if
    
    ipar(3) = 2
    !  -2 == || dx(i) || <= rtol * || s || + atol
    !  -1 == || dx(i) || <= rtol * || dx(1) || + atol
    !   0 == solver will choose test 1 (next)
    !   1 == || residual || <= rtol * || initial residual || + atol
    !   2 == || residual || <= rtol * || s || + atol
    ! 999 == caller will perform the test
    ! rtol = fpar(1) and atol = fpar(2)
    ipar(4) = wsize
    ipar(5) = m
    if (PRESENT(max_iter)) then
        ipar(6) = max_iter ! User input?
    else
        ipar(6) = -1 ! This makes it run until convergence. Probably use this in general.
    end if
    ipar(7) = 0 ! Number of matvec and precond applications, after terminating it gets written. So initialize to 0
                ! Maybe gets overwritten to 0 by bisinit. But fpar(11) doesn't?
    ! ipar(8) and ipar(9) are set by the program to inform what parts of w are used
    ! in the requested matrix operation.
    ! ipar(10) -- Set and used by program to know where to start after returning from a matrix operation request.
    ! ipar(11) -- Needed if caller will perform convergence, to inform the SPARSKIT subroutine
    ! to clean up and terminate. Not relevant to me.
    
    ! MHY - try changing these to 10^-5 and 10^-10 respectively
    fpar(1) = 1.0e-8_KREAL  ! Relative tolerance
    fpar(2) = 1.0e-12_KREAL ! Absolute tolerance ! Consider softening to 1.0e-10
    ! Meanings of ipar(3:7) depend on ipar(3)
    ! This will be good for printing though.
    !fpar(3) -- initial residual/error norm
    !fpar(4) -- target residual/error norm
    !fpar(5) -- current residual norm (if available) !!! THIS IS THE ONE I WANT.
    !fpar(6) -- current residual/error norm
    !fpar(7) -- convergence rate
    !fpar(8:10) -- Used by invidual subroutines. Do I need to set some? Maybe, but not for gmres.
    fpar(11) = 0 ! Number of FLOPS. So initialize to 0
    
    iter = 0
    write (iuout,'(10X,A)') '-------------------------------------'
    write (iuout,'(10X,A)') ' ITERATE           RESIDUAL NORM     '
    write (iuout,'(10X,A)') '---------     -----------------------'
    do
        iter = iter + 1
        
        call gmres (n, s%v, x%v, ipar, fpar, w)
        
        write (iuout,'(13X,I3,7X,ES24.16)') iter, fpar(5)
        
        if (ISNAN(fpar(5))) then
            call stophere ('RESIDUAL NORM is NaN.')
        end if
        
        ! Check ipar. No need for separate subroutine because this subroutine only exists once. That's the whole point of GMRESWrapper
        ! BUT : If I ever wanted to use more of sparskit's iter.f subroutines, then I would definitely want to
        ! consolidate the below and SOME of the above in its own subroutine. Actually would be a really good idea. 
        select case (ipar(1))
        case (1)
            ! IS THIS SLOW, REDUNDANT, ETC.?
            tx%v = w(ipar(8):ipar(8) + n - 1)
            call Op%MatVec (tx)
            w(ipar(9):ipar(9) + n - 1) = tx%v
        case (2)
            ! NOT USED IN GMRES
            ! Request matvec with A^T.
        case (3)
            ! Left precondition
            if (doprecond) then
                tx%v = w(ipar(8):ipar(8) + n - 1)
                call Precond%MatInv (tx)
                w(ipar(9):ipar(9) + n - 1) = tx%v
            end if
        ! case (4) --> Left precond transposed
        ! case (5) --> Right precond
        ! case (6) --> Right precond transposed
        case (10)
            ! CHECK TOLERANCE. Exit or set ipar(1) = 0?
        case (0)
            ! Normal termination
            write (iuout,*)                 ''
            write (iuout,'(10X,A)')         'CONVERGED.'
            write (iuout,'(10X,A,I3)')      'Total number of iterates : ', iter
            write (iuout,'(10X,A,ES24.16)') 'Final Residual Norm      : ', fpar(5)
            write (iuout,*)                 ''
            return
        case (-1)
            ! Too many iterates have been carried out
            errmess = 'GMRES : Iteration did not converge in specified iteration limits.'
            call stophere (errmess)
        case (-2)
            errmess = 'GMRES : Insufficient work space.' ! Will user ever encounter this?
            call stophere (errmess)
        case (-3)
            errmess = 'GMRES : Anticipated breakdown / division by zero.'
            ! CHECK ipar(12) NOW???
            call stophere (errmess)
        case (-4)
            ! SOME STUFF TO DO WITH fpar, REVISE LATER.
        case (-9)
            ! "while trying to detect a break-down, an abnormal number is detected."
            errmess = 'GMRES : '
        case (-10)
            ! "return due to some non-numerical reasons, e.g. invalid floating-point numbers etc."
            errmess = 'GMRES : '
        end select
    end do
    
End Subroutine

End Module SPARSKITWrappers