Module BoltzmannOperator
    !  ========================================================================
    !    This module implements the Boltzmann operator as an extension of the
    !    abstract LinearOperator class.
    !    Read the docs at the OpenRPS website if you are interested in 
    !    learning more.
    !  ========================================================================
    use constants
    use profiler
    use NittanyAPI
    use Geometry
    use AngularSpace
    use LinearOperatorClass
    use TransportOperator
    use TransportAlgebra
    use ScatteringOperator
    use ScatteringAlgebra
    Implicit None
    Private
    Public :: BoltzmannOp
    
    ! Operator: L = T - K OR T^{-1}L = 1 - T^{-1}K (if Tprecond is true, which will be set) ! MHY LATER - must make Tprecond user-accessible...
    Type, Extends (LinearOperator) :: BoltzmannOp
        ! We note, this operator's phase space objects are never associated
        ! General
        Logical                        :: Tprecond ! Whether to use L = T - K or T^{-1}L = 1 - T^{-1}K
        Character (LEN=:), Allocatable :: solver   ! Which solver to use
        Type (SpaceAngleVector)        :: xw       ! Workspace vector for applying T and K separately
        Type (SpaceAngleVector)        :: xg       ! Guess vector
        ! Transport
        Type (TransportOp)             :: T
        ! Scattering
        Type (ScatteringOp)            :: K
    Contains
        Procedure :: MatVec     => MatVec
        Procedure :: MatInv     => MatInv
        Procedure :: Full       => Full
        Procedure :: Space      => Space
        Procedure :: SpaceAngle => SpaceAngle
        Procedure :: Physics    => Physics
        Procedure :: Destroy    => Destroy
        Generic   :: Build      => Full, Space, SpaceAngle, Physics
    End Type
    
    ! ! COARSE ! Operator: L = T - K OR T^{-1}L = 1 - T^{-1}K (if SN)
    ! Type, Extends (OpL) :: OpLc
    !     ! MAKE L HERE. Or perhaps make a coarseoperator type which takes any operator
    !     ! and also includes the coarse2fine stuff
    ! Contains
    !     Procedure :: MatVec => MatVec_c
    !     Procedure :: MatInv => MatInv_c
    ! End Type
    ! Then, in the MatVec and MatInv, use restriction operator (fine2coarse) and prolongation operator (coarse2fine)
    
    ! Could end up revisiting the Space subroutines if I ever do adaptive mesh refinement. May have to consider working in single-element operators though. Or something like that
    
    Contains
    
    Subroutine MatVec (Op, x)
        use constants
        use LinearOperatorClass
        use TransportAlgebra,  ONLY : ApplyTransport, Sweep
        use ScatteringAlgebra, ONLY : ApplyScattering
        Implicit None
        ! WRAPPER SUBROUTINE
        ! This subroutine applies L (PN) or T^{-1}L (SN) to a vector, x, in collapsed notation.
        ! That is, x is a simultaneous space-angle array, with every linear entry mapping a unique
        ! triplet of (e,k,i).
        Class (BoltzmannOp),     Intent (InOut) :: Op
        
        Type (SpaceAngleVector), Intent (InOut) :: x
        
        if (.not. (Op%SPACEREADY .AND. Op%ANGLEREADY .AND. Op%PHYSREADY)) &
            call stophere ('BoltzmannOp%MatVec was called, but the operator was not ready.')
        
        Op%xw = x
        if (Op%Tprecond) then
            ! Feed x
            if (gtiming) call Timestamp ('Scattering')
            call Op%K%MatVec (x)
            if (gtiming) call Timestamp ('Scattering')
            ! Returned Kx
            
            ! Feed Kx
            if (gtiming) call Timestamp ('Sweeping')
            call Op%T%MatInv (x)
            if (gtiming) call Timestamp ('Sweeping')
            ! Returned T^{-1}Kx
            
            ! Now form x - T^{-1}Kx
        else
            ! Feed xw
            call Op%T%MatVec (Op%xw)
            ! Returned Tx
            
            ! Feed x
            call Op%K%MatVec (x)
            ! Returned Kx
            
            ! Now form (T - K)x
        end if
        x%v = Op%xw%v - x%v
        
    End Subroutine
    
    Subroutine MatInv (Op, s)
        use constants
        use SPARSKITWrappers
        use SolversInterface
        Implicit None
        ! NOT A WRAPPER SUBROUTINE.
        ! This subroutine determines L^{-1}s, returning s, which was initialized
        ! as the source, as the solution.
        Class (BoltzmannOp),     Intent (InOut) :: Op
        
        Type (SpaceAngleVector), Intent (InOut) :: s
        
        if (.not. (Op%SPACEREADY .AND. Op%ANGLEREADY .AND. Op%PHYSREADY)) &
            call stophere ('BoltzmannOp%MatInv was called, but the operator was not ready.')
        
        if (Op%Tprecond) then
            call Op%T%MatInv (s)
        end if
        
        select case (Op%solver)
        case ('GMRES')
            ! GMRES
            
            ! For now, I'll just take the guess to be the source divided by
            ! the attenuation of material 1 (arbitrary choice of material)
            ! Maybe later devise a fast way to assign s / attn depending on the
            ! element's material
            Op%xg%v = s%v / Op%T%XS%t(1)
            
            call GMRESWrapper (10, s, Op, Op%xg)
            
            s = Op%xg
            
        case ('SI')
            ! SOURCE ITERATION
            
            if (.not. Op%Tprecond) then
                ! If you didn't do T precond, you must do it
                call Op%T%MatInv (s)
            end if
            
            call SourceIteration (Op%T, Op%K, s)
            
        end select
        
    End Subroutine
    
    Subroutine Full (Op, mesh, angular, XS)
        Implicit None
        Class (BoltzmannOp),         Intent (InOut) :: Op
        
        Type (MeshClass),    Target, Intent (In)    :: mesh
        Type (AngularClass), Target, Intent (In)    :: angular
        Type (XSType),       Target, Intent (In)    :: XS
        
        ! Give the operator its name
        Op%name = 'Boltzmann Operator'
        
        ! Now construct the dimensional parts of L
        call Space      ( Op, mesh          )
        call SpaceAngle ( Op, mesh, angular )
        call Physics    ( Op, angular, XS   )
        
    End Subroutine
    
    Subroutine Space (Op, mesh)
        Implicit None
        Class (BoltzmannOp),      Intent (InOut) :: Op
        
        Type (MeshClass), Target, Intent (In)    :: mesh
        
        if (.not. mesh%READY) then
            call stophere ('BoltzmannOp.f90: Space: Mesh used to build a Boltzmann operator is not READY. ' // &
                           'Make sure to initialize & post-process before use.')
        end if
        
        ! Give the operator its name
        Op%name = 'Boltzmann Operator'
        
        !  ======================
        !    Transport Operator  
        !  ======================
        if (gmemprofiling) call RSSLogger (' --- IN BOLTZMANN OPERATOR --- Before T space')
        call Op%T%Build (mesh)
        if (gmemprofiling) call RSSLogger (' --- IN BOLTZMANN OPERATOR --- After T space')
        !  =======================
        !    Scattering Operator  
        !  =======================
        if (gmemprofiling) call RSSLogger (' --- IN BOLTZMANN OPERATOR --- Before K space')
        call Op%K%Build (mesh)
        if (gmemprofiling) call RSSLogger (' --- IN BOLTZMANN OPERATOR --- After K space')
        !  ===========================================
        !    Declare the space portion of K as ready  
        !  ===========================================
        
        Op%SPACEREADY = .TRUE.
        
    End Subroutine
    
    Subroutine SpaceAngle (Op, mesh, angular)
        Implicit None
        Class (BoltzmannOp),         Intent (InOut) :: Op
        
        Type (MeshClass),    Target, Intent (In)    :: mesh
        Type (AngularClass), Target, Intent (In)    :: angular
        
        Integer                                     :: NENK
        Integer                                     :: NI
        
        NENK = mesh%NENK
        NI   = angular%NI
        
        ! Give the operator its name
        Op%name = 'Boltzmann Operator'
        
        !  ================================================
        !    Form the space part if it hasn't been formed  
        !  ================================================
        
        if (.not. Op%SPACEREADY) then
            call Space (Op, mesh)
        end if
        
        !  ===========================================
        !    Allocate or resize the workspace vector  
        !  ===========================================
        
        if (ALLOCATED(Op%xw%v)) then
            ! If xw is already allocated, check if you need to change its size
            if (Op%xw%NENK /= NENK .or. Op%xw%NI /= NI) then
                
                Op%xw = SpaceAngleVector (mesh, angular)
                
            end if
        else
            if (gmemprofiling) call RSSLogger (' --- IN BOLTZMANN OPERATOR --- Before allocating xw')
            ! If xw is not allocated, allocate it
            Op%xw = SpaceAngleVector (mesh, angular)
            if (gmemprofiling) call RSSLogger (' --- IN BOLTZMANN OPERATOR --- After allocating xw')
        end if
        
        !  =======================================
        !    Allocate or resize the guess vector  
        !  =======================================
        
        if (ALLOCATED(Op%xg%v)) then
            ! If xg is already allocated, check if you need to change its size
            if (Op%xg%NENK /= NENK .or. Op%xg%NI /= NI) then
                
                Op%xg = SpaceAngleVector (mesh, angular)
                
            end if
        else
            if (gmemprofiling) call RSSLogger (' --- IN BOLTZMANN OPERATOR --- Before allocating xg')
            ! If xg is not allocated, allocate it
            Op%xg = SpaceAngleVector (mesh, angular)
            if (gmemprofiling) call RSSLogger (' --- IN BOLTZMANN OPERATOR --- After allocating xg')
        end if
        
        !  ======================
        !    Transport Operator  
        !  ======================
        if (gmemprofiling) call RSSLogger (' --- IN BOLTZMANN OPERATOR --- Before building T spaceangle')
        call Op%T%Build (mesh, angular)
        if (gmemprofiling) call RSSLogger (' --- IN BOLTZMANN OPERATOR --- After building T spaceangle')
        !  =======================
        !    Scattering Operator  
        !  =======================
        if (gmemprofiling) call RSSLogger (' --- IN BOLTZMANN OPERATOR --- Before building K spaceangle')
        call Op%K%Build (mesh, angular)
        if (gmemprofiling) call RSSLogger (' --- IN BOLTZMANN OPERATOR --- After building K spaceangle')
        ! As constructed here, energy groups all will use the same angular disc and solver. Could revisit this later.
        !  ==================
        !    Set SN logical  
        !  ==================
        
        Op%Tprecond = angular%disc == 'SN'
        
        !  ==============
        !    Set solver  
        !  ==============
        
        Op%solver = angular%solver
        
        !  =================================================
        !    Declare the space-angle portion of L as ready  
        !  =================================================
        
        Op%ANGLEREADY = .TRUE.
        
    End Subroutine
    
    Subroutine Physics (Op, angular, XS)
        Implicit None
        Class (BoltzmannOp),         Intent (InOut) :: Op
        
        Type (AngularClass), Target, Intent (In)    :: angular
        Type (XSType),       Target, Intent (In)    :: XS
            
        ! Give the operator its name
        Op%name = 'Boltzmann Operator'
        
        !  ======================
        !    Transport Operator  
        !  ======================
        
        call Op%T%Build (angular, XS)
        
        !  =======================
        !    Scattering Operator  
        !  =======================
        
        call Op%K%Build (angular, XS)
        
        !  =============================================
        !    Declare the physics portion of L as ready  
        !  =============================================
        
        Op%PHYSREADY = .TRUE.
        
    End Subroutine
    
    Subroutine Destroy (Op)
        Implicit None
        Class (BoltzmannOp), Intent (InOut) :: Op
        
        if (ALLOCATED(Op%solver)) DEALLOCATE(Op%solver)
        call Op%xw%Destroy ()
        call Op%xg%Destroy ()
        call Op%T%Destroy  ()
        call Op%K%Destroy  ()
        
    End Subroutine
    
End Module BoltzmannOperator