Module TransportOperator
    !  ========================================================================
    !    This module implements the transport operator as an extension of the
    !    abstract LinearOperator class.
    !    Read the docs at the OpenRPS website if you are interested in 
    !    learning more.
    !  ========================================================================
    use constants
    use profiler
    use types
    use NittanyAPI
    use Geometry
    use AngularSpace
    use SpaceAngleInterface
    use LinearOperatorClass
    use TransportAlgebra
    Implicit None
    Private
    Public :: TransportOp
    
    ! At the moment this operator assumes vacuum boundary condition (solution is zero at boundary for incident angles)
    ! While for convex external beam problems we can (and do) transform our problem to reach this point,
    ! this will likely need to be generalized somehow for problems with cavities, vacua, and internal sources.
    ! There will be a more generalized framework to specify and work with boundary conditions in the future.
    
    ! Operator: T (Transport)
    Type, Extends (LinearOperator)           :: TransportOp
        Integer,                 Allocatable :: NK       (:)   ! Number of nodes in element e
        Integer,                 Allocatable :: el2idx   (:)
        Integer,                 Allocatable :: el2mat   (:)
        Integer,                 Allocatable :: sl       (:,:) ! Sweep list                              ! (es, i)
        ! Integer,                 Allocatable :: isl      (:,:) ! Inverse sweep list                      ! (e, i)
        Character (LEN=:),       Allocatable :: solver         ! Solver for T^{-1}
        Type (CRealM),           Allocatable :: n        (:)   ! Normal vectors                          ! (e)%m(dir, f)
        Type (UpstreamSys)                   :: IAF
        Type (SpaceAngleVector)              :: xw             ! Workspace vector for sweep encoding
        Type (SpaceAngleVector)              :: xg             ! Guess vector
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
    
    Contains
    
    Subroutine MatVec (Op, x)
        use constants
        use LinearOperatorClass
        use TransportAlgebra, ONLY : ApplyTransport
        Implicit None
        ! WRAPPER SUBROUTINE
        ! The raw subroutine is ApplyTransport.
        ! This subroutine applies T to a vector, x, in collapsed notation.
        Class (TransportOp),     Intent (InOut) :: Op
        
        Type (SpaceAngleVector), Intent (InOut) :: x
        
        if (.not. (Op%SPACEREADY .AND. Op%ANGLEREADY .AND. Op%PHYSREADY)) &
            call stophere ('TransportOp%MatVec was called, but the operator was not ready.')
        
        ! MHY LATER - Could probably easily implement.
        
        call stophere ('TransportOp%MatVec was called, but the operator has yet to be implemented.')
        
    End Subroutine
    
    Subroutine MatInv (Op, s)
        use constants
        use LinearOperatorClass
        use SPARSKITWrappers
        use SpaceAngleInterface, ONLY : EncodeSweep
        use TransportAlgebra, ONLY : Sweep
        Implicit None
        ! WRAPPER SUBROUTINE
        ! This subroutine
        Class (TransportOp),     Intent (InOut) :: Op
        
        Type (SpaceAngleVector), Intent (InOut) :: s
        
        Integer                                 :: d
        
        if (.not. (Op%SPACEREADY .AND. Op%ANGLEREADY .AND. Op%PHYSREADY)) &
            call stophere ('TransportOp%MatInv was called, but the operator was not ready.')
        
        d = MERGE(1, 3, Op%mesh%slab) ! Dimensionality of space
        
        select case (Op%solver)
        case ('GMRES')
            ! GMRES
            Op%xg%v = s%v / Op%XS%t(1)
            
            call GMRESWrapper (10, s, Op, Op%xg) ! CHANGE 10? MAKE AN OPTION???
            
            s = Op%xg
            
        case ('SWEEP')
            ! SWEEP
            
            ! Encode sweep
            call EncodeSweep (  &
                .TRUE.,         &
                Op%mesh%NE,     &
                Op%NK,          &
                Op%mesh%NENK,   &
                Op%angular%NI,  &
                Op%mesh%offset, &
                Op%IAF%SAoffset, &
                Op%sl,          &
                Op%xw%v,        &
                s%v             &
            )
            
            ! Do the sweep
            call Sweep (       &
                d,             &
                Op%mesh%NE,    &
                Op%angular%NI, &
                Op%NK,         &
                Op%el2idx,     &
                Op%el2mat,     &
                Op%sl,         &
                Op%XS%t,       &
                Op%angular%k,  &
                Op%n,          &
                Op%mesh%IP,    &
                Op%IAF,        &
                s%v            &
            )
            
            ! Decode sweep
            call EncodeSweep (   &
                .FALSE.,         &
                Op%mesh%NE,      &
                Op%NK,           &
                Op%mesh%NENK,    &
                Op%angular%NI,   &
                Op%mesh%offset,  &
                Op%IAF%SAoffset, &
                Op%sl,           &
                Op%xw%v,         &
                s%v              &
            )
            
        end select
        
    End Subroutine
    
    Subroutine Full (Op, mesh, angular, XS)
        ! Creates, or re-writes, the transport operator, whether or not new objects are being passed.
        Implicit None
        Class (TransportOp),         Intent (InOut) :: Op
        
        Type (MeshClass),    Target, Intent (In)    :: mesh
        Type (AngularClass), Target, Intent (In)    :: angular
        Type (XSType),       Target, Intent (In)    :: XS
            
        ! Give the operator its name
        Op%name = 'Transport Operator'
        
        ! Now construct the dimensional parts of T
        call Space      ( Op, mesh          )
        call SpaceAngle ( Op, mesh, angular )
        call Physics    ( Op, angular, XS   )
        
    End Subroutine
    
    Subroutine Space (Op, mesh)
        Implicit None
        Class (TransportOp),      Intent (InOut) :: Op
        
        Type (MeshClass), Target, Intent (In)    :: mesh
        
        Integer                                  :: e
        Integer                                  :: Nm
        Integer                                  :: NE
        Integer                                  :: NENK
        
        if (.not. mesh%READY) then
            call stophere ('TransportOperator.f90: Space: Mesh used to build a transport operator is not READY. ' // &
                           'Make sure to initialize & post-process before use.')
        end if
        
        ! The space part of an operator will never be re-written, so if it's already built, just return.
        if (Op%SPACEREADY) return
        
        ! Give the operator its name
        Op%name = 'Transport Operator'
        
        !  =================
        !    Point to mesh  
        !  =================
        
        if (ASSOCIATED(Op%mesh)) NULLIFY(Op%mesh)
        Op%mesh => mesh
        
        !  ====================
        !    Assign constants  
        !  ====================
        
        Nm   = mesh%Nm
        NE   = mesh%NE
        NENK = mesh%NENK
        
        !  =======================================
        !    Number of local nodes in an element  
        !  =======================================
        
        if (ALLOCATED(Op%NK)) DEALLOCATE(Op%NK)
        ALLOCATE(Op%NK(NE))
        do e = 1, NE
            Op%NK(e) = mesh%el(e)%NK
        end do
        
        !  =======================================
        !    The mapping from element to FEM idx  
        !  =======================================
        
        if (ALLOCATED(Op%el2idx)) DEALLOCATE(Op%el2idx)
        ALLOCATE(Op%el2idx(NE))
        do e = 1, NE
            Op%el2idx(e) = mesh%el(e)%idx
        end do
        
        !  ========================================
        !    The mapping from element to material  
        !  ========================================
        
        if (ALLOCATED(Op%el2mat)) DEALLOCATE(Op%el2mat)
        ALLOCATE(Op%el2mat(NE))
        do e = 1, NE
            Op%el2mat(e) = mesh%el(e)%mat
        end do
        
        !  ==================
        !    Normal vectors  
        !  ==================
        
        if (ALLOCATED(Op%n)) DEALLOCATE(Op%n)
        ALLOCATE(Op%n(NE))
        do e = 1, NE
            ALLOCATE(Op%n(e)%m, source = mesh%el(e)%n) ! MHY LATER - WASTE OF SPACE... FIGURE OUT
        end do
        
        !  ===========================================
        !    Declare the space portion of T as ready  
        !  ===========================================
        
        Op%SPACEREADY = .TRUE.
        
    End Subroutine
    
    Subroutine SpaceAngle (Op, mesh, angular)
        Implicit None
        Class (TransportOp),         Intent (InOut) :: Op
        
        Type (MeshClass),    Target, Intent (In)    :: mesh
        Type (AngularClass), Target, Intent (In)    :: angular
        
        Integer                                     :: e, es, i, j
        Logical                                     :: buildSA
        Integer                                     :: Nm
        Integer                                     :: NE
        Integer                                     :: NENK
        Integer                                     :: NI
        Integer                                     :: NKe
        Integer,                     Allocatable    :: sl     (:,:)
        Integer,                     Allocatable    :: isl    (:,:)
        
        ! Initialize buildSA (for some reason this variable appears to be saved if I declare it as false above)
        buildSA = .FALSE.
        
        ! Give the operator its name
        Op%name = 'Transport Operator'
        
        !  =======================
        !    Form the space part  
        !  =======================
        if (gmemprofiling) call RSSLogger (' --- IN TRANSPORT OPERATOR --- Before building T space')
        call Space (Op, mesh) ! Is it slow to always do this?
        if (gmemprofiling) call RSSLogger (' --- IN TRANSPORT OPERATOR --- Before building T space')
        !  ====================
        !    Point to angular  
        !  ====================
        
        if (ASSOCIATED(Op%angular)) NULLIFY(Op%angular)
        Op%angular => angular
        
        !  ====================
        !    Assign constants  
        !  ====================
        
        Nm   = mesh%Nm
        NE   = mesh%NE
        NENK = mesh%NENK
        NI   = angular%NI
        
        !  ===========================================
        !    Allocate or resize the workspace vector  
        !  ===========================================
        
        if (ALLOCATED(Op%xw%v)) then
            ! If xw is already allocated, check if you need to change its size
            if (Op%xw%NENK /= NENK .or. Op%xw%NI /= NI) then
                
                Op%xw = SpaceAngleVector (mesh, angular)
                
            end if
        else
            ! If xw is not allocated, allocate it
            if (gmemprofiling) call RSSLogger (' --- IN TRANSPORT OPERATOR --- Before allocating xw')
            Op%xw = SpaceAngleVector (mesh, angular)
            if (gmemprofiling) call RSSLogger (' --- IN TRANSPORT OPERATOR --- After allocating xw')
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
            ! If xg is not allocated, allocate it
            if (gmemprofiling) call RSSLogger (' --- IN TRANSPORT OPERATOR --- Before allocating xg')
            Op%xg = SpaceAngleVector (mesh, angular)
            if (gmemprofiling) call RSSLogger (' --- IN TRANSPORT OPERATOR --- After allocating xg')
        end if
        
        ! Determine if space-angle arrays must be built/rebuilt
        if (.not. buildSA) then
            ! First, check if space-angle arrays were ever built
            buildSA = buildSA .or. (.not. ALLOCATED(Op%sl))
            ! buildSA = buildSA .or. (.not. ALLOCATED(Op%isl))
            
            ! If they were built already, check if they need to be rebuilt
            if (.not. buildSA) then
                buildSA = buildSA .or. SIZE(Op%sl,  DIM = 1) /= NE
                buildSA = buildSA .or. SIZE(Op%sl,  DIM = 2) /= NI
                ! buildSA = buildSA .or. SIZE(Op%isl, DIM = 1) /= NE
                ! buildSA = buildSA .or. SIZE(Op%isl, DIM = 2) /= NI
            end if
        end if
        
        if (buildSA) then
            
            !  ============================
            !    Build space-angle arrays  
            !  ============================
            if (gmemprofiling) call RSSLogger (' --- IN TRANSPORT OPERATOR --- Before building spaceangle arrays')
            call BuildSpaceAngle (mesh, angular, sl, isl, Op%IAF)
            if (gmemprofiling) call RSSLogger (' --- IN TRANSPORT OPERATOR --- After building spaceangle arrays')
            
            !  ==============
            !    Sweep list  
            !  ==============
            
            if (ALLOCATED(Op%sl)) DEALLOCATE(Op%sl)
            call move_alloc (sl, Op%sl)
            
            ! Remove this if its not needed for matvec
            ! !  ======================
            ! !    Inverse sweep list  
            ! !  ======================
            ! 
            ! if (ALLOCATED(Op%isl)) DEALLOCATE(Op%isl)
            ! ALLOCATE(Op%isl(NE, NI))
            ! do i = 1, NI
            !     Op%isl(1:NE,i) = sweep(i)%invlist
            ! end do
            
        end if
        
        !  ===========================================
        !    Declare the angle portion of T as ready  
        !  ===========================================
        
        Op%ANGLEREADY = .TRUE.
        
    End Subroutine
    
    Subroutine Physics (Op, angular, XS)
        Implicit None
        Class (TransportOp),         Intent (InOut) :: Op
        
        Type (AngularClass), Target, Intent (In)    :: angular
        Type (XSType),       Target, Intent (In)    :: XS
        
        ! Give the operator its name
        Op%name = 'Transport Operator'
        
        ! We do not point to angular here because this is the physics routine.
        ! To update the angular part you must provide the mesh as well as the angular.
        
        !  ===============
        !    Point to XS  
        !  ===============
        
        if (ASSOCIATED(Op%XS)) NULLIFY(Op%XS)
        Op%XS => XS
        
        !  ==============
        !    Set solver  
        !  ==============
        
        ! MHY LATER - Currently only gonna allow sweep for solver. Must generalize, not sure how...
        !             Of course the other option is to use a solver like GMRES, but I don't know how users should be
        !             able to control this.
        !             Maybe the T solver is solely determined by discretization method
        ! IDEA - have a separate 'configure' routine for each operator, that takes its own
        ! arguments and thus doesn't have an abstract interface. Allow defaults to be usable
        ! if configure is not called. This also allows me to remove solver info from angular.
        Op%solver = 'SWEEP'
        
        !  =============================================
        !    Declare the physics portion of T as ready  
        !  =============================================
        
        Op%PHYSREADY = .TRUE.
        
    End Subroutine
    
    Subroutine Destroy (Op)
        Implicit None
        Class (TransportOp), Intent (InOut) :: Op
        
        if (ALLOCATED(Op%NK))       DEALLOCATE(Op%NK)
        if (ALLOCATED(Op%el2mat))   DEALLOCATE(Op%el2mat)
        if (ALLOCATED(Op%sl))       DEALLOCATE(Op%sl)
        ! if (ALLOCATED(Op%isl)) DEALLOCATE(Op%isl)
        if (ALLOCATED(Op%solver))   DEALLOCATE(Op%solver)
        if (ALLOCATED(Op%n))        DEALLOCATE(Op%n)
        call Op%IAF%Destroy ()
        call Op%xw%Destroy  ()
        call Op%xg%Destroy  ()
        
    End Subroutine
    
End Module TransportOperator