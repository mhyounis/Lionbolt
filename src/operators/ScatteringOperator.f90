Module ScatteringOperator
    !  =========================================================================
    !    This module implements the scattering operator as an extension of the
    !    abstract LinearOperator class.
    !    Read the docs at the OpenRPS website if you are interested in 
    !    learning more.
    !  =========================================================================
    use constants
    use NittanyAPI
    use Geometry
    use AngularSpace
    use LinearOperatorClass
    use ScatteringAlgebra
    Implicit None
    Private
    Public :: ScatteringOp
    
    ! Operator: K (Scattering)
    Type, Extends (LinearOperator) :: ScatteringOp
        Logical                    :: deltadown = .FALSE. ! This is only relevant if the operator is used for delta-down. Generally avoid tampering with it.
        Real (KREAL), Allocatable  :: XSdd   (:)          ! Delta-down cross sections (if present). NOT set by Build, set by BuildDeltaDown
        Real (KREAL), Allocatable  :: xw     (:,:)        ! Workspace vector for matmul ! (e x k, i)
        Real (KREAL), Allocatable  :: YYS    (:,:,:)      ! Single group transfer scattering kernel ! (i', i, mat)
        Type (CIntV), Allocatable  :: mat2sd (:)          ! (mat)%(list of elements having material mat)
    Contains
        Procedure :: MatVec     => MatVec
        Procedure :: MatInv     => MatInv ! The scattering operator is not invertible.
        Procedure :: Full       => Full
        Procedure :: Space      => Space
        Procedure :: SpaceAngle => SpaceAngle
        Procedure :: Physics    => Physics
        Procedure :: Destroy    => Destroy
        Procedure :: BuildDeltaDown
        Generic   :: Build      => Full, Space, SpaceAngle, Physics
    End Type
    
    Contains
    
    Subroutine MatVec (Op, x)
        use constants
        use ScatteringAlgebra, ONLY : ApplyScattering
        Implicit None
        ! WRAPPER SUBROUTINE
        ! The raw subroutine is ApplyScattering.
        ! This subroutine applies K to a vector, x, in collapsed notation.
        Class (ScatteringOp),    Intent (InOut) :: Op
        
        Type (SpaceAngleVector), Intent (InOut) :: x
        
        if (.not. (Op%SPACEREADY .AND. Op%ANGLEREADY .AND. Op%PHYSREADY)) &
            call stophere ('ScatteringOp%MatVec was called, but the operator was not ready.')
        
        if (.not. Op%deltadown) then
            call ApplyScattering (Op%mesh%Nm,    &
                                  Op%angular%NI, &
                                  Op%mesh%NENK,  &
                                  Op%YYS,        &
                                  Op%mat2sd,     &
                                  Op%xw,         &
                                  x%v             )
        else
            call ApplyDeltaDown  (Op%mesh%Nm,    &
                                  Op%angular%NI, &
                                  Op%mesh%NENK,  &
                                  Op%XSdd,       &
                                  Op%mat2sd,     &
                                  x%v             )
        end if
        
    End Subroutine
    
    Subroutine MatInv (Op, s)
        Implicit None
        Class (ScatteringOp),    Intent (InOut) :: Op
        
        Type (SpaceAngleVector), Intent (InOut) :: s
        
        call stophere ('ScatteringOp%MatInv was called, but this operation is undefined.')
        
    End Subroutine
    
    Subroutine Full (Op, mesh, angular, XS)
        Implicit None
        Class (ScatteringOp),        Intent (InOut) :: Op
        
        Type (MeshClass),    Target, Intent (In)    :: mesh
        Type (AngularClass), Target, Intent (In)    :: angular
        Type (XSType),       Target, Intent (In)    :: XS
        
        ! Give the operator its name
        Op%name = 'Scattering Operator'
        
        ! Now construct the dimensional parts of K
        call Space      ( Op, mesh          )
        call SpaceAngle ( Op, mesh, angular )
        call Physics    ( Op, angular, XS   )
        
    End Subroutine
    
    Subroutine Space (Op, mesh)
        Implicit None
        Class (ScatteringOp),     Intent (InOut) :: Op
        
        Type (MeshClass), Target, Intent (In)    :: mesh
        
        Integer                                  :: e, m, jp
        Logical,                  Allocatable    :: elmask (:)
        Integer                                  :: Nm
        Integer                                  :: NE
        Integer                                  :: NKe
        Integer                                  :: NENK
        
        if (.not. mesh%READY) then
            call stophere ('ScatteringOperator.f90: Space: Mesh used to build a scattering operator is not READY. ' // &
                           'Make sure to initialize & post-process before use.')
        end if
        
        ! The space part of an operator will never be re-written, so if it's already built, just return.
        if (Op%SPACEREADY) return
        
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
        
        ! Give the operator its name
        Op%name = 'Scattering Operator'
        
        !  ===========================================
        !    Gives the material of an element x node  ! Is there a better way to do this kind of thing than using mat2sd?
        !  ===========================================
        
        if (ALLOCATED(Op%mat2sd)) DEALLOCATE(Op%mat2sd)
        ALLOCATE(Op%mat2sd(Nm), source = mesh%mat2sd) ! Should I use move alloc?
        
        !  ===========================================
        !    Declare the space portion of K as ready  
        !  ===========================================
        
        Op%SPACEREADY = .TRUE.
        
    End Subroutine
    
    Subroutine SpaceAngle (Op, mesh, angular)
        Implicit None
        Class (ScatteringOp),        Intent (InOut) :: Op
        
        Type (MeshClass),    Target, Intent (In)    :: mesh
        Type (AngularClass), Target, Intent (In)    :: angular
            
        Integer                                     :: e, g, i, jp, m
        Logical,                     Allocatable    :: elmask (:)
        Integer                                     :: NENK
        Integer                                     :: NI
        Integer                                     :: NK
        Integer                                     :: total
        
        ! Give the operator its name
        Op%name = 'Scattering Operator'
        
        !  ================================================
        !    Form the space part if it hasn't been formed  
        !  ================================================
        
        call Space (Op, mesh)
        
        !  ====================
        !    Point to angular  
        !  ====================
        
        if (ASSOCIATED(Op%angular)) NULLIFY(Op%angular)
        Op%angular => angular
        
        !  ====================
        !    Assign constants  
        !  ====================
        
        NENK = mesh%NENK
        NI   = angular%NI
        
        !  =======================================
        !    Allocate or resize the guess vector  
        !  =======================================
        
        if (ALLOCATED(Op%xw)) then
            ! If xw is already allocated, check if you need to change its size
            if (SIZE(Op%xw, DIM=1) /= NENK .AND. SIZE(Op%xw, DIM=2) /= NI) then
                DEALLOCATE(Op%xw)
                
                ALLOCATE(Op%xw(NENK, NI))
            end if
        else
            ! If xw is not allocated, allocate it
            ALLOCATE(Op%xw(NENK, NI))
        end if
        
        !  ===========================================
        !    Declare the angle portion of K as ready  
        !  ===========================================
        
        Op%ANGLEREADY = .TRUE.
        
    End Subroutine
    
    Subroutine Physics (Op, angular, XS)
        Implicit None
        Class (ScatteringOp),        Intent (InOut) :: Op
        
        Type (AngularClass), Target, Intent (In)    :: angular
        Type (XSType),       Target, Intent (In)    :: XS
        
        Integer                                     :: ell, i, m
        Integer                                     :: Nm
        Integer                                     :: NI
        
        ! Give the operator its name
        Op%name = 'Scattering Operator'
        
        ! Ensure that the scattering operator doesn't use delta-down logic
        ! (for a delta-down scattering operator use Op%BuildDeltaDown procedure)
        Op%deltadown = .FALSE.
        
        ! We do not point to angular here because this is the physics routine.
        ! To update the angular part you must provide the mesh as well as the angular.
        
        !  ===============
        !    Point to XS  
        !  ===============
        
        if (ASSOCIATED(Op%XS)) NULLIFY(Op%XS)
        Op%XS => XS
        
        !  ====================
        !    Assign constants  
        !  ====================
        
        Nm = SIZE(XS%s, dim=2)
        NI = angular%NI
        
        !  =======================
        !    Populate the kernel  
        !  =======================
        
        if (ALLOCATED(Op%YYS)) then
            if (SIZE(Op%YYS, DIM=1) /= NI) then
                DEALLOCATE(Op%YYS)
                ALLOCATE(Op%YYS(NI, NI, Nm))
            end if
        else
            ALLOCATE(Op%YYS(NI, NI, Nm))
        end if
        
        Op%YYS = ZERO
        do m = 1, Nm
            do ell = 0, angular%L
                Op%YYS(1:NI,1:NI,m) = Op%YYS(1:NI,1:NI,m) &
                                    + angular%YY(1:NI,1:NI,ell) * XS%s(ell,m)
            end do
        end do
        
        !  =============================================
        !    Declare the physics portion of K as ready  
        !  =============================================
        
        Op%PHYSREADY = .TRUE.
        
    End Subroutine
    
    Subroutine Destroy (Op)
        Implicit None
        Class (ScatteringOp), Intent (InOut) :: Op
        
        Op%deltadown = .FALSE.
        if (ALLOCATED(Op%XSdd))   DEALLOCATE(Op%XSdd)
        if (ALLOCATED(Op%xw))     DEALLOCATE(Op%xw)
        if (ALLOCATED(Op%YYS))    DEALLOCATE(Op%YYS)
        if (ALLOCATED(Op%mat2sd)) DEALLOCATE(Op%mat2sd)
        
    End Subroutine
    
    Subroutine BuildDeltaDown (Op, XSdd)
        ! This routine builds a delta-down scattering kernel specifically
        ! Perhaps I can include this as part of the generic interface Build. Would that be more user-friendly?
        Implicit None
        Class (ScatteringOp), Intent (InOut) :: Op
        
        Real (KREAL),         Intent (In)    :: XSdd (:)
        
        ! Give the operator its name
        Op%name = 'Scattering Operator'
        
        Op%deltadown = .TRUE.
        
        !  =======================
        !    Populate the kernel  
        !  =======================
        
        ! Deallocating and allocating this size is not enough of a burden for the below to be inefficient
        if (ALLOCATED(Op%XSdd)) DEALLOCATE(Op%XSdd)
        ALLOCATE(Op%XSdd, source=XSdd)
        
        !  =============================================
        !    Declare the physics portion of K as ready  
        !  =============================================
        
        Op%PHYSREADY = .TRUE.
        
    End Subroutine
    
End Module ScatteringOperator