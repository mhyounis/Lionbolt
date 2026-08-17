Module LinearOperatorClass
    ! =================================================================================================== 
    !   Emphasize that these operators are specifically space-angle operators
    !   
    !   
    !   
    !   
    !   
    !   
    !   
    !   
    !   
    !   
    !   
    !   
    ! =================================================================================================== 
    use constants
    use Geometry
    use AngularSpace
    use NittanyAPI
    Implicit None
    
    ! HUGE NOTE : Keeping track of allocation status of these objects is going to be a pain for users, since
    !             so many routines demand pre-allocated vectors. Thus, enforce an initialization routine like psi = SpaceAngleVector(NENK, NI)
    !             Allow for deallocation but users will have to reallocate it.
    !             Maybe also have a type-bound routine that accepts mesh and angular and then confirms that the vector belongs to them.
    !             It would be called by internal routines, rather than the user.
    
    ! SHOULD I MAKE THIS A PARAMETRIZED DERIVED TYPE? This way NENK and NI are length parameters.
    Type :: SpaceAngleVector
        Integer                        :: NENK
        Integer                        :: NI
        Real (KREAL),      Allocatable :: v (:)
        ! Character (LEN=:), Allocatable :: disc ! Not for users. 'SN' or 'PN' depending on v, not necessarily depending on the angular used to build.
    Contains
        Procedure :: Destroy
        Procedure :: Fluence
        ! Procedure :: MapToPN
        ! Procedure :: MapToSN
        ! Procedure :: AngularDist ! Returns the angular distribution
    End Type
    
    Type, Abstract :: LinearOperator
        Logical                          :: SPACEREADY = .FALSE.
        Logical                          :: ANGLEREADY = .FALSE.
        Logical                          :: PHYSREADY  = .FALSE.
        Character (LEN=:),   Allocatable :: name
        Type (MeshClass),    Pointer     :: mesh
        Type (AngularClass), Pointer     :: angular
        Type (XSType),       Pointer     :: XS
    Contains
        Procedure (MatVecInterface),          Deferred :: MatVec
        Procedure (MatInvInterface),          Deferred :: MatInv
        Procedure (BuildFullInterface),       Deferred :: Full
        Procedure (BuildSpaceInterface),      Deferred :: Space
        Procedure (BuildSpaceAngleInterface), Deferred :: SpaceAngle
        Procedure (BuildPhysicsInterface),    Deferred :: Physics
        Procedure (DestroyInterface),         Deferred :: Destroy
    End Type
    
    Abstract Interface
        
        ! MHY LATER - could *possibly* speed things up by figuring out a way to make Op not InOut.
        ! Also in general would be great to get rid of the workspace and guess vectors somehow.
        ! Storage wise.
        
        Subroutine MatVecInterface (Op, x)
            Import :: KREAL, SpaceAngleVector, LinearOperator
            Class (LinearOperator),  Intent (InOut) :: Op
            Type (SpaceAngleVector), Intent (InOut) :: x
        End Subroutine
        
        Subroutine MatInvInterface (Op, s)
            Import :: KREAL, SpaceAngleVector, LinearOperator
            Class (LinearOperator),  Intent (InOut) :: Op
            Type (SpaceAngleVector), Intent (InOut) :: s
        End Subroutine
        
        Subroutine BuildFullInterface (Op, mesh, angular, XS)
            Import :: LinearOperator, MeshClass, AngularClass, XSType
            Class (LinearOperator),         Intent (InOut) :: Op
            Type (MeshClass),       Target, Intent (In)    :: mesh
            Type (AngularClass),    Target, Intent (In)    :: angular
            Type (XSType),          Target, Intent (In)    :: XS
        End Subroutine
        
        Subroutine BuildSpaceInterface (Op, mesh)
            Import :: LinearOperator, MeshClass
            Class (LinearOperator),         Intent (InOut) :: Op
            Type (MeshClass),       Target, Intent (In)    :: mesh
        End Subroutine
        
        Subroutine BuildSpaceAngleInterface (Op, mesh, angular)
            Import :: LinearOperator, MeshClass, AngularClass
            Class (LinearOperator),         Intent (InOut) :: Op
            Type (MeshClass),       Target, Intent (In)    :: mesh
            Type (AngularClass),    Target, Intent (In)    :: angular
        End Subroutine
        
        Subroutine BuildPhysicsInterface (Op, angular, XS)
            Import :: LinearOperator, AngularClass, XSType
            Class (LinearOperator),         Intent (InOut) :: Op
            Type (AngularClass),    Target, Intent (In)    :: angular
            Type (XSType),          Target, Intent (In)    :: XS
        End Subroutine
        
        Subroutine DestroyInterface (Op)
            Import :: LinearOperator
            Class (LinearOperator), Intent (InOut) :: Op
        End Subroutine
        
    End Interface
    
    Interface SpaceAngleVector
        Module Procedure init_SpaceAngleVector
    End Interface
    
    Private :: init_SpaceAngleVector, Destroy
    
Contains
    
    Function init_SpaceAngleVector (mesh, angular) Result (x)
        ! Should I make x point to mesh and angular?
        ! If not should I just make this routine take in NENK and NI?
        Implicit None
        Type (MeshClass),    Intent (In) :: mesh
        Type (AngularClass), Intent (In) :: angular
        
        Type (SpaceAngleVector)          :: x
        
        ! Initialize the vector, IF NEEDED.
        if (ALLOCATED(x%v)) then
            if (x%NENK == mesh%NENK .and. x%NI == angular%NI) return
            
            DEALLOCATE(x%v)
        end if
        
        x%NENK = mesh%NENK
        x%NI   = angular%NI
        ALLOCATE(x%v(mesh%NENK * angular%NI))
        ! x%disc = angular%disc
        
    End Function
    
    Subroutine Destroy (self)
        Class (SpaceAngleVector), Intent (InOut) :: self
        
        if (ALLOCATED(self%v)) DEALLOCATE(self%v)
        self%NENK = 0
        self%NI   = 0
        
    End Subroutine
    
    Function Fluence (self, w) Result (fl) ! This is one of the routines that might make me have SpaceAngleVector store pointers to mesh and angular... good idea or not???
        Implicit None
        Class (SpaceAngleVector), Intent (In) :: self
        Real (KREAL),             Intent (In) :: w  (:)
        
        Real (KREAL), Allocatable             :: fl (:)
        
        Integer                               :: i
        Integer                               :: start
        
        ALLOCATE(fl(self%NENK))
        
        fl = ZERO
        do i = 1, self%NI
            start = (i - 1) * self%NENK
            
            fl(1:self%NENK) = fl(1:self%NENK) &
                            + w(i) * self%v(1 + start : self%NENK + start)
        end do
        
    End Function
    
    ! Subroutine MapToPN (self)
    !     Implicit None
    !     Class (SpaceAngleVector), Intent (In) :: self
        
    !     Real (KREAL), Allocatable
        
        
    ! End Subroutine
    
End Module