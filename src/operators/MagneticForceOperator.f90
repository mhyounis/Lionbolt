!  ===============================================
!    THIS OPERATOR IS CURRENTLY WORK-IN-PROGRESS  
!  ===============================================
!    Plan is to allow transport in a magnetic field. Useful for MR-guided therapy
!    A prelude to more generalized Vlasov transport capabilities. Vlasov, however,
!    includes an electric field as well as nonlinear transport due to the EM field
!    being sourced by the motion of particles.
!    Even for just the inclusion of an E-field, constant or not, there will need
!    to be a lot more writing as well as an entirely new driver file, because
!    electric fields do work, and thus, can give energy to a particle inbetween
!    scattering events, meaning we must allow for upscatter and therefore we 
!    can not use energy iteration.

Module MagneticForceOperator
    use constants
    use NittanyAPI
    use Geometry
    use AngularSpace
    use LinearOperatorClass
    use MagneticForceAlgebra
    Implicit None
    Private
    ! Public :: MagneticForceOp ! This will be un-commented when this operator is ready for use
    
    ! Operator: F (Velocity-space advection)
    Type, Extends (LinearOperator) :: MagneticForceOp
        Logical                   :: luse = .FALSE.               ! Logical that tells whether or not a Boltzmann operator will actually use this term
        Real (KREAL)              :: invp                         ! Energy-grouped form of 1/|p|, where p is the particle momentum
        Real (KREAL)              :: B (3) = [ ZERO, ZERO, ZERO ] ! Magnetic field vector ! Currently only constant
        Real (KREAL), Allocatable :: ddmu  (:,:)                  ! The angular matrix which applies a derivative with respect to polar angle cosine mu
        Real (KREAL), Allocatable :: ddphi (:,:)                  ! The angular matrix which applies a derivative with respect to azimuth phi
    Contains
        Procedure :: MatVec     => MatVec
        Procedure :: MatInv     => MatInv
        Procedure :: Full       => Full
        Procedure :: Space      => Space
        Procedure :: SpaceAngle => SpaceAngle
        Procedure :: Physics    => Physics
        Procedure :: Destroy    => Destroy
        Procedure :: BuildFields
        Generic   :: Build      => Full, Space, SpaceAngle, Physics
    End Type
    
    Contains
    
    Subroutine MatVec (Op, x)
        use constants
        use LinearOperatorClass
        Class (MagneticForceOp), Intent (InOut) :: Op
        
        Type (SpaceAngleVector), Intent (InOut) :: x
        
        if (.not. (Op%SPACEREADY .AND. Op%ANGLEREADY .AND. Op%PHYSREADY)) &
            call stophere ('MagneticForceOp%MatVec was called, but the operator was not ready.')
        
        
        
    End Subroutine
    
    Subroutine MatInv (Op, s)
        use constants
        use LinearOperatorClass
        Class (MagneticForceOp), Intent (InOut) :: Op
        
        Type (SpaceAngleVector), Intent (InOut) :: s
        
        if (.not. (Op%SPACEREADY .AND. Op%ANGLEREADY .AND. Op%PHYSREADY)) &
            call stophere ('MagneticForceOp%MatInv was called, but the operator was not ready.')
        
        
        
    End Subroutine
    
    Subroutine Full (Op, mesh, angular, XS)
        Implicit None
        Class (MagneticForceOp),     Intent (InOut) :: Op
        
        Type (MeshClass),    Target, Intent (In)    :: mesh
        Type (AngularClass), Target, Intent (In)    :: angular
        Type (XSType),       Target, Intent (In)    :: XS
        
        ! Give the operator its name
        Op%name = 'Magnetic Force Operator'
        
        ! Now construct the dimensional parts of F
        call Space      ( Op, mesh          )
        call SpaceAngle ( Op, mesh, angular )
        call Physics    ( Op, angular, XS   )
        
    End Subroutine
    
    Subroutine Space (Op, mesh)
        Implicit None
        Class (MagneticForceOp),  Intent (InOut) :: Op
        
        Type (MeshClass), Target, Intent (In)    :: mesh
        
        if (.not. mesh%READY) then
            call stophere ('MagneticForceOp.f90: Space: Mesh used to build a force operator is not READY. ' // &
                           'Make sure to initialize & post-process before use.')
        end if
        
        ! Give the operator its name
        Op%name = 'Magnetic Force Operator'
        
        Op%SPACEREADY = .TRUE.
        
    End Subroutine
    
    Subroutine SpaceAngle (Op, mesh, angular)
        Implicit None
        Class (MagneticForceOp),     Intent (InOut) :: Op
        
        Type (MeshClass),    Target, Intent (In)    :: mesh
        Type (AngularClass), Target, Intent (In)    :: angular
        
        Integer                                     :: NENK
        Integer                                     :: NI
        Real (KREAL),                Allocatable    :: P (:,:) ! Associated Legendre polynomials
        Real (KREAL),                Allocatable    :: y (:,:,:) ! Real spherical harmonics
        
        NENK = mesh%NENK
        NI   = angular%NI
        
        ! Give the operator its name
        Op%name = 'Magnetic Force Operator'
        
        !  ================================================
        !    Form the space part if it hasn't been formed  
        !  ================================================
        
        if (.not. Op%SPACEREADY) then
            call Space (Op, mesh)
        end if
        
        ALLOCATE(Op%ddmu (NI, NI))
        ALLOCATE(Op%ddphi(NI, NI))
        
        
        
        Op%ANGLEREADY = .TRUE.
        
    End Subroutine
    
    Subroutine Physics (Op, angular, XS)
        Implicit None
        Class (MagneticForceOp),     Intent (InOut) :: Op
        
        Type (AngularClass), Target, Intent (In)    :: angular
        Type (XSType),       Target, Intent (In)    :: XS
            
        ! Give the operator its name
        Op%name = 'Magnetic Force Operator'
        
        Op%PHYSREADY = .TRUE.
        
    End Subroutine
    
    Subroutine Destroy (Op)
        Implicit None
        Class (MagneticForceOp), Intent (InOut) :: Op
        
    End Subroutine
    
    Subroutine BuildFields (Op, B)
        Class (MagneticForceOp), Intent (InOut) :: Op
        
        Real (KREAL),    Intent (In)    :: B (3)
        
        Op%B = B
        
    End Subroutine
    
End Module