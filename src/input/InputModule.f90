Module InputModule
    use constants
    Implicit None
    
    !  =====================================================================
    !    The input type is used to read and retain all input parameters.    
    !    It is split into blocks corresponding to the input file's format.  
    !    Furthermore, the input blocks each come with their own procedures  
    !    for generating objects based on user inputs.                       
    !                                                                       
    !    See $LIONBOLT/src/main/Main.f90 for how this ends up being used.   
    !  =====================================================================
    
    Integer, Parameter :: inunit = 1
    
    Type ProblemInput
        Logical :: General
        Logical :: Slab
    Contains
        Procedure :: Read => ProblemInputRead
    End Type
    
    Type MeshInput
        Logical                        :: dotrans   = .FALSE.
        Logical                        :: doscale   = .FALSE.
        Real (KREAL)                   :: trans (3) = [ ZERO, ZERO, ZERO ]
        Real (KREAL)                   :: scale (3) = [ ONE , ONE , ONE  ]
        Character (LEN=:), Allocatable :: fname
        Type (SlabInput),  Allocatable :: slabs (:)
    Contains
        Procedure :: Read  => MeshInputRead
        Procedure :: Build => MeshInputBuild
    End Type
    
    Type SlabInput
        Integer                        :: material
        Integer                        :: elements
        Real (KREAL)                   :: thickness
        Character (LEN=:), Allocatable :: structure
    Contains
        Procedure :: Read => SlabInputRead
    End Type
    
    Type MaterialInput ! ALSO DEFINED IN NittanyPhysics peripheral routines. There will be no conflict if this alone is changed, but consistency is better
        Logical                        :: norm294formula ! Will treat Zeff stuff later
        Integer,           Allocatable :: Z (:)
        Integer,           Allocatable :: N (:)
        Real (KREAL)                   :: Dens = ZERO
        Real (KREAL)                   :: Zeff           ! Will treat Zeff stuff later
        Character (LEN=:), Allocatable :: name
    Contains
        Procedure :: Read  => MaterialInputRead
    End Type
    
    Type MaterialsInput
        Type (MaterialInput), Allocatable :: v (:)
    Contains
        Procedure :: Read  => MaterialsInputRead
        Procedure :: Build => MaterialsInputBuild
    End Type
    
    ! LATER :
    ! SPECIFY COARSE STRUCTURE HERE.
    ! ALSO PHYSICS HERE. LIKE MGCS, NITTANYPHYSICS/CUSTOMPHYSICS, SOME OF THE OPTIONS, ETC.
    Type ParticleInput
        Logical                        :: noscat  = .FALSE.
        Logical                        :: slab            ! Just a copy of the slab logical. Set manually in ReadInput.f90. All particles will have the same value
        Integer                        :: L               ! Just a copy of the angular discretization parameter. Set manually in ReadInput.f90. All particles will have the same value
        Integer                        :: coarse          ! Number of coarse sub-systems for preconditioning (WIP)
        Integer,           Allocatable :: LC       (:)    ! Legendre orders used in coarse solves (WIP)
        Integer,           Allocatable :: angularC (:)    ! Angular discretization in the respective coarse solves (WIP)
        ! In input, ask for coarse systems to be specified like:
        ! COARSE
        !   SN 8
        !   PN 3
        ! END
        Integer                        :: G    = -1     ! User-specified number of energy groups (if provided)
        Real (KREAL)                   :: Emin = - ONE  ! User-specified minimum energy (if provided)
        Real (KREAL)                   :: Emax = - ONE  ! User-specified maximum energy (if provided)
        Character (LEN=:), Allocatable :: name          ! User-specified name of the particle
        Character (LEN=:), Allocatable :: identity      ! Actual physical identity of the particle (required if generating physics via NittanyPhysics)
        Character (LEN=:), Allocatable :: gridstruct    ! Energy group structure (if provided)
        Character (LEN=:), Allocatable :: fname
        Character (LEN=:), Allocatable :: angular       ! Angular discretization method for this particle
        Character (LEN=:), Allocatable :: solver        ! Solver for this particle
    Contains
        Procedure :: Read         => ParticleInputRead
        Procedure :: BuildAngular => ParticleInputBuildAngular
        Procedure :: BuildEnergy  => ParticleInputBuildEnergy
    End Type
    
    Type ParticlesInput
        Type (ParticleInput), Allocatable :: v (:)
    Contains
        Procedure :: Read => ParticlesInputRead
        ! You never build all particles at once, so let the individual particles
        ! build
        ! MHY LATER - wait why did I do it like this? Could I just make input carry an allocatable array???
    End Type
    
    Type BeamInput
        Logical                        :: slab                                        ! Just a copy of the slab logical. Set depending on the number of inputs for axis
        Logical                        :: used = .FALSE.                              ! Just a logical to help me keep track of whether a beam has already been used to source particles
        Integer                        :: particle                                    ! Index of particle (according to user-definition) which this beam sources
        Integer,           Allocatable :: mcgroups (:)                                ! List of groups to populate (starting from the highest energies, g = 1) in a MONOCHROMATIC spectrum
        Real (KREAL)                   :: E0     = -ONE                               ! Forced width of a single group in a MONOCHROMATIC spectrum.
        Real (KREAL)                   :: dE     = -ONE                               ! Forced width of a single group in a MONOCHROMATIC spectrum.
        Real (KREAL)                   :: weight = ONE                                ! Relative weight of a given beam.
        Real (KREAL)                   :: axis         (3) = [ZERO, ZERO, ZERO]       ! Beam axis vector (direction of inflow as well). In the slab case, only the first component is used, and it corresponds to the beam angle IN DEGREES.
        Real (KREAL)                   :: origin       (3) = [ZERO, ZERO, ZERO]       ! Beam origin
        Real (KREAL)                   :: cutoutparams (4) = [ZERO, ZERO, ZERO, ZERO] ! Cutout parameters (see examples/annotated_input.in)
        Character (LEN=:), Allocatable :: fname                                       ! Linac spectrum file name
        Character (LEN=:), Allocatable :: angdist                                     ! Beam angular
        Character (LEN=:), Allocatable :: cutout                                      ! Beam cutout
        Character (LEN=:), Allocatable :: spectrum                                    ! Spectrum type
    Contains
        Procedure :: Build => BeamInputBuild ! Strictly speaking this builds a FieldGeometry... naming conventions are a bit inconsistent since I split the original BeamType into Beam and FieldGeometry
    End Type
    
    Type BeamsInput
        Type (BeamInput), Allocatable :: v (:)
    Contains
        Procedure :: Read => BeamsInputRead
        ! You never build all beams at once, so let the individual beams build
        ! Furthermore, because beams are the only repeated header/subheader
        ! that is NOT at all dependent on ordering, the read routine is structured
        ! somewhat uniquely.
    End Type
    
    Type PostProcInput
        Logical :: fluence   = .FALSE.
        Logical :: outputunc = .FALSE.
        Logical :: energy    = .FALSE.
        Logical :: dose      = .FALSE.
        Logical :: charge    = .FALSE.
    Contains
        Procedure :: Read => PostProcInputRead
    End Type
    
    Type OptionsInput
        Logical :: storlim       = .FALSE.
        Logical :: memlim        = .FALSE.
        Logical :: debug         = .FALSE.
        Logical :: NoFCS         = .FALSE.
    Contains
        Procedure :: AssignGlobals
    End Type
    
    Type InputFile
        Integer               :: L
        Type (ProblemInput)   :: problem
        Type (MeshInput)      :: mesh
        Type (MaterialsInput) :: materials
        Type (ParticlesInput) :: particles
        Type (BeamsInput)     :: beams
        Type (PostProcInput)  :: postproc
        Type (OptionsInput)   :: options
    Contains
        Procedure :: ReadPNScattering ! Rather special treatment for a rather poorly conceptualized input block...
        Procedure :: GenerateXSLibrary
        Procedure :: Validate => ValidateInput
    End Type
    
    Interface InputFile
        Module Procedure ReadInput
    End Interface
    
    Interface
        
        ! #########################################################################################
        !   Implementation: ReadInput.f90
        
        Module Function ReadInput (fname) Result (self)
            use IO
            Implicit None
            Character (*),   Optional, Intent (In) :: fname
            Type (InputFile)                       :: self
        End Function
        
        ! #########################################################################################
        !   Implementations: InputReadModule.f90
        
        Module Subroutine ProblemInputRead (self)
            use IO
            Implicit None
            Class (ProblemInput), Intent (Out) :: self
        End Subroutine
        
        Module Subroutine MeshInputRead (self)
            use constants
            use IO
            Implicit None
            Class (MeshInput), Intent (Out) :: self
        End Subroutine
        
        Module Subroutine SlabInputRead (self)
            use constants
            use IO
            Implicit None
            Class (SlabInput), Intent (Out) :: self
        End Subroutine
        
        Module Subroutine MaterialsInputRead (self)
            use constants
            use IO
            Implicit None
            Class (MaterialsInput), Intent (Out) :: self
        End Subroutine
        
        Module Subroutine MaterialInputRead (self)
            use constants
            use IO
            Implicit None
            Class (MaterialInput), Intent (Out) :: self
        End Subroutine
        
        Module Subroutine ReadPNScattering (self)
            use constants
            use IO
            Implicit None
            Class (InputFile), Intent (InOut) :: self
        End Subroutine
        
        Module Subroutine ParticlesInputRead (self)
            use constants
            use IO
            Implicit None
            Class (ParticlesInput), Intent (Out) :: self
        End Subroutine
        
        Module Subroutine ParticleInputRead (self)
            use constants
            use IO
            Implicit None
            Class (ParticleInput), Intent (Out) :: self
        End Subroutine
        
        Module Subroutine BeamsInputRead (self)
            use constants
            use IO
            Implicit None
            Class (BeamsInput), Intent (InOut) :: self
        End Subroutine
        
        Module Subroutine PostProcInputRead (self)
            use constants
            use IO
            Implicit None
            Class (PostProcInput), Intent (InOut) :: self
        End Subroutine
        
        ! #########################################################################################
        !   Implementation: ValidateInput.f90
        
        Module Subroutine ValidateInput (self)
            Implicit None
            Class (InputFile), Intent (In) :: self
        End Subroutine
        
        ! #########################################################################################
        !   Implementations: InputBuildModule.f90
        
        Module Subroutine AssignGlobals (self)
            use globals
            Implicit None
            Class (OptionsInput), Intent (In) :: self
        End Subroutine
        
        Module Function MeshInputBuild (self) Result (mesh)
            use constants
            use LionboltAPI
            Implicit None
            Class (MeshInput), Intent (In) :: self
            Type (MeshClass)               :: mesh
        End Function
        
        Module Function MaterialsInputBuild (self) Result (mats)
            use NittanyAPI
            Implicit None
            Class (MaterialsInput), Intent (In) :: self
            Type (MaterialSet)                  :: mats
        End Function
        
        Module Function ParticleInputBuildAngular (self) Result (angular)
            use constants
            use LionboltAPI
            Implicit None
            Class (ParticleInput), Intent (In) :: self
            Type (AngularClass)                :: angular
        End Function
        
        Module Function ParticleInputBuildEnergy (self, beams) Result (energy)
            use constants
            use NittanyAPI
            use LionboltAPI
            Implicit None
            Class (ParticleInput),       Intent (In) :: self
            Type (BeamInput),  Optional, Intent (In) :: beams (:)
            Type (EnergyGrid)                        :: energy
        End Function
        
        Module Function BeamInputBuild (self) Result (fldgeo)
            use constants
            use LionboltAPI
            Implicit None
            Class (BeamInput),   Intent (In)  :: self
            Type (FieldGeometry)              :: fldgeo
        End Function
        
        Module Function GenerateXSLibrary (self, energy) Result (XSLib)
            use NittanyAPI
            use LionboltAPI
            Implicit None
            Class (InputFile), Intent (In) :: self
            Type (EnergyGrid), Intent (In) :: energy (:)
            Type (XSLibrary),  Allocatable :: XSLib  (:,:)
        End Function
        
        ! #########################################################################################
        !   Implementation: NittanyPhysicsLibrary.f90
        
        Module Subroutine NittanyPhysicsLibrary (L, mats, particles, energy, XSLib)
            use constants
            use NittanyAPI
            use LionboltAPI
            Implicit None
            Integer,                           Intent (In)  :: L
            Type (MaterialSet),                Intent (In)  :: mats
            Type (ParticleInput),              Intent (In)  :: particles (:)
            Type (EnergyGrid),                 Intent (In)  :: energy    (:)
            Type (XSLibrary),     Allocatable, Intent (Out) :: XSLib     (:,:)
        End Subroutine
        
        !   Implementation: CustomPhysicsLibrary.f90
        Module Subroutine CustomPhysicsLibrary ()
            Implicit None
        End Subroutine
        
    End Interface
    
End Module InputModule