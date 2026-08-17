Submodule (InputModule) submodInputBuildModule
Contains
    
    Module Subroutine AssignGlobals (self)
        use globals
        Implicit None
        Class (OptionsInput), Intent (In) :: self
        
        gstorlim    = self%storlim
        gmemlim     = self%memlim
        gdebug      = self%debug
        
    End Subroutine
    
    Module Function MeshInputBuild (self) Result (mesh)
        use constants
        use LionboltAPI
        use InputModule
        Implicit None
        Class (MeshInput), Intent (In) :: self
        
        Type (MeshClass)               :: mesh
        
        Logical                        :: slab
        Integer                        :: s
        Integer                        :: NE
        Integer                        :: mat
        Real (KREAL)                   :: T
        Character (LEN=:), Allocatable :: structure
        
        ! This is a good condition as long as :
        !       self%slabs is indeed allocated regardless in MeshInputRead
        !       self%slabs is NEVER appended to in Problem General
        slab = SIZE(self%slabs) > 0
        
        if (.not. slab) then
            
            mesh = MeshClass (self%fname)
            
            if (self%dotrans) call mesh%Translate (self%trans)
            if (self%doscale) call mesh%Scale     (self%scale)
            
        else
            
            do s = 1, SIZE(self%slabs)
                NE        = self%slabs(s)%elements
                T         = self%slabs(s)%thickness
                mat       = self%slabs(s)%material
                structure = self%slabs(s)%structure
                call mesh%AddSlab (NE, T, mat, structure)
            end do
            
        end if
        
        call mesh%PostProcess ()
        
    End Function
    
    ! Should this just be in NittanyPhysics somehow???
    ! Sort this out later...
    Module Function MaterialsInputBuild (self) Result (mats)
        use NittanyAPI
        use InputModule
        Implicit None
        Class (MaterialsInput), Intent (In) :: self
        
        Type (MaterialSet)                  :: mats
        
        Integer                             :: m
        
        do m = 1, SIZE(self%v)
            call mats%AddMaterial (self%v(m)%Z, self%v(m)%N, self%v(m)%Dens, self%v(m)%name)
        end do
        
    End Function
    
    Module Function ParticleInputBuildAngular (self) Result (angular)
        use constants
        use LionboltAPI
        use InputModule
        Implicit None
        Class (ParticleInput), Intent (In) :: self
        
        Type (AngularClass)                :: angular
        
        Integer                            :: i
        Logical                            :: DoSN
        Logical                            :: DoPN
        Character (LEN=:),     Allocatable :: solverstr
        
        DoSN = self%angular == 'SN'
        DoPN = self%angular == 'PN'
        
        solverstr = self%solver
        if (solverstr /= 'SI' .and. solverstr /= 'GMRES') then
            call stophere ('InputBuildModule.f90: ParticleInputBuildAngular: Unknown solver given to BuildAngular.')
        end if
        
        if (DoSN) then
            angular = AngularClass (self%L, self%solver, self%slab)
        else if (DoPN) then
            ! WIP of course
            angular = AngularClass ()
        end if
        
    End Function
    
    Module Function ParticleInputBuildEnergy (self, beams) Result (energy)
        use constants
        use IO
        use NittanyAPI
        use LionboltAPI
        use InputModule
        Implicit None
        Class (ParticleInput),       Intent (In) :: self
        Type (BeamInput),  Optional, Intent (In) :: beams (:)
        
        Type (EnergyGrid)                        :: energy
        
        Integer                                  :: pb
        Integer                                  :: useG
        
        useG = self%G
        if (PRESENT(beams)) then
            ! Determine the number of groups to actually construct with.
            ! That is, if a user is specifying custom monochromatic widths,
            ! we need to insert those in a background grid, in the end
            ! having the user-specified number of groups
            do pb = 1, SIZE(beams)
                if (beams(pb)%E0 <= ZERO) cycle
                useG = useG - beams(pb)%mcgroups(1)
            end do
        end if
        
        if (useG == self%G) then
            select case (self%gridstruct)
            case ('linear', 'logarithmic', 'exponential')
                energy = EnergyGrid (useG, self%Emin, self%Emax, self%gridstruct)
            case ('null')
                energy = EnergyGrid (self%fname)
            end select
        else
            if (PRESENT(beams)) then
                ! Now insert energy groups if a beam for this particle
                ! requested that
                do pb = 1, SIZE(beams)
                    if (beams(pb)%E0 <= ZERO) cycle
                    
                    call energy%InsertGroups (beams(pb)%mcgroups(1), beams(pb)%E0, beams(pb)%dE, spread=.TRUE.)
                end do
            end if
        end if
        
    End Function
    
    Module Function BeamInputBuild (self) Result (fldgeo)
        use constants
        use LionboltAPI
        use InputModule
        Implicit None
        Class (BeamInput),   Intent (In) :: self
        
        Type (FieldGeometry)             :: fldgeo
        
        if (self%slab) then
            fldgeo = FieldGeometry (self%weight, self%spectrum, self%axis(1))
        else
            fldgeo = FieldGeometry (self%weight, self%spectrum, self%origin, &
                                    self%axis, self%angdist, self%cutout, self%cutoutparams)
        end if
        
    End Function
    
    Module Function GenerateXSLibrary (self, energy) Result (XSLib)
        use NittanyAPI
        use LionboltAPI
        use InputModule
        Implicit None
        Class (InputFile), Intent (In) :: self
        Type (EnergyGrid), Intent (In) :: energy (:)
        
        Type (XSLibrary),  Allocatable :: XSLib  (:,:)
        
        Integer                        :: p
        Type (MaterialSet)             :: mats
        
        if (.TRUE.) then
            ! CURRENTLY ONLY NITTANYPHYSICS IS ALLOWED.
            
            mats = self%materials%Build ()
            
            call NittanyPhysicsLibrary (self%L, mats, self%particles%v, energy, XSLib)
            
        else
            ! CustomPhysics
            ! First discretize energy according to input
            
            ! call CustomPhysicsLibrary ()
            
        end if
        
        call PrintXSLibrary (iuout, mats, XSLib)
        
        call WriteXSLibraryLionbolt (XSLib)
        
    End Function
    
End Submodule