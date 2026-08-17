Submodule (InputModule) submodInputReadModule
Contains
    
    Module Subroutine ProblemInputRead (self)
        use IO
        use InputModule
        Implicit None
        Class (ProblemInput), Intent (Out) :: self
        
        Integer                            :: ierr
        Character (LEN=:),    Allocatable  :: line
        Type (WordsType),     Allocatable  :: words (:)
        
        ! Read first line
        1 continue
        call ReadSplitSimplify (inunit, line, words, ierr, ic=.TRUE.)
        if (ierr < 0) then
            call stophere ('Input read error: Blank input file.')
        end if
        
        ! Allow for comments or line skips
        if (line == '') go to 1
        
        if (SIZE(words) /= 2) then
            call stophere ('Input read error: Either too many or too few words in first input line. ' // &
                           'Use two words - "PROBLEM" followed by the desired problem type. Your input was: ' // line)
        end if
        
        if (ierr < 0) then
            call stophere ('Input read error: First line has too few words. Use key "PROBLEM" followed by ' // &
                           'desired problem type. Your input was: ' // line)
        end if
        
        if (words(1)%str /= 'problem') then
            call stophere ('Input read error: First line must denote solve type. Use key "PROBLEM". '// &
                           'Key used was: ' // words(1)%str)
        end if
        
        self%General = .FALSE.
        self%Slab    = .FALSE.
        
        select case (words(2)%str)
        case ('general'); self%General = .TRUE.
        case ('slab');    self%Slab = .TRUE.
        case default;     call stophere ('Input read error: Invalid problem type : ' // line)
        end select
        
    End Subroutine
    
    Module Subroutine MeshInputRead (self)
        use constants
        use IO
        use InputModule
        Implicit None
        Class (MeshInput), Intent (Out) :: self
        
        Integer                        :: ierr
        Character (LEN=:), Allocatable :: line
        Type (SlabInput)               :: blank
        Type (WordsType),  Allocatable :: words (:)
        
        ! Initialize the slab subheaders
        ALLOCATE(self%slabs(0))
        
        do
            call ReadSplitSimplify (inunit, line, words, ierr, ic=.TRUE.)
            
            select case (words(1)%str)
            case ('gmsh')
                
                if (ierr < 0) then
                    call stophere ('Input read error: Missing file name for MESH GMSH key.')
                end if
                
                ! Read file name in the true case (not lower case)
                call GoBack (inunit, 1)
                call ReadAndSplit (inunit, line, words, ierr, ic=.TRUE.)
                
                ! Use full line after 'gmsh'
                self%fname = TRIM(ADJUSTL(line(5:)))
                
            case ('translate')
                if (ierr < 0) then
                    call stophere ('Input read error: Missing parameters for MESH Translate key.')
                end if
                
                read(words(2)%str,*,iostat=ierr) self%trans(1)
                read(words(3)%str,*,iostat=ierr) self%trans(2)
                read(words(4)%str,*,iostat=ierr) self%trans(3)
                
                self%dotrans = .TRUE.
                
                if (ierr < 0) then
                    call stophere ('Input read error: Unable to read MESH translation vector. Make sure this ' // &
                                   'is given as three separate real numbers.')
                end if
                
            case ('scale')
                if (ierr < 0) then
                    call stophere ('Input read error: Missing parameters for MESH Scale key.')
                end if
                
                read(words(2)%str,*,iostat=ierr) self%scale(1)
                read(words(3)%str,*,iostat=ierr) self%scale(2)
                read(words(4)%str,*,iostat=ierr) self%scale(3)
                
                self%doscale = .TRUE.
                
                if (ierr < 0) then
                    call stophere ('Input read error: Unable to read MESH scale vector. Make sure this ' // &
                                   'is given as three separate real numbers.')
                end if
                
            case ('slab')
                ! DO I NEED TO DESTROY THE OLD BLANK?
                call blank%Read ()
                
                self%slabs = [self%slabs, blank]
                
            case ('')
                cycle
                
            case ('end')
                return
                
            case default
                ! Keep? This prevents mesh from trying to read other headers when end is not present.
                call stophere ('Input read error: Unknown key in MESH block : ' // words(1)%str)
                
            end select
        end do
        
    End Subroutine
    
    Module Subroutine SlabInputRead (self)
        use constants
        use IO
        use InputModule
        Implicit None
        Class (SlabInput), Intent (Out) :: self
        
        Integer                         :: ierr
        Logical                         :: matgiven = .FALSE.
        Integer                         :: NE
        Integer                         :: mat
        Real (KREAL)                    :: T
        Character (LEN=:), Allocatable  :: line
        Type (WordsType),  Allocatable  :: words (:)
        
        do
            call ReadSplitSimplify (inunit, line, words, ierr, ic=.TRUE.)
            
            select case (words(1)%str)
            case ('structure')
                if (ierr < 0) then
                    call stophere ('Input read error: Missing value for SLAB Structure key.')
                end if
                
                select case (words(2)%str)
                case ('linear')
                    self%structure = 'linear'
                case ('logarithmic')
                    self%structure = 'logarithmic'
                case default
                    call stophere ('Input read error: Unknown structure for SLAB Structure key:' // words(2)%str)
                end select
                
            case ('thickness')
                if (ierr < 0) then
                    call stophere ('Input read error: Missing value for SlSLABaSLABb Thickness key.')
                end if
                
                read(words(2)%str,*) T
                
                self%thickness = T
                
            case ('elements')
                if (ierr < 0) then
                    call stophere ('Input read error: Missing value for SLAB Elements key.')
                end if
                
                read(words(2)%str,*) NE
                
                self%elements = NE
                
            case ('material')
                if (ierr < 0) then
                    call stophere ('Input read error: Missing value for SLAB Material key.')
                end if
                
                matgiven = .TRUE.
                
                read(words(2)%str,*) mat
                
                self%material = mat
                
            case ('')
                print *, 'CYCLING'
                cycle
                
            case ('end')
                if (.not. matgiven) self%material = 1
                
                return
            end select
        end do
        
    End Subroutine
    
    Module Subroutine MaterialsInputRead (self)
        use constants
        use IO
        use InputModule
        Implicit None
        Class (MaterialsInput), Intent (Out) :: self
        
        Integer                              :: ierr
        Logical                              :: namegiven = .FALSE.
        Character (LEN=:),      Allocatable  :: line
        Type (WordsType),       Allocatable  :: words (:)
        Type (MaterialInput)                 :: blank
        
        ALLOCATE(self%v(0))
        
        do
            call ReadSplitSimplify (inunit, line, words, ierr, ic=.TRUE.)
            
            select case (words(1)%str)
            case default
                ! It is assumed the first line is the name of the material.
                
                ! Read material name in the true case (not lower case)
                call GoBack (inunit, 1)
                call ReadAndSplit (inunit, line, words, ierr, ic=.TRUE.)
                
                call blank%Read ()
                
                blank%name = line
                
                self%v = [self%v, blank]
                
            case ('')
                cycle
                
            case ('end')
                ! Wrap up and return
                if (SIZE(self%v) == 0) then
                    call stophere ('Input read error: No materials listed in MATERIALS block.')
                end if
                
                return
                
            end select
        end do
    
    End Subroutine
    
    Module Subroutine MaterialInputRead (self)
        use constants
        use IO
        use InputModule
        Implicit None
        Class (MaterialInput), Intent (Out) :: self
        
        Logical                             :: norm294formula
        Integer                             :: ierr
        Integer                             :: Z
        Integer                             :: N
        Character (LEN=:),     Allocatable  :: line
        Type (WordsType),      Allocatable  :: words (:)
        Interface
            Function PeriodicTable (name) Result (Z)
                Character (*), Intent (In) :: name
                Integer                    :: Z
            End Function PeriodicTable
        End Interface
        
        ALLOCATE(self%Z(0))
        ALLOCATE(self%N(0))
        
        do
            call ReadSplitSimplify (inunit, line, words, ierr, ic=.TRUE.)
            
            select case (words(1)%str)
            case ('density')
                read(words(2)%str,*) self%Dens
                
            ! case ('zeff')
            !     
            !     if (word2 == 'formula') then
            !         norm294formula = .TRUE.
            !     else
            !         read(word2,*) Zeff
            !     end if
            !     
            case default
                
                ! Going to try reading an atom. First, go back and read in the true case (not lowercase)
                ! Should I make it case-insensitive? Would be easy
                call GoBack (inunit, 1)
                call ReadAndSplit (inunit, line, words, ierr, ic=.TRUE.)
                
                Z = PeriodicTable (words(1)%str)
                
                if (Z == 0) then
                    ! words(1)%str was NOT an atom. It needs to be
                    call stophere ('Input read error: Unknown key in a MATERIAL block : ' // line)
                    
                else
                    ! words(1)%str WAS an atom, so process it
                    
                    read(words(2)%str,*) N
                    
                    self%Z = [self%Z, Z]
                    self%N = [self%N, N]
                    
                end if
                
            case ('')
                cycle
                
            case ('end')
                ! ! Wrap up and return
                ! if (Zeff < ZERO) then
                !     if (norm294formula) then
                !         ! Use the L2.94-norm-based expression
                !         ALLOCATE(felec(Ne))
                !         felec = Na * Z
                !         felec = felec / sum(felec)
                !         Zeff = sum(felec * (Z**(2.94_KREAL)))**(ONE / 2.94)
                !     else
                !         ! Take Zeff as the total number of electrons per molecule
                !         Zeff = sum(Na * Z)
                !     end if
                ! end if
                ! ! Otherwise, Zeff was already written from user input
                
                return
                
            end select
        end do
        
    End Subroutine
    
    Module Subroutine ReadPNScattering (self)
        use constants
        use IO
        use InputModule
        Implicit None
        Class (InputFile), Intent (InOut) :: self
        
        Integer                           :: ierr
        Character (LEN=:), Allocatable    :: line
        Type (WordsType),  Allocatable    :: words (:)
        
        do
            call ReadSplitSimplify (inunit, line, words, ierr, ic=.TRUE.)
            
            select case (words(1)%str)
            case ('pnscattering')
                if (ierr < 0) then
                    call stophere ('Input read error: Missing value for PNSCATTERING key.')
                end if
                
                read(words(2)%str,*) self%L
                
            case ('')
                cycle
                
            case ('end')
                return
                
            case default
                call stophere ('Input read error: Unknown key in ANGULAR block : ' // words(1)%str)
                
            end select
        end do
        
    End Subroutine
    
    Module Subroutine ParticlesInputRead (self)
        use constants
        use IO
        use InputModule
        Implicit None
        Class (ParticlesInput), Intent (Out) :: self
        
        Integer                              :: p
        Integer                              :: ierr
        Character (LEN=:),      Allocatable  :: line
        Type (WordsType),       Allocatable  :: words (:)
        Type (ParticleInput)                 :: blank
        
        ! Initialize
        ALLOCATE(self%v(0))
        
        do
            call ReadSplitSimplify (inunit, line, words, ierr, ic=.TRUE.)
            select case (words(1)%str)
            case default
                
                ! Read particle name in the true case (not lower case)
                call GoBack (inunit, 1)
                call ReadAndSplit (inunit, line, words, ierr, ic=.TRUE.)
                
                call blank%Read ()
                
                ! Name gets set after reading because ParticleInputRead uses Intent (Out)
                blank%name = line
                
                ! Names must be unique (although PHYSICS not necessarily so)
                do p = 1, SIZE(self%v)
                    if (self%v(p)%name == blank%name) then
                        call stophere ('Input read error: Duplicate particle name ' // blank%name // '. '         // & 
                                       'If you wish to transport multiple particles of the same kind, name them ' // &
                                       'differently, but give them the same PHYSICS key.')
                    end if
                end do
                
                ! If the (case-insensitive) name is a known particle identity, then 
                ! its physics MUST belong to that particle identity.
                ! Otherwise its physics MUST have been provided in the PHYSICS key.
                if (streq(blank%name, 'electrons')) then
                    blank%identity = 'electrons'
                else if (streq(blank%name, 'photons')) then
                    blank%identity = 'photons'
                end if
                
                if (.not. ALLOCATED(blank%identity)) then
                    call stophere ('Input read error: Particle ' // blank%name // ' was given without a PHYSICS key.')
                end if
                
                self%v = [self%v, blank]
                
            case ('')
                cycle
                
            case ('end')
                exit
                
            end select
            
        end do
        
    End Subroutine
    
    Module Subroutine ParticleInputRead (self)
        use constants
        use IO
        use InputModule
        Implicit None
        Class (ParticleInput), Intent (Out) :: self
        
        Integer                             :: i, s
        Logical                             :: form1
        Logical                             :: form2
        Integer                             :: pos1
        Integer                             :: pos2
        Integer                             :: ierr
        Integer                             :: G
        Real (KREAL)                        :: Emin
        Real (KREAL)                        :: Emax
        Character (LEN=:),     Allocatable  :: part1
        Character (LEN=:),     Allocatable  :: part2
        Character (LEN=:),     Allocatable  :: part3
        Character (LEN=:),     Allocatable  :: wtmp
        Character (LEN=:),     Allocatable  :: line
        Type (WordsType),      Allocatable  :: words (:)
        
        do
            call ReadSplitSimplify (inunit, line, words, ierr, ic=.TRUE.)
            
            select case (words(1)%str)
            case ('physics')
                select case (words(2)%str)
                case ('electrons', 'photons'); self%identity = words(2)%str
                case default
                    ! MHY LATER - when I implement custom physics, I'll allow this block to point
                    ! to a directory/file.
                    call stophere ('Input read error: Unknown particle in PHYSICS key : ' // words(2)%str)
                end select
            case ('grid')
                ! THREE OPTIONS - Specify lin/log/exp, give filename, or give like Emin:Emax:G
                ! If option 1, must specify min max and groups
                
                select case (words(2)%str)
                case ('linear', 'logarithmic', 'exponential'); self%gridstruct = words(2)%str
                case default
                    ! It is not any of the above, but is it Emin:Emax:G or a file name?
                    ! Consider the following forms
                    
                    ! Form 1: Emin:Emax:G. Indicated by two semi-colons and no commas
                    s = 0
                    wtmp = words(2)%str
                    do i = 1, LEN_TRIM(wtmp)
                        if (wtmp(i:i) == ':') s = s + 1
                    end do
                    form1 = s == 2 .AND. INDEX(wtmp, ',') == 0
                    
                    ! Form 2: File name. Indicated by no semi-colons, no commas.
                    form2 = INDEX(wtmp, ':') == 0 .AND. INDEX(wtmp, ',') == 0
                    
                    if (form1) then
                        self%gridstruct = 'linear'
                        
                        pos1 = INDEX(wtmp, ':')
                        pos2 = INDEX(wtmp(pos1+1:), ':') + pos1
                        
                        part1 = wtmp(1:pos1-1)
                        part2 = wtmp(pos1+1:pos2-1)
                        part3 = wtmp(pos2+1:)
                        
                        read(part1,*) Emin
                        read(part2,*) Emax
                        read(part3,*) G
                        
                        self%Emin = Emin
                        self%Emax = Emax
                        self%G    = G
                        
                    else if (form2) then
                        
                        self%gridstruct = 'null'
                        
                        ! Read file name in the true case (not lower case)
                        call GoBack (inunit, 1)
                        call ReadAndSplit (inunit, line, words, ierr, ic=.TRUE.)
                        self%fname = words(2)%str
                        
                    else
                        call stophere ('Input read error: Unknown value for PARTICLE GRID key  : ' // words(2)%str)
                        
                    end if
                    
                end select
                
            case ('min')
                read(words(2)%str,*) Emin
                
                self%Emin = Emin
                
            case ('max')
                read(words(2)%str,*) Emax
                
                self%Emax = Emax
                
            case ('groups')
                read(words(2)%str,*) G
                
                self%G = G
                
            case ('solver')
                
                select case (words(2)%str)
                case ('gmres')
                    self%solver = 'GMRES'
                    
                case ('si')
                    self%solver = 'SI'
                    
                case default
                    call stophere ('Input read error: Unknown value for PARTICLE solver key  : ' // words(2)%str)
                    
                end select
                
            case ('angular')
                select case (words(2)%str)
                case ('sn')
                    self%angular = 'SN'
                    
                case ('pn')
                    self%angular = 'PN'
                    
                case default
                    call stophere ('Input read error: Unknown value for PARTICLE angular key : ' // words(2)%str)
                    
                end select
                
            case ('noscatter')
                self%noscat = .TRUE.
                
            case ('')
                cycle
                
            case ('end')
                exit
                
            case default
                call stophere ('Input read error: Unknown key in PARTICLE block : ' // words(1)%str)
                
            end select
            
        end do
        
    End Subroutine
    
    Module Subroutine BeamsInputRead (self)
        ! Because the Beam header is the only repeated header/subheader that
        ! has no enforcement of ordering, this routine is structured uniquely.
        ! That is, it acts as an InOut, where whenever a beam is read in, it's
        ! just appended directly to the array inside of the BeamsInput object.
        ! Furthermore, unlike other repeated headers/subheaders, it is permitted
        ! for the beam headers to be non-consecutive, so in init_InputFile,
        ! this routine is called as many times as there are headers.
        use constants
        use IO
        use InputModule
        Implicit None
        Class (BeamsInput), Intent (InOut) :: self
        
        Integer                            :: i, j
        Logical                            :: axisset = .FALSE.
        Integer                            :: ierr
        Integer                            :: n
        Integer                            :: ni
        Integer                            :: nr
        Real (KREAL)                       :: axis  (3)
        Character (LEN=:), Allocatable     :: line
        Type (WordsType),  Allocatable     :: words (:)
        Type (BeamInput)                   :: blank
        
        if (.not. ALLOCATED(self%v)) ALLOCATE(self%v(0))
        
        do
            call ReadSplitSimplify (inunit, line, words, ierr, ic=.TRUE.)
            
            select case (words(1)%str)
            case ('particle')
                read(words(2)%str,*) blank%particle
                
            case ('polychromatic')
                if (ierr < 0) then
                    call stophere ('Input read error: Missing file name for BEAM Polychromatic key.')
                end if
                
                ! Read file name in the true case (not lower case)
                call GoBack (inunit, 1)
                call ReadAndSplit (inunit, line, words, ierr, ic=.TRUE.)
                blank%spectrum = words(2)%str
                
            case ('monochromatic')
                blank%spectrum = 'constant' ! While monochromatic =/= constant, this just indicates to a given %Populate call that you want constant bounds between Ei and Ef.
                
                ! User can either provide a list of integers, which gives the groups to populate, or they can
                ! provide a number of groups, and E0, and a dE, and we force this into the energy group spectrum.
                ! Note that the 2nd option cannot be done if a particle's energy grid is given from a file.
                ! In that instance, a user should just use the first option + customize their energy grid
                ! to yield whatever desired widths/centers
                
                ! Determine number of integers and number of reals in words(:)
                ni = 0
                nr = 0
                do j = 2, SIZE(words)
                    i = INDEX(words(j)%str,'.')
                    if (i == 0) then
                        ni = ni + 1
                    else
                        nr = nr + 1
                    end if
                end do
                
                if (nr == 2) then
                    ! If two reals were given then this is an E0 specification. You must find 1 integer and it must be the first
                    if (ni /= 1) then
                        call stophere ('Input read error: Improper input after BEAM MONOCHROMATIC key. ' // &
                                       'Please consult documentation.')
                    end if
                    ! Verify the first entry is the integer
                    i = INDEX(words(2)%str,'.')
                    if (i /= 0) then
                        call stophere ('Input read error: Improper input after BEAM MONOCHROMATIC key. ' // &
                                       'Please consult documentation.')
                    end if
                    
                    ALLOCATE(blank%mcgroups(1)) ! In this case, mcgroups isn't the set of groups to use, it's the number of groups that are spanned by dE
                    read(words(2)%str,*,iostat=ierr) blank%mcgroups (1)
                    read(words(3)%str,*) blank%E0
                    read(words(4)%str,*) blank%dE
                else
                    ! If two reals were not given, then there ought to be NO reals
                    if (nr /= 0) then
                        call stophere ('Input read error: Improper input after BEAM MONOCHROMATIC key. ' // &
                                       'Please consult documentation.')
                    end if
                    ! We assume the remaining integers are the groups
                    ALLOCATE(blank%mcgroups(ni))
                    do j = 2, ni + 1
                        read(words(j)%str,*) blank%mcgroups(j - 1)
                    end do
                end if
                
            case ('weight')
                
                read(words(2)%str,*) blank%weight
                
            case ('origin')
                if (ierr < 0) then
                    call stophere ('Input read error: Missing parameters for BEAM Origin key.')
                end if
                
                read(words(2)%str,*) blank%origin(1)
                read(words(3)%str,*) blank%origin(2)
                read(words(4)%str,*) blank%origin(3)
                
            case ('axis')
                if (ierr < 0) then
                    call stophere ('Input read error: Missing parameters for BEAM Axis key.')
                end if
                
                ! First determine if one or three parameters were given
                ! In the slab case, zero or one parameters are needed (i.e., this is optional). 
                ! In the general case, exactly three are needed
                n = SIZE(words) - 1
                
                select case (n)
                case (1)
                    read(words(2)%str,*) axis(1)
                    blank%slab = .TRUE.
                    
                case (3)
                    read(words(2)%str,*) axis(1)
                    read(words(3)%str,*) axis(2)
                    read(words(4)%str,*) axis(3)
                    
                    ! Normalize the axis
                    blank%axis = axis / NORM2(axis)
                    
                    blank%slab = .FALSE.
                case default
                    call stophere ('Input read error: Incorrect number of parameters for BEAM Axis key.')
                end select
                
                axisset = .TRUE.
                
            case ('cutout')
                if (ierr < 0) then
                    call stophere ('Input read error: Missing all parameters for BEAM Cutout key.')
                end if
                
                select case (words(2)%str)
                case ('none')
                    blank%cutout = 'none'
                    
                    ! No parameters to read
                case ('circle')
                    blank%cutout = 'circle'
                    
                    read(words(3)%str,*,iostat=ierr) blank%cutoutparams(1)
                    read(words(4)%str,*,iostat=ierr) blank%cutoutparams(2)
                    
                    if (ierr < 0) then
                        call stophere ('Input read error: Missing Cutout Circle parameters. ' // &
                                       'You must enter two.')
                    end if
                    
                case ('rectangle')
                    blank%cutout = 'rectangle'
                    
                    read(words(3)%str,*,iostat=ierr) blank%cutoutparams(1)
                    read(words(4)%str,*,iostat=ierr) blank%cutoutparams(2)
                    read(words(5)%str,*,iostat=ierr) blank%cutoutparams(3)
                    read(words(6)%str,*,iostat=ierr) blank%cutoutparams(4)
                    
                    if (ierr < 0) then
                        call stophere ('Input read error: Missing Cutout Rectangle parameters. ' // &
                                       'You must enter four.')
                    end if
                    
                case default
                    call stophere ('Input read error: Unknown cutout type for BEAM Cutout key: ' // words(2)%str)
                end select
                
            case ('spherical', 'planar')
                blank%angdist = words(1)%str
                
            case ('')
                cycle
                
            case ('end')
                
                ! If the axis was not provided it is assumed this corresponds to a slab solve,
                ! which is the only case where the axis is optional.
                ! Consistency with other input blocks is checked during ValidateInput subroutine
                if (.not. axisset) blank%slab = .TRUE.
                
                self%v = [self%v, blank]
                
                return
                
            case default
                call stophere ('Input read error: Unknown key in BEAM block : ' // words(1)%str)
            end select
        end do
        
    End Subroutine
    
    Module Subroutine PostProcInputRead (self)
        use constants
        use IO
        Implicit None
        Class (PostProcInput), Intent (InOut) :: self
        
        Integer                            :: i
        Integer                            :: ierr
        Character (LEN=:), Allocatable     :: line
        Type (WordsType),  Allocatable     :: words (:)
        
        do
            call ReadSplitSimplify (inunit, line, words, ierr, ic=.TRUE.)
            
            select case (words(1)%str)
            case ('fluence');    self%fluence   = .TRUE.
            case ('uncollided'); self%outputunc = .TRUE.
            case ('energy');     self%energy    = .TRUE.
            case ('dose');       self%dose      = .TRUE.
            case ('charge');     self%charge    = .TRUE.
            case ('');           cycle
            case ('end');        return
            case default;        call stophere ('Input read error: Unknown key in POSTPROCESSING block : ' // words(1)%str)
            end select
        end do
        
    End Subroutine
    
End Submodule