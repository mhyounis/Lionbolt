Submodule (InputModule) submodReadInput
Contains
Module Function ReadInput (fname) Result (self)
    use globals
    use IO
    use InputModule
    Implicit None
    Character (*),     Optional, Intent (In) :: fname
    
    Type (InputFile)                         :: self
    
    Integer                                  :: p
    Logical                                  :: donemesh      = .FALSE.
    Logical                                  :: donematerials = .FALSE.
    Logical                                  :: doneangular   = .FALSE.
    Logical                                  :: doneparticles = .FALSE.
    Logical                                  :: doneabeam     = .FALSE.
    Logical                                  :: donepostproc  = .FALSE.
    Integer                                  :: ierr
    Character (LEN=:),          Allocatable  :: line
    Type (WordsType),           Allocatable  :: words (:)
    
    if (PRESENT(fname)) then
        open(unit=inunit, file=fname, action='read')
    else
        call GetInputFname ()
        
        open(unit=inunit, file=InputFname, action='read')
    end if
    
    call self%problem%Read ()
    
    ! Loop through headers
    do
        ! call ReadLine (inunit, line, ierr, ic=.TRUE.)
        call ReadSplitSimplify (inunit, line, words, ierr, ic=.TRUE.)
        
        if (ierr < 0) exit
        
        select case (line)
        case ('mesh')
            
            if (donemesh) call stophere ('Input read error: Duplicate MESH block.')
            
            call self%mesh%Read ()
            
            donemesh = .TRUE.
            
        case ('materials')
            
            if (donematerials) call stophere ('Input read error: Duplicate MATERIALS block.')
            
            call self%materials%Read ()
            
            donematerials = .TRUE.
            
        case ('angular')
            
            if (doneangular) call stophere ('Input read error: Duplicate ANGULAR block.')
            
            call self%ReadPNScattering ()
            
            doneangular = .TRUE.
            
        case ('particles')
            
            if (doneparticles) call stophere ('Input read error: Duplicate PARTICLES block.')
            
            call self%particles%Read ()
            
            doneparticles = .TRUE.
            
        case ('beam')
            
            call self%beams%Read ()
            
            doneabeam = .TRUE. ! Duplicates are allowed, so this key really just helps me know if the user provided at least one beam
            
        case ('postprocessing')
            
            if (donepostproc) call stophere ('Input read error: Duplicate POSTPROCESSING block.')
            
            call self%postproc%Read ()
            
            donepostproc = .TRUE.
        
        case ('memorylimited')
            self%options%memlim = .TRUE.
            
        case ('storagelimited')
            self%options%storlim = .TRUE.
            
        case ('debug')
            self%options%debug = .TRUE.
            
        case ('')
            cycle
            
        case default
            call stophere ('Input read error: Invalid header : ' // line)
            
        end select
    end do
    
    close(unit=inunit)
    
    ! Now check that every required block has been visited
    if (.not. donemesh)      call stophere ('Input read error: MESH block was not found.')
    if (.not. donematerials) call stophere ('Input read error: MATERIALS block was not found.') ! MHY LATER - this will only be when customphysics is not used
    if (.not. doneangular)   call stophere ('Input read error: ANGULAR block was not found.')
    if (.not. doneparticles) call stophere ('Input read error: PARTICLES block was not found.')
    if (.not. doneabeam)     call stophere ('Input read error: Must provide at least on BEAM block.')
    
    ! Now assign the Legendre order and slab logical to all particles
    do p = 1, SIZE(self%particles%v)
        self%particles%v(p)%L    = self%L
        self%particles%v(p)%slab = self%problem%slab
    end do
    
    ! Verify all necessary parameters have been set and that there are no conflicts.
    call self%Validate ()
    
End Function
End Submodule