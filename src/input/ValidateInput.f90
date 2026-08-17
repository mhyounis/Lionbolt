Submodule (InputModule) submodValidateInput
Contains
Module Subroutine ValidateInput (self)
    use constants
    use IO
    use InputModule
    Implicit None
    Class (InputFile), Intent (In) :: self
    
    Integer                        :: s, m
    Integer,           Parameter   :: unit = 10
    Integer                        :: ierr
    Logical                        :: exists
    Integer                        :: Nm1
    Integer                        :: Nm2
    Character (LEN=:), Allocatable :: line
    
    !  ========
    !    MESH  
    !  ========
    
    if (self%problem%General) then
        ! Check that no slabs have been defined
        if (SIZE(self%mesh%slabs) > 0) then
            call stophere ('Input validation error: Using PROBLEM GENERAL, but slabs were specified.')
        end if
        
        ! Check that mesh file has been given and exists
        if (.not. ALLOCATED(self%mesh%fname)) then
            call stophere ('Input validation error: Using PROBLEM GENERAL, but mesh file was not specified.')
        end if
        
        inquire (file=trim(self%mesh%fname), exist=exists)
        
        if (.not. exists) then
            call stophere ('Input validation error: Provided mesh file does not exist: ' // self%mesh%fname)
        end if
    else
        ! Check that no mesh file has been given
        if (ALLOCATED(self%mesh%fname)) then
            call stophere ('Input validation error: Using PROBLEM SLAB, but mesh file was specified.')
        end if
        
        ! Check that slabs have been defined
        if (SIZE(self%mesh%slabs) == 0) then
            call stophere ('Input validation error: Using PROBLEM SLAB, but no slabs were specified.')
        end if
    end if
    
    !  =============
    !    MATERIALS  
    !  =============
    
    ! Check that all materials have densities, Z, N
    do m = 1, SIZE(self%materials%v)
        if (self%materials%v(m)%DENS <= ZERO) then
            call stophere ('Input validation error: Density negative, zero, or not provided for material ' &
                           // self%materials%v(m)%name // '.')
        end if
        if (SIZE(self%materials%v(m)%Z) == 0 .or. SIZE(self%materials%v(m)%N) == 0) then
            call stophere ('Input validation error: Atomic composition not provided for material ' &
                           // self%materials%v(m)%name // '.')
        end if
    end do
    
    ! Now check that the mesh has an agreeing number of materials, if there is more than one.
    ! If there isn't more than one, user doesn't need to specify physical tags anyway.
    Nm1 = SIZE(self%materials%v)
    if (self%problem%General) then
        if (Nm1 > 1) then
            open(unit=unit, file=self%mesh%fname, action='read')
            do
                call ReadLine (unit, line, ierr)
                
                if (ierr < 0) then
                    call stophere ('Input validation error: Could not find the number of materials in the mesh. ' // &
                                   'Are you sure you defined physical volumes?')
                end if
                
                if (line == '$PhysicalNames') then
                    read(unit,*) Nm2
                    exit
                end if
                
            end do
            close(unit=unit)
            
            if (Nm1 /= Nm2) then
                call stophere ('Input validation error: ' // &
                               'Number of materials specified does not match with number of materials in mesh.')
            end if
            
        end if
    else
        do s = 1, SIZE(self%mesh%slabs)
            if (self%mesh%slabs(s)%material > Nm1) then
                call stophere ('Input validation error: ' // &
                                               'Number of materials specified does not match with number of materials in mesh.')
            end if
        end do
    end if
    
    !  =============
    !    PARTICLES  
    !  =============
    ! This and in NittanyPhysics are still WIP
    ! Check that all required particle properties are given
    
    ! Check that each particle has a valid grid
    ! If fname is given use inquire
    
    !  =========
    !    BEAMS  
    !  =========
    
    ! Enforce that if a user provided an energy grid file, they can't use
    ! monochromatic E0 and dE. They can still specify the energy groups though.
    
    ! Also make sure that Emax is less than or equal to E0 + dE/2
    
    if (self%problem%General) then
        ! Check that all required beam properties are given
        
    else
        ! Check that all required beam properties are given.
        ! For a slab there are things you don't need.
        
    end if
    
End Subroutine
End Submodule