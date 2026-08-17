Submodule (Sources) submodGenerateENSRCFile
Contains
Module Subroutine GenerateENSRCFile (self, energy, fname)
    use constants
    use Interpolation
    use NittanyAPI
    
    ! This routine will generate a spectrum file suitable for EGS calculations
    ! directly from the MULTIGROUP discretization of the user's spectrum
    ! (not from the user's spectrum itself).
    ! This is ideal for validation. It is not currently accessible via
    ! input but that may be arranged soon.
    
    Implicit None
    Class (ExternalBeam), Intent (In) :: self
    Type (EnergyGrid),    Intent (In) :: energy
    Character (*),        Intent (In) :: fname
    
    Integer                           :: g, i, j
    Integer,              Parameter   :: unit  = 10
    Integer,              Parameter   :: IMODE = 0
    Integer                           :: ierr
    Real (KREAL)                      :: norm
    Real (KREAL)                      :: binstart
    Real (KREAL)                      :: tmp2   (2)
    Real (KREAL)                      :: binend (energy%G) ! Bin ENDPOINTS (requested by EGS)
    Real (KREAL)                      :: N      (energy%G) ! Bin counts (we use IMODE = 0)
    Real (KREAL),         Allocatable :: E      (:)
    Real (KREAL),         Allocatable :: f      (:)
    
    if (.not. streq(self%fldgeo%spectrum, 'constant')) then
        ! Open the file and do lin-lin interp
        open(unit=unit, file=self%fldgeo%spectrum, action='read', iostat=ierr)
        if (ierr /= 0) then
            call stophere ('GenerateENSRCFile.f90: GenerateENSRCFile: ' // &
                           "Spectrum could not be opened. Are you sure it's a valid file?")
        end if
        j = 0
        do
            read(unit,*,iostat=ierr) tmp2
            if (ierr /= 0) exit
            j = j + 1
        end do
        close(unit)
        
        ALLOCATE(E(j))
        ALLOCATE(f(j))
        open(unit=unit, file=self%fldgeo%spectrum, action='read')
        do i = 1, j
            read(unit,*,iostat=ierr) E(i), f(i)
        end do
        close(unit)
        
        norm = Trapezoidal1D (E, f, energy%E(energy%G + 1), energy%E(1), LINLIN)
        
        ! Bin start is just the minimum energy
        binstart = energy%E(energy%G + 1)
        
        ! Bin ENDPOINTS are given like this
        binend = energy%E(1 : energy%G)
        
        ! Now determine bin counts
        do g = 1, energy%G
            N(g) = Trapezoidal1D (E, f, energy%E(g + 1), energy%E(g), LINLIN) / norm ! Normalization is not required by EGS though
        end do
        
    else
        ! Still working on this. Should be simple but may require some logic or inputs
        ! I don't have here
    end if
    
    ! Now write
    open(unit=unit, file=fname, status='replace', action='write')
    write (unit,'(A)') 'Lionbolt_generated_spectrum' ! Default name, can be changed by user
    write (unit,'(I0,", ",F12.6,", ",I0)') energy%G, binstart, IMODE
    do g = energy%G, 1, -1
        write (unit, '(F7.5,", ",ES13.6)') binend(g), N(g)
    end do
    close(unit)
    
End Subroutine
End Submodule