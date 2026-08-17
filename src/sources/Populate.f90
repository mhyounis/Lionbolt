Submodule (Sources) submodGenerateParticles
Contains
Module Subroutine Populate (self, Emin, Emax, Ei, Ef)
    use constants
    use IO
    use Interpolation
    use Sources
    Implicit None
    Class (ExternalBeam), Intent (InOut) :: self
    Real (KREAL),         Intent (In)    :: Emin ! Lower energy bound for normalization (required sadly)
    Real (KREAL),         Intent (In)    :: Emax ! Upper energy bound for normalization
    Real (KREAL),         Intent (In)    :: Ei   ! Lower energy bound for integration
    Real (KREAL),         Intent (In)    :: Ef   ! Upper energy bound for integration
    
    Integer                              :: i
    Integer,              Parameter      :: unit = 10
    Integer                              :: ierr
    Integer                              :: N
    Real (KREAL)                         :: Emin
    Real (KREAL)                         :: Emax
    Real (KREAL)                         :: norm
    Real (KREAL)                         :: tmp2 (2)
    Real (KREAL),         Allocatable    :: E    (:)
    Real (KREAL),         Allocatable    :: f    (:)
    
    if (Ef < Ei) then
        call stophere ('Populate.f90: Populate: Provided final energy is below initial energy.')
    end if
    
    if (streq(self%fldgeo%spectrum, 'constant')) then
        call stophere ('Populate.f90: Populate: Should not be called if the field spectrum is constant. ' // &
                       'Just assign self%n manually')
    end if
    
    ! Open the file and do interpolation
    
    open(unit=unit, file=self%fldgeo%spectrum, action='read', iostat=ierr)
    if (ierr /= 0) then
        call stophere ('Populate.f90: Populate: ' // &
                       "Spectrum could not be opened. Are you sure it's a valid file?")
    end if
    N = 0
    do
        read(unit,*,iostat=ierr) tmp2
        if (ierr /= 0) exit
        N = N + 1
    end do
    close(unit)
    
    ALLOCATE(E(N))
    ALLOCATE(f(N))
    open(unit=unit, file=self%fldgeo%spectrum, action='read')
    do i = 1, N
        read(unit,*,iostat=ierr) E(i), f(i)
    end do
    close(unit)
    
    ! Normalize the user's spectrum
    norm = Trapezoidal1D (E, f, Emin, Emax, LINLIN)
    
    ! Now form the integral, normalize, and weight it
    self%n = self%fldgeo%weight * Trapezoidal1D (E, f, Ei, Ef, LINLIN) / norm
    
End Subroutine
End Submodule