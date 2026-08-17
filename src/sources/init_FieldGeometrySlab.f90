Submodule (Sources) submodinit_FieldGeometrySlab
Contains
Module Function init_FieldGeometrySlab (weight, spectrum, angle) Result (self)
    use constants
    use IO
    use Sources
    Implicit None
    Real (KREAL),        Intent (In) :: weight
    Character (*),       Intent (In) :: spectrum
    Real (KREAL),        Intent (In) :: angle
    
    Type (FieldGeometry)             :: self
    
    ! I must make it so that you can check if a given string (like spectrum) is:
    !    an absolute file
    !    a relative file
    !    neither
    ! For now I'm going to just assume it's a valid file name
    ! Should also have a function that can turn relative files into absolute ones
    ! and I should process all user input through that as well as make it visible to users of LionboltAPI (HOW? core is technically internal)
    
    self%slab      = .TRUE.
    self%weight    = weight
    self%spectrum  = spectrum
    self%axis(1)   = angle
    self%axis(2:3) = ZERO
    self%angdist   = 'planar'
    self%span      = ONE
    
    if (.not. streq(self%spectrum, 'constant')) then
        ! Verify this spectrum is a valid file name
        call VerifyFile (fname=self%spectrum, errmess='init_FieldGeometrySlab.f90: init_FieldGeometrySlab: ' // &
                                                      'Spectrum file does not exist --- ' // self%spectrum)
    end if
    
End Function
End Submodule