Submodule (Sources) submodinit_FieldGeometryGen
Contains
Module Function init_FieldGeometryGen (weight, spectrum, origin, axis, angdist, cutout, cutoutparams) Result (self)
    use constants
    use IO
    use BasicMathFunctions
    use Sources
    Implicit None
    Real (KREAL),        Intent (In) :: weight
    Character (*),       Intent (In) :: spectrum
    Real (KREAL),        Intent (In) :: origin       (3)
    Real (KREAL),        Intent (In) :: axis         (3)
    Character (*),       Intent (In) :: angdist
    Character (*),       Intent (In) :: cutout
    Real (KREAL),        Intent (In) :: cutoutparams (:)
    
    Type (FieldGeometry)             :: self
    
    Integer                          :: nparam
    Real (KREAL)                     :: a
    Real (KREAL)                     :: b
    Real (KREAL)                     :: c
    Real (KREAL)                     :: tmp (3)
    
    if (streq(cutout, 'none')) then
        ! Don't care what params they input
        
        nparam = 0
    else if (streq(cutout, 'circle')) then
        ! 2 params
        if (SIZE(cutoutparams) < 2) then
            call stophere ('init_FieldGeometryGen.f90: init_FieldGeometryGen: ' // &
                           'For a circle cutout, at least 2 cutout parameters are needed.')
        end if
        
        nparam = 2
    else if (streq(cutout, 'rectangle')) then
        ! 4 params
        if (SIZE(cutoutparams) < 4) then
            call stophere ('init_FieldGeometryGen.f90: init_FieldGeometryGen: ' // &
                           'For a circle cutout, at least 4 cutout parameters are needed.')
        end if
        
        nparam = 4
    else
        call stophere ('init_FieldGeometryGen.f90: init_FieldGeometryGen: ' // &
                       'Unknown cutout used : ' // TRIM(ADJUSTL(cutout)))
    end if
    
    if ((.not. streq(angdist, 'planar')) .and. (.not. streq(angdist, 'spherical'))) then
        call stophere ('init_FieldGeometryGen.f90: init_FieldGeometryGen: ' // &
                       'Unknown angdist used : ' // TRIM(ADJUSTL(angdist)))
    end if
    
    self%slab                   = .FALSE.
    self%weight                 = weight
    self%spectrum               = spectrum
    self%origin                 = origin
    self%axis                   = axis / NORM2(axis) ! Normalize the axis for the user
    self%angdist                = angdist
    self%cutout                 = cutout
    self%cutoutparams(1:nparam) = cutoutparams(1:nparam)
    
    if (.not. streq(self%spectrum, 'constant')) then
        ! Verify this spectrum is a valid file name
        call VerifyFile (fname=self%spectrum, errmess='init_FieldGeometryGen.f90: init_FieldGeometryGen: ' // &
                                                      'Spectrum file does not exist. --- ' // self%spectrum)
    end if
    
    ! Now assign the beam span ! I DON'T THINK I NEED THIS ANYMORE.
    if (streq(self%angdist, 'planar')) then
        ! Gives area
        select case (SimplifyString(self%cutout))
        case ('none');      self%span = ONE
        case ('circle');    self%span = TWOPI * self%cutoutparams(2) * self%cutoutparams(2)
        case ('rectangle'); self%span = self%cutoutparams(3) * self%cutoutparams(4)
        end select
    else if (streq(self%angdist, 'spherical')) then
        select case (SimplifyString(self%cutout))
        case ('none')
            self%span = FOURPI
        case ('circle')
            a = self%cutoutparams(2)
            c = self%cutoutparams(1)
            self%span = TWOPI * (ONE - c / SQRT(c * c + a * a))
        case ('rectangle')
            a = self%cutoutparams(3)
            b = self%cutoutparams(4)
            c = self%cutoutparams(1)
            self%span = FOUR * ATAN(a * b / (TWO * SQRT(FOUR * c * c + a * a + b * b)))
        end select
    end if
    
    ! Form the beam-frame x and y vectors (e1 and e2 respectively)
    ! The beam's x and y axes are given by rotating a standard xyz frame
    ! using polar and azimuthal angles such that the beam axis corresponds to z.
    ! The beam rotation angle is thus a final intrinsic rotation with respect to this new z axis.
    ! Thus obtain from axis, alpha and beta, then use cutoutparams(2) as gamma.
    if (streq(self%cutout, 'rectangle')) then
        
        a = ATAN2(self%axis(2), self%axis(1))
        b = ACOS (self%axis(3))
        c = cutoutparams(2) * PI / 180_KREAL ! Convert degrees to radians
        
        call ZYZBasis (a, b, c, self%e1, self%e2, tmp) ! Throw away the third vector because we have axis already
        
    end if
    
End Function
End Submodule