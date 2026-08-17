Submodule (Sources) submodBoundaryFluence
Contains
Module Function BoundaryFluence (fldgeo, r, k, override) Result (phi)
    ! DEPRECATED 
    ! Private function which facilitates the construction of the boundary source
    ! or the ray traced angular fluence. Should not be visible in Lionbolt API,
    ! is not user friendly.
    ! Note for developers, phi is returned WITH an (unnormalized) angular distribution ONLY if the 
    ! beam angular distribution is spherical. Otherwise you are expected to apply
    ! that factor externally
    use constants
    use IO
    use BasicMathFunctions
    use Polynomials
    use Sources
    use, INTRINSIC :: ieee_arithmetic
    Implicit None
    Class (FieldGeometry), Intent (In) :: fldgeo
    Real (KREAL),          Intent (In) :: r (3) ! For slab case, just use the first entry
    Real (KREAL),          Intent (In) :: k (3) ! For slab case, just use the first entry
    Logical, Optional,     Intent (In) :: override
    
    Real (KREAL)                       :: phi
    
    Logical                            :: proceed
    Logical                            :: doangular
    Real (KREAL)                       :: A
    Real (KREAL)                       :: G
    Real (KREAL)                       :: k0 (3)
    Real (KREAL)                       :: R0 (3)
    
    ! First determine if this point is on the boundary.
    ! Or override this check if the user desires
    if (PRESENT(override)) then
        if (override) then
            proceed = .TRUE.
        else
            proceed = fldgeo%InField (r)
        end if
    else
        proceed = fldgeo%InField (r)
    end if
    
    ! If the point is not in the field, return with NaN.
    ! NOTE - Returning with NaN is a very motivated decision, but perhaps not user-friendly.
    ! In PrepSource I use NaN to tell me that a point is not in the field, so I can cycle.
    ! It is a stricter condition than ZERO. However I think perhaps ZERO may be preferable for
    ! user-friendliness as well as sufficiently safe for the later functions called in PrepSource.
    if (.not. proceed) then
        phi = IEEE_VALUE(phi, ieee_quiet_nan) ! Gives NaN
        return
    end if
    
    if (fldgeo%slab) then
        
        !  =============
        !    SLAB CASE  
        !  =============
        
        ! Angular factor must be applied externally
        ! ! Convert angle of incidence to radians then take the cosine
        ! k0(1) = fldgeo%axis(1) * PI / 180.0_KREAL
        ! k0(1) = COS(k0(1))
        
        ! A = EvalBeamAngDist (k(1), k0(1))
        phi = ONE
        
    else if (streq(fldgeo%angdist, 'planar')) then
        
        !  ====================
        !    PLANAR BEAM CASE  
        !  ====================
        
        ! ! Angular factor must be applied externally
        ! k0 = fldgeo%axis
        ! 
        ! A = EvalBeamAngDist (k, k0)
        
        phi = ONE / fldgeo%span
        
    else if (streq(fldgeo%angdist, 'spherical')) then
        
        !  =======================
        !    SPHERICAL BEAM CASE  
        !  =======================
        
        R0 = fldgeo%origin
        
        k0 = (r - R0) / NORM2(r - R0)
        
        A = EvalBeamAngDist (k, k0)
        G = ONE / DOT_PRODUCT(r - R0, r - R0)
        
        phi = A * G
        
    end if
    
End Function
End Submodule