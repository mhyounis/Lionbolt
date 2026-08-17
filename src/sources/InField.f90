Submodule (Sources) submodInField
Contains
Module Function InField (self, r) Result (l)
    ! This function returns TRUE if r is within the field geometry,
    ! FALSE otherwise.
    use constants
    use IO
    use BasicMathFunctions
    use Sources
    Implicit None
    Class (FieldGeometry), Intent (In) :: self
    Real (KREAL),          Intent (In) :: r (3)
    
    Logical                            :: l
    
    Real (KREAL),          Parameter   :: TOL = 1.0e-12_KREAL
    Real (KREAL)                       :: nu
    Real (KREAL)                       :: h
    Real (KREAL)                       :: dx
    Real (KREAL)                       :: dy
    Real (KREAL)                       :: x
    Real (KREAL)                       :: u (3)
    Real (KREAL)                       :: v (3)
    
    if (streq(self%cutout, 'none')) then
        ! The point is definitely within an open collimator field.
        ! The only edge case is if you have a planar field and then somehow geometry
        ! on the other side of the field. Which is mathematically possible, I suppose, when I
        ! implement non-convex geometries, but it is unphysical so I won't treat it.
        l = .TRUE.
        return
    end if
    
    if (streq(self%angdist, 'planar')) then
        ! Both cases need (r - R0) perp k0
        
        if (streq(self%cutout, 'circle')) then
            ! Check if the point is within a cylinder
            
            ! ALGORITHM : Check the magnitude of (r - R0) perp k0, to see if it is less than the circle radius
            
            u = r - self%origin
            v = u - DOT_PRODUCT(u, self%axis) * self%axis ! k0 is normalized so we do not need the |k0|^{2} denominator
            
            l = NORM2(v) <= self%cutoutparams(2) + TOL
            
        else if (streq(self%cutout, 'rectangle')) then
            ! Check if the point is within a rectangular prism
            
            ! ALGORITHM : Check the components of (r - R0) perp k0, to see if they are less than the rotated 
            !             field cutout. Thus, we further project v onto e1 and e2 and check the magnitude
            
            u = r - self%origin
            v = u - DOT_PRODUCT(u, self%axis) * self%axis
            
            l = .TRUE.
            
            l = l .and. ABS(DOT_PRODUCT(v, self%e1)) <= HALF * self%cutoutparams(3) + TOL
            
            if (.not. l) return
            
            l = l .and. ABS(DOT_PRODUCT(v, self%e2)) <= HALF * self%cutoutparams(4) + TOL
            
        end if
    else if (streq(self%angdist, 'spherical')) then
        if (streq(self%cutout, 'circle')) then
            ! Check if the point is within a cone
            
            ! ALGORITHM : Compare the angle cosine of (r - R0) with respect to k0, and see if it exceeds
            !             what is allowed by the radius.
            
            ! The radius rho, defined at S (distance from R0) allows cos(atan(rho / S)), thus:
            x = self%cutoutparams(2) / self%cutoutparams(1)
            nu = ONE / SQRT(ONE + x * x)
            
            ! The angle cosine spanned by (r - R0) is (r - R0) * k0 / |r - R0|
            u = r - self%origin
            v = u / NORM2(u) ! This is (r - R0) / |r - R0|
            
            ! The criterion is that the angle spanned by (r - R0) must be less than nu
            l = DOT_PRODUCT(v, self%axis) <= nu + TOL
            
        else if (streq(self%cutout, 'rectangle')) then
            ! Check if the point is within a pyramid
            
            ! ALGORITHM : Here we do something similar to the planar case, except we compare with
            !             a rectangle scaled to the depth of the point.
            
            ! First determine the depth along the axis of the beam
            u = r - self%origin
            
            h = DOT_PRODUCT(u, self%axis) ! Keep in mind this is the SIGNED distance. Deliberate choice, because points behind the beam should not be counted (unrealistic as they may be)
            
            v = u - h * self%axis
            
            ! The field sizes at this depth can be given like d_2 = h * d_1 / S, where
            ! d1 is the field size at depth S1.
            ! But we actually want half of the field size on one side, for our check later on.
            
            dx = HALF * h * self%cutoutparams(3) / self%cutoutparams(1)
            dy = HALF * h * self%cutoutparams(4) / self%cutoutparams(1)
            ! print *, dx, dy
            ! Now we take a vector to originate at the point along the beam axis which intersects
            ! this depth. This is still just (r - R0) perp k0. Process this vector through the same
            ! criterion as the planar case but with dx and dy instead of cutoutparams
            
            l = .TRUE.
            
            l = l .and. ABS(DOT_PRODUCT(v, self%e1)) <= dx + TOL
            
            if (.not. l) return
            
            l = l .and. ABS(DOT_PRODUCT(v, self%e2)) <= dy + TOL
            
        end if
    end if
    
End Function
End Submodule