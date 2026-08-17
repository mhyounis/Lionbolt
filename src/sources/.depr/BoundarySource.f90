Submodule (Sources) submodBoundarySource
Contains
Module Subroutine BoundarySource (self, s)
    ! DEPRECATED
    use LinearOperatorClass
    use Sources
    Implicit None
    Class (ExternalBeam),    Intent (InOut) :: self
    
    Type (SpaceAngleVector), Intent (InOut) :: s
    
    if (.not. ALLOCATED(self%rows)) call self%PrepSource ()
    
    s%v(self%rows) = s%v(self%rows) + self%n * self%vals
    
End Subroutine
End Submodule