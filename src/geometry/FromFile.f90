Submodule (Geometry) submodFromFile
Contains
Module Subroutine FromFile (self, fname)
    use Geometry
    Implicit None
    Class (MeshClass), Intent (InOut) :: self
    Character (*),     Intent (In)    :: fname
    
    if (self%init) then
        call stophere ('FromFile.f90: FromFile: Attempting to re-initialize an already existing mesh. Use destroy () first.')
    end if
    
    self%init = .TRUE.
    self%slab = .FALSE.
    
    ! Currently only GMSH reader is implemented
    call ReadGMSH (fname, self)
    
End Subroutine
End Submodule