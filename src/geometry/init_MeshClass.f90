Submodule (Geometry) submodinit_MeshClass
Contains
Module Function init_MeshClass (fname) Result (self)
    use Geometry
    Implicit None
    Character (*), Intent (In) :: fname
    
    Type (MeshClass)           :: self
    
    self%slab = .FALSE.
    
    call self%FromFile (fname)
    
End Function
End Submodule