Submodule (PostProcInterface) submodTallyDeposition
Contains
Module Subroutine TallyDeposition (mesh, XSDEP, fl, DEP)
    use constants
    use types
    use IO
    use Geometry
    Implicit None
    Type (MeshClass), Intent (In)    :: mesh
    Real (KREAL),     Intent (In)    :: XSDEP  (:) ! Deposition cross section for this material     ! (m)
    Real (KREAL),     Intent (In)    :: fl     (:) ! Fluence                                        ! (sdof)
    
    Real (KREAL),     Intent (InOut) :: DEP    (:) ! Deposition                                     ! (sdof)
    
    Integer                          :: m
    
    ! mesh%mat2sd(m)%v gives spatial dofs corresponding to material m
    do m = 1, mesh%Nm
        DEP(mesh%mat2sd(m)%v) = DEP(mesh%mat2sd(m)%v) + XSDEP(m) * fl(mesh%mat2sd(m)%v)
    end do
    
End Subroutine
End Submodule