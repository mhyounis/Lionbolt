Module PostProcInterface
    Implicit None
    
    Interface
        
        Module Subroutine CalculateFluence (s, w, fl)
            use constants
            Implicit None
            Real (KREAL), Intent (In)  :: s  (:)
            Real (KREAL), Intent (In)  :: w  (:)
            Real (KREAL), Intent (Out) :: fl (:)
        End Subroutine
        
        Module Subroutine TallyDeposition (mesh, XSDEP, fl, DEP)
            use constants
            use types
            use IO
            use Geometry
            Implicit None
            Type (MeshClass), Intent (In)    :: mesh
            Real (KREAL),     Intent (In)    :: XSDEP (:)
            Real (KREAL),     Intent (In)    :: fl    (:)
            Real (KREAL),     Intent (InOut) :: DEP   (:)
        End Subroutine
        
    End Interface
    
End Module PostProcInterface