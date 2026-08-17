Module SpaceAngleInterface
    use constants
    Implicit None
    
    ! This type just holds addressing info. for the application of the upstream term
    ! due to the sweep.
    ! (Name contradicts wiscobolt notes, where I use F_{\downarrow} for this, but its just because of
    ! silly perspective convention flip-flopping.)
    ! Point being, this is the off-diagonal elements of the transport array that do not form the
    ! boundary source
    Type UpstreamSys
        Integer, Allocatable :: SAoffset (:,:) ! General space-angle offset in the sweep. Gives the value of the first sweep-space-angle index of a given (es,i) block in a vector using collapsed notation
        Integer, Allocatable :: UFoffset (:,:) ! Offset array for the indexing of any upstream-face properties (see below)
        Integer, Allocatable :: upfaces  (:)   ! Flattened list of local face indices
        Integer, Allocatable :: ipsf     (:)   ! Flattened list of starting sweep-space-angle indices for the given local faces ! (See TransportAlgebra.f90)
        Integer, Allocatable :: ipef     (:)   ! Flattened list of ending sweep-space-angle indices for the given local faces   ! (See TransportAlgebra.f90)
    Contains
        Procedure :: Destroy
    End Type
    
    Private :: Destroy
    
    Interface
        
        Module Subroutine BuildSpaceAngle (mesh, angular, sl, isl, IAF)
            use constants
            use types
            use Geometry
            use AngularSpace
            Implicit None
            Type (MeshClass),                 Intent (In)  :: mesh
            Type (AngularClass),              Intent (In)  :: angular
            Integer,             Allocatable, Intent (Out) :: sl  (:,:)
            Integer,             Allocatable, Intent (Out) :: isl (:,:)
            Type (UpstreamSys),               Intent (Out) :: IAF
        End Subroutine
        
        Module Subroutine SweepOrder (mesh, angular, IsUpstream, IUoffset, sl, isl, estart)
            use constants
            use Parallelism
            use profiler
            use types
            use Geometry
            use AngularSpace
            Implicit None
            Type (MeshClass),                Intent (In)  :: mesh
            Type (AngularClass),             Intent (In)  :: angular
            Logical,            Allocatable, Intent (Out) :: IsUpstream (:)   ! Tells you if the given face in the given element for the given angle is upstream
            Integer,            Allocatable, Intent (Out) :: IUoffset   (:,:) ! Offset array for IsUpstream
            Integer,            Allocatable, Intent (Out) :: sl         (:,:) ! Sweep list
            Integer,            Allocatable, Intent (Out) :: isl        (:,:) ! Inverse sweep list
            Integer,            Allocatable, Intent (Out) :: estart     (:)   ! Index of the last element that can be solved without sweep (name is bad)
        End Subroutine
        
        Module Subroutine SweepArrays (mesh, angular, IsUpstream, IUoffset, sl, isl, IAF)
            use constants
            use Parallelism
            use types
            use Geometry
            use AngularSpace
            Implicit None
            Type (MeshClass),    Intent (In)  :: mesh
            Type (AngularClass), Intent (In)  :: angular
            Logical,             Intent (In)  :: IsUpstream (:)
            Integer,             Intent (In)  :: IUoffset   (:,:)
            Integer,             Intent (In)  :: sl         (:,:)
            Integer,             Intent (In)  :: isl        (:,:)
            Type (UpstreamSys),  Intent (Out) :: IAF
        End Subroutine
        
        Module Subroutine EncodeSweep (encode, NE, NK, NENK, NI, offset, SAoffset, sl, t, s)
            use constants
            Implicit None
            Logical,      Intent (In)    :: encode
            Integer,      Intent (In)    :: NE
            Integer,      Intent (In)    :: NK       (NE)
            Integer,      Intent (In)    :: NENK
            Integer,      Intent (In)    :: NI
            Integer,      Intent (In)    :: offset   (NE + 1)
            Integer,      Intent (In)    :: SAoffset (NE + 1, NI)
            Integer,      Intent (In)    :: sl       (NE, NI)
            Real (KREAL), Intent (InOut) :: t        (:)
            Real (KREAL), Intent (InOut) :: s        (:)
        End Subroutine
        
    End Interface
    
    Contains
    
    Subroutine Destroy (self)
        Implicit None
        Class (UpstreamSys), Intent (InOut) :: self
        
        if (ALLOCATED(self%SAoffset)) DEALLOCATE(self%SAoffset)
        if (ALLOCATED(self%UFoffset)) DEALLOCATE(self%UFoffset)
        if (ALLOCATED(self%upfaces )) DEALLOCATE(self%upfaces )
        if (ALLOCATED(self%ipsf    )) DEALLOCATE(self%ipsf    )
        if (ALLOCATED(self%ipef    )) DEALLOCATE(self%ipef    )
        
    End Subroutine
    
End Module SpaceAngleInterface