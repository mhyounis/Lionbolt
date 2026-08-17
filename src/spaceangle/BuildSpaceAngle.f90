Submodule (SpaceAngleInterface) submodBuildSpaceAngle
Contains
Module Subroutine BuildSpaceAngle (mesh, angular, sl, isl, IAF)
    use constants
    use profiler
    use types
    use Geometry
    use AngularSpace
    use SpaceAngleInterface
    Implicit None
    Type (MeshClass),                 Intent (In)  :: mesh
    Type (AngularClass),              Intent (In)  :: angular
    
    Integer,             Allocatable, Intent (Out) :: sl  (:,:)
    Integer,             Allocatable, Intent (Out) :: isl (:,:)
    Type (UpstreamSys),               Intent (Out) :: IAF
    
    Logical,             Allocatable               :: IsUpstream (:)
    Integer,             Allocatable               :: IUoffset   (:,:)
    Integer,             Allocatable               :: estart     (:)
    
    ! if (.not. mesh%slab) call profile ('Time to determine sweep order ')
    if (gmemprofiling) call RSSLogger (' --- IN BUILDSPACEANGLE --- Before sweep')
    call SweepOrder (mesh, angular, IsUpstream, IUoffset, sl, isl, estart)
    if (gmemprofiling) call RSSLogger (' --- IN BUILDSPACEANGLE --- After sweep')
    ! if (.not. mesh%slab) call profile ('Time to determine sweep order ')
    if (gmemprofiling) call RSSLogger (' --- IN BUILDSPACEANGLE --- Before sweeparrays')
    call SweepArrays (mesh, angular, IsUpstream, IUoffset, sl, isl, IAF)
    if (gmemprofiling) call RSSLogger (' --- IN BUILDSPACEANGLE --- After sweeparrays')
    
End Subroutine
End Submodule