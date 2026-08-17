Module Geometry
    use constants
    use types
    Implicit None
    
    !  ==================================================================================================
    !    A set of types, most importantly, the "MeshClass" which stores all of the relevant information  
    !    about a user-defined mesh                                                                       
    !  ==================================================================================================
    
    Type InnerProductType
        ! All inner products listed here are mass-eliminated.
        ! i.e., they have the inverse of the mass integral applied
        ! Also note, we still must figure out generalization for when meshes of multiple element types
        ! are introduced. We will likely make the IFM and IIM matrices be face-shape dependent (i.e., tri, rect, etc.),
        ! however we need to then figure out permutations to map these to full element matrices.
        Integer ,     Allocatable :: elmap (:)       ! Takes the mesh element and gives the index 'le.' Returns 0 if the element isn't of this type
        Real (KREAL), Allocatable :: IC    (:,:,:,:) ! Tetrahedral, convective            ! (k, k', dir, le)
        Real (KREAL), Allocatable :: IFM   (:,:,:,:) ! Tet-Tet Face mass                  ! (k, k', f, le)
        Real (KREAL), Allocatable :: IIM   (:,:,:,:) ! Tet-Tet Interface mass, transposed ! (k, k', f, le)
    Contains
        Procedure :: Destroy => DestroyIP
    End Type
    
    Type MeshClass
        Logical                              :: init  = .FALSE.    ! Just a small boolean to prevent a mesh from being formed multiple times, or from adding a slab to a 3D mesh or vice versa
        Logical                              :: READY = .FALSE.    ! For the case of a 3D mesh, because one can translate and scale their mesh after reading, they need to manually post-process when they are ready to essentially publish the mesh. That sets READY = .TRUE.
        Logical                              :: slab               ! True if slab, false if full 3D geometry
        Logical                              :: homog              ! Whether or not all elements have the same structure (i.e., all tetrahedra, all prisms, etc.)
        Integer                              :: Nm = 0             ! Number of materials
        Integer                              :: NK                 ! Number of global nodes
        Integer                              :: NF                 ! Number of global faces
        Integer                              :: NBF                ! Number of boundary global faces
        Integer                              :: NE = 0             ! Number of elements
        Integer                              :: NENK               ! The total number of element-nodes (i.e., \sum_{e}mesh%el(e)%NK)
        Integer,                 Allocatable :: b2reg        (:)   ! For multi-material meshes this gives the 'region' index from the block index (block is recorded with elements)
        Integer,                 Allocatable :: offset       (:)   ! Offset array for flattened spatial dof index
        Integer,                 Allocatable :: connectivity (:)   ! Connectivity array for flattened spatial dof index
        Real (KREAL),            Allocatable :: rg           (:,:) ! Global mesh      ! (dir, kg)
        Type (InnerProductType)              :: IP           (7)   ! (idx)
        Type (CRealM),           Allocatable :: IMr          (:)   ! The inverse finite element mass matrix in a reference element of index idx ! (idx)%m(k, k') 
        Type (CIntV),            Allocatable :: mat2sd       (:)   ! The list of sdofs corresponding to material m ! (m)
        Type (MeshFace),         Allocatable :: face         (:)   ! (f)
        Type (MeshElement),      Allocatable :: el           (:)   ! (e)
        Type (MeshRegionType),   Allocatable :: region       (:)   ! For multi-material meshes this is used in raytracing.
        ! Type (MeshOctree)                    :: octree ! WIP but will be very important in the near future
    Contains
        Procedure :: FromFile
        Procedure :: Translate
        Procedure :: Scale
        Procedure :: AddSlab
        Procedure :: PostProcess
        Procedure :: Destroy => DestroyMesh
    End Type
    
    Type MeshElement
        Integer                   :: block           ! In the mesh, the 'block' or volume tag this element belongs to
        Integer                   :: NK              ! Number of nodes in the element
        Integer                   :: NF              ! Number of faces in the element
        Integer                   :: order           ! Element order (linear = 1, quadratic = 2, etc.) NOT CURRENTLY IN USE
        Integer                   :: mat             ! Index of the element's material
        Integer                   :: idx             ! Index giving the type of finite element (according to GMSH's indexing system)
        Integer,      Allocatable :: node      (:)   ! Index of a global node that composes the element
        Integer,      Allocatable :: face      (:)   ! Index of a global face that composes the element                ! (f)
        Integer,      Allocatable :: neighbors (:)   ! List of elements neighboring this one
        Real (KREAL)              :: vol             ! Volume of element
        Real (KREAL)              :: detJv           ! Determinant of Jacobian ! Sure it's 6 * V for a tetrahedron but not all elements are tetrahedra
        Real (KREAL), Allocatable :: n         (:,:) ! Normal vectors of the faces                                     ! (dir, f)
    End Type
    
    Type MeshFace
        Logical                  :: bdy = .FALSE. ! Whether or not the face is on the boundary
        Integer                  :: NK            ! Number of nodes on this face
        Integer                  :: idx           ! Index giving the type of face (the shape, according to GMSH's indexing system)
        Integer                  :: sharedby (2)  ! The elements on either side of the face
        Integer,     Allocatable :: node     (:)  ! Index of a global node that composes the face
        Real (KREAL)             :: area          ! Surface area of face
    End Type
    
    ! Type InnerProductType
    !     ! All inner products listed here are mass-eliminated.
    !     ! i.e., they have the inverse of the mass integral applied
    !     Real (KREAL)               :: detJv        ! Determinant of volume Jacobian
    !     Real (KREAL),  Allocatable :: IC   (:,:,:) ! Convective     ! (k, k', dir)
    !     Type (CRealM), Allocatable :: IFM  (:)     ! Face mass      ! (f)%(k, k')
    !     Type (CRealM), Allocatable :: IIM  (:)     ! Interface mass ! (f)%(k, k')
    ! End Type
    
    Type MeshRegionType
        Integer              :: NF     ! Number of global faces enclosing the region
        Integer              :: mat    ! Material of this region
        Integer, Allocatable :: b  (:) ! List of blocks in this region
        Integer, Allocatable :: fg (:) ! List of global faces enclosing the region
    End Type
    
    ! Type LeafType
    !     Integer,           Allocatable :: elements (:) ! List of elements in this leaf
    !     Real (KREAL)                   :: xmin         ! Min x value of leaf box
    !     Real (KREAL)                   :: xmax         ! Max x value of leaf box
    !     Real (KREAL)                   :: ymin         ! Min y value of leaf box
    !     Real (KREAL)                   :: ymax         ! Max y value of leaf box
    !     Real (KREAL)                   :: zmin         ! Min z value of leaf box
    !     Real (KREAL)                   :: zmax         ! Max z value of leaf box
    !     ! Does this type need a parent? If so, use Type (LeafType), Pointer :: Parent => NULL ()
    !     Type (LeafTypePW), Allocatable :: Children (:)
    ! Contains
    !     Procedure :: Split
    !     Procedure :: Populate
    !     Procedure :: PointIsInLeaf
    ! End Type
    ! 
    ! Type LeafTypePW
    !     Type (LeafType), Pointer :: p => NULL()
    ! End Type
    ! 
    ! Type MeshOctree
    !     ! WAIT - if root is the only interesting part of this then I really dont need this type.
    !     ! Because xmin xmax etc should agree with root. Maybe this is useful for type-bound procedures?
    !     ! Also how will I easily traverse the leaf? I can't use like getattr. But what I can do is have, e.g.,
    !     ! a function where I feed like 'n' to go down n levels... figure it out when I have a specific use case.
    !     Real (KREAL)    :: xmin ! Min x value of mesh
    !     Real (KREAL)    :: xmax ! Max x value of mesh
    !     Real (KREAL)    :: ymin ! Min y value of mesh
    !     Real (KREAL)    :: ymax ! Max y value of mesh
    !     Real (KREAL)    :: zmin ! Min z value of mesh
    !     Real (KREAL)    :: zmax ! Max z value of mesh
    !     Type (LeafType) :: root ! Highest level leaf
    ! End Type
    
    Private :: FromFile, ReadGMSH, PostProcess, AddSlab, init_MeshClass, Translate, Scale, DestroyMesh, DestroyIP
    
    Interface MeshClass
        Module Procedure init_MeshClass
    End Interface
    
    Interface
        
        Module Function init_MeshClass (fname) Result (self)
            Implicit None
            Character (*), Intent (In) :: fname
            Type (MeshClass)           :: self
        End Function
        
        Module Subroutine FromFile (self, fname)
            Implicit None
            Class (MeshClass), Intent (InOut) :: self
            Character (*),     Intent (In)    :: fname
        End Subroutine
        
        Module Subroutine ReadGMSH (fname, mesh)
            use constants
            use IO
            Implicit None
            Character (*),    Intent (In)    :: fname
            Type (MeshClass), Intent (InOut) :: mesh
        End Subroutine
        
        Module Subroutine PostProcess (self)
            use constants
            use IO
            Implicit None
            Class (MeshClass), Intent (InOut) :: self
        End Subroutine
        
        Module Subroutine AddSlab (self, NE, T, mat, structure)
            use constants
            Implicit None
            Class (MeshClass),       Intent (InOut) :: self
            Integer,                 Intent (In)    :: NE
            Real (KREAL),            Intent (In)    :: T
            Integer,       Optional, Intent (In)    :: mat
            Character (*), Optional, Intent (In)    :: structure
        End Subroutine
        
        Module Subroutine FiniteElementAnalysis (mesh)
            use constants
            use types
            use ShapeFunctions
            use FEAInnerProducts
            use Jacobians
            Implicit None
            Type (MeshClass), Intent (InOut) :: mesh
        End Subroutine
        
    End Interface
    
    Contains
    
    ! Translate and scale are such simple procedures that they really don't warrant
    ! entirely new files.
    
    ! MHY LATER - for main program make it so that depending on which key the user inputs first,
    ! the translate or the scale happens first. And they can be chained.
    
    Subroutine Translate (self, t)
        use constants
        Implicit None
        Class (MeshClass), Intent (InOut) :: self
        
        Real (KREAL),      Intent (In)    :: t (3)
        
        Integer                           :: i
        
        if (.not. self%init) then
            ! write (iuout,*) 
            ! Figure out what I should do here. Return? Stop? Write to iuout? Print?
            return
        end if
        
        if (self%READY) then
            call stophere ('Geometry.f90: Translate: Mesh translation was requested but this mesh is already READY. ' // &
                           'Destroy this object before reusing.')
        end if
        
        do i = 1, 3
            self%rg(i,:) = self%rg(i,:) + t(i)
        end do
        
    End Subroutine
    
    Subroutine Scale (self, s)
        use constants
        Implicit None
        Class (MeshClass), Intent (InOut) :: self
        
        Real (KREAL),      Intent (In)    :: s (3)
        
        Integer                           :: i
        
        if (.not. self%init) then
            ! write (iuout,*) 
            ! Figure out what I should do here. Return? Stop? Write to iuout? Print?
            return
        end if
        
        if (self%READY) then
            call stophere ('Geometry.f90: Scale: Mesh scaling was requested but this mesh is already READY. ' // &
                           'Destroy this object before reusing.')
        end if
        
        do i = 1, 3
            self%rg(i,:) = self%rg(i,:) * s(i)
        end do
        
    End Subroutine
    
    Subroutine DestroyMesh (self)
        Implicit None
        Class (MeshClass), Intent (InOut) :: self
        
        Integer                           :: idx
        
        ! Reset to defaults and deallocate
        
        self%init = .FALSE.
        self%Nm   = 0
        self%NE   = 0
        if (ALLOCATED(self%offset))       DEALLOCATE (self%offset)
        if (ALLOCATED(self%connectivity)) DEALLOCATE (self%connectivity)
        if (ALLOCATED(self%rg))           DEALLOCATE (self%rg)
        if (ALLOCATED(self%IMr))          DEALLOCATE (self%IMr)
        if (ALLOCATED(self%mat2sd))       DEALLOCATE (self%mat2sd)
        if (ALLOCATED(self%face))         DEALLOCATE (self%face)
        if (ALLOCATED(self%el))           DEALLOCATE (self%el)
        do idx = 1, 7
            call self%IP(idx)%Destroy ()
        end do
        self%READY = .FALSE.
        
    End Subroutine
    
    Subroutine DestroyIP (self)
        Implicit None
        Class (InnerProductType), Intent (InOut) :: self
        
        if (ALLOCATED(self%elmap)) DEALLOCATE(self%elmap)
        if (ALLOCATED(self%IC))    DEALLOCATE(self%IC)
        if (ALLOCATED(self%IFM))   DEALLOCATE(self%IFM)
        if (ALLOCATED(self%IIM))   DEALLOCATE(self%IIM)
        
    End Subroutine
    
End Module Geometry