Module LionboltHDF5
    use constants
    use IO
    use HDF5Tree
    use HDF5Attributes
    use HDF5Datasets
    use NittanyAPI
    use LionboltHDF5Schema
    use Geometry
    use HDF5
    Implicit None
    
Contains

Subroutine OpenParticleH5 (particle, nolc, noroot)
    ! By default, opens a particle group in the root.
    ! Providing nolc   = TRUE makes the group name NOT lowercase.
    ! Providing noroot = TRUE creates/finds the particle group in the current active group
    Implicit None
    Character (*),     Intent (In) :: particle
    Logical, Optional, Intent (In) :: nolc
    Logical, Optional, Intent (In) :: noroot
    
    Character (LEN=:), Allocatable :: Name
    Logical                        :: stat
    
    ! Return active node to root if not requested otherwise
    if (.not. PRESENT(noroot) .or. .not. noroot) call ReturnToRoot ()
    
    ! Determine the group name
    if (PRESENT(nolc) .AND. nolc) then
        ALLOCATE(Name, source = TRIM(ADJUSTL(particle)))
    else
        ALLOCATE(Name, source = TRIM(ADJUSTL(Lowercase(particle))))
    end if
    
    ! Create/find the group in the active group
    call OpenGroup (Name, stat)
    
    ! If this group was just now created, add an attribute
    ! explicitly specifying it as a particle group
    if (stat) then
        call AddAttribute (ATTR_OBJTYPE, IS_PARTICLE)
    end if
    
End Subroutine

Subroutine CloseParticleH5 (particle, nolc)
    Implicit None
    Character (*),     Intent (In) :: particle
    Logical, Optional, Intent (In) :: nolc
    
    Character (LEN=:), Allocatable :: Name
    
    ! Determine the group name
    if (PRESENT(nolc) .AND. nolc) then
        ALLOCATE(Name, source = TRIM(ADJUSTL(particle)))
    else
        ALLOCATE(Name, source = TRIM(ADJUSTL(Lowercase(particle))))
    end if
    
    ! CAN I THROW AN ERROR IF YOU ARE NEVER GOING TO ENCOUNTER PARTICLE???
    ! HOW???
    ! This may lead to infinite loops...
    do while (ActiveNode%Name /= Name)
        call RewindActiveNode (1)
    end do
    call CloseGroup (Name)
    
End Subroutine

Subroutine WriteProblemTypeH5 (probtype)
    ! Writes the problem type to the active group
    Implicit None
    Character (*), Intent (In) :: probtype
    
    call AddAttribute (ATTR_PROBTYPE, probtype)
    
End Subroutine

Subroutine WriteMeshH5 (mesh, noroot)
    ! Writes a mesh to the HDF5 file.
    ! If noroot = TRUE or noroot isn't given, it writes the mesh to current location
    Implicit None
    Type (MeshClass),          Intent (In) :: mesh
    
    Logical,         Optional, Intent (In) :: noroot
    
    Integer                                :: e, m
    Integer                                :: Nm
    Integer                                :: NE
    Integer                                :: NENK
    Integer                                :: NK
    Real (KREAL),              Allocatable :: vol (:) ! Temporary array for the volumes
    Character (512)                        :: matstrtmp
    Character (LEN=:),         Allocatable :: matstr
    
    Nm   = mesh%Nm
    NE   = mesh%NE
    NENK = mesh%NENK
    NK   = mesh%NK
    
    if (.not. PRESENT(noroot) .or. .not. noroot) call ReturnToRoot ()
    
    !  =========================
    !    Create the mesh group
    !  =========================
    
    call OpenGroup (GROUP_MESH)
    
    !  ========================================
    !    Add the problem type as an attribute  
    !  ========================================
    
    if (mesh%slab) then
        call AddAttribute (ATTR_PROBTYPE, 'slab')
    else
        call AddAttribute (ATTR_PROBTYPE, 'general')
    end if
    
    !  ===========================================
    !    Add number of materials as an attribute  
    !  ===========================================
    
    call AddAttribute (ATTR_NUMMATS_M, Nm)
    
    !  ==========================================
    !    Add number of elements as an attribute  
    !  ==========================================
    
    call AddAttribute (ATTR_NUMELS, NE)
    
    !  ==============================================
    !    Add number of spatial dofs as an attribute  
    !  ==============================================
    
    call AddAttribute (ATTR_NUMSD, NENK)
    
    !  ==============================================
    !    Add number of global nodes as an attribute  
    !  ==============================================
    
    call AddAttribute (ATTR_NUMKG, NK)
    
    !  ================================
    !    Write the connectivity array  
    !  ================================
    
    ! NOTES BECAUSE I KEEP FORGETTING:
    !   connectivity(sdof) gives the global node index of spatial dof sdof.
    !   offset(e) + k gives the spatial dof of element e and local node k.
    !   Thus, to get from (e,k) to kg, you would use connectivity(offset(e) + k).
    !   
    !   Therefore, I originally made connectivity as a flattened array corresponding to
    !   wiscobolt's 'Cekk' array that I carried around everywhere. This is standard practice I believe.
    !   Why haven't I needed connectivity anywhere in Lionbolt though? Is it because I flattened 
    !   my space-angle index as ip?
    
    ! Using one-based indexing (e.g. Fortran), connectivity(offset(e) + k) for e, k, starting from 1, gives the global index of node (e, k)
    ! Using zero-based indexing (e.g. Python), connectivity[offset[e] + k] - 1 for e, k, starting from 0, gives the global index of node (e, k)
    ! Here, we write using the one-based indexing convention. It is the responsibility of the user / Terpdose
    ! to adjust to the zero-based, or any other, convention
    
    call AddDataset (DATASET_CONNECTIVITY, mesh%connectivity)
    
    !  ==========================
    !    Write the offset array  
    !  ==========================
    
    call AddDataset (DATASET_OFFSET, mesh%offset)
    
    !  ===============================
    !    Write the global mesh nodes  
    !  ===============================
    
    call AddDataset (DATASET_R, mesh%rg)
    
    !  ==================
    !    Material nodes  
    !  ==================
    
    call OpenGroup (GROUP_MAT2SD)
    do m = 1, mesh%Nm
        
        write (matstrtmp, '(A,I0)') PATTERN_MATERIALS, m
        ALLOCATE(matstr, source = TRIM(ADJUSTL(matstrtmp)))
        
        call AddDataset (matstr, mesh%mat2sd(m)%v)
        
        DEALLOCATE(matstr)
        
    end do
    call CloseGroup (GROUP_MAT2SD)
    
    !  ===================
    !    Element Volumes  
    !  ===================
    
    ALLOCATE(vol(NE))
    do e = 1, NE
        vol(e) = mesh%el(e)%vol
    end do
    
    call AddDataset (DATASET_VOLS, vol)
    
    call CloseGroup (GROUP_MESH)
    
End Subroutine

Subroutine AppendOrdinatesH5 (w, k)
    ! Writes the abscissae and weights to the active group (expected to be a particle group)
    ! Also gives the number of angular dofs as an attribute
    Implicit None
    Real (KREAL), Intent (In) :: w (:)
    Real (KREAL), Intent (In) :: k (:,:)
    
    Integer                   :: NI
    
    NI = SIZE(w)
    
    call AddAttribute (ATTR_NUMANGLES,     NI)
    call AddDataset   (DATASET_ANGWEIGHTS, w)
    call AddDataset   (DATASET_ABSCISSAE,  k)
    
End Subroutine

Subroutine AppendEnergyGridH5 (E)
    Implicit None
    Real (KREAL), Intent (In) :: E (:)
    
    Integer                   :: G
    
    G = SIZE(E) - 1
    
    call AddAttribute (ATTR_NUMENERGIES, G)
    call AddDataset   (DATASET_EGRID,    E)
    
End Subroutine

Subroutine AppendSpatialArray (Name, x, g)
    ! Deprecated but could be useful to me in the future. 
    ! Thus it doesn't conform to typical form of AppendAngularFluence
    ! If you need some other similar functionality or functionalities, make your own
    ! routine using the HDF5Tree and HDF5 attribute/dataset routines.
    Implicit None
    Character (*),     Intent (In) :: Name
    Real (KREAL),      Intent (In) :: x (:)
    Integer, Optional, Intent (In) :: g          ! Energy group
    
    Character (32)                 :: tmp
    Character (LEN=:), Allocatable :: gstr
    Character (LEN=:), Allocatable :: grouplabel
    
    if (PRESENT(g)) then
        ! If g is present this implies that the Name should be a new group
        ! in which g indexes datasets.
        
        ! Open / create the group
        call OpenGroup (Name)
        
        ! Write the energy group to a character
        write (tmp, '(I0)') g
        ALLOCATE(gstr, source = TRIM(tmp))
        
        ! Get the group label by combining the energy pattern and gstr
        ALLOCATE(grouplabel, source = PATTERN_ENERGY // gstr)
        
        call AddDataset (grouplabel, x)
        
        call CloseGroup (Name)
        
    else
        ! If g is not present, Name is just the name of the dataset added to the active group.
        ! User should probably just use AddDataset... Must rethink a lot of the user-friendliness of HDF5 especially considering some of it lies in core.
        
        call AddDataset (Name, x)
        
    end if
    
End Subroutine

Subroutine AppendAngularFluence (x, NENK, uncollided, g)
    Implicit None
    Real (KREAL),      Intent (In) :: x (:)
    Integer,           Intent (In) :: NENK       ! Size of an individual spatial array (any other way to get this?)
    Logical,           Intent (In) :: uncollided ! Whether or not this is an uncollided angular fluence
    Integer,           Intent (In) :: g          ! Energy group
    
    Integer                        :: i
    Integer                        :: NI
    Real (KREAL)                   :: t (NENK)   ! Temporary space array for writing
    Character (32)                 :: tmp
    Character (LEN=:), Allocatable :: Name
    Character (LEN=:), Allocatable :: gstr
    Character (LEN=:), Allocatable :: grouplabel
    Character (LEN=:), Allocatable :: istr
    Character (LEN=:), Allocatable :: dsetlabel
    
    ! Determine the type of angular fluence group
    if (uncollided) then
        ALLOCATE(Name, source = GROUP_ANG_FL_UNC)
    else
        ALLOCATE(Name, source = GROUP_ANG_FL)
    end if
    
    ! Open / create the angular fluence group
    call OpenGroup (Name)
    
    ! Write the energy group to a character
    write (tmp, '(I0)') g
    ALLOCATE(gstr, source = TRIM(tmp))
    
    ! Get the group label by combining the energy pattern and gstr
    ALLOCATE(grouplabel, source = PATTERN_ENERGY // gstr)
    
    ! Create the group. (there's no circumstance where this should already exist though)
    call OpenGroup (grouplabel)
    
    ! Loop over angles and append the space arrays
    NI = SIZE(x) / NENK
    do i = 1, NI
        ! Write the angle index to a character
        write (tmp, '(I0)') i
        ALLOCATE(istr, source = TRIM(tmp))
        
        ! Get the group label by combining the angle pattern and istr
        ALLOCATE(dsetlabel, source = PATTERN_ANGLE // istr)
        
        ! Add a dataset for this angle
        t(1:NENK) = x(1 + (i - 1) * NENK : i * NENK)
        call AddDataset (dsetlabel, t)
        
        DEALLOCATE(istr)
        DEALLOCATE(dsetlabel)
        
        ! I could add attributes here
        
    end do
    
    ! Close the groups
    call CloseGroup (grouplabel)
    call CloseGroup (Name)
    
End Subroutine

Subroutine WriteXSLibraryLionbolt (XSLib)
    ! Essentially a wrapper for NittanyPhysics' HDF5 format, but it will first
    ! create a group so that everything that NittanyPhysics does is in that group.
    Implicit None
    Type (XSLibrary), Intent (In) :: XSLib (:,:)
    
    call OpenGroup (GROUP_XSLIBRARY)
    call WriteXSLibrary (XSLib)
    
End Subroutine

End Module