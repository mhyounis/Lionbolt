Submodule (Geometry) submodReadGMSH
Contains
Module Subroutine ReadGMSH (fname, mesh)
    use LionboltBasePath
    use constants
    use IO
    use Geometry
    Implicit None
    Character (*),     Intent (In)    :: fname
    
    Type (MeshClass),  Intent (InOut) :: mesh
    
    Integer                           :: b, e, f, i, j, k
    Logical                           :: needhomog
    Logical                           :: deletescratch = .FALSE.
    Integer,           Parameter      :: unit = 10
    Integer                           :: iPyerr
    Integer                           :: ierr
    Integer                           :: lNF
    Integer                           :: lNK
    Integer                           :: NvB
    Integer                           :: Niter
    Integer                           :: e1
    Integer                           :: e2
    Integer                           :: b1
    Integer                           :: b2
    Integer                           :: m1
    Integer                           :: m2
    Integer                           :: reg1
    Integer                           :: reg2
    Integer                           :: ints    (5) ! MAY HAVE TO GENERALIZE THIS
    Integer,           Allocatable    :: volinfo (:)
    Real (KREAL)                      :: reals   (6) ! MAY HAVE TO GENERALIZE THIS
    Character (LEN=:), Allocatable    :: line
    Character (LEN=:), Allocatable    :: shfname
    
    ! https://gmsh.info//doc/texinfo/gmsh.html#MSH-ASCII-file-format
    ! https://gmsh.info/dev/doc/texinfo/gmsh.pdf
    
    ! Get the shell base path and file name
    call GetShellPath (LBpath)
    shfname = WrapFolders (fname)
    
    if (.not. ALLOCATED(ScratchName)) then
        ! Default scratch is just gonna be in the Lionbolt folder
        ! Probably a bad construct. Just know that gmshprocessing.py is on the chopping block
        deletescratch = .TRUE.
        ScratchName = LBpath
    end if
    
    ! Pre-process the mesh using the Python API
    ! write (iuout,*) 'Running script gmshprocessing.py'
    ! MHY LATER - Must end all dependence on executing a python script externally.
    ! AND/OR handle scratch stuff... driver users don't set a scratch...
    ! This will be more important during MPI update.
    call system ('python ' // shpath // '/src/geometry/gmshprocessing.py -b ' &
                 // shpath // ' -f ' // shfname // ' -s ' // ScratchName // '> /dev/null 2>&1', iPyerr)
    if (iPyerr /= 0) then
        call stophere ('ReadGMSH: gmshprocessing.py failed')
    end if
    ! write (iuout,*) 'Script terminated successfully'
    
    ! NOW READ IN ALL QUANTITIES.
    ! Read in nodes
    open(unit=unit, file=ScratchName // 'TAPE%0_raw_nodes.txt', action='read')
    read(unit,*) ints(1)
    
    mesh%NK = ints(1)
    ALLOCATE(mesh%rg(3, mesh%NK))
    do k = 1, mesh%NK
        read(unit,*) ints(1), mesh%rg(1:3,k)
    end do
    close(unit=unit, status='delete')
    
    ! Read in elements. Must first go through the raw GMSH file, to abstract materials/beam information.
    mesh%Nm = 0
    
    open(unit=unit, file=fname, action='read')
    do
        call ReadLine (unit, line, ierr)
        
        if (ierr < 0) exit
        
        if (line == '$Entities') then
            ! Number of points in geo, Number of curves, Number of surfaces, Number of volumes 
            read(unit,*) ints(1:4)
            
            NvB = ints(4)
            
            ! MAKE ARRAY(S) WITH VOL INFO, FOR USE IN ELEMENTS b.
            ALLOCATE(volinfo(NvB))
            volinfo   = 0
            needhomog = .FALSE.
            
            ! volinfo should never need more than 1. If it does, that means a b has multiple materials, so call an error.
            
            ! Skip to volumes (for now, may have use for surface info later)
            do i = 1, SUM(ints(1:3))
                call ReadLine (unit, line, ierr)
            end do
            
            do b = 1, NvB
                call ReadLine (unit, line, ierr)
                
                ! Read the number of physical tags in this volume
                read(line,*) ints(1), reals(1:6), ints(2)
                
                ! The number of physical tags for a single volume should never exceed two
                ! (One for material, one for beam specification)
                if (ints(2) > 2) then
                    call stophere ('ReadGMSH: A volume in the user-provided mesh file has more than' // &
                                   ' two physical tags. There should never be more than two (one for a material,' // &
                                   ' one for the beam).')
                end if
                
                ! Once you know the number, you can determine the actual physical tags present
                read(line,*) ints(1), reals(1:6), ints(2), ints(3:2 + ints(2))
                
                ! Assign the physical tags to volinfo
                if (ints(2) == 0) then
                    ! If there are no physical tags in this material, you assume homogeneous transport
                    volinfo(b) = 1
                    
                    needhomog = .TRUE.
                else if (ints(2) == 1) then
                    volinfo(b) = ints(3)
                else
                    call stophere ('ReadGMSH: A volume in the user-provided mesh has more than '//     &
                                   'one physical tag. There should never be more than two (as this '// &
                                   'suggests two materials in the same volume).')
                end if
            end do
            
            ! If materials were not provided for any of the volumes, then check that they were
            ! not provided for ALL of the volumes
            if (needhomog) then
                if (ANY(volinfo /= 1)) then
                    call stophere ('ReadGMSH: User-provided mesh does not specify material for ALL volumes.')
                end if
            end if
            
        end if
    end do
    close(unit=unit)
    
    ! Read raw elements file. First correlate the bs to materials.
    open(unit=unit, file=ScratchName // 'TAPE%0_elements.txt', action='read')
    read(unit,*) ints(1)
    
    mesh%NE = ints(1)
    
    ALLOCATE(mesh%el(mesh%NE))
    
    e = 0
    mesh%Nm = 0
    do b = 1, NvB
        read(unit,*) ints(1:2)
        Niter = ints(2)
        
        select case (ints(1))
        case (4)
            lNK = 4
            lNF = 4
        case default
            call stophere ('Error reading mesh: Currently only 4-node tetrahedra are implemented.')
        end select
        
        do i = 1, Niter
            e = e + 1
            
            mesh%el(e)%block = b
            mesh%el(e)%idx   = ints(1)
            mesh%el(e)%NK    = lNK
            mesh%el(e)%NF    = lNF
            mesh%el(e)%mat   = volinfo(b) ! Make sure that user-defined physical tags don't skip any numbers
            mesh%Nm          = MAX(mesh%Nm, mesh%el(e)%mat)
            
            ALLOCATE(mesh%el(e)%node(1:lNK))
            
            read(unit,*) ints(2), mesh%el(e)%node(1:lNK) ! ints(2) is trash here
            
        end do
        
    end do
    close(unit=unit, status='delete')
    
    ! Some quick allocations
    do e = 1, mesh%NE
        lNF = mesh%el(e)%NF
        ALLOCATE(mesh%el(e)%neighbors(1:lNF))
    end do
    do e = 1, mesh%NE
        lNF = mesh%el(e)%NF
        ALLOCATE(mesh%el(e)%face(1:lNF))
    end do
    do e = 1, mesh%NE
        lNF = mesh%el(e)%NF
        ALLOCATE(mesh%el(e)%n(3, 1:lNF))
    end do
    
    ! Read in local face information
    open(unit=unit, file=ScratchName // 'TAPE%0_local_faces.txt', action='read')
    do e = 1, mesh%NE
        read(unit,*) ints(1), mesh%el(e)%face(1:mesh%el(e)%NF)
    end do
    close(unit=unit, status='delete')
    
    ! Read in face elements and nodes
    open(unit=unit, file=ScratchName // 'TAPE%0_global_faces.txt', action='read')
    read(unit,*) ints(1:2)
    
    mesh%NF  = ints(1)
    mesh%NBF = ints(2)
    
    ALLOCATE(mesh%face(mesh%NF))
    do f = 1, mesh%NBF
        mesh%face(f)%bdy = .TRUE.
    end do
    
    do f = 1, mesh%NF
        read(unit,*) ints(1:3)
        mesh%face(f)%NK  = ints(1)
        mesh%face(f)%idx = ints(1) - 1 ! Index is just NK minus 1. So if NK = 3 (triangle), index is 2.
        
        mesh%face(f)%sharedby(1:2) = ints(2:3)
        
        ALLOCATE(mesh%face(f)%node(mesh%face(f)%NK))
        
        read(unit,*) mesh%face(f)%node
    end do
    close(unit=unit, status='delete')
    
    ! Read in element neighbors
    open(unit=unit, file=ScratchName // 'TAPE%0_neighbors.txt', action='read')
    do e = 1, mesh%NE
        read(unit,*) ints(1), mesh%el(e)%neighbors(1:mesh%el(e)%NF)
    end do
    close(unit=unit, status='delete')
    
    ! Determine if all elements have the same structure
    mesh%homog = .TRUE.
    do e = 1, mesh%NE - 1
        mesh%homog = mesh%homog .AND. mesh%el(e)%idx == mesh%el(e + 1)%idx
    end do
    
    ALLOCATE(mesh%offset(mesh%NE + 1))
    ! mesh%offset(1) = 1
    ! do e = 1, mesh%NE
    !     mesh%offset(e + 1) = mesh%offset(e) + mesh%el(e)%NK
    ! end do
    
    j = 0
    do e = 1, mesh%NE
        mesh%offset(e) = j
        j = j + mesh%el(e)%NK
    end do
    mesh%offset(mesh%NE + 1) = j
    
    ! Determine region info
    if (mesh%Nm == 1) then
        ! For one material there is just one region, properties are all known easily
        ALLOCATE(mesh%region(1))
        ALLOCATE(mesh%region(1)%fg, source=PACK([(f, f=1, mesh%NBF)], mask=.TRUE.))
        mesh%region(1)%mat = 1
        mesh%region(1)%NF  = mesh%NF
        ! Block info is not needed
    else
        ! First determine block connectivity. This will allow you to map from block to region later
        ! You start with b = reg, and then just re-define regions as you connect things
        ALLOCATE(mesh%b2reg, source=PACK([(b, b=1, NvB)], mask=.TRUE.))
        
        ! Boundary faces give no useful info
        do f = mesh%NBF + 1, mesh%NF
            e1 = mesh%face(f)%sharedby(1)
            e2 = mesh%face(f)%sharedby(2)
            b1 = mesh%el(e1)%block
            b2 = mesh%el(e2)%block
            m1 = mesh%el(e1)%mat
            m2 = mesh%el(e2)%mat
            
            ! If b1 and b2 are different, these two blocks are connected
            if (b1 /= b2) then
                ! They may be in the same region, as long as they have the same material
                if (m1 == m2) then
                    ! These two blocks are in the same region.
                    ! First, make the new union region have the index of the lower-indexed region
                    reg1 = mesh%b2reg(b1)
                    reg2 = mesh%b2reg(b2)
                    
                    if (reg2 > reg1) then
                        ! Replace reg2 with reg1
                        mesh%b2reg(b2) = reg1
                    else
                        ! Replace reg1 with reg2
                        mesh%b2reg(b1) = reg2
                    end if
                    
                    ! This convention makes the region indexes go from 1 to however many regions there are naturally
                    
                end if
            end if
            
        end do
        
        ALLOCATE(mesh%region(MAXVAL(mesh%b2reg)))
        do reg1 = 1, MAXVAL(mesh%b2reg)
            ALLOCATE(mesh%region(reg1)%fg(0)) ! May be concerned about data storage location and stuff
            
            b = FINDLOC(mesh%b2reg, reg1, DIM=1) ! Find a block to assign the material
            mesh%region(reg1)%mat = volinfo(b)
        end do
        
    end if
    
    if (deletescratch) DEALLOCATE(ScratchName)
    
End Subroutine
End Submodule


! Algorithm plans ---
!   Make local faces while reading. Follow the GMSH convention for local node indexing. May also have to have a sorted set of faces that is only used for connectivity
!   connectivity - Use sorted faces. Have a separate array of faces for each possible number of nodes on a face (because only these will match). Match indices successively.
!   