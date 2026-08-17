Submodule (Geometry) submodFiniteElementAnalysis
Contains
Module Subroutine FiniteElementAnalysis (mesh)
    
    ! This all will need to be restructured when mesh elements can be mixed with eachother
    ! Specifically, the IFM storage corresponding to the full element idx (i.e., not a face idx or anything)
    ! is not really logical, but it's being done because I don't want to rearrange matrices using the permutation
    ! matrices that would be required. You can define an inner product based only on the face, but to apply it
    ! in matrix form you need to actually commit to the two element idx's on either side of the face, which
    ! is an issue because face shapes don't correspond one:one with full element shapes.
    ! In general just need a more unified treatment of ALL element shape and dimension types.
    ! 
    ! Anyway changes are coming
    
    use constants
    use types
    use ShapeFunctions
    use FEAInnerProducts
    use Jacobians
    use Geometry
    Implicit None
    Type (MeshClass),        Intent (InOut) :: mesh
    
    Integer                                 :: e, ep, f, fg, idx, idxf, le
    Logical                                 :: idxp    (7)   ! Whether or not certain mesh element types are present
    Integer                                 :: NE
    Integer                                 :: d
    Integer                                 :: NK
    Integer                                 :: NKp
    Integer                                 :: NF
    Integer                                 :: info
    Integer                                 :: nwork
    Integer                                 :: NKidx   (7)   ! Number of nodes in an element of GMSH index idx
    Integer                                 :: NFidx   (7)   ! Number of faces in an element of GMSH index idx
    Real (KREAL)                            :: detJv         ! Determinant of Jacobian matrix
    Real (KREAL)                            :: normJs
    Real (KREAL)                            :: Js      (3)
    Real (KREAL),            Allocatable    :: work    (:)
    Real (KREAL)                            :: Jv      (3,3) ! Jacobian matrix
    Real (KREAL)                            :: JvinvT  (3,3) ! Inverse-Transpose of Jacobian matrix
    Type (CIntV)                            :: ipiv    (7)   ! ipiv for BLAS routine
    Type (CIntVV)                           :: f2k     (7)   ! Gives the element node index of a face and ref node (ref node is a node labeled in the face)
    Type (CRealM)                           :: IMr     (7)   ! Reference mass integral / inverse info. Entry idx gives for mesh element type idx
    Type (CRealM)                           :: IM      (7)   ! IMr / detJv for the given element
    Type (CRealM)                           :: IFMr    (0:2) ! Reference face mass integral / inverse info. idx = 0 : line "faces", idx = 1 : triangle, idx = 2 : rectangle
    Type (CRealT3)                          :: ICr     (7)   ! Reference convective integral / inverse info. Entry idx gives for mesh element type idx
    
    !  ======================================
    !    Briefly, initialize FEA quadrature
    !  ======================================
    
    call InitializeFEAQuadrature (4)
    
    !  =============================================
    !    Determine which element types are present
    !  =============================================
    
    NE = mesh%NE
    
    idxp = .FALSE.
    do e = 1, NE
        idxp(mesh%el(e)%idx) = .TRUE.
    end do
    
    !  ===============================================================================
    !    Allocations and assignments. Here, mainly temporaries and addressing arrays
    !  ===============================================================================
    
    NKidx = 0
    NKidx(1) = 2
    NKidx(4) = 4
    NKidx(5) = 8
    NKidx(6) = 6
    NKidx(7) = 5
    
    NFidx = 0
    NFidx(1) = 2
    NFidx(4) = 4
    NFidx(5) = 6
    NFidx(6) = 5
    NFidx(7) = 5
    
    do idx = 1, 7
        if (.not. idxp(idx)) cycle
        ALLOCATE(ipiv(idx)%v(NKidx(idx)))
    end do
    
    ! Linear 1D elements - all faces (the two edges) have one node
    ALLOCATE(f2k(1)%v(NFidx(1)))
    do f = 1, NFidx(1)
        ALLOCATE(f2k(1)%v(f)%v(1))
    end do
    
    ! Tetrahedra - all faces have same number of nodes (3)
    ALLOCATE(f2k(4)%v(NFidx(4)))
    
    do f = 1, NFidx(4)
        ALLOCATE(f2k(4)%v(f)%v(3))
    end do
    
    ! Hexahedra - all faces have same number of nodes (4)
    ALLOCATE(f2k(5)%v(NFidx(5)))
    
    do f = 1, NFidx(5)
        ALLOCATE(f2k(5)%v(f)%v(4))
    end do
    
    ! Prism - first two faces (by GMSH convention) have 3 nodes, last three have 4
    ALLOCATE(f2k(6)%v(NFidx(6)))
    
    do f = 1, 2
        ALLOCATE(f2k(6)%v(f)%v(3))
    end do
    do f = 3, NFidx(6)
        ALLOCATE(f2k(6)%v(f)%v(4))
    end do
    
    ! Pyramid - first face (by GMSH convention) has 4 nodes, last four have 3
    ALLOCATE(f2k(7)%v(NFidx(7)))
    
    ALLOCATE(f2k(7)%v(1)%v(4))
    do f = 2, NFidx(7)
        ALLOCATE(f2k(7)%v(f)%v(3))
    end do
    
    f2k(1)%v(1)%v = [1]
    f2k(1)%v(2)%v = [2]
    
    f2k(4)%v(1)%v = [1, 2, 3]
    f2k(4)%v(2)%v = [1, 2, 4]
    f2k(4)%v(3)%v = [1, 3, 4]
    f2k(4)%v(4)%v = [2, 3, 4]
    
    f2k(5)%v(1)%v = [1, 2, 3, 4]
    f2k(5)%v(2)%v = [5, 6, 7, 8]
    f2k(5)%v(3)%v = [1, 2, 6, 5]
    f2k(5)%v(4)%v = [3, 4, 8, 7]
    f2k(5)%v(5)%v = [1, 4, 8, 5]
    f2k(5)%v(6)%v = [2, 3, 7, 6]
    
    f2k(6)%v(1)%v = [1, 2, 3]
    f2k(6)%v(2)%v = [4, 5, 6]
    f2k(6)%v(3)%v = [1, 2, 5, 4]
    f2k(6)%v(4)%v = [2, 3, 6, 5]
    f2k(6)%v(5)%v = [3, 1, 4, 6]
    
    f2k(7)%v(1)%v = [1, 2, 3, 4]
    f2k(7)%v(2)%v = [1, 2, 5]
    f2k(7)%v(3)%v = [2, 3, 5]
    f2k(7)%v(4)%v = [3, 4, 5]
    f2k(7)%v(5)%v = [4, 1, 5]
    
    !  ===================================================================
    !    Form the reference inner products for all present element types
    !  ===================================================================
    
    ! Make all face mass ref IP's, even if they won't be used
    call FEAFaceMassRef (1, IFMr(0)%m)
    call FEAFaceMassRef (3, IFMr(1)%m)
    call FEAFaceMassRef (4, IFMr(2)%m)
    
    do idx = 1, 7
        if (.not. idxp(idx)) cycle
        
        NK = NKidx(idx)
        
        call FEAMassRef       (idx, IMr(idx)%m)
        call FEAConvectionRef (idx, ICr(idx)%t)
        
        ! Invert the mass reference inner product
        call dgetrf (NK, NK, IMr(idx)%m, NK, ipiv(idx)%v, info)
        
        nwork = -1
        ALLOCATE(work(1))
        call dgetri (NK, IMr(idx)%m, NK, ipiv(idx)%v, work, nwork, info)
        nwork = INT(work(1))
        DEALLOCATE(work)
        ALLOCATE(work(nwork))
        
        call dgetri (NK, IMr(idx)%m, NK, ipiv(idx)%v, work, nwork, info)
        
        DEALLOCATE(work)
        
        ALLOCATE(IM(idx)%m(NK, NK))
        
    end do
    
    ALLOCATE(mesh%IMr, source = IMr)
    
    !  ===========================================================
    !    Now visit each element and determine its inner products  
    !  ===========================================================
    
    ! First, allocation loops structured to make storage of I arrays optimal
    ! Also generally just planning out what needs to be allocated. Visit each idx
    do idx = 1, 7
        if (.not. idxp(idx)) cycle ! Only really need to visit present idx's
        le = 0
        do e = 1, NE
            if (mesh%el(e)%idx /= idx) cycle ! Make sure you cycle for not-present elements
            le = le + 1
        end do
        
        ! If idx is 1, set the spatial dimension to 1
        if (idx == 1) then
            d = 1
        else
            d = 3
        end if
        
        ! Make allocations based on how many elements were present
        ALLOCATE(mesh%IP(idx)%elmap(le), source=0) ! Will need to re-loop to actually get the values. Better than appending I think
        ALLOCATE(mesh%IP(idx)%IC (NKidx(idx), NKidx(idx), d,          le))
        ALLOCATE(mesh%IP(idx)%IFM(NKidx(idx), NKidx(idx), NFidx(idx), le))
        ALLOCATE(mesh%IP(idx)%IIM(NKidx(idx), NKidx(idx), NFidx(idx), le))
        
        le = 0
        do e = 1, NE
            if (mesh%el(e)%idx /= idx) cycle
            le = le + 1
            mesh%IP(idx)%elmap(e) = le
            ! mesh%IP(idx)%elmap(le) = e ! SWITCH THIS BACK AND FIGURE OUT WHY IT SEEMED TO WORK, AND BE FASTER...
        end do
        
    end do
    
    ! MHY LATER - Should parallelize.
    do e = 1, NE
        
        idx = mesh%el(e)%idx
        le  = mesh%IP(idx)%elmap(e)
        ! le  = mesh%IP(idx)%elmap(le) ! SWITCH THIS BACK AND FIGURE OUT WHY IT SEEMED TO WORK, AND BE FASTER...
        
        NK = NKidx(idx)
        NF = NFidx(idx)
        
        !  ----------------------
        !    Determine Jacobian
        !  ----------------------
        
        call Jacobian (idx, mesh%el(e)%node, mesh%rg, Jv, JvinvT, detJv)
        
        mesh%el(e)%detJv = detJv
        
        IM(idx)%m = IMr(idx)%m / detJv
        
        !  -----------------------------------------------
        !    Store the inverse of the mass inner product  
        !  -----------------------------------------------
        
        !  ------------------------------------------------------
        !    Determine mass-eliminated convective inner product  
        !  ------------------------------------------------------
        
        call MEFEAConvection (IM(idx)%m, ICr(idx)%t, JvinvT, detJv, mesh%IP(idx)%IC(:,:,:,le))
        
        do f = 1, NF
            
            fg = mesh%el(e)%face(f)
            
            if (idx == 1) then
                idxf = 0
            else if (idx == 4 .or. (idx == 6 .AND. f <= 2) .or. (idx == 7 .AND. f >= 2)) then
                ! Face is triangular
                idxf = 1
            else
                ! Face is rectangular
                idxf = 2
            end if
            
            !  ------------------------------
            !    Determine surface Jacobian  
            !  ------------------------------
            
            ! call SurfaceJacobian (idxf, mesh%face(fg)%node, mesh%rg, Js, normJs) ! THIS MIGHT BE BROKEN DUE TO mesh%face(fg)%node NOT BEING ORDERED ANYMORE. Not needed unless I use nonlinear elements at some point
            if (idxf == 0) then
                normJs = ONE ! May not be strictly correct, but cancels out with whatever's in IFMr anyway
            else if (idxf == 1) then
                normJs = TWO * mesh%face(fg)%area
            else
                ! MHY LATER - WHEN DOING HEXAHEDRAL CALCS., VERIFY THIS. IT MAY DEPEND ON THE BOUNDS OF THE INTEGRAL I'M TRANSFORMING...
                normJs = mesh%face(fg)%area
            end if
            
            !  -----------------------------------------------------
            !    Determine mass-eliminated face mass inner product  
            !  -----------------------------------------------------
            
            call MEFEAFaceMass (IM(idx)%m, f2k(idx)%v(f)%v, IFMr(idxf)%m, normJs, mesh%IP(idx)%IFM(:,:,f,le))
            
            !  ----------------------------------------------------------
            !    Determine mass-eliminated interface mass inner product  
            !  ----------------------------------------------------------
            
            ! First, determine the properties of the neighboring face. If there is no neighboring
            ! face (i.e., element is on the boundary), take IIM = IFM and cycle
            ep = mesh%el(e)%neighbors(f)
            if (ep == 0) cycle
            
            call MEFEAInterfaceMass &
                (IM(idx)%m, f2k(idx)%v(f)%v, mesh%el(e)%node, mesh%el(ep)%node, IFMr(idxf)%m, normJs, mesh%IP(idx)%IIM(:,:,f,le))
            
        end do
    end do
    
    call DestroyFEAQuadrature ()
    
    ! Deprecated slab code. Merging the slab and tetrahedral case assists in validation, even though it makes the slab case slower (more than worth it).
    !if (idxp(1)) then
    !    
    !    NK = NKidx(1)
    !    NF = NFidx(1)
    !    
    !    DEALLOCATE(IM(1)%m)
    !    
    !    ALLOCATE(IM(1)%m(NK, NK))
    !    
    !    ALLOCATE(mesh%I(NE))
    !    
    !    do e = 1, NE
    !        
    !        detJv = HALF * mesh%el(e)%vol
    !        
    !        IM(1)%m(1,1) =   TWO / detJv
    !        IM(1)%m(1,2) = - ONE / detJv
    !        IM(1)%m(2,1) = IM(1)%m(1,2)
    !        IM(1)%m(2,2) = IM(1)%m(1,1)
    !        
    !        ALLOCATE(mesh%I(e)%IC(NK, NK, 1))
    !        
    !        mesh%I(e)%IC(1,1,1) = - HALF
    !        mesh%I(e)%IC(1,2,1) = - HALF
    !        mesh%I(e)%IC(2,1,1) = HALF
    !        mesh%I(e)%IC(2,2,1) = HALF
    !        
    !        mesh%I(e)%IC(:,:,1) = MATMUL(IM(1)%m, mesh%I(e)%IC(:,:,1))
    !        
    !        ALLOCATE(mesh%I(e)%IFM(NF))
    !        ALLOCATE(mesh%I(e)%IIM(NF))
    !        
    !        do f = 1, NF
    !            ALLOCATE(mesh%I(e)%IFM(f)%m(NK, NK))
    !            ALLOCATE(mesh%I(e)%IIM(f)%m(NK, NK))
    !            
    !            mesh%I(e)%IFM(f)%m(1, 1) = MERGE(ONE,  ZERO, f == 1)
    !            mesh%I(e)%IFM(f)%m(1, 2) = ZERO
    !            mesh%I(e)%IFM(f)%m(2, 1) = ZERO
    !            mesh%I(e)%IFM(f)%m(2, 2) = MERGE(ZERO, ONE,  f == 1)
    !            
    !            if ((e == 1 .AND. f == 1) .or. (e == NE .AND. f == NF)) then
    !                mesh%I(e)%IIM(f)%m = mesh%I(e)%IFM(f)%m
    !            else
    !                mesh%I(e)%IIM(f)%m(1, 1) = ZERO
    !                mesh%I(e)%IIM(f)%m(1, 2) = MERGE(ONE,  ZERO, f == 1)
    !                mesh%I(e)%IIM(f)%m(2, 1) = MERGE(ZERO, ONE,  f == 1)
    !                mesh%I(e)%IIM(f)%m(2, 2) = ZERO
    !            end if
    !            
    !            mesh%I(e)%IFM(f)%m = MATMUL(IM(1)%m, mesh%I(e)%IFM(f)%m)
    !            mesh%I(e)%IIM(f)%m = MATMUL(IM(1)%m, mesh%I(e)%IIM(f)%m)
    !        end do
    !        
    !    end do
    !    
    !end if
    
End Subroutine
End Submodule

        ! Implemented Elements:
        !
        ! Tetrahedron (idx = 4)
        ! Faces:            Reference indices:
        ! 1 2 3             1 2 3
        ! 1 2 4             1 2 3
        ! 1 3 4             1 2 3
        ! 2 3 4             1 2 3
        !
        !                      3           
        !                    ,/|`\         
        !                  ,/  |  `\       
        !                ,/    '.   `\     
        !              ,/       |     `\   
        !            ,/         |       `\ 
        !         1 <-----------'.--------> 2
        !            `\.         |      ,/ 
        !               `\.      |    ,/   
        !                  `\.   '. ,/     
        !                     `\. |/       
        !                        `' 4
        !
        ! Hexahedron (idx = 5)
        ! Faces:            Reference incides:
        ! 1 2 3 4           1 2 3 4
        ! 5 6 7 8           1 2 3 4
        ! 1 2 6 5           1 2 3 4
        ! 3 4 8 7           1 2 3 4
        ! 1 4 8 5           1 2 3 4
        ! 2 3 7 6           1 2 3 4
        !
        !            4----------3    
        !            |\         |\   
        !            | \        | \  
        !            |  \       |  \ 
        !            |   8------+---7
        !            |   |      |   |
        !            1---+------2   |
        !             \  |       \  |
        !              \ |        \ |
        !               \|         \|
        !                5----------6
        !
        ! Prism (idx = 6)
        ! Faces:            Reference indices:
        ! 1 2 3             1 2 3
        ! 4 5 6             1 2 3
        ! 1 2 5 4           1 2 3 4
        ! 2 3 6 5           1 2 3 4
        ! 3 1 4 6           1 2 3 4
        !                        
        !                  4        
        !                ,/|`\      
        !              ,/  |  `\    
        !            ,/    |    `\  
        !           5------+------6 
        !           |      |      | 
        !           |      |      | 
        !           |      |      | 
        !           |      |      | 
        !           |      |      | 
        !           |      1      | 
        !           |    ,/ `\    | 
        !           |  ,/     `\  | 
        !           |,/         `\| 
        !           2-------------3 
        !
        ! 5-node pyramid (idx = 7)
        ! Faces:            Reference indices
        ! 1 2 3 4           1 2 3 4
        ! 1 2 5             1 2 3
        ! 2 3 5             1 2 3
        ! 3 4 5             1 2 3
        ! 4 1 5             1 2 3
        !
        !                      5            
        !                    ,/|\           
        !                  ,/ .'|\          
        !                ,/   | | \         
        !              ,/    .' | `.        
        !            ,/      |  '.  \       
        !          ,/       .'   |   \      
        !        ,/         |    |    \     
        !       1----------.'----4    `.    
        !        `\        |      `\    \   
        !          `\     .'        `\   \ 
        !            `\   |           `\  \ 
        !              `\.'             `\ \ 
        !                 2-----------------3 