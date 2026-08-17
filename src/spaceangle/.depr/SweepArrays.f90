Submodule (SpaceAngleInterface) submodSweepArrays
Contains
Module Subroutine SweepArrays (mesh, angular, upfaces, upfaces_os, sl, isl, estart, nup, IAF)
    
    ! This deprecated file contains code that was relevant when I was storing the off-diagonal entries
    ! of the entire transport matrix (which is a huge albeit sparse matrix, it is (NENK * NI)^2)
    ! to apply in the sweep. This is a very very bad construct, because it demands an insane amount of
    ! memory, while not really speeding up much.
    ! 
    ! As far as speed is concerned what do you actually save doing this?
    ! 1. You don't need to multiply k dot n by the NK x NK' matrix during the sweep (negligible)
    ! 2. There are several zero entries in the interface mass matrix that you would be able to avoid 
    !    multiplying by
    ! Both of these are absolutely negligible and well worth saving the significant amount of space.
    ! Furthermore, it is not a guarantee that a SBR multiplication of a small matrix with ~1/2 of its
    ! entries being zero is cheaper than a Fortran MATMUL of the same matrix. Also, this routine itself
    ! does take quite some time, which even for problems requiring many iterates may not be worth it.
    ! In all, I had not thoroughly profiled the sweep before and after getting rid of this, but
    ! I do see calculations taking about the same amount of time. So this routine will likely never
    ! make a comeback, even in the (.not. storagelimited) case.
    ! 
    ! Later on I may also define and store permutation matrices for the interface mass matrix so that
    ! it doesn't need to be stored, only the face mass matrix will need to be stored. However this will
    ! require a revised matrix multiplication construct in TransportAlgebra, since MATMUL will not work
    ! anymore.
    
    use constants
    use Parallelism
    use types
    use Geometry
    use AngularSpace
    use SpaceAngleInterface
    Implicit None
    Type (MeshClass),    Intent (In)  :: mesh
    Type (AngularClass), Intent (In)  :: angular
    Logical,             Intent (In)  :: upfaces    (:)
    Integer,             Intent (In)  :: upfaces_os (:,:)
    Integer,             Intent (In)  :: sl         (:,:)
    Integer,             Intent (In)  :: isl        (:,:)
    Integer,             Intent (In)  :: estart     (:)
    Integer,             Intent (In)  :: nup        (:,:)
    
    Type (UpstreamSys), Intent (Out) :: IAF
    
    Integer                           :: e, es, esp, ep, f, fp, fup, fupp, i, j, k, kp, row
    Real (KREAL),        Parameter    :: TOL = EPSILON(ONE)
    Integer                           :: d
    Integer                           :: NE
    Integer                           :: NI
    Integer                           :: NK
    Integer                           :: NKf
    Integer                           :: NF
    Integer                           :: NK1
    Integer                           :: NK2
    Integer                           :: ipcol
    Integer                           :: iprow
    Integer                           :: lnup
    Integer                           :: lclupfaces (5) ! Because there's up to 6 faces. You will never have all faces be dependent.
    Integer                           :: upels      (5)
    Integer,             Allocatable  :: NZoffset   (:)
    Integer,             Allocatable  :: trowptr    (:)
    Integer (KINT)                    :: NZ
    Integer (KINT),      Allocatable  :: NZi        (:)
    Real (KREAL)                      :: kdotn      (5)
    Real (KREAL),        Allocatable  :: kv         (:)
    Real (KREAL),        Allocatable  :: n          (:)
    Real (KREAL)                      :: tmp8       (8) ! Because there's up to 8 nodes in a general mesh element
    !Type (CIntV)                      :: NZf        (7)    ! For a given index and face, gives the number of nonzero entries in the IIM matrix ! (idx)%v(f) ! NOTE - could also use face idx instead, to make a 1D array...
    
    ! HOW TO DETERMINE NZ WITHOUT VISITING ELEMENTS:
    ! In MATMUL(A, B), assuming A is fully dense (which it is for tetrahedra, although maybe for some shapes it's not),
    ! the number of nonzero entries will be determined ONLY by the nonzero columns in B. If B has m nonzero columns,
    ! and each column has n rows, then the number of nonzero entries is m * n.
    ! So for every interface type (i.e., given the idx's of elements on both sides of the face), we can determine
    ! the NZ.
    !
    ! Rows are NK1. Columns are NK2.
    ! The number of nonzero columns is NKf. Thus, m = NKf
    ! NK1 is the number of rows. Thus:
    ! NZf = NK1 * NKf
    ! Great, we only need to visit the element.
    ! BIG QUESTION : For meshes with various types of elements, will NZ be the same for every i??? I am actually not sure.
    ! The number of upstream interior faces is ALWAYS NF - NBF... But since for different angles, NK1 changes, I think it might actually be i
    ! dependent after all...
    
    !ALLOCATE(NZf(1)%v(2))
    !ALLOCATE(NZf(4)%v(4))
    !ALLOCATE(NZf(5)%v(6))
    !ALLOCATE(NZf(6)%v(5))
    !ALLOCATE(NZf(7)%v(5))
    
    !NZf(1)%v(1) = 1
    !NZf(1)%v(1) = 1
    
    ! MHY LATER - Do this with MPI/HPC compat. - I found that for a mesh of ~250,000 elements and for 512 angles,
    !       the size of IAF exceeds the max for single precision integers
    !       (2,147,483,647)
    ! Need to generalize to higher precision integers... should not be complicated at all...
    
    d  = SIZE(angular%k, DIM=1)
    NE = mesh%NE
    NI = angular%NI
    
    ALLOCATE(kv(d))
    ALLOCATE(n(d))
    
    ! Create offset array.
    ! IAF%offset(es + (i - 1) * NE) + 1 gives the value of the first index of a given (es,i) block in a vector using collapsed notation
    ALLOCATE(IAF%offset(NE * NI + 1))
    j = 0
    do i = 1, NI
        do es = 1, NE
            e  = sl(es,i)
            NK = mesh%el(e)%NK
            
            IAF%offset(es + (i - 1) * NE) = j
            j = j + NK
            
        end do
    end do
    
    IAF%offset(NE * NI + 1) = j
    
    ! Now you must visit every nonzero value in Fdown, and determine:
    ! IAF%vals   ---> A list of every value in sequential order, running along the rows for a given i.        SIZE : NZ    (Number of nonzero entries)
    ! IAF%cols   ---> A list of the corresponding cols                                                        SIZE : NZ
    ! IAF%rowptr ---> Defined such that rowptr(col) + 1 : rowptr(col + 1) gives the set of j values in a row. SIZE : N + 1 (Dimension of the system + 1)
    
    ! First, determine the number of nonzero entries per i.
    ! Each global interior face is shared by two elements. For ANY given i, one of these faces will be downstream.
    ! Thus, we need only to visit every face, determine the number of nodes it has, determine which of the shared
    ! elements is downstream, and from that we can determine the number of nonzero elements it contributes
    ALLOCATE(NZi(NI))
    if (mesh%homog) then
        ! All angles have the same number of dependent faces
        NZi = 0
        do es = estart(1) + 1, NE
            e    = sl(es,1)
            NF   = mesh%el(e)%NF
            NK1  = mesh%el(e)%NK
            lnup = nup(es,1)
            
            ! WHAT I SHOULD DO IS NOT keep upfaces the mask, keep upfaces the actual list, and do indirect addressing on THAT.
            ! SHOULD I REPLACE upfaces WITH THE ACTUAL FACE INDICES???
            lclupfaces(1:lnup) = PACK([(f,f=1,NF)], mask = upfaces(upfaces_os(e,1) : upfaces_os(e,1) + NF - 1))
            ! The number of dependent faces in total is NF - NBF
            ! The number of NZ elements should be easy to determine... visit each face and depending on its number of nodes, determine?
            ! As for upfaces, make it work with face els... in sweeporder, visit a global face, and then determine isdownstream == k dot n < ZERO
            ! using n in element 1 (of 2). Then, here, we visit the faces of element e such that fg(e,f) is downstream according to isdownstream(fg). 
            ! IS that slower? 
            ! SHOULD I DO THIS FOR ALL NORMAL VECTORS? MAKE THEM A PROPERTY OF fg, AND THEN CARRY SIGNS(e) WHICH IS +1 IF THAT e USES n DIRECTLY,
            ! -1 OTHERWISE???
            ! I don't know the best way to do it... for now just move forward
            
            do fup = 1, lnup
                f   = lclupfaces(fup)
                NKf = mesh%face(mesh%el(e)%face(f))%NK
                ep  = mesh%el(e)%neighbors(f)
                
                NZi = NZi + NK1 * NKf
            end do
        end do
    else
        ! Loop over i for this one. Use the NZf rules? Remove lclupfaces and use something else? Treat this when
        ! it comes time to look at meshes with differently shaped elements... not soon
        print *, 'STOPPING - wdwpdksod' ! This means WIP. Currently only tetrahedral meshes are allowed.
        stop
    end if
    
    ALLOCATE(NZoffset(NI))
    NZoffset(1) = 0
    do i = 2, NI
        NZoffset(i) = NZoffset(i - 1) + NZi(i - 1)
    end do
    
    NZ = NZoffset(NI) + NZi(NI)
    
    ! Initialize IAF%cols and IAF%vals
    ALLOCATE(IAF%vals(NZ))
    ALLOCATE(IAF%cols(NZ))
    
    IAF%cols = 0
    IAF%vals = ZERO
    
    ! Initialize trowptr, which will be used to create the number of nonzero entries before a given row (IAF%rowptr)
    ALLOCATE(trowptr(mesh%NENK * NI), source = 0)
    
    ! Begin loop to construct IAF
    do i = 1, NI ! THE ONLY PARALLELIZABLE INDEX
        kv(1:d) = angular%k(1:d,i)
        
        j = 0
        do es = estart(i) + 1, NE
            e    = sl(es,i)
            NK1  = mesh%el(e)%NK
            lnup = nup(es,i)
            
            ! Starting row MINUS 1 of a given (es, i) block
            iprow = IAF%offset(es + (i - 1) * NE)
            
            ! Determine the list of upstream faces
            NF = mesh%el(e)%NF
            lclupfaces(1:lnup) = PACK([(f,f=1,NF)], mask = upfaces(upfaces_os(e,i) : upfaces_os(e,i) + NF - 1))
            
            ! Now determine their corresponding element indices in the sweep
            upels = NE + 1
            upels(1:lnup) = isl(mesh%el(e)%neighbors(lclupfaces(1:lnup)), i)
            
            ! Do a quick insertion sort so that the first values in upels are the smallest esp values
            ! The purpose of this is to match up lclupfaces so you can visit the correct faces.
            do fup = 2, lnup
                esp  = upels(fup)
                fp   = lclupfaces(fup)
                fupp = fup - 1
                do
                    if (fupp < 1) exit
                    if (upels(fupp) > esp) then
                        upels(fupp + 1)   = upels(fupp)
                        lclupfaces(fupp + 1) = lclupfaces(fupp)
                        fupp = fupp - 1
                    else
                        exit
                    end if
                end do
                upels(fupp + 1)      = esp
                lclupfaces(fupp + 1) = fp
            end do
            
            ! Now populate the off-diagonal entries of the transport matrix
            ! Since we're going row-by-row, but the faces indicate columns in the element block structure,
            ! you visit every face and write down the dot product of k and n. Later on, faces will be
            ! iterated through, in each individual row, so this avoids calculating it multiple times.
            ! We also update trowptr here
            do fup = 1, lnup
                f      = lclupfaces(fup)
                ep     = mesh%el(e)%neighbors(f)
                NKf    = mesh%face(mesh%el(e)%face(f))%NK
                NK2    = mesh%el(ep)%NK
                n(1:d) = mesh%el(e)%n(1:d,f)
                
                ! Now determine the column and row where the NK1 x NK2 block begins.
                esp = upels(fup)
                
                ! Starting column MINUS 1 of a given (esp, i) block
                ipcol = IAF%offset(esp + (i - 1) * NE)
                
                ! Update trowptr with the number of nonzero entries in this block's row, which is NKf.
                trowptr(iprow + 1 : iprow + NK1) = trowptr(iprow + 1 : iprow + NK1) + NKf
                
                ! Determine and temporarily store the needed object
                kdotn(fup) = DOT_PRODUCT(n(1:d), kv(1:d))
            end do
            
            do k = 1, NK1
                do fup = 1, lnup
                    ! MHY LATER - Since this subroutine is only called once and since it's the only place using IIM, should I just restructure IFM here on-the-fly?
                    ! Use/store permutation matrices for this purpose?
                    f     = lclupfaces(fup)
                    ep    = mesh%el(e)%neighbors(f)
                    esp   = upels(fup)
                    ipcol = IAF%offset(esp + (i - 1) * NE)
                    NK2   = mesh%el(ep)%NK
                    
                    tmp8(1:NK2) = kdotn(fup) * mesh%I(e)%IIM(f)%m(1:NK2, k)
                    do kp = 1, NK2
                        if (ABS(tmp8(kp)) <= TOL) cycle ! Restructure IIM to avoid this? In general IIM can be foregone too
                        
                        j = j + 1
                        
                        IAF%cols(j + NZoffset(i)) = ipcol + kp
                        IAF%vals(j + NZoffset(i)) = tmp8(kp)
                    end do
                end do
            end do
        end do
    end do
    
    ! Finally, create rowptr
    ALLOCATE(IAF%rowptr(1:SIZE(trowptr) + 1))
    ! BEFORE the first row, there are zero nonzero elements
    IAF%rowptr(1) = 0
    do row = 2, SIZE(trowptr) + 1
        ! At the nth row, you sum trowptr from 1 to row - 1. This can be done with:
        ! IAF%rowptr(row) = SUM(trowptr(1:row - 1))
        ! OR, more efficiently, iterate:
        IAF%rowptr(row) = IAF%rowptr(row - 1) + trowptr(row - 1)
    end do
    
End Subroutine
End Submodule