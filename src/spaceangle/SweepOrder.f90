Submodule (SpaceAngleInterface) submodSweepOrder
Contains
Module Subroutine SweepOrder (mesh, angular, IsUpstream, IUoffset, sl, isl, estart)
    use constants
    use Parallelism
    use profiler
    use types
    use Geometry
    use AngularSpace
    use SpaceAngleInterface
    Implicit None
    Type (MeshClass),                Intent (In)  :: mesh
    Type (AngularClass),             Intent (In)  :: angular
    
    Logical,            Allocatable, Intent (Out) :: IsUpstream (:)   ! Tells you if the given face in the given element for the given angle is upstream
    Integer,            Allocatable, Intent (Out) :: IUoffset   (:,:) ! Offset array for IsUpstream
    Integer,            Allocatable, Intent (Out) :: sl         (:,:) ! Sweep list
    Integer,            Allocatable, Intent (Out) :: isl        (:,:) ! Inverse sweep list
    Integer,            Allocatable, Intent (Out) :: estart     (:)   ! Index of the last element that can be solved without sweep (name is bad)
    
    Integer                                       :: b, e, es, esp, ep, f, i, eshypo
    Real (KREAL)                                  :: TOL = EPSILON(ONE)
    Logical                                       :: fwd
    Logical,            Allocatable               :: rem        (:)
    Integer                                       :: NENINF
    Integer                                       :: NI
    Integer                                       :: NE
    Integer                                       :: NF
    Integer                                       :: lsss
    Integer                                       :: floc
    Integer                                       :: d
    Integer,            Allocatable               :: nup        (:)
    Integer,            Allocatable               :: lss        (:)
    Real (KREAL),       Allocatable               :: k          (:)
    Real (KREAL),       Allocatable               :: rtemp      (:)
    Real (KREAL),       Allocatable               :: n          (:)
    Character (512)                               :: errmess
    Type (CLogV),       Allocatable               :: sldepend   (:)
    
    d = SIZE(angular%k, DIM=1)
    ALLOCATE(k(d))
    ALLOCATE(rtemp(d))
    ALLOCATE(n(d))
    NE = mesh%NE
    NI = angular%NI
    
    !  ===============================
    !    Determining the sweep order
    !  ===============================
    !  
    !  -----------------------------------------------------------------------------------------------------
    !    ALGORITHM :
    !    1. Choose an angular ordinate k = angular%k(i).
    !    2. Loop through all elements and determine their dependencies, 
    !       i.e., neighbors which are upstream.
    !       When you find an element with no dependencies, add this element to the sweeplist sl(:,i)
    !    3. Loop through elements not present on the sweeplist, and if an element is found with no 
    !       dependency then it is added to the sweeplist sl(:,i).
    !    4. The number of elements added to sl(:,i) is recorded in ss(:,i)
    !    5. Now visit each new addition to the sweeplist, and inform its neighbors that this element is
    !       on the sweeplist, so they no longer depend on it.
    !    6. Repeat from 3 until all elements are on the sweeplist.
    !    7. Repeat from 1 for all angular ordinates
    !    
    !    Step 1 - i loop, 2 - First e loop, 3 - e loop within unbounded es loop (this unbounded loop is 6)
    !    
    !    ORGANIZE THESE NOTES LATER:
    !    IsUpstream(i)%v(e)%v(f) tells you if element e depends on face f (true or false), in general
    !    sldepend(e)%v(f) tells you this information but updates it as the sweeplist is being made
    !  -----------------------------------------------------------------------------------------------------
    
    ! Using these doubly concatenated arrays in this context is dumb and bad and causes crazy storage issues.
    ! So i'll just convert them. see this:
    ! IsUpstream(i)%v(e)%v(f) = IsUpstream(IUoffset(e,i) + f)
    
    ! Pre-allocate rem
    ALLOCATE(rem(NE))
    if (gmemprofiling) call RSSLogger (' --- IN SWEEPORDER --- 1')
    ! Pre-allocate nup
    ALLOCATE(nup(NE))
    if (gmemprofiling) call RSSLogger (' --- IN SWEEPORDER --- 2')
    ! Pre-allocate dependency arrays
    ALLOCATE(sldepend(NE))
    do e = 1, NE
        NF = mesh%el(e)%NF
        ALLOCATE(sldepend(e)%v(NF))
    end do
    if (gmemprofiling) call RSSLogger (' --- IN SWEEPORDER --- 3')
    ALLOCATE(IUoffset(NE + 1, NI))
    NENINF = 0
    do i = 1, NI
        do e = 1, NE
            IUoffset(e,i) = NENINF
            NENINF = NENINF + mesh%el(e)%NF
        end do
        IUoffset(NE + 1, i) = NENINF
    end do
    ALLOCATE(IsUpstream(NENINF))
    if (gmemprofiling) call RSSLogger (' --- IN SWEEPORDER --- 4')
    ! Pre-allocate sweep list arrays
    ALLOCATE(sl (NE, NI), source=0) ! Must be initialized as zero because sl(NE,i) /= 0 is exit condition
    ALLOCATE(isl(NE, NI))
    ALLOCATE(estart(NI))
    if (gmemprofiling) call RSSLogger (' --- IN SWEEPORDER --- 5')
    !$OMP PARALLEL DO DEFAULT (SHARED) PRIVATE (rem, i, k, sldepend, e, ep, esp, f, NF, lss, lsss, nup, es, n, floc)
    do i = 1, NI
        ALLOCATE(lss(0))
        
        k = angular%k(:,i)
        lsss = 0
        ! Initialize depend
        do e = 1, NE
            NF = mesh%el(e)%NF
            do f = 1, NF
                ! If it's a boundary face, then you don't depend on it ever
                if (mesh%face(mesh%el(e)%face(f))%bdy) then
                    IsUpstream(IUoffset(e,i) + f) = .FALSE.
                    sldepend(e)%v(f) = .FALSE.
                    cycle
                end if
                
                ! upstr is true if n dot k is negative.
                ! So, we take .not. n dot k
                n = mesh%el(e)%n(:,f)
                
                IsUpstream(IUoffset(e,i) + f) = DOT_PRODUCT(k, n) < - TOL
                sldepend(e)%v(f) = IsUpstream(IUoffset(e,i) + f)
            end do
            
            nup(e) = COUNT(sldepend(e)%v)
            
        end do
        
        rem = .TRUE.
        es  = 0
        ! Begin sweep set iteration
        do
            lsss = 0
            ! Visit each remaining element and determine its upstream faces
            do e = 1, NE
                if (.not. rem(e)) cycle
                
                if (.not. ANY(sldepend(e)%v)) then
                    ! This is an upstream face. Record it
                    es       = es + 1
                    lsss     = lsss + 1
                    sl(es,i) = e
                    isl(e,i) = es
                    rem(e)   = .FALSE.
                    ! Here I could do an f loop where I set neighbor_with_prev_ss(mesh%el(e)%face(f)) = .TRUE.
                    ! in use with if ((.not. rem(e)) .or. (.not. neighbor_with_prev_ss(e))) cycle
                    ! Note, I'd have to set the entire thing to true at every sweep set iteration (es loop)
                end if
            end do
            
            ! Record the total number of elements in this sweep set
            lss = [lss, lsss]
            
            ! Visit the swept elements' neighbors and set their dependency to false. MUST BE DONE HERE, rather than in e loop, because then you get possibly jagged sweep blocks
            do esp = es - lsss + 1, es
                e  = sl(esp,i)
                NF = mesh%el(e)%NF
                
                do f = 1, NF
                    ! Visit a neighboring element
                    ep = mesh%el(e)%neighbors(f)
                    if (ep == 0) cycle
                    ! This asks, in element ep, what is the face index corresponding to e?
                    floc = findloc(mesh%el(ep)%neighbors, e, DIM=1)
                    
                    ! Set upstream of this face to FALSE
                    sldepend(ep)%v(floc) = .FALSE.
                end do
            end do
            
            ! Exit condition
            if (sl(NE,i) /= 0) exit
            if (lsss == 0) then
                !$OMP CRITICAL
                errmess = 'SweepOrder.f90: Sweep is stagnant. Check your mesh file for errors, ' // &
                          'as this may imply discontinuity in face/element indexing.'
                call stophere (errmess)
                !$OMP END CRITICAL
            end if
        end do
        
        ! Fix 0 as the first index of lss so you can use ss(s) + 1 to ss(s+1) as bounds
        ! THIS IS A HOLDOVER FROM WHEN I WAS DOING ss STUFF.
        ! I still need estart(i), which is technically lss(1) right now, but after the line
        ! below it becomes lss(2). I'm keeping that because I will eventually use ss stuff
        ! and I don't want to confuse myself later
        lss = [0, lss]
        
        ! ! Not doing the ss stuff at the moment, but will return to it when I do 
        ! ! MPI. It really is crucial for that. But I'll have to figure out how to structure the array
        ! ! because I'm getting rid of concatenated arrays.
        ! ! Append lss to ss
        ! ALLOCATE(sweep(i)%ss(SIZE(lss)))
        ! sweep(i)%ss(1) = 0
        ! do b = 2, SIZE(lss)
        !     sweep(i)%ss(b) = sweep(i)%ss(b - 1) + lss(b)
        ! end do
        estart(i) = lss(2) ! The starting index for the forward iteration is the second index of lss. Because lss(b) gives the b'th sweep block, and the first sweep block is all boundary faces.
        DEALLOCATE(lss)
        
    end do
    !$OMP END PARALLEL DO
    if (gmemprofiling) call RSSLogger (' --- IN SWEEPORDER --- 6')
    
End Subroutine
End Submodule