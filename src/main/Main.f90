Program Main
    ! Core
    use constants
    use Parallelism
    use profiler
    use IO
    ! NittanyAPI
    use NittanyAPI
    ! LionboltAPI
    use LionboltAPI
    ! Peripheral routines
    use globals
    use Utilities
    use InputModule
    Implicit None
    
    !                               ####################
    !                                     Lionbolt      
    !                                    Written by     
    !                                 Muhsin H. Younis  
    !                               ####################
    
    !  ==============================================================================
    !    License is found at $LIONBOLT/LICENSE/gpl-3.0.txt                           
    !    where $LIONBOLT is the main Lionbolt directory                              
    !  ==============================================================================
    
    !  ==============================================================================
    !    Lionbolt is part of the Open Radiation Physics Suite (OpenRPS). For         
    !    thorough documentation, including quickstart, installation, and other       
    !    important info., see the OpenRPS website, at:                               
    !        https://mhyounis.github.io/OpenRPS/                                     
    !  ==============================================================================
    
    !  ==============================================================================
    !    This file uses the Lionbolt API to carry out the standard Lionbolt program  
    !    solution procedure, based on a user-provided input file. Users have the     
    !    option to alternatively:                                                    
    !        1) Create their own solution procedure ('driver'), not relying on       
    !           input.                                                               
    !        2) Directly link Lionbolt to their own program and rely on the API      
    !           there.                                                               
    !  ==============================================================================
    
! Objects
Logical                                  :: deltadown
Type (MeshClass),    Target              :: mesh
Type (AngularClass), Target              :: angular
Type (XSType),       Target, Allocatable :: XS     (:)
Type (XSType),               Allocatable :: CXS    (:,:) ! Must come up with better way to store and use coupling XS
Type (XSLibrary),            Allocatable :: XSLib  (:,:)
Type (EnergyGrid),           Allocatable :: energy (:)
Type (FieldGeometry)                     :: fldgeo
Type (ExternalBeam),         Allocatable :: beams  (:)
Type (ExternalBeam)                      :: beam
Type (InputFile)                         :: input
Type (ParticleInput)                     :: particle
Type (BoltzmannOp)                       :: OpL
Type (ScatteringOp)                      :: OpKc         ! Coupling operator (an instance of the scattering operator that couples energy groups as well as particles to eachother)

! Vectors
Type (SpaceAngleVector)   :: s      ! Flattened source array (in) / solution array (out)
Type (SpaceAngleVector)   :: w1     ! Flattened workspace vector
Type (SpaceAngleVector)   :: w2     ! Flattened workspace vector ! May be able to forego this if I refactor PIEISource
Type (SpaceAngleVector)   :: unc    ! Flattened 'uncollided' fluence (includes contributions from delta-down scattering if applicable)
Real (KREAL), Allocatable :: fl (:) ! Fluence (if requested)

! Iteration variables and other relevant objects
Integer                   :: b, e, g, gp, p, pb, pin
Logical                   :: unciszero
Integer                   :: NP
Integer                   :: nfmt
Integer                   :: mfmt
Integer                   :: Gfmt
Integer,      Allocatable :: Gs     (:)
Real (KREAL)              :: wnorm
Real (KREAL)              :: dEtot
Real (KREAL), Allocatable :: DEP    (:)
Character (STDMAXLEN)     :: fmt
Type (CIntV), Allocatable :: p2bind (:) ! Maps particle index to beam index

!  ==============================================================================================
!    Open the output file, HDF5, determine base path and write header, then start tracking time  
!  ==============================================================================================

call OpenOutput ()
call Header     ()
call StartTime  ()
call OpenHDF5   (ScratchName // 'results.tmp.h5')

if (gtiming)       call InitializeTimer     ()
if (gmemprofiling) call InitializeRSSLogger ()

!  ===================
!    Read input file  
!  ===================

write (iuout,*) 'READING INPUT FILE ---'
input = InputFile ()
call input%options%AssignGlobals ()
NP = SIZE(input%particles%v)
write (iuout,*) 'INPUT FILE SUCCESSFULLY READ.'
write (iuout,*) ''

if (gmemprofiling) call RSSLogger ('After reading input file')

!  ===================================================
!    From the problem type, determine how to proceed  
!  ===================================================

if (input%problem%General) then
    
    write (iuout,*) '#############################'
    write (iuout,*) ' PROBLEM : GENERAL TRANSPORT '
    write (iuout,*) '#############################'
    write (iuout,*) ''
    
    call WriteProblemTypeH5 ('general') ! Should I make the HDF5 file an object and then have this (and other routines) be type-bound?
    
else if (input%problem%Slab) then
    
    write (iuout,*) '##########################'
    write (iuout,*) ' PROBLEM : SLAB TRANSPORT '
    write (iuout,*) '##########################'
    write (iuout,*) ''
    
    call WriteProblemTypeH5 ('slab')
    
end if

write (iuout,*) '=========================='
write (iuout,*) '  PRE-PROCESSING : BEGIN  '
write (iuout,*) '=========================='
write (iuout,*) ''

!  ==============================
!    Plan for beam construction  
!  ==============================

write (iuout,*) 'PLANNING BEAM CONSTRUCTION ---'
ALLOCATE(p2bind(NP))
do p = 1, NP
    ALLOCATE(p2bind(p)%v(0))
end do
do b = 1, SIZE(input%beams%v)
    p2bind(input%beams%v(b)%particle)%v = [p2bind(input%beams%v(b)%particle)%v, b]
end do
write (iuout,*) 'DONE.'
write (iuout,*) ''

!  ==========================================
!    Discretize all particles' energy grids  
!  ==========================================

write (iuout,*) 'DISCRETIZING ENERGY ---'
ALLOCATE(energy(NP))
ALLOCATE(Gs(NP))
do p = 1, NP
    energy(p) = input%particles%v(p)%BuildEnergy (input%beams%v(p2bind(p)%v))
    Gs(p) = energy(p)%G
end do
write (iuout,*) 'DONE.'
write (iuout,*) ''

!  ================================
!    Generate the physics library  
!  ================================

write (iuout,*) 'GENERATING A PHYSICS LIBRARY ---'
XSLib = input%GenerateXSLibrary (energy)
write (iuout,*) 'DONE.'
write (iuout,*) ''

if (gmemprofiling) call RSSLogger ('After creating XS library')

write (iuout,*) 'DISCRETIZING GEOMETRY ---'
mesh = input%mesh%Build ()
call WriteMeshH5 (mesh)
write (iuout,*) 'DONE.'
write (iuout,*) ''

if (gmemprofiling) call RSSLogger ('After discretizing mesh')

write (iuout,*) '========================'
write (iuout,*) '  PRE-PROCESSING : END  '
write (iuout,*) '========================'
write (iuout,*) ''

write (iuout,*) '=============================='
write (iuout,*) '  PARTICLE ITERATION : BEGIN  '
write (iuout,*) '=============================='

! Admittedly, this part of the code has gotten very messy.
! 
! The initial construction was a series of subroutine calls
! that take minimal inputs, like mesh, angular, etc., but the
! severing of Lionbolt the program from Lionbolt the library has forced
! a lot of that subroutine code to hide in the input object, and this
! has notably led to a lot of pre-processing for input attributes that
! are arrays (like beams). Also, some stray objects like the particle
! object is not an API feature yet is important here.

do p = 1, NP
    
    particle = input%particles%v(p)
    
    unciszero = .TRUE.
    
    !  ===========================================
    !    Open a new HDF5 group for this particle  
    !  ===========================================
    
    call OpenParticleH5 (particle%name)
    
    write (iuout,*) REPEAT('-', 29 + LEN(particle%name) + 3)
    write (iuout,'(1X,A,I0,A)') '   TRANSPORTING PARTICLE ', p, ' ('// particle%name //')'
    write (iuout,*) REPEAT('-', 29 + LEN(particle%name) + 3)
    write (iuout,'(4X,A)') ''
    
    ! Briefly write the energy
    call energy(p)%Print (iuout)
    call AppendEnergyGridH5 (energy(p)%E)
    
    !  =================================================
    !    Discretize the particle's angular coordinates  
    !  =================================================
    
    write (iuout,'(4X,A)') 'DISCRETIZING ANGLE ---'
    angular = particle%BuildAngular ()
    call angular%Print (iuout)
    call AppendOrdinatesH5 (angular%w, angular%k)
    write (iuout,'(4X,A)') 'DONE.'
    write (iuout,'(4X,A)') ''
    
    if (gmemprofiling) call RSSLogger ('    After creating angular for particle : ' // particle%name)
    
    !  ==================================
    !    Initialize the solution vector  
    !  ==================================
    
    write (iuout,'(4X,A)') 'INITIALIZING SOLUTION VECTORS ---'
    s   = SpaceAngleVector (mesh, angular)
    w1  = SpaceAngleVector (mesh, angular)
    w2  = SpaceAngleVector (mesh, angular)
    unc = SpaceAngleVector (mesh, angular)
    if (input%postproc%fluence) then
        if (ALLOCATED(fl)) DEALLOCATE(fl)
        ALLOCATE(fl(mesh%NENK))
    end if
    write (iuout,'(4X,A)') 'DONE.'
    write (iuout,'(4X,A)') ''
    
    if (gmemprofiling) call RSSLogger ('    After initializing solution vectors')
    
    ! Prepare the coupling XS array. It is shaped to
    ! the largest used particle energy, even though many entries may
    ! not be used. This is fine considering most attributes of XS are
    ! themselves allocatable.
    
    if (ALLOCATED(CXS)) DEALLOCATE(CXS)
    ALLOCATE(CXS(MAXVAL(Gs), p - 1))
    
    !  ================================================================
    !    If this particle is sourced by an external beam, generate it  
    !  ================================================================
    
    if (.not. ALLOCATED(beams)) ALLOCATE(beams(SIZE(input%beams%v)))
    
    ! Now normalize weights if this is particle 1
    if (p == 1) then
        wnorm = ZERO
        do pb = 1, SIZE(p2bind(1)%v)
            b = p2bind(1)%v(pb)
            
            wnorm = wnorm + input%beams%v(b)%weight
        end do
        do pb = 1, SIZE(p2bind(1)%v)
            b = p2bind(1)%v(pb)
            
            input%beams%v(b)%weight = input%beams%v(b)%weight / wnorm
        end do
    end if
    
    ! Now finally construct the beam
    do pb = 1, SIZE(p2bind(p)%v)
        b = p2bind(p)%v(pb)
        
        write (iuout,'(4X,A)') 'CREATING AN EXTERNAL BEAM ---'
        fldgeo = input%beams%v(b)%Build ()
        beam   = ExternalBeam (mesh, angular, fldgeo, particle%name)
        
        beams(b)  = beam
        
        input%beams%v(b)%used = .TRUE.
        
        unciszero = .FALSE.
        
        write (iuout,'(4X,A)') 'DONE.'
        write (iuout,'(4X,A)') ''
    end do
    
    if (gmemprofiling) call RSSLogger ('    After constructing beams')
    
    !  ==================================================
    !    Build a Boltzmann operator's space-angle parts  
    !  ==================================================
    
    write (iuout,'(4X,A)') 'CREATING A BOLTZMANN OPERATOR ---'
    call OpL%Build (mesh, angular)
    write (iuout,'(4X,A)') 'DONE.'
    write (iuout,'(4X,A)') ''
    
    if (gmemprofiling) call RSSLogger ('    After creating Boltzmann operator')
    
    !  ==============================
    !    Create a coupling operator  
    !  ==============================
    
    write (iuout,'(4X,A)') 'CREATING A COUPLING OPERATOR ---'
    call OpKc%Build (mesh, angular)
    write (iuout,'(4X,A)') 'DONE.'
    write (iuout,'(4X,A)') ''
    
    if (gmemprofiling) call RSSLogger ('    After creating coupling operator')
    
    do g = 1, energy(p)%G
        
        nfmt = MAX(INT(LOG10(energy(p)%E(g + 1))) + 1, 1) ! MAX(.,1) is to retain the leading zero digit in the case that nfmt <= 1
        mfmt = MAX(INT(LOG10(energy(p)%E(g)))     + 1, 1)
        Gfmt = INT(LOG10(REAL(energy(p)%G))) + 1
        
        write (fmt,'( "(7X,A,I" , I0 , ",A,F" , I0 , ".5,A,F" , I0 , ".5,A)" )') Gfmt, nfmt + 6, mfmt + 6
        
        write (iuout,'(4X,A)') REPEAT('-', 43 + Gfmt + nfmt + mfmt)
        write (iuout, fmt)    'ENERGY GROUP ', g, ' ( ' , energy(p)%E(g+1), ' - ', energy(p)%E(g), ' MeV )   '
        write (iuout,'(4X,A)') REPEAT('-', 43 + Gfmt + nfmt + mfmt)
        
        !  ===========================================
        !    Read the physics library to generate XS  
        !  ===========================================
        
        write (iuout,'(7X,A)') 'Reading cross sections from library'
        if (ALLOCATED(XS)) DEALLOCATE(XS)
        ALLOCATE(XS(g))
        do gp = 1, g
            call ReadXSLibrary (p, p, gp, g, XSLib, XS(gp))
        end do
        
        !  ===================================================
        !    Read the physics library to generate coupled XS  
        !  ===================================================
        
        write (iuout,'(7X,A)') 'Reading coupling cross sections from library'
        do pin = 1, p - 1
            do gp = 1, Gs(pin)
                call ReadXSLibrary (pin, p, gp, g, XSLib, CXS(gp, pin))
            end do
        end do
        
        !  ================================================
        !    Update the Boltzmann operator for this group  
        !  ================================================
        
        write (iuout,'(7X,A)') 'Updating Boltzmann Operator'
        call OpL%Build (angular, XS(g))
        
        if (gmemprofiling) call RSSLogger ('        After updating Boltzmann operator')
        
        !  =======================================
        !    Construct the source for this group  
        !  =======================================
        
        write (iuout,'(7X,A)') 'Constructing Source'
        deltadown = SIZE(XS(g)%s, dim=1) == angular%L + 2
        unciszero = unciszero .and. (.not. deltadown)
        ! Now populate each beam's particles
        do pb = 1, SIZE(p2bind(p)%v)
            b = p2bind(p)%v(pb)
            
            if (streq(beams(b)%fldgeo%spectrum, 'constant')) then
                ! First determine if the user specified groups or E0 and dE
                if (input%beams%v(b)%E0 <= ZERO) then
                    ! User specified groups
                    if (.not. ANY(input%beams%v(b)%mcgroups == g)) then
                        ! This group has no particles populated
                        beams(b)%n = ZERO
                        cycle
                    else
                        ! The constant value is given by the fraction of particles in this group times the beam's weight
                        ! The fraction is determined by putting the same number of particles per unit energy. That is,
                        ! if you plotted the spectrum, you'd find the same height over each group, so the n values are
                        ! generally different.
                        
                        ! Here you determine the normalization by just adding up the energy widths
                        ! of each group. That accomplishes what's stated above.
                        dEtot = ZERO
                        do gp = 1, energy(p)%G
                            if (.not. ANY(input%beams%v(b)%mcgroups == gp)) cycle
                            dEtot = dEtot + energy(p)%E(gp) - energy(p)%E(gp + 1)
                        end do
                        
                        beams(b)%n = input%beams%v(b)%weight &
                                   * (energy(p)%E(g) - energy(p)%E(g + 1)) &
                                   / dEtot
                        write (iuout,'(10X,A,I0,A,F9.6)') 'Population of beam ', b, ' : ', beams(b)%n
                    end if
                else
                    if (g > input%beams%v(b)%mcgroups(1)) then
                        ! No beam
                        beams(b)%n = ZERO
                    else
                        ! User specified width
                        beams(b)%n = input%beams%v(b)%weight &
                                   * (energy(p)%E(g) - energy(p)%E(g + 1)) &
                                   / input%beams%v(b)%dE
                        write (iuout,'(10X,A,I0,A,F9.6)') 'Population of beam ', b, ' : ', beams(b)%n
                    end if
                end if
            else
                call beams(b)%Populate (energy(p)%E(energy(p)%G + 1), energy(p)%E(1), energy(p)%E(g + 1), energy(p)%E(g))
                write (iuout,'(10X,A,I0,A,F9.6)') 'Population of beam ', b, ' : ', beams(b)%n
            end if
        end do
        
        call PIEISource (p, g, deltadown, &
                         p2bind, beams, mesh, angular, XS, CXS, OpL, OpKc, w1, w2, unc, s)
        
        if (gmemprofiling) call RSSLogger ('        After creating the source')
        
        !  =========
        !    SOLVE  
        !  =========
        
        call OpL%MatInv (s)
        
        if (gmemprofiling) call RSSLogger ('        After solving')
        
        !  =====================================
        !    Recording : Total Angular Fluence  
        !  =====================================
        
        ! Add the uncollided angular fluence if present (s is the solution at this stage, so s <-- s + unc is the total angular fluence)
        if (.not. unciszero) s%v = s%v + unc%v
        
        ! Write angular fluence to tape (for coupling to later groups/particles)
        call TapeWriter (p, g, s%v)
        
        ! Write the angular fluence to the output HDF5 file
        if (.not. gstorlim) call AppendAngularFluence (s%v, mesh%NENK, .FALSE., g)
        
        ! Write the fluence to the output HDF5 file if requested
        if (input%postproc%fluence) then
            fl = s%Fluence (angular%w)
            call AppendSpatialArray (GROUP_FL, fl, g)
        end if
        
        !  ==========================================
        !    Recording : Uncollided Angular Fluence  
        !  ==========================================
        
        if (.not. unciszero .and. input%postproc%outputunc) then
            
            ! Write the primary angular fluence to the output HDF5 file
            if (.not. gstorlim) call AppendAngularFluence (unc%v, mesh%NENK, .TRUE., g)
            
            ! Write the primary fluence to the output HDF5 file if requested
            if (input%postproc%fluence) then
                fl = unc%Fluence (angular%w)
                call AppendSpatialArray (GROUP_FL_UNC, fl, g)
            end if
            
        end if
        
        if (gmemprofiling) call RSSLogger ('        After writing everything')
        
    end do
    
    !  ====================================
    !    Close this particle's HDF5 group  
    !  ====================================
    
    call CloseParticleH5 (particle%name)
    
    !  ==================================================
    !    Destroy particle-dependent phase space objects  
    !  ==================================================
    
    ! MHY LATER - issue with this is that postproc needs angular. Easy fix, just make postproc regenerate it. So it'll need particles.
    ! For now though not destroying angular is really not an issue at all.
    ! This is because I redefine the variable as a whole using input build routines
    ! call angular%Destroy ()
    
end do

write (iuout,*) '============================'
write (iuout,*) '  PARTICLE ITERATION : END  '
write (iuout,*) '============================'
write (iuout,*) ''

if (ANY([input%postproc%energy, input%postproc%dose, input%postproc%charge])) then
    write (iuout,*) '==========================='
    write (iuout,*) '  POST-PROCESSING : BEGIN  '
    write (iuout,*) '==========================='
    write (iuout,*) ''
    
    ! Free up some space
    call w1%Destroy  ()
    call w2%Destroy  ()
    call unc%Destroy ()
    
    call CalculateDeposition (input%postproc, NP, Gs, mesh, angular%w, s)
    
    write (iuout,*) '========================='
    write (iuout,*) '  POST-PROCESSING : END  '
    write (iuout,*) '========================='
    write (iuout,*) ''
    
end if

!  ===============
!    TERMINATION  
!  ===============

if (gmemprofiling) call PublishRSSLogger ()

call CloseHDF5 ()
call CleanUp   ()

if (gtiming) call Timestamp ('Publishing HDF5')
call PublishHDF5 ()
if (gtiming) call Timestamp ('Publishing HDF5')
if (gtiming) call PublishTimer     ()

call StopTime    ()
call CloseOutput ()

    Contains
    
    ! A subroutine for forming the source within the PI/EI loop. Just to keep the main program code clean.
    ! Very very ugly.
    ! This is NOT meant to be highly generalizable nor user/developer-friendly. 
    ! It IS meant to be somewhat optimized however. That will be addressed more fully in a later update.
    ! Anyway the source construction is the slowest part of the program right now.
    Subroutine PIEISource (p, g, deltadown, p2bind, beams, mesh, angular, XS, CXS, OpL, OpKc, w1, w2, unc, s)
        Implicit None
        Integer,                  Intent (In)    :: p
        Integer,                  Intent (In)    :: g
        Logical,                  Intent (In)    :: deltadown
        Type (CIntV),             Intent (In)    :: p2bind (:)
        Type (ExternalBeam),      Intent (InOut) :: beams  (:)
        Type (MeshClass),         Intent (In)    :: mesh
        Type (AngularClass),      Intent (In)    :: angular
        Type (XSType),            Intent (In)    :: XS     (:)
        Type (XSType),            Intent (In)    :: CXS    (:,:)
        
        Type (BoltzmannOp),       Intent (InOut) :: OpL  ! MHY LATER - REMOVE INOUT IF I CAN REMOVE IT FROM THE MATVECS
        Type (ScatteringOp),      Intent (InOut) :: OpKc ! Keep InOut because of updating via Coupling
        Type (SpaceAngleVector),  Intent (InOut) :: w1
        Type (SpaceAngleVector),  Intent (InOut) :: w2
        Type (SpaceAngleVector),  Intent (InOut) :: unc
        Type (SpaceAngleVector),  Intent (InOut) :: s
        
        Integer                                  :: gp, pb, pin
        Logical                                  :: unciszero = .TRUE.
        Real (KREAL)                             :: attn (mesh%Nm)
        Real (KREAL)                             :: XSdd (mesh%Nm)
        
        !  ===================
        !    Initializations  
        !  ===================
        
        attn = XS(g)%t
        
        w1%v  = ZERO
        w2%v  = ZERO
        unc%v = ZERO
        s%v   = ZERO
        
        !  ======================================
        !    Form the generic PI and EI sources  
        !  ======================================
        
        ! EI
        do gp = 1, g - 1
            if (gtiming) call Timestamp ('EI')
            if (ALL(XS(gp)%s == ZERO)) go to 1
            
            ! Prepare the coupling operator with the transfer moments for (gp -> g)
            call OpKc%Build (angular, XS(gp))
            
            ! Read in the solution at energy group gp
            call TapeReader (p, gp, w1%v)
            
            ! If doing delta-down, copy w1 for later
            if (deltadown) w2%v = w1%v
            
            ! Scatter the solution to produce the source
            call OpKc%MatVec (w1)
            
            ! Append to the source
            s%v = s%v + w1%v
            
            1 continue
            ! Add delta down terms to s
            ! NOTE - here we do not lump delta down into the uncollided source,
            ! which would otherwise improve convergence. I probably will refactor this
            ! here eventually. The only issue would be I'd end up either introducing
            ! more memory requirement or reading from tape twice for the same array.
            ! Perhaps I should be far more willing to do repeated allocation/deallocation
            ! MHY - Memory is no longer such a big concern so I should just read once and copy to two
            ! arrays
            if (deltadown) then
                XSdd = XS(gp)%s(angular%L + 1,:)
                
                if (ALL(XSdd == ZERO)) go to 2
                
                ! Prepare the coupling operator with the delta down transfer for (gp -> g)
                call OpKc%BuildDeltaDown (XSdd)
                
                ! Scatter the solution to produce the source
                call OpKc%MatVec (w2)
                
                ! Append to the source
                s%v = s%v + w2%v
                
            end if
            2 continue
            if (gtiming) call Timestamp ('EI')
        end do
        
        ! PI
        do pin = 1, p - 1
            if (gtiming) call Timestamp ('PI')
            do gp = 1, SIZE(CXS, dim=1)
                if (.not. ALLOCATED(CXS(gp, pin)%s)) cycle
                
                if (ALL(CXS(gp,pin)%s == ZERO)) cycle
                
                ! Prepare the coupling operator with the transfer moments for (pin, gp) -> (p, g)
                call OpKc%Build (angular, CXS(gp,pin))
                
                ! Read in the solution for particle pin at energy group gp
                call TapeReader (pin, gp, w1%v)
                
                ! Scatter the solution to produce the source
                call OpKc%MatVec (w1)
                
                ! Append to the source
                s%v = s%v + w1%v
            end do
            if (gtiming) call Timestamp ('PI')
        end do
        
        !  =======================
        !    Form the FCS source  
        !  =======================
        
        do pb = 1, SIZE(p2bind(p)%v)
            b = p2bind(p)%v(pb)
            
            if (beams(b)%n <= ZERO) cycle
            
            if (gtiming) call Timestamp ('Creating external beam source')
            call beams(b)%Source (OpL%K, s)
            if (gtiming) call Timestamp ('Creating external beam source')
            
            if (gtiming) call Timestamp ('Creating uncollided angular fluence')
            call beams(b)%AddAngularFluence (OpL%T%XS%t, unc) ! SLOW b/c FluenceQuadrature IS CALLED AGAIN HERE. What to do?
            if (gtiming) call Timestamp ('Creating uncollided angular fluence')
            
        end do
        
        ! !  =========================================
        ! !    Form the boundary + delta down source  
        ! !  =========================================
        
        ! if (input%options%NoFCS) then
        !     !  ===================================================
        !     !    No FCS Case : Form the boundary source directly  
        !     !  ===================================================
            
        !     print *, 'NoFCS WIP'
        !     stop
            
        !     ! Form bdy source
        !     do b = 1, SIZE(beams)
        !         if (.not. validb(b)) cycle
                
        !         if (beams(b)%n <= ZERO) cycle
                
        !         ! call beams(b)%source ()
                
        !     end do
            
        ! else
            
        !     !  ==============================================================
        !     !    FCS Case : Form the primary ('uncollided') angular fluence  
        !     !  ==============================================================
            
        !     !  ----------------------------------------------------------------
        !     !    If raytracing     : This loop gives unc = T^{-1}s^{bdy}
        !     !    If not raytracing : This loop gives unc = s^{bdy}
        !     !  ----------------------------------------------------------------
            
        !     do b = 1, SIZE(beams)
        !         if (.not. validb(b)) cycle
                
        !         if (beams(b)%n <= ZERO) cycle
                
        !         if (raytrace) then
                    
        !             call beams(b)%AngularFluence (attn, w1)
                    
        !             unc%v = unc%v + w1%v
                    
        !         else
                    
        !             call beams(b)%source (unc)
                    
        !         end if
                
        !         unciszero = .FALSE.
                
        !     end do
            
        !     !  ------------------------------------------------------------------
        !     !    If not raytracing : You have unc ~ s, complete unc = T^{-1}s    
        !     !                        If unc was not written to, nothing happens  
        !     !  ------------------------------------------------------------------
            
        !     if (.not. unciszero .and. .not. raytrace) call OpL%T%MatInv (unc)
            
        !     !  --------------------------------------------------------------
        !     !    At this stage you have the uncollided fluence.              
        !     !    Now scatter it and add it to the generic PI and EI sources  
        !     !  --------------------------------------------------------------
            
        !     w1%v = unc%v
            
        !     if (.not. unciszero) call OpL%K%MatVec (w1)
            
        !     s%v = s%v + w1%v
            
        ! end if
        
    End Subroutine
    
    ! A brief subroutine to calculate deposition maps as requested
    Subroutine CalculateDeposition (request, NP, Gs, mesh, w, s)
        Implicit None
        Type (PostProcInput),      Intent (In)    :: request
        Integer,                   Intent (In)    :: NP
        Integer,                   Intent (In)    :: Gs  (:)
        Type (MeshClass),          Intent (In)    :: mesh
        Real (KREAL),              Intent (In)    :: w   (:)
        Type (SpaceAngleVector),   Intent (InOut) :: s
        
        Integer                                   :: g, p
        Type (XSType)                             :: XS
        Real (KREAL), Allocatable                 :: fl    (:)
        Real (KREAL), Allocatable                 :: XSDEP (:)
        Real (KREAL), Allocatable                 :: EDEP  (:)
        Real (KREAL), Allocatable                 :: DDEP  (:)
        Real (KREAL), Allocatable                 :: CDEP  (:)
        
        if (request%energy) then
            ALLOCATE(EDEP(mesh%NENK))
            EDEP = ZERO
        end if
        if (request%dose) then
            ALLOCATE(DDEP(mesh%NENK))
            DDEP = ZERO
        end if
        if (request%charge) then
            ALLOCATE(CDEP(mesh%NENK))
            CDEP = ZERO
        end if
        ALLOCATE(XSDEP(mesh%Nm))
        
        do p = 1, NP
            
            do g = 1, Gs(p)
                call ReadXSLibrary (p, p, g, g, XSLib, XS)
                
                call TapeReader (p, g, s%v)
                
                fl = s%Fluence (w)
                
                if (request%energy) then
                    XSDEP = XS%ED
                    call TallyDeposition (mesh, XSDEP, fl, EDEP)
                end if
                if (request%dose) then
                    XSDEP = XS%DD
                    call TallyDeposition (mesh, XSDEP, fl, DDEP)
                end if
                if (request%charge) then
                    XSDEP = XS%CD
                    call TallyDeposition (mesh, XSDEP, fl, CDEP)
                end if
                
            end do
            
        end do
        
        if (request%energy) then
            write (iuout,*) 'WRITING ENERGY DEPOSITION ---'
            call AddDataset (DATASET_EDEP // PATTERN_DEP, EDEP)
            write (iuout,*) 'DONE.'
        end if
        if (request%dose) then
            write (iuout,*) 'WRITING DOSE DEPOSITION ---'
            call AddDataset (DATASET_DDEP // PATTERN_DEP, DDEP)
            write (iuout,*) 'DONE.'
        end if
        if (request%charge) then
            write (iuout,*) 'WRITING CHARGE DEPOSITION ---'
            call AddDataset (DATASET_CDEP // PATTERN_DEP, CDEP)
            write (iuout,*) 'DONE.'
        end if
        write (iuout,*) ''
        
    End Subroutine

End Program