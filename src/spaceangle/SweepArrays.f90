Submodule (SpaceAngleInterface) submodSweepArrays
Contains
Module Subroutine SweepArrays (mesh, angular, IsUpstream, IUoffset, sl, isl, IAF)
    use constants
    use Parallelism
    use types
    use Geometry
    use AngularSpace
    use SpaceAngleInterface
    Implicit None
    Type (MeshClass),    Intent (In)  :: mesh
    Type (AngularClass), Intent (In)  :: angular
    Logical,             Intent (In)  :: IsUpstream (:)
    Integer,             Intent (In)  :: IUoffset   (:,:)
    Integer,             Intent (In)  :: sl         (:,:)
    Integer,             Intent (In)  :: isl        (:,:)
    
    Type (UpstreamSys),  Intent (Out) :: IAF
    
    Integer                           :: e, es, esp, eup, f, i, iupf, j
    Integer                           :: NE
    Integer                           :: NI
    Integer                           :: NKe
    Integer                           :: NFe
    Integer                           :: nup
    Integer (KINT)                    :: NZ
    
    NE = mesh%NE
    NI = angular%NI
    
    ! Create space-angle offset array.
    ! IAF%SAoffset(es,i) + 1 gives the value of the first index of a given (es,i) block in a vector using collapsed notation
    ! Inner index gets the extra because we don't jump i values during indexing in the transport
    ALLOCATE(IAF%SAoffset(NE + 1, NI), source=0)
    j = 0
    do i = 1, NI
        do es = 1, NE
            e   = sl(es,i)
            NKe = mesh%el(e)%NK
            
            IAF%SAoffset(es,i) = j
            j = j + NKe
            
        end do
        IAF%SAoffset(NE + 1, i) = j
    end do
    
    ! Create the rest of IAF
    
    ALLOCATE(IAF%UFoffset(NE + 1, NI), source=0)
    NZ = COUNT(IsUpstream)
    ALLOCATE(IAF%upfaces(NZ), source=0)
    ALLOCATE(IAF%ipsf   (NZ), source=0)
    ALLOCATE(IAF%ipef   (NZ), source=0)
    
    ! You visit the valid faces according to IsUpstream and record their values as well.
    iupf = 0 ! Index for NZ
    j    = 0 ! Index for offset construction
    do i = 1, NI
        do es = 1, NE
            e = sl(es,i)
            
            NFe = mesh%el(e)%NF
            
            nup = 0 ! Keep track of how many upstream faces there are here
            do f = 1, NFe
                if (.not. IsUpstream(IUoffset(e,i) + f)) cycle ! If this face is not upstream, just cycle
                nup  = nup + 1
                iupf = iupf + 1
                
                ! All the info is readily present
                ! Local face index of this upstream face (used for normal vector assignment later)
                IAF%upfaces(iupf) = f
                
                ! Now the starting and ending indices of the upstream neighbor in the sweep-space-angle indexing system
                ! Only need the number of nodes in the upstream element
                eup = mesh%el(e)%neighbors(f)
                esp = isl(eup,i)
                NKe = mesh%el(eup)%NK
                
                IAF%ipsf(iupf) = IAF%SAoffset(esp,i) + 1
                IAF%ipef(iupf) = IAF%SAoffset(esp,i) + NKe
                
            end do
            
            ! Construct the offset itself by counting the number of such upstream faces
            IAF%UFoffset(es,i) = j
            j = j + nup
            
        end do
        IAF%UFoffset(NE + 1,i) = j
    end do
    
End Subroutine
End Submodule