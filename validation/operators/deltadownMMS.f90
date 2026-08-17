Program deltadownMMS
    use constants
    use HDF5Tree
    use NittanyAPI
    use LionboltAPI
    Implicit None
    
    Type (MeshClass),    Target :: mesh
    Type (AngularClass), Target :: angular
    Type (ScatteringOp)         :: OpK
    
    Type (SpaceAngleVector)     :: s
    
    Integer                     :: e, ell, i, k, sdof
    Integer                     :: NE
    Integer                     :: L
    Real (KREAL),   Allocatable :: XSdd (:)
    
    ! Discretization parameters
    NE = 6
    
    ! Initialize geometry
    call mesh%AddSlab (NE=NE / 2, T=HALF * NE, mat=1)
    call mesh%AddSlab (NE=NE / 2, T=HALF * NE, mat=2)
    call mesh%PostProcess ()
    
    L = 1
    angular = AngularClass (L=L, solver='GMRES', slab=.FALSE.)
    
    ALLOCATE(XSdd(2))
    
    XSdd(1) = ONE
    XSdd(2) = TEN
    
    call OpK%Build (mesh)
    call OpK%Build (mesh, angular)
    call OpK%BuildDeltaDown (XSdd)
    
    s = SpaceAngleVector (mesh, angular)
    s%v = ONE
    
    call OpK%MatVec (s)
    
    do i = 1, angular%NI
        do sdof = 1, mesh%NENK
            print *, i, sdof, s%v(sdof + (i - 1) * mesh%NENK)
        end do
    end do
    
End Program