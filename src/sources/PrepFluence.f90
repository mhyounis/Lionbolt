Submodule (Sources) submodPrepFluence
Contains
Module Subroutine PrepFluence (self)
    ! Stores the spatial profile of the uncollided fluence.
    ! This can be re-used if energy iteration is being performed,
    ! preventing the neeed to continually reconstruct these values.
    use constants
    use IO
    use profiler
    use Quadrature
    use Sources
    Implicit None
    Class (ExternalBeam), Intent (InOut) :: self
    
    Integer                              :: e, eb, kq, q
    Integer                              :: NE
    Integer                              :: NKe
    Integer,              Allocatable    :: einb (:)
    Real (KREAL)                         :: R0   (3)
    Real (KREAL)                         :: wq   (4)
    Real (KREAL)                         :: re   (3, 4)
    Real (KREAL)                         :: cq   (3, 4)
    Real (KREAL)                         :: v    (3, 4)
    Real (KREAL) :: tmp
    
    if (ALLOCATED(self%phi0s)) then
        call stophere ('PrepFluence.f90: PrepFluence: PrepFluence was called, but this ExternalBeam ' // &
                       'has already done fluence preparation.')
    end if
    
    ALLOCATE(self%phi0s(self%mesh%NENK), source=ZERO)
    
    NE = self%mesh%NE
    
    if (self%fldgeo%slab) then
        
        !  =============
        !    SLAB CASE  
        !  =============
        
        self%phi0s = ONE
        
    else if (streq(self%fldgeo%angdist, 'planar')) then
        
        !  ====================
        !    PLANAR BEAM CASE  
        !  ====================
        
        ALLOCATE(einb, source=PACK([(e,e=1,NE)], mask=self%elemmask))
        
        do eb = 1, SIZE(einb)
            e = einb(eb)
            
            do kq = 1, 4 ! ----- HERE IS WHERE THE TETRAHEDRAL CASE IS ENFORCED -----
                q = self%mesh%offset(e) + kq
                if (.not. self%quadmask(q)) cycle
                
                self%phi0s(q) = ONE
            end do
        end do
        
    else if (streq(self%fldgeo%angdist, 'spherical')) then
        
        !  =======================
        !    SPHERICAL BEAM CASE  
        !  =======================
        
        ALLOCATE(einb, source=PACK([(e,e=1,NE)], mask=self%elemmask))
        
        R0 = self%fldgeo%origin
        ! ----- HERE IS WHERE THE TETRAHEDRAL CASE IS ENFORCED ----- (AND QUADRATURE ORDER...)
        call PrepTetrahedralQuadrature (cq, wq)
        
        do eb = 1, SIZE(einb)
            e = einb(eb)
            
            re(1:3, 1:4) = self%mesh%rg(1:3, self%mesh%el(e)%node(1:4))
            v = MapToTetrahedralQuad (cq, re)
            
            do kq = 1, 4
                q = self%mesh%offset(e) + kq
                if (.not. self%quadmask(q)) cycle
                
                self%phi0s(q) = ONE / DOT_PRODUCT(v(:,kq) - R0, v(:,kq) - R0)
            end do
        end do
        
        ! Normalize spherical case by multiplying by S^2
        self%phi0s = self%phi0s * (self%fldgeo%cutoutparams(1)**2)
        
    end if
    
End Subroutine
End Submodule