Submodule (Sources) submodinit_ExternalBeam
Contains
Module Function init_ExternalBeam (mesh, angular, fldgeo, particle) Result (self)
    use constants
    use IO
    use BasicMathFunctions
    use Quadrature
    use Geometry
    use AngularSpace
    use Sources
    Implicit None
    Type (MeshClass),     Target, Intent (In) :: mesh
    Type (AngularClass),  Target, Intent (In) :: angular
    Type (FieldGeometry),         Intent (In) :: fldgeo
    Character (*),                Intent (In) :: particle
        
    Type (ExternalBeam)                       :: self
    
    Integer                                   :: e, i, ia, j, kq
    Logical                                   :: l
    Logical                                   :: elemisinbeam
    Logical                                   :: NONE      ! Stored logicals for cutout
    Logical                                   :: CIRCLE    ! Stored logicals for cutout
    Logical                                   :: RECTANGLE ! Stored logicals for cutout
    Integer                                   :: NE
    Integer                                   :: NENK
    Integer                                   :: NI
    Integer                                   :: n
    Integer                                   :: m
    Real (KREAL)                              :: dmu
    Real (KREAL)                              :: mui
    Real (KREAL)                              :: muf
    Real (KREAL)                              :: xi
    Real (KREAL)                              :: xf
    Real (KREAL)                              :: dx
    Real (KREAL)                              :: yi
    Real (KREAL)                              :: yf
    Real (KREAL)                              :: dy
    Real (KREAL)                              :: x
    Real (KREAL)                              :: k0   (3)
    Real (KREAL)                              :: wq   (4)
    Real (KREAL)                              :: cq   (3, 4)
    Real (KREAL)                              :: r    (3, 4)
    Real (KREAL)                              :: v    (3, 4)
    Real (KREAL),                 Allocatable :: grid (:,:,:)
    
    ! Make assignments
    self%particle =  particle
    self%mesh     => mesh
    self%angular  => angular
    self%fldgeo   =  fldgeo
    self%n        =  ONE ! Populate with one particle. User can later change using self%Populate and providing energetic information
    
    NE   = mesh%NE
    NENK = mesh%NENK
    NI   = angular%NI
    
    ! Determine the quadrature points inside of the beam
    ! Note, all of this would be facilitated by an octree. Just more reason to implement an octree in Lionbolt too
    if (mesh%slab) then
        ! In the slab case, the quadrature points are determined from order 2 GL quadrature
        ! Everything is in the beam
        ALLOCATE(self%quadmask(NENK), source = .TRUE.)
        ALLOCATE(self%elemmask(NE),   source = .TRUE.)
    else
        ! Loop through the elements and check that quadrature points of an element
        ! are within the beam, this will classify a sdof's quadrature point as
        ! being within the beam
        
        ! We note, this whole block assumes that the number of quadrature points is equal
        ! to the number of nodes, as well as that we are using tetrahedra. 
        ! Will not always be the case. But this will be updated when different finite elements
        ! are implemented
        ALLOCATE(self%quadmask(NENK))
        ALLOCATE(self%elemmask(NE))
        
        call PrepTetrahedralQuadrature (cq, wq)
        
        do e = 1, NE
            
            r(1:3, 1:4) = mesh%rg(1:3, mesh%el(e)%node(1:4))
            v = MapToTetrahedralQuad (cq, r)
            
            elemisinbeam = .FALSE.
            
            ! MHY LATER - generalize max kq of course
            do kq = 1, 4
                l = self%fldgeo%InField (v(1:3,kq))
                self%quadmask(mesh%offset(e) + kq) = l
                ! print *, e, kq, l
                if (l) elemisinbeam = .TRUE. ! If a quadrature node is in the beam, this gets flipped to true. Otherwise, don't touch it
            end do
            
            self%elemmask(e) = elemisinbeam
            
        end do
        
    end if
    
    ! IGNORE THE BELOW. IT IS LARGELY WIP.
    
    !  ==================================================================
    !    Set up source normalization factor on a grid for interpolation  
    !  ==================================================================
    
    ! What is this? --- Here, Lionbolt does something that is perhaps non-standard but should still
    ! be reasonable according to basic analysis. Thus, it requires quite a bit of elaboration.
    ! 
    ! Now, ideally, the sources thus implemented in Lionbolt are described by singular angular
    ! distributions. That means, particles are only populated for angles evaluated PRECISELY at some
    ! angle. In the planar case, this angle is fixed, but in the spherical case, this angle is
    ! position-dependent, and is just the ray traced from the origin of the beam to the physical
    ! point at which one is evaluating the source/any derived quantities (like the fluence of the source).
    ! Nevertheless, we use the discrete ordinates approximation, in which we do not choose our ordinates,
    ! and especially we can not choose them in a position-dependent manner. Thus, we approximate our
    ! source angular distribution (in Lionbolt, we rely on forward-peaked Gaussians).
    ! 
    ! However, though we may normalize these Gaussians, when we evaluate them on a quadrature grid, they
    ! are indeed NOT necessarily normalized consistent with the quadrature grid. Indeed, for too sharp
    ! a Gaussian, you can easily find an equivalent (relative) distribution of fluence across all angular
    ! ordiantes. Thus, your Gaussian becomes a polynomial in the space of the abscissae, and THIS needs to
    ! be normalized. So, in principle, EVERY evaluation of the angular distribution requires ALL discrete 
    ! ordinates in order to perform normalization.
    ! 
    ! When the angular distribution is position-dependent, this becomes quite costly. My solution, which
    ! can be argued to be physically reasonable, is to actually evaluate the position dependence
    ! on a grid on the unit sphere (its position dependence only ever depends on the unit vector
    ! connecting a mesh point and the source origin), and then when one requires evaluation of the 
    ! angular distribution, the normalization is actually provided by interpolating on this grid. This
    ! way, rather than evaluate the normalizaiton factor by evaluating the angular distribution on ALL 
    ! angular ordinates for EVERY mesh point, an expensive task, we just get the normalization factor
    ! by interpolation. 
    ! 
    ! The question is then - does our function, a function on the unit sphere, behave nicely? My argument
    ! is this:
    !   The function we are evaluating is:
    !   
    !   Q(\mathbf{r}) = \sum_{i=1}^{N_{I}} \ w_{i}f(\mathbf{\hat{k}}_{i}, \mathbf{\hat{k}}_{0}(\mathbf{r}))
    !   
    !   where $f(\mathbf{\hat{k}}, \mathbf{\hat{k}}_{0})$ is a Gaussian distribution on a unit sphere, 
    !   centered about the unit vector $\mathbf{\hat{k}}_{0}$. We are simply, therefore, sweeping around
    !   $\mathbf{\hat{k}}_{0}$ as a function of $\mathbf{r}$. 
    !   
    !   Imagine instead that we are sweeping around the a quadrature basis $\mathbf{\hat{k}}_{i}$. Does
    !   Changing the sampling of our function sharply change the value of $Q(\mathbf{r})$? In fact,
    !   $Q(\mathbf{r})$ will change sharply if ordinates are arranged in such a way that, at one evaluation,
    !   a discrete ordinate lands close to $\mathbf{\hat{k}}_{0}$, and then at another evaluation, no 
    !   ordinates land so close to $\mathbf{\hat{k}}_{0}$. However, this can be avoided if the spacing between
    !   ordinates does not exceed the width of the Gaussian --- which is already a requirement, because
    !   otherwise, the Gaussian gets evaluated exclusively in relatively flat regions, and then upon
    !   normalization all ordinates are roughly equally populated.
    
    ! if (mesh%slab .or. streq(self%fldgeo%angdist, 'planar')) then
    !     ! This is position-independent. So, you actually know the normalization factor immediately
    !     ALLOCATE(self%Q0(1,1))
        
    !     self%Q0 = ZERO
        
    !     if (mesh%slab) then
    !         k0(1) = COS(self%fldgeo%axis(1) * PI / 180.0_KREAL)
            
    !         ! FOR LATER WHEN I TRY TO CONTROL s DYNAMICALLY
    !         ! dmu = ABS(ACOS(angular%k(1,1)) - ACOS(angular%k(1,2))) ! This is actually dtheta
            
    !         do ia = 1, NI
    !             self%Q0(1,1) = self%Q0(1,1) + angular%w(ia) * EvalBeamAngDist (angular%k(1,ia), k0(1))
    !         end do
    !     else if (streq(self%fldgeo%angdist, 'planar')) then
    !         k0 = self%fldgeo%axis
            
    !         ! FOR LATER WHEN I TRY TO CONTROL s DYNAMICALLY
    !         ! dmu = ABS(angular%k(3,1) - angular%k(3,2))
            
    !         do ia = 1, NI
    !             self%Q0(1,1) = self%Q0(1,1) + angular%w(ia) * EvalBeamAngDist (angular%k(1:3,ia), k0)
    !         end do
    !     end if
        
    !     ! Return Q0 as the actual inverse of the normalization, to avoid having to take
    !     ! one over this a whole bunch of times (this is NOT done in the spherical case)
    !     self%Q0(1,1) = ONE / self%Q0(1,1)
        
    ! else
    !     ! This is position-dependent, so we set up a grid on the unit sphere on which we can interpolate
        
    !     NONE      = streq(self%fldgeo%cutout, 'none')
    !     CIRCLE    = streq(self%fldgeo%cutout, 'circle')
    !     RECTANGLE = streq(self%fldgeo%cutout, 'rectangle')
        
    !     if (NONE .or. CIRCLE) then
    !         ! Discretize the grid with 8 times as many points in each coordinate as the user's chosen discretization
    !         ! We know that NI is divisible by 2, and that NI / 2 is a perfect square, equaling n = Nmu
    !         n = NINT(SQRT(NI / 2))
            
    !         n = 8 * n
    !         m = 2 * n
            
    !         ALLOCATE(self%Qx(n))
    !         ALLOCATE(self%Qy(m))
    !         ALLOCATE(grid(3,n,m))
            
    !         mui = ONE
    !         if (NONE) then
    !             ! We take the grid to have poles on the z-axis
    !             muf = - ONE
    !         else if (CIRCLE) then
    !             ! For now, we work in the coordinate system where mu is the polar angle
    !             ! with respect to the beam axis. So when we interpolate we must keep this
    !             ! in mind, however, for the grid, we also need to rotate z to the beam
    !             ! axis
    !             x   = self%fldgeo%cutoutparams(2) / self%fldgeo%cutoutparams(1)
    !             muf = ONE / SQRT(ONE + x * x)
    !         end if
    !         dmu = muf - mui ! Keep in mind this is negative, so the below loop uses + dmu
            
    !         do i = 1, n
    !             self%Qx(i) = mui + (2 * i - 1) * HALF * dmu / n
    !         end do
    !         do j = 1, m
    !             self%Qy(j) = (j - 1) * TWOPI / (m - 1)
    !         end do
            
    !         ! Briefly, we flip Qx so that our values are ascending
    !         self%Qx(1:n) = self%Qx(n:1:-1)
            
    !         do j = 1, m
    !             do i = 1, n
    !                 grid(1,i,j) = COS(self%Qy(j)) * SQRT(ONE - self%Qx(i)**2)
    !                 grid(2,i,j) = SIN(self%Qy(j)) * SQRT(ONE - self%Qx(i)**2)
    !                 grid(3,i,j) = self%Qx(i)
    !             end do
    !         end do
            
    !         if (CIRCLE) then
    !             ! Rotate the grid to have the correct beam axis
    !             call RotateGrid (n, m, self%fldgeo%e1, self%fldgeo%e2, self%fldgeo%axis, grid)
    !         end if
            
    !     else if (RECTANGLE) then
    !         ! We create a simple lattice on the cutout
            
    !         ! Discretize the grid with 8 times the angular quadrature in mu, for both axes
    !         n = NINT(SQRT(NI / 2))
            
    !         n = 8 * n
    !         m = n
            
    !         ALLOCATE(self%Qx(n))
    !         ALLOCATE(self%Qy(m))
    !         ALLOCATE(grid(3,n,m))
            
    !         dx = self%fldgeo%cutoutparams(2)
    !         xi = - HALF * dx
    !         xf = - xi
    !         dy = self%fldgeo%cutoutparams(3)
    !         yi = - HALF * dy
    !         yf = - yi
            
    !         do i = 1, n
    !             self%Qx(i) = xi + (i - 1) * dx / (n - 1)
    !         end do
    !         do j = 1, m
    !             self%Qy(i) = xi + (i - 1) * dx / (n - 1)
    !         end do
            
    !         ! We take, prior to transformation, the beam as pointing up. So the 3rd component of each grid point is the source-to-surface distance
    !         grid(3,1:n,1:m) = self%cutoutparams(1)
    !         do i = 1, n
    !             grid(1,i,1:m) = self%Qx(i)
    !         end do
    !         do j = 1, m
    !             grid(2,1:n,j) = self%Qy(j) 
    !         end do
            
    !         ! Finally transform.
    !         call RotateGrid (n, m, self%fldgeo%e1, self%fldgeo%e2, self%fldgeo%axis, grid)
            
    !     end if
        
    !     ! For later
    !     ! dmu = ABS(angular%k(3,1) - angular%k(3,2))
        
    !     ALLOCATE(self%Q0(n, m))
    !     self%Q0 = ZERO
    !     do j = 1, m
    !         do i = 1, n
    !             k0 = (grid(1:3,i,j) - self%fldgeo%origin) / NORM2(grid(1:3,i,j) - self%fldgeo%origin)
    !             do ia = 1, NI
    !                 self%Q0(i,j) = self%Q0(i,j) + angular%w(ia) * EvalBeamAngDist (angular%k(1:3,ia), k0)
    !             end do
    !         end do
    !     end do
        
    ! end if
    
    Contains
    
    Subroutine RotateGrid (n, m, e1, e2, e3, grid)
        use constants
        Implicit None
        Integer,      Intent (In)    :: n
        Integer,      Intent (In)    :: m
        Real (KREAL), Intent (In)    :: e1 (3)
        Real (KREAL), Intent (In)    :: e2 (3)
        Real (KREAL), Intent (In)    :: e3 (3)
        
        Real (KREAL), Intent (InOut) :: grid (3,n,m)
        
        Integer                      :: i, j
        Real (KREAL)                 :: tmp (3)
        
        ! R(:,1) = e1
        ! R(:,2) = e2
        ! R(:,3) = e3
        
        do j = 1, m
            do i = 1, n
                tmp(1:3) = grid(1:3,i,j)
                
                grid(1,i,j) = e1(1) * tmp(1) + e2(1) * tmp(2) + e3(1) * tmp(3)
                grid(2,i,j) = e1(2) * tmp(1) + e2(2) * tmp(2) + e3(2) * tmp(3)
                grid(3,i,j) = e1(3) * tmp(1) + e2(3) * tmp(2) + e3(3) * tmp(3)
            end do
        end do
        
    End Subroutine
    
End Function
End Submodule