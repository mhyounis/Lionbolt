Submodule (AngularSpace) submodSNPrep
Contains
Module Function SNPrep (L, solver, slab) Result (self)
    use constants
    use IO
    use AngularSpace
    use Quadrature
    use Polynomials
    Implicit None
    Integer,            Intent (In) :: L
    Character (*),      Intent (In) :: solver
    Logical,            Intent (In) :: slab
    
    Type (AngularClass)             :: self
    
    Integer                         :: ell, i, ii, ip, j
    Integer                         :: Ng
    Integer                         :: Nc
    Integer                         :: N
    Real (KREAL),       Allocatable :: x     (:)
    Real (KREAL),       Allocatable :: w     (:)
    Real (KREAL),       Allocatable :: xdoty (:)
    Real (KREAL),       Allocatable :: P     (:,:)
    
    !  ==============
    !    Discretize  
    !  ==============
    
    ! Make assignments
    self%disc = 'SN'
    self%L    = L
    
    ! Set the solver
    self%solver = solver
    if (.not. streq(solver, 'SI') .and. .not. streq(solver, 'GMRES')) then
        call stophere ('SNPrep.f90 : Unknown solver given to SNPrep.')
    end if
    
    Ng = L + 1
    
    ! Polar angle
    call PrepGaussLegendreQuadrature (Ng, x, w)
    
    call move_alloc (x, self%mu)
    call move_alloc (w, self%wG)
    
    if (.not. slab) then
        !   Azimuth
        Nc = 2 * Ng
        
        self%NI = Nc * Ng
        
        self%wC = TWOPI / Nc
        
        ALLOCATE(self%phi(Nc))
        
        do j = 1, Nc
            self%phi(j) = (2 * j - 1) * PI / Nc
        end do
        
        !   Propagation direction
        ALLOCATE(self%k(3, self%NI))
        ALLOCATE(self%w(self%NI))
        
        ii = 0
        do j = 1, Nc
            do i = 1, Ng
                ii = ii + 1
                
                self%k(1, ii) = COS(self%phi(j)) * SQRT(ONE - self%mu(i)**2)
                self%k(2, ii) = SIN(self%phi(j)) * SQRT(ONE - self%mu(i)**2)
                self%k(3, ii) = self%mu(i)
                
                self%w(ii) = self%wG(i) * self%wC
            end do
        end do
        
    else
        self%NI = Ng
        
        ALLOCATE(self%w, source = self%wG)
        
        ALLOCATE(self%k(1, self%NI))
        
        self%k(1,1:Ng) = self%mu(1:Ng)
    end if
    
    !  ===============================
    !    Prepare scattering operator  
    !  ===============================
    ! 
    ! Note so I can stop thinking about this
    ! I decided to put this here rather than in the scattering build routines
    ! because those use this as a kernel for the full scattering operator.
    ! So it's really important to have this stored outside of the actual
    ! operator whose cross section can then easily be swapped out.
    ! 
    ! Note again
    ! No, I should move this to K and just let K store the kernel.
    ! The only issue with that is that Lionbolt uses two instances of a scattering operator
    ! (coupling and the one in L), but that could possibly be changed anyway. I don't want to
    ! store this twice. Anyway this is important because of how I'm planning to do MPI.
    
    N = self%NI
    
    if (.not. slab) then
        ALLOCATE(xdoty(N * (N + 1) / 2))
        
        ! Form k \dot k'
        do ip = 1, N
            xdoty(ip + (ip - 1) * ip / 2) = ONE ! k vectors are unit length, so k \dot k' = 1
            do i = 1, ip - 1
                xdoty(i + (ip - 1) * ip / 2) = DOT_PRODUCT(self%k(:,i), self%k(:,ip))
            end do
        end do
        
        call LegendrePolynomials (L, xdoty, P)
        
        ALLOCATE(self%YY(N, N, 0:L))
        
        do ell = 0, L
            do ip = 1, N
                self%YY(ip,ip,ell) = (2 * ell + 1) / FOURPI ! Diagonal has k \dot k' = 1, so P(1) = 1
                do i = 1, ip - 1
                    self%YY(ip,i,ell) = (2 * ell + 1) * P(i + (ip - 1) * ip / 2, ell) / FOURPI
                    self%YY(i,ip,ell) = self%YY(ip,i,ell)
                end do
            end do
        end do
        
        do ip = 1, N
            self%YY(ip,:,:) = self%YY(ip,:,:) * self%w(ip)
        end do
        
    else
        call LegendrePolynomials (L, self%mu, P)
        
        ALLOCATE(self%YY(N, N, 0:L))
        
        do ell = 0, L
            do ip = 1, N
                do i = 1, N
                    self%YY(ip,i,ell) = (2 * ell + 1) * self%w(ip) * P(i,ell) * P(ip,ell) / 2
                end do
            end do
        end do
        
    end if
    
End Function
End Submodule