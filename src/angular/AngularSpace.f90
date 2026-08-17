Module AngularSpace
    use constants
    Implicit None
    
    Type AngularClass
        Integer                        :: L = -1      ! Number of Gauss-Legendre ordinates. This determines SN Chebyshev ordinates (2N) as well as PN/PN scattering Legendre order (L = N-1).
        Integer                        :: NI          ! Total number of discrete ordinates, regardless of method, regardless of problem type (SN). OR Total number of spherical harmonic moments (PN)
        Real (KREAL)                   :: wC          ! Chebyshev quadrature weight (Always TWOPI / 2 * N)
        Real (KREAL),      Allocatable :: mu  (:)     ! Polar angle cosine abscissae
        Real (KREAL),      Allocatable :: wG  (:)     ! Polar angle cosine weights
        Real (KREAL),      Allocatable :: w   (:)     ! Full angular space weights
        Real (KREAL),      Allocatable :: phi (:)     ! Azimuth
        Real (KREAL),      Allocatable :: k   (:,:)   ! Propagation direction vector ! (dir, i)
        Real (KREAL),      Allocatable :: YY  (:,:,:) ! Scattering operator angular kernel
        Character (LEN=:), Allocatable :: disc        ! Discretization method (SN or PN) (PN is WIP)
        Character (LEN=:), Allocatable :: solver      ! Solver to use for particle p ! (p) ! NOTE - this is included here because coarse solves may use different solvers
    Contains
        Procedure :: Destroy
        Procedure :: Print
    End Type
    
    Private :: SNPrep, PNPrep, Destroy
    
    Interface AngularClass
        Module Procedure SNPrep
        Module Procedure PNPrep
    End Interface
    
    Interface
        
        Module Function SNPrep (L, solver, slab) Result (self)
            use constants
            use Quadrature
            use Polynomials
            Implicit None
            Integer,            Intent (In) :: L
            Character (*),      Intent (In) :: solver
            Logical,            Intent (In) :: slab
            Type (AngularClass)             :: self
        End Function
        
        Module Function PNPrep () Result (self)
            use constants
            Implicit None
            Type (AngularClass) :: self
        End Function
        
    End Interface
    
    Contains
    
    Subroutine Print (self, iuf)
        Implicit None
        Class (AngularClass), Intent (In) :: self
        Integer,              Intent (In) :: iuf ! Unit to print to
        
        Integer                           :: i, iDir
        Logical                           :: DoSN
        Logical                           :: DoPN
        Logical                           :: slab
        Character (5)                     :: istr
        Character (30)                    :: kstr (3)
        Character (300)                   :: line
        
        DoSN = self%disc == 'SN'
        DoPN = self%disc == 'PN'
        slab = SIZE(self%k, DIM=1) == 1
        
        if (DoSN) then
            
            write (iuf,'(4X,A)') 'Ordinates : '
            if (.not. slab) then
                
                write (iuf,'(4X,A)') 'index )' // REPEAT(' ', 31) //  '(𝑘𝑥, 𝑘𝑦, 𝑘𝑧)'
                write (iuf,'(4X,A)')  REPEAT('‾', 78)
                
                do i = 1, self%NI
                    
                    do iDir = 1, 3
                        write (kstr(iDir), '(F19.16)') self%k(iDir,i)
                    end do
                    
                    write(istr,'(I5)') i
                    
                    line = TRIM(istr) // ' )    ( ' // TRIM(kstr(1)) // ', ' // TRIM(kstr(2)) // ', ' // TRIM(kstr(3)) // ' )'
                    
                    write (iuf,'(4X,A)') TRIM(line)
                    
                end do
                
            else
                
                write (iuf,'(4X,A)') 'index )  𝜇 = cos𝜃'
                write (iuf,'(4X,A)') REPEAT('‾', 29)
                
                do i = 1, self%NI
                    
                    write(kstr(1),'(F19.16)') self%k(1,i)
                    
                    write(istr,'(I5)') i
                    
                    line = TRIM(istr) // ' ) ' // TRIM(kstr(1))
                    
                    write (iuf,'(4X,A)') TRIM(line)
                    
                end do
                
            end if
            
            write (iuf,*) ''
            
        else if (DoPN) then
            
            ! WIP of course
            
        end if
        
    End Subroutine
    
    Subroutine Destroy (self)
        Class (AngularClass), Intent (InOut) :: self
        
        self%L = -1
        if (ALLOCATED(self%mu))     DEALLOCATE(self%mu)
        if (ALLOCATED(self%wG))     DEALLOCATE(self%wG)
        if (ALLOCATED(self%w))      DEALLOCATE(self%w)
        if (ALLOCATED(self%phi))    DEALLOCATE(self%phi)
        if (ALLOCATED(self%k))      DEALLOCATE(self%k)
        if (ALLOCATED(self%YY))     DEALLOCATE(self%YY)
        if (ALLOCATED(self%disc))   DEALLOCATE(self%disc)
        if (ALLOCATED(self%solver)) DEALLOCATE(self%solver)
        
    End Subroutine
    
End Module