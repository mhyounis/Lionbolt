Submodule (Geometry) submodPostProcess
Contains
Module Subroutine PostProcess (self)
    use constants
    use BasicMathFunctions
    use Geometry
    use Parallelism
    Implicit None
    Class (MeshClass), Intent (InOut) :: self
    
    Integer                           :: e, f, jp, k, m
    Logical                           :: muststop = .FALSE.
    Logical,           Allocatable    :: done   (:)
    Logical,           Allocatable    :: elmask (:)
    Integer                           :: NE
    Integer                           :: NF
    Integer                           :: fg
    Integer                           :: k1
    Integer                           :: k2
    Integer                           :: k3
    Integer                           :: ep
    Integer                           :: fp
    Integer                           :: reg1
    Integer                           :: reg2
    Integer                           :: NKe
    Integer                           :: nodes     (8)
    Real (KREAL)                      :: normfac
    Real (KREAL)                      :: tmp
    Real (KREAL)                      :: r1        (3)
    Real (KREAL)                      :: r2        (3)
    Real (KREAL)                      :: r3        (3)
    Real (KREAL)                      :: n         (3)
    Real (KREAL)                      :: CoM       (3)
    Real (KREAL)                      :: r         (3, 8)
    
    if (self%READY) then
        call stophere ('PostProcess.f90: PostProcess: Mesh post-processing was requested but this mesh is already READY. ' // &
                       'Destroy this object before reusing.')
    end if
    
    if (self%slab) go to 100 ! For slabs you get the below information while adding slabs, so just skip to the FEM call.
    
    NE = self%NE
    
    ALLOCATE(done(self%NF), source = .FALSE.) ! Used to tell if a local face's normal has already been constructed by its neighbor
    
    !   Face-element connectivity
    !$OMP PARALLEL DO DEFAULT (SHARED) PRIVATE (e, f, fg, NF, k1, k2, k3, ep, fp, CoM, r1, r2, n, normfac) ! Note for parallelizing in the future: Be sure that you're not making too many large arrays private. You can assign the master array inside the loop as shared, but to write independent values you can use temporary arrays
    do e = 1, NE
        NF = self%el(e)%NF
        
        ! Construct the center of the element
        CoM = ZERO
        do k = 1, self%el(e)%NK
            k1  = self%el(e)%node(k)
            CoM = CoM + self%rg(:,k1)
        end do
        CoM = CoM / self%el(e)%NK
        
        do f = 1, NF
            ! Identify the global face
            fg = self%el(e)%face(f)
            
            ! If this face has been treated, cycle
            if (done(fg)) cycle
            
            ! Construct normal. Only need three nodes, any three nodes, ever
            k1  = self%face(fg)%node(1)
            k2  = self%face(fg)%node(2)
            k3  = self%face(fg)%node(3)
            
            r1 = self%rg(:,k2) - self%rg(:,k1)
            r2 = self%rg(:,k3) - self%rg(:,k1)
            
            n = CROSS_PRODUCT(r2, r1)
            
            ! Check sign and normalize
            ! r1 now points from the center of the element TO one of the nodes on the face
            r1 = self%rg(:,k1) - CoM
            
            ! Thus, we want to keep the sign of n if n dot r1 is positive, but
            ! flip it if n dot r1 is negative. So we multiply by sign(n dot r1).
            
            normfac = SIGN(ONE, DOT_PRODUCT(n, r1)) / NORM2(n)
            
            self%el(e)%n(:,f) = normfac * n
            
            ! Now prepare to assign normal opposite to neighbor
            ep = self%el(e)%neighbors(f)
            
            if (ep == 0) cycle
            
            fp = FINDLOC(self%el(ep)%neighbors, e, DIM=1)
            
            self%el(ep)%n(:,fp) = - self%el(e)%n(:,f)
            
            ! Set this global face to done
            done(fg) = .TRUE.
        end do
        
    end do
    !$OMP END PARALLEL DO
    
    if (self%Nm > 1) then
        do fg = 1, self%NF
            ! Assign elements
            e  = self%face(fg)%sharedby(1)
            ep = self%face(fg)%sharedby(2)
            
            ! First quickly determine if this face belongs to a region interface.
            if (ep == 0) then
                ! If this is on the boundary assign to a region and cycle
                reg1 = self%b2reg(self%el(e)%block)
                self%region(reg1)%fg = [self%region(reg1)%fg, fg]
                cycle
            end if
            
            ! If not on the boundary determine the regions this face lies between
            reg1 = self%b2reg(self%el(e)%block)
            reg2 = self%b2reg(self%el(ep)%block)
            
            ! If the regions are the same then just pass
            if (reg1 == reg2) cycle
            
            ! Otherwise record this face in both regions
            self%region(reg1)%fg = [self%region(reg1)%fg, fg]
            self%region(reg2)%fg = [self%region(reg2)%fg, fg]
            
        end do
    end if
    
    !   Face area
    ! IS THIS SELECT CASE LOOP SLOW??? If so what could I even do given that idx changes with f???
    ! SHOULD I USE GMSH API FOR THIS STUFF TOO?
    do f = 1, self%NF
        select case (self%face(f)%idx)
        case (2)
            nodes(1:3) = self%face(f)%node
            
            r1 = self%rg(:,nodes(2)) - self%rg(:,nodes(1))
            r2 = self%rg(:,nodes(3)) - self%rg(:,nodes(1))
            
            self%face(f)%area = HALF * norm2(CROSS_PRODUCT(r2, r1))
            
        !case (3)
            ! ASSUMPTION IS THAT NODES ARE ORDERED IN TERMS OF CONNECTIVITY:
            ! 1 - 2
            ! |   |
            ! 4 - 3
            ! This assumption is satisfied in GMSH, BUT NOT FOR GLOBAL FACES. So, figure this out
            
            !nodes(1:4) = self%face(f)%node
            !
            !! We do the triangle formula, but with 1->2, 1->3, and add it to that with 1->4, 1->3
            !r1 = self%rg(:,nodes(3)) - self%rg(:,nodes(1))
            !r2 = self%rg(:,nodes(2)) - self%rg(:,nodes(1))
            !
            !self%face(f)%area = HALF * norm2(CROSS_PRODUCT(r2, r1))
            !
            !r2 = self%rg(:,nodes(4)) - self%rg(:,nodes(1))
            !self%face(f)%area = self%face(f)%area + HALF * norm2(CROSS_PRODUCT(r2, r1))
            
        case default
            print *, 'Developer error: PostProcess.f90: self face with other than 3 nodes given. WIP'
            stop
        end select
        
    end do
    
    ! ! ACCURACY CHECK:
    ! ! Could be rather slow, maybe a bad construct in general, but has saved me
    ! do e = 1, NE
    !     r1 = ZERO
    !     do f = 1, self%el(e)%NF
    !         r1 = r1 + self%face(self%el(e)%face(f))%area * self%el(e)%n(:,f)
    !     end do
    !     if (NORM2(r1) > 1.0e-8_KREAL * SUM(self%face(self%el(e)%face)%area)) then ! Could change the tolerance level to be smarter...
    !         !print *, '------------------'
    !         !print *, e
    !         !
    !         !NF = self%el(e)%NF
    !         !
    !         !! Construct the center of the element
    !         !CoM = ZERO
    !         !do k = 1, self%el(e)%NK
    !         !    k1  = self%el(e)%node(k)
    !         !    CoM = CoM + self%rg(:,k1)
    !         !end do
    !         !CoM = CoM / self%el(e)%NK
    !         !
    !         !print *, 'Center of mass:', CoM
    !         !
    !         !do f = 1, NF
    !         !    print *, 'f:', f, self%face(self%el(e)%face(f))%area, '|', self%el(e)%n(:,f)
    !         !    
    !         !    ! Identify the global face
    !         !    fg = self%el(e)%face(f)
    !         !    print *, 'Global f:', fg
    !         !    ! Construct normal. Only need three nodes, ever
    !         !    k1  = self%face(fg)%node(1)
    !         !    k2  = self%face(fg)%node(2)
    !         !    k3  = self%face(fg)%node(3)
    !         !    print *, 'k values:', k1, k2, k3
    !         !    r1 = self%rg(:,k2) - self%rg(:,k1)
    !         !    r2 = self%rg(:,k3) - self%rg(:,k1)
    !         !    print *, 'r1:', r1
    !         !    print *, 'r2:', r2
    !         !    n = CROSS_PRODUCT(r2, r1)
    !         !    print *, 'n:', n
    !         !    ! Check sign and normalize
    !         !    ! r1 now points from the center of the element TO one of the nodes on the face
    !         !    r1 = self%rg(:,k1) - CoM
    !         !    print *, 'r1 to CoM:', r1
    !         !    ! Thus, we want to keep the sign of n if n dot r1 is positive, but
    !         !    ! flip it if n dot r1 is negative. So we multiply by sign(n dot r1).
    !         !    
    !         !    normfac = SIGN(ONE, DOT_PRODUCT(n, r1)) / norm2(n)
    !         !    print *, 'normfac:', normfac
    !         !    self%el(e)%n(:,f) = normfac * n
    !         !    print *, 'final:', self%el(e)%n(:,f)
    !         !    print *, 'norm of final:', norm2(self%el(e)%n(:,f))
    !         !end do
    !         !
    !         !do f = 1, self%el(e)%NF
    !         !    print *, 'f:', f, self%face(self%el(e)%face(f))%area, '|', self%el(e)%n(:,f)
    !         !end do
    !         !print *, '------------------'
    !         if (.not. muststop) write (iuout,*) ''
    !         write (iuout,*) 'PostProcess.f90: Element area vector found exceeding tolerance: ', e
    !         write (iuout,*) 'Norm of sum (should be zero, or less than 1.0e-8): ', norm2(r1)
    !         muststop = .TRUE.
    !     end if
    ! end do
    
    if (muststop) then
        write (iuout,*) ''
        call stophere ('PostProcess.f90: Some element area vectors were found exceeding tolerance. ' &
                        // 'Verify your mesh.')
    end if
    
    !   Element volume
    ! INSTEAD OF THIS SELECT CASE LOOP SHOULD I HAVE idx2el?
    do e = 1, NE
        select case (self%el(e)%idx)
        case (4)
            ! Tetrahedron
            ! Faces:
            ! 1 2 3
            ! 1 2 4
            ! 1 3 4
            ! 2 3 4
            !
            ! ! From GMSH website:
            ! !
            ! !                      3           
            ! !                    ,/|`\         
            ! !                  ,/  |  `\       
            ! !                ,/    '.   `\     
            ! !              ,/       |     `\   
            ! !            ,/         |       `\ 
            ! !         1 <-----------'.--------> 2
            ! !            `\.         |      ,/ 
            ! !               `\.      |    ,/   
            ! !                  `\.   '. ,/     
            ! !                     `\. |/       
            ! !                        `' 4
            
            ! One sixth the triple product of the vectors originating at 1
            NKe = 4
            
            nodes(1:4) = self%el(e)%node
            
            r(1:3,1:NKe) = self%rg(:,nodes(1:4))
            
            r1 = r(1:3,2) - r(1:3,1)
            r2 = r(1:3,3) - r(1:3,1)
            r3 = r(1:3,4) - r(1:3,1)
            
            self%el(e)%vol = SIXTH * ABS(TRIPLE_PRODUCT(r3, r2, r1))
            
        case (5)
            ! Hexahedron
            ! Faces:
            ! 1 2 3 4
            ! 5 6 7 8
            ! 1 2 6 5
            ! 3 4 8 7
            ! 1 4 8 5
            ! 2 3 7 6
            !
            ! ! From GMSH website:
            ! !
            ! !            4----------3    
            ! !            |\         |\   
            ! !            | \        | \  
            ! !            |  \       |  \ 
            ! !            |   8------+---7
            ! !            |   |      |   |
            ! !            1---+------2   |
            ! !             \  |       \  |
            ! !              \ |        \ |
            ! !               \|         \|
            ! !                5----------6
            !
            ! Tetrahedra are: 1 2 5 4
            !                 6 5 2 7
            !                 8 7 4 5
            !                 3 7 4 2
            !                 4 7 5 2
            
            NKe = 8
            
            nodes(1:8) = self%el(e)%node
            
            r(1:3,1:NKe) = self%rg(:,nodes(1:8))
            
            r1 = r(:,2) - r(:,1)
            r2 = r(:,5) - r(:,1)
            r3 = r(:,4) - r(:,1)
            
            tmp = SIXTH * ABS(TRIPLE_PRODUCT(r3, r2, r1))
            
            r1 = r(:,5) - r(:,6)
            r2 = r(:,2) - r(:,6)
            r3 = r(:,7) - r(:,6)
            
            tmp = tmp + SIXTH * ABS(TRIPLE_PRODUCT(r3, r2, r1))
            
            r1 = r(:,7) - r(:,8)
            r2 = r(:,4) - r(:,8)
            r3 = r(:,5) - r(:,8)
            
            tmp = tmp + SIXTH * ABS(TRIPLE_PRODUCT(r3, r2, r1))
            
            r1 = r(:,7) - r(:,3)
            r2 = r(:,4) - r(:,3)
            r3 = r(:,2) - r(:,3)
            
            tmp = tmp + SIXTH * ABS(TRIPLE_PRODUCT(r3, r2, r1))
            
            r1 = r(:,7) - r(:,4)
            r2 = r(:,5) - r(:,4)
            r3 = r(:,2) - r(:,4)
            
            tmp = tmp + SIXTH * ABS(TRIPLE_PRODUCT(r3, r2, r1))
            
            self%el(e)%vol = tmp
            
        case (6)
            ! Prism
            ! Faces:
            ! 1 2 3
            ! 4 5 6
            ! 1 2 5 4
            ! 2 3 6 5
            ! 3 1 4 6
            !
            ! ! From GMSH website:
            ! !
            ! !                   
            ! !                   
            ! !                  4        
            ! !                ,/|`\      
            ! !              ,/  |  `\    
            ! !            ,/    |    `\  
            ! !           5------+------6 
            ! !           |      |      | 
            ! !           |      |      | 
            ! !           |      |      | 
            ! !           |      |      | 
            ! !           |      |      | 
            ! !           |      1      | 
            ! !           |    ,/ `\    | 
            ! !           |  ,/     `\  | 
            ! !           |,/         `\| 
            ! !           2-------------3 
            
            
            
        case (7)
            ! 5-node pyramid
            ! Faces:
            ! 1 2 3 4
            ! 1 2 5
            ! 2 3 5
            ! 3 4 5
            ! 4 1 5
            !
            ! ! From GMSH website:
            ! !
            ! !                      5            
            ! !                    ,/|\           
            ! !                  ,/ .'|\          
            ! !                ,/   | | \         
            ! !              ,/    .' | `.        
            ! !            ,/      |  '.  \       
            ! !          ,/       .'   |   \      
            ! !        ,/         |    |    \     
            ! !       1----------.'----4    `.    
            ! !        `\        |      `\    \   
            ! !          `\     .'        `\   \
            ! !            `\   |           `\  \ 
            ! !              `\.'             `\ \
            ! !                 2-----------------3
            
            
        case default
            call stophere ('PostProcess.f90: PostProcess: self element contains an unknown index.')
        end select
    
    end do
    
    !! ACCURACY CHECK: (developer only since I don't think it's worth getting the self's volume another way to check against)
    !r1 = ZERO
    !do e = 1, NE
    !    r1(1) = r1(1) + self%el(e)%vol
    !end do
    !print *, r1(1)
    !stop
    
    ! Determine NENK
    self%NENK = 0
    do e = 1, NE
        self%NENK = self%NENK + self%el(e)%NK
    end do
    
    ! Determine connectivity array
    ALLOCATE(self%connectivity(self%NENK))
    do e = 1, NE
        self%connectivity(self%offset(e) + 1 : self%offset(e) + self%el(e)%NK) = self%el(e)%node
    end do
    
    ! Build the material 2 spatial DOF map
    ALLOCATE(self%mat2sd(self%Nm))
    if (self%Nm > 1) then
        
        ALLOCATE(elmask(self%NENK))
        
        do m = 1, self%Nm
            
            elmask = .FALSE.
            jp = 0
            do e = 1, self%NE
                NKe = self%el(e)%NK
                elmask(jp + 1 : jp + NKe) = self%el(e)%mat == m
                jp = jp + NKe
            end do
            
            ALLOCATE(self%mat2sd(m)%v, source = PACK([(e,e=1,self%NENK)], mask = elmask))
        end do
    else
        ALLOCATE(self%mat2sd(1)%v(self%NENK), source = PACK([(e,e=1,self%NENK)], mask = .TRUE.))
    end if
    
    ! Finally, perform finite element analysis and set the mesh to ready
    100 continue
    call FiniteElementAnalysis (self)
    
    self%READY = .TRUE.
    
End Subroutine
End Submodule