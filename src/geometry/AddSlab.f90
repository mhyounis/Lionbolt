Submodule (Geometry) submodAddSlab
Contains
Module Subroutine AddSlab (self, NE, T, mat, structure)
    use constants
    use IO
    use Geometry
    Implicit None
    Class (MeshClass),       Intent (InOut) :: self
    
    Integer,                 Intent (In)    :: NE
    Real (KREAL),            Intent (In)    :: T
    Integer,       Optional, Intent (In)    :: mat
    Character (*), Optional, Intent (In)    :: structure
    
    Integer                                 :: e, fg, j, jp, k, m
    Logical,                 Allocatable    :: elmask (:)
    Integer                                 :: umat
    Integer                                 :: NK
    Integer                                 :: NKe
    Integer                                 :: NF
    Integer                                 :: NEold
    Integer                                 :: NKold
    Integer,                 Allocatable    :: matslist (:)
    Real (KREAL),            Allocatable    :: dz       (:)
    Real (KREAL),            Allocatable    :: tmpr     (:,:)
    Character (LEN=:),       Allocatable    :: ustructure
    Type (MeshClass)                        :: blank
    Type (MeshFace),         Allocatable    :: tmpfc    (:)
    Type (MeshElement),      Allocatable    :: tmpel    (:)
    
    ! SLAB MESH: First element is at the bottom. (Beam is coming in from there too)
    
    ! Could forego this by making Slab and Mesh an extension of some other abstract type.
    ! Could be useful.
    if (self%init .and. .not. self%slab) then
        call stophere ('AddSlab.f90: AddSlab: Given mesh object does not correspond to a slab.')
    end if
    
    self%slab = .TRUE.
    
    blank%NE   = NE
    blank%NK   = NE + 1
    blank%NENK = 2 * NE
    
    ALLOCATE(blank%rg(1, blank%NK))
    ALLOCATE(blank%el(NE))
    
    if (PRESENT(mat)) then
        umat = mat
    else
        umat = 1
    end if
    
    if (PRESENT(structure)) then
        ustructure = TRIM(ADJUSTL(Lowercase (structure)))
    else
        ustructure = 'linear'
    end if
    
    !  ==================================
    !    Create this slab's global mesh  
    !  ==================================
    
    ALLOCATE(dz(NE))
    
    select case (ustructure)
    case ('linear')
        dz = T / NE
    case ('logarithmic')
        do e = 1, NE
            dz(e) = ONE / (NE - e + 1)
        end do
        
        dz = dz * T / SUM(dz)
    case default
        call stophere ('AddSlab.f90: AddSlab: Unknown structure given : ' // ustructure)
    end select
    
    do e = 1, NE
        if (e == 1) then
            blank%rg(1,1) = ZERO
        end if
        blank%rg(1,e + 1) = blank%rg(1,e) + dz(e)
    end do
    
    !  ==================================
    !    Generate this slab's face info  
    !  ==================================
    
    blank%NF = blank%NK
    ALLOCATE(blank%face(blank%NF))
    
    do fg = 1, blank%NF
        ALLOCATE(blank%face(fg)%node(1))
        
        blank%face(fg)%bdy     = fg == 1 .or. fg == NE + 1
        blank%face(fg)%node(1) = fg
        blank%face(fg)%NK      = 1
    end do
    
    ! ----- Rewrite this using the MERGE stuff I used below
    blank%face(1)%sharedby(1) = 1
    blank%face(1)%sharedby(2) = 0
    do fg = 2, blank%NF - 1
        blank%face(fg)%sharedby(1) = fg - 1
        blank%face(fg)%sharedby(2) = fg
    end do
    blank%face(blank%NF)%sharedby(1) = blank%NF - 1
    blank%face(blank%NF)%sharedby(2) = 0
    
    !  =====================================
    !    Generate this slab's element info  
    !  =====================================
    
    do e = 1, NE
        blank%el(e)%NK   = 2
        blank%el(e)%mat  = umat
        blank%el(e)%vol  = dz(e)
        blank%el(e)%idx  = 1
        blank%el(e)%NF   = 2
        
        NF = blank%el(e)%NF
        NK = blank%el(e)%NK
        
        ALLOCATE(blank%el(e)%node(NK))
        ALLOCATE(blank%el(e)%face(NF))
        ALLOCATE(blank%el(e)%neighbors(NF))
        ALLOCATE(blank%el(e)%n(1, NF))
        do k = 1, 2
            blank%el(e)%node(k) = e + k - 1
        end do
        
        blank%el(e)%neighbors(1) = e - 1 ! MERGE(e - 1, 0, e > 1)
        blank%el(e)%neighbors(2) = MERGE(e + 1, 0, e < NE)
        
        blank%el(e)%face(1) = e     ! fg = kg
        blank%el(e)%face(2) = e + 1 ! fg = kg
        
        blank%el(e)%n(1,1) = - ONE
        blank%el(e)%n(1,2) =   ONE
    end do
    
    if (self%init) then
        
        ! Code here is awful but works. Slab case is mostly for validation anyway (at least for me it is)
        
        !  =============================
        !    Add this slab to the mesh  
        !  =============================
        
        NKold = SIZE(self%rg, dim=2)
        NEold = SIZE(self%el)
        
        ! First, shift this slab's global mesh by the thickness of the full mesh up to this point
        blank%rg(1,:) = blank%rg(1,:) + self%rg(1, NKold)
        ! Now, directly reallocate self%rg such that blank%rg sits on top of it, but DON'T include
        ! the very first node
        ALLOCATE(tmpr(1, blank%NK + NKold - 1))
        tmpr(1, 1:NKold - 1)                 = self%rg(1,1:NKold - 1)
        tmpr(1, NKold: blank%NK + NKold - 1) = blank%rg(1,1:blank%NK)
        
        call move_alloc (tmpr, self%rg)
        
        self%NK = SIZE(self%rg, dim=2)
        
        ! Next, add the faces. Again, shifts must be made
        blank%face(1)%bdy = .FALSE. ! The first face of the new slab is no longer on the boundary
        do fg = 1, blank%NF
            blank%face(fg)%sharedby(:) = blank%face(fg)%sharedby(:) + NEold
            blank%face(fg)%node(:)     = blank%face(fg)%node(:) + NKold - 1
        end do
        blank%face(blank%NF)%sharedby(2) = 0 ! Very last face sharedby is still zero
        
        ALLOCATE(tmpfc(NKold - 1 + blank%NF)) ! The -1 accounts for the fact that the boundary face of the old slab is getting replaced
        tmpfc(1:NKold - 1)         = self%face(1:NKold - 1)
        tmpfc(NKold : SIZE(tmpfc)) = blank%face
        
        call move_alloc (tmpfc, self%face)
        
        self%NF = SIZE(self%face)
        
        ! Finally, add the elements. Several indices within el need to be shifted first.
        do e = 1, NE
            blank%el(e)%node(:)      = blank%el(e)%node(:) + NKold - 1
            blank%el(e)%face(:)      = blank%el(e)%face(:) + NKold - 1
            blank%el(e)%neighbors(:) = blank%el(e)%neighbors(:) + NEold
        end do
        blank%el(NE)%neighbors(2)   = 0 ! Very last neighbor should still be zero, not 0 + NEold
        self%el(NEold)%neighbors(2) = NEold + 1 ! Very last neighbor should be first neighbor in new slab
        
        ALLOCATE(tmpel(NE + NEold))
        tmpel(1:NEold)                = self%el(1:NEold)
        tmpel(NEold + 1 : NEold + NE) = blank%el(1:NE)
        
        call move_alloc (tmpel, self%el)
        
        self%NE = SIZE(self%el)
        
        ! Lastly, reset NENK. Relationship is always NE * 2
        self%NENK = 2 * self%NE
        
    else
        
        self%homog = .TRUE.
        self%NK    = blank%NK
        self%NF    = blank%NF
        self%NE    = blank%NE
        self%NENK  = blank%NENK
        call move_alloc (blank%rg,   self%rg)
        call move_alloc (blank%face, self%face)
        call move_alloc (blank%el,   self%el)
        
    end if
    
    !  ================================
    !    Some general post-processing  
    !  ================================
    
    ALLOCATE(matslist, source = PACK([(self%el(e)%mat, e=1, self%NE)], mask=.TRUE.))
    self%Nm = MAXVAL(matslist)
    
    ! Rebuild the mat2sd array (again, wasteful, but fine)
    if (ALLOCATED(self%mat2sd)) DEALLOCATE(self%mat2sd)
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
    
    ! The below two blocks are very inefficient but whatever, slab is just for validation anyway
    if (ALLOCATED(self%offset)) DEALLOCATE(self%offset)
    ALLOCATE(self%offset(self%NE + 1))
    j = 0
    do e = 1, self%NE
        self%offset(e) = j
        j = j + self%el(e)%NK
    end do
    self%offset(self%NE + 1) = j
    
    ! Determine connectivity array
    if (ALLOCATED(self%connectivity)) DEALLOCATE(self%connectivity)
    ALLOCATE(self%connectivity(self%NENK))
    do e = 1, self%NE
        self%connectivity(self%offset(e) + 1 : self%offset(e) + self%el(e)%NK) = self%el(e)%node
    end do
    
    self%init  = .TRUE.
    
End Subroutine
End Submodule