Submodule (Sources) submodFluenceQuadrature
Contains
Module Function FluenceQuadrature (self, attn) Result (phi0q)
    use constants
    use profiler
    use Sources
    Implicit None
    Class (ExternalBeam), Intent (InOut) :: self
    Real (KREAL),         Intent (In)    :: attn  (:)
    
    Real (KREAL),         Allocatable    :: phi0q (:)
    
    Integer                              :: iRay, q
    Integer                              :: NR
    
    if (.not. ALLOCATED(self%rays))  call self%RayTrace    ()
    if (.not. ALLOCATED(self%phi0s)) call self%PrepFluence ()
    ! call profile ('FluenceQuadrature')
    NR = SIZE(self%rays)
    
    ALLOCATE(phi0q(self%mesh%NENK), source=self%phi0s)
    
    ! Now just go through and add Beer's law factor
    do iRay = 1, NR
        ! Visit a ray (mesh sdof)
        q = self%rays(iRay)%q
        
        phi0q(q) = phi0q(q) * EXP( - DOT_PRODUCT(attn, self%rays(iRay)%pathlen) )
    end do
    
    ! Apply population factor
    phi0q = self%n * phi0q
    
End Function
End Submodule