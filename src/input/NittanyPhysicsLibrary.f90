Submodule (InputModule) submodNittanyPhysicsLibrary
Contains
Module Subroutine NittanyPhysicsLibrary (L, mats, particles, energy, XSLib)
    !  =====================================================================
    !    This subroutine generates a physics library using NittanyPhysics.  
    !  =====================================================================
    use constants
    use NittanyAPI
    use LionboltAPI
    use InputModule
    Implicit None
    Integer,                            Intent (In)  :: L
    Type (MaterialSet),                 Intent (In)  :: mats
    Type (ParticleInput),               Intent (In)  :: particles (:) ! MHY LATER - Hate this variable and class naming... I don't like how many different particle classes there are, when none of them are actually important to Lionbolt itself... very bad design
    Type (EnergyGrid),                  Intent (In)  :: energy    (:)
    
    Type (XSLibrary),      Allocatable, Intent (Out) :: XSLib     (:,:)
    
    Integer                                         :: g, gp, m, pin, pout
    Integer                                         :: Nm
    Integer                                         :: NP
    Class (ParticleClass), Allocatable              :: p
    Type ParticleSet
        Class (ParticleClass), Allocatable :: p
    End Type
    Type (ParticleSet),    Allocatable              :: particleObjs (:)
    
    Nm = mats%Nm
    NP = SIZE(particles)
    
    ALLOCATE(XSLib(NP, NP))
    ALLOCATE(particleObjs(NP))
    
    do pout = 1, NP
        ! Do pin = pout case first
        XSLib(pout,pout)%names(1)%str = particles(pout)%name
        XSLib(pout,pout)%names(2)%str = particles(pout)%name
        
        ! Allocate the particleObj to the given particle and then generate the physics
        select case (particles(pout)%identity)
        case ('electrons')
            XSLib(pout,pout)%hasddown = .TRUE.
            ALLOCATE(Electrons :: p) ! Just used as a short-hand object
        case ('photons')
            ALLOCATE(Photons   :: p)
        end select
        
        call p%Build (L, mats, energy(pout), XSLib(pout,pout)%hasddown)
        
        ! Apparently must separately visit a select type block if I want to
        ! access the type-bound procedures of an abstract type object that
        ! was allocated to a concrete extension
        select type (p)
        type is (Electrons)
            
            call p%Elastic        ()
            call p%Moller         ()
            ! Bremsstrahlung is not working at the moment
            ! call p%Bremsstrahlung ()
            
            ALLOCATE(particleObjs(pout)%p, source=p)
            ALLOCATE(XSLib(pout,pout)%XS, source = p%XS)
            
        type is (Photons)
            
            call p%Compton         ()
            call p%Photoabsorption ()
            call p%PairAbsorption  ()
            
            ALLOCATE(particleObjs(pout)%p, source=p)
            ALLOCATE(XSLib(pout,pout)%XS, source = p%XS)
            
        end select
        
        ! Do coupling cross sections
        do pin = 1, pout - 1
            
            XSLib(pin,pout)%names(1)%str = particles(pin)%name
            XSLib(pin,pout)%names(2)%str = particles(pout)%name
            
            call Couple (particleObjs(pin)%p, p, XSLib(pin,pout)%XS)
            
        end do
        
        call p%Destroy ()
        DEALLOCATE(p)
        
    end do
    
    ! If noscat is true for any particles, take the scattering cross section to zero.
    do pout = 1, NP
        if (particles(pout)%noscat) then
            do g = 1, SIZE(XSLib(pout, pout)%XS, dim=2)
                do gp = 1, g
                    XSLib(pout, pout)%XS(gp,g)%s = ZERO
                end do
            end do
        end if
    end do
    
    Contains
    
    Subroutine Couple (p1, p2, XS)
        Implicit None
        Class (ParticleClass),      Intent (In)    :: p1       ! Incident particle
        
        Class (ParticleClass),      Intent (InOut) :: p2       ! Outgoing particle
        Type (XSType), Allocatable, Intent (Out)   :: XS (:,:)
        
        select type (p1)
        type is (Electrons)
            
            select type (p2)
            type is (Electrons)
                ! Electrons in electrons out - do not generate XS object
            type is (Photons)
                ! Electrons in photons out - these mechanisms are still WIP
            end select
            
        type is (Photons)
            
            select type (p2)
            type is (Electrons)
                ! Photons in electrons out
                
                if (ALLOCATED(p2%PhotonInt)) DEALLOCATE(p2%PhotonInt)
                call p2%Compton        (p1)
                call p2%Photoelectrons (p1)
                call p2%PairProduction (p1)
                
                ALLOCATE(XS, source=p2%PhotonInt)
                
            type is (Photons)
                ! Photons in photons out - do not generate XS object
            end select
            
        end select
        
    End Subroutine
    
End Subroutine
End Submodule