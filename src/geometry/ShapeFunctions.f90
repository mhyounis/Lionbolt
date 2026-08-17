Module ShapeFunctions
    use constants
    use IO
    Implicit None
    
    !  =================================================
    !    Linear shape functions for reference elements
    !  =================================================
    
Contains

Function LineSF (k, z) Result (u)
    Implicit None
    Integer,      Intent (In) :: k
    Real (KREAL), Intent (In) :: z
    
    Real (KREAL)              :: u
    
    ! Node 1 is -1, node 2 is +1
    select case (k)
    case (1); u = HALF * (ONE - z)
    case (2); u = HALF * (ONE + z)
    end select
    
End Function

Function GradLineSF (k, z) Result (Du)
    Implicit None
    Integer,      Intent (In) :: k
    Real (KREAL), Intent (In) :: z
    
    Real (KREAL)              :: Du
    
    ! Node 1 is -1, node 2 is +1
    select case (k)
    case (1); Du = - HALF
    case (2); Du = HALF
    end select
    
End Function

Function TetrahedralSF (k, r) Result (u)
    Implicit None
    
    ! INCLUDE SHAPE LABELS WITH IMAGE LIKE GMSH HAS
    
    Integer,      Intent (In) :: k     ! Node index
    Real (KREAL), Intent (In) :: r (3) ! Point
    
    Real (KREAL)              :: u
    
    select case (k)
    case (1); u = ONE - SUM(r)
    case (2); u = r(1)
    case (3); u = r(2)
    case (4); u = r(3)
    end select
    
End Function

Function GradTetrahedralSF (k, r) Result (Du)
    Implicit None
    Integer,      Intent (In) :: k      ! Node index
    Real (KREAL), Intent (In) :: r  (3) ! Point
    
    Real (KREAL)              :: Du (3)
    
    select case (k)
    case (1); Du = [-ONE,  -ONE,  -ONE  ]
    case (2); Du = [ ONE,   ZERO,  ZERO ]
    case (3); Du = [ ZERO,  ONE,   ZERO ]
    case (4); Du = [ ZERO,  ZERO,  ONE  ]
    end select
    
End Function

Function HexahedralSF (k, r) Result (u)
    Implicit None
    Integer,      Intent (In) :: k
    Real (KREAL), Intent (In) :: r (3)
    
    Real (KREAL)              :: u
    ! MUST REWRITE THE CONVENTIONS HERE
    print *, 'wip3409235'
    stop
    select case (k)
    case (1); u = (1 - r(1)) * (1 - r(2)) * (1 - r(3)) / EIGHT
    case (2); u = (1 - r(1)) * (1 + r(2)) * (1 - r(3)) / EIGHT
    case (3); u = (1 - r(1)) * (1 + r(2)) * (1 + r(3)) / EIGHT
    case (4); u = (1 - r(1)) * (1 - r(2)) * (1 + r(3)) / EIGHT
    case (5); u = (1 + r(1)) * (1 - r(2)) * (1 - r(3)) / EIGHT
    case (6); u = (1 + r(1)) * (1 + r(2)) * (1 - r(3)) / EIGHT
    case (7); u = (1 + r(1)) * (1 + r(2)) * (1 + r(3)) / EIGHT
    case (8); u = (1 + r(1)) * (1 - r(2)) * (1 + r(3)) / EIGHT
    end select
    
End Function

Function GradHexahedralSF (k, r) Result (Du)
    Implicit None
    Integer,      Intent (In) :: k
    Real (KREAL), Intent (In) :: r  (3)
    
    Real (KREAL)              :: Du (3)
    ! MUST REWRITE THE CONVENTIONS HERE
    print *, 'wip34092355'
    stop
    select case (k)
    case (1); Du = [ - (1 - r(2)) * (1 - r(3)) / EIGHT, &
                     - (1 - r(1)) * (1 - r(3)) / EIGHT, &
                     - (1 - r(1)) * (1 - r(2)) / EIGHT    ]
    case (2); Du = [ - (1 + r(2)) * (1 - r(3)) / EIGHT, &
                       (1 - r(1)) * (1 - r(3)) / EIGHT, &
                     - (1 - r(1)) * (1 + r(2)) / EIGHT    ]
    case (3); Du = [ - (1 + r(2)) * (1 + r(3)) / EIGHT, &
                       (1 - r(1)) * (1 + r(3)) / EIGHT, &
                       (1 - r(1)) * (1 + r(2)) / EIGHT    ]
    case (4); Du = [ - (1 - r(2)) * (1 + r(3)) / EIGHT, &
                     - (1 - r(1)) * (1 + r(3)) / EIGHT, &
                       (1 - r(1)) * (1 - r(2)) / EIGHT    ]
    case (5); Du = [   (1 - r(2)) * (1 - r(3)) / EIGHT, &
                     - (1 + r(1)) * (1 - r(3)) / EIGHT, &
                     - (1 + r(1)) * (1 - r(2)) / EIGHT    ]
    case (6); Du = [   (1 + r(2)) * (1 - r(3)) / EIGHT, &
                       (1 + r(1)) * (1 - r(3)) / EIGHT, &
                     - (1 + r(1)) * (1 + r(2)) / EIGHT    ]
    case (7); Du = [   (1 + r(2)) * (1 + r(3)) / EIGHT, &
                       (1 + r(1)) * (1 + r(3)) / EIGHT, &
                       (1 + r(1)) * (1 + r(2)) / EIGHT    ]
    case (8); Du = [   (1 - r(2)) * (1 + r(3)) / EIGHT, &
                     - (1 + r(1)) * (1 + r(3)) / EIGHT, &
                       (1 + r(1)) * (1 - r(2)) / EIGHT    ]
    end select
    
End Function

Function TriangularSF (k, r) Result (u)
    Implicit None
    Integer,      Intent (In) :: k     ! Node index
    Real (KREAL), Intent (In) :: r (2) ! Point
    
    Real (KREAL)              :: u
    
    select case (k)
    case (1); u = ONE - SUM(r)
    case (2); u = r(1)
    case (3); u = r(2)
    end select
    
End Function

Function GradTriangularSF (k, r) Result (Du)
    Implicit None
    Integer,      Intent (In) :: k      ! Node index
    Real (KREAL), Intent (In) :: r  (2) ! Point
    
    Real (KREAL)              :: Du (2)
    
    select case (k)
    case (1); Du = [-ONE, -ONE]
    case (2); Du = [ ONE, ZERO]
    case (3); Du = [ZERO, ONE ]
    end select
    
End Function

Function RectangularSF (k, r) Result (u)
    Implicit None
    Integer,      Intent (In) :: k     ! Node index
    Real (KREAL), Intent (In) :: r (2) ! Point
    
    Real (KREAL)              :: u
    ! MUST REWRITE THE CONVENTIONS HERE
    select case (k)
    case (1); u = (1 - r(1)) * (1 - r(2)) / FOUR
    case (2); u = (1 - r(1)) * (1 + r(2)) / FOUR
    case (3); u = (1 + r(1)) * (1 - r(2)) / FOUR
    case (4); u = (1 + r(1)) * (1 + r(2)) / FOUR
    end select
    
End Function

Function GradRectangularSF (k, r) Result (Du)
    Implicit None
    Integer,      Intent (In) :: k      ! Node index
    Real (KREAL), Intent (In) :: r  (2) ! Point
    
    Real (KREAL)              :: Du (2)
    ! MUST REWRITE THE CONVENTIONS HERE
    select case (k)
    case (1); Du = [- (1 - r(2)) / FOUR, - (1 - r(1)) / FOUR]
    case (2); Du = [- (1 + r(2)) / FOUR,   (1 - r(1)) / FOUR]
    case (3); Du = [  (1 - r(2)) / FOUR, - (1 + r(1)) / FOUR]
    case (4); Du = [  (1 + r(2)) / FOUR,   (1 + r(1)) / FOUR]
    end select
    
End Function

End Module ShapeFunctions