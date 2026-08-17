Module Jacobians
    use constants
    use BasicMathFunctions
    use ShapeFunctions
    Implicit None
    
    !  =======================================================
    !    Jacobians for transformations to reference elements
    !  =======================================================
    
Contains

Subroutine Jacobian (idx, nodes, rg, Jv, JvinvT, detJv)
    Implicit None
    Integer,      Intent (In) :: idx           ! Element type index
    Integer,      Intent (In) :: nodes   (:)   ! Node addresses
    Real (KREAL), Intent (In) :: rg      (:,:) ! Global mesh
    
    Real (KREAL), Intent (Out) :: Jv     (3,3)
    Real (KREAL), Intent (Out) :: JvinvT (3,3)
    Real (KREAL), Intent (Out) :: detJv
    
    Integer                    :: j, k, k1, kg
    Real (KREAL)               :: rp  (3)
    Real (KREAL)               :: Du  (3)
    
    rp = ZERO
    
    select case (idx)
    case (1)
        ! For a line element, only the (1,1) component of Jv is meaningful.
        
        Jv(1,1) = HALF * (rg(1,nodes(2)) - rg(1,nodes(1)))
        JvinvT  = ONE / Jv(1,1)
        detJv   = Jv(1,1)
        
        return
    case (4)
        ! Formula for Jacobian is:
        ! J_{ij} = r(i,j+1) - r(i,1)
        ! for local coordinates r(1:3,1:NK)
        k1 = nodes(1)
        do j = 1, 3
            k = nodes(j + 1)
            
            Jv(1:3,j) = rg(1:3,k) - rg(1:3,k1)
        end do
        
    case (5)
        ! More general formula, involving gradient of shape functions
        do k = 1, SIZE(nodes)
            kg = nodes(k)
            Du = GradHexahedralSF (kg, rp)
            do j = 1, 3
                Jv(1:3,j) = Jv(1:3,j) + rg(1:3,kg) * Du(j)
            end do
        end do
        
    case (6)
        ! WIP
    case (7)
        ! WIP
    end select
    
    !  -------------------------------------
    !    Invert the Jacobian and transpose
    !  -------------------------------------
    ! Crude but for well-conditioned Jacobians (which come from evenly shaped elements)
    ! this should be sufficiently numerically stable.
    ! Will certainly consider adding a LU inversion as a fallback
    
    detJv = Jv(1,1) * (Jv(2,2) * Jv(3,3) - Jv(2,3) * Jv(3,2)) &
          + Jv(1,2) * (Jv(2,3) * Jv(3,1) - Jv(2,1) * Jv(3,3)) &
          + Jv(1,3) * (Jv(2,1) * Jv(3,2) - Jv(2,2) * Jv(3,1))
    
    JvinvT(1,1) = Jv(2,2) * Jv(3,3) - Jv(2,3) * Jv(3,2)
    JvinvT(1,2) = Jv(1,3) * Jv(3,2) - Jv(1,2) * Jv(3,3)
    JvinvT(1,3) = Jv(1,2) * Jv(2,3) - Jv(1,3) * Jv(2,2)
    
    JvinvT(2,1) = Jv(2,3) * Jv(3,1) - Jv(2,1) * Jv(3,3)
    JvinvT(2,2) = Jv(1,1) * Jv(3,3) - Jv(1,3) * Jv(3,1)
    JvinvT(2,3) = Jv(1,3) * Jv(2,1) - Jv(1,1) * Jv(2,3)
    
    JvinvT(3,1) = Jv(2,1) * Jv(3,2) - Jv(2,2) * Jv(3,1)
    JvinvT(3,2) = Jv(1,2) * Jv(3,1) - Jv(1,1) * Jv(3,2)
    JvinvT(3,3) = Jv(1,1) * Jv(2,2) - Jv(1,2) * Jv(2,1)
    
    JvinvT = TRANSPOSE(JvinvT) / detJv
    
    ! ! Check that this Jacobian is sufficiently well-conditioned.
    ! ! If it's not, use LAPACK routine to get the proper Jacobian.
    ! if (MAXVAL(ABS(Jv)) * MAXVAL(ABS(JinvT)) / ABS(detJv) > TOL) then
    !     ! Do an LU inversion for numerical stability
    !     ! Maybe not - if an element is ill-conditioned LU inversion likely won't save it. User should just provide a better mesh.
    ! end if
    
End Subroutine

Subroutine SurfaceJacobian (idxf, nodes, rg, Js, normJs)
    Implicit None
    Integer,      Intent (In) :: idxf          ! Face element type index
    Integer,      Intent (In) :: nodes   (:)   ! Node addresses
    Real (KREAL), Intent (In) :: rg      (:,:) ! Global mesh
    
    Real (KREAL), Intent (Out) :: Js     (3)
    Real (KREAL), Intent (Out) :: normJs
    
    Integer                    :: k, kg
    Integer                    :: NKf
    Real (KREAL)               :: rp (2)
    Real (KREAL)               :: Du (2)
    Real (KREAL)               :: t1 (3)
    Real (KREAL)               :: t2 (3)
    !Real (KREAL)               :: rTEST (3,3)
    
    !rTEST(:,1) = [ZERO, ZERO, ZERO]
    !rTEST(:,2) = [ONE, ZERO, ZERO]
    !rTEST(:,3) = [TWO, ZERO, ZERO]
    
    rp  = ZERO
    t1  = ZERO
    t2  = ZERO
    NKf = SIZE(nodes)
    
    do k = 1, NKf
        select case (idxf)
        case (0)
            print *, 'wip13243 not sure if Ill ever need this'
            stop
            !Du = 
        case (1)
            Du = GradTriangularSF  (k, rp)
        case (2)
            Du = GradRectangularSF (k, rp)
        end select
        
        kg = nodes(k)
        
        t1 = t1 + rg(1:3,kg) * Du(1)
        t2 = t2 + rg(1:3,kg) * Du(2)
        
        ! t1 = t1 + rTEST(:,k) * Du(1)
        ! t2 = t2 + rTEST(:,k) * Du(2)
        
    end do
    
    Js = CROSS_PRODUCT(t1, t2)
    
    normJs = norm2(Js)
    
End Subroutine

End Module Jacobians