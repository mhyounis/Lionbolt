Module SolversInterface
    Implicit None
    
    Interface
        
        Module Subroutine SourceIteration (OpT, OpK, s)
            use constants
            use profiler
            use LinearOperatorClass
            Implicit None
            Class (LinearOperator),  Intent (InOut) :: OpT
            Class (LinearOperator),  Intent (InOut) :: OpK
            Type (SpaceAngleVector), Intent (InOut) :: s
        End Subroutine
        
        ! Should eventually have more here
        
    End Interface
    
End Module SolversInterface