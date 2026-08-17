Module globals
    Implicit None
    
    ! Globally-used objects
    Logical :: gstorlim    = .FALSE. ! Whether or not storage is a concern. TRUE = angular fluences are NOT saved, but fluences still are.
    Logical :: gmemlim     = .FALSE. ! Whether or not memory is a concern. TRUE = solves are temporarily saved to disk and then read in during PI/EI
    Logical :: gdebug      = .FALSE.
    
End Module globals