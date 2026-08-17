Submodule (Sources) submodRayTrace
Contains
Module Subroutine RayTrace (self)
    ! I want to keep abstract RayTrace routines purely internal, while
    ! this subroutine is type-bound to the ExternalBeam type. Thus, this
    ! subroutine is essentially just a wrapper for a routine in RayTracing
    use RayTracing
    use profiler
    Implicit None
    Class (ExternalBeam), Intent (InOut) :: self
    ! call profile ('RayTrace')
    if (ALLOCATED(self%rays)) then
        call stophere ('RayTrace.f90: RayTrace: RayTrace was called, but this ExternalBeam has already done ray tracing')
    end if
    call RayTraceBeam (self)
    ! call profile ('RayTrace')
End Subroutine
End Submodule