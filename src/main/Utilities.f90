Module Utilities
    use LionboltAPI, ONLY : LBpath
    use constants
    use globals
    use IO
    Implicit None
    
    ! NEW PLAN FOR THOSE USING gmemlim = .FALSE.
    ! Just use a huge real array and flatten it. Store offsets. Is this bad/slow?
    
    ! Generic interfaces for tape writing and reading
    Interface TapeWriter
        Module Procedure TWR1 ! Real, rank one
        Module Procedure TWR2 ! Real, rank two
    End Interface
    
    Interface TapeReader
        Module Procedure TRR1 ! Real, rank one
        Module Procedure TRR2 ! Real, rank two
    End Interface
    
    Private :: TWR1, TWR2, TRR1, TRR2
    
Contains
    
    Function TapeName (tape, tag) Result (fname)
        Implicit None
        Integer,           Intent (In) :: tape
        Integer,           Intent (In) :: tag
        
        Character (LEN=:), Allocatable :: fname
        
        Character (STDMAXLEN)          :: tmp ! 512 is definitely safe for tapes
        
        ! Tape to type of data:
        ! 0 : Geometry-related data generated my mesh readers
        ! p : Fully solved fluence of particle p
        
        write(tmp, '(I0, "%", I0, ".dat")') tape, tag
        fname = ScratchName // tmp
        
    End Function
    
    Subroutine TWR1 (tape, tag, obj)
        Implicit None
        Integer,           Intent (In) :: tape
        Integer,           Intent (In) :: tag
        Real (KREAL),      Intent (In) :: obj (:)
        
        Character (LEN=:), Allocatable :: fname
        
        Integer                        :: i
        
        fname = TapeName (tape, tag)
        
        open(unit=TAPEUNIT, file=fname, form='unformatted', access='sequential', status='replace')
        write(TAPEUNIT) obj
        close(TAPEUNIT)
        
        call RecordTapeName (fname)
        
    End Subroutine
    
    Subroutine TWR2 (tape, tag, obj)
        Implicit None
        Integer,           Intent (In) :: tape
        Integer,           Intent (In) :: tag
        Real (KREAL),      Intent (In) :: obj (:,:)
        
        Character (LEN=:), Allocatable :: fname
        
        Integer                        :: i
        
        fname = TapeName (tape, tag)
        
        open(unit=TAPEUNIT, file=fname, form='unformatted', access='sequential', status='replace')
        write(TAPEUNIT) obj
        close(TAPEUNIT)
        
        call RecordTapeName (fname)
        
    End Subroutine
    
    Subroutine RecordTapeName (fname)
        Implicit None
        Character (*),    Intent (In) :: fname
        
        Type (WordsType), Allocatable :: tmp (:)
        
        if (ALLOCATED(activetapes)) then
            ALLOCATE(tmp(SIZE(activetapes) + 1))
            tmp(1:SIZE(activetapes))       = activetapes
            tmp(SIZE(activetapes) + 1)%str = fname
            call move_alloc (tmp, activetapes)
        else
            ALLOCATE(activetapes(1))
            activetapes(1)%str = fname
        end if
        
    End Subroutine
    
    Subroutine TRR1 (tape, tag, obj)
        Implicit None
        Integer,        Intent (In)    :: tape
        Integer,        Intent (In)    :: tag
        
        Real (KREAL),   Intent (InOut) :: obj (:)
        
        Logical                        :: exists
        Character (LEN=:), Allocatable :: fname
        
        fname = TapeName (tape, tag)
        
        inquire(file=fname, exist=exists)
        
        if (.not. exists) then
            call stophere ("TapeReader : Attempting to read a tape that doesn't exist : " // fname)
        end if
        
        open(unit=TAPEUNIT, file=fname, form='unformatted', access='sequential')
        read(TAPEUNIT) obj
        close(TAPEUNIT)
        
    End Subroutine
    
    Subroutine TRR2 (tape, tag, obj)
        Implicit None
        Integer,        Intent (In)    :: tape
        Integer,        Intent (In)    :: tag
        
        Real (KREAL),   Intent (InOut) :: obj (:,:)
        
        Logical                        :: exists
        Character (LEN=:), Allocatable :: fname
        
        fname = TapeName (tape, tag)
        
        inquire(file=TRIM(ADJUSTL(fname)), exist=exists)
        
        if (.not. exists) then
            call stophere ("TapeReader : Attempting to read a tape that doesn't exist : " // fname)
        end if
        
        open(unit=TAPEUNIT, file=TRIM(ADJUSTL(fname)), form='unformatted', access='sequential')
        read(TAPEUNIT) obj
        close(TAPEUNIT)
        
    End Subroutine
    
End Module