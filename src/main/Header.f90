Subroutine Header ()
    use IO
    Implicit None
    
    ! Version
    Character (5) :: v = '1.0.0'
    
    ! Announcements
    write (iuout,*) '#############################################################################################'
    write (iuout,*) '   Lionbolt --- Developed by Muhsin H. Younis (myounis@psu.edu / myounis@umd.edu)            '
    write (iuout,*) '   VERSION: ',v,' (Released 8-17-26)                                                         '
    write (iuout,*) '      Licensed under the GNU General Public License 3.0.                                     '
    write (iuout,*) '      See LICENSE/gpl-3.0.txt or visit https://www.gnu.org/licenses/ for more information.   '
    write (iuout,*) '#############################################################################################'
    
    ! State the input and output file locations
    call GetCurrentWorkingDirectory ()
    write (iuout,*) ''
    write (iuout,*) 'Input File  : ', cwd // '/' // InputFname
    write (iuout,*) 'Output File : ', cwd // '/' // OutputFname
    write (iuout,*) 'Scratch     : ', ScratchName
    write (iuout,*) ''
    
    ! PRINT THE HEADER AS WELL
    print *, ''
    print *, '#############################################################################################'
    print *, '   Lionbolt --- Developed by Muhsin H. Younis (myounis@psu.edu / myounis@umd.edu)            '
    print *, '   VERSION: ',v,' (Released 8-17-26)                                                         '
    print *, '      Licensed under the GNU General Public License 3.0.                                     '
    print *, '      See LICENSE/gpl-3.0.txt or visit https://www.gnu.org/licenses/ for more information.   '
    print *, '#############################################################################################'
    
    print *, ''
    print *, 'Input File  : ', cwd // '/' // InputFname
    print *, 'Output File : ', cwd // '/' // OutputFname
    print *, 'Scratch     : ', ScratchName
    print *, ''
    
End Subroutine