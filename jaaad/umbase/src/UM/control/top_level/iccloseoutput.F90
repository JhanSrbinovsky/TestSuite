#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! History
! -------
! Version    Date                 Comment
! -------    ----                 -------
!  6.2     11/01/06    Original code  T.Edwards

!       The deck, icdiag, routines that are used to manage a
!       fortran unit file stream which records a list of
!       teration counts.


!         Subroutine ICOpenOutput
!         -----------------------
!       Opens the file stream constructing a filename from
!       the the ICFILNAME variables that has been set in a
!       namelist. If this is blank or just spaces, then
!       construct a name from the $DATAW and $RUNID
!       environment variables.


!       Subroutine ICCloseOutput
!       ------------------------
!       Closes the output stream if it is open

      subroutine ICCloseOutput()

          implicit none
#include "cntlatm.h"

!       Local Variables
        integer :: icode ! Error return code

        if(l_icount) then
          close(unit=152,iostat=icode)
          l_icount = .false.
          if(icode > 0) then
            write(6,*) "ITERCOUNT : Failed to close output file"
          end if
        end if
      end subroutine ICCloseOutput

!      Subroutine ICWriteEntry(value)
!      ------------------------------
!      This routine is passed the number of iterations in its
!      argument and it then writes this to the output file if
!      it is open.

#endif
