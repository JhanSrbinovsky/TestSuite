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


!      Subroutine ICWriteEntry(value)
!      ------------------------------
!      This routine is passed the number of iterations in its
!      argument and it then writes this to the output file if
!      it is open.

      subroutine ICwriteEntry(value)

        implicit none
#include "cntlatm.h"
!      Arguments
        integer, intent(in) :: value ! The number of iterations

!      If the unit is open the write the value to the file stream
        if(l_icount) then
          write(152,*) value
        end if
      end subroutine ICwriteEntry
#endif
