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

      subroutine ICOpenOutput(runtype)

          implicit none

#include "chsunits.h"
#include "clfhist.h"
#include "cntlatm.h"
#include "parvars.h"

!       Arguments
      character(len=4)  :: runtype

!       Local variables

      integer           :: icode       ! Error return strings
      character(len=64) :: datawString ! To hold DATAW env var
      character(len=5)  :: runidString ! To hold RUNID env var
      character(len=80) :: fileName    ! Name of filename
      logical           :: fileExists  ! Check if file already exists

!     Checks if output is required, and if unit is already open
      if(l_icount .and. (mype == 0)) then
        fileName = adjustl(icfile(10:))
        if(len_trim(fileName) == 0) then
          call fort_get_env('DATAW',5,datawString,64,icode)
          if (icode /= 0) then
            write(6,*) 'ITERCOUNT : Failed to get value of $DATAW'
          end if
          call fort_get_env('RUNID',5,runidString,5,icode)
          if (icode /= 0) then
            write(6,*) 'ITERCOUNT : Failed to get value of $RUNID'
          end if
          fileName = trim(datawString) // "/" // trim(runidString) //   &
     &        ".itercount"
        end if

        if(runtype /= "CRUN") then
          open(unit=152,file=fileName,status='replace', iostat=icode)
        else
!         Check that the file exists if this is a CRUN
          inquire(file=fileName,exist=fileExists)
          if(fileExists) then
            open(unit=152,file=fileName,status='old', iostat=icode,     &
     &          position="append")
!         Otherwise complain loudly but overwrite
          else
            write(6,*) "ITERCOUNT : CRUN but can't find old output " // &
     &                  "file! Starting from scratch"
            open(unit=152,file=fileName,status='new', iostat=icode)
          end if
        end if
        if(icode > 0) then
          write(6,*) "ITERCOUNT : Failed to open ", fileName
          l_icount = .false.
        end if
      else
!       If we are not PE0 then ensure l_icount is set to false
        l_icount = .false.
      end if

      end subroutine ICOpenOutput

!       Subroutine ICCloseOutput
!       ------------------------
!       Closes the output stream if it is open


!      Subroutine ICWriteEntry(value)
!      ------------------------------
!      This routine is passed the number of iterations in its
!      argument and it then writes this to the output file if
!      it is open.

#endif
