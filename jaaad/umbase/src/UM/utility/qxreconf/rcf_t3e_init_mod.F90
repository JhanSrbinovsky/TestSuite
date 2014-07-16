#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ T3E specific initialisation

Module Rcf_T3E_Init_Mod

!  Subroutine Rcf_T3E_Init - T3E initialisation tasks
!
! Description:
! This module does T3E specific initialisation - specifically
! specifying resolution for the timer and splitting the output
! streams - 1 per pe
!
! Method:
!   Similar to code in UM_SHELL from the UM.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   6.0   07/01/04   (re)Move section to send unit 6 output to a unique
!                    file for each PE. R.Sharp
!   6.2   11/08/05   Replace quotes in includes with angle brackets.
!                    P.Selwood.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_T3E_Init()

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Parvars_Mod, Only : &
    nproc,              &
    mype

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

Implicit None

#if defined(T3E)
! Local Variables
Integer                        :: iclktck
Integer                        :: ierr
Integer                        :: ErrorStatus
Integer                        :: len_dataw       ! length of $DATAW
Integer                        :: len_runid       ! length of $RUNID
Integer                        :: len_basename    ! length of file base
Integer                        :: get_char_len    ! function
Integer                        :: um_nam_max_seconds
                                  ! timeout value for shmem nam
Integer                        :: um_rnl_skip     ! namelist skipping
Character (Len=*), Parameter   :: RoutineName = 'Rcf_T3E_Init'
Character (Len=80)             :: Cmessage
Character (Len=200)            :: stdout_filename ! file for stdout
Character (Len=180)            :: stdout_basename ! base of filename
Character (Len=170)            :: dataw_char      ! value of $DATAW
Character (Len=5)              :: runid_char      ! value of $RUNID
Character (Len=8)              :: c_nam_max_seconds
                                          ! value of $UM_NAM_MAX_SECONDS
Character (Len=8)              :: c_um_rnl_skip
                                          ! value of $UM_RNL_SKIP

! Comdecks
#include "t3eclktk.h"
#include "gccom.h"

External PxfConst, Pxfsysconf, Fort_Get_Env, Get_Char_Len

!---------------------------------------------------------------------
! Find out the number of clock ticks per second on this machine.
! This information is required to calculate the wallclock times
! in TIMER
!---------------------------------------------------------------------

iclktck=0
ticks_per_second=0

CALL PXFCONST('CLK_TCK',iclktck,ierr)

IF (ierr .NE. 0) THEN
  Write (Cmessage,*)  'Failure in PXFCONST, err= ',ierr
  ErrorStatus = 10
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
ENDIF

CALL PXFSYSCONF(iclktck,ticks_per_second,ierr)

IF (ierr .NE. 0) THEN
  Write (Cmessage, *) 'Failure in PXFSYSCONF, err= ',ierr
  ErrorStatus = 20
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
ENDIF

!---------------------------------------------------------------------
!  See if the SHMEM/NAM Timeout value has been set in a shell
!  variable.
!---------------------------------------------------------------------
Call fort_get_env('UM_NAM_MAX_SECONDS', 18, c_nam_max_seconds,  8, ierr)
If (ierr .ne. 0) Then
  um_nam_max_seconds=300.
  If (mype == 0) then
    write(6,*) 'Warning: Environment variable UM_NAM_MAX_SECONDS ',  &
               'has not been set.'
    write(6,*) 'Setting UM_NAM_MAX_SECONDS to ', um_nam_max_seconds
    write(0,*) 'Warning: Environment variable UM_NAM_MAX_SECONDS ',  &
               'has not been set.'
    write(0,*) 'Setting UM_NAM_MAX_SECONDS to ', um_nam_max_seconds
  End If
Else
  Read(c_nam_max_seconds,'(i8)') um_nam_max_seconds
  If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
    write(6,*) 'Setting UM_NAM_MAX_SECONDS to ',  um_nam_max_seconds
    write(0,*) 'Setting UM_NAM_MAX_SECONDS to ',  um_nam_max_seconds
  Endif
Endif

!--now set the value
Call gc_setopt(GC_NAM_TIMEOUT, um_nam_max_seconds, ierr)
If (ierr /= GC_OK) Then
  Write (Cmessage,*) 'Response from GC_SETOPT was ', ierr
  ErrorStatus = 50
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

!---------------------------------------------------------------------
! Do we wish to see namelist skip messages?
!---------------------------------------------------------------------
Call Fort_Get_Env( 'UM_RNL_SKIP', 11, c_um_rnl_skip, 8, ierr )

If ( ierr /= 0 ) Then
  Write(6,*) 'Warning: Environment variable UM_RNL_SKIP has ', &
             'not been set.'
  Write(6,*) 'Setting um_rnl_skip to 0 - Omit Skipped Messages'
  um_rnl_skip=0
Else
  Read(c_um_rnl_skip,'(I4)') um_rnl_skip
End If

Call rnlskip(um_rnl_skip)

#endif
Return
End Subroutine Rcf_T3E_Init
End Module Rcf_T3E_Init_Mod
#endif
