#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine READSIZE --------------------------------------------
!LL
!LL    Purpose:
!LL  This routine reads the length required for the dimensioning
!LL  of all main data arrays that have dimensions set via questions
!LL  in the USER INTERFACE.
!LL  It is called by the shell program UM_SHELL and passes the
!LL  data via common block to that routine. Data is then passed by
!LL  argument lists to allow for proper dynamic allocation in the
!LL  portable model.
!LL
!LL Programming standard:
!LL
!LL Logical components covered:
!LL
!LL Project task:
!LL
!LL External documentation:
!LL
!LL ------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE READSIZE()

      IMPLICIT NONE

!  Local variables

      INTEGER                       :: ERRORSTATUS
      CHARACTER (LEN=*), PARAMETER  :: ROUTINENAME='READSIZE'
      CHARACTER (LEN=80)            :: CMESSAGE

#include "parparm.h"
#include "cnlsizes.h"
#include "amaxsize.h"

!L 1.1 Read size information from namelists

!     Initialise TPPS_OZONE_LEVELS from OZONE_LEVELS
!     Prevents failures when using 'preset=nan'
!     Currently not set through UMUI
      TPPS_OZONE_LEVELS = 0
! Similarly, JO_NMAX_OBS_ICE and JO_MAX_OBS_VAL
! which are only conditionally set in UMUI
!     Currently not set through UMUI
!     Initialise river routing sizes
      RIVER_ROW_LENGTH  = 0
      RIVER_ROWS        = 0

!     FILE ON UNIT 5 DOES NOT NEED TO BE OPENED.
!     READ THE UI PROVIDED NAMELIST OF ARRAY SIZES.
      READ(5,NLSIZES)

!      Copy OCEAN values from NAMELIST to COMMON

! Check that dynamically provided sizes are not larger than their
! statically defined maximums.
      IF (GLOBAL_ROW_LENGTH > ROW_LENGTH_MAX) THEN
        ERRORSTATUS = 10
        CMESSAGE = 'Row length is larger than maximum defined in'       &
     &           //' AMAXSIZE'
! DEPENDS ON: ereport
        CALL EREPORT( ROUTINENAME, ERRORSTATUS, CMESSAGE )
      END IF

      IF (GLOBAL_ROWS > ROWS_MAX) THEN
        ERRORSTATUS = 20
        CMESSAGE = 'Number of rows is larger than maximum defined in'   &
     &           //' AMAXSIZE'
! DEPENDS ON: ereport
        CALL EREPORT( ROUTINENAME, ERRORSTATUS, CMESSAGE )
      END IF

      IF (MODEL_LEVELS > MODEL_LEVELS_MAX) THEN
        ERRORSTATUS = 30
        CMESSAGE = 'Number of levels is larger than maximum defined in' &
     &           //' AMAXSIZE'
! DEPENDS ON: ereport
        CALL EREPORT( ROUTINENAME, ERRORSTATUS, CMESSAGE )
      END IF

!  -------------------------------------------------------------------
!  Namelist STSIZES used to be read here, from the former STASHC file.
!  The STSIZES common block is now set within STASH_PROC.
!  -------------------------------------------------------------------
! Call to DERVSIZE is moved into UM_SHELL

      RETURN
      END SUBROUTINE READSIZE
#endif
