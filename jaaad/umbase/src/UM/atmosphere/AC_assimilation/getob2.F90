#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE GETOBS,GETOB2,GETOB3----------------------------------
!LL
!LL  Purpose : Extract appropriate observations from the COMOBS common
!LL            block and the OBS array (passed from RDOBS via argument
!LL            list)
!LL            GETOBS called from AC gets list of obs in time window
!LL            GETOB2 called from AC2 gets lat, long, time and
!LL                   model obs type for an obs type.
!LL            GETOB3 called from VERTANL gets data specific to a type
!LL                   (eg data values and assoc error ratios)
!LL
!LL  For Cray - Global  ; Enable defs GLOBAL
!LL
!LL  Transfer from Cyber by Dave Robinson.
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL  3.1     18/11/92  : Option to thin GLOSS      S Bell
!LL  3.2     14/6/93   : use OBTHIN in place of THINxxx; cater for 306;
!LL                    : Eliminate QA FORTRAN complaints    S Bell
!LL  3.3     2/12/93   : New ob type 506/407   Bruce M/ Nigel R
!LL  3.4     9/09/94   : New ob type 901 Pete Clark
!LL  4.0     11/8/95   : deleting redundant comments relating to
!LL                    : OBS_FORMAT=1  SB
!LL  4.1   31/05/96     The number of v points to be processed on a
!LL                     C grid differs from u by row_length. u,v
!LL                     dimensioned separately in call to WLLTOEQ.
!LL                     Requirement for VAR.
!LL                     Author I.Edmond       Reviewer D. Goddard
!LL  4.2   25/11/96     Mods for T3E + adjust time window Stuart Bell
!LL  5.2   18/12/00     remove refs to daco-related diagnostic routines
!LL                                                B Macpherson
!LL  6.0   10/10/03     Replace SHMEM with GCOM for SX6. Clive JOnes
!LL  6.2   21/10/05     Replace GSYNC with SSYNC. P.Selwood.
!LL
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
!
      SUBROUTINE GETOB2 (KACT,OBS,INOBS,OBS_LAT,OBS_LONG,OBS_TIME,      &
     &                   MODEL_OBS_TYPE,OBS_NO,LENOBT,                  &
     &                   ICODE,CMESSAGE)
!
!
      IMPLICIT NONE
!     UM AC comdecks
#include "acparm.h"
#include "comobs.h"
#include "comacp.h"
#include "comacdg.h"
#if !defined(GLOBAL)
#include "commg.h"
#endif
!     UM comdeck
#include "c_pi.h"
!-----------------------------------------------------------------------
      INTEGER KACT                     ! IN  Obs type index
      INTEGER INOBS                    ! IN  No of obs for this type
      INTEGER LENOBT                   ! IN  No of chosen obs
      INTEGER OBS_NO (LENOBT)          ! IN  pointers to chosen obs
      INTEGER MODEL_OBS_TYPE (LENOBT)  ! OUT model obs types

      REAL    OBS(INOBS,*)             ! IN  obs from RDOBS
      REAL    OBS_LAT (LENOBT)         ! OUT latitudes
      REAL    OBS_LONG(LENOBT)         ! OUT longitudes
      REAL    OBS_TIME(LENOBT)         ! OUT times

      INTEGER ICODE                    ! OUT error code
      CHARACTER*256 CMESSAGE           ! OUT error message
!-----------------------------------------------------------------------
!L    Local work arrays
      REAL WLT (LENOBT)
      REAL WLN (LENOBT)
!-----------------------------------------------------------------------
!     Local variables
      INTEGER JOB       !  Loop conter over obs
      INTEGER IP_LAT    !  Pointer to obs latitude
      INTEGER IP_LONG   !  Pointer to obs longitude
      INTEGER IP_TIME   !  Pointer to obs times
      INTEGER IP_TYPE   !  Pointer to model obs type numbers
!-----------------------------------------------------------------------
      EXTERNAL TIMER
#if !defined(GLOBAL)
      EXTERNAL EQTOLL
#endif
!-----------------------------------------------------------------------
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('GETOB2  ',3)

!     Set up pointers to data in OBS array
      IP_LAT  = 1   !  Latitude
      IP_LONG = 2   !  Longitude
      IP_TIME = 3   !  Time
      IP_TYPE = 4   !  Model observation type (**)

!     Get latitudes, longitudes, time and model obs type numbers
      DO JOB=1,LENOBT
       OBS_LAT(JOB)        = OBS(OBS_NO(JOB) - OBS_NO_ST(KACT),IP_LAT)
       OBS_LONG(JOB)       = OBS(OBS_NO(JOB) - OBS_NO_ST(KACT),IP_LONG)
       OBS_TIME(JOB)       = OBS(OBS_NO(JOB) - OBS_NO_ST(KACT),IP_TIME)
       MODEL_OBS_TYPE(JOB) =                                            &
     &                NINT(OBS(OBS_NO(JOB) - OBS_NO_ST(KACT),IP_TYPE))
      ENDDO


! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('GETOB2  ',4)
      RETURN
      END SUBROUTINE GETOB2
#endif
