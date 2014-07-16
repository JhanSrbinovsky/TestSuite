#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set pointers to types of clouds
!
! Method:
!       The types of condensate included are examined. Their phases
!       are set and depending on the representation of clouds adopted
!       it is determined to which type of cloud they contribute.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SET_CLOUD_POINTER(IERR                                 &
     &   , N_CONDENSED, TYPE_CONDENSED, I_CLOUD_REPRESENTATION          &
     &   , L_DROP, L_ICE                                                &
     &   , I_PHASE_CMP, I_CLOUD_TYPE, L_CLOUD_CMP                       &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!
!     INCLUDE COMDECKS
#include "stdio3a.h"
#include "error3a.h"
#include "dimfix3a.h"
#include "phase3a.h"
#include "cldcmp3a.h"
#include "clrepp3a.h"
#include "cldtyp3a.h"
!
!     DUMMY VARIABLES.
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_CONDENSED                                                  &
!             NUMBER OF CONDENSED COMPONENTS
     &   , TYPE_CONDENSED(NPD_CLOUD_COMPONENT)                          &
!             TYPES OF COMPONENTS
     &   , I_CLOUD_REPRESENTATION
!             REPRESENTATION OF CLOUDS USED
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_DROP                                                       &
!             FLAG FOR INCLUSION OF DROPLETS
     &   , L_ICE
!             FLAG FOR INCLUSION OF ICE CRYSTALS
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     I_PHASE_CMP(NPD_CLOUD_COMPONENT)                             &
!             PHASES OF COMPONENTS
     &   , I_CLOUD_TYPE(NPD_CLOUD_COMPONENT)
!             TYPES OF CLOUD TO WHICH EACH COMPONENT CONTRIBUTES
      LOGICAL                                                           &
                !, INTENT(OUT)
     &     L_CLOUD_CMP(NPD_CLOUD_COMPONENT)
!             LOGICAL SWITCHES TO INCLUDE COMPONENTS
!
!
!     LOCAL VARIABLES
      INTEGER                                                           &
     &     K
!            LOOP VARIABLE
!
#include "clrepd3a.h"
!
!
!
      DO K=1, N_CONDENSED
!
         I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_MAP(TYPE_CONDENSED(K)            &
     &      , I_CLOUD_REPRESENTATION)
!
!        CHECK FOR 0 FLAGGING ILLEGAL TYPES.
         IF (I_CLOUD_TYPE(K) == 0) THEN
            WRITE(IU_ERR, '(/A)')                                       &
     &         '*** ERROR: A COMPONENT IS NOT COMPATIBLE WITH THE'      &
     &         //'REPRESENTATION OF CLOUDS SELECTED.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
!
         IF (TYPE_CONDENSED(K) == IP_CLCMP_ST_WATER) THEN
!
            I_PHASE_CMP(K)=IP_PHASE_WATER
            L_CLOUD_CMP(K)=L_DROP
!
         ELSE IF (TYPE_CONDENSED(K) == IP_CLCMP_ST_ICE) THEN
!
            I_PHASE_CMP(K)=IP_PHASE_ICE
            L_CLOUD_CMP(K)=L_ICE
!
         ELSE IF (TYPE_CONDENSED(K) == IP_CLCMP_CNV_WATER) THEN
!
            I_PHASE_CMP(K)=IP_PHASE_WATER
            L_CLOUD_CMP(K)=L_DROP
!
         ELSE IF (TYPE_CONDENSED(K) == IP_CLCMP_CNV_ICE) THEN
!
            I_PHASE_CMP(K)=IP_PHASE_ICE
            L_CLOUD_CMP(K)=L_ICE
!
         ENDIF
!
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SET_CLOUD_POINTER
#endif
#endif
