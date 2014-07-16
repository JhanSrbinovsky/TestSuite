#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to set number of cloudy parameters.
!
! Method:
!       Straightforward
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.5             18-05-98                Code for new parametr-
!                                               ization of droplets
!                                               included.
!                                               (J. M. Edwards)
!       5.5             24-02-03                Code for aggregate
!                                               parametrization of
!                                               ice crystals included.
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77 with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      FUNCTION SET_N_CLOUD_PARAMETER(I_SCHEME, I_COMPONENT              &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     INCLUDE COMDECKS
#include "cldcmp3a.h"
#include "wclprm3a.h"
#include "iclprm3a.h"
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_SCHEME                                                     &
!             PARAMETRIZATION SCHEME
     &   , I_COMPONENT
!             COMPONENT IN CLOUD
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     SET_N_CLOUD_PARAMETER
!             RETURNED NUMBER OF COEFFICIENTS IN PARAMETRIZATION
!
!
!
      IF ( (I_COMPONENT == IP_CLCMP_ST_WATER).OR.                       &
     &     (I_COMPONENT == IP_CLCMP_CNV_WATER) ) THEN
!
         IF (I_SCHEME == IP_SLINGO_SCHRECKER) THEN
           SET_N_CLOUD_PARAMETER =6
         ELSE IF (I_SCHEME == IP_ACKERMAN_STEPHENS) THEN
            SET_N_CLOUD_PARAMETER=9
         ELSE IF (I_SCHEME == IP_DROP_PADE_2) THEN
            SET_N_CLOUD_PARAMETER=16
         ENDIF
!
      ELSE IF ( (I_COMPONENT == IP_CLCMP_ST_ICE).OR.                    &
     &          (I_COMPONENT == IP_CLCMP_CNV_ICE) ) THEN
!
         IF (I_SCHEME == IP_SLINGO_SCHRECKER_ICE) THEN
            SET_N_CLOUD_PARAMETER=6
         ELSE IF (I_SCHEME == IP_ICE_ADT) THEN
            SET_N_CLOUD_PARAMETER=30
         ELSE IF (I_SCHEME == IP_SUN_SHINE_VN2_VIS) THEN
            SET_N_CLOUD_PARAMETER=6
         ELSE IF (I_SCHEME == IP_SUN_SHINE_VN2_IR) THEN
            SET_N_CLOUD_PARAMETER=0
         ELSE IF (I_SCHEME == IP_ICE_AGG_DE) THEN
            SET_N_CLOUD_PARAMETER=14
         ENDIF
!
      ENDIF
!
!
!
      RETURN
      END FUNCTION SET_N_CLOUD_PARAMETER
#endif
#endif
