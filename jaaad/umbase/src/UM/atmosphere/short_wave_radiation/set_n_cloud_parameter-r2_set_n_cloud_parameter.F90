#if defined(A70_1Z)
#if defined(A01_3Z) || defined(A02_3Z)
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
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                  Comment
!       6.2             13-02-06              Include code into UM
!                                             (J.-C. Thelen)
!
! Description of code:
!   FORTRAN 77 with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      FUNCTION SET_N_CLOUD_PARAMETER(I_SCHEME, I_COMPONENT, N_PHASE_TERM&
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     Include header files.
#include "c_kinds.h"
#include "cloud_component_pcf3z.h"
#include "cloud_parametrization_pcf3z.h"
#include "ice_cloud_param_pcf3z.h"
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    I_SCHEME                                                      &
!           Parametrization scheme
     &  , I_COMPONENT                                                   &
!           Component in cloud
     &  , N_PHASE_TERM
!           Number of terms in the phase function
!
      INTEGER ::                                                        &
     &    SET_N_CLOUD_PARAMETER
!           Returned number of coefficients in parametrization
!
!
!
      IF ( (I_COMPONENT == IP_CLCMP_ST_WATER).OR.                       &
     &     (I_COMPONENT == IP_CLCMP_CNV_WATER) ) THEN
!
        IF (I_SCHEME == IP_SLINGO_SCHRECKER) THEN
          SET_N_CLOUD_PARAMETER=6
        ELSE IF (I_SCHEME == IP_ACKERMAN_STEPHENS) THEN
          SET_N_CLOUD_PARAMETER=9
        ELSE IF (I_SCHEME == IP_DROP_PADE_2) THEN
          SET_N_CLOUD_PARAMETER=16
        ELSE IF (I_SCHEME == IP_SLINGO_SCHR_PHF) THEN
          SET_N_CLOUD_PARAMETER=4+2*N_PHASE_TERM
        ENDIF
!
      ELSE IF ( (I_COMPONENT == IP_CLCMP_ST_ICE).OR.                    &
     &          (I_COMPONENT == IP_CLCMP_CNV_ICE) ) THEN
!
        IF (I_SCHEME == IP_SLINGO_SCHRECKER_ICE) THEN
          SET_N_CLOUD_PARAMETER=6
        ELSE IF (I_SCHEME == IP_ICE_ADT) THEN
          SET_N_CLOUD_PARAMETER=30
        ELSE IF (I_SCHEME == IP_ICE_ADT_10) THEN
          SET_N_CLOUD_PARAMETER=36
        ELSE IF (I_SCHEME == IP_SUN_SHINE_VN2_VIS) THEN
          SET_N_CLOUD_PARAMETER=6
        ELSE IF (I_SCHEME == IP_SUN_SHINE_VN2_IR) THEN
          SET_N_CLOUD_PARAMETER=0
        ELSE IF (I_SCHEME == IP_ICE_FU_SOLAR) THEN
          SET_N_CLOUD_PARAMETER=14
        ELSE IF (I_SCHEME == IP_ICE_FU_IR) THEN
          SET_N_CLOUD_PARAMETER=10
        ELSE IF (I_SCHEME == IP_SLINGO_SCHR_ICE_PHF) THEN
          SET_N_CLOUD_PARAMETER=4+2*N_PHASE_TERM
        ELSE IF (I_SCHEME == IP_ICE_FU_PHF) THEN
          SET_N_CLOUD_PARAMETER=9+5*N_PHASE_TERM
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
