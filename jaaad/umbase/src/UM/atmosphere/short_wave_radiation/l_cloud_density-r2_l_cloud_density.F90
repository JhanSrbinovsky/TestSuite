#if defined(A70_1Z)
#if defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to determine whether densities are required for clouds.
!
! Method:
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
      FUNCTION L_CLOUD_DENSITY(N_CONDENSED, I_PHASE_CMP, L_CLOUD_CMP    &
     &   , I_CONDENSED_PARAM                                            &
     &   , ND_CLOUD_COMPONENT                                           &
     &   )
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays
      INTEGER, INTENT(IN) ::                                            &
     &    ND_CLOUD_COMPONENT
!           Size allocated for components of clouds
!
!     Include header files.
#include "c_kinds.h"
#include "cloud_parametrization_pcf3z.h"
#include "ice_cloud_param_pcf3z.h"
#include "phase_pcf3z.h"
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N_CONDENSED                                                   &
!           Number of types of condensate
     &  , I_PHASE_CMP(ND_CLOUD_COMPONENT)                               &
!           Phases of components
     &  , I_CONDENSED_PARAM(ND_CLOUD_COMPONENT)
!           Parametrizations of components
      LOGICAL, INTENT(IN) ::                                            &
     &    L_CLOUD_CMP(ND_CLOUD_COMPONENT)
!           Flags for enabled components
      LOGICAL ::                                                        &
     &    L_CLOUD_DENSITY
!           Returned flag for calculating density
!
!
!     Local variables.
      INTEGER                                                           &
     &    K
!           Loop variable
!
!
      L_CLOUD_DENSITY=.FALSE.
!
!     Densities must be calculated if Sun & Shine's parametrizations
!     are used.
      DO K=1, N_CONDENSED
        L_CLOUD_DENSITY=L_CLOUD_DENSITY                                 &
#if defined(STANDARD)
     &    .OR.( (I_PHASE_CMP(K) == IP_PHASE_WATER).AND.                 &
     &            (I_CONDENSED_PARAM(K) == IP_DROP_UNPARAMETRIZED))     &
     &    .OR.( (I_PHASE_CMP(K) == IP_PHASE_ICE).AND.                   &
     &            (I_CONDENSED_PARAM(K) == IP_ICE_UNPARAMETRIZED))      &
#endif
     &    .OR.( L_CLOUD_CMP(K).AND.(I_PHASE_CMP(K) == IP_PHASE_ICE)     &
     &          .AND.((I_CONDENSED_PARAM(K) == IP_SUN_SHINE_VN2_VIS)    &
     &            .OR.(I_CONDENSED_PARAM(K) == IP_SUN_SHINE_VN2_IR)) )
      ENDDO
!
!
!
      RETURN
      END FUNCTION L_CLOUD_DENSITY
#endif
#endif
