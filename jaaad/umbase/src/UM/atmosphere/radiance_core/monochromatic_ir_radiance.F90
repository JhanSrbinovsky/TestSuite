#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the infra-red radiance ignoring scattering.
!
! Method:
!       Using the secant of the ray transmission coefficients for
!       each layer may be defined and source terms may be calculated.
!       The upward and downward radiances are integrated along
!       their paths.
!
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE MONOCHROMATIC_IR_RADIANCE(N_PROFILE, N_LAYER           &
     &   , TAU                                                          &
     &   , RAD_INC_DOWN                                                 &
     &   , DIFF_PLANCK, SOURCE_GROUND, ALBEDO_SURFACE_DIFF              &
     &   , SECANT_RAY                                                   &
     &   , RADIANCE                                                     &
     &   , ND_PROFILE, ND_LAYER                                         &
     &   )
!
!
      IMPLICIT NONE
!
! Include Header Files
#include "c_kinds.h"
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Maximum number of profiles
     &  , ND_LAYER
!           Maximum number of layers
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER
!           Number of layers
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU(ND_PROFILE, ND_LAYER)                                     &
!           Optical depths of layers
     &  , RAD_INC_DOWN(ND_PROFILE)                                      &
!           Incident downward radiance
     &  , SOURCE_GROUND(ND_PROFILE)                                     &
!           Source function of ground
     &  , ALBEDO_SURFACE_DIFF(ND_PROFILE)                               &
!           Diffuse albedo
     &  , DIFF_PLANCK(ND_PROFILE, ND_LAYER)                             &
!           Difference in Planckian function
     &  , SECANT_RAY
!           Secant of angle with vertical
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    RADIANCE(ND_PROFILE, 2*ND_LAYER+2)
!           Diffuse radiance
!
!     Local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L
!           Loop variable
      REAL  (Real64) ::                                                 &
     &    TRANS(ND_PROFILE, ND_LAYER)                                   &
!           Transmissivities
     &  , SOURCE_UP(ND_PROFILE, ND_LAYER)                               &
!           Upward source function
     &  , SOURCE_DOWN(ND_PROFILE, ND_LAYER)
!           Downward source function
!
!     Variables related to the treatment of ill-conditioning
      REAL  (Real64) ::                                                 &
     &    EPS_R                                                         &
!           The smallest real number such that 1.0-EPS_R is not 1
!           to the computer's precision
     &  , SQ_EPS_R
!           The square root of the above
!
!
!
!     Set the tolerances used in avoiding ill-conditioning, testing
!     on any variable.
      EPS_R=EPSILON(TAU(1, 1))
      SQ_EPS_R=SQRT(EPS_R)
!
      DO I=1, N_LAYER
        DO L=1, N_PROFILE
          TRANS(L, I)=EXP(-SECANT_RAY*TAU(L, I))
        ENDDO
      ENDDO
!
      DO I=1, N_LAYER
        DO L=1, N_PROFILE
          SOURCE_UP(L, I)=(1.0E+00_Real64-TRANS(L, I)+SQ_EPS_R)         &
     &      *DIFF_PLANCK(L, I)                                          &
     &      /(SECANT_RAY*TAU(L, I)+SQ_EPS_R)
          SOURCE_DOWN(L, I)=-SOURCE_UP(L, I)
        ENDDO
      ENDDO
!
!     Downward radiance.
      DO L=1, N_PROFILE
        RADIANCE(L, 2)=RAD_INC_DOWN(L)
      ENDDO
      DO I=1, N_LAYER
        DO L=1, N_PROFILE
          RADIANCE(L, 2*I+2)=TRANS(L, I)*RADIANCE(L, 2*I)               &
     &      +SOURCE_DOWN(L, I)
        ENDDO
      ENDDO
!
!     Upward radiance.
      DO L=1, N_PROFILE
        RADIANCE(L, 2*N_LAYER+1)=SOURCE_GROUND(L)                       &
     &    +ALBEDO_SURFACE_DIFF(L)*RADIANCE(L, 2*N_LAYER+2)
      ENDDO
      DO I=N_LAYER, 1, -1
        DO L=1, N_PROFILE
          RADIANCE(L, 2*I-1)=TRANS(L, I)*RADIANCE(L, 2*I+1)             &
     &      +SOURCE_UP(L, I)
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE MONOCHROMATIC_IR_RADIANCE
#endif
