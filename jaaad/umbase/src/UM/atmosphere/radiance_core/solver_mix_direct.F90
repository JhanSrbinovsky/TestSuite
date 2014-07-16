#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to solve for mixed fluxes scattering without a matrix.
!
! Method:
!       Gaussian elimination in an upward direction is employed to
!       determine effective albedos for lower levels of the atmosphere.
!       This allows a downward pass of back-substitution to be carried
!       out to determine the upward and downward fluxes.
!
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOLVER_MIX_DIRECT(N_PROFILE, N_LAYER, N_CLOUD_TOP      &
     &   , T, R, S_DOWN, S_UP                                           &
     &   , T_CLOUD, R_CLOUD, S_DOWN_CLOUD, S_UP_CLOUD                   &
     &   , V11, V21, V12, V22                                           &
     &   , U11, U12, U21, U22                                           &
     &   , FLUX_INC_DOWN                                                &
     &   , SOURCE_GROUND_FREE, SOURCE_GROUND_CLOUD, ALBEDO_SURFACE_DIFF &
     &   , FLUX_TOTAL                                                   &
     &   , ND_PROFILE, ND_LAYER, ID_CT                                  &
     &   )
!
!
!
      IMPLICIT NONE
!
! Include header files
#include "c_kinds.h"
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ID_CT
!           Topmost declared cloudy layer
!
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , N_CLOUD_TOP
!           Topmost cloudy layer
      REAL  (Real64), INTENT(IN) ::                                     &
     &    T(ND_PROFILE, ND_LAYER)                                       &
!           Clear-sky transmission
     &  , R(ND_PROFILE, ND_LAYER)                                       &
!           Clear-sky reflection
     &  , S_DOWN(ND_PROFILE, ND_LAYER)                                  &
!           Clear-sky downward source function
     &  , S_UP(ND_PROFILE, ND_LAYER)                                    &
!           Clear-sky upward source function
     &  , T_CLOUD(ND_PROFILE, ND_LAYER)                                 &
!           Cloudy transmission
     &  , R_CLOUD(ND_PROFILE, ND_LAYER)                                 &
!           Cloudy reflection
     &  , S_DOWN_CLOUD(ND_PROFILE, ND_LAYER)                            &
!           Downward cloudy source function
     &  , S_UP_CLOUD(ND_PROFILE, ND_LAYER)
!           Upward cloudy source function
      REAL  (Real64), INTENT(IN) ::                                     &
     &    V11(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient
     &  , V21(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient
     &  , V12(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient
     &  , V22(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient
     &  , U11(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient
     &  , U12(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient
     &  , U21(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient
     &  , U22(ND_PROFILE, ID_CT-1: ND_LAYER)
!           Energy transfer coefficient
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FLUX_INC_DOWN(ND_PROFILE)                                     &
!           Incident flux
     &  , SOURCE_GROUND_FREE(ND_PROFILE)                                &
!           Source from ground (clear sky)
     &  , SOURCE_GROUND_CLOUD(ND_PROFILE)                               &
!           Source from ground (cloudy region)
     &  , ALBEDO_SURFACE_DIFF(ND_PROFILE)
!           Diffuse albedo
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_TOTAL(ND_PROFILE, 2*ND_LAYER+2)
!           Total flux
!
!     Local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L
!           Loop variable
!
!     Effective coupling albedos and source functions:
      REAL  (Real64) ::                                                 &
     &    ALPHA11(ND_PROFILE, ND_LAYER+1)                               &
     &  , ALPHA22(ND_PROFILE, ND_LAYER+1)                               &
     &  , ALPHA21(ND_PROFILE, ND_LAYER+1)                               &
     &  , ALPHA12(ND_PROFILE, ND_LAYER+1)                               &
     &  , G1(ND_PROFILE, ND_LAYER+1)                                    &
     &  , G2(ND_PROFILE, ND_LAYER+1)
!     Terms for downward propagation:
      REAL  (Real64) ::                                                 &
     &    GAMMA11(ND_PROFILE, ND_LAYER)                                 &
     &  , GAMMA12(ND_PROFILE, ND_LAYER)                                 &
     &  , GAMMA21(ND_PROFILE, ND_LAYER)                                 &
     &  , GAMMA22(ND_PROFILE, ND_LAYER)                                 &
     &  , BETA11_INV(ND_PROFILE, ND_LAYER)                              &
     &  , BETA21(ND_PROFILE, ND_LAYER)                                  &
     &  , BETA22_INV(ND_PROFILE, ND_LAYER)                              &
     &  , H1(ND_PROFILE, ND_LAYER)                                      &
     &  , H2(ND_PROFILE, ND_LAYER)
!
!     Auxilairy numerical variables required only in the current layer:
      REAL  (Real64) ::                                                 &
     &    THETA11                                                       &
     &  , THETA12                                                       &
     &  , THETA21                                                       &
     &  , THETA22                                                       &
     &  , LAMBDA                                                        &
     &  , LAMBDA_BAR
!
!     Temporary fluxes
      REAL  (Real64) ::                                                 &
     &    FLUX_DOWN_1(ND_PROFILE, 0: ND_LAYER)                          &
!           Downward fluxes outside clouds just below I'th level
     &  , FLUX_DOWN_2(ND_PROFILE, 0: ND_LAYER)                          &
!           Downward fluxes inside clouds just below I'th level
     &  , FLUX_UP_1(ND_PROFILE, 0: ND_LAYER)                            &
!           Upward fluxes outside clouds just above I'th level
     &  , FLUX_UP_2(ND_PROFILE, 0: ND_LAYER)
!           Upward fluxes inside clouds just above I'th level
!
!
!
!     Initialize at the bottom of the column for upward elimination.
      DO L=1, N_PROFILE
        ALPHA11(L, N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
        ALPHA22(L, N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
        ALPHA21(L, N_LAYER+1)=0.0E+00_Real64
        ALPHA12(L, N_LAYER+1)=0.0E+00_Real64
        G1(L, N_LAYER+1)=SOURCE_GROUND_FREE(L)
        G2(L, N_LAYER+1)=SOURCE_GROUND_CLOUD(L)
      ENDDO
!
!     Upward elimination through the cloudy layers.
      DO I=N_LAYER, N_CLOUD_TOP, -1
        DO L=1, N_PROFILE
!
           THETA11=ALPHA11(L, I+1)*V11(L, I)+ALPHA12(L, I+1)*V21(L, I)
           THETA12=ALPHA11(L, I+1)*V12(L, I)+ALPHA12(L, I+1)*V22(L, I)
           THETA21=ALPHA21(L, I+1)*V11(L, I)+ALPHA22(L, I+1)*V21(L, I)
           THETA22=ALPHA21(L, I+1)*V12(L, I)+ALPHA22(L, I+1)*V22(L, I)
           BETA21(L, I)=-THETA21*R(L, I)
           BETA22_INV(L, I)=1.0E+00_Real64                              &
     &       /(1.0E+00_Real64-THETA22*R_CLOUD(L, I))
           GAMMA21(L, I)=THETA21*T(L, I)
           GAMMA22(L, I)=THETA22*T_CLOUD(L, I)
           H2(L, I)=G2(L, I+1)+THETA21*S_DOWN(L, I)                     &
     &       +THETA22*S_DOWN_CLOUD(L, I)
           LAMBDA=THETA12*R_CLOUD(L, I)*BETA22_INV(L, I)
           BETA11_INV(L, I)=1.0E+00_Real64                              &
     &       /(1.0E+00_Real64-THETA11*R(L, I)+LAMBDA*BETA21(L, I))
           GAMMA11(L, I)=THETA11*T(L, I)+LAMBDA*GAMMA21(L, I)
           GAMMA12(L, I)=THETA12*T_CLOUD(L, I)+LAMBDA*GAMMA22(L, I)
           H1(L, I)=G1(L, I+1)+THETA11*S_DOWN(L, I)                     &
     &       +THETA12*S_DOWN_CLOUD(L, I)+LAMBDA*H2(L, I)
           LAMBDA=U22(L, I-1)*T_CLOUD(L, I)*BETA22_INV(L, I)
           LAMBDA_BAR=(U21(L, I-1)*T(L, I)+LAMBDA*BETA21(L, I))         &
     &       *BETA11_INV(L, I)
           ALPHA21(L, I)=U21(L, I-1)*R(L, I)+LAMBDA*GAMMA21(L, I)       &
     &       +LAMBDA_BAR*GAMMA11(L, I)
           ALPHA22(L, I)=U22(L, I-1)*R_CLOUD(L, I)                      &
     &       +LAMBDA*GAMMA22(L, I)+LAMBDA_BAR*GAMMA12(L, I)
           G2(L, I)=U21(L, I-1)*S_UP(L, I)+U22(L, I-1)*S_UP_CLOUD(L, I) &
     &       +LAMBDA*H2(L, I)+LAMBDA_BAR*H1(L, I)
!
           LAMBDA=U12(L, I-1)*T_CLOUD(L, I)*BETA22_INV(L, I)
           LAMBDA_BAR=(U11(L, I-1)*T(L, I)+LAMBDA*BETA21(L, I))         &
     &       *BETA11_INV(L, I)
           ALPHA11(L, I)=U11(L, I-1)*R(L, I)+LAMBDA*GAMMA21(L, I)       &
     &       +LAMBDA_BAR*GAMMA11(L, I)
           ALPHA12(L, I)=U12(L, I-1)*R_CLOUD(L, I)                      &
     &       +LAMBDA*GAMMA22(L, I)+LAMBDA_BAR*GAMMA12(L, I)
           G1(L, I)=U11(L, I-1)*S_UP(L, I)+U12(L, I-1)*S_UP_CLOUD(L, I) &
     &       +LAMBDA*H2(L, I)+LAMBDA_BAR*H1(L, I)
!
        ENDDO
      ENDDO
!
!     The layer above the cloud: only one set of alphas is now needed.
!     This will not be presented if there is cloud in the top layer.
!
      IF (N_CLOUD_TOP >  1) THEN
!
        I=N_CLOUD_TOP-1
        DO L=1, N_PROFILE
!
          IF (N_CLOUD_TOP <  N_LAYER) THEN
!           If there is no cloud in the column the V's will not be
!           assigned so an if test is required.
            THETA11=ALPHA11(L, I+1)*V11(L, I)+ALPHA12(L, I+1)*V21(L, I)
          ELSE
            THETA11=ALPHA11(L, I+1)
          ENDIF
!
          BETA11_INV(L,I)=1.0E+00_Real64/(1.0E+00_Real64-THETA11*R(L,I))
          GAMMA11(L, I)=THETA11*T(L, I)
          H1(L, I)=G1(L, I+1)+THETA11*S_DOWN(L, I)
!
          LAMBDA=T(L, I)*BETA11_INV(L, I)
          ALPHA11(L, I)=R(L, I)+LAMBDA*GAMMA11(L, I)
          G1(L, I)=S_UP(L, I)+LAMBDA*H1(L, I)
!
        ENDDO
!
      ENDIF
!
!
      DO I=N_CLOUD_TOP-2, 1, -1
        DO L=1, N_PROFILE
!
          BETA11_INV(L, I)=1.0E+00_Real64                               &
     &      /(1.0E+00_Real64-ALPHA11(L, I+1)*R(L, I))
          GAMMA11(L, I)=ALPHA11(L, I+1)*T(L, I)
          H1(L, I)=G1(L, I+1)+ALPHA11(L, I+1)*S_DOWN(L, I)
!
          LAMBDA=T(L, I)*BETA11_INV(L, I)
          ALPHA11(L, I)=R(L, I)+LAMBDA*GAMMA11(L, I)
          G1(L, I)=S_UP(L, I)+LAMBDA*H1(L, I)
!
        ENDDO
      ENDDO
!
!
!     Initialize for downward back-substitution.
      DO L=1, N_PROFILE
        FLUX_TOTAL(L, 2)=FLUX_INC_DOWN(L)
      ENDDO
      IF (N_CLOUD_TOP >  1) THEN
        DO L=1, N_PROFILE
          FLUX_TOTAL(L, 1)=ALPHA11(L, 1)*FLUX_TOTAL(L, 2)+G1(L, 1)
        ENDDO
      ELSE
        DO L=1, N_PROFILE
          FLUX_TOTAL(L, 1)=G1(L, 1)+FLUX_INC_DOWN(L)                    &
     &      *(V11(L, 0)*ALPHA11(L, 1)+V21(L, 0)*ALPHA12(L, 1))
        ENDDO
      ENDIF
!
!     Sweep downward through the clear-sky region, finding the downward
!     flux at the top of the layer and the upward flux at the bottom.
      DO I=1, N_CLOUD_TOP-1
        DO L=1, N_PROFILE
          FLUX_TOTAL(L, 2*I+1)=(GAMMA11(L, I)*FLUX_TOTAL(L, 2*I)        &
     &      +H1(L, I))*BETA11_INV(L, I)
          FLUX_TOTAL(L, 2*I+2)=T(L, I)*FLUX_TOTAL(L, 2*I)               &
     &      +R(L, I)*FLUX_TOTAL(L, 2*I+1)+S_DOWN(L, I)
        ENDDO
      ENDDO
!
!     Pass into the top cloudy layer. Use FLUX_DOWN_[1,2] to hold,
!     provisionally, the downward fluxes just below the top of the
!     layer, then calculate the upward fluxes at the bottom and
!     finally the downward fluxes at the bottom of the layer.
      I=N_CLOUD_TOP
      DO L=1, N_PROFILE
        FLUX_DOWN_1(L, I)=V11(L, I-1)*FLUX_TOTAL(L, 2*I)
        FLUX_DOWN_2(L, I)=V21(L, I-1)*FLUX_TOTAL(L, 2*I)
        FLUX_UP_1(L, I)=(GAMMA11(L, I)*FLUX_DOWN_1(L, I)                &
     &    +GAMMA12(L, I)*FLUX_DOWN_2(L, I)+H1(L, I))*BETA11_INV(L, I)
        FLUX_UP_2(L, I)=(GAMMA21(L, I)*FLUX_DOWN_1(L, I)                &
     &    +GAMMA22(L, I)*FLUX_DOWN_2(L, I)+H2(L, I)                     &
     &    -BETA21(L, I)*FLUX_UP_1(L, I))*BETA22_INV(L, I)
        FLUX_DOWN_1(L, I)=T(L, I)*FLUX_DOWN_1(L, I)                     &
     &    +R(L, I)*FLUX_UP_1(L, I)+S_DOWN(L, I)
        FLUX_DOWN_2(L, I)=T_CLOUD(L, I)*FLUX_DOWN_2(L, I)               &
     &    +R_CLOUD(L, I)*FLUX_UP_2(L, I)+S_DOWN_CLOUD(L, I)
      ENDDO
!
!     The main loop of back-substitution. The provisional use of the
!     downward fluxes is as above.
      DO I=N_CLOUD_TOP+1, N_LAYER
        DO L=1, N_PROFILE
          FLUX_DOWN_1(L, I)=V11(L, I-1)*FLUX_DOWN_1(L, I-1)             &
     &      +V12(L, I-1)*FLUX_DOWN_2(L, I-1)
          FLUX_DOWN_2(L, I)=V21(L, I-1)*FLUX_DOWN_1(L, I-1)             &
     &      +V22(L, I-1)*FLUX_DOWN_2(L, I-1)
          FLUX_UP_1(L, I)=(GAMMA11(L, I)*FLUX_DOWN_1(L, I)              &
     &      +GAMMA12(L, I)*FLUX_DOWN_2(L, I)+H1(L, I))                  &
     &      *BETA11_INV(L, I)
          FLUX_UP_2(L, I)=(GAMMA21(L, I)*FLUX_DOWN_1(L, I)              &
     &      +GAMMA22(L, I)*FLUX_DOWN_2(L, I)                            &
     &      -BETA21(L, I)*FLUX_UP_1(L, I)+H2(L, I))                     &
     &      *BETA22_INV(L, I)
          FLUX_DOWN_1(L, I)=T(L, I)*FLUX_DOWN_1(L, I)                   &
     &      +R(L, I)*FLUX_UP_1(L, I)+S_DOWN(L, I)
          FLUX_DOWN_2(L, I)=T_CLOUD(L, I)*FLUX_DOWN_2(L, I)             &
     &      +R_CLOUD(L, I)*FLUX_UP_2(L, I)+S_DOWN_CLOUD(L, I)
        ENDDO
      ENDDO
!
!
!     Calculate the overall flux.
      DO I=N_CLOUD_TOP, N_LAYER
        DO L=1, N_PROFILE
          FLUX_TOTAL(L, 2*I+1)=FLUX_UP_1(L, I)+FLUX_UP_2(L, I)
          FLUX_TOTAL(L, 2*I+2)=FLUX_DOWN_1(L, I)+FLUX_DOWN_2(L, I)
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SOLVER_MIX_DIRECT
#endif
