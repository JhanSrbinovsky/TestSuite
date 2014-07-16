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
      SUBROUTINE SOLVER_TRIPLE(N_PROFILE, N_LAYER, N_CLOUD_TOP          &
     &   , T, R, S_DOWN, S_UP                                           &
     &   , T_STRAT, R_STRAT, S_DOWN_STRAT, S_UP_STRAT                   &
     &   , T_CONV, R_CONV, S_DOWN_CONV, S_UP_CONV                       &
     &   , V11, V12, V13, V21, V22, V23, V31, V32, V33                  &
     &   , U11, U12, U13, U21, U22, U23, U31, U32, U33                  &
     &   , FLUX_INC_DOWN                                                &
     &   , SOURCE_GROUND_FREE, SOURCE_GROUND_STRAT                      &
     &   , SOURCE_GROUND_CONV, ALBEDO_SURFACE_DIFF                      &
     &   , FLUX_TOTAL                                                   &
     &   , ND_PROFILE, ND_LAYER, ID_CT                                  &
     &   )
!
!
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
!           Size allocated for atmospheric profiles
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ID_CT
!           Topmost declared cloudy layer
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
     &  , T_STRAT(ND_PROFILE, ND_LAYER)                                 &
!           Stratfiform transmission
     &  , R_STRAT(ND_PROFILE, ND_LAYER)                                 &
!           Stratfiform reflection
     &  , S_DOWN_STRAT(ND_PROFILE, ND_LAYER)                            &
!           Downward stratfiform source function
     &  , S_UP_STRAT(ND_PROFILE, ND_LAYER)                              &
!           Upward stratfiform source function
     &  , T_CONV(ND_PROFILE, ND_LAYER)                                  &
!           Convective transmission
     &  , R_CONV(ND_PROFILE, ND_LAYER)                                  &
!           Convective reflection
     &  , S_DOWN_CONV(ND_PROFILE, ND_LAYER)                             &
!           Downward convective source function
     &  , S_UP_CONV(ND_PROFILE, ND_LAYER)
!           Upward convective source function
      REAL  (Real64), INTENT(IN) ::                                     &
     &    V11(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for downward radiation
     &  , V12(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for downward radiation
     &  , V13(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for downward radiation
     &  , V21(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for downward radiation
     &  , V22(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for downward radiation
     &  , V23(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for downward radiation
     &  , V31(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for downward radiation
     &  , V32(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for downward radiation
     &  , V33(ND_PROFILE, ID_CT-1: ND_LAYER)
!           Energy transfer coefficient for downward radiation
      REAL  (Real64) ::                                                 &
     &    U11(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for upward radiation
     &  , U12(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for upward radiation
     &  , U13(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for upward radiation
     &  , U21(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for upward radiation
     &  , U22(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for upward radiation
     &  , U23(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for upward radiation
     &  , U31(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for upward radiation
     &  , U32(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for upward radiation
     &  , U33(ND_PROFILE, ID_CT-1: ND_LAYER)
!           Energy transfer coefficient for upward radiation
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FLUX_INC_DOWN(ND_PROFILE)                                     &
!           Incident flux
     &  , SOURCE_GROUND_FREE(ND_PROFILE)                                &
!           Source from ground (clear sky)
     &  , SOURCE_GROUND_STRAT(ND_PROFILE)                               &
!           Source from ground (cloudy region)
     &  , SOURCE_GROUND_CONV(ND_PROFILE)                                &
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
     &  , ALPHA12(ND_PROFILE, ND_LAYER+1)                               &
     &  , ALPHA13(ND_PROFILE, ND_LAYER+1)                               &
     &  , ALPHA21(ND_PROFILE, ND_LAYER+1)                               &
     &  , ALPHA22(ND_PROFILE, ND_LAYER+1)                               &
     &  , ALPHA23(ND_PROFILE, ND_LAYER+1)                               &
     &  , ALPHA31(ND_PROFILE, ND_LAYER+1)                               &
     &  , ALPHA32(ND_PROFILE, ND_LAYER+1)                               &
     &  , ALPHA33(ND_PROFILE, ND_LAYER+1)                               &
     &  , G1(ND_PROFILE, ND_LAYER+1)                                    &
     &  , G2(ND_PROFILE, ND_LAYER+1)                                    &
     &  , G3(ND_PROFILE, ND_LAYER+1)
!     Terms for downward propagation:
      REAL  (Real64) ::                                                 &
     &    GAMMA11(ND_PROFILE, ND_LAYER)                                 &
     &  , GAMMA12(ND_PROFILE, ND_LAYER)                                 &
     &  , GAMMA13(ND_PROFILE, ND_LAYER)                                 &
     &  , GAMMA21(ND_PROFILE, ND_LAYER)                                 &
     &  , GAMMA22(ND_PROFILE, ND_LAYER)                                 &
     &  , GAMMA23(ND_PROFILE, ND_LAYER)                                 &
     &  , GAMMA31(ND_PROFILE, ND_LAYER)                                 &
     &  , GAMMA32(ND_PROFILE, ND_LAYER)                                 &
     &  , GAMMA33(ND_PROFILE, ND_LAYER)                                 &
     &  , BETA11_INV(ND_PROFILE, ND_LAYER)                              &
     &  , BETA21(ND_PROFILE, ND_LAYER)                                  &
     &  , BETA22_INV(ND_PROFILE, ND_LAYER)                              &
     &  , BETA31(ND_PROFILE, ND_LAYER)                                  &
     &  , BETA32(ND_PROFILE, ND_LAYER)                                  &
     &  , BETA33_INV(ND_PROFILE, ND_LAYER)                              &
     &  , H1(ND_PROFILE, ND_LAYER)                                      &
     &  , H2(ND_PROFILE, ND_LAYER)                                      &
     &  , H3(ND_PROFILE, ND_LAYER)
!
!     Auxilairy numerical variables required only in the current layer:
      REAL  (Real64) ::                                                 &
     &    THETA11                                                       &
     &  , THETA12                                                       &
     &  , THETA13                                                       &
     &  , THETA21                                                       &
     &  , THETA22                                                       &
     &  , THETA23                                                       &
     &  , THETA31                                                       &
     &  , THETA32                                                       &
     &  , THETA33                                                       &
     &  , LAMBDA3                                                       &
     &  , LAMBDA2                                                       &
     &  , LAMBDA1                                                       &
     &  , LAMBDA
!
!     Temporary fluxes
      REAL  (Real64) ::                                                 &
     &    FLUX_DOWN_1(ND_PROFILE, 0: ND_LAYER)                          &
!           Downward fluxes outside clouds just below i'th level
     &  , FLUX_DOWN_2(ND_PROFILE, 0: ND_LAYER)                          &
!           Downward fluxes inside clouds just below i'th level
     &  , FLUX_DOWN_3(ND_PROFILE, 0: ND_LAYER)                          &
!           Downward fluxes inside clouds just below i'th level
     &  , FLUX_UP_1(ND_PROFILE, 0: ND_LAYER)                            &
!           Upward fluxes outside clouds just above i'th level
     &  , FLUX_UP_2(ND_PROFILE, 0: ND_LAYER)                            &
!           Upward fluxes inside clouds just above i'th level
     &  , FLUX_UP_3(ND_PROFILE, 0: ND_LAYER)
!           Upward fluxes inside clouds just above i'th level
!
!
!
!     This routine is specific to cases of three regions and it is
!     assumed that 1 represents clear skies, 2 represents startiform
!     clouds and 3 represents convective cloud.
!
!     Initialize at the bottom of the column for upward elimination.
      DO L=1, N_PROFILE
        ALPHA11(L, N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
        ALPHA12(L, N_LAYER+1)=0.0E+00_Real64
        ALPHA13(L, N_LAYER+1)=0.0E+00_Real64
        ALPHA21(L, N_LAYER+1)=0.0E+00_Real64
        ALPHA22(L, N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
        ALPHA23(L, N_LAYER+1)=0.0E+00_Real64
        ALPHA31(L, N_LAYER+1)=0.0E+00_Real64
        ALPHA32(L, N_LAYER+1)=0.0E+00_Real64
        ALPHA33(L, N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
        G1(L, N_LAYER+1)=SOURCE_GROUND_FREE(L)
        G2(L, N_LAYER+1)=SOURCE_GROUND_STRAT(L)
        G3(L, N_LAYER+1)=SOURCE_GROUND_CONV(L)
      ENDDO
!
!     Upward elimination through the cloudy layers.
      DO I=N_LAYER, N_CLOUD_TOP, -1
        DO L=1, N_PROFILE
!
          THETA11=ALPHA11(L, I+1)*V11(L, I)+ALPHA12(L, I+1)*V21(L, I)   &
     &      +ALPHA13(L, I+1)*V31(L, I)
          THETA12=ALPHA11(L, I+1)*V12(L, I)+ALPHA12(L, I+1)*V22(L, I)   &
     &      +ALPHA13(L, I+1)*V32(L, I)
          THETA13=ALPHA11(L, I+1)*V13(L, I)+ALPHA12(L, I+1)*V23(L, I)   &
     &      +ALPHA13(L, I+1)*V33(L, I)
          THETA21=ALPHA21(L, I+1)*V11(L, I)+ALPHA22(L, I+1)*V21(L, I)   &
     &      +ALPHA23(L, I+1)*V31(L, I)
          THETA22=ALPHA21(L, I+1)*V12(L, I)+ALPHA22(L, I+1)*V22(L, I)   &
     &      +ALPHA23(L, I+1)*V32(L, I)
          THETA23=ALPHA21(L, I+1)*V13(L, I)+ALPHA22(L, I+1)*V23(L, I)   &
     &      +ALPHA23(L, I+1)*V33(L, I)
          THETA31=ALPHA31(L, I+1)*V11(L, I)+ALPHA32(L, I+1)*V21(L, I)   &
     &      +ALPHA33(L, I+1)*V31(L, I)
          THETA32=ALPHA31(L, I+1)*V12(L, I)+ALPHA32(L, I+1)*V22(L, I)   &
     &      +ALPHA33(L, I+1)*V32(L, I)
          THETA33=ALPHA31(L, I+1)*V13(L, I)+ALPHA32(L, I+1)*V23(L, I)   &
     &      +ALPHA33(L, I+1)*V33(L, I)
          BETA31(L, I)=-THETA31*R(L, I)
          BETA32(L, I)=-THETA32*R_STRAT(L, I)
          BETA33_INV(L, I)=1.0E+00_Real64                               &
     &      /(1.0E+00_Real64-THETA33*R_CONV(L, I))
          GAMMA31(L, I)=THETA31*T(L, I)
          GAMMA32(L, I)=THETA32*T_STRAT(L, I)
          GAMMA33(L, I)=THETA33*T_CONV(L, I)
          H3(L, I)=G3(L, I+1)+THETA31*S_DOWN(L, I)                      &
     &      +THETA32*S_DOWN_STRAT(L, I)                                 &
     &      +THETA33*S_DOWN_CONV(L, I)
!
          LAMBDA3=THETA23*R_CONV(L, I)*BETA33_INV(L, I)
          BETA22_INV(L, I)=1.0E+00_Real64                               &
     &      /(1.0E+00_Real64-THETA22*R_STRAT(L, I)+LAMBDA3*BETA32(L, I))
          BETA21(L, I)=-THETA21*R(L, I)+LAMBDA3*BETA31(L, I)
          GAMMA21(L, I)=THETA21*T(L, I)+LAMBDA3*GAMMA31(L, I)
          GAMMA22(L, I)=THETA22*T_STRAT(L, I)+LAMBDA3*GAMMA32(L, I)
          GAMMA23(L, I)=THETA23*T_CONV(L, I)+LAMBDA3*GAMMA33(L, I)
          H2(L, I)=G2(L, I+1)+THETA21*S_DOWN(L, I)                      &
     &      +THETA22*S_DOWN_STRAT(L, I)+THETA23*S_DOWN_CONV(L, I)       &
     &      +LAMBDA3*H3(L, I)
!
          LAMBDA3=THETA13*R_CONV(L, I)*BETA33_INV(L, I)
          LAMBDA2=(THETA12*R_STRAT(L, I)-LAMBDA3*BETA32(L, I))          &
     &      *BETA22_INV(L, I)
          BETA11_INV(L, I)=1.0E+00_Real64                               &
     &      /(1.0E+00_Real64-THETA11*R(L, I)+LAMBDA3*BETA31(L, I)       &
     &      +LAMBDA2*BETA21(L, I))
          GAMMA11(L, I)=THETA11*T(L, I)+LAMBDA3*GAMMA31(L, I)           &
     &      +LAMBDA2*GAMMA21(L, I)
          GAMMA12(L, I)=THETA12*T_STRAT(L, I)+LAMBDA3*GAMMA32(L, I)     &
     &      +LAMBDA2*GAMMA22(L, I)
          GAMMA13(L, I)=THETA13*T_CONV(L, I)+LAMBDA3*GAMMA33(L, I)      &
     &      +LAMBDA2*GAMMA23(L, I)
          H1(L, I)=G1(L, I+1)+THETA11*S_DOWN(L, I)                      &
     &      +THETA12*S_DOWN_STRAT(L, I)+THETA13*S_DOWN_CONV(L, I)       &
     &      +LAMBDA3*H3(L, I)+LAMBDA2*H2(L, I)
!
          LAMBDA3=U33(L, I-1)*T_CONV(L, I)*BETA33_INV(L, I)
          LAMBDA2=(U32(L, I-1)*T_STRAT(L, I)+LAMBDA3*BETA32(L, I))      &
     &      *BETA22_INV(L, I)
          LAMBDA1=(U31(L, I-1)*T(L, I)+LAMBDA3*BETA31(L, I)             &
     &      +LAMBDA2*BETA21(L, I))*BETA11_INV(L, I)
          ALPHA31(L, I)=U31(L, I-1)*R(L, I)+LAMBDA3*GAMMA31(L, I)       &
     &      +LAMBDA2*GAMMA21(L, I)+LAMBDA1*GAMMA11(L, I)
          ALPHA32(L, I)=U32(L, I-1)*R_STRAT(L, I)                       &
     &      +LAMBDA3*GAMMA32(L, I)+LAMBDA2*GAMMA22(L, I)                &
     &      +LAMBDA1*GAMMA12(L, I)
          ALPHA33(L, I)=U33(L, I-1)*R_CONV(L, I)                        &
     &      +LAMBDA3*GAMMA33(L, I)+LAMBDA2*GAMMA23(L, I)                &
     &      +LAMBDA1*GAMMA13(L, I)
          G3(L, I)=U31(L, I-1)*S_UP(L, I)+U32(L, I-1)*S_UP_STRAT(L, I)  &
     &      +U33(L, I-1)*S_UP_CONV(L, I)                                &
     &      +LAMBDA3*H3(L, I)+LAMBDA2*H2(L, I)+LAMBDA1*H1(L, I)
!
          LAMBDA3=U23(L, I-1)*T_CONV(L, I)*BETA33_INV(L, I)
          LAMBDA2=(U22(L, I-1)*T_STRAT(L, I)+LAMBDA3*BETA32(L, I))      &
     &      *BETA22_INV(L, I)
          LAMBDA1=(U21(L, I-1)*T(L, I)+LAMBDA3*BETA31(L, I)             &
     &      +LAMBDA2*BETA21(L, I))*BETA11_INV(L, I)
          ALPHA21(L, I)=U21(L, I-1)*R(L, I)+LAMBDA3*GAMMA31(L, I)       &
     &      +LAMBDA2*GAMMA21(L, I)+LAMBDA1*GAMMA11(L, I)
          ALPHA22(L, I)=U22(L, I-1)*R_STRAT(L, I)                       &
     &      +LAMBDA3*GAMMA32(L, I)+LAMBDA2*GAMMA22(L, I)                &
     &      +LAMBDA1*GAMMA12(L, I)
          ALPHA23(L, I)=U23(L, I-1)*R_CONV(L, I)                        &
     &      +LAMBDA3*GAMMA33(L, I)+LAMBDA2*GAMMA23(L, I)                &
     &      +LAMBDA1*GAMMA13(L, I)
          G2(L, I)=U21(L, I-1)*S_UP(L, I)+U22(L, I-1)*S_UP_STRAT(L, I)  &
     &      +U23(L, I-1)*S_UP_CONV(L, I)                                &
     &      +LAMBDA3*H3(L, I)+LAMBDA2*H2(L, I)+LAMBDA1*H1(L, I)
!
          LAMBDA3=U13(L, I-1)*T_CONV(L, I)*BETA33_INV(L, I)
          LAMBDA2=(U12(L, I-1)*T_STRAT(L, I)+LAMBDA3*BETA32(L, I))      &
     &      *BETA22_INV(L, I)
          LAMBDA1=(U11(L, I-1)*T(L, I)+LAMBDA3*BETA31(L, I)             &
     &      +LAMBDA2*BETA21(L, I))*BETA11_INV(L, I)
          ALPHA11(L, I)=U11(L, I-1)*R(L, I)+LAMBDA3*GAMMA31(L, I)       &
     &      +LAMBDA2*GAMMA21(L, I)+LAMBDA1*GAMMA11(L, I)
          ALPHA12(L, I)=U12(L, I-1)*R_STRAT(L, I)                       &
     &      +LAMBDA3*GAMMA32(L, I)+LAMBDA2*GAMMA22(L, I)                &
     &      +LAMBDA1*GAMMA12(L, I)
          ALPHA13(L, I)=U13(L, I-1)*R_CONV(L, I)                        &
     &      +LAMBDA3*GAMMA33(L, I)+LAMBDA2*GAMMA23(L, I)                &
     &      +LAMBDA1*GAMMA13(L, I)
          G1(L, I)=U11(L, I-1)*S_UP(L, I)+U12(L, I-1)*S_UP_STRAT(L, I)  &
     &      +U13(L, I-1)*S_UP_CONV(L, I)                                &
     &      +LAMBDA3*H3(L, I)+LAMBDA2*H2(L, I)+LAMBDA1*H1(L, I)
!
        ENDDO
      ENDDO
!
!     The layer above the cloud, if the cloud does not reach to the top
!     of the column: only one set of alphas is now needed.
!
      IF (N_CLOUD_TOP >  1) THEN
!
        I=N_CLOUD_TOP-1
        DO L=1, N_PROFILE
!
          IF (N_CLOUD_TOP <  N_LAYER) THEN
!           If there is no cloud in the column the V's will not be
!           assigned so an IF-test is required.
            THETA11=ALPHA11(L, I+1)*V11(L, I)+ALPHA12(L, I+1)*V21(L, I) &
     &        +ALPHA13(L, I+1)*V31(L, I)
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
!     Continue through the cloud free region: if there is no such
!     region the DO-loop will not be executed.
      DO I=N_CLOUD_TOP-2, 1, -1
        DO L=1, N_PROFILE
!
          BETA11_INV(L, I)=1.0E+00_Real64                               &
     &      /(1.0E+00_Real64-ALPHA11(L, I+1)*R(L, I))
          GAMMA11(L, I)=ALPHA11(L, I+1)*T(L, I)
          H1(L, I)=G1(L, I+1)+ALPHA11(L, I+1)*S_DOWN(L, I)
!
          LAMBDA1=T(L, I)*BETA11_INV(L, I)
          ALPHA11(L, I)=R(L, I)+LAMBDA1*GAMMA11(L, I)
          G1(L, I)=S_UP(L, I)+LAMBDA1*H1(L, I)
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
     &      *(V11(L, 0)*ALPHA11(L, 1)+V21(L, 0)*ALPHA12(L, 1)           &
     &      + V31(L, 0)*ALPHA13(L, 1))
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
!     Pass into the top cloudy layer. use flux_down_[1,2,3] to hold,
!     provisionally, the downward fluxes just below the top of the
!     layer, then calculate the upward fluxes at the bottom and
!     finally the downward fluxes at the bottom of the layer.
      I=N_CLOUD_TOP
      DO L=1, N_PROFILE
        FLUX_DOWN_1(L, I)=V11(L, I-1)*FLUX_TOTAL(L, 2*I)
        FLUX_DOWN_2(L, I)=V21(L, I-1)*FLUX_TOTAL(L, 2*I)
        FLUX_DOWN_3(L, I)=V31(L, I-1)*FLUX_TOTAL(L, 2*I)
        FLUX_UP_1(L, I)=(GAMMA11(L, I)*FLUX_DOWN_1(L, I)                &
     &    +GAMMA12(L, I)*FLUX_DOWN_2(L, I)                              &
     &    +GAMMA13(L, I)*FLUX_DOWN_3(L, I)                              &
     &    +H1(L, I))*BETA11_INV(L, I)
        FLUX_UP_2(L, I)=(GAMMA21(L, I)*FLUX_DOWN_1(L, I)                &
     &    +GAMMA22(L, I)*FLUX_DOWN_2(L, I)                              &
     &    +GAMMA23(L, I)*FLUX_DOWN_3(L, I)+H2(L, I)                     &
     &    -BETA21(L, I)*FLUX_UP_1(L, I))*BETA22_INV(L, I)
        FLUX_UP_3(L, I)=(GAMMA31(L, I)*FLUX_DOWN_1(L, I)                &
     &    +GAMMA32(L, I)*FLUX_DOWN_2(L, I)                              &
     &    +GAMMA33(L, I)*FLUX_DOWN_3(L, I)+H3(L, I)                     &
     &    -BETA31(L, I)*FLUX_UP_1(L, I)-BETA32(L, I)*FLUX_UP_2(L, I))   &
     &    *BETA33_INV(L, I)
        FLUX_DOWN_1(L, I)=T(L, I)*FLUX_DOWN_1(L, I)                     &
     &    +R(L, I)*FLUX_UP_1(L, I)+S_DOWN(L, I)
        FLUX_DOWN_2(L, I)=T_STRAT(L, I)*FLUX_DOWN_2(L, I)               &
     &    +R_STRAT(L, I)*FLUX_UP_2(L, I)+S_DOWN_STRAT(L, I)
        FLUX_DOWN_3(L, I)=T_CONV(L, I)*FLUX_DOWN_3(L, I)                &
     &    +R_CONV(L, I)*FLUX_UP_3(L, I)+S_DOWN_CONV(L, I)
      ENDDO
!
!     The main loop of back-substitution. the provisional use of the
!     downward fluxes is as above.
      DO I=N_CLOUD_TOP+1, N_LAYER
        DO L=1, N_PROFILE
          FLUX_DOWN_1(L, I)=V11(L, I-1)*FLUX_DOWN_1(L, I-1)             &
     &      +V12(L, I-1)*FLUX_DOWN_2(L, I-1)                            &
     &      +V13(L, I-1)*FLUX_DOWN_3(L, I-1)
          FLUX_DOWN_2(L, I)=V21(L, I-1)*FLUX_DOWN_1(L, I-1)             &
     &      +V22(L, I-1)*FLUX_DOWN_2(L, I-1)                            &
     &      +V23(L, I-1)*FLUX_DOWN_3(L, I-1)
          FLUX_DOWN_3(L, I)=V31(L, I-1)*FLUX_DOWN_1(L, I-1)             &
     &      +V32(L, I-1)*FLUX_DOWN_2(L, I-1)                            &
     &      +V33(L, I-1)*FLUX_DOWN_3(L, I-1)
          FLUX_UP_1(L, I)=(GAMMA11(L, I)*FLUX_DOWN_1(L, I)              &
     &      +GAMMA12(L, I)*FLUX_DOWN_2(L, I)                            &
     &      +GAMMA13(L, I)*FLUX_DOWN_3(L, I)                            &
     &      +H1(L, I))*BETA11_INV(L, I)
          FLUX_UP_2(L, I)=(GAMMA21(L, I)*FLUX_DOWN_1(L, I)              &
     &      +GAMMA22(L, I)*FLUX_DOWN_2(L, I)                            &
     &      +GAMMA23(L, I)*FLUX_DOWN_3(L, I)+H2(L, I)                   &
     &      -BETA21(L, I)*FLUX_UP_1(L, I))*BETA22_INV(L, I)
          FLUX_UP_3(L, I)=(GAMMA31(L, I)*FLUX_DOWN_1(L, I)              &
     &      +GAMMA32(L, I)*FLUX_DOWN_2(L, I)                            &
     &      +GAMMA33(L, I)*FLUX_DOWN_3(L, I)+H3(L, I)                   &
     &      -BETA31(L, I)*FLUX_UP_1(L, I)                               &
     &      -BETA32(L, I)*FLUX_UP_2(L, I))                              &
     &      *BETA33_INV(L, I)
          FLUX_DOWN_1(L, I)=T(L, I)*FLUX_DOWN_1(L, I)                   &
     &      +R(L, I)*FLUX_UP_1(L, I)+S_DOWN(L, I)
          FLUX_DOWN_2(L, I)=T_STRAT(L, I)*FLUX_DOWN_2(L, I)             &
     &      +R_STRAT(L, I)*FLUX_UP_2(L, I)+S_DOWN_STRAT(L, I)
          FLUX_DOWN_3(L, I)=T_CONV(L, I)*FLUX_DOWN_3(L, I)              &
     &      +R_CONV(L, I)*FLUX_UP_3(L, I)+S_DOWN_CONV(L, I)
        ENDDO
      ENDDO
!
!
!     Calculate the overall flux.
      DO I=N_CLOUD_TOP, N_LAYER
        DO L=1, N_PROFILE
          FLUX_TOTAL(L, 2*I+1)=FLUX_UP_1(L, I)+FLUX_UP_2(L, I)          &
     &      +FLUX_UP_3(L, I)
          FLUX_TOTAL(L, 2*I+2)=FLUX_DOWN_1(L, I)+FLUX_DOWN_2(L, I)      &
     &      +FLUX_DOWN_3(L, I)
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SOLVER_TRIPLE
#endif
