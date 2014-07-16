#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
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
!       This version has been modified by Robin Hogan to allow 
!       shadowing, as documented in Shonk & Hogan, 2007, J. Climate.
!
! Current Owner of Code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOLVER_MIX_DIRECT_HOGAN(N_PROFILE, N_LAYER, N_CLOUD_TOP&
     &   , T, R, S_DOWN, S_UP                                           &
     &   , T_CLOUD, R_CLOUD, S_DOWN_CLOUD, S_UP_CLOUD                   &
     &   , V11, V21, V12, V22                                           &
     &   , U11, U12, U21, U22                                           &
     &   , L_NET                                                        &
     &   , FLUX_INC_DOWN                                                &
     &   , SOURCE_GROUND_FREE, SOURCE_GROUND_CLOUD, ALBEDO_SURFACE_DIFF &
     &   , FLUX_TOTAL                                                   &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!
      IMPLICIT NONE
!
!
!     SIZES OF DUMMY ARRAYS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , N_CLOUD_TOP
!             TOPMOST CLOUDY LAYER
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_NET
!             FLAG FOR CALCULATION OF NET FLUXES
      REAL                                                              &
                !, INTENT(IN)
     &     T(NPD_PROFILE, NPD_LAYER)                                    &
!             CLEAR-SKY TRANSMISSION
     &   , R(NPD_PROFILE, NPD_LAYER)                                    &
!             CLEAR-SKY REFLECTION
     &   , S_DOWN(NPD_PROFILE, NPD_LAYER)                               &
!             CLEAR-SKY DOWNWARD SOURCE FUNCTION
     &   , S_UP(NPD_PROFILE, NPD_LAYER)                                 &
!             CLEAR-SKY UPWARD SOURCE FUNCTION
     &   , T_CLOUD(NPD_PROFILE, NPD_LAYER)                              &
!             CLOUDY TRANSMISSION
     &   , R_CLOUD(NPD_PROFILE, NPD_LAYER)                              &
!             CLOUDY REFLECTION
     &   , S_DOWN_CLOUD(NPD_PROFILE, NPD_LAYER)                         &
!             DOWNWARD CLOUDY SOURCE FUNCTION
     &   , S_UP_CLOUD(NPD_PROFILE, NPD_LAYER)
!             UPWARD CLOUDY SOURCE FUNCTION
      REAL                                                              &
                !, INTENT(IN)
     &     V11(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT
     &   , V21(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT
     &   , V12(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT
     &   , V22(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT
     &   , U11(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT
     &   , U12(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT
     &   , U21(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT
     &   , U22(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
      REAL                                                              &
                !, INTENT(IN)
     &     FLUX_INC_DOWN(NPD_PROFILE)                                   &
!             INCIDENT FLUX
     &   , SOURCE_GROUND_FREE(NPD_PROFILE)                              &
!             SOURCE FROM GROUND (CLEAR SKY)
     &   , SOURCE_GROUND_CLOUD(NPD_PROFILE)                             &
!             SOURCE FROM GROUND (CLOUDY REGION)
     &   , ALBEDO_SURFACE_DIFF(NPD_PROFILE)
!             DIFFUSE ALBEDO
      REAL                                                              &
                !, INTENT(OUT)
     &     FLUX_TOTAL(NPD_PROFILE, 2*NPD_LAYER+2)
!             TOTAL FLUX
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
!     EFFECTIVE COUPLING ALBEDOS AND SOURCE FUNCTIONS:
      REAL                                                              &
     &     ALPHA11(NPD_PROFILE, NPD_LAYER+1)                            &
     &   , ALPHA22(NPD_PROFILE, NPD_LAYER+1)                            &
     &   , G1(NPD_PROFILE, NPD_LAYER+1)                                 &
     &   , G2(NPD_PROFILE, NPD_LAYER+1)
!     TERMS FOR DOWNWARD PROPAGATION:
      REAL                                                              &
     &     GAMMA11(NPD_PROFILE, NPD_LAYER)                              &
     &   , GAMMA22(NPD_PROFILE, NPD_LAYER)                              &
     &   , BETA11_INV(NPD_PROFILE, NPD_LAYER)                           &
     &   , BETA22_INV(NPD_PROFILE, NPD_LAYER)                           &
     &   , H1(NPD_PROFILE, NPD_LAYER)                                   &
     &   , H2(NPD_PROFILE, NPD_LAYER)
!
!     AUXILAIRY NUMERICAL VARIABLES REQUIRED ONLY IN THE CURRENT LAYER:
      REAL                                                              &
     &     THETA11                                                      &
     &   , THETA22                                                      &
     &   , LAMBDA1                                                      &
     &   , LAMBDA2                                                      &
     &   , LAMBDA
!
!     TEMPORARY FLUXES
      REAL                                                              &
     &     FLUX_DOWN_1(NPD_PROFILE, 0: NPD_LAYER)                       &
!             DOWNWARD FLUXES OUTSIDE CLOUDS JUST BELOW I'TH LEVEL
     &   , FLUX_DOWN_2(NPD_PROFILE, 0: NPD_LAYER)                       &
!             DOWNWARD FLUXES INSIDE CLOUDS JUST BELOW I'TH LEVEL
     &   , FLUX_UP_1(NPD_PROFILE, 0: NPD_LAYER)                         &
!             UPWARD FLUXES OUTSIDE CLOUDS JUST ABOVE I'TH LEVEL
     &   , FLUX_UP_2(NPD_PROFILE, 0: NPD_LAYER)
!             UPWARD FLUXES INSIDE CLOUDS JUST ABOVE I'TH LEVEL
!
!
!
!     INITIALIZE AT THE BOTTOM OF THE COLUMN FOR UPWARD ELIMINATION.
      DO L=1, N_PROFILE
         ALPHA11(L, N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
         ALPHA22(L, N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
         G1(L, N_LAYER+1)=SOURCE_GROUND_FREE(L)
         G2(L, N_LAYER+1)=SOURCE_GROUND_CLOUD(L)
      ENDDO
!
!     UPWARD ELIMINATION THROUGH THE CLOUDY LAYERS.
      DO I=N_LAYER, N_CLOUD_TOP, -1
         DO L=1, N_PROFILE
!
            THETA11=ALPHA11(L, I+1)*V11(L, I)+ALPHA22(L, I+1)*V21(L, I)
            THETA22=ALPHA11(L, I+1)*V12(L, I)+ALPHA22(L, I+1)*V22(L, I)

            BETA11_INV(L, I)=1.0E+00/(1.0E+00-THETA11*R(L, I))
            GAMMA11(L, I)=THETA11*T(L, I)
            H1(L, I)=G1(L, I+1)+THETA11*S_DOWN(L, I)

            BETA22_INV(L, I)=1.0E+00/(1.0E+00-THETA22*R_CLOUD(L, I))
            GAMMA22(L, I)=THETA22*T_CLOUD(L, I)
            H2(L, I)=G2(L, I+1)+THETA22*S_DOWN_CLOUD(L, I)

            LAMBDA1 = S_UP(L, I)+H1(L, I)*T(L, I)*BETA11_INV(L, I)
            LAMBDA2 = S_UP_CLOUD(L, I)+H2(L, I)*T_CLOUD(L, I)           &
     &           *BETA22_INV(L, I)

            ALPHA11(L, I)=R(L, I)                                       &
     &           + THETA11*T(L, I)*T(L, I)*BETA11_INV(L, I)
            G1(L, I)=U11(L, I-1)*LAMBDA1 + U12(L, I-1)*LAMBDA2

            ALPHA22(L, I)=R_CLOUD(L, I)                                 &
     &           + THETA22*T_CLOUD(L, I)*T_CLOUD(L, I)*BETA22_INV(L, I)
            G2(L, I)=U21(L, I-1)*LAMBDA1 + U22(L, I-1)*LAMBDA2
         ENDDO
      ENDDO
!
      IF (N_CLOUD_TOP >  1) THEN
!
!     THE LAYER ABOVE THE CLOUD: ONLY ONE SET OF ALPHAS IS NOW NEEDED.
!
         DO L=1, N_PROFILE
!
            IF (N_CLOUD_TOP <  N_LAYER) THEN
!              If all columns are clear down to the surface, the
!              coefficients V11 etc. will not be set, so the case when
!              N_CLOUD_TOP equals N_LAYER must be treated separately.
               THETA11=ALPHA11(L, I+1)*V11(L, I)                        &
     &            +ALPHA22(L, I+1)*V21(L, I)
            ELSE
               THETA11=ALPHA11(L, I+1)
            ENDIF
!
            BETA11_INV(L, I)=1.0E+00/(1.0E+00-THETA11*R(L, I))
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
      DO I=N_CLOUD_TOP-2, 1, -1
         DO L=1, N_PROFILE
!
            BETA11_INV(L, I)=1.0E+00/(1.0E+00-ALPHA11(L, I+1)*R(L, I))
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
!     Initialize for downward back-substitution. If there is cloud
!     in the top layer the upward flux must be calculated allowing
!     for reflection from clear air and the cloud.
      DO L=1, N_PROFILE
         FLUX_TOTAL(L, 2)=FLUX_INC_DOWN(L)
      ENDDO
      IF (N_CLOUD_TOP >  1) THEN
         DO L=1, N_PROFILE
            FLUX_TOTAL(L, 1)=G1(L, 1)+ALPHA11(L, 1)*FLUX_INC_DOWN(L)
         ENDDO
      ELSE
         DO L=1, N_PROFILE
            FLUX_TOTAL(L, 1)=G1(L, 1)+FLUX_INC_DOWN(L)                  &
     &         *(V11(L, 0)*ALPHA11(L, 1)+V21(L, 0)*ALPHA22(L, 1))
         ENDDO
      ENDIF
!
!     SWEEP DOWNWARD THROUGH THE CLEAR-SKY REGION, FINDING THE DOWNWARD
!     FLUX AT THE TOP OF THE LAYER AND THE UPWARD FLUX AT THE BOTTOM.
      DO I=1, N_CLOUD_TOP-1
         DO L=1, N_PROFILE
            FLUX_TOTAL(L, 2*I+1)=(GAMMA11(L, I)*FLUX_TOTAL(L, 2*I)      &
     &         +H1(L, I))*BETA11_INV(L, I)
            FLUX_TOTAL(L, 2*I+2)=T(L, I)*FLUX_TOTAL(L, 2*I)             &
     &         +R(L, I)*FLUX_TOTAL(L, 2*I+1)+S_DOWN(L, I)
         ENDDO
      ENDDO
!
!     PASS INTO THE TOP CLOUDY LAYER. USE FLUX_DOWN_[1,2] TO HOLD,
!     PROVISIONALLY, THE DOWNWARD FLUXES JUST BELOW THE TOP OF THE
!     LAYER, THEN CALCULATE THE UPWARD FLUXES AT THE BOTTOM AND
!     FINALLY THE DOWNWARD FLUXES AT THE BOTTOM OF THE LAYER.
      I=N_CLOUD_TOP
      DO L=1, N_PROFILE
         FLUX_DOWN_1(L, I)=V11(L, I-1)*FLUX_TOTAL(L, 2*I)
         FLUX_DOWN_2(L, I)=V21(L, I-1)*FLUX_TOTAL(L, 2*I)
         FLUX_UP_1(L, I)=(GAMMA11(L, I)*FLUX_DOWN_1(L, I)               &
     &      +H1(L, I))*BETA11_INV(L, I)
         FLUX_UP_2(L, I)=(GAMMA22(L, I)*FLUX_DOWN_2(L, I)               &
     &        +H2(L, I))*BETA22_INV(L, I)
         FLUX_DOWN_1(L, I)=T(L, I)*FLUX_DOWN_1(L, I)                    &
     &      +R(L, I)*FLUX_UP_1(L, I)+S_DOWN(L, I)
         FLUX_DOWN_2(L, I)=T_CLOUD(L, I)*FLUX_DOWN_2(L, I)              &
     &      +R_CLOUD(L, I)*FLUX_UP_2(L, I)+S_DOWN_CLOUD(L, I)
      ENDDO
!
!     THE MAIN LOOP OF BACK-SUBSTITUTION. THE PROVISIONAL USE OF THE
!     DOWNWARD FLUXES IS AS ABOVE.
      DO I=N_CLOUD_TOP+1, N_LAYER
         DO L=1, N_PROFILE
            FLUX_DOWN_1(L, I)=V11(L, I-1)*FLUX_DOWN_1(L, I-1)           &
     &         +V12(L, I-1)*FLUX_DOWN_2(L, I-1)
            FLUX_DOWN_2(L, I)=V21(L, I-1)*FLUX_DOWN_1(L, I-1)           &
     &         +V22(L, I-1)*FLUX_DOWN_2(L, I-1)
            FLUX_UP_1(L, I)=(GAMMA11(L, I)*FLUX_DOWN_1(L, I)            &
     &         +H1(L, I))*BETA11_INV(L, I)
            FLUX_UP_2(L, I)=(GAMMA22(L, I)*FLUX_DOWN_2(L, I)            &
     &         +H2(L, I))*BETA22_INV(L, I)
            FLUX_DOWN_1(L, I)=T(L, I)*FLUX_DOWN_1(L, I)                 &
     &         +R(L, I)*FLUX_UP_1(L, I)+S_DOWN(L, I)
            FLUX_DOWN_2(L, I)=T_CLOUD(L, I)*FLUX_DOWN_2(L, I)           &
     &         +R_CLOUD(L, I)*FLUX_UP_2(L, I)+S_DOWN_CLOUD(L, I)
         ENDDO
      ENDDO
!
!
!     CALCULATE THE OVERALL FLUX.
      DO I=N_CLOUD_TOP, N_LAYER
         DO L=1, N_PROFILE
            FLUX_TOTAL(L, 2*I+1)=FLUX_UP_1(L, I)+FLUX_UP_2(L, I)
            FLUX_TOTAL(L, 2*I+2)=FLUX_DOWN_1(L, I)+FLUX_DOWN_2(L, I)
         ENDDO
      ENDDO
!
!     REDUCE TO NET FLUXES IF REQUIRED.
      IF (L_NET) THEN
         DO I=0, N_LAYER
            DO L=1, N_PROFILE
               FLUX_TOTAL(L, I+1)                                       &
     &            =FLUX_TOTAL(L, 2*I+2)-FLUX_TOTAL(L, 2*I+1)
            ENDDO
         ENDDO
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE SOLVER_MIX_DIRECT_HOGAN
#endif
#endif
