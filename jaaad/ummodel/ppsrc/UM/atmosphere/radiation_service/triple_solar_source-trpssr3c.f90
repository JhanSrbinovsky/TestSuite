
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the solar solar terms in a triple column.
!
! Method:
!       The Direct beam is calculated by propagating down through
!       the column. These direct fluxes are used to  the
!       source terms in each layer.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.2             24-05-96                Original Code
!                                               (J. M. Edwards)
!       5.3             04-10-01                Number of regions
!                                               passed explicitly.
!                                               (J. M. Edwards)
!       6.2             15-02-05                Direct flux at
!                                               ground corrected
!                                               for sloping terrain.
!                                               (James Manners)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE TRIPLE_SOLAR_SOURCE(N_PROFILE, N_LAYER, N_CLOUD_TOP    &
     &   , N_REGION, FLUX_INC_DIRECT                                    &
     &   , L_SCALE_SOLAR, ADJUST_SOLAR_KE                               &
     &   , TRANS_0, SOURCE_COEFF                                        &
     &   , V11, V12, V13, V21, V22, V23, V31, V32, V33                  &
     &   , FLUX_DIRECT                                                  &
     &   , FLUX_DIRECT_GROUND                                           &
     &   , S_UP, S_DOWN                                                 &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!
!
      USE solinc_data, ONLY: lg_orog_corr, L_orog
      IMPLICIT NONE
!
!
!     SIZES OF DUMMY ARRAYS
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     COMDECKS INCLUDED
! DIMFIX3A defines internal dimensions tied to algorithms for
! two-stream radiation code, mostly for clouds

      ! number of components of clouds
      INTEGER,PARAMETER:: NPD_CLOUD_COMPONENT=4

      ! number of permitted types of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_TYPE=4

      ! number of permitted representations of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_REPRESENTATION=4

      ! number of overlap coefficients for clouds
      INTEGER,PARAMETER:: NPD_OVERLAP_COEFF=18

      ! number of coefficients for two-stream sources
      INTEGER,PARAMETER:: NPD_SOURCE_COEFF=2

      ! number of regions in a layer
      INTEGER,PARAMETER:: NPD_REGION=3

! DIMFIX3A end
! SCFPT3A defines pointers to source coefficients in two-stream
! radiation code.

      ! pointer to source coeficient for upward solar beam
      INTEGER,PARAMETER::IP_SCF_SOLAR_UP=1

      ! pointer to source coeficient for downward solar beam
      INTEGER,PARAMETER:: IP_SCF_SOLAR_DOWN=2

      ! pointer to source coeficient for 1st difference of planckian
      INTEGER,PARAMETER:: IP_SCF_IR_1D=1

      ! pointer to source coeficient for 2nd difference of planckian
      INTEGER,PARAMETER:: IP_SCF_IR_2D=2

! SCFPT3A end
! CLDREG3A defines reference numbers for regions of clouds.in two-stream
! radiation code.

      INTEGER,PARAMETER:: IP_REGION_CLEAR=1 ! clear-sky region
      INTEGER,PARAMETER:: IP_REGION_STRAT=2 ! stratiform cloudy region
      INTEGER,PARAMETER:: IP_REGION_CONV=3  ! convective cloudy region

! CLDREG3A end
!
!
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , N_CLOUD_TOP                                                  &
!             TOP CLOUDY LAYER
     &   , N_REGION
!             Number of cloudy regions
!
!     SPECIAL ARRAYS FOR EQUIVALENT EXTINCTION:
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_SCALE_SOLAR
!             SCALING APPLIED TO SOLAR FLUX
      REAL                                                              &
                !, INTENT(IN)
     &     ADJUST_SOLAR_KE(NPD_PROFILE, NPD_LAYER)
!             ADJUSTMENT TO SOLAR FLUXES WITH EQUIVALENT EXTINCTION
!
      REAL                                                              &
                !, INTENT(IN)
     &     FLUX_INC_DIRECT(NPD_PROFILE)
!             INCIDENT DIRECT SOLAR FLUX
!
!     OPTICAL PROPERTIES:
      REAL                                                              &
                !, INTENT(IN)
     &     TRANS_0(NPD_PROFILE, NPD_LAYER, NPD_REGION)                  &
!             DIRECT TRANSMISSION
     &   , SOURCE_COEFF(NPD_PROFILE, NPD_LAYER                          &
     &      , NPD_SOURCE_COEFF, NPD_REGION)
!             SOURCE COEFFICIENTS
!
!     ENERGY TRANSFER COEFFICIENTS:
      REAL                                                              &
                !, INTENT(IN)
     &     V11(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT
     &   , V12(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT
     &   , V13(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT
     &   , V21(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT
     &   , V22(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT
     &   , V23(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT
     &   , V31(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT
     &   , V32(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT
     &   , V33(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
!
!     CALCULATED DIRECT FLUX AND SOURCE TERMS:
      REAL                                                              &
                !, INTENT(OUT)
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)                       &
!             OVERALL DIRECT FLUX
     &   , FLUX_DIRECT_GROUND(NPD_PROFILE, NPD_REGION)                  &
!             DIRECT FLUXES AT GROUND BENEATH EACH REGION
     &   , S_UP(NPD_PROFILE, NPD_LAYER, NPD_REGION)                     &
!             UPWARD SOURCE FUNCTIONS
     &   , S_DOWN(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             DOWNWARD SOURCE FUNCTIONS
!
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L                                                            &
!             LOOP VARIABLE
     &   , K
!             LOOP VARIABLE
!
      REAL                                                              &
     &     SOLAR_TOP(NPD_PROFILE, NPD_REGION)                           &
!             SOLAR FLUXES AT TOP OF LAYER
     &   , SOLAR_BASE(NPD_PROFILE, NPD_REGION)
!             SOLAR FLUXES AT BASE OF LAYER
!
!
!
!     THE CLEAR AND CLOUDY DIRECT FLUXES ARE CALCULATED SEPARATELY
!     AND ADDED TOGETHER TO FORM THE TOTAL DIRECT FLUX.
!
!     SET INCIDENT FLUXES.
      DO L=1, N_PROFILE
         FLUX_DIRECT(L, 0)=FLUX_INC_DIRECT(L)
      ENDDO
!
!     WITH EQUIVALENT EXTINCTION THE DIRECT SOLAR FLUX MUST BE
!     CORRECTED.
!
      IF (L_SCALE_SOLAR) THEN
!
         DO I=1, N_CLOUD_TOP-1
            DO L=1, N_PROFILE
               FLUX_DIRECT(L, I)                                        &
     &            =FLUX_DIRECT(L, I-1)*TRANS_0(L, I, IP_REGION_CLEAR)   &
     &            *ADJUST_SOLAR_KE(L, I)
               S_UP(L, I, IP_REGION_CLEAR)                              &
     &            =SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP, IP_REGION_CLEAR) &
     &            *FLUX_DIRECT(L, I-1)
               S_DOWN(L, I, IP_REGION_CLEAR)                            &
     &            =(SOURCE_COEFF(L, I                                   &
     &            , IP_SCF_SOLAR_DOWN, IP_REGION_CLEAR)                 &
     &            -TRANS_0(L, I, IP_REGION_CLEAR))*FLUX_DIRECT(L, I-1)  &
     &            +FLUX_DIRECT(L, I)
            ENDDO
         ENDDO
!
      ELSE
!
         DO I=1, N_CLOUD_TOP-1
            DO L=1, N_PROFILE
               FLUX_DIRECT(L, I)                                        &
     &            =FLUX_DIRECT(L, I-1)*TRANS_0(L, I, IP_REGION_CLEAR)
               S_UP(L, I, IP_REGION_CLEAR)                              &
     &            =SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP, IP_REGION_CLEAR) &
     &            *FLUX_DIRECT(L, I-1)
               S_DOWN(L, I, IP_REGION_CLEAR)                            &
     &            =SOURCE_COEFF(L, I                                    &
     &            , IP_SCF_SOLAR_DOWN, IP_REGION_CLEAR)                 &
     &            *FLUX_DIRECT(L, I-1)
            ENDDO
         ENDDO
!
      ENDIF
!
!
!
!     CLEAR AND CLOUDY REGION.
!     INITIALIZE PARTIAL FLUXES:
      DO L=1, N_PROFILE
         SOLAR_BASE(L, IP_REGION_CLEAR)=FLUX_DIRECT(L, N_CLOUD_TOP-1)
         SOLAR_BASE(L, IP_REGION_STRAT)=0.0E+00
         SOLAR_BASE(L, IP_REGION_CONV)=0.0E+00
      ENDDO
!
!
      DO I=N_CLOUD_TOP, N_LAYER
!
!        TRANSFER FLUXES ACROSS THE INTERFACE.
!
         DO L=1, N_PROFILE
            SOLAR_TOP(L, IP_REGION_CLEAR)                               &
     &         =V11(L, I-1)*SOLAR_BASE(L, IP_REGION_CLEAR)              &
     &         +V12(L, I-1)*SOLAR_BASE(L, IP_REGION_STRAT)              &
     &         +V13(L, I-1)*SOLAR_BASE(L, IP_REGION_CONV)
            SOLAR_TOP(L, IP_REGION_STRAT)                               &
     &         =V21(L, I-1)*SOLAR_BASE(L, IP_REGION_CLEAR)              &
     &         +V22(L, I-1)*SOLAR_BASE(L, IP_REGION_STRAT)              &
     &         +V23(L, I-1)*SOLAR_BASE(L, IP_REGION_CONV)
            SOLAR_TOP(L, IP_REGION_CONV)                                &
     &         =V31(L, I-1)*SOLAR_BASE(L, IP_REGION_CLEAR)              &
     &         +V32(L, I-1)*SOLAR_BASE(L, IP_REGION_STRAT)              &
     &         +V33(L, I-1)*SOLAR_BASE(L, IP_REGION_CONV)
         ENDDO
!
!
!        PROPAGATE THE FLUXES THROUGH THE LAYER:
!
         IF (L_SCALE_SOLAR) THEN
!
            DO K=1, N_REGION
               DO L=1, N_PROFILE
                  SOLAR_BASE(L, K)                                      &
     &               =SOLAR_TOP(L, K)                                   &
     &               *TRANS_0(L, I, K)*ADJUST_SOLAR_KE(L, I)
                  S_UP(L, I, K)                                         &
     &               =SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP, K)            &
     &               *SOLAR_TOP(L, K)
                  S_DOWN(L, I, K)                                       &
     &               =(SOURCE_COEFF(L, I, IP_SCF_SOLAR_DOWN, K)         &
     &               -TRANS_0(L, I, K))*SOLAR_TOP(L, K)                 &
     &               +SOLAR_BASE(L, K)
               ENDDO
            ENDDO
!
         ELSE
!
            DO K=1, N_REGION
               DO L=1, N_PROFILE
                  SOLAR_BASE(L, K)=SOLAR_TOP(L, K)                      &
     &               *TRANS_0(L, I, K)
                  S_UP(L, I, K)                                         &
     &               =SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP, K)            &
     &               *SOLAR_TOP(L, K)
                  S_DOWN(L, I, K)                                       &
     &               =SOURCE_COEFF(L, I, IP_SCF_SOLAR_DOWN, K)          &
     &               *SOLAR_TOP(L, K)
               ENDDO
            ENDDO
!
         ENDIF
!
!
!        CALCULATE THE TOTAL DIRECT FLUX.
!
         DO L=1, N_PROFILE
            FLUX_DIRECT(L, I)=SOLAR_BASE(L, IP_REGION_CLEAR)            &
     &         +SOLAR_BASE(L, IP_REGION_STRAT)                          &
     &         +SOLAR_BASE(L, IP_REGION_CONV)
         ENDDO
!
      ENDDO
!
!     PASS THE LAST VALUE AT THE BASE OF THE CLOUD OUT.
      DO K=1, N_REGION
         DO L=1, N_PROFILE
            FLUX_DIRECT_GROUND(L, K)=SOLAR_BASE(L, K)
         ENDDO
      ENDDO
!
!
!     CORRECT THE DIRECT FLUX AT THE GROUND FOR SLOPING TERRAIN

      IF (L_orog) THEN
         FLUX_DIRECT(1:N_PROFILE, N_LAYER) =                            &
     &      FLUX_DIRECT(1:N_PROFILE, N_LAYER) *                         &
     &      lg_orog_corr(1:N_PROFILE)

         DO K=1, N_REGION
            FLUX_DIRECT_GROUND(1:N_PROFILE, K) =                        &
     &         FLUX_DIRECT_GROUND(1:N_PROFILE, K) *                     &
     &         lg_orog_corr(1:N_PROFILE)

            S_DOWN(1:N_PROFILE, N_LAYER, K) =                           &
     &         S_DOWN(1:N_PROFILE, N_LAYER, K) +                        &
     &         SOLAR_BASE(1:N_PROFILE, K) *                             &
     &         (lg_orog_corr(1:N_PROFILE) - 1.0)

         ENDDO
      ENDIF

!
      RETURN
      END SUBROUTINE TRIPLE_SOLAR_SOURCE
