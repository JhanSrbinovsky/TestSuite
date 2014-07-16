


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate transmission and reflection coefficients.
!
! Method:
!        Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.5             11-06-98                Optimised Code
!                                               (P. Burton)
!       5.1             04-04-00                Tolerances replaced
!                                               by F90 intrinsics.
!                                               (J. M. Edwards)
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE TRANS_SOURCE_COEFF(N_PROFILE                           &
     &   , I_LAYER_FIRST, I_LAYER_LAST                                  &
     &   , ISOLIR, L_IR_SOURCE_QUAD                                     &
     &   , TAU, SUM, DIFF, LAMBDA, SEC_0                                &
     &   , GAMMA_UP, GAMMA_DOWN                                         &
     &   , TRANS, REFLECT, TRANS_0, SOURCE_COEFF                        &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!
      USE rad_switches_mod, ONLY: lrad_quad_src_fix
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
! SPCRG3A defines flags for different portions of the spectrum in
! two-stream radiation code.
      INTEGER,PARAMETER:: IP_SOLAR=1
      INTEGER,PARAMETER:: IP_INFRA_RED=2
! SPCRG3A end
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
!
!     DUMMY VARIABLES.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , I_LAYER_FIRST                                                &
!             FIRST LAYER TO CONSIDER
     &   , I_LAYER_LAST
!             LAST LAYER TO CONSIDER
!
!     ALGORITHMIC control
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_IR_SOURCE_QUAD
!             QUADRATIC SOURCE IN INFRA-RED
      INTEGER                                                           &
                !, INTENT(IN)
     &     ISOLIR
!             SPECTRAL REGION
!
!     OPTICAL PROPERTIES OF THE LAYER
      REAL                                                              &
                !, INTENT(IN)
     &     TAU(NPD_PROFILE, NPD_LAYER)                                  &
!             OPTICAL DEPTHS OF LAYERS
     &   , SUM(NPD_PROFILE, NPD_LAYER)                                  &
!             SUM OF ALPHA_1 AND ALPHA_2
     &   , DIFF(NPD_PROFILE, NPD_LAYER)                                 &
!             DIFFERENCE OF ALPHA_1 AND ALPHA_2
     &   , LAMBDA(NPD_PROFILE, NPD_LAYER)                               &
!             LAMBDA
     &   , SEC_0(NPD_PROFILE)                                           &
!             SECANT OF SOLAR ZENITH ANGLE
     &   , GAMMA_UP(NPD_PROFILE, NPD_LAYER)                             &
!             BASIC SOLAR COEFFICIENT FOR UPWARD RADIATION
     &   , GAMMA_DOWN(NPD_PROFILE, NPD_LAYER)
!             BASIC SOLAR COEFFICIENT FOR DOWNWARD RADIATION
!
!     TRANSMISSION AND REFLECTION COEFFICIENTS AND COEFFICIENTS FOR
!     SOURCE TERMS.
      REAL                                                              &
                !, INTENT(OUT)
     &     TRANS(NPD_PROFILE, NPD_LAYER)                                &
!             DIFFUSE TRANSMISSION COEFFICIENT
     &   , REFLECT(NPD_PROFILE, NPD_LAYER)                              &
!             DIFFUSE REFLECTION COEFFICIENT
     &   , TRANS_0(NPD_PROFILE, NPD_LAYER)                              &
!             DIRECT TRANSMISSION COEFFICIENT
     &   , SOURCE_COEFF(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)
!             SOURCE COEFFICIENTS
!
!
!     LOCAL VARIABLES
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
      REAL                                                              &
     &     GAMMA                                                        &
!             GAMMA
     &   , EXPONENTIAL
!             EXPONENTIAL OF SCALED OPTICAL DEPTH
!     REAL  XLAMTAU(N_PROFILE,I_LAYER_LAST-I_LAYER_FIRST+1) !Workspace
!     INTEGER n_input     ! No. of inputs for exp_v
      REAL TMP_INV, TEMP_NUM2,TEMP_DEN2
      REAL TEMP_NUM1,TEMP_DEN1
      REAL TEMP(NPD_PROFILE)
      real ttt
!
!     Variables related to the treatment of ill-conditioning
      REAL ::                                                           &
     &    TOL
!           The tolerance used for switching to the asymptotic form 
!           of the quadratic source function term.
!
!
!
!
!     Set the tolerances used in avoiding ill-conditioning.
      IF (lrad_quad_src_fix) THEN
        TOL=SQRT(SQRT(EPSILON(TOL)))
      ELSE
        TOL=EXP(3.3E-01*LOG(EPSILON(TOL)))
      ENDIF
!
!     DETERMINE THE DIFFUSE TRANSMISSION AND REFLECTION COEFFICIENTS.
!
!     DO I=I_LAYER_FIRST, I_LAYER_LAST
!        DO L=1, N_PROFILE
!           XLAMTAU(L,I-I_LAYER_FIRST+1)=-LAMBDA(L,I)*TAU(L,I)
!        ENDDO
!     ENDDO
!     n_input=(I_LAYER_LAST-I_LAYER_FIRST+1)*N_PROFILE
!     call exp_v(n_input,xlamtau,xlamtau)
!
      DO I=I_LAYER_FIRST, I_LAYER_LAST
         DO L=1, N_PROFILE
            EXPONENTIAL=exp(-LAMBDA(L,I)*TAU(L,I))
            GAMMA=(SUM(L, I)-LAMBDA(L, I))                              &
     &         /(SUM(L, I)+LAMBDA(L, I))
            TMP_INV = 1.0E+00                                           &
     &         /(1.0E+00-(EXPONENTIAL*GAMMA)**2)
            TRANS(L, I)=EXPONENTIAL*(1.0E+00-GAMMA**2)                  &
     &         *TMP_INV
            REFLECT(L, I)=GAMMA*(1.0E+00-EXPONENTIAL**2)                &
     &         *TMP_INV
         ENDDO
      ENDDO
!
!
!
      IF (ISOLIR == IP_SOLAR) THEN
!
!        CALCULATE THE DIRECT TRANSMISSION AND THE SOURCE COEFFICIENTS
!        FOR THE SOLAR BEAM: IN THE SOLAR CASE THESE ARE
!        THE COEFFICIENTS WHICH WILL MULTIPLY THE DIRECT FLUX AT THE
!        TOP OF THE LAYER TO GIVE THE SOURCE TERMS FOR THE UPWARD
!        DIFFUSE FLUX AND THE TOTAL DOWNWARD FLUX.
!
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
              TRANS_0(L, I)=EXP(-TAU(L, I)*SEC_0(L))
            ENDDO
            DO L=1,N_PROFILE
               SOURCE_COEFF(L,I,IP_SCF_SOLAR_UP)                        &
     &              = (GAMMA_UP(L,I) - REFLECT(L,I)                     &
     &              *(1.0E+00+GAMMA_DOWN(L,I)))
               SOURCE_COEFF(L,I,IP_SCF_SOLAR_DOWN)                      &
     &              =(1.0E+00+GAMMA_DOWN(L,I)                           &
     &              -GAMMA_UP(L,I)*REFLECT(L,I))
            END DO
            DO L=1,N_PROFILE
               SOURCE_COEFF(L,I,IP_SCF_SOLAR_UP)                        &
     &              = SOURCE_COEFF(L,I,IP_SCF_SOLAR_UP)                 &
     &              -GAMMA_UP(L,I)*TRANS(L,I)*TRANS_0(L,I)
            END DO
            DO L=1,N_PROFILE
               SOURCE_COEFF(L,I,IP_SCF_SOLAR_DOWN)                      &
     &              = TRANS_0(L,I)*SOURCE_COEFF(L,I,IP_SCF_SOLAR_DOWN)  &
     &              -(1.0E+00+GAMMA_DOWN(L,I))*TRANS(L,I)
            END DO
         ENDDO
!
!
      ELSE IF (ISOLIR == IP_INFRA_RED) THEN
!
!        IN THE CASE OF INFRA-RED RADIATION, THE FIRST SOURCE
!        COEFFICIENT HOLDS THE MULTIPLIER FOR THE FIRST DIFFERENCE
!        OF THE PLANCK FUNCTION ACROSS THE LAYER, AND THE SECOND
!        THAT FOR THE SECOND DIFFERENCE.
!
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE

!              A TOLERANCE IS ADDED TO THE NUMERATOR AND THE DENOMIATOR
!              TO AVOID ILL-CONDITIONING AT SMALL OPTICAL DEPTHS.
!
               SOURCE_COEFF(L, I, IP_SCF_IR_1D)=(1.0E+00-TRANS(L, I)    &
     &            +REFLECT(L, I)+SQRT(EPSILON(TRANS)))                  &
     &            /(SQRT(EPSILON(TRANS))+TAU(L, I)*SUM(L, I))

            ENDDO
         ENDDO
!
!
         IF (L_IR_SOURCE_QUAD) THEN
!
!           QUADRATIC CORRECTION TO SOURCE FUNCTION.
!           THIS CORRECTION IS VERY ILL-CONDITIONED FOR
!           SMALL OPTICAL DEPTHS SO THE ASYMPTOTIC FORM IS THEN USED.
!
            DO I=I_LAYER_FIRST, I_LAYER_LAST
               DO L=1, N_PROFILE
!
!                 USE A SEPARATE ASYMPTOTIC FORM WHEN THE OPTICAL
!                 DEPTH IS SMALL, MAKING THE TRANSITION WHEN THE
!                 OPTICAL DEPTH IS ROUGHLY EQUAL TO THE CUBE ROOT
!                 OF THE MACHINE'S PRECISION.
!
                  IF (TAU(L, I) > TOL) THEN
                     ttt                                                &
     &                  =-2.0E+00*(1.0E+00-TRANS(L, I)-REFLECT(L, I)    &
     &                  +SQRT(EPSILON(TRANS)))                          &
     &                  /(DIFF(L, I)*TAU(L, I)+SQRT(EPSILON(TRANS)))
                  ELSE
                     ttt                                                &
     &                  =-2.0E+00+DIFF(L, I)*TAU(L, I)
                  ENDIF
!
                  SOURCE_COEFF(L, I, IP_SCF_IR_2D)                      &
     &               =-(1.0E+00+REFLECT(L, I)+TRANS(L, I)               &
     &               +ttt)                                              &
     &               /(SUM(L, I)*TAU(L, I)+SQRT(EPSILON(TRANS)))
               ENDDO
            ENDDO
!
         ENDIF
!
      ENDIF
!
      RETURN
      END SUBROUTINE TRANS_SOURCE_COEFF
