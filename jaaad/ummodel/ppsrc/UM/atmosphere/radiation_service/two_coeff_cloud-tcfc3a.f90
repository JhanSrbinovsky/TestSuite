


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate cloudy two-stream coefficients.
!
! Method:
!       The coeffients for each type of cloud are determined and
!       averaged.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             04-03-96                Gathering and scattering
!                                               of types of clouds
!                                               introduced for speed.
!                                               This change has no
!                                               physical effect.
!       4.2             Nov. 96   T3E migration: CALL WHENFGT replaced
!                                  by portable fortran code.
!                                                S.J.Swarbrick
!LL  4.5  27/04/98  Add Fujitsu vectorization directive.
!LL                                           RBarnes@ecmwf.int
!  6.0  21/08/03  NEC SX-6 optimisation - add vectorisation
!                 !CDIR NODEP directives.  R Barnes & J-C Rioual.
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
! Fujitsu directive to encourage vectorization for whole routine
!OCL NOVREC
      SUBROUTINE TWO_COEFF_CLOUD(IERR                                   &
     &   , N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST                       &
     &   , I_2STREAM, L_IR_SOURCE_QUAD, N_SOURCE_COEFF                  &
     &   , N_CLOUD_TYPE, FRAC_CLOUD                                     &
     &   , ASYMMETRY_CLOUD, OMEGA_CLOUD, TAU_CLOUD                      &
     &   , ISOLIR, SEC_0                                                &
     &   , TRANS_CLOUD, REFLECT_CLOUD, TRANS_0_CLOUD                    &
     &   , SOURCE_COEFF_CLOUD                                           &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
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
!     INCLUDE COMDECKS.
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
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET ERROR FLAGS IN THE RADIATION CODE.
!
      INTEGER                                                           &
     &     I_NORMAL                                                     &
!             ERROR FREE CONDITION
     &   , I_ERR_FATAL                                                  &
!             FATAL ERROR: IMMEDIATE RETURN
     &   , I_ABORT_CALCULATION                                          &
!             CALCULATION ABORTED
     &   , I_MISSING_DATA                                               &
!             MISSING DATA ERROR: CONDITIONAL
     &   , I_ERR_IO                                                     &
!             I/O ERROR
     &   , I_ERR_RANGE                                                  &
!             INTERPOLATION RANGE ERROR
     &   , I_ERR_EXIST
!             EXISTENCE ERROR
!
      PARAMETER(                                                        &
     &     I_NORMAL=0                                                   &
     &   , I_ERR_FATAL=1                                                &
     &   , I_ABORT_CALCULATION=2                                        &
     &   , I_MISSING_DATA=3                                             &
     &   , I_ERR_IO=4                                                   &
     &   , I_ERR_RANGE=5                                                &
     &   , I_ERR_EXIST=6                                                &
     &   )
!
!     ------------------------------------------------------------------
!
!
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , I_LAYER_FIRST                                                &
!             FIRST LAYER TO CONSIDER
     &   , I_LAYER_LAST                                                 &
!             LAST LAYER TO CONSIDER
     &   , ISOLIR                                                       &
!             SPECTRAL REGION
     &   , N_CLOUD_TYPE                                                 &
!             NUMBER OF TYPES OF CLOUDS
     &   , I_2STREAM                                                    &
!             TWO STREAM SCHEME
     &   , N_SOURCE_COEFF
!             NUMBER OF SOURCE COEFFICIENTS
!
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_IR_SOURCE_QUAD
!             USE A QUADRATIC SOURCE IN THE INFRA-RED
!
!     OPTICAL PROPERTIES OF LAYER:
      REAL                                                              &
                !, INTENT(IN)
     &     FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)           &
!             FRACTIONS OF DIFFERENT TYPES OF CLOUDS
     &   , ASYMMETRY_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)      &
!             ASYMMETRY FACTOR
     &   , OMEGA_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)          &
!             ALBEDO OF SINGLE SCATTERING
     &   , TAU_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             OPTICAL DEPTH
!
!     SOLAR BEAM
      REAL                                                              &
                !, INTENT(IN)
     &     SEC_0(NPD_PROFILE)
!             SECANT OF ZENITH ANGLE
!
!
!     COEFFICIENTS IN THE TWO-STREAM EQUATIONS:
      REAL                                                              &
                !, INTENT(OUT)
     &     TRANS_CLOUD(NPD_PROFILE, NPD_LAYER)                          &
!             MEAN DIFFUSE TRANSMISSION COEFFICIENT
     &   , REFLECT_CLOUD(NPD_PROFILE, NPD_LAYER)                        &
!             MEAN DIFFUSE REFLECTION COEFFICIENT
     &   , TRANS_0_CLOUD(NPD_PROFILE, NPD_LAYER)                        &
!             MEAN DIRECT TRANSMISSION COEFFICIENT
     &   , SOURCE_COEFF_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)
!             MEAN SOURCE COEFFICIENTS IN TWO-STREAM EQUATIONS
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , J                                                            &
!             LOOP VARIABLE
     &   , K                                                            &
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
!     COEFFICIENTS IN THE TWO-STREAM EQUATIONS:
      REAL                                                              &
                !, INTENT(OUT)
     &     TRANS_TEMP(NPD_PROFILE, NPD_LAYER)                           &
!             TEMPORARY DIFFUSE TRANSMISSION COEFFICIENT
     &   , REFLECT_TEMP(NPD_PROFILE, NPD_LAYER)                         &
!             TEMPORARY DIFFUSE REFLECTION COEFFICIENT
     &   , TRANS_0_TEMP(NPD_PROFILE, NPD_LAYER)                         &
!             TEMPORARY DIRECT TRANSMISSION COEFFICIENT
     &   , SOURCE_COEFF_TEMP(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)
!             TEMPORARY SOURCE COEFFICIENTS IN TWO-STREAM EQUATIONS
!
!     VARIABLES FOR GATHERING:
      INTEGER                                                           &
     &     N_LIST                                                       &
!             NUMBER OF POINTS IN LIST
     &   , L_LIST(NPD_PROFILE)                                          &
!             LIST OF COLLECTED POINTS
     &   , LL
!             LOOP VARIABLE
      REAL                                                              &
     &     TARGET
!             TARGET FOR SEARCHING
      REAL                                                              &
     &     TAU_GATHERED(NPD_PROFILE, NPD_LAYER)                         &
!             GATHERED OPTICAL DEPTH
     &   , OMEGA_GATHERED(NPD_PROFILE, NPD_LAYER)                       &
!             GATHERED ALEBDO OF SINGLE SCATTERING
     &   , ASYMMETRY_GATHERED(NPD_PROFILE, NPD_LAYER)                   &
!             GATHERED ASYMMETRY
     &   , SEC_0_GATHERED(NPD_PROFILE)
!             GATHERED ASYMMETRY
!
!     SUBROUTINES CALLED:
      EXTERNAL                                                          &
     &     TWO_COEFF
!
!     CRAY DIRECTIVES FOR THE WHOLE ROUTINE:
!     POINTS ARE NOT REPEATED IN THE INDEXING ARRAY, SO IT IS SAFE
!     TO VECTORIZE OVER INDIRECTLY ADDRESSED ARRAYS.
!fpp$ NODEPCHK R
!
!
!
!     INITIALIZE THE FULL ARRAYS.
!
      DO I=I_LAYER_FIRST, I_LAYER_LAST
         DO L=1, N_PROFILE
            TRANS_CLOUD(L, I)=0.0E+00
            REFLECT_CLOUD(L, I)=0.0E+00
         ENDDO
      ENDDO
      DO J=1, N_SOURCE_COEFF
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               SOURCE_COEFF_CLOUD(L, I, J)=0.0E+00
            ENDDO
         ENDDO
      ENDDO
!
      IF (ISOLIR == IP_SOLAR) THEN
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               TRANS_0_CLOUD(L, I)=0.0E+00
            ENDDO
         ENDDO
      ENDIF
!
!
!     CALCULATE THE TRANSMISSION AND REFLECTION COEFFICIENTS FOR
!     EACH TYPE OF CLOUD AND INCREMENT THE TOTALS, WEIGHTING WITH
!     THE CLOUD FRACTION.
!
      DO K=1, N_CLOUD_TYPE
!
!        GATHER POINTS WHERE THERE IS CLOUD OF THE PRESENT TYPE.
         TARGET=0.0E+00
!
         DO I=I_LAYER_FIRST, I_LAYER_LAST
!
!           DETERMINE WHETHER THIS TYPE OF CLOUD EXISTS IN THIS ROW.
!
            N_LIST =0
            DO L   =1,N_PROFILE
              IF (FRAC_CLOUD(L,I,K) >  TARGET) THEN
                N_LIST   =N_LIST+1
                L_LIST(N_LIST)=L
              END IF
            END DO
!
            IF (N_LIST >  0) THEN
!
!              GATHER THE OPTICAL PROPERTIES. THOUGH WE CONSIDER ONLY
!              ONE LAYER AT A TIME THE LOWER ROUTINES WILL OPERATE ON
!              ARRAYS WITH VERTICAL STRUCTURE, SO THE GATHERED ARRAYS
!              ARE TWO-DIMENSIONAL.
!
               DO L=1, N_LIST
                  TAU_GATHERED(L, I)                                    &
     &              =TAU_CLOUD(L_LIST(L), I, K)
                  OMEGA_GATHERED(L, I)                                  &
     &              =OMEGA_CLOUD(L_LIST(L), I, K)
                  ASYMMETRY_GATHERED(L, I)                              &
     &              =ASYMMETRY_CLOUD(L_LIST(L), I, K)
               ENDDO
               IF (ISOLIR == IP_SOLAR) THEN
                  DO L=1, N_LIST
                     SEC_0_GATHERED(L)=SEC_0(L_LIST(L))
                  ENDDO
               ENDIF
!
!
! DEPENDS ON: two_coeff
               CALL TWO_COEFF(IERR                                      &
     &            , N_LIST, I, I                                        &
     &            , I_2STREAM, L_IR_SOURCE_QUAD                         &
     &            , ASYMMETRY_GATHERED, OMEGA_GATHERED                  &
     &            , TAU_GATHERED                                        &
     &            , ISOLIR, SEC_0_GATHERED                              &
     &            , TRANS_TEMP, REFLECT_TEMP, TRANS_0_TEMP              &
     &            , SOURCE_COEFF_TEMP                                   &
     &            , NPD_PROFILE, NPD_LAYER                              &
     &            )
               IF (IERR /= I_NORMAL) RETURN
!
!CDIR NODEP
               DO L=1, N_LIST
                  LL=L_LIST(L)
                  TRANS_CLOUD(LL, I)=TRANS_CLOUD(LL, I)                 &
     &               +FRAC_CLOUD(LL, I, K)*TRANS_TEMP(L, I)
                  REFLECT_CLOUD(LL, I)=REFLECT_CLOUD(LL, I)             &
     &               +FRAC_CLOUD(LL, I, K)*REFLECT_TEMP(L, I)
               ENDDO
               DO J=1, N_SOURCE_COEFF
!CDIR NODEP
                  DO L=1, N_LIST
                     LL=L_LIST(L)
                     SOURCE_COEFF_CLOUD(LL, I, J)                       &
     &                  =SOURCE_COEFF_CLOUD(LL, I, J)                   &
     &                  +FRAC_CLOUD(LL, I, K)*SOURCE_COEFF_TEMP(L, I, J)
                  ENDDO
               ENDDO
               IF (ISOLIR == IP_SOLAR) THEN
!CDIR NODEP
                  DO L=1, N_LIST
                     LL=L_LIST(L)
                     TRANS_0_CLOUD(LL, I)=TRANS_0_CLOUD(LL, I)          &
     &                  +FRAC_CLOUD(LL, I, K)*TRANS_0_TEMP(L, I)
                  ENDDO
               ENDIF
            ENDIF
!
         ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE TWO_COEFF_CLOUD
