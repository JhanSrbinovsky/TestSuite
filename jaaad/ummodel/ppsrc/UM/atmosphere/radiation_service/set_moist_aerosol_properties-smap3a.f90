


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set moist aerosol properties independent of bands.
!
! Method:
!       The mean relative humidities are calculated and pointers to
!       the lookup tables are set.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.3             17-12-96                Code extended to permit
!                                               use with both moist
!                                               and dry aerosols.
!                                               (J. M. Edwards)
!       6.2             23-11-05                Added possibility to
!                                               use the clear-sky mean
!                                               relative humidity.
!                                               (N. Bellouin)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SET_MOIST_AEROSOL_PROPERTIES(IERR                      &
     &   , N_PROFILE, N_LAYER, L_LAYER                                  &
     &   , L_USE_CLEARRH                                                &
     &   , N_AEROSOL, I_AEROSOL_PARAMETRIZATION, NHUMIDITY              &
     &   , WATER_MIX_RATIO, T, P, W_CLOUD                               &
     &   , DELTA_HUMIDITY, MEAN_REL_HUMIDITY, I_HUMIDITY_POINTER        &
     &   , NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES                  &
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
     &   , NPD_LAYER                                                    &
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_AEROSOL_SPECIES
!             MAXIMUM NUMBER OF AEROSOLS
!
!     INCLUDE COMDECKS.
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
! AERPRM3A aerosol parametrizations for two-stream radiation code.

      INTEGER,PARAMETER:: IP_AEROSOL_PARAM_DRY=1
      INTEGER,PARAMETER:: IP_AEROSOL_PARAM_MOIST=2
      INTEGER,PARAMETER:: IP_AEROSOL_UNPARAMETRIZED=3 ! Observational

! AERPRM3A end
! STDIO3A defines unit numbers for standard i/o in two-stream radiation
! code.
      INTEGER,PARAMETER:: IU_STDIN=5
      INTEGER,PARAMETER:: IU_STDOUT=6
      INTEGER,PARAMETER:: IU_ERR=6
! STDIO3A end
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
!
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , N_AEROSOL                                                    &
!             NUMBER OF AEROSOL SPECIES
     &   , I_AEROSOL_PARAMETRIZATION(NPD_AEROSOL_SPECIES)               &
!             PARAMETRIZATIONS OF AEROSOL
!             SPECIES
     &   , NHUMIDITY(NPD_AEROSOL_SPECIES)
!             NUMBER OF HUMIDITY VALUES
      INTEGER                                                           &
                !, INTENT(OUT)
     &     I_HUMIDITY_POINTER(NPD_PROFILE, NPD_LAYER)
!             POINTERS TO LOOK-UP TABLES
      LOGICAL                                                           &
     &     L_USE_CLEARRH
!             control OF THE TYPE OF RELATIVE HUMIDITY
!             TO COMPUTE (CLEAR-SKY MEAN IF TRUE,
!             GRID-BOX MEAN IF FALSE)
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_LAYER
!             LAYER FLAG
      REAL                                                              &
                !, INTENT(IN)
     &     WATER_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER)                   &
!             MIXING RATIO OF WATER VAPOUR
     &   , T(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             TEMPERATURES
     &   , P(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             PRESSURES
     &   , W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUD FRACTION
      REAL                                                              &
                !, INTENT(OUT)
     &     MEAN_REL_HUMIDITY(NPD_PROFILE, NPD_LAYER)                    &
!             MEAN HUMIDITIES OF LAYERS
     &   , DELTA_HUMIDITY
!             INCREMENT IN HUMIDITY
!
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , J                                                            &
!             LOOP VARIABLE
     &   , L                                                            &
!             LOOP VARIABLE
     &   , NHUMIDITY_COMMON
!             COMMON NUMBER OF HUMIDITIES FOR MOIST AEROSOLS
      REAL                                                              &
     &     MIX_RATIO_SAT(NPD_PROFILE, 0: NPD_LAYER)
!             SATURATED HUMIDITY MIXING RATIO
      REAL                                                              &
     &     MIX_RATIO_CLR(NPD_PROFILE, 0: NPD_LAYER)
!             HUMIDITY MIXING RATIO OF THE CLEAR-SKY PART OF
!             THE GRID-BOX
!
!     SUBROUTINES CALLED:
      EXTERNAL                                                          &
     &     QSAT_WAT
!
!
!
!     SET UP ARRAY OF POINTERS TO INCLUDE THE EFFECTS OF HUMIDITY.
!     CALCULATE THE SATURATED MIXING RATIO.
      DO I=1, N_LAYER
! DEPENDS ON: qsat_wat
         CALL QSAT_WAT(MIX_RATIO_SAT(1, I), T(1, I), P(1, I)            &
     &      , N_PROFILE)
      ENDDO
!     CALCULATE THE CLEAR-SKY MIXING RATIO IF NEEDED
      IF (L_USE_CLEARRH) THEN
        DO I=1, N_LAYER
          DO L = 1, N_PROFILE
            MIX_RATIO_CLR(L, I) =                                       &
     &        (WATER_MIX_RATIO(L, I) - (W_CLOUD(L, I) *                 &
     &         MIX_RATIO_SAT(L, I))) /                                  &
     &         MAX( (1.00E+00 - W_CLOUD(L, I)), 0.001)
            ! MAKE SURE THAT 0 <= MIX_RATIO_CLR <= WATER_MIX_RATIO
            MIX_RATIO_CLR(L, I) = MIN(MIX_RATIO_CLR(L, I),              &
     &                                WATER_MIX_RATIO(L, I))
            MIX_RATIO_CLR(L, I) = MAX(MIX_RATIO_CLR(L, I),              &
     &                                0.00E+00)
          ENDDO
        ENDDO
      ENDIF
!
!     DETERMINE THE NUMBER OF HUMIDITIES TO BE USED FOR MOIST
!     AEROSOLS. THIS MUST BE THE SAME FOR ALL MOIST AEROSOLS
!     IN THE CURRENT VERSION OF THE CODE.
      NHUMIDITY_COMMON=0
      DO J=1, N_AEROSOL
         IF (I_AEROSOL_PARAMETRIZATION(J) == IP_AEROSOL_PARAM_MOIST)    &
     &         THEN
            IF (NHUMIDITY(J) >  0) THEN
!              SET THE ACTUAL COMMON VALUE.
               IF (NHUMIDITY_COMMON == 0) THEN
                  NHUMIDITY_COMMON=NHUMIDITY(J)
               ELSE IF (NHUMIDITY(J) /= NHUMIDITY_COMMON) THEN
!                 THERE IS AN INCONSISTENCY.
                  WRITE(IU_ERR, '(/A)')                                 &
     &               '***ERROR: THE LOOK-UP TABLES FOR MOIST AEROSOLS ' &
     &               , 'ARE OF DIFFERENT SIZES. THIS IS NOT PERMITTED.'
                  IERR=I_ERR_FATAL
                  RETURN
               ENDIF
            ENDIF
         ENDIF
      ENDDO
!     THE LOOK-UP TABLE IS ASSUMED TO BE UNIFORM IN HUMIDITY.
      DELTA_HUMIDITY=1.0E+00/(REAL(NHUMIDITY_COMMON)-1.0E+00)

!     COMPUTE THE CLEAR-SKY MEAN RELATIVE HUMIDITY OR THE
!     GRID-BOX MEAN RELATIVE HUMIDITY
      IF(L_USE_CLEARRH) THEN
        DO I = 1, N_LAYER
          DO L = 1, N_PROFILE
            MEAN_REL_HUMIDITY(L,I)                                      &
     &         = MIX_RATIO_CLR(L, I)*(1.0E+00-MIX_RATIO_SAT(L, I))      &
     &         /((1.0E+00-MIX_RATIO_CLR(L, I))*MIX_RATIO_SAT(L, I))
!           CHECK THAT THE CLEAR-SKY MEAN RELATIVE HUMIDITY
!           LIES BETWEEN 0.0 (INCL.) AND 1.0 (EXCL.)
            MEAN_REL_HUMIDITY(L, I)=MIN(MEAN_REL_HUMIDITY(L, I)         &
     &        , 0.99999)
            MEAN_REL_HUMIDITY(L, I)=MAX(MEAN_REL_HUMIDITY(L, I)         &
     &        , 0.00000)
            I_HUMIDITY_POINTER(L, I)=1                                  &
     &         +INT(MEAN_REL_HUMIDITY(L, I)*(NHUMIDITY_COMMON-1))
          ENDDO
        ENDDO
      ELSE
        DO I=1, N_LAYER
          DO L=1, N_PROFILE
            MEAN_REL_HUMIDITY(L, I)                                     &
     &         =WATER_MIX_RATIO(L, I)*(1.0E+00-MIX_RATIO_SAT(L, I))     &
     &         /((1.0E+00-WATER_MIX_RATIO(L, I))*MIX_RATIO_SAT(L, I))
!           CHECK THAT THE MEAN RELATIVE HUMIDITY
!           DOES NOT EXCEED OR EQUAL 1.0.
            MEAN_REL_HUMIDITY(L, I)=MIN(MEAN_REL_HUMIDITY(L, I)         &
     &        , 0.99999)
            I_HUMIDITY_POINTER(L, I)=1                                  &
     &         +INT(MEAN_REL_HUMIDITY(L, I)*(NHUMIDITY_COMMON-1))
          ENDDO
        ENDDO
      ENDIF
!
!
!

      RETURN
      END SUBROUTINE SET_MOIST_AEROSOL_PROPERTIES
