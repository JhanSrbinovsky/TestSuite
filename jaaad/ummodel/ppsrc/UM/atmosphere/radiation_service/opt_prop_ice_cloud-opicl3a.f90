


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate optical properties of ice clouds.
!
! Method:
!       If the optical properties come from an observational
!       distribution a separte subroutine is called. Otherwise
!       appropriate mean quantities in the layer are calculated
!       as the parametrization requires and these values are
!       substituted into the parametrization to give the optical
!       properties.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.4             30-09-96                Old scheme for
!                                               Cirrus removed.
!                                               New scheme based on
!                                               anomalous diffraction
!                                               theory introduced.
!                                               Effective radius is
!                                               relabelled as the
!                                               characteristic
!                                               dimension for more
!                                               general formulations.
!                                               (J. M. Edwards)
!       5.5             24-02-03                Code for aggregate
!                                               parametrization of
!                                               ice crystals included.
!                                               (J. M. Edwards)
!  6.0  21/08/03  NEC SX-6 optimisation - R Barnes & J-C Rioual.
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE OPT_PROP_ICE_CLOUD(IERR                                &
     &   , N_PROFILE, N_LAYER, N_CLOUD_TOP                              &
     &   , N_CLOUD_PROFILE, I_CLOUD_PROFILE                             &
     &   , L_RESCALE, L_LAYER, L_CLOUD_LAYER                            &
     &   , I_PARAMETRIZATION_ICE, ICE_CLOUD_PARAMETER                   &
     &   , ICE_MASS_FRAC, DIM_CHAR_ICE                                  &
     &   , T, DENSITY                                                   &
     &   , K_EXT_TOT_CLOUD, K_EXT_SCAT_CLOUD                            &
     &   , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD                       &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   , NPD_CLOUD_PARAMETER                                          &
     &   )
!
!
      IMPLICIT NONE
!
!
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER                                                    &
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_CLOUD_PARAMETER
!             MAXIMUM NUMBER OF CLOUD PARAMETERS
!
!     INCLUDE COMDECKS
! STDIO3A defines unit numbers for standard i/o in two-stream radiation
! code.
      INTEGER,PARAMETER:: IU_STDIN=5
      INTEGER,PARAMETER:: IU_STDOUT=6
      INTEGER,PARAMETER:: IU_ERR=6
! STDIO3A end
! ICLPRM3A defines numbers for ice cloud schemes in two-stream radiation
! code.
!
!   5.5   Feb 2003     Addition of the aggregate parametrization.
!                                                 John Edwards
!   6.2   Jan 2006     Various options for radiation code
!                      3Z added.   (J.-C. Thelen)
!

      ! number of cloud fitting schemes
      INTEGER,PARAMETER:: NPD_ICE_CLOUD_FIT=6

      ! parametrization of slingo and schrecker.
      INTEGER,PARAMETER:: IP_SLINGO_SCHRECKER_ICE=1

      ! unparametrized ice crystal data
       INTEGER,PARAMETER:: IP_ICE_UNPARAMETRIZED=3

      ! sun and shine's parametrization in the visible (version 2)
      INTEGER,PARAMETER:: IP_SUN_SHINE_VN2_VIS=4

      ! sun and shine's parametrization in the ir (version 2)
      INTEGER,PARAMETER:: IP_SUN_SHINE_VN2_IR=5

      ! scheme based on anomalous diffraction theory for ice crystals
      INTEGER,PARAMETER:: IP_ICE_ADT=6

      ! Provisional agregate parametrization.
      INTEGER,PARAMETER:: IP_ICE_AGG_DE=12

! ICLPRM3A end
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
!     DUMMY VARIABLES.
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , N_CLOUD_TOP                                                  &
!             TOPMOST CLOUDY LAYER
     &   , I_PARAMETRIZATION_ICE                                        &
!             TREATMENT OF ICE CRYSTALS
     &   , N_CLOUD_PROFILE(NPD_LAYER)                                   &
!             NUMBER OF CLOUDY PROFILES
     &   , I_CLOUD_PROFILE(NPD_PROFILE, NPD_LAYER)
!             PROFILES CONTAINING CLOUDS
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_LAYER                                                      &
!             VARIABLES GIVEN IN LAYERS
     &   , L_CLOUD_LAYER                                                &
!             CLOUD VARIABLES GIVEN IN LAYERS
     &   , L_RESCALE
!             DELTA-RESCALING REQUIRED
      REAL                                                              &
                !, INTENT(IN)
     &     ICE_CLOUD_PARAMETER(NPD_CLOUD_PARAMETER)                     &
!             ICE CLOUD PARAMETERS
     &   , ICE_MASS_FRAC(NPD_PROFILE, 0: NPD_LAYER)                     &
!             ICE MASS FRACTION
     &   , DIM_CHAR_ICE(NPD_PROFILE, 0: NPD_LAYER)                      &
!             CHARACTERISTIC DIMENSION FOR CRYSTALS
     &   , T(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             TEMPERATURE
     &   , DENSITY(NPD_PROFILE, 0: NPD_LAYER)
!             DENSITY AT LEVELS
      REAL                                                              &
                !, INTENT(OUT)
     &     K_EXT_SCAT_CLOUD(NPD_PROFILE, NPD_LAYER)                     &
!             SCATTERING EXTINCTION
     &   , K_EXT_TOT_CLOUD(NPD_PROFILE, NPD_LAYER)                      &
!             TOTAL EXTINCTION
     &   , ASYMMETRY_CLOUD(NPD_PROFILE, NPD_LAYER)                      &
!             CLOUDY ASYMMETRIES
     &   , FORWARD_SCATTER_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY FORWARD SCATTERING
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     L                                                            &
!             LOOP VARIABLE
     &   , LL                                                           &
!             LOOP VARIABLE
     &   , I
!             LOOP VARIABLE
      REAL                                                              &
     &     ASYMMETRY_PROCESS(NPD_PROFILE)                               &
!             ASYMMETRY FACTOR FOR CURRENT PROC.
     &   , DIM_CHAR_AVE                                                 &
!             AVERAGE CHARACTERISTIC DIMENSION IN LAYER
     &   , X                                                            &
!             TEMPORARY ALGEBRAIC VARIABLE
     &   , ICE_MASS_FRAC_AVE                                            &
!             AVERAGE ICE MASS FRACTION
     &   , DENSITY_AVE                                                  &
!             AVERAGE DENSITY
     &   , T_CELSIUS                                                    &
!             TEMPERATURE IN CELSIUS
     &   , TEMP_CORRECTION
!             TEMPERATURE CORRECTION
!
!
!
!
      IF (I_PARAMETRIZATION_ICE == IP_SLINGO_SCHRECKER_ICE) THEN
!
!        WE CALCULATE AVERAGE PROPERTIES FOR THE LAYER AND PUT THESE
!        INTO THE PARAMETRIZATION, RATHER THAN CALCULATING THE
!        PARAMETRIZATION AT EACH LEVEL: USUALLY THIS IS MORE ACCURATE.
!        IT ALSO FITS MORE NATURALLY WITH CASES WHERE DATA ARE GIVEN
!        IN LAYERS.
!
!        THE TOTAL EXTINCTIONS SHOULD BE INCREMENTED BY THE TOTAL
!        CONTRIBUTIONS FROM CLOUDS, NOT JUST BY THE ABSORPTIVE
!        EXTINCTIONS.
!
         DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
            DO LL=1, N_CLOUD_PROFILE(I)
               L=I_CLOUD_PROFILE(LL, I)
               ICE_MASS_FRAC_AVE=ICE_MASS_FRAC(L, I)
               DIM_CHAR_AVE=DIM_CHAR_ICE(L, I)
               K_EXT_TOT_CLOUD(L, I)                                    &
     &            =ICE_MASS_FRAC_AVE*(ICE_CLOUD_PARAMETER(1)            &
     &            +ICE_CLOUD_PARAMETER(2)/DIM_CHAR_AVE)
               K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)             &
     &            *(1.0E+00-ICE_CLOUD_PARAMETER(3)                      &
     &            -ICE_CLOUD_PARAMETER(4)*DIM_CHAR_AVE)
               ASYMMETRY_PROCESS(L)=                                    &
     &            ICE_CLOUD_PARAMETER(5)+ICE_CLOUD_PARAMETER(6)         &
     &            *DIM_CHAR_AVE
               ASYMMETRY_CLOUD(L, I)=                                   &
     &            K_EXT_SCAT_CLOUD(L, I)*ASYMMETRY_PROCESS(L)
            ENDDO
            IF (L_RESCALE) THEN
!CDIR NODEP
               DO LL=1, N_CLOUD_PROFILE(I)
                  L=I_CLOUD_PROFILE(LL, I)
                  FORWARD_SCATTER_CLOUD(L, I)                           &
     &               =K_EXT_SCAT_CLOUD(L, I)                            &
     &               *ASYMMETRY_PROCESS(L)**2
               ENDDO
            ENDIF
         ENDDO
!
      ELSE IF (I_PARAMETRIZATION_ICE == IP_ICE_ADT) THEN
!
         DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
            DO LL=1, N_CLOUD_PROFILE(I)
               L=I_CLOUD_PROFILE(LL, I)
               ICE_MASS_FRAC_AVE=ICE_MASS_FRAC(L, I)
               DIM_CHAR_AVE=DIM_CHAR_ICE(L, I)
               X=LOG(DIM_CHAR_AVE/ICE_CLOUD_PARAMETER(10))
               IF (X >  0.0E+00) THEN
!                 LARGE MODE.
                  K_EXT_TOT_CLOUD(L, I)=ICE_MASS_FRAC_AVE               &
     &               *EXP(ICE_CLOUD_PARAMETER(1)                        &
     &               +X*(ICE_CLOUD_PARAMETER(2)                         &
     &               +X*(ICE_CLOUD_PARAMETER(3)                         &
     &               +X*(ICE_CLOUD_PARAMETER(4)                         &
     &               +X*ICE_CLOUD_PARAMETER(5)))))
               ELSE
!                 SMALL MODE.
                  K_EXT_TOT_CLOUD(L, I)=ICE_MASS_FRAC_AVE               &
     &               *EXP(ICE_CLOUD_PARAMETER(1)                        &
     &               +X*(ICE_CLOUD_PARAMETER(6)                         &
     &               +X*(ICE_CLOUD_PARAMETER(7)                         &
     &               +X*(ICE_CLOUD_PARAMETER(8)                         &
     &               +X*ICE_CLOUD_PARAMETER(9)))))
               ENDIF
               X=LOG(DIM_CHAR_AVE/ICE_CLOUD_PARAMETER(20))
               IF (X >  0.0E+00) THEN
!                 LARGE MODE.
                  K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)          &
     &               *(1.0E+00-(ICE_CLOUD_PARAMETER(11)                 &
     &               +X*(ICE_CLOUD_PARAMETER(12)                        &
     &               +X*(ICE_CLOUD_PARAMETER(13)                        &
     &               +X*(ICE_CLOUD_PARAMETER(14)                        &
     &               +X*ICE_CLOUD_PARAMETER(15))))))
               ELSE
!                 SMALL MODE.
                  K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)          &
     &               *(1.0E+00-(ICE_CLOUD_PARAMETER(11)                 &
     &               +X*(ICE_CLOUD_PARAMETER(16)                        &
     &               +X*(ICE_CLOUD_PARAMETER(17)                        &
     &               +X*(ICE_CLOUD_PARAMETER(18)                        &
     &               +X*ICE_CLOUD_PARAMETER(19))))))
               ENDIF
               X=LOG(DIM_CHAR_AVE/ICE_CLOUD_PARAMETER(30))
               IF (X >  0.0E+00) THEN
!                 LARGE MODE.
                  ASYMMETRY_PROCESS(L)=ICE_CLOUD_PARAMETER(21)          &
     &               +X*(ICE_CLOUD_PARAMETER(22)                        &
     &               +X*(ICE_CLOUD_PARAMETER(23)                        &
     &               +X*(ICE_CLOUD_PARAMETER(24)                        &
     &               +X*ICE_CLOUD_PARAMETER(25))))
               ELSE
!                 SMALL MODE.
                  ASYMMETRY_PROCESS(L)=ICE_CLOUD_PARAMETER(21)          &
     &               +X*(ICE_CLOUD_PARAMETER(26)                        &
     &               +X*(ICE_CLOUD_PARAMETER(27)                        &
     &               +X*(ICE_CLOUD_PARAMETER(28)                        &
     &               +X*ICE_CLOUD_PARAMETER(29))))
               ENDIF
               ASYMMETRY_CLOUD(L, I)=                                   &
     &            K_EXT_SCAT_CLOUD(L, I)*ASYMMETRY_PROCESS(L)
            ENDDO
            IF (L_RESCALE) THEN
!CDIR NODEP
               DO LL=1, N_CLOUD_PROFILE(I)
                  L=I_CLOUD_PROFILE(LL, I)
                  FORWARD_SCATTER_CLOUD(L, I)                           &
     &               =K_EXT_SCAT_CLOUD(L, I)                            &
     &               *ASYMMETRY_PROCESS(L)**2
               ENDDO
            ENDIF
         ENDDO
!
!
      ELSE IF (I_PARAMETRIZATION_ICE == IP_ICE_AGG_DE) THEN
!
        DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
          DO LL=1, N_CLOUD_PROFILE(I)
            L=I_CLOUD_PROFILE(LL, I)
            X=DIM_CHAR_ICE(L, I)/ICE_CLOUD_PARAMETER(4)
            K_EXT_TOT_CLOUD(L, I)=ICE_MASS_FRAC(L, I)                   &
     &        *((ICE_CLOUD_PARAMETER(3)/X                               &
     &        +ICE_CLOUD_PARAMETER(2))/X                                &
     &        +ICE_CLOUD_PARAMETER(1))
            X=DIM_CHAR_ICE(L, I)/ICE_CLOUD_PARAMETER(9)
            K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)                &
     &        *(1.0                                                     &
     &        -(ICE_CLOUD_PARAMETER(5)+X                                &
     &        *(ICE_CLOUD_PARAMETER(6)+X                                &
     &        *(ICE_CLOUD_PARAMETER(7)+X                                &
     &        *ICE_CLOUD_PARAMETER(8)))))
            X=DIM_CHAR_ICE(L, I)/ICE_CLOUD_PARAMETER(14)
            ASYMMETRY_PROCESS(L)=ICE_CLOUD_PARAMETER(10)                &
     &        +X*(ICE_CLOUD_PARAMETER(11)                               &
     &        +X*(ICE_CLOUD_PARAMETER(12)                               &
     &        +X*ICE_CLOUD_PARAMETER(13)))
            ASYMMETRY_CLOUD(L, I)                                       &
     &        =K_EXT_SCAT_CLOUD(L, I)*ASYMMETRY_PROCESS(L)
          ENDDO
          IF (L_RESCALE) THEN
!CDIR NODEP
            DO LL=1, N_CLOUD_PROFILE(I)
              L=I_CLOUD_PROFILE(LL, I)
              FORWARD_SCATTER_CLOUD(L, I)                               &
     &          =K_EXT_SCAT_CLOUD(L, I)                                 &
     &          *ASYMMETRY_PROCESS(L)**2
            ENDDO
          ENDIF
        ENDDO
!
      ELSE IF (I_PARAMETRIZATION_ICE == IP_SUN_SHINE_VN2_VIS) THEN
!
         DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
            DO LL=1, N_CLOUD_PROFILE(I)
               L=I_CLOUD_PROFILE(LL, I)
               ICE_MASS_FRAC_AVE=ICE_MASS_FRAC(L, I)
               DENSITY_AVE=DENSITY(L, I)
               T_CELSIUS=T(L, I)-2.7316E+02
               TEMP_CORRECTION=1.047E+00+T_CELSIUS*(-9.13E-05+T_CELSIUS &
     &            *(2.026E-04-1.056E-05*T_CELSIUS))
               K_EXT_TOT_CLOUD(L, I)=TEMP_CORRECTION*ICE_MASS_FRAC_AVE  &
     &            /(3.05548E-02                                         &
     &            +2.54802E+02*DENSITY_AVE*ICE_MASS_FRAC_AVE)
               K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)             &
     &            *(1.0E+00-ICE_CLOUD_PARAMETER(1)                      &
     &            *EXP(ICE_CLOUD_PARAMETER(2)                           &
     &            *LOG(DENSITY_AVE*ICE_MASS_FRAC_AVE+1.0E-12)))         &
     &            *(1.0E+00+ICE_CLOUD_PARAMETER(5)                      &
     &            *(TEMP_CORRECTION-1.0E+00)/TEMP_CORRECTION)
               ASYMMETRY_PROCESS(L)=                                    &
     &            ICE_CLOUD_PARAMETER(3)*EXP(ICE_CLOUD_PARAMETER(4)     &
     &            *LOG(DENSITY_AVE*ICE_MASS_FRAC_AVE+1.0E-12))          &
     &            *(1.0E+00+ICE_CLOUD_PARAMETER(6)                      &
     &            *(TEMP_CORRECTION-1.0E+00)/TEMP_CORRECTION)
               ASYMMETRY_CLOUD(L, I)=                                   &
     &            K_EXT_SCAT_CLOUD(L, I)*ASYMMETRY_PROCESS(L)
            ENDDO
            IF (L_RESCALE) THEN
!CDIR NODEP
               DO LL=1, N_CLOUD_PROFILE(I)
                  L=I_CLOUD_PROFILE(LL, I)
                  FORWARD_SCATTER_CLOUD(L, I)                           &
     &               =K_EXT_SCAT_CLOUD(L, I)                            &
     &               *ASYMMETRY_PROCESS(L)**2
               ENDDO
            ENDIF
         ENDDO
!
      ELSE IF (I_PARAMETRIZATION_ICE == IP_SUN_SHINE_VN2_IR) THEN
!
         DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
            DO LL=1, N_CLOUD_PROFILE(I)
               L=I_CLOUD_PROFILE(LL, I)
               ICE_MASS_FRAC_AVE=ICE_MASS_FRAC(L, I)
               DENSITY_AVE=DENSITY(L, I)
               T_CELSIUS=T(L, I)-2.7316E+02
               TEMP_CORRECTION=1.047E+00+T_CELSIUS*(-9.13E-05+T_CELSIUS &
     &            *(2.026E-04-1.056E-05*T_CELSIUS))
               K_EXT_TOT_CLOUD(L, I)=TEMP_CORRECTION*ICE_MASS_FRAC_AVE  &
     &            /(6.30689E-02                                         &
     &            +2.65874E+02*DENSITY_AVE*ICE_MASS_FRAC_AVE)
            ENDDO
         ENDDO
!
      ELSE
         WRITE(IU_ERR, '(/A)') '*** ERROR: AN INVALID PARAMETRIZATION ' &
     &      //'OF ICE CRYSTALS HAS BEEN USED..'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE OPT_PROP_ICE_CLOUD
