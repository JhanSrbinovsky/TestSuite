

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the mixing ratios of gases.
!
! Purpose:
!   The full array of mass mixing ratios of gases is filled.
!
! Method:
!   The arrays of supplied mixing ratios are inverted and fed
!   into the array to pass to the radiation code. For well-mixed
!   gases the constant mixing ratios are fed into this array.
!
! Current Owner of Code: James Manners
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_GAS_MIX_RATIO(IERR                              &
     &   , N_PROFILE, NLEVS, N_LAYER, NWET, NOZONE                      &
     &   , I_GATHER, L_EXTRA_TOP                                        &
     &   , N_ABSORB, TYPE_ABSORB                                        &
     &   , L_N2O, L_CH4, L_CFC11, L_CFC12, L_O2                         &
     &   , L_CFC113, L_HCFC22, L_HFC125, L_HFC134A                      &
     &   , H2O, CO2, O3, N2O_MIX_RATIO, CH4_MIX_RATIO                   &
     &   , C11_MIX_RATIO, C12_MIX_RATIO, O2_MIX_RATIO                   &
     &   , C113_MIX_RATIO, HCFC22_MIX_RATIO                             &
     &   , HFC125_MIX_RATIO, HFC134A_MIX_RATIO                          &
     &   , GAS_MIX_RATIO                                                &
     &   , CO2_DIM1, CO2_DIM2, CO2_3D, L_CO2_3D                         &
     &   , stochem_CH4, L_use_stochem_CH4                               &
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_SPECIES               &
! chemical greenhouse gas fields
     &   , ngrgas, grgas_field                                          &
     &   )
!
!
!     COMDECKS INCLUDED
! GASID3A defines indexing numbers of gaseous absorbing species for
! two-stream radiation code.
! the numbering 1-12 corresponds to lowtran 7.

      INTEGER,PARAMETER:: NPD_GASES=19 ! Number of indexed gases

      INTEGER,PARAMETER:: IP_H2O=1
      INTEGER,PARAMETER:: IP_CO2=2
      INTEGER,PARAMETER:: IP_O3=3
      INTEGER,PARAMETER:: IP_N2O=4
      INTEGER,PARAMETER:: IP_CO=5
      INTEGER,PARAMETER:: IP_CH4=6
      INTEGER,PARAMETER:: IP_O2=7
      INTEGER,PARAMETER:: IP_NO=8
      INTEGER,PARAMETER:: IP_SO2=9
      INTEGER,PARAMETER:: IP_NO2=10
      INTEGER,PARAMETER:: IP_NH3=11
      INTEGER,PARAMETER:: IP_HNO3=12
      INTEGER,PARAMETER:: IP_N2=13
      INTEGER,PARAMETER:: IP_CFC11=14
      INTEGER,PARAMETER:: IP_CFC12=15
      INTEGER,PARAMETER:: IP_CFC113=16
      INTEGER,PARAMETER:: IP_HCFC22=17
      INTEGER,PARAMETER:: IP_HFC125=18
      INTEGER,PARAMETER:: IP_HFC134A=19

! GASID3A end
! STDIO3A defines unit numbers for standard i/o in two-stream radiation
! code.
      INTEGER,PARAMETER:: IU_STDIN=5
      INTEGER,PARAMETER:: IU_STDOUT=6
      INTEGER,PARAMETER:: IU_ERR=6
! STDIO3A end
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
! UKCA_FEEDBACK start
!
! Purpose: define positions of individual greenhouse gases in
!          greenhouse gas array, grgas_field.
!
      INTEGER, PARAMETER :: p_o3 = 1, p_ch4 = 2, p_n2o = 3, p_f11 = 4,  &
     &                      p_f12 = 5, p_f113 = 6, p_f22 = 7, p_h2os=8

!
!     DUMMY ARGUMENTS.
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
!
!     SIZES OF ARRAYS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_FIELD                                                    &
!             SIZE OF ARRAY FROM UM
     &   , NPD_PROFILE                                                  &
!             SIZE OF ARRAY
     &   , NPD_LAYER                                                    &
!             SIZE OF ARRAY
     &   , NPD_SPECIES
!             SIZE OF ARRAY
!
!     SIZES USED:
!       5.1             04-04-00                Tolerances replaced by
!                                               F90 intrinsics.
!                                               (J. M. Edwards)
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , NLEVS                                                        &
!             Number of layers in the main model
     &   , N_LAYER                                                      &
!             Number of radiative layers
     &   , NWET                                                         &
!             NUMBER OF WET LEVELS
     &   , NOZONE
!             NUMBER OF OZONE LEVELS
!
!     GATHERING ARRAY:
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_GATHER(NPD_FIELD)
!             LIST OF POINTS TO BE GATHERED
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_EXTRA_TOP
!             Flag to use an extra top layer in radiative calculations
!
!     TYPES OF GASES:
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_ABSORB                                                     &
!             NUMBER OF ABSORBERS
     &   , TYPE_ABSORB(NPD_SPECIES)
!             TYPES OF ABSORBERS
!
!     FLAGS FOR MINOR GASES:
      LOGICAL                                                           &
                !,INTENT(IN)
     &     L_N2O                                                        &
!             FLAG FOR NITROUS OXIDE
     &   , L_CH4                                                        &
!             FLAG FOR METHANE
     &   , L_CFC11                                                      &
!             FLAG FOR CFC11
     &   , L_CFC12                                                      &
!             FLAG FOR CFC12
     &   , L_O2                                                         &
!             FLAG FOR O2
     &   , L_CFC113                                                     &
!             FLAG FOR CFC113
     &   , L_HCFC22                                                     &
!             FLAG FOR HCFC22
     &   , L_HFC125                                                     &
!             FLAG FOR HFC125
     &   , L_HFC134A
!             FLAG FOR HFC134A
!
!     MIXING RATIOS SUPPLIED:
      INTEGER  CO2_DIM1, CO2_DIM2   ! dimensions of CO2_3D field
      LOGICAL  L_CO2_3D    !  controls use of 3D co2 field
      LOGICAL  L_use_stochem_CH4 !  controls use of STOCHEM CH4
      REAL                                                              &
                !, INTENT(IN)
     &     H2O(NPD_FIELD, NWET)                                         &
!             MASS MIXING RATIO OF WATER VAPOUR
     &   , CO2                                                          &
!             MASS MIXING RATIO OF CARBON DIOXIDE
     &   , CO2_3D(CO2_DIM1, CO2_DIM2)                                   &
!             3D MASS MIXING RATIO OF CO2 (full field)
     &   , stochem_CH4(NPD_FIELD,NLEVS)                                 &
!             Mass mixing ratio of CH4 from STOCHEM
     &   , O3(NPD_FIELD, NOZONE)                                        &
!             MASS MIXING RATIO OF OZONE
     &   , N2O_MIX_RATIO                                                &
!             MASS MIXING RATIO OF NITROUS OXIDE
     &   , CH4_MIX_RATIO                                                &
!             MASS MIXING RATIO OF METHANE
     &   , C11_MIX_RATIO                                                &
!             MASS MIXING RATIO OF CFC11
     &   , C12_MIX_RATIO                                                &
!             MASS MIXING RATIO OF CFC12
     &   , O2_MIX_RATIO                                                 &
!             MASS MIXING RATIO OF O2
     &   , C113_MIX_RATIO                                               &
!             MASS MIXING RATIO OF CFC113
     &   , HCFC22_MIX_RATIO                                             &
!             MASS MIXING RATIO OF HCFC22
     &   , HFC125_MIX_RATIO                                             &
!             MASS MIXING RATIO OF HFC125
     &   , HFC134A_MIX_RATIO
!             MASS MIXING RATIO OF HFC134A

!  ngrgas is either non-zero, if called from LWRAD, or 0
!  if called from swrad. In the latter case, the supplied
!  fields are ignored.
      INTEGER, INTENT(IN) :: ngrgas
      REAL, INTENT(IN) :: grgas_field(npd_field, nlevs, ngrgas) 
!
!     ARRAY OF ASSIGNED MXING RATIOS:
      REAL                                                              &
                !, INTENT(OUT)
     &     GAS_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER, NPD_SPECIES)
!             MIXING RATIOS
!
!     LOCAL VARIABLES.
!
!     POINTERS TO GASES:
      INTEGER                                                           &
     &     IUMP_H2O                                                     &
!             POINTER TO WATER VAPOUR
     &   , IUMP_CO2                                                     &
!             POINTER TO CARBON DIOXIDE
     &   , IUMP_O3                                                      &
!             POINTER TO OZONE
     &   , IUMP_N2O                                                     &
!             POINTER TO NITOUS OXIDE
     &   , IUMP_CH4                                                     &
!             POINTER TO METHANE
     &   , IUMP_CFC11                                                   &
!             POINTER TO CFC11
     &   , IUMP_CFC12                                                   &
!             POINTER TO CFC12
     &   , IUMP_O2                                                      &
!             POINTER TO O2
     &   , IUMP_CFC113                                                  &
!             POINTER TO CFC113
     &   , IUMP_HCFC22                                                  &
!             POINTER TO HCFC22
     &   , IUMP_HFC125                                                  &
!             POINTER TO HFC125
     &   , IUMP_HFC134A
!             POINTER TO HFC134A
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L                                                            &
!             LOOP VARIABLE
     &   , LG                                                           &
!             CORRESPONDING UNGATHERED INDEX
     &   , I_TOP_COPY
!             Topmost layer where properties are set by copying the
!             input fields.
!
!
!
!
!     MATCH THE INDEXING NUMBERS OF GASEOUS SPECIES IN THE SPECTRAL
!     FILE WITH ACTUAL TYPES OF GASES KNOWN TO THE UM.
!
!     SET ALL POINTERS TO 0 INITIALLY TO FLAG MISSING GASES.
      IUMP_H2O=0
      IUMP_CO2=0
      IUMP_O3=0
      IUMP_N2O=0
      IUMP_CH4=0
      IUMP_CFC11=0
      IUMP_CFC12=0
      IUMP_O2=0
      IUMP_CFC113=0
      IUMP_HCFC22=0
      IUMP_HFC125=0
      IUMP_HFC134A=0
!
!
      DO I=1, N_ABSORB
!
         IF (TYPE_ABSORB(I) == IP_H2O) THEN
            IUMP_H2O=I
         ELSE IF (TYPE_ABSORB(I) == IP_CO2) THEN
            IUMP_CO2=I
         ELSE IF (TYPE_ABSORB(I) == IP_O3) THEN
            IUMP_O3=I
         ELSE IF (TYPE_ABSORB(I) == IP_N2O) THEN
            IUMP_N2O=I
         ELSE IF (TYPE_ABSORB(I) == IP_CH4) THEN
            IUMP_CH4=I
         ELSE IF (TYPE_ABSORB(I) == IP_CFC11) THEN
            IUMP_CFC11=I
         ELSE IF (TYPE_ABSORB(I) == IP_CFC12) THEN
            IUMP_CFC12=I
         ELSE IF (TYPE_ABSORB(I) == IP_O2) THEN
            IUMP_O2=I
         ELSE IF (TYPE_ABSORB(I) == IP_CFC113) THEN
            IUMP_CFC113=I
         ELSE IF (TYPE_ABSORB(I) == IP_HCFC22) THEN
            IUMP_HCFC22=I
         ELSE IF (TYPE_ABSORB(I) == IP_HFC125) THEN
            IUMP_HFC125=I
         ELSE IF (TYPE_ABSORB(I) == IP_HFC134A) THEN
            IUMP_HFC134A=I
         ENDIF
!
      ENDDO
!
!
      IF (L_EXTRA_TOP) THEN
!       The second radiative layer will be the first to have properties
!       set by copying input fields.
        I_TOP_COPY=2
      ELSE
!       The first radiative layer will be the first to have properties
!       set by copying input fields.
        I_TOP_COPY=1
      ENDIF
!
!
!     ASSIGN MIXING RATIOS OF THE GASES TO THE MAIN ARRAYS.
!
!     WATER VAPOUR:
!
      IF (IUMP_H2O >  0) THEN
!        No water exists above the wet levels.
         DO I=1, N_LAYER-NWET
            DO L=1, N_PROFILE
               GAS_MIX_RATIO(L, I, IUMP_H2O)=0.0E+00
            ENDDO
         ENDDO
         DO I=N_LAYER-NWET+1, N_LAYER
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               GAS_MIX_RATIO(L, I, IUMP_H2O)                            &
     &           =MAX(H2O(LG, N_LAYER-I+1), 0.0E+00)
            ENDDO
         ENDDO
      ELSE
         WRITE(IU_ERR, '(/A)')                                          &
     &      '*** ERROR: WATER VAPOUR IS NOT IN THE SPECTRAL FILE.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!     CARBON DIOXIDE:
!
      IF (IUMP_CO2 >  0) THEN
         DO I=1, N_LAYER
           IF (L_CO2_3D) THEN
             DO L=1, N_PROFILE
               LG=I_GATHER(L)
               GAS_MIX_RATIO(L, I, IUMP_CO2)=CO2_3D(LG, N_LAYER-I+1)
             ENDDO
           ELSE
             DO L=1, N_PROFILE
               GAS_MIX_RATIO(L, I, IUMP_CO2)=CO2
             ENDDO
           ENDIF
         ENDDO
      ELSE
         WRITE(IU_ERR, '(/A)')                                          &
     &      '*** ERROR: CARBON DIOXIDE IS NOT IN THE SPECTRAL FILE.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!     OZONE:
!
      IF (IUMP_O3 >  0) THEN
!        The input field of ozone is supplied on NOZONE levels.
!        These values apply to the upper layers used by the full UM.
!        If NOZONE is smaller than NLEVS, the mixing ratio on the
!        bottom level supplied is copied to lower levels. If an
!        extra top level is used its mixing ratio is set by copying
!        the value for the top non-radiative level.
         IF (L_EXTRA_TOP) THEN
           DO L=1, N_PROFILE
             LG=I_GATHER(L)
             GAS_MIX_RATIO(L, 1, IUMP_O3)=O3(LG, NOZONE)
           ENDDO
         ENDIF
         DO I=I_TOP_COPY, NOZONE+I_TOP_COPY-1
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               GAS_MIX_RATIO(L, I, IUMP_O3)=O3(LG, NOZONE+I_TOP_COPY-I)
            ENDDO
         ENDDO
         DO I=NOZONE+I_TOP_COPY, N_LAYER
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               GAS_MIX_RATIO(L, I, IUMP_O3)=O3(LG, 1)
            ENDDO
         ENDDO
      ELSE
         WRITE(IU_ERR, '(/A)')                                          &
     &      '*** ERROR: OZONE IS NOT IN THE SPECTRAL FILE.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!
!
!     OTHER TRACE GASES:
!
!     THESE GASES ARE NOT ALWAYS INCLUDED IN THE CALCULATION.
!     TESTING IS THEREFORE MORE INTRICATE.
!
      IF (IUMP_N2O >  0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (L_N2O) THEN
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  IF (ngrgas >= p_n2o) THEN
                    LG=I_GATHER(L)
                    GAS_MIX_RATIO(L, I, IUMP_N2O)=                      &
     &                    grgas_field(LG, N_LAYER-I+1,p_n2o)
                  ELSE
                    gas_mix_ratio(L, I, IUMP_N2O)=N2O_MIX_RATIO
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_N2O)=0.0E+00
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (L_N2O) THEN
            WRITE(IU_ERR, '(/A)')                                       &
     &         '*** ERROR: NITROUS OXIDE IS NOT IN THE SPECTRAL FILE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (IUMP_CH4 >  0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (L_CH4) THEN
           IF (L_use_stochem_CH4) THEN
             DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  LG=I_GATHER(L)
                  GAS_MIX_RATIO(L, I, IUMP_CH4)=                        &
     &                        stochem_CH4(LG,N_LAYER-I+1)
               ENDDO
             ENDDO
           ELSE
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  IF (ngrgas >= p_ch4) THEN
                    LG=I_GATHER(L)
                    GAS_MIX_RATIO(L, I, IUMP_CH4)=                      &
     &                    grgas_field(LG, N_LAYER-I+1,p_ch4)
                  ELSE
                    gas_mix_ratio(L, I, IUMP_CH4)=CH4_MIX_RATIO
                  ENDIF
               ENDDO
            ENDDO
           ENDIF
         ELSE
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_CH4)=0.0E+00
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (L_CH4) THEN
            WRITE(IU_ERR, '(/A)')                                       &
     &         '*** ERROR: METHANE IS NOT IN THE SPECTRAL FILE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (IUMP_CFC11 >  0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (L_CFC11) THEN
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  IF (ngrgas >= p_f11) THEN
                    LG=I_GATHER(L)
                    GAS_MIX_RATIO(L, I, IUMP_CFC11)=                    &
     &                    grgas_field(LG, N_LAYER-I+1,p_f11)
                  ELSE
                    gas_mix_ratio(L, I, IUMP_CFC11)=C11_MIX_RATIO
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_CFC11)=0.0E+00
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (L_CFC11) THEN
            WRITE(IU_ERR, '(/A)')                                       &
     &         '*** ERROR: CFC11 IS NOT IN THE SPECTRAL FILE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (IUMP_CFC12 >  0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (L_CFC12) THEN
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  IF (ngrgas >= p_f12) THEN
                    LG=I_GATHER(L)
                    GAS_MIX_RATIO(L, I, IUMP_CFC12)=                    &
     &                    grgas_field(LG, N_LAYER-I+1,p_f12)
                  ELSE
                    gas_mix_ratio(L, I, IUMP_CFC12)=C12_MIX_RATIO
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_CFC12)=0.0E+00
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (L_CFC12) THEN
            WRITE(IU_ERR, '(/A)')                                       &
     &         '*** ERROR: CFC12 IS NOT IN THE SPECTRAL FILE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (IUMP_O2 >  0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (L_O2) THEN
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_O2)=O2_MIX_RATIO
               ENDDO
            ENDDO
         ELSE
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_O2)=0.0E+00
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (L_O2) THEN
            WRITE(IU_ERR, '(/A)')                                       &
     &         '*** ERROR: O2 IS NOT IN THE SPECTRAL FILE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (IUMP_CFC113 >  0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (L_CFC113) THEN
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  IF (ngrgas >= p_f113) THEN
                    LG=I_GATHER(L)
                    GAS_MIX_RATIO(L, I, IUMP_CFC113)=                   &
     &                    grgas_field(LG, N_LAYER-I+1,p_f113)
                  ELSE
                    gas_mix_ratio(L, I, IUMP_CFC113)=C113_MIX_RATIO
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_CFC113)=0.0E+00
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (L_CFC113) THEN
            WRITE(IU_ERR, '(/A)')                                       &
     &         '*** ERROR: CFC113 IS NOT IN THE SPECTRAL FILE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (IUMP_HCFC22 >  0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (L_HCFC22) THEN
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  IF (ngrgas >= p_f22) THEN
                    LG=I_GATHER(L)
                    GAS_MIX_RATIO(L, I, IUMP_HCFC22)=                   &
     &                    grgas_field(LG, N_LAYER-I+1,p_f22)
                  ELSE
                    gas_mix_ratio(L, I, IUMP_HCFC22)=HCFC22_MIX_RATIO
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_HCFC22)=0.0E+00
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (L_HCFC22) THEN
            WRITE(IU_ERR, '(/A)')                                       &
     &         '*** ERROR: HCFC22 IS NOT IN THE SPECTRAL FILE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (IUMP_HFC125 >  0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (L_HFC125) THEN
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_HFC125)=HFC125_MIX_RATIO
               ENDDO
            ENDDO
         ELSE
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_HFC125)=0.0E+00
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (L_HFC125) THEN
            WRITE(IU_ERR, '(/A)')                                       &
     &         '*** ERROR: HFC125 IS NOT IN THE SPECTRAL FILE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
      IF (IUMP_HFC134A >  0) THEN
!        THE GAS IS IN THE SPECTRAL FILE. IF IT HAS BEEN SELECTED
!        FROM THE UI ITS MIXING RATIO MUST BE SET. IF IT IS IN THE
!        FILE BUT NOT SELECTED THE MIXING RATIO IS SET TO 0.
         IF (L_HFC134A) THEN
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_HFC134A)=HFC134A_MIX_RATIO
               ENDDO
            ENDDO
         ELSE
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  GAS_MIX_RATIO(L, I, IUMP_HFC134A)=0.0E+00
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        THE GAS IS NOT IN THE SPECTRAL FILE. AN ERROR RESULTS IF
!        IT WAS TO BE INCLUDED IN THE CALCULATION.
         IF (L_HFC134A) THEN
            WRITE(IU_ERR, '(/A)')                                       &
     &         '*** ERROR: HFC134A IS NOT IN THE SPECTRAL FILE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE R2_SET_GAS_MIX_RATIO
