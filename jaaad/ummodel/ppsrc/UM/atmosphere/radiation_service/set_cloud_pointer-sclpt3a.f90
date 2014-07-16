


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set pointers to types of clouds
!
! Method:
!       The types of condensate included are examined. Their phases
!       are set and depending on the representation of clouds adopted
!       it is determined to which type of cloud they contribute.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SET_CLOUD_POINTER(IERR                                 &
     &   , N_CONDENSED, TYPE_CONDENSED, I_CLOUD_REPRESENTATION          &
     &   , L_DROP, L_ICE                                                &
     &   , I_PHASE_CMP, I_CLOUD_TYPE, L_CLOUD_CMP                       &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!
!     INCLUDE COMDECKS
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
! PHASE3A defines indices for phases in two-stream radiation code.
      INTEGER,PARAMETER:: IP_PHASE_WATER = 1
      INTEGER,PARAMETER:: IP_PHASE_ICE   = 2
! PHASE3A end
! CLDCMP3A sets components of clouds for two-stream radiation code.

      ! stratiform water droplets
      INTEGER,PARAMETER:: IP_CLCMP_ST_WATER=1

      ! stratiform ice crystals
      INTEGER,PARAMETER:: IP_CLCMP_ST_ICE=2

      ! convective water droplets
      INTEGER,PARAMETER:: IP_CLCMP_CNV_WATER=3

      ! convective ice crystals
      INTEGER,PARAMETER:: IP_CLCMP_CNV_ICE=4

! CLDCMP3A end
! CLREPP3A defines representations of clouds in two-stream radiation
! code.

      ! all components are mixed homogeneously
      INTEGER,PARAMETER:: IP_CLOUD_HOMOGEN     = 1

      ! ice and water clouds are treated separately
      INTEGER,PARAMETER:: IP_CLOUD_ICE_WATER   = 2

      ! clouds are divided into homogeneously mixed stratiform and
      ! convective parts
      INTEGER,PARAMETER:: IP_CLOUD_CONV_STRAT  = 3

      ! clouds divided into ice and water phases and into stratiform and
      ! convective components.
      INTEGER,PARAMETER:: IP_CLOUD_CSIW        = 4

! Types of clouds (values in CLREPD3A)

      ! number of type of clouds in representation
      INTEGER :: NP_CLOUD_TYPE(NPD_CLOUD_REPRESENTATION)

      ! map of components contributing to types of clouds
      INTEGER :: IP_CLOUD_TYPE_MAP(NPD_CLOUD_COMPONENT,                 &
     &  NPD_CLOUD_REPRESENTATION)

! CLREPP3A end
! CLDTYP3A defines cloud types for TWO-STREAM RADIATION CODE.

      INTEGER,PARAMETER:: IP_CLOUD_TYPE_HOMOGEN=1 ! water and ice
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_WATER=1   ! Water only
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_ICE=2     ! Ice only

      ! mixed-phase stratiform cloud
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_STRAT=1

      ! mixed-phase convective cloud
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_CONV=2

      INTEGER,PARAMETER:: IP_CLOUD_TYPE_SW=1 ! stratiform water cloud
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_SI=2 ! stratiform ice cloud
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_CW=3 ! convective water cloud
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_CI=4 ! convective ice cloud

! CLDTYP3A end
!
!     DUMMY VARIABLES.
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_CONDENSED                                                  &
!             NUMBER OF CONDENSED COMPONENTS
     &   , TYPE_CONDENSED(NPD_CLOUD_COMPONENT)                          &
!             TYPES OF COMPONENTS
     &   , I_CLOUD_REPRESENTATION
!             REPRESENTATION OF CLOUDS USED
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_DROP                                                       &
!             FLAG FOR INCLUSION OF DROPLETS
     &   , L_ICE
!             FLAG FOR INCLUSION OF ICE CRYSTALS
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     I_PHASE_CMP(NPD_CLOUD_COMPONENT)                             &
!             PHASES OF COMPONENTS
     &   , I_CLOUD_TYPE(NPD_CLOUD_COMPONENT)
!             TYPES OF CLOUD TO WHICH EACH COMPONENT CONTRIBUTES
      LOGICAL                                                           &
                !, INTENT(OUT)
     &     L_CLOUD_CMP(NPD_CLOUD_COMPONENT)
!             LOGICAL SWITCHES TO INCLUDE COMPONENTS
!
!
!     LOCAL VARIABLES
      INTEGER                                                           &
     &     K
!            LOOP VARIABLE
!
! CLREPD3A defines representations of clouds in two-stream radiation
! code.

      DATA NP_CLOUD_TYPE(IP_CLOUD_HOMOGEN)/1/
      DATA NP_CLOUD_TYPE(IP_CLOUD_ICE_WATER)/2/
      DATA NP_CLOUD_TYPE(IP_CLOUD_CONV_STRAT)/2/
      DATA NP_CLOUD_TYPE(IP_CLOUD_CSIW)/4/

      ! the array ip_cloud_type_map indicates to which type of cloud
      ! each component belongs in a particular representation. an
      ! entry of 0 indicates that that component should not be
      ! present in the representation.

      DATA                                                              &
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_HOMOGEN)        &
     &       /IP_CLOUD_TYPE_HOMOGEN/                                    &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_HOMOGEN)          &
     &       /IP_CLOUD_TYPE_HOMOGEN/                                    &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_HOMOGEN)       &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_HOMOGEN)         &
     &       /0/
      DATA                                                              &
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_ICE_WATER)      &
     &       /IP_CLOUD_TYPE_WATER/                                      &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_ICE_WATER)        &
     &       /IP_CLOUD_TYPE_ICE/                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_ICE_WATER)     &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_ICE_WATER)       &
     &       /0/
      DATA                                                              &
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_CONV_STRAT)     &
     &       /IP_CLOUD_TYPE_STRAT/                                      &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_CONV_STRAT)       &
     &       /IP_CLOUD_TYPE_STRAT/                                      &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_CONV_STRAT)    &
     &       /IP_CLOUD_TYPE_CONV/                                       &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_CONV_STRAT)      &
     &       /IP_CLOUD_TYPE_CONV/
      DATA                                                              &
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_CSIW)           &
     &       /IP_CLOUD_TYPE_SW/                                         &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_CSIW)             &
     &       /IP_CLOUD_TYPE_SI/                                         &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_CSIW)          &
     &       /IP_CLOUD_TYPE_CW/                                         &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_CSIW)            &
     &       /IP_CLOUD_TYPE_CI/

! CLREPD3A end
!
!
!
      DO K=1, N_CONDENSED
!
         I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_MAP(TYPE_CONDENSED(K)            &
     &      , I_CLOUD_REPRESENTATION)
!
!        CHECK FOR 0 FLAGGING ILLEGAL TYPES.
         IF (I_CLOUD_TYPE(K) == 0) THEN
            WRITE(IU_ERR, '(/A)')                                       &
     &         '*** ERROR: A COMPONENT IS NOT COMPATIBLE WITH THE'      &
     &         //'REPRESENTATION OF CLOUDS SELECTED.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
!
         IF (TYPE_CONDENSED(K) == IP_CLCMP_ST_WATER) THEN
!
            I_PHASE_CMP(K)=IP_PHASE_WATER
            L_CLOUD_CMP(K)=L_DROP
!
         ELSE IF (TYPE_CONDENSED(K) == IP_CLCMP_ST_ICE) THEN
!
            I_PHASE_CMP(K)=IP_PHASE_ICE
            L_CLOUD_CMP(K)=L_ICE
!
         ELSE IF (TYPE_CONDENSED(K) == IP_CLCMP_CNV_WATER) THEN
!
            I_PHASE_CMP(K)=IP_PHASE_WATER
            L_CLOUD_CMP(K)=L_DROP
!
         ELSE IF (TYPE_CONDENSED(K) == IP_CLCMP_CNV_ICE) THEN
!
            I_PHASE_CMP(K)=IP_PHASE_ICE
            L_CLOUD_CMP(K)=L_ICE
!
         ENDIF
!
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SET_CLOUD_POINTER
