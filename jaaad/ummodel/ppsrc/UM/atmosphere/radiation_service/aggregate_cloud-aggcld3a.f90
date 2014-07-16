


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to aggregate clouds into regions.
!
! Method:
!       The clouds in a layer are combined in groups to form regions
!       which will be considered as bulk entities in the solution of the
!       equation of transfer. The extents of these regions are also
!       determined.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       HADAM3          05-06-96                Original Code
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE AGGREGATE_CLOUD(IERR                                   &
     &   , N_PROFILE, N_LAYER, N_CLOUD_TOP                              &
     &   , I_CLOUD, I_CLOUD_REPRESENTATION, N_CLOUD_TYPE                &
     &   , FRAC_CLOUD                                                   &
     &   , I_REGION_CLOUD, FRAC_REGION                                  &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     DUMMY ARRAY SIZES
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
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
! CLSCHM3A defines reference numbers for cloud schemes in two-stream
! radiation code.

      ! maximum/random overlap in a mixed column
      INTEGER,PARAMETER:: IP_CLOUD_MIX_MAX=2

      ! random overlap in a mixed column
      INTEGER,PARAMETER:: IP_CLOUD_MIX_RANDOM=4

      ! maximum overlap in a column model
      INTEGER,PARAMETER:: IP_CLOUD_COLUMN_MAX=3

      ! clear column
      INTEGER,PARAMETER:: IP_CLOUD_CLEAR=5

      ! mixed column with split between  convective and layer cloud.
      INTEGER,PARAMETER:: IP_CLOUD_TRIPLE=6

      ! Coupled overlap with partial correlation of cloud
      INTEGER,Parameter:: IP_cloud_part_corr=7

      ! Coupled overlap with partial correlation of cloud
      ! with a separate treatment of convective cloud
      INTEGER,Parameter:: IP_cloud_part_corr_cnv=8

! CLSCHM3A end
! CLDREG3A defines reference numbers for regions of clouds.in two-stream
! radiation code.

      INTEGER,PARAMETER:: IP_REGION_CLEAR=1 ! clear-sky region
      INTEGER,PARAMETER:: IP_REGION_STRAT=2 ! stratiform cloudy region
      INTEGER,PARAMETER:: IP_REGION_CONV=3  ! convective cloudy region

! CLDREG3A end
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
     &   , N_CLOUD_TOP
!             TOPMOST CLOUDY LAYER
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_CLOUD                                                      &
!             CLOUD SCHEME USED
     &   , I_CLOUD_REPRESENTATION                                       &
!             REPRESENTATION OF CLOUDS USED
     &   , N_CLOUD_TYPE
!             NUMBER OF TYPES OF CLOUD
!
      REAL                                                              &
                !, INTENT(OUT)
     &     FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             FRACTIONS OF EACH TYPE OF CLOUD
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     I_REGION_CLOUD(NPD_CLOUD_TYPE)
!             REGIONS IN WHICH PARTICULAR TYPES OF CLOUD FALL
      REAL                                                              &
                !, INTENT(OUT)
     &     FRAC_REGION(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             FRACTIONS OF TOTAL CLOUD OCCUPIED BY EACH REGION
!
!
!     LOCAL VARIABLES
      INTEGER                                                           &
     &     I                                                            &
!            LOOP VARIABLE
     &   , L                                                            &
!            LOOP VARIABLE
     &   , K
!            LOOP VARIABLE
!
!
!
      IF ( (I_CLOUD == IP_CLOUD_TRIPLE).OR.                             &
     &     (I_CLOUD == IP_CLOUD_PART_CORR_CNV) ) THEN
!
         IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CSIW) THEN
!
            DO K=1, N_CLOUD_TYPE
               IF (K == IP_CLOUD_TYPE_SW) THEN
                  I_REGION_CLOUD(K)=IP_REGION_STRAT
               ELSE IF (K == IP_CLOUD_TYPE_SI) THEN
                  I_REGION_CLOUD(K)=IP_REGION_STRAT
               ELSE IF (K == IP_CLOUD_TYPE_CW) THEN
                  I_REGION_CLOUD(K)=IP_REGION_CONV
               ELSE IF (K == IP_CLOUD_TYPE_CI) THEN
                  I_REGION_CLOUD(K)=IP_REGION_CONV
               ENDIF
            ENDDO
!
            DO I=N_CLOUD_TOP, N_LAYER
               DO L=1, N_PROFILE
                  FRAC_REGION(L, I, IP_REGION_STRAT)                    &
     &               =FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SW)                &
     &               +FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SI)
                  FRAC_REGION(L, I, IP_REGION_CONV)                     &
     &               =FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CW)                &
     &               +FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CI)
               ENDDO
            ENDDO
!
         ELSE IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CONV_STRAT) THEN
!
            DO K=1, N_CLOUD_TYPE
               IF (K == IP_CLOUD_TYPE_STRAT) THEN
                  I_REGION_CLOUD(K)=IP_REGION_STRAT
               ELSE IF (K == IP_CLOUD_TYPE_CONV) THEN
                  I_REGION_CLOUD(K)=IP_REGION_CONV
               ENDIF
            ENDDO
!
            DO I=N_CLOUD_TOP, N_LAYER
               DO L=1, N_PROFILE
                  FRAC_REGION(L, I, IP_REGION_STRAT)                    &
     &               =FRAC_CLOUD(L, I, IP_CLOUD_TYPE_STRAT)
                  FRAC_REGION(L, I, IP_REGION_CONV)                     &
     &               =FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CONV)
               ENDDO
            ENDDO
!
!
         ELSE
            WRITE(IU_ERR, '(/A)')                                       &
     &         '*** ERROR: THIS REPRESENTATION OF CLOUDS IS NOT '       &
     &         //'COMPATIBLE WITH THE TRIPLE OVERLAP.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE AGGREGATE_CLOUD
