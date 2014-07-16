


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to find energy transfer coefficients for triple overlap.
!
! Method:
!       Energy transfer coefficients for upward and downward radiation
!       at the edges of the layers are calculated assuming maximal
!       overlap of regions of the same nature and random overlap of
!       regions of a different nature.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.2             24-05-96                Original Code
!       4.3             20-02-97                Vector searching
!                                               routine WHENFGT
!                                               replaced by IF-tests.
!                                               (J. M. Edwards)
!       5.1             04-04-00                Tolerances replaced
!                                               by F90 intrinsics.
!                                               (J. M. Edwards)
!       5.3             04-10-01                Number of regions
!                                               passed explicitly.
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)

!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE OVERLAP_TRIPLE(N_PROFILE, N_LAYER, N_CLOUD_TOP         &
     &   , W_CLOUD, W_FREE, N_REGION, FRAC_REGION                       &
     &   , CLOUD_OVERLAP                                                &
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
! CLCFPT3A defines pointers in CLOUD_OVERLAP for two-stream radiation
! code.
!
! note that several pointers are identical since only certain
! groups of coefficients are relevant to a particular scheme.

      INTEGER,PARAMETER:: IP_CLOVLP_GFF=1 ! GAMMA-FREE-FREE
      INTEGER,PARAMETER:: IP_CLOVLP_GFC=2 ! GAMMA-FREE-CLOUD
      INTEGER,PARAMETER:: IP_CLOVLP_GCF=3 ! GAMMA-CLOUD-FREE
      INTEGER,PARAMETER:: IP_CLOVLP_GCC=4 ! GAMMA-CLOUD-CLOUD
      INTEGER,PARAMETER:: IP_CLOVLP_BFF=5 ! BETA-FREE-FREE
      INTEGER,PARAMETER:: IP_CLOVLP_BFC=6 ! BETA-FREE-CLOUD
      INTEGER,PARAMETER:: IP_CLOVLP_BCF=7 ! BETA-CLOUD-FREE
      INTEGER,PARAMETER:: IP_CLOVLP_BCC=8 ! BETA-CLOUD-CLOUD
      INTEGER,PARAMETER:: IP_CLOVLP_GFM=5 ! GAMMA_F-
      INTEGER,PARAMETER:: IP_CLOVLP_GFP=6 ! GAMMA_F+
      INTEGER,PARAMETER:: IP_CLOVLP_BFM=7 ! BETA_F-
      INTEGER,PARAMETER:: IP_CLOVLP_BFP=8 ! BETA_F+
      INTEGER,PARAMETER:: IP_CLOVLP_GM=5  ! GAMMA_-
      INTEGER,PARAMETER:: IP_CLOVLP_GP=6  ! GAMMA_+
      INTEGER,PARAMETER:: IP_CLOVLP_BM=7  ! BETA_-
      INTEGER,PARAMETER:: IP_CLOVLP_BP=8  ! BETA_+

      ! pointers for triple overlaps:

      INTEGER,PARAMETER:: IP_CLOVLP_V11=1
      INTEGER,PARAMETER:: IP_CLOVLP_V12=2
      INTEGER,PARAMETER:: IP_CLOVLP_V13=3
      INTEGER,PARAMETER:: IP_CLOVLP_V21=4
      INTEGER,PARAMETER:: IP_CLOVLP_V22=5
      INTEGER,PARAMETER:: IP_CLOVLP_V23=6
      INTEGER,PARAMETER:: IP_CLOVLP_V31=7
      INTEGER,PARAMETER:: IP_CLOVLP_V32=8
      INTEGER,PARAMETER:: IP_CLOVLP_V33=9
      INTEGER,PARAMETER:: IP_CLOVLP_U11=10
      INTEGER,PARAMETER:: IP_CLOVLP_U12=11
      INTEGER,PARAMETER:: IP_CLOVLP_U13=12
      INTEGER,PARAMETER:: IP_CLOVLP_U21=13
      INTEGER,PARAMETER:: IP_CLOVLP_U22=14
      INTEGER,PARAMETER:: IP_CLOVLP_U23=15
      INTEGER,PARAMETER:: IP_CLOVLP_U31=16
      INTEGER,PARAMETER:: IP_CLOVLP_U32=17
      INTEGER,PARAMETER:: IP_CLOVLP_U33=18

! CLCFPT3A end
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
! CLDREG3A defines reference numbers for regions of clouds.in two-stream
! radiation code.

      INTEGER,PARAMETER:: IP_REGION_CLEAR=1 ! clear-sky region
      INTEGER,PARAMETER:: IP_REGION_STRAT=2 ! stratiform cloudy region
      INTEGER,PARAMETER:: IP_REGION_CONV=3  ! convective cloudy region

! CLDREG3A end
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , N_CLOUD_TOP                                                  &
!             TOPMOST CLOUDY LAYER
     &   , N_REGION
!             Number of cloudy regions
      REAL                                                              &
                !, INTENT(IN)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)                              &
!             CLOUD AMOUNTS
     &   , FRAC_REGION(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             FRACTIONS OF TOTAL CLOUD AMOUNT OCCUPIED BY
!             DIFFERENT REGIONS
!
      REAL                                                              &
                !, INTENT(OUT)
     &     W_FREE(NPD_PROFILE, NPD_LAYER)                               &
!             CLOUD-FREE AMOUNTS
     &   , CLOUD_OVERLAP(NPD_PROFILE, 0: NPD_LAYER, NPD_OVERLAP_COEFF)
!             COEFFICIENTS FOR TRANSFER OF ENERGY AT INTERFACE
!
!
!     LOCAL ARGUMENTS.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , J                                                            &
!             LOOP VARIABLE
     &   , L                                                            &
!             LOOP VARIABLE
     &   , K
!             LOOP VARIABLE
!
!
!     FIXED LOCAL VALUES:
      REAL                                                              &
     &     AREA_LOWER(NPD_PROFILE, NPD_REGION)                          &
!             AREAS OF REGIONS IN LOWER LAYER
     &   , AREA_UPPER(NPD_PROFILE, NPD_REGION)                          &
!             AREAS OF REGIONS IN LOWER LAYER
     &   , AREA_OVERLAP(NPD_PROFILE, NPD_REGION, NPD_REGION)            &
!             AREAS OF REGIONS IN LOWER LAYER
     &   , AREA_RANDOM(NPD_PROFILE)
!             AREAS OF RANDOM OVERLAP
!
!
!
!     SET THE FREE FRACTIONS IN EACH LAYER.
      DO I=1, N_CLOUD_TOP-1
         DO L=1, N_PROFILE
            W_FREE(L, I)=1.0E+00
         ENDDO
      ENDDO
      DO I=N_CLOUD_TOP, N_LAYER
         DO L=1, N_PROFILE
            W_FREE(L, I)=1.0E+00-W_CLOUD(L, I)
         ENDDO
      ENDDO
!
!
!
!     WE CONSIDER EACH BOUNDARY IN TURN, COMPARING THE FRACTIONS
!     OF EACH REGION IN THE LAYERS ABOVE AND BELOW THE BOUNDARY.
!
!     INITIALIZE FOR THE LAYER ABOVE THE CLOUDS.
      DO L=1, N_PROFILE
         AREA_UPPER(L, IP_REGION_CLEAR)=1.0E+00
         AREA_UPPER(L, IP_REGION_STRAT)=0.0E+00
         AREA_UPPER(L, IP_REGION_CONV)=0.0E+00
      ENDDO

      DO I=N_CLOUD_TOP-1, N_LAYER
!
!        SET AREAS OF THE REGIONS IN THE LOWER LAYER.
         IF (I <  N_LAYER) THEN
            DO L=1, N_PROFILE
               AREA_LOWER(L, IP_REGION_CLEAR)=W_FREE(L, I+1)
               AREA_LOWER(L, IP_REGION_STRAT)=W_CLOUD(L, I+1)           &
     &            *FRAC_REGION(L, I+1, IP_REGION_STRAT)
               AREA_LOWER(L, IP_REGION_CONV)=W_CLOUD(L, I+1)            &
     &            *FRAC_REGION(L, I+1, IP_REGION_CONV)
            ENDDO
         ELSE
            DO L=1, N_PROFILE
               AREA_LOWER(L, IP_REGION_CLEAR)=1.0E+00
               AREA_LOWER(L, IP_REGION_STRAT)=0.0E+00
               AREA_LOWER(L, IP_REGION_CONV)=0.0E+00
            ENDDO
         ENDIF
!
!        SET THE AREAS OF OVERLAP BETWEEN LIKE REGIONS.
         DO K=1, N_REGION
            DO L=1, N_PROFILE
               AREA_OVERLAP(L, K, K)=MIN(AREA_LOWER(L, K)               &
     &            , AREA_UPPER(L, K))
            ENDDO
         ENDDO
!
!        FIND THE AREAS OF OVERLAP BETWEEN UNLIKE REGIONS. THE OVERLAP
!        BETWEEN UNLIKE REGIONS IS ASSUMED TO BE RANDOM. THESE AREAS
!        ARE SET EQUAL TO 0 FOR THE CASE WHERE THERE IS NO SUCH AREA
!        AND ARE RESET WHEN SUCH AN AREA IS PRESENT.
         DO K=1, N_REGION
            DO J=1, K-1
               DO L=1, N_PROFILE
                  AREA_OVERLAP(L, K, J)=0.0E+00
                  AREA_OVERLAP(L, J, K)=0.0E+00
               ENDDO
            ENDDO
         ENDDO
         DO L=1, N_PROFILE
            AREA_RANDOM(L)=1.0E+00
         ENDDO
         DO K=1, N_REGION
            DO L=1, N_PROFILE
               AREA_RANDOM(L)=AREA_RANDOM(L)-AREA_OVERLAP(L, K, K)
            ENDDO
         ENDDO
         DO K=1, N_REGION
            DO J=1, K-1
               DO L=1, N_PROFILE
                  IF (AREA_RANDOM(L) >  EPSILON(AREA_RANDOM)) THEN
                     AREA_OVERLAP(L, K, J)                              &
     &                  =(AREA_UPPER(L, K)-AREA_OVERLAP(L, K, K))       &
     &                  *(AREA_LOWER(L, J)-AREA_OVERLAP(L, J, J))       &
     &                  /AREA_RANDOM(L)
                     AREA_OVERLAP(L, J, K)                              &
     &                  =(AREA_UPPER(L, J)-AREA_OVERLAP(L, J, J))       &
     &                  *(AREA_LOWER(L, K)-AREA_OVERLAP(L, K, K))       &
     &                  /AREA_RANDOM(L)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
!
!        NOW PROCEED TO FIND THE ENERGY TRANSFER COEFFICIENTS
!        BETWEEN THE VARIOUS REGIONS.
!
!        COEFFICIENTS FOR THE DOWNWARD TRANSFER OF ENERGY:
!
!        TO AVOID DIVISION BY 0 WE INITIALIZE TO DEFAULT VALUES
!        AND RESET.
         DO L=1, N_PROFILE
            CLOUD_OVERLAP(L, I, IP_CLOVLP_V11)=1.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_V21)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_V31)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_V12)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_V22)=1.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_V32)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_V13)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_V23)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_V33)=1.0E+00
         ENDDO
!
!        TRANSFER FROM CLEAR-SKY REGION:
         DO L=1, N_PROFILE
            IF (AREA_UPPER(L, IP_REGION_CLEAR) >                        &
     &          EPSILON(AREA_UPPER)) THEN
               CLOUD_OVERLAP(L, I, IP_CLOVLP_V11)                       &
     &            =AREA_OVERLAP(L, IP_REGION_CLEAR, IP_REGION_CLEAR)    &
     &            /AREA_UPPER(L, IP_REGION_CLEAR)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_V21)                       &
     &            =AREA_OVERLAP(L, IP_REGION_CLEAR, IP_REGION_STRAT)    &
     &            /AREA_UPPER(L, IP_REGION_CLEAR)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_V31)                       &
     &            =AREA_OVERLAP(L, IP_REGION_CLEAR, IP_REGION_CONV)     &
     &            /AREA_UPPER(L, IP_REGION_CLEAR)
            ENDIF
         ENDDO
!
!        TRANSFER FROM STRATIFORM REGION:
         DO L=1, N_PROFILE
            IF (AREA_UPPER(L, IP_REGION_STRAT) >                        &
     &          EPSILON(AREA_UPPER)) THEN
               CLOUD_OVERLAP(L, I, IP_CLOVLP_V12)                       &
     &            =AREA_OVERLAP(L, IP_REGION_STRAT, IP_REGION_CLEAR)    &
     &            /AREA_UPPER(L, IP_REGION_STRAT)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_V22)                       &
     &            =AREA_OVERLAP(L, IP_REGION_STRAT, IP_REGION_STRAT)    &
     &            /AREA_UPPER(L, IP_REGION_STRAT)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_V32)                       &
     &            =AREA_OVERLAP(L, IP_REGION_STRAT, IP_REGION_CONV)     &
     &            /AREA_UPPER(L, IP_REGION_STRAT)
            ENDIF
         ENDDO
!
!        TRANSFER FROM CONVECTIVE REGION:
         DO L=1, N_PROFILE
            IF (AREA_UPPER(L, IP_REGION_CONV) >                         &
     &          EPSILON(AREA_UPPER)) THEN
               CLOUD_OVERLAP(L, I, IP_CLOVLP_V13)                       &
     &            =AREA_OVERLAP(L, IP_REGION_CONV, IP_REGION_CLEAR)     &
     &            /AREA_UPPER(L, IP_REGION_CONV)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_V23)                       &
     &            =AREA_OVERLAP(L, IP_REGION_CONV, IP_REGION_STRAT)     &
     &            /AREA_UPPER(L, IP_REGION_CONV)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_V33)                       &
     &            =AREA_OVERLAP(L, IP_REGION_CONV, IP_REGION_CONV)      &
     &            /AREA_UPPER(L, IP_REGION_CONV)
            ENDIF
         ENDDO
!
!
!        TRANSFER COEFFICIENTS FOR UPWARD FLOW OF ENERGY:
!
!        TO AVOID DIVISION BY 0 WE INITIALIZE TO DEFAULT VALUES
!        AND RESET.
         DO L=1, N_PROFILE
            CLOUD_OVERLAP(L, I, IP_CLOVLP_U11)=1.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_U21)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_U31)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_U12)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_U22)=1.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_U32)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_U13)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_U23)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_U33)=1.0E+00
         ENDDO
!
!        TRANSFER FROM CLEAR-SKY REGION:
         DO L=1, N_PROFILE
            IF (AREA_LOWER(L, IP_REGION_CLEAR) >                        &
     &          EPSILON(AREA_LOWER)) THEN
               CLOUD_OVERLAP(L, I, IP_CLOVLP_U11)                       &
     &            =AREA_OVERLAP(L, IP_REGION_CLEAR, IP_REGION_CLEAR)    &
     &            /AREA_LOWER(L, IP_REGION_CLEAR)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_U21)                       &
     &            =AREA_OVERLAP(L, IP_REGION_STRAT, IP_REGION_CLEAR)    &
     &            /AREA_LOWER(L, IP_REGION_CLEAR)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_U31)                       &
     &            =AREA_OVERLAP(L, IP_REGION_CONV, IP_REGION_CLEAR)     &
     &            /AREA_LOWER(L, IP_REGION_CLEAR)
            ENDIF
         ENDDO
!
!        TRANSFER FROM STRATIFORM REGION:
         DO L=1, N_PROFILE
            IF (AREA_LOWER(L, IP_REGION_STRAT) >                        &
     &          EPSILON(AREA_LOWER)) THEN
               CLOUD_OVERLAP(L, I, IP_CLOVLP_U12)                       &
     &            =AREA_OVERLAP(L, IP_REGION_CLEAR, IP_REGION_STRAT)    &
     &            /AREA_LOWER(L, IP_REGION_STRAT)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_U22)                       &
     &            =AREA_OVERLAP(L, IP_REGION_STRAT, IP_REGION_STRAT)    &
     &            /AREA_LOWER(L, IP_REGION_STRAT)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_U32)                       &
     &            =AREA_OVERLAP(L, IP_REGION_CONV, IP_REGION_STRAT)     &
     &            /AREA_LOWER(L, IP_REGION_STRAT)
            ENDIF
         ENDDO
!
!        TRANSFER FROM CONVECTIVE REGION:
         DO L=1, N_PROFILE
            IF (AREA_LOWER(L, IP_REGION_CONV) >                         &
     &          EPSILON(AREA_LOWER)) THEN
               CLOUD_OVERLAP(L, I, IP_CLOVLP_U13)                       &
     &            =AREA_OVERLAP(L, IP_REGION_CLEAR, IP_REGION_CONV)     &
     &            /AREA_LOWER(L, IP_REGION_CONV)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_U23)                       &
     &            =AREA_OVERLAP(L, IP_REGION_STRAT, IP_REGION_CONV)     &
     &            /AREA_LOWER(L, IP_REGION_CONV)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_U33)                       &
     &            =AREA_OVERLAP(L, IP_REGION_CONV, IP_REGION_CONV)      &
     &            /AREA_LOWER(L, IP_REGION_CONV)
            ENDIF
         ENDDO
!
!        REASSIGN THE FRACTIONS IN THE UPPER LAYER TO STEP DOWN
!        THROUGH THE ATMOSPHERE.
         IF (I <  N_LAYER) THEN
            DO K=1, N_REGION
               DO L=1, N_PROFILE
                  AREA_UPPER(L, K)=AREA_LOWER(L, K)
               ENDDO
            ENDDO
         ENDIF
!
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE OVERLAP_TRIPLE
