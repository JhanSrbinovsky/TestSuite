
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to find energy transfer coefficients for coupled overlap.
!
! Method:
!       Energy transfer coefficients for upward and downward radiation
!       at the edges of the layers are calculated assuming maximal
!       overlap of regions of the same nature and random overlap of
!       regions of a different nature.
!
!       Storage and Indexing: Now that solvers for the net flux are no
!       longer supported, the overlap coefficients can be stored more
!       easily. The coefficient referring to downward transfer from the
!       kth to the jth region is stored with a third index of
!       K+N_REGION*(J-1): the coeffieint for upward transfer is stored
!       with an index of N_REGION*(N_REGION+J-1)+K, so that in both
!       cases the originating region changes most frequently with the
!       index.
!
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE OVERLAP_COUPLED(N_PROFILE, N_LAYER, N_CLOUD_TOP        &
     &  , W_CLOUD, W_FREE, N_REGION, TYPE_REGION, FRAC_REGION, P        &
     &  , I_CLOUD                                                       &
     &  , CLOUD_OVERLAP                                                 &
     &  , ND_PROFILE, ND_LAYER, ND_OVERLAP_COEFF, ND_REGION, ND_FIELD   &
     &  , ID_CT, DP_CORR_STRAT, DP_CORR_CONV, TOT_CLOUD_COVER           &
     &  , LATITUDE)
!
!
      Use rad_switches_mod, ONLY:                                       &
          lrad_exprand
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_FIELD                                                      &
!            FIELD SIZE IN CALLING PROGRAM
     &  , ND_PROFILE                                                    &
!           Maximum number of profiles
     &  , ND_LAYER                                                      &
!           Maximum number of layers
     &  , ND_OVERLAP_COEFF                                              &
!           Maximum number of overlap coefficients
     &  , ND_REGION                                                     &
!           Maximum number of regions
     &  , ID_CT
!           Topmost declared cloudy layer
!
!     Include header files.
! Start C_KINDS

! Description:
!   Contains parameters defining kinds for 32 and 64 integers
!   and reals.
!
! Current Code Owner: Paul Selwood
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   6.1  10/03/04   Original code. Paul Selwood

! Parameters for 32 and 64 bit kinds

! Precision and range for 64 bit real
      Integer, Parameter :: prec64  = 15
      Integer, Parameter :: range64 = 307

! Precision and range for 32 bit real
      Integer, Parameter :: prec32  = 6
      Integer, Parameter :: range32 = 37

! Range for integers
      Integer, Parameter :: irange64=15
      Integer, Parameter :: irange32=9

! Kind for 64 bit real
      Integer, Parameter :: real64  = selected_real_kind(prec64,range64)
! Kind for 32 bit real
      Integer, Parameter :: real32  = selected_real_kind(prec32,range32)
! Kind for 64 bit integer
      Integer, Parameter :: integer64 = selected_int_kind(irange64)
! Kind for 32 bit integer
      Integer, Parameter :: integer32 = selected_int_kind(irange32)

! Kinds for 64 and 32 bit logicals. Note that there is no
! "selected_logical_kind", but using the equivalent integer kind is a
! workaround that works on every platform we've tested.
      Integer, Parameter :: logical64 = integer64
      Integer, Parameter :: logical32 = integer32

! End C_KINDS
!     ------------------------------------------------------------------
!     Module to define reference numbers for regions of clouds.
!
      INTEGER                                                           &
     &    IP_REGION_CLEAR                                               &
!           Reference number for clear-sky region
     &  , IP_REGION_STRAT                                               &
!           Reference number for stratiform cloudy region
     &  , IP_REGION_CONV
!           Reference number for convective cloudy region
!
      PARAMETER(                                                        &
     &    IP_REGION_CLEAR=1                                             &
     &  , IP_REGION_STRAT=2                                             &
     &  , IP_REGION_CONV=3                                              &
     &  )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     Module to define reference numbers for cloud schemes.
!
      INTEGER, Parameter :: IP_cloud_mix_max         = 2
!           Maximum/random overlap in a mixed column
      INTEGER, Parameter :: IP_cloud_mix_random      = 4
!           Random overlap in a mixed column
      INTEGER, Parameter :: IP_cloud_column_max      = 3
!           Maximum overlap in a column model
      INTEGER, Parameter :: IP_cloud_clear           = 5
!           Clear column
      INTEGER, Parameter :: IP_cloud_triple          = 6
!           Mixed column with split between convective and layer cloud
      INTEGER, Parameter :: IP_cloud_part_corr       = 7
!           Coupled overlap with partial correlation of cloud
      INTEGER, Parameter :: IP_cloud_part_corr_cnv   = 8
!           Coupled overlap with partial correlation of cloud
!           with a separate treatment of convective cloud
!
!     ------------------------------------------------------------------
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , N_CLOUD_TOP                                                   &
!           Topmost cloudy layer
     &  , N_REGION                                                      &
!           Number of cloudy regions
     &  , TYPE_REGION(ND_REGION)                                        &
!           Array holding the type of each region
     &  , I_CLOUD
!           Cloud scheme selected
      REAL  (Real64), INTENT(IN) ::                                     &
     &    W_CLOUD(ND_PROFILE, ID_CT: ND_LAYER)                          &
!           Cloud amounts
     &  , FRAC_REGION(ND_PROFILE, ID_CT: ND_LAYER, ND_REGION)           &
!           Fractions of total cloud amount occupied by
!           different regions
     &  , P(ND_PROFILE, ND_LAYER)                                       &
!           Pressures at the middles of layers
     &  , DP_CORR_STRAT                                                 &
!           Decorrelation pressure scale for large scale cloud
     &  , DP_CORR_CONV                                                  &
!           Decorrelation pressure scale for convective cloud
     &  , LATITUDE(ND_FIELD)
!           Latitude distribution
!
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    W_FREE(ND_PROFILE, ID_CT: ND_LAYER)                           &
!           Cloud-free amounts
     &  , CLOUD_OVERLAP(ND_PROFILE, ID_CT-1: ND_LAYER                   &
     &      , ND_OVERLAP_COEFF)
!           Coefficients for transfer of energy at interface
!
      REAL, INTENT(OUT) ::                                              &
     &     TOT_CLOUD_COVER(ND_PROFILE)
!             Total cloud cover
!
!     Local arguments.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , J                                                             &
!           Loop variable
     &  , L                                                             &
!           Loop variable
     &  , K
!           Loop variable
!
!
!     Fixed local values:
      REAL  (Real64) ::                                                 &
     &    DP_CORR                                                       &
!           Pressure scale over which correlation between cloudy
!           layers is lost
     &  , DZ_CORR(ND_PROFILE)                                           &
!           Height scale over which correlation between cloudy
!           layers is lost
     &  , DP_CORR_FULL(ND_PROFILE, ND_REGION)                           &
!           Field of pressure scales calculated from the decorrelation
!           height 
     &  , CORR_FACTOR(ND_PROFILE, ND_REGION)
!           Correlation factors for each region across the boundary
!           between layers: this represents the fraction of the
!           potentially maximally overlapped region that is actually
!           maximally overlapped.
      REAL  (Real64) ::                                                 &
     &    AREA_LOWER(ND_PROFILE, ND_REGION)                             &
!           Areas of regions in lower layer
     &  , AREA_UPPER(ND_PROFILE, ND_REGION)                             &
!           Areas of regions in lower layer
     &  , AREA_OVERLAP(ND_PROFILE, ND_REGION, ND_REGION)                &
!           Areas of overlap between the different regions:
!           the first index refers to the upper layer
     &  , AREA_RANDOM_UPPER(ND_PROFILE, ND_REGION)                      &
!           Areas of each region in the upper layer
!           to be overlapped randomly
     &  , AREA_RANDOM_LOWER(ND_PROFILE, ND_REGION)                      &
!           Areas of each region in the lower layer
!           to be overlapped randomly
     &  , AREA_RANDOM_TOT(ND_PROFILE)                                   &
!           Total randomly overlapped area
     &  , TOL_CLOUD
!           Tolerance used to detect cloud amounts of 0
!
!
      TOL_CLOUD=1.0E+02_Real64*EPSILON(TOL_CLOUD)
!
!     Set the free fractions in each layer.
      DO I=N_CLOUD_TOP, N_LAYER
        DO L=1, N_PROFILE
          W_FREE(L, I)=1.0E+00_Real64-W_CLOUD(L, I)
        ENDDO
      ENDDO
!
!     Use the total cloud cover temporarily to hold the clear-sky
!     fraction and convert back to cloud cover later.
      DO L=1, N_PROFILE
         TOT_CLOUD_COVER(L)=1.0E+00_Real64 - W_CLOUD(L, N_CLOUD_TOP)
      ENDDO
!
!
!     We consider each boundary in turn, comparing the fractions
!     of each region in the layers above and below the boundary.
!
!     Initialize for the layer above the clouds: here the clear
!     region will cover the grid-box.
      DO K=1, N_REGION
        IF (TYPE_REGION(K) == IP_REGION_CLEAR) THEN
          DO L=1, N_PROFILE
            AREA_UPPER(L, K)=1.0E+00_Real64
          ENDDO
        ELSE
          DO L=1, N_PROFILE
            AREA_UPPER(L, K)=0.0E+00_Real64
          ENDDO
        ENDIF
      ENDDO
!
      DO I=N_CLOUD_TOP-1, N_LAYER
!
!       Set the correlations between like regions at each interface.
!
        IF ( (I_CLOUD == IP_CLOUD_TRIPLE).OR.                           &
     &       (I_CLOUD == IP_CLOUD_MIX_MAX) ) THEN
!
          DO K=1, N_REGION
            DO L=1, N_PROFILE
              CORR_FACTOR(L, K)=1.0E+00_Real64
            ENDDO
          ENDDO
!
        ELSE IF (I_CLOUD == IP_CLOUD_MIX_RANDOM) THEN
!
          DO K=1, N_REGION
            DO L=1, N_PROFILE
              CORR_FACTOR(L, K)=0.0E+00_Real64
            ENDDO
          ENDDO
!
        ELSE IF ( (I_CLOUD == IP_CLOUD_PART_CORR).OR.                   &
     &            (I_CLOUD == IP_CLOUD_PART_CORR_CNV) ) THEN
!
!
!        If exponential-random overlap is selected, we 
!        calculate decorrelation according to a constant height scale that
!        varies with latitude; if not, we perform decorrelation overlap as in
!        the original code.
!
          IF (LRAD_EXPRAND) THEN
!
            DO K=1, N_REGION
!
!             Calculate the constant decorrelation height scales from the 
!             local latitude.
              IF (TYPE_REGION(K) == IP_REGION_CLEAR) THEN
                DO L=1, N_PROFILE
                  DZ_CORR(L)=(2.174-0.02069*ABS(LATITUDE(L)))*1000.0
!                  DZ_CORR(L)=(2.174-0.02069*45.0)*1000.0
                ENDDO
              ELSE IF (TYPE_REGION(K) == IP_REGION_STRAT) THEN
                DO L=1, N_PROFILE
                  IF (I_CLOUD == IP_CLOUD_PART_CORR) THEN
                    DZ_CORR(L)=(2.174-0.02069*ABS(LATITUDE(L)))*1000.0
!                    DZ_CORR(L)=(2.174-0.02069*45.0)*1000.0
                  ELSE IF (I_CLOUD == IP_CLOUD_PART_CORR_CNV) THEN
                    DZ_CORR(L)=(2.174-0.02069*ABS(LATITUDE(L)))*500.0
!                    DZ_CORR(L)=(2.174-0.02069*45.0)*500.0
                  ENDIF
                ENDDO
              ELSE IF (TYPE_REGION(K) == IP_REGION_CONV) THEN
                DO L=1, N_PROFILE
                  DZ_CORR(L)=(2.174-0.02069*ABS(LATITUDE(L)))*500.0
!                  DZ_CORR(L)=(2.174-0.02069*45.0)*500.0
                ENDDO
              ENDIF
!
!             Convert the decorrelation heights to decorrelation pressures
!             assuming hydrostatic balance.
              DO L=1, N_PROFILE
                DP_CORR_FULL(L, K) = DZ_CORR(L)*P(L, I)*9.81/(287.0*250.0)
              ENDDO
!
!             Calculate the overlap parameter from the decorrelation
!             pressures.
              IF ( (I <  N_LAYER).AND.(I >  1) ) THEN
                DO L=1, N_PROFILE
                  CORR_FACTOR(L, K)=EXP((P(L, I)-P(L, I+1))              &
                                       /DP_CORR_FULL(L, K))
                  IF (CORR_FACTOR(L, K) .GT. 1) THEN
                     CORR_FACTOR(L, K) = 1
                  ELSEIF (CORR_FACTOR(L, K) .LT. 0) THEN
                     CORR_FACTOR(L, K) = 0
                  ENDIF
                ENDDO
              ELSE
!               At the surface and the top of the atmosphere
!               the correlation factor is irrelevant.
                DO L=1, N_PROFILE
                  CORR_FACTOR(L, K)=1.0E+00_Real64

                  IF (CORR_FACTOR(L, K) .GT. 1) THEN
                     CORR_FACTOR(L, K) = 1
                  ELSEIF (CORR_FACTOR(L, K) .LT. 0) THEN
                     CORR_FACTOR(L, K) = 0
                  ENDIF
                ENDDO
              ENDIF
!
            ENDDO
!
          ELSE
!
            DO K=1, N_REGION
!
!             Experimental version: set the pressure scales over
!             which decorrelation occurs.
              IF (TYPE_REGION(K) == IP_REGION_CLEAR) THEN
                DP_CORR=1.0E+00_Real64 
              ELSE IF (TYPE_REGION(K) == IP_REGION_STRAT) THEN
                DP_CORR=DP_CORR_STRAT
              ELSE IF (TYPE_REGION(K) == IP_REGION_CONV) THEN
                DP_CORR=DP_CORR_CONV
              ENDIF
!
              IF ( (I <  N_LAYER).AND.(I >  1) ) THEN
                DO L=1, N_PROFILE
                  CORR_FACTOR(L, K)=EXP((P(L, I)-P(L, I+1))/DP_CORR)
                ENDDO
              ELSE
!               At the surface and the top of the atmosphere
!               the correlation factor is irrelevant.
                DO L=1, N_PROFILE
                  CORR_FACTOR(L, K)=1.0E+00_Real64
                ENDDO
              ENDIF
!
            ENDDO
!
          ENDIF  ! LRAD_EXPRAND
!
        ENDIF
!
!       Set areas of the regions in the lower layer.
        DO K=1, N_REGION
          IF (I <  N_LAYER) THEN
            IF (TYPE_REGION(K) == IP_REGION_CLEAR) THEN
              DO L=1, N_PROFILE
                AREA_LOWER(L, K)=W_FREE(L, I+1)
              ENDDO
            ELSE
              DO L=1, N_PROFILE
                AREA_LOWER(L, K)=W_CLOUD(L, I+1)                        &
     &            *FRAC_REGION(L, I+1, K)
              ENDDO
            ENDIF
          ELSE
!           At the very bottom of the column we imagine a notional
!           clear layer below the ground surface.
            IF (TYPE_REGION(K) == IP_REGION_CLEAR) THEN
              DO L=1, N_PROFILE
                AREA_LOWER(L, K)=1.0E+00_Real64
              ENDDO
            ELSE
              DO L=1, N_PROFILE
                AREA_LOWER(L, K)=0.0E+00_Real64
              ENDDO
            ENDIF
          ENDIF
!
!         Begin by setting the maximally overlapped parts of the
!         atmospheric column. The area of common overlap betwen
!         like regions may be incremented by randomly overlapped
!         fractions later.
!
          DO L=1, N_PROFILE
            AREA_OVERLAP(L, K, K)=CORR_FACTOR(L, K)                     &
     &        *MIN(AREA_LOWER(L, K), AREA_UPPER(L, K))
          ENDDO
!
        ENDDO
!
!       Find the remaining areas of overlap on the assumption that
!       the overlap is random. We initialize the areas of overlap to
!       0 and reset later when such an area is present.
        DO K=1, N_REGION
          DO J=1, K-1
            DO L=1, N_PROFILE
              AREA_OVERLAP(L, K, J)=0.0E+00_Real64
              AREA_OVERLAP(L, J, K)=0.0E+00_Real64
            ENDDO
          ENDDO
        ENDDO
!
        DO L=1, N_PROFILE
          AREA_RANDOM_TOT(L)=1.0E+00_Real64-AREA_OVERLAP(L, 1, 1)
        ENDDO
        DO K=2, N_REGION
          DO L=1, N_PROFILE
            AREA_RANDOM_TOT(L)=AREA_RANDOM_TOT(L)-AREA_OVERLAP(L, K, K)
          ENDDO
        ENDDO
        DO K=1, N_REGION
          DO L=1, N_PROFILE
            AREA_RANDOM_UPPER(L, K)                                     &
     &        =AREA_UPPER(L, K)-AREA_OVERLAP(L, K, K)
            AREA_RANDOM_LOWER(L, K)                                     &
     &        =AREA_LOWER(L, K)-AREA_OVERLAP(L, K, K)
          ENDDO
        ENDDO
!       To calculate the contributions of random overlap to the
!       areas of overlap we take the randomly overlapped portion
!       of the kth region in the upper layer. The probability that
!       this is overalpped with the randomly overlapped portion of
!       the jth region in the lower layer will be equal to
!       the randomly overlapped area of the lower jth region divided
!       by the total randomly overalpped area. The ratio might become
!       ill-conditioned for small amounts of cloud, the but this
!       should not be an issue as the randomly overalpped area would
!       then be small.
        DO K=1, N_REGION
          DO J=1, N_REGION
            DO L=1, N_PROFILE
              IF (AREA_RANDOM_TOT(L) >  TOL_CLOUD) THEN
                AREA_OVERLAP(L, K, J)=AREA_OVERLAP(L, K, J)             &
     &            +AREA_RANDOM_UPPER(L, K)                              &
     &            *AREA_RANDOM_LOWER(L, J)/AREA_RANDOM_TOT(L)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
!
!       Now proceed to find the energy transfer coefficients
!       between the various regions.
!
!       Coefficients for the downward transfer of energy:
!
!       To avoid division by 0 we initialize to default values
!       and reset.
        DO K=1, N_REGION
          DO L=1, N_PROFILE
            CLOUD_OVERLAP(L, I, N_REGION*(K-1)+K)=1.0E+00_Real64
          ENDDO
          DO J=1, K-1
            DO L=1, N_PROFILE
              CLOUD_OVERLAP(L, I, N_REGION*(J-1)+K)=0.0E+00_Real64
              CLOUD_OVERLAP(L, I, N_REGION*(K-1)+J)=0.0E+00_Real64
            ENDDO
          ENDDO
        ENDDO
!
        DO K=1, N_REGION
          DO L=1, N_PROFILE
            IF (AREA_UPPER(L, K) >  TOL_CLOUD) THEN
              DO J=1, N_REGION
                CLOUD_OVERLAP(L, I, N_REGION*(J-1)+K)                   &
     &            =AREA_OVERLAP(L, K, J)/AREA_UPPER(L, K)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
!
!
!       Coefficients for upward flow of energy:
!
!       To avoid division by 0 we initialize to default values
!       and reset.
        DO K=1, N_REGION
          DO L=1, N_PROFILE
            CLOUD_OVERLAP(L, I, N_REGION*(N_REGION+K-1)+K)              &
     &        =1.0E+00_Real64
          ENDDO
          DO J=1, K-1
            DO L=1, N_PROFILE
              CLOUD_OVERLAP(L, I, N_REGION*(N_REGION+J-1)+K)            &
     &          =0.0E+00_Real64
              CLOUD_OVERLAP(L, I, N_REGION*(N_REGION+K-1)+J)            &
     &          =0.0E+00_Real64
            ENDDO
          ENDDO
        ENDDO
!
        DO K=1, N_REGION
          DO L=1, N_PROFILE
            IF (AREA_LOWER(L, K) >  TOL_CLOUD) THEN
              DO J=1, N_REGION
                CLOUD_OVERLAP(L, I, N_REGION*(N_REGION+J-1)+K)          &
     &            =AREA_OVERLAP(L, J, K)/AREA_LOWER(L, K)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
!
!
!       Reassign the fractions in the upper layer to step down
!       through the atmosphere.
        IF (I <  N_LAYER) THEN
          DO K=1, N_REGION
            DO L=1, N_PROFILE
              AREA_UPPER(L, K)=AREA_LOWER(L, K)
            ENDDO
          ENDDO
        ENDIF


!       Now calculate total column cloud cover for use as a diagnostic.
!       We calculate this quantity by imagining a totally transparent
!       atmosphere containing totally opaque clouds and finding the
!       transmission.

        IF ( (I_CLOUD == IP_CLOUD_MIX_MAX)    .OR.                      &
     &       (I_CLOUD == IP_CLOUD_MIX_RANDOM) .OR.                      &
     &       (I_CLOUD == IP_CLOUD_PART_CORR) ) THEN
!
! 
           IF ( (I >= N_CLOUD_TOP).AND.(I < N_LAYER) ) THEN
              DO L=1, N_PROFILE
                 TOT_CLOUD_COVER(L)=TOT_CLOUD_COVER(L)                  &
     &              *( 1.0E+00_Real64-W_CLOUD(L, I+1) ) /               &
     &               MAX( 1.0E+00_Real64-CORR_FACTOR(L, 2)*             &
     &                      MIN(W_CLOUD(L, I),W_CLOUD(L, I+1)),         &
     &                          EPSILON(1.0_Real64) )
              ENDDO
           ENDIF
       
        ELSE IF ( (I_CLOUD == IP_CLOUD_TRIPLE)         .OR.             &
     &            (I_CLOUD == IP_CLOUD_PART_CORR_CNV) ) THEN

           IF ( (I >= N_CLOUD_TOP).AND.(I < N_LAYER) ) THEN

              DO L=1, N_PROFILE
                 TOT_CLOUD_COVER(L)=TOT_CLOUD_COVER(L)                  &
     &              *( 1.0E+00_Real64-W_CLOUD(L, I+1) ) /               &
     &               MAX( 1.0E+00_Real64-CORR_FACTOR(L, 2)*             &
     &                      MIN( W_CLOUD(L,I  )*FRAC_REGION(L,I,  2),   &
     &                           W_CLOUD(L,I+1)*FRAC_REGION(L,I+1,2) )  &
     &                         -CORR_FACTOR(L, 3)*                      &
     &                      MIN( W_CLOUD(L,I  )*FRAC_REGION(L,I,  3),   &
     &                           W_CLOUD(L,I+1)*FRAC_REGION(L,I+1,3) ), &
     &                    EPSILON(1.0_Real64) )
              ENDDO
           ENDIF        
        
        ENDIF

      ENDDO

      DO L=1, N_PROFILE
         TOT_CLOUD_COVER(L)=1.0E+00_Real64-TOT_CLOUD_COVER(L)
      ENDDO


      RETURN
      END SUBROUTINE OVERLAP_COUPLED
