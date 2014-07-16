

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the parametrization schemes for clouds.
!
! Purpose:
!   The parametrization schemes for each component within a cloud
!   are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: James Manners
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_CLOUD_PARAMETRIZATION(IERR, N_BAND              &
     &   , I_ST_WATER, I_CNV_WATER, I_ST_ICE, I_CNV_ICE                 &
     &   , L_DROP_TYPE, I_DROP_PARAMETRIZATION, DROP_PARAMETER_LIST     &
     &   , DROP_PARM_MIN_DIM, DROP_PARM_MAX_DIM                         &
     &   , L_ICE_TYPE, I_ICE_PARAMETRIZATION, ICE_PARAMETER_LIST        &
     &   , ICE_PARM_MIN_DIM, ICE_PARM_MAX_DIM                           &
     &   , I_CONDENSED_PARAM, CONDENSED_PARAM_LIST                      &
     &   , CONDENSED_MIN_DIM, CONDENSED_MAX_DIM                         &
     &   , NPD_BAND, NPD_DROP_TYPE, NPD_ICE_TYPE, NPD_CLOUD_PARAMETER   &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED.
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
!
!
!     DUMMY ARGUMENTS:
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
!
!     SIZES OF ARRAYS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_BAND                                                     &
!             MAXIMUM NUMBER OF SPECTRAL BANDS
     &   , NPD_DROP_TYPE                                                &
!             MAXIMUM NUMBER OF TYPES OF DROPLETS
     &   , NPD_ICE_TYPE                                                 &
!             MAXIMUM NUMBER OF TYPES OF ICE CRYSTALS
     &   , NPD_CLOUD_PARAMETER
!             MAXIMUM NUMBER OF PARAMETERS FOR CLOUDS
!
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_BAND
!             NUMBER OF SPECTRAL BANDS
!
!     TYPES OF DROPLETS AND CRYSTALS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_ST_WATER                                                   &
!             TYPE OF WATER DROPLETS IN STRATIFORM CLOUDS
     &   , I_CNV_WATER                                                  &
!             TYPE OF WATER DROPLETS IN CONVECTIVE CLOUDS
     &   , I_ST_ICE                                                     &
!             TYPE OF ICE CRYSTALS IN STRATIFORM CLOUDS
     &   , I_CNV_ICE
!             TYPE OF ICE CRYSTALS IN CONVECTIVE CLOUDS
!
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_DROP_TYPE(NPD_DROP_TYPE)                                   &
!             FLAGS FOR TYPES OF DROPLET PRESENT
     &   , L_ICE_TYPE(NPD_ICE_TYPE)
!             FLAGS FOR TYPES OF ICE CRYSTAL PRESENT
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_DROP_PARAMETRIZATION(NPD_DROP_TYPE)                        &
!             PARAMETRIZATIONS OF TYPES OF DROPLETS
     &   , I_ICE_PARAMETRIZATION(NPD_ICE_TYPE)
!             PARAMETRIZATIONS OF TYPES OF ICE CRYSTALS
      REAL                                                              &
                !, INTENT(IN)
     &     DROP_PARAMETER_LIST(NPD_CLOUD_PARAMETER                      &
     &        , NPD_BAND, NPD_DROP_TYPE)                                &
!             PARAMETERS FOR OPTICAL PARAMETRIZATIONS OF DROPLETS
     &   , DROP_PARM_MIN_DIM(NPD_DROP_TYPE)                             &
!             MINIMUM SIZE OF DROPLETS PERMITTED IN PARAMETRIZATIONS
     &   , DROP_PARM_MAX_DIM(NPD_DROP_TYPE)                             &
!             MAXIMUM SIZE OF DROPLETS PERMITTED IN PARAMETRIZATIONS
     &   , ICE_PARAMETER_LIST(NPD_CLOUD_PARAMETER                       &
     &        , NPD_BAND, NPD_ICE_TYPE)                                 &
!             PARAMETERS FOR OPTICAL PARAMETRIZATIONS OF ICE CRYSTALS
     &   , ICE_PARM_MIN_DIM(NPD_ICE_TYPE)                               &
!             MINIMUM SIZE OF ICE CRYSTALS PERMITTED IN PARAMETRIZATIONS
     &   , ICE_PARM_MAX_DIM(NPD_ICE_TYPE)
!             MAXIMUM SIZE OF ICE CRYSTALS PERMITTED IN PARAMETRIZATIONS
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     I_CONDENSED_PARAM(NPD_CLOUD_COMPONENT)
!             TYPES OF PARAMETRIZATION USED FOR CONDENSED
!             COMPONENTS IN CLOUDS
      REAL                                                              &
                !, INTENT(OUT)
     &     CONDENSED_PARAM_LIST(NPD_CLOUD_PARAMETER                     &
     &        , NPD_CLOUD_COMPONENT, NPD_BAND)                          &
!             COEFFICIENTS FOR PARAMETRIZATION OF CONDENSED PHASES
     &   , CONDENSED_MIN_DIM(NPD_CLOUD_COMPONENT)                       &
!             MINIMUM DIMENSION OF EACH CONDENSED COMPONENT
     &   , CONDENSED_MAX_DIM(NPD_CLOUD_COMPONENT)
!             MAXIMUM DIMENSION OF EACH CONDENSED COMPONENT
!
!
!     LOCAL VARIABLES:
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , J                                                            &
!             LOOP VARIABLE
     &   , I_SCHEME
!             PARAMETRIZATION SCHEME
!
!     FUNCTIONS CALLED:
      INTEGER                                                           &
     &     SET_N_CLOUD_PARAMETER
!             FUNCTION TO FIND NUMBER OF PARAMETERS FOR CLOUDS
      EXTERNAL                                                          &
     &     SET_N_CLOUD_PARAMETER
!
!
!
!     SELECT PARAMETRIZATION FOR WATER IN STRATIFORM CLOUDS:
!
      IF ( (I_ST_WATER <= NPD_DROP_TYPE).AND.                           &
     &     (L_DROP_TYPE(I_ST_WATER)) ) THEN
         I_SCHEME=I_DROP_PARAMETRIZATION(I_ST_WATER)
         I_CONDENSED_PARAM(IP_CLCMP_ST_WATER)=I_SCHEME
         CONDENSED_MIN_DIM(IP_CLCMP_ST_WATER)                           &
     &      =DROP_PARM_MIN_DIM(I_ST_WATER)
         CONDENSED_MAX_DIM(IP_CLCMP_ST_WATER)                           &
     &      =DROP_PARM_MAX_DIM(I_ST_WATER)
      ELSE
         WRITE(IU_ERR, '(/A, /A)') '*** ERROR: NO DATA EXIST FOR TYPE ' &
     &     , 'OF DROPLET SELECTED IN STRATIFORM WATER CLOUDS.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      DO I=1, N_BAND
! DEPENDS ON: set_n_cloud_parameter
         DO J=1, SET_N_CLOUD_PARAMETER(I_SCHEME, IP_CLCMP_ST_WATER)
            CONDENSED_PARAM_LIST(J, IP_CLCMP_ST_WATER, I)               &
     &         =DROP_PARAMETER_LIST(J, I, I_ST_WATER)
         ENDDO
      ENDDO
!
!
!     SELECT PARAMETRIZATION FOR WATER IN CONVECTIVE CLOUDS:
!
      IF ( (I_CNV_WATER <= NPD_DROP_TYPE).AND.                          &
     &     (L_DROP_TYPE(I_CNV_WATER)) ) THEN
         I_SCHEME=I_DROP_PARAMETRIZATION(I_CNV_WATER)
         I_CONDENSED_PARAM(IP_CLCMP_CNV_WATER)=I_SCHEME
         CONDENSED_MIN_DIM(IP_CLCMP_CNV_WATER)                          &
     &      =DROP_PARM_MIN_DIM(I_CNV_WATER)
         CONDENSED_MAX_DIM(IP_CLCMP_CNV_WATER)                          &
     &      =DROP_PARM_MAX_DIM(I_CNV_WATER)
      ELSE
         WRITE(IU_ERR, '(/A, /A)') '*** ERROR: NO DATA EXIST FOR TYPE ' &
     &     , 'OF CRYSTAL SELECTED IN CONVECTIVE WATER CLOUDS.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      DO I=1, N_BAND
! DEPENDS ON: set_n_cloud_parameter
         DO J=1, SET_N_CLOUD_PARAMETER(I_SCHEME, IP_CLCMP_CNV_WATER)
            CONDENSED_PARAM_LIST(J, IP_CLCMP_CNV_WATER, I)              &
     &         =DROP_PARAMETER_LIST(J, I, I_CNV_WATER)
         ENDDO
      ENDDO
!
!
!     SELECT PARAMETRIZATION FOR ICE IN STRATIFORM CLOUDS:
!
      IF ( (I_ST_ICE <= NPD_ICE_TYPE).AND.                              &
     &     (L_ICE_TYPE(I_ST_ICE)) ) THEN
         I_SCHEME=I_ICE_PARAMETRIZATION(I_ST_ICE)
         I_CONDENSED_PARAM(IP_CLCMP_ST_ICE)=I_SCHEME
         CONDENSED_MIN_DIM(IP_CLCMP_ST_ICE)                             &
     &      =ICE_PARM_MIN_DIM(I_ST_ICE)
         CONDENSED_MAX_DIM(IP_CLCMP_ST_ICE)                             &
     &      =ICE_PARM_MAX_DIM(I_ST_ICE)
      ELSE
         WRITE(IU_ERR, '(/A, /A)') '*** ERROR: NO DATA EXIST FOR TYPE ' &
     &      , 'OF CRYSTAL SELECTED IN STRATIFORM ICE CLOUDS.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      DO I=1, N_BAND
! DEPENDS ON: set_n_cloud_parameter
         DO J=1, SET_N_CLOUD_PARAMETER(I_SCHEME, IP_CLCMP_ST_ICE)
            CONDENSED_PARAM_LIST(J, IP_CLCMP_ST_ICE, I)                 &
     &         =ICE_PARAMETER_LIST(J, I, I_ST_ICE)
         ENDDO
      ENDDO
!
!
!     SELECT PARAMETRIZATION FOR ICE IN CONVECTIVE CLOUDS:
!
      IF ( (I_CNV_ICE <= NPD_ICE_TYPE).AND.                             &
     &     (L_ICE_TYPE(I_CNV_ICE)) ) THEN
         I_SCHEME=I_ICE_PARAMETRIZATION(I_CNV_ICE)
         I_CONDENSED_PARAM(IP_CLCMP_CNV_ICE)=I_SCHEME
         CONDENSED_MIN_DIM(IP_CLCMP_CNV_ICE)                            &
     &      =ICE_PARM_MIN_DIM(I_CNV_ICE)
         CONDENSED_MAX_DIM(IP_CLCMP_CNV_ICE)                            &
     &      =ICE_PARM_MAX_DIM(I_CNV_ICE)
      ELSE
         WRITE(IU_ERR, '(/A, /A)') '*** ERROR: NO DATA EXIST FOR TYPE ' &
     &      , 'OF CRYSTAL SELECTED IN CONVECTIVE ICE CLOUDS.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      DO I=1, N_BAND
! DEPENDS ON: set_n_cloud_parameter
         DO J=1, SET_N_CLOUD_PARAMETER(I_SCHEME, IP_CLCMP_CNV_ICE)
            CONDENSED_PARAM_LIST(J, IP_CLCMP_CNV_ICE, I)                &
     &         =ICE_PARAMETER_LIST(J, I, I_CNV_ICE)
         ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE R2_SET_CLOUD_PARAMETRIZATION
