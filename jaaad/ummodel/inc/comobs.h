!-----------------------------------------------------------------------
!LCOMDECK COMOBS
! Description
!   This comdeck provides parameters relating to observations in
!   atmospheric AC assimilation (see DOCOBS for detail)
!
!   History:
!   Model    Date     Modification history
!  version
!   4.4      3/7/97   Add MAX_NUM_ACOB_FILES,PER_FILE_TNDVMAX S. Bell
!   4.5      5/8/98   Increase USED_FILES size S. Bell
!   5.3    07/06/01   Move data statements to blkdata.  A van der Wal
!L-------------------------------------------------------------------
      INTEGER NDATAVMX, NOBLEVMX
      INTEGER MAX_NUM_ACOB_FILES
      INTEGER NUM_OB_FILE_TYPES
      PARAMETER (NDATAVMX = 6+3*P_LEVELS_MAX)
      PARAMETER (NOBLEVMX = P_LEVELS_MAX+1)
      PARAMETER (MAX_NUM_ACOB_FILES=100)
      PARAMETER (NUM_OB_FILE_TYPES = 10)
      INTEGER NOBTYP,NDVHDR,MAXNLEV1,                                   &
     & OBSTYP(NOBTYPMX),NOBLEV(NOBTYPMX),NDATAV(NOBTYPMX),              &
     & NERLEV1(NOBTYPMX),NOBS(NOBTYPMX),OBLEVTYP(NOBTYPMX),             &
     & MDISPOBT(NOBTYPMX),OBS_NO_ST(NOBTYPMX),                          &
     & OBS_REF_YY, OBS_REF_MM, OBS_REF_DD, OBS_REF_HH, OBS_REF_MIN

      REAL           MISSD,                                             &
     &               OBLEVELS(NOBLEVMX,NOBTYPMX),                       &
     &               OBLAYERB(NOBLEVMX+1,NOBTYPMX),                     &
     &               TIMEINT,TIMENEXT,                                  &
     &               OBS_LAT_N, OBS_LAT_S, OBS_LONG_E, OBS_LONG_W

      INTEGER PER_FILE_TNDVMAX(MAX_NUM_ACOB_FILES)

      CHARACTER*30 OB_FILE_TYPE(NUM_OB_FILE_TYPES)

      CHARACTER*256 USED_FILES(MAX_NUM_ACOB_FILES)
      INTEGER FILENAME_LEN(MAX_NUM_ACOB_FILES)
      INTEGER NUM_USED_FILES

      COMMON /COMOBS/ NOBTYP,NDVHDR,MAXNLEV1,                           &
     & OBSTYP,NOBLEV,NDATAV,NERLEV1,NOBS,OBLEVTYP,                      &
     & MDISPOBT,OBS_NO_ST,MISSD,OBLEVELS,OBLAYERB,                      &
     & OBS_REF_YY, OBS_REF_MM, OBS_REF_DD, OBS_REF_HH, OBS_REF_MIN,     &
     & TIMEINT, TIMENEXT,                                               &
     & OBS_LAT_N, OBS_LAT_S, OBS_LONG_E, OBS_LONG_W, PER_FILE_TNDVMAX,  &
     & FILENAME_LEN,NUM_USED_FILES,                                     &
     ! Character variables at the end
     & OB_FILE_TYPE,USED_FILES
!-----------------------------------------------------------------------
