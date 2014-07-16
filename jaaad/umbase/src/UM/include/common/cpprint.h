!LL  Comdeck: CPPRINT -------------------------------------------------
!LL
!LL  Purpose: Hold control variables for point print and field
!LL           increment diagnostics.
!LL
!LL  Model            Modification history:
!LL version  Date
!CL  4.4  28/08/97  Extend to include extra field increment diagnostic
!LL                 variables (not passed by UMUI). R.Rawlins
!LL  ------------------------------------------------------------------


! Switches
      LOGICAL :: LPRVXN      ! min/max monitor or print

      ! min/max print everytime (at PRVXN_STEP intervals)
      LOGICAL :: LPRVXNP

      LOGICAL :: LMOISTP     ! moisture prints in max/min print
      LOGICAL :: LPPRINT     ! patch print around point
      LOGICAL :: LPPRINT_A   ! patch print atmosphere(LPPRINT must be T)
      LOGICAL :: LPPRINT_S   ! patch print surface (LPPRINT must be T)
      LOGICAL :: LVPRINT     ! patch print vertical format
      LOGICAL :: LHPRINT     ! patch print horizontal format
      LOGICAL :: LCLOUDP     ! cloud in min/max print or patch print

      ! radiative heat prints in max/min print  or patch print
      LOGICAL :: LRADP

      LOGICAL :: LTHETAP     ! THETA/T in D1 array
      LOGICAL :: LGLOBALP    ! global/LAM model
      LOGICAL :: LPRFLD      ! model field increments

      REAL :: DUMMY
      REAL :: PPRINT_LAT    !latitude of point at centre of patch print
      REAL :: PPRINT_LONG   !longitude of point at centre of patch print


      INTEGER :: PPRINT_STEP   !timestep interval between patch prints
      INTEGER :: PPRINT_FIRST  !first timestep for patch prints

      ! last timestep for patch prints(0=unlimited)
      INTEGER :: PPRINT_LAST

      INTEGER :: PPRINT_POINT  !point at centre of patch print

      ! tolerance of patch print i.e.no of points each way about centre
      ! point default=1, PPRINT_TOL <=4
      INTEGER :: PPRINT_TOL

      INTEGER :: PRVXN_STEP    !timestep interval between min/max print
      INTEGER :: PRVXN_FIRST   !first timestep for min/max prints

      ! last timestep for min/max prints(0=unlimited)
      INTEGER :: PRVXN_LAST

      INTEGER :: PRVXN_LEVEL ! level to test max/min (<=0 = all levels)
      INTEGER :: PRFLD_STEP  ! timestep interval between increment print
      INTEGER :: PRFLD_FIRST ! first timestep for increment prints
      INTEGER :: PRFLD_LAST  ! last  timestep for increment prints

      COMMON/CPPRINT/LPRVXN,LPRVXNP,LPPRINT,LPPRINT_A,LPPRINT_S,LVPRINT,&
     &  LHPRINT,PPRINT_STEP,PPRINT_POINT,PPRINT_TOL, PPRINT_LAT,        &
     &  PPRINT_LONG,PPRINT_FIRST,PPRINT_LAST,PRVXN_FIRST,PRVXN_LAST,    &
     &  PRVXN_LEVEL,PRVXN_STEP,DUMMY,LCLOUDP,LMOISTP,LRADP,LTHETAP,     &
     &  LGLOBALP,LPRFLD,PRFLD_STEP,PRFLD_FIRST,PRFLD_LAST

      ! Extra field increment diagnostic variables, not passed by UMUI

      ! LOGICAL DEVICE NO.
      INTEGER,PARAMETER:: NDEV_FLD=138 ! for 'CACHED' env. variable

      INTEGER LEN_FLD_FILENAME       ! Filename length of NDEV_FLD
      CHARACTER*80 FLD_FILENAME      ! Filename of NDEV_FLD file
      COMMON/CPPRINT_FLD/LEN_FLD_FILENAME,FLD_FILENAME

! CPPRINT end
