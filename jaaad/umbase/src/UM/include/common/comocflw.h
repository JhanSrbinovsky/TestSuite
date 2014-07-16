! COMOCFLW start
!LL 5.2  23/08/00   MXSCAN and MSCAN moved to UMSCALAR.  D.Storkey
      REAL ::  AREA            ! Surface area of ocean
      REAL ::  VOLUME          ! Volume of ocean
      REAL ::  PNU             ! Time smoothing parameter
      REAL ::  PNU2M           ! (1-2*PNU)
      REAL ::  BBUD
      REAL ::  CCUD
      REAL ::  DDUD
      REAL ::  GGUD
      REAL ::  HHUD
      REAL ::  BBUB
      REAL ::  CCUB
      REAL ::  DDUB
      REAL ::  GGUB
      REAL ::  HHUB
      REAL ::  BBTD
      REAL ::  CCTD
      REAL ::  DDTD
      REAL ::  BBTB
      REAL ::  CCTB
      REAL ::  DDTB

      INTEGER:: NMIX      ! Number of timestems between forward steps
      INTEGER:: NNERGY    ! ----"------------- energy diagnostics
      INTEGER:: NDISKB    ! Pointer to logical "disk" of previous tstep
      INTEGER:: NDISK     ! ----------"----------------- current -"-
      INTEGER:: NDISKA    ! ------------"--------------  next    -"-
      INTEGER:: NDISKB_SAVE
      INTEGER:: LABS(3)   ! "Unit numbers" for logical disks
      INTEGER:: MIX       ! Indicator for "mixing" timestep
      INTEGER:: NERGY     ! Indicator for energy diagnostic timestep
      INTEGER:: NTSI      ! Timesteps between single line output
      INTEGER:: JRPRT(20) ! Rows to be printed during ROWCALC
      INTEGER:: O_ADVECT_SCHEME(2,20)

      COMMON /FULLWD/ AREA,VOLUME,NMIX,NNERGY,                          &
     &  NDISKB,NDISK,NDISKA,NDISKB_SAVE,LABS,MIX,NERGY,                 &
     &  NTSI,PNU,PNU2M,JRPRT,                                           &
     &  O_ADVECT_SCHEME,                                                &
     &  BBUD,CCUD,DDUD,GGUD,HHUD,BBUB,CCUB,DDUB,GGUB,HHUB               &
     &,BBTD,CCTD,DDTD,BBTB,CCTB,DDTB
! COMOCFLW end
