#if defined(OCEAN) || defined(S40_1A) || defined(S40_2B) \
 || defined(C82_1A)
!========================== COMDECK COCNINDX ==========================
!
! Description:
!
!       This comdeck contains all the indices and row-wise loop
!       control variables required by the ocean MPP code.
!
      ! Note: All variables prefixed "J_" contain values which
      ! take account of halo sizes. Eg: for the 3 row domain defined
      ! by JST = 10 and JFIN = 12, with a halo of 2 rows, J_1
      ! will be 3, J_JMT will be 5.

      INTEGER :: J_1     ! Local value of loop control for J = 1, n
      INTEGER :: J_2     !   "     "    "   "     "     "  J = 2, n
      INTEGER :: J_3     !   "     "    "   "     "     "  J = 3, n
      INTEGER :: J_JMT   !   "     "    "   "     "     "  J = n, JMT
      INTEGER :: J_JMTM1 !   "     "    "   "     "     "  J = n, JMTM1
      INTEGER :: J_JMTM2 !   "     "    "   "     "     "  J = n, JMTM2
      INTEGER :: J_JMTP1 !   "     "    "   "     "     "  J = n, JMTP1
      INTEGER :: JMT_NOHALO  ! New variable needed for SWAP_BOUNDS!
      INTEGER :: JMTM1_NOHALO
      INTEGER :: JST     ! First row this process considers (no halo)
      INTEGER :: JFIN    ! Last   "    "     "     "        "
      INTEGER :: J_FROM_LOC       ! Local value of start index
      INTEGER :: J_TO_LOC         ! Local value of end index
      INTEGER :: JMT_GLOBAL       ! Global value of JMT
      INTEGER :: JMTM1_GLOBAL     ! Global value of JMT - 1
      INTEGER :: JMTM2_GLOBAL     ! Global value of JMT - 2
      INTEGER :: JMTP1_GLOBAL     ! Global value of JMT + 1
      INTEGER :: J_OFFSET         ! Start row - 1
      INTEGER :: O_MYPE           ! MYPE value in arg lists for ocean
      INTEGER :: O_EW_HALO        ! EW_HALO for ocean arg lists
      INTEGER :: O_NS_HALO        ! NS_HALO for ocean arg lists
      INTEGER :: J_PE_JSTM1
      INTEGER :: J_PE_JSTM2
      INTEGER :: J_PE_JFINP1
      INTEGER :: J_PE_JFINP2
      INTEGER :: O_NPROC          ! NPROC for ocean
      INTEGER :: imout(4),jmout(4)! i,j indices for pts in Med outflow
      INTEGER :: J_PE_IND_MED(4)  ! no for each PE in Med outflow
      INTEGER :: NMEDLEV          ! no of levels for Med outflow
      INTEGER :: lev_med  ! level at which deep flow from Med occurs
      INTEGER :: lev_hud  ! level at which deep flow into Hudson Bay
      INTEGER :: imout_hud(4),jmout_hud(4)  ! Hudson Bay i,j
      INTEGER :: J_PE_IND_HUD(4)  ! PE's involved in Hudson Bay outflow

      ! last level for which there is inflow to Mediterranean
      INTEGER :: med_topflow
      INTEGER :: JMT_MAX ! Highest value of JMT on any PE
      REAL :: LAT_CHECK  ! Save NSSPACEO for consistency checks.

      COMMON /COCNINDX/                                                 &
     &  J_1, J_2, J_3,                                                  &
     &  J_JMT, J_JMTM1, J_JMTM2, J_JMTP1,                               &
     &  JST, JFIN, JMT_NOHALO,JMTM1_NOHALO, J_FROM_LOC, J_TO_LOC,       &
     &  JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL,                         &
     &  JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO,           &
     &  J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2,               &
     &  O_NPROC,JMT_MAX,                                                &
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,                               &
     &  lev_med,lev_hud,imout_hud,jmout_hud,J_PE_IND_HUD,med_topflow    &
     &  ,LAT_CHECK

! COCNINDX end
#endif
