#if defined(OCEAN) || defined(S40_1A) || defined(S40_2B) \
 || defined(C82_1A)
! TYPOINDX contains all the indices and row-wise loop control variables
! required by the ocean MPP code.
!LL   5.4   19.06.02   Added JMT_MAX for easier dimensioning
!                      of memory aligned arrays. R. Hill
!     6.2  14.02.06   Add to section S40_2B. R. Hill


      INTEGER :: J_1     ! Local value of loop control for J = 1, n
      INTEGER :: J_2     !   "     "    "   "     "     "  J = 2, n
      INTEGER :: J_3     !   "     "    "   "     "     "  J = 3, n
      INTEGER :: J_JMT   !   "     "    "   "     "     "  J = n, JMT
      INTEGER :: J_JMTM1 !   "     "    "   "     "     "  J = n, JMTM1
      INTEGER :: J_JMTM2 !   "     "    "   "     "     "  J = n, JMTM2
      INTEGER :: J_JMTP1 !   "     "    "   "     "     "  J = n, JMTP1
      INTEGER :: JST     ! First row this process considers (no halo)
      INTEGER :: JFIN    ! Last   "    "     "     "        "
      INTEGER :: JMT_NOHALO
      INTEGER :: JMTM1_NOHALO
      INTEGER :: J_FROM_LOC
      INTEGER :: J_TO_LOC
      INTEGER :: JMT_GLOBAL       ! Global value of JMT
      INTEGER :: JMTM1_GLOBAL     ! Global value of JMT - 1
      INTEGER :: JMTM2_GLOBAL     ! Global value of JMT - 2
      INTEGER :: JMTP1_GLOBAL     ! Global value of JMT + 1
      INTEGER :: J_OFFSET         ! Global value of JST - 1
      INTEGER :: O_MYPE           ! MYPE for ocean arg lists
      INTEGER :: O_EW_HALO        ! EW HALO for ocean arg lists
      INTEGER :: O_NS_HALO        ! NS HALO for ocean arg lists
      INTEGER :: J_PE_JSTM1
      INTEGER :: J_PE_JSTM2
      INTEGER :: J_PE_JFINP1
      INTEGER :: J_PE_JFINP2
      INTEGER :: O_NPROC
      INTEGER :: JMT_MAX ! Highest value of JMT on any PE
      INTEGER :: imout(4)
      INTEGER :: jmout(4)! i,j indices for pts in Med outflow
      INTEGER :: J_PE_IND_MED(4)  ! no for each PE in Med outflow
      INTEGER :: NMEDLEV          ! no of levels for Med outflow

      !level at which deep advective Med outflow exits the Mediterranean
      INTEGER :: lev_med

      ! level at which deep advective flow enters the Hudson Bay
      INTEGER :: lev_hud

      INTEGER :: imout_hud(4)  ! zonal index for Hudson Bay outflow
      INTEGER :: jmout_hud(4)  ! merid index for Hudson Bay outflow
      INTEGER :: J_PE_IND_HUD(4)  ! PE's involved in HB outflow

      ! last level for which there is inflow toMediterranean
      INTEGER :: med_topflow
! TYPOINDX end
#endif
