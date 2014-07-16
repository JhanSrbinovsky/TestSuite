! TYPOC3DG Define 3D variables which are to be passed through all
! subroutine argument lists from ROW_CTL to ROWCALC. Arguments being
! placed in this file must also be included in the call to BLOKCALC
! from ROW_CTL. (Typically, these will take the form of pointers to the
! STASHWORK array).
! History:
! Version  Date     Comment
! -------  ----     -------
!   6.2  14.07.05   Added items for MLD and ILD diagnostics

      INTEGER :: imt_stash           ! Number of distinct columns
      REAL :: dtrc(imt_stash,jmt,km) ! Total heating rate in TRACER (K/s)
      LOGICAL :: sfrc                ! Stash flag for dtrc.
      REAL::UISOP_OUT(imt_stash,jmt,km)  ! u* isopycnal velocity
      REAL::VISOPN_OUT(imt_stash,jmt,km) ! v* at north face of T gridbox
      REAL::WISOP_OUT(imt_stash,jmt,km-1)  ! w* at top face of T gridbox
      REAL::DTGM_OUT(imt_stash,jmt,km)  ! Heating rate due to GM scheme
      REAL::DSGM_OUT(imt_stash,jmt,km)  ! dsalinity/dt due to GM scheme
!     Stash flags for the above
      LOGICAL :: SF_UISOP
      LOGICAL :: SF_VISOP
      LOGICAL :: SF_WISOP
      LOGICAL :: SF_DTGM
      LOGICAL :: SF_DSGM
      REAL :: gnum_OUT(imt_stash,jmt,kmm1)
      REAL :: gnuT_OUT(imt_stash,jmt,kmm1)
      REAL :: Rim_OUT(imt_stash,jmt,kmm1)
      REAL :: RiT_OUT(imt_stash,jmt,kmm1)
      REAL :: HM_OUT(IMT_STASH,JMT)
      REAL :: HT_OUT(IMT_STASH,JMT)
      REAL :: LM_OUT(IMT_STASH,JMT)
      REAL :: LT_OUT(IMT_STASH,JMT)
      REAL :: RIMLDCALC_OUT(IMT_STASH,JMT,KMM1)
      LOGICAL :: SF_LM
      LOGICAL :: SF_LT
      LOGICAL :: SF_MLDCALC
      LOGICAL :: SF_gnum
      LOGICAL :: SF_gnuT
      LOGICAL :: SF_Rim
      LOGICAL :: SF_RiT
      LOGICAL :: SF_hm
      LOGICAL :: SF_hT
      LOGICAL :: SFUTOT
      LOGICAL :: SFVTOT
      LOGICAL :: SFTEMP
      LOGICAL :: SFMLD_DENS
      LOGICAL :: SFMLD_TEMP

      REAL :: UTOT(IMT_STASH,JMT,KM)
      REAL :: VTOT(IMT_STASH,JMT,KM)
      REAL :: TEMPERATURE(IMT_STASH,JMT,KM)
      REAL :: MLD_DENS(IMT_STASH,JMT)
      REAL :: MLD_TEMP(IMT_STASH,JMT)
! TYPOC3DG end
