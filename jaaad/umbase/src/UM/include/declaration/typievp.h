! Arguments passed around within EVP module
! Current Code Owner: R Hill
! History:
! Version  Date     Comment
! -------  ----     -------
!   5.5   21/02/03  Added this header + Changes necessitated
!                   by removal of SWAP_1D.     R. Hill
!   6.2   30/03/06  Added L_ONOPOLO for use in EVP_SETUP. R. Hill
! TYPIEVP start
      LOGICAL :: L_OCYCLIC
      LOGICAL :: L_ICESSTILT
      LOGICAL :: L_ICYNPOL
      LOGICAL :: L_ONOPOLO
      INTEGER :: IMT
      INTEGER :: IMTM1
      INTEGER :: IMTM2  !No. points in row, minus 1, and minus 2
      INTEGER :: JMT
      INTEGER :: JMTM1        !No. Tracer (Velocity) rows on processor
      REAL :: ndte             !No. EVP subcycles per timestep
      REAL :: dtts               !Timestep of forcing
      REAL :: dte
      REAL :: dter         !Subcycling timestep, and reciprocal
      REAL :: t_dx(IMT,JMT) !Grid lengths between tracer points & recips
      REAL :: t_dxr(IMT,JMT)
      REAL :: t_dy(IMT,JMT)
      REAL :: t_dyr(IMT,JMT)
      REAL :: cs(JMT)
      REAL :: cst(JMT)
      REAL :: csr(jmt)  ! 1/cs
      REAL :: tng(jmt) !Tan(lat) on U grid (dimension jmt?) see TYPOCONE
      REAL :: radius_si         !Radius of the earth
      REAL :: HT_E(IMT,JMT)     !Length of eastern edge of T-cell
      REAL :: HT_N(IMT,JMT)     !Length of northern edge of T-cell
      REAL :: xymin
      LOGICAL :: ocean(imt,jmt)    !Land-sea mask (Tracer grid)
      LOGICAL :: icy(imt,jmt)      !Ice mask
      LOGICAL :: icy_uv(imt,jmtm1)
      LOGICAL :: icy_uvp(imt,jmtm1)
      REAL :: ocean_t_mask(IMT,JMT)
      REAL :: ocean_uv_mask(IMT,JMTM1)
      REAL :: uice(IMT,JMTM1)
      REAL :: vice(IMT,JMTM1)
      REAL :: aice(IMT,JMT)
      REAL :: hice(IMT,JMT)
      REAL :: hsnow(IMT,JMT)
      REAL :: rhoice
      REAL :: rhosnow
      REAL :: rho_water_si
      REAL :: coriolis(IMT,JMT)
      REAL :: ucurrent(IMT,JMTM1)
      REAL :: vcurrent(IMT,JMTM1)
      REAL :: wsx_ice(IMT,JMTM1)
      REAL :: wsy_ice(IMT,JMTM1)
      REAL :: isx(IMT,JMTM1)
      REAL :: isy(IMT,JMTM1)
      REAL :: sstiltx(IMT,JMTM1)
      REAL :: sstilty(IMT,JMTM1)
      REAL :: waterx(IMT,JMTM1)
      REAL :: watery(IMT,JMTM1)
      REAL :: quad_ice_drag
      REAL :: ecc2              ! squared eccentricity of yield curve
      REAL :: ecci              ! 4./ecc2
      REAL :: eccm             ! 2.-0.5*ecci)
      REAL :: eccp             ! 1.+ecci*0.25
      REAL :: cw             ! cosine of turn. angle in water (phiwater)
      REAL :: sw               ! sine of turn. angle in water (phiwater)
      REAL :: eyc              ! coefficient for calc. the parameter E
      REAL :: cstar            ! c*, parameter in the defn of prss
      REAL :: pstar            ! P*, param in  defn of prss (dyn/cm^2)
      REAL :: etan(imt,jmt)    ! shear viscosity in n triangle of T-cell
      REAL :: etas(imt,jmt)    ! shear viscosity in s triangle of T-cell
      REAL :: etae(imt,jmt)    ! shear viscosity in e triangle of T-cell
      REAL :: etaw(imt,jmt)    ! shear viscosity in w triangle of T-cell
      REAL :: prss(imt,jmt)   ! pressure P (centered in T-cell) (dyn/cm)
      REAL :: umass(imt,jmtm1) ! mass of U-cell (kg)
      REAL :: fm(imt,jmtm1)     ! mass*coriolis param. on U grid
      REAL :: sig11ne(imt,jmt) ! int. ice stress tensor, sigma_11, north
      REAL :: sig11se(imt,jmt)   !   sigma_11, east       (g/s^2)
      REAL :: sig11sw(imt,jmt)   !   sigma_11, south
      REAL :: sig11nw(imt,jmt)   !   sigma_11, west
      REAL :: sig12ne(imt,jmt)   !   sigma_12, north
      REAL :: sig12se(imt,jmt)   !   sigma_12, east
      REAL :: sig12sw(imt,jmt)   !   sigma_12, south
      REAL :: sig12nw(imt,jmt)   !   sigma_12, west
      REAL :: sig22ne(imt,jmt)   !   sigma_22, north
      REAL :: sig22se(imt,jmt)   !   sigma_22, east
      REAL :: sig22sw(imt,jmt)   !   sigma_22, south
      REAL :: sig22nw(imt,jmt)   !   sigma_22  west
      REAL :: a2na(imt,jmt)
      REAL :: a2sa(imt,jmt)
      REAL :: a2ea(imt,jmt)
      REAL :: a2wa(imt,jmt)
      REAL :: b2n(imt,jmt)
      REAL :: b2s(imt,jmt)
      REAL :: b2e(imt,jmt)
      REAL :: b2w(imt,jmt)
      REAL :: t_dx8(imt,jmt)
      REAL :: t_dy8(imt,jmt)
      REAL :: edy(imt,jmt)
      REAL :: edx(imt,jmt)
      REAL :: eHN(imt,jmt)
      REAL :: eHE(imt,jmt)
      REAL :: eHNm(imt,jmt)
      REAL :: eHEm(imt,jmt)
      REAL :: HT_N4(imt,jmt)
      REAL :: HT_E4(imt,jmt)
      REAL :: h2n(imt,jmt)
      REAL :: h2s(imt,jmt)
      REAL :: h2e(imt,jmt)
      REAL :: h2w(imt,jmt)
      REAL :: prssn(imt,jmt)
      REAL :: prsss(imt,jmt)
      REAL :: prsse(imt,jmt)
      REAL :: prssw(imt,jmt)
      REAL :: Tdter
      REAL :: tdamp2
      REAL :: c2n
      REAL :: d2n
      REAL :: acoef
      REAL :: bcoef
      REAL :: ccoef
      REAL :: dcoef
      REAL :: ecoef
      REAL :: dxtr4(imt,jmt)
      REAL :: dytr4(imt,jmt)
! TYPIEVP end
