#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Read run-time control information from namelists for atmos model
!
! Subroutine Interface:
      SUBROUTINE Readlsta(                                              &
#include "arglndm.h"
     &Dummyarg)
!
! If new version of radiation code is required then we
! need to use the modules for improved time-stepping
!
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)

      USE CORADOCA
      USE SW_CONTROL_STRUCT
      USE LW_CONTROL_STRUCT
#endif
      USE CSENARIO_MOD

      Use cv_run_mod                         ! Access to all variables


! module for RUN_LAND namelist
      USE LAND_SURF_MOD, ONLY :                                         &
     & FRAC_SNOW_SUBL_MELT                                              &
     &,MASKD                                                            &
     &,SOILHC_METHOD                                                    &
     &,ALL_TILES                                                        &
     &,L_VG_SOIL                                                        &
     &,RUN_LAND

!
      USE BL_OPTION_MOD, ONLY : WeightLouisToLong
!
      IMPLICIT NONE
!
! Description: Read run-time control information passed as namelists
!   from the UMUI, as required to set up parametrization constants and
!   logical switches needed by physics and dynamics schemes for the
!   Atmosphere model.
! Method:  Sequential read of namelists. Note that defaults are set to
!   missing data indicators. Namelist variables/declarations are held in
!   cruntimc.h file.
!
! Current Code Owner: Rick Rawlins
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component: Control Atmos
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):
#include "cmaxsize.h"
#include "parvars.h"
#include "typsize.h"
#include "typcona.h" 
#include "typlndm.h"
#include "nstypes.h"
#include "natforce.h"

#include "c_gwave.h" 

#include "csubmodl.h"
#include "chsunits.h"
#include "ccontrol.h"
#include "cruntimc.h"
#include "cpprint.h"
#include "c_mdi.h"
#include "c_pi.h"
#include "swopt3a.h"
#include "lwopt3a.h"
#if defined(A01_3A) || defined(A02_3A)
#include "swcopt3a.h"
#include "lwcopt3a.h"
#endif
#include "ctlnl3a.h"
#include "blopt8a.h"
#include "cprintst.h"
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
#include "satopt.h"
#endif
!
! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER                                                           &
     & Dummyarg    ! Not used - required for legal routine call
!   Array  arguments with intent(in):
!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
!   Scalar arguments with intent(out):
!   Array  arguments with intent(out):

! Local parameters:
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='Readlsta')

! Local scalars:
      INTEGER                                                           &
     & ErrorStatus      ! Return code : 0 Normal Exit : >0 Error
      CHARACTER*256                                                     &
     & CMessage         ! Error message if Errorstatus >0


      INTEGER                                                           &
     & LEVEL                                                            &
     &,J                                                                &
     &,JJ
! Local dynamic arrays:
! Function & Subroutine calls:
      External SET_H_SECT,Ereport

! Extra namelists:
#include "cv_run_nml.h"

      NAMELIST/PPRINTXN/                                                &
     &               LPRVXN,LPRVXNP,LPPRINT,LPPRINT_A,LPPRINT_S,LVPRINT,&
     &               LHPRINT,PPRINT_STEP,PPRINT_POINT,PPRINT_TOL,       &
     &               PPRINT_FIRST,PPRINT_LAST,                          &
     &               PRVXN_FIRST,PRVXN_LAST,PRVXN_LEVEL,                &
     &               PRVXN_STEP,LMOISTP,                                &
     &               LPRFLD,PRFLD_STEP,PRFLD_FIRST,PRFLD_LAST

      NAMELIST / CLMCHFCG /  CLIM_FCG_NYEARS, CLIM_FCG_YEARS,           &
     &                       CLIM_FCG_LEVLS,  CLIM_FCG_RATES,           &
                             L_CLMCHFCG
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      NAMELIST / RADFCDIA / C2C_O2, C2C_O3, C2C_CO2, C2C_N2O, C2C_CH4,  &
     &  C2C_CFC11, C2C_CFC12, C2C_C113, C2C_HCFC22, C2C_HFC125,         &
     &  C2C_HFC134, C2C_AEROSOL, C2C_SULPC_D, C2C_SEAS_D, C2C_SOOT_D,   &
     &  C2C_BMB_D, C2C_OCFF_D, C2C_SUN, C2C_VOL, C2C_LAND_S, C2C_ALL,   &
     &  C2C_WMG
#endif

      ErrorStatus = 0

!  Set assimilation logicals from character variable.
      L_AC        = A_ASSIM_MODE  ==  "AC"

!  Set H_SECT to values found in ATMOS_SR in COMDECK MODEL
!  Needs tidied at a later version to be done in a less messy way.
! DEPENDS ON: set_h_sect
      CALL SET_H_SECT(H_SECT,MAXSECTS)

! Read atmosphere run consts. control data into COMMON

! Boundary layer physics
      DO level=1,max_number_alpha_cds
         alpha_Cd (level) = RMDI
      ENDDO
      L_SBLeq             = .FALSE.
      L_SBLco             = .TRUE.
      Muw_SBL             = 1.5
      Mwt_SBL             = 1.0
      SBL_OP              = Long_tails
      ISHEAR_BL           = ON
      NG_STRESS           = OFF
      DECFIX              = ON
      STOPWE_SBL          = ON
      TRWEIGHTS1          = ON
      FLUX_GRAD           = Locketal2000
      SeaSalinityFactor   = 1.0
      ISeaZ0T             = Fixed_Z0T
      ISeaDynDiag         = OFF
      COR_UST             = ON
      COR_MO_ITER         = OFF
      NON_LOCAL_BL        = ON
      Buddy_sea           = OFF
      FORMDRAG            = Explicit_stress
      OROG_DRAG_PARAM     = 0.3
      FD_stab_dep         = ON
      l_use_bl_diag_term  = .false.
      LOCAL_FA            = OFF
      PRANDTL             = Constant_SBL
      Keep_Ri_FA          = OFF
      NL_BL_LEVELS        = OFF
      DO J = 1, 20
        BL_OPTIONS(J)     = OFF
      ENDDO
! STPH_RP Boundary Layer physics
      par_mezcla_max = 0.5
      par_mezcla     = 0.15
      par_mezcla_min = 0.05
      G0_max = 20.0
      G0_RP  = 10.0
      G0_min = 5.0
      Charnock_max = 0.026
      Charnock_min = 0.01
! Large scale precipitation physics
      DO level=1,max_model_levels
        RHCRIT(level) = RMDI
      ENDDO
      l_droplet_settle=.false.
! STPH_RP Large scale precipitation physics
       RHCRIT_max = 0.910
       RHCRIT_min = 0.875
       CI_max = 33.0
       CI_min = 17.0
       M_CI_max = 1.4
       M_CI     = 1.0
       M_CI_min = 0.6

! Radiation physics
      l_ovrlap            = .FALSE.  ! Allow Convective and LS Fractions
                                     ! to overlap
      l_ccw_scav          = .FALSE.  ! Maintain condensate in overlapped
                                     ! Cloud fractions (if l_ovrlap=T)
      l_emis_land_gen     = .FALSE.  ! Aggregate Land surface emissivity
      emis_land_gen       = 1.0      ! Aggregate Land surface emissivity
      L_climat_aerosol    = .FALSE.
      L_HadGEM1_Clim_Aero = .FALSE.
      L_SEC_VAR           = .FALSE.
      L_EqT               = .FALSE.
      L_cldtop_t_fix      = .FALSE.
      A_SW_SEGMENTS       = IMDI
      A_SW_SEG_SIZE       = -99     ! -ive is disabled
      A_LW_SEGMENTS       = IMDI
      A_LW_SEG_SIZE       = -99     ! -ive is disabled
      CO2_MMR             = RMDI
      INHOM_CLOUD_SW      = RMDI
      INHOM_CLOUD_LW      = RMDI
      DP_CORR_STRAT       = RMDI
      DP_CORR_CONV        = RMDI
!   (Note: some radiation variables are defined in rad_com.h file:)
      ALPHAC = RMDI
      ALPHAM = RMDI
      DTICE  = RMDI
      DT_BARE = RMDI
      SW_BETA = RMDI
      DALB_BARE_WET = RMDI
      PEN_RAD_FRAC  = RMDI
      N2OMMR = RMDI
      CH4MMR = RMDI
      C11MMR = RMDI
      C12MMR = RMDI
      O2MMR  = RMDI
      C113MMR    = RMDI
      HCFC22MMR  = RMDI
      HFC125MMR  = RMDI
      HFC134AMMR = RMDI
      IS_NCOL = IMDI

      CAN_RAD_MOD=1
      ILAYERS=10

!     In common block LWOPT3A
      L_EXTRA_TOP_LW=.FALSE.
      L_LOCAL_CNV_PARTITION_LW = .FALSE.

!     In common block SWOPT3A
      L_LOCAL_CNV_PARTITION_SW = .FALSE.
      L_EXTRA_TOP_SW=.FALSE.

! Gravity wave drag physics
      KAY_GWAVE           = RMDI
      GWD_FRC             = RMDI
      GWD_FSAT            = RMDI
      L_TAUS_SCALE        = .FALSE.
      L_FIX_GWSATN        = .FALSE.
      L_gwd_40km          = .FALSE.
      L_USSP_OPAQUE       = .FALSE.
      SAT_SCHEME          = Stress_saturation
! STPH_RP Gravity Wave Drag
       GWD_FRC_max = 6.0
       GWD_FRC_min = 2.0
       KAY_GWAVE_max = 4.4E+03
       KAY_GWAVE_min = 2.5E+03

! River routing defaults
      RIVER_VEL    = RMDI
      RIVER_MCOEF  = RMDI

! FILE natural forcing defaults
      FILE_SCVARY  = ''
      FILE_VOLCTS  = ''

!     **********    Boundary layer solver defaults   *********
      L_us_blsol = .false.
      Puns = 0.5
      Pstb = 2.0
      WeightLouisToLong_in = 0.0

!     **********    Dynamics defaults   *********

      L_regular = .true.
      var_ratio = 1.0
      lam_ratio = 1.0
      phi_ratio = 1.0
      lam_var = 0
      lam_frac = 0.5
      phi_var = 0
      phi_frac = 0.5

      n_rims_to_do = 1
      L_LBC_balance = .FALSE.
      L_lbc_new = .FALSE.

      L_conserve_tracers=.TRUE.

! Dynamics MDI defaults
      ramp_lat_radians = RMDI
      L_Physics = .true.
      L_Backwards = .false.
      L_Run_With_Physics2 = .true.
      L_perturb_IC_theta = .false.
      IntRand_seed=0

! Defaults for iterative semi-Lagrangian scheme (cycling)
      NumCycles = 1
      L_new_tdisc = .false.
      extrp_weight = 1.5
      GCR_tol_abs2 = 1.0E-02
      GCR_tol_res2 = 1.0E-07
      L_GCR_cycle_opt = .false.
      alpha_1_2 = 0.6
      alpha_2_2 = 1.0
      alpha_3_2 = 0.6
      alpha_4_2 = 1.0
! Defaults for theta advection method
      L_fint_theta = .false.
! Diffusion defaults
      L_filter = .false.
      L_filter_incs = .false.
      L_diff_auto = .false.
      max_sweeps = 8
      ref_lat_deg = 0.0
      diff_coeff_ref = 0.25
      vdiffuv_test = 100.0
      L_vdiff_uv = .false.
      vdiffuv_start = model_levels
      vdiffuv_end = model_levels
      L_adjust_theta = .false.
      adjust_theta_start = model_levels
      adjust_theta_end = model_levels
      vdiffuv_timescale = 1
      L_upper_ramp = .false.
      top_filt_start = model_levels + 1
      top_filt_end = model_levels
      top_diff = 0.1
      up_diff_scale = 0.5
      L_pofil_new = .false.
      L_sponge = .false.
      sponge_ew = 0
      sponge_ns = 0
      sponge_power = 1

      tar_horizontal = 0

!   thmono limiter default
      thmono_height = 0.0
      thmono_levels = 0


!     **********    Dynamics defaults   *********

!   Diagnostic printing defaults
      L_print_pe = .false.
      L_flush6 = .false.
      L_diag_print = .false.
      L_diag_print_ops = .false.
      L_print_w = .false.
      L_print_wmax = .false.
      L_print_div = .false.
      L_print_lapse = .false.
      L_print_theta1 = .false.
      L_print_max_wind = .false.
      L_print_shear = .false.
      L_diag_L2norms = .false.
      L_diag_L2helm = .false.
      print_step = 1
      diag_interval = 1
      first_norm_print = 1
      w_print_limit = 1000.0
      norm_lev_start = 1
      norm_lev_end = model_levels
      L_print_pe = .false.
      L_diag_wind = .false.
      L_diag_noise = .false.

!     In common block CDERIVED (in comdeck CCONSTS called in TYPCONA)

! Default settings for h_split (can be overridden by UMUI-supplied
!  values if this is required at a future release, when rmdi defaults
!  should be reinstated):
!  low:middle:high cloud model levels =(1)->(2):(2)->(3):(3)->(4)
      h_split(1) =   111.        ! ICAO 1000mb height (m)
      h_split(2) =  1949.        ! ICAO  800mb height (m)
      h_split(3) =  5574.        ! ICAO  500mb height (m)
      h_split(4) = 13608.        ! ICAO  150mb height (m)
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! Defaults for things to change in diagnostic 2nd call to radiation.
      DO J=1, C2C_size
        C2C_ALL(J) = .FALSE.
        C2C_WMG(J) = .FALSE.
        C2C_O2(J) = .FALSE.
        C2C_O3(J) = .FALSE.
        C2C_CO2(J) = .FALSE.
        C2C_N2O(J) = .FALSE.
        C2C_CH4(J) = .FALSE.
        C2C_CFC11(J) = .FALSE.
        C2C_CFC12(J) = .FALSE.
        C2C_C113(J) = .FALSE.
        C2C_HCFC22(J) = .FALSE.
        C2C_HFC125(J) = .FALSE.
        C2C_HFC134(J) = .FALSE.
        C2C_AEROSOL(J) = .FALSE.
        C2C_SULPC_D(J) = .FALSE.
        C2C_SEAS_D(J) = .FALSE.
        C2C_SOOT_D(J) = .FALSE.
        C2C_BMB_D(J) = .FALSE.
        C2C_OCFF_D(J) = .FALSE.
        C2C_SUN(J) = .FALSE.
        C2C_VOL(J) = .FALSE.
        C2C_LAND_S(J) = .FALSE.
      ENDDO
#endif

      DO J = 1,MAX_REQ_THPV_LEVS
        REQ_THETA_PV_LEVS(J) = RMDI
      ENDDO

!=====================================================================
!     Set defaults for RUN_Convection namelist
!=====================================================================
#include "cv_run_nml_data.h"


!  Read in Physics/Dynamics namelists

      READ (5,RUN_BL)           ! Boundary layer physics
      READ (5,RUN_BLICE)      ! Sea ice roughness characteristics
      READ (5,RUN_BLVEG)    ! Surface type characteristics
!     ! Store BL options in integer array for ease of passing around
!     ! Changes here must be duplicated in S_SHELL for the SCM
      BL_OPTIONS(1) = ISHEAR_BL
      BL_OPTIONS(2) = NG_STRESS
      BL_OPTIONS(3) = DECFIX
      BL_OPTIONS(4) = STOPWE_SBL
      BL_OPTIONS(5) = TRWEIGHTS1
      BL_OPTIONS(6) = FLUX_GRAD
      BL_OPTIONS(7) = SBL_OP
      BL_OPTIONS(8) = ISeaZ0T
      BL_OPTIONS(9) = ISeaDynDiag
      BL_OPTIONS(10) = COR_UST
      BL_OPTIONS(11) = NON_LOCAL_BL
      BL_OPTIONS(12) = LOCAL_FA
      BL_OPTIONS(13) = PRANDTL
      BL_OPTIONS(14) = Buddy_sea
      BL_OPTIONS(15) = I_SCRN_T_DIAG
      BL_OPTIONS(16) = COR_MO_ITER
      BL_OPTIONS(17) = Keep_Ri_FA
      BL_OPTIONS(18) = NL_BL_LEVELS
      BL_OPTIONS(19) = FD_stab_dep
!     Copy WeightLouisToLong from the namelist to the module
!     variable. This is a temporary measure to minimize the
!     number of changes required at this release. At a future
!     release everything should be based on the module.
      WeightLouisToLong = WeightLouisToLong_in

      READ (5,RUN_PFT)           ! Surface parameters

      READ (5,RUN_LAND)         ! land-surface parameters

      READ (5,RUN_Precip)       ! Large scale precipitation physics

      READ (5,RUN_Cloud)        ! Large scale cloud physics

      READ (5,RUN_Convection)   ! Convection physics

      READ (5,RUN_Radiation)    ! Radiation physics

      READ (5,RUN_GWD)          ! Gravity wave drag physics
!
      READ (5,RUN_Aerosol)       ! Aerosol Modelling

      READ (5,RUN_UKCA)       ! UKCA Sub-model

      READ (5,RUN_Dyn)          ! Generalised integration and GCR
                                ! dynamics
      READ (5,RUN_SL)           ! Semi-Lagrangian advection dynamics

      READ (5,RUN_Diffusion)    ! Diffusion, divergence damping and
                                ! filtering dynamics

      READ (5,RUN_RIVERS)       ! River routing
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      READ (5,RADFCDIA)         ! Diagnostic double call to radiation
#endif

      IF ( L_VOLCTS .or. L_SCVARY ) then
         READ (5,FILENATFORCE)     ! filenames for solar/volcanic forcing
      ENDIF

!------------------------------------------------------- 
! Check that if L_CCRAD is selected certain other
! switches are not set so that they conflict.
!------------------------------------------------------- 
! Error capture 
!------------------------------------------------------- 
                 
      If (L_ccrad) Then 

        If (.NOT. L_3D_CCA) Then 
          ErrorStatus = 100 
          CMessage    = '**ERROR**: CCRad is not yet available without'// &
                                  ' the anvil scheme (L_3D_CCA = .True.)'
! DEPENDS ON: ereport 
          CALL ereport(RoutineName, ErrorStatus, CMessage) 
        End If

        If (l_convcld_hadgem1) Then 
          ErrorStatus = 100 
          CMessage    = '**ERROR**: L_CCRad and l_convcld_hadgem1'//      &
                                  ' should not be both set to true.'
! DEPENDS ON: ereport 
          CALL ereport(RoutineName, ErrorStatus, CMessage) 
        End If

        If (l_fix_udfactor) Then 
          ErrorStatus = 100 
          CMessage    = '**ERROR**: L_CCRad and l_fix_udfactor'//         &
                                  ' should not be both set to true.'
! DEPENDS ON: ereport 
          CALL ereport(RoutineName, ErrorStatus, CMessage) 
        End If

        If (l_pc2_diag_sh) Then 
          ErrorStatus = 100 
          CMessage    = '**ERROR**: L_CCRad and l_pc2_diag_sh'//          &
                                  ' should not be both set to true.'
! DEPENDS ON: ereport 
          CALL ereport(RoutineName, ErrorStatus, CMessage) 
        End If

      End If       ! l_ccrad 

!------------------------------------------------------- 
! End error capture 
!-------------------------------------------------------

      IF (L_Backwards) L_Physics = .false.

      IF(PrintStatus >= PrStatus_Normal) THEN

        WRITE(6,*) RoutineName,':Atmosphere run-time constants:-'
        WRITE(6,RUN_BLICE)
        WRITE(6,RUN_BLVEG)
        WRITE(6,RUN_BL)
        WRITE(6,RUN_PFT)
        WRITE(6,RUN_LAND)

        WRITE(6,RUN_Precip)
        WRITE(6,RUN_Cloud)
        WRITE(6,RUN_Convection)
        WRITE(6,RUN_Radiation)
        WRITE(6,RUN_GWD)
        WRITE(6,RUN_Aerosol)
        WRITE(6,RUN_UKCA)
        WRITE(6,RUN_Dyn)
        WRITE(6,RUN_SL)
        WRITE(6,RUN_Diffusion)
        WRITE(6,RUN_RIVERS)
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      WRITE(6,RADFCDIA)
#endif
        WRITE(6,FILENATFORCE)

      ENDIF ! PrintStatus



! Convert from degrees to radians
      polar_filter_north_lat_limit= polar_filter_north_lat_limit*pi/180.
      polar_filter_south_lat_limit= polar_filter_south_lat_limit*pi/180.
      polar_filter_lat_limit      = polar_filter_lat_limit      *pi/180.
      polar_filter_step_per_sweep = polar_filter_step_per_sweep *pi/180.
      ramp_lat_radians = ramp_lat_radians *pi/180.



!  Multiply input req_theta_pv_levs by 100. to convert mb to pascals.
      DO LEVEL = 1,THETA_PV_P_LEVS
        REQ_THETA_PV_LEVS(LEVEL) = REQ_THETA_PV_LEVS(LEVEL)*100.
      END DO

!  Read atmosphere patch print and min/max controls into COMMON/CPPRINT/

      READ(5,PPRINTXN)


!
!     Read controlling information for version 3A of the SW or LW
!     radiation schemes.
!
#if defined(A01_3A) || defined(A02_3A)
      IF (H_SECT(1) == '03A') THEN
!        Options for the shortwave
         READ(5, R2SWCLNL)
      ENDIF
!
      IF (H_SECT(2) == '03A') THEN
!        Options for the longwave
         READ(5, R2LWCLNL)
      ENDIF
#endif
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      IF ((H_SECT(1) == '03C').OR.(H_SECT(1) == '03Z')) THEN
        READ(5, R2SWNCAL)

        DO J=1,N_SWCALL
! DEPENDS ON: sw_control_default
          CALL SW_CONTROL_DEFAULT(SW_CONTROL(J))

!         Initialize data to safe values

          SPECTRAL_FILE_SW=SW_CONTROL(J)%SPECTRAL_FILE
          FIRST_BAND_SW=SW_CONTROL(J)%FIRST_BAND
          LAST_BAND_SW=SW_CONTROL(J)%LAST_BAND
          I_2STREAM_SW=SW_CONTROL(J)%I_2STREAM
          I_GAS_OVERLAP_SW=SW_CONTROL(J)%I_GAS_OVERLAP
          I_CLOUD_SW=SW_CONTROL(J)%I_CLOUD
          I_CLOUD_REPRESENTATION_SW=SW_CONTROL(J)%I_CLOUD_REPRESENTATION
          I_SOLVER_SW=SW_CONTROL(J)%I_SOLVER
          L_O2_SW=SW_CONTROL(J)%L_O2
          I_ST_WATER_SW=SW_CONTROL(J)%I_ST_WATER
          I_CNV_WATER_SW=SW_CONTROL(J)%I_CNV_WATER
          I_ST_ICE_SW=SW_CONTROL(J)%I_ST_ICE
          I_CNV_ICE_SW=SW_CONTROL(J)%I_CNV_ICE
          L_LOCAL_CNV_PARTITION_SW=SW_CONTROL(J)%L_LOCAL_CNV_PARTITION
          I_ANGULAR_INTEGRATION_SW=SW_CONTROL(J)%I_ANGULAR_INTEGRATION
          I_SPH_ALGORITHM_SW=SW_CONTROL(J)%I_SPH_ALGORITHM
          N_ORDER_PHASE_SOLAR_SW=SW_CONTROL(J)%N_ORDER_PHASE_SOLAR
          L_EULER_TRNF_SW=SW_CONTROL(J)%L_EULER_TRNF
          I_TRUNCATION_SW=SW_CONTROL(J)%I_TRUNCATION
          LS_GLOBAL_TRUNC_SW=SW_CONTROL(J)%LS_GLOBAL_TRUNC
          MS_MIN_SW=SW_CONTROL(J)%MS_MIN
          MS_MAX_SW=SW_CONTROL(J)%MS_MAX
          ACCURACY_ADAPTIVE_SW=SW_CONTROL(J)%ACCURACY_ADAPTIVE
          LS_BRDF_TRUNC_SW=SW_CONTROL(J)%LS_BRDF_TRUNC
          L_HENYEY_GREENSTEIN_PF_SW=SW_CONTROL(J)%L_HENYEY_GREENSTEIN_PF
          I_SPH_MODE_SW=SW_CONTROL(J)%I_SPH_MODE

!         Satellite Data:
          L_SUBSAMPLE=SW_CONTROL(J)%L_SUBSAMPLE
          L_GEOSTATIONARY=SW_CONTROL(J)%L_GEOSTATIONARY
          SAT_DESC=SW_CONTROL(J)%SAT_DESC
          SAT_HGT=SW_CONTROL(J)%SAT_HGT
          SAT_LON=SW_CONTROL(J)%SAT_LON
          SAT_LAT=SW_CONTROL(J)%SAT_LAT
          MAX_VIEW_LON=SW_CONTROL(J)%MAX_VIEW_LON
          MIN_VIEW_LON=SW_CONTROL(J)%MIN_VIEW_LON
          MAX_VIEW_LAT=SW_CONTROL(J)%MAX_VIEW_LAT
          MIN_VIEW_LAT=SW_CONTROL(J)%MIN_VIEW_LAT

!         Options for the shortwave
          READ(5, R2SWCLNL)

!         Transfer values which we allow the user to change
!         from the namelist to the structure.

          SW_CONTROL(J)%SPECTRAL_FILE=SPECTRAL_FILE_SW
          SW_CONTROL(J)%FIRST_BAND=FIRST_BAND_SW
          SW_CONTROL(J)%LAST_BAND=LAST_BAND_SW
          SW_CONTROL(J)%I_2STREAM=I_2STREAM_SW
          SW_CONTROL(J)%I_GAS_OVERLAP=I_GAS_OVERLAP_SW
          SW_CONTROL(J)%I_CLOUD=I_CLOUD_SW
          SW_CONTROL(J)%I_CLOUD_REPRESENTATION=I_CLOUD_REPRESENTATION_SW
          SW_CONTROL(J)%I_SOLVER=I_SOLVER_SW
          SW_CONTROL(J)%L_O2=L_O2_SW
          SW_CONTROL(J)%I_ST_WATER=I_ST_WATER_SW
          SW_CONTROL(J)%I_CNV_WATER=I_CNV_WATER_SW
          SW_CONTROL(J)%I_ST_ICE=I_ST_ICE_SW
          SW_CONTROL(J)%I_CNV_ICE=I_CNV_ICE_SW
          SW_CONTROL(J)%L_LOCAL_CNV_PARTITION=L_LOCAL_CNV_PARTITION_SW
          SW_CONTROL(J)%I_ANGULAR_INTEGRATION=I_ANGULAR_INTEGRATION_SW
          SW_CONTROL(J)%I_SPH_ALGORITHM=I_SPH_ALGORITHM_SW
          SW_CONTROL(J)%N_ORDER_PHASE_SOLAR=N_ORDER_PHASE_SOLAR_SW
          SW_CONTROL(J)%L_EULER_TRNF=L_EULER_TRNF_SW
          SW_CONTROL(J)%I_TRUNCATION=I_TRUNCATION_SW
          SW_CONTROL(J)%LS_GLOBAL_TRUNC=LS_GLOBAL_TRUNC_SW
          SW_CONTROL(J)%MS_MIN=MS_MIN_SW
          SW_CONTROL(J)%MS_MAX=MS_MAX_SW
          SW_CONTROL(J)%ACCURACY_ADAPTIVE=ACCURACY_ADAPTIVE_SW
          SW_CONTROL(J)%LS_BRDF_TRUNC=LS_BRDF_TRUNC_SW
          SW_CONTROL(J)%L_HENYEY_GREENSTEIN_PF=L_HENYEY_GREENSTEIN_PF_SW
          SW_CONTROL(J)%I_SPH_MODE=I_SPH_MODE_SW

!         Satellite Data:
          SW_CONTROL(J)%L_SUBSAMPLE=L_SUBSAMPLE
          SW_CONTROL(J)%L_GEOSTATIONARY=L_GEOSTATIONARY
          SW_CONTROL(J)%SAT_DESC=SAT_DESC
          SW_CONTROL(J)%SAT_HGT=SAT_HGT
          SW_CONTROL(J)%SAT_LON=SAT_LON
          SW_CONTROL(J)%SAT_LAT=SAT_LAT
          SW_CONTROL(J)%MAX_VIEW_LON=MAX_VIEW_LON
          SW_CONTROL(J)%MIN_VIEW_LON=MIN_VIEW_LON
          SW_CONTROL(J)%MAX_VIEW_LAT=MAX_VIEW_LAT
          SW_CONTROL(J)%MIN_VIEW_LAT=MIN_VIEW_LAT
        ENDDO
      ENDIF
      IF ((H_SECT(2) == '03C').OR.(H_SECT(2) == '03Z')) THEN
        READ(5, R2LWNCAL)
        DO J=1, N_LWCALL
! DEPENDS ON: lw_control_default
          CALL LW_CONTROL_DEFAULT(LW_CONTROL(J))

!         Initialize data to safe values

          SPECTRAL_FILE_LW=LW_CONTROL(J)%SPECTRAL_FILE
          FIRST_BAND_LW=LW_CONTROL(J)%FIRST_BAND
          LAST_BAND_LW=LW_CONTROL(J)%LAST_BAND
          I_2STREAM_LW=LW_CONTROL(J)%I_2STREAM
          L_IR_SOURCE_QUAD_LW=LW_CONTROL(J)%L_IR_SOURCE_QUAD
          I_GAS_OVERLAP_LW=LW_CONTROL(J)%I_GAS_OVERLAP
          I_CLOUD_LW=LW_CONTROL(J)%I_CLOUD
          I_CLOUD_REPRESENTATION_LW=LW_CONTROL(J)%I_CLOUD_REPRESENTATION
          I_SOLVER_LW=LW_CONTROL(J)%I_SOLVER
          L_N2O_LW=LW_CONTROL(J)%L_N2O
          L_CH4_LW=LW_CONTROL(J)%L_CH4
          L_CFC11_LW=LW_CONTROL(J)%L_CFC11
          L_CFC12_LW=LW_CONTROL(J)%L_CFC12
          L_CFC113_LW=LW_CONTROL(J)%L_CFC113
          L_HCFC22_LW=LW_CONTROL(J)%L_HCFC22
          L_HFC125_LW=LW_CONTROL(J)%L_HFC125
          L_HFC134A_LW=LW_CONTROL(J)%L_HFC134A
          I_ST_WATER_LW=LW_CONTROL(J)%I_ST_WATER
          I_CNV_WATER_LW=LW_CONTROL(J)%I_CNV_WATER
          I_ST_ICE_LW=LW_CONTROL(J)%I_ST_ICE
          I_CNV_ICE_LW=LW_CONTROL(J)%I_CNV_ICE
          L_MICROPHYSICS_LW=LW_CONTROL(J)%L_MICROPHYSICS
          L_LOCAL_CNV_PARTITION_LW=LW_CONTROL(J)%L_LOCAL_CNV_PARTITION
          I_ANGULAR_INTEGRATION_LW=LW_CONTROL(J)%I_ANGULAR_INTEGRATION
          I_SPH_ALGORITHM_LW=LW_CONTROL(J)%I_SPH_ALGORITHM
          L_EULER_TRNF_LW=LW_CONTROL(J)%L_EULER_TRNF
          I_TRUNCATION_LW=LW_CONTROL(J)%I_TRUNCATION
          LS_GLOBAL_TRUNC_LW=LW_CONTROL(J)%LS_GLOBAL_TRUNC
          MS_MIN_LW=LW_CONTROL(J)%MS_MIN
          MS_MAX_LW=LW_CONTROL(J)%MS_MAX
          ACCURACY_ADAPTIVE_LW=LW_CONTROL(J)%ACCURACY_ADAPTIVE
          LS_BRDF_TRUNC_LW=LW_CONTROL(J)%LS_BRDF_TRUNC
          L_HENYEY_GREENSTEIN_PF_LW=LW_CONTROL(J)%L_HENYEY_GREENSTEIN_PF
          I_SPH_MODE_LW=LW_CONTROL(J)%I_SPH_MODE

!         Satellite Data:
          L_SUBSAMPLE=LW_CONTROL(J)%L_SUBSAMPLE
          L_GEOSTATIONARY=LW_CONTROL(J)%L_GEOSTATIONARY
          SAT_DESC=LW_CONTROL(J)%SAT_DESC
          SAT_HGT=LW_CONTROL(J)%SAT_HGT
          SAT_LON=LW_CONTROL(J)%SAT_LON
          SAT_LAT=LW_CONTROL(J)%SAT_LAT
          MAX_VIEW_LON=LW_CONTROL(J)%MAX_VIEW_LON
          MIN_VIEW_LON=LW_CONTROL(J)%MIN_VIEW_LON
          MAX_VIEW_LAT=LW_CONTROL(J)%MAX_VIEW_LAT
          MIN_VIEW_LAT=LW_CONTROL(J)%MIN_VIEW_LAT

!         Options for the longwave
          READ(5, R2LWCLNL)

!         Transfer values which we allow the user to change
!         from the namelist to the structure.

          LW_CONTROL(J)%SPECTRAL_FILE=SPECTRAL_FILE_LW
          LW_CONTROL(J)%FIRST_BAND=FIRST_BAND_LW
          LW_CONTROL(J)%LAST_BAND=LAST_BAND_LW
          LW_CONTROL(J)%I_2STREAM=I_2STREAM_LW
          LW_CONTROL(J)%L_IR_SOURCE_QUAD=L_IR_SOURCE_QUAD_LW
          LW_CONTROL(J)%I_GAS_OVERLAP=I_GAS_OVERLAP_LW
          LW_CONTROL(J)%I_CLOUD=I_CLOUD_LW
          LW_CONTROL(J)%I_CLOUD_REPRESENTATION=I_CLOUD_REPRESENTATION_LW
          LW_CONTROL(J)%I_SOLVER=I_SOLVER_LW
          LW_CONTROL(J)%L_N2O=L_N2O_LW
          LW_CONTROL(J)%L_CH4=L_CH4_LW
          LW_CONTROL(J)%L_CFC11=L_CFC11_LW
          LW_CONTROL(J)%L_CFC12=L_CFC12_LW
          LW_CONTROL(J)%L_CFC113=L_CFC113_LW
          LW_CONTROL(J)%L_HCFC22=L_HCFC22_LW
          LW_CONTROL(J)%L_HFC125=L_HFC125_LW
          LW_CONTROL(J)%L_HFC134A=L_HFC134A_LW
          LW_CONTROL(J)%I_ST_WATER=I_ST_WATER_LW
          LW_CONTROL(J)%I_CNV_WATER=I_CNV_WATER_LW
          LW_CONTROL(J)%I_ST_ICE=I_ST_ICE_LW
          LW_CONTROL(J)%I_CNV_ICE=I_CNV_ICE_LW
          LW_CONTROL(J)%L_MICROPHYSICS=L_MICROPHYSICS_LW
          LW_CONTROL(J)%L_LOCAL_CNV_PARTITION=L_LOCAL_CNV_PARTITION_LW
          LW_CONTROL(J)%I_ANGULAR_INTEGRATION=I_ANGULAR_INTEGRATION_LW
          LW_CONTROL(J)%I_SPH_ALGORITHM=I_SPH_ALGORITHM_LW
          LW_CONTROL(J)%L_EULER_TRNF=L_EULER_TRNF_LW
          LW_CONTROL(J)%I_TRUNCATION=I_TRUNCATION_LW
          LW_CONTROL(J)%LS_GLOBAL_TRUNC=LS_GLOBAL_TRUNC_LW
          LW_CONTROL(J)%MS_MIN=MS_MIN_LW
          LW_CONTROL(J)%MS_MAX=MS_MAX_LW
          LW_CONTROL(J)%ACCURACY_ADAPTIVE=ACCURACY_ADAPTIVE_LW
          LW_CONTROL(J)%LS_BRDF_TRUNC=LS_BRDF_TRUNC_LW
          LW_CONTROL(J)%L_HENYEY_GREENSTEIN_PF=L_HENYEY_GREENSTEIN_PF_LW
          LW_CONTROL(J)%I_SPH_MODE=I_SPH_MODE_LW

!         Satellite Data:
          LW_CONTROL(J)%L_SUBSAMPLE=L_SUBSAMPLE
          LW_CONTROL(J)%L_GEOSTATIONARY=L_GEOSTATIONARY
          LW_CONTROL(J)%SAT_DESC=SAT_DESC
          LW_CONTROL(J)%SAT_HGT=SAT_HGT
          LW_CONTROL(J)%SAT_LON=SAT_LON
          LW_CONTROL(J)%SAT_LAT=SAT_LAT
          LW_CONTROL(J)%MAX_VIEW_LON=MAX_VIEW_LON
          LW_CONTROL(J)%MIN_VIEW_LON=MIN_VIEW_LON
          LW_CONTROL(J)%MAX_VIEW_LAT=MAX_VIEW_LAT
          LW_CONTROL(J)%MIN_VIEW_LAT=MIN_VIEW_LAT
        ENDDO
      ENDIF
#endif

      READ (5, CLMCHFCG)
      IF ( L_CLMCHFCG ) THEN
         !  Convert rates from percent to multiplicative factors:
         DO J=1, NCLMFCGS
!          !  This is a null loop, as it should be, if CLIM_FCG_NYEARS=0
           DO JJ=1, CLIM_FCG_NYEARS(J)
             IF ( CLIM_FCG_RATES(JJ,J)  >   -100. )                     &
     &           CLIM_FCG_RATES(JJ,J) = 1. + 0.01 * CLIM_FCG_RATES(JJ,J)
           ENDDO
         ENDDO
       ELSE
!        ! If the namelist is not to be read, set number of designated
!        !   years to zero for all possible forcings, as this may be
!        !   used to test if this system is being used for each forcing.
         DO J=1, NCLMFCGS
           CLIM_FCG_NYEARS(J) = 0
         ENDDO
      ENDIF !  L_CLMCHFCG
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      DO J=1, C2C_size
        IF ( .NOT. ( C2C_O2(J) .OR. C2C_O3(J) .OR. C2C_CO2(J) .OR.      &
     &     C2C_N2O(J) .OR. C2C_CH4(J) .OR. C2C_CFC11(J) .OR.            &
     &     C2C_CFC12(J) .OR. C2C_C113(J) .OR. C2C_HCFC22(J) .OR.        &
     &     C2C_HFC125(J) .OR. C2C_HFC134(J) .OR. C2C_AEROSOL(J) .OR.    &
     &     C2C_SULPC_D(J) .OR. C2C_SEAS_D(J) .OR. C2C_SOOT_D(J) .OR.    &
     &     C2C_BMB_D(J) .OR. C2C_OCFF_D(J) .OR. C2C_SUN(J) .OR.         &
     &     C2C_VOL(J) .OR. C2C_LAND_S(J) .OR. C2C_ALL(J) .OR.           &
     &     C2C_WMG(J) ) ) THEN
!          ErrorStatus = 614
          Write(6,*) 'Warning: Forcing switched on but nothing,         &
     &                set to change in diagnostic call!'
!          CALL Ereport(RoutineName,ErrorStatus,Cmessage)
        ENDIF
        IF ( C2C_ALL(J) ) THEN
! C_O3(J) = .TRUE.  ! Commented out till we have the ancillary
! Best made conditional on the alternate ozone file being specified
! if we don't hardwire the latter's name (in which case it'd always
! be safe to read from that file, & add nothing much to the cost of
! a 2nd call which is going to be made anyway).
          C2C_WMG(J) = .TRUE.
          C2C_AEROSOL(J) = .TRUE.
          C2C_LAND_S(J) = .TRUE.
! This is used for surface albedo forcing
        ENDIF
        IF ( C2C_WMG(J) ) THEN
          C2C_O2(J) = .TRUE.
          C2C_CO2(J) = .TRUE.
          C2C_N2O(J) = .TRUE.
          C2C_CH4(J) = .TRUE.
          C2C_CFC11(J) = .TRUE.
          C2C_CFC12(J) = .TRUE.
          C2C_C113(J) = .TRUE.
          C2C_HCFC22(J) = .TRUE.
          C2C_HFC125(J) = .TRUE.
          C2C_HFC134(J) = .TRUE.
        ENDIF
        IF ( C2C_AEROSOL(J) ) THEN
          C2C_SULPC_D(J) = .TRUE.
          C2C_SEAS_D(J) = .TRUE.
          C2C_SOOT_D(J) = .TRUE.
          C2C_BMB_D(J) = .TRUE.
          C2C_OCFF_D(J) = .TRUE.
        ENDIF
      ENDDO
#endif
!

! Check error condition
      IF(ErrorStatus >  0) THEN
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

      RETURN
      END SUBROUTINE Readlsta

#endif
