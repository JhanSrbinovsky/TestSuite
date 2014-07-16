#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! SCMSHELL is the main calling program for the Single Column Model.
! It sets up the vertical level information read in from the UMUI and
! passes the information down to S_MAIN which then performs the run.
!
! Program scmshell
!=====================================================================
!                     SCM
!           Single Column Unified Model
!                  Master Deck
!
! This is a new deck created at 5.3 by Z. Gardner
!
!=====================================================================
      PROGRAM SCM_SHELL
 
      Use cv_run_mod                         ! Access to all variables

! module for RUN_LAND namelist
      USE LAND_SURF_MOD, ONLY :                                         &
     & FRAC_SNOW_SUBL_MELT                                              &
     &,MASKD                                                            &
     &,SOILHC_METHOD                                                    &
     &,ALL_TILES                                                        &
     &,L_VG_SOIL                                                        &
     &,RUN_LAND


      USE BL_OPTION_MOD, ONLY : WeightLouisToLong


      Implicit none

#include "nstypes.h"
#include "chsunits.h"
#include "s_dims.h"
#include "cntlatm.h"
#include "cntlall.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "swopt3a.h"
#include "swcopt3a.h"
#include "lwopt3a.h"
#include "lwcopt3a.h"
#include "blopt8a.h"
#include "c_mdi.h"
#include "ctlnl3a.h"

#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
#include "satopt.h"
#endif
! for N_INTERNAL_MODEL in typsts.h
#include "csubmodl.h"

#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      Use SW_CONTROL_STRUCT
      Use LW_CONTROL_STRUCT
#endif
!     Local Variables

      Integer row_length        ! Number of points the model is
      Parameter(row_length = 1) !   running on (1 now) = row_length
                                !                      * rows
      Integer rows
      Parameter(rows       = 1)

      Integer ntrop             ! Max number of levels in the
      Parameter(ntrop      = 20)!    troposphere STATS forcing

      Integer sec_day
      Parameter(sec_day    = 86400)

      Integer nsprog            ! no. of single level prognostics
      Parameter(nsprog     = 10)

      Integer ntab              ! Dimension of array used in random
      Parameter(ntab       = 32)!  generator (Do not change this
                                !  value as it is hard coded into
                                !  the S_RANDOM deck)
      Integer                                                           &
     & CO2_DIM_LEN                                                      &
                                ! Length of a CO2 field row.
     &,CO2_DIM_ROW              ! Number of CO2 field rows.

      Parameter (co2_dim_len = 1, co2_dim_row = 1)

      Real                                                              &
     &  dummy

      Integer                                                           &
     &  ISTATUS                                                         &
     &, ICODE                                                           &
     &, first_blank                                                     &
     &, level                                                           &
     &, j

      Logical :: l_ts_log ! Option for timestep information

      Character*200                                                     &
     &  sdum0,sdum1                                                     &
                          ! Dummy strings
     &, dname                                                           &
                          ! Directory where the sizes file is kept
     &, filename                                                        &
                          ! Sizes filename to read in basic model
                          ! dimensions
     &, vert_lev
                          ! Vertical level file


      CHARACTER*(256) CMESSAGE  ! Error message if ErrorStatus > 0
      CHARACTER (Len=*), Parameter :: RoutineName = 'scm_shell'

      !-----------------------------------------------------------------
      ! Variables read in which refer to the SCM forcing namelist 
      !-----------------------------------------------------------------
      Integer :: land_points      = 0    ! Default number of land points
      Integer :: nfor             = -999 ! Number terms for observational
                                         ! forcing, requires setting in
                                         ! namelist
      Integer :: model_levels_nml = -999 ! Number of model levels
                                         ! specified in supplied
                                         ! namelist. Must be set in
                                         ! namelist.

      Namelist/NLCFILES/ vert_lev
      Namelist/CNTLSCM/ nfor, model_levels_nml, l_ts_log, land_points
#include "cv_run_nml.h"

!=====================================================================
!     First read in directory where UMUI files
!=====================================================================
      Open(10, File='dir_name', Iostat=Istatus)

      Read(10,*) dname

      Close(10)

      first_blank=len_trim(dname)+1
      dname = trim(dname)

!=====================================================================
!     Read SCM runtime namelist CNTLSCM
!=====================================================================

      Filename = dname(1:first_blank-1)//'/namelist.scm'

      Open(10, File=Filename, Iostat=IstatuS, Status='old')

      If (Istatus /= 0) Then
        Icode=500
        Write(Cmessage,*)  " Error opening " //TRIM(ADJUSTL(Filename))//      &
                           " file on unit 10"

! DEPENDS ON: ereport
        Call Ereport(RoutineName,Icode,Cmessage)

      End if

      Read(10,CNTLSCM)

      Close(10)

      If (model_levels_nml == -999) Then
        Icode=501
        Write(6,*) " "
        Write(6,*) "============================================="
        Write(6,*) " Number of model levels (model_levels_nml) in"
        Write(6,*) " SCM namelist (&CNTLSCM) has not been set.   "
        Write(6,*) "============================================="
        Write(6,*) " "
        Write(6,*) "Run ABORTED"
        Write(6,*) " "

        Write(Cmessage,*) ' Variable MODEL_LEVELS_NML has not been set'

! DEPENDS ON: ereport
        Call Ereport(RoutineName,Icode,Cmessage)
      End if

      If (land_points > row_length*rows) Then
        Icode=502
        Write(6,*) " "
        Write(6,*) "============================================="
        Write(6,*) " Specified number of land points greater than"
        Write(6,*) " row_length*rows.                            "
        Write(6,*) "============================================="
        Write(6,*) " "
        Write(6,*) "Run ABORTED"
        Write(6,*) " "

        Write(Cmessage,*) ' Too many land points specified'

! DEPENDS ON: ereport
        Call Ereport(RoutineName,Icode,Cmessage)
      End if
      
!=====================================================================
! Read in model dimensions from SIZES file
!=====================================================================

      Filename = dname(1:first_blank-1)//'/SIZES'

      Open(10, File=Filename, Iostat=Istatus, Status='old')
      If(Istatus  /=  0) Then
        Icode=500
        Write(Cmessage,*) " Error opening " //TRIM(ADJUSTL(Filename))//       &
                          " on unit 10"

! DEPENDS ON: ereport
        Call Ereport(RoutineName,Icode,Cmessage)
      End If

      Read(10,NLSIZES)

      Close(10)

      !---------------------------------------------------------------
      ! Check that number of Wet and Cloud levels do not exceed the
      ! number of model levels
      !---------------------------------------------------------------
      If (wet_levels  >   model_levels) Then

        write (*,'(4(TR1,A52/))')                                       &
     &        '===================================================='    &
     &,       '| Warning : Wet_levels > Model_levels              |'    &
     &,       '| Setting : Wet_levels = Model_levels              |'    &
     &,       '===================================================='

        wet_levels = model_levels

      End If

      If (cloud_levels  >   wet_levels) Then
        write (*,'(4(TR1,A52/))')                                       &
     &        '===================================================='    &
     &,       '| Warning : Cloud_levels > Wet_levels              |'    &
     &,       '| Setting : Cloud_levels = Wet_levels              |'    &
     &,       '===================================================='

        cloud_levels = wet_levels
      End If

      !---------------------------------------------------------------
      ! Check to see if model_levels for job in SIZES namelist matches
      ! with the number of model levels (model_levels_nml) specified
      ! by the input namelist.
      !---------------------------------------------------------------

      If (model_levels_nml /= model_levels) Then
        Icode=502
        Write(sdum0,*) model_levels_nml
        Write(sdum1,*) model_levels

        Write(6,*) " "
        Write(6,*) "============================================="
        Write(6,*) " Number of model levels in SCM namelist(" //      &
     &              TRIM(ADJUSTL(sdum0))                      // ")"
        Write(6,*) " is different from that expected by SCM "        
        Write(6,*) " executable(" // TRIM(ADJUSTL(sdum1))     // ")"
        Write(6,*) "============================================="
        Write(6,*) " "
        Write(6,*) "Run ABORTED"
        Write(6,*) " "

        Write(Cmessage,*)                                                       &
                   ' Model level mismatch between namelist and executable'

! DEPENDS ON: ereport
        Call Ereport(RoutineName,Icode,Cmessage)
      End if
   

!=====================================================================
!     Now read in the name of the file containing the vertical level
!     information
!=====================================================================
      Filename = dname(1:first_blank-1)//'/INITHIS'

      Open(10, File=Filename, Iostat=IstatuS, Status='old')
      If(Istatus  /=  0) Then
        Icode=500
        Write(Cmessage,*) " Error opening " //TRIM(ADJUSTL(Filename))//       &
                          " on unit 10"
! DEPENDS ON: ereport
        Call Ereport(RoutineName,Icode,Cmessage)
      End If

      Read(10,NLCFILES)

      Close(10)





!=====================================================================
!     Read in the namelists from the CNTLALL file
!=====================================================================
      Filename = dname(1:first_blank-1)//'/CNTLALL'

      Open(10, File=Filename, Iostat=Istatus, Status='old')
      If(Istatus  /=  0) Then
        Icode=500
        Write(Cmessage,*) " Error opening " //TRIM(ADJUSTL(Filename))//       &
                          " on unit 10"

! DEPENDS ON: ereport
        Call Ereport(RoutineName,Icode,Cmessage)
      End If

      Read(10,NLSTCALL)

      Close(10)

!=====================================================================
!     Set defaults for RUN_BL namelist
!=====================================================================
      DO level=1,max_number_alpha_cds
         alpha_Cd (level) = RMDI
      ENDDO
      L_SBLeq             = .FALSE.
      L_SBLco             = .TRUE.
      Muw_SBL             = 1.5
      Mwt_SBL             = 1.0
      ISHEAR_BL           = ON
      NG_STRESS           = OFF
      DECFIX              = ON
      STOPWE_SBL          = ON
      TRWEIGHTS1          = ON
      FLUX_GRAD           = Locketal2000
      SBL_OP              = Long_tails
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
      LOCAL_FA            = OFF
      PRANDTL             = Constant_SBL
      Keep_Ri_FA          = OFF
      NL_BL_LEVELS        = OFF
      WeightLouisToLong_in = 0.0
      l_use_bl_diag_term  = .false.
      l_emis_land_gen     = .false.
      emis_land_gen       = 1.0
      DO J = 1, 20
        BL_OPTIONS(J)     = OFF
      ENDDO

!=====================================================================
!     Set defaults for RUN_Convection namelist
!=====================================================================
#include "cv_run_nml_data.h"

!=====================================================================
!     Read in the namelists from the CNTLATM file
!=====================================================================
      Filename = dname(1:first_blank-1)//'/CNTLATM'

      Open(10, File=Filename, Iostat=Istatus, Status='old')
      If(Istatus  /=  0) Then
        Icode=500
        Write(Cmessage,*) " Error opening " //TRIM(ADJUSTL(Filename))//       &
                          " on unit 10"

! DEPENDS ON: ereport
        Call Ereport(RoutineName,Icode,Cmessage)
      End If

!     L_MR_PHYSICS2 default
      L_MR_PHYSICS2 = .FALSE.

!     Default value of CAN_RAD_MOD which is defined in RUN_BLVEG
      CAN_RAD_MOD = 1     

      Read(10,NLSTCATM)
      Read(10,RUN_BL)

      If (Buddy_sea /= OFF) Then 
         
        Write (*,'(8(TR1,A52/))')                                       & 
              '===================================================='    & 
      ,       '| Invalid Boundary Layer options for SCM           |'    & 
      ,       '| Altering the following RUN_BL namelist variable  |'    & 
      ,       '|                                                  |'    & 
      ,       '|  Buddy_sea                                       |'    & 
      ,       '|                                                  |'    & 
      ,       '====================================================' 

        Buddy_sea = OFF 
 
      End If ! Test for invalid RUN_BL options 

      Read(10,RUN_BLICE)
      Read(10,RUN_BLVEG)
!     ! Store BL options in integer array for ease of passing around
!     ! Changes here must be duplicated in READLSA2
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


      READ (10,RUN_PFT)           ! Surface parameters
      Read(10,RUN_LAND)           ! land-surface parameters
      Read(10,RUN_Precip)
      Read(10,RUN_Cloud)
      Read(10,RUN_Convection)
      Read(10,RUN_Radiation)
      Read(10,RUN_GWD)
      Read(10,RUN_Aerosol)
      Read(10,RUN_UKCA)
      Read(10,RUN_Dyn)
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
!
! Read in Namelists R2SWCLNL and R2LWCLNL and transfer data to
! the structure SW_CONTROL and LW_CONTROL. This part of the
! code (under 3C and 3Z) uses modules to pass arguments around.
!
      Read(10, R2SWNCAL)

      Do j=1,n_swcall
! DEPENDS ON: sw_control_default
        Call sw_control_default(sw_control(j))
!
!         Initialize satellite data to safe values
!
        l_subsample    =.FALSE.
        l_geostationary=.FALSE.
        sat_desc="                    " //                              &
     &           "                    " //                              &
     &           "                    " //                              &
     &           "                    "
        sat_hgt=0.0
        sat_lon=0.0
        sat_lat=0.0
        max_view_lon=0.0
        min_view_lon=0.0
        max_view_lat=0.0
        min_view_lat=0.0
!
!         Options for the shortwave
!
        Read(5, R2SWCLNL)
!
!         Transfer values which we allow the user to change
!         from the namelist to the structure.
!
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
!
!         Satellite Data:
!
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

      Enddo       !  Loop over SW Calls

      Read(10, R2LWNCAL)
      Do j=1, n_lwcall
! DEPENDS ON: lw_control_default
        Call lw_control_default(lw_control(j))
!
!         Initialize satellite data to safe values
!
        l_subsample    =.FALSE.
        l_geostationary=.FALSE.
        sat_desc="                    " //                              &
     &           "                    " //                              &
     &           "                    " //                              &
     &           "                    "
        sat_hgt=0.0
        sat_lon=0.0
        sat_lat=0.0
        max_view_lon=0.0
        min_view_lon=0.0
        max_view_lat=0.0
        min_view_lat=0.0
!
!         Options for the longwave
!
        Read(5, R2LWCLNL)
!
!         Transfer values which we allow the user to change
!         from the namelist to the structure.
!
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
!
!         Satellite Data:
!
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
!
      Enddo       !  Loop over LW calls
#else
      Read(10,R2SWCLNL)
      Read(10,R2LWCLNL)
#endif

      Close(10)

      n_internal_model = 1

!------------------------------------------------------- 
! Check that if L_CCRAD is selected certain other
! switches are not set so that they conflict.
!------------------------------------------------------- 
! Error capture 
!------------------------------------------------------- 
                 
      If (L_ccrad) Then 
 
       If (.NOT. l_3d_cca) Then 
          Icode       = 100 
          CMessage    = '**ERROR**: CCRad is not yet available without'// &
                                  ' the anvil scheme (L_3D_CCA = .True.)'
! DEPENDS ON: ereport 
          CALL ereport(RoutineName, Icode, CMessage) 
        End If

        If (l_convcld_hadgem1) Then 
          Icode       = 100 
          CMessage    = '**ERROR**: L_CCRad and l_convcld_hadgem1'//      &
                                  ' should not be both set to true.'
! DEPENDS ON: ereport 
          CALL ereport(RoutineName, Icode, CMessage) 
        End If

        If (l_fix_udfactor) Then 
          Icode       = 100 
          CMessage    = '**ERROR**: L_CCRad and l_fix_udfactor'//         &
                                  ' should not be both set to true.'
! DEPENDS ON: ereport 
          CALL ereport(RoutineName, Icode, CMessage) 
        End If

        If (l_pc2_diag_sh) Then 
          Icode       = 100 
          CMessage    = '**ERROR**: L_CCRad and l_pc2_diag_sh'//          &
                                  ' should not be both set to true.'
! DEPENDS ON: ereport 
          CALL ereport(RoutineName, Icode, CMessage) 
        End If

      End If      ! l_ccrad 

!------------------------------------------------------- 
! End error capture 
!-------------------------------------------------------

!=====================================================================
!     Call the main part of the model
!=====================================================================
! DEPENDS ON: scm_main
      Call SCM_MAIN(vert_lev, nfor, l_ts_log, ntrop, sec_day                  &
         , land_points, nsprog, ntab, co2_dim_len, co2_dim_row                &
         , cloud_levels, model_levels, wet_levels, tr_levels                  &
         , tr_vars, tr_ukca, st_levels, sm_levels, bl_levels                  &
         , ozone_levels, ntiles                                               &
#include "swcarg3a.h"
#include "lwcarg3a.h"
         , atmos_sr, nsectp                                                   &
         )

      END PROGRAM SCM_SHELL

#endif
