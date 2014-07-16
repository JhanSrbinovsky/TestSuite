!include file: s_maind.h
! Description:
!   Sets up data for the SCM.
!
! Current Code Owner: R Wong
!

#if defined(ATMOS)

! &INDATA

      Data                                                              &
     &  year_init, month_init, day_init, hour_init, min_init,sec_init,  &
     &  tapeyear_init, tapemonth_init, tapeday_init, tapehour_init,     &
     &  tapemin_init, tapesec_init /                                    &
     &  1998,      1,          1,        0,          0,       0,        &
     &  1998 ,         1,              1,             0,                &
     &  0,            0 /

      Data salt_dim1, salt_dim2, salt_dim3 /3*1/
      Data gather  /.false./

! &INOBSFOR

      Data ichgf /1/
! Initialise new variables for forcing
      Data L_windrlx, L_vertadv /.False., .False./
      Data tau_rlx /3600.0/
      Data tau_T_eurocs /3600.0/
      Data tau_q_eurocs /3600.0/
      Data tau_uv_eurocs /3600.0/

! &LOGIC

      Data ancyc, land_ice_mask, land_sea_mask, soil_mask,              &
     &            altdat, local_time/ .true., .false.,4*.true./
      Data obs, prindump_obs, stats, noforce/ 4*.false./
      Data geoforce, geoinit / 2*.false./
      Data test, prinstat, radcloud_fixed/ 3*.false./
      Data prindump_day, prindump_days, prindump_step/ 3*.false./
      Data grafdump_day, grafdump_days, grafdump_step/ 3*.false./
      Data tapein, tapeout/ 2*.false./
      Data cumulus /.false./
      Data l_do_t, l_do_inc_vels, l_do_inc_q /3*.false./
      Data L_flux_bc,L_spec_z0 /2*.false./
      Data relaxT_eurocs /.false./
      Data relaxq_eurocs /.false./
      Data relaxuv_eurocs /.false./


! &RUNDATA

      Data cort, cord, corvn, corw / 0.9, 0.9, 0.5, 0.5 /
      Data dump_step /0/
      Data dump_days / 4*1/
      Data resdump_days, change_clim / 1, 10/
      Data exname_in, exname_out/ 'XXXXXXXX', 'XXXXXXXX'/
      Data ndayin, timestep, ntrad, ntrad1/ 1, 1800.0, 6, 1/
      Data nminin, nsecin/2*0/
      Data runno_in, runno_out/ 0, 999/
      Data tstar_land, tstar_sea, tstar_sice /3*-999.0/
      Data nbdsc, ntdsc /0,0/
      Data co2start /  4.9E-4 /
      Data co2end   /  4.9E-4 /
      Data co2rate  /  0.0    /
      Data min_trop_level, max_trop_level /0, 0/


! &PHYSWITCH

      Data conv_mode /0/

!---------------------------------------------------------------------
!     Initialise the array giving the unit nos. for output of the
!     initial data
!---------------------------------------------------------------------

      Data nout / 19*0 /

!---------------------------------------------------------------------
#endif
