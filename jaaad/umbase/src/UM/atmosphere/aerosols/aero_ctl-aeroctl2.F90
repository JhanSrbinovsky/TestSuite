#if defined(A17_2A) || defined(A17_2B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      SUBROUTINE AERO_CTL(                                              &
!
! Parallel variables
     &  halo_i, halo_j, off_x, off_y, global_row_length, global_rows    &
     &, proc_row_group, proc_col_group, at_extremity, n_proc, n_procx   &
     &, n_procy, neighbour, g_rows, g_row_length, g_datastart, me       &
! model dimensions
     &, row_length, rows, n_rows, land_points                           &
     &, model_levels, wet_model_levels, bl_levels, n_cca_levels         &
     &, theta_field_size                                                &
     &, salt_dim1, salt_dim2, salt_dim3                                 &
     &, aero_dim1, aero_dim2, aero_dim3                                 &
! Model switches
     &, model_domain, L_CAL360, L_SEC_VAR, L_EqT, Ltimer                &
! Model parameters
     &, Ntot_land, Ntot_sea                                             &
! Physical constants
     &, lc, lf, cp, two_Omega, p_zero, kappa                            &
     &, R, g, Lapse_Rate, earth_radius, Pi                              &
! Co-ordinate information
     &, r_rho_levels, r_theta_levels                                    &
     &, eta_theta_levels, eta_rho_levels                                &
     &, delta_lambda, delta_phi                                         &
     &, lat_rot_NP, long_rot_NP                                         &
! Time stepping information
     &, timestep                                                        &
     &, val_year, val_day_number, val_hour, val_minute                  &
     &, val_second, timestep_number                                     &
     &, PREVIOUS_TIME                                                   &
     &, CALL_CHEM_FREQ                                                  &
! Trig arrays
     &, sin_theta_longitude, cos_theta_longitude                        &
     &, FV_cos_theta_latitude                                           &
! Grid-dependent arrays
     &, f3_at_u, true_longitude, true_latitude                          &
!
! Data fields IN
     &, u, v, Tstar, Tstar_sea                                          &
     &, theta, q, qcl, qcf                                              &
     &, rho, land_mask, fland_ctile                                     &
     &, p_theta_levels, exner_rho_levels, exner_theta_levels            &
     &, ice_fract, snow_depth                                           &
     &, CLOUDF_halos                                                    &
     &, OH_CONC, H2O2_LMT, HO2_CONC, O3                                 &
     &, SO2_SURFEM, SO2_HILEM, SO2_NATEM                                &
     &, DMS_em_ancil, DMS_conc, NH3_EM                                  &
     &, DMS_Ointer                                                      &
     &, SOOT_SUREM, SOOT_HILEM, BMASS_SUREM, BMASS_HILEM, OCFF_SUREM    &
     &, OCFF_HILEM, SO2_HIGH_LEVEL, SOOT_HIGH_LEVEL, BMASS_HIGH_LEVEL_1 &
     &, BMASS_HIGH_LEVEL_2, OCFF_HIGH_LEVEL, land_index                 &
! Logicals IN
     &, L_SULPC_SO2, L_SULPC_DMS, L_SULPC_OZONE                         &
     &, L_SULPC_SO2_O3_NONBUFFERED, L_SULPC_NH3                         &
     &, L_sulphate_CCN, L_seasalt_CCN, L_SOOT, L_soot_CCN               &
     &, L_biomass, L_biomass_CCN, L_ocff, L_ocff_CCN                    &
     &, L_SO2_SURFEM, L_SO2_HILEM, L_SO2_NATEM, L_DMS_EM                &
     &, L_DMS_em_inter, L_DMS_Liss_Merlivat                             &
     &, L_DMS_Wanninkhof, L_DMS_Nightingale                             &
     &, L_DMS_Ointer                                                    &
     &, L_NH3_EM, L_ctile                                               &
     &, L_SOOT_SUREM, L_SOOT_HILEM, L_BMASS_SUREM, L_BMASS_HILEM        &
     &, L_OCFF_SUREM, L_OCFF_HILEM, L_USE_BIOGENIC                      &
     &, L_USE_SEASALT_DIRECT, L_USE_SEASALT_INDIRECT                    &
     &, L_USE_SEASALT_AUTOCONV, L_DUST                                  &
!
! Data fields IN/OUT
     &, SO2, DMS                                                        &
     &, SO4_AIT, SO4_ACC, SO4_DIS                                       &
     &, H2O2_MXR, NH3                                                   &
     &, SOOT_NEW, SOOT_AGD, SOOT_CLD                                    &
     &, BMASS_NEW, BMASS_AGD, BMASS_CLD                                 &
     &, OCFF_NEW, OCFF_AGD, OCFF_CLD                                    &
     &, BIOGENIC                                                        &
!
! Data fields IN
     &, DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5           &
!
! Data fields OUT
! Diagnostic info
     &,                                                                 &
#include "argsts.h"
     &  STASHwork17                                                     &
     &, Error_code                                                      &
     & )
!
!---------------------------------------------------------------------
! Purpose: Interface for Aerosol Modelling, to include Sulphur Cycle
!          and SOOT modelling
!
! Level 2 control routine
!
! Current owners of code:                C E Johnson
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
! System components covered:
!
! System task:
!
! Documentation: Not yet available
!
!-----------------------------------------------------------------
!
      IMPLICIT NONE
!
! Arguments with intent IN:
!
! Parallel setup variables
      Integer                                                           &
     &  halo_i                                                          &
                   ! Size of halo in i direction.
     &, halo_j                                                          &
                   ! Size of halo in j direction.
     &, off_x                                                           &
                   ! Size of small halo in i
     &, off_y                                                           &
                   ! Size of small halo in j.
     &, global_row_length                                               &
                           ! number of points on a row
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, proc_col_group                                                  &
                       ! Group id for processors on the same col
     &, global_rows                                                     &
                           ! NUMBER OF global rows
     &, n_proc                                                          &
                   ! Total number of processors
     &, n_procx                                                         &
                   ! Number of processors in longitude
     &, n_procy                                                         &
                   ! Number of processors in latitude
     &, neighbour(4)                                                    &
                             ! Array with the Ids of the four neighbour
                             ! in the horizontal plane
     &, g_rows (0:n_proc-1)                                             &
     &, g_row_length (0:n_proc-1)                                       &
     &, g_datastart (3,0:n_proc-1)                                      &
     &, me         ! My processor number
!
      Logical                                                           &
     &  at_extremity(4) ! Indicates if this processor is at north, sout
                        ! east or west of the processor grid
!
#include "parparm.h"
!
! Model dimensions
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, n_rows                                                          &
     &, land_points                                                     &
     &, model_levels                                                    &
     &, wet_model_levels                                                &
     &, bl_levels                                                       &
     &, n_cca_levels                                                    &
                        ! No. conv cloud levels (1 if 2D, nlevs if 3D)
     &, theta_field_size                                                &
     &,salt_dim1                                                        &
                                           !dimensns of seasalt array
     &,salt_dim2                                                        &
     &,salt_dim3                                                        &
     &,aero_dim1                                                        &
                                           !row_length or 1
     &,aero_dim2                                                        &
                                           !rows or 1
     &,aero_dim3                           !model_levels or 1
! Model switches
      Integer                                                           &
     &  model_domain
      Logical                                                           &
     &  L_CAL360                                                        &
                        ! T if using 360 day calendar
     &, L_SEC_VAR                                                       & 
                        ! if T include secular varn of earth's orbit,
     &, L_EqT                                                           &
                        ! if T include equn of time         in SOLPOS
     &, Ltimer                                                           
                        ! if T then output some timing information

! Physical constants
      Real                                                              &
     &  lc, lf, cp                                                      &
     &, two_Omega                                                       & 
                        ! twice Earth's rotation rate
     &, p_zero                                                          &
     &, kappa                                                           &
     &, R, g, Lapse_Rate, earth_radius, Pi
! Co-ordinate arrays
      Real                                                              &
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:model_levels)           &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, model_levels)              &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)                                    &
     &, delta_lambda                                                    &
     &, delta_phi                                                       &
     &, lat_rot_NP                                                      &
     &, long_rot_NP

! Trig arrays
      Real                                                              &
     &  cos_theta_longitude (row_length, rows)                          &
     &, sin_theta_longitude (row_length, rows)
! Grid-dependent arrays
      Real                                                              &
     &  f3_at_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y)          &
     &, true_longitude(row_length, rows)                                &
     &, true_latitude(row_length, rows)                                 &

     &, FV_cos_theta_latitude (1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y)
! Time stepping information
      Real                                                              &
     &  timestep                       !atmosphere model timetsep
      Integer                                                           &
     &  val_year                                                        &
     &, val_day_number                                                  &
     &, val_hour                                                        &
     &, val_minute                                                      &
     &, val_second                                                      &
     &, timestep_number                                                 &
     &, PREVIOUS_TIME(7)                                                &
     &, CALL_CHEM_FREQ                  !frequency of calling chemistry
!                                         per atmos phys timestep
!
#include "csubmodl.h"
#include "typsts.h"
!
! Diagnostics info
      Real                                                              &
     &  STASHwork17(*)  ! STASH workspace for section 17 (Aero_Ctl)
!
      Real                                                              &
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &                                                 model_levels)    &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &                                                 model_levels)    &
     &, theta(1-off_x:row_length+off_x,1-off_y:rows+off_y,              &
     &                                                 model_levels)    &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &                                             wet_model_levels)    &
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &                                             wet_model_levels)    &
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &                                             wet_model_levels)    &
     &, rho(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &                                                 model_levels)    &
     &, p_theta_levels(1-off_x:row_length+off_x,                        &
     &                 1-off_y:rows+off_y, model_levels)                &
     &, exner_rho_levels(1-off_x:row_length+off_x,                      &
     &                           1-off_y:rows+off_y, model_levels+1)    &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                             1-off_y:rows+off_y, model_levels)    &
     &, ice_fract(row_length, rows)                                     &
     &, snow_depth(row_length, rows)                                    &
     &, land_fract(row_length, rows)                                    &
     &, Tstar(row_length, rows)                                         &
     &, Tstar_sea(row_length, rows)                                     &
     &, fland_ctile(land_points)                                        &
     &, CLOUDF_halos(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,  &
     &                                        wet_model_levels)         &
!
     &,OH_CONC(row_length,rows,model_levels)                            &
     &,HO2_CONC(row_length,rows,model_levels)                           &
     &,H2O2_LMT(row_length,rows,model_levels)                           &
     &,O3(row_length,rows,model_levels)                                 &
     &,SO2_SURFEM(row_length,rows)                                      &
                                          !SO2 emiss at surface
     &,SO2_HILEM(row_length,rows)                                       &
                                          !SO2 emiss at chimney level
     &,SO2_NATEM(row_length,rows,model_levels)                          &
                                                !Volcanic SO2 emiss
     &,DMS_em_ancil(row_length,rows)                                    &
                                          !Ancillary DMS emiss (surf)
     &,DMS_conc(row_length,rows)                                        &
                                          !Seawater DMS conc (n mol l-1)
     &,NH3_EM(row_length,rows)                                          &
                                         !NH3 emiss (surface)
     &,SOOT_SUREM(row_length,rows)                                      &
                                        !SOOT emiss at surface
     &,SOOT_HILEM(row_length,rows)                                      &
                                        !SOOT emiss at chimney lev
     &,BMASS_SUREM(row_length,rows)                                     &
                                        !Biomass surface emissions
     &,BMASS_HILEM(row_length,rows)                                     &
                                        !Biomass high level emissions
     &,OCFF_SUREM(row_length,rows)                                      &
                                        !OCFF emiss at surface
     &,OCFF_HILEM(row_length,rows)
                                        !OCFF emiss at chimney lev

!
      Real                                                              &
     & Ntot_land                                                        &
                           ! Number of droplets over land / m-3
     &,Ntot_sea            ! Number of droplets over sea / m-3

      Integer                                                           &
     & SO2_HIGH_LEVEL                                                   &
                                         !Model level for chimney emiss
     &,SOOT_HIGH_LEVEL                                                  &
                                         !Model level for chimney emiss
     &,BMASS_HIGH_LEVEL_1                                               &
                                         !Lowest and highest model
     &,BMASS_HIGH_LEVEL_2                                               &
                                         !levels with biomass emiss
     &,OCFF_HIGH_LEVEL                                                  &
                                         !Model level for chimney emiss
     &,land_index(land_points)
!
      LOGICAL                                                           &
     & L_SULPC_SO2                                                      &
                                         !T if Sulphur Cycle required
     &,L_SULPC_DMS                                                      &
                                         !T if DMS chemistry required
     &,L_SULPC_OZONE                                                    &
                                         !T if O3 oxidn required
     &,L_SULPC_SO2_O3_NONBUFFERED                                       &
                                         !T if SO2+O3 reaction is NOT
                                         !  to be buffered by NH3.
     &,L_SULPC_NH3                                                      &
                                         !T if NH3 buffering required
                                         !(always T if L_SULPC_OZONE
                                         ! is T)
     &,L_sulphate_CCN                                                   &
                                         !T if sulphate used for CCN
     &,L_seasalt_CCN                                                    &
                                         !T if sea-salt used for CCN
     &,L_soot_CCN                                                       &
                                         !T if soot used for CCN
     &,L_biomass_CCN                                                    &
                                         !T if biomass used for CCN
     &,L_ocff_CCN                                                       &
                                         !T if OCFF used for CCN
     &,L_ctile                                                          &
                                         !T if coastal tiling is on
     &,L_SOOT                                                           &
                                         !T if SOOT modelling required
     &,L_BIOMASS                                                        &
                                         !T if biomass modelling reqd
     &,L_OCFF                                                           &
                                         !T if OCFF modelling required
     &,L_SO2_SURFEM                                                     &
                                         !T if surface SO2 ems present
     &,L_SO2_HILEM                                                      &
                                         !T if high lev SO2 em present
     &,L_SO2_NATEM                                                      &
                                         !T if volcanic SO2 ems present
     &,L_DMS_EM                                                         &
                                         !T if DMS emiss present (surf)
     &,L_DMS_em_inter                                                   &
                                         !T if interactive DMS emiss
     &,L_DMS_Ointer                                                     &
                                         !T if ocean DMS emiss model
     &,L_DMS_Liss_Merlivat                                              &
                                         !Switches to determine which
     &,L_DMS_Wanninkhof                                                 &
                                         !  scheme to use for inter-
     &,L_DMS_Nightingale                                                &
                                         !  active DMS emissions
     &,L_NH3_EM                                                         &
                                         !T if NH3 emiss present (surf)
     &,L_SOOT_SUREM                                                     &
                                         !T if surface SOOT ems present
     &,L_SOOT_HILEM                                                     &
                                         !T if high lev SOOT ems presnt
     &,L_BMASS_SUREM                                                    &
                                         !T if sfc biomass ems present
     &,L_BMASS_HILEM                                                    &
                                         !T if hi lev bmass ems present
     &,L_OCFF_SUREM                                                     &
                                         !T if surface OCFF ems present
     &,L_OCFF_HILEM                                                     &
                                         !T if high lev OCFF ems presnt
     &,L_USE_BIOGENIC                                                   &
                                         !T if using biogenics for CCN
     &,LAND_MASK(row_length,rows)                                       &    
                                         !T IF LAND, F IF SEA
     &,L_DUST                                                           &
                                         !T if mineral dust used
     &,L_USE_SEASALT_DIRECT                                             &
                                         !T if SS dir. rad. effect.
     &,L_USE_SEASALT_INDIRECT                                           &
                                         !T if SS 1st indir. effect
     &,L_USE_SEASALT_AUTOCONV
                                         !T if SS 2nd indir. effect
!
! Arguments with intent IN/OUT:
      REAL                                                              &
     & SO2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &                                model_levels)                     &
                                                         !mmr S in SO2
     &,DMS(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &                                model_levels)                     &
                                                         !mmr S in DMS
     &,SO4_AIT(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &                                model_levels)                     &
                                                         !mmr S in AIT
     &,SO4_ACC(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &                                model_levels)                     &
                                                         !mmr S in ACC
     &,SO4_DIS(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &                                model_levels)                     &
                                                         !mmr S in DIS
     &,NH3(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &                                model_levels)                     &
                                                         !mmr N in NH3
     &,H2O2_MXR(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &                                model_levels)                     &
                                                         !mmr H2O2
     &,SOOT_NEW(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &                                model_levels)                     &
                                                      !mmr fresh soot
     &,SOOT_AGD(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &                                model_levels)                     &
                                                      !mmr aged soot
     &,SOOT_CLD(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &                                model_levels)                     &
                                                      !mmr sootincloud
     &,BMASS_NEW(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                model_levels)                     &
                                                      !mmr fresh smoke
     &,BMASS_AGD(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                model_levels)                     &
                                                      !mmr aged smoke
     &,BMASS_CLD(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                model_levels)                     &
                                                      !mmr cloud smoke
     &,OCFF_NEW(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &                                model_levels)                     &
                                                      !mmr fresh OCFF
     &,OCFF_AGD(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &                                model_levels)                     &
                                                      !mmr aged OCFF
     &,OCFF_CLD(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &                                model_levels)                     &
                                                      !mmr OCFF incloud
     &,BIOGENIC(row_length, rows, model_levels)       !mmr biogenics
!
! Arguments with intent IN:
      REAL, INTENT(IN) ::                                               &
     & DUST_DIV1(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                model_levels)                     &   
                                                      !mmr Dust div 1
     &,DUST_DIV2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                model_levels)                     &      
                                                      !mmr Dust div 2
     &,DUST_DIV3(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &         
     &                                model_levels)                     &   
                                                      !mmr Dust div 3
     &,DUST_DIV4(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &         
     &                                model_levels)                     &
                                                      !mmr Dust div 4
     &,DUST_DIV5(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &         
     &                                model_levels)                        
                                                      !mmr Dust div 5
! arguments with intent OUT:
      Integer                                                           &
     &  Error_code
!
#include "c_sulchm.h"
!
! Variables with intent OUT (diagnostics):
      REAL                                                              &
     & MSA(row_length,rows,model_levels)                                &
                                                  !mmr S in MSA
     &,F_DMS_TO_SO2(row_length,rows,model_levels)                       &
                                                  !frac oxid DMS to SO2
     &,F_DMS_TO_SO4(row_length,rows,model_levels)                       &
                                                  !frac oxid DMS to SO4
     &,F_DMS_TO_MSA(row_length,rows,model_levels)                       &
                                                  !frac oxid DMS to MSA
     &,NH3_DEP(row_length,rows,model_levels)                            &
                                                  !NH3 depleted
     &,DELTAS_DRY(row_length,rows,model_levels)                         &
                                                   !SO2 dry ox per ts
     &,DELTAS_WET(row_length,rows,model_levels)                         &
                                                   !SO2 wet ox by H2O2
     &,DELTAS_WET_O3(row_length,rows,model_levels)                      &
                                                   !SO2 wet ox by O3
     &,DELTAS_TOT(row_length,rows,model_levels)                         &
                                                   !total SO2 ox per ts
     &,DELTAS_DMS(row_length,rows,model_levels)                         &
                                                   !DMS dry ox per ts
     &,DELTAS_EVAP(row_length,rows,model_levels)                        &
                                                   !SO4_DIS released by
!                             evapn of cloud droplets to SO4_ACC per ts
     &,DELTAS_NUCL(row_length,rows,model_levels)                        &
                                                   !SO4_ACC transfd by
!                                          nucleation to SO4_DIS per ts
     &,DELTAS_DIFFUSE(row_length,rows,model_levels)                     &
                                                    !SO4_AIT transfd to
!                                           SO4_DIS by diffusion per ts
#if defined(A17_2B)
     &,DELTAS_MERGE(row_length,rows,model_levels)                       &
#endif
     &,DELTAS_COAG(row_length,rows,model_levels)                        &
     &,PSI(row_length,rows,model_levels)
!
! Local variables
!
      Integer i,j,k,n                                                   &
     &, i_start                                                         &
                        ! Row start point for polar row tidying
     &, istat           ! Status (error code) indicator
!
      REAL CHEMSTEP                         ! chemistry timestep
! Diagnostic increments generated by each call of SULPHR
      REAL                                                              &
     & MSA_inc(row_length,rows,model_levels)                            &
     &,NH3_DEP_inc(row_length,rows,model_levels)
!
!
      REAL T(row_length,rows,model_levels)                              &
                                            ! Temp (K) on theta levels
     &,    T_surf(row_length,rows)          ! Surface temperature (K)
!
! For calculation of cos zenith angle
      REAL COSZA2D(row_length,rows)      !cos zenith angle
      Real                                                              &
     &  Sindec                                                          &
                                         ! sin(solar declination)
     &, hour_angle(row_length,rows)                                     &
     &, SCS                                                             &
                                         ! solar constant scaling factor
     &, seconds_since_midnight                                          &
     &, Eq_Time                                                         &
                                         ! The Equation of time
     &, day_fraction(row_length,rows)                                   &
     &, sin_true_latitude(row_length, rows)
!
! For sea-salt aerosol and DMS emissions calculations
      Real                                                              &
     &  u_1(aero_dim1,aero_dim2)                                        &
                                              ! surface wind u
     &, v_1(aero_dim1,aero_dim2)                                        &
                                              ! surface wind v
     &, u_1_mean                                                        &
                                              ! mean of u_1 for N Pole
     &, v_1_mean                                                        &
                                              ! mean of v_1 for N Pole
     &, sea_salt_film(salt_dim1,salt_dim2,salt_dim3)                    &
                                                     ! Sea-salt
     &, sea_salt_jet(salt_dim1,salt_dim2,salt_dim3)                     &
                                                     !    aerosol.
     &, height(aero_dim1,aero_dim2,aero_dim3)                           &
                                              ! Layer centre heights
     &, DMS_em_inter(aero_dim1,aero_dim2)                               &
                                              ! Interactive DMS emissns
     &, DMS_Ointer(aero_dim1,aero_dim2)                                 &
                                            ! Interactive DMS emissns
     &, DMS_Ointer_mean                                                 &
                                            ! N Pole mean DMS_Ointer
     &, DMS_em_combined(row_length,rows)      ! Combined DMS emissns
!
! Increments from soot ageing and diffusional scavenging
      Real                                                              &
     &  delta_agesoot(row_length,rows,model_levels)                     &
                                   ! Increment from soot ageing
     & ,delta_sootdiffscav(row_length,rows,model_levels)
                                   ! Increment from soot diff. scav.
!
! Increments from biomass smoke ageing and nucleation scavenging
      Real                                                              &
     &  delta_agebmass(row_length,rows,model_levels)                    &
                                    ! Increment from smoke ageing
     & ,delta_bmassnuclscav(row_length,rows,model_levels)
                                    ! Increment from smoke nucl. scav.
! For emissions of biomass smoke
      Integer biomass_level         ! Model level of biomass emissions
      Real biomass_per_level(row_length,rows)
                                    ! Quantity of biomass smoke emitted
                                    ! per model level
!
! Increments from fossil-fuel organic carbon ageing and nucl. scavenging
      Real                                                              &
     &  delta_ageocff(row_length, rows, model_levels)                   &
                                    ! Increment from OCFF ageing
     & ,delta_ocffnuclscav(row_length, rows, model_levels)
                                    ! Increment from OCFF nucl. scav.

!
! For cloud fraction without halos
      Real                                                              &
     &  cloudf(row_length,rows,wet_model_levels)
!
      LOGICAL L_SULPC_NEWDMS        ! TRUE if new DMS scheme required
!
      PARAMETER (L_SULPC_NEWDMS = .TRUE.)
!
! For droplet number calculation
      Real n_droplet(row_length, rows, wet_model_levels)                &
     &,    rho_air                                                      &
     &,    number_droplet                                               &
     &,    ss_film, ss_jet                                              &
     &,    sulp_ait, sulp_acc, sulp_dis                                 &
     &,    agd_bmass, cld_bmass                                         &
     &,    agd_ocff, cld_ocff                                           &
     &,    biogenic_ccn
!
! For PM10 & PM2.5 calculation
      Real                                                              &
     &  PM10(row_length, rows, model_levels)                            &
     &, PM2p5(row_length, rows, model_levels)                           &
     &, PM10_SO4 (row_length, rows, model_levels)                       &
     &, PM2p5_SO4(row_length, rows, model_levels)                       &
     &, PM10_BC (row_length, rows, model_levels)                        &
     &, PM2p5_BC(row_length, rows, model_levels)                        &
     &, PM10_BB (row_length, rows, model_levels)                        &
     &, PM2p5_BB(row_length, rows, model_levels)                        &
     &, PM10_OCFF (row_length, rows, model_levels)                      &
     &, PM2p5_OCFF(row_length, rows, model_levels)                      &
     &, PM10_SOA (row_length, rows, model_levels)                       &
     &, PM2p5_SOA(row_length, rows, model_levels)                       &
     &, PM10_SS (row_length, rows, model_levels)                        &
     &, PM2p5_SS(row_length, rows, model_levels)                        &
     &, PM10_DUST (row_length, rows, model_levels)                      &
     &, PM2p5_DUST(row_length, rows, model_levels)
!
! External functions and subroutines
      External TRSRCE, SOLPOS, SOLANG, SET_SEASALT, SULPHR              &
     &, AGESOOT, SOOTDIFFSCAV                                           &
     &, DMS_FLUX                                                        &
     &, AGEBMASS, BMASSNUCLSCAV                                         &
     &, AGEOCFF, OCFFNUCLSCAV                                           &
     &, gcg_rvecsumr, gcg_rvecsumf                                      &
     &, diagnostics_aero                                                &
     &, number_droplet                                                  & 
     &, pm10_pm2p5
!
! Set up global land fraction field
!
      Do j = 1, rows
        Do i = 1, row_length
          land_fract(i,j) = 0.0
        End Do
      End Do

      If (L_ctile) Then
        Do k = 1,land_points
          j = (land_index(k)-1)/row_length + 1
          i = land_index(k) - (j-1)*row_length
          land_fract(i,j) = fland_ctile(k)
        End Do
      Else
        Do j = 1, rows
          Do i = 1, row_length
            If (LAND_MASK(i,j)) Then
              land_fract(i,j) = 1.0
            Else
              land_fract(i,j) = 0.0
            End If
          End Do
        End Do
      End If

!
! Copy cloud fraction into an array without halos and set to zero if
! smaller than machine precision
          Do k = 1,wet_model_levels
            Do j = 1,rows
              Do i = 1,row_length
                cloudf(i,j,k)=cloudf_halos(i,j,k)
                if (cloudf(i,j,k) <  epsilon(cloudf(i,j,k)))            &
     &            cloudf(i,j,k)=0.0
                if (cloudf(i,j,k) >  1.0) then
                  cloudf(i,j,k) = 1.0
                end if
              End Do
            End Do
          End Do
!
! If sea-salt aerosol or interactive DMS emissions are to be
! used then calculate surface windspeed and height information
!
      If (L_seasalt_CCN .OR. L_DMS_em_inter) Then
!
        Do j = 1,rows
          Do i = 1,row_length
            u_1(i,j) = u(i,j,1)
            v_1(i,j) = v(i,j,1)
          End Do
        End Do
!
! Tidy up at North Pole; not done at South Pole as it's land.
! Polar values are mean of 1st interior row:
        If (at_extremity(PNorth)) Then
! Start point of first interior (i.e. non-polar) row:
          i_start  = (rows - 2) * row_length + 1
! Sum over points on PEs in order along first interior row:
#if defined(REPROD)
          Call gcg_rvecsumr(row_length*rows, row_length, i_start, 1,    &
     &                      u_1, proc_row_group, istat, u_1_mean)
          Call gcg_rvecsumr(row_length*rows, row_length, i_start, 1,    &
     &                      v_1, proc_row_group, istat, v_1_mean)
#else
          Call gcg_rvecsumf(row_length*rows, row_length, i_start, 1,    &
     &                      u_1, proc_row_group, istat, u_1_mean)
          Call gcg_rvecsumf(row_length*rows, row_length, i_start, 1,    &
     &                      v_1, proc_row_group, istat, v_1_mean)
#endif
          u_1_mean=u_1_mean/global_row_length
          v_1_mean=v_1_mean/global_row_length
          Do i = 1, row_length
            u_1(i, rows) = u_1_mean
            v_1(i, rows) = v_1_mean
          End Do
        Endif
!
          Do k = 1, aero_dim3
            Do j = 1, aero_dim2
              Do i = 1, aero_dim1
                height(i,j,k) = r_theta_levels(i,j,k)                   &
     &                                          -r_theta_levels(i,j,0)
              End Do
            End Do
          End Do
!
      End If
      If (L_DMS_Ointer) Then
!
! Tidy up ocean DMS flux at North Pole, as done above for u_1,v_1.
! Polar values are mean of 1st interior row:
        If (at_extremity(PNorth)) Then
! Start point of first interior (i.e. non-polar) row:
          i_start  = (rows - 2) * row_length + 1
! Sum over points on PEs in order along first interior row:
#if defined(REPROD)
          Call gcg_rvecsumr(row_length*rows, row_length, i_start, 1,    &
     &         dms_Ointer, proc_row_group, istat, dms_Ointer_mean)
#else
          Call gcg_rvecsumf(row_length*rows, row_length, i_start, 1,    &
     &         dms_Ointer, proc_row_group, istat, dms_Ointer_mean)
#endif
          dms_Ointer_mean=dms_Ointer_mean/global_row_length
          Do i = 1, row_length
            dms_Ointer(i, rows) = dms_Ointer_mean
          End Do
        Endif
      Endif ! L_DMS_Ointer
!

!
      If (L_seasalt_CCN) Then
! DEPENDS ON: set_seasalt
          CALL SET_SEASALT(u_1, v_1, height, land_fract, ice_fract      &
     &                    , row_length, rows, model_levels              &
     &                    , bl_levels, sea_salt_film, sea_salt_jet)
      Else
        sea_salt_film(1,1,1) = 0.0
        sea_salt_jet(1,1,1) = 0.0
      Endif
!
! Calculate T from theta if soot, S cycle, biomass or OCFF switched on
! and droplet numbers for diffusional scavenging.
!
      IF (L_SULPC_SO2 .OR. L_SOOT .OR. L_BIOMASS .OR. L_OCFF) THEN
        Do k=1,model_levels
          Do j=1,rows
            Do i=1,row_length
              T(i,j,k)=theta(i,j,k)*exner_theta_levels(i,j,k)
            End Do
          End Do
        End Do
        Do k=1,wet_model_levels
          Do j=1,rows
            Do i=1,row_length
              rho_air=p_theta_levels(i,j,k)/                            &
     &                (R * T(i,j,k))
!
              If (L_seasalt_CCN) Then
                ss_film=sea_salt_film(i,j,k)
                ss_jet=sea_salt_jet(i,j,k)
              Else
                ss_film=0.0
                ss_jet=0.0
              End If
!
              If (L_sulphate_CCN) Then
                sulp_ait=so4_ait(i,j,k)
                sulp_acc=so4_acc(i,j,k)
                sulp_dis=so4_dis(i,j,k)
              Else
                sulp_ait=0.0
                sulp_acc=0.0
                sulp_dis=0.0
              End If
!
              If (L_biomass_CCN) Then
                agd_bmass=bmass_agd(i,j,k)
                cld_bmass=bmass_cld(i,j,k)
              Else
                agd_bmass=0.0
                cld_bmass=0.0
              End If
!
              If (L_ocff_CCN) Then
                agd_ocff=ocff_agd(i,j,k)
                cld_ocff=ocff_cld(i,j,k)
              Else
                agd_ocff=0.0
                cld_ocff=0.0
              End If
!
              If (L_USE_BIOGENIC) Then
                biogenic_ccn=BIOGENIC(i,j,k)
              Else
                biogenic_ccn=0.0
              End If
!
! DEPENDS ON: number_droplet
              n_droplet(i,j,k)=NUMBER_DROPLET(                          &
     &                 L_sulphate_CCN,.FALSE.                           &
     &                ,sulp_ait,sulp_acc,sulp_dis                       &
     &                ,L_seasalt_CCN,ss_film,ss_jet                     &
     &                ,L_USE_BIOGENIC,biogenic_ccn                      &
     &                ,L_biomass_CCN,agd_bmass,cld_bmass                &
     &                ,L_ocff_CCN,agd_ocff,cld_ocff                     &
     &                ,rho_air,snow_depth(i,j),land_fract(i,j)          &
     &                ,Ntot_land,Ntot_sea                               &
     &                 )
!
            End Do
          End Do
        End Do
      END IF
!
! Soot cycle code.
! Calculate emissions of soot from surface and chimney level, ageing
! of soot and diffusional scavenging of aged soot to soot-in-cloud.
!
      IF (L_SOOT) THEN  ! If soot modelling is included

        If (L_SOOT_SUREM) then  ! If surface soot emissions are included
! DEPENDS ON: trsrce
          CALL TRSRCE(                                                  &
     &    rows, row_length, off_x, off_y, halo_i, halo_j                &
     &,   model_levels, wet_model_levels                                &
     &,   halo_i, halo_j                                                &
     &,   r_rho_levels, r_theta_levels                                  &
     &,   theta, q , qcl , qcf , exner_rho_levels, rho                  &
     &,   SOOT_NEW(:,:,1), SOOT_SUREM, 1                                &
     &,   TIMESTEP, 1, 1, 0.0                                           &
     &    )
        End If  ! L_SOOT_SUREM
!
        If (L_SOOT_HILEM) then  ! If chimney level soot emissions
                                ! are included
! DEPENDS ON: trsrce
          CALL TRSRCE(                                                  &
     &    rows, row_length, off_x, off_y, halo_i, halo_j                &
     &,   model_levels, wet_model_levels                                &
     &,   halo_i, halo_j                                                &
     &,   r_rho_levels, r_theta_levels                                  &
     &,   theta, q , qcl , qcf , exner_rho_levels, rho                  &
     &,   SOOT_NEW(:,:,SOOT_HIGH_LEVEL), SOOT_HILEM                     &
     &,   SOOT_HIGH_LEVEL, TIMESTEP, 1, 1, 0.0                          &
     &    )
        End If  ! L_SOOT_HILEM
!
! Calculate quantity of fresh soot converted to aged soot
!
! DEPENDS ON: agesoot
        CALL AGESOOT(                                                   &
     &  row_length, rows, off_x, off_y,                                 &
     &  model_levels, timestep,                                         &
     &  soot_new,                                                       &
     &  delta_agesoot                                                   &
     &  )
!
! Calculate quantity of aged soot scavenged to cloud soot
!
! DEPENDS ON: sootdiffscav
        CALL SOOTDIFFSCAV(                                              &
     &  rows, row_length, off_x, off_y, halo_i, halo_j,                 &
     &  model_levels, wet_model_levels, timestep,                       &
     &  cloudf, qcl, qcf, p_theta_levels, t,                            &
     &  n_droplet,                                                      &
     &  soot_agd, soot_cld,                                             &
     &  delta_sootdiffscav                                              &
     &  )
!
! Update soot arrays with increments from ageing and
! diffusional scavenging
!
        Do k=1,model_levels
          Do j=1,rows
            Do i=1,row_length
              soot_new(i,j,k)=soot_new(i,j,k) - delta_agesoot(i,j,k)
              soot_agd(i,j,k)=soot_agd(i,j,k) + delta_agesoot(i,j,k)    &
     &                               - delta_sootdiffscav(i,j,k)
              soot_cld(i,j,k)=soot_cld(i,j,k)                           &
     &                               + delta_sootdiffscav(i,j,k)
            End Do
          End Do
        End Do
!
      End If  ! L_SOOT
!
! End of soot cycle code
!
! Biomass aerosol code.
! Calculate emissions of smoke from surface and high level, ageing
! of smoke and nucleation scavenging of aged smoke to smoke-in-cloud.
!
      IF (L_BIOMASS) THEN  ! If biomass smoke modelling is included

        If (L_BMASS_SUREM) then  ! If surface smoke emissions included
! DEPENDS ON: trsrce
          CALL TRSRCE(                                                  &
     &    rows, row_length, off_x, off_y, halo_i, halo_j                &
     &,   model_levels, wet_model_levels                                &
     &,   halo_i, halo_j                                                &
     &,   r_rho_levels, r_theta_levels                                  &
     &,   theta, q , qcl , qcf , exner_rho_levels, rho                  &
     &,   BMASS_NEW(:,:,1), BMASS_SUREM, 1                              &
     &,   TIMESTEP, 1, 1, 0.0                                           &
     &    )
        End If  ! L_BMASS_SUREM
!
        If (L_BMASS_HILEM) then  ! If high level smoke emissions
                                 ! are included
!
! Check that the range of emission levels is correctly specified
! and, if so, emit equal quantities of smoke on all model levels
! between bmass_high_level_1 and bmass_high_level_2.
!
          If (bmass_high_level_1  >   bmass_high_level_2 .or.           &
     &        bmass_high_level_1  >   model_levels .or.                 &
     &        bmass_high_level_2  >   model_levels) then
            Write(6,*) 'Aero_Ctl: Invalid range of biomass emission '
            Write(6,*) 'levels specified.'
            Write(6,*) 'Lowest level: ',bmass_high_level_1
            Write(6,*) 'Highest level: ',bmass_high_level_2
          Else
!
            Do j=1,rows
              Do i=1,row_length
                biomass_per_level(i,j) = bmass_hilem(i,j)/              &
     &          ((bmass_high_level_2 - bmass_high_level_1) + 1)
              End Do
            End Do
!
            Do biomass_level = bmass_high_level_1, bmass_high_level_2
!
! DEPENDS ON: trsrce
              CALL TRSRCE(                                              &
     &        rows, row_length, off_x, off_y, halo_i, halo_j            &
     &,       model_levels, wet_model_levels                            &
     &,       halo_i, halo_j                                            &
     &,       r_rho_levels, r_theta_levels                              &
     &,       theta, q , qcl , qcf , exner_rho_levels, rho              &
     &,       bmass_new(:,:,biomass_level), biomass_per_level           &
     &,       biomass_level, timestep, 1, 1, 0.0                        &
     &        )
!
            End Do  ! biomass_level
!
          End If  ! Test on emissions levels
!
        End If  ! L_BMASS_HILEM
!
! Calculate quantity of fresh smoke converted to aged smoke
!
! DEPENDS ON: agebmass
        CALL AGEBMASS(                                                  &
     &  row_length, rows, off_x, off_y,                                 &
     &  model_levels, timestep,                                         &
     &  bmass_new,                                                      &
     &  delta_agebmass                                                  &
     &  )
!
! Calculate quantity of aged smoke scavenged to cloud smoke
!
! DEPENDS ON: bmassnuclscav
        CALL BMASSNUCLSCAV(                                             &
     &  rows, row_length, off_x, off_y, halo_i, halo_j,                 &
     &  model_levels, wet_model_levels, timestep,                       &
     &  cloudf, qcl, qcf,                                               &
     &  bmass_agd, bmass_cld,                                           &
     &  delta_bmassnuclscav                                             &
     &  )
!
! Update smoke arrays with increments from ageing and
! nucleation scavenging
!
        Do k=1,model_levels
          Do j=1,rows
            Do i=1,row_length
              bmass_new(i,j,k)=bmass_new(i,j,k) - delta_agebmass(i,j,k)
#if defined(A17_2A)
              bmass_agd(i,j,k)=bmass_agd(i,j,k) + delta_agebmass(i,j,k) &
     &                               - delta_bmassnuclscav(i,j,k)
#endif
#if defined(A17_2B)
!    Simulate the condensation of VOCs onto aged biomass by
!    increasing the mass transferred upon ageing by (8.75/5.4),
!    the ratio of the fraction of BC in fresh (8.75%) and aged
!    (5.4%) biomass aerosol:
              bmass_agd(i,j,k)=bmass_agd(i,j,k)                         &
     &                        + ((8.75/5.4)*delta_agebmass(i,j,k))      &
     &                               - delta_bmassnuclscav(i,j,k)
#endif
              bmass_cld(i,j,k)=bmass_cld(i,j,k)                         &
     &                               + delta_bmassnuclscav(i,j,k)
            End Do
          End Do
        End Do
!
      End If  ! L_BIOMASS
!
! End of biomass aerosol code
!
!
      IF ( L_SULPC_SO2 ) THEN
!
! Calculate cos zenith angle for SULPH2 by calling SOLPOS and SOLANG
! Calculate number of seconds since midnight to the beginning of
! the timetsep.
!
      seconds_since_midnight = Real( PREVIOUS_TIME(4) * 3600            &
     &         + PREVIOUS_TIME(5) * 60  + PREVIOUS_TIME(6))

! Calculate true latitude at each grid-point.
        Do j = 1, rows
          Do i = 1, row_length
            sin_true_latitude(i,j) = f3_at_u(i,j) / two_Omega
          End Do
        End Do
!
! DEPENDS ON: solpos
        CALL SOLPOS (PREVIOUS_TIME(7), PREVIOUS_TIME(1),                &
     &               L_CAL360, L_SEC_VAR, L_EqT, Eq_Time, Sindec, SCS)
!
! DEPENDS ON: solang
        CALL SOLANG(                                                    &
! input constants
     &                Sindec, seconds_since_midnight,                   &
     &                timestep, Eq_Time,                                &
! row and column dependent constants
     &                sin_true_latitude,                                &
     &                true_longitude,                                   &
! size variables
     &                row_length*rows,                                  &
! output fields
     &                day_fraction, COSZA2D, hour_angle )
!
!
! Calculate interactive DMS emissions if required
!
       If (L_DMS_em_inter) Then

         Do j = 1, rows
           Do i = 1, row_length
             If (L_ctile) Then
               T_surf(i, j) = Tstar_sea(i, j)
             Else
               T_surf(i, j) = Tstar(i, j)
             End If
           End Do
         End Do

! DEPENDS ON: dms_flux
         CALL DMS_FLUX (                                                &
! size variables
     &                row_length, rows                                  &
! input fields
     &,               u_1, v_1                                          &
     &,               height(:,:,1)                                     &
     &,               T_surf                                            &
     &,               land_fract                                        &
     &,               DMS_conc                                          &
! logical switches
     &,               L_DMS_Liss_Merlivat                               &
     &,               L_DMS_Wanninkhof                                  &
     &,               L_DMS_Nightingale                                 &
! output field
     &,               DMS_em_inter )

       End If
!
!
! Zero output diagnostics for full model timestep
      Do k=1,model_levels
        Do j=1,rows
          Do i=1,row_length
            MSA(i,j,k) = 0.0
            NH3_DEP(i,j,k) = 0.0
          End Do
        End Do
      End Do
!
! Calculate length of chemistry timestep and use it to control input
! of emissions and S Chemistry
!
      CHEMSTEP=timestep/CALL_CHEM_FREQ
!
        Do n=1,CALL_CHEM_FREQ
!
! Call TRSRCE to insert emissions
!
      If (L_SO2_SURFEM) Then         ! Insert surface SO2 emiss
! DEPENDS ON: trsrce
        CALL TRSRCE(                                                    &
     &  rows, row_length, off_x, off_y, halo_i, halo_j                  &
     &, model_levels, wet_model_levels                                  &
     &, halo_i, halo_j                                                  &
     &, r_rho_levels, r_theta_levels                                    &
     &, theta, q , qcl , qcf , exner_rho_levels, rho                    &
     &, SO2(:,:,1), SO2_SURFEM, 1                                       &
     &, CHEMSTEP, 1, 1, 0.0                                             &
     & )
      End If
!
      If (L_SO2_HILEM) Then          ! Insert chimney SO2 emiss
! DEPENDS ON: trsrce
        CALL TRSRCE(                                                    &
     &  rows, row_length, off_x, off_y, halo_i, halo_j                  &
     &, model_levels, wet_model_levels                                  &
     &, halo_i, halo_j                                                  &
     &, r_rho_levels, r_theta_levels                                    &
     &, theta, q , qcl , qcf , exner_rho_levels, rho                    &
     &, SO2(:,:,SO2_HIGH_LEVEL), SO2_HILEM, SO2_HIGH_LEVEL              &
     &, CHEMSTEP, 1, 1, 0.0                                             &
     & )
      End If
!
      If (L_SO2_NATEM) Then          ! Insert volcanic SO2 emiss
        Do k= 1,model_levels
! DEPENDS ON: trsrce
          CALL TRSRCE(                                                  &
     &  rows, row_length, off_x, off_y, halo_i, halo_j                  &
     &, model_levels, wet_model_levels                                  &
     &, halo_i, halo_j                                                  &
     &, r_rho_levels, r_theta_levels                                    &
     &, theta, q , qcl , qcf , exner_rho_levels, rho                    &
     &, SO2(:,:,k), SO2_NATEM(:,:,k), k                                 &
     &, CHEMSTEP, 1, 1, 0.0                                             &
     & )
        End Do
      End If
!
! Merge ocean DMS emissions with land DMS emissions
      If (L_DMS_Ointer) Then
        DMS_em_inter(:,:) = DMS_Ointer(:,:)
      End If

      If (L_DMS_em_inter .OR. L_DMS_Ointer) Then
        Do j = 1, rows
          Do i = 1, row_length
            DMS_em_combined(i, j)                                       &
     &          = (land_fract(i, j) * DMS_em_ancil(i, j))               &
     &          + ( ((1.0 - land_fract(i, j)) * DMS_em_inter(i, j))     &
     &             * (1.0 - ice_fract(i, j)) )
          End Do
        End Do
      Else
        If (L_DMS_em) Then     ! Just copy over the standard ancil
          Do j = 1, rows
            Do i = 1, row_length
              DMS_em_combined(i, j) = DMS_em_ancil(i, j)
            End Do
          End Do
        End If
      End If
!
      If (L_DMS_EM) Then             ! Insert DMS emiss (surface)
! DEPENDS ON: trsrce
        CALL TRSRCE(                                                    &
     &  rows, row_length, off_x, off_y, halo_i, halo_j                  &
     &, model_levels, wet_model_levels                                  &
     &, halo_i, halo_j                                                  &
     &, r_rho_levels, r_theta_levels                                    &
     &, theta, q , qcl , qcf , exner_rho_levels, rho                    &
     &, DMS(:,:,1), DMS_em_combined, 1                                  &
     &, CHEMSTEP, 1, 1, 0.0                                             &
     & )
      End If
!
      If (L_NH3_EM) Then             ! Insert NH3 emiss (surface)
! DEPENDS ON: trsrce
        CALL TRSRCE(                                                    &
     &  rows, row_length, off_x, off_y, halo_i, halo_j                  &
     &, model_levels, wet_model_levels                                  &
     &, halo_i, halo_j                                                  &
     &, r_rho_levels, r_theta_levels                                    &
     &, theta, q , qcl , qcf , exner_rho_levels, rho                    &
     &, NH3(:,:,1), NH3_EM, 1                                           &
     &, CHEMSTEP, 1, 1, 0.0                                             &


     & )
      End If
!
!
! DEPENDS ON: sulphr
      CALL SULPHR(                                                      &
! Arguments IN
     &             halo_i, halo_j, off_x, off_y                         &
     &,            row_length, rows                                     &
     &,            model_levels, wet_model_levels                       &
     &, theta_field_size                                                &
     &,            CHEMSTEP                                             &
     &,            CLOUDF, COSZA2D                                      &
     &,            p_theta_levels, T, q, qcl, qcf                       &
     &,            OH_CONC, H2O2_LMT, HO2_CONC, O3                      &
     &,            n_droplet                                            &
     &,            L_SULPC_DMS, L_SULPC_NEWDMS                          &
     &,            L_SULPC_OZONE                                        &
#if defined(A17_2B)
     &,            L_SULPC_SO2_O3_NONBUFFERED                           &
#endif
     &,            L_SULPC_NH3                                          &
! Arguments IN/OUT
     &,            SO2, DMS, SO4_AIT, SO4_ACC, SO4_DIS                  &
     &,            NH3, H2O2_MXR                                        &
! Arguments OUT (diagnostics)
     &,            MSA_inc                                              &
     &,            NH3_DEP_inc                                          &
     &,            F_DMS_TO_SO2, F_DMS_TO_SO4, F_DMS_TO_MSA             &
     &,            DELTAS_DRY, DELTAS_WET, DELTAS_WET_O3                &
     &,            DELTAS_TOT, DELTAS_DMS                               &
     &,            DELTAS_EVAP, DELTAS_NUCL, DELTAS_DIFFUSE             &
#if defined(A17_2B)
     &,            DELTAS_MERGE                                         &
#endif
     &,            DELTAS_COAG, PSI                                     &
     &            )
!
! Add diagnostic increments from SULPHR to total for model timestep
!
      Do k=1,model_levels
        Do j=1,rows
          Do i=1,row_length
            MSA(i,j,k)=MSA(i,j,k)+MSA_inc(i,j,k)
          NH3_DEP(i,j,k)=NH3_DEP(i,j,k)+NH3_DEP_inc(i,j,k)
          End Do
        End Do
      End Do
!
!
        End Do            ! End CALL_CHEM_FREQ loop
!
      ENDIF               ! End L_SULPC_SO2 test

!
! OCFF cycle code.
! Calculate emissions of ocff from surface and chimney level, ageing
! of ocff and nucleation scavenging of aged ocff to ocff-in-cloud.
!
      IF (L_OCFF) THEN  ! If ocff modelling is included

        If (L_OCFF_SUREM) then  ! If surface ocff emissions are included
! DEPENDS ON: trsrce
          CALL TRSRCE(                                                  &
     &    rows, row_length, off_x, off_y, halo_i, halo_j                &
     &,   model_levels, wet_model_levels                                &
     &,   halo_i, halo_j                                                &
     &,   r_rho_levels, r_theta_levels                                  &
     &,   theta, q , qcl , qcf , exner_rho_levels, rho                  &
     &,   OCFF_NEW(:,:,1), OCFF_SUREM, 1                                &
     &,   TIMESTEP, 1, 1, 0.0                                           &
     &    )
        End If  ! L_OCFF_SUREM
!
        If (L_OCFF_HILEM) then  ! If chimney level ocff emissions
                                ! are included
! DEPENDS ON: trsrce
          CALL TRSRCE(                                                  &
     &    rows, row_length, off_x, off_y, halo_i, halo_j                &
     &,   model_levels, wet_model_levels                                &
     &,   halo_i, halo_j                                                &
     &,   r_rho_levels, r_theta_levels                                  &
     &,   theta, q , qcl , qcf , exner_rho_levels, rho                  &
     &,   OCFF_NEW(:,:,OCFF_HIGH_LEVEL), OCFF_HILEM                     &
     &,   OCFF_HIGH_LEVEL, TIMESTEP, 1, 1, 0.0                          &
     &    )
        End If  ! L_OCFF_HILEM
!
! Calculate quantity of fresh ocff converted to aged ocff
!
! DEPENDS ON: ageocff
        CALL AGEOCFF(                                                   &
     &  row_length, rows, off_x, off_y,                                 &
     &  model_levels, timestep,                                         &
     &  ocff_new,                                                       &
     &  delta_ageocff                                                   &
     &  )
!
! Calculate quantity of aged ocff scavenged to cloud ocff
!
! DEPENDS ON: ocffnuclscav
        CALL OCFFNUCLSCAV(                                              &
     &  rows, row_length, off_x, off_y, halo_i, halo_j,                 &
     &  model_levels, wet_model_levels, timestep,                       &
     &  cloudf, qcl, qcf,                                               &
     &  ocff_agd, ocff_cld,                                             &
     &  delta_ocffnuclscav                                              &
     &  )
!
! Update ocff arrays with increments from ageing and
! nucleation scavenging
!
        Do k=1,model_levels
          Do j=1,rows
            Do i=1,row_length
              ocff_new(i,j,k)=ocff_new(i,j,k) - delta_ageocff(i,j,k)
              ocff_agd(i,j,k)=ocff_agd(i,j,k) + delta_ageocff(i,j,k)    &
     &                               - delta_ocffnuclscav(i,j,k)
              ocff_cld(i,j,k)=ocff_cld(i,j,k)                           &
     &                               + delta_ocffnuclscav(i,j,k)
            End Do
          End Do
        End Do
!
      End If  ! L_OCFF
!
! End of ocff cycle code
!

!
! PM10 / PM2.5 code: If any of the diagnostics on PM10 and PM2.5
! mass concentrations (item numbers 220-235, Sect. 17) is requested
! then call the subroutine that considers all aerosol species and modes, 
! to calculate their contributions to PM concentrations

      If (sf(220,17) .OR. sf(221,17) .OR. sf(222,17) .OR.               & 
          sf(223,17) .OR. sf(224,17) .OR. sf(225,17) .OR.               &  
          sf(226,17) .OR. sf(227,17) .OR. sf(228,17) .OR.               & 
          sf(229,17) .OR. sf(230,17) .OR. sf(231,17) .OR.               & 
          sf(232,17) .OR. sf(233,17) .OR. sf(234,17) .OR.               & 
          sf(235,17)) then

! DEPENDS ON: pm10_pm2p5
        Call pm10_pm2p5 (                                               &
          off_x, off_y,                                                 &
          row_length, rows,                                             &
          model_levels,                                                 &
          salt_dim1, salt_dim2, salt_dim3,                              &
          L_SULPC_SO2, L_SOOT, L_BIOMASS,                               & 
          L_OCFF, L_USE_BIOGENIC, L_DUST,                               &   
          L_seasalt_CCN, L_USE_SEASALT_DIRECT,                          & 
          L_USE_SEASALT_INDIRECT, L_USE_SEASALT_AUTOCONV,               &   
          p_theta_levels, T,                                            &
          SO4_AIT, SO4_ACC, SOOT_NEW, SOOT_AGD, BMASS_NEW, BMASS_AGD,   &
          OCFF_NEW, OCFF_AGD, BIOGENIC, sea_salt_film, sea_salt_jet,    &
          DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5,        &
          PM10, PM2p5,                                                  &
          PM10_SO4, PM2p5_SO4,                                          &
          PM10_BC, PM2p5_BC,                                            &
          PM10_BB, PM2p5_BB,                                            &
          PM10_OCFF, PM2p5_OCFF,                                        &
          PM10_SOA, PM2p5_SOA,                                          &
          PM10_SS, PM2p5_SS,                                            &
          PM10_DUST, PM2p5_DUST)
      
      End if

!
! Call diagnostic routine
!
! DEPENDS ON: diagnostics_aero
      Call diagnostics_aero(                                            &
     &                       row_length, rows, model_levels             &
     &,                      wet_model_levels                           &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      timestep                                   &
     &,                      at_extremity                               &
!
     &,                      L_SULPC_SO2, L_DMS_em                      &
     &,                      L_SULPC_DMS, L_SULPC_NEWDMS                &
     &,                      L_SULPC_OZONE, L_SULPC_NH3                 &
     &,                      L_SOOT                                     &
     &,                      MSA, NH3_DEP                               &
     &,                      DMS_em_combined                            &
     &,                      DELTAS_DMS                                 &
     &,                      F_DMS_TO_SO2                               &
     &,                      F_DMS_TO_SO4                               &
     &,                      F_DMS_TO_MSA                               &
     &,                      DELTAS_DRY                                 &
     &,                      DELTAS_WET                                 &
     &,                      DELTAS_WET_O3                              &
     &,                      DELTAS_EVAP                                &
     &,                      DELTAS_NUCL                                &
     &,                      DELTAS_DIFFUSE                             &
     &,                      DELTAS_COAG                                &
#if defined(A17_2B)
     &,                      DELTAS_MERGE                               &
#endif
     &,                      PSI                                        &
     &,                      PM10, PM2p5                                &
     &,                      PM10_SO4, PM2p5_SO4                        &
     &,                      PM10_BC, PM2p5_BC                          &
     &,                      PM10_BB, PM2p5_BB                          &
     &,                      PM10_OCFF, PM2p5_OCFF                      &
     &,                      PM10_SOA, PM2p5_SOA                        &
     &,                      PM10_SS, PM2p5_SS                          &
     &,                      PM10_DUST, PM2p5_DUST                      &
     &     ,                                                            &
#include "argsts.h"
     & STASHwork17                                                      &
     &  )
      RETURN
      END SUBROUTINE AERO_CTL
!
#endif
