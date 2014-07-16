#if defined(ATMOS)
#if defined(CONTROL) || defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine NI_gwd_ctl

      Subroutine NI_gwd_ctl (                                           &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y, global_row_length, global_rows    &
     &, proc_row_group, at_extremity, n_proc, n_procx, n_procy          &
     &, neighbour, g_rows, g_row_length, g_datastart, me                &

! model dimensions.
     &, row_length, rows, n_rows, land_points                           &
     &, model_levels                                                    &

! Model switches
     &, model_domain, gw_kay, gwd_frc                                   &
     &, Ltimer, l_gwd, L_use_ussp                                       &
     &, l_taus_scale, l_fix_gwsatn, l_gwd_40km, l_ussp_opaque           &
     &, sat_scheme, gwd_fsat                                            &

! trig arrays
     &, sin_theta_longitude, sin_theta_latitude                         &

! in coordinate information
     &, r_rho_levels, r_theta_levels                                    &
     &, eta_theta_levels, eta_rho_levels                                &

! in time stepping information.
     &, timestep, timestep_number                                       &

! diagnostic info
     &     ,                                                            &
#include "argsts.h"
     & STASHwork                                                        &

! in data fields.
     &, u, v                                                            &
     &, land_sea_mask, p_layer_boundaries                               &
     &, rho, theta_latest, sd_orog_land, orog_grad_xx_land              &
     &, orog_grad_xy_land, orog_grad_yy_land, land_index                &

! in/out
     &, R_u, R_v                                                        &

! error information
     &, Error_code  )

! purpose: Interface to Atmospheric Physics GWD Schemes.
!         Scheme 1: Orographic flow blocking and gravity wave scheme.
!         Scheme 2: Non-orographic ultra-simple spectral gw scheme.
!
!
! current code owner: S. Webster
!
! history:
! Version   Date     Comment
! ----   -------   -------
! 5.0  30/11/99 Original version. J-C Thil.
! 5.1  09/12/99 Call diagnostics routine only when diagnostics
!               requested + comment changes. Rick Rawlins
! 5.1 15/03/00 Correct diagnostic dimensioning and remove now
!              superfluous lists. R Rawlins
! 5.1  28/02/00 Provide explicit model increments as STASH output
!               diagnostics. R Rawlins
! 5.1  10/04/00 Remove setting GWD gw_p_ref in namelist. S. Webster
!
! 5.2  15/11/00 Replace call to gwave with call to gwd_glue to allow
!               3B or 4A schemes to be called. Update call to daggwd
!               to include all new diagnostics. Allocate memory
!               for all diagnostics only when required. Pass 4A
!               switches down from cruntimc.
!                                                         S. Webster
! 5.3  17/10/01 Changes required for Single Column Model
!                                             Z. Gardner
! 5.3  16/10/01 Pass l_use_ussp thru' routine rather than l_gwd_wake
!                                                         S. Webster
! 5.3  11/10/01 Add arguments to diagnostic_gwd to calculate
!               orographic standard deviation. D.M. Goddard
! 5.4  28/08/02 Add arrays for numerical limiter diagnostics. S.Webster
!
!  5.4   05/09/02   Add diagnostics for spectral (non-orographic)
!               gravity wave forcing scheme (GW_USSP).      Adam Scaife
!
! 5.5  25/02/03 Remove 3B GWD code and hence call 4A code from this
!               deck rather than via glue routines. Also remove other
!               redundant arrays.                            S. Webster
! 6.2  21/02/06 Pass l_taus_scale and l_fix_gwsatn through to GWAVE4A
!               and l_ussp_opaque thru' to GW_USSP         S. Webster
!
! 6.2  22/06/06 Fix dependency on diagnostics_gwd for SCM.
!               P.Selwood
! 6.2  21/02/06 Add orographic surface stress and
!               sigma_xx, xy and yy diagnostics.             S. Webster
! 6.3  11/10/06 Fix unmatched endif in SCMA extraction.  R Barnes
!
! 6.4  19/05/06 Correction of USSP diagnostic flux output.
!                                                   A.C. Bushell
!
! code description:
!   language: fortran 77 + cray extensions
!   this code is written to umdp3 programming standards.

      Implicit None

! arguments with intent in. ie: input variables.

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
     &, global_rows                                                     &
                           ! NUMBER OF global rows
     &, n_proc                                                          &
                   ! Total number of processors
     &, n_procx                                                         &
                   ! Number of processors in longitude
     &, n_procy                                                         &
                   ! Number of processors in latitude
     &, neighbour(4)                                                    &
                             ! Array with the Ids of the four neighbours
                             ! in the horizontal plane
     &, g_rows (0:n_proc-1)                                             &
     &, g_row_length (0:n_proc-1)                                       &
     &, g_datastart (3,0:n_proc-1)                                      &
     &, me         ! My processor number

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north, sout
                         ! east or west of the processor grid


! Model dimensions
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, n_rows                                                          &
     &, model_levels                                                    &
     &, land_points

! Model switches
      Integer                                                           &
     &  model_domain                                                    &
     &, sat_scheme

      Logical                                                           &
     &  Ltimer                                                          &
                 ! true then output some timing information
     &, L_gwd                                                           &
                    ! switch for the orographic GWD scheme
     &, L_use_ussp                                                      &
                    ! switch for the non-orographic USSP scheme
     &, L_taus_scale                                                    &
                      ! switch to scale orog surf stress by Froude No.
     &, L_fix_gwsatn                                                    &
                      ! switch to invoke minor bug fixes in gwsatn4a
     &, L_gwd_40km                                                      &
                      ! switch to turn off orographic GWD above 40km     
     &, L_ussp_opaque ! switch to change lid condition in gw_ussp

! model parameters
      Real                                                              &
     &  timestep                                                        &
     &, GW_kay                                                          &
     &, gwd_frc                                                         &
     &, gwd_fsat

#include "csubmodl.h"
#include "typsts.h"

! Diagnostics info
       Real                                                             &
     & STASHwork(*) ! STASH workspace for section 6 (GW Drag)


! Data arrays

      Integer                                                           &
     &  land_index (land_points)      ! set from land_sea_mask

      Real                                                              &
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)                                               &
     &, rho(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &      model_levels)    ! density*r*r

      Real                                                              &
     &  theta_latest(row_length, rows, model_levels)

      Real                                                              &
     &  p_layer_boundaries(row_length, rows, 0:model_levels)            &
              ! pressure at layer boundaries. Same as p except at
              ! bottom level = pstar, and at top = 0.
     &, sd_orog_land (land_points)                                      &
                                   ! orog/qrparm.orog.stdev
     &, orog_grad_xx_land(land_points)                                  &
                                       ! orog/qrparm.orog.sigmaxx
     &, orog_grad_xy_land(land_points)                                  &
                                       ! orog/qrparm.orog.sigmaxy
     &, orog_grad_yy_land(land_points) ! orog/qrparm.orog.sigmayy

      logical                                                           &
     &  land_sea_mask(row_length, rows)

! Co-ordinate arrays
      Real                                                              &
           ! local vertical co-ordinate information
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:model_levels)           &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, model_levels)              &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

      Real                                                              &
     &  sin_theta_longitude (row_length, rows)                          &
     &, sin_theta_latitude  (row_length, rows)

! time information for current timestep
      Integer                                                           &
     &  timestep_number


! arguments with intent in/out. ie: input variables changed on output.

! arguments with intent out. ie: output variables.

      Real                                                              &
     &  R_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &        model_levels)                                             &
     &, R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,             &
     &        model_levels)

      Integer                                                           &
     &  Error_code

! local variables
      Logical                                                           &
     &  stress_ud_on                                                    &
     &, stress_vd_on                                                    &
     &, stress_ud_satn_on                                               &
     &, stress_vd_satn_on                                               &
     &, stress_ud_wake_on                                               &
     &, stress_vd_wake_on                                               &
     &, du_dt_satn_on                                                   &
     &, dv_dt_satn_on                                                   &
     &, du_dt_wake_on                                                   &
     &, dv_dt_wake_on                                                   &
     &, GWSPEC_EFLUX_ON                                                 &
     &, GWSPEC_SFLUX_ON                                                 &
     &, GWSPEC_WFLUX_ON                                                 &
     &, GWSPEC_NFLUX_ON                                                 &
     &, GWSPEC_EWACC_ON                                                 &
     &, GWSPEC_NSACC_ON                                                 &
     &, u_s_d_on                                                        &
     &, v_s_d_on                                                        &
     &, nsq_s_d_on                                                      &
     &, fr_d_on                                                         &
     &, bld_d_on                                                        &
     &, bldt_d_on                                                       &
     &, num_lim_d_on                                                    &
     &, num_fac_d_on                                                    &
     &, tausx_d_on                                                      &
     &, tausy_d_on                                                      &
     &, taus_scale_d_on                                                 &
     &, L_u_incr_gwd                                                    &
     &, L_v_incr_gwd

      Integer                                                           &
     & i,j,k       ! loop counters


! Local data arrays

! Allocatable arrays for diagnostic variables - required to save memory
! when diagnostic not requested
      Real,Dimension(:,:,:),Allocatable::                               &
     &  u_incr_diagnostic                                               &
                              ! u wind increment for STASH
     &, v_incr_diagnostic                                               &
                              ! v wind increment for STASH
     &, stress_ud                                                       &
     &, stress_vd                                                       &
     &, stress_ud_satn                                                  &
     &, stress_vd_satn                                                  &
     &, stress_ud_wake                                                  &
     &, stress_vd_wake                                                  &
     &, du_dt_satn                                                      &
     &, dv_dt_satn                                                      &
     &, du_dt_wake                                                      &
     &, dv_dt_wake                                                      &
     &, GWSPEC_EFLUX                                                    &
     &, GWSPEC_SFLUX                                                    &
     &, GWSPEC_WFLUX                                                    &
     &, GWSPEC_NFLUX                                                    &
     &, GWSPEC_EWACC                                                    &
     &, GWSPEC_NSACC

      Real,Dimension(:,:),Allocatable::                                 &
     &  u_s_d                                                           &
     &, v_s_d                                                           &
     &, nsq_s_d                                                         &
     &, fr_d                                                            &
     &, bld_d                                                           &
     &, bldt_d                                                          &
     &, num_lim_d                                                       &
     &, num_fac_d                                                       &
     &, tausx_d                                                         &
     &, tausy_d                                                         &
     &, taus_scale_d


! Diagnostic land_point array sizes
      Integer                                                           &
     & points_stress_ud                                                 &
     &,points_stress_vd                                                 &
     &,points_stress_ud_satn                                            &
     &,points_stress_vd_satn                                            &
     &,points_stress_ud_wake                                            &
     &,points_stress_vd_wake                                            &
     &,points_du_dt_satn                                                &
     &,points_dv_dt_satn                                                &
     &,points_du_dt_wake                                                &
     &,points_dv_dt_wake                                                &
     &,points_u_s_d                                                     &
     &,points_v_s_d                                                     &
     &,points_nsq_s_d                                                   &
     &,points_fr_d                                                      &
     &,points_bld_d                                                     &
     &,points_bldt_d                                                    &
     &,points_num_lim_d                                                 &
     &,points_num_fac_d                                                 &
     &,points_tausx_d                                                   &
     &,points_tausy_d                                                   &
     &,points_taus_scale_d

! Local arrays holding information to be passed between physics
! routines.

! Diagnostics controlled by Diagnostic switches


! External Routines:
      External timer, g_wave
#if !defined(SCMA)
      External diagnostics_gwd
#endif

! ----------------------------------------------------------------------
! Section GWD.1 Set stash diagnostic switches
! ----------------------------------------------------------------------

! General case of the atmosphere model, ie : with stash.
#if !defined(SCMA)
      GWSPEC_EFLUX_ON = sf(101,6)
      GWSPEC_SFLUX_ON = sf(102,6)
      GWSPEC_WFLUX_ON = sf(103,6)
      GWSPEC_NFLUX_ON = sf(104,6)
      GWSPEC_EWACC_ON = sf(105,6)
      GWSPEC_NSACC_ON = sf(106,6)
      L_u_incr_gwd      = sf(185,6)
      L_v_incr_gwd      = sf(186,6)
      stress_ud_on      = sf(201,6)
      stress_vd_on      = sf(202,6)
      du_dt_satn_on     = sf(207,6)
      dv_dt_satn_on     = sf(208,6)
      u_s_d_on          = sf(214,6)
      v_s_d_on          = sf(215,6)
      nsq_s_d_on        = sf(216,6)
      fr_d_on           = sf(217,6)
      bld_d_on          = sf(218,6)
      bldt_d_on         = sf(222,6)
      stress_ud_satn_on = sf(223,6)
      stress_vd_satn_on = sf(224,6)
      stress_ud_wake_on = sf(227,6)
      stress_vd_wake_on = sf(228,6)
      du_dt_wake_on     = sf(231,6)
      dv_dt_wake_on     = sf(232,6)
      num_lim_d_on      = sf(233,6)
      num_fac_d_on      = sf(234,6)
      tausx_d_on        = sf(235,6)
      tausy_d_on        = sf(236,6)
      taus_scale_d_on   = sf(237,6)
#else
      GWSPEC_EFLUX_ON = .false.
      GWSPEC_SFLUX_ON = .false.
      GWSPEC_WFLUX_ON = .false.
      GWSPEC_NFLUX_ON = .false.
      GWSPEC_EWACC_ON = .false.
      GWSPEC_NSACC_ON = .false.
      L_u_incr_gwd      = .false.
      L_v_incr_gwd      = .false.
      stress_ud_on      = .false.
      stress_vd_on      = .false.
      stress_ud_satn_on = .false.
      stress_vd_satn_on = .false.
      stress_ud_wake_on = .false.
      stress_vd_wake_on = .false.
      du_dt_satn_on     = .false.
      dv_dt_satn_on     = .false.
      du_dt_wake_on     = .false.
      dv_dt_wake_on     = .false.
      u_s_d_on          = .false.
      v_s_d_on          = .false.
      nsq_s_d_on        = .false.
      fr_d_on           = .false.
      bld_d_on          = .false.
      bldt_d_on         = .false.
      num_lim_d_on      = .false.
      num_fac_d_on      = .false.
      tausx_d_on        = .false.
      tausy_d_on        = .false.
      taus_scale_d_on   = .false.
#endif


! ----------------------------------------------------------------------
! Section GWD.1
! ----------------------------------------------------------------------

      If (error_code  ==  0 ) Then
! DEPENDS ON: timer
        If (Ltimer) Call timer ('GW Drag ',3)
! Save R_u    before updating
      If ( L_u_incr_gwd) Then  ! STASHflag set
        Allocate ( u_incr_diagnostic(row_length,rows,model_levels) )

        Do k=1,model_levels
         Do j=1,rows
          Do i=1,row_length
            u_incr_diagnostic(i,j,k) = R_u(i,j,k)
          Enddo ! i
         Enddo ! j
        Enddo ! k

      Endif                    ! on STASHflag

! Save R_v    before updating
      If ( L_v_incr_gwd) Then  ! STASHflag set
        Allocate ( v_incr_diagnostic(row_length,n_rows,model_levels) )

        Do k=1,model_levels
         Do j=1,n_rows
          Do i=1,row_length
            v_incr_diagnostic(i,j,k) = R_v(i,j,k)
          Enddo ! i
         Enddo ! j
        Enddo ! k

      Endif                    ! on STASHflag

      If ( stress_ud_on ) Then  ! STASHflag set
        Allocate ( stress_ud(row_length,rows,0:model_levels) )
        points_stress_ud = land_points
      Else
        points_stress_ud = 1
      Endif

      If ( stress_vd_on ) Then  ! STASHflag set
        Allocate ( stress_vd(row_length,n_rows,0:model_levels) )
        points_stress_vd = land_points
      Else
        points_stress_vd = 1
      Endif

      If ( du_dt_satn_on ) Then  ! STASHflag set
        Allocate ( du_dt_satn(row_length,rows,model_levels) )
        points_du_dt_satn = land_points
      Else
        points_du_dt_satn = 1
      Endif

      If ( dv_dt_satn_on ) Then  ! STASHflag set
        Allocate ( dv_dt_satn(row_length,n_rows,model_levels) )
        points_dv_dt_satn = land_points
      Else
        points_dv_dt_satn = 1
      Endif

      If ( u_s_d_on ) Then  ! STASHflag set
        Allocate ( u_s_d(row_length,rows) )
        points_u_s_d = land_points
      Else
        points_u_s_d = 1
      Endif

      If ( v_s_d_on ) Then  ! STASHflag set
        Allocate ( v_s_d(row_length,rows) )
        points_v_s_d = land_points
      Else
        points_v_s_d = 1
      Endif

      If ( nsq_s_d_on ) Then  ! STASHflag set
        Allocate ( nsq_s_d(row_length,rows) )
        points_nsq_s_d = land_points
      Else
        points_nsq_s_d = 1
      Endif

      If ( fr_d_on ) Then  ! STASHflag set
        Allocate ( fr_d(row_length,rows) )
        points_fr_d = land_points
      Else
        points_fr_d = 1
      Endif

      If ( bld_d_on ) Then  ! STASHflag set
        Allocate ( bld_d(row_length,rows) )
        points_bld_d = land_points
      Else
        points_bld_d = 1
      Endif

      If ( bldt_d_on ) Then  ! STASHflag set
        Allocate ( bldt_d(row_length,rows) )
        points_bldt_d = land_points
      Else
        points_bldt_d = 1
      Endif

      If ( num_lim_d_on ) Then  ! STASHflag set
        Allocate ( num_lim_d(row_length,rows) )
        points_num_lim_d = land_points
      Else
        points_num_lim_d = 1
      Endif

      If ( num_fac_d_on ) Then  ! STASHflag set
        Allocate ( num_fac_d(row_length,rows) )
        points_num_fac_d = land_points
      Else
        points_num_fac_d = 1
      Endif

      If ( stress_ud_satn_on ) Then  ! STASHflag set
        Allocate ( stress_ud_satn(row_length,rows,0:model_levels) )
        points_stress_ud_satn = land_points
      Else
        points_stress_ud_satn = 1
      Endif

      If ( stress_vd_satn_on ) Then  ! STASHflag set
        Allocate ( stress_vd_satn(row_length,n_rows,0:model_levels) )
        points_stress_vd_satn = land_points
      Else
        points_stress_vd_satn = 1
      Endif

      If ( stress_ud_wake_on ) Then  ! STASHflag set
        Allocate ( stress_ud_wake(row_length,rows,0:model_levels) )
        points_stress_ud_wake = land_points
      Else
        points_stress_ud_wake = 1
      Endif

      If ( stress_vd_wake_on ) Then  ! STASHflag set
        Allocate ( stress_vd_wake(row_length,n_rows,0:model_levels) )
        points_stress_vd_wake = land_points
      Else
        points_stress_vd_wake = 1
      Endif

      If ( du_dt_wake_on ) Then  ! STASHflag set
        Allocate ( du_dt_wake(row_length,rows,model_levels) )
        points_du_dt_wake = land_points
      Else
        points_du_dt_wake = 1
      Endif

      If ( dv_dt_wake_on ) Then  ! STASHflag set
        Allocate ( dv_dt_wake(row_length,n_rows,model_levels) )
        points_dv_dt_wake = land_points
      Else
        points_dv_dt_wake = 1
      Endif

      If ( tausx_d_on ) Then  ! STASHflag set
        Allocate ( tausx_d(row_length,rows) )
        points_tausx_d = land_points
      Else
        points_tausx_d = 1
      Endif

      If ( tausy_d_on ) Then  ! STASHflag set
        Allocate ( tausy_d(row_length,n_rows) )
        points_tausy_d = land_points
      Else
        points_tausy_d = 1
      Endif

      If ( taus_scale_d_on ) Then  ! STASHflag set
        Allocate ( taus_scale_d(row_length,rows) )
        points_taus_scale_d = land_points
      Else
        points_taus_scale_d = 1
      Endif

      IF ( GWSPEC_EFLUX_ON ) Then  ! STASHflag set
        Allocate (GWSPEC_EFLUX(row_length,rows,model_levels))
      ENDIF
      IF ( GWSPEC_SFLUX_ON ) Then  ! STASHflag set
        Allocate (GWSPEC_SFLUX(row_length,n_rows,model_levels))
      ENDIF
      IF ( GWSPEC_WFLUX_ON ) Then  ! STASHflag set
        Allocate (GWSPEC_WFLUX(row_length,rows,model_levels))
      ENDIF
      IF ( GWSPEC_NFLUX_ON ) Then  ! STASHflag set
        Allocate (GWSPEC_NFLUX(row_length,n_rows,model_levels))
      ENDIF
      IF ( GWSPEC_EWACC_ON ) Then  ! STASHflag set
       Allocate (GWSPEC_EWACC(row_length,rows,model_levels))
      ENDIF
      IF ( GWSPEC_NSACC_ON ) Then  ! STASHflag set
       Allocate (GWSPEC_NSACC(row_length,n_rows,model_levels))
      ENDIF


!-------------------------------------------------------------------
!L Section GWD.2  Call Gravity Wave Drag Scheme version 4A
!-------------------------------------------------------------------

       If ( L_gwd ) Then

! DEPENDS ON: g_wave
        CALL G_WAVE(                                                    &
     &   theta_latest, u, v, row_length, rows, n_rows,                  &
     &   off_x, off_y, halo_i, halo_j,                                  &
     &   model_domain, at_extremity, model_levels,                      &
     &   rho, r_rho_levels, r_theta_levels, sd_orog_land,               &
     &   orog_grad_xx_land, orog_grad_xy_land, orog_grad_yy_land,       &
     &   land_index, land_points, timestep, gw_kay, gwd_frc,            &
     &   r_u, r_v, l_taus_scale, l_fix_gwsatn, l_gwd_40km,              &
     &   sat_scheme, gwd_fsat,                                          &
     &   stress_ud      , stress_ud_on     , points_stress_ud     ,     &
     &   stress_vd      , stress_vd_on     , points_stress_vd     ,     &
     &   stress_ud_satn , stress_ud_satn_on, points_stress_ud_satn,     &
     &   stress_vd_satn , stress_vd_satn_on, points_stress_vd_satn,     &
     &   stress_ud_wake , stress_ud_wake_on, points_stress_ud_wake,     &
     &   stress_vd_wake , stress_vd_wake_on, points_stress_vd_wake,     &
     &   du_dt_satn     , du_dt_satn_on    , points_du_dt_satn    ,     &
     &   dv_dt_satn     , dv_dt_satn_on    , points_dv_dt_satn    ,     &
     &   du_dt_wake     , du_dt_wake_on    , points_du_dt_wake    ,     &
     &   dv_dt_wake     , dv_dt_wake_on    , points_dv_dt_wake    ,     &
     &   u_s_d          , u_s_d_on         , points_u_s_d         ,     &
     &   v_s_d          , v_s_d_on         , points_v_s_d         ,     &
     &   nsq_s_d        , nsq_s_d_on       , points_nsq_s_d       ,     &
     &   fr_d           , fr_d_on          , points_fr_d          ,     &
     &   bld_d          , bld_d_on         , points_bld_d         ,     &
     &   bldt_d         , bldt_d_on        , points_bldt_d        ,     &
     &   num_lim_d      , num_lim_d_on     , points_num_lim_d     ,     &
     &   num_fac_d      , num_fac_d_on     , points_num_fac_d     ,     &
     &   tausx_d        , tausx_d_on       , points_tausx_d       ,     &
     &   tausy_d        , tausy_d_on       , points_tausy_d       ,     &
     &   taus_scale_d   , taus_scale_d_on  , points_taus_scale_d  ,     &
     &   error_code)

      End If

      If ( L_use_ussp ) Then

! DEPENDS ON: gw_ussp
      CALL GW_USSP(MODEL_LEVELS, MODEL_DOMAIN, ROWS, N_ROWS,            &
     &    OFF_X, OFF_Y, HALO_I, HALO_J, ROW_LENGTH,                     &
     &    R_RHO_LEVELS, R_THETA_LEVELS, P_LAYER_BOUNDARIES,             &
     &    R_U,R_V, L_USSP_OPAQUE,                                       &
     &    SIN_THETA_LONGITUDE, SIN_THETA_LATITUDE,                      &
     &    THETA_LATEST, RHO, TIMESTEP, U, V, AT_EXTREMITY,              &
     &    GWSPEC_EFLUX,GWSPEC_SFLUX,GWSPEC_WFLUX,GWSPEC_NFLUX,          &
     &    GWSPEC_EWACC,GWSPEC_NSACC,ETA_THETA_LEVELS,                   &
     &    GWSPEC_EFLUX_ON,GWSPEC_SFLUX_ON,GWSPEC_WFLUX_ON,              &
     &    GWSPEC_NFLUX_ON,GWSPEC_EWACC_ON,GWSPEC_NSACC_ON)

      End If


! DEPENDS ON: timer
        If (Ltimer) Call timer ('GW Drag ',4)

! ----------------------------------------------------------------------
! Section GWD.3 Call GWD diagnostics
! ----------------------------------------------------------------------

#if !defined(SCMA)
        If(sf(0,6)) Then ! diagnostics requested this timestep
! DEPENDS ON: timer
        If (Ltimer) Call timer ('Diags   ',3)

! DEPENDS ON: diagnostics_gwd
          Call diagnostics_gwd(                                         &
     &                       row_length, rows, model_levels             &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      at_extremity                               &
     &,                      u, v, R_u, R_v                             &
     &,                      u_incr_diagnostic,v_incr_diagnostic        &
     &,                      stress_ud     ,  stress_vd                 &
     &,                      stress_ud_satn,  stress_vd_satn            &
     &,                      stress_ud_wake,  stress_vd_wake            &
     &,                      du_dt_satn    ,  dv_dt_satn                &
     &,                      du_dt_wake   ,  dv_dt_wake                 &
     &,                      u_s_d, v_s_d, nsq_s_d                      &
     &,                      num_lim_d, num_fac_d                       &
     &,                      fr_d, bld_d, bldt_d                        &
     &, tausx_d, tausy_d, taus_scale_d                                  &
     &, sd_orog_land, land_sea_mask, land_points                        &
     &, orog_grad_xx_land, orog_grad_xy_land, orog_grad_yy_land         &
     &, GWSPEC_EFLUX, GWSPEC_SFLUX, GWSPEC_WFLUX                        &
     &, GWSPEC_NFLUX, GWSPEC_EWACC, GWSPEC_NSACC,                       &
#include "argsts.h"
     & STASHwork                                                        &
     & )
#endif

      If ( L_u_incr_gwd      ) Deallocate ( u_incr_diagnostic )
      If ( L_v_incr_gwd      ) Deallocate ( v_incr_diagnostic )
      If ( stress_ud_on      ) Deallocate ( stress_ud         )
      If ( stress_vd_on      ) Deallocate ( stress_vd         )
      If ( du_dt_satn_on     ) Deallocate ( du_dt_satn        )
      If ( dv_dt_satn_on     ) Deallocate ( dv_dt_satn        )
      If ( u_s_d_on          ) Deallocate ( u_s_d             )
      If ( v_s_d_on          ) Deallocate ( v_s_d             )
      If ( nsq_s_d_on        ) Deallocate ( nsq_s_d           )
      If ( fr_d_on           ) Deallocate ( fr_d              )
      If ( bld_d_on          ) Deallocate ( bld_d             )
      If ( bldt_d_on         ) Deallocate ( bldt_d            )
      If ( num_lim_d_on      ) Deallocate ( num_lim_d         )
      If ( num_fac_d_on      ) Deallocate ( num_fac_d         )
      If ( stress_ud_satn_on ) Deallocate ( stress_ud_satn    )
      If ( stress_vd_satn_on ) Deallocate ( stress_vd_satn    )
      If ( stress_ud_wake_on ) Deallocate ( stress_ud_wake    )
      If ( stress_vd_wake_on ) Deallocate ( stress_vd_wake    )
      If ( du_dt_wake_on     ) Deallocate ( du_dt_wake        )
      If ( dv_dt_wake_on     ) Deallocate ( dv_dt_wake        )
      If ( tausx_d_on        ) Deallocate ( tausx_d           )
      If ( tausy_d_on        ) Deallocate ( tausy_d           )
      If ( taus_scale_d_on   ) Deallocate ( taus_scale_d      )
      IF ( GWSPEC_EFLUX_ON   ) Deallocate ( GWSPEC_EFLUX      )
      IF ( GWSPEC_SFLUX_ON   ) Deallocate ( GWSPEC_SFLUX      )
      IF ( GWSPEC_WFLUX_ON   ) Deallocate ( GWSPEC_WFLUX      )
      IF ( GWSPEC_NFLUX_ON   ) Deallocate ( GWSPEC_NFLUX      )
      IF ( GWSPEC_EWACC_ON   ) Deallocate ( GWSPEC_EWACC      )
      IF ( GWSPEC_NSACC_ON   ) Deallocate ( GWSPEC_NSACC      )

#if !defined(SCMA)
! DEPENDS ON: timer
        If (Ltimer) Call timer ('Diags   ',4)
        Endif            ! on sf(0,6)
#endif

      End If ! on error code equal to zero

! end of routine NI_gwd_ctl
      Return
      END SUBROUTINE NI_gwd_ctl
#endif
#endif
