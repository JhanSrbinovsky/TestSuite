#if defined(A70_1C) || defined(A70_1Z)
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to fill in missing data
!
! Method:
!
!   Radiative fluxes may not have been calculated at all
!   points: we now fill in as required. This part of the
!   code was originally located in RAD_CTL2 (v6.1 and below)
!   but has been move into a subroutine in order to make
!   RAD_CTL2 more readable.

! Current Code Owner:  J.-C. Thelen
!
! History
!
! Version   Data        Comment
! 6.2       19/01/06    Original code
!                       (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77/90  with extensions listed in documentation.
!
!-----------------------------------------------------------------------
!
      Subroutine fill_missing_data_sw(                                  &
     &  off_x, off_y, row_length, rows,model_levels,                    &
     &  salt_dim1, salt_dim2, salt_dim3,                                &
     &  cloud_levels,first_row,last_row,                                &
     &  first_data_interp, ES_space_interp,                             &
     &  L_complete_North, L_complete_South, L_complete_deg,             &
     &  L_flux_below_690nm_surf, n_channel, j_sw, l_extra_top,          &
     &  SW_incs, netsw, SWsea,flux_below_690nm_surf,                    &
     &  top_absorption, surf_down_SW, dirpar_inc,                       &
     &  sea_salt_film, sea_salt_jet)

      Use solinc_data, Only:                                            &
        f_orog, L_orog
      Use sw_diag_mod, Only:                                            &
        SW_diag

      Implicit None

!
! VARIABLES WITH INTENT IN
!
      Integer                                                           &
     &  off_x                                                           &
     &, off_y                                                           &
     &, row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, cloud_levels                                                    &
     &, salt_dim1                                                       &
     &, salt_dim2                                                       &
     &, salt_dim3

      Integer                                                           &
     &  first_row                                                       &
     &, last_row                                                        &
     &, first_data_interp                                               &
     &, n_channel                                                       &
     &, j_sw

      Real                                                              &
     &   ES_SPACE_INTERP(4, row_length, rows)

      Logical                                                           &
     &  L_complete_North                                                &
     &, L_complete_South                                                &
     &, L_complete_deg                                                  &
     &, L_flux_below_690nm_surf                                         &
     &, l_extra_top
!
! VARIABLES WITH INTENT IN/OUT
!
      Real                                                              &
     &  SW_incs(row_length, rows, 0:model_levels+1)                     &
     &, netSW(row_length, rows)                                         &
     &, SWsea(row_length, rows)                                         &
     &, surf_down_sw(row_length,rows,4)                                 &
     &, top_absorption(row_length, rows)                                &
     &, flux_below_690nm_surf(row_length,rows)                          &
     &, sea_salt_film(salt_dim1, salt_dim2, salt_dim3)                  &
     &, sea_salt_jet(salt_dim1, salt_dim2, salt_dim3)                   &
     &, dirpar_inc(row_length, rows)

!
!           Primary Fields:
!
! DEPENDS ON: rad3d_inp
        Call rad3d_inp(                                                 &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp,                        &
     &       model_levels+2,                                            &
     &       SW_incs                                                    &
     &       )
! DEPENDS ON: rad3d_inp
        Call rad3d_inp(                                                 &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       netSW                                                      &
     &       )
! DEPENDS ON: rad3d_inp
        Call rad3d_inp(                                                 &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       SWsea                                                      &
     &       )
!
! The next field is purely diagnostic in some versions of
! the model, but is always required with MOSES.
!
        If ( L_flux_below_690nm_surf ) Then

! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       flux_below_690nm_surf                                      &
     &       )
        Endif

        If (SW_diag(j_sw)%L_direct_par) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, 1,                  &
     &          dirpar_inc                                              &
     &          )
        Endif

        If (l_extra_top) then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       top_absorption                                             &
     &       )
        Endif

! DEPENDS ON: rad3d_inp
        Call rad3d_inp(                                                 &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp,                        &
     &       4,                                                         &
     &       surf_down_SW                                               &
     &       )

        If (L_orog) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, 1,                    &
     &        f_orog                                                    &
     &        )
        Endif
!
!           Complete the diagnostic fields as required.
!
        If ( SW_diag(j_sw)%L_solar_out_toa ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       SW_diag(j_sw)%solar_out_toa                                &
     &       )
        Endif
        If ( SW_diag(j_sw)%L_solar_out_clear ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       SW_diag(j_sw)%solar_out_clear                              &
     &       )
        Endif
        If ( SW_diag(j_sw)%L_surface_down_flux ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       SW_diag(j_sw)%surface_down_flux                            &
     &       )
        Endif
        If ( SW_diag(j_sw)%L_surf_down_clr ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       SW_diag(j_sw)%surf_down_clr                                &
     &       )
        Endif
        If ( SW_diag(j_sw)%L_surf_up_clr ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       SW_diag(j_sw)%surf_up_clr                                  &
     &       )
        Endif
        If ( SW_diag(j_sw)%L_clear_hr ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp,                        &
     &       model_levels,                                              &
     &       SW_diag(j_sw)%clear_hr                                     &
     &       )
        Endif
        If ( SW_diag(j_sw)%L_net_flux_trop ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       SW_diag(j_sw)%net_flux_trop                                &
     &       )
        Endif
        If ( SW_diag(j_sw)%L_up_flux_trop ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       SW_diag(j_sw)%up_flux_trop                                 &
     &       )
        Endif
!
!           Microphysical diagnostics
!
        If ( SW_diag(j_sw)%re_conv_flag ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp,                        &
     &       cloud_levels,                                              &
     &       SW_diag(j_sw)%re_conv                                      &
     &       )
        Endif
        If ( SW_diag(j_sw)%re_strat_flag ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp,                        &
     &       cloud_levels,                                              &
     &       SW_diag(j_sw)%re_strat                                     &
     &       )
        Endif
        If ( SW_diag(j_sw)%wgt_conv_flag ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp,                        &
     &       cloud_levels,                                              &
     &       SW_diag(j_sw)%wgt_conv                                     &
     &       )
        Endif
        If ( SW_diag(j_sw)%wgt_strat_flag ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp,                        &
     &       cloud_levels,                                              &
     &       SW_diag(j_sw)%wgt_strat                                    &
     &       )
        Endif
        If ( SW_diag(j_sw)%lwp_strat_flag ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp,                        &
     &       cloud_levels,                                              &
     &       SW_diag(j_sw)%lwp_strat                                    &
     &       )
        Endif
        If ( SW_diag(j_sw)%weighted_re_flag ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       SW_diag(j_sw)%weighted_re                                  &
     &       )
        Endif
        If ( SW_diag(j_sw)%sum_weight_re_flag ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       SW_diag(j_sw)%sum_weight_re                                &
     &       )
        Endif
        If ( SW_diag(j_sw)%wgtd_warm_re_flag ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       SW_diag(j_sw)%weighted_warm_re                             &
     &       )
        Endif
        If ( SW_diag(j_sw)%sum_wgt_warm_re_flag ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       SW_diag(j_sw)%sum_weight_warm_re                           &
     &       )
        Endif
        If ( SW_diag(j_sw)%Nc_diag_flag ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       SW_diag(j_sw)%Nc_diag                                      &
     &       )
        Endif
        If ( SW_diag(j_sw)%Nc_weight_flag ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       SW_diag(j_sw)%Nc_weight                                    &
     &       )
        Endif
        If ( SW_diag(j_sw)%ntot_diag_flag ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp,                        &
     &       cloud_levels,                                              &
     &       SW_diag(j_sw)%ntot_diag                                    &
     &       )
        Endif
        If ( SW_diag(j_sw)%strat_lwc_diag_flag ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp,                        &
     &       cloud_levels,                                              &
     &       SW_diag(j_sw)%strat_lwc_diag                               &
     &       )
        Endif
        If ( SW_diag(j_sw)%so4_ccn_diag_flag ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp,                        &
     &       cloud_levels,                                              &
     &       SW_diag(j_sw)%so4_ccn_diag                                 &
     &       )
        Endif
        If ( SW_diag(j_sw)%cond_samp_wgt_flag ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp,                        &
     &       cloud_levels,                                              &
     &       SW_diag(j_sw)%cond_samp_wgt                                &
     &       )
        Endif
        If ( SW_diag(j_sw)%seasalt_film_flag ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp,                        &
     &       salt_dim3,                                                 &
     &       sea_salt_film                                              &
     &       )
        Endif
        If ( SW_diag(j_sw)%seasalt_jet_flag ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp,                        &
     &       salt_dim3,                                                 &
     &       sea_salt_jet                                               &
     &       )
        Endif
!
! Diagnostics for MOSES
!
        If ( SW_diag(j_sw)%L_FlxSolBelow690nmSurf ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       SW_diag(j_sw)%FlxSolBelow690nmSurf                         &
     &       )
        Endif
        If ( SW_diag(j_sw)%L_FlxSeaBelow690nmSurf ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       SW_diag(j_sw)%FlxSeaBelow690nmSurf                         &
     &       )
        Endif
!
! Diagnostic for Radiance
!
        If ( SW_diag(j_sw)%L_toa_radiance ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, n_channel,          &
     &          SW_diag(j_sw)%toa_radiance                              &
     &          )
        Endif
!
! Diagnostic for Fluxes
!
        If ( SW_diag(j_sw)%L_uvflux_direct ) Then
! DEPENDS ON: rad3d_inp
           Call rad3d_inp(                                              &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp,                       &
     &        model_levels+1,                                           &
     &        SW_diag(j_sw)%uvflux_direct                               &
     &       )
        Endif
        If ( SW_diag(j_sw)%L_uvflux_up ) Then
! DEPENDS ON: rad3d_inp
           Call rad3d_inp(                                              &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp,                       &
     &        model_levels+1,                                           &
     &        SW_diag(j_sw)%uvflux_up                                   &
     &       )
        Endif
        If ( SW_diag(j_sw)%L_uvflux_net ) Then
! DEPENDS ON: rad3d_inp
           Call rad3d_inp(                                              &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp,                       &
     &        model_levels+1,                                           &
     &        SW_diag(j_sw)%uvflux_net                                  &
     &       )
        Endif
        If ( SW_diag(j_sw)%L_flux_direct ) Then
! DEPENDS ON: rad3d_inp
           Call rad3d_inp(                                              &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp,                       &
     &        model_levels+1,                                           &
     &        SW_diag(j_sw)%flux_direct                                 &
     &       )
        Endif
        If ( SW_diag(j_sw)%L_flux_diffuse ) Then
! DEPENDS ON: rad3d_inp
           Call rad3d_inp(                                              &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp,                       &
     &        model_levels+1,                                           &
     &        SW_diag(j_sw)%flux_diffuse                                &
     &       )
        Endif
!
! Diagnostics for orography correction
!
        If ( SW_diag(j_sw)%L_orog_corr ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 1,                     &
     &       SW_diag(j_sw)%orog_corr                                    &
     &       )
        Endif
!
! Extinction diagnostics:
!
        If ( SW_diag(j_sw)%L_cloud_extinction ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, cloud_levels,          &
     &       SW_diag(j_sw)%cloud_extinction                             &
     &       )
        Endif
        If ( SW_diag(j_sw)%L_cloud_weight_extinction ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, cloud_levels,          &
     &       SW_diag(j_sw)%cloud_weight_extinction                      &
     &       )
        Endif
        If ( SW_diag(j_sw)%L_ls_cloud_extinction ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, cloud_levels,          &
     &       SW_diag(j_sw)%ls_cloud_extinction                          &
     &       )
        Endif
        If ( SW_diag(j_sw)%L_ls_cloud_weight_extinction ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, cloud_levels,          &
     &       SW_diag(j_sw)%ls_cloud_weight_extinction                   &
     &       )
        Endif
        If ( SW_diag(j_sw)%L_cnv_cloud_extinction ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, cloud_levels,          &
     &       SW_diag(j_sw)%cnv_cloud_extinction                         &
     &       )
        Endif
        If ( SW_diag(j_sw)%L_cnv_cloud_weight_extinction ) Then
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, cloud_levels,          &
     &       SW_diag(j_sw)%cnv_cloud_weight_extinction                  &
     &       )
        Endif

      return
      END SUBROUTINE fill_missing_data_sw
#endif
#endif
