#if defined(ATMOS)
! ARGNORM start

      if( norm_lev_start == norm_lev_end ) then
        do k =1, model_levels
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   k, k,                                          &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
        end do !  k =1, model_levels
      else ! norms from norm_lev_start to norm_lev_end
! DEPENDS ON: print_l2norms
        call print_l2norms(                                             &
     &                   exner_rho_levels,                              &
     &                   u, v, w,                                       &
     &                   u_adv, v_adv, w_adv,                           &
     &                   theta,q,qcl,qcf,                               &
     &                   R_u, R_v, R_w, theta_star,                     &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   row_length, rows, n_rows,                      &
     &                   model_levels, wet_levels,                      &
     &                   norm_lev_start, norm_lev_end,                  &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   mype, nproc_x, nproc_y,                        &
     &                   gc_proc_row_group, gc_proc_col_group,          &
     &                   global_row_length, global_rows,                &
     &                   at_extremity, datastart, model_domain,         &
     &                   .false., .false., L_print_pe )
      endif !  norm_lev_start == norm_lev_end
! DEPENDS ON: um_fort_flush
        if ( L_flush6 ) call UM_FORT_FLUSH(6,info)
! ARGNORM end
#endif
