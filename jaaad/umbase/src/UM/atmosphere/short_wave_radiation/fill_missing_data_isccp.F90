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
      Subroutine fill_missing_data_isccp(                               &
     &  off_x, off_y, row_length, rows,                                 &
     &  first_row,last_row,                                             &
     &  first_data_interp, j_lw, ES_space_interp,                       &
     &  L_complete_North, L_complete_South, L_complete_deg)

      Use lw_diag_mod, Only: LW_diag

      Implicit None

!
! VARIABLES WITH INTENT IN
!
      Integer                                                           &
     &  off_x                                                           &
     &, off_y                                                           &
     &, row_length                                                      &
     &, rows                                                            &
     &, model_levels

      Integer                                                           &
     &  first_row                                                       &
     &, last_row                                                        &
     &, first_data_interp                                               &
     &, j_lw

      Real                                                              &
     &   ES_SPACE_INTERP(4, row_length, rows)

      Logical                                                           &
     &  L_complete_North                                                &
     &, L_complete_South                                                &
     &, L_complete_deg                                                  &
     &, l_extra_top

!
!   Isccp diagnostics
!
        If ( LW_diag(j_lw)%L_isccp_weights ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, 1,                      &
     &      LW_diag(j_lw)%isccp_weights                                 &
     &      )
        Endif

        If ( LW_diag(j_lw)%L_isccp_cf ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, 7,                      &
     &      LW_diag(j_lw)%isccp_cf                                      &
     &      )
        Endif

        If ( LW_diag(j_lw)%L_isccp_cf_tau_0_to_p3 ) THEN
! DEPENDS ON: rad3d_inp
           Call rad3d_inp(                                              &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 7,                     &
     &       LW_diag(j_lw)%isccp_cf_tau_0_to_p3                         &
     &       )
        Endif

        If ( LW_diag(j_lw)%L_isccp_cf_tau_p3_to_1p3 ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 7,                     &
     &       LW_diag(j_lw)%isccp_cf_tau_p3_to_1p3                       &
     &       )
        Endif

        If ( LW_diag(j_lw)%L_isccp_cf_tau_1p3_to_3p6 ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 7,                     &
     &       LW_diag(j_lw)%isccp_cf_tau_1p3_to_3p6                      &
     &       )
        Endif

        If ( LW_diag(j_lw)%L_isccp_cf_tau_3p6_to_9p4 ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 7,                     &
     &       LW_diag(j_lw)%isccp_cf_tau_3p6_to_9p4                      &
     &       )
        Endif

        If ( LW_diag(j_lw)%L_isccp_cf_tau_9p4_to_23 ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, 7,                      &
     &      LW_diag(j_lw)%isccp_cf_tau_9p4_to_23                        &
     &      )
        Endif

        If ( LW_diag(j_lw)%L_isccp_cf_tau_23_to_60 ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &       L_complete_North, L_complete_South, L_complete_deg,        &
     &       row_length, rows, off_x, off_y, first_row, last_row,       &
     &       first_data_interp, ES_space_interp, 7,                     &
     &       LW_diag(j_lw)%isccp_cf_tau_23_to_60                        &
     &       )
        Endif

        If ( LW_diag(j_lw)%L_isccp_cf_tau_ge_60 ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, 7,                      &
     &      LW_diag(j_lw)%isccp_cf_tau_ge_60                            &
     &      )
        Endif

      return
      END SUBROUTINE fill_missing_data_isccp
#endif
#endif
