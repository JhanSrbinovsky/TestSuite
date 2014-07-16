#if defined(A15_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculate PV at all theta points.

      SUBROUTINE Calc_PV_at_theta ( u, v, theta, rho,                   &
                                                                  ! in
     &                              r_theta_levels, r_rho_levels,       &
                                                                  ! in
     &                              r_at_u, r_at_v,                     &
                                                                  ! in
     &                              sec_v_latitude,                     &
                                                                  ! in
     &                              tan_v_latitude,                     &
                                                                  ! in
     &                              sec_theta_latitude,                 &
                                                                  ! in
     &                              f3_at_v,                            &
                                                                  ! in
     &                              delta_lambda, delta_phi,            &
                                                                  ! in
     &                              Model_domain,                       &
                                                                  ! in
     &                              pv_at_theta )                 ! out

! Description:
!
!   Calculate PV at all theta points.
!
! Method:
!
!   1. Call Calc_PV to obtain PV midway E-W between v points on rho
!      levels.
!   2. Add haloes to resulting field.
!   3. Interpolate horizontally and vertically to theta points.
!      (PV at top theta level is set to PV at top rho level.)
!   4. Reset polar rows to their mean values.
!
! Current Code Owner: Adam Clayton.
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 5.2       14/12/00 Original code. Adam Clayton
! 5.3       19/10/01 Use appropriate gcg routines.   S. Cusack
! 6.0       10/11/03 Fix South Pole means on 1x1 PEs. P.Selwood.
! 6.2       21/10/05 Replace Gsync with Ssync. P.Selwood.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! Declarations:

      IMPLICIT NONE

! Common blocks:

#include "parvars.h"
#include "typsize.h"
#include "domtyp.h"

! Subroutine arguments:

      REAL, INTENT(IN) ::                                               &
      
     &  u                  ( 1 - offx   : row_length + offx,            &
     &                       1 - offy   : rows       + offy,            &
     &                       1          : model_levels ),               &
      
     &  v                  ( 1 - offx   : row_length + offx,            &
     &                       1 - offy   : n_rows     + offy,            &
     &                       1          : model_levels ),               &
      
     &  theta              ( 1 - offx   : row_length + offx,            &
     &                       1 - offy   : rows       + offy,            &
     &                       1          : model_levels ),               &
      
     &  rho                ( 1 - offx   : row_length + offx,            &
     &                       1 - offy   : rows       + offy,            &
     &                       1          : model_levels ),               &
      
     &  r_theta_levels     ( 1 - halo_i : row_length + halo_i,          &
     &                       1 - halo_j : rows       + halo_j,          &
     &                       0          : model_levels ),               &
      
     &  r_rho_levels       ( 1 - halo_i : row_length + halo_i,          &
     &                       1 - halo_j : rows       + halo_j,          &
     &                       1          : model_levels ),               &
      
     &  r_at_u             ( 1 - halo_i : row_length + halo_i,          &
     &                       1 - halo_j : rows       + halo_j,          &
     &                       1          : model_levels ),               &
      
     &  r_at_v             ( 1 - halo_i : row_length + halo_i,          &
     &                       1 - halo_j : n_rows     + halo_j,          &
     &                       1          : model_levels ),               &
      
     &  sec_v_latitude     ( 1 - offx   : row_length + offx,            &
     &                       1 - offy   : n_rows     + offy ),          &
      
     &  tan_v_latitude     ( row_length : n_rows ),                     &
      
     &  sec_theta_latitude ( 1 - offx   : row_length + offx,            &
     &                       1 - offy   : rows       + offy ),          &
      
     &  f3_at_v            ( 1 - offx   : row_length + offx,            &
     &                       1 - offy   : n_rows     + offy ),          &
      
     &  delta_lambda,                                                   &
     &  delta_phi

      INTEGER, INTENT(IN) ::                                            &
      
     &  Model_domain

      REAL, INTENT(OUT) ::                                              &
      
     &  pv_at_theta (row_length, rows, model_levels)

! Local variables:

      INTEGER :: i, j, k
      INTEGER :: ICode

      REAL :: pv (row_length, n_rows, model_levels)

      REAL :: pv_plus_haloes ( 1 - offx   : row_length + offx,          &
     &                         1 - offy   : n_rows     + offy,          &
     &                         1          : model_levels )

      REAL :: polar_sums (model_levels)
      REAL :: polar_means(model_levels)

! External subroutines called:

      EXTERNAL Calc_PV, Swap_Bounds

!- End of header ------------------------------------------------------

!----------------------------------------------------------------------
! [1]: Calculate PV midway E-W between v points on rho levels.
!----------------------------------------------------------------------

! DEPENDS ON: calc_pv
      CALL Calc_PV  ( u, v, theta, rho,                                 &
                                                      ! in
     &                r_theta_levels, r_rho_levels,                     &
                                                      ! in
     &                r_at_u, r_at_v,                                   &
                                                      ! in
     &                sec_v_latitude, tan_v_latitude,                   &
                                                      ! in
     &                sec_theta_latitude, f3_at_v,                      &
                                                      ! in
     &                delta_lambda, delta_phi,                          &
                                                      ! in
     &                row_length, rows, n_rows,                         &
                                                      ! in
     &                model_levels,                                     &
                                                      ! in
     &                offx, offy, halo_i, halo_j,                       &
                                                      ! in
     &                at_extremity,                                     &
                                                      ! in
     &                pv )                            ! out

!----------------------------------------------------------------------
! [2]: Add haloes.
!----------------------------------------------------------------------

      pv_plus_haloes(:,:,:) = 0.0

      DO k = 1, model_levels
        DO j = 1, n_rows
          DO i = 1, row_length
            pv_plus_haloes(i,j,k) = pv(i,j,k)
          END DO
        END DO
      END DO

! DEPENDS ON: swap_bounds
      CALL Swap_Bounds ( pv_plus_haloes,                                &
                                                             ! inout
     &                   row_length, n_rows, model_levels,              &
                                                             ! in
     &                   offx, offy, fld_type_v, .FALSE. )   ! in

! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(pv_plus_haloes,row_length, n_rows,       &
     &                     model_levels,offx,offy)

!----------------------------------------------------------------------
! [3]: Interpolate to theta points.
!----------------------------------------------------------------------

      DO k = 1, model_levels - 1
        DO j = 1, rows
          DO i = 1, row_length
            pv_at_theta(i,j,k) = 0.25                                   &
     &                           * ( ( pv_plus_haloes(i,  j,  k)        &
     &                               + pv_plus_haloes(i,  j-1,k)        &
     &                               + pv_plus_haloes(i-1,j,  k)        &
     &                               + pv_plus_haloes(i-1,j-1,k)        &
     &                               )                                  &
     &                             * ( r_rho_levels  (i,j,k+1)          &
     &                               - r_theta_levels(i,j,k)            &
     &                               )                                  &
     &                             + ( pv_plus_haloes(i,  j,  k+1)      &
     &                               + pv_plus_haloes(i,  j-1,k+1)      &
     &                               + pv_plus_haloes(i-1,j,  k+1)      &
     &                               + pv_plus_haloes(i-1,j-1,k+1)      &
     &                               )                                  &
     &                             * ( r_theta_levels(i,j,k)            &
     &                               - r_rho_levels  (i,j,k)            &
     &                             ) )                                  &
     &                           / ( r_rho_levels(i,j,k+1)              &
     &                             - r_rho_levels(i,j,k)                &
     &                             )
          END DO
        END DO
      END DO

      ! Set PV at top theta level equal to PV at top rho level.
      k = model_levels
      DO j = 1, rows
        DO i = 1, row_length
          pv_at_theta(i,j,k) = 0.25                                     &
     &                       * ( pv_plus_haloes(i,  j,  k)              &
     &                         + pv_plus_haloes(i,  j-1,k)              &
     &                         + pv_plus_haloes(i-1,j,  k)              &
     &                         + pv_plus_haloes(i-1,j-1,k)              &
     &                         )
        END DO
      END DO

!----------------------------------------------------------------------
! [4]: Reset polar rows to their mean values.
!----------------------------------------------------------------------

      IF (Model_domain == mt_global) THEN

        CALL GC_SSYNC ( nproc, ICode )

        IF (at_extremity(PNorth)) THEN

#if defined(REPROD)
          CALL GCG_RVECSUMR (row_length,                                &
                                                    ! in
     &                       row_length,                                &
                                                    ! in
     &                       1,                                         &
                                                    ! in
     &                       model_levels,                              &
                                                    ! in
     &                       pv_at_theta(:,rows,:),                     &
                                                    ! in
     &                       gc_proc_row_group,                         &
                                                    ! in
     &                       ICode,                                     &
                                                    ! out
     &                       polar_sums)            ! out
#else
          CALL GCG_RVECSUMF (row_length,                                &
                                                    ! in
     &                       row_length,                                &
                                                    ! in
     &                       1,                                         &
                                                    ! in
     &                       model_levels,                              &
                                                    ! in
     &                       pv_at_theta(:,rows,:),                     &
                                                    ! in
     &                       gc_proc_row_group,                         &
                                                    ! in
     &                       ICode,                                     &
                                                    ! out
     &                       polar_sums)            ! out
#endif

           polar_means(:) = polar_sums(:) / REAL(global_row_length)

           DO k = 1, model_levels
             pv_at_theta(:,rows,k) = polar_means(k)
           END DO

        END IF   ! PNorth

        IF (at_extremity(PSouth)) THEN

#if defined(REPROD)
          CALL GCG_RVECSUMR (row_length,                                &
                                                 ! in
     &                       row_length,                                &
                                                 ! in
     &                       1,                                         &
                                                 ! in
     &                       model_levels,                              &
                                                 ! in
     &                       pv_at_theta(:,1,:),                        &
                                                 ! in
     &                       gc_proc_row_group,                         &
                                                 ! in
     &                       ICode,                                     &
                                                 ! out
     &                       polar_sums)         ! out
#else
          CALL GCG_RVECSUMF (row_length,                                &
                                                 ! in
     &                       row_length,                                &
                                                 ! in
     &                       1,                                         &
                                                 ! in
     &                       model_levels,                              &
                                                 ! in
     &                       pv_at_theta(:,1,:),                        &
                                                 ! in
     &                       gc_proc_row_group,                         &
                                                 ! in
     &                       ICode,                                     &
                                                 ! out
     &                       polar_sums)         ! out
#endif

           polar_means(:) = polar_sums(:) / REAL(global_row_length)

           DO k = 1, model_levels
             pv_at_theta(:,1,k) = polar_means(k)
           END DO

        END IF

      END IF ! (Model_domain == mt_global)


      RETURN
      END SUBROUTINE Calc_PV_at_theta
#endif
