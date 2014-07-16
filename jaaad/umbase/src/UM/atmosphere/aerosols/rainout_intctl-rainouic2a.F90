#if defined(A17_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
      subroutine rainout_intctl(                                        &
     &               row_length, rows,                                  &
     &               off_x, off_y, halo_i, halo_j,                      &
     &               model_levels, wet_model_levels,                    &
     &               r_rho_levels, r_theta_levels,                      &
     &               rho_r2, Q,                                         &
     &               QCF_REMAIN, QCL_REMAIN,                            &
     &               QCF_PREVIOUS, QCL_PREVIOUS,                        &
     &               LS_RAIN3D, LS_SNOW3D,                              &
     &               timestep,                                          &
     &               AERO_INCLOUD, AERO_ACCUM,                          &
     &               RNOUT_AERO)
!
!---------------------------------------------------------------------
! Purpose: Wrapper to version 2A of the aerosol rain-out routine.
!
!          Called by microphys_ctl (deck mcr_ctl2)
!
! Current code owners: N Bellouin
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.2   20/01/06   New deck                               N Bellouin
!
! Code description:
!   Language: FORTRAN 77 + common extensions
!
! Documentation: UMDP 20
!
!---------------------------------------------------------------------
!
      IMPLICIT NONE
!
! Arguments with intent in:
!
      INTEGER                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, off_x                                                           &
                                        !EW size of std. halo
     &, off_y                                                           &
                                        !NS size of std. halo
     &, halo_i                                                          &
                                        !EW extended halo
     &, halo_j                                                          &
                                        !NS extended halo
     &, model_levels                                                    &
     &, wet_model_levels
!
      REAL                                                              &
     &  r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j,     model_levels)            &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j, 0:model_levels)            &
     &, rho_r2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &                                      model_levels)  !density*r*r
!
      REAL                                                              &
     &  QCL_REMAIN(row_length,rows,wet_model_levels)                    &
     &, QCF_REMAIN(row_length,rows,wet_model_levels)                    &
     &, QCF_PREVIOUS(row_length,rows,wet_model_levels)                  &
     &, QCL_PREVIOUS(row_length,rows,wet_model_levels)                  &
     &, Q(row_length,rows,wet_model_levels)                             &
     &, LS_RAIN3D(row_length,rows,wet_model_levels)                     &
     &, LS_SNOW3D(row_length,rows,wet_model_levels)
!
      REAL timestep ! timestep in seconds
!
! Arguments with intent IN/OUT:
!   mass mixing ratio of in-cloud and accumulation/aged aerosol
      REAL                                                              &
     & AERO_INCLOUD(1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &                                         model_levels),           &
     & AERO_ACCUM(1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &                                         model_levels)

!
! Arguments with intent OUT (diagnostics):
      REAL                                                              &
     & RNOUT_AERO(row_length,rows)      !tracer removed kg/m2/ts
!
! No local variables
!

!
! External routines
!

      EXTERNAL RAINOUT

! DEPENDS ON: rainout
      CALL RAINOUT(                                                     &
     &             row_length, rows,                                    &
     &             off_x, off_y, halo_i, halo_j,                        &
     &             model_levels, wet_model_levels,                      &
     &             r_rho_levels, r_theta_levels,                        &
     &             rho_r2, Q,                                           &
     &             QCF_REMAIN, QCL_REMAIN,                              &
     &             QCF_PREVIOUS, QCL_PREVIOUS,                          &
     &             LS_RAIN3D, LS_SNOW3D,                                &
     &             AERO_INCLOUD,                                        &
     &             RNOUT_AERO)

      RETURN
      END SUBROUTINE rainout_intctl

#endif
