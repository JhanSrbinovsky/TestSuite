#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Restore convection variables.
!
! Subroutine Interface:
      SUBROUTINE RESTORE_CONV                                           &
! Input data
     &    (resdump, nvars, row_length, rows, model_levels,              &
     &     wet_levels, tr_levels, tr_vars,                              &
     &     OFFX,OFFY,                                                   &
! Output data
     &     theta_conv, q_conv, qcl_conv, qcf_conv,                      &
     &     cf_liquid_conv, cf_frozen_conv, bulk_cf_conv,                &
     &     theta_inc, q_inc, qcl_inc, qcf_inc,                          &
     &     cf_liquid_inc, cf_frozen_inc, bulk_cf_inc,                   &
     &     R_U, R_V, AEROSOL,                                           &
     &     DUST_DIV1,DUST_DIV2,DUST_DIV3,                               &
     &     DUST_DIV4,DUST_DIV5,DUST_DIV6,                               &
     &     SO2, SO4_AITKEN,                                             &
     &     so4_accu, so4_diss, dms, nh3, soot_new, soot_agd,            &
     &     soot_cld, ocff_new, ocff_agd, ocff_cld,                      &
     &     co2, free_tracers, ozone_tracer,                             &
     &     cclwp, conv_rain, conv_snow)

      IMPLICIT NONE

!
! Description: To restore convection variables after the call to
!              convection in the case where convection is called
!              to get the diagnostics out only.
!
! Method:
!
! Owner: UM System team
!
! History:
! Version  Date     Comment
! =======  ====     =======
! 5.3      08/01 Original code. (Zoe Gardner)
!
! 5.4   24/10/02 Passing layer cloud increments in manner that mimics
!                how q_inc is handled. (A.C.Bushell)
!
!    5.5  12/02/03  Include code for mineral dust scheme. S Woodward
!
! Code description:
!   FORTRAN 77 + common extensions also in fortran 90.
!   This code is written to UM programming standards version 7.4

!     INCLUDED COMDECKS

!     Inputs

      INTEGER, Intent(In) ::                                            &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, wet_levels                                                      &
     &, tr_levels                                                       &
     &, tr_vars                                                         &
     &, OFFX                                                            &
     &, OFFY                                                            &
     &, nvars

       REAL, Intent(Out) ::                                             &
     &  theta_conv(row_length, rows, model_levels)                      &
     &, q_conv(row_length, rows, wet_levels)                            &
     &, qcl_conv(row_length, rows, wet_levels)                          &
     &, qcf_conv(row_length, rows, wet_levels)                          &
     &, cf_liquid_conv(row_length, rows, wet_levels)                    &
     &, cf_frozen_conv(row_length, rows, wet_levels)                    &
     &, bulk_cf_conv(row_length, rows, wet_levels)                      &
     &, theta_inc(row_length, rows, model_levels)                       &
     &, q_inc(row_length, rows, wet_levels)                             &
     &, qcl_inc(row_length, rows, wet_levels)                           &
     &, qcf_inc(row_length, rows, wet_levels)                           &
     &, cf_liquid_inc(row_length, rows, wet_levels)                     &
     &, cf_frozen_inc(row_length, rows, wet_levels)                     &
     &, bulk_cf_inc(row_length, rows, wet_levels)                       &
     &, R_u(row_length, rows, model_levels)                             &
     &, R_v(row_length, rows, model_levels)                             &
     &, aerosol(row_length, rows, model_levels)                         &
     &, DUST_DIV1(1-OFFX:ROW_LENGTH+OFFX, 1-OFFY:ROWS+OFFY,             &
     &       MODEL_LEVELS)                                              &
                           !dust mmr in div1
     &, DUST_DIV2(1-OFFX:ROW_LENGTH+OFFX, 1-OFFY:ROWS+OFFY,             &
     &       MODEL_LEVELS)                                              &
                           !dust mmr in div2
     &, DUST_DIV3(1-OFFX:ROW_LENGTH+OFFX, 1-OFFY:ROWS+OFFY,             &
     &       MODEL_LEVELS)                                              &
                           !dust mmr in div3
     &, DUST_DIV4(1-OFFX:ROW_LENGTH+OFFX, 1-OFFY:ROWS+OFFY,             &
     &       MODEL_LEVELS)                                              &
                           !dust mmr in div4
     &, DUST_DIV5(1-OFFX:ROW_LENGTH+OFFX, 1-OFFY:ROWS+OFFY,             &
     &       MODEL_LEVELS)                                              &
                           !dust mmr in div5
     &, DUST_DIV6(1-OFFX:ROW_LENGTH+OFFX, 1-OFFY:ROWS+OFFY,             &
     &       MODEL_LEVELS)                                              &
                           !dust mmr in div6
     &, so2(row_length, rows, model_levels)                             &
     &, so4_aitken(row_length, rows, model_levels)                      &
     &, so4_accu(row_length, rows, model_levels)                        &
     &, so4_diss(row_length, rows, model_levels)                        &
     &, dms(row_length, rows, model_levels)                             &
     &, nh3(row_length, rows, model_levels)                             &
     &, soot_new(row_length, rows, model_levels)                        &
     &, soot_agd(row_length, rows, model_levels)                        &
     &, soot_cld(row_length, rows, model_levels)                        &
     &, ocff_new(row_length, rows, model_levels)                        &
     &, ocff_agd(row_length, rows, model_levels)                        &
     &, ocff_cld(row_length, rows, model_levels)                        &
     &, co2(row_length, rows, model_levels)                             &
     &, free_tracers(row_length, rows, tr_levels, tr_vars)              &
     &, cclwp(row_length, rows)                                         &
     &, conv_rain(row_length, rows)                                     &
     &, conv_snow(row_length, rows)                                     &
     &, ozone_tracer(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &         model_levels)                                        
                           ! Cariolle ozone tracer

      Real, Intent(In) ::                                               &
     &  resdump(row_length, rows, nvars)

! Local Variables

      Integer                                                           &
     &  i,j,k,l, kcount

!----------------------------------------------------------------------

      Do i = 1, row_length
        Do j = 1, rows

          Do k = 1, model_levels
            theta_conv(i,j,k) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + wet_levels - 1
            q_conv(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + wet_levels - 1
            qcl_conv(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + wet_levels - 1
            qcf_conv(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + wet_levels - 1
            cf_liquid_conv(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + wet_levels - 1
            cf_frozen_conv(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + wet_levels - 1
            bulk_cf_conv(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + model_levels - 1
            theta_inc(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + wet_levels - 1
            q_inc(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + wet_levels - 1
            qcl_inc(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + wet_levels - 1
            qcf_inc(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + wet_levels - 1
            cf_liquid_inc(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + wet_levels - 1
            cf_frozen_inc(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + wet_levels - 1
            bulk_cf_inc(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + model_levels - 1
            R_u(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + model_levels - 1
            R_v(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + model_levels - 1
            aerosol(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          KCOUNT = K
          DO K = KCOUNT, KCOUNT + MODEL_LEVELS - 1
            DUST_DIV1(I,J,K-KCOUNT+1) = RESDUMP(I,J,K)
          END DO
          KCOUNT = K
          DO K = KCOUNT, KCOUNT + MODEL_LEVELS - 1
            DUST_DIV2(I,J,K-KCOUNT+1) = RESDUMP(I,J,K)
          END DO
          KCOUNT = K
          DO K = KCOUNT, KCOUNT + MODEL_LEVELS - 1
            DUST_DIV3(I,J,K-KCOUNT+1) = RESDUMP(I,J,K)
          END DO
          KCOUNT = K
          DO K = KCOUNT, KCOUNT + MODEL_LEVELS - 1
            DUST_DIV4(I,J,K-KCOUNT+1) = RESDUMP(I,J,K)
          END DO
          KCOUNT = K
          DO K = KCOUNT, KCOUNT + MODEL_LEVELS - 1
            DUST_DIV5(I,J,K-KCOUNT+1) = RESDUMP(I,J,K)
          END DO
          KCOUNT = K
          DO K = KCOUNT, KCOUNT + MODEL_LEVELS - 1
            DUST_DIV6(I,J,K-KCOUNT+1) = RESDUMP(I,J,K)
          END DO
          kcount = k
          Do k = kcount, kcount + model_levels - 1
            so2(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + model_levels - 1
            so4_aitken(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + model_levels - 1
            so4_accu(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + model_levels - 1
            so4_diss(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + model_levels - 1
            dms(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + model_levels - 1
            nh3(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + model_levels - 1
            soot_new(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + model_levels - 1
            soot_agd(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + model_levels - 1
            soot_cld(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + model_levels - 1
            ocff_new(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + model_levels - 1
            ocff_agd(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + model_levels - 1
            ocff_cld(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + model_levels - 1
            co2(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          Do k = kcount, kcount + model_levels - 1
            ozone_tracer(i,j,k-kcount+1) = resdump(i,j,k)
          End Do
          kcount = k
          cclwp(i,j) = resdump(i,j,k)
          kcount = kcount+1
          conv_snow(i,j) = resdump(i,j,k)
          kcount = kcount+1
          conv_rain(i,j) = resdump(i,j,k)
          kcount = kcount+1
          Do l = 1, tr_vars
            Do k = kcount, kcount + tr_levels - 1
              free_tracers(i,j,k-kcount+1,l) = resdump(i,j,k)
            End Do
          End Do
        End Do
      End Do

      Return
      END SUBROUTINE RESTORE_CONV

#endif
!----------------------------------------------------------------------
