#if defined(C92_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! SUBROUTINE P_TO_T

      SUBROUTINE P_TO_T (                                               &
                           ! For temperature.
! IN Field dimensions and pointers.
     &  row_length,rows, halo_i, halo_j                                 &
     &, halo_i_data, halo_j_data, LEVELS                                &
! IN Vertical coordinate levels.
     &, R_THETA_LEVELS, R_RHO_LEVELS                                    &
! IN field to be interpolated.
     &, FIELD_IN                                                        &
! OUT Interpolated field.
     &, FIELD_OUT                                                       &
     & )

! PURPOSE:
!
!
! HISTORY:
! DATE   VERSION   COMMENT
! ----   -------   -------
! 03/02/06  6.2   Move to C92_2A. P.Selwood
!
! CODE DESCRIPTION:
!   LANGUAGE: FORTRAN 77 + CRAY EXTENSIONS
!   THIS CODE IS WRITTEN TO UMDP 3 PROGRAMMING STANDARDS.
!

      IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.

      INTEGER                                                           &
     &  row_length,rows                                                 &
     &, halo_i,halo_j                                                   &
     &, halo_i_data,halo_j_data                                         &
     &, LEVELS

      REAL                                                              &
     &  R_THETA_LEVELS (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j,      0:LEVELS+1)          &
     &, R_RHO_LEVELS (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, LEVELS+1)                   &
     &, FIELD_IN (1-halo_i_data:row_length+halo_i_data,                 &
     &            1-halo_j_data:rows+halo_j_data, LEVELS+1)


! ARGUMENTS WITH INTENT IN/OUT. IE: INPUT VARIABLES CHANGED ON OUTPUT.


! ARGUMENTS WITH INTENT OUT. IE: OUTPUT VARIABLES.

      REAL                                                              &
     &  FIELD_OUT (1-halo_i_data:row_length+halo_i_data,                &
     &             1-halo_j_data:rows+halo_j_data, LEVELS)


! LOCAL VARIABLES.

      INTEGER                                                           &
     &  i,j,K

      REAL                                                              &
     &  WEIGHT_1                                                        &
     &, WEIGHT_2                                                        &
     &, WEIGHT_3

! ----------------------------------------------------------------------
! Interpolate field_in (on P grid) to T grid (field_out).
! ----------------------------------------------------------------------

      Do k = 1, levels
        Do j = 1-halo_j_data, rows+halo_j_data
          Do i = 1-halo_i_data, row_length+halo_i_data
            WEIGHT_1 = R_RHO_LEVELS(i,j, K+1) -                         &
     &                 R_RHO_LEVELS(i,j, K)
            WEIGHT_2 = R_THETA_LEVELS(i,j, K) -                         &
     &                 R_RHO_LEVELS(i,j, K)
            WEIGHT_3 = R_RHO_LEVELS(i,j, K+1) -                         &
     &                 R_THETA_LEVELS(i,j, K)
            FIELD_OUT (i,j, K) =                                        &
     &               WEIGHT_3/WEIGHT_1 * FIELD_IN (i,j,K+1)             &
     &             + WEIGHT_2/WEIGHT_1 * FIELD_IN (i,j,K)
          END DO
        End do
      END DO

      RETURN
      END SUBROUTINE P_TO_T
#endif
