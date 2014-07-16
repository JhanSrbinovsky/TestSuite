#if defined(A12_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine calc_mix_star

      Subroutine calc_mix_star                                          &
     &                  (row_length, rows, wet_levels,                  &
     &                   halo_i, halo_j, offx, offy,                    &
     &                   q, qcl, qcf,                                   &
     &                   qcf2, qrain, qgraup,                           &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   qcf2_star, qrain_star, qgraup_star,            &
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   mix_v_star, mix_cl_star, mix_cf_star           &
     &                  ,mix_cf2_star, mix_rain_star, mix_graup_star    &
     &                   )

! Purpose:
!          calculate mixing ratio increments
!
! Method:
!
! Original Progammer: Andrew J. Malcolm
! Current code owner: Andrew J. Malcolm
!
! History:
! Version   Date     Comment
! ----     -------     -------
!  5.4     17/08/01  This deck introduced           Andy Malcolm
!  5.5     25/02/03  Add code for extra moisture variables Andy Malcolm
!  6.2     05/04/05  Fix code to not use unset values      Andy Malcolm
!
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

      Integer                                                           &
     &  row_length                                                      &
                        ! number of points on a row
     &, rows                                                            &
                        ! number of rows of data
     &, wet_levels                                                      &
     &, halo_i                                                          &
     &, halo_j                                                          &
     &, offx                                                            &
     &, offy

      Logical                                                           &
     &  L_mcr_qcf2                                                      &
                          ! true if second cloud ice active
     &, L_mcr_qrain                                                     &
                          ! true if rain active
     &, L_mcr_qgraup      ! true if graupel active

! Arguments with Intent IN. ie: Input

      Real                                                              &
     &  q    (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        wet_levels)                                               &
     &, qcl  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        wet_levels)                                               &
     &, qcf  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        wet_levels)                                               &
     &, qcf2    (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &           wet_levels)                                            &
     &, qrain   (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &           wet_levels)                                            &
     &, qgraup  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &           wet_levels)

      Real                                                              &
     &  q_star    (1-offx:row_length+offx, 1-offy:rows+offy,            &
     &        wet_levels)                                               &
     &, qcl_star  (1-offx:row_length+offx, 1-offy:rows+offy,            &
     &        wet_levels)                                               &
     &, qcf_star  (1-offx:row_length+offx, 1-offy:rows+offy,            &
     &        wet_levels)                                               &
     &, qcf2_star    (1-offx:row_length+offx, 1-offy:rows+offy,         &
     &                wet_levels)                                       &
     &, qrain_star   (1-offx:row_length+offx, 1-offy:rows+offy,         &
     &                wet_levels)                                       &
     &, qgraup_star  (1-offx:row_length+offx, 1-offy:rows+offy,         &
     &                wet_levels)

! Arguments with Intent OUT. ie: Output

      Real                                                              &
     &  mix_v_star  (1-offx:row_length+offx, 1-offy:rows+offy,          &
     &          wet_levels)                                             &
     &, mix_cl_star (1-offx:row_length+offx, 1-offy:rows+offy,          &
     &          wet_levels)                                             &
     &, mix_cf_star (1-offx:row_length+offx, 1-offy:rows+offy,          &
     &          wet_levels)                                             &
     &, mix_cf2_star   (1-offx:row_length+offx, 1-offy:rows+offy,       &
     &                  wet_levels)                                     &
     &, mix_rain_star  (1-offx:row_length+offx, 1-offy:rows+offy,       &
     &                  wet_levels)                                     &
     &, mix_graup_star (1-offx:row_length+offx, 1-offy:rows+offy,       &
     &                  wet_levels)

! local variables
      Real                                                              &
     & conv,sum_q,sum_q_star                                            &
     &, moist(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &             wet_levels)                                          &
     &, moist_star(1-offx:row_length+offx, 1-offy:rows+offy,            &
     &             wet_levels)

      Integer                                                           &
     & i, j, k

! ----------------------------------------------------------------------
! Section 1. calculate mix_v_star, mix_cl_star, mix_cf_star
! ----------------------------------------------------------------------

      do k=1,wet_levels
        do j=1,rows
          do i=1,row_length
            moist(i,j,k) = 1 - q(i,j,k) - qcl(i,j,k) - qcf(i,j,k)
            moist_star(i,j,k) = q_star(i,j,k) + qcl_star(i,j,k) +       &
     &                          qcf_star(i,j,k)
          end do
        end do
      end do
      If (L_mcr_qcf2) Then
        do k=1,wet_levels
          do j=1,rows
            do i=1,row_length
              moist(i,j,k)      = moist(i,j,k) - qcf2(i,j,k)
              moist_star(i,j,k) = moist_star(i,j,k) + qcf2_star(i,j,k)
            end do
          end do
        end do
      End If
      If (L_mcr_qrain) Then
        do k=1,wet_levels
          do j=1,rows
            do i=1,row_length
              moist(i,j,k)      = moist(i,j,k) - qrain(i,j,k)
              moist_star(i,j,k) = moist_star(i,j,k) + qrain_star(i,j,k)
            end do
          end do
        end do
      End If
      If (L_mcr_qgraup) Then
        do k=1,wet_levels
          do j=1,rows
            do i=1,row_length
             moist(i,j,k)      = moist(i,j,k) - qgraup(i,j,k)
              moist_star(i,j,k) = moist_star(i,j,k) + qgraup_star(i,j,k)
            end do
          end do
        end do
      End If


        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              conv= 1./(  moist(i,j,k) *                                &
     &                   (moist(i,j,k) - moist_star(i,j,k)) )
              mix_v_star(i,j,k)  =                                      &
     &               ( q_star(i,j,k)*moist(i,j,k) +                     &
     &                 q(i,j,k)*moist_star(i,j,k) )  *conv
              mix_cl_star(i,j,k) =                                      &
     &               ( qcl_star(i,j,k)*moist(i,j,k) +                   &
     &                 qcl(i,j,k)*moist_star(i,j,k) )  *conv
              mix_cf_star(i,j,k) =                                      &
     &               ( qcf_star(i,j,k)*moist(i,j,k) +                   &
     &                 qcf(i,j,k)*moist_star(i,j,k) )  *conv
            End Do
          End Do
        End Do

      If (L_mcr_qcf2) Then
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              conv= 1./(  moist(i,j,k) *                                &
     &                   (moist(i,j,k) - moist_star(i,j,k)) )
              mix_cf2_star(i,j,k)  =                                    &
     &               ( qcf2_star(i,j,k)*moist(i,j,k) +                  &
     &                 qcf2(i,j,k)*moist_star(i,j,k) )  *conv
            End Do
          End Do
        End Do
      End If
      If (L_mcr_qrain) Then
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              conv= 1./(  moist(i,j,k) *                                &
     &                   (moist(i,j,k) - moist_star(i,j,k)) )
              mix_rain_star(i,j,k)  =                                   &
     &               ( qrain_star(i,j,k)*moist(i,j,k) +                 &
     &                 qrain(i,j,k)*moist_star(i,j,k) )  *conv
            End Do
          End Do
        End Do
      End If
      If (L_mcr_qgraup) Then
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              conv= 1./(  moist(i,j,k) *                                &
     &                   (moist(i,j,k) - moist_star(i,j,k)) )
              mix_graup_star(i,j,k)  =                                  &
     &               ( qgraup_star(i,j,k)*moist(i,j,k) +                &
     &                 qgraup(i,j,k)*moist_star(i,j,k) )  *conv
            End Do
          End Do
        End Do
      End If

! end of routine

      Return
      END SUBROUTINE calc_mix_star

#endif
