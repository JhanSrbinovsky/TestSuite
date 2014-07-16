#if defined(A12_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine calc_q_star

      Subroutine calc_q_star                                            &
     &                  (row_length, rows, wet_levels,                  &
     &                   halo_i, halo_j, offx, offy,                    &
     &                   mix_v, mix_cl, mix_cf,                         &
     &                   mix_cf2, mix_rain, mix_graup,                  &
     &                   mix_v_star, mix_cl_star, mix_cf_star,          &
     &                   mix_cf2_star, mix_rain_star, mix_graup_star,   &
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   q, qcl, qcf,                                   &
     &                   qcf2, qrain, qgraup,                           &
     &                   q_star, qcl_star, qcf_star                     &
     &                  ,qcf2_star, qrain_star, qgraup_star             &
     &                   )

! Purpose:
!          calculate specific humidity increments
!
! Method:
!
! Original Progammer: Andrew J. Malcolm
! Current code owner: Andrew J. Malcolm
!
! History:
! Version   Date     Comment
! ----     -------     -------
!  5.4     17/08/01  This deck introduced                 Andy Malcolm
!  5.5     25/02/03  Add code for extra moisture variables Andy Malcolm
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
     &  mix_v  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          wet_levels)                                             &
     &, mix_cl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          wet_levels)                                             &
     &, mix_cf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          wet_levels)                                             &
     &, mix_cf2   (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,    &
     &             wet_levels)                                          &
     &, mix_rain  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,    &
     &             wet_levels)                                          &
     &, mix_graup (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,    &
     &             wet_levels)

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

! Arguments with Intent OUT. ie: Output

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

! local variables
      Real                                                              &
     & conv,sum_mix,sum_mix_star                                        &
     &, moist(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &             wet_levels)                                          &
     &, moist_star(1-offx:row_length+offx, 1-offy:rows+offy,            &
     &             wet_levels)

      Integer                                                           &
     & i, j, k

! ----------------------------------------------------------------------
! Section 1. calculate q_star, qcl_star, qcf_star
! ----------------------------------------------------------------------
      moist = 1 + mix_v + mix_cl + mix_cf
      moist_star = mix_v_star + mix_cl_star + mix_cf_star
      If (L_mcr_qcf2) Then
        moist      = moist + mix_cf2
        moist_star = moist_star + mix_cf2_star
      End If
      If (L_mcr_qrain) Then
        moist      = moist + mix_rain
        moist_star = moist_star + mix_rain_star
      End If
      If (L_mcr_qgraup) Then
        moist      = moist + mix_graup
        moist_star = moist_star + mix_graup_star
      End If

        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              conv= 1./(  moist(i,j,k) *                                &
     &                   (moist(i,j,k) + moist_star(i,j,k)) )
              q_star(i,j,k)  =                                          &
     &               ( mix_v_star(i,j,k)*moist(i,j,k) -                 &
     &                 mix_v(i,j,k)*moist_star(i,j,k) )  * conv
              qcl_star(i,j,k)  =                                        &
     &               ( mix_cl_star(i,j,k)*moist(i,j,k) -                &
     &                 mix_cl(i,j,k)*moist_star(i,j,k) )  * conv
              qcf_star(i,j,k)  =                                        &
     &               ( mix_cf_star(i,j,k)*moist(i,j,k) -                &
     &                 mix_cf(i,j,k)*moist_star(i,j,k) )  * conv

            End Do
          End Do
        End Do

      If (L_mcr_qcf2) Then
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              conv= 1./(  moist(i,j,k) *                                &
     &                   (moist(i,j,k) + moist_star(i,j,k)) )
              qcf2_star(i,j,k)  =                                       &
     &               ( mix_cf2_star(i,j,k)*moist(i,j,k) -               &
     &                 mix_cf2(i,j,k)*moist_star(i,j,k) )  * conv
            End Do
          End Do
        End Do
      End If
      If (L_mcr_qrain) Then
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              conv= 1./(  moist(i,j,k) *                                &
     &                   (moist(i,j,k) + moist_star(i,j,k)) )
              qrain_star(i,j,k)  =                                      &
     &               ( mix_rain_star(i,j,k)*moist(i,j,k) -              &
     &                 mix_rain(i,j,k)*moist_star(i,j,k) ) * conv
            End Do
          End Do
        End Do
      End If
      If (L_mcr_qgraup) Then
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              conv= 1./(  moist(i,j,k) *                                &
     &                   (moist(i,j,k) + moist_star(i,j,k)) )
              qgraup_star(i,j,k)  =                                     &
     &               ( mix_graup_star(i,j,k)*moist(i,j,k) -             &
     &                 mix_graup(i,j,k)*moist_star(i,j,k) ) * conv
            End Do
          End Do
        End Do
      End If

! end of routine

      Return
      END SUBROUTINE calc_q_star

#endif
