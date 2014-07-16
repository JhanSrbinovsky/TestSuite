#if defined(A12_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine mix_to_q

      Subroutine mix_to_q                                               &
     &                  (row_length, rows, wet_levels,                  &
     &                   halo_i, halo_j,                                &
     &                   mix_v, mix_cl, mix_cf,                         &
     &                   mix_cf2, mix_rain, mix_graup,                  &
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   q, qcl, qcf                                    &
     &                  ,qcf2,qrain,qgraup                              &
     &                   )

! Purpose:
!          Convert from mixing ratios to specific humidities
!
! Method:
!
! Original Progammer: Andrew J. Malcolm
! Current code owner: Andrew J. Malcolm
!
! History:
! Version   Date     Comment
! ----     -------     -------
!  5.4     17/08/01  This deck introduced                Andy Malcolm
!  5.5     25/02/03  Add code for extra moisture variables Andy Malcolm
!  6.2     21/07/05  Fix to not use unset values          Andy Malcolm
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
     &, halo_j

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

! Arguments with Intent OUT. ie: Output

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

! local variables
      Real                                                              &
     & conv                                                             &
     &, moist   (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &           wet_levels)

      Integer                                                           &
     & i, j, k

! ----------------------------------------------------------------------
! Section 1. convert mix_v, mix_cl,mix_cf to q, qcl,qcf
! ----------------------------------------------------------------------

        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
             moist(i,j,k)=1.0+ mix_v (i,j,k)+mix_cl(i,j,k)+mix_cf(i,j,k)
            End Do
          End Do
        End Do

      If (L_mcr_qcf2) Then
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
             moist(i,j,k)= moist(i,j,k)+mix_cf2(i,j,k)
            End Do
          End Do
        End Do
      End If
      If (L_mcr_qrain) Then
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
             moist(i,j,k)= moist(i,j,k)+mix_rain(i,j,k)
            End Do
          End Do
        End Do
      End If
      If (L_mcr_qgraup) Then
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
             moist(i,j,k)= moist(i,j,k)+mix_graup(i,j,k)
            End Do
          End Do
        End Do
      End If

        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              conv= 1./ moist(i,j,k)
              q  (i,j,k) = mix_v (i,j,k)*conv
              qcl(i,j,k) = mix_cl(i,j,k)*conv
              qcf(i,j,k) = mix_cf(i,j,k)*conv
            End Do
          End Do
        End Do

      If (L_mcr_qcf2) Then
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              conv= 1./moist(i,j,k)
              qcf2(i,j,k)  = mix_cf2(i,j,k) * conv
            End Do
          End Do
        End Do
      End If
      If (L_mcr_qrain) Then
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              conv= 1./moist(i,j,k)
              qrain(i,j,k)  = mix_rain(i,j,k) * conv
            End Do
          End Do
        End Do
      End If
      If (L_mcr_qgraup) Then
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              conv= 1./moist(i,j,k)
              qgraup(i,j,k)  = mix_graup(i,j,k) * conv
            End Do
          End Do
        End Do
      End If

! end of routine

      Return
      END SUBROUTINE mix_to_q

#endif
