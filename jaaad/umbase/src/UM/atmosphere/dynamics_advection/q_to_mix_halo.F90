#if defined(A12_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine q_to_mix


! Subroutine q_to_mix_halo

      Subroutine q_to_mix_halo                                          &
     &                  (row_length, rows, wet_levels,                  &
     &                   halo_i, halo_j,                                &
     &                   q, qcl, qcf,                                   &
     &                   qcf2, qrain, qgraup,                           &
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   mix_v, mix_cl, mix_cf                          &
     &                  ,mix_cf2, mix_rain, mix_graup                   &
     &                   )

! Purpose:
!          Convert from specific humidities to mixing ratios
!
! Method:
!
! Original Progammer: Andrew J. Malcolm
! Current code owner: Andrew J. Malcolm
!
! History:
! Version   Date     Comment
! ----     -------     -------
!  6.2     21/07/05  This subroutine added                 Andy Malcolm
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

! local variables
      Real                                                              &
     & conv                                                             &
     &, moist(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &             wet_levels)

      Integer                                                           &
     & i, j, k

! ----------------------------------------------------------------------
! Section 1. convert q, qcl,qcf to mix_v, mix_cl,mix_cf
! ----------------------------------------------------------------------
      do k=1,wet_levels
        do j=1-halo_j,rows+halo_j
          do i=1-halo_i,row_length+halo_i
            moist(i,j,k) = 1 - q(i,j,k) - qcl(i,j,k) - qcf(i,j,k)
          end do
        end do
      end do
      If (L_mcr_qcf2) Then
        do k=1,wet_levels
        do j=1-halo_j,rows+halo_j
          do i=1-halo_i,row_length+halo_i
              moist(i,j,k)      = moist(i,j,k) - qcf2(i,j,k)
            end do
          end do
        end do
      End If
      If (L_mcr_qrain) Then
        do k=1,wet_levels
        do j=1-halo_j,rows+halo_j
          do i=1-halo_i,row_length+halo_i
              moist(i,j,k)      = moist(i,j,k) - qrain(i,j,k)
            end do
          end do
        end do
      End If
      If (L_mcr_qgraup) Then
        do k=1,wet_levels
        do j=1-halo_j,rows+halo_j
          do i=1-halo_i,row_length+halo_i
              moist(i,j,k)      = moist(i,j,k) - qgraup(i,j,k)
            end do
          end do
        end do
      End If


        Do k = 1, wet_levels
        do j=1-halo_j,rows+halo_j
          do i=1-halo_i,row_length+halo_i
              conv= 1./moist(i,j,k)
              mix_v (i,j,k) = q  (i,j,k)*conv
              mix_cl(i,j,k) = qcl(i,j,k)*conv
              mix_cf(i,j,k) = qcf(i,j,k)*conv
            End Do
          End Do
        End Do

      If (L_mcr_qcf2) Then
        Do k = 1, wet_levels
        do j=1-halo_j,rows+halo_j
          do i=1-halo_i,row_length+halo_i
              conv= 1./moist(i,j,k)
              mix_cf2(i,j,k) = qcf2(i,j,k) * conv
            End Do
          End Do
        End Do
      End If
      If (L_mcr_qrain) Then
        Do k = 1, wet_levels
        do j=1-halo_j,rows+halo_j
          do i=1-halo_i,row_length+halo_i
              conv= 1./moist(i,j,k)
              mix_rain(i,j,k) = qrain(i,j,k) * conv
            End Do
          End Do
        End Do
      End If
      If (L_mcr_qgraup) Then
        Do k = 1, wet_levels
        do j=1-halo_j,rows+halo_j
          do i=1-halo_i,row_length+halo_i
              conv= 1./moist(i,j,k)
              mix_graup(i,j,k) = qgraup(i,j,k) * conv
            End Do
          End Do
        End Do
      End If

! end of routine

      Return
      END SUBROUTINE q_to_mix_halo

#endif
