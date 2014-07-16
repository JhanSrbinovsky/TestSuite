#if defined(A17_2A) || defined(A17_2B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Convert a proportion of fresh ocff to aged ocff.
!
!
      SUBROUTINE AGEOCFF(                                               &
! Arguments IN
     & row_length, rows, off_x, off_y,                                  &
     & model_levels, timestep,                                          &
! Arguments INOUT
     & ocff_new,                                                        &
! Arguments OUT
     & delta_ageocff                                                    &
     & )
!
! Purpose:
!   To convert a proportion of the fresh ocff to aged ocff. This
!   conversion takes place as an exponential decay with an e-folding
!   time of 1.0 days.
!
!   Called by Aero_ctl
!
! Current owner of code:  Nicolas Bellouin
!
! Code Description:
!  Language: Fortran 90
!  This code is written to UMDP3 v8 programming standards
!
! Documentation: UMDP20
!
!
      IMPLICIT NONE
!
! Arguments with intent IN:
!
      integer ::                                                        &
       row_length,                                                      &
                     !no. of pts along a row
       rows,                                                            &
                     !no. of rows
       off_x,                                                           &
                     !size of small halo in i
       off_y,                                                           &
                     !size of small halo in j.
       model_levels  !no. of model levels
      real :: timestep      !timestep
!
! Arguments with intent IN:
!
      real ::                                                           &
        ocff_new(1-off_x:row_length+off_x,1-off_y:rows+off_y,           &
                model_levels)    !mmr of fresh ocff
!
! Arguments with intent OUT:
!
      real ::                                                           &
        delta_ageocff(row_length, rows, model_levels)
                                 !ocff increment due to ageing
!
! Local variables:
!
      integer :: i,j,k  ! loop counters
      real, parameter :: rate = 1.157E-5 ! conversion rate
                                         ! equivalent to 1/(1.0 day
                                         ! expressed in seconds)

      real :: A ! Local workspace, equal to 1.0-exp(-rate*timestep)
!
!
! Cycle through all points in the field, calculating the amount
! of ocff converted on each point on this timestep.
!
      A = 1.0 - exp(-rate*timestep)
!
      Do k=1,model_levels
        Do j=1,rows
          Do i=1,row_length
            delta_ageocff(i,j,k) = A * ocff_new(i,j,k)
          End Do
        End Do
      End Do

      Return
      END SUBROUTINE AGEOCFF
#endif
