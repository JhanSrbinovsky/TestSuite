
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Convert a proportion of fresh soot to aged soot
!
      SUBROUTINE AGESOOT(                                               &
! Arguments IN
     & row_length, rows, off_x, off_y,                                  &
     & model_levels, timestep,                                          &
! Arguments INOUT
     & soot_new,                                                        &
! Arguments OUT
     & delta_agesoot                                                    &
     & )
!
! Purpose:
!   To convert a proportion of the fresh soot to aged soot. This
!   conversion takes place as an exponential decay with an e-folding
!   time of 1.0 days.
!
!   Called by Aero_ctl
!
! Current owners of code: P. Davison
!
! History:
! Version     Date     Comment
! -------     ----     -------
!
!   5.4      10/07/02  Original code, based on a re-write of
!                                     v4.5 deck SOOT1A    P. Davison
!   6.2      20/10/05  DEF for a17_2b added.              A. Jones
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
! Documentation: UMDP20
!
!
      IMPLICIT NONE
!
! Arguments with intent IN:
!
      INTEGER                                                           &
     & row_length,                                                      &
                     !no. of pts along a row
     & rows,                                                            &
                     !no. of rows
     & off_x,                                                           &
                     !size of small halo in i
     & off_y,                                                           &
                     !size of small halo in j.
     & model_levels  !no. of model levels
      REAL                                                              &
     & timestep      !timestep
!
! Arguments with intent IN:
!
      REAL                                                              &
     &  soot_new(1-off_x:row_length+off_x,1-off_y:rows+off_y,           &
     &          model_levels)    !mmr of fresh soot
!
! Arguments with intent OUT:
!
      REAL                                                              &
     &  delta_agesoot(row_length, rows, model_levels)
                                 !soot increment due to ageing
!
! Local variables:
!
      INTEGER                                                           &
     & i,j,k  ! loop counters
      REAL                                                              &
     & rate   ! conversion rate
!
      parameter(rate=1.157E-5) ! equivalent to 1/(1.0 days expressed
                               ! in seconds)
      REAL                                                              &
     & A      ! Local workspace, equal to 1.0-exp(-rate*timestep)
!
!
! Cycle through all points in the field, calculating the amount
! of soot converted on each point on this timestep.
!
      A = 1.0 - exp(-rate*timestep)
!
      Do k=1,model_levels
        Do j=1,rows
          Do i=1,row_length
            delta_agesoot(i,j,k) = A * soot_new(i,j,k)
          End Do
        End Do
      End Do

      Return
      END SUBROUTINE AGESOOT
!
!-----------------------------------------------------------------------
!
!+ Perform diffusional scavenging of aged soot to cloud soot
!
