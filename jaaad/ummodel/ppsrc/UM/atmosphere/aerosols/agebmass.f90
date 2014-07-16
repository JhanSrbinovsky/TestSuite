
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Convert a proportion of fresh smoke to aged smoke
!
      SUBROUTINE AGEBMASS(                                              &
! Arguments IN
     & row_length, rows, off_x, off_y,                                  &
     & model_levels, timestep,                                          &
     & bmass_new,                                                       &
! Arguments OUT
     & delta_agebmass                                                   &
     & )
!
! Purpose:
!   To convert a proportion of the fresh smoke to aged smoke. This
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
!   5.5      05/02/03  Original code, based on SOOT2A.     P. Davison
!   6.2      20/10/05  DEF for a17_2b added.               A. Jones
!   6.2      01/11/05  Ageing timescale decreased to 6 hours if version
!                      2B of the aerosol scheme is used.
!                                            A. Jones & J. M. Haywood
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
     &  bmass_new(1-off_x:row_length+off_x,1-off_y:rows+off_y,          &
     &            model_levels)    !mmr of fresh smoke
!
! Arguments with intent OUT:
!
      REAL                                                              &
     &  delta_agebmass(row_length, rows, model_levels)
                                 !smoke increment due to ageing
!
! Local variables:
!
      INTEGER                                                           &
     & i,j,k  ! loop counters
      REAL                                                              &
     & rate   ! conversion rate
!





      parameter(rate=4.6296E-5) ! equivalent to 1/(6.0 hours expressed
                                ! in seconds)

      REAL                                                              &
     & A      ! Local workspace, equal to 1.0-exp(-rate*timestep)
!
!
! Cycle through all points in the field, calculating the amount
! of smoke converted on each point on this timestep.
!
      A = 1.0 - exp(-rate*timestep)
!
      Do k=1,model_levels
        Do j=1,rows
          Do i=1,row_length
            delta_agebmass(i,j,k) = A * bmass_new(i,j,k)
          End Do
        End Do
      End Do

      Return
      END SUBROUTINE AGEBMASS
!
!-----------------------------------------------------------------------
!
!+ Perform nucleation scavenging of aged smoke to cloud smoke.
!
