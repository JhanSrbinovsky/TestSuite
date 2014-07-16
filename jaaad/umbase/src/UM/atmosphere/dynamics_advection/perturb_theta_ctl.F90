#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine perturb_theta_ctl
      subroutine perturb_theta_ctl(                                     &
#include "argd1.h"
#include "argptra.h"
     &                         row_length_in, rows_in, model_levels_in, &
     &                         global_row_length_in, global_rows_in,    &
     &                         model_domain, at_extremity,              &
     &                   offx_in, offy_in, IntRand_Seed, l_datastart    &
     &                        )

! Purpose:
!          Perturb theta at the bit level using a random number  
!
! Method:
!          Is described in ;
!
! Original Programmer: Andrew J. Malcolm
! Current code owner: Andrew J. Malcolm
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

      Implicit None
#include "cmaxsize.h"
#include "parparm.h"
#include "typsize.h"
#include "typd1.h"
#include "typptra.h"

! Arguments with Intent IN. ie: Input variables.
      Integer, Intent(In) :: ROW_LENGTH_in     ! No of points per local row
      Integer, Intent(In) :: ROWS_in           ! No of local (theta) rows
      Integer, Intent(In) :: global_ROW_LENGTH_in 
      Integer, Intent(In) :: global_ROWS_in       
      Integer, Intent(In) :: MODEL_LEVELS_in   ! No of model levels
      Integer, Intent(In) :: Offx_in    ! standard halo size in East-West
      Integer, Intent(In) :: Offy_in    ! standard halo size in North-South
      Logical, Intent(In) :: At_extremity(4)
      Integer, Intent(In) :: model_domain
      Integer :: IntRand_Seed
      Integer, Intent(In) :: l_datastart(2)

!----------------------------------------------------------------------
! DEPENDS ON:perturb_theta
      call perturb_theta( d1(jtheta(1)),                                &
     &                         row_length_in, rows_in, model_levels_in, &
     &                         global_row_length_in, global_rows_in,    &
     &                         model_domain, at_extremity,              &
     &                   offx_in, offy_in, IntRand_Seed, l_datastart    &
     &                        )

      RETURN
      END SUBROUTINE perturb_theta_ctl

#endif
