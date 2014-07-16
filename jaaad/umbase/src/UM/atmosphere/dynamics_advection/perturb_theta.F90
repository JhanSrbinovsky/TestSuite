#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine perturb_theta
      subroutine perturb_theta(                                         &
     &                         theta, row_length, rows, model_levels,   &
     &                         global_row_length, global_rows,          &
     &                         model_domain, at_extremity,              &
     &                         offx, offy, IntRand_Seed, l_datastart    &
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

! Arguments with Intent IN. ie: Input variables.
      Integer, Intent(In) :: ROW_LENGTH     ! No of points per local row
      Integer, Intent(In) :: ROWS           ! No of local (theta) rows
      Integer, Intent(In) :: MODEL_LEVELS   ! No of model levels
      Integer, Intent(In) :: Offx    ! standard halo size in East-West
      Integer, Intent(In) :: Offy    ! standard halo size in North-South
      Logical, Intent(In) :: At_extremity(4)
      Integer, Intent(In) :: model_domain
      Integer :: IntRand_Seed
      Integer, Intent(In) :: l_datastart(2)
      Integer, Intent(In) :: global_ROW_LENGTH    
      Integer, Intent(In) :: global_ROWs

      Real, Intent (InOut) ::                                           &
     &  theta(1-offx:row_length+offx, 1-offy:rows+offy,                 &
     &        model_levels) 

! Local Variables
      Integer :: i,j,k    ! loop variables
      Integer :: j_start,j_end
      REAL :: RandomNumber(global_row_length,global_rows,model_levels)
      REAL :: eps_machine
      Integer :: gi,gj

#include "domtyp.h"
#include "parparm.h"

!----------------------------------------------------------------------

      if(model_domain==mt_global .and. at_extremity(Psouth))then
        j_start=2
      else 
        j_start=1
      endif
      if(model_domain==mt_global .and. at_extremity(PNorth))then
        j_end=rows-1
      else 
        j_end=rows
      endif

      write(6,*)'PERTURBING THETA FIELD'
      eps_machine=epsilon(1.0)
      write(6,*)' MACHINE epsilon ', eps_machine

! get random number in range 0-1:
! DEPENDS ON: var_randomNumber
      call var_randomNumber(IntRand_seed,RandomNumber,global_row_length,&
     &                      global_rows,model_levels)

      Do k=1,model_levels
        Do j=j_start,j_end
          gj=j+l_datastart(2)-1
          Do i=1, row_length
            gi=i+l_datastart(1)-1
            theta(i,j,k) = theta(i,j,k) +                               &
     &           theta(i,j,k) * (2.0*Randomnumber(gi,gj,k) - 1.0) *     &
     &             eps_machine
          End Do
        End Do
      End Do

      RETURN
      END SUBROUTINE perturb_theta

#endif
