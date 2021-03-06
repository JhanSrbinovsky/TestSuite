#if defined(A10_2B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine GCR_precon_1_exec

      Subroutine GCR_precon_1_exec(                                     &
     &                     RHS,model_domain,                            &
     &                     row_length, rows, model_levels,              &
     &                     FV_sec_theta_latitude,                       &
     &                     a0, a1, factor, Soln,                        &
     &                     offx, offy,                                  &
     &                     i_start, i_stop, j_start, j_stop             &
     &                     )

! Purpose:
!          Calculates pre-conditioning operator applied to field.
!          Set-up of matrix pre-computed.
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
! Original Progammer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
! Version  Date        Comment
! ----     -------     -------
!   6.2   25/12/05  Pass loop bounds variables as arguments
!                                                     Terry Davies
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit none

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, model_domain

!  Parallel variables

      Integer                                                           &
     &  offx, offy                                                      &
     &, i_start, i_stop                                                 &
                                   ! loop bounds set in PE_Helmholtz
     &, j_start, j_stop            ! loop bounds set in PE_Helmholtz

      Real                                                              &
     &  RHS(row_length,rows,model_levels)                               &
                                          ! RIGHT-HAND-SIDE OF EQUATION.
     &, FV_sec_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy)

      Real                                                              &
     &  a0(row_length,rows,model_levels)                                &
     &, a1(row_length,rows,model_levels)                                &
     &, factor(row_length,rows,model_levels)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.
      Real                                                              &
     &  Soln(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
                ! SOLUTION.

! Local Variables.

      Integer                                                           &
     &  i,j,k

! Local arrays.

#include "parparm.h"
#include "domtyp.h"
!    No External routines.

!-----------------------------------------------------------------------
!     Section 1.
!-----------------------------------------------------------------------

      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            Soln(i,j,k) = RHS(i,j,k) * FV_sec_theta_latitude(i,j)
          End Do
        End Do
      End Do

! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

      Do k= 2, model_levels
        Do j = j_start, j_stop
          Do i = i_start, i_stop
            Soln(i,j,k) = Soln(i,j,k) -                                 &
     &                          factor(i,j,k)*Soln(i,j,k-1)
          End Do
        End Do
      End Do

! Back substitute to get solution.

      Do j = j_start, j_stop
        Do i = i_start, i_stop
          Soln(i,j,model_levels)  = a0(i,j,model_levels) *              &
     &                              Soln(i,j,model_levels)
        End Do
      End Do
      Do k= model_levels-1, 1, -1
        Do j = j_start, j_stop
          Do i = i_start, i_stop
            Soln(i,j,k)  = a0(i,j,k) *                                  &
     &    ( Soln(i,j,k) - a1(i,j,k) * Soln(i,j,k+1) )
          End Do
        End Do
      End Do

!     end of routine GCR_precon_1_exec

      return
      END SUBROUTINE GCR_precon_1_exec

#endif
