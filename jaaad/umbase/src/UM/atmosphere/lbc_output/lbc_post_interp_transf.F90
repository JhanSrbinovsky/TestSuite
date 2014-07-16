#if defined(A32_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Performs transforms after a LBC field has been interpolated
!
! Subroutine Interface:

      Subroutine LBC_Post_Interp_Transf (                               &
     &           lbc_size                                               &
     &,          lbc_first_level                                        &
     &,          lbc_last_level                                         &
     &,          lbc_stashcode                                          &
     &,          lbc_data                                               &
     & )

      IMPLICIT NONE
!
! Description:
!   Perform tranformations/processing on a LBC field after it has
!   been interpolated.
!
! Method:
!   Choice of transform/processing is based on stashcode. Data passes
!   through unchanged if no processing done.
!
! Current Code Owner: Dave Robinson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2    13/11/00   Original code. Dave Robinson
!   5.3    22/10/01   Remove write statements. Dave Robinson
!   5.5    03/02/03   Include check on qcf2,qrain,qgraup. R.M.Forbes
!   6.0    29/07/03   Include bounds check for cloud fraction lbcs
!                                              Damian Wilson
!   6.2    01/10/04   Reset murk if < 0.1  R.M.Forbes
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!

! Arguments
      Integer  ::  lbc_size
      Integer  ::  lbc_first_level
      Integer  ::  lbc_last_level
      Integer  ::  lbc_stashcode

      Real  ::  lbc_data (lbc_size, lbc_first_level:lbc_last_level)

#include "cmaxsize.h"
#include "cintfa.h"
#include "parvars.h"
#include "parlbcs.h"

! Local
      Integer :: Level
      Integer :: i

! -------------------------
! Do transforms as required
! -------------------------

      Select Case ( lbc_stashcode )

        Case (lbc_stashcode_w, lbc_stashcode_w_adv)

          lbc_data(:,lbc_first_level) = 0.0
          lbc_data(:,lbc_last_level)  = 0.0

        Case (lbc_stashcode_q)

          Do level = lbc_first_level, lbc_last_level
            Do i = 1, lbc_size
              If ( lbc_data(i,level) < lbc_q_min ) Then
                lbc_data(i,level) = lbc_q_min
              End If
            End Do
          End Do

        Case (lbc_stashcode_qcf, lbc_stashcode_qcl,                     &
     &        lbc_stashcode_qcf2, lbc_stashcode_qrain,                  &
     &        lbc_stashcode_qgraup)

          Do level = lbc_first_level, lbc_last_level
            Do i = 1, lbc_size
              If ( lbc_data(i,level) < 0.0 ) Then
                lbc_data(i,level) = 0.0
              End If
            End Do
          End Do

        Case (lbc_stashcode_cf_bulk, lbc_stashcode_cf_liquid,           &
     &        lbc_stashcode_cf_frozen)

          Do level = lbc_first_level, lbc_last_level
            Do i = 1, lbc_size
              If (lbc_data(i,level) < 0.0 ) Then
                lbc_data(i,level) = 0.0
              End If
              If (lbc_data(i,level) > 1.0 ) Then
                lbc_data(i,level) = 1.0
              End If
            End Do
          End Do

        Case (lbc_stashcode_murk)

        ! Reset murk aerosol if less than 0.1 after interpolation

          Do level = lbc_first_level, lbc_last_level
            Do i = 1, lbc_size
              If ( lbc_data(i,level) < 0.1 ) Then
                lbc_data(i,level) = 0.1
              End If
            End Do
          End Do

      End Select

      Return
      END SUBROUTINE LBC_Post_Interp_Transf
#endif
