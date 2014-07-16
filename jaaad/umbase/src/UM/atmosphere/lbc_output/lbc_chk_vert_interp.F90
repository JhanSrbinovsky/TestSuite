#if defined(A32_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Check if Vertical Interpolation is required for LBC data
!
! Subroutine Interface:

      Subroutine LBC_Chk_Vert_Interp (                                  &
#include "arginfa.h"
#include "argduma.h"
#include "argptra.h"
     & idummy )

      IMPLICIT NONE
!
! Description:
!   Determine whether vertical interpolation of Lateral Boundary
!   conditions is required between model and LBC levels.
!
! Method:
!   Vertical interpolation is switched on if any of the following are
!   different.
!   - No of model or wet levels.
!   - First rho level at which height is constant.
!   - Height of top of model.
!   - Eta values for theta or rho levels.
!
! Current Code Owner: Dave Robinson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2    13/11/00   Original code. Dave Robinson
!   6.1    18/08/04   Do checking for VERTLEVS namelist only.
!                     D Robinson.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Arguments

#include "parvars.h"
#include "typsize.h"
#include "typinfa.h"
#include "typduma.h"
#include "typptra.h"

      Integer idummy

#include "cmaxsize.h"
#include "cintfa.h"
#include "cprintst.h"

! Local variables

      Character(Len=*), Parameter :: RoutineName= 'LBC_Chk_Vert_Interp'

      Integer :: jlev
      Integer :: jintf


      Do jintf = 1, n_intf_a

        If (lbc_nd(jintf) == 1) Then

! ------------------------
! Check no of model levels
! ------------------------

        If ( model_levels /= intf_p_levels(jintf) ) Then
          Intf_Vert_Interp(jintf) = .true.
          If (PrintStatus >= PrStatus_Oper ) Then
            Write (6,*) 'LBC_CHK_VI : Difference in no of model levels'
          End If
        End If

! ----------------------
! Check no of wet levels
! ----------------------

        If ( wet_levels /= intf_q_levels(jintf) ) Then
          Intf_Vert_Interp(jintf) = .true.
          If (PrintStatus >= PrStatus_Oper ) Then
            Write (6,*) 'LBC_CHK_VI : Difference in no of wet levels'
          End If
        End If

! -------------------------------------------------
! Check first rho level at which height is constant
! -------------------------------------------------

        If ( a_inthd(24) /= lbc_first_r_rho(jintf) ) Then
          Intf_Vert_Interp(jintf) = .true.
          If (PrintStatus >= PrStatus_Oper ) Then
            Write (6,*) 'LBC_CHK_VI : Difference in first constant ',   &
     &                  'rho level.'
          End If
        End If

! ----------------------------
! Check height of top of model
! ----------------------------

        If ( abs ( a_realhd(16) - lbc_z_top_model(jintf) ) >            &
     &       Epsilon (1.0) ) Then
          Intf_Vert_Interp(jintf) = .true.
          If (PrintStatus >= PrStatus_Oper ) Then
            Write (6,*) 'LBC_CHK_VI : Difference in hts of top level'
          End If
        End If

! -----------------------------------------
! Check eta values for theta and rho levels
! -----------------------------------------

        If ( .not. Intf_Vert_Interp(jintf) ) Then

          Do jlev = 1, model_levels+1
            If ( abs ( a_levdepc(jetatheta+jlev-1) -                    &
     &           lbc_eta_theta(jlev,jintf) ) > Epsilon (1.0) ) Then
              Intf_Vert_Interp(jintf) = .true.
              If (PrintStatus >= PrStatus_Oper ) Then
                Write (6,*) 'LBC_CHK_VI : Difference in eta_theta'
              End If
              Exit
            End If
          End Do

          Do jlev = 1, model_levels
            If ( abs ( a_levdepc(jetarho+jlev-1) -                      &
     &           lbc_eta_rho(jlev,jintf) ) > Epsilon (1.0) ) Then
              Intf_Vert_Interp(jintf) = .true.
              If (PrintStatus >= PrStatus_Oper ) Then
                Write (6,*) 'LBC_CHK_VI : Difference in eta_rho'
              End If
              Exit
            End If
          End Do

        End If

        If (PrintStatus >= PrStatus_Normal ) Then
          If ( Intf_Vert_Interp(jintf) ) Then
            write (6,*) 'LBC_CHK_VI : Vertical Interpolation required ',&
     &                  'for LBC Area ',lbc_stream_a(jintf)
          End If
        End If

        End If

      End Do

      Return
      END SUBROUTINE LBC_Chk_Vert_Interp
#endif
