#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine: RESETATM ----------------------------------------------
!LL
!LL  Purpose: Control routine to perform recalculations of prognostic
!LL           quantities following dump.  For reproducibility, various
!LL           fields need to be recalculated using rounded values.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Author:   T.C.Johns
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.1    8/02/93 : added comdeck CHSUNITS to define NUNITS for
!LL                    comdeck CCONTROL.
!LL   3.2    27/03/93 Dynamic allocation of main data arrays. R. Rawlins
!LL   5.1    27/06/00 Subroutine reintroduced & adapted for the ND code
!LL                   (It had been removed at 5.0). JC Thil
!     5.3    27/07/01 Add calc_Exner_at_theta and remove calc_Exner
!                     and calc_P_at_theta              Clive Wilson
!     6.0    21/11/03 Do pressure calculations throughout halos (as
!                     in SETCONA). P.Selwood.
!     6.1    18/08/04 Remove repeated declaration of CMESSAGE. P.Dando
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: C0
!LL
!LL  Project task: C0
!LL
!LL  External documentation: On-line UM document C0 - The top-level
!LL                          control system
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE RESETATM(                                              &
#include "argd1.h"
#include "argptra.h"
#include "argcona.h"
     &                  ICODE,CMESSAGE)

      Use level_heights_mod, Only : r_theta_levels, r_rho_levels

      IMPLICIT NONE

#include "csubmodl.h"
#include "cmaxsize.h"
#include "parvars.h"

#include "typsize.h"
#include "typd1.h"
#include "typptra.h"
#include "typcona.h"
!L       Addresses of arrays within super arrays

!*----------------------------------------------------------------------
!  Common blocks
#include "chsunits.h"
#include "ccontrol.h"
#include "cphyscon.h"
! Local scalars:
      REAL                                                              &
     & constant

      Integer                                                           &
     &  i, j, k, l                                                      &
                          ! loop counters
     &, ij                ! 2D array index for offx/y variables
!  Subroutines called
      EXTERNAL Calc_P_star, Calc_Exner_at_theta, Swap_Bounds,           &
     &         Calc_P_from_Exner
!L
! Error reporting
      INTEGER       ICODE       ! =0 normal exit; >0 error exit
      CHARACTER*256 Cmessage    ! Error message
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='RESETATM')

      ICODE=0
      CMESSAGE=''
!
!L----------------------------------------------------------------------
!L
       constant = 1./ kappa
       Do k = 1, model_levels+1
         Do ij = 1, (rows + 2*offy) * (row_length + 2* offx)
           d1(jp(k) +ij-1) =                                            &
     &       (d1(jexner_rho_levels(k) +ij-1)**constant) * p_zero
         End Do
       End Do

! DEPENDS ON: calc_exner_at_theta
        Call Calc_Exner_at_theta(                                       &
     &  r_theta_levels, r_rho_levels,                                   &
     & d1(jexner_rho_levels(1)),                                        &
     &  row_length, rows, model_levels,                                 &
     &  offx, offy, halo_i, halo_j,                                     &
     &  d1(jexner_theta_levels(1)), .TRUE.)

! Calculate pressure from Exner at theta levels.
! DEPENDS ON: calc_p_from_exner
       call Calc_P_from_Exner(                                          &
     &                      d1(jp_theta_levels(1)), kappa, p_zero,      &
     &                      row_length, rows, model_levels,             &
     &                      offx, offy,                                 &
     &                      d1(jexner_theta_levels(1)),.TRUE.)

! DEPENDS ON: calc_p_star
       Call Calc_P_star(                                                &
     &  r_theta_levels, r_rho_levels,                                   &
     &  d1(jp(1)), d1(jrho(1)),                                         &
     &  g, row_length, rows, model_levels,                              &
     &  offx, offy, halo_i, halo_j,                                     &
     &  d1(jpstar))


      END SUBROUTINE RESETATM
#endif
