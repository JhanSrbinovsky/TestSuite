#if defined(A14_1B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine: INIT_EMCORR--------------------------------------------
!LL
!LL  Purpose: Interface routine required to pass super arrays down into
!LL           ENG_MASS_CAL which initialises the energy correction.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 6.1.5A
!LL
!LL  Model            Modification history:
!LL version  date
!LL  3.2   30/03/93  Redefine INIT_EMCORR to become a new routine
!LL                  INIT_EMCORR2 and introduce a control interface
!LL                  routine INIT_EMCORR for dynamic allocation of
!LL                  main data arrays.
!LL  3.4   23/06/94  Arguments LLINTS,LWHITBROM added and passed to
!LL                           INIT_EMCORR2        S.J.Swarbrick
!LL  4.1   23/04/96  Added TYPFLDPT variables to pass down to
!                    ENG_MASS_DIAG                 P.Burton
!LL  Model            Modification history:
!LL version  Date
!LL  5.1  27/03/00  Split include file argcona.h into argcona.h &
!LL                 arglndm.h JC Thil
!LL  5.1   20/01/00  Correct for new dynamics R.A.Stratton
!LL  5.3   08/06/01  Remove duplicate declarations.  A van der Wal
!  6.2   25/12/05  Variable resolution changes            Yongming Tang
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: C0
!LL
!LL  Project task: C0
!LL
!LL  External documentation: On-line UM document C1 - The top-level
!LL                          dynamic allocation
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE INIT_EMCORR(                                           &
#include "argd1.h"
#include "argduma.h"
#include "argptra.h"
#include "argcona.h"
#include "arglndm.h"
     &              ICODE,CMESSAGE)
!
!*----------------------------------------------------------------------

use level_heights_mod
use dyn_var_res_mod
use trignometric_mod
      IMPLICIT NONE
!
!  Subroutines called
!
      EXTERNAL ENG_MASS_DIAG
!
!  Arguments
!
!  Configuration-dependent sizes and arrays

#include "cmaxsize.h"
#include "csubmodl.h"
#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
#include "typduma.h"
#include "typptra.h"
#include "typcona.h"
#include "typlndm.h"
#include "chsunits.h"
#include "ccontrol.h"
#include "c_mdi.h"
#include "c_pi.h"
#include "ctime.h"

      INTEGER ICODE             ! Work - Internal return code
      CHARACTER*256 CMESSAGE    ! Work - Internal error message
!  NOTE icode and cmessage are not altered by this subroutine
!       returned unchanged.
!-------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------
      Logical                                                           &
     &  lmass_cor0                                                      &
                        ! switch for mass correction at T+0
     &, lmoist_cor0     ! switch for moist correction at T+0

      Real                                                              &
     &  dummy           ! dummy initial mass
      INTEGER i,j,k,ij
!-------------------------------------------------------------------
! Only set initial energy and mass if the header values in the start
! dump are set to missing data otherwise the model will start from
! the values in the dump.
! Values in the dump will come from a previous climate integration.

      IF(A_REALHD(rh_tot_mass_init) == RMDI .OR.                        &
     &   A_REALHD(rh_tot_energy_init) == RMDI .OR.                      &
     &   A_REALHD(rh_energy_corr) == RMDI) THEN
!

        A_REALHD(rh_tot_mass_init) = 0.0
        A_REALHD(rh_tot_energy_init) = 0.0

        lmass_cor0=.false.        ! do not attempt mass correction
        lmoist_cor0=.false.       ! do not attempt moist correction
        dummy = 0.0
        delta_lambda = a_realhd(rh_deltaEW)*PI_over_180
        delta_phi    = a_realhd(rh_deltaNS)*PI_over_180

! DEPENDS ON: eng_mass_diag
          Call eng_mass_diag (                                          &
     &                      halo_i, halo_j, offx, offy                  &
     &,                     global_row_length, gc_proc_row_group        &
     &,                     gc_proc_col_group                           &
     &,                     at_extremity, nproc, nproc_x, nproc_y       &
     &,                     neighbour                                   &
     &,                     mype                                        &
     &,                     row_length, rows, n_rows                    &
     &,                     model_domain                                &
     &,                     model_levels,wet_levels                     &
     &,                     r_theta_levels,r_rho_levels                 &
     &,                     delta_lambda,delta_phi                      &
     &,                     FV_cos_theta_latitude,cos_v_latitude        &
     &,                     cos_theta_longitude,sin_theta_longitude     &
     &,                   D1(JTHETA(1)),D1(JU(1)),D1(JV(1)),D1(JW(0))   &
     &,                   D1(Jrho(1))                                   &
     &,                   D1(JQ(1)),D1(JQCL(1)),D1(JQCF(1))             &
     &,                   D1(JEXNER_THETA_LEVELS(1))                    &
!  passing D1(jpstar) as a dummy array
     &,                   D1(jpstar)                                    &
     &,                     dummy,dummy,lmass_cor0,Lmoist_cor0          &
     &,                     LEMQ_print                                  &
     &,                     A_energysteps,secs_per_stepim(atmos_im)     &
     &,                     A_REALHD(rh_tot_energy_init)                &
     &,                     A_REALHD(rh_tot_mass_init)                  &
     &,                     A_REALHD(rh_tot_m_init) )


! INITIAL RATE OF ENERGY CORRECTION TO ZERO

        A_REALHD(rh_energy_corr) = 0.0


      ENDIF
!
      RETURN
      END SUBROUTINE INIT_EMCORR

#endif
