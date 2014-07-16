#if defined(A14_1B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE CAL_ENG_MASS_CORR--------------------------------------
!LL
!LL  PURPOSE : PART OF ENERGY CORRECTION SUITE OF ROUTINES
!LL            - TO CALCUATE THE NECESSARY CORRECTION TO
!LL              TEMPERATURE TO CONSERVE TOTAL ENERGY
!LL
!LL  NOT SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  MODEL            MODIFICATION HISTORY:
!LL VERSION  DATE
!LL   5.1  3/02/00  : Altered for new dynamics, removing assumption
!LL                   energy correction is done daily
!  6.2   25/12/05  Variable resolution changes          Yongming Tang
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL
!LL  DOCUMENTATION :
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE CAL_ENG_MASS_CORR (row_length,rows,me                  &
     &,                             A_ENERGYSTEPS,timestep              &
     &,                             delta_lambda,delta_phi              &
     &,                             halo_i, halo_j                      &
     &,                             Sum_energy_fluxes                   &
     &,                             TOT_MASS_INIT,TOT_ENERGY_INIT       &
     &,                             ENERGY_CORR,tot_energy_final)

!
      IMPLICIT NONE
!
#include "c_r_cp.h"
#include "c_a.h"
#include "c_pi.h"

!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!----------------------------------------------------------------------
!
      Integer                                                           &
     &   row_length                                                     &
                          ! row length
     &,  rows                                                           &
                          ! umber of model rows
     &,  me                                                             &
                          ! processor number
     &,  A_ENERGYSTEPS                                                  &
                          ! number of timesteps per energy period
     &,  halo_i                                                         &
     &,  halo_j
      Real                                                              &
     &   timestep                                                       &
                          ! model timestep in seconds
     &,  delta_lambda                                                   &
                          ! EW GRID SPACING IN radians
     &,  delta_phi                                                      &
                          ! NS GRID SPACING IN radians
     &, SUM_ENERGY_FLUXES(row_length,rows) ! total energy flux into
!                                          ! atmosphere
!
      REAL TOT_ENERGY_INIT      ! IN TOTAL ENERGY OF ATMOSPHERE
                                ! AT START OF energy correction period
!
      REAL TOT_MASS_INIT        ! IN TOTAL dry MASS OF ATMOSPHERE
                                ! AT START OF MODEL
!
      REAL TOT_ENERGY_FINAL     ! IN TOTAL ENERGY OF ATMOSPHERE
                                ! AT END OF energy corr period
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE IN AND OUT
!----------------------------------------------------------------------
!
      REAL ENERGY_CORR          ! INOUT ENERGY CORRECTION FACTOR
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!----------------------------------------------------------------------
!
      REAL CHG_ENERGY          ! CHANGE IN ENERGY
!
      REAL ERROR_ENERGY                                                 &
                               ! ERROR IN ENERGY CALCULATION
     &,    factor                                                       &
                               ! grid resolution factor
     &,    tot_fluxes                                                   &
                               ! global total energy flux
     &,    flux_corr           ! flux correction
!
!----------------------------------------------------------------------
! INTERNAL LOOP COUNTERS
!----------------------------------------------------------------------
!
      INTEGER I,J                ! LOOP COUNTER
!
!
!----------------------------------------------------------------------
! EXTERNAL SUBROUTINE CALLS  -
!----------------------------------------------------------------------
       External do_sums
!
!*---------------------------------------------------------------------
!
!======================================================================
! CALCULATION OF TEMPERATURE CHANGE TO BE ADDED OVER THE NEXT DAY
! TO CORRECT FOR ERROR IN ENERGY BUDGET DURING PRESENT DAY
!======================================================================
!

! Do global sum of heat fluxes added to the atmosphere during period
! This sum of fluxes is accumulated in D1(jnet_flux)
! Note no halo on input array, grid_type 1, halo_type 3

      tot_fluxes = 0.0

! DEPENDS ON: do_sums
      Call do_sums(sum_energy_fluxes, row_length, rows, 0, 0, 1, 3      &
     &,                1,tot_fluxes)

! add on energy correction added during last energy correction period

      factor=delta_lambda*delta_phi   ! already in radians
      tot_fluxes = tot_fluxes*factor

      tot_fluxes = tot_fluxes +                                         &
     &     (CP-R)*ENERGY_CORR*TOT_MASS_INIT*a_energysteps*timestep
!
!----------------------------------------------------------------------
! CALCULATE ENERGY CHANGE DURING THE energy period
!----------------------------------------------------------------------
!
      CHG_ENERGY = TOT_ENERGY_FINAL - TOT_ENERGY_INIT
!
!
!----------------------------------------------------------------------
! CALCULATE ERROR = DIFFERENCE BETWEEN CHANGE IN TOTAL ENERGY
! DURING THE period AND THAT EXPECTED DUE TO THE FLUXES OF ENERGY INTO
! THE ATMOSPHERE
!-----------------------------------------------------------------------
!
      ERROR_ENERGY = TOT_FLUXES - CHG_ENERGY
!
!
!-----------------------------------------------------------------------
! CALCULATE TEMPERATURE CHANGE TO BE APPLIED OVER THE
! NEXT period TO CORRECT FOR THE ERROR
! Now use Cv instead of Cp
!-----------------------------------------------------------------------

!      ENERGY_CORR = ERROR_ENERGY / (CP*TOT_MASS_INIT)
      ENERGY_CORR = ERROR_ENERGY / ((Cp-R)*TOT_MASS_INIT)
      FLUX_corr = ERROR_ENERGY/                                         &
     &(4.*pi*Earth_Radius*Earth_Radius*a_energysteps*timestep)
!
!----------------------------------------------------------------------
! CALCULATE RATE OF TEMPERATURE CHANGE WHICH NEEDS TO BE
! APPLIED OVER NEXT PERIOD TO CORRECT FOR ERROR
!----------------------------------------------------------------------
!
      ENERGY_CORR = ENERGY_CORR / (a_energysteps*timestep)
!
!
!----------------------------------------------------------------------
! DIAGNOSTICS
!----------------------------------------------------------------------
!


      If (me  ==  0)                                                    &
     &             WRITE(6,100) TOT_ENERGY_FINAL,TOT_ENERGY_INIT,       &
     &             CHG_ENERGY,TOT_FLUXES,ERROR_ENERGY,                  &
     &             ENERGY_CORR*86400.0,ENERGY_CORR,                     &
     &             flux_corr


 100  FORMAT (1X,'FINAL TOTAL ENERGY           = ',E13.5,' J/ '/        &
     &        1X,'INITIAL TOTAL ENERGY         = ',E13.5,' J/ '/        &
     &        1X,'CHG IN TOTAL ENERGY OVER Per = ',E13.5,' J/ '/        &
     &        1X,'FLUXES INTO ATM OVER PERIOD  = ',E13.5,' J/ '/        &
     &        1X,'ERROR IN ENERGY BUDGET       = ',E13.5,' J/ '/        &
     &        1X,'TEMP CORRECTION OVER A DAY   = ',E13.5,' K    '/      &
     &        1X,'TEMPERATURE CORRECTION RATE  = ',E13.5,' K/S  '/      &
     &        1X,'FLUX CORRECTION (ATM)        = ',E13.5,' W/M2 ')
!
!----------------------------------------------------------------------
! Note dry mass correction applied in eng_mass_cal routine
!----------------------------------------------------------------------

      RETURN
      END SUBROUTINE CAL_ENG_MASS_CORR
#endif
