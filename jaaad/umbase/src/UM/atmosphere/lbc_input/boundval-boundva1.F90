#if defined(A31_1A)
#if defined(ATMOS) && !defined(GLOBAL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine BOUNDVAL
!
! Purpose : Checks whether a boundary incrementing step and increments
!           boundary values if required.
!
! Author :  Z. Gardner
!
! Model             Modification history from model version 5.2
! version    Date
! 5.2      04/12/00 Deck rewritten for 5.x format LBCs      Z. Gardner
! 5.3      11/02/02 Bug fix for rim_stepsa = 2     Z. Gardner
! 5.3   12/10/01  extended haloes for exner in lbcs    Terry Davies
! 5.5      02/03/02  Include qcf2,qrain,qgraup lbcs
!                    and associated logicals           R.M.Forbes
! 6.0      29/07/03 Cloud fractions in lbcs            Damian Wilson
! 6.2      30/03/05 Use BNDARY_OFFSETim to calculate
!                   steps_to_next_update.              Dave Robinson
! 6.2      01/10/04 Include murk aerosol lbc.     R.M.Forbes
!
! ---------------------------------------------------------
!
      SUBROUTINE BOUNDVAL(                                              &
     &  LENRIM                                                          &
     &, L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, L_pc2_lbc    &
     &, L_murk_lbc, L_int_uvw_lbc                                       &
     &, U_LBC, U_LBC_TEND                                               &
     &, V_LBC, V_LBC_TEND                                               &
     &, W_LBC, W_LBC_TEND                                               &
     &, RHO_LBC, RHO_LBC_TEND                                           &
     &, THETA_LBC, THETA_LBC_TEND                                       &
     &, Q_LBC, Q_LBC_TEND                                               &
     &, QCL_LBC, QCL_LBC_TEND                                           &
     &, QCF_LBC, QCF_LBC_TEND                                           &
     &, QCF2_LBC, QCF2_LBC_TEND                                         &
     &, QRAIN_LBC, QRAIN_LBC_TEND                                       &
     &, QGRAUP_LBC, QGRAUP_LBC_TEND                                     &
     &, CF_BULK_LBC, CF_BULK_LBC_TEND                                   &
     &, CF_LIQUID_LBC, CF_LIQUID_LBC_TEND                               &
     &, CF_FROZEN_LBC, CF_FROZEN_LBC_TEND                               &
     &, EXNER_LBC, EXNER_LBC_TEND                                       &
     &, U_ADV_LBC, U_ADV_LBC_TEND                                       &
     &, V_ADV_LBC, V_ADV_LBC_TEND                                       &
     &, W_ADV_LBC, W_ADV_LBC_TEND                                       &
     &, MURK_LBC, MURK_LBC_TEND                                         &
     &, TRACER_LBC, TRACER_LBC_TEND                                     &
     &, IO1, IO2                                                        &
     &, ICODE, CMESSAGE)

      IMPLICIT NONE

! Parameters required for argument declarations
#include "cmaxsize.h"
#include "parparm.h"
#include "typsize.h"
#include "typbnd.h"
#include "csubmodl.h"
#include "clookadd.h"
#include "ctime.h"

! Arguments:

      INTEGER                                                           &
     &  LENRIM(Nfld_max,NHalo_Max)                                      &
                                  ! IN : Size of a level of LBC
     &, IO1, IO2                  ! Offsets to allow for old or new lbcs

      Logical, Intent (In) ::                                           &
     &  L_mcr_qcf2_lbc                                                  &
                         ! true if using second cloud ice lbcs
     &, L_mcr_qrain_lbc                                                 &
                         ! true if using rain lbcs
     &, L_mcr_qgraup_lbc                                                &
                         ! true if using graupel lbcs
     &, L_pc2_lbc                                                       &
                         ! true if cloud fractions in lbcs
     &, L_int_uvw_lbc                                                   &
                          ! true if using interpolated advecting winds 
                          !  in lateral boundaries
     &, L_murk_lbc       ! true if murk aerosol in lbcs

! Note: LBCs are the current value of the LBC (and will be updated)
!       LBC_TENDs are the value of the LBC at the end of the LBC
!       period, towards which the LBC will tend.

      Real, Intent (InOut) ::                                           &
     &  U_LBC(LENRIM(fld_type_u,halo_type_extended),                    &
     &        MODEL_LEVELS)                                             &
                                    ! IN/OUT : U LBC
     &, U_LBC_TEND(LENRIM(fld_type_u,halo_type_extended),               &
     &             MODEL_LEVELS)                                        &
                                    ! IN : U LBC tendency
     &, V_LBC(LENRIM(fld_type_v,halo_type_extended),                    &
     &        MODEL_LEVELS)                                             &
                                    ! IN/OUT : V LBC
     &, V_LBC_TEND(LENRIM(fld_type_v,halo_type_extended),               &
     &             MODEL_LEVELS)                                        &
                                    ! IN : V LBC tendency
     &, W_LBC(LENRIM(fld_type_p,halo_type_extended),                    &
     &        0:MODEL_LEVELS)                                           &
                                    ! IN/OUT : V LBC
     &, W_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),               &
     &             0:MODEL_LEVELS)                                      &
                                    ! IN : V LBC tendency
     &, RHO_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
     &          MODEL_LEVELS)                                           &
                                    ! IN/OUT : Rho LBC
     &, RHO_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),             &
     &               MODEL_LEVELS)                                      &
                                    ! IN : Rho LBC tendency
     &, THETA_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &            MODEL_LEVELS)                                         &
                                    ! IN/OUT : Theta LBC
     &, THETA_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),           &
     &                 MODEL_LEVELS)                                    &
                                    ! IN : Theta LBC tendency
     &, Q_LBC(LENRIM(fld_type_p,halo_type_extended),                    &
     &        WET_LEVELS)                                               &
                                    ! IN/OUT : Q LBC
     &, Q_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),               &
     &             WET_LEVELS)                                          &
                                    ! IN : Q LBC tendency
     &, QCL_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
     &          WET_LEVELS)                                             &
                                    ! IN/OUT : QCL LBC
     &, QCL_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),             &
     &               WET_LEVELS)                                        &
                                    ! IN : QCL LBC tendency
     &, QCF_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
     &          WET_LEVELS)                                             &
                                    ! IN/OUT : QCL LBC
     &, QCF_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),             &
     &               WET_LEVELS)                                        &
                                    ! IN : QCL LBC tendency
     &, QCF2_LBC(LENRIM(fld_type_p,halo_type_extended),                 &
     &          WET_LEVELS)                                             &
                                    ! IN/OUT : QCF2 LBC
     &, QCF2_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),            &
     &               WET_LEVELS)                                        &
                                    ! IN : QCF2 LBC tendency
     &, QRAIN_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &          WET_LEVELS)                                             &
                                    ! IN/OUT : QRAIN LBC
     &, QRAIN_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),           &
     &               WET_LEVELS)                                        &
                                    ! IN : QRAIN LBC tendency
     &, QGRAUP_LBC(LENRIM(fld_type_p,halo_type_extended),               &
     &          WET_LEVELS)                                             &
                                    ! IN/OUT : QGRAUP LBC
     &, QGRAUP_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),          &
     &               WET_LEVELS)                                        &
                                    ! IN : QGRAUP LBC tendency
     &, CF_BULK_LBC(LENRIM(fld_type_p,halo_type_extended),              &
     &          WET_LEVELS)                                             &
                                    ! IN/OUT : CF_BULK LBC
     &, CF_BULK_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),         &
     &               WET_LEVELS)                                        &
                                    ! IN : CF_BULK LBC tendency
     &, CF_LIQUID_LBC(LENRIM(fld_type_p,halo_type_extended),            &
     &          WET_LEVELS)                                             &
                                    ! IN/OUT : CF_LIQUID LBC
     &, CF_LIQUID_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),       &
     &               WET_LEVELS)                                        &
                                    ! IN : CF_LIQUID LBC tendency
     &, CF_FROZEN_LBC(LENRIM(fld_type_p,halo_type_extended),            &
     &          WET_LEVELS)                                             &
                                    ! IN/OUT : CF_FROZEN LBC
     &, CF_FROZEN_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),       &
     &               WET_LEVELS)                                        &
                                    ! IN : CF_FROZEN LBC tendency
     &, EXNER_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &            MODEL_LEVELS+1)                                       &
                                    ! IN/OUT : Exner LBC
     &, EXNER_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),           &
     &                 MODEL_LEVELS+1)                                  &
                                         ! IN : Exner LBC tendency
     &, U_ADV_LBC(LENRIM(fld_type_u,halo_type_extended),                &
     &            MODEL_LEVELS)                                         &
                                    ! IN/OUT : U_ADV LBC
     &, U_ADV_LBC_TEND(LENRIM(fld_type_u,halo_type_extended),           &
     &                 MODEL_LEVELS)                                    &
                                    ! IN : U_ADV LBC tendency
     &, V_ADV_LBC(LENRIM(fld_type_v,halo_type_extended),                &
     &            MODEL_LEVELS)                                         &
                                    ! IN/OUT : V_ADV LBC
     &, V_ADV_LBC_TEND(LENRIM(fld_type_v,halo_type_extended),           &
     &                 MODEL_LEVELS)                                    &
                                    ! IN : V_ADV LBC tendency
     &, W_ADV_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &            0:MODEL_LEVELS)                                       &
                                    ! IN/OUT : W LBC
     &, W_ADV_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),           &
     &                 0:MODEL_LEVELS)                                  &
                                       ! IN : W LBC tendency
     &, MURK_LBC(LENRIM(fld_type_p,halo_type_single),                   &
     &           MODEL_LEVELS)                                          &
                                    ! IN/OUT : MURK LBC
     &, MURK_LBC_TEND(LENRIM(fld_type_p,halo_type_single),              &
     &                MODEL_LEVELS)                                     &
                                    ! IN : MURK LBC tendency
     &, TRACER_LBC(LENRIM(fld_type_p,halo_type_extended),               &
     &              TR_LEVELS,TR_VARS)                                  &
                                            ! IN/OUT : Tracer LBCs
     &, TRACER_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),          &
     &                   TR_LEVELS,TR_VARS) ! IN : Tracer LBCs tendency

      INTEGER                                                           &
     &  ICODE                           ! Error code

      CHARACTER*(80)                                                    &
     &  CMESSAGE                        ! Error message

! Local variables

      Integer       :: steps_to_next_update

! --------------------------------------------------------------------

      CMESSAGE=' '

! Check that the lbcs haven't lost synchronisation and then call
! increment_lbcs to increment them

      If (current_lbc_step  /=  stepim(atmos_im) + IO1) Then

        IF ((current_LBC_STEP /= STEPim(atmos_im) - IO2) .and.          &
     &         (MOD(BNDARY_OFFSETim(atmos_im)+STEPim(atmos_im) + IO1,   &
     &              RIM_STEPSA) /= 2) .and.                             &
     &         (MOD(BNDARY_OFFSETim(atmos_im)+STEPim(atmos_im) + IO1,   &
     &              RIM_STEPSA) /= 0)) Then

          ICODE=10
! DEPENDS ON: Ereport
          Call Ereport("BOUNDVA1 ", ICODE,                              &
     &                 "LBC values have lost synchronisation" )
        ENDIF

        steps_to_next_update = RIM_STEPSA -                             &
     &               MOD(BNDARY_OFFSETim(atmos_im)+STEPim(atmos_im)-1,  &
     &                   RIM_STEPSA) + IO2

! DEPENDS ON: increment_atmos_lbcs
        CALL INCREMENT_ATMOS_LBCS(                                      &
     &        steps_to_next_update,                                     &
     &        LENRIM,                                                   &
     &        MODEL_LEVELS,WET_LEVELS,TR_VARS,TR_LEVELS,                &
     &        L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc,        &
     &        l_pc2_lbc, L_murk_lbc, L_int_uvw_lbc,                     &
     &        U_LBC,U_LBC_TEND,                                         &
     &        V_LBC,V_LBC_TEND,                                         &
     &        W_LBC,W_LBC_TEND,                                         &
     &        RHO_LBC,RHO_LBC_TEND,                                     &
     &        THETA_LBC,THETA_LBC_TEND,                                 &
     &        Q_LBC,Q_LBC_TEND,                                         &
     &        QCL_LBC,QCL_LBC_TEND,                                     &
     &        QCF_LBC,QCF_LBC_TEND,                                     &
     &        QCF2_LBC,QCF2_LBC_TEND,                                   &
     &        QRAIN_LBC,QRAIN_LBC_TEND,                                 &
     &        QGRAUP_LBC,QGRAUP_LBC_TEND,                               &
     &        CF_BULK_LBC,CF_BULK_LBC_TEND,                             &
     &        CF_LIQUID_LBC,CF_LIQUID_LBC_TEND,                         &
     &        CF_FROZEN_LBC,CF_FROZEN_LBC_TEND,                         &
     &        EXNER_LBC,EXNER_LBC_TEND,                                 &
     &        U_ADV_LBC,U_ADV_LBC_TEND,                                 &
     &        V_ADV_LBC,V_ADV_LBC_TEND,                                 &
     &        W_ADV_LBC,W_ADV_LBC_TEND,                                 &
     &        MURK_LBC,MURK_LBC_TEND,                                   &
     &        TRACER_LBC,TRACER_LBC_TEND,                               &
     &        RIM_STEPSA,ICODE,CMESSAGE)

!       IF (ICODE  >   0) THEN
!             WRITE(6,*) 'Failure in INCREMENT_ATMOS_LBCS while ',      &
!    &                   'attempting to update LBCS for timestep ',     &
!    &                   STEPim(a_im)
!         GOTO 9999
!       ENDIF

        current_LBC_STEP = STEPim(atmos_im) + IO1

      ENDIF ! IF (current_LBC_STEP  /=  STEPim(atmos_im))

 9999 CONTINUE
      RETURN
      END SUBROUTINE BOUNDVAL
#endif
#endif
