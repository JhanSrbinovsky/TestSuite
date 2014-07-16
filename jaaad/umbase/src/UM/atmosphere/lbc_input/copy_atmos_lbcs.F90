#if defined(A31_1A)
#if defined(ATMOS) && !defined(GLOBAL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine COPY_ATMOS_LBCS
!
! Purpose : Copies LBC_TEND array to LBC array, ready for new
!           tendencies to be read into LBC_TEND
!
! Author :  P.Burton
!
! Model                Modification history from model version 5.2
! version    Date
! 5.2        02/10/00  New deck at 5.2                    P.Burton
! 5.3   12/10/01  extended haloes for exner in lbcs    Terry Davies
! 5.5   03/02/03  Include qcf2,qrain,qgraup lbcs. R.M.Forbes
! 6.0   29/07/03  Include cloud fractions in lbcs. Damian Wilson
! 6.2   01/10/04  Include murk aerosol lbcs.  R.M.Forbes
!
! ---------------------------------------------------------

      SUBROUTINE COPY_ATMOS_LBCS(                                       &
     &  LENRIM,                                                         &
     &  MODEL_LEVELS,WET_LEVELS,                                        &
     &  TR_VARS,TR_LEVELS,                                              &
     &  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, l_pc2_lbc,   &
     &  L_murk_lbc,                                                     &
     &  U_LBC, U_LBC_TEND,                                              &
     &  V_LBC, V_LBC_TEND,                                              &
     &  W_LBC, W_LBC_TEND,                                              &
     &  RHO_LBC, RHO_LBC_TEND,                                          &
     &  THETA_LBC, THETA_LBC_TEND,                                      &
     &  Q_LBC, Q_LBC_TEND,                                              &
     &  QCL_LBC, QCL_LBC_TEND,                                          &
     &  QCF_LBC, QCF_LBC_TEND,                                          &
     &  QCF2_LBC, QCF2_LBC_TEND,                                        &
     &  QRAIN_LBC, QRAIN_LBC_TEND,                                      &
     &  QGRAUP_LBC, QGRAUP_LBC_TEND,                                    &
     &  CF_BULK_LBC, CF_BULK_LBC_TEND,                                  &
     &  CF_LIQUID_LBC, CF_LIQUID_LBC_TEND,                              &
     &  CF_FROZEN_LBC, CF_FROZEN_LBC_TEND,                              &
     &  EXNER_LBC, EXNER_LBC_TEND,                                      &
     &  U_ADV_LBC, U_ADV_LBC_TEND,                                      &
     &  V_ADV_LBC, V_ADV_LBC_TEND,                                      &
     &  W_ADV_LBC, W_ADV_LBC_TEND,                                      &
     &  MURK_LBC, MURK_LBC_TEND,                                        &
     &  TRACER_LBCS, TRACER_LBCS_TEND                                   &
     &  )

      IMPLICIT NONE

! Parameters required for argument declarations
#include "parparm.h"

! Arguments:

      INTEGER                                                           &
     &  LENRIM(Nfld_max,NHalo_Max)                                      &
                                    ! IN : Size of a level of LBC
     &, MODEL_LEVELS                                                    &
                                    ! IN : Number of model levels
     &, WET_LEVELS                                                      &
                                    ! IN : Number of wet model levels
     &, TR_VARS                                                         &
                                    ! IN : Number of tracers
     &, TR_LEVELS                   ! IN : Number of tracer levels

      Logical, Intent (In) ::                                           &
     &  L_mcr_qcf2_lbc                                                  &
                          ! true if second cloud ice lbcs active
     &, L_mcr_qrain_lbc                                                 &
                          ! true if rain lbcs active
     &, L_mcr_qgraup_lbc                                                &
                          ! true if graupel lbcs active
     &, L_pc2_lbc                                                       &
                          ! true if cloud fractions in lbcs
     &, L_murk_lbc        ! true if murk aerosol in lbcs

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
     &          WET_LEVELS)                                             &
                                    ! IN : CF_BULK LBC tendency
     &, CF_LIQUID_LBC(LENRIM(fld_type_p,halo_type_extended),            &
     &          WET_LEVELS)                                             &
                                    ! IN/OUT : CF_LIQUID LBC
     &, CF_LIQUID_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),       &
     &          WET_LEVELS)                                             &
                                    ! IN : CF_LIQUID LBC tendency
     &, CF_FROZEN_LBC(LENRIM(fld_type_p,halo_type_extended),            &
     &          WET_LEVELS)                                             &
                                    ! IN/OUT : CF_FROZEN LBC
     &, CF_FROZEN_LBC_TEND(LENRIM(fld_type_p,halo_type_extended),       &
     &          WET_LEVELS)                                             &
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
     &, TRACER_LBCS(LENRIM(fld_type_p,halo_type_extended),              &
     &              TR_LEVELS,TR_VARS)                                  &
                                            ! IN/OUT : Tracer LBCs
     &, TRACER_LBCS_TEND(LENRIM(fld_type_p,halo_type_extended),         &
     &                   TR_LEVELS,TR_VARS) ! IN : Tracer LBCs tendency


! Local variables

      INTEGER                                                           &
     &  tracer                                                          &
                  ! loop counter over tracer variables
     &, k                                                               &
                  ! loop counter over levels
     &, i         ! loop counter over horizontal dimension

! ---------------------------------------------------------

      DO k=1,MODEL_LEVELS
        DO i=1,LENRIM(fld_type_u,halo_type_extended)
          U_LBC(i,k)=U_LBC_TEND(i,k)
          U_ADV_LBC(i,k)=U_ADV_LBC_TEND(i,k)
        ENDDO ! i

        DO i=1,LENRIM(fld_type_v,halo_type_extended)
          V_LBC(i,k)=V_LBC_TEND(i,k)
          V_ADV_LBC(i,k)=V_ADV_LBC_TEND(i,k)
        ENDDO ! i

        DO i=1,LENRIM(fld_type_p,halo_type_extended)
          RHO_LBC(i,k)=RHO_LBC_TEND(i,k)
          THETA_LBC(i,k)=THETA_LBC_TEND(i,k)
        ENDDO ! i
      ENDDO ! k

      DO k=1,WET_LEVELS
        DO i=1,LENRIM(fld_type_p,halo_type_extended)
          Q_LBC(i,k)=Q_LBC_TEND(i,k)
          QCL_LBC(i,k)=QCL_LBC_TEND(i,k)
          QCF_LBC(i,k)=QCF_LBC_TEND(i,k)
        ENDDO ! i
      ENDDO ! k

      IF (L_mcr_qcf2_lbc) THEN  ! qcf2 lbcs active
        DO k=1,WET_LEVELS
          DO i=1,LENRIM(fld_type_p,halo_type_extended)
            QCF2_LBC(i,k) = QCF2_LBC_TEND(i,k)
          ENDDO ! i
        ENDDO ! k
      ENDIF

      IF (L_mcr_qrain_lbc) THEN  ! qrain lbcs active
        DO k=1,WET_LEVELS
          DO i=1,LENRIM(fld_type_p,halo_type_extended)
            QRAIN_LBC(i,k) = QRAIN_LBC_TEND(i,k)
          ENDDO ! i
        ENDDO ! k
      ENDIF

      IF (L_mcr_qgraup_lbc) THEN  ! qgraup lbcs active
        DO k=1,WET_LEVELS
          DO i=1,LENRIM(fld_type_p,halo_type_extended)
            QGRAUP_LBC(i,k) = QGRAUP_LBC_TEND(i,k)
          ENDDO ! i
        ENDDO ! k
      ENDIF

      IF (L_pc2_lbc) THEN  ! cloud fraction lbc's active
        DO k=1,WET_LEVELS
          DO i=1,LENRIM(fld_type_p,halo_type_extended)
            CF_BULK_LBC(i,k)   = CF_BULK_LBC_TEND(i,k)
            CF_LIQUID_LBC(i,k) = CF_LIQUID_LBC_TEND(i,k)
            CF_FROZEN_LBC(i,k) = CF_FROZEN_LBC_TEND(i,k)
          ENDDO ! i
        ENDDO ! k
      ENDIF

      DO k=0,MODEL_LEVELS
        DO i=1,LENRIM(fld_type_p,halo_type_extended)
          W_LBC(i,k)=W_LBC_TEND(i,k)
          W_ADV_LBC(i,k)=W_ADV_LBC_TEND(i,k)
        ENDDO ! i
      ENDDO ! k

      DO k=1,MODEL_LEVELS+1
        DO i=1,LENRIM(fld_type_p,halo_type_extended)
           EXNER_LBC(i,k)=EXNER_LBC_TEND(i,k)
        ENDDO ! i
      ENDDO ! k

      IF (L_murk_lbc) THEN
      ! murk aerosol turned on and murk lbcs active
        DO k=1,WET_LEVELS
          DO i=1,LENRIM(fld_type_p,halo_type_single)
            MURK_LBC(i,k) = MURK_LBC_TEND(i,k)
          ENDDO ! i
        ENDDO ! k
      ENDIF

      DO tracer=1,TR_VARS
        DO k=1,TR_LEVELS
          DO i=1,LENRIM(fld_type_p,halo_type_extended)
            TRACER_LBCS(i,k,tracer)=TRACER_LBCS_TEND(i,k,tracer)
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE COPY_ATMOS_LBCS
#endif
#endif
