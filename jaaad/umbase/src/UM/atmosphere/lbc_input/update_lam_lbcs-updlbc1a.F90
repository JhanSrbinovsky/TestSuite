#if defined(A31_1A)
#if defined(ATMOS) && !defined(GLOBAL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Update the Lateral Boundary Conditions (LBCs) of LAM fields

      SUBROUTINE UPDATE_LAM_LBCs(                                       &
     &  r_rho_levels, r_theta_levels,                                   &
     &  ROW_LENGTH,ROWS,N_ROWS,MODEL_LEVELS,WET_LEVELS,                 &
     &  OFFX,OFFY,HALO_I,HALO_J,AT_EXTREMITY,                           &
     &  L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_pc2,                   &
     &  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, L_pc2_lbc,   &
     &  L_murk, L_murk_lbc,                                             &
     &  L_LBC_balance, L_int_uvw_lbc,                                   &
     &  RIMWIDTH,RIMWEIGHTS,LENRIM,LBC_SIZE,LBC_START,                  &
     &  THETA_LBC, Q_LBC, QCL_LBC, QCF_LBC,                             &
     &  QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,                                &
     &  CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC,                      &
     &  RHO_LBC, EXNER_LBC,                                             &
     &  U_LBC, V_LBC, W_LBC, U_ADV_LBC, V_ADV_LBC, W_ADV_LBC,           &
     &  MURK_LBC,                                                       &
     &  THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP,                        &
     &  CF_BULK, CF_LIQUID, CF_FROZEN,                                  &
     &  RHO, EXNER, U, V, W, U_ADV, V_ADV, W_ADV, MURK,                 &
     &  DELTA_PHI, DELTA_LAMBDA,                                        &
     &  BASE_PHI, BASE_LAMBDA,                                          &
     &  TWO_OMEGA, DATASTART,                                           &
     &  lat_rot_NP,                                                     &
     &  GLOBAL_ROW_LENGTH, GLOBAL_ROWS                                  &
     &     )

      IMPLICIT NONE
! Updates all the input fields with their LBCs
!
! Author : Paul Burton
!
! History
! Date       Version    Comment
! ----       -------    -------
!            5.1        New deck introduced             P.Burton
! 25/8/00    5.2        Corrected number of exner levels
!                       Corrected update of QCF         P.Burton
! 5.3   12/10/01  extended haloes for exner in lbcs    Terry Davies
! 5.5   03/02/03  Include code for qcf2,qrain,qgraup lbcs.  R.M.Forbes
! 6.0   30/07/03  Include cloud fraction in lbcs.  Damian Wilson
! 6.1   30/06/04  Change logicals L_mcr_q* to L_mcr_q*_lbc. R.M.Forbes
! 6.2   05/01/06  Include both L_mcr_q* and L_mcr_q*_lbc.  R.Forbes
! 6.2   01/10/04  Include code for murk aerosol lbcs.  R.M.Forbes
! 6.2   04/10/05  Changes for iterative semi-Lagrangian scheme.
!                                                        M. Diamantakis
! 6.4   12/12/06  Remove changes introduced at 6.2 for iterative
!                 semi-Lagrangian scheme (use a separate 
!                 subroutine for this instead)           M. Diamantakis

! Parameters required for dimensioning some of the arguments
#include "parparm.h"

! Arguments

      INTEGER                                                           &
     &  ROW_LENGTH                                                      &
                          ! IN : Length of a model row
     &, ROWS                                                            &
                          ! IN : Number of rows for theta,u fields
     &, N_ROWS                                                          &
                          ! IN : Number of rows for v fields
     &, MODEL_LEVELS                                                    &
                          ! IN : Number of model levels
     &, WET_LEVELS                                                      &
                          ! IN : Number of moist model levels
     &, OFFX                                                            &
                          ! IN : Size of "single" halo (EW direction)
     &, OFFY                                                            &
                          ! IN : Size of "single" halo (NS direction)
     &, HALO_I                                                          &
                          ! IN : Size of extended halo (EW direction)
     &, HALO_J                                                          &
                          ! IN : Size of extended halo (NS direction)
     &, RIMWIDTH                                                        &
                          ! IN : Size of boundary region
     &, LENRIM(Nfld_max,NHalo_max)                                      &
                          ! IN : Size of single level of LBC
     &, LBC_SIZE(4,Nfld_max,NHalo_max)                                  &
                          ! IN : Size of a side of a LBC
     &, LBC_START(4,Nfld_max,NHalo_max)                                 &
                          ! IN : Start of a side in a LBC
     &, DATASTART                                                       &
                          ! IN : position of first data point for PE
     &, GLOBAL_ROW_LENGTH                                               &
                          ! IN : no. of columns in LAM
     &, GLOBAL_ROWS                                                     
                          ! IN : no. or rows in LAM

      Logical, Intent (In) ::                                           &
     &  AT_EXTREMITY(4)                                                 &
                      ! IN : At an edge?
     &, L_mcr_qcf2                                                      &
                      ! true if prognostic 2nd cloud ice active
     &, L_mcr_qrain                                                     &
                      ! true if prognostic rain active
     &, L_mcr_qgraup                                                    &
                      ! true if prognostic graupel active
     &, L_mcr_qcf2_lbc                                                  &
                          ! true if prognostic 2nd cloud ice in lbcs
     &, L_mcr_qrain_lbc                                                 &
                          ! true if prognostic rain in lbcs
     &, L_mcr_qgraup_lbc                                                &
                          ! true if prognostic graupel in lbcs
     &, L_pc2                                                           &
                          ! true if prognostic cloud fracs active
     &, L_pc2_lbc                                                       &
                          ! true if prognostic cloud fracs in lbcs
     &, L_murk                                                          &
                          ! true if murk aerosol active
     &, L_murk_lbc                                                      &
                          ! true if murk aerosol lbcs active
     &, L_int_uvw_lbc                                                   &
                          ! true if using interpolated advecting winds 
                          !  in lateral boundaries
     &, L_LBC_balance     ! true if imposing balance in vertical 
                          !      momentum equation in LBC regions


      Real, Intent(In) ::                                               &
           ! vertical co-ordinate information
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:model_levels)             &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j, model_levels)

      Real, Intent (In) ::                                              &
     &  TWO_OMEGA                                                       &
                    ! twice angular velocity of Earth
     &, DELTA_PHI                                                       &
                    ! grid spacing (latitude)
     &, DELTA_LAMBDA                                                    &
                    ! grid spacing (longitude)
     &, BASE_PHI                                                        &
                    ! first latitude
     &, BASE_LAMBDA                                                     &
                    ! first longitude
     &, lat_rot_NP                                                      &
                    ! lat of rotated pole
     &, RIMWEIGHTS(RIMWIDTH)  ! weight to apply to LBC

      Real, Intent (InOut) ::                                           &
     &  THETA_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &            MODEL_LEVELS)                                         &
     &, Q_LBC(LENRIM(fld_type_p,halo_type_extended),                    &
     &        WET_LEVELS)                                               &
     &, QCL_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
     &          WET_LEVELS)                                             &
     &, QCF_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
     &          WET_LEVELS)                                             &
     &, QCF2_LBC(LENRIM(fld_type_p,halo_type_extended),                 &
     &          WET_LEVELS)                                             &
     &, QRAIN_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &          WET_LEVELS)                                             &
     &, QGRAUP_LBC(LENRIM(fld_type_p,halo_type_extended),               &
     &          WET_LEVELS)                                             &
     &, CF_BULK_LBC(LENRIM(fld_type_p,halo_type_extended),              &
     &          WET_LEVELS)                                             &
     &, CF_LIQUID_LBC(LENRIM(fld_type_p,halo_type_extended),            &
     &          WET_LEVELS)                                             &
     &, CF_FROZEN_LBC(LENRIM(fld_type_p,halo_type_extended),            &
     &          WET_LEVELS)                                             &
     &, RHO_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
     &          MODEL_LEVELS)                                           &
     &, EXNER_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &            MODEL_LEVELS+1)                                       &
     &, U_LBC(LENRIM(fld_type_u,halo_type_extended),                    &
     &        MODEL_LEVELS)                                             &
     &, V_LBC(LENRIM(fld_type_v,halo_type_extended),                    &
     &        MODEL_LEVELS)                                             &
     &, W_LBC(LENRIM(fld_type_p,halo_type_extended),                    &
     &        0:MODEL_LEVELS)                                           &
     &, U_ADV_LBC(LENRIM(fld_type_u,halo_type_extended),                &
     &            MODEL_LEVELS)                                         &
     &, V_ADV_LBC(LENRIM(fld_type_v,halo_type_extended),                &
     &            MODEL_LEVELS)                                         &
     &, W_ADV_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &            0:MODEL_LEVELS)                                       &
     &, MURK_LBC(LENRIM(fld_type_p,halo_type_single),                   &
     &            MODEL_LEVELS)

      Real, Intent (InOut) ::                                           &
     &  THETA(1-OFFX:ROW_LENGTH+OFFX,                                   &
     &        1-OFFY:ROWS+OFFY,                                         &
     &        MODEL_LEVELS)                                             &
     &, Q(1-HALO_I:ROW_LENGTH+HALO_I,                                   &
     &    1-HALO_J:ROWS+HALO_J,                                         &
     &    WET_LEVELS)                                                   &
     &, QCL(1-HALO_I:ROW_LENGTH+HALO_I,                                 &
     &      1-HALO_J:ROWS+HALO_J,                                       &
     &      WET_LEVELS)                                                 &
     &, QCF(1-HALO_I:ROW_LENGTH+HALO_I,                                 &
     &      1-HALO_J:ROWS+HALO_J,                                       &
     &      WET_LEVELS)                                                 &
     &, QCF2(1-HALO_I:ROW_LENGTH+HALO_I,                                &
     &      1-HALO_J:ROWS+HALO_J,                                       &
     &      WET_LEVELS)                                                 &
     &, QRAIN(1-HALO_I:ROW_LENGTH+HALO_I,                               &
     &      1-HALO_J:ROWS+HALO_J,                                       &
     &      WET_LEVELS)                                                 &
     &, QGRAUP(1-HALO_I:ROW_LENGTH+HALO_I,                              &
     &      1-HALO_J:ROWS+HALO_J,                                       &
     &      WET_LEVELS)                                                 &
     &, CF_BULK(1-HALO_I:ROW_LENGTH+HALO_I,                             &
     &      1-HALO_J:ROWS+HALO_J,                                       &
     &      WET_LEVELS)                                                 &
     &, CF_LIQUID(1-HALO_I:ROW_LENGTH+HALO_I,                           &
     &      1-HALO_J:ROWS+HALO_J,                                       &
     &      WET_LEVELS)                                                 &
     &, CF_FROZEN(1-HALO_I:ROW_LENGTH+HALO_I,                           &
     &      1-HALO_J:ROWS+HALO_J,                                       &
     &      WET_LEVELS)                                                 &
     &, RHO(1-OFFX:ROW_LENGTH+OFFX,                                     &
     &      1-OFFY:ROWS+OFFY,                                           &
     &      MODEL_LEVELS)                                               &
     &, EXNER(1-OFFX:ROW_LENGTH+OFFX,                                   &
     &        1-OFFY:ROWS+OFFY,                                         &
     &            MODEL_LEVELS+1)                                       &
     &, U(1-OFFX:ROW_LENGTH+OFFX,                                       &
     &    1-OFFY:ROWS+OFFY,                                             &
     &    MODEL_LEVELS)                                                 &
     &, V(1-OFFX:ROW_LENGTH+OFFX,                                       &
     &    1-OFFY:N_ROWS+OFFY,                                           &
     &    MODEL_LEVELS)                                                 &
     &, W(1-OFFX:ROW_LENGTH+OFFX,                                       &
     &    1-OFFY:ROWS+OFFY,                                             &
     &    0:MODEL_LEVELS)                                               &
     &, U_ADV(1-HALO_I:ROW_LENGTH+HALO_I,                               &
     &        1-HALO_J:ROWS+HALO_J,                                     &
     &        MODEL_LEVELS)                                             &
     &, V_ADV(1-HALO_I:ROW_LENGTH+HALO_I,                               &
     &        1-HALO_J:N_ROWS+HALO_J,                                   &
     &        MODEL_LEVELS)                                             &
     &, W_ADV(1-HALO_I:ROW_LENGTH+HALO_I,                               &
     &        1-HALO_J:ROWS+HALO_J,                                     &
     &        0:MODEL_LEVELS)                                           &
     &, MURK(1-OFFX:ROW_LENGTH+OFFX,                                    &
     &       1-OFFY:ROWS+OFFY,                                          &
     &       MODEL_LEVELS)

! Local variables

      LOGICAL                                                           &
     &  L_Do_Boundaries                                                 &
     &, L_Do_Halos
      INTEGER                                                           &
     &  N_RIMS_TO_DO

! ----------------------------------------------------------------------

       L_Do_Boundaries = .TRUE.
       L_Do_Halos = .TRUE.
       N_RIMS_TO_DO = RIMWIDTH


       IF(L_LBC_balance)THEN

! Reset LBC array values for Exner, rho and w 
! DEPENDS ON: balance_lbc_values
         CALL BALANCE_LBC_VALUES(                                       &
     &    EXNER_LBC, RHO_LBC, THETA_LBC,Q_LBC, W_LBC, W_ADV_LBC,        &
     &    U_LBC, V_LBC,                                                 &
     &    R_RHO_LEVELS, R_THETA_LEVELS,                                 &
     &    ROW_LENGTH, ROWS, model_levels, wet_levels, HALO_I, HALO_J,   &
     &    LENRIM(fld_type_p,halo_type_extended),                        &
     &    LENRIM(fld_type_u,halo_type_extended),                        &
     &    LENRIM(fld_type_v,halo_type_extended),                        &
     &    LBC_START(1,fld_type_p,halo_type_extended),                   &
     &    LBC_START(1,fld_type_u,halo_type_extended),                   &
     &    LBC_START(1,fld_type_v,halo_type_extended),                   &
     &    RIMWIDTH, N_RIMS_TO_DO, RIMWEIGHTS, AT_EXTREMITY,             &
     &    DELTA_PHI, DELTA_LAMBDA,                                      &
     &    BASE_PHI, BASE_LAMBDA,                                        &
     &    TWO_OMEGA, DATASTART,                                         &
     &    LAT_ROT_NP,                                                   &
     &    GLOBAL_ROW_LENGTH, GLOBAL_ROWS, L_int_uvw_lbc                 &
     &    )

       ENDIF
       
! Theta
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,OFFX,OFFY,MODEL_LEVELS,fld_type_p,THETA,        &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,THETA_LBC,                                        &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! Q
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,Q,          &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,Q_LBC,                                            &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! QCL
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,QCL,        &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,QCL_LBC,                                          &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! QCF
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,QCF,        &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,QCF_LBC,                                          &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! QCF2 - Update if prognostic is active and data is in lbcs
      IF (L_mcr_qcf2 .AND. L_mcr_qcf2_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,QCF2,       &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,QCF2_LBC,                                         &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)
      End If

! QRAIN - Update if prognostic is active and data is in lbcs
      IF (L_mcr_qrain .AND. L_mcr_qrain_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,QRAIN,      &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,QRAIN_LBC,                                        &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)
      End If

! QGRAUP - Update if prognostic is active and data is in lbcs
      IF (L_mcr_qgraup .AND. L_mcr_qgraup_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,QGRAUP,     &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,QGRAUP_LBC,                                       &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)
      End If

! Cloud fractions - Update if prognostics are active and data in lbcs
      IF (L_pc2 .AND. L_pc2_lbc) THEN

! CF_BULK
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,CF_BULK,    &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,CF_BULK_LBC,                                      &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! CF_LIQUID
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,CF_LIQUID,  &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,CF_LIQUID_LBC,                                    &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! CF_FROZEN
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,CF_FROZEN,  &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,CF_FROZEN_LBC,                                    &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

      End If ! L_pc2 .AND. L_pc2_lbc

! RHO
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,OFFX,OFFY,MODEL_LEVELS,fld_type_p,RHO,          &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,RHO_LBC,                                          &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! EXNER
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,OFFX,OFFY,MODEL_LEVELS+1,fld_type_p,EXNER,      &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  halo_i,halo_j,EXNER_LBC,                                        &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! U
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,OFFX,OFFY,MODEL_LEVELS,fld_type_u,U,            &
     &  LENRIM(fld_type_u,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_u,halo_type_extended),                      &
     &  LBC_START(1,fld_type_u,halo_type_extended),                     &
     &  HALO_I, HALO_J,U_LBC,                                           &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! V
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,N_ROWS,OFFX,OFFY,MODEL_LEVELS,fld_type_v,V,          &
     &  LENRIM(fld_type_v,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_v,halo_type_extended),                      &
     &  LBC_START(1,fld_type_v,halo_type_extended),                     &
     &  HALO_I, HALO_J,V_LBC,                                           &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)
! W
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,OFFX,OFFY,MODEL_LEVELS+1,fld_type_p,W,          &
     &  LENRIM(fld_type_p,halo_type_extended),                          &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I, HALO_J,W_LBC,                                           &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

       IF ( .NOT. L_int_uvw_lbc ) THEN

! U_ADV
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,MODEL_LEVELS,fld_type_u,          &
     &  U_ADV,LENRIM(fld_type_u,halo_type_extended),                    &
     &  LBC_SIZE(1,fld_type_u,halo_type_extended),                      &
     &  LBC_START(1,fld_type_u,halo_type_extended),                     &
     &  HALO_I,HALO_J,U_ADV_LBC,                                        &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! V_ADV
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,N_ROWS,HALO_I,HALO_J,MODEL_LEVELS,fld_type_v,        &
     &  V_ADV,LENRIM(fld_type_v,halo_type_extended),                    &
     &  LBC_SIZE(1,fld_type_v,halo_type_extended),                      &
     &  LBC_START(1,fld_type_v,halo_type_extended),                     &
     &  HALO_I,HALO_J,V_ADV_LBC,                                        &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

! W_ADV
! DEPENDS ON: set_lateral_boundaries
       CALL SET_LATERAL_BOUNDARIES(                                     &
     &  ROW_LENGTH,ROWS,HALO_I,HALO_J,MODEL_LEVELS+1,fld_type_p,        &
     &  W_ADV,LENRIM(fld_type_p,halo_type_extended),                    &
     &  LBC_SIZE(1,fld_type_p,halo_type_extended),                      &
     &  LBC_START(1,fld_type_p,halo_type_extended),                     &
     &  HALO_I,HALO_J,W_ADV_LBC,                                        &
     &  RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                  &
     &  L_Do_Boundaries,L_Do_Halos)

       ENDIF ! .NOT. L_int_uvw_lbc

! MURK Aerosol
      ! If murk aerosol is in use and murk lbcs are active,
      IF (L_murk .AND. L_murk_lbc) THEN
! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
     &   ROW_LENGTH,ROWS,OFFX,OFFY,MODEL_LEVELS,fld_type_p,MURK,        &
     &   LENRIM(fld_type_p,halo_type_single),                           &
     &   LBC_SIZE(1,fld_type_p,halo_type_single),                       &
     &   LBC_START(1,fld_type_p,halo_type_single),                      &
     &   OFFX,OFFY,MURK_LBC,                                            &
     &   RIMWIDTH,N_RIMS_TO_DO,RIMWEIGHTS,AT_EXTREMITY,                 &
     &   L_Do_Boundaries,L_Do_Halos)
      End If

       RETURN
       END SUBROUTINE UPDATE_LAM_LBCs
#endif
#endif
