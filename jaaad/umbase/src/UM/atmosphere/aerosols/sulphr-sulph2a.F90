#if defined(A17_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      SUBROUTINE SULPHR(                                                &
!
! Arguments IN
     &             halo_i, halo_j, off_x, off_y                         &
     &,            row_length, rows                                     &
     &,            model_levels, wet_model_levels                       &
     &, theta_field_size                                                &
     &,            TSTEP                                                &
     &,            CLOUDF, COSZA2D                                      &
     &,            P, T, Q, QCL, QCF                                    &
     &,            OH_CONC, H2O2_LMT, HO2_CONC, O3                      &
     &,            n_droplet                                            &
     &,            L_SULPC_DMS, L_SULPC_NEWDMS                          &
     &,            L_SULPC_OZONE, L_SULPC_NH3                           &
! Arguments IN/OUT
     &,            SO2, DMS, SO4_AIT, SO4_ACC, SO4_DIS                  &
     &,            NH3, H2O2_MXR                                        &
! Arguments OUT
     &,            MSA, NH3_DEP                                         &
     &,            F_DMS_TO_SO2, F_DMS_TO_SO4, F_DMS_TO_MSA             &
     &,            DELTAS_DRY, DELTAS_WET, DELTAS_WET_O3                &
     &,            DELTAS_TOT, DELTAS_DMS                               &
     &,            DELTAS_EVAP, DELTAS_NUCL, DELTAS_DIFFUSE             &
     &,            DELTAS_COAG, PSI                                     &
     &            )
!
!---------------------------------------------------------------------
! Purpose: To perform oxidation chemistry of sulphur dioxide to 3 modes
!          of Sulphate aerosol (Aitken, accumulation and dissolved),
!          and dimethyl sulphide to sulphur dioxide and methyl sulphonic
!          acid.
!          There is also exchange between the 3 modes of sulphate
!          aerosol due to nucleation, diffusion and evaporation
!          processes.
!          Called by Aero_Ctl
!
! Current owners of code:          C.E. Johnson and J.G.L. Rae
!
! History:
! Version     Date     Comment
! -------     ----     -------
!  5.2    2/11/00   Change array dimensions and DO loops to 3D to
!                   conform to ND model style.                 MJW
!
!  5.2    8/11/00   Add NUMBER_DROPLET function to calculate
!                   diffusional scavenging parameter for SO4_AIT
!                                                  A Jones, M Woodage
!  5.2    8/11/00   Add new DMS oxidation scheme; uses O3 as oxidant
!                   and produces SO2, MSA, SO4_AIT, SO4_ACC
!                                                C Johnson, M Woodage
!  5.3    29/8/01   Introduce switch to control calculation of droplet
!                   number.                        A Jones
!  5.3  09/10/01  Improve diffusional scavenging of SO4_AIT.
!                                                D Roberts, M Woodage
!  5.4  25/04/02  Replace land/sea mask with land fraction in call
!                 to NUMBER_DROPLET.               A Jones
!  5.4  05/09/02  Add code to model coagulation of Aitken to
!                    accumulation mode SO4.       D Roberts, M Woodage
!  5.5  05/02/03  Calculation of n_droplet moved up to Aero_Ctl.
!                                                 P Davison.
!  6.1  17/08/04  NEC Opimisations    S.S.Wilson
!  6.2  21/12/05 Add code to output as diagnostics the sulphur fluxes
!               associated with the various processes.
!                   N. Bellouin and J. Rae.
!  6.2  21/11/05  make powers integers                   A.Malcolm
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
! System components covered:
!
! System task:
!
! Documentation: UM Doc 20
!
!---------------------------------------------------------------------
!
      IMPLICIT NONE
!
! CONTAINS GAS CONST R=287.05 J/KG/K
#include "c_r_cp.h"
! CONTAINS DENSITY OF WATER, RHO_WATER.
#include "c_densty.h"
! CONTAINS OTHER REQUIRED PARAMETERS
#include "c_sulchm.h"
#include "c_pi.h"
#include "c_g.h"
!
! Arguments with intent IN:
!
      INTEGER                                                           &
     & row_length                                                       &
                                  !no. of pts along a row
     &,rows                                                             &
                                  !no. of rows
     &,model_levels                                                     &
                                  !no. of model levels
     &,wet_model_levels                                                 &
                                  !no. of wet model levels
     &, halo_i                                                          &
                                  !EW halo size
     &, halo_j                                                          &
                                  !NS halo size
     &, off_x                                                           &
     &, off_y                                                           &
     &, theta_field_size
!
      REAL TSTEP                !chemistry tstep: LE physics tstep
!
      REAL                                                              &
     & CLOUDF(row_length,rows,wet_model_levels)                         &
                                                  !cloud fraction (0-1)
     &,QCL(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
     &                        wet_model_levels)                         &
                                                  !cloud liquid water
     &,QCF(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
     &                        wet_model_levels)                         &
                                                  !cloud frozen water
     &,Q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,              &
     &                        wet_model_levels)                         &
                                                  !specific humidity
     &,COSZA2D(row_length,rows)                                         &
                                                  !cos zenith angle
     &,P(1-off_x:row_length+off_x,                                      &
                                                  !pressure
     &                 1-off_y:rows+off_y, model_levels)                &
     &,T(row_length,rows,model_levels)                                  &
                                                  !temperature
     &,OH_CONC(row_length,rows,model_levels)                            &
                                                  !OH concn
     &,HO2_CONC(row_length,rows,model_levels)                           &
                                                  !HO2 concn
     &,H2O2_LMT(row_length,rows,model_levels)                           &
                                                  !H2O2 max limit
     &,O3(row_length,rows,model_levels)                                 &
                                                  !O3 mmr
     &,n_droplet(row_length, rows, wet_model_levels)
                                                  !Drop concentration
!
      LOGICAL                                                           &
     & L_SULPC_DMS                                                      &
                                         !T if DMS chemistry required
     &,L_SULPC_NEWDMS                                                   & 
                                         !T if new DMS scheme req'd,
!                                            else old scheme used
     &,L_SULPC_OZONE                                                    &
                                         !T if O3 oxidn required
     &,L_SULPC_NH3                       !T if NH3 buffering used
!
! Arguments with intent IN/OUT:
!
      REAL                                                              &
     & SO2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &                                                 model_levels)    &
     &,DMS(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &                                                 model_levels)    &
     &,SO4_AIT(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &                                                 model_levels)    &
     &,SO4_ACC(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &                                                 model_levels)    &
     &,SO4_DIS(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &                                                 model_levels)    &
     &,NH3(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &                                                 model_levels)    &
     &,H2O2_MXR(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &                                                 model_levels)
!
! Arguments with intent OUT (diagnostics):
!
      REAL                                                              &
     & MSA(row_length,rows,model_levels)                                &
                                                  !mmr S in MSA
     &,F_DMS_TO_SO2(row_length,rows,model_levels)                       &
                                                  !frac oxid DMS to SO2
     &,F_DMS_TO_SO4(row_length,rows,model_levels)                       &
                                                  !frac oxid DMS to SO4
     &,F_DMS_TO_MSA(row_length,rows,model_levels)                       &
                                                  !frac oxid DMS to MSA
     &,NH3_DEP(row_length,rows,model_levels)                            &
                                                  !NH3 depleted
     &,DELTAS_DRY(row_length,rows,model_levels)                         &
                                                   !SO2 dry ox per ts
     &,DELTAS_WET(row_length,rows,model_levels)                         &
                                                   !SO2 wet ox by H2O2
     &,DELTAS_WET_O3(row_length,rows,model_levels)                      &
                                                   !SO2 wet ox by O3
     &,DELTAS_TOT(row_length,rows,model_levels)                         &
                                                   !total SO2 ox per ts
     &,DELTAS_DMS(row_length,rows,model_levels)                         &
                                                   !DMS dry ox per ts
     &,DELTAS_EVAP(row_length,rows,model_levels)                        &
                                                   !SO4_DIS released by
!                             evapn of cloud droplets to SO4_ACC per ts
     &,DELTAS_NUCL(row_length,rows,model_levels)                        &
                                                   !SO4_ACC transfd by
!                                          nucleation to SO4_DIS per ts
     &,DELTAS_DIFFUSE(row_length,rows,model_levels)                     &
                                                    !SO4_AIT transfd to
!                                           SO4_DIS by diffusion per ts
     &,DELTAS_COAG(row_length,rows,model_levels)                        &
                                                 !SO4_AIT transfd by
!                                        coagulation to SO4_ACC per ts
     &,DELTAS_TOT_AIT(row_length,rows,model_levels)                     &
                                                    !total SO4_AIT
!                                                       transfd per ts
     &,PSI(row_length,rows,model_levels)           !fraction of dry ox
!                                                   SO2 becoming SO4_AIT
!
! Local variables:
!
      INTEGER i,j,k               !loop variables
!
      REAL                                                              &
     &     DRYRATE                                                      &
                                  ! rate of dry oxidn SO2 in 1/s
     &,    WETRATE                                                      &
                                  ! rate of wet oxidn SO2 in 1/s
     &,    DMSRATE                ! rate of dry oxidn DMS in 1/s
!
      REAL                                                              &
     & CLEARF(row_length,rows,wet_model_levels)                         &
                                                   !clear air fraction
     &,QCTOTAL(row_length,rows,wet_model_levels)                        &
                                                   !tot condensed water
     &,RHO_AIR(row_length,rows,model_levels)                            &
                                                   !air density kg/m3
     &,AIR_CONC(row_length,rows,model_levels)                           &
                                                   !air concn  mcls/cc
     &,VISCOS_AIR(row_length,rows,model_levels)                         &
                                                !air viscosity (kg/m/s)
     &,MEAN_FREEP(row_length,rows,model_levels)                         &
                                                !mean free path of air
!                                                   molecules (m)
     &,REL_HUM(row_length,rows,model_levels)                            &
                                                !relative humidity(0-1)
     &,P_ARRAY(row_length,rows,model_levels)                            &
                                                !pressure, no halo
     &,SO2_ARRAY(row_length,rows,model_levels)                          &
                                                !SO2, no halo
     &,SO4_ACC_ARRAY(row_length,rows,model_levels)                      &
                                                    !SO4_ACC,no halo
     &,SO4_AIT_ARRAY(row_length,rows,model_levels)                      &
                                                    !SO4_AIT,no halo
     &,SO4_DIS_ARRAY(row_length,rows,model_levels)                      &
                                                    !SO4_DIS,no halo
     &,QCL_ARRAY(row_length,rows,wet_model_levels)                      &
                                                    !QCL, no halo
     &,QCF_ARRAY(row_length,rows,wet_model_levels)                      &
                                                    !QCF, no halo
     &,Q_ARRAY(row_length,rows,wet_model_levels)                        &
                                                  !Q, no halo
     &,H2O2_MXR_ARRAY(row_length,rows,model_levels)                     &
                                                    !H2O2_MXR,no halo
     &,NH3_ARRAY(row_length,rows,model_levels)                          &
                                                !NH3, no halo
     &,DMS_ARRAY(row_length,rows,model_levels)  !DMS, no halo

      REAL                                                              &
     &     EVAPTIME,                                                    &
                           ! timescale for cloud droplets to evaporate
     &     NUCTIME,                                                     &
                           ! timescale for particles to enter a cloud
!                            and nucleate.
     &     DIFFUSE_TAU,                                                 &
                           ! diffusive lifetime of Aitken mode particles
!                            once they enter a cloud
     &     RHO_CUBEROOT,                                                &
                           ! cube root of density of water.
     &     DIFFUSE_TAURATIO,                                            &
                             ! CLOUDTAU/DIFFUSE_TAU
     &     PROBDIFF_INV,                                                &
                           ! inverse of probability of a particle being
!                            diffusionallly scavenged in a cloud.
     &     PROBDIFF_FN1,                                                &
                           ! PROBDIFF_INV - 0.5
     &     PROBDIFF_FN2,                                                &
                           ! PROBDIFF_INV*EXP(DIFFUSE_TAURATIO*0.5)
     &     PROBDIFF_CLOUD,                                              &
                           ! probability of an Aitken mode particle
!                            being in cloud at the start of a step.
     &     PROBDIFF_CLEAR,                                              &
                           ! probability of an Aitken mode particle
!                            being in clear air at the start of a step.
     &     LAMBDA_AITKEN,                                               &
                           ! ratio of concentrations of Aitken mode
!                            particles in cloudy to clear air.
     &     DIFFRATE ! rate of diffusive capture of Aitken mode particles
!
      REAL                                                              &
     &     W_CONC,                                                      &
                                  ! H2O concentration in molecules/cc
     &     K_SO2_OH,                                                    &
                                  ! reaction rate for SO2+OH  cc/mcl/s
     &     K_SO2OH_LO,                                                  &
                                  ! low_pressure reaction rate limit
     &     K_LIMRAT,                                                    &
                                  ! ratio  reaction rate limits (LO/HI)
     &     K_HO2_HO2,                                                   &
                                  ! reaction rate for HO2 to H2O2
     &     H2O2_CON_TO_MXR,                                             &
                            ! factor = RMM_H2O2/AVOGADRO/AIR DENSITY
     &     O2_CONC,                                                     &
                            ! O2 concentration in molecules/cc
     &     O3_CONC,                                                     &
                            ! O3 concentration in molecules/cc
     &     K1_DMS_OH,                                                   &
                            ! reaction rate for DMS+OH via path 1
     &     K2_DMS_OH,                                                   &
                            ! reaction rate for DMS+OH via path 2
     &     K3_CH3SO2,                                                   &
                            ! thermal decomposition rate for CH3SO2
     &     K6_CH3SO3,                                                   &
                            ! thermal decomposition rate for CH3SO3
     &     F_DMS_TO_CH3SO2,                                             &
                            ! fraction of oxidised DMS becoming CH3SO2
     &     F_DMS_TO_CH3SO3,                                             &
                            ! fraction of oxidised DMS becoming CH3SO3
     &     F_CH3SO2_TO_SO2,                                             &
                            ! fraction of CH3SO2 becoming SO2
     &     F_CH3SO3_TO_SO4  ! fraction of CH3SO3 becoming SO4
!
      REAL FAC2,FAC3              !for interp between LO and HI K_SO2OH
!                                 !    reaction rate limits
!
      REAL TERM1,                                                       &
                                  ! dummy variables to assist
     &     TERM2,                                                       &
                                  !  wet oxidn calculation
     &     DENOM,                                                       &
                                  !
     &     TERM3,                                                       &
                                  ! dummy variables to assist
     &     TERM4,                                                       &
                                  ! calculation of diffusive capture.
     &     TAURATIO,                                                    &
                                  ! CLOUDTAU/CHEMTAU
     &     HFTAURAT,                                                    &
                                  ! CLOUDTAU/2*CHEMTAU
     &     PROBOX_INV,                                                  &
                                  ! 1/probability of oxidn in cloud
     &     PROBOX_FN1,                                                  &
                                  ! PROBOX_INV - 0.5
     &     PROBOX_FN2,                                                  &
                                  ! PROBOX_INV * EXP(-HFTAURAT)
     &     PROB_CLOUD,                                                  &
                                  ! probability SO2 starts in cloud
     &     PROB_CLEAR,                                                  &
                                  ! probability SO2 starts in clear air
     &     LAMDA,                                                       &
                                  ! ratio SO2 in cloud/clear air
     &     H2O2_REP               ! replenished H2O2
!
!     Extra variables for improved diffusional scavenging.
!
      REAL DIFFUSIVITY   !Mean diffusion coefficent of Aitken particles
      REAL VISCOSITY_AIR !Dynamic viscosity of air (kg m-1 s-1).
      REAL MEAN_FREE_PATH  !Mean free path of air molecules.
      REAL KNUDSEN_WEBER   !Expression related to the Cunningham slip
!                          flow correction to the diffusion coefficient
      REAL LOG_SIGMA_AIT    !Natural log of SIGMA_AIT (in C_SULCHM).
      REAL SQ_LOG_SIGMA_AIT !Square of the previous parameter.
      REAL DIFF_CON1        !Term in diffusion coefficent formula.
      REAL DIFF_CON2        !Term in diffusion coefficent formula.
      REAL PEC          !Quantity associated with (not =) Peclet number
      REAL WORK_RADIUS  !Variable linked to average droplet radius.
      REAL SCAVCD       !Scavenging coefficient for in-cloud
!                  advective-diffusive removal of Aitken mode particles
!
! Extra variables for coagulation code
!
      REAL COAGRATE          !Rate of transfer of SO4_AIT to SO4_ACC
!                                    due to coagulation (kg kg-1 s-1)
      REAL                                                              &
     &     ALPHA                                                        &
                             !Hygroscopic growth parameter in
                             ! Fitzgerald's scheme (assume BETA=1)
     &,    A                                                            &
                             !Parameters required for COAG_QI function
     &,    B                                                            &
                             !          "
     &,    SUM_IY                                                       &
                             !Sum of COAG_QI_ functions in coag calc
     &,    SUM_IZ                                                       &
                             !              "
     &,    YFAC                                                         &
                             !Coefficient used in coag calcn
     &,    ZFAC              !         "
!
! COAG_QI_ functions
      REAL                                                              &
     &     COAG_QI_1                                                    &
                             !COAG_QI for A=1,    B=0
     &,    COAG_QI_2                                                    &
                             !   "        A=4/3,  B=-1/3
     &,    COAG_QI_3                                                    &
                             !   "        A=2/3,  B=1/3
     &,    COAG_QI_4                                                    &
                             !   "        A=2/3,  B=0
     &,    COAG_QI_5                                                    &
                             !   "        A=4/3,  B=-2/3
     &,    COAG_QI_6                                                    &
                             !   "        A=1/3,  B=1/3
     &,    COAG_QI_7         !   "        A=1,    B=-1/3
!
      REAL ACUNN                                                        &
                             !Cunningham correction factor
     &,    LOG_SIGMA_ACC                                                &
                             !nat log of geom std dev of acc mode distn
     &,    SQ_LOG_SIGMA_ACC  !square of   "
!
      PARAMETER (ACUNN = 1.591)
!
!
      REAL MMR_STAR        ! threshold mmr of SO4_ACC below which PSI=1
      REAL CON_1,                                                       &
                           ! constants required for PSI calcn
     &     CON_2
      REAL SCALE_FACTOR                                                 &
                                 ! to prevent negative SO2
     &,    SCALE_FACT_AIT           !to prevent negative SO4_AIT
!
! External routines:
!
      EXTERNAL QSAT_WAT                                                 &
     &,        HYGRO_FACT
      EXTERNAL COAG_QI
!
      REAL COAG_QI           !Function for coagulation calcn
!
!--------------------------------------------------------------------
! 0. Set up constants and initialise OUT arrays to 0
!--------------------------------------------------------------------
!
! The next constant cannot be declared as PARAMETERs because
! they involve non-integer exponents.
!
      RHO_CUBEROOT = RHO_WATER**0.333333
!
!     Extra parameters for improved diffusional scavenging.
!
      LOG_SIGMA_AIT = LOG(SIGMA_AIT)
      SQ_LOG_SIGMA_AIT = LOG_SIGMA_AIT * LOG_SIGMA_AIT
      DIFF_CON1 = EXP(-2.5*SQ_LOG_SIGMA_AIT)/RAD_AIT
      DIFF_CON2 = EXP(-4.0*SQ_LOG_SIGMA_AIT)/(RAD_AIT*RAD_AIT)
!
!
!     Extra parameters for coagulation
!
      LOG_SIGMA_ACC=LOG(SIGMA)
      SQ_LOG_SIGMA_ACC=LOG_SIGMA_ACC*LOG_SIGMA_ACC
!
! Calculate COAG_QI_ functions (not dependent on ALPHA)
!
! DEPENDS ON: coag_qi
      COAG_QI_1=COAG_QI(1.0, 0.0,                                       &
     &             RAD_AIT,RAD_ACC,SQ_LOG_SIGMA_AIT,SQ_LOG_SIGMA_ACC)
!
! DEPENDS ON: coag_qi
      COAG_QI_2=COAG_QI(1.333333, -0.333333,                            &
     &             RAD_AIT,RAD_ACC,SQ_LOG_SIGMA_AIT,SQ_LOG_SIGMA_ACC)

! DEPENDS ON: coag_qi
      COAG_QI_3=COAG_QI(0.666667, 0.333333,                             &
     &             RAD_AIT,RAD_ACC,SQ_LOG_SIGMA_AIT,SQ_LOG_SIGMA_ACC)
!
! DEPENDS ON: coag_qi
      COAG_QI_4=COAG_QI(0.666667, 0.0,                                  &
     &             RAD_AIT,RAD_ACC,SQ_LOG_SIGMA_AIT,SQ_LOG_SIGMA_ACC)
!
! DEPENDS ON: coag_qi
      COAG_QI_5=COAG_QI(1.333333, -0.666667,                            &
     &             RAD_AIT,RAD_ACC,SQ_LOG_SIGMA_AIT,SQ_LOG_SIGMA_ACC)
!
! DEPENDS ON: coag_qi
      COAG_QI_6=COAG_QI(0.333333, 0.333333,                             &
     &             RAD_AIT,RAD_ACC,SQ_LOG_SIGMA_AIT,SQ_LOG_SIGMA_ACC)
!
! DEPENDS ON: coag_qi
      COAG_QI_7=COAG_QI(1.0, -0.333333,                                 &
     &             RAD_AIT,RAD_ACC,SQ_LOG_SIGMA_AIT,SQ_LOG_SIGMA_ACC)
!
! Calculate sums of COAG_QI_ functions
!
      SUM_IZ = 2.0*COAG_QI_1 + COAG_QI_2 + COAG_QI_3
!
      SUM_IY = COAG_QI_4 + COAG_QI_5 + COAG_QI_6 + COAG_QI_7
!
! Calculate parameters which are same for each grid box (for wet oxidn)
! and set up other arrays
!
      TAURATIO=CLOUDTAU/CHEMTAU
      HFTAURAT=TAURATIO/2.0
      PROBOX_INV=1.0/(1.0-EXP(-TAURATIO))
      PROBOX_FN1=PROBOX_INV-0.5
      PROBOX_FN2=PROBOX_INV*EXP(-HFTAURAT)
!
      P_ARRAY(:,:,:)=P(1:row_length,1:rows,1:model_levels)
      SO2_ARRAY(:,:,:)=SO2(1:row_length,1:rows,1:model_levels)
      SO4_ACC_ARRAY(:,:,:)=SO4_ACC(1:row_length,1:rows,1:model_levels)
      SO4_AIT_ARRAY(:,:,:)=SO4_AIT(1:row_length,1:rows,1:model_levels)
!cdir collapse
      Do k=1,model_levels
        Do j=1,rows
          Do i=1,row_length
            MSA(i,j,k)  = 0.0
            F_DMS_TO_SO2(i,j,k)  = 0.0
            F_DMS_TO_SO4(i,j,k)  = 0.0
            F_DMS_TO_MSA(i,j,k)  = 0.0
            NH3_DEP(i,j,k)  = 0.0
            DELTAS_DRY(i,j,k)  = 0.0
            DELTAS_WET(i,j,k)  = 0.0
            DELTAS_WET_O3(i,j,k) = 0.0
            DELTAS_TOT(i,j,k)  = 0.0
            DELTAS_DMS(i,j,k)  = 0.0
            DELTAS_EVAP(i,j,k) = 0.0
            DELTAS_NUCL(i,j,k) = 0.0
            DELTAS_DIFFUSE(i,j,k) = 0.0
           DELTAS_COAG(i,j,k) = 0.0
           DELTAS_TOT_AIT(i,j,k) = 0.0
!
        RHO_AIR(i,j,k)=P_ARRAY(i,j,k)/(R*T(i,j,k))
        AIR_CONC(i,j,k)=                                                &
     &       (AVOGADRO*P_ARRAY(i,j,k))/(R*T(i,j,k)*RMM_AIR*1.0E6)
!
! Compute dynamic viscosity of air and mean free path of air molecules
! using the formulae on P.417 of Pruppacher and Klett, 2nd edition.
!
      VISCOS_AIR(i,j,k)=(1.718 + (T(i,j,k)-273.15)*0.0049 )*1.0E-5
      MEAN_FREEP(i,j,k)=                                                &
     &       (MFP_REF*PREF_MFP*T(i,j,k))/(TREF_MFP*P_ARRAY(i,j,k))
!
          End Do
        End Do
      End Do
!
!---------------------------------------------------------------------
! 1.  This section calculates amount of SO2 dry oxidised using a
!     3_D OH concentration field to control the oxidn rate when
!     cos(zenith angle) is G.T. 10-6  (else rate is zero).
!  Reaction rates are dependent on temperature and pressure
!---------------------------------------------------------------------
!
      CON_1 = EXP(4.5*LOG(SIGMA)*LOG(SIGMA))
      CON_2 = 4.0*PI*RHO_SO4*CHI*NUM_STAR*(RAD_ACC)**3
!
      Do k=1,model_levels
!cdir collapse
        Do j=1,rows
          Do i=1,row_length
!
          K_SO2OH_LO=AIR_CONC(i,j,k)*K1*(T1/T(i,j,k))**PARH
!
          K_LIMRAT = K_SO2OH_LO / K_SO2OH_HI
              FAC2 = K_SO2OH_LO /(1.0+ K_LIMRAT)
              FAC3 = 1.0 + (LOG10(K_LIMRAT)/FAC1)**2
!
          K_SO2_OH = FAC2 * FC**(1.0/FAC3)         ! Variable K_SO2_OH
!
!  Only calculate the dry oxidation rate in the daytime.
!
         IF (COSZA2D(i,j) >  1.0E-06) THEN
           DRYRATE=K_SO2_OH * OH_CONC(i,j,k)
         ELSE                     ! Zero rate if nightime
           DRYRATE=0.0
         ENDIF
!
         DELTAS_DRY(i,j,k)=DRYRATE*TSTEP*SO2_ARRAY(i,j,k) !AmntSO2dryox
!
!  Calculate fraction becoming Aitken mode SO4
      MMR_STAR = CON_1*CON_2/(3.0*RHO_AIR(i,j,k))
!
      IF (SO4_ACC_ARRAY(i,j,k) <  MMR_STAR)   THEN
        PSI(i,j,k)=1.0
      ELSE
        PSI(i,j,k)=RAD_ACC*E_PARM*SO4_AIT_ARRAY(i,j,k)/                 &
     &          (RAD_ACC*E_PARM*SO4_AIT_ARRAY(i,j,k)+                   &
     &           RAD_AIT*SO4_ACC_ARRAY(i,j,k))
      END IF
!
          End Do
        End Do
      End Do
!
!--------------------------------------------------------------
! 2. This section calculates amount of SO2 wet oxidised using a
!    3_D H2O2 field to control the oxidn rate. As H2O2 is depleted, it
!    is replenished from the HO2 concn field (via a reaction rate), but
!    it may not exceed the H2O2 LIMIT field .
!    The reaction rate is a function of pressure, temp and humidity.
!    It uses D.Roberts' formula for the wet oxidn rate, and R.Colvile's
!    H2O2 and HO2 chemistry.
!    Depletion and replenishment of H2O2 is done at the end of the
!    routine, after scaling to avoid SO2 becoming negative.
!
!    If H2O2 is limiting, further oxidation by O3 occurs if there is
!    sufficient NH3 available to buffer the reaction; Choularton's
!    suggested procedure is:
!    a) Do oxidn of SO2 by H2O2
!    b) Neutralise sulphate produced with NH3
!    c) If SO2 remains do oxidn by O3 until SO2 or NH3 required to
!       neutralise it is exhausted.
!
!--------------------------------------------------------------
!
      QCL_ARRAY(:,:,:)=QCL(1:row_length,1:rows,1:wet_model_levels)
      QCF_ARRAY(:,:,:)=QCF(1:row_length,1:rows,1:wet_model_levels)
      H2O2_MXR_ARRAY(:,:,:)=H2O2_MXR(1:row_length,1:rows,1:model_levels)
      NH3_ARRAY(:,:,:)=NH3(1:row_length,1:rows,1:model_levels)
      Do k=1,wet_model_levels
!cdir collapse
        Do j=1,rows
          Do i=1,row_length
!
! Calculate clear air fraction
          CLEARF(i,j,k)=1.0-CLOUDF(i,j,k)
!
      IF ((QCL_ARRAY(i,j,k) >= THOLD) .AND. (CLOUDF(i,j,k) >  0.0)) THEN
!
!   CALCULATE LAMDA
           TERM1=(CLEARF(i,j,k)*TAURATIO)**2
           TERM1=TERM1+2.0*TAURATIO*CLEARF(i,j,k)*                      &
     &                             (CLEARF(i,j,k)-CLOUDF(i,j,k))
          TERM1=SQRT(1.0+TERM1)
           TERM2=2.0*CLOUDF(i,j,k)-TAURATIO*CLEARF(i,j,k)-1.0
          TERM2=TERM2+TERM1
           LAMDA=TERM2/(2.0*CLOUDF(i,j,k))
!
!   CALCULATE PROB_CLEAR AND PROB_CLOUD
           DENOM=CLEARF(i,j,k)+CLOUDF(i,j,k)*LAMDA
           PROB_CLEAR=CLEARF(i,j,k)/DENOM
           PROB_CLOUD=CLOUDF(i,j,k)*LAMDA/DENOM
!
!   CALCULATE EXPECTED LIFETIME OF SO2 (WHICH IS 1/WETRATE)
           TERM1=PROBOX_FN1*PROB_CLEAR
           TERM2=PROBOX_FN2*PROB_CLOUD
           TERM2=(TERM1+TERM2)*CLEARF(i,j,k)/CLOUDF(i,j,k)
        DENOM=TERM2*CLOUDTAU+CHEMTAU
        WETRATE=1.0/DENOM
!
!!  RC'S H2O2 CHEMISTRY: OXIDATION OF SO2 TO SO4 IS CONTROLLED BY THE
!                        SMALLER OF THE SO2 AND H2O2 MIX RAT FIELDS
!
       IF (SO2_ARRAY(i,j,k) <= (H2O2_MXR_ARRAY(i,j,k)*RELM_S_H2O2)) THEN
!
         DELTAS_WET(i,j,k)=(1.0-EXP(-WETRATE*TSTEP))*SO2_ARRAY(i,j,k)
!
       ELSE
!
       DELTAS_WET(i,j,k)=(1.0-EXP(-WETRATE*TSTEP))*                     &
     &                    H2O2_MXR_ARRAY(i,j,k)*RELM_S_H2O2
!
!   H2O2 field is limiting, so oxidise SO2 further with O3, using NH3
!   as buffer, if sufficient O3 and NH3 available.
!
         If (L_SULPC_OZONE) Then
 !
         IF ( (O3(i,j,k) >= O3_MIN) .AND.                               &
     &    (NH3_ARRAY(i,j,k)  >   DELTAS_WET(i,j,k)/RELM_S_2N) ) THEN
!
! O3 oxidation is controlled by smaller of SO2 and NH3 fields.
! 2 atoms of N required for each S atom in ammonium sulphate.
! Sufficient NH3 must be available to neutralise all sulphate produced
! in grid box for O3 reaction to continue (otherwise PH too low)
! Note that all SO2 or NH3 is used in DELTAS_WET_O3 calcn; adjustment by
! SCALE_FACTOR at end of routine prevents removing too much.

!SO2 limiting
          IF (SO2_ARRAY(i,j,k) <= (NH3_ARRAY(i,j,k)*RELM_S_2N)) THEN
!
            DELTAS_WET_O3(i,j,k)=                                       &
     &             (1.0-EXP(-WETRATE*TSTEP))*SO2_ARRAY(i,j,k)
!
          ELSE                                        ! NH3 limiting
!
           DELTAS_WET_O3(i,j,k)=(1.0-EXP(-WETRATE*TSTEP))*              &
     &                                   NH3_ARRAY(i,j,k)*RELM_S_2N
!
          END IF
!
         END IF               ! End ozone oxidation block
!
         End If               ! End L_SULPC_OZONE Test
!
       END IF                 ! End peroxide oxidation block
!
       END IF                 ! End cloud present condition
!
!
          End Do
        End Do
      End Do
!
!------------------------------------------------------------------
! 3.  This section calculates amount of DMS dry oxidised to SO2
!     and MSA using a 3D OH concentration field to control the oxidn
!     rate when cos zenith angle is G.T. 10-6  (else rate is zero).
!     MSA is not saved, and K_DMS_OH is constant in this version.
!    OR, if O3 available and new DMS chemistry required -
!     This section calculates the amount of DMS dry oxidised to SO2,
!     SO4 and MSA using OH, HO2, O3 and O2 as oxidants. OH and HO2
!     are available in daylight only. Intermediate species have
!     short lifetimes so are not stored between timesteps.
!------------------------------------------------------------------
!
      IF (L_SULPC_DMS) THEN
!
        DMS_ARRAY(:,:,:)=DMS(1:row_length,1:rows,1:model_levels)

        If (L_SULPC_NEWDMS .AND. L_SULPC_OZONE) Then
!
      Do k=1,model_levels
!cdir collapse
        Do j=1,rows
          Do i=1,row_length
!
! Calculate rate coefficients K1_DMS_OH and K2_DMS_OH for the two
!    DMS+OH -> CH3SO2 reactions, which are daylight dependent
!
          IF (COSZA2D(i,j) >  1.0E-06) THEN     ! daylight
!
          K1_DMS_OH=1.13E-11*EXP(-254.0/T(i,j,k))
!
!  Calculate concentration of O2 as 20.95% of air concentration

          O2_CONC=0.2095*AIR_CONC(i,j,k)
!
          K2_DMS_OH=     ( 1.7E-42*O2_CONC*EXP(7810.0/T(i,j,k)) ) /     &
     &         ( 1.0 + ( 5.5E-31*O2_CONC*EXP(7460.0/T(i,j,k)) ) )
!
!  Calculate thermal decomposition rate coefficient K3_CH3SO2 for
!         CH3SO2 -> CH3+SO2
!
          K3_CH3SO2=2.6E11*EXP(-9056.0/T(i,j,k))
!
!  Rate coefficient K4_CH3SO2_O3 for CH3SO2+O3 -> CH3SO3+O2 is const
!         K4_CH3SO2_O3=1.0E-14  in C_SULCHM
!
!  Rate coefficient K5_CH3SO3_HO2 for CH3SO3+HO2 -> MSA+O2 is const
!         K5_CH3SO3_HO2=4.0E-11  in C_SULCHM
!
!  Calculate thermal decomposition rate coefficient K6_CH3SO3 for
!         CH3SO3 -> CH3+SO3
!
          K6_CH3SO3=1.1E17*EXP(-12057.0/T(i,j,k))
!
!  Calculate DMS oxidation rate
!
          DMSRATE = (K1_DMS_OH + K2_DMS_OH) * OH_CONC(i,j,k)
!
!  Calculate fraction of DELTAS_DMS becoming CH3SO2
          F_DMS_TO_CH3SO2=(K1_DMS_OH+0.9*K2_DMS_OH)/                    &
     &                    (K1_DMS_OH+K2_DMS_OH)
!
!  Calculate fraction of DELTAS_DMS becoming SO2 (requires O3_CONC)
          O3_CONC= O3(i,j,k) * AIR_CONC(i,j,k) * RMM_AIR/RMM_O3
          F_CH3SO2_TO_SO2=K3_CH3SO2/(K3_CH3SO2+K4_CH3SO2_O3*O3_CONC)
          F_DMS_TO_SO2(i,j,k)=F_DMS_TO_CH3SO2 * F_CH3SO2_TO_SO2
!
!  Calculate fraction of DELTAS_DMS becoming sulphate
          F_CH3SO3_TO_SO4= K6_CH3SO3 /                                  &
     &                    (K6_CH3SO3+K5_CH3SO3_HO2*HO2_CONC(i,j,k))
          F_DMS_TO_SO4(i,j,k)=F_DMS_TO_CH3SO2 * (1.0-F_CH3SO2_TO_SO2)   &
     &                    * F_CH3SO3_TO_SO4
!
!  Calculate fraction of DELTAS_DMS becoming MSA
          F_DMS_TO_MSA(i,j,k)=1.0 - F_DMS_TO_SO2(i,j,k)                 &
     &                            - F_DMS_TO_SO4(i,j,k)
!
!  Calculate amount of DMS oxidised in TSTEP (in daylight only)
!
          DELTAS_DMS(i,j,k)=DMSRATE*TSTEP*DMS_ARRAY(i,j,k) !DELTA_DMS
!
          END IF             ! End daylight condition
!
          End Do
        End Do
      End Do
!
        Else                 ! implement simple DMS scheme
!
      Do k=1,model_levels
!cdir collapse
        Do j=1,rows
          Do i=1,row_length
!
!  Calculate DMS oxidation rate if daytime.
          IF (COSZA2D(i,j) >  1.0E-06) THEN
            DMSRATE = K_DMS_OH * OH_CONC(i,j,k)
          ELSE                 ! ZERO RATE IF NIGHTIME
            DMSRATE=0.0
          ENDIF                   ! END COSZA2D IF
!
!  CALCULATE FRACTION OF DMS OXIDISED
!
          DELTAS_DMS(i,j,k)=DMSRATE*TSTEP*DMS_ARRAY(i,j,k)!AmtDMSdryox
!
          End Do
        End Do
      End Do
!
        End If                 ! End L_SULPC_NEWDMS Test
!
      END IF                   ! END L_SULPC_DMS IF

!---------------------------------------------------------------------
! 4. Release of aerosol from evaporating cloud droplets:
!    if no condensed water (liquid + ice) in grid box, release
!    dissolved sulphate as accumulation mode aerosol.
!     If cloud fraction less than 0.95, release some in clear  air.
!--------------------------------------------------------------------
!
      SO4_DIS_ARRAY(:,:,:)=SO4_DIS(1:row_length,1:rows,1:model_levels)
!cdir collapse
      Do k=1,wet_model_levels
        Do j=1,rows
          Do i=1,row_length
!
            QCTOTAL(i,j,k)= QCL_ARRAY(i,j,k) + QCF_ARRAY(i,j,k)
          IF (QCTOTAL(i,j,k)  <   THOLD) THEN
            DELTAS_EVAP(i,j,k)=SO4_DIS_ARRAY(i,j,k)
          ELSE IF  (CLOUDF(i,j,k) <  0.95) THEN
            EVAPTIME = EVAPTAU + 0.5*CLOUDTAU
            DELTAS_EVAP(i,j,k)=                                         &
     &             (1.0-EXP(-TSTEP/EVAPTIME))*SO4_DIS_ARRAY(i,j,k)
          ELSE
            DELTAS_EVAP(i,j,k)=0.0
          ENDIF
!
          End Do
        End Do
      End Do
!
!     Also release any dissolved aerosol in a non-wet level,
!     where it should not be.
!
      IF (wet_model_levels  <   model_levels)  THEN
!
        Do k=wet_model_levels+1,model_levels
!cdir collapse
          Do j=1,rows
            Do i=1,row_length
              DELTAS_EVAP(i,j,k)= SO4_DIS_ARRAY(i,j,k)
            End Do
          End Do
        End Do
!
      ENDIF
!
!-------------------------------------------------------------------
! 5. Nucleation of accumulation mode aerosol forming dissolved SO4
!-------------------------------------------------------------------
!
!    THIS CODE ASSUMES THAT THE PARAMETER NUCTAU, WHICH IS THE
!    TIMESCALE FOR NUCLEATION ONCE A PARTICLE ENTERS A CLOUD, IS
!    VERY SHORT COMPARED WITH CLOUDTAU.
!
!cdir collapse
      Do k=1,wet_model_levels
        Do j=1,rows
          Do i=1,row_length
!
          IF ((QCTOTAL(i,j,k) >= THOLD) .AND.                           &
     &                                (CLOUDF(i,j,k) >  0.0))  THEN
          NUCTIME=NUCTAU+((CLEARF(i,j,k)*CLOUDTAU)/(2.0*CLOUDF(i,j,k)))
          DELTAS_NUCL(i,j,k)=                                           &
     &       (1.0-EXP(-TSTEP/NUCTIME))*SO4_ACC_ARRAY(i,j,k)
          ENDIF
!
          End Do
        End Do
      End Do
!
!-------------------------------------------------------------------
! 6. Diffusional scavenging of Aitken mode SO4 forming dissolved SO4
! Improved code (introduced at vn5.3) allows for:
! (1) The enhancement of diffusional scavenging due to the
! relative motion of cloud droplets and aerosol is modelled.
! (2) The variation of particle diffusivity with ambient
! temperature and pressure is taken into account.
! (3) Averaging over the size distributions of both the cloud
! droplets and the aerosol particles has been done, assuming
! a Khrgian-Mazin distribution for the droplets and a lognormal
! distribution for the aerosol.
!
!-------------------------------------------------------------------
!
!    THIS IS A MUCH SLOWER PROCESS THAN NUCLEATION AND THEREFORE WE
!    CANNOT MAKE THE SAME APPROXIMATIONS USED IN THAT CASE.
!
!cdir collapse
      Do k=1,wet_model_levels
        Do j=1,rows
          Do i=1,row_length
!
          IF ((QCTOTAL(i,j,k) >= THOLD) .AND.                           &
     &                                 (CLOUDF(i,j,k) >  0.0)) THEN
!
!    FIRST COMPUTE IN-CLOUD TIMESCALE FOR DIFFUSIONAL CAPTURE,
!    USING TOTAL CONDENSED WATER WITHIN THE CLOUDY PORTION OF THE BOX.
!    THE DIFFERENCE BETWEEN LIQUID WATER AND ICE IS NEGLECTED HERE.
!    THIS SHOULD BE IMPROVED ON EVENTUALLY.
!
!
! Compute mean free path of air molecules.
!  (See P.417 of Pruppacher and Klett, 2nd edition.)
!
      MEAN_FREE_PATH=                                                   &
     &       (MFP_REF*PREF_MFP*T(i,j,k))/(TREF_MFP*P_ARRAY(i,j,k))
!
! Compute the Knudsen and Weber term for a particle of the median
!  radius (note approximation here: we do not average over the size
!  distribution). See P.450 of Pruppacher and Klett, 2nd edition.
!
      KNUDSEN_WEBER=1.257 + 0.4*EXP(-((1.1*RAD_AIT)/MEAN_FREE_PATH))
!
! Temporarily use DIFFUSIVITY to store working value.
!
      DIFFUSIVITY=DIFF_CON1 + DIFF_CON2*MEAN_FREE_PATH*KNUDSEN_WEBER
!
! Compute dynamic viscosity of air, using an approximate version of
!  the formula on P.417 of Pruppacher and Klett, 2nd edition.
!
      VISCOSITY_AIR=(1.718 + (T(i,j,k)-273.15)*0.0049 )*1.0E-5
!
! Now compute mean diffusion coefficient.
!
      DIFFUSIVITY=(BOLTZMANN*T(i,j,k)*DIFFUSIVITY)/                     &
     &                               (6.0*PI*VISCOSITY_AIR)
!
! Now compute the term PEC related to (but not equal to) the cube
!  root of the Peclet Number.
!
      PEC=((4.0*G*RHO_WATER)/(9.0*DIFFUSIVITY*VISCOSITY_AIR))**0.333333
!
      WORK_RADIUS=QCTOTAL(i,j,k)/                                       &
     &            (CLOUDF(i,j,k)*10.0*PI*RHO_WATER*N_DROPLET(i,j,k))
      WORK_RADIUS=WORK_RADIUS**0.333333
!
! We can finally compute the timescale for diffusive
!  scavenging once inside cloud, DIFFUSE_TAU.
!
      SCAVCD=6.0*PI*DIFFUSIVITY*N_DROPLET(i,j,k)*WORK_RADIUS*           &
     &         (1.0 + PEC*WORK_RADIUS)
!
      DIFFUSE_TAU=1.0/SCAVCD
!
!
            DIFFUSE_TAURATIO = CLOUDTAU/DIFFUSE_TAU
            PROBDIFF_INV = 1.0/( 1.0 - EXP(-DIFFUSE_TAURATIO) )
            PROBDIFF_FN1 = PROBDIFF_INV - 0.5
            PROBDIFF_FN2 = PROBDIFF_INV*EXP(-(0.5*DIFFUSE_TAURATIO) )
!
!     CALCULATE LAMBDA_AITKEN.
!
            TERM3 = (CLEARF(i,j,k)*DIFFUSE_TAURATIO)**2
            TERM3 = TERM3 + ( 2.0*DIFFUSE_TAURATIO                      &
     &                 *CLEARF(i,j,k)*(CLEARF(i,j,k)-CLOUDF(i,j,k)) )
            TERM3 = SQRT(1.0+TERM3)
            TERM4=2.0*CLOUDF(i,j,k)-DIFFUSE_TAURATIO*CLEARF(i,j,k)-1.0
            TERM4 = TERM4 + TERM3
            LAMBDA_AITKEN = TERM4/(2.0*CLOUDF(i,j,k))
!
!   CALCULATE PROBDIFF_CLEAR AND PROBDIFF_CLOUD
!
            DENOM = CLEARF(i,j,k)+CLOUDF(i,j,k)*LAMBDA_AITKEN
            PROBDIFF_CLEAR = CLEARF(i,j,k)/DENOM
            PROBDIFF_CLOUD = (CLOUDF(i,j,k)*LAMBDA_AITKEN)/DENOM
!
!   CALCULATE EXPECTED LIFETIME OF AN AITKEN MODE PARTICLE WITH
!   RESPECT TO DIFFUSIVE CAPTURE BY CLOUD DROPLETS.
!
            TERM3 = PROBDIFF_FN1*PROBDIFF_CLEAR
            TERM4 = PROBDIFF_FN2*PROBDIFF_CLOUD
            TERM4 = (TERM3+TERM4)*CLEARF(i,j,k)/CLOUDF(i,j,k)
            DENOM = TERM4*CLOUDTAU + DIFFUSE_TAU
            DIFFRATE = 1.0/DENOM
!
!   NOW COMPUTE THE AMOUNT OF AITKEN MODE SO4 CONVERTED TO DISSOLVED
!   SO4 BY DIFFUSIVE CAPTURE DURING THE STEP.
!
            DELTAS_DIFFUSE(i,j,k)=(1.0-EXP(-DIFFRATE*TSTEP))*           &
     &                                    SO4_AIT_ARRAY(i,j,k)
            ENDIF
         End Do
        End Do
      End Do
!
!---------------------------------------------------------------------
! 6.5. Coagulation of Aitken mode with accumulation mode SO4 to form
!      more of the accumulation mode. Allowance is made for hygroscopic
!      growth of particles for values of RH up to 0.97 using
!      Fitzgerald's scheme. For RH > 0.97, coagulation is neglected.
!---------------------------------------------------------------------
!
! First calculate relative humidity to determine hygroscopic
! growth factor ALPHA for particles; RH is obtained from q and
! qsat using RH=q(1-qsat)/(qsat(1-q))
!
! Use REL_HUM array to store returned QSAT values in call to QSAT_WAT
! to save space; then calculate fractional REL_HUM
!
      Q_ARRAY(:,:,:)=Q(1:row_length,1:rows,1:wet_model_levels)

      Do k=1,wet_model_levels
!
! DEPENDS ON: qsat_wat
        CALL QSAT_WAT(REL_HUM(1,1,k)                                    &
     &,               T(1,1,k), P_ARRAY(1,1,k), theta_field_size)
!cdir collapse
        Do j=1,rows
          Do i=1,row_length
!
            REL_HUM(i,j,k) = Q_ARRAY(i,j,k)*(1.0-REL_HUM(i,j,k)) /      &
     &                     ( REL_HUM(i,j,k)*(1.0-Q_ARRAY(i,j,k)) )
!
            If (REL_HUM(i,j,k)  <   0.3) Then
!
! Assume no hygroscopic growth occurs, so ALPHA=1
! Calculate factors multiplying SUM_IZ, SUM_IY (simplified for ALPHA=1)
!
              ZFAC = 2.0 * BOLTZMANN * T(i,j,k) *                       &
     &               SO4_AIT_ARRAY(i,j,k)*SO4_ACC_ARRAY(i,j,k)
              ZFAC = ZFAC * (RHO_AIR(i,j,k)/RHO_SO4)
              ZFAC = ZFAC/( 3.0*VISCOS_AIR(i,j,k) )
!
              YFAC = (1.333333*PI)**0.333333
              YFAC = YFAC*ACUNN*MEAN_FREEP(i,j,k)
!
              COAGRATE = ZFAC*(SUM_IZ + YFAC*SUM_IY)
!
            Else If (REL_HUM(i,j,k)  <=  0.97) Then
!
! Calculate hygroscopic growth factor ALPHA
!
! DEPENDS ON: hygro_fact
              CALL HYGRO_FACT(REL_HUM(i,j,k), ALPHA)
!
! Calculate factors multiplying SUM_IZ, SUM_IY (for ALPHA /= 1)
!
              ZFAC = 2.0 * BOLTZMANN * T(i,j,k) *                       &
     &               SO4_AIT_ARRAY(i,j,k)*SO4_ACC_ARRAY(i,j,k)
              ZFAC = ZFAC * (RHO_AIR(i,j,k)/RHO_SO4)
              ZFAC = ZFAC/( 3.0*VISCOS_AIR(i,j,k) )
              ZFAC = ZFAC/( ALPHA**6 )
!
              YFAC = (1.333333*PI)**0.333333
              YFAC = YFAC*ACUNN*MEAN_FREEP(i,j,k)
              YFAC = YFAC/ALPHA
!
              COAGRATE = ZFAC*(SUM_IZ + YFAC*SUM_IY)
!
            Else
!
! Neglect coagulation for RH > 0.97
!
              COAGRATE = 0.0

            End If         !End RH condition
!
! Calculate amount of SO4 tranferred from Ait to acc mode so4:
!
            DELTAS_COAG(i,j,k) = COAGRATE*TSTEP
!
          End Do           !End i loop
        End Do             !End j loop
!
      End Do               !End k loop
!
!-------------------------------------------------------------------
! 7. UPDATE SO2, SO4 modes and DMS, assuming:
!    Dry oxidation SO2 produces SO4 Aitken mode aerosol
!    Wet oxidation SO2 produces dissolved SO4
!    Dry oxidation DMS produces SO2 and MSA,
!                       and SO4 if L_SULPC_NEWDMS=T
!---------------------------------------------------------------------
!
! CHECK THAT AMOUNTS OF SO2 & DMS OXIDISED NOT GREATER THAN SO2 & DMS
! AVAILABLE AND SCALE INCREMENTS IF NECESSARY
!
!cdir collapse
      Do k=1,model_levels
        Do j=1,rows
          Do i=1,row_length
!
          DELTAS_TOT(i,j,k)=DELTAS_DRY(i,j,k)+DELTAS_WET(i,j,k)         &
     &                              + DELTAS_WET_O3(i,j,k)
!
          IF (DELTAS_TOT(i,j,k)  >   SO2_ARRAY(i,j,k))  THEN
               SCALE_FACTOR=SO2_ARRAY(i,j,k)/DELTAS_TOT(i,j,k)
              DELTAS_DRY(i,j,k)=DELTAS_DRY(i,j,k)*SCALE_FACTOR
              DELTAS_WET(i,j,k)=DELTAS_WET(i,j,k)*SCALE_FACTOR
!
            IF (L_SULPC_OZONE) THEN
              DELTAS_WET_O3(i,j,k)=DELTAS_WET_O3(i,j,k)*SCALE_FACTOR
           END IF
!
            DELTAS_TOT(i,j,k) = SO2_ARRAY(i,j,k)
          END IF
!
          DELTAS_TOT_AIT(i,j,k)=DELTAS_DIFFUSE(i,j,k)+DELTAS_COAG(i,j,k)
!
          IF ( DELTAS_TOT_AIT(i,j,k)  >   SO4_AIT_ARRAY(i,j,k) ) THEN
            SCALE_FACT_AIT=SO4_AIT_ARRAY(i,j,k)/DELTAS_TOT_AIT(i,j,k)
            DELTAS_DIFFUSE(i,j,k)=DELTAS_DIFFUSE(i,j,k)*SCALE_FACT_AIT
            DELTAS_COAG(i,j,k)=DELTAS_COAG(i,j,k)*SCALE_FACT_AIT
            DELTAS_TOT_AIT(i,j,k)=SO4_AIT_ARRAY(i,j,k)
          END IF
!
          If ( L_SULPC_DMS) Then
            IF (DELTAS_DMS(i,j,k)  >   DMS_ARRAY(i,j,k)) THEN
              DELTAS_DMS(i,j,k) = DMS_ARRAY(i,j,k)
            END IF
          End If
!
!  UPDATE SO2, SO4 MODES AND DMS
!
      If (L_SULPC_DMS) Then
!
        IF (L_SULPC_NEWDMS .AND. L_SULPC_OZONE) THEN
!
          SO2_ARRAY(i,j,k)=                                             &
     &       SO2_ARRAY(i,j,k) +DELTAS_DMS(i,j,k)*F_DMS_TO_SO2(i,j,k)    &
     &                       - DELTAS_TOT(i,j,k)
!
          SO4_AIT_ARRAY(i,j,k)=                                         &
     &                 SO4_AIT_ARRAY(i,j,k) - DELTAS_TOT_AIT(i,j,k)     &
     &                                  + PSI(i,j,k)*                   &
     &      (DELTAS_DRY(i,j,k)+DELTAS_DMS(i,j,k)*F_DMS_TO_SO4(i,j,k) )
!
          SO4_ACC_ARRAY(i,j,k)=SO4_ACC_ARRAY(i,j,k)                     &
     &                                  + DELTAS_EVAP(i,j,k)            &
     &                                  - DELTAS_NUCL(i,j,k)            &
     &                                  + (1.0-PSI(i,j,k))*             &
     &      (DELTAS_DRY(i,j,k)+DELTAS_DMS(i,j,k)*F_DMS_TO_SO4(i,j,k) )  &
     &                                  + DELTAS_COAG(i,j,k)
!
          SO4_DIS_ARRAY(i,j,k)=SO4_DIS_ARRAY(i,j,k)                     &
     &                                  + DELTAS_WET(i,j,k)             &
     &                                  - DELTAS_EVAP(i,j,k)            &
     &                                  + DELTAS_NUCL(i,j,k)            &
     &                                  + DELTAS_DIFFUSE(i,j,k)         &
     &                                  + DELTAS_WET_O3(i,j,k)
!
          DMS_ARRAY(i,j,k)=DMS_ARRAY(i,j,k) - DELTAS_DMS(i,j,k)
!
          MSA(i,j,k)=DELTAS_DMS(i,j,k)*F_DMS_TO_MSA(i,j,k)
!
!
        ELSE                 ! update for simple DMS scheme
!
          SO2_ARRAY(i,j,k)=SO2_ARRAY(i,j,k)                             &
     &                          + DELTAS_DMS(i,j,k)*BRAT_SO2            &
     &                          - DELTAS_TOT(i,j,k)
!
!
          SO4_AIT_ARRAY(i,j,k)=SO4_AIT_ARRAY(i,j,k)                     &
     &                                  + PSI(i,j,k)*DELTAS_DRY(i,j,k)  &
     &                                  - DELTAS_TOT_AIT(i,j,k)
!
          SO4_ACC_ARRAY(i,j,k)=SO4_ACC_ARRAY(i,j,k)                     &
     &                                  + DELTAS_EVAP(i,j,k)            &
     &                                  - DELTAS_NUCL(i,j,k)            &
     &                 + (1.0-PSI(i,j,k))*DELTAS_DRY(i,j,k)             &
     &                                  + DELTAS_COAG(i,j,k)
!
          SO4_DIS_ARRAY(i,j,k)=SO4_DIS_ARRAY(i,j,k)                     &
     &                                  + DELTAS_WET(i,j,k)             &
     &                                  - DELTAS_EVAP(i,j,k)            &
     &                                  + DELTAS_NUCL(i,j,k)            &
     &                                  + DELTAS_DIFFUSE(i,j,k)         &
     &                                  + DELTAS_WET_O3(i,j,k)
!
            DMS_ARRAY(i,j,k) = DMS_ARRAY(i,j,k) - DELTAS_DMS(i,j,k)
!
            MSA(i,j,k) = DELTAS_DMS(i,j,k)*BRAT_MSA
!
!
        END IF            ! END L_SULPC_NEWDMS TEST
!
      Else                   ! no DMS
!
          SO2_ARRAY(i,j,k)=SO2_ARRAY(i,j,k) - DELTAS_TOT(i,j,k)
!
!
          SO4_AIT_ARRAY(i,j,k)=SO4_AIT_ARRAY(i,j,k)                     &
     &                                  + PSI(i,j,k)*DELTAS_DRY(i,j,k)  &
     &                                  - DELTAS_TOT_AIT(i,j,k)
!
          SO4_ACC_ARRAY(i,j,k)=SO4_ACC_ARRAY(i,j,k)                     &
     &                                  + DELTAS_EVAP(i,j,k)            &
     &                                  - DELTAS_NUCL(i,j,k)            &
     &                 + (1.0-PSI(i,j,k))*DELTAS_DRY(i,j,k)             &
     &                                  + DELTAS_COAG(i,j,k)
!
          SO4_DIS_ARRAY(i,j,k)=SO4_DIS_ARRAY(i,j,k)                     &
     &                                  + DELTAS_WET(i,j,k)             &
     &                                  - DELTAS_EVAP(i,j,k)            &
     &                                  + DELTAS_NUCL(i,j,k)            &
     &                                  + DELTAS_DIFFUSE(i,j,k)         &
     &                                  + DELTAS_WET_O3(i,j,k)
!
!
      End If                 ! End L_SULPC_DMS condition
!
         End Do
        End Do
      End Do

      If (L_SULPC_DMS) Then
        DMS(1:row_length,1:rows,1:model_levels)=DMS_ARRAY(:,:,:)
      End If
      SO2(1:row_length,1:rows,1:model_levels)=SO2_ARRAY(:,:,:)
      SO4_ACC(1:row_length,1:rows,1:model_levels)=SO4_ACC_ARRAY(:,:,:)
      SO4_AIT(1:row_length,1:rows,1:model_levels)=SO4_AIT_ARRAY(:,:,:)
      SO4_DIS(1:row_length,1:rows,1:model_levels)=SO4_DIS_ARRAY(:,:,:)
!
!--------------------------------------------------------------------
! 8.  DEPLETE  NH3 and H2O2, REPLENISH H2O2 from HO2
!--------------------------------------------------------------------
!
!cdir collapse
      Do k=1,wet_model_levels
        Do j=1,rows
          Do i=1,row_length
!
!  Deplete H2O2
!
        IF (DELTAS_WET(i,j,k)  >   0.0)  THEN
          H2O2_MXR_ARRAY(i,j,k)=                                        &
     &          H2O2_MXR_ARRAY(i,j,k)-DELTAS_WET(i,j,k)/RELM_S_H2O2
        END IF
!
!   Replenish H2O2_MXR field from HO2 field in dry part of grid box
!   (whether or not oxidn has occurred) in daylight only.
!
        IF ( (H2O2_MXR_ARRAY(i,j,k)  <   H2O2_LMT(i,j,k)) .AND.         &
     &          (COSZA2D(i,j)  >   1.0E-06) )  THEN
!
!   Calculate replenishment in concn, then convert to mix ratio
!   Reaction rate K_HO2_HO2 is a fn of T, P, AND Q
!
          W_CONC=Q_ARRAY(i,j,k)*AIR_CONC(i,j,k)*RMM_AIR/RMM_W !MCL/CC
!
          H2O2_CON_TO_MXR=RMM_H2O2/(RMM_AIR*AIR_CONC(i,j,k))   !CC/MCL
!
          K_HO2_HO2=( K2*EXP(T2/T(i,j,k)) +                             &
     &                          AIR_CONC(i,j,k)*K3*EXP(T3/T(i,j,k)) )   &
     &                * (1.0 + W_CONC*K4*EXP(T4/T(i,j,k)) )
!
          H2O2_REP = HO2_CONC(i,j,k)*HO2_CONC(i,j,k)*K_HO2_HO2          &
     &                                        * TSTEP*CLEARF(i,j,k)
!
          H2O2_REP = H2O2_CON_TO_MXR * H2O2_REP
!
!  Increment H2O2_MXR up to H2O2_LMT value
!
          H2O2_MXR_ARRAY(i,j,k) = H2O2_MXR_ARRAY(i,j,k) + H2O2_REP
!
          IF ( H2O2_MXR_ARRAY(i,j,k)  >   H2O2_LMT(i,j,k) ) THEN
            H2O2_MXR_ARRAY(i,j,k) = H2O2_LMT(i,j,k)
          END IF
!
      END IF                     !End replenishment condition
!
! Deplete NH3 if ozone oxidation included and wet oxidn has occurred
!
      If (L_SULPC_OZONE) Then
!
      IF ( (DELTAS_WET(i,j,k)  >   0.0) .OR.                            &
     &                        (DELTAS_WET_O3(i,j,k)  >   0.0) )  THEN
!
        NH3_DEP(i,j,k) = (DELTAS_WET(i,j,k)+DELTAS_WET_O3(i,j,k))       &
     &                          / RELM_S_2N
!
          IF ( NH3_DEP(i,j,k)  >   NH3_ARRAY(i,j,k) ) THEN
            NH3_DEP(i,j,k) = NH3_ARRAY(i,j,k)
          END IF
!
        NH3_ARRAY(i,j,k) = NH3_ARRAY(i,j,k) - NH3_DEP(i,j,k)
!
      END IF       ! end wet oxidn test
!
      End If       ! End L_SULPC_OZONE Test
!
         End Do
        End Do
      End Do

      H2O2_MXR(1:row_length,1:rows,1:model_levels)=H2O2_MXR_ARRAY(:,:,:)
      NH3(1:row_length,1:rows,1:model_levels)=NH3_ARRAY(:,:,:)
!
!-----------------------------------------------------------------------
!
      RETURN
      END SUBROUTINE SULPHR
!
#endif
