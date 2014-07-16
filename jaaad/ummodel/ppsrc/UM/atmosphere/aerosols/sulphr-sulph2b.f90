
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
     &,            L_SULPC_OZONE, L_SULPC_SO2_O3_NONBUFFERED            &
     &,            L_SULPC_NH3                                          &
! Arguments IN/OUT
     &,            SO2, DMS, SO4_AIT, SO4_ACC, SO4_DIS                  &
     &,            NH3, H2O2_MXR                                        &
! Arguments OUT
     &,            MSA, NH3_DEP                                         &
     &,            F_DMS_TO_SO2, F_DMS_TO_SO4, F_DMS_TO_MSA             &
     &,            DELTAS_DRY, DELTAS_WET, DELTAS_WET_O3                &
     &,            DELTAS_TOT, DELTAS_DMS                               &
     &,            DELTAS_EVAP, DELTAS_NUCL, DELTAS_DIFFUSE             &
     &,            DELTAS_MERGE, DELTAS_COAG, PSI                       &
     &            )
!
!---------------------------------------------------------------------
! Purpose: To perform oxidation chemistry of sulphur dioxide to 3 modes
!          of Sulphate aerosol (Aitken, accumulation and dissolved),
!          and dimethyl sulphide to sulphur dioxide and methyl sulphonic
!          acid.
!          There is also exchange between the 3 modes of sulphate
!          aerosol due to nucleation, diffusion and evaporation
!          processes, and Aitken-accumulation mode merging caused by
!          particle growth due to the condensation of newly-formed
!          sulphuric acid onto Aitken mode sulphate particles.
!          Called by Aero_Ctl
!
! Current owners of code:    C.E. Johnson and J.G.L. Rae
!
! History:
! Version     Date     Comment
! -------     ----     -------
!
!  6.2  20/10/05  New deck.  Copy of SULPH2A with modifications for
!                 changes to condensation of sulphuric acid onto
!                 sulphate aerosol, and for inclusion of
!                 Aitken-accumulaition mode merging.
!                 N. Bellouin, O. Boucher, J. Haywood, C. Johnson,
!                 A. Jones and J. Rae.
!
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
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------
! CONTAINS DENSITY OF WATER, RHO_WATER.
! C_DENSTY for subroutine SF_EXCH
      REAL,PARAMETER:: RHOSEA = 1026.0 ! density of sea water (kg/m3)
      REAL,PARAMETER:: RHO_WATER = 1000.0! density of pure water (kg/m3)
! C_DENSTY end
! CONTAINS OTHER REQUIRED PARAMETERS
!-------------------COMDECK C_SULCHM--------------------------------
! Parameters for Sulphur Cycle Chemistry
      REAL                                                              &
     &     EVAPTAU,                                                     &
                          ! timescale for dissolved SO4 to evaporate
     &     NUCTAU,                                                      &
                          ! timescale for accumulation mode particles
!                           to nucleate once they enter a cloud.
     &     DIFFUSE_AIT,                                                 &
                          ! diffusion coefficient of Aitken particles
     &     K_SO2OH_HI,                                                  &
                                  ! high pressure reaction rate limit
     &     K_DMS_OH,                                                    &
                                  ! reaction rate for DMS+OH  cc/mcl/s
     &      K4_CH3SO2_O3,                                               &
                             ! Rate coeff for CH3SO2+O3 -> CH3SO3+O2
     &      K5_CH3SO3_HO2,                                              &
                             ! Rate coeff for CH3SO3+HO2 -> MSA+O2
     &      RMM_O3,                                                     &
                             ! relative molecular mass O3
     &     BRAT_SO2,                                                    &
                                  ! branching ratio for SO2 in DMS oxidn
     &     BRAT_MSA,                                                    &
                                  ! branching ratio for MSA in DMS oxidn
     &     AVOGADRO,                                                    &
                                 ! no. of molecules in 1 mole
     &     RMM_H2O2,                                                    &
                                 ! relative molecular mass H2O2 kg/mole
     &     RMM_AIR,                                                     &
                                 ! relative molecular mass dry air
     &     RMM_W,                                                       &
                                 ! relative molecular mass water
     &     RELM_S_H2O2,                                                 &
                                 ! rel atomic mass sulphur/RMM_H2O2
     &     RELM_S_2N,                                                   &
                              ! rel atomic mass Sulphur/2*Nitrogen
     &     PARH,                                                        &
                                ! power of temp dependence of K_SO2OH_LO
     &     K1,                                                          &
                                ! parameters for calcn of K_SO2OH_LO
     &     T1,                                                          &
                                !
     &     FC,                                                          &
                                ! parameters for interpolation between
     &     FAC1,                                                        &
                                !   LO and HI reaction rate limits
     &     K2,K3,K4,                                                    &
                                ! parameters for calcn of K_HO2_HO2
     &     T2,T3,T4,                                                    &
                                !
     &     CLOUDTAU,                                                    &
                                  ! air parcel lifetime in cloud
     &     CHEMTAU,                                                     &
                                  ! chem lifetime in cloud before oxidn
     &     O3_MIN,                                                      &
                              ! min mmr of O3 required for oxidn
     &     THOLD                  ! threshold for cloud liquid water
!
!
      PARAMETER (                                                       &
     &           EVAPTAU = 300.0,                                       &
                                              ! secs  (=5 mins)
     &             NUCTAU = 30.0,                                       &
                                          ! secs
     &       DIFFUSE_AIT = 1.7134E-9,                                   &
                                             ! sq m/s
     &        K_SO2OH_HI = 2.0E-12,                                     &
                                       ! cc/mcl/s from STOCHEM model
     &           K_DMS_OH = 9.1E-12,                                    &
                                          ! cc/mcl/s
     &       K4_CH3SO2_O3 = 1.0E-14,                                    &
                                        ! cc/mcl/s
     &      K5_CH3SO3_HO2 = 4.0E-11,                                    &
     &             RMM_O3 = 4.8E-2,                                     &
                                        ! kg/mole
     &          BRAT_SO2 = 0.9,                                         &
     &           BRAT_MSA = 1.0-BRAT_SO2,                               &
     &           AVOGADRO = 6.022E23,                                   &
                                          ! per mole
     &           RMM_H2O2 = 3.40E-2,                                    &
                                          ! kg/mole
     &            RMM_AIR = 2.896E-2,                                   &
                                          ! kg/mole
     &              RMM_W = 1.8E-2,                                     &
                                          ! kg/mole
     &        RELM_S_H2O2 = 3.206/3.40,                                 &
     &           RELM_S_2N = 3.206/2.80,                                &
     &               PARH = 3.3,                                        &
     &                K1 = 4.0E-31,                                     &
                                       ! (cc/mcl)2/s from STOCHEM
     &                 T1 = 300.0,                                      &
                                          ! K
     &                FC = 0.45,                                        &
                                        ! from STOCHEM model
     &              FAC1 = 1.1904,                                      &
                                    ! 0.75-1.27*LOG10(FC) from STOCHEM
     &                 K2 = 2.2E-13,                                    &
                                          ! cc/mcl/s
     &                 K3 = 1.9E-33,                                    &
                                          ! (cc/mcl)2/s
     &                 K4 = 1.4E-21,                                    &
                                          ! cc/mcl
     &                 T2 = 600.0,                                      &
                                          ! K
     &                 T3 = 890.0,                                      &
                                          ! K
     &                 T4 = 2200.0,                                     &
                                          ! K
     &           CLOUDTAU = 1.08E4,                                     &
                                          ! secs (=3 hours)
     &            CHEMTAU = 9.0E2,                                      &
                                          ! secs (=15 mins)
     &              O3_MIN = 1.6E-8,                                    &
                                        !(kg/kg, equiv. 10ppbv)
     &              THOLD = 1.0E-8                                      &
                                          ! kg/kg
     &          )
!
      REAL RAD_AIT,                                                     &
                            ! median radius of Aitken mode particles
     &     DIAM_AIT,                                                    &
                            !   "    diameter    "
     &     RAD_ACC,                                                     &
                            ! median radius of acccumulation mode
     &     DIAM_ACC,                                                    &
                            !   "    diameter    "
     &     CHI,                                                         &
                            ! mole fraction of S in particle
     &     RHO_SO4,                                                     &
                            ! density of  SO4 particle
     &     SIGMA,                                                       &
                            ! standard devn of particle size distn
!                                 for accumulation mode
     &     E_PARM,                                                      &
                            ! param relating size distns of Ait & Acc
     &     NUM_STAR         ! threshold concn of accu mode particles
                            !  below which PSI=1
!
      PARAMETER (                                                       &
     &           RAD_AIT = 6.5E-9,                                      &
                                             ! m
     &          DIAM_AIT = 2.0*RAD_AIT,                                 &
     &           RAD_ACC = 95.0E-9,                                     &
                                             ! m
     &          DIAM_ACC = 2.0*RAD_ACC,                                 &
     &               CHI = 32.0/132.0,                                  &
     &           RHO_SO4 = 1769.0,                                      &
                                              ! kg/m3
     &             SIGMA = 1.4,                                         &
     &            E_PARM = 0.9398,                                      &
     &          NUM_STAR = 1.0E6                                        &
                                             ! m-3
     &          )
!
      REAL BOLTZMANN       !Boltzmanns constant.
      REAL MFP_REF         !Reference value of mean free path.
      REAL TREF_MFP        !Reference temperature for mean free path.
      REAL PREF_MFP        !Reference pressure for mean free path.
      REAL SIGMA_AIT       !Geometric standard deviation of the Aitken
!                             mode distribution.
!
      PARAMETER (BOLTZMANN = 1.3804E-23)  ! J K-1
      PARAMETER (MFP_REF = 6.6E-8                                       &
                                          ! m
     &        ,  TREF_MFP = 293.15                                      &
                                          ! K
     &        ,  PREF_MFP = 1.01325E5)    ! Pa
      PARAMETER (SIGMA_AIT = 1.30)
!
!*---------------------------------------------------------------------
!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

!*----------------------------------------------------------------------
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
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
     &,L_SULPC_SO2_O3_NONBUFFERED                                       &
                                         !T if SO2+O3 reaction is NOT
                                         !  to be buffered by NH3.
     &,L_SULPC_NH3                       !T if NH3 buffering used
                                         !(always T if L_SULPC_OZONE
                                         ! is T)
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
!----------------COMDECK C_SULCOND---------------------------------
! Arrays and variables to be used in interpolation to determine
! dependence of partitioning of dry oxidised SO2 to Aitken and
! accumulation modes:
      REAL                                                              &
     &       REL_HUM_TABLE(101),                                        &
                                 ! Relative humidities between 0 and 1
     &       A_TABLE(101),                                              &
                              ! CdotN_Ait/Mtot_Ait at corresponding RH
     &       B_TABLE(101)     ! CdotN_Acc/Mtot_Acc at corresponding RH
      PARAMETER(REL_HUM_TABLE = (/                                      &
     &        0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,           &
     &        0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15,           &
     &        0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23,           &
     &        0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31,           &
     &        0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39,           &
     &        0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47,           &
     &        0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55,           &
     &        0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63,           &
     &        0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71,           &
     &        0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79,           &
     &        0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87,           &
     &        0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95,           &
     &        0.96, 0.97, 0.98, 0.99, 1.00 /) )
      PARAMETER(a_table = (/                                            &
     &   3.73408E+07, 3.73408E+07, 3.73408E+07, 3.73408E+07,            &
     &   3.73408E+07, 3.73408E+07, 3.73408E+07, 3.73408E+07,            &
     &   3.73408E+07, 3.73408E+07, 3.73408E+07, 3.73408E+07,            &
     &   3.73408E+07, 3.73408E+07, 3.73408E+07, 3.73408E+07,            &
     &   3.73408E+07, 3.73408E+07, 3.73408E+07, 3.73408E+07,            &
     &   3.73408E+07, 3.73408E+07, 3.73408E+07, 3.73408E+07,            &
     &   3.73408E+07, 3.73408E+07, 3.73408E+07, 3.73408E+07,            &
     &   3.73408E+07, 3.73408E+07, 3.73408E+07, 3.80313E+07,            &
     &   3.87275E+07, 3.94294E+07, 4.01369E+07, 4.08502E+07,            &
     &   4.15693E+07, 4.22940E+07, 4.30243E+07, 4.37602E+07,            &
     &   4.45019E+07, 4.52490E+07, 4.60017E+07, 4.67600E+07,            &
     &   4.75240E+07, 4.82933E+07, 4.90682E+07, 4.98487E+07,            &
     &   5.06346E+07, 5.14260E+07, 5.22229E+07, 5.30252E+07,            &
     &   5.38329E+07, 5.46460E+07, 5.54646E+07, 5.62885E+07,            &
     &   5.71177E+07, 5.79524E+07, 5.87925E+07, 5.96377E+07,            &
     &   6.04882E+07, 6.13442E+07, 6.22053E+07, 6.30717E+07,            &
     &   6.39433E+07, 6.48201E+07, 6.57019E+07, 6.65891E+07,            &
     &   6.74816E+07, 6.83793E+07, 6.92821E+07, 7.01899E+07,            &
     &   7.11028E+07, 7.20209E+07, 7.29444E+07, 7.38725E+07,            &
     &   7.48059E+07, 7.57442E+07, 7.66875E+07, 7.76360E+07,            &
     &   7.85893E+07, 7.95478E+07, 8.13327E+07, 8.33181E+07,            &
     &   8.55388E+07, 8.80381E+07, 9.08706E+07, 9.41064E+07,            &
     &   9.78351E+07, 1.02175E+08, 1.07284E+08, 1.13379E+08,            &
     &   1.20761E+08, 1.29866E+08, 1.41339E+08, 1.56177E+08,            &
     &   1.75991E+08, 2.03530E+08, 2.43830E+08, 3.06957E+08,            &
     &   4.15500E+08  /) )
      PARAMETER(b_table = (/                                            &
     &   8.49410e+05, 8.49410e+05, 8.49410e+05, 8.49410e+05,            &
     &   8.49410e+05, 8.49410e+05, 8.49410e+05, 8.49410e+05,            &
     &   8.49410e+05, 8.49410e+05, 8.49410e+05, 8.49410e+05,            &
     &   8.49410e+05, 8.49410e+05, 8.49410e+05, 8.49410e+05,            &
     &   8.49410e+05, 8.49410e+05, 8.49410e+05, 8.49410e+05,            &
     &   8.49410e+05, 8.49410e+05, 8.49410e+05, 8.49410e+05,            &
     &   8.49410e+05, 8.49410e+05, 8.49410e+05, 8.49410e+05,            &
     &   8.49410e+05, 8.49410e+05, 8.49410e+05, 8.59822e+05,            &
     &   8.70240e+05, 8.80669e+05, 8.91104e+05, 9.01551e+05,            &
     &   9.12005e+05, 9.22465e+05, 9.32935e+05, 9.43409e+05,            &
     &   9.53896e+05, 9.64389e+05, 9.74885e+05, 9.85393e+05,            &
     &   9.95906e+05, 1.00643e+06, 1.01695e+06, 1.02748e+06,            &
     &   1.03802e+06, 1.04857e+06, 1.05912e+06, 1.06967e+06,            &
     &   1.08023e+06, 1.09080e+06, 1.10137e+06, 1.11195e+06,            &
     &   1.12253e+06, 1.13312e+06, 1.14371e+06, 1.15431e+06,            &
     &   1.16491e+06, 1.17552e+06, 1.18613e+06, 1.19675e+06,            &
     &   1.20737e+06, 1.21799e+06, 1.22862e+06, 1.23925e+06,            &
     &   1.24989e+06, 1.26053e+06, 1.27117e+06, 1.28182e+06,            &
     &   1.29247e+06, 1.30312e+06, 1.31378e+06, 1.32445e+06,            &
     &   1.33511e+06, 1.34578e+06, 1.35645e+06, 1.36712e+06,            &
     &   1.37780e+06, 1.38848e+06, 1.40824e+06, 1.43001e+06,            &
     &   1.45410e+06, 1.48093e+06, 1.51096e+06, 1.54481e+06,            &
     &   1.58324e+06, 1.62724e+06, 1.67810e+06, 1.73752e+06,            &
     &   1.80784e+06, 1.89228e+06, 1.99547e+06, 2.12424e+06,            &
     &   2.28913e+06, 2.50720e+06, 2.80779e+06, 3.24555e+06,            &
     &   3.93362e+06 /) )
!
      INTEGER                                                           &
     &        nearest_index(row_length, rows, model_levels)
!  Index of element of REL_HUM_TABLE
!  that is closest to actual RH(i,j,k).
!
      REAL                                                              &
     &       dRH,                                                       &
                   ! Interval between values in REL_HUM_TABLE
     &       A_ARRAY(row_length, rows, model_levels),                   &
! A_TABLE(nearest_index) * MMR of Aitken sulphate
     &       B_ARRAY(row_length, rows, model_levels)
! B_TABLE(nearest_index) * MMR of accumulation sulphate
!
! c_sulcond.h contains arrays and variables used in calculation of PSI.
!
! Arrays and variables required for mode-merging calculation:
      REAL AIT_PRODUCTION(row_length, rows, model_levels)
      REAL FRAC_TRANS(row_length, rows, model_levels)
      REAL K_MERGE            ! FRAC_TRANS = K_MERGE * AIT_PRODUCTION
      PARAMETER(K_MERGE = 3.068E+09)
      REAL DELTAS_MERGE(row_length, rows, model_levels)
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
         IF (L_SULPC_OZONE) THEN
           IF (L_SULPC_SO2_O3_NONBUFFERED) Then
!
! SO2+O3 NOT buffered by NH3.
! When L_SULPC_SO2_O3_NONBUFFERED = T, O3 oxidation proceeds only if 
! O3 MMR is above a threshold (O3_MIN). The condition on NH3, which 
! applies when L_SULPC_SO2_O3_NONBUFFERED = F, does not apply when 
! L_SULPC_SO2_O3_NONBUFFERED = T.
!
             IF (O3(i,j,k) >= O3_MIN) THEN
!
! When L_SULPC_SO2_O3_NONBUFFERED = T, O3 oxidation is controlled by SO2 
! field only.  Note that all SO2 is used in DELTAS_WET_O3 calcn; 
! adjustment by SCALE_FACTOR at end of routine prevents removing too much.
!
                  DELTAS_WET_O3(i,j,k)=                                 &
     &                    (1.0-EXP(-WETRATE*TSTEP))*SO2_ARRAY(i,j,k)
             END IF               ! End ozone threshold condtion.
           ELSE
!
! SO2+O3 buffered by NH3.
!
! When L_SULPC_SO2_O3_NONBUFFERED = F, O3 oxidation proceeds only if O3 
! MMR is above a threshold (O3_MIN) and NH3 MMR is large enough that there 
! is still NH3 left in the gridbox after production of ammonium sulphate 
! in the aqueous reaction of SO2 with H2O2 (see above).
!
             IF ((O3(i,j,k) >= O3_MIN) .AND.                            &
     &           (NH3_ARRAY(i,j,k) > DELTAS_WET(i,j,k)/RELM_S_2N)) THEN
!
! When L_SULPC_SO2_O3_NONBUFFERED = F, O3 oxidation is controlled by s
! smaller of SO2 and NH3 fields. 2 atoms of N required for each S atom 
! in ammonium sulphate.  Sufficient NH3 must be available to neutralise 
! all sulphate produced in grid box for O3 reaction to continue (otherwise 
! PH too low) Note that all SO2 or NH3 is used in DELTAS_WET_O3 calcn; 
! adjustment by SCALE_FACTOR at end of routine prevents removing too much.
!
               IF(SO2_ARRAY(i,j,k)<=(NH3_ARRAY(i,j,k)*RELM_S_2N)) THEN
! SO2 limiting
                  DELTAS_WET_O3(i,j,k)=                                 &
     &                    (1.0-EXP(-WETRATE*TSTEP))*SO2_ARRAY(i,j,k)
!
               ELSE                                   
! NH3 limiting
                  DELTAS_WET_O3(i,j,k)=(1.0-EXP(-WETRATE*TSTEP))*       &
     &                                 NH3_ARRAY(i,j,k)*RELM_S_2N
!
               END IF
!
             END IF    ! End ozone and NH3 thresholds condition.




           END IF   ! End of IF(L_SULPC_SO2_O3_NONBUFFERED)-THEN-ELSE.
!
         END IF     ! End L_SULPC_OZONE test
!



       END IF           ! End peroxide oxidation block
!
       END IF           ! End cloud present condition
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
!---------------------------------------------------------------------
! 6a.          Calculation of PSI, fraction of dry-oxidied SO2 going
!              into Aitken mode (as opposed to accumulation mode).  The
!              method of calculating PSI has been changed, and the
!              calculation has been moved here from Section 1, as it
!              now depends on relative humidity, which is calculated in
!              Section 6.
!----------------------------------------------------------------------
      dRH = REL_HUM_TABLE(2) - REL_HUM_TABLE(1)
      DO k=1,model_levels
         DO j=1,rows
            DO i=1,row_length
               nearest_index(i,j,k) = (NINT((rel_hum(i,j,k)             &
     &                 - rel_hum_table(1)) / dRH)) + 1
               IF(nearest_index(i,j,k)  >   101) THEN
                  nearest_index(i,j,k) = 101
               ENDIF
               IF(nearest_index(i,j,k)  <   1) THEN
                  nearest_index(i,j,k) = 1
               ENDIF
               A_ARRAY(i,j,k) = A_TABLE(nearest_index(i,j,k))           &
     &                               * SO4_AIT_ARRAY(i,j,k)
               B_ARRAY(i,j,k) = B_TABLE(nearest_index(i,j,k))           &
     &                               * SO4_ACC_ARRAY(i,j,k)
               IF (A_ARRAY(i,j,k)  ==  0.0) THEN
                  PSI(i,j,k) = 1.0
               ELSE
                  PSI(i,j,k) = A_ARRAY(i,j,k) /                         &
     &                          (A_ARRAY(i,j,k) + B_ARRAY(i,j,k))
               ENDIF
            END DO        !End i loop
         END DO           !End j loop
      END DO              !End k loop
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
          If ( L_SULPC_DMS) Then
            IF (DELTAS_DMS(i,j,k)  >   DMS_ARRAY(i,j,k)) THEN
              DELTAS_DMS(i,j,k) = DMS_ARRAY(i,j,k)
            END IF
          End If
! Calculation of AIT_PRODUCTION, the amount of Aitken mode sulphate
! produced in a timestep:
               IF(L_SULPC_DMS) THEN
                  IF(L_SULPC_NEWDMS .AND. L_SULPC_OZONE) THEN
                    AIT_PRODUCTION(i,j,k) = PSI(i,j,k)                  &
     &                    * (DELTAS_DRY(i,j,k)                          &
     &                  + DELTAS_DMS(i,j,k) * F_DMS_TO_SO4(i,j,k))
                  ELSE
                     AIT_PRODUCTION(i,j,k) = PSI(i,j,k)                 &
     &                                        * DELTAS_DRY(i,j,k)
                  END IF
               ELSE
                  AIT_PRODUCTION(i,j,k) = PSI(i,j,k)                    &
     &                                        * DELTAS_DRY(i,j,k)
               END IF
!
! Fraction of Aitken mode transferred to accumulation mode in a
! timestep by mode merging:
               FRAC_TRANS(i,j,k) = K_MERGE * AIT_PRODUCTION(i,j,k)
! Change in MMR of Aitken mode sulphate as a result of this transfer:
               DELTAS_MERGE(i,j,k) = FRAC_TRANS(i,j,k)                  &
     &                                   * SO4_AIT_ARRAY(i,j,k)
!
          DELTAS_TOT_AIT(i,j,k)=DELTAS_DIFFUSE(i,j,k)+DELTAS_COAG(i,j,k)&
     &     + DELTAS_MERGE(i,j,k)
!
          IF ( DELTAS_TOT_AIT(i,j,k)  >   SO4_AIT_ARRAY(i,j,k) ) THEN
            SCALE_FACT_AIT=SO4_AIT_ARRAY(i,j,k)/DELTAS_TOT_AIT(i,j,k)
            DELTAS_DIFFUSE(i,j,k)=DELTAS_DIFFUSE(i,j,k)*SCALE_FACT_AIT
            DELTAS_COAG(i,j,k)=DELTAS_COAG(i,j,k)*SCALE_FACT_AIT
            DELTAS_MERGE(i,j,k) = DELTAS_MERGE(i,j,k) * SCALE_FACT_AIT
            DELTAS_TOT_AIT(i,j,k)=SO4_AIT_ARRAY(i,j,k)
          END IF
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
     &                                  + DELTAS_COAG(i,j,k)            &
     &     + DELTAS_MERGE(i,j,k)
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
     &                                  + DELTAS_COAG(i,j,k)            &
     &     + DELTAS_MERGE(i,j,k)
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
     &                                  + DELTAS_COAG(i,j,k)            &
     &     + DELTAS_MERGE(i,j,k)
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
