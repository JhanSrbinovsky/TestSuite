
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
     SUBROUTINE PM10_PM2P5 (                                               &
!
! Arguments IN
!
! Parallel variables
        off_x, off_y,                                                      &
!   Model dimensions
        row_length, rows,                                                  &
        model_levels,                                                      &
        salt_dim1, salt_dim2, salt_dim3,                                   &
!   Logicals IN
        L_SULPC_SO2, L_SOOT, L_BIOMASS, L_OCFF, L_USE_BIOGENIC, L_DUST,    &   
        L_seasalt_CCN, L_USE_SEASALT_DIRECT,                               & 
        L_USE_SEASALT_INDIRECT, L_USE_SEASALT_AUTOCONV,                    &
!   Data fields IN        
        p_theta_levels, T,                                                 &
        SO4_AIT, SO4_ACC, SOOT_NEW, SOOT_AGD, BMASS_NEW, BMASS_AGD,        &
        OCFF_NEW, OCFF_AGD, BIOGENIC, sea_salt_film, sea_salt_jet,         &
        DUST_1, DUST_2, DUST_3, DUST_4, DUST_5,                            &
!
! Arguments OUT
        PM10, PM2p5,                                                       &
        PM10_SO4, PM2p5_SO4,                                               &
        PM10_BC, PM2p5_BC,                                                 &
        PM10_BB, PM2p5_BB,                                                 &
        PM10_OCFF, PM2p5_OCFF,                                             &
        PM10_SOA, PM2p5_SOA,                                               &
        PM10_SS, PM2p5_SS,                                                 &
        PM10_DUST, PM2p5_DUST                                              &
     )
   
!
!--------------------------------------------------------------------------
! Purpose: Calculation of the PM10 and PM2.5 mass concentrations (in 
!          microgram/m3) by using the cumulative volume (from 0 to 
!          cut-off diameter) of a log-normal, with median radius 
!          (or gemetric mean radius) r_bar, geometric standard deviation 
!          sigma and total number Ntot. 
!          The contributions by the different aerosol species and modes 
!          are added up to get the total PM10 & PM2.5 concs (microgram/m3).
!           
!          Formulation derived from eqs. (7.51) & (7.51a) of Seinfeld
!          and Pandis (1998) using method of (7.39). Also, checked with
!          "The Mechanics of Inhaled Pharmaceutical Aerosols: An
!          Introduction", by Warren A. Finlay, Academic Press, 2001
!
!          Called by Aero_Ctl
!
! Current owner of code:    C. Ordonez
!
! Language:  FORTRAN 90
!
! Documentation:  Will be available in UMDP 20
!
!--------------------------------------------------------------------------
!
    
      IMPLICIT NONE      
! CONTAINS GAS CONST FOR DRY AIR: R=287.05 J/KG/K
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
! CONTAINS VALUE OF PI
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
!
! Model dimensions
      INTEGER, INTENT(IN) ::                                               &
        off_x,                                                             &
                   ! Size of small halo in i
        off_y,                                                             &
                   ! Size of small halo in j.
        row_length,                                                        &
        rows,                                                              &
        model_levels,                                                      &
        salt_dim1,                                                         &
        salt_dim2,                                                         &
        salt_dim3                                                           
                                                 !dimens of seasalt array
! 
      LOGICAL, INTENT(IN) ::                                               &
        L_SULPC_SO2,                                                       &
                                         !T if Sulphur Cycle required
        L_SOOT,                                                            &
                                         !T if SOOT modelling required
        L_BIOMASS,                                                         &
                                         !T if biomass modelling reqd
        L_OCFF,                                                            &
                                         !T if OCFF modelling required
        L_USE_BIOGENIC,                                                    &
                                         !T if clim. of biog. aeros. used
        L_DUST,                                                            &
                                         !T if mineral dust used
        L_seasalt_CCN,                                                     &
                                         !T if sea-salt used for CCN 
        L_USE_SEASALT_DIRECT,                                              &
                                         !T if SS direct rad. effect      
        L_USE_SEASALT_INDIRECT,                                            &
                                         !T if SS 1st indir. effect  
        L_USE_SEASALT_AUTOCONV
                                         !T if SS 2nd indir. effect   

! Data fields IN
      REAL, INTENT(IN) ::                                                  &
        p_theta_levels(1-off_x:row_length+off_x,                           &
                       1-off_y:rows+off_y, model_levels),                  &
                                           ! pressure on theta levels
      
        T(row_length,rows,model_levels),                                   &                                
                                            ! Temp (K) on theta levels
        SO4_AIT(1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
                                      model_levels),                       &
                                                         !mmr S in AIT
        SO4_ACC(1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
                                      model_levels),                       &
                                                         !mmr S in ACC
        SOOT_NEW(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
                                      model_levels),                       &
                                                      !mmr fresh soot
        SOOT_AGD(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
                                      model_levels),                       &
                                                      !mmr aged soot
        BMASS_NEW(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
                                      model_levels),                       &
                                                      !mmr fresh smoke
        BMASS_AGD(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
                                      model_levels),                       &
                                                      !mmr aged smoke
        OCFF_NEW(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
                                      model_levels),                       &
                                                      !mmr fresh OCFF
        OCFF_AGD(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
                                      model_levels),                       &
                                                      !mmr aged OCFF
        BIOGENIC(row_length, rows, model_levels),                          &
                                                      !mmr biogenics
        sea_salt_film(salt_dim1,salt_dim2,salt_dim3),                      &
        sea_salt_jet(salt_dim1,salt_dim2,salt_dim3),                       &
                                                      ! Sea-salt
                                                      ! aerosol nrs.
        DUST_1(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &         
                                      model_levels),                       &      
                                                      !mmr Dust div 1
        DUST_2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &         
                                      model_levels),                       &      
                                                      !mmr Dust div 2
        DUST_3(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &         
                                      model_levels),                       &      
                                                      !mmr Dust div 3
        DUST_4(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &        
                                      model_levels),                       &    
                                                      !mmr Dust div 4
        DUST_5(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &        
                                      model_levels)     
                                                      !mmr Dust div 5
     
! Output arguments      
      REAL, INTENT(OUT) ::                                                 &
        PM10 (row_length, rows, model_levels),                             &
                                                      !PM10  (ug m-3)
        PM2p5(row_length, rows, model_levels),                             &                            
                                                      !PM2.5 (ug m-3)
        PM10_SO4 (row_length, rows, model_levels),                         &
        PM2p5_SO4(row_length, rows, model_levels),                         &
                                      !PM concentrations due to sulphate
        PM10_BC (row_length, rows, model_levels),                          &
        PM2p5_BC(row_length, rows, model_levels),                          &
                                      !PM concs. due to black carbon
        PM10_BB (row_length, rows, model_levels),                          &
        PM2p5_BB(row_length, rows, model_levels),                          &
                             !PM concs. due to biomass burning aerosol 
        PM10_OCFF (row_length, rows, model_levels),                        &
        PM2p5_OCFF(row_length, rows, model_levels),                        &
                                      !PM concs. due to OCFF 
        PM10_SOA (row_length, rows, model_levels),                         &
        PM2p5_SOA(row_length, rows, model_levels),                         &
                                      !PM concs. due to SOA 
        PM10_SS (row_length, rows, model_levels),                          &
        PM2p5_SS(row_length, rows, model_levels),                          &
                                      !PM concs. due to sea-salt aerosol 
        PM10_DUST (row_length, rows, model_levels),                        &
        PM2p5_DUST(row_length, rows, model_levels)
                                      !PM concs. due to mineral dust 
 

!Local variables for PM10 and PM2.5 calculations        
      REAL ::                                                              &
          rho_air, N_tot, d_cut1, d_cut2,                                  &
                     ! air density, total nr. aerosol conc. & cut-off diams.
          rho_SU,                                                          &
          r_bar_SU_Ait, sigma_SU_Ait, denom_SU_Ait,                        &
          A_param_SU_Ait, B_param_SU_Ait,                                  &
          r_bar_SU_acc, sigma_SU_acc, denom_SU_acc,                        &
          A_param_SU_acc, B_param_SU_acc,                                  &
                              ! Variables for sulphate
          rho_BC,                                                          &
          r_bar_BC_fr, sigma_BC_fr, denom_BC_fr,                           &
          A_param_BC_fr, B_param_BC_fr,                                    &
          r_bar_BC_ag, sigma_BC_ag, denom_BC_ag,                           &
          A_param_BC_ag, B_param_BC_ag,                                    &
                              ! Variables for black carbon
          rho_BB,                                                          &
          r_bar_BB_fr, sigma_BB_fr, denom_BB_fr,                           &
          A_param_BB_fr, B_param_BB_fr,                                    &
          r_bar_BB_ag, sigma_BB_ag, denom_BB_ag,                           &
          A_param_BB_ag, B_param_BB_ag,                                    &
                              ! Variables for biomass burning aerosol
          rho_OCFF,                                                        &
          r_bar_OCFF_fr, sigma_OCFF_fr, denom_OCFF_fr,                     &
          A_param_OCFF_fr, B_param_OCFF_fr,                                &
          r_bar_OCFF_ag, sigma_OCFF_ag, denom_OCFF_ag,                     &
          A_param_OCFF_ag, B_param_OCFF_ag,                                &
                              ! Variables for OCFF
          rho_SOA,                                                         &
          r_bar_SOA, sigma_SOA, denom_SOA,                                 &
          A_param_SOA, B_param_SOA,                                        &
                              ! Variables for SOA
          rho_SS,                                                          &
          r_bar_SS_fi, sigma_SS_fi,                                        &
          A_param_SS_fi, B_param_SS_fi,                                    &
          r_bar_SS_je, sigma_SS_je,                                        &
          A_param_SS_je, B_param_SS_je
                              ! Variables for sea salt aerosol
      INTEGER :: i, j         ! loop counters (horizontal field indices)
      INTEGER :: k            ! loop counter  (vertical index)
      REAL    :: ERF          ! error function                     
      CHARACTER(LEN=80)  :: cmessage
!           
!     Method: 
!
!     (1) Cumulative_volume for particles with diameter between
!     0 and d_cut (assuming spherical particles): 
!            V(d_cut) = (A*Ntot/2.0) + (A*Ntot/2.0) *
!                   erf((alog(d_cut)-B)/(sqrt(2.0)*alog(sigma)))
!     with
!            A = (pi/6.0)*exp(3.0*alog(dbar)+4.5*(alog(sigma))**2)
!            B = alog (volume median diameter) =
!              = alog(d_bar)+3.0*(alog(sigma))**2
!     See also note (3) for the calculation of Ntot   
!
!     (2) Known V(d_cut) in m3, the mass concentration (ug / m3) of
!         particles with diameter lower than d_cut is:
!            V(d_cut) * rho_particle * 1e9
!            (need to multiply by 1e9 because rho is in kg/m3)
!
!     (3) Sea salt diagnostics are aerosol number (N) while the
!     diagnostics for the other aerosol species are m.m.r.
!     In those cases one can calculate Ntot with the relationship:
! 
!                            rho_air           1   
!      N (m-3) = m.m.r. * ------------  * -------------, with  
!                         rho_particle     V_particle
!
!      V_particle = (4*PI/3) * (r_bar**3) * exp (4.5*(ln (sigma))**2),
!      because 2*r_bar is the geometric diameter
!
!                                          3 * rho_air
!      We finally use:  N (m-3) = m.m.r. * ------------- 
!                                            denom
!                      (where denom = 3*rho_particle*V_particle)
!
!     (4) The size distribution of mineral dust is currently not lognormal. 
!         Six size bins are used instead, covering 
!         3.162e-8 to 1.000e-7 m (bin 1), 1.000e-8 to 3.162e-7 m (bin 2),
!         3.162e-7 to 1.000e-6 m (bin 3), 1.000e-6 to 3.162e-6 m (bin 4), 
!         3.162e-6 to 1.000e-5 (bin 5), and 1.000e-5 to 3.162e-5 m (bin 6). 
!         Consequently, PM2.5 concentrations are derived from m.m.r. of 
!         bins 1 to 3 and 19.3941% of bin 4, and PM10 concentrations are 
!         derived from m.m.r. of bins 1 to 4 and 39.8317% of bin 5: 
!            m.m.r. * rho_air (kg/m3) * 1e9 --> ug / m3
!
!     (5) All contributions are finally added to get the total
!         PM10 and PM2.5 mass concentrations
!   


!     Sulphate

      rho_SU         = 1769.0                !all densities in kg/m3
      r_bar_SU_Ait   = 0.0065e-6             !all radii in metres
      sigma_SU_Ait   = 1.3
      r_bar_SU_acc   = 0.095e-6
      sigma_SU_acc   = 1.4
      denom_SU_Ait   = 4.0*Pi*rho_SU*(r_bar_SU_Ait**3)                     &
                       *exp(4.5*(alog(sigma_SU_Ait))**2)
      A_param_SU_Ait = Pi*(1.0/6.0)                                        &
            *exp(3.0*alog(r_bar_SU_Ait*2.0)+4.5*(alog(sigma_SU_Ait))**2)
      B_param_SU_Ait = alog(r_bar_SU_Ait*2.0)+3.0*(alog(sigma_SU_Ait))**2
      denom_SU_acc   = 4.0*Pi*rho_SU*(r_bar_SU_acc**3)                     &
                       *exp(4.5*(alog(sigma_SU_acc))**2)
      A_param_SU_acc = Pi*(1.0/6.0)                                        &
            *exp(3.0*alog(r_bar_SU_acc*2.0)+4.5*(alog(sigma_SU_acc))**2)
      B_param_SU_acc = alog(r_bar_SU_acc*2.0)+3.0*(alog(sigma_SU_acc))**2


!     Black carbon

      rho_BC        = 1900.0
      r_bar_BC_fr   = 0.04e-6
      sigma_BC_fr   = 2.0
      r_bar_BC_ag   = 0.04e-6
      sigma_BC_ag   = 2.0
      denom_BC_fr   = 4.0*Pi*rho_BC*(r_bar_BC_fr**3)                       &
                      *exp(4.5*(alog(sigma_BC_fr))**2)
      A_param_BC_fr = Pi*(1.0/6.0)                                         &
            *exp(3.0*alog(r_bar_BC_fr*2.0)+4.5*(alog(sigma_BC_fr))**2)
      B_param_BC_fr = alog(r_bar_BC_fr*2.0)+3.0*(alog(sigma_BC_fr))**2
      denom_BC_ag   = 4.0*Pi*rho_BC*(r_bar_BC_ag**3)                       &
                      *exp(4.5*(alog(sigma_BC_ag))**2)
      A_param_BC_ag = Pi*(1.0/6.0)                                         &
            *exp(3.0*alog(r_bar_BC_ag*2.0)+4.5*(alog(sigma_BC_ag))**2)
      B_param_BC_ag = alog(r_bar_BC_ag*2.0)+3.0*(alog(sigma_BC_ag))**2
        

!     Biomass-burning aerosol

      rho_BB        = 1350.0
      r_bar_BB_fr   = 0.1e-6
      sigma_BB_fr   = 1.3
      r_bar_BB_ag   = 0.12e-6
      sigma_BB_ag   = 1.3
      denom_BB_fr   = 4.0*Pi*rho_BB*(r_bar_BB_fr**3)                       &
                      *exp(4.5*(alog(sigma_BB_fr))**2)
      A_param_BB_fr = Pi*(1.0/6.0)                                         &
            *exp(3.0*alog(r_bar_BB_fr*2.0)+4.5*(alog(sigma_BB_fr))**2)
      B_param_BB_fr = alog(r_bar_BB_fr*2.0)+3.0*(alog(sigma_BB_fr))**2
      denom_BB_ag   = 4.0*Pi*rho_BB*(r_bar_BB_ag**3)                       &
                      *exp(4.5*(alog(sigma_BB_ag))**2)
      A_param_BB_ag = Pi*(1.0/6.0)                                         &
            *exp(3.0*alog(r_bar_BB_ag*2.0)+4.5*(alog(sigma_BB_ag))**2)
      B_param_BB_ag = alog(r_bar_BB_ag*2.0)+3.0*(alog(sigma_BB_ag))**2
        

!     Organic Carbon (from fossil fuels) aerosol

      rho_OCFF        = 1350.0
      r_bar_OCFF_fr   = 0.1e-6
      sigma_OCFF_fr   = 1.3
      r_bar_OCFF_ag   = 0.12e-6
      sigma_OCFF_ag   = 1.3
      denom_OCFF_fr   = 4.0*Pi*rho_OCFF*(r_bar_OCFF_fr**3)                 &
                        *exp(4.5*(alog(sigma_OCFF_fr))**2)
      A_param_OCFF_fr = Pi*(1.0/6.0)                                       &
         *exp(3.0*alog(r_bar_OCFF_fr*2.0)+4.5*(alog(sigma_OCFF_fr))**2)
      B_param_OCFF_fr = alog(r_bar_OCFF_fr*2.0)                            &
                        +3.0*(alog(sigma_OCFF_fr))**2
      denom_OCFF_ag   = 4.0*Pi*rho_OCFF*(r_bar_OCFF_ag**3)                 &
                        *exp(4.5*(alog(sigma_OCFF_ag))**2)
      A_param_OCFF_ag = Pi*(1.0/6.0)                                       &
         *exp(3.0*alog(r_bar_OCFF_ag*2.0)+4.5*(alog(sigma_OCFF_ag))**2)
      B_param_OCFF_ag = alog(r_bar_OCFF_ag*2.0)                            &
                        +3.0*(alog(sigma_OCFF_ag))**2
        

!    Biogenic Secondary Organic Aerosol

      rho_SOA     = 1300.0
      r_bar_SOA   = 0.095e-6
      sigma_SOA   = 1.5
      denom_SOA   = 4.0*Pi*rho_SOA*(r_bar_SOA**3)                          &
                    *exp(4.5*(alog(sigma_SOA))**2)
      A_param_SOA = Pi*(1.0/6.0)                                           &
                    *exp(3.0*alog(r_bar_SOA*2.0)+4.5*(alog(sigma_SOA))**2)
      B_param_SOA = alog(r_bar_SOA*2.0)+3.0*(alog(sigma_SOA))**2


!     Sea-salt aerosol

      rho_SS        = 2165.0
      r_bar_SS_fi   = 0.1e-6
      sigma_SS_fi   = 1.9
      r_bar_SS_je   = 1.0e-6
      sigma_SS_je   = 2.0
      A_param_SS_fi = Pi*(1.0/6.0)                                         &
            *exp(3.0*alog(r_bar_SS_fi*2.0)+4.5*(alog(sigma_SS_fi))**2)
      B_param_SS_fi = alog(r_bar_SS_fi*2.0)+3.0*(alog(sigma_SS_fi))**2
      A_param_SS_je = Pi*(1.0/6.0)                                         &
            *exp(3.0*alog(r_bar_SS_je*2.0)+4.5*(alog(sigma_SS_je))**2)
      B_param_SS_je = alog(r_bar_SS_je*2.0)+3.0*(alog(sigma_SS_je))**2

      
      
!     Add contributions from the different aerosol species and modes 
!     to calculate pm10 & pm2.5 concentrations
      
      d_cut1 = 10.0e-06      !cut-off aerodyn. diameter for PM10  (metres)
      d_cut2 =  2.5e-06      !cut-off aerodyn. diameter for PM2.5 (metres)
      
      !initialise pm10 & pm2.5 concentrations
      PM10       (:,:,:) = 0.0
      PM2p5      (:,:,:) = 0.0
      PM10_SO4   (:,:,:) = 0.0
      PM2p5_SO4  (:,:,:) = 0.0
      PM10_BC    (:,:,:) = 0.0
      PM2p5_BC   (:,:,:) = 0.0
      PM10_BB    (:,:,:) = 0.0
      PM2p5_BB   (:,:,:) = 0.0
      PM10_OCFF  (:,:,:) = 0.0
      PM2p5_OCFF (:,:,:) = 0.0
      PM10_SOA   (:,:,:) = 0.0
      PM2p5_SOA  (:,:,:) = 0.0
      PM10_SS    (:,:,:) = 0.0
      PM2p5_SS   (:,:,:) = 0.0
      PM10_DUST  (:,:,:) = 0.0
      PM2p5_DUST (:,:,:) = 0.0
      

      DO k = 1, model_levels 
        DO j = 1, rows      
          DO i = 1, row_length

            rho_air = p_theta_levels(i, j, k)/(R * T(i, j, k))

            IF (L_SULPC_SO2) THEN
              !What is advected is sulphur. Therefore,      
              !necessary to convert m.m.r. of sulphur to ammonium
              !sulphate by multiplying by ratio of molecular weights: 
              !(NH4)2SO4 / S = 132 / 32 =  4.125
              N_tot = 3.0*(SO4_AIT(i, j, k)*4.125)*rho_air / denom_SU_Ait
              
              PM10_SO4(i, j, k) =                                          &
                                  ((A_param_SU_Ait * N_tot / 2.0)          &
                                  +(A_param_SU_Ait * N_tot / 2.0)          &
                                  * erf ((alog(d_cut1) - B_param_SU_Ait)   &
                                  / (sqrt(2.0) * alog(sigma_SU_Ait))))     &
                                  * rho_SU * 1.0e9
              
              PM2p5_SO4(i, j, k) =                                         &
                                  ((A_param_SU_Ait * N_tot/2.0)            &
                                  +(A_param_SU_Ait * N_tot/2.0)            &
                                  * erf((alog(d_cut2) - B_param_SU_Ait)    &
                                  / (sqrt(2.0) * alog(sigma_SU_Ait))))     &
                                  * rho_SU * 1.0e9

              N_tot = 3.0*(SO4_ACC(i, j, k)*4.125)*rho_air / denom_SU_acc
      
              PM10_SO4(i, j, k) = PM10_SO4(i, j, k) +                      &
                                  (((A_param_SU_acc * N_tot / 2.0)         &
                                  + (A_param_SU_acc * N_tot / 2.0)         &
                                  * erf((alog(d_cut1) - B_param_SU_acc)    &
                                  / (sqrt(2.0) * alog(sigma_SU_acc))))     &
                                  * rho_SU * 1.0e9)

              PM2p5_SO4(i, j, k) = PM2p5_SO4(i, j, k) +                    &
                                   (((A_param_SU_acc * N_tot/2.0)          &
                                   + (A_param_SU_acc * N_tot/2.0)          &
                                   * erf((alog(d_cut2) - B_param_SU_acc)   &
                                   / (sqrt(2.0) * alog(sigma_SU_acc))))    &
                                   * rho_SU * 1.0e9)
     
              PM10  (i, j, k) = PM10_SO4 (i, j, k) 
              PM2p5 (i, j, k) = PM2p5_SO4(i, j, k) 
            ELSE
              cmessage='Sulphur cycle must be selected to calculate PM diags.'
! DEPENDS ON: ereport
              CALL EREPORT('PM10_PM2P5',1,cmessage)
            END IF

                 
            IF (L_SOOT) THEN
              N_tot = 3.0 * SOOT_NEW(i, j, k) * rho_air / denom_BC_fr
      
              PM10_BC(i, j, k) =                                           &
                                  (((A_param_BC_fr * N_tot / 2.0)          &
                                  + (A_param_BC_fr * N_tot / 2.0)          &
                                  * erf((alog(d_cut1) - B_param_BC_fr)     &
                                  / (sqrt(2.0) * alog(sigma_BC_fr))))      &
                                  * rho_BC * 1.0e9)

              PM2p5_BC(i, j, k) =                                          &
                                  (((A_param_BC_fr * N_tot / 2.0)          &
                                  + (A_param_BC_fr * N_tot / 2.0)          &
                                  *erf((alog(d_cut2) - B_param_BC_fr)      &
                                  / (sqrt(2.0) * alog(sigma_BC_fr))))      &
                                  * rho_BC * 1.0e9)
     
              N_tot = 3.0 * SOOT_AGD (i, j, k) * rho_air / denom_BC_ag
      
              PM10_BC(i, j, k) =  PM10_BC(i, j, k) +                       &
                                  (((A_param_BC_ag * N_tot / 2.0)          &
                                  + (A_param_BC_ag * N_tot / 2.0)          &
                                  * erf((alog(d_cut1) - B_param_BC_ag)     &
                                  / (sqrt(2.0) * alog(sigma_BC_ag))))      &
                                  * rho_BC * 1.0e9)
     
              PM2p5_BC(i, j, k) = PM2p5_BC(i, j, k) +                      &
                                  (((A_param_BC_ag * N_tot / 2.0)          &
                                  + (A_param_BC_ag * N_tot / 2.0)          &
                                  * erf((alog(d_cut2) - B_param_BC_ag)     &
                                  / (sqrt(2.0) * alog(sigma_BC_ag))))      &
                                  * rho_BC * 1.0e9)
     
              PM10  (i, j, k) = PM10  (i, j, k) + PM10_BC (i, j, k)
              PM2p5 (i, j, k) = PM2p5 (i, j, k) + PM2p5_BC(i, j, k)
            END IF

         
            IF (L_BIOMASS) THEN
              N_tot = 3.0 * BMASS_NEW(i, j, k) * rho_air / denom_BB_fr
      
              PM10_BB(i, j, k) =                                           &
                                  (((A_param_BB_fr * N_tot / 2.0)          &
                                  + (A_param_BB_fr * N_tot / 2.0)          &
                                  * erf((alog(d_cut1) - B_param_BB_fr)     &
                                  / (sqrt(2.0) * alog(sigma_BB_fr))))      &
                                  * rho_BB * 1.0e9)

              PM2p5_BB(i, j, k) =                                          &
                                  (((A_param_BB_fr * N_tot / 2.0)          &
                                  + (A_param_BB_fr * N_tot / 2.0)          &
                                  * erf((alog(d_cut2) - B_param_BB_fr)     &
                                  / (sqrt(2.0) * alog(sigma_BB_fr))))      &
                                  * rho_BB * 1.0e9)

     
              N_tot = 3.0 * BMASS_AGD(i, j, k) * rho_air / denom_BB_ag
      
              PM10_BB(i, j, k) =  PM10_BB(i, j, k) +                       &
                                  (((A_param_BB_ag * N_tot / 2.0)          &
                                  + (A_param_BB_ag * N_tot / 2.0)          &
                                  * erf((alog(d_cut1) - B_param_BB_ag)     &
                                  / (sqrt(2.0) * alog(sigma_BB_ag))))      &
                                  * rho_BB * 1.0e9)
            
              PM2p5_BB(i, j, k) = PM2p5_BB(i, j, k) +                      &
                                  (((A_param_BB_ag * N_tot / 2.0)          &
                                  + (A_param_BB_ag * N_tot / 2.0)          &
                                  * erf((alog(d_cut2) - B_param_BB_ag)     &
                                  / (sqrt(2.0) * alog(sigma_BB_ag))))      &
                                  * rho_BB * 1.0e9)
     
              PM10  (i, j, k) = PM10  (i, j, k) + PM10_BB (i, j, k)
              PM2p5 (i, j, k) = PM2p5 (i, j, k) + PM2p5_BB(i, j, k)
            END IF

            
            IF (L_OCFF) THEN
              N_tot = 3.0 * OCFF_NEW(i, j, k) * rho_air / denom_OCFF_fr
      
              PM10_OCFF(i, j, k) =                                         &
                                   (((A_param_OCFF_fr * N_tot / 2.0)       &
                                    +(A_param_OCFF_fr * N_tot / 2.0)       &
                                    *erf((alog(d_cut1) - B_param_OCFF_fr)  &
                                    /(sqrt(2.0) * alog(sigma_OCFF_fr))))   &
                                    * rho_OCFF * 1.0e9)
      
              PM2p5_OCFF(i, j, k) =                                        &
                                    (((A_param_OCFF_fr * N_tot/2.0)        &
                                    + (A_param_OCFF_fr * N_tot/2.0)        &
                                    * erf((alog(d_cut2) - B_param_OCFF_fr) &
                                    / (sqrt(2.0) * alog(sigma_OCFF_fr))))  &
                                    * rho_OCFF * 1.0e9)

              N_tot = 3.0 * OCFF_AGD(i, j, k) * rho_air / denom_OCFF_ag
     
              PM10_OCFF(i, j, k) =  PM10_OCFF(i, j, k) +                   &
                                    (((A_param_OCFF_ag * N_tot / 2.0)      &
                                    + (A_param_OCFF_ag * N_tot / 2.0)      &
                                    * erf((alog(d_cut1) - B_param_OCFF_ag) &
                                    / (sqrt(2.0) * alog(sigma_OCFF_ag))))  &
                                    * rho_OCFF * 1.0e9)
      
              PM2p5_OCFF(i, j, k) = PM2p5_OCFF(i, j, k) +                  &
                                    (((A_param_OCFF_ag * N_tot/2.0)        &
                                    + (A_param_OCFF_ag * N_tot/2.0)        &
                                    * erf((alog(d_cut2) - B_param_OCFF_ag) &
                                    / (sqrt(2.0) * alog(sigma_OCFF_ag))))  &
                                    * rho_OCFF * 1.0e9)
     
              PM10  (i, j, k) = PM10  (i, j, k) + PM10_OCFF (i, j, k)
              PM2p5 (i, j, k) = PM2p5 (i, j, k) + PM2p5_OCFF(i, j, k)
            END IF
     
            IF (L_USE_BIOGENIC) THEN
              N_tot = 3.0 * BIOGENIC(i, j, k) * rho_air / denom_SOA
      
              PM10_SOA(i, j, k) =                                          &
                                   (((A_param_SOA * N_tot / 2.0)           &
                                   + (A_param_SOA * N_tot / 2.0)           &
                                   * erf((alog(d_cut1) - B_param_SOA)      &
                                   / (sqrt(2.0) * alog(sigma_SOA))))       &
                                   * rho_SOA * 1.0e9)
            
              PM2p5_SOA(i, j, k) =                                         &
                                   (((A_param_SOA * N_tot / 2.0)           &
                                   + (A_param_SOA * N_tot / 2.0)           &
                                   * erf((alog(d_cut2) - B_param_SOA)      &
                                   / (sqrt(2.0) * alog(sigma_SOA))))       &
                                   * rho_SOA * 1.0e9)
     
              PM10  (i, j, k) = PM10  (i, j, k) + PM10_SOA (i, j, k)
              PM2p5 (i, j, k) = PM2p5 (i, j, k) + PM2p5_SOA(i, j, k)
            END IF

            IF (L_seasalt_CCN .OR. L_USE_SEASALT_DIRECT .OR.               &
                L_USE_SEASALT_INDIRECT .OR. L_USE_SEASALT_AUTOCONV) THEN
              N_tot = sea_salt_film(i, j, k)

              PM10_SS(i, j, k) =                                           &
                                  (((A_param_SS_fi * N_tot / 2.0)          &
                                  + (A_param_SS_fi * N_tot / 2.0)          &
                                  * erf((alog(d_cut1) - B_param_SS_fi)     &
                                  / (sqrt(2.0) * alog(sigma_SS_fi))))      &
                                  * rho_SS * 1.0e9)

     
              PM2p5_SS(i, j, k) =                                          &
                                  (((A_param_SS_fi * N_tot / 2.0)          &
                                  + (A_param_SS_fi * N_tot / 2.0)          &
                                  * erf((alog(d_cut2) - B_param_SS_fi)     &
                                  / (sqrt(2.0) * alog(sigma_SS_fi))))      &
                                  * rho_SS * 1.0e9)

              N_tot = sea_salt_jet(i, j, k)
     
              PM10_SS(i, j, k) =  PM10_SS(i, j, k) +                       &
                                  (((A_param_SS_je * N_tot / 2.0)          &
                                  + (A_param_SS_je * N_tot / 2.0)          &
                                  * erf((alog(d_cut1) - B_param_SS_je)     &
                                  / (sqrt(2.0) * alog(sigma_SS_je))))      &
                                  * rho_SS * 1.0e9)
    
              PM2p5_SS(i, j, k) = PM2p5_SS(i, j, k) +                      &
                                  (((A_param_SS_je * N_tot / 2.0)          &
                                  + (A_param_SS_je * N_tot / 2.0)          &
                                  * erf((alog(d_cut2) - B_param_SS_je)     &
                                  / (sqrt(2.0) * alog(sigma_SS_je))))      &
                                  * rho_SS * 1.0e9)
     
              PM10  (i, j, k) = PM10  (i, j, k) + PM10_SS (i, j, k)
              PM2p5 (i, j, k) = PM2p5 (i, j, k) + PM2p5_SS(i, j, k)
            END IF
          
            IF (L_DUST) THEN    
              PM10_DUST(i, j, k) =                                         &
                                    ((DUST_1(i, j, k)                      &
                                    + DUST_2(i, j, k)                      &
                                    + DUST_3(i, j, k)                      &
                                    + DUST_4(i, j, k)                      &
                                    + (0.398317 * DUST_5(i, j, k)))        &
                                    * rho_air * 1.0e9)
     
              PM2p5_DUST(i, j, k) =                                        &
                                    ((DUST_1(i, j, k)                      &
                                    + DUST_2(i, j, k)                      &
                                    + DUST_3(i, j, k)                      &
                                    + (0.193941 * DUST_4(i, j, k)))        &
                                    * rho_air * 1.0e9)
     
              PM10  (i, j, k) = PM10  (i, j, k) + PM10_DUST (i, j, k)
              PM2p5 (i, j, k) = PM2p5 (i, j, k) + PM2p5_DUST(i, j, k)
            END IF

          END DO
        END DO
      END DO
      
     RETURN
     END SUBROUTINE PM10_PM2P5
!
