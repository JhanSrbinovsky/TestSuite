

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine to calculate the spectral albedo of the land surface using
! the two stream approach of Sellers, 1995.
!
!**********************************************************************
      SUBROUTINE ALBPFT (P_FIELD,LAND_FIELD,                            &
     &                       LAND_INDEX,TILE_INDEX,TILE_PTS,ILAYERS,    &
     &                       ALBSOIL,COSZ,LAI,ALB_TYPE,                 &
     &                       FAPAR_DIR,FAPAR_DIF,CAN_RAD_MOD)


      IMPLICIT NONE

!------------------------ nstypes.h ----------------------------------
!jhan:further renovation of ths file may be necessary params are dependent on dataset
!jhan: ALSO nstypes_cable.h should be unecessary nsoil/soil is only difference
      !--- Number of non-vegetation surface types
      Integer, Parameter :: NNVG  = 4

      !--- Number of plant functional types.
      Integer, Parameter :: NPFT  = 13
      
      !--- Number of surface types.
      Integer, Parameter :: NTYPE =17 
      
      !--- Index of the surface type 'Soil'
      !Integer, Parameter :: SOIL  = 16 
      !dhb599, 20110615: change made as per Peter Vohralik, item 1:
      Integer, Parameter :: SOIL  = 14

!--- Land surface types :
!--- original veg. tiles 
!     1 - Broadleaf Tree
!     2 - Needleleaf Tree
!     3 - C3 Grass
!     4 - C4 Grass
!     5 - Shrub
!--- for testing these tiles are set = 1:5 
!     6 - Broadleaf Tree
!     7 - Needleleaf Tree
!     8 - C3 Grass
!     9 - C4 Grass
!    10 - Shrub
!--- for testing these tiles are set = 0
!    11 - 0 
!    11 - 0
!    11 - 0
!--- original non-veg tiles moved to these indices
!     14 - Urban
!     15 - Water
!     16 - Soil
!     17 - Ice









! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER                                                           &
     & P_FIELD                                                          &
                                   ! Total number of grid points.
     &,LAND_FIELD                  ! Number of land points.

!   Array arguments with intent(in):
      INTEGER                                                           &
     & LAND_INDEX(LAND_FIELD)                                           &
                                   ! Index of land points.
     &,TILE_PTS(NTYPE)                                                  &
                                   ! Number of land points which
!                                  ! include the surface type.
     &,TILE_INDEX(LAND_FIELD,NTYPE)                                     &
                                   ! Indices of land points which
!                                  ! include the surface type.

     &,ILAYERS                                                          &
                                   ! IN Number of layers over which
                                   !    the PAR absorption profile is to
                                   !    be calculated.
     &,CAN_RAD_MOD                 ! IN Which canopy radiation model is 
                                   ! used (1/2)
      REAL                                                              &
     & ALBSOIL(LAND_FIELD)                                              &
                                   ! Soil albedo.
     &,COSZ(P_FIELD)                                                    &
                                   ! Cosine of the zenith angle.
     &,LAI(LAND_FIELD,NPFT)        ! Leaf area index.

!   Array arguments with intent(out):
      REAL                                                              &
     & ALB_TYPE(LAND_FIELD,NTYPE,4)! Albedos for surface types.
!                                  !   (*,*,1) - Direct beam visible
!                                  !   (*,*,2) - Diffuse visible
!                                  !   (*,*,3) - Direct beam near-IR
!                                  !   (*,*,4) - Diffuse near-IR

! Local arrays:
      REAL                                                              &
     & ALBUDIF(LAND_FIELD,2)                                            &
                                  ! Diffuse albedo of the underlying
!                                 ! surface.
     &,ALBUDIR(LAND_FIELD,2)                                            &
                                  ! Direct albedo of the underlying
!                                 ! surface.
     &,FAPAR_DIR(LAND_FIELD,NPFT,ILAYERS)                               &
!                                  ! OUT Profile of absorbed PAR
!                                  !     - Direct (fraction/LAI)
     &,FAPAR_DIF(LAND_FIELD,NPFT,ILAYERS)                               &
!                                  ! OUT Profile of absorbed PAR
!                                  !     - Diffuse (fraction/LAI)
     &,ALPL(2)                                                          &
                                  ! Leaf reflection coefficient.
     &,OM(2)                                                            &
                                  ! Leaf scattering coefficient.
     &,TAUL(2)                    ! Leaf transmission coefficient.

! Local scalars:
      REAL                                                              &
     & BETADIR                                                          &
                                  ! Upscatter parameter for direct beam.
     &,BETADIF                                                          &
                                  ! Upscatter parameter for diffuse beam
     &,COSZM                                                            &
                                  ! Mean value of COSZ.
     &,K                                                                &
                                  ! Optical depth per unit leaf area.
     &,G                                                                &
                                  ! Relative projected leaf area in
                                  ! direction cosz.
     &,SALB                                                             &
                                  ! Single scattering albedo.
     &,SQCOST                                                           &
                                  ! Cosine squared of the mean leaf
                                  ! angle to the horizontal.
     &,B,C,CA,D,F,H,U1                                                  &
                                  ! Miscellaneous variables from
     &,P1,P2,P3,P4,D1                                                   &
                                  ! Sellers (1985).
     &,H1,H2,H3,H7,H8                                                   &
                                  !
     &,S1,S2,SIG                  !

!-----------------------------------------------------------------------
! Additional work variables to calculate profiles of absorbed PAR
!-----------------------------------------------------------------------
      REAL                                                              &
     & DLAI                                                             &
                                      ! LAI increment.
     &,LA                                                               &
                                      ! Cumulative LAI through canopy.
     &,DRDIRD_DLAI,DRDIRU_DLAI                                          &
     &,DRDIFD_DLAI,DRDIFU_DLAI                                          &
     &,U2,U3,D2,H4,H5,H6,H9,H10
                                      ! Rate of change of fluxes with LAI
                                      ! (W/m2/LAI). 

      INTEGER                                                           &
     & BAND,I,J,L,N               ! Loop counters.

!-----------------------------------------------------------------------
!Functional Type dependent parameters
!-----------------------------------------------------------------------
      INTEGER                                                           &
     & C3(NPFT)                                                         &
                                  ! 1 for C3 Plants, 0 for C4 Plants.
     &,CROP(NPFT)                                                       &
                                  ! 1 for crop type, 0 for non-crop.
     &,ORIENT(NPFT)               ! 1 for horizontal, 0 for spherical.

      REAL                                                              &
     & ALPHA(NPFT)                                                      &
                                  ! Quantum efficiency
!                                ! (mol CO2/mol PAR photons).
     &,ALNIR(NPFT)                                                      &
                                  ! Leaf reflection coefficient for
!                                ! near infra-red.
     &,ALPAR(NPFT)                                                      &
                                  ! Leaf reflection coefficient for
!                                ! PAR.
     &,A_WL(NPFT)                                                       &
                                  ! Allometric coefficient relating
!                                ! the target woody biomass to
!                                ! the leaf area index (kg C/m2).
     &,A_WS(NPFT)                                                       &
                                  ! Woody biomass as a multiple of
!                                ! live stem biomass.
     &,B_WL(NPFT)                                                       &
                                  ! Allometric exponent relating
!                                ! the target woody biomass to
!                                ! the leaf area index.
     &,DGL_DM(NPFT)                                                     &
                                  ! Rate of change of leaf turnover
!                                ! rate with moisture availability.
     &,DGL_DT(NPFT)                                                     &
                                  ! Rate of change of leaf turnover
!                                ! rate with temperature (/K)
     &,DQCRIT(NPFT)                                                     &
                                  ! Critical humidity deficit
!                                ! (kg H2O/kg air).
     &,ETA_SL(NPFT)                                                     &
                                  ! Live stemwood coefficient
!                                ! (kg C/m/LAI).
     &,FSMC_OF(NPFT)                                                    &
                                  ! Moisture availability below
!                                ! which leaves are dropped.
     &,F0(NPFT)                                                         &
                                  ! CI/CA for DQ = 0.
     &,GLMIN(NPFT)                                                      &
                                  ! Minimum leaf conductance for H2O
     &,G_AREA(NPFT)                                                     &
                                  ! Disturbance rate (/360days).
     &,G_GROW(NPFT)                                                     &
                                  ! Rate of leaf growth (/360days).
     &,G_LEAF_0(NPFT)                                                   &
                                  ! Minimum turnover rate for leaves
!                                 ! (/360days).
     &,G_ROOT(NPFT)                                                     &
                                  ! Turnover rate for root biomass
!                                 ! (/360days).
     &,G_WOOD(NPFT)                                                     &
                                  ! Turnover rate for woody biomass
!                                 ! (/360days).
     &,KPAR(NPFT)                                                       &
                                  ! PAR Extinction coefficient
!                                ! (m2 leaf/m2 ground).
     &,LAI_MAX(NPFT)                                                    &
                                  ! Maximum projected LAI.
     &,LAI_MIN(NPFT)                                                    &
                                  ! Minimum projected LAI.
     &,NL0(NPFT)                                                        &
                                  ! Top leaf nitrogen concentration
!                                ! (kg N/kg C).
     &,NR_NL(NPFT)                                                      &
                                  ! Ratio of root nitrogen
!                                ! concentration to leaf
!                                ! nitrogen concentration.
     &,NS_NL(NPFT)                                                      &
                                  ! Ratio of stem nitrogen
!                                ! concentration to leaf
!                                ! nitrogen concentration.
     &,OMEGA(NPFT)                                                      &
                                  ! Leaf scattering coefficient
!                                ! for PAR.
     &,OMNIR(NPFT)                                                      &
                                  ! Leaf scattering coefficient for
!                                ! near infra-red.
     &,R_GROW(NPFT)                                                     &
                                  ! Growth respiration fraction.
     &,SIGL(NPFT)                                                       &
                                  ! Specific density of leaf carbon
!                                ! (kg C/m2 leaf).
     &,TLEAF_OF(NPFT)                                                   &
                                  ! Temperature below which leaves are
!                                ! dropped.
     &,TLOW(NPFT)                                                       &
                                  ! Lower temperature for
!                                ! photosynthesis (deg C)
     &,TUPP(NPFT)                 ! Upper temperature for
!                                ! photosynthesis (deg C)

!----------------------------------------------------------------------
!                       BT     NT    C3G    C4G     S
!----------------------------------------------------------------------
      DATA C3      /      1,     1,     1,     0,     1 /
      DATA CROP    /      0,     0,     1,     1,     0 /
      DATA ORIENT  /      0,     0,     0,     0,     0 /
      DATA ALPHA   /   0.08,  0.08,  0.08, 0.040,  0.08 /
      DATA ALNIR   /   0.45,  0.35,  0.58,  0.58,  0.58 /
      DATA ALPAR   /   0.10,  0.07,  0.10,  0.10,  0.10 /
      DATA A_WL    /   0.65,  0.65, 0.005, 0.005,  0.10 /
      DATA A_WS    /  10.00, 10.00,  1.00,  1.00, 10.00 /
      DATA B_WL    /  1.667, 1.667, 1.667, 1.667, 1.667 /
      DATA DGL_DM  /    0.0,   0.0,   0.0,   0.0,   0.0 /
      DATA DGL_DT  /    9.0,   9.0,   0.0,   0.0,   9.0 /
      DATA DQCRIT  /  0.090, 0.060, 0.100, 0.075, 0.100 /
      DATA ETA_SL  /   0.01,  0.01,  0.01,  0.01,  0.01 /
      DATA F0      /  0.875, 0.875, 0.900, 0.800, 0.900 /
      DATA FSMC_OF /   0.00,  0.00,  0.00,  0.00,  0.00 /
      DATA GLMIN   / 1.0E-6,1.0E-6,1.0E-6,1.0E-6,1.0E-6 /
      DATA G_AREA  /  0.005, 0.004,  0.25,  0.25,  0.05 /
      DATA G_GROW  /  20.00, 20.00, 20.00, 20.00, 20.00 /
      DATA G_LEAF_0/   0.25,  0.25,  0.25,  0.25,  0.25 /
      DATA G_ROOT  /   0.25,  0.25,  0.25,  0.25,  0.25 /
      DATA G_WOOD  /   0.01,  0.01,  0.20,  0.20,  0.05 /
      DATA KPAR    /   0.50,  0.50,  0.50,  0.50,  0.50 /
      DATA LAI_MAX /   9.00,  9.00,  4.00,  4.00,  4.00 /
      DATA LAI_MIN /   3.00,  3.00,  1.00,  1.00,  1.00 /
      DATA NL0     /  0.040, 0.030, 0.060, 0.030, 0.030 /
      DATA NR_NL   /   1.00,  1.00,  1.00,  1.00,  1.00 /
      DATA NS_NL   /   0.10,  0.10,  1.00,  1.00,  0.10 /
      DATA OMEGA   /   0.15,  0.15,  0.15,  0.17,  0.15 /
      DATA OMNIR   /   0.70,  0.45,  0.83,  0.83,  0.83 /
      DATA R_GROW  /   0.25,  0.25,  0.25,  0.25,  0.25 /
      DATA SIGL    / 0.0375,0.1000,0.0250,0.0500,0.0500 /
      DATA TLEAF_OF/ 273.15,243.15,258.15,258.15,243.15 /
      DATA TLOW    /    0.0,  -5.0,   0.0,  13.0,   0.0 /
      DATA TUPP    /   36.0,  31.0,  36.0,  45.0,  36.0 /

      DO L=1,LAND_FIELD
        ALBUDIF(L,1) = ALBSOIL(L)
        ALBUDIF(L,2) = ALBSOIL(L)
        ALBUDIR(L,1) = ALBSOIL(L)
        ALBUDIR(L,2) = ALBSOIL(L)
      ENDDO

      DO N=1,NPFT

        OM(1) = OMEGA(N)
        OM(2) = OMNIR(N)
        ALPL(1) = ALPAR(N)
        ALPL(2) = ALNIR(N)

        DO BAND=1,2  ! Visible and near-IR bands
          TAUL(BAND) = OM(BAND) - ALPL(BAND)
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            I = LAND_INDEX(L)
            IF (ORIENT(N) == 0) THEN
              SQCOST = 1./3.
              G = 0.5
              COSZM = 1.0
              IF (CAN_RAD_MOD == 2) THEN
                K = G / 0.01
              ENDIF
              SALB = 0.5*OM(BAND)
              IF (COSZ(I) > 0.01) THEN
                IF (CAN_RAD_MOD == 2) THEN
                  K = G / COSZ(I)
              ENDIF
              SALB = 0.5*OM(BAND) *                                     &
     &                 ( 1. - COSZ(I)*LOG((COSZ(I)+1.)/COSZ(I)) )
          ENDIF
          ELSEIF (ORIENT(N) == 1) THEN
            SQCOST = 1.
            IF (CAN_RAD_MOD == 1) THEN
              G = COSZ(I)
            ENDIF
            COSZM = 1.
            IF (CAN_RAD_MOD == 2) THEN
!             K = G / COSZ(I) = 1.0 because G = COSZ(I)
!                                  for horizontal leaves
              K = 1.0
            ENDIF

            SALB = OM(BAND)/4.0
          ENDIF
          IF (CAN_RAD_MOD == 1) THEN
            K = G / 0.01
            IF (COSZ(I) >  0.01) K = G / COSZ(I)
          ENDIF

          BETADIR = (1. + COSZM*K)/(OM(BAND)*COSZM*K)*SALB
          C = 0.5*( ALPL(BAND) + TAUL(BAND) +                           &
     &               (ALPL(BAND) - TAUL(BAND))*SQCOST )
          BETADIF = C / OM(BAND)
          B = 1. - (1. - BETADIF)*OM(BAND)
          D = OM(BAND)*COSZM*K*BETADIR
          F = OM(BAND)*COSZM*K*(1. - BETADIR)
          H = SQRT(B*B - C*C) / COSZM
          SIG = (COSZM*K)**2 + C*C - B*B
          U1 = B - C/ALBUDIF(L,BAND)
          CA = C*ALBUDIR(L,BAND)/ALBUDIF(L,BAND)
          S1 = EXP(-H*LAI(L,N))
          S2 = EXP(-K*LAI(L,N))
          P1 = B + COSZM*H
          P2 = B - COSZM*H
          P3 = B + COSZM*K
          P4 = B - COSZM*K
          D1 = P1*(U1 - COSZM*H)/S1 - P2*(U1 + COSZM*H)*S1
          H1 = -D*P4 - C*F
          H2 = ( (D - P3*H1/SIG) * (U1 - COSZM*H) / S1 -                &
     &             (D - CA - (U1 + COSZM*K)*H1/SIG)*P2*S2 ) / D1
          H3 = - ( (D - P3*H1/SIG) * (U1 + COSZM*H)*S1 -                &
     &               (D - CA - (U1 + COSZM*K)*H1/SIG)*P1*S2 ) / D1
          H7 = (C/D1)*(U1 - COSZM*H) / S1
          H8 = - (C/D1)*(U1 + COSZM*H) * S1
          ALB_TYPE(L,N,2*BAND-1) = H1/SIG + H2 + H3   ! Direct beam
          ALB_TYPE(L,N,2*BAND) = H7 + H8              ! Diffuse
          IF (CAN_RAD_MOD == 2) THEN
!-----------------------------------------------------------------------
! If required calculate the profile of absorbed PAR through the canopy
! of direct and diffuse beams (BEWARE: assumes PAR is band 1)
!-----------------------------------------------------------------------
            IF (BAND==1) THEN
              U2 = B-C*ALBUDIF(L,BAND)
              U3 = F+C*ALBUDIF(L,BAND)
              D2 = (U2+COSZM*H)/S1-(U2-COSZM*H)*S1
              H4 = -F*P3-C*D
              H5 = -1./D2*(H4/SIG*(U2+COSZM*H)/S1                       &
     &           +(U3-H4/SIG*(U2-COSZM*K))*S2)
              H6 = 1./D2*(H4/SIG*(U2-COSZM*H)*S1                        &
     &          +(U3-H4/SIG*(U2-COSZM*K))*S2)
              H9 = 1./D2*(U2+COSZM*H)/S1
              H10 = -1./D2*(U2-COSZM*H)*S1
!-----------------------------------------------------------------------
! Two-stream equations for direct and diffuse upward and downward beams:
!             RDIRU(I)=H1/SIG*EXP(-K*LA)+H2*EXP(-H*LA)+H3*EXP(H*LA)
!             RDIRD(I)=(H4/SIG+1)*EXP(-K*LA)+H5*EXP(-H*LA)+H6*EXP(H*LA)
!             RDIFU(I)=H7*EXP(-H*LA)+H8*EXP(H*LA)
!             RDIFD(I)=H9*EXP(-H*LA)+H10*EXP(H*LA)
! Differentiate these equations to calculate PAR absorption per unit
! LAI down through the canopy. Centre derivatives in the centre of each
! LAI layer.
!-----------------------------------------------------------------------
              DLAI=LAI(L,N)/REAL(ILAYERS)
              LA=0.5*DLAI
              DO I=1,ILAYERS

                DRDIRU_DLAI=-K*H1/SIG*EXP(-K*LA)-H*H2*EXP(-H*LA)        &
     &                    +H*H3*EXP(H*LA)
                DRDIRD_DLAI=-K*(H4/SIG+1)*EXP(-K*LA)-H*H5*EXP(-H*LA)    &
     &                    +H*H6*EXP(H*LA)
                DRDIFU_DLAI=-H*H7*EXP(-H*LA)+H*H8*EXP(H*LA)
                DRDIFD_DLAI=-H*H9*EXP(-H*LA)+H*H10*EXP(H*LA)
  
                FAPAR_DIR(L,N,I)=-DRDIRD_DLAI+DRDIRU_DLAI
                FAPAR_DIF(L,N,I)=-DRDIFD_DLAI+DRDIFU_DLAI
                LA=LA+DLAI


              ENDDO  !layers
            ENDIF    !abs par
           ENDIF !Can_rad_mod=2 loop

          ENDDO
        ENDDO

      ENDDO

      RETURN
      END SUBROUTINE ALBPFT
