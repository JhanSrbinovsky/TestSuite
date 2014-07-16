
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!**********************************************************************
! Subroutine to calculate gridbox mean values of surface conductance
! and carbon fluxes. Also returns net primary productivity, leaf
! turnover and wood respiration of each plant functional type for
! driving TRIFFID.
!
!
!
!    Programming standard:
!
!**********************************************************************
      SUBROUTINE PHYSIOL (ROW_LENGTH,ROWS,LAND_PTS,LAND_INDEX           &
     &,                   NSHYD,NTILES,TILE_PTS,TILE_INDEX              &
     &,                   DIM_CS1, DIM_CS2                              &
     &,                   CO2,CO2_3D,CO2_DIM_LEN                        &
     &,                   CO2_DIM_ROW,L_CO2_INTERACTIVE                 &
     &,                   L_TRIFFID, L_Q10                              &
     &,                   CAN_MODEL,CS,FRAC,HT,IPAR,LAI,PSTAR,Q1        &
     &,                   STHU,TIMESTEP,TSOIL,TSTAR_TILE                &
     &,                   V_CRIT,V_SAT,V_WILT,WIND,Z0_TILE,Z1           &
     &,                   CANHC_TILE,VFRAC_TILE,FLAKE,G_LEAF,GS,GS_TILE &
     &,                   GPP,GPP_FT,NPP,NPP_FT,RESP_P,RESP_P_FT        &
     &,                   RESP_S,RESP_W_FT,SMCT,WT_EXT_TILE,FSMC,WT_EXT &
     &,                   RA,F_ROOT,ALBSOIL,COS_ZENITH_ANGLE            &
     &,                   CAN_RAD_MOD, ILAYERS)

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



      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                                  ! IN Number of points on a row
     &,ROWS                                                             &
                                  ! IN Number of rows in a theta field
     &,LAND_PTS                                                         &
                                  ! IN Number of land points to be
!                                 !    processed.
     &,LAND_INDEX(LAND_PTS)                                             &
                                  ! IN Index of land points on the
!                                 !    P-grid.
     &,CO2_DIM_LEN                                                      &
                                  ! IN Length of a CO2 field row.
     &,CO2_DIM_ROW                                                      &
                                  ! IN Number of CO2 field rows.
     &,NSHYD                                                            &
                                  ! IN Number of soil moisture
!                                 !    levels.
     &,NTILES                                                           &
                                  ! IN Number of surface tiles.
     &,TILE_PTS(NTYPE)                                                  &
                                  ! IN Number of land points which
!                                 !    include the nth surface type.
     &,TILE_INDEX(LAND_PTS,NTYPE)                                       &
                                  ! IN Indices of land points which
!                                 !    include the nth surface type.
     &,CAN_MODEL                                                        &
                                  ! IN Swith for thermal vegetation
!                                 !    canopy
     &,DIM_CS1, DIM_CS2           ! soil carbon dimensions
     
      LOGICAL                                                           &
     &        L_CO2_INTERACTIVE                                         &
                                  ! switch for 3D CO2 field
     &,       L_TRIFFID                                                 &
                                  ! TRUE if using TRIFFID
     &,       L_Q10               ! TRUE if using Q10 for soil resp
     
      INTEGER                                                           &
     & CAN_RAD_MOD                                                      &
!                                  !Switch for canopy radiation model
     & ,ILAYERS
!                                  !No of layers in canopy radiation model
      REAL                                                              &
     & ALB_TYPE_DUMMY(LAND_PTS,NTYPE,4)                                 &
!                                 ! WORK Dummy argument for albedo
!                                 ! subroutine
     &,FAPAR_DIR(LAND_PTS,NPFT,ILAYERS)                                 &
!                                 ! WORK Profile of absorbed PAR -
!                                 ! Direct beam
     &,FAPAR_DIF(LAND_PTS,NPFT,ILAYERS)                                 &
!                                 ! WORK Profile of absorbed PAR -
!                                 ! Diffuse beam
     &,FAPAR_DIR_TOT(LAND_PTS,NPFT)                                     &
!                                 ! WORK Total absorbed PAR -
!                                 ! Direct beam
     &,FAPAR_DIF_TOT(LAND_PTS,NPFT)                                     &
!                                 ! WORK Total absorbed PAR -
!                                 ! Diffuse beam
     &,FAPARV(LAND_PTS,ILAYERS) ! WORK Profile of absorbed PAR.



      REAL                                                              &
     & CO2                                                              &
                                  ! IN Atmospheric CO2 concentration
     &,CO2_3D(CO2_DIM_LEN,CO2_DIM_ROW)                                  &
!                                 ! IN 3D atmos CO2 concentration
!                                 !    (kg CO2/kg air).
     &,CS(LAND_PTS,DIM_CS1)                                             &
                                 ! IN Soil carbon (kg C/m2).
     &,VEG_FRAC(DIM_CS2)                                                &
                                 ! WORK vegetated fraction of gridbox
     &,FRAC(LAND_PTS,NTYPE)                                             &
                                  ! IN Surface type fractions.
     &,HT(LAND_PTS,NPFT)                                                &
                                  ! IN Canopy height (m).
     &,IPAR(ROW_LENGTH,ROWS)                                            &
                                  ! IN Incident PAR (W/m2).
     &,LAI(LAND_PTS,NPFT)                                               &
                                  ! IN Leaf area index.
     &,PSTAR(ROW_LENGTH,ROWS)                                           &
                                  ! IN Surface pressure (Pa).
     &,Q1(ROW_LENGTH,ROWS)                                              &
                                  ! IN Specific humidity at level 1
!                                 !    (kg H2O/kg air).
     &,STHU(LAND_PTS,NSHYD)                                             &
                                  ! IN Soil moisture content in each
!                                 !    layer as a fraction of saturation
     &,TIMESTEP                                                         &
                                  ! IN Timestep (s).
     &,TSOIL(LAND_PTS)                                                  &
                                  ! IN Soil temperature (K).
     &,TSTAR_TILE(LAND_PTS,NTILES)                                      &
                                  ! IN Tile surface temperatures (K).
     &,V_CRIT(LAND_PTS)                                                 &
                                  ! IN Volumetric soil moisture
!                                 !    concentration above which
!                                 !    stomata are not sensitive
!                                 !    to soil water (m3 H2O/m3 soil).
     &,V_SAT(LAND_PTS)                                                  &
                                  ! IN Volumetric soil moisture
!                                 !    concentration at saturation
!                                 !    (m3 H2O/m3 soil).
     &,V_WILT(LAND_PTS)                                                 &
                                  ! IN Volumetric soil moisture
!                                 !    concentration below which
!                                 !    stomata close (m3 H2O/m3 soil).

     &,WIND(ROW_LENGTH,ROWS)                                            &
                                  ! IN Windspeed (m/s).
     &,Z0_TILE(LAND_PTS,NTILES)                                         &
                                  ! IN Tile roughness lengths (m).
     &,Z1(ROW_LENGTH,ROWS)                                              &
                                  ! IN Windspeed reference height(m).
     &,GS(LAND_PTS)                                                     &
                                  ! INOUT Gridbox mean surface
!                                 !       conductance (m/s).
     &,WT_EXT(LAND_PTS,NSHYD)                                           &
                                  ! OUT Gridbox-mean WT_EXT.
     &,RA(LAND_PTS)                                                     &
                                  ! OUT Aerodynamic resistance (s/m).
     &,F_ROOT(NSHYD)                                                    &
                                  ! OUT Fraction of roots in each soil
!                                 !     layer.
     &,ALBSOIL(LAND_PTS)                                                &
!                                  ! Soil albedo.
     &, COS_ZENITH_ANGLE(ROW_LENGTH, ROWS)
!                                  ! Cosine of the zenith angle

      REAL                                                              &
     & CANHC_TILE(LAND_PTS,NTILES)                                      &
                                  ! OUT Areal heat capacity of canopy
!                                 !     for land tiles (J/K/m2).
     &,FLAKE(LAND_PTS,NTILES)                                           &
                                  ! OUT Lake fraction.
     !kdcorbin, 11/10 - changed from NPFT
     &,G_LEAF(LAND_PTS,NTILES)                                          &
                                  ! OUT Leaf turnover rate (/360days).
     &,GS_TILE(LAND_PTS,NTILES)                                         &
                                  ! OUT Surface conductance for
!                                 !     land tiles (m/s).
     &,GPP(LAND_PTS)                                                    &
                                  ! OUT Gridbox mean gross primary
!                                 !     productivity (kg C/m2/s).
     !kdcorbin, 11/10 - changed from NPFT
     &,GPP_FT(LAND_PTS,NTILES)                                          &
                                  ! OUT Gross primary productivity
!                                 !     (kg C/m2/s).
     &,NPP(LAND_PTS)                                                    &
                                  ! OUT Gridbox mean net primary
!                                 !     productivity (kg C/m2/s).
     !kdcorbin, 11/10 - changed from NPFT
     &,NPP_FT(LAND_PTS,NTILES)                                          &
                                  ! OUT Net primary productivity
!                                 !     (kg C/m2/s).
     &,RESP_P(LAND_PTS)                                                 &
                                  ! OUT Gridbox mean plant respiration
!                                 !     (kg C/m2/s).
     !kdcorbin, 11/10 - changed from NPFT
     &,RESP_P_FT(LAND_PTS,NTILES)                                       &
                                  ! OUT Plant respiration (kg C/m2/s).
     &,RESP_S(LAND_PTS,DIM_CS1)                                         &
                                 ! OUT Soil respiration (kg C/m2/s).
     &,RESP_W_FT(LAND_PTS,NPFT)                                         &
                                  ! OUT Wood maintenance respiration
!                                 !     (kg C/m2/s).
     &,SMCT(LAND_PTS)                                                   &
                                  ! OUT Available moisture in the
!                                 !     soil profile (mm).
     &,VFRAC_TILE(LAND_PTS,NTILES)                                      &
                                  ! OUT Fractional canopy coverage for
!                                 !     land tiles.
     &,WT_EXT_TILE(LAND_PTS,NSHYD,NTILES)                               &
!                                 ! OUT Fraction of evapotranspiration
!                                 !     which is extracted from each
!                                 !     soil layer by each tile.
     &,FSMC(LAND_PTS,NPFT)        ! OUT Moisture availability factor.


!  External routines called :-
      EXTERNAL ROOT_FRAC,SMC_EXT,RAERO,SF_STOM,SOIL_EVAP,               &
     & LEAF_LIT,CANCAP,MICROBE


      REAL                                                              &
     & CANHC(LAND_PTS)                                                  &
                                  ! WORK Canopy heat capacity (J/K/m2).
     &,CH_TYPE(LAND_PTS,NTYPE)                                          &
                                  ! WORK CANHC for surface types.
     &,GSOIL(LAND_PTS)                                                  &
                                  ! WORK Bare soil conductance.
     &,GS_TYPE(LAND_PTS,NTYPE)                                          &
                                  ! WORK Conductance for surface types.
     &,PSTAR_LAND(LAND_PTS)                                             &
                                  ! WORK Surface pressure (Pa).
     &,RIB(ROW_LENGTH,ROWS)                                             &
                                  ! WORK Bulk Richardson Number.
     &,TSTAR(LAND_PTS)                                                  &
                                  ! WORK Surface temperature (K).
     &,VFRAC(LAND_PTS)                                                  &
                                  ! WORK Fractional canopy coverage.
     &,VF_TYPE(LAND_PTS,NTYPE)                                          &
                                  ! WORK VFRAC for surface types.
     &,WT_EXT_TYPE(LAND_PTS,NSHYD,NTYPE)                                &
!                                 ! WORK WT_EXT for surface types.
     &,Z0(LAND_PTS)               ! WORK Roughness length (m).

      INTEGER                                                           &
     & I,J,K,L,M,N                ! Loop indices

! PFTPARM defines Surface parameters for each Plant Functional Type

      REAL:: ALBSNC_MAX(NPFT)  ! Snow-covered albedo for large LAI.
      REAL:: ALBSNC_MIN(NPFT)  ! Snow-covered albedo for zero LAI.
      REAL:: ALBSNF_MAX(NPFT)  ! Snow-free albedo for large LAI.
      REAL:: DZ0V_DH(NPFT)     ! Rate of change of vegetation
                               ! roughness length with height.
      REAL:: CATCH0(NPFT)      ! Minimum canopy capacity (kg/m2).
      REAL:: DCATCH_DLAI(NPFT) ! Rate of change of canopy capacity
                               ! with LAI.
      REAL:: INFIL_F(NPFT)     ! Infiltration enhancement factor.
      REAL:: KEXT(NPFT)        ! Light extinction coefficient.
      REAL:: ROOTD_FT(NPFT)    ! Rootdepth (m).
      !----------------------------------------------------------------
      !                     BT    NT   C3G   C4G    S
      !----------------------------------------------------------------
      COMMON  /RUN_PFT/ALBSNC_MAX,ALBSNC_MIN,ALBSNF_MAX,DZ0V_DH,        &
     &  CATCH0,DCATCH_DLAI,INFIL_F,KEXT,ROOTD_FT

! PFTPARM end
! NVEGPARM start

! Surface and vegetation parameters
      REAL :: ALBSNC_NVG(NNVG)  ! Snow-covered albedo.
      REAL :: ALBSNF_NVG(NNVG)  ! Snow-free albedo.
      REAL :: CATCH_NVG(NNVG)   ! Canopy capacity (kg/m2).
      REAL :: GS_NVG(NNVG)          ! Surface conductance (m/s).
      REAL :: INFIL_NVG(NNVG)       ! Infiltration enhancement factor.
      REAL :: ROOTD_NVG(NNVG)   ! Rootdepth (m).
      REAL :: Z0_NVG(NNVG)          ! Roughness length (m).
      REAL :: CH_NVG(NNVG)         ! "Canopy" heat capacity (J/K/m2)
      REAL :: VF_NVG(NNVG)         ! Fractional "canopy" coverage

      COMMON  /RUN_BLVEG/ALBSNC_NVG,ALBSNF_NVG,CATCH_NVG,GS_NVG,        &
     &  INFIL_NVG,ROOTD_NVG,Z0_NVG,CH_NVG,VF_NVG

! NVEGPARM end
! C_DENSTY for subroutine SF_EXCH
      REAL,PARAMETER:: RHOSEA = 1026.0 ! density of sea water (kg/m3)
      REAL,PARAMETER:: RHO_WATER = 1000.0! density of pure water (kg/m3)
! C_DENSTY end
!---Soil layer thicknesses (m)
!---6 layers => CABLE else revert to MOSES
   REAL,PARAMETER:: DZSOIL(6) =(/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/)

      REAL                                                              &
     & DIFF_FRAC(ROW_LENGTH*ROWS)
!                             ! Fraction of diffuse radiation
!
!-----------------------------------------------------------------------
! Parameters
!-----------------------------------------------------------------------
      REAL, PARAMETER          ::   LAI_MAX = 30.0  ! Max allowable LAI
      REAL, PARAMETER          ::   GMIN = 1.0E-8   ! Min allowable GSOIL

!-----------------------------------------------------------------------
! Initialisations
!-----------------------------------------------------------------------
      IF (CAN_RAD_MOD == 2) THEN
      DO L=1,ROW_LENGTH*ROWS
         DIFF_FRAC(L) = 0.0
      ENDDO
      ENDIF
      DO K=1,NSHYD
        F_ROOT(K)=0.0
        DO L=1,LAND_PTS
          WT_EXT(L,K)=0.0
          DO N=1,NTYPE
            WT_EXT_TYPE(L,K,N)=0.0
          ENDDO
        ENDDO
      ENDDO
      IF (CAN_RAD_MOD == 2) THEN
      DO L=1,LAND_PTS
       DO N=1,NTYPE
        DO K=1,4
           ALB_TYPE_DUMMY(L,N,K)=0.0
        ENDDO
       ENDDO
      ENDDO

      DO L=1,LAND_PTS
       DO N=1,NPFT

! Bodge to stop TRIFFID giving a huge LAI at a small
! number of points
        LAI(L,N) = MIN(LAI(L,N),LAI_MAX)

        DO K=1,ILAYERS
          FAPAR_DIR(L,N,K)=0.0
          FAPAR_DIF(L,N,K)=0.0
        ENDDO
       ENDDO
      ENDDO
      ENDIF


      DO L=1,LAND_PTS
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        PSTAR_LAND(L) = PSTAR(I,J)
      ENDDO

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        RIB(I,J)=0.0
        ENDDO
      ENDDO

      DO N=1,NPFT
        DO L=1,LAND_PTS
          G_LEAF(L,N)=0.0
          GPP_FT(L,N)=0.0
          NPP_FT(L,N)=0.0
          RESP_P_FT(L,N)=0.0
          RESP_W_FT(L,N)=0.0
          FSMC(L,N)=0.0
        ENDDO
      ENDDO

      DO N=1,NTYPE
        DO L=1,LAND_PTS
          GS_TYPE(L,N)=GS(L)
        ENDDO
      ENDDO

      DO L=1,LAND_PTS
        GPP(L)=0.0
        NPP(L)=0.0
        RESP_P(L)=0.0
        SMCT(L)=0.0
        GS(L)=0.0
        RA(L)=0.0
        CANHC(L)=0.0
        VFRAC(L)=0.0
      ENDDO
      DO N=1,DIM_CS1
        DO L=1,LAND_PTS
          RESP_S(L,N)=0.0
        ENDDO
      ENDDO

      DO L=1,LAND_PTS
        GSOIL(L) = 0.
        IF (V_CRIT(L) > 0.) THEN
          GSOIL(L) = GS_NVG(SOIL-NPFT)*(STHU(L,1)*V_SAT(L)/V_CRIT(L))**2
          IF (CAN_RAD_MOD == 2) THEN
            GSOIL(L) = MAX(GSOIL(L),GMIN)
          ENDIF
        ENDIF
      ENDDO

!-----------------------------------------------------------------------
! Calculate light absorption by the plant canopy
!-----------------------------------------------------------------------
      IF (CAN_RAD_MOD == 2) THEN
        CALL ALBPFT       (ROW_LENGTH*ROWS,LAND_PTS,                    &
     &                     LAND_INDEX,TILE_INDEX,TILE_PTS,ILAYERS,      &
     &                     ALBSOIL,COS_ZENITH_ANGLE,LAI,ALB_TYPE_DUMMY, &
     &                     FAPAR_DIR,FAPAR_DIF,CAN_RAD_MOD)

      ENDIF


!-----------------------------------------------------------------------
! Loop over Plant Functional Types to calculate the available moisture
! and the values of canopy conductance, the carbon fluxes and the leaf
! turnover rate
!-----------------------------------------------------------------------
      DO N=1,NPFT

        IF (NTILES == 1) THEN
          DO L=1,LAND_PTS
            TSTAR(L) = TSTAR_TILE(L,1)
            Z0(L) = Z0_TILE(L,1)
          ENDDO
        ELSE
          DO L=1,LAND_PTS
            TSTAR(L) = TSTAR_TILE(L,N)
            Z0(L) = Z0_TILE(L,N)
          ENDDO
        ENDIF

! DEPENDS ON: root_frac
        CALL ROOT_FRAC(NSHYD,DZSOIL,ROOTD_FT(N),F_ROOT)

! DEPENDS ON: smc_ext
        CALL SMC_EXT (LAND_PTS,NSHYD,TILE_PTS(N),TILE_INDEX(1,N)        &
     &,               F_ROOT ,STHU,V_CRIT,V_SAT,V_WILT                  &
     &,               WT_EXT_TYPE(1,1,N),FSMC(1,N))

! DEPENDS ON: raero
        CALL RAERO (ROW_LENGTH,ROWS,LAND_PTS,LAND_INDEX                 &
     &,             TILE_PTS(N),TILE_INDEX(1,N)                         &
     &,             RIB,WIND,Z0,Z0,Z1,RA)
!-----------------------------------------------------------------------
! Calculate light absorption by the plant canopy
!-----------------------------------------------------------------------

        IF (CAN_RAD_MOD == 2) THEN

        DO M=1,ILAYERS
          DO L=1,LAND_PTS
            FAPARV(L,M)=0.0
          ENDDO
        ENDDO

          DO L=1,LAND_PTS
            FAPAR_DIF_TOT(L,N)=0.0
            FAPAR_DIR_TOT(L,N)=0.0
          ENDDO

          DO M=1,ILAYERS
            DO K=1,TILE_PTS(N)
              L = TILE_INDEX(K,N)
              I = LAND_INDEX(L)

              FAPAR_DIR_TOT(L,N)=FAPAR_DIR_TOT(L,N)                     &
     &                           +FAPAR_DIR(L,N,M)
              FAPAR_DIF_TOT(L,N)=FAPAR_DIF_TOT(L,N)                     &
     &                           +FAPAR_DIF(L,N,M)

              FAPARV(L,M) = (1-DIFF_FRAC(I))*FAPAR_DIR(L,N,M)           &
     &                    + DIFF_FRAC(I)*FAPAR_DIF(L,N,M)
            ENDDO
          ENDDO
        ENDIF



! DEPENDS ON: sf_stom
        CALL SF_STOM (ROW_LENGTH,ROWS,LAND_PTS,LAND_INDEX               &
     &,               TILE_PTS(N),TILE_INDEX(1,N),N                     &
     &,               CO2,CO2_3D,CO2_DIM_LEN                            &
     &,               CO2_DIM_ROW,L_CO2_INTERACTIVE                     &
     &,               FSMC(1,N),HT(1,N),IPAR,LAI(1,N),PSTAR_LAND        &
     &,               Q1,RA,TSTAR                                       &
     &,               CAN_RAD_MOD, ILAYERS, FAPARV                      &
     &,               GPP_FT(1,N),NPP_FT(1,N),RESP_P_FT(1,N)            &
     &,               RESP_W_FT(1,N),GS_TYPE(1,N))

! DEPENDS ON: soil_evap
        CALL SOIL_EVAP (LAND_PTS,NSHYD,TILE_PTS(N),TILE_INDEX(1,N)      &
     &,                 GSOIL,LAI(1,N),GS_TYPE(1,N),WT_EXT_TYPE(1,1,N))

! DEPENDS ON: leaf_lit
        CALL LEAF_LIT (LAND_PTS,TILE_PTS(N),TILE_INDEX(1,N)             &
     &,                N,FSMC(1,N),TSTAR,G_LEAF(1,N))

! DEPENDS ON: cancap
        CALL CANCAP (LAND_PTS,TILE_PTS(N),TILE_INDEX(1,N),CAN_MODEL,N   &
     &,              HT(1,N),LAI(1,N),CH_TYPE(1,N),VF_TYPE(1,N))

      ENDDO

!----------------------------------------------------------------------
! Non-vegetated surface types
!----------------------------------------------------------------------
      DO N=NPFT+1,NTYPE
        DO M=1,TILE_PTS(N)
          L=TILE_INDEX(M,N)
          GS_TYPE(L,N) = GS_NVG(N-NPFT)
          DO K=1,NSHYD
            WT_EXT_TYPE(L,K,N) = 0.
          ENDDO
        ENDDO
      ENDDO

! Copy soil conductance and add bare soil fraction to extraction from
! surface layer
      N = SOIL
      DO M=1,TILE_PTS(N)
        L=TILE_INDEX(M,N)
        GS_TYPE(L,N) = GSOIL(L)
        WT_EXT_TYPE(L,1,N) = 1.
      ENDDO

!----------------------------------------------------------------------
! Canopy heat capacity and coverage for non-vegetated surfaces
!----------------------------------------------------------------------
      DO N=NPFT+1,NTYPE
        DO M=1,TILE_PTS(N)
          L=TILE_INDEX(M,N)
          CH_TYPE(L,N) = CH_NVG(N-NPFT)
          VF_TYPE(L,N) = VF_NVG(N-NPFT)
        ENDDO
      ENDDO

!----------------------------------------------------------------------
! Calculate the rate of soil respiration
!----------------------------------------------------------------------
! set VEG_FRAC according to whether it is full or dummy field
      IF (L_TRIFFID) THEN
        DO L=1,LAND_PTS
          VEG_FRAC(L) = FRAC(L,1) + FRAC(L,2) + FRAC(L,3)               &
     &                            + FRAC(L,4) + FRAC(L,5)
        ENDDO
      ELSE
        VEG_FRAC(1)=0.0
      ENDIF
! DEPENDS ON: microbe
      CALL MICROBE (LAND_PTS,DIM_CS1,DIM_CS2,L_TRIFFID,L_Q10,CS,        &
     &              STHU,V_SAT,V_WILT,TSOIL,RESP_S,VEG_FRAC)

!----------------------------------------------------------------------
! Form gridbox mean values
!----------------------------------------------------------------------

      DO N=1,NTYPE
        DO M=1,TILE_PTS(N)
          L=TILE_INDEX(M,N)
          GS(L) = GS(L) + FRAC(L,N)*GS_TYPE(L,N)
        ENDDO
      ENDDO

      IF (NTILES == 1) THEN
        DO N=1,NTYPE
          DO M=1,TILE_PTS(N)
            L=TILE_INDEX(M,N)
            CANHC(L) = CANHC(L) + FRAC(L,N)*CH_TYPE(L,N)
            VFRAC(L) = VFRAC(L) + FRAC(L,N)*VF_TYPE(L,N)
            DO K=1,NSHYD
              WT_EXT(L,K) = WT_EXT(L,K) + FRAC(L,N)*WT_EXT_TYPE(L,K,N)
            ENDDO
          ENDDO
        ENDDO
        DO L=1,LAND_PTS
          FLAKE(L,1) = FRAC(L,7)
          GS_TILE(L,1) = 0.
          IF (FLAKE(L,1) <  1.)                                         &
     &      GS_TILE(L,1) = GS(L) / (1. - FLAKE(L,1))
          CANHC_TILE(L,1) = CANHC(L)
          VFRAC_TILE(L,1) = VFRAC(L)
          DO K=1,NSHYD
            WT_EXT_TILE(L,K,1) = WT_EXT(L,K)
          ENDDO
        ENDDO
      ELSE
        GS_TILE(:,:)=0.0
        CANHC_TILE(:,:)=0.0
        VFRAC_TILE(:,:)=0.0
        DO N=1,NTYPE
          DO M=1,TILE_PTS(N)
            L=TILE_INDEX(M,N)
            FLAKE(L,N) = 0.
            GS_TILE(L,N) = GS_TYPE(L,N)
            CANHC_TILE(L,N) = CH_TYPE(L,N)
            VFRAC_TILE(L,N) = VF_TYPE(L,N)
            DO K=1,NSHYD
              WT_EXT_TILE(L,K,N) = WT_EXT_TYPE(L,K,N)
            ENDDO
          ENDDO
        ENDDO
        N = 7    ! Lake tile
        DO M=1,TILE_PTS(N)
          L=TILE_INDEX(M,N)
          FLAKE(L,N) = 1.
        ENDDO
      ENDIF

      DO N=1,NPFT
        DO M=1,TILE_PTS(N)
          L=TILE_INDEX(M,N)

          GPP(L)=GPP(L)+FRAC(L,N)*GPP_FT(L,N)
          NPP(L)=NPP(L)+FRAC(L,N)*NPP_FT(L,N)
          RESP_P(L)=RESP_P(L)+FRAC(L,N)*RESP_P_FT(L,N)

        ENDDO
      ENDDO

!----------------------------------------------------------------------
! Diagnose the available moisture in the soil profile
!----------------------------------------------------------------------
      DO N=1,NSHYD
        DO L=1,LAND_PTS
          SMCT(L) = SMCT(L) + MAX( 0. ,                                 &
     &             RHO_WATER*DZSOIL(N)*(STHU(L,N)*V_SAT(L)-V_WILT(L)))
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE PHYSIOL
