
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE HYDROL-------------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE HYDROL (                                               &
     &                   LICE_PTS,LICE_INDEX,SOIL_PTS,SOIL_INDEX,       &
     &                   NPNTS,NSHYD,B,CAN_CPY,CON_RAIN,CON_SNOW,       &
     &                   E_CANOPY,EXT,HCAP,HCON,LS_RAIN,                &
     &                   LS_SNOW,SATCON,SATHH,                          &
     &                   SURF_HT_FLUX,TIMESTEP,                         &
     &                   V_SAT,V_WILT,RGRAIN,L_SNOW_ALBEDO,             &
     &                   CAN_WCNT,STF_HF_SNOW_MELT,                     &
     &                   STF_SUB_SURF_ROFF,SMCL,STHF,STHU,TSOIL,        &
     &                   CAN_WCNT_GB,HF_SNOW_MELT,INFIL,SMC,            &
     &                   SNOW_MELT,SNOMLT_SUB_HTF,                      &
     &                   SUB_SURF_ROFF,SURF_ROFF,TOT_TFALL,             &

! add new inland basin variable
     &                   inlandout_atm,L_INLAND,                        &

! Additional variables for MOSES II
     &                   NTILES,TILE_PTS,TILE_INDEX,CAN_MODEL,          &
     &                   INFIL_TILE,                                    &
     &                   MELT_TILE,TSTAR_TILE,FRAC,                     &
     &                   SNOW_TILE,LYING_SNOW,                          &
     &                   CATCH_SNOW,SNOW_GRND,                          &

! Additional variables required for large-scale hydrology:
     &                   L_TOP,L_PDM,FEXP,GAMTOT,TI_MEAN,TI_SIG,CS,     &
     &                   DUN_ROFF,DRAIN,FSAT,FWETL,QBASE,QBASE_ZW,      &
     &                   ZW,STHZW,A_FSAT,C_FSAT,A_FWET,C_FWET,          &
     &                   FCH4_WETL,L_SOIL_SAT_DOWN,                     &
! sxy
       L_TILE_PTS,ROW_LENGTH,ROWS,N_ROWS,                               &
!       ROW_LENGTH,ROWS,N_ROWS,                               &
       LAND_INDEX,SMVCCL,SMVCST,SMVCWT,TSTAR,                           &
       FQW_1,FTL_1,FTL_TILE,LE_TILE,RADNET_TILE,FQW_TILE,               &
       SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L,                   &
       SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE,             &
       TSOIL_TILE,T_SURF_TILE,HCONS,                                    &
       SOIL_TYPE,VEG_TYPE,ISNOW_FLG3L,total_nsteps,                     &
       FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB,                         &
       TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB,                     &
       SURF_HTF_CAB,SURF_HT_FLUX_LAND,                                  &
       ESOIL_TILE,EI_TILE,                                              &
       SNAGE_TILE,GFLUX_TILE,SGFLUX_TILE,                               &
! EAK
     &    l_cable, WB_LAKE,                                             &

     & LTIMER                                                           &
     & )
! sxy        

      IMPLICIT NONE
!
! Description:
!     Surface hydrology module which also updates the
!     sub-surface temperatures. Includes soil water phase
!     changes and the effect of soil water and ice on the
!     thermal and hydraulic characteristics of the soil.
!     This version is for use with MOSES II land surface
!     scheme.
!
! Documentation : UM Documentation Paper 25
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.2    15/11/00   New Deck         M. Best
!  5.3    08/06/01   Declare ntiles before it is used.  A van der Wal
!  5.4    11/04/02   Canopy snow model for needleleaf tree tile.
!                    R. Essery
!  5.4    04/09/02   If there are no soil points, set sub-surface
!                    runoff to zero. C. Durman.
!  5.5  12/02/03   Code added for large-scale hydrology.
!                                          Nic Gedney.
!  6.1  17/08/04  Add SSFM code                           Ian Pearman
!  6.2  02/02/06  Remove SSFM. P.Selwood
!   6.2   21/2/06  Re-route outflow from inland basins to soil moisture
!                  P. Falloon

! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! System component covered: P25
! System Task: P25
!

! Global variables:
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
!---Soil layer thicknesses (m)
!---6 layers => CABLE else revert to MOSES
   REAL,PARAMETER:: DZSOIL(6) =(/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/)


! Subroutine arguments
!   Scalar arguments with intent(IN) :
      INTEGER                                                           &
     & LICE_PTS                                                         &
                           ! IN Number of land ice points.
     &,NPNTS                                                            &
                           ! IN Number of gridpoints.
     &,NSHYD                                                            &
                           ! IN Number of soil moisture levels.
     &,SOIL_PTS                                                         &
                           ! IN Number of soil points.
     &, NTILES                                                          &
     &,CAN_MODEL           ! IN Switch for thermal vegetation canopy.

      REAL                                                              &
     & TIMESTEP            ! IN Model timestep (s).

      LOGICAL LTIMER       ! Logical switch for TIMER diags

! EAK
      LOGICAL  l_cable

      LOGICAL                                                           &
     & STF_HF_SNOW_MELT                                                 &
                           ! IN Stash flag for snowmelt heat flux.
     &,STF_SUB_SURF_ROFF                                                &
                           ! IN Stash flag for sub-surface runoff.
     &,L_SNOW_ALBEDO                                                    &
                           ! IN Flag for prognostic snow albedo
     &,L_TOP                                                            &
                   ! IN Flag for TOPMODEL-based hydrology.
     &,L_PDM                                                            &
                   ! IN Flag for PDM hydrology.
     &,L_SOIL_SAT_DOWN
                   ! IN Direction of super_saturated soil moisture


!   Array arguments with intent(IN) :
      INTEGER                                                           &
     & LICE_INDEX(NPNTS)                                                &
                           ! IN Array of land ice points.
     &,SOIL_INDEX(NPNTS)   ! IN Array of soil points.

      REAL                                                              &
     & B(NPNTS)                                                         &
                           ! IN Clapp-Hornberger exponent.
     &,CAN_CPY(NPNTS,NTILES)                                            &
                            !IN Canopy/surface capacity of
!                          !    land tiles (kg/m2).
     &,CATCH_SNOW(NPNTS)                                                &
                           ! IN Coniferous canopy snow capacity (kg/m2)
     &,CON_RAIN(NPNTS)                                                  &
                           ! IN Convective rain (kg/m2/s).
     &,CON_SNOW(NPNTS)                                                  &
                           ! IN Convective snowfall (kg/m2/s).
     &,E_CANOPY(NPNTS,NTILES)                                           &
!                          ! IN Canopy evaporation from
!                          !    land tiles (kg/m2/s).
     &,EXT(NPNTS,NSHYD)                                                 &
                           ! IN Extraction of water from each soil
!                          !    layer (kg/m2/s).
     &,HCAP(NPNTS)                                                      &
                           ! IN Soil heat capacity (J/K/m3).
     &,HCON(NPNTS)                                                      &
                           ! IN Soil thermal conductivity (W/m/K).
     &,LS_RAIN(NPNTS)                                                   &
                           ! IN Large-scale rain (kg/m2/s).
     &,LS_SNOW(NPNTS)                                                   &
                           ! IN Large-scale snowfall (kg/m2/s).
     &,SATCON(NPNTS)                                                    &
                           ! IN Saturated hydraulic conductivity
!                          !    (kg/m2/s).
     &,SATHH(NPNTS)                                                     &
                           ! IN Saturated soil water pressure (m).
     &,SURF_HT_FLUX(NPNTS)                                              &
                           ! IN Net downward surface heat flux (W/m2)
     &,V_SAT(NPNTS)                                                     &
                           ! IN Volumetric soil moisture
!                          !    concentration at saturation
!                          !    (m3 H2O/m3 soil).
     &,V_WILT(NPNTS)                                                    &
                           ! IN Volumetric soil moisture
!                          !    concentration below which
!                          !    stomata close (m3 H2O/m3 soil).
     &,FEXP(NPNTS)                                                      &
                           ! IN Decay factor in Sat. Conductivity
!                          !    in water table layer.
     &,GAMTOT(NPNTS)                                                    &
                           ! IN Integrated complete Gamma function.
     &,TI_MEAN(NPNTS)                                                   &
                           ! IN Mean topographic index.
     &,TI_SIG(NPNTS)                                                    &
                           ! IN Standard dev. of topographic index.
     &,CS(NPNTS)                                                        &
                           ! IN Soil carbon (kg C/m2).
     &,A_FSAT(NPNTS)                                                    &
                           ! IN Fitting parameter for Fsat in LSH model
     &,C_FSAT(NPNTS)                                                    &
                           ! IN Fitting parameter for Fsat in LSH model
     &,A_FWET(NPNTS)                                                    &
                           ! IN Fitting parameter for Fwet in LSH model
     &,C_FWET(NPNTS)                                                    &
                           ! IN Fitting parameter for Fwet in LSH model

     &,RGRAIN(NPNTS,NTILES)                                             &
                           ! INOUT Snow grain size (microns).
     &,CAN_WCNT(NPNTS,NTILES)
!                          ! INOUT Canopy water content for
!                          !       land tiles (kg/m2).
!
!   Array arguments with intent(INOUT) :
!
      REAL                                                              &
     & SMCL(NPNTS,NSHYD)                                                &
                           ! INOUT Soil moisture content of each
!                          !       layer (kg/m2).
     &,STHF(NPNTS,NSHYD)                                                &
                           ! INOUT Frozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
     &,STHU(NPNTS,NSHYD)                                                &
                           ! INOUT Unfrozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
     &,TSOIL(NPNTS,NSHYD)                                               &
                           ! INOUT Sub-surface temperatures (K).
     &,SNOW_GRND(NPNTS)                                                 &
                           ! INOUT Snow below canopy (kg/m2).
     &,FSAT(NPNTS)                                                      &
                            ! INOUT Surface saturation fraction.
     &,FWETL(NPNTS)                                                     &
                            ! INOUT Wetland fraction.
     &,ZW(NPNTS)                                                        &
                            ! INOUT Water table depth (m).
     &,STHZW(NPNTS)         ! INOUT soil moist fract. in deep-zw layer.



!   Array arguments with intent(OUT) :
      REAL                                                              &
     & CAN_WCNT_GB(NPNTS)                                               &
                            ! OUT Gridbox canopy water content (kg/m2).
     &,HF_SNOW_MELT(NPNTS)                                              &
                            ! OUT Gridbox snowmelt heat flux (W/m2).
     &,INFIL(NPNTS)                                                     &
                            ! OUT Maximum surface infiltration
!                           !     rate (kg/m2/s).
     &,SMC(NPNTS)                                                       &
                            ! OUT Soil moisture in the top metre (kg/m2)
     &,SNOW_MELT(NPNTS)                                                 &
                            ! OUT Snowmelt (kg/m2/s).
     &,SNOMLT_SUB_HTF(NPNTS)                                            &
                            ! OUT Sub-surface snowmelt heat
!                           !     flux (W/m2).
     &,WB_LAKE(NPNTS,NTILES)                                            &
     &,SUB_SURF_ROFF(NPNTS)                                             &
                            ! OUT Sub-surface runoff (kg/m2/s).
     &,SURF_ROFF(NPNTS)                                                 &
                            ! OUT Surface runoff (kg/m2/s).
     &,TOT_TFALL(NPNTS)                                                 &
                            ! OUT Total throughfall (kg/m2/s).
     &,DUN_ROFF(NPNTS)                                                  &
                            ! OUT Dunne part of sfc runoff (kg/m2/s).
     &,QBASE(NPNTS)                                                     &
                            ! OUT Base flow (kg/m2/s).
     &,QBASE_ZW(NPNTS)                                                  &
                            ! OUT Base flow from ZW layer (kg/m2/s).
     &,QBASE_L(NPNTS,NSHYD+1)                                           &
!                           ! OUT Base flow from each level (kg/m2/s).
     &,DRAIN(NPNTS)                                                     & 
                            ! OUT Drainage out of nshyd'th level (kg/m2/s).
     &,FCH4_WETL(NPNTS)     ! OUT Wetland methane flux. (kg C/m2/s).


! Additional variables for MOSES II
      INTEGER                                                           &
     & TILE_PTS(NTILES)                                                 &
                           ! IN Number of tile points.
     &,TILE_INDEX(NPNTS,NTILES)
!                          ! IN Index of tile points.

      REAL                                                              &
     & INFIL_TILE(NPNTS,NTILES)                                         &
!                          ! IN Maximum surface infiltration
     &,MELT_TILE(NPNTS,NTILES)                                          &
!                          ! IN Snowmelt on tiles (kg/m2/s).
     &,TSTAR_TILE(NPNTS,NTILES)                                         &
!                          ! IN Surface temperature (K).
     &,FRAC(NPNTS,NTILES)                                               &
                           ! IN Tile fractions.
     &,SNOW_TILE(NPNTS,NTILES)                                          &
!                          ! INOUT Snowmass on tiles (kg/m2).
     &,LYING_SNOW(NPNTS)                                                &
                           ! OUT Gridbox snowmass (kg/m2).
! Declare variable for inland basin outflow
     &,INLANDOUT_ATM(NPNTS)            ! IN TRIP INLAND BASIN
!                       OUTFLOW FOR LAND POINTS ONLY,kg/m2/s=mm
      LOGICAL                                                           &
     & L_INLAND                   ! IN True if re-routing inland
                                  !   basin flow to soil moisture       
!sxy
      Integer                                                           &
     &  ROW_LENGTH                                                      &
     &, ROWS                                                            &
     &, N_ROWS                                                          &
     &, LAND_INDEX (NPNTS)                                              &
     &, total_nsteps                                                    &
                                ! Total number of steps in run
     &, SOIL_TYPE(row_length,rows)                                      &
     &, VEG_TYPE(row_length,rows)                                       &
     &, ISNOW_FLG3L(NPNTS,NTILES)

      Real                                                              &
     &  SMVCCL (NPNTS)                                                  &
                             ! soil/qrparm.soil.crit
     &, SMVCWT (NPNTS)                                                  &
                             ! soil/qrparm.soil.wilt
     &, SMVCST (NPNTS)                                                  &
                            ! soil/qrparm.soil.satn
     &, TSTAR(row_length, rows)                                         &
                             ! Surface temperature (K)
     &, FQW_1(row_length, rows)                                         &
     &, FTL_1(row_length, rows)                                         &
     &, FQW_TILE(NPNTS,NTILES)                                          &
     &, FTL_TILE(NPNTS,NTILES)                                          &
     &, LE_TILE(NPNTS,NTILES)                                           &
     &, RADNET_TILE(NPNTS,NTILES)                                       &
     &, SNOW_DEPTH3L(NPNTS,NTILES,3)                          &
     &, SNOW_MASS3L(NPNTS,NTILES,3)                           &
     &, SNOW_COND(NPNTS,NTILES,3)                             &
     &, SNOW_TMP3L(NPNTS,NTILES,3)                            &
     &, SNOW_RHO3L(NPNTS,NTILES,3)                            &
     &, SNOW_RHO1L(NPNTS,NTILES)                              &
     &, SNAGE_TILE(NPNTS,NTILES)                              &
     &, SMCL_TILE(NPNTS,NTILES,NSHYD)                         &
     &, STHU_TILE(NPNTS,NTILES,NSHYD)                         &
     &, STHF_TILE(NPNTS,NTILES,NSHYD)                         &
     &, TSOIL_TILE(NPNTS,NTILES,NSHYD)                        &
     &, T_SURF_TILE(NPNTS,NTILES)                             &
     &, GFLUX_TILE(NPNTS,NTILES)                              &
     &, SGFLUX_TILE(NPNTS,NTILES)                             &
     &, HCONS(NPNTS)                                          &
     &, FTL_TILE_CAB(NPNTS,NTILES)                                      &
     &, FTL_CAB(NPNTS)                                                  &
     &, LE_TILE_CAB(NPNTS,NTILES)                                       &
     &, LE_CAB(NPNTS)                                                   &
     &, TSTAR_TILE_CAB(NPNTS,NTILES)                                    &
     &, TSTAR_CAB(NPNTS)                                                &
     &, SMCL_CAB(NPNTS,NSHYD)                                           &
     &, TSOIL_CAB(NPNTS,NSHYD)                                          &
     &, SURF_HTF_CAB(NPNTS)                                             &
     &, SURF_HT_FLUX_LAND(row_length,rows)                              &
     &, ESOIL_TILE(NPNTS,NTILES)                                        &
                                ! Evaporation from bare soil (kg/m2)
     &, EI_TILE(NPNTS,NTILES)                                           

     Logical                                                            &
     &  L_TILE_PTS(NPNTS,NTILES)



! Local scalars:
      INTEGER                                                           &
     & I,J                                                              &
                            ! WORK Loop counters.
     &,N                    ! WORK Tile loop counter.

! Local arrays:

      REAL                                                              &
     & DSMC_DT(NPNTS)                                                   &
                            ! WORK Rate of change of soil moisture
!                           !      due to water falling onto the
!                           !      surface after surface runoff
!                           !      (kg/m2/s).
     &,W_FLUX(NPNTS,0:NSHYD)                                            &
                            ! WORK Fluxes of water between layers
!                           !      (kg/m2/s).
     &,KSZ(NPNTS,0:NSHYD)                                               &
                            ! WORK Saturated hydraulic
!                           !      conductivity in layer (kg/m2/s).
     &,TOP_CRIT(NPNTS)                                                  &
                            ! WORK Critical TI when ZW <=0.0
     &,ZDEPTH(0:NSHYD)                                                  &
                            ! WORK Lower soil layer boundary depth (m).
     &,TSOIL_D(NPNTS)                                                   &
                            ! WORK Soil temperature in the top metre
     &,WUTOT(NPNTS)         ! WORK Ratio of unfrozen to total soil
!                                             !    moisture at ZW.

! C_TOPOG start
! 5.5 17/02/03    Required for large-scale hydrology L_TOP code.
!
! Topographic index increment:
      REAL,PARAMETER :: DTI = 4.0
! Maximum topographic index considered:
      REAL,PARAMETER :: TI_MAX = 10.0
! Maximum allowed water table depth (m):
      REAL,PARAMETER :: ZW_MAX = 5.5
! Standard deviation of LOG(Ksat(0)):
      REAL,PARAMETER :: SIGMA_LOGK = 0.0
! Parameter to remove very high water tables
! from the calculated wetland fraction:
      REAL,PARAMETER :: TI_WETL = 1.5
!
! C_TOPOG end
! C_DENSTY for subroutine SF_EXCH
      REAL,PARAMETER:: RHOSEA = 1026.0 ! density of sea water (kg/m3)
      REAL,PARAMETER:: RHO_WATER = 1000.0! density of pure water (kg/m3)
! C_DENSTY end
! Function & Subroutine calls:
      EXTERNAL                                                          &
     & SFSNOW,SURF_HYD,SOIL_HYD,SOIL_HTC,ICE_HTC,SOILMC                 &
     &,CALC_BASEFLOW,SOILT,CH4_WETL
      EXTERNAL                                                          &
     & TIMER

! End of header--------------------------------------------------------

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('HYDROL ',103)
      ENDIF

!sxy
      if( l_cable ) then
! DEPENDS ON: cable_hyd_driver.o
         call cable_hyd_driver( &
            SNOW_TILE, LYING_SNOW, SURF_ROFF, SUB_SURF_ROFF, &
            TOT_TFALL, WB_LAKE )
      endif

!----------------------------------------------------------------------
! Set up variables required for LSH scheme:
!----------------------------------------------------------------------
      ZDEPTH(0)=0.0
      DO N=1,NSHYD
         ZDEPTH(N)=ZDEPTH(N-1)+DZSOIL(N)
      ENDDO
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Add snowfall to snow mass and update snow grain size
! Calculate interception of snow by coniferous canopy and melt of snow
! beneath canopy if canopy snow model selected
!----------------------------------------------------------------------

      IF ( .NOT. l_cable ) THEN
! DEPENDS ON: sfsnow
        CALL SFSNOW(NPNTS,NTILES,TILE_PTS,TILE_INDEX,CAN_MODEL,         &
     &            CON_SNOW,LS_SNOW,CATCH_SNOW,DZSOIL(1),HCAP,MELT_TILE, &
     &            SMCL,STHF,FRAC,TSTAR_TILE,TIMESTEP,V_SAT,             &
     &            RGRAIN,SNOW_TILE,SNOW_GRND,TSOIL,                     &
     &            LYING_SNOW,SNOW_MELT,SNOMLT_SUB_HTF,                  &
     &            L_SNOW_ALBEDO,LTIMER)
      ENDIF ! .NOT. l_cable
!----------------------------------------------------------------------
! Update the total snowmelt heat flux
!----------------------------------------------------------------------
      IF (STF_HF_SNOW_MELT) THEN
        DO I=1,NPNTS
          HF_SNOW_MELT(I)=LF*SNOW_MELT(I)
          SNOMLT_SUB_HTF(I) = 0.
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
! Calculate throughfall and surface runoff, and update the canopy water
! content
!-----------------------------------------------------------------------
      IF ( .NOT. l_cable ) THEN
! DEPENDS ON: surf_hyd
         CALL SURF_HYD (NPNTS,NTILES,TILE_PTS,TILE_INDEX,               &
     &               CAN_CPY,E_CANOPY,FRAC,INFIL_TILE,CON_RAIN,LS_RAIN, &
     &               MELT_TILE,SNOW_MELT,TIMESTEP,                      &
     &               CAN_WCNT,CAN_WCNT_GB,DSMC_DT,                      &
     &               L_TOP,L_PDM,NSHYD,SOIL_PTS,SOIL_INDEX,             &
     &               SURF_ROFF,TOT_TFALL,                               &
     &               DUN_ROFF,FSAT,V_SAT,STHU,STHF)

      ENDIF ! .NOT. l_cable
!-----------------------------------------------------------------------
! Specify the reduction of hydraulic conductivity with depth:
! Initial base flow to zero:
!-----------------------------------------------------------------------

      DO N=0,NSHYD
!CDIR NODEP
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          KSZ(I,N)=SATCON(I)
        ENDDO
      ENDDO  
      
      DO N=1,NSHYD
!CDIR NODEP
        DO J=1,SOIL_PTS
          QBASE_L(SOIL_INDEX(J),N)=0.0
        END DO
      END DO

      DO I=1,NPNTS
        QBASE(I)=0.0
        QBASE_ZW(I)=0.0
        WUTOT(I)=0.0
        DRAIN(I)=0.0
      ENDDO

      IF ( .NOT. l_cable ) THEN
       IF(L_TOP)THEN
        IF (SOIL_PTS /= 0) THEN
! DEPENDS ON: calc_baseflow
          CALL CALC_BASEFLOW(                                           &
     &        SOIL_PTS,SOIL_INDEX,NPNTS,NSHYD                           &
     &       ,ZDEPTH,KSZ                                                &
     &       ,B,FEXP,TI_MEAN,ZW,STHF,STHU                               &
     &       ,WUTOT,TOP_CRIT,QBASE,QBASE_L                              &
     & )
        ENDIF
       ENDIF
      ENDIF ! .NOT. l_cable


      IF(L_INLAND)THEN

       do i=1,NPNTS

! Add inland basin outflow to change in soil moisture store

         DSMC_DT(I)=DSMC_DT(I)+INLANDOUT_ATM(I)

       Enddo
      ENDIF


!-----------------------------------------------------------------------
! Update the layer soil moisture contents and calculate the
! gravitational drainage.
!-----------------------------------------------------------------------

      IF ( .NOT. l_cable ) THEN
       IF ( SOIL_PTS /= 0 ) THEN
! DEPENDS ON: soil_hyd
       CALL SOIL_HYD (NPNTS,NSHYD,SOIL_PTS,SOIL_INDEX,B,DZSOIL,          &
     &               EXT,DSMC_DT,SATCON,KSZ,SATHH,TIMESTEP,V_SAT,       &
     &               SUB_SURF_ROFF,SMCL,STHU,SURF_ROFF,W_FLUX,          &
     &               STF_SUB_SURF_ROFF,                                 &
     &               ZW,STHZW,ZDEPTH,QBASE,QBASE_L,                     &
     &               DUN_ROFF,DRAIN,L_TOP,L_SOIL_SAT_DOWN,              &
     &               LTIMER)
!-----------------------------------------------------------------------
! Calculate surface saturation and wetland fractions:
!-----------------------------------------------------------------------
       IF(L_TOP)THEN
        DO I=1,NPNTS
           FSAT(I)=0.0
           FWETL(I)=0.0
! Zero soil porosity over land ice:
          IF(V_SAT(I) <= 0.0)ZW(I)=ZW_MAX
        ENDDO
        IF (SOIL_PTS /= 0) THEN
          DO J=1,SOIL_PTS
            I=SOIL_INDEX(J)
            QBASE_ZW(I)=QBASE_L(I,NSHYD+1)
!Now use fit for fsat and fwet:
            FSAT(I)=A_FSAT(I)*EXP(-C_FSAT(I)*TOP_CRIT(I))
            FWETL(I)=A_FWET(I)*EXP(-C_FWET(I)*TOP_CRIT(I))
            IF(TOP_CRIT(I) >= TI_MAX)THEN
              FSAT(I)=0.0
              FWETL(I)=0.0
            ENDIF
          ENDDO
        ENDIF
       ENDIF ! IF(L_TOP)

       ELSE   ! no soil points

!---------------------------------------------------------------------
! If required by STASH flag and there are no soil points,
! set sub-surface runoff to zero.
!---------------------------------------------------------------------

        IF(STF_SUB_SURF_ROFF) THEN
          DO I=1,NPNTS
            SUB_SURF_ROFF(I)=0.0
          ENDDO
        ENDIF

       ENDIF ! IF ( SOIL_PTS /= 0
      ENDIF ! IF NOT l_cable

!-----------------------------------------------------------------------
! Update the soil temperatures and the frozen moisture fractions
!-----------------------------------------------------------------------

      IF (.NOT. l_cable ) THEN
       IF (SOIL_PTS /= 0) THEN
! DEPENDS ON: soil_htc
        CALL SOIL_HTC (NPNTS,NSHYD,NTILES,SOIL_PTS,SOIL_INDEX,          &
     &                 TILE_PTS,TILE_INDEX,B,DZSOIL,FRAC,HCAP,HCON,     &
     &                 SATHH,SNOW_TILE,SURF_HT_FLUX,TIMESTEP,V_SAT,     &
     &                 W_FLUX,SMCL,STHU,STHF,TSOIL,LTIMER)
       ENDIF

!-----------------------------------------------------------------------
! Update the sub-surface temperatures for land ice
!-----------------------------------------------------------------------
       IF (LICE_PTS /= 0) THEN
! DEPENDS ON: ice_htc
        CALL ICE_HTC (NPNTS,NSHYD,LICE_PTS,LICE_INDEX,DZSOIL,           &
     &                SURF_HT_FLUX,TIMESTEP,                            &
     &                TSOIL,LTIMER)
       ENDIF
      ENDIF  ! .NOT. l_cable

!-----------------------------------------------------------------------
! Diagnose the soil moisture in the top metre.
!-----------------------------------------------------------------------
! DEPENDS ON: soilmc
      CALL SOILMC ( NPNTS,NSHYD,SOIL_PTS,SOIL_INDEX,                    &
     &              DZSOIL,STHU,V_SAT,V_WILT,SMC )

!-----------------------------------------------------------------------
! Calculate mean soil temperature and scaled CH4 flux:
!-----------------------------------------------------------------------

      DO I=1,NPNTS
        FCH4_WETL(I)=0.0
      ENDDO
      IF (.NOT. l_cable ) THEN
      IF(L_TOP)THEN
        IF (SOIL_PTS /= 0) THEN
! DEPENDS ON: soilt
          CALL SOILT(NPNTS,NSHYD,SOIL_PTS,SOIL_INDEX                    &
     &             ,DZSOIL,TSOIL,TSOIL_D)
! DEPENDS ON: ch4_wetl
          CALL CH4_WETL(NPNTS,NSHYD,SOIL_PTS,SOIL_INDEX                 &
     &      ,TSOIL_D,CS,FWETL,FCH4_WETL)
        ENDIF
      ENDIF
      ENDIF ! .NOT. l_cable

!     Diagnostic for CABLE
!      print 88,SURF_HT_FLUX(100),TSOIL(100,:),SMCL(100,:), &
!               STHU(100,:),STHF(100,:),SNOW_TILE(100,:)
!88    format(1x,'HYDR100',f5.0,4f4.0,x,2f4.0,2f5.0,x,4f3.2,x,4f3.2,9f6.0)
!      print 89,SURF_HT_FLUX(1000),TSOIL(1000,:),SMCL(1000,:), &
!               STHU(1000,:),STHF(1000,:)
!!               STHU(1000,:),STHF(1000,:),SNOW_TILE(1000,:)
!89    format(1x,'HYD1000',f5.0,4f4.0,x,2f4.0,2f5.0,x,4f3.2,x,4f3.2,9f6.0)
!      print 91,SURF_HT_FLUX(1500),TSOIL(1500,:),SMCL(1500,:), &
!               STHU(1500,:),STHF(1500,:)
!!               STHU(1500,:),STHF(1500,:),SNOW_TILE(1500,:)
!91    format(1x,'HYD1500',f5.0,4f4.0,x,2f4.0,2f5.0,x,4f3.2,x,4f3.2,9f6.0)

!-----------------------------------------------------------------------
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('HYDROL ',104)
      ENDIF

      RETURN
      END SUBROUTINE HYDROL
