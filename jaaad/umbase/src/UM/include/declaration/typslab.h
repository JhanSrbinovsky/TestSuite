! TYPSLAB Declaration of fields passed around within the slab model.
! Requires TYPOCPAR called from TYPSIZE to be included first.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.4   28/06/02   Original list.                     K.Williams
!   5.5   14/02/03   Consistent changes with ARGSLAB
!                    to allow compatibilty with the ITD ice model
!                                                       M. Crucifix
!

      LOGICAL ICY(IMT,JMT)                                              &
                                   ! TRUE if gridbox contains ice
     & ,OPENSEA(IMT,JMT)                                                &
                                   ! TRUE if gridbox is ice-free ocean
     & ,NEWICE(IMT,JMT)                                                 &
                                   ! TRUE if ice is forming this tstep
     & ,LAND(IMT,JMT)              ! TRUE if land. If tiled, this is the
                                   !  only true when gridbox is all land

      REAL FLAND(IMT,JMT)                                               &
                                   ! Fraction of grid box which is land
     & ,TSTAR_SICE(IMT,JMT)                                             &
                                   ! Surface temp at ice points (atmos)
     & ,TSTAR_LAND(IMT,JMT)                                             &
                                   ! Surface temp at land points (atmos)
     & ,COS_LAT_T(IMT,JMT)                                              &
                                   ! COS latitude on theta grid
     & ,COS_LAT_V(IMT,JMTM1)                                            &
                                   ! COS latitude on velocity grid
     & ,SIN_LAT_T(IMT,JMT)                                              &
                                   ! SIN latitude on theta grid
     & ,SIN_LAT_V(IMT,JMTM1)                                            &
                                   ! SIN latitude on velocity grid
     & ,SEC_LAT_T(IMT,JMT)                                              &
                                   ! 1.0/COS latitude on theta grid
     & ,SEC_LAT_V(IMT,JMTM1)                                            &
                                   ! 1.0/COS latitude on velocity grid
     & ,TAN_LAT_V(IMT,JMTM1)                                            &
                                   ! TAN latitude on velocity grid
     & ,CORIOLIS(IMT,JMT)                                               &
                                   ! 2*OMEGA*SIN(lat) on theta grid
     & ,AMX(IMT,JMT)                                                    &
                                   ! Maximum permitted ice conc.
     & ,SNOWRATE(IMT,JMT)                                               &
                                   ! Rate of snowfall (kg m-2 s-1)
     & ,SNOWSLAB(IMT,JMT)                                               &
                                   ! Snowfall rate melting in ocean
     & ,SNOWLEAD(IMT,JMT)                                               &
                                   ! Snowfall rate melting in leads
     & ,ATMSFLUX(IMT,JMT)                                               &
                                   ! Net heat into ocean thro leads
     & ,LEADFLUX(IMT,JMT)                                               &
                                   ! Net heat into ice thro leads
     & ,HSNOWATMN(IMT,JMT,NICE)                                         &
                                    ! Snow depth (m) (av ovr ice-atmos)
     & ,HSNOWSLBN(IMT,JMT,NICE)                                         &
                                    ! Snow depth (m) (ice-slab by catego
     & ,AICEATMN(IMT,JMT,NICE)                                          &
                                    ! Ice concentration (atmos model)
     & ,AICESLBN(IMT,JMT,NICE)                                          &
                                    ! Ice concentration (slab model by c
     & ,HICESLBN(IMT,JMT,NICE)                                          &
                                    ! Ice depth (m) (slb my category)
     & ,HSNOWSLB(IMT,JMT)                                               &
                                    ! Snow depth (m) (aggreted)
     & ,AICESLB(IMT,JMT)                                                &
                                    ! Ice concentration (aggregated)
     & ,HICESLB(IMT,JMT)                                                &
                                    ! Ice depth (m) (aggregated)
     & ,SOLARIN(IMT,JMT)                                                &
                                   ! Net down SW from atmos (all freq)
     & ,BLUEIN(IMT,JMT)                                                 &
                                   ! Net down SW from atmos (band 1)
     & ,EVAP(IMT,JMT)                                                   &
                                   ! Surf evap (weighted by lead area)
     & ,LONGWAVE(IMT,JMT)                                               &
                                   ! Net down LW heat flux
     & ,SENSIBLE(IMT,JMT)                                               &
                                   ! Snsble (+ve up) (wtd by lead area)
     & ,SNOWLS(IMT,JMT)                                                 &
                                   ! Large-scale snow rate (kg m-2 s-1)
     & ,SNOWCONV(IMT,JMT)                                               &
                                   ! Convective snow rate (kg m-2 s-1)
     & ,TCLIMC(IMT,JMT)                                                 &
                                   ! Climatological SSTs (C)
     & ,TSTARATM(IMT,JMT)                                               &
                                   ! SST in atmos model (K)
     & ,SLABTEMP(IMT,JMT)                                               &
                                   ! Temperature of slab ocean (C)
     & ,ADJHCONV(IMT,JMT)                                               &
                                   ! Heat convergence (W m-2)
     & ,SUBLIMA(IMT,JMT)                                                &
                                   ! Accumulated sublimation (kg m-2)
     & ,SUBLIMZ(IMT,JMT)           ! Sublimation rate (kg m-2 s-1)

      REAL                                                              &
     &  TOPMELTZN(IMT,JMT,NICE)                                         &
                                    ! Rate of snow melt (W m-2) by categ
     & ,BOTMELTZN(IMT,JMT,NICE)                                         &
                                    ! Diffus flux thro ice (+ve for melt
     & ,OIFLUX(IMT,JMT)                                                 &
                                   ! Ocean to ice heat flux
     & ,CARYHEAT(IMT,JMT)                                               &
                                   ! Misc heat flux from ice to ocean
     & ,CARYHEAT_O2I(IMT,JMT)                                           &
                                   ! Misc heat flux from ocean to ice
     & ,UCURRENT(IMT,JMTM1)                                             &
                                   ! X component of surf current (m s-1)
     & ,VCURRENT(IMT,JMTM1)                                             &
                                   ! Y component of surf current (m s-1)
     & ,UICE(IMT,JMTM1)                                                 &
                                   ! X component of ice velocity (m s-1)
     & ,VICE(IMT,JMTM1)                                                 &
                                   ! Y component of ice velocity (m s-1)
     & ,WSX(IMT,JMTM1)                                                  &
                                   ! X component of surface wind stress
     & ,WSY(IMT,JMTM1)                                                  &
                                   ! Y component of surface wind stress
     & ,AINC_THERM(IMT,JMT)                                             &
                                   ! ice fraction inc (thermo)
     & ,HINC_THERM(IMT,JMT)                                             &
                                   ! ice depth inc (thermo)
     & ,HSINC_THERM(IMT,JMT)                                            &
                                   ! snow depth inc *ice frac (thermo)
     & ,AINC_DYN(IMT,JMT)                                               &
                                   ! ice fraction inc (dynamics)
     & ,HINC_DYN(IMT,JMT)                                               &
                                   ! ice depth inc (dynamics)
     & ,HSINC_DYN(IMT,JMT)                                              &
                                   ! snow depth inc *ice fract (dyn)
     & ,HINC_DIFF(IMT,JMT)                                              &
                                   ! ice depth inc (diffusion)
     & ,HINC_ADV(IMT,JMT)                                               &
                                   ! ice depth inc (advection)
     & ,HSINC_ADV(IMT,JMT)                                              &
                                   ! snow depth inc *ice fract (adv)
     & ,SIG11NE(IMT,JMT)                                                &
                                   ! Ice internal stresses for EVP
     & ,SIG11SE(IMT,JMT)                                                &
                                   !            ice dynamics (m s-1)
     & ,SIG11SW(IMT,JMT)                                                &
     & ,SIG11NW(IMT,JMT)                                                &
     & ,SIG12NE(IMT,JMT)                                                &
     & ,SIG12SE(IMT,JMT)                                                &
     & ,SIG12SW(IMT,JMT)                                                &
     & ,SIG12NW(IMT,JMT)                                                &
     & ,SIG22NE(IMT,JMT)                                                &
     & ,SIG22SE(IMT,JMT)                                                &
     & ,SIG22SW(IMT,JMT)                                                &
     & ,SIG22NW(IMT,JMT)                                                &
   !kdcorbin, 10/10 - tracer mass
     & ,TMASS(IMT,JMT)                                 &
                                   ! Mass of ice and snow
     & ,PRSS(IMT,JMT)                                                   &
                                   ! Ice Pressure
     & ,DSIG11DX(IMT,JMTM1)                                             &
                                   ! EVP diagnostics
     & ,DSIG12DY(IMT,JMTM1)                                             &
     & ,DSIG21DX(IMT,JMTM1)                                             &
     & ,DSIG22DY(IMT,JMTM1)                                             &
     & ,DDT_AICEN_THERM(IMT,JMT,NICE)                                   &
                                      ! d/dt AICEN due to thermodynamics
     & ,DDT_HICEN_THERM(IMT,JMT,NICE)                                   &
                                      ! d/dt HICEN due to thermodynamics
     & ,DDT_HSNOWN_THERM(IMT,JMT,NICE)                                  &
                                      ! d/dt HSNOWNdue to thermodynamics
     & ,DDT_AICE_THERM(IMT,JMT)                                         &
                                      ! d/dt AICE  due to thermodynamics
     & ,DDT_HICE_THERM(IMT,JMT)                                         &
                                      ! d/dt HICE  due to thermodynamics
     & ,DDT_HSNOW_THERM(IMT,JMT)                                        &
                                      ! d/dt HSNOW due to thermodynamics
     & ,DDT_AICEN_DYN(IMT,JMT,NICE)                                     &
                                      ! d/dt AICEN due to dynamics
     & ,DDT_HICEN_DYN(IMT,JMT,NICE)                                     &
                                      ! d/dt HICEN due to dynamics
     & ,DDT_HSNOWN_DYN(IMT,JMT,NICE)                                    &
                                           ! d/dt HSNOWNdue to dynamics
     & ,DDT_AICE_DYN(IMT,JMT)                                           &
                                      ! d/dt AICE  due to dynamics
     & ,DDT_HICE_DYN(IMT,JMT)                                           &
                                      ! d/dt HICE  due to dynamics
     & ,DDT_HSNOW_DYN(IMT,JMT)        ! d/dt HSNOW due to dynamics
