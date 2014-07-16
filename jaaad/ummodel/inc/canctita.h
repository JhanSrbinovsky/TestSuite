! CANCTITA Store Titles for Atmos/Slab Ancillary Files
!
! History :
! Version   Date   Comment
! -------   ----   -------
!  4.1   20/05/96  New comdeck. New titles for sulphur ancillary
!                  files. Rename SOURCE TERMS to SULPEMIS. D. Robinson
!  4.3    18/3/97  Add aerosol forcing of climate change. William Ingram
!  4.4    12/9/97  Add titles for initial surface type fracs, initial
!                  vegetation state and vegetation disturbance.
!                                                   R.A.Betts
!  4.5   02/09/98  Add title for CO2 emissions ancillary file -
!                  CO2EMITS. Chris Jones.
!  4.5   22/04/98  Add title for soot emissions ancillary file -
!                  SOOTEMIS. R.Rawlins
!  5.3   19/06/01  Add title for Tropopause-Ozone ancillary file
!                  TPPSOZON (tropopause-based ozone) Dave Tan
!  5.5   30/12/02  Add title for Land Fraction, Dust and Biomass
!                  ancillary files. Dave Robinson.
!  6.1   07/04/04  Add title for DMS concentration.  Andy Jones
!  6.1   08/11/04  Add titles for River Routing Files.  R.Sharp
!  7.3   31/05/10  Add titles for tracer fluxes - kdcorbin
!
! -------------------------------------------------------------------
      CHARACTER(LEN=80) :: TITLE_ANCIL_A (NANCIL_DATASETSA)

      DATA TITLE_ANCIL_A /                                              &
     &  'OZONE',                   'SOIL MOISTURE/SNOWDEPTH',           &
     &  'SOIL TEMPERATURES',       'SOIL TYPES',                        &
     &  'VEGETATION TYPE',         'SEA SURFACE TEMPERATURES',          &
     &  'SEA ICE',                 'OCEAN CURRENTS',                    &
     &  'LAND SEA MASK',           'OROGRAPHY',                         &
     &  'HEAT CONVERGENCE (SLAB)', 'Single Level Sulphur Emissions',    &
     &  'AEROSOL TRACERS',         'MURKINESS',                         &
     &  'USER ANCIL.SINGLE',       'USER ANCIL.MULTI',                  &
     &  'Natural SO2 Emissions',   'Chemistry Oxidants',                &
     &  'Sulphate aerosol forcing','Initial fractions of surface types',&
     &  'Initial vegetation state','Disturbed fraction of vegetation',  &
     &  'Soot Emissions'          ,'Surface CO2 emissions'              &
     & ,'Tropopause Ozone'        ,'Land Fraction',                     &
     &  'Dust Soil Properties'    ,'Biomass Emissions'                  &
     & ,'DMS Conc in Seawater'    ,'River Water Storage'                &
     & ,'Riv Channel Dir & Seqnce','River Routing 2A fields'            &
     & ,''                        ,''                                   &
     & ,''                        ,''                                   &
     & ,''                        ,'Clim biogenic aerosol'              &
     &,'Clim biomass burning '    ,'Clim black carbon aero'             &
     &,'Clim sea salt aerosol'    ,'Clim sulphate aerosol'              &
     &,'Clim dust aerosol'        ,'Clim Org C Fossil Fuel'             &
     &,'Clim delta aerosol'       ,'Cariolle ozone tracer'              &
     &,'Fossil-fuel OC aer emiss.','Tracer Flux 1'                      &
     &,'Tracer Flux 2'            ,'Tracer Flux 3'                      &
     &,'Tracer Flux 4'            ,'Tracer Flux 5'                      &
     &,'Tracer Flux 6'            ,'Tracer Flux 7'                      &
     &,'Tracer Flux 8'            ,'Tracer Flux 9'                      &
     &,'Tracer Flux 10'           ,'Tracer Flux 11'                     &
     &,'Tracer Flux 12'           ,'Tracer Flux 13'                     &
     &,'Tracer Flux 14'           ,'Tracer Flux 15'                     &
     &,'Tracer Flux 16'           ,'Tracer Flux 17'                     &
     &,'Tracer Flux 18'           ,'Tracer Flux 19'                     &
     &,'Tracer Flux 20'/
! CANCTITA end
