#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!     UKCA-MODE aerosol code: interface routine called from
!     UKCA_MAIN1 to perform a 1-timestep integration.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Current Code Owner:       Graham Mann
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
! Subroutine Interface:
      SUBROUTINE UKCA_AERO_CTL(i_month, i_day_number,                   &
                      i_hour,i_minute, DTC,                             &
                      p_levelsda, rows, row_length,                     &
                      wet_levels,                                       &
                      global_row_length,global_rows,                    &
                      n_sulf_tracers,                                   &
                      n_mode_tracers,                                   &
                      area,                                             &
                      pres,                                             &
                      temp,                                             &
                      q,                                                &
                      rh3d,                                             &
                      p_layer_boundaries,                               &
                      sulf_tracers,                                     &
                      mode_tracers,                                     &
                      bl_levels,                                        &
                      t_surf,                                           &
                      sea_ice_frac,                                     &
                      z0m,                                              &
                      u_s,                                              &
                      u_10m,                                            &
                      drain, crain,                                     &
                      land_fraction,                                    &
                      NBOX,                                             &
                      delso2_wet_h2o2,                                  &
                      delso2_wet_o3,                                    &
                      mode_diags,                                       &
                      n_emissions,                                      &
                      em_spec,                                          &
                      emissions,                                        &
                      SO2emiss_3D,                                      &
                      cloud_frac)

      USE UKCA_D1_DEFS,     ONLY: N_mode_diags,                         &
                                  ukca_sect, Nukca_D1items,             &
                                  ukcaD1codes, mode_diag_sect,          &
                                  n_chem_emissions
      USE UKCA_CONSTANTS,   ONLY: PPI, AVC, ZBOLTZ, VKARMN, RA, RR,     &
                                  GG, RAD_E, MM_DA, MMSUL, NMOL,        &
                                  EMS_EPS, CONC_EPS, DN_EPS,            &
                                  RHOSUL
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES
!
      IMPLICIT NONE

#include "cmaxsize.h"
#include "nstypes.h"
#include "cruntimc.h"
! molecular weights are in c_v_m in 6.6
#include "c_v_m.h"

! Inputs
      INTEGER, INTENT(IN) :: i_month           ! month
      INTEGER, INTENT(IN) :: i_day_number      ! day
      INTEGER, INTENT(IN) :: i_hour            ! hour
      INTEGER, INTENT(IN) :: i_minute          ! minute
      INTEGER, INTENT(IN) :: p_levelsda        ! # of model levels
      INTEGER, INTENT(IN) :: wet_levels        ! # of wet levels
      INTEGER, INTENT(IN) :: rows              ! # of rows in patch
      INTEGER, INTENT(IN) :: row_length        ! # of pts in a patch row
      INTEGER, INTENT(IN) :: global_rows       ! # of rows (global)
      INTEGER, INTENT(IN) :: global_row_length ! # of pts in a row (global)
      INTEGER, INTENT(IN) :: n_sulf_tracers    ! # of sulphur tracers
      INTEGER, INTENT(IN) :: n_mode_tracers    ! # of mode tracers
      INTEGER, INTENT(IN) :: n_emissions       ! dimension of emission
      INTEGER, INTENT(IN) :: bl_levels         ! # of levels in BL
      INTEGER, INTENT(IN) :: nbox              ! dimension of slice

      REAL, INTENT(IN) :: dtc                              ! timestep(s)
      REAL, INTENT(IN) :: area(row_length,rows,p_levelsda) ! area m^2
      REAL, INTENT(IN) :: pres(row_length,rows,p_levelsda) ! pressure
      REAL, INTENT(IN) :: temp(row_length,rows,p_levelsda) ! temperature
      REAL, INTENT(IN) :: q(row_length,rows,wet_levels)    ! sp humidity
      REAL, INTENT(IN) :: rh3d(row_length,rows,wet_levels) ! RH (frac)
      REAL, INTENT(IN) :: p_layer_boundaries(row_length,                &
                                        rows,0:p_levelsda) ! pressure
      REAL, INTENT(IN)    :: emissions(row_length,rows,                 &
                                  n_chem_emissions)   ! 2D emissions fields
! 3-D volcanic SO2 emissions
      REAL, INTENT(IN) :: SO2emiss_3D(row_length,rows,p_levelsda)
      REAL, INTENT(IN) :: t_surf(row_length, rows)           ! surface temp
      REAL, INTENT(IN) :: sea_ice_frac(row_length, rows)     ! sea ice
      REAL, INTENT(IN) :: u_s(row_length, rows)              ! friction velocity
      REAL, INTENT(IN) :: u_10m(row_length, rows)            ! wind at 10m
      REAL, INTENT(IN) :: z0m(row_length, rows)              ! roughness length
      REAL, INTENT(IN) :: drain(row_length,rows, p_levelsda) ! 3-D LS rain rate
      REAL, INTENT(IN) :: crain(row_length,rows, p_levelsda) ! 3-D conv rain
      REAL, INTENT(IN) :: land_fraction(row_length,rows)     ! land_fraction
! cloud oxidation rates (molecules/cc/DTC) from h2o2 and ozone:
      REAL, INTENT(IN) :: delso2_wet_h2o2(row_length,rows,p_levelsda)
      REAL, INTENT(IN) :: delso2_wet_o3  (row_length,rows,p_levelsda)
! cloud fraction
      REAL, INTENT(IN) :: cloud_frac(row_length, rows, wet_levels)

! sulphur tracers mass mixing ratio
      REAL, INTENT(INOUT) :: sulf_tracers(row_length,rows,              &
                                  p_levelsda,n_sulf_tracers)
! aerosol tracers mass mixing ratio
      REAL, INTENT(INOUT) :: mode_tracers(row_length,rows,              &
                                  p_levelsda,n_mode_tracers)
! 3-D diagnostic array
      REAL, INTENT(INOUT) :: mode_diags(row_length,rows,                &
                                  p_levelsda,n_mode_diags)

! Emission name array
      CHARACTER(len=10), INTENT(IN) :: em_spec(n_emissions)

! Local variables
      INTEGER, PARAMETER :: NMTS=1
! No. of microphysical sub-steps per DTC
      INTEGER :: NZTS
! No. of condensation-nucleation competition sub-steps per DTM
      INTEGER :: PRIMSU_ON
! Switch for whether primary sulfate particle emissions are on/off
      INTEGER :: PRIMBCOC_ON
! Switch for whether primary carbonaceous particle emissions are on/off
      INTEGER, PARAMETER :: PRIMDU_ON=0
! Switch for whether primary dust particle emissions are on/off
      INTEGER :: PRIMSS_ON
! Switch for whether primary sea-salt ptcl emissions are on/off
      INTEGER, PARAMETER :: RAINOUT_ON=1
! Switch for whether rainout (nucl. scav.) is on/off
      INTEGER, PARAMETER :: IMSCAV_ON=1
! Switch for whether impaction scavenging is on/off
      INTEGER, PARAMETER :: WETOX_ON=1
! Switch for whether wet oxidation (cloud processing) is on/off
      INTEGER, PARAMETER :: DDEPAER_ON=1
! Switch for whether aerosol dry deposition is on/off
      INTEGER :: SEDI_ON
! Switch for whether aerosol sedimentation is on/off
      INTEGER, PARAMETER :: COND_ON=1
! Switch for whether vapour condensation is  on/off
      INTEGER :: NUCL_ON
! Switch for whether binary nucleation is on/off
      INTEGER, PARAMETER :: COAG_ON=1
! Switch for whether coagulation is on/off
      INTEGER, PARAMETER :: ICOAG=1
! Switch for KIJ method (1:GLOMAP, 2: M7, 3: UMorig, 4:UMorig MFPP)
!   =3 Cunnigham scheme as in UM, =4 as in UM but computing values)
      INTEGER, PARAMETER :: IMERGE=2
! Switch to use mid-pts (=1), edges (2) or dynamic (=3) in remode
      INTEGER, PARAMETER :: IFUCHS=2
! Switch for Fuchs(1964) (=1) or Fuchs-Sutugin(1971) for CC (=2)
      INTEGER, PARAMETER :: IACTMETHOD=1
! Switch for activation method (0=off,1=fixed ract,2=NSO3 scheme)
      INTEGER, PARAMETER :: IWVOLMETHOD=2
! Switch for wet volume method (1=behave-as-H2SO4,2=multi-cpt)
      INTEGER :: INUCSCAV
! Switch for nucl scav method (1=as GLOMAP Spr05, 2=use scav coeffs)
      INTEGER :: IDDEPAER
! Switch for dry dep method (1=as GLOMAP Spr05, 2=incl. sedi)
      INTEGER, PARAMETER :: VERBOSE=0
! Switch to determine level of debug o/p (0=none, 1, 2)
      INTEGER, PARAMETER :: CHECKMD_ND=0
! Switch for whether to check for bad values of MD and ND
      INTEGER, PARAMETER :: INTRAOFF=0
! Switch to turn off intra-modal coagulation
      INTEGER, PARAMETER :: INTEROFF=0
! Switch to turn off inter-modal coagulation
      INTEGER, PARAMETER :: IDUSTEMS=0
! Switch for using Pringle scheme (=1) or AEROCOMdaily (=2)
      INTEGER, PARAMETER :: DRYOX_IN_AER=0
! Switch : update SU & SO by SO2 & terp dryox in aerosol code? (1/0)
      INTEGER, PARAMETER :: WETOX_IN_AER=0
! Switch : calculate wet oxidation of SO2 in aerosol code? (1/0)

      INTEGER  :: NEMVOLCONSO2(NBOX)
! No. of volcanic SO2 ems sources in gridbox (continuously eruptive)
      INTEGER  :: NEMVOLEXPSO2(NBOX)
! No. of volcanic SO2 ems sources in gridbox (explosively eruptive)
      INTEGER  :: LDAY(NBOX)
! Switch for day/night (1/0)
      INTEGER  :: JLABOVE(NBOX)
! Index of box directly above this grid box
      REAL     :: EMANSO2(NBOX,6)
! Anthrop. SO2 ems rates, low sources (kgSO2/box/s)
      REAL     :: EMVOLCONSO2(NBOX,10)
! Volcanic SO2 ems rates (cont. src) (kgSO2/box/s)
      REAL     :: EMVOLEXPSO2(NBOX,10)
! Volcanic SO2 ems rates (expl. src) (kgSO2/box/s)
      REAL     :: EMBIOMSO2(NBOX,6)
! Biomass SO2 ems rates (kgSO2/box/s) [6 alt ranges]
      REAL     :: EMC(NBOX,4)
! BC/OC emission rates from bio- & fossil-fuels (kgC/box/s)
      REAL     :: EMCBM(NBOX,2,6)
! BC/OC emission rates from biomass burning (kgC/box/s)
      REAL     :: ND(NBOX,NMODES)
! Aerosol ptcl number density for mode (cm^-3)
      REAL     :: MDT(NBOX,NMODES)
! Avg tot mass of aerosol ptcl in mode (particle^-1)
      REAL     :: MDWAT(NBOX,NMODES)
! Molecular concentration of water (molecules per particle)
      REAL     :: DRYDP(NBOX,NMODES)
! Geometric mean dry diameter of particles in each mode (m)
      REAL     :: WETDP(NBOX,NMODES)
! Geometric mean wet diameter of particles in each mode (m)
      REAL     :: RHOPAR(NBOX,NMODES)
! Total particle density [incl. H2O & insoluble cpts] (kgm^-3)
      REAL     :: DVOL(NBOX,NMODES)
! Geometric mean dry volume of particles in each mode (m^3)
      REAL     :: WVOL(NBOX,NMODES)
! Geometric mean wet volume of particles in each mode (m^3)
      REAL     :: MD(NBOX,NMODES,NCP)
! Avg cpt mass of aerosol particle in mode (particle^-1)
      REAL     :: S0(NBOX,NADVG)
! Partial masses of gas phase species (kg per gridbox)
      REAL     :: S0_DOT_CONDENSABLE(NBOX,NCHEMG)
! ASAD tendencies for condensable gas phase species (vmr per s)
      REAL     :: SM(NBOX)     ! Mass of air in gridbox (kg)
! Grid box mass of air (kg)
      REAL     :: AIRD(NBOX)
! Number density of air (per cm3)
      REAL     :: AIRDM3(NBOX)
! Number density of air (per m3)
      REAL     :: RHOA(NBOX)
! Air density (kg/m3)
      REAL     :: VBA(NBOX)
! Mean free speed of air molecules (m/s)
      REAL     :: TSQRT(NBOX)
! Square-root of centre level temperature (K)
      REAL     :: DVISC(NBOX)
! Dynamic viscosity of air (kg m^-1 s^-1)
      REAL     :: MFPA(NBOX)
! Mean free path of air (m)
      REAL     :: T(NBOX)
! Air temperature at mid-point (K)
      REAL     :: RH(NBOX)
! Relative humidity (fraction)
      REAL     :: S(NBOX)
! Specific humidity (kg/kg)
      REAL     :: PMID(NBOX)
! Air pressure at mid-point (Pa)
      REAL     :: PUPPER(NBOX)
! Air pressure at upper interface (Pa)
      REAL     :: PLOWER(NBOX)
! Air pressure at lower interface (Pa)
      REAL     :: ZO3(NBOX)
! Background vmr of O3 (dimensionless)
      REAL     :: ZHO2(NBOX)
! Background conc. of HO2 (molecules per cc)
      REAL     :: ZH2O2(NBOX)
! Background conc. of H2O2 (molecules per cc)
      REAL     :: USTR(NBOX)
! Surface friction velocity (m/s)
      REAL     :: US10M(NBOX)
! Scalar wind at 10m (ms-1)
      REAL     :: ZNOT(NBOX)
! Roughness length (m)
      REAL     :: SURTP(NBOX)
! Surface type: 0=seasurf,1=landsurf,2=above-seasurf,3=above-landsurf
      REAL     :: SURF(NBOX)
! Surface area of box (horizontal) (m^2)
      REAL     :: LAND_FRAC(NBOX)
! Fraction of horizontal gridbox area covered by land
      REAL     :: SEAICE(NBOX)
! Fraction of horizontal gridbox area containing seaice
      REAL     :: CRAING(NBOX)
! Rain rate for conv precip. in box (kgm^-2s^-1)
      REAL     :: DRAING(NBOX)
! Rain rate for dyn. precip. in box (kgm^-2s^-1)
      REAL     :: CRAING_UP(NBOX)
! Rain rate for conv precip. in box above (kgm^-2s^-1)
      REAL     :: DRAING_UP(NBOX)
! Rain rate for dyn. precip. in box above (kgm^-2s^-1)
      REAL     :: FCONV_CONV(NBOX)
! Fraction of box condensate --> rain in 6 hours (conv)
      REAL     :: LOWCLOUD(NBOX)
! Horizontal low cloud fraction
      REAL     :: VFAC(NBOX)
! Vertical low cloud fraction
      REAL     :: LWC(NBOX)
! Cloud liquid water content [kg/m3]
      REAL     :: DTM
! Microphysics time step (s)
      REAL     :: DTZ
! Competition (cond/nucl) time step (s)
      REAL     :: RMOIS
! Month of year
      REAL     :: DLON
! Delta longitude
      REAL     :: DLAT
! Delta latitude

! Outputs for budget calculations
      REAL :: BUD_AER_MAS(NBOX,NBUDAER)

! Outputs for budget calculations
      REAL :: MASPRIMSUAITSOL(NBOX),MASPRIMSUACCSOL(NBOX)
      REAL :: MASPRIMSUCORSOL(NBOX),MASPRIMSSACCSOL(NBOX)
      REAL :: MASPRIMSSCORSOL(NBOX),MASPRIMBCAITINS(NBOX)
      REAL :: MASPRIMOCAITSOL(NBOX),MASPRIMOCAITINS(NBOX)
      REAL :: MASDDEPSUNUCSOL(NBOX),MASDDEPSUAITSOL(NBOX)
      REAL :: MASDDEPSUACCSOL(NBOX),MASDDEPSUCORSOL(NBOX)
      REAL :: MASDDEPSSACCSOL(NBOX),MASDDEPSSCORSOL(NBOX)
      REAL :: MASDDEPBCAITSOL(NBOX),MASDDEPBCACCSOL(NBOX)
      REAL :: MASDDEPBCCORSOL(NBOX),MASDDEPBCAITINS(NBOX)
      REAL :: MASDDEPOCAITSOL(NBOX),MASDDEPOCACCSOL(NBOX)
      REAL :: MASDDEPOCCORSOL(NBOX),MASDDEPOCAITINS(NBOX)
      REAL :: MASDDEPSONUCSOL(NBOX),MASDDEPSOAITSOL(NBOX)
      REAL :: MASDDEPSOACCSOL(NBOX),MASDDEPSOCORSOL(NBOX)
      REAL :: MASNUSCSUNUCSOL(NBOX),MASNUSCSUAITSOL(NBOX)
      REAL :: MASNUSCSUACCSOL(NBOX),MASNUSCSUCORSOL(NBOX)
      REAL :: MASNUSCSSACCSOL(NBOX),MASNUSCSSCORSOL(NBOX)
      REAL :: MASNUSCBCAITSOL(NBOX),MASNUSCBCACCSOL(NBOX)
      REAL :: MASNUSCBCCORSOL(NBOX),MASNUSCBCAITINS(NBOX)
      REAL :: MASNUSCOCAITSOL(NBOX),MASNUSCOCACCSOL(NBOX)
      REAL :: MASNUSCOCCORSOL(NBOX),MASNUSCOCAITINS(NBOX)
      REAL :: MASNUSCSONUCSOL(NBOX),MASNUSCSOAITSOL(NBOX)
      REAL :: MASNUSCSOACCSOL(NBOX),MASNUSCSOCORSOL(NBOX)
      REAL :: MASIMSCSUNUCSOL(NBOX),MASIMSCSUAITSOL(NBOX)
      REAL :: MASIMSCSUACCSOL(NBOX),MASIMSCSUCORSOL(NBOX)
      REAL :: MASIMSCSSACCSOL(NBOX),MASIMSCSSCORSOL(NBOX)
      REAL :: MASIMSCBCAITSOL(NBOX),MASIMSCBCACCSOL(NBOX)
      REAL :: MASIMSCBCCORSOL(NBOX),MASIMSCBCAITINS(NBOX)
      REAL :: MASIMSCOCAITSOL(NBOX),MASIMSCOCACCSOL(NBOX)
      REAL :: MASIMSCOCCORSOL(NBOX),MASIMSCOCAITINS(NBOX)
      REAL :: MASIMSCSONUCSOL(NBOX),MASIMSCSOAITSOL(NBOX)
      REAL :: MASIMSCSOACCSOL(NBOX),MASIMSCSOCORSOL(NBOX)
      REAL :: MASCLPRSUAITSOL1(NBOX),MASCLPRSUACCSOL1(NBOX)
      REAL :: MASCLPRSUCORSOL1(NBOX)
      REAL :: MASCLPRSUAITSOL2(NBOX),MASCLPRSUACCSOL2(NBOX)
      REAL :: MASCLPRSUCORSOL2(NBOX)
      REAL :: MASCONDSUNUCSOL(NBOX),MASCONDSUAITSOL(NBOX)
      REAL :: MASCONDSUACCSOL(NBOX),MASCONDSUCORSOL(NBOX)
      REAL :: MASCONDSUAITINS(NBOX)
      REAL :: MASCONDSONUCSOL(NBOX),MASCONDSOAITSOL(NBOX)
      REAL :: MASCONDSOACCSOL(NBOX),MASCONDSOCORSOL(NBOX)
      REAL :: MASCONDSOAITINS(NBOX)
      REAL :: MASNUCLSUNUCSOL(NBOX)
      REAL :: MASCOAGSUINTR12(NBOX),MASCOAGSUINTR13(NBOX)
      REAL :: MASCOAGSUINTR14(NBOX),MASCOAGSUINTR15(NBOX)
      REAL :: MASCOAGSOINTR12(NBOX),MASCOAGSOINTR13(NBOX)
      REAL :: MASCOAGSOINTR14(NBOX),MASCOAGSOINTR15(NBOX)
      REAL :: MASCOAGSUINTR23(NBOX),MASCOAGBCINTR23(NBOX)
      REAL :: MASCOAGOCINTR23(NBOX),MASCOAGSOINTR23(NBOX)
      REAL :: MASCOAGSUINTR24(NBOX),MASCOAGBCINTR24(NBOX)
      REAL :: MASCOAGOCINTR24(NBOX),MASCOAGSOINTR24(NBOX)
      REAL :: MASCOAGSUINTR34(NBOX),MASCOAGBCINTR34(NBOX)
      REAL :: MASCOAGOCINTR34(NBOX),MASCOAGSSINTR34(NBOX)
      REAL :: MASCOAGSOINTR34(NBOX)
      REAL :: MASCOAGBCINTR53(NBOX),MASCOAGOCINTR53(NBOX)
      REAL :: MASCOAGBCINTR54(NBOX),MASCOAGOCINTR54(NBOX)
      REAL :: MASAGEDSUINTR52(NBOX),MASAGEDBCINTR52(NBOX)
      REAL :: MASAGEDOCINTR52(NBOX),MASAGEDSOINTR52(NBOX)
      REAL :: MASMERGSUINTR12(NBOX),MASMERGSOINTR12(NBOX)
      REAL :: MASMERGSUINTR23(NBOX),MASMERGBCINTR23(NBOX)
      REAL :: MASMERGOCINTR23(NBOX),MASMERGSOINTR23(NBOX)
      REAL :: MASMERGSUINTR34(NBOX),MASMERGSSINTR34(NBOX)
      REAL :: MASMERGBCINTR34(NBOX),MASMERGOCINTR34(NBOX)
      REAL :: MASMERGSOINTR34(NBOX)
      REAL :: MASPROCSUINTR23(NBOX),MASPROCBCINTR23(NBOX)
      REAL :: MASPROCOCINTR23(NBOX),MASPROCSOINTR23(NBOX)

      INTEGER :: I,J,K,L,N,JL,IMODE,ICP
      INTEGER :: JLAT
      INTEGER :: n_reqd_tracers           ! No of tracers required
      INTEGER :: lso2emlo,lso2emhi
      INTEGER :: ITRA,ITRA2
      INTEGER :: II_ND(NMODES)
      INTEGER :: II_MD(NMODES,NCP)

      INTEGER :: ISO2EMS ! setting for emissions (now from namelist)

! ISO2EMS == 1 : 3.0% of SO2 ems --> primary SO4 ems --> 15%/85% to
!                10/70 nm g.m.diam. modes as in Spracklen (2005)
!                and Binkowski & Shankar (1995).
!
! ISO2EMS == 2 : 2.5% of SO2 ems --> primary SO4 ems
!                road/off-road/domestic      all -->   30nm gm.diam mode
!                industrial/power-plant/ship all --> 1000nm gm.diam mode
!                as for original AEROCOM size recommendations.
!
! ISO2EMS >= 3 : 2.5% of SO2 ems --> primary SO4 ems --> 50%/50% to
!                150/1500nm  g.m.diam. modes as for Stier et al (2005)
!                modified AEROCOM sizdis recommendations.
!
! Note also that currently, ISO2EMS also controls size assumptions for
! primary carbonaceous aerosol (this needs to be changed):
!
! ISO2EMS /= 3 : biofuel & biomass BC/OC emissions -->  80nm g.m.diam.
!                fossil-fuel BC/OC emissions --> 30nm g.m.diam.
!
! ISO2EMS == 3 : biofuel & biomass BC/OC emissions --> 150nm g.m.diam.
!                fossil-fuel BC/OC emissions --> 60nm g.m.diam.
!
      REAL    :: PARFRAC
! PARFRAC is fraction of mass of SO2 emitted as primary SO4
      REAL    :: DELSO2(NBOX)
!  S(IV) --> S(VI) by H2O2 (molecules per cc) [input if WETOX_IN_AER=0]
      REAL    :: DELSO2_2(NBOX)
!  S(IV) --> S(VI) by O3   (molecules per cc) [input if WETOX_IN_AER=0]

      REAL, PARAMETER :: MA=4.78E-26 ! mass of air molecule (kg)
      REAL :: ACT                    ! radius for activation (m)

      CHARACTER(LEN=36)             :: sulf_tracer_names(3)=(/          &
       'DMS MMR                             ',                          &
       'SO2 MMR                             ',                          &
       'H2SO4 MMR                           '/)

      CHARACTER(LEN=36)            :: mode_tracer_names(25)=(/          &
       'NUCLEATION MODE (SOLUBLE) NUMBER    ',                          &
       'NUCLEATION MODE (SOLUBLE) H2SO4 MMR ',                          &
       'AITKEN MODE (SOLUBLE) NUMBER        ',                          &
       'AITKEN MODE (SOLUBLE) H2SO4 MMR     ',                          &
       'AITKEN MODE (SOL) BLACK CARBON MMR  ',                          &
       'AITKEN MODE (SOL) ORGANIC CARBON MMR',                          &
       'ACCUMULATION MODE (SOLUBLE) NUMBER  ',                          &
       'ACCUMULATION MODE (SOL) H2SO4 MMR   ',                          &
       'ACCUM MODE (SOL) BLACK CARBON MMR   ',                          &
       'ACCUM MODE (SOL) ORGANIC CARBON MMR ',                          &
       'ACCUMULATION MODE (SOL) SEA SALT MMR',                          &
       'ACCUMULATION MODE (SOLUBLE) DUST MMR',                          &
       'COARSE MODE (SOLUBLE) NUMBER        ',                          &
       'COARSE MODE (SOLUBLE) SO4 MMR       ',                          &
       'COARSE MODE (SOLUBLE) BLACK CARBON  ',                          &
       'COARSE MODE (SOLUBLE) ORGANIC CARBON',                          &
       'COARSE MODE (SOLUBLE) SEA SALT MMR  ',                          &
       'COARSE MODE (SOLUBLE) DUST MMR      ',                          &
       'AITKEN MODE (INSOLUBLE) NUMBER      ',                          &
       'AITKEN MODE (INSOLUBLE) BLACK CARBON',                          &
       'AITKEN MODE (INSOL) ORGANIC CARBON  ',                          &
       'ACCUMULATION MODE (INSOLUBLE) NUMBER',                          &
       'ACCUMULATION MODE (INSOL) DUST MMR  ',                          &
       'COARSE MODE (INSOLUBLE) NUMBER      ',                          &
       'COARSE MODE (INSOLUBLE) DUST MMR    '/)

! used for debug output
      LOGICAL            :: mode_tracer_debug(25)=(/                    &
       .true., .true., .true., .true., .false.,                         &
       .false., .true., .true., .false., .false.,                       &
       .true., .false., .true., .true., .false.,                        &
       .false., .true., .false., .false., .false.,                      &
       .false., .false., .false., .false., .false./)

      CHARACTER(LEN=72) :: cmessage     ! Error message

! Controls output, to be replaced by info from umui eventually
      LOGICAL, PARAMETER :: L_ukca_mode_diags=.false.

      LOGICAL :: TOOHI ! switch to check if dry radius out of bounds
      LOGICAL :: TOOLO ! switch to check if dry radius out of bounds
      INTEGER :: N_DPFIX_1D(NBOX,NMODES) ! counter: fixes to drydp
      INTEGER :: N_MERGE_1D(NBOX,NMODES) ! counter: mode-merges applied
      INTEGER :: N_DPFIX_3D(row_length,rows,p_levelsda,NMODES)
      INTEGER :: N_MERGE_3D(row_length,rows,p_levelsda,NMODES)
      INTEGER :: IDX(NBOX)
      INTEGER :: NFIX

! Set local variables to values in run_ukca namelist
      NZTS = i_mode_nzts
      IDDEPAER = i_mode_ddepaer
      INUCSCAV = i_mode_nucscav
      ISO2EMS = i_mode_sizeprim
      ACT = mode_actdryr
      IF (L_UKCA_PRIMSU) THEN
        PRIMSU_ON = 1
      ELSE
        PRIMSU_ON = 0
      ENDIF
      IF (L_UKCA_PRIMSS) THEN
        PRIMSS_ON = 1
      ELSE
        PRIMSS_ON = 0
      ENDIF
      IF (L_UKCA_PRIMBCOC) THEN
        PRIMBCOC_ON = 1
      ELSE
        PRIMBCOC_ON = 0
      ENDIF
      IF (L_UKCA_SEDI) THEN
        SEDI_ON = 1
      ELSE
        SEDI_ON = 0
      ENDIF
      IF (L_UKCA_NUCL) THEN
        NUCL_ON = 1
      ELSE
        NUCL_ON = 0
      ENDIF

      DTM=DTC/REAL(NMTS)
      DTZ=DTM/REAL(NZTS)
      RMOIS=REAL(i_month)
      DLAT=PPI/REAL(global_rows)
      DLON=2.0*PPI/REAL(global_row_length)

      IF(VERBOSE.GE.1) THEN

      write(6,*) 'UKCA_MODE INPUTS: '
      write(6,*) 'i_month: ',i_month
      write(6,*) 'i_day_number: ',i_day_number
      write(6,*) 'i_hour: ',i_hour
      write(6,*) 'i_minute: ',i_minute
      write(6,*) 'DTC: ',DTC
      write(6,*) 'p_levelsda: ',p_levelsda
      write(6,*) 'rows: ',rows
      write(6,*) 'row_length: ',row_length
      write(6,*) 'global_rows: ',global_rows
      write(6,*) 'global_row_length: ',global_row_length
      write(6,*) 'n_sulf_tracers: ',n_sulf_tracers
      write(6,*) 'n_mode_tracers: ',n_mode_tracers

      write(6,*) 'Array:     MIN        MAX         MEAN'
      write(6,*) 'area: ',minval(area),maxval(area),                    &
                 sum(area)/real(size(area))
      i=0
      write(6,*) 'Level: ',i
      write(6,*) 'p_layer_boundaries: ',                                &
                 minval(p_layer_boundaries(:,:,i)),                     &
                 sum(p_layer_boundaries(:,:,i))/                        &
                 real(size(p_layer_boundaries(:,:,i)))
      do i=1,2            ! p_levelsda
      write(6,*) 'Level: ',i
      write(6,*) 'pres: ',minval(pres(:,:,i)),maxval(pres(:,:,i)),      &
                 sum(pres(:,:,i))/real(size(pres(:,:,i)))
      write(6,*) 'temp: ',minval(temp(:,:,i)),maxval(temp(:,:,i)),      &
                 sum(temp(:,:,i))/real(size(temp(:,:,i)))
      write(6,*) 'q: ',minval(q(:,:,i)),maxval(q(:,:,i)),               &
                 sum(q(:,:,i))/real(size(q(:,:,i)))
      write(6,*) 'rh3d: ',minval(rh3d(:,:,i)),maxval(rh3d(:,:,i)),      &
                 sum(rh3d(:,:,i))/real(size(rh3d(:,:,i)))
      write(6,*) 'p_layer_boundaries: ',                                &
                 minval(p_layer_boundaries(:,:,i)),                     &
                 sum(p_layer_boundaries(:,:,i))/                        &
                 real(size(p_layer_boundaries(:,:,i)))
      write(6,*) 'delso2_wet_h2o2: ',                                   &
                 minval(delso2_wet_h2o2(:,:,i)),                        &
                 sum(delso2_wet_h2o2(:,:,i))/                           &
                 real(size(delso2_wet_h2o2(:,:,i)))
      enddo

      do i=1,2       ! p_levelsda
      do j=1,n_sulf_tracers
      write(6,*) 'Level: ',i,' Tracer: ',j,sulf_tracer_names(j)
      write(6,*) 'sulf_tracers: ',minval(sulf_tracers(:,:,i,j)),        &
                 maxval(sulf_tracers(:,:,i,j)),                         &
                 sum(sulf_tracers(:,:,i,j))/                            &
                 real(size(sulf_tracers(:,:,i,j)))
      enddo
      do j=1,n_mode_tracers
        if (mode_tracer_debug(j)) then
        write(6,*) 'Level: ',i,' Tracer: ',j,mode_tracer_names(j)
        write(6,*) 'mode_tracers: ',minval(mode_tracers(:,:,i,j)),      &
                   maxval(mode_tracers(:,:,i,j)),                       &
                   sum(mode_tracers(:,:,i,j))/                          &
                   real(size(mode_tracers(:,:,i,j)))
        endif
      enddo
      enddo     ! p_levelsda

      ENDIF ! IF(VERBOSE.GE.1)

! Calculate number of aerosol tracers required
      n_reqd_tracers = 0 
      DO i=1,nmodes 
        DO j=1,ncp 
          IF (component(i,j)) n_reqd_tracers = n_reqd_tracers + 1 
        ENDDO 
      ENDDO 
      n_reqd_tracers = n_reqd_tracers + sum(mode_choice) 

! Check the number of tracers, warn if too many, stop if too few
      IF (n_mode_tracers > n_reqd_tracers) THEN
       cmessage=' Too many tracers input to UKCA_MODE'
       write(6,*) cmessage,n_mode_tracers,n_reqd_tracers
! DEPENDS ON: ereport
       CALL EREPORT('UKCA_AERO_CTL',-1,cmessage)
      ENDIF
      IF (n_mode_tracers < n_reqd_tracers) THEN
       cmessage=' Warning: not all advected aerosol tracers being used'
       write(6,*) cmessage,n_mode_tracers,n_reqd_tracers
! DEPENDS ON: ereport
       CALL EREPORT('UKCA_AERO_CTL',1,cmessage)
      ENDIF

      N_DPFIX_3D=0
      N_MERGE_3D=0

      DO l=1,n_chem_emissions
       IF (em_spec(l) == 'SO2_low ' ) lso2emlo=l
       IF (em_spec(l) == 'SO2_high' ) lso2emhi=l
      ENDDO

      IF(PRIMSU_ON == 1) THEN
       IF(ISO2EMS == 0) PARFRAC=0.0
       IF(ISO2EMS == 1) PARFRAC=0.03
       IF(ISO2EMS > 1) PARFRAC=0.025
      ELSE
       PARFRAC=0.0
      ENDIF

!! .. below is the main template which matches the
!! .. 25 prognostic tracers as indexed in mode_tracers to
!! .. the ND and MD values (as set using variables ITRA,ITRA2)
      II_ND                =(/ 1, 3, 7,13,19,22,24/)
      II_MD(1:NMODES,CP_SU)=(/ 2, 4, 8,14, 0, 0, 0/)
      II_MD(1:NMODES,CP_BC)=(/ 0, 5, 9,15,20, 0, 0/)
      II_MD(1:NMODES,CP_OC)=(/ 0, 6,10,16,21, 0, 0/)
      II_MD(1:NMODES,CP_CL)=(/ 0, 0,11,17, 0, 0, 0/)
      II_MD(1:NMODES,CP_DU)=(/ 0, 0,12,18, 0,23,25/)
      II_MD(1:NMODES,CP_SO)=(/ 0, 0, 0, 0, 0, 0, 0/) ! not used yet

      DO l=1,p_levelsda

       IF(VERBOSE.GE.2) THEN
        write(6,*) 'In AERO_CTL before AERO_STEP call, level ',l
       ENDIF
       DO i=1,rows
        DO k=1,row_length
         JL=(i-1)*row_length+k
         T(JL)=temp(k,i,l)
        IF(l <= wet_levels) THEN
          RH(JL)=rh3d(k,i,l)
          S(JL)=q(k,i,l)
        ELSE
         RH(JL)=0.0
         S(JL)=0.0
        ENDIF
         PMID(JL)=pres(k,i,l)
         PUPPER(JL)=p_layer_boundaries(k,i,l)
         PLOWER(JL)=p_layer_boundaries(k,i,l-1)
         USTR(JL)=u_s(k,i) ! still set to surface value even above surf.
         US10M(JL)=u_10m(k,i) ! still set to surface value.....
         ZNOT(JL)=z0m(k,i) ! still set to surface value even above surf.
        SEAICE(JL)=sea_ice_frac(k,i) ! fraction of sea_ice at surface
         LAND_FRAC(JL)=land_fraction(k,i) ! fraction of land at surface
         IF(land_frac(JL) < 1.0) THEN    ! some sea
          IF(l.EQ.1) SURTP(JL)=0.0 ! 0 if at seasurf
          IF(l.NE.1) SURTP(JL)=2.0 ! 2 if above sea
         ELSE
          IF(l.EQ.1) SURTP(JL)=1.0 ! 1 if at landsurf
          IF(l.NE.1) SURTP(JL)=3.0 ! 3 if above land
         ENDIF
         SURF(JL)=area(k,i,l)
         CRAING(JL)=crain(k,i,l) ! convective rain
         DRAING(JL)=drain(k,i,l) ! dynamic rain
         IF (l == p_levelsda) THEN
           CRAING_UP(JL)=0.0
           DRAING_UP(JL)=0.0
         ELSE
           CRAING_UP(JL)=crain(k,i,l+1) ! convective rain in box above
           DRAING_UP(JL)=drain(k,i,l+1) ! dynamic rain in box above
         ENDIF
         FCONV_CONV(JL)=0.99 ! (fraction of condensate-->rain in 6 hrs)
        IF(l <= wet_levels) THEN
          LOWCLOUD(JL)=cloud_frac(k,i,l)
         ELSE
         LOWCLOUD(JL)=0.0
        ENDIF
         VFAC(JL)=1.0 ! set to 1 so that VFAC*LOWCLOUD=cloud_frac
         DELSO2(JL)=delso2_wet_h2o2(k,i,l)
         DELSO2_2(JL)=delso2_wet_o3(k,i,l)
         AIRD(JL)=PMID(JL)/(T(JL)*ZBOLTZ*1.0E6) ! no conc of air (/cm3)
         AIRDM3(JL)=AIRD(JL)*1.0e6              ! no conc of air (/m3)
         RHOA(JL)=PMID(JL)/(T(JL)*RA)
         VBA(JL)=SQRT(8.0*ZBOLTZ*T(JL)/(PPI*MA))
         TSQRT(JL)=SQRT(T(JL))
         DVISC(JL)=1.83E-5*(416.16/(T(JL)+120.0))*                      &
                   (SQRT(T(JL)/296.16)**3)
         MFPA(JL)=2.0*DVISC(JL)/(RHOA(JL)*VBA(JL))
         SM(JL)=SURF(JL)*(PLOWER(JL)-PUPPER(JL))/GG ! mass air in box
! .. set H2SO4 entry in S0 from H2SO4 mmr in gas phase tracer array
         S0(JL,:)=0.0
! .. sulf_tracers(:,:,:,3) is H2SO4 mmr in kgH2SO4/kgdryair
         IF(sulf_tracers(k,i,l,3) > 0.0) THEN
          S0(JL,MH2SO4)=(MM_DA/MMSUL)*sulf_tracers(k,i,l,3)*SM(JL)
        ENDIF
         S0_DOT_CONDENSABLE(JL,:)=0.0
! Set this to zero since do condensable update in CHEM_CTL (DRYOX_AER=0)
         EMANSO2(JL,:)=0.0

         IF(PRIMSU_ON == 1) THEN      ! do primary S emissions
          IF(L == 1) THEN ! if gridbox is at surface

! .. multiplying by SURF below converts from kgSO2/m2/s to kgSO2/box/s
           EMANSO2(JL,1)=emissions(k,i,lso2emlo)*SURF(JL) ! kgSO2/box/s

!!! ** NOTE HAVE COMMENTED THIS OUT AS CAN PRODUCE NEGATIVES **
!!! ** SINCE PBL MIXING HAD BEEN CARRIED OUT AFTER EMISSIONS **
!!! ** Allow extra % of SO2 emissions --> to sulfate for now **
!!
!!!! .. below reduces SO2 mmr according to amount emitted as primary SO4
!!!! .. emissions are in mmrS whereas sulf_traces are in mmrSO2, etc.
!!!! .. sulf_tracers(:,:,:,2) is SO2 mmr in kgSO2/kgdryair
!!!           sulf_tracers(k,i,l,2)=sulf_tracers(k,i,l,2)-                 &
!!!     &   PARFRAC*(m_so2/m_s)*emissions(k,i,lso2emlo)*DTC ! kgSO2/kgair

!!
!! order of sulf_tracers is i) DMS, ii) SO2, iii) H2SO4
!!
          ENDIF

          IF(L == SO2_HIGH_LEVEL) THEN

! .. Add high-level primary SO4 emissions into level specified in UMUI
! .. multiplying by SURF below converts from kgSO2/m2/s to kgSO2/box/s
           EMANSO2(JL,2)=emissions(k,i,lso2emhi)*SURF(JL) ! kgSO2/box/s

!!! ** NOTE HAVE COMMENTED THIS OUT AS CAN PRODUCE NEGATIVES **
!!! ** SINCE PBL MIXING HAD BEEN CARRIED OUT AFTER EMISSIONS **
!!! ** Allow extra % of SO2 emissions --> to sulfate for now **
!!
!!!! .. below reduces SO2 mmr according to amount emitted as primary SO4
!!!! .. emissions are in mmrS whereas sulf_traces are in mmrSO2, etc.
!!!! .. sulf_tracers(:,:,:,2) is SO2 mmr in kgSO2/kgdryair
!!!           sulf_tracers(k,i,l,2)=sulf_tracers(k,i,l,2)-                 &
!!!     &      PARFRAC*(m_so2/m_s)*emissions(k,i,lso2emhi)*DTC

          ENDIF
         ENDIF ! if PRIMSU_ON = 1

         EMBIOMSO2(JL,:)=0.0
         EMVOLCONSO2(JL,:)=0.0
         NEMVOLCONSO2(JL)=0

         IF(PRIMSU_ON == 1) THEN

! .. Add volcanic SO4 emissions,
! ..  multiplying by SURF below converts from kgSO2/m2/s to kgSO2/box/s
           EMVOLCONSO2(JL,1)=SO2emiss_3D(k,i,l)*SURF(JL) ! kgSO2/box/s

!!! ** NOTE HAVE COMMENTED THIS OUT AS CAN PRODUCE NEGATIVES **
!!! ** SINCE PBL MIXING HAD BEEN CARRIED OUT AFTER EMISSIONS **
!!! ** Allow extra % of SO2 emissions --> to sulfate for now **
!!
!!!! .. below reduces SO2 mmr according to amount emitted as primary SO4
!!!! .. emissions are in mmrS whereas sulf_traces are in mmrSO2, etc.
!!!! .. sulf_tracers(:,:,:,2) is SO2 mmr in kgSO2/kgdryair
!!!          sulf_tracers(k,i,l,2)=sulf_tracers(k,i,l,2)-                  &
!!!     &     PARFRAC*(m_so2/m_s)*SO2emiss_3D(k,i,l)*DTC ! kgSO2/kgair

          EMVOLCONSO2(JL,2:10)=0.0
          NEMVOLCONSO2(JL)=1
         ENDIF ! if PRIMSU_ON = 1

         EMVOLEXPSO2(JL,:)=0.0
         NEMVOLEXPSO2(JL)=0
         EMC(JL,:)=0.0
         EMCBM(JL,:,:)=0.0
         ZO3(JL)=0.0     ! currently do wet ox separately in UM
         ZHO2(JL)=0.0    ! currently do wet ox separately in UM
         ZH2O2(JL)=0.0   ! currently do wet ox separately in UM
         LWC(JL)=0.0     ! currently do wet ox separately in UM.
         LDAY(JL)=0      ! currently do wet ox separately in UM
         JLABOVE(JL)=0   ! for now dont do sedimentation (IDDEPAER=1)
        ENDDO ! loop over longitudes (i=1,row_length)
       ENDDO ! loop over latitudes (k=1,rows)

       DO IMODE=1,NMODES
        IF(MODE(IMODE)) THEN
         ITRA=II_ND(IMODE)
         DO i=1,rows
         JLAT=i
          DO k=1,row_length
           JL=(i-1)*row_length+k
! .. set ND (particles per cc) from aerosol tracer array
           ND(JL,IMODE)=mode_tracers(k,i,l,ITRA)*AIRD(JL)
           MDT(JL,IMODE)=0.0
          ENDDO ! loop over longitudes (i=1,row_length)
         ENDDO ! loop over latitudes (k=1,rows)
         DO ICP=1,NCP
          IF(COMPONENT(IMODE,ICP)) THEN
           ITRA=II_MD(IMODE,ICP)
          IF(ITRA == 0) THEN
           write(6,*) '***** ERROR: MD_CP specified by COMPONENT'
           write(6,*) '*****        doesnt exist in II_MD'
           write(6,*) 'IMODE,ICP,II_MD=',IMODE,ICP,II_MD(IMODE,ICP)
           write(6,*) 'COMPONENT(IMODE,ICP)=',COMPONENT(IMODE,ICP)
           write(6,*) 'II_MD(IMODE,ICP)=',II_MD(IMODE,ICP)
           cmessage='MD_CP specified by COMPONENT not in II_MD'
! DEPENDS ON: ereport
            CALL EREPORT('UKCA_AERO_CTL',1,cmessage)
          ENDIF
           DO i=1,rows
            DO k=1,row_length
             JL=(i-1)*row_length+k
             IF (ND(JL,IMODE) > NUM_EPS(IMODE)) THEN
! .. set MD (molecules per particle) from aerosol tracer array
              MD(JL,IMODE,ICP)=(MM_DA/MM(ICP))*AIRD(JL)*                &
                 mode_tracers(k,i,l,ITRA)/ND(JL,IMODE)
             ELSE
              MD(JL,IMODE,ICP)=MMID(IMODE)*MFRAC_0(IMODE,ICP)
             ENDIF
             MDT(JL,IMODE)=MDT(JL,IMODE)+MD(JL,IMODE,ICP)
            ENDDO ! loop over longitudes (i=1,row_length)
           ENDDO ! loop over latitudes (k=1,rows)
          ELSE
           MD(:,IMODE,ICP)=0.0
          ENDIF
         ENDDO ! loop over cpts
        ELSE
         ND(:,IMODE)=0.0
         MDT(:,IMODE)=MMID(IMODE)
         DO ICP=1,NCP
          IF(COMPONENT(IMODE,ICP)) THEN
           MD(:,IMODE,ICP)=MMID(IMODE)*MFRAC_0(IMODE,ICP)
          ELSE
           MD(:,IMODE,ICP)=0.0
          ENDIF
         ENDDO
        ENDIF
       ENDDO ! loop over modes

       IF(VERBOSE.GE.2) THEN
        write(6,*) 'n_sulf_tracers=',n_sulf_tracers
        do itra=1,n_sulf_tracers
         write(6,*) itra,sulf_tracer_names(itra),                       &
           minval(sulf_tracers(:,:,l,itra)),                            &
           maxval(sulf_tracers(:,:,l,itra))
       enddo
        do itra=1,nadvg
         write(6,*) 'S0 : ',itra,minval(s0(:,itra)),maxval(s0(:,itra))
        enddo
        do imode=1,nmodes
         if(mode(imode)) then
          write(6,*) 'ND : ',IMODE,minval(ND (:,IMODE)),                &
                                   maxval(ND (:,IMODE))
          write(6,*) 'MDT: ',IMODE,minval(MDT(:,IMODE)),                &
                                   maxval(MDT(:,IMODE))
          do icp=1,ncp
           if(component(imode,icp)) then
            write(6,*) 'MD : ',IMODE,minval(MD(:,IMODE,ICP))            &
                                    ,maxval(MD(:,IMODE,ICP))
           endif
          enddo
         endif
        enddo

       ENDIF

! DEPENDS ON: ukca_calc_drydiam
       CALL UKCA_CALC_DRYDIAM(NBOX,ND,MD,MDT,DRYDP,DVOL,VERBOSE)

       IF(VERBOSE.GE.2) THEN
        do imode=1,nmodes
         if(mode(imode)) then
          write(6,*) 'DRYDP',IMODE,minval(DRYDP(:,IMODE)),              &
                                   maxval(DRYDP(:,IMODE))
         endif
        enddo
       ENDIF

       N_DPFIX_1D=0
       N_MERGE_1D=0

       DO IMODE=1,NMODES
        IF(MODE(IMODE)) THEN
         NFIX = 0
         DO JL=1,NBOX
          TOOHI=DRYDP(JL,IMODE) > DDPLIM1(IMODE)
          TOOLO=DRYDP(JL,IMODE) < DDPLIM0(IMODE)
          IF(TOOHI.OR.TOOLO) THEN
           NFIX=NFIX+1
           IDX(NFIX)=JL
          ENDIF
         ENDDO
         DO I=1,NFIX
          JL = IDX(I)
          IF(VERBOSE >= 2) THEN
            write(6,*) '*********************************'
            write(6,*) 'MD,ND imply DRPDP out of bounds'//              &
        ', re-setting MD,MDT-->MMID for JL,IMODE=',JL,IMODE
            write(6,*) 'Before re-set DRYDP=',DRYDP(JL,IMODE)
          ENDIF
          MDT(JL,IMODE)=MMID(IMODE)
          DVOL(JL,IMODE)=MMID(IMODE)*MMSUL/(AVC*RHOSUL)
          DO ICP=1,NCP
           IF(COMPONENT(IMODE,ICP)) THEN
            MD(JL,IMODE,ICP)=MMID(IMODE)*MFRAC_0(IMODE,ICP)
           ENDIF
          ENDDO
          DRYDP(JL,IMODE)=                                              &
         ((6.0*DVOL(JL,IMODE))/(PPI*X(IMODE)))**(1.0/3.0)
          IF(VERBOSE.GE.2) THEN
           write(6,*) 'After reset DRYDP=',DRYDP(JL,IMODE)
           write(6,*) '*********************************'
          ENDIF
          N_DPFIX_1D(JL,IMODE)=N_DPFIX_1D(JL,IMODE)+1
! .. increment counter to track number of DRYDP fixes
         ENDDO
        ENDIF
       ENDDO

       IF(VERBOSE.GE.2) THEN
        DO imode=1,nmodes
         IF(mode(imode)) THEN
          write(6,*) 'DRYDP',IMODE,minval(DRYDP(:,IMODE)),              &
                                   maxval(DRYDP(:,IMODE))
         ENDIF
        ENDDO

       ENDIF

!! zero aerosol budget terms before calling UKCA_AERO_STEP
       BUD_AER_MAS(:,:)=0.0

!! below is call to UKCA_AERO_STEP as at ukca_mode_v12_e17_gm20_new3.f (Jul

! DEPENDS ON: ukca_aero_step
       CALL UKCA_AERO_STEP(NBOX,                                        &
       ND,MDT,MD,MDWAT,S0,DRYDP,WETDP,RHOPAR,DVOL,WVOL,SM,              &
       AIRD,AIRDM3,RHOA,MFPA,DVISC,T,TSQRT,RH,S,PMID,PUPPER,PLOWER,     &
       EMC,EMCBM,ZO3,ZHO2,ZH2O2,USTR,US10M,ZNOT,                        &
       SURTP,LAND_FRAC,SURF,SEAICE,                                     &
       CRAING,DRAING,CRAING_UP,DRAING_UP,FCONV_CONV,LOWCLOUD,VFAC,      &
       EMANSO2,EMVOLCONSO2,NEMVOLCONSO2,EMVOLEXPSO2,NEMVOLEXPSO2,       &
       EMBIOMSO2,ISO2EMS,                                               &
       DTC,DTM,DTZ,NMTS,NZTS,LDAY,ACT,BUD_AER_MAS,                      &
       PRIMSU_ON,PRIMBCOC_ON,PRIMSS_ON,PRIMDU_ON,RAINOUT_ON,            &
       IMSCAV_ON,WETOX_ON,DDEPAER_ON,SEDI_ON,                           &
       DRYOX_IN_AER,WETOX_IN_AER,DELSO2,DELSO2_2,                       &
       COND_ON,NUCL_ON,COAG_ON,ICOAG,IMERGE,IFUCHS,IWVOLMETHOD,         &
       IACTMETHOD,IDDEPAER,INUCSCAV,VERBOSE,CHECKMD_ND,JLAT,INTRAOFF,   &
       INTEROFF,IDUSTEMS,S0_DOT_CONDENSABLE,LWC,JLABOVE,N_MERGE_1D)

! Update tracers
       DO i=1,rows
        DO k=1,row_length
         JL=(i-1)*row_length+k
! .. update gas phase H2SO4 mmr following H2SO4 condensation/nucleation
         sulf_tracers(k,i,l,3)=S0(JL,MH2SO4)*(MMSUL/MM_DA)/SM(JL) ! if mmrH
        ENDDO ! loop over longitudes (k=1,row_length)
       ENDDO ! loop over latitudes (i=1,nrows)
       DO IMODE=1,NMODES
        IF(MODE(IMODE)) THEN
         ITRA2=II_ND(IMODE)
         DO i=1,rows
          DO k=1,row_length
           JL=(i-1)*row_length+k
! .. update aerosol no. conc. following aerosol microphysics
           mode_tracers(k,i,l,ITRA2)=ND(JL,IMODE)/AIRD(JL)
          ENDDO ! loop over longitudes (k=1,row_length)
         ENDDO ! loop over latitudes (i=1,nrows)
         DO ICP=1,NCP
          IF(COMPONENT(IMODE,ICP)) THEN
           ITRA2=II_MD(IMODE,ICP)
           DO i=1,rows
            DO k=1,row_length
             JL=(i-1)*row_length+k
             IF (ND(JL,IMODE).LE.NUM_EPS(IMODE)) THEN
              MD(JL,IMODE,ICP)=MMID(IMODE)*MFRAC_0(IMODE,ICP)
             ENDIF
! .. update aerosol mmr following aerosol microphysics
             mode_tracers(k,i,l,ITRA2)=(MM(ICP)/MM_DA)*                 &
                       MD(JL,IMODE,ICP)*ND(JL,IMODE)/AIRD(JL)
            ENDDO ! loop over longitudes (k=1,row_length)
           ENDDO ! loop over latitudes (i=1,nrows)
          ENDIF
         ENDDO ! loop over cpts
        ENDIF
        DO i=1,rows
         DO k=1,row_length
          JL=(i-1)*row_length+k
          N_DPFIX_3D(k,i,l,IMODE)=N_DPFIX_1D(JL,IMODE)
          N_MERGE_3D(k,i,l,IMODE)=N_MERGE_1D(JL,IMODE)
         ENDDO ! loop over longitudes (k=1,row_length)
        ENDDO ! loop over latitudes (i=1,nrows)
       ENDDO ! loop over modes

!!!     Calculate rbar_i, S_i, V_i and set to mode_diags diagnostics array
!!!     so 7*3=21 diagnostics (leave out for now)
!!!
!!!     At some stage also pass back deposition arrays for input to earth s
!!!
!!!     Eventually also pass back budget variables.

        IF(VERBOSE.GE.2) THEN

        write(6,*) 'AFTER CALL TO UKCA_AERO_STEP_SUSS'
        do itra=1,nadvg
         write(6,*) 'S0 : ',itra,minval(s0(:,itra)),maxval(s0(:,itra))
        enddo
        do imode=1,nmodes
         if(mode(imode)) then
          write(6,*) 'ND : ',IMODE,minval(ND (:,IMODE)),                &
                                   maxval(ND (:,IMODE))
          write(6,*) 'MDT: ',IMODE,minval(MDT(:,IMODE)),                &
                                   maxval(MDT(:,IMODE))
          do icp=1,ncp
           if(component(imode,icp)) then
            write(6,*) 'MD : ',IMODE,minval(MD(:,IMODE,ICP))            &
                                    ,maxval(MD(:,IMODE,ICP))
           endif
          enddo
         endif
        enddo

! DEPENDS ON: ukca_calc_drydiam
        CALL UKCA_CALC_DRYDIAM(NBOX,ND,MD,MDT,DRYDP,DVOL,VERBOSE)

        do imode=1,nmodes
         if(mode(imode)) then
          write(6,*) 'DRYDP',IMODE,minval(DRYDP(:,IMODE)),              &
                                   maxval(DRYDP(:,IMODE))
         endif
        enddo

        ENDIF ! if VERBOSE.GE.2

!! For now set MODE diagnostics from BUD_AER_MAS
!! to those currently in MODE diagnostics list:

        MASPRIMSUAITSOL(:)=BUD_AER_MAS(:,NMASPRIMSUAITSOL)
        MASPRIMSUACCSOL(:)=BUD_AER_MAS(:,NMASPRIMSUACCSOL)
        MASPRIMSUCORSOL(:)=BUD_AER_MAS(:,NMASPRIMSUCORSOL)
        MASPRIMSSACCSOL(:)=BUD_AER_MAS(:,NMASPRIMSSACCSOL)
        MASPRIMSSCORSOL(:)=BUD_AER_MAS(:,NMASPRIMSSCORSOL)
       MASDDEPSUNUCSOL(:)=BUD_AER_MAS(:,NMASDDEPSUNUCSOL)
       MASDDEPSUAITSOL(:)=BUD_AER_MAS(:,NMASDDEPSUAITSOL)
       MASDDEPSUACCSOL(:)=BUD_AER_MAS(:,NMASDDEPSUACCSOL)
       MASDDEPSUCORSOL(:)=BUD_AER_MAS(:,NMASDDEPSUCORSOL)
       MASDDEPSSACCSOL(:)=BUD_AER_MAS(:,NMASDDEPSSACCSOL)
       MASDDEPSSCORSOL(:)=BUD_AER_MAS(:,NMASDDEPSSCORSOL)
       MASNUSCSUNUCSOL(:)=BUD_AER_MAS(:,NMASNUSCSUNUCSOL)
       MASNUSCSUAITSOL(:)=BUD_AER_MAS(:,NMASNUSCSUAITSOL)
       MASNUSCSUACCSOL(:)=BUD_AER_MAS(:,NMASNUSCSUACCSOL)
       MASNUSCSUCORSOL(:)=BUD_AER_MAS(:,NMASNUSCSUCORSOL)
       MASNUSCSSACCSOL(:)=BUD_AER_MAS(:,NMASNUSCSSACCSOL)
       MASNUSCSSCORSOL(:)=BUD_AER_MAS(:,NMASNUSCSSCORSOL)
       MASIMSCSUNUCSOL(:)=BUD_AER_MAS(:,NMASIMSCSUNUCSOL)
       MASIMSCSUAITSOL(:)=BUD_AER_MAS(:,NMASIMSCSUAITSOL)
       MASIMSCSUACCSOL(:)=BUD_AER_MAS(:,NMASIMSCSUACCSOL)
       MASIMSCSUCORSOL(:)=BUD_AER_MAS(:,NMASIMSCSUCORSOL)
       MASIMSCSSACCSOL(:)=BUD_AER_MAS(:,NMASIMSCSSACCSOL)
       MASIMSCSSCORSOL(:)=BUD_AER_MAS(:,NMASIMSCSSCORSOL)
       MASCLPRSUAITSOL1(:)=BUD_AER_MAS(:,NMASCLPRSUAITSOL1)
       MASCLPRSUACCSOL1(:)=BUD_AER_MAS(:,NMASCLPRSUACCSOL1)
       MASCLPRSUCORSOL1(:)=BUD_AER_MAS(:,NMASCLPRSUCORSOL1)
       MASCLPRSUAITSOL2(:)=BUD_AER_MAS(:,NMASCLPRSUAITSOL2)
       MASCLPRSUACCSOL2(:)=BUD_AER_MAS(:,NMASCLPRSUACCSOL2)
       MASCLPRSUCORSOL2(:)=BUD_AER_MAS(:,NMASCLPRSUCORSOL2)
       MASCONDSUNUCSOL(:)=BUD_AER_MAS(:,NMASCONDSUNUCSOL)
       MASCONDSUAITSOL(:)=BUD_AER_MAS(:,NMASCONDSUAITSOL)
       MASCONDSUACCSOL(:)=BUD_AER_MAS(:,NMASCONDSUACCSOL)
       MASCONDSUCORSOL(:)=BUD_AER_MAS(:,NMASCONDSUCORSOL)
       MASNUCLSUNUCSOL(:)=BUD_AER_MAS(:,NMASNUCLSUNUCSOL)
       MASCOAGSUINTR12(:)=BUD_AER_MAS(:,NMASCOAGSUINTR12)
       MASCOAGSUINTR13(:)=BUD_AER_MAS(:,NMASCOAGSUINTR13)
       MASCOAGSUINTR14(:)=BUD_AER_MAS(:,NMASCOAGSUINTR14)
       MASCOAGSUINTR23(:)=BUD_AER_MAS(:,NMASCOAGSUINTR23)
       MASCOAGSUINTR24(:)=BUD_AER_MAS(:,NMASCOAGSUINTR24)
       MASCOAGSUINTR34(:)=BUD_AER_MAS(:,NMASCOAGSUINTR34)
       MASCOAGSSINTR34(:)=BUD_AER_MAS(:,NMASCOAGSSINTR34)
       MASMERGSUINTR12(:)=BUD_AER_MAS(:,NMASMERGSUINTR12)
       MASMERGSUINTR23(:)=BUD_AER_MAS(:,NMASMERGSUINTR23)
       MASMERGSUINTR34(:)=BUD_AER_MAS(:,NMASMERGSUINTR34)
       MASMERGSSINTR34(:)=BUD_AER_MAS(:,NMASMERGSSINTR34)
       MASPROCSUINTR23(:)=BUD_AER_MAS(:,NMASPROCSUINTR23)

        IF (VERBOSE > 0 .AND. l == 1) THEN
       write(6,*) 'MODE Fluxes (SUM, MAX, MIN, MEAN) for Level: ',l
       write(6,'(A18,4E14.4)') 'MPRIMSUAITSOL: ',SUM(MASPRIMSUAITSOL),  &
          MAXVAL(MASPRIMSUAITSOL), MINVAL(MASPRIMSUAITSOL),             &
          SUM(MASPRIMSUAITSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MPRIMSUACCSOL: ',SUM(MASPRIMSUACCSOL),  &
          MAXVAL(MASPRIMSUACCSOL), MINVAL(MASPRIMSUACCSOL),             &
          SUM(MASPRIMSUACCSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MPRIMSUCORSOL: ',SUM(MASPRIMSUCORSOL),  &
          MAXVAL(MASPRIMSUCORSOL), MINVAL(MASPRIMSUCORSOL),             &
          SUM(MASPRIMSUCORSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MPRIMSSACCSOL: ',SUM(MASPRIMSSACCSOL),  &
          MAXVAL(MASPRIMSSACCSOL), MINVAL(MASPRIMSSACCSOL),             &
          SUM(MASPRIMSSACCSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MPRIMSSCORSOL: ',SUM(MASPRIMSSCORSOL),  &
          MAXVAL(MASPRIMSSCORSOL), MINVAL(MASPRIMSSCORSOL),             &
          SUM(MASPRIMSSCORSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MDDEPSUNUCSOL: ',SUM(MASDDEPSUNUCSOL),  &
          MAXVAL(MASDDEPSUNUCSOL), MINVAL(MASDDEPSUNUCSOL),             &
          SUM(MASDDEPSUNUCSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MDDEPSUAITSOL: ',SUM(MASDDEPSUAITSOL),  &
          MAXVAL(MASDDEPSUAITSOL), MINVAL(MASDDEPSUAITSOL),             &
          SUM(MASDDEPSUAITSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MDDEPSUACCSOL: ',SUM(MASDDEPSUACCSOL),  &
          MAXVAL(MASDDEPSUACCSOL), MINVAL(MASDDEPSUACCSOL),             &
          SUM(MASDDEPSUACCSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MDDEPSUCORSOL: ',SUM(MASDDEPSUCORSOL),  &
          MAXVAL(MASDDEPSUCORSOL), MINVAL(MASDDEPSUCORSOL),             &
          SUM(MASDDEPSUCORSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MDDEPSSACCSOL: ',SUM(MASDDEPSSACCSOL),  &
          MAXVAL(MASDDEPSSACCSOL), MINVAL(MASDDEPSSACCSOL),             &
          SUM(MASDDEPSSACCSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MDDEPSSCORSOL: ',SUM(MASDDEPSSCORSOL),  &
          MAXVAL(MASDDEPSSCORSOL), MINVAL(MASDDEPSSCORSOL),             &
          SUM(MASDDEPSSCORSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MNUSCSUNUCSOL: ',SUM(MASNUSCSUNUCSOL),  &
          MAXVAL(MASNUSCSUNUCSOL), MINVAL(MASNUSCSUNUCSOL),             &
          SUM(MASNUSCSUNUCSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MNUSCSUAITSOL: ',SUM(MASNUSCSUAITSOL),  &
          MAXVAL(MASNUSCSUAITSOL), MINVAL(MASNUSCSUAITSOL),             &
          SUM(MASNUSCSUAITSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MNUSCSUACCSOL: ',SUM(MASNUSCSUACCSOL),  &
          MAXVAL(MASNUSCSUACCSOL), MINVAL(MASNUSCSUACCSOL),             &
          SUM(MASNUSCSUACCSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MNUSCSUCORSOL: ',SUM(MASNUSCSUCORSOL),  &
          MAXVAL(MASNUSCSUCORSOL), MINVAL(MASNUSCSUCORSOL),             &
          SUM(MASNUSCSUCORSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MNUSCSSACCSOL: ',SUM(MASNUSCSSACCSOL),  &
          MAXVAL(MASNUSCSSACCSOL), MINVAL(MASNUSCSSACCSOL),             &
          SUM(MASNUSCSSACCSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MNUSCSSCORSOL: ',SUM(MASNUSCSSCORSOL),  &
          MAXVAL(MASNUSCSSCORSOL), MINVAL(MASNUSCSSCORSOL),             &
          SUM(MASNUSCSSCORSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MIMSCSUNUCSOL: ',SUM(MASIMSCSUNUCSOL),  &
          MAXVAL(MASIMSCSUNUCSOL), MINVAL(MASIMSCSUNUCSOL),             &
          SUM(MASIMSCSUNUCSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MIMSCSUAITSOL: ',SUM(MASIMSCSUAITSOL),  &
          MAXVAL(MASIMSCSUAITSOL), MINVAL(MASIMSCSUAITSOL),             &
          SUM(MASIMSCSUAITSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MIMSCSUACCSOL: ',SUM(MASIMSCSUACCSOL),  &
          MAXVAL(MASIMSCSUACCSOL), MINVAL(MASIMSCSUACCSOL),             &
          SUM(MASIMSCSUACCSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MIMSCSUCORSOL: ',SUM(MASIMSCSUCORSOL),  &
          MAXVAL(MASIMSCSUCORSOL), MINVAL(MASIMSCSUCORSOL),             &
          SUM(MASIMSCSUCORSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MIMSCSSACCSOL: ',SUM(MASIMSCSSACCSOL),  &
          MAXVAL(MASIMSCSSACCSOL), MINVAL(MASIMSCSSACCSOL),             &
          SUM(MASIMSCSSACCSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MIMSCSSCORSOL: ',SUM(MASIMSCSSCORSOL),  &
          MAXVAL(MASIMSCSSCORSOL), MINVAL(MASIMSCSSCORSOL),             &
          SUM(MASIMSCSSCORSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MCLPRSUAITSOL1: ',                      &
          SUM(MASCLPRSUAITSOL1),                                        &
          MAXVAL(MASCLPRSUAITSOL1), MINVAL(MASCLPRSUAITSOL1),           &
          SUM(MASCLPRSUAITSOL1)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MCLPRSUAITSOL2: ',                      &
          SUM(MASCLPRSUAITSOL2),                                        &
          MAXVAL(MASCLPRSUAITSOL2), MINVAL(MASCLPRSUAITSOL2),           &
          SUM(MASCLPRSUAITSOL2)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MCLPRSUACCSOL1: ',                      &
          SUM(MASCLPRSUACCSOL1),                                        &
          MAXVAL(MASCLPRSUACCSOL1), MINVAL(MASCLPRSUACCSOL1),           &
          SUM(MASCLPRSUACCSOL1)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MCLPRSUACCSOL2: ',                      &
          SUM(MASCLPRSUACCSOL2),                                        &
          MAXVAL(MASCLPRSUACCSOL2), MINVAL(MASCLPRSUACCSOL2),           &
          SUM(MASCLPRSUACCSOL2)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MCLPRSUCORSOL1: ',                      &
          SUM(MASCLPRSUCORSOL1),                                        &
          MAXVAL(MASCLPRSUCORSOL1), MINVAL(MASCLPRSUCORSOL1),           &
          SUM(MASCLPRSUCORSOL1)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MCLPRSUCORSOL2: ',                      &
          SUM(MASCLPRSUCORSOL2),                                        &
          MAXVAL(MASCLPRSUCORSOL2), MINVAL(MASCLPRSUCORSOL2),           &
          SUM(MASCLPRSUCORSOL2)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MCONDSUNUCSOL: ',SUM(MASCONDSUNUCSOL),  &
          MAXVAL(MASCONDSUNUCSOL), MINVAL(MASCONDSUNUCSOL),             &
          SUM(MASCONDSUNUCSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MCONDSUAITSOL: ',SUM(MASCONDSUAITSOL),  &
          MAXVAL(MASCONDSUAITSOL), MINVAL(MASCONDSUAITSOL),             &
          SUM(MASCONDSUAITSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MCONDSUACCSOL: ',SUM(MASCONDSUACCSOL),  &
          MAXVAL(MASCONDSUACCSOL), MINVAL(MASCONDSUACCSOL),             &
          SUM(MASCONDSUACCSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MCONDSUCORSOL: ',SUM(MASCONDSUCORSOL),  &
          MAXVAL(MASCONDSUCORSOL), MINVAL(MASCONDSUCORSOL),             &
          SUM(MASCONDSUCORSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MNUCLSUNUCSOL: ',SUM(MASNUCLSUNUCSOL),  &
          MAXVAL(MASNUCLSUNUCSOL), MINVAL(MASNUCLSUNUCSOL),             &
          SUM(MASNUCLSUNUCSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MNUSCSUACCSOL: ',SUM(MASNUSCSUACCSOL),  &
          MAXVAL(MASNUSCSUACCSOL), MINVAL(MASNUSCSUACCSOL),             &
          SUM(MASNUSCSUACCSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MNUSCSSACCSOL: ',SUM(MASNUSCSSACCSOL),  &
          MAXVAL(MASNUSCSSACCSOL), MINVAL(MASNUSCSSACCSOL),             &
          SUM(MASNUSCSSACCSOL)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MCOAGSUINTR12: ',SUM(MASCOAGSUINTR12),  &
          MAXVAL(MASCOAGSUINTR12), MINVAL(MASCOAGSUINTR12),             &
          SUM(MASCOAGSUINTR12)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MCOAGSUINTR13: ',SUM(MASCOAGSUINTR13),  &
          MAXVAL(MASCOAGSUINTR13), MINVAL(MASCOAGSUINTR13),             &
          SUM(MASCOAGSUINTR13)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MCOAGSUINTR14: ',SUM(MASCOAGSUINTR14),  &
          MAXVAL(MASCOAGSUINTR14), MINVAL(MASCOAGSUINTR14),             &
          SUM(MASCOAGSUINTR14)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MCOAGSUINTR23: ',SUM(MASCOAGSUINTR23),  &
          MAXVAL(MASCOAGSUINTR23), MINVAL(MASCOAGSUINTR23),             &
          SUM(MASCOAGSUINTR23)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MCOAGSUINTR24: ',SUM(MASCOAGSUINTR24),  &
          MAXVAL(MASCOAGSUINTR24), MINVAL(MASCOAGSUINTR24),             &
          SUM(MASCOAGSUINTR24)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MCOAGSUINTR34: ',SUM(MASCOAGSUINTR34),  &
          MAXVAL(MASCOAGSUINTR34), MINVAL(MASCOAGSUINTR34),             &
          SUM(MASCOAGSUINTR34)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MCOAGSSINTR34: ',SUM(MASCOAGSSINTR34),  &
          MAXVAL(MASCOAGSSINTR34), MINVAL(MASCOAGSSINTR34),             &
          SUM(MASCOAGSSINTR34)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MMERGSUINTR12: ',SUM(MASMERGSUINTR12),  &
          MAXVAL(MASMERGSUINTR12), MINVAL(MASMERGSUINTR12),             &
          SUM(MASMERGSUINTR12)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MMERGSUINTR23: ',SUM(MASMERGSUINTR23),  &
          MAXVAL(MASMERGSUINTR23), MINVAL(MASMERGSUINTR23),             &
          SUM(MASMERGSUINTR23)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MMERGSUINTR34: ',SUM(MASMERGSUINTR34),  &
          MAXVAL(MASMERGSUINTR34), MINVAL(MASMERGSUINTR34),             &
          SUM(MASMERGSUINTR34)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MMERGSSINTR34: ',SUM(MASMERGSSINTR34),  &
          MAXVAL(MASMERGSSINTR34), MINVAL(MASMERGSSINTR34),             &
          SUM(MASMERGSSINTR34)/REAL(NBOX)
       write(6,'(A18,4E14.4)') 'MPROCSUINTR23: ',SUM(MASPROCSUINTR23),  &
          MAXVAL(MASPROCSUINTR23), MINVAL(MASPROCSUINTR23),             &
          SUM(MASPROCSUINTR23)/REAL(NBOX)
        ENDIF   ! verbose

! Write 3_D diagnostics to mode_diags array
!
        IF (L_UKCA_MODE_diags) THEN  ! fill 3D array
          k=0
          DO N=1,nukca_D1items
            IF (ukcaD1codes(N)%section == MODE_diag_sect .AND.          &
                ukcaD1codes(N)%item >= 301 .AND.                        &
                ukcaD1codes(N)%item <= 405 .AND.                        &
                ukcaD1codes(N)%required) THEN
                k=k+1

!! BUD_AER_MAS now used rather than MASPRIMSUAITSOL etc.
!!
                SELECT CASE(UkcaD1codes(N)%item)
                  CASE(301)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASPRIMSUAITSOL,(/row_length,rows/))
                  CASE(302)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASPRIMSUACCSOL,(/row_length,rows/))
                  CASE(303)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASPRIMSUCORSOL,(/row_length,rows/))
                  CASE(304)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASPRIMSSACCSOL,(/row_length,rows/))
                  CASE(305)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASPRIMSSCORSOL,(/row_length,rows/))
                  CASE(306)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASDDEPSUNUCSOL,(/row_length,rows/))
                  CASE(307)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASDDEPSUAITSOL,(/row_length,rows/))
                  CASE(308)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASDDEPSUACCSOL,(/row_length,rows/))
                  CASE(309)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASDDEPSUCORSOL,(/row_length,rows/))
                  CASE(310)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASDDEPSSACCSOL,(/row_length,rows/))
                  CASE(311)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASDDEPSSCORSOL,(/row_length,rows/))
                  CASE(324)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASNUSCSUNUCSOL,(/row_length,rows/))
                  CASE(325)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASNUSCSUAITSOL,(/row_length,rows/))
                  CASE(326)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASNUSCSUACCSOL,(/row_length,rows/))
                  CASE(327)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASNUSCSUCORSOL,(/row_length,rows/))
                  CASE(328)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASNUSCSSACCSOL,(/row_length,rows/))
                  CASE(329)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASNUSCSSCORSOL,(/row_length,rows/))
                  CASE(342)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASCLPRSUAITSOL1,(/row_length,rows/))
                  CASE(343)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASCLPRSUACCSOL1,(/row_length,rows/))
                  CASE(344)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASCLPRSUCORSOL1,(/row_length,rows/))
                  CASE(345)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASCLPRSUAITSOL2,(/row_length,rows/))
                  CASE(346)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASCLPRSUACCSOL2,(/row_length,rows/))
                  CASE(347)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASCLPRSUCORSOL2,(/row_length,rows/))
                  CASE(348)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASCONDSUNUCSOL,(/row_length,rows/))
                  CASE(349)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASCONDSUAITSOL,(/row_length,rows/))
                  CASE(350)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASCONDSUACCSOL,(/row_length,rows/))
                  CASE(351)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASCONDSUCORSOL,(/row_length,rows/))
                  CASE(358)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASNUCLSUNUCSOL,(/row_length,rows/))
                  CASE(359)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASCOAGSUINTR12,(/row_length,rows/))
                  CASE(360)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASCOAGSUINTR13,(/row_length,rows/))
                  CASE(361)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASCOAGSUINTR14,(/row_length,rows/))
                  CASE(367)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASCOAGSUINTR23,(/row_length,rows/))
                  CASE(371)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASCOAGSUINTR24,(/row_length,rows/))
                  CASE(375)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASCOAGSUINTR34,(/row_length,rows/))
                  CASE(378)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASCOAGSSINTR34,(/row_length,rows/))
                  CASE(384)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASMERGSUINTR12,(/row_length,rows/))
                  CASE(386)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASMERGSUINTR23,(/row_length,rows/))
                  CASE(390)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASMERGSUINTR34,(/row_length,rows/))
                  CASE(391)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(MASMERGSSINTR34,(/row_length,rows/))
                  CASE(395)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(DRYDP(:,1),(/row_length,rows/))
                 CASE(396)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(DRYDP(:,2),(/row_length,rows/))
                  CASE(397)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(DRYDP(:,3),(/row_length,rows/))
                  CASE(398)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(DRYDP(:,4),(/row_length,rows/))
                  CASE(402)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(WETDP(:,1),(/row_length,rows/))
                  CASE(403)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(WETDP(:,2),(/row_length,rows/))
                  CASE(404)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(WETDP(:,3),(/row_length,rows/))
                  CASE(405)
                  mode_diags(:,:,l,k)=                                  &
                     RESHAPE(WETDP(:,4),(/row_length,rows/))
                  CASE DEFAULT
                   cmessage=' Item not found in CASE statement'
! DEPENDS ON: ereport
                   CALL EREPORT('UKCA_AERO_CTL',                        &
                                UkcaD1codes(N)%item,cmessage)
                END SELECT
                if (verbose > 0 .AND. l == 1) then      ! DEBUG
                write(6,'(A16,3i4,3e14.4)') 'UKCA_MODE diag: ',         &
                              k,N,ukcaD1codes(N)%item,                  &
                              sum(mode_diags(:,:,l,k)),                 &
                           maxval(mode_diags(:,:,l,k)),                 &
                           minval(mode_diags(:,:,l,k))
                endif
            ENDIF    ! item >=301 etc
          ENDDO      ! nukca_D1items
        ENDIF        ! L_UKCA_MODE_diags

      ENDDO ! loop over p_levelsda

      IF (VERBOSE.GE.1) THEN

        WRITE(6,*) ' Tracers at end of UKCA_MODE:'
        DO i=1,2          !p_levelsda
         DO j=1,n_sulf_tracers
         WRITE(6,*) 'Level: ',i,' Tracer: ',j
         WRITE(6,*) 'sulf_tracers: ',minval(sulf_tracers(:,:,i,j)),     &
                 maxval(sulf_tracers(:,:,i,j)),                         &
                 sum(sulf_tracers(:,:,i,j))/                            &
                 real(size(sulf_tracers(:,:,i,j)))
         ENDDO
         DO j=1,n_mode_tracers
          IF (mode_tracer_debug(j)) then
            WRITE(6,*) 'Level: ',i,' Tracer: ',j,mode_tracer_names(j)
            WRITE(6,*) 'mode_tracers: ',minval(mode_tracers(:,:,i,j)),  &
                   maxval(mode_tracers(:,:,i,j)),                       &
                   sum(mode_tracers(:,:,i,j))/                          &
                   real(size(mode_tracers(:,:,i,j)))
          ENDIF
         ENDDO
         WRITE(6,*) 'Number of merges and dpfixes for Level: ',i
         DO j=1,nmodes
          WRITE(6,*) j,sum(N_MERGE_3D(:,:,i,j)),                        &
                     sum(N_DPFIX_3D(:,:,i,j))
         ENDDO
       ENDDO      ! i

      ENDIF

      write(6,'(1a36,2i6,2e12.3)')                                      &
           'Total Number of merges and dpfixes=:',                      &
            sum(N_MERGE_3D),sum(N_DPFIX_3D),                            &
       real(sum(N_MERGE_3D))/real(size(N_MERGE_3D)),                    &
       real(sum(N_DPFIX_3D))/real(size(N_DPFIX_3D))
      DO j=1,nmodes
       write(6,'(3i6,2e12.3)') j,                                       &
            sum(N_MERGE_3D(:,:,:,j)),sum(N_DPFIX_3D(:,:,:,j)),          &
       real(sum(N_MERGE_3D(:,:,:,j)))/real(size(N_MERGE_3D(:,:,:,j))),  &
       real(sum(N_DPFIX_3D(:,:,:,j)))/real(size(N_DPFIX_3D(:,:,:,j)))
      ENDDO

      RETURN
      END SUBROUTINE UKCA_AERO_CTL
#endif
