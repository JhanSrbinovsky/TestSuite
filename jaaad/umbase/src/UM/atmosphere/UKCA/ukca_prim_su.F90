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
!    Calculates primary particulate sulfate emissions.
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
      SUBROUTINE UKCA_PRIM_SU(NBOX,                                     &
       ND,MDT,MD,EMANSO2,EMVOLCONSO2,NEMVOLCONSO2,                      &
       EMVOLEXPSO2,NEMVOLEXPSO2,EMBIOMSO2G,                             &
       DTC,SM,ISO2EMS,AIRD,BUD_AER_MAS)

!--------------------------------------------------------
!
!     Calculates primary particulate sulfate emissions
!
!     Assumes particulate emissions emitted as lognormal modes
!     with geometric number mean diameters of GMDIAM and geometric
!     standard deviations of GSD.
!
!     Inputs
!     -------
!     NBOX       : Number of grid boxes
!     ND         : Aerosol ptcl number density for mode (cm^-3)
!     MDT        : Avg tot mass of aerosol ptcl in mode (particle^-1)
!     MD         : Avg cpt mass of aerosol ptcl in mode (particle^-1)
!     EMANSO2    : Anthrop. SO2 ems rates (N_EM_TYPE) (kgSO2/box/s)
!     EMVOLCONSO2(NBOX,10) Natural SO2 ems rates (from cont. volc. src)
!                                                       (in kgSO2/box/s)
!     NEMVOLCONSO2(NBOX)   # of 1x1 cont. volc. SO2 ems sources in box
!     EMVOLEXPSO2(NBOX,10) Natural SO2 ems rates (from expl. volc. src)
!                                                       (in kgSO2/box/s)
!     NEMVOLEXPSO2(NBOX)   # of 1x1 expl. volc. SO2 ems sources in box
!     EMBIOMSO2G(NBOX,6) Biomass burning SO2 ems rates (kgSO2/box/s)
!                                                    [6 altitude ranges]
!     DTC        : Chemical time step (s)
!     SM         : Grid box mass of air (kg)
!     ISO2EMS    : Switch for size distribution of primary sulfate ems
!     AIRD       : Number density of air (per cc)
!
!     Outputs
!     -------
!     ND,MD,MDT  : Updated no. conc, total avg mass, cpt avg mass.
!     BUD_AER_MAS:
!
!     Local variables
!     ---------------
!     MMAIR               : Molar mass of dry air (kg/mol)
!     PARFRAC             : Fraction of anthrop. S mass emitted as ptcls
!     N_EM_TYPES          : Number of emission types for ISO2EMS setting
!     N_EM_MODES          : Number of emission modes for ISO2EMS setting
!     EM_TYPE_MODE_GMDIAM : Geom. mean diam. of ems type & mode (nm)
!     EM_TYPE_MODE_GSD    : Geom. st. dev. of ems type & mode
!     EM_TYPE_MODE_FRAC   : Fraction of ems type which goes to this mode
!     EM_TYPE_MODE_MODE   : Index of mode to emit ems type & mode into
!     N_EM_BIOMSO2_MODES  : Number of emission modes for biomass SO4
!     EM_BIOMSO2_MODE_MODE: Index of mode to emit biomass  SO4 into
!     EM_BIOMSO2_MODE_GMDIAM: Geom. mean diam. of each biomass SO4 mode
!     EM_BIOMSO2_MODE_GSD : Geom. std. dev. of each biomass SO4 mode
!     N_EM_VOLSO2_MODES   : Number of emission modes for volcanic SO4
!     EM_VOLSO2_MODE_GMDIAM: Geom. mean diam. of each volcanic SO4 mode
!     EM_VOLSO2_MODE_GSD  : Geom. std. dev. of each volcanic SO4 mode
!     EM_VOLSO2_MODE_MODE : Index of mode to emit volcanic SO4 into
!     EM_TYPE_MASS        : Mass emitted in ems type (kgH2SO4/box/s)
!     MODE_EM        : Mode into which emit into
!     LGSD           : Natural log of geometric standard deviation
!     NEWN,DELN      : New and change in no. conc (per cc)
!     PARMASS        : Mass emitted to box (type) (kgH2SO4/box/s)
!     MODEMASS       : Mass emitted to box (type & mode) (kgH2SO4/box/s)
!     MODEVOL        : Vol. conc. emitted (type & mode) (nm3/box/s)
!     TOTNUMMODE     : No. conc. of ptcls emitted (type & mode) (/box/s)
!     FACTOR         : Converts from /gridbox/s to /cc/tstep
!
!     Inputted by module UKCA_CONSTANTS
!     ---------------------------------
!     PPI        : 3.1415927
!     AVC        : Avogadro's constant (molecules per mole)
!     RA         : Dry air gas constant = 287.05 Jkg^-1 K^-1
!     ZBOLTZ     : Stefan-Boltzmann constant (kg m2 s-2 K-1 molec-1)
!     DN_EPS     : Value of DELN below which do not carry out process
!     EMS_EPS    : Ems flux value below which don't emit (kg/gridbox/s)
!
!     Inputted by module UKCA_MODE_SETUP
!     ----------------------------------
!     NMODES     : Number of modes set
!     NCP        : Number of components set
!     MODE       : Logical defining which modes are set.
!     COMPONENT  : Logical defining which cpts are in which mode
!     MM         : Molar masses of condensable components (kg/mol)
!     RHOCOMP    : Density of each of the aerosol components (kg/m3)
!     NUM_EPS    : Value of NEWN below which don't recalculate MD
!                                                  or carry out process
!     CP_SU      : index of cpt which sulfate ptcl mass is stored
!
!     Inputted by module UKCA_SETUP_INDICES
!     -------------------------------------
!     NCHEMG     : Number of species for gas phase chemistry scheme
!     MSOTWO     : Index of MOLWT, WTRATC and S0G for SO2
!     MM_GAS     : Array of molar masses for gas phase species (kg/mol)
!     Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES

      IMPLICIT NONE
!
!     Subroutine interface
      INTEGER :: NBOX
      INTEGER :: ISO2EMS
      INTEGER :: NEMVOLCONSO2(NBOX)
      INTEGER :: NEMVOLEXPSO2(NBOX)
      REAL    :: MD(NBOX,NMODES,NCP)
      REAL    :: ND(NBOX,NMODES)
      REAL    :: MDT(NBOX,NMODES)
      REAL    :: DTC
      REAL    :: SM(NBOX)
      REAL    :: AIRD(NBOX)
      REAL    :: EMANSO2(NBOX,6)
      REAL    :: EMVOLCONSO2(NBOX,10)
      REAL    :: EMVOLEXPSO2(NBOX,10)
      REAL    :: EMBIOMSO2G(NBOX,6)
      REAL    :: BUD_AER_MAS(NBOX,NBUDAER)
!
!     Local variables
      INTEGER :: JL
      INTEGER :: IMODE
      INTEGER :: ICP
      INTEGER :: IA
      INTEGER :: IB
      INTEGER :: EM_TYPE_MODE_MODE(3,2)
      INTEGER :: EM_BIOMSO2_MODE_MODE(2)
      INTEGER :: EM_VOLSO2_MODE_MODE(2)
      INTEGER :: N_EM_TYPES
      INTEGER :: N_EM_MODES
      INTEGER :: N_EM_BIOMSO2_MODES
      INTEGER :: N_EM_VOLSO2_MODES
      INTEGER :: MODE_EM
      REAL    :: MMAIR
      REAL    :: PARFRAC
      REAL    :: TOTNUMMODE
      REAL    :: PARMASS
      REAL    :: MODEMASS
      REAL    :: MODEVOL
      REAL    :: LGSD
      REAL    :: FACTOR
      REAL    :: NEWN
      REAL    :: DELN
      REAL    :: EM_TYPE_MODE_GMDIAM(3,2)
      REAL    :: EM_TYPE_MODE_GSD(3,2)
      REAL    :: EM_TYPE_MODE_FRAC(3,2)
      REAL    :: EM_TYPE_MASS(3)
      REAL    :: EM_BIOMSO2_MODE_GMDIAM(2)
      REAL    :: EM_BIOMSO2_MODE_GSD(2)
      REAL    :: EM_BIOMSO2_MODE_FRAC(2)
      REAL    :: EM_VOLSO2_MODE_GMDIAM(2)
      REAL    :: EM_VOLSO2_MODE_GSD(2)
      REAL    :: EM_VOLSO2_MODE_FRAC(2)
!
      MMAIR=AVC*ZBOLTZ/RA
!
      DO JL=1,NBOX
       IF(ISO2EMS == 1) THEN
! Follow BS95 as in GLOMAP version from Spracklen et al (2005)
!
! .. type 1 is GEIA SO2 ems from surface  src, size:BS95, alt:0-100m
! .. type 2 is GEIA SO2 ems from elevated src, size:BS95, alt:level 2

        PARFRAC=0.03
        N_EM_TYPES=2
        N_EM_MODES=2
        EM_TYPE_MODE_GMDIAM(1,1)=10.0
        EM_TYPE_MODE_GMDIAM(1,2)=70.0
        EM_TYPE_MODE_GMDIAM(2,1)=10.0
        EM_TYPE_MODE_GMDIAM(2,2)=70.0
        EM_TYPE_MODE_GSD(1,1)=1.6
        EM_TYPE_MODE_GSD(1,2)=2.0
        EM_TYPE_MODE_GSD(2,1)=1.6
        EM_TYPE_MODE_GSD(2,2)=2.0
        EM_TYPE_MODE_FRAC(1,1)=0.15
        EM_TYPE_MODE_FRAC(1,2)=0.85
        EM_TYPE_MODE_FRAC(2,1)=0.15
        EM_TYPE_MODE_FRAC(2,2)=0.85
        EM_TYPE_MODE_MODE(1,1)=2
        EM_TYPE_MODE_MODE(1,2)=2
        EM_TYPE_MODE_MODE(2,1)=2
        EM_TYPE_MODE_MODE(2,2)=2
        EM_TYPE_MASS(1)=EMANSO2(JL,1)*PARFRAC*MM(CP_SU)/MM_GAS(MSOTWO)
        EM_TYPE_MASS(2)=EMANSO2(JL,2)*PARFRAC*MM(CP_SU)/MM_GAS(MSOTWO)
! .. biomass burning not done for ISO2EMS=1 (NEMBIOMSO2G=0)
        N_EM_VOLSO2_MODES=2
        EM_VOLSO2_MODE_GMDIAM(1)=30.0
        EM_VOLSO2_MODE_GMDIAM(2)=80.0
        EM_VOLSO2_MODE_GSD(1)=1.8
        EM_VOLSO2_MODE_GSD(2)=1.8
        EM_VOLSO2_MODE_FRAC(1)=0.5
        EM_VOLSO2_MODE_FRAC(2)=0.5
        EM_VOLSO2_MODE_MODE(1)=2
        EM_VOLSO2_MODE_MODE(2)=2
       ENDIF
       IF(ISO2EMS == 2) THEN
! Follow AEROCOM recommendations for ACB
!
! .. type 1: roads, off-road & domestic, size:TRAFFIC, alt:lowest layer
! .. type 2: industrial, power-plants, size:INDUSTR, alt:100-300m
! .. type 3: shipping, size:INDUSTR, alt:lowest layer

        PARFRAC=0.025
        N_EM_TYPES=3
        N_EM_MODES=1
        EM_TYPE_MODE_GMDIAM(1,1)=30.0
        EM_TYPE_MODE_GMDIAM(2,1)=1000.0
        EM_TYPE_MODE_GMDIAM(3,1)=1000.0
        EM_TYPE_MODE_GSD(1,1)=1.8
        EM_TYPE_MODE_GSD(2,1)=2.0
        EM_TYPE_MODE_GSD(3,1)=2.0
        EM_TYPE_MODE_FRAC(1,1)=1.0
        EM_TYPE_MODE_FRAC(2,1)=1.0
        EM_TYPE_MODE_FRAC(3,1)=1.0
        EM_TYPE_MODE_MODE(1,1)=2
        EM_TYPE_MODE_MODE(2,1)=3
        EM_TYPE_MODE_MODE(3,1)=3
        EM_TYPE_MASS(1)=(EMANSO2(JL,3)+EMANSO2(JL,5)+EMANSO2(JL,6))*    &
                         PARFRAC*MM(CP_SU)/MM_GAS(MSOTWO)
        EM_TYPE_MASS(2)=(EMANSO2(JL,1)+EMANSO2(JL,2))*                  &
                         PARFRAC*MM(CP_SU)/MM_GAS(MSOTWO)
        EM_TYPE_MASS(3)=EMANSO2(JL,4)*PARFRAC*MM(CP_SU)/MM_GAS(MSOTWO)
        N_EM_BIOMSO2_MODES=1
        EM_BIOMSO2_MODE_GMDIAM(1)=80.0
        EM_BIOMSO2_MODE_GSD(1)=1.8
        EM_BIOMSO2_MODE_FRAC(1)=1.0
        EM_BIOMSO2_MODE_MODE(1)=2
        N_EM_VOLSO2_MODES=2
        EM_VOLSO2_MODE_GMDIAM(1)=30.0
        EM_VOLSO2_MODE_GMDIAM(2)=80.0
        EM_VOLSO2_MODE_GSD(1)=1.8
        EM_VOLSO2_MODE_GSD(2)=1.8
        EM_VOLSO2_MODE_FRAC(1)=0.5
        EM_VOLSO2_MODE_FRAC(2)=0.5
        EM_VOLSO2_MODE_MODE(1)=2
        EM_VOLSO2_MODE_MODE(2)=2
       ENDIF
       IF((ISO2EMS == 3).OR.(ISO2EMS == 5).OR.(ISO2EMS == 6)) THEN
! Follow Stier05 modifications to AEROCOM ACB
!
! .. type 1 is roads, off-road & domestic, size:modifiedTRAFFIC,
! ..                                       alt:lowest layer
! .. type 2 is industrial, power-plants  , size:modifiedINDUSTR,
! ..                                       alt:100-300m
! .. type 3 is shipping, size:modifiedINDUSTR, alt:lowest layer

        PARFRAC=0.025
        N_EM_TYPES=3
        N_EM_MODES=2
        EM_TYPE_MODE_GMDIAM(1,1)=60.0
        EM_TYPE_MODE_GMDIAM(1,2)=150.0
        EM_TYPE_MODE_GMDIAM(2,1)=150.0
        EM_TYPE_MODE_GMDIAM(2,2)=1500.0
        EM_TYPE_MODE_GMDIAM(3,1)=150.0
        EM_TYPE_MODE_GMDIAM(3,2)=1500.0
        EM_TYPE_MODE_GSD(1,1)=1.59
        EM_TYPE_MODE_GSD(1,2)=1.59
        EM_TYPE_MODE_GSD(2,1)=1.59
        EM_TYPE_MODE_GSD(2,2)=2.0
        EM_TYPE_MODE_GSD(3,1)=1.59
        EM_TYPE_MODE_GSD(3,2)=2.0
        EM_TYPE_MODE_FRAC(1,1)=0.5
        EM_TYPE_MODE_FRAC(1,2)=0.5
        EM_TYPE_MODE_FRAC(2,1)=0.5
        EM_TYPE_MODE_FRAC(2,2)=0.5
        EM_TYPE_MODE_FRAC(3,1)=0.5
        EM_TYPE_MODE_FRAC(3,2)=0.5
        EM_TYPE_MODE_MODE(1,1)=2
        EM_TYPE_MODE_MODE(1,2)=3
        EM_TYPE_MODE_MODE(2,1)=3
        EM_TYPE_MODE_MODE(2,2)=4
        EM_TYPE_MODE_MODE(3,1)=3
        EM_TYPE_MODE_MODE(3,2)=4
        EM_TYPE_MASS(1)=(EMANSO2(JL,3)+EMANSO2(JL,5)+EMANSO2(JL,6))*    &
                         PARFRAC*MM(CP_SU)/MM_GAS(MSOTWO)
        EM_TYPE_MASS(2)=(EMANSO2(JL,1)+EMANSO2(JL,2))*                  &
                         PARFRAC*MM(CP_SU)/MM_GAS(MSOTWO)
        EM_TYPE_MASS(3)=EMANSO2(JL,4)*PARFRAC*MM(CP_SU)/MM_GAS(MSOTWO)
        N_EM_BIOMSO2_MODES=2
        EM_BIOMSO2_MODE_GMDIAM(1)=60.0
        EM_BIOMSO2_MODE_GMDIAM(2)=150.0
        EM_BIOMSO2_MODE_GSD(1)=1.59
        EM_BIOMSO2_MODE_GSD(2)=1.59
        EM_BIOMSO2_MODE_FRAC(1)=0.5
        EM_BIOMSO2_MODE_FRAC(2)=0.5
        EM_BIOMSO2_MODE_MODE(1)=2
        EM_BIOMSO2_MODE_MODE(2)=3
        N_EM_VOLSO2_MODES=2
        EM_VOLSO2_MODE_GMDIAM(1)=60.0
        EM_VOLSO2_MODE_GMDIAM(2)=150.0
        EM_VOLSO2_MODE_GSD(1)=1.59
        EM_VOLSO2_MODE_GSD(2)=1.59
        EM_VOLSO2_MODE_FRAC(1)=0.5
        EM_VOLSO2_MODE_FRAC(2)=0.5
        EM_VOLSO2_MODE_MODE(1)=2
        EM_VOLSO2_MODE_MODE(2)=3
       ENDIF
       IF(ISO2EMS == 4) THEN
! GEIA Emissions with AEROCOM recommendations for ACB
!
! .. type 1 is low level emissions, size:TRAFFIC, alt:lowest layer
! .. type 2 is high level emissions, size:INDUSTR, alt:100-300m

        PARFRAC=0.025
        N_EM_TYPES=2
        N_EM_MODES=2
        EM_TYPE_MODE_GMDIAM(1,1)=60.0
        EM_TYPE_MODE_GMDIAM(1,2)=150.0
        EM_TYPE_MODE_GMDIAM(2,1)=150.0
        EM_TYPE_MODE_GMDIAM(2,2)=1500.0
        EM_TYPE_MODE_GSD(1,1)=1.59
        EM_TYPE_MODE_GSD(1,2)=1.59
        EM_TYPE_MODE_GSD(2,1)=1.59
        EM_TYPE_MODE_GSD(2,2)=2.0
        EM_TYPE_MODE_FRAC(1,1)=0.5
        EM_TYPE_MODE_FRAC(1,2)=0.5
        EM_TYPE_MODE_FRAC(2,1)=0.5
        EM_TYPE_MODE_FRAC(2,2)=0.5
        EM_TYPE_MODE_MODE(1,1)=2
        EM_TYPE_MODE_MODE(1,2)=3
        EM_TYPE_MODE_MODE(2,1)=3
        EM_TYPE_MODE_MODE(2,2)=4
        EM_TYPE_MASS(1)=EMANSO2(JL,1)*PARFRAC*MM(CP_SU)/MM_GAS(MSOTWO)
        EM_TYPE_MASS(2)=EMANSO2(JL,2)*PARFRAC*MM(CP_SU)/MM_GAS(MSOTWO)
        N_EM_BIOMSO2_MODES=1
        EM_BIOMSO2_MODE_GMDIAM(1)=80.0
        EM_BIOMSO2_MODE_GSD(1)=1.8
        EM_BIOMSO2_MODE_FRAC(1)=1.0
        EM_BIOMSO2_MODE_MODE(1)=2
        N_EM_VOLSO2_MODES=2
        EM_VOLSO2_MODE_GMDIAM(1)=30.0
        EM_VOLSO2_MODE_GMDIAM(2)=80.0
        EM_VOLSO2_MODE_GSD(1)=1.8
        EM_VOLSO2_MODE_GSD(2)=1.8
        EM_VOLSO2_MODE_FRAC(1)=0.5
        EM_VOLSO2_MODE_FRAC(2)=0.5
        EM_VOLSO2_MODE_MODE(1)=2
        EM_VOLSO2_MODE_MODE(2)=2
       ENDIF
!
!
       FACTOR=DTC*AIRD(JL)*MMAIR/(SM(JL)*AVC)  ! converts from
                                               ! per gridbox/s to
                                               ! per cc per timestep

!----------------------------------------------------------------------
!      This section does primary SO4 from SO2 from anthropogenic src

       DO IA=1,N_EM_TYPES

        PARMASS=EM_TYPE_MASS(IA) ! kgH2SO4/box/s

        IF(PARMASS > EMS_EPS) THEN

         DO IB=1,N_EM_MODES

! .. Calulate natural logs of standard deviation
          LGSD=LOG(EM_TYPE_MODE_GSD(IA,IB))
!
! .. Particulate mass emissions (kg_H2SO4 per gridbox per s) in mode
          MODEMASS=PARMASS*EM_TYPE_MODE_FRAC(IA,IB)
!
! .. Calculate total particle volume (nm3 per gridbox per s)
          MODEVOL=1E27*MODEMASS/RHOCOMP(CP_SU)
!
! .. Calculate total particle number (per gridbox per s)
          TOTNUMMODE=MODEVOL/((PPI/6.0)*                                &
             (EM_TYPE_MODE_GMDIAM(IA,IB)**3.0)*EXP(4.5*LGSD*LGSD))
!
! .. Calculate change in number conc for emission type and mode
! .. FACTOR converts from number/gridbox/s to number/cc/tstep
          DELN=FACTOR*TOTNUMMODE

! .. Store which mode to emit primary SO4 into
          MODE_EM=EM_TYPE_MODE_MODE(IA,IB)

! .. Calculate new particle number concetration
          NEWN=ND(JL,MODE_EM)+DELN

          IF((DELN > DN_EPS).AND.(NEWN > NUM_EPS(MODE_EM))) THEN

           DO ICP=1,NCP
            IF(COMPONENT(MODE_EM,ICP)) THEN

! .. Update cpt masses per particle from SO4 from anthropogenic src
             IF(ICP == CP_SU) THEN
              MD(JL,MODE_EM,ICP)=(MD(JL,MODE_EM,ICP)*ND(JL,MODE_EM)+    &
                MODEMASS*FACTOR*AVC/MM(ICP))/NEWN
! .. Add fraction of ems type ems for mode to ptcl mass in chosen mode
! .. (FACTOR*AVC/MM(ICP)) converts kg/gridbox/s to molecules/cc/tstep
              IF((MODE_EM == 2).AND.(NMASPRIMSUAITSOL > 0))             &
                               BUD_AER_MAS(JL,NMASPRIMSUAITSOL)=        &
                               BUD_AER_MAS(JL,NMASPRIMSUAITSOL)+        &
                                  MODEMASS*FACTOR*AVC/MM(ICP)
              IF((MODE_EM == 3).AND.(NMASPRIMSUACCSOL > 0))             &
                               BUD_AER_MAS(JL,NMASPRIMSUACCSOL)=        &
                               BUD_AER_MAS(JL,NMASPRIMSUACCSOL)+        &
                                  MODEMASS*FACTOR*AVC/MM(ICP)
              IF((MODE_EM == 4).AND.(NMASPRIMSUCORSOL > 0))             &
                               BUD_AER_MAS(JL,NMASPRIMSUCORSOL)=        &
                               BUD_AER_MAS(JL,NMASPRIMSUCORSOL)+        &
                                  MODEMASS*FACTOR*AVC/MM(ICP)
             ELSE
              MD(JL,MODE_EM,ICP)=MD(JL,MODE_EM,ICP)*ND(JL,MODE_EM)/NEWN
             ENDIF

            ENDIF
           ENDDO ! end loop over cpts

! .. Update number concentration from SO4 from anthropogenic src
           ND(JL,MODE_EM)=NEWN

          ENDIF ! if DELN>DN_EPS and NEWN>NUM_EPS
!
         ENDDO
!
        ENDIF ! end of if PARMASS>EMS_EPS
!
       ENDDO ! end loop over types
!
!----------------------------------------------------------------------
!      This section does primary SO4 from SO2 from biomass burning src

       DO IA=1,6 ! loop over 6 emission levels

        PARMASS=EMBIOMSO2G(JL,IA)*PARFRAC*                              &
                    MM(CP_SU)/MM_GAS(MSOTWO) ! in kgH2SO4/box/s

        IF(PARMASS > EMS_EPS) THEN

         DO IB=1,N_EM_BIOMSO2_MODES ! loop over # of emission modes

! .. Calulate natural logs of standard deviation
          LGSD=LOG(EM_BIOMSO2_MODE_GSD(IB))
!
! .. Particulate mass emissions (kg_H2SO4 per gridbox per s) in mode
          MODEMASS=PARMASS*EM_BIOMSO2_MODE_FRAC(IB)
!
! .. Calculate total particle volume (nm3 per gridbox per s)
          MODEVOL=1E27*MODEMASS/RHOCOMP(CP_SU)
!
! .. Calculate total particle number (per gridbox per s)
          TOTNUMMODE=MODEVOL/((PPI/6.0)*                                &
             (EM_BIOMSO2_MODE_GMDIAM(IB)**3.0)*EXP(4.5*LGSD*LGSD))
!
! .. Calculate change in number conc for emission type and mode
! .. FACTOR converts from number/gridbox/s to number/cc/tstep
          DELN=FACTOR*TOTNUMMODE

! .. Store which mode to emit primary SO4 into
          MODE_EM=EM_BIOMSO2_MODE_MODE(IB)

! .. Calculate new particle number concentration
          NEWN=ND(JL,MODE_EM)+DELN

          IF((DELN > DN_EPS).AND.(NEWN > NUM_EPS(MODE_EM))) THEN

           DO ICP=1,NCP
            IF(COMPONENT(MODE_EM,ICP)) THEN

! .. Update cpt masses per particle from SO4 from biomass burning src
             IF(ICP == CP_SU) THEN
              MD(JL,MODE_EM,ICP)=(MD(JL,MODE_EM,ICP)*ND(JL,MODE_EM)+    &
                MODEMASS*FACTOR*AVC/MM(ICP))/NEWN
! .. Add fraction of ems type ems for mode to ptcl mass in chosen mode
! .. (FACTOR*AVC/MM(ICP)) converts kg/gridbox/s to molecules/cc/tstep
              IF((MODE_EM == 2).AND.(NMASPRIMSUAITSOL > 0))             &
                               BUD_AER_MAS(JL,NMASPRIMSUAITSOL)=        &
                               BUD_AER_MAS(JL,NMASPRIMSUAITSOL)+        &
                                  MODEMASS*FACTOR*AVC/MM(ICP)
              IF((MODE_EM == 3).AND.(NMASPRIMSUACCSOL > 0))             &
                               BUD_AER_MAS(JL,NMASPRIMSUACCSOL)=        &
                               BUD_AER_MAS(JL,NMASPRIMSUACCSOL)+        &
                                  MODEMASS*FACTOR*AVC/MM(ICP)
              IF((MODE_EM == 4).AND.(NMASPRIMSUCORSOL > 0))             &
                               BUD_AER_MAS(JL,NMASPRIMSUCORSOL)=        &
                               BUD_AER_MAS(JL,NMASPRIMSUCORSOL)+        &
                                  MODEMASS*FACTOR*AVC/MM(ICP)
             ELSE
              MD(JL,MODE_EM,ICP)=MD(JL,MODE_EM,ICP)*ND(JL,MODE_EM)/NEWN
             ENDIF
            ENDIF
           ENDDO ! end loop over cpts

! .. Update number concentration from SO4 from biomass burning src
           ND(JL,MODE_EM)=NEWN

          ENDIF ! if DELN>DN_EPS and NEWN>NUM_EPS
!
         ENDDO
!
        ENDIF ! if PARMASS>CONC_EPS
       ENDDO
!
!----------------------------------------------------------------------
!      This section does primary SO4 from SO2 from volcanic (cont.) src

       DO IA=1,NEMVOLCONSO2(JL)

        PARMASS=EMVOLCONSO2(JL,IA)*PARFRAC*                             &
                     MM(CP_SU)/MM_GAS(MSOTWO) ! in kgH2SO4/box/s

        IF(PARMASS > EMS_EPS) THEN
!
         DO IB=1,N_EM_VOLSO2_MODES

! .. Calulate natural logs of standard deviation
          LGSD=LOG(EM_VOLSO2_MODE_GSD(IB))
!
! .. Particulate mass emissions (kg_H2SO4 per gridbox per s) in mode
          MODEMASS=PARMASS*EM_VOLSO2_MODE_FRAC(IB)
!
! .. Calculate total particle volume (nm3 per gridbox per s)
          MODEVOL=1E27*MODEMASS/RHOCOMP(CP_SU)
!
! .. Calculate total particle number (per gridbox per s)
          TOTNUMMODE=MODEVOL/((PPI/6.0)*                                &
             (EM_VOLSO2_MODE_GMDIAM(IB)**3.0)*EXP(4.5*LGSD*LGSD))
!
! .. Calculate change in number conc for emission type and mode
! .. FACTOR converts from number/gridbox/s to number/cc/tstep
          DELN=FACTOR*TOTNUMMODE

! .. Store which mode to emit primary SO4 into
          MODE_EM=EM_VOLSO2_MODE_MODE(IB)

! .. Calculate new particle number concetration
          NEWN=ND(JL,MODE_EM)+DELN

          IF((DELN > DN_EPS).AND.(NEWN > NUM_EPS(MODE_EM))) THEN

           DO ICP=1,NCP
            IF(COMPONENT(MODE_EM,ICP)) THEN

! .. Update cpt masses per particle from SO4 from volcanic (cont.) src
             IF(ICP == CP_SU) THEN
              MD(JL,MODE_EM,ICP)=(MD(JL,MODE_EM,ICP)*ND(JL,MODE_EM)+    &
                MODEMASS*FACTOR*AVC/MM(ICP))/NEWN
! .. Add fraction of ems type ems for mode to ptcl mass in chosen mode
! .. (FACTOR*AVC/MM(ICP)) converts kg/gridbox/s to molecules/cc/tstep
              IF((MODE_EM == 2).AND.(NMASPRIMSUAITSOL > 0))             &
                               BUD_AER_MAS(JL,NMASPRIMSUAITSOL)=        &
                               BUD_AER_MAS(JL,NMASPRIMSUAITSOL)+        &
                                  MODEMASS*FACTOR*AVC/MM(ICP)
              IF((MODE_EM == 3).AND.(NMASPRIMSUACCSOL > 0))             &
                               BUD_AER_MAS(JL,NMASPRIMSUACCSOL)=        &
                               BUD_AER_MAS(JL,NMASPRIMSUACCSOL)+        &
                                  MODEMASS*FACTOR*AVC/MM(ICP)
              IF((MODE_EM == 4).AND.(NMASPRIMSUCORSOL > 0))             &
                               BUD_AER_MAS(JL,NMASPRIMSUCORSOL)=        &
                               BUD_AER_MAS(JL,NMASPRIMSUCORSOL)+        &
                                  MODEMASS*FACTOR*AVC/MM(ICP)
             ELSE
              MD(JL,MODE_EM,ICP)=MD(JL,MODE_EM,ICP)*ND(JL,MODE_EM)/NEWN
             ENDIF
            ENDIF
           ENDDO ! end loop over cpts

! .. Update number concentration from SO4 from volcanic (cont.) src
           ND(JL,MODE_EM)=NEWN

          ENDIF ! if DELN>DN_EPS and NEWN>NUM_EPS
!
         ENDDO
!
        ENDIF ! if PARMASS>CONC_EPS
       ENDDO
!
!----------------------------------------------------------------------
!      This section does primary SO4 from SO2 from volcanic (expl.) src

       DO IA=1,NEMVOLEXPSO2(JL)

        PARMASS=EMVOLEXPSO2(JL,IA)*PARFRAC*                             &
                     MM(CP_SU)/MM_GAS(MSOTWO) ! in kgH2SO4/box/s

        IF(PARMASS > EMS_EPS) THEN
!
         DO IB=1,N_EM_VOLSO2_MODES

! .. Calulate natural logs of standard deviation
          LGSD=LOG(EM_VOLSO2_MODE_GSD(IB))
!
! .. Particulate mass emissions (kg_H2SO4 per gridbox per s) in mode
          MODEMASS=PARMASS*EM_VOLSO2_MODE_FRAC(IB)
!
! .. Calculate total particle volume (nm3 per gridbox per s)
          MODEVOL=1E27*MODEMASS/RHOCOMP(CP_SU)
!
! .. Calculate total particle number (per gridbox per s)
          TOTNUMMODE=MODEVOL/((PPI/6.0)*                                &
             (EM_VOLSO2_MODE_GMDIAM(IB)**3.0)*EXP(4.5*LGSD*LGSD))
!
! .. Calculate change in number conc for emission type and mode
! .. FACTOR converts from number/gridbox/s to number/cc/tstep
          DELN=FACTOR*TOTNUMMODE

! .. Store which mode to emit primary SO4 into
          MODE_EM=EM_VOLSO2_MODE_MODE(IB)

! .. Calculate new particle number concetration
          NEWN=ND(JL,MODE_EM)+DELN

          IF((DELN > DN_EPS).AND.(NEWN > NUM_EPS(MODE_EM))) THEN

           DO ICP=1,NCP
            IF(COMPONENT(MODE_EM,ICP)) THEN

! .. Update cpt masses per particle from SO4 from volcanic (expl.) src
             IF(ICP == CP_SU) THEN
              MD(JL,MODE_EM,ICP)=(MD(JL,MODE_EM,ICP)*ND(JL,MODE_EM)+    &
                MODEMASS*FACTOR*AVC/MM(ICP))/NEWN
! .. Add fraction of ems type ems for mode to ptcl mass in chosen mode
! .. (FACTOR*AVC/MM(ICP)) converts kg/gridbox/s to molecules/cc/tstep
              IF((MODE_EM == 2).AND.(NMASPRIMSUAITSOL > 0))             &
                               BUD_AER_MAS(JL,NMASPRIMSUAITSOL)=        &
                               BUD_AER_MAS(JL,NMASPRIMSUAITSOL)+        &
                                  MODEMASS*FACTOR*AVC/MM(ICP)
              IF((MODE_EM == 3).AND.(NMASPRIMSUACCSOL > 0))             &
                               BUD_AER_MAS(JL,NMASPRIMSUACCSOL)=        &
                               BUD_AER_MAS(JL,NMASPRIMSUACCSOL)+        &
                                  MODEMASS*FACTOR*AVC/MM(ICP)
              IF((MODE_EM == 4).AND.(NMASPRIMSUCORSOL > 0))             &
                               BUD_AER_MAS(JL,NMASPRIMSUCORSOL)=        &
                               BUD_AER_MAS(JL,NMASPRIMSUCORSOL)+        &
                                  MODEMASS*FACTOR*AVC/MM(ICP)
             ELSE
              MD(JL,MODE_EM,ICP)=MD(JL,MODE_EM,ICP)*ND(JL,MODE_EM)/NEWN
             ENDIF
            ENDIF
           ENDDO ! end loop over cpts

! .. Update number concentration from SO4 from volcanic (expl.) src
           ND(JL,MODE_EM)=NEWN

          ENDIF ! if DELN>DN_EPS and NEWN>NUM_EPS
!
         ENDDO
!
        ENDIF ! if PARMASS>CONC_EPS
       ENDDO
!----------------------------------------------------------------------
!
! .. Update total mode mass per particle MDT
       DO IMODE=1,NMODES
        IF(MODE(IMODE)) THEN
         MDT(JL,IMODE) = 0.0
         DO ICP=1,NCP
          IF(COMPONENT(IMODE,ICP)) MDT(JL,IMODE)=                       &
                   MDT(JL,IMODE) + MD(JL,IMODE,ICP)
         ENDDO ! end of DO ICP=1,NCP
        ENDIF ! end of IF(MODE(IMODE))
       ENDDO ! end of DO IMODE=1,NMODES
!
      ENDDO ! end of DO JL=1,NBOX
!
      RETURN
      END SUBROUTINE UKCA_PRIM_SU
#endif
