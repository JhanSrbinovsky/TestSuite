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
!    Calculates primary carbonaceous aerosol emissions.
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
      SUBROUTINE UKCA_PRIM_CAR(NBOX,ND,MDT,MD,                          &
       EMC,EMCBM,DTC,SM,AIRD,SURTP,ISO2EMS,BUD_AER_MAS)
!-------------------------------------------------------------
!
!     Calculates primary carbonaceous aerosol emissions
!
!     Assumes Bond et al (2004) emissions (which give mass fluxes
!     of BC and OC for bio-fuel and fossil-fuel emissions on a
!     1x1 degree grid) are for particles consisting of an
!     internal mixture of BC and OC (rather than pure
!     BC particles and pure OC particles).
!
!     EMC  : Carbon Emissions (kgC/box/s)
!      1: Black carbon - Bio fuel
!      2: Black carbon - Fossil fuel
!      3: Organic carbon - Bio fuel
!      4: Organic carbon - Fossil fuel
!
!     EMCBM : Carbon emissions from biomass burning (into heights) (kgC/box/s)
!    1,1: Black carbon - biomass burning in 1st height interval (0-100m)
!    1,2: Black carbon - biomass burning in 2nd height interval (100-500m)
!    1,3: Black carbon - biomass burning in 3rd height interval (500-1000m)
!    1,4: Black carbon - biomass burning in 4th height interval (1000-2000m)
!    1,5: Black carbon - biomass burning in 5th height interval (2000-3000m)
!    1,6: Black carbon - biomass burning in 6th height interval (3000-6000m)
!    2,1: Organic carbon - biomass burning in 1st height interval (0-100m)
!    2,2: Organic carbon - biomass burning in 2nd height interval (100-500m)
!    2,3: Organic carbon - biomass burning in 3rd height interval (500-1000m)
!    2,4: Organic carbon - biomass burning in 4th height interval (1000-2000m)
!    2,5: Organic carbon - biomass burning in 5th height interval (2000-3000m)
!    2,6: Organic carbon - biomass burning in 6th height interval (3000-6000m)
!
!     Currently assume all BC and OC insoluble at emission
!
!     Parameters
!     ----------
!     MODE_DIAM  : Geometric mean diameter of modes into which emit bf/ff ems
!     STDEV      : Standard deviation of modes into which emit bf/ff ems
!
!     Inputs
!     -----
!     NBOX       : Number of grid boxes
!     ND         : Aerosol ptcl number density for mode (cm^-3)
!     MDT        : Avg tot mass of aerosol ptcl in mode (particle^-1)
!     MD         : Avg cpt mass of aerosol ptcl in mode (particle^-1)
!     EMC        : BC/OC ems rates [bio- & fossil-fuel srcs] (kgC/box/s)
!     EMCBM      : BC/OC ems rates [biomass burning srcs] (kgC/box/s)
!     DTC        : Chemical time step (s)
!     SM         : Grid box mass of air (kg)
!     AIRD       : Number density of air (per cc)
!     SURTP      : Surface type [0/1=at sea/land-surf,2/3=above-sea/land]
!     ISO2EMS    : Switch for choice of SO2 emissions
!
!     Outputs
!     -------
!     ND,MD,MDT  : Updated no. conc, total avg mass, cpt avg mass.
!     BUD_AER_MAS: Updated budget mass fluxes including BC/OC emissions
!
!     Local variables
!     ---------------
!     MMAIR      : Molar mass of dry air
!     MODE_DIAM  : Geom. mean diam. of modes into which emit bf/ff ems
!     STDEV      : Geom. std dev. of modes into which emit bf/ff ems
!     LSTDEV     : Natural logarithm of STDEV
!     MDNDCP     : Latest cpt mass conc. (molecules/cm3) as updated
!     DELMDNDCP  : Change in cpt mass conc. for particular emission
!     EMCBC      : BC ems rate into this particular mode (kgC/box/s)
!     EMCOC      : OC ems rate into this particular mode (kgC/box/s)
!     TOTEMC     : BC+OC ems rate into this particular mode (kgC/box/s)
!     EMCVOL     : Total (BC+OC) emitted particle volume (nm3/box/s)
!     TOTNUM     : Total emitted particle number (ptcls/box/s)
!     A          : Loops over 1 (bio-fuel & fire) and 2 (fossil-fuel)
!     FACTOR     : Converts from /gridbox/s to /cc/timestep
!     NEWN       : Updated particle number concentration (/cm3)
!     DELN       : Change in particle number concentration (/cm3)
!
!     Inputted by module UKCA_CONSTANTS
!     ---------------------------------
!     PPI        : 3.1415927
!     AVC        : Avogadro's constant (molecules per mole)
!     RA         : Dry air gas constant = 287.05 Jkg^-1 K^-1
!     ZBOLTZ     : Stefan-Boltzmann constant (kg m2 s-2 K-1 molec-1)
!     EMS_EPS    : Value of ems flux below which dont do ems (kg/box/s)
!     DN_EPS     : Value of DELN below which do not carry out process
!
!     Inputted by module UKCA_MODE_SETUP
!     ----------------------------------
!     NMODES     : Number of modes set
!     NCP        : Number of components set
!     MODE       : Logical variable defining which modes are set.
!     COMPONENT  : Logical variable defining which cpt are in which dsts
!     MM         : Molar masses of condensable components (kg/mole)
!     RHOCOMP    : Component mass densities (kg/m3)
!     NUM_EPS    : Value of NEWN below which don't recalculate MD
!                                                  or carry out process
!     CP_BC      : index of cpt into which emit BC ptcl mass
!     CP_OC      : index of cpt into which emit OC ptcl mass
!     FRACBCEM   : Fraction of BC ems to go into each mode
!     FRACOCEM   : Fraction of OC ems to go into each mode
!
!     Inputted by module UKCA_SETUP_INDICES
!     -------------------------------------
!     Various indices for budget terms in BUD_AER_MAS
!
!     References
!     ----------
!     Bond et al (2004), "A technology-based global inventory of black &
!       organic carbon emissions from combustion",
!       JGR Atmos, 109 (D14): Art. No. D14203.
!
!     Stier et al (2005), "The aerosol-climate model ECHAM5-HAM",
!       ACP, 5,  pp. 1125-1156.
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
      REAL    :: MD(NBOX,NMODES,NCP)
      REAL    :: ND(NBOX,NMODES)
      REAL    :: MDT(NBOX,NMODES)
      REAL    :: DTC
      REAL    :: SM(NBOX)
      REAL    :: AIRD(NBOX)
      REAL    :: SURTP(NBOX)
      REAL    :: EMC(NBOX,4)
      REAL    :: EMCBM(NBOX,2,6)
      REAL    :: BUD_AER_MAS(NBOX,NBUDAER)
!
!     Local variables
      INTEGER :: JL
      INTEGER :: IMODE
      INTEGER :: ICP
      INTEGER :: A
      INTEGER :: I
      INTEGER :: IA
      REAL    :: MMAIR
      REAL    :: EMCVOL
      REAL    :: TOTNUM
      REAL    :: NUM
      REAL    :: FACTOR
      REAL    :: TOTEMC
      REAL    :: EMCBC
      REAL    :: EMCOC
      REAL    :: NEWN
      REAL    :: DELN
      REAL    :: MODE_DIAM(2)
      REAL    :: STDEV(2)
      REAL    :: LSTDEV(2)
      REAL    :: MDNDCP(NCP)
      REAL    :: DELMDNDCP

      IF(ISO2EMS /= 3) THEN
! If ISO2EMS /= 3 then use AEROCOM ACB recommendations for BC/OC size
       MODE_DIAM=(/80.0,30.0/)
       STDEV=(/1.8,1.8/)
      ENDIF
      IF(ISO2EMS == 3) THEN
! If ISO2EMS == 3 then use AEROCOM ACB recommendations for BC/OC size
!                          as modified by Stier et al (2005)
       MODE_DIAM=(/150.0,60.0/)
       STDEV=(/1.59,1.59/)
      ENDIF
!
      MMAIR=AVC*ZBOLTZ/RA
!
      DO JL=1,NBOX
       FACTOR=DTC*AIRD(JL)*MMAIR/(SM(JL)*AVC)   ! converts from
                                                ! per gridbox/s to
                                                ! per cc per timestep

       DO IMODE=1,NMODES
        IF(MODE(IMODE)) THEN

! If some BC/OC is emitted into this mode
         IF((FRACBCEM(IMODE)+FRACOCEM(IMODE)) > 0.0) THEN

          DO A=1,2 ! loop over bio-fuel/fires(1) and fossil-fuel(2)
!
! A=1 accesses EMC(1) and EMC(3) (BC and OC) for bio-fuel & wildfires
! A=2 accesses EMC(2) and EMC(4) (BC and OC) for fossil-fuel

!-----------------------------------------------------------------------
! .. This section updates due to bio-fuel and fossil-fuel emissions

           LSTDEV(A)=LOG(STDEV(A))

! If gridbox is at surf then inject emissions there (bio & fossil fuel)
           IF(SURTP(JL) < 1.5) THEN
!           Calculate natural logs of standard deviations
!
            EMCBC=EMC(JL,A)*FRACBCEM(IMODE)
            EMCOC=EMC(JL,A+2)*FRACOCEM(IMODE)
            TOTEMC=EMCBC+EMCOC

! If amount to emit is significant
            IF (TOTEMC > EMS_EPS) THEN
!
! Calculate total emitted particle volume (nm3 per gridbox per s)
! (assumes carbonaceous ptcls are internal mixture of BC & OC).
             EMCVOL=1E27*(EMCBC/RHOCOMP(CP_BC)+EMCOC/RHOCOMP(CP_OC))
!
! Calculate total emitted particle number (per gridbox per s)
             TOTNUM=EMCVOL/( (PPI/6.0)*((MODE_DIAM(A))**3.0)*           &
                EXP(4.5*LSTDEV(A)*LSTDEV(A)) )
!
             DELN=TOTNUM*FACTOR
             NEWN=ND(JL,IMODE)+DELN

             IF((DELN > DN_EPS).AND.(NEWN > NUM_EPS(IMODE))) THEN

              DO ICP=1,NCP
               IF(COMPONENT(IMODE,ICP)) THEN
                MDNDCP(ICP)=MD(JL,IMODE,ICP)*ND(JL,IMODE)
               ENDIF
              ENDDO

              DO ICP=1,NCP
               IF(COMPONENT(IMODE,ICP)) THEN

! .. Update BC component mass
                IF(ICP == CP_BC) THEN

                 DELMDNDCP=(FACTOR*AVC/MM(ICP))*EMCBC
! .. EMCBC is in kgC/gridbox/s
! .. FACTOR converts from /gridbox/s to /cc/tstep
! .. so EMCBC*FACTOR is in kgC/cc/tstep
! .. so EMCBC*FACTOR/MM(ICP) is in moles of CP_BC/cc/tstep
! .. so DELMDNDCP=EMCBC*FACTOR*AVC/MM(ICP) is in
! ..                                       molecules of CP_BC/cc/tstep

                 MDNDCP(ICP)=MDNDCP(ICP)+DELMDNDCP

                 IF(IMODE == 5) THEN ! if insoluble Aitken mode
                  IF(NMASPRIMBCAITINS > 0)                              &
                   BUD_AER_MAS(JL,NMASPRIMBCAITINS)=                    &
                   BUD_AER_MAS(JL,NMASPRIMBCAITINS)+DELMDNDCP
                 ENDIF
                 IF(IMODE == 2) THEN ! if   soluble Aitken mode
                  IF(NMASPRIMBCAITSOL > 0)                              &
                   BUD_AER_MAS(JL,NMASPRIMBCAITSOL)=                    &
                   BUD_AER_MAS(JL,NMASPRIMBCAITSOL)+DELMDNDCP
                 ENDIF

                ENDIF ! if ICP=CP_BC

! .. Update OC component mass
                IF(ICP == CP_OC) THEN

                 DELMDNDCP=(FACTOR*AVC/MM(ICP))*EMCOC
! .. EMCOC is in kgC/gridbox/s
! .. FACTOR converts from /gridbox/s to /cc/tstep
! .. so EMCOC*FACTOR is in kgC/cc/tstep
! .. so EMCOC*FACTOR/MM(ICP) is in moles of CP_OC/cc/tstep
! .. so DELMDNDCP=EMCOC*FACTOR*AVC/MM(ICP) is in
! ..                                       molecules of CP_OC/cc/tstep

                 MDNDCP(ICP)=MDNDCP(ICP)+DELMDNDCP

                 IF(IMODE == 5) THEN ! if insoluble Aitken mode
                  IF(NMASPRIMOCAITINS > 0)                              &
                    BUD_AER_MAS(JL,NMASPRIMOCAITINS)=                   &
                    BUD_AER_MAS(JL,NMASPRIMOCAITINS)+DELMDNDCP
                 ENDIF
                 IF(IMODE == 2) THEN ! if   soluble Aitken mode
                  IF(NMASPRIMOCAITSOL > 0)                              &
                    BUD_AER_MAS(JL,NMASPRIMOCAITSOL)=                   &
                    BUD_AER_MAS(JL,NMASPRIMOCAITSOL)+DELMDNDCP
                 ENDIF

                ENDIF ! if ICP=CP_OC

               ENDIF ! if COMPONENT(IMODE,ICP)
              ENDDO ! loop over cpts

! .. Update MD,MDT and ND due to BC/OC emissions from ff and bf sources
              MDT(JL,IMODE) = 0.0
              DO ICP=1,NCP
               IF(COMPONENT(IMODE,ICP)) THEN
                MD(JL,IMODE,ICP)=MDNDCP(ICP)/NEWN
                MDT(JL,IMODE)=MDT(JL,IMODE)+MD(JL,IMODE,ICP)
               ENDIF
              ENDDO
              ND(JL,IMODE)=NEWN

             ENDIF ! if DELN>DN_EPS and NEWN>NUM_EPS
            ENDIF ! if TOTEMC>EMS_EPS
           ENDIF ! if at surface (above is for bio-fuel and fossil-fuel)
!
!-----------------------------------------------------------------------
! .. This section updates due to biomass burning BC/OC emissions

           IF(A == 1) THEN

! .. If A=1 then do biomass burning BC/OC

            DO IA=1,6 ! loop over each injection height range

! .. If some biomass burning emissions are injected into this level
             IF((EMCBM(JL,1,IA)+EMCBM(JL,2,IA)) > 0.0) THEN

              EMCBC=EMCBM(JL,1,IA)*FRACBCEM(IMODE) ! in kgC/gridbox/s
              EMCOC=EMCBM(JL,2,IA)*FRACOCEM(IMODE) ! in kgC/gridbox/s
              TOTEMC=EMCBC+EMCOC

! If amount to emit is significant
              IF (TOTEMC > EMS_EPS) THEN

! .. Calculate total emitted particle volume (nm3 per gridbox per s)
! .. (assumes carbonaceous ptcls are internal mixture of BC & OC).
               EMCVOL=1E27*(EMCBC/RHOCOMP(CP_BC)+EMCOC/RHOCOMP(CP_OC))
!
! .. Calculate total emitted particle number (per gridbox per s)
               TOTNUM=EMCVOL/((PPI/6.0)*((MODE_DIAM(A))**3.0)*          &
                      EXP(4.5*LSTDEV(A)*LSTDEV(A))        )
!
               DELN=TOTNUM*FACTOR ! in particles per cc per timestep
               NEWN=ND(JL,IMODE)+DELN ! in particles per cc

               IF((DELN > DN_EPS).AND.(NEWN > NUM_EPS(IMODE))) THEN

                DO ICP=1,NCP
                 IF(COMPONENT(IMODE,ICP)) THEN
                  MDNDCP(ICP)=MD(JL,IMODE,ICP)*ND(JL,IMODE)
                 ENDIF
                ENDDO

                DO ICP=1,NCP
                 IF(COMPONENT(IMODE,ICP)) THEN

! .. Update BC component mass
                  IF(ICP == CP_BC) THEN

                   DELMDNDCP=(FACTOR*AVC/MM(ICP))*EMCBC
! .. EMCBC is in kgC/gridbox/s
! .. FACTOR converts from /gridbox/s to /cc/tstep
! .. so EMCBC*FACTOR is in kgC/cc/tstep
! .. so EMCBC*FACTOR/MM(ICP) is in moles of CP_BC/cc/tstep
! .. so DELMDNDCP=EMCBC*FACTOR*AVC/MM(ICP) is in
! ..                                       molecules of CP_BC/cc/tstep

                   MDNDCP(ICP)=MDNDCP(ICP)+DELMDNDCP

                   IF(IMODE == 5) THEN ! if insoluble Aitken mode
                    IF(NMASPRIMBCAITINS > 0)                            &
                     BUD_AER_MAS(JL,NMASPRIMBCAITINS)=                  &
                     BUD_AER_MAS(JL,NMASPRIMBCAITINS)+DELMDNDCP
                   ENDIF
                   IF(IMODE == 2) THEN ! if   soluble Aitken mode
                    IF(NMASPRIMBCAITSOL > 0)                            &
                     BUD_AER_MAS(JL,NMASPRIMBCAITSOL)=                  &
                     BUD_AER_MAS(JL,NMASPRIMBCAITSOL)+DELMDNDCP
                   ENDIF

                  ENDIF ! if ICP=CP_BC

! .. Update BC component mass
                  IF(ICP == CP_OC) THEN

                   DELMDNDCP=(FACTOR*AVC/MM(ICP))*EMCOC
! .. EMCOC is in kgC/gridbox/s
! .. FACTOR converts from /gridbox/s to /cc/tstep
! .. so EMCOC*FACTOR is in kgC/cc/tstep
! .. so EMCOC*FACTOR/MM(ICP) is in moles of CP_OC/cc/tstep
! .. so DELMDNDCP=EMCOC*FACTOR*AVC/MM(ICP) is in
! ..                                       molecules of CP_OC/cc/tstep

                   MDNDCP(ICP)=MDNDCP(ICP)+DELMDNDCP

                   IF(IMODE == 5) THEN ! if insoluble Aitken mode
                    IF(NMASPRIMOCAITINS > 0)                            &
                      BUD_AER_MAS(JL,NMASPRIMOCAITINS)=                 &
                      BUD_AER_MAS(JL,NMASPRIMOCAITINS)+DELMDNDCP
                   ENDIF
                   IF(IMODE == 2) THEN ! if   soluble Aitken mode
                    IF(NMASPRIMOCAITSOL > 0)                            &
                      BUD_AER_MAS(JL,NMASPRIMOCAITSOL)=                 &
                      BUD_AER_MAS(JL,NMASPRIMOCAITSOL)+DELMDNDCP
                   ENDIF

                  ENDIF ! if ICP=CP_OC

                 ENDIF ! if COMPONENT(IMODE,ICP)
                ENDDO ! loop over components

! .. Update MD,MDT and ND due to BC/OC emissions from wildfire sources
                MDT(JL,IMODE) = 0.0
                DO ICP=1,NCP
                 IF(COMPONENT(IMODE,ICP)) THEN
                  MD(JL,IMODE,ICP)=MDNDCP(ICP)/NEWN
                  MDT(JL,IMODE)=MDT(JL,IMODE)+MD(JL,IMODE,ICP)
                 ENDIF
                ENDDO
                ND(JL,IMODE)=NEWN

               ENDIF ! if DELN>DN_EPS and NEWN>NUM_EPS
              ENDIF ! if TOTEMC>EMS_EPS
             ENDIF ! IF total biomass BC + OC > 0.0
            ENDDO ! loop over each range
           ENDIF ! if A=1

!-----------------------------------------------------------------------

          ENDDO ! loop over A (1=biofuel/biomass, 2=fossilfuel)
         ENDIF ! if mode has BC/OC emitted into it
        ENDIF ! if mode is defined
       ENDDO ! loop over NMODES
      ENDDO ! loop over NBOX
!
      RETURN
      END SUBROUTINE UKCA_PRIM_CAR
#endif
