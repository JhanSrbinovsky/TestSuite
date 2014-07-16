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
!    Calculate emissions of sea salt aerosol as sectional representation
!    of Gong-Monahan and add to soluble accum and coarse modes,
!    changing ND, MDT and MD accordingly (only sea surface boxes).
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
      SUBROUTINE UKCA_PRIM_SS(NBOX,ND,MDT,MD,SURTP,LAND_FRAC,           &
            SURF,SEAICE,DTC,SM,US10M,ZNOT,AIRD,                         &
            BUD_AER_MAS,VERBOSE)
!----------------------------------------------------------------------
!
! Calculate emissions of sea salt aerosol as sectional representation
! of Gong-Monahan and add to soluble accum and coarse modes.
! changing ND, MDT and MD accordingly (only sea surface boxes)
!
! Purpose
! -------
! Use flux parametrisation (Smith,93; Smith,98; Monahan, 86;
! Gong-Monahan, 03) combined with
! wind stress to calculate flux of sea salt aerosol.
!
! Inputs
! ------
! NBOX      : Number of grid boxes
! ND        : Aerosol ptcl number density for mode (cm^-3)
! MDT       : Avg tot mass of aerosol ptcl in mode (particle^-1)
! MD        : Avg cpt mass of aerosol ptcl in mode (particle^-1)
! SURTP     : Surface type [0/1=at sea/land-surf,2/3=above-sea/land]
! LAND_FRAC : Fraction of horizontal gridbox area covered by land
! SURF      : Surface area of box (horizontal) (m^2)
! SEAICE    : Fraction of horizontal gridbox area containing seaice
! DTC       : Chemical time step (s)
! SM        : Grid box mass of air (kg)
! US10M     : Scalar wind at 10m (ms-1)     !CEJ used instead of USTAR
! ZNOT      : Surface roughness length (m)
! AIRD      : Number density of air (per cc)
!
! Outputs
! -------
! ND,MDT,MD   : Updated ptcl no. conc., total avg mass and cpt avg mass
! BUD_AER_MAS : Updated mass budget terms
!
! Local Variables
! ---------------
! MMAIR     : Molar mass of dry air (kg/mole)
! S98FLUX   : Smith sea salt aerosol ems rate [dF/dr] (m-2 s-1 um-1)
! M86FLUX   : Gong/Monahan ssalt aerosol ems rate [dF/dr] (m-2 s-1 um-1)
! COMBFLUX  : Combined Smith-Monahan [dF/dr] (m-2 s-1 um-1)
! DFDR      : Chosen (from above 3) [dF/dr] (m-2 s-1 um-1)
! FLUX      : Sea salt aerosol emissions rate [(dF/dr)*deltar] (m-2 s-1)
! BOXFLUX   : Grid box sea salt aerosol flux [(dF/dr)*deltar*SURF] (s-1)
! DELN      : Change in aerosol bin no. dens. due to ssalt ems (cm^-3)
! BET       : Parameter in Monahan formulation of ssalt aerosol ems
! AGONG     : Parameter in Gong extension of Monahan flux formulation
! MODEMT    : Tmpry varble -- which mode to emit bin-seasalt into
! CUTOFF    : Smallest dry radius that gong-monahan scheme is applicable
!             use 35nm cut-off for r at rh=80 (~17.5 nm cutoff for dryr)
! MSMALL    : Lower bin edge for smallest sea-salt emission bin
! MLARGE    : Upper bin edge for largest  sea-salt emission bin
! NBINS     : Number of bins for sea-salt distribution
! DELTAR    : Difference in dry radii edges for each bin (microns)
! MBLO/MBMID/MBHI : Dry mass of lower/mid/upper particle bin edge (m)
! DRLO/DRMID/DRHI : Dry radius of lower/mid/upper particle bin edge (m)
! MODEMT    : Index of mode to emit this bin-resolved ems flux into
! NEWN      : Updated particle number conc. after ssalt ems (/cm3)
! DUM       : Dummy variable in calculation of DRLO,DRHI,DELTAR etc.
!
! References
! ----------
! Smith 1993, QJRMS, 119, 809-824
! Smith 1998, JAS, 29, S189-190
! Monahan 1986, Oceanic Whitecaps, 167-174
! Gong, 2003, Global Biogeochemical Cycles, 17(4), 8-1,6.
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! AVC       : Avogadro's constant
! VKARMN    : Von Karman's constant = 0.4
! PPI       : 3.1415927
! RA        : Dry air gas constant = 287.05 Jkg^-1 K^-1
! ZBOLTZ    : Stefan-Boltzmann constant (kg m2 s-2 K-1 molec-1)
! DN_EPS    : Value of DELN below which do not carry out process
! MM_DA     : Molar mass of dry air (kg/mol)
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES    : Number of possible aerosol modes
! NCP       : Number of possible aerosol components
! MODE      : Logical variable defining which modes are set.
! COMPONENT : Logical variable defining which cpt are in which modes
! MM        : Molar mass of each of the aerosol components (kg/mol)
! RHOCOMP   : Mass density of each of the aerosol components (kgm^-3)
! DDPLIM0   : Lower limit for dry diameter in mode (m)
! DDPLIM1   : Upper limit for dry diameter in mode (m)
! NUM_EPS   : Value of NEWN below which do not recalculate MD
!                                              or carry out process
! CP_CL     : index of cpt in which sea-salt mass is stored
!
! Inputted by module UKCA_SETUP_INDICES
! -------------------------------------
! Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES

      IMPLICIT NONE
!
!     Subroutine interface
      INTEGER :: NBOX
      INTEGER :: VERBOSE
      REAL    :: ND(NBOX,NMODES)
      REAL    :: MDT(NBOX,NMODES)
      REAL    :: MD(NBOX,NMODES,NCP)
      REAL    :: SURTP(NBOX)
      REAL    :: DTC
      REAL    :: SM(NBOX)
      REAL    :: LAND_FRAC(NBOX)
      REAL    :: SURF(NBOX)
      REAL    :: SEAICE(NBOX)
      REAL    :: US10M(NBOX)
      REAL    :: ZNOT(NBOX)
      REAL    :: AIRD(NBOX)
      REAL    :: BUD_AER_MAS(NBOX,NBUDAER)
      LOGICAL :: LOGICLO
      LOGICAL :: LOGICHI
!
!     Local variables
      INTEGER :: IMODE
      INTEGER :: ICP
      INTEGER :: JL
      INTEGER :: JV
      INTEGER :: MODEMT
      INTEGER, PARAMETER :: NBINS=20
      REAL    :: MBMID(NBINS)
      REAL    :: MBLO(NBINS)
      REAL    :: MBHI(NBINS)
      REAL    :: DRMID(NBINS)
      REAL    :: DELTAR(NBINS)
      REAL    :: DRHI(NBINS)
      REAL    :: DRLO(NBINS)
      REAL    :: FLUX
      REAL    :: BOXFLUX
      REAL    :: DELN
      REAL    :: MMAIR
      REAL    :: BET
      REAL    :: AGONG
      REAL    :: S98FLUX
      REAL    :: M86FLUX
      REAL    :: COMBFLUX
      REAL    :: DFDR
      REAL    :: DUM
      REAL    :: NEWN
      REAL, PARAMETER :: MSMALL=1.0e3
      REAL, PARAMETER :: MLARGE=4.0e14
! .. use 35 nm cut-off for r at rh=80 (approx 17.5 nm cutoff for dryr)
      REAL, PARAMETER :: CUTOFF=0.0175e-6
!
      MMAIR=AVC*ZBOLTZ/RA
!
! DEPENDS ON: ukca_ingridg
      CALL UKCA_INGRIDG(NBINS,MSMALL,MLARGE,MBMID,MBLO,MBHI)

      DUM=3.0*MM(CP_CL)/(AVC*RHOCOMP(CP_CL)*4.0*PPI)
      DO JV=1,NBINS
        DELTAR(JV)=((DUM*MBHI(JV))**(1.0/3.0)-                          &
                    (DUM*MBLO(JV))**(1.0/3.0))*1.0E6 ! in microns
        DRMID(JV)=(DUM*MBMID(JV))**(1.0/3.0) ! in m
        DRLO(JV) =(DUM*MBLO(JV ))**(1.0/3.0) ! in m
        DRHI(JV) =(DUM*MBHI(JV ))**(1.0/3.0) ! in m
      ENDDO
!
      DO JL=1,NBOX
! .. Only emit sea-salt for ocean surface boxes not covered by sea-ice
        IF((SURTP(JL) < 0.5).AND.(SEAICE(JL) < 1.0)) THEN

         DO JV=1,NBINS ! loop over sea-salt emission bins
          IF (DRHI(JV) > CUTOFF) THEN ! if bin size is > cutoff for flux

           MODEMT=0
           IF(DRMID(JV) < 0.5e-6) MODEMT=3 ! soluble accum.
           IF(DRMID(JV) >= 0.5e-6) MODEMT=4 ! soluble coarse

           IF(MODEMT > 0) THEN
!           CUTOFF is lowest dry radius for the Gong-Monahan scheme

!           Smith et al. (1998)
            S98FLUX=1.4367* 0.2*(US10M(JL)**3.5)                        &
                *EXP(-1.5*(LOG(2.0*DRMID(JV)/2.088E-6 ))**2) +          &
                 1.4367 * 6.8E-3*(US10M(JL)**3  )                       &
               *EXP(-1.0*(LOG(2.0*DRMID(JV)/20.88E-6))**2)
!
!           Gong (2003) ---- extension from Monahan et al. (1986)
            BET=(0.433-LOG10(DRMID(JV)*2.0E6))/0.433
! .. checked that this should be log10
            AGONG=4.7*(1.0+30.0*(DRMID(JV)*2.0E6))**(-0.017*            &
                   (DRMID(JV)*2.0E6)**(-1.44))
!           Factor of 2.0 to convert dF/dr(80) to dF/dr(dry)
            M86FLUX=2.0*1.373*(US10M(JL)**3.41)*                        &
                   (DRMID(JV)*2.0E6)**(-AGONG)*                         &
                   (1.0+0.057*((2.0E6*DRMID(JV))**(3.45)))*             &
                   10.0**(1.607*EXP(-BET*BET))
!
!           If using a combination of Smith and Monahan, then
!           the recommendation is to use the larger of the
!           two fluxes at any radius
            IF (S98FLUX > M86FLUX) THEN
              COMBFLUX=S98FLUX
            ELSE
              COMBFLUX=M86FLUX
            ENDIF

            DFDR=M86FLUX ! particles m-2 s-1 um-1

!           Calculate ssalt ems flux in particles m-2 s-1
            IF (DRLO(JV) > CUTOFF) THEN
              FLUX=DFDR*DELTAR(JV)
            ELSE
              FLUX=DFDR*(DRHI(JV)*1.0E6-CUTOFF*1.0E6)
            ENDIF

!           Calculate grid box ssalt aerosol source (particles s-1)
            BOXFLUX=FLUX*SURF(JL)*(1.-SEAICE(JL))*(1.0-LAND_FRAC(JL))
!
!           Calculate change in grid box aerosol number density
            DELN=BOXFLUX*DTC*AIRD(JL)*MMAIR/(SM(JL)*AVC)
            IF((MODEMT == 3).AND.(NMASPRIMSSACCSOL > 0))                &
                  BUD_AER_MAS(JL,NMASPRIMSSACCSOL)=                     &
                  BUD_AER_MAS(JL,NMASPRIMSSACCSOL)+DELN*MBMID(JV)
            IF((MODEMT == 4).AND.(NMASPRIMSSCORSOL > 0))                &
                  BUD_AER_MAS(JL,NMASPRIMSSCORSOL)=                     &
                  BUD_AER_MAS(JL,NMASPRIMSSCORSOL)+DELN*MBMID(JV)

! .. Calculate new mode number conc. due to ssalt ptcls in this ems bin
            NEWN=ND(JL,MODEMT)+DELN

! .. Update all cpt mass concentrations due to ssalt in this ems bin
            IF((DELN > DN_EPS).AND.(NEWN > NUM_EPS(MODEMT))) THEN
             DO ICP=1,NCP

! .. Check whether component icp is in mode MODEMT
              IF(COMPONENT(MODEMT,ICP)) THEN

               IF(ICP == CP_CL) MD(JL,MODEMT,ICP)=                      &
                    (MD(JL,MODEMT,ICP)*ND(JL,MODEMT)                    &
                    +DELN*MBMID(JV))/NEWN
               IF(ICP /= CP_CL) MD(JL,MODEMT,ICP)=                      &
                    MD(JL,MODEMT,ICP)*ND(JL,MODEMT)/NEWN

              ENDIF ! if COMPONENT(MODEMT,ICP)

             ENDDO ! loop over NCP

! .. Update particle number concentration
             ND(JL,MODEMT)=NEWN

! .. Update total (over all cpts) mass per particle MDT
             MDT(JL,MODEMT) = 0.0
             DO ICP=1,NCP
              IF(COMPONENT(MODEMT,ICP)) MDT(JL,MODEMT)=                 &
                    MDT(JL,MODEMT) + MD(JL,MODEMT,ICP)
             ENDDO

            ENDIF ! if DELN>DN_EPS and NEWN>NUM_EPS

           ENDIF ! IF MODEMT > 0
          ENDIF ! IF DRHI(JV) > CUTOFF
         ENDDO
        ENDIF ! IF AT OCEAN SURFACE AND NOT ICE
      ENDDO

      RETURN
      END SUBROUTINE UKCA_PRIM_SS
#endif
