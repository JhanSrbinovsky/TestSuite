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
!    Calculate (wet) volume corresponding to mid-pt particles in each
!    aerosol mode.
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
      SUBROUTINE UKCA_VOLUME_MODE(NBOX,ND,MD,MDT,                       &
       RH,PMID,T,WVOL,WETDP,RHOPAR,IWVOLMETHOD,                         &
       DVOL,DRYDP,MDWAT,VERBOSE)
!----------------------------------------------------------------------
!
! Calculate (wet) volume corresponding to mid-pt particles in each mode.
! Two methods can be followed (as specified by IWVOLMETHOD:
!
! 1) Quick method --- soluble particle mass assumed to behave as SO4
!                     and wet volume calculated from Kohler equation
!                     ln(rh)=A/Dp - B/Dp^3 --- in sub-saturated
!                     conditions can be approximated as
!                     ln(rh)=-B/Dp^3, [Dp is wet diameter of particle]
!                     so  wetvol=(pi/6)Dp^3=-(pi/6)B/ln(rh)
!                     Seinfeld & Pandis have B=3.44e13.nu.(ms/Ms)
!                     which we set as B=3.44e13.nu.(MDSOL/AVC) where
!                     MDSOL is all soluble mass (per particle),
!                     # of dissociating ions, nu=3 and then
!                     wetvol=-8.974e-29*MDSOL/ln(rh)z
!                     Also add on dry volume of insoluble
!                     molecules to give overall (wet) ptcl volume.
!
! 2) Better method -- calculate the water-content and density of the
!                     aerosol using ZSR and water activity coeffs for
!                     H+,SO42-,Na+,Cl- from Jacobsen (pg 610).
!                     Complete dissociation of solute ions is assumed
!                     to give the electrolyte concentrations and
!                     associated water content associated with each.
!                     The dry volume of insoluble molecules is then
!                     added to give overall (wet) ptcl volume.
!                     Note that OC is assumed water-insoluble in the
!                     insoluble mode but is assumed to have aged
!                     chemically in the aerosol to become hygroscopic
!                     once transferred to the soluble distribution.
!                     To respresent this, in the ZSR calculation,
!                     the concentration of SO4 ions is incremented
!                     by FHYG_AOM*MD(:,:,CP_OC)/F_AO -- i.e. the aged
!                     OC is assumed to take up water at a fraction
!                     FHYG_AOM (set at 0.65) of SO4.
!
! In each case sphericity is assumed to give the ptcl wet radius.
!
! Purpose
! -------
! Calculate avg wet volume WVOL & wet diameter WETDP for each aerosol
! particle size mode from the relative humidity and avg. number of
! molecules per particle (MD) of each cpt in each mode.
!
! (N.b. Local rel. hum. values are corrected to lie between 0.1 & 0.9)
!
! Inputs
! ------
! NBOX     : Number of grid boxes
! ND       : Aerosol ptcl no. concentration (ptcls per cc)
! MD       : Component median aerosol mass (molecules per ptcl)
! MDT      : Total median aerosol mass (molecules per ptcl)
! RH       : Relative humidity (corrected to lie within range 0.1-0.9)
! PMID     : Centre level air pressure (Pa)
! T        : Centre level air temperature (K)
! IWVOLMETHOD: Chosen method for calculation of wet volume
! DVOL     : Dry volume of particle (m^3)
! DRYDP    : Dry diameter of particle (m)
! VERBOSE  : Switch for level of verbosity
!
! Outputs
! ------
! WVOL     : Avg wet volume of size mode (m3)
! WETDP    : Avg wet diameter of size mode (m)
! MDWAT    : Molecular concentration of water (molecules per particle)
! RHOPAR   : Particle density [incl. H2O & insoluble cpts] (kgm^-3)
!
! Local Variables
! ---------------
! CORRH    : Locally corrected RH (to lie between 0.1 and 0.9)
! FHYG_AOM : Hygroscopicity of aged organic (fraction relative to SO4)
! F_AO     : No. of "moles" of POM in 1 mole of aged organic species
! MM_AGE_ORG: Molar mass of aged organic species (kg/mol)
! MM_POM   : Molar mass of particulate organic matter [for CP_OC,CP_SO]
! IONS     : Logical indicating presence of ions
! CL       : Ion concentrations (moles/cm3 of air)
! MASK     : Mask where in domain to calculate values
! MDSOL    : Mass per particle (total over soluble cpts) (mlcls/ptcl)
! RHOSOL   : Density of particle solution [excl. insoluble cpts] (kg/m3)
! B        : Factor in solute term in Kohler equation =
!              BCONST*(no. of ions)*(solute mass)/(solute molar mass)
! DENOM    : Temporary variable calculating denominator of expression
! DENOM2   : Temporary variable calculating denominator of expression
! RHOTMP   : Temporary variable in calculation of particle density
! RHOTMP2  : Temporary variable in calculation of particle density
! WVOLTMP  : Temporary variable in calculation of wet volume
! WVOL_INS : Contribution to wet volume from insoluble components (m^3)
! WVOL_SOL : Contribution to wet volume from soluble components (m^3)
! WC       : Water content for aerosol (moles/cm3 of air)
! WDPCUB   : Cube of particle wet diameter (m^3)
! SIXOVRPIX: (6.0/pi)*{ 1.0/EXP((9/2)*LOG^2(SIGMA_G)) }
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! PPI      : 3.1415927....
! AVC      : Avogadro's constant (molecules per mole)
! CONC_EPS : Value of soluble material mass conc. below which
!            assume no soluble mass
! RHOW     : Density of water (=1000.0 kgm^-3)
! MMW      : Molecular mass of water (=0.018 kg/mol)
! BCONST   : Value of constant in B term of Kohler equation
! NU_H2SO4 : Number of dissociated ions for H2SO4
! CONVERT  : Conversion from micron^3 to m^3
! RHOSUL   : Mass density of a pure H2SO4 aerosol (kg per m^3)
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES   : Number of possible aerosol modes
! NCP      : Number of possible aerosol components
! MODE     : Logical variable denoting where mode is defined
! COMPONENT: Logical variable denoting where cpt is defined
! SOLUBLE  : Logical variable defining which cpts are soluble
! MM       : Molar masses of components (kg per mole)
! NO_IONS  : Number of dissociating ions in solute (H2SO4=3,NaCl=2)
! MODESOL  : Defines which modes are soluble (integer)
! RHOCOMP  : Densities (dry) of each component (kg/m^3)
! MMID     : Avg mass of mode when rmed_g=exp(0.5*(lnr0+lnr1)) (ptcl^-1)
! X        : EXP((9/2)*LOG^2(SIGMA_G))
! NUM_EPS  : Value of NEWN below which don't recalculate MD (per cc)
!                                            or carry out process
! CP_SU    : Index of component containing SO4    component
! CP_OC    : Index of component containing 1st OC component
! CP_CL    : Index of component containing NaCl   component
! CP_SO    : Index of component containing 2nd OC component
!
! Inputted by module UKCA_SETUP_INDICES
! -------------------------------------
! NCHEMG   : Number of gas phase tracers for gas phase chemistry scheme
! MM_GAS   : Array of molar masses for gas phase species (kg/mol)
! MSEC_ORG : Index of MM_GAS, WTRATC and S0G for SEC_ORG
!
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES

      IMPLICIT NONE
!
!     Subroutine interface
      INTEGER :: NBOX
      INTEGER :: IWVOLMETHOD
      INTEGER :: VERBOSE
      REAL    :: ND(NBOX,NMODES)
      REAL    :: MD(NBOX,NMODES,NCP)
      REAL    :: MDT(NBOX,NMODES)
      REAL    :: RH(NBOX)
      REAL    :: PMID(NBOX)
      REAL    :: T(NBOX)
      REAL    :: WVOL(NBOX,NMODES)
      REAL    :: WETDP(NBOX,NMODES)
      REAL    :: RHOPAR(NBOX,NMODES)
      REAL    :: DVOL(NBOX,NMODES)
      REAL    :: DRYDP(NBOX,NMODES)
      REAL    :: MDWAT(NBOX,NMODES)
!
!     Local variables
      INTEGER :: I
      INTEGER :: JL
      INTEGER :: IMODE
      INTEGER :: ICP
      LOGICAL :: MASK(NBOX)
      LOGICAL :: IONS(NBOX,-NANION:NCATION) !ION PRESENCE SWITCHES
      REAL    :: CORRH(NBOX)
      REAL    :: B(NBOX)
      REAL    :: MDSOL(NBOX)
      REAL    :: RHOSOL(NBOX)
      REAL    :: WC(NBOX)
      REAL    :: WVOLTMP(NBOX)
      REAL    :: WVOL_SOL(NBOX)
      REAL    :: WVOL_INS(NBOX)
      REAL    :: RHOTMP(NBOX)
      REAL    :: RHOTMP2(NBOX)
      REAL    :: DENOM(NBOX)
      REAL    :: DENOM2(NBOX)
      REAL    :: CBRT
      REAL    :: WDPCUB(NBOX)
      REAL    :: SIXOVRPIX(NMODES)
      REAL    :: F_AO
      REAL    :: CL(NBOX,-NANION:NCATION) !ION CONCS (MOL/CC OF AIR)
      REAL, PARAMETER :: FHYG_AOM=0.65
      REAL, PARAMETER :: MM_AGE_ORG=0.150
      REAL, PARAMETER :: MM_POM=0.0168
!
      IF((IWVOLMETHOD /= 1).AND.(IWVOLMETHOD /= 2)) THEN
! DEPENDS ON: ereport
       CALL EREPORT('UKCA_VOLUME_MODE',imode,'IWVOLMETHOD not 1 or 2')
      ENDIF
!
      SIXOVRPIX(:)=6.0/(PPI*X(:))
!
!     Correct relative humidities to lie within the range of 10-90%
      CORRH(:)=RH(:)
      WHERE (CORRH(:) > 0.9) CORRH(:)=0.9
      WHERE (CORRH(:) < 0.1) CORRH(:)=0.1
      DO IMODE=1,NMODES
       IF(MODE(IMODE)) THEN
        IF(MODESOL(IMODE) == 1) THEN

         IF(IWVOLMETHOD == 1) THEN
! .. assume all soluble components are H2SO4 & uses Kohler theory
          MASK(:) = ND(:,IMODE) > NUM_EPS(IMODE)
          MDSOL(:)=0.0
          WVOL_INS(:)=0.0
          DO ICP=1,NCP
           IF(COMPONENT(IMODE,ICP)) THEN
            IF(SOLUBLE(ICP)) THEN
             WHERE(MASK(:)) MDSOL(:)=MDSOL(:)+MD(:,IMODE,ICP)
            ELSE
! .. store contribution to volume from insoluble components (m^3)
             WHERE(MASK(:))                                             &
              WVOL_INS(:)=WVOL_INS(:)+                                  &
                MD(:,IMODE,ICP)*MM(ICP)/RHOCOMP(ICP)/AVC
            ENDIF
           ENDIF
          ENDDO
          WHERE(MASK(:))
!
!     Kohler equation is lnS=A/Dp - B/Dp^3
!     In sub-saturated environment, B/Dp^3 > A/Dp
!     Reasonable assumption B/Dp^3 >> A/Dp for all but smallest ptcls
!     Then,   Dp=(-B/ln(rh))^(1/3)
!       or  volp=-(pi/6)*B/ln(rh)
!
!     (CONVERT converts from micron^3 to m^3 in B term, see S&P pg 787).
!
           B(:)=BCONST*NU_H2SO4*MDSOL(:)/AVC
           WVOL(:,IMODE)=-(PPI/6.0)*B(:)/LOG(CORRH(:))*CONVERT          &
                         +WVOL_INS(:)
           MDWAT(:,IMODE)=(WVOL(:,IMODE)-DVOL(:,IMODE))*RHOW*AVC/MMW
           RHOTMP(:)=RHOW*MDWAT(:,IMODE)*MMW
           DENOM(:)=MDWAT(:,IMODE)*MMW
          ENDWHERE
          DO ICP=1,NCP
           IF(COMPONENT(IMODE,ICP)) THEN
            WHERE(MASK(:))
             RHOTMP(:)=RHOTMP(:)+RHOCOMP(ICP)*MD(:,IMODE,ICP)           &
                        *MM(ICP)
             DENOM(:)=DENOM(:)+MD(:,IMODE,ICP)*MM(ICP)
            ENDWHERE
           ENDIF
          ENDDO
! .. RHOPAR total particle density [including insoluble cpts] (kgm^-3)
          WHERE(MASK(:))
           RHOPAR(:,IMODE)=RHOTMP(:)/DENOM(:)
           WVOL(:,IMODE)=MAX(WVOL(:,IMODE),DVOL(:,IMODE))
          ENDWHERE

          WHERE(.NOT.MASK(:))
           WVOL(:,IMODE)=DVOL(:,IMODE)
           MDWAT(:,IMODE)=0.0
           RHOPAR(:,IMODE)=RHOSUL
          ENDWHERE

         ENDIF ! if IWVOLMETHOD = 1

         F_AO=MM_AGE_ORG/MM_POM

         IF(IWVOLMETHOD == 2) THEN

! .. Use composition information to calculate water uptake by
! .. each component according to ZSR using water activity data from
! .. Jacobsen page 610 (Table B.10) for binary electrolyte molalities

!***********************************************************************
!**  Liquid Phase Species:
!**
!**  **Cations**           **Anions**             **Neutrals**
!**  1: H                  -1: HSO4               0: H2O
!**  2: NH4                -2: SO4
!**  3: Na                 -3: NO3
!**                        -4: Cl
!***********************************************************************

          MASK(:) = ND(:,IMODE) > NUM_EPS(IMODE)

          DO I=-NANION,NCATION
           CL(:,I)=0.0 ! set all concentrations to zero initially
          ENDDO

          IF(COMPONENT(IMODE,CP_SU)) THEN ! assume all H2SO4 --> SO4
           CL(:,-2)=MD(:,IMODE,CP_SU)/AVC   ! [SO4] in moles/cc (air)
          ENDIF

          IF(MSEC_ORG > 0) THEN ! if secondary organic species included
           IF(COMPONENT(IMODE,CP_SO)) THEN
            CL(:,-2)=CL(:,-2)+(FHYG_AOM/AVC)*(MD(:,IMODE,CP_SO)/F_AO)
! .. Increment concentration of SO4 ions to represent the
! .. presence of hygroscopic aged organic aerosol mass in CP_SO.
! .. Assume it has uptake behaviour at fraction FHYG_AOM of SO4.
! .. Need to divide by F_AO because MD of CP_SO is in "moles" of POM
! .. whereas CL needs to be in moles of aged organic (MM=0.15kg/mol).
           ENDIF
          ENDIF

          IF(COMPONENT(IMODE,CP_OC)) THEN
           CL(:,-2)=CL(:,-2)+(FHYG_AOM/AVC)*(MD(:,IMODE,CP_OC)/F_AO)
! .. Increment concentration of SO4 ions to represent the
! .. presence of hygroscopic aged organic aerosol mass in CP_OC.
! .. Assume it has uptake behaviour at fraction FHYG_AOM of SO4.
! .. Need to divide by F_AO because MD of CP_OC is in "moles" of POM
! .. whereas CL needs to be in moles of aged organic (MM=0.15kg/mol).
!
! .. This effectively says that by the time the primary carbonaceous
! .. aerosol has been microphysically aged to the soluble mode, the
! .. organic component has been chemically aged to become hygroscopic.
! .. (it is assumed to be water-insoluble in the insoluble mode)
!
          ENDIF

          IF(COMPONENT(IMODE,CP_CL)) THEN ! assume complete dissociation
           CL(:,3)=MD(:,IMODE,CP_CL)/AVC  ! [Na] in moles per cc (air)
           CL(:,-4)=MD(:,IMODE,CP_CL)/AVC ! [Cl] in moles per cc (air)
          ENDIF

! SET H+ FOR CHARGE BALANCE  -- CL(1) is [H] in moles per cc(air)
          CL(:,1)=MAX((2.0*CL(:,-2)+CL(:,-1)+CL(:,-3)+CL(:,-4)          &
                       -CL(:,2)-CL(:,3)),0.0)

          DO I=-NANION,NCATION
           IONS(:,I)=(CL(:,I) > 0.)
          ENDDO

! DEPENDS ON: ukca_water_content_v
          CALL UKCA_WATER_CONTENT_V(NBOX,MASK,CL,CORRH,IONS,WC)

          WHERE(MASK(:))
           MDWAT(:,IMODE)=WC(:)*AVC
! .. calculate solution density (avg over each cpt mass contribution)
           RHOTMP(:)=RHOW*MDWAT(:,IMODE)*MMW
           RHOTMP2(:)=RHOTMP(:)
           DENOM(:)=MDWAT(:,IMODE)*MMW
           DENOM2(:)=DENOM(:)
          ENDWHERE

          DO ICP=1,NCP
           IF(COMPONENT(IMODE,ICP)) THEN
            IF(SOLUBLE(ICP)) THEN
             WHERE(MASK(:))
              RHOTMP(:)=RHOTMP(:)+RHOCOMP(ICP)*MD(:,IMODE,ICP)*MM(ICP)
              DENOM(:)=DENOM(:)+MD(:,IMODE,ICP)*MM(ICP)
             ENDWHERE
            ENDIF
            WHERE(MASK(:))
             RHOTMP2(:)=RHOTMP2(:)+RHOCOMP(ICP)*MD(:,IMODE,ICP)*MM(ICP)
             DENOM2(:)=DENOM2(:)+MD(:,IMODE,ICP)*MM(ICP)
            ENDWHERE
           ENDIF
          ENDDO

! .. RHOSOL is density of ptcl solution [exclud. insoluble cpts] (kgm/3)
!
! .. where some soluble material
          WHERE(MASK(:) .AND. DENOM(:) > CONC_EPS)
           RHOSOL(:)=RHOTMP(:)/DENOM(:)
          ELSEWHERE
           RHOSOL(:)=RHOSUL
          END WHERE

! .. RHOPAR total particle density [incl H2O & insoluble cpts] (kgm^-3)
          WHERE(MASK(:) .AND. DENOM2(:) > 1E-20)
            RHOPAR(:,IMODE)=RHOTMP2(:)/DENOM2(:)
          ELSEWHERE
            RHOPAR(:,IMODE)=RHOSUL
          ENDWHERE
! .. calculate ptcl wet volume (including that from insoluble cpts)
          WHERE(MASK(:))
           WVOLTMP(:)=MDWAT(:,IMODE)*MMW
           WVOL_INS(:)=0.0
          ENDWHERE
          DO ICP=1,NCP
           IF(COMPONENT(IMODE,ICP)) THEN
            IF(SOLUBLE(ICP)) THEN
             WHERE(MASK(:))                                             &
              WVOLTMP(:)=WVOLTMP(:)+MD(:,IMODE,ICP)*MM(ICP)
            ELSE
! .. WVOL_INS contribution to volume from insoluble components (m^3)
             WHERE(MASK(:))                                             &
              WVOL_INS(:)=WVOL_INS(:)+                                  &
                 MD(:,IMODE,ICP)*MM(ICP)/RHOCOMP(ICP)/AVC
            ENDIF
           ENDIF
          ENDDO

! .. WVOL_SOL is contribution to volume from soluble components (m^3)

! .. where some soluble material
          WHERE(MASK(:) .AND. DENOM(:) > CONC_EPS)
           WVOL_SOL(:)=WVOLTMP(:)/RHOSOL(:)/AVC
          ELSEWHERE
           WVOL_SOL(:)=0.0
          ENDWHERE
          WHERE(MASK(:)) WVOL(:,IMODE)=WVOL_SOL(:)+WVOL_INS(:)

          WHERE(.NOT.MASK(:))
           WVOL(:,IMODE)=DVOL(:,IMODE)
           MDWAT(:,IMODE)=0.0
           RHOPAR(:,IMODE)=RHOSUL
          ENDWHERE

         ENDIF ! if IWVOLMETHOD=2
        ELSE  ! not soluble
         WVOL(:,IMODE)=DVOL(:,IMODE)
         MDWAT(:,IMODE)=0.0
         RHOPAR(:,IMODE)=RHOSUL
        ENDIF ! if mode is soluble
       ENDIF ! if(mode(imode))
      ENDDO
!
      DO IMODE=1,NMODES
       IF(MODE(IMODE)) THEN
        WDPCUB(:)=SIXOVRPIX(IMODE)*WVOL(:,IMODE)
        WETDP(:,IMODE)=CBRT(WDPCUB(:))
!!        WETDP(:,IMODE)=WDPCUB(:)**(1.0/3.0)
       ENDIF  ! if(mode(imode))
      ENDDO ! loop over modes
!
      RETURN
      END SUBROUTINE UKCA_VOLUME_MODE
#endif
