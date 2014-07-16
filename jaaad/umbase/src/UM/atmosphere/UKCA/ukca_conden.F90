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
!    Calculates condensation of condensable cpt vapours onto pre-existing
!    aerosol particles. Includes switch for using either Fuchs (1964) or
!    modified Fuchs and Sutugin (1971) calculation of CC.
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
      SUBROUTINE UKCA_CONDEN(NBOX,GC,ND,MD,MDT,                         &
       DTZ,DRYDP,WETDP,TSQRT,RHOA,AIRDM3,                               &
       DELTAGC,IFUCHS,AGETERM1,BUD_AER_MAS)
!----------------------------------------------------------------------
!
! Purpose
! -------
! Calculates condensation of condensable cpt vapours onto pre-existing
! aerosol particles. Includes switch for using either Fuchs (1964) or
! modified Fuchs and Sutugin (1971) calculation of CC.
!
! Parameters
! ----------
! SE_SOL : Sticking efficiency for   soluble modes [set to 1.0 as in M7]
! SE_INS : Sticking efficiency for insoluble modes [set to 0.3 as in M7]
!
! Inputs
! ------
! NBOX     : Number of grid boxes
! GC       : Condensable cpt number density (molecules cm-3)
! ND       : Aerosol ptcl no. concentration (ptcls per cc)
! MD       : Component median aerosol mass (molecules per ptcl)
! MDT      : Total median aerosol mass (molecules per ptcl)
! DTZ      : Time Step for nucl/cond competition (s)
! DRYDP    : Avg dry diameter for each aerosol mode (m)
! WETDP    : Avg wet diameter for each aerosol mode (m)
! TSQRT    : Square-root of mid-level temperature (K)
! RHOA     : Air density (kg/m3)
! AIRDM3   : Number density of air (per m3)
! IFUCHS   : Switch for Fuchs (1964) or Fuchs-Sutugin (1971) for CC
! BUD_AER_MAS : Aerosol mass fluxes (molecules/cc/DTC)
!
! Outputs:
! -------
! MD    : Avg aerosol ptcl mass in mode (molecules per ptcl)
! GC    : Condensable cpt number density (molecules per cc)
! AGETERM1: stores mass of soluble material which has condensed onto
!           each of the insoluble modes for use in UKCA_AGEING
!           (molecules per cc)
!
! Local variables:
! ---------------
! MMAIR   : Molar mass of air (kg mol-1)
! DMOL    : Molecular diameter of condensable cpt (m)
! CC      : Conden. coeff. for condensable cpt onto particle (m^3s^-1)
! RP      : Radius of aerosol particle (m)
! SE      : Sticking efficiency (accomodation coeff)
! SE_SOL  : Sticking efficiency (accomodation coeff) for soluble mode
! SE_INS  : Sticking efficiency (accomodation coeff) for insoluble mode
! MMCG    : Molar mass of condensing gas (kg/mole)
! NC      : Product of number conc and condensation coefficient
! SUMNC   : Sum of NC over all modes
! DELTAMS : Mass of condensing gas taken up by this   soluble mode
! DELTAMI : Mass of condensing gas taken up by this insoluble mode
! DELTAM  : Mass of condensing gas taken up by both modes (-->soluble)
!    n.b. DELTAMS,DELTAMI,DELTAM all in molecules per cc.
! MASK1-3 : Logical array to define regions of domain to work on
!
! References
! ----------
! Gong et al, JGR, 108(D1), 4007, doi:10.1029/2001JD002002, 2003.
! Raes et al, J. Aerosol Sci., 23 (7), pp. 759--771, 1992.
! Fuchs & Sutugin, Topics in aerosol research, 1971.
! Fuchs, "Mechanics of aerosols", Pergamon Press, 1964.
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! RA       : Dry air gas constant = 287.05 Jkg^-1 K^-1
! RR       : Universal gas constant = 8.314 Jmol^-1 K^-1
! PPI      : 3.1415927...........
! AVC      : Avogadros constant (mol-1)
! CONC_EPS : Threshold for condensable conc (molecules per cc)
! ZBOLTZ   : Boltzmann's constant (kg m2 s-2 K-1 molec-1)
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES   : Number of possible aerosol modes
! NCP      : Number of possible aerosol components
! MODE     : Defines which modes are set
! COMPONENT: Defines which cpts are allowed in each mode
! CONDENSABLE : Logical variable defining which cpts are condensable
! MODESOL  : Defines whether the mode is soluble or not (=1 or 0)
! MM       : Molar masses of components (kg/mole)
! DIMEN    : Molecular diamters of condensable components (m)
! NUM_EPS  : Value of NEWN below which do not recalculate MD (per cc)
!                                             or carry out process
! CP_SU    : Index of component in which sulfate is stored
! CP_OC    : Index of component in which 1st OC cpt is stored
! CP_SO    : Index of component in which 2nd OC cpt is stored
!
! Inputted by module UKCA_SETUP_INDICES
! -------------------------------------
! NCHEMG   : Number of gas phase tracers for gas phase chemistry scheme
! MM_GAS   : Array of molar masses for gas phase species (kg/mol)
! MH2SO4   : Index of MM_GAS, WTRATC and S0G for H2SO4
! MSEC_ORG : Index of MM_GAS, WTRATC and S0G for SEC_ORG
! Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES

      IMPLICIT NONE
!
! Subroutine interface:
      INTEGER :: NBOX
      INTEGER :: IFUCHS
      REAL    :: ND(NBOX,NMODES)
      REAL    :: MD(NBOX,NMODES,NCP)
      REAL    :: MDT(NBOX,NMODES)
      REAL    :: TSQRT(NBOX)
      REAL    :: RHOA(NBOX)
      REAL    :: AIRDM3(NBOX)
      REAL    :: DTZ
      REAL    :: GC(NBOX,NCHEMG)
      REAL    :: DELTAGC(NBOX,NCHEMG)
      REAL    :: WETDP(NBOX,NMODES)
      REAL    :: DRYDP(NBOX,NMODES)
      REAL    :: AGETERM1(NBOX,3,NCHEMG)
      REAL    :: BUD_AER_MAS(NBOX,NBUDAER)
!
! Local variables
      INTEGER :: ICP
      INTEGER :: IMODE
      INTEGER :: JV
      LOGICAL :: MASK1(NBOX)
      LOGICAL :: MASK2(NBOX)
      LOGICAL :: MASK3(NBOX)
      REAL    :: DMOL
      REAL    :: MMCG
      REAL    :: MMAIR
      REAL    :: CC(NBOX)
      REAL    :: RP(NBOX)
      REAL    :: SUMNC(NBOX)
      REAL    :: NC(NBOX,NMODES)
      REAL    :: DELTAM(NBOX)
      REAL    :: DELTAMS(NBOX)
      REAL    :: DELTAMI(NBOX)
      REAL    :: SE
      REAL, PARAMETER :: SE_SOL=1.0
      REAL, PARAMETER :: SE_INS=0.3
!
      MMAIR=AVC*ZBOLTZ/RA
!
      AGETERM1(:,:,:)=0.0
!
      DO JV=1,NCHEMG
       IF(CONDENSABLE(JV)) THEN
!
!       Set component into which component will condense
        ICP=CONDENSABLE_CHOICE(JV)
!
        DMOL=DIMEN(JV)
        MMCG=MM_GAS(JV)
        DELTAGC(:,JV)=0.0
!
        MASK1(:) = GC(:,JV) > CONC_EPS
!
        SUMNC(:)=0.0
        DO IMODE=1,NMODES
         IF(MODE(IMODE)) THEN
!
          NC(:,IMODE)=0.0
          MASK2(:) = MASK1(:) .AND. ND(:,IMODE) > NUM_EPS(IMODE)
!
          RP(:)=WETDP(:,IMODE)*0.5
!
          IF(MODESOL(IMODE) == 1) SE=SE_SOL
          IF(MODESOL(IMODE) == 0) SE=SE_INS
!
!         Calculate change in condensable cpt conc (molecules cm^-3)
! DEPENDS ON: ukca_cond_coff_v
          CALL UKCA_COND_COFF_V(NBOX,MASK2,RP,TSQRT,AIRDM3,RHOA,        &
                                MMCG,MMAIR,SE,DMOL,IFUCHS,CC)
          WHERE(MASK2(:))
           NC(:,IMODE)=ND(:,IMODE)*CC(:)
           SUMNC(:)=SUMNC(:)+NC(:,IMODE)
          ENDWHERE
!
         ENDIF ! if mode is present
        ENDDO ! Over modes
!
        WHERE(MASK1(:))                                                 &
         DELTAGC(:,JV)=GC(:,JV)*(1.0-EXP(-SUMNC(:)*DTZ))
!
!       Update condensable cpt concentration (molecules cm^-3)
!
        MASK2(:) = MASK1(:) .AND. DELTAGC(:,JV) > CONC_EPS
        WHERE(MASK2(:) .AND. DELTAGC(:,JV) > GC(:,JV))
         DELTAGC(:,JV)=DELTAGC(:,JV)/GC(:,JV) ! make sure no -ves
        ENDWHERE
        WHERE(MASK2(:)) GC(:,JV)=GC(:,JV)-DELTAGC(:,JV)
!
        DO IMODE=1,4 ! loop over sol modes (do cond sol -> ins here too)
         IF(MODE(IMODE)) THEN
!
!         Calculate increase in total & cpt masses in each soluble mode
!
          DELTAMS(:)=0.0
          DELTAMI(:)=0.0
!
          MASK3(:) = MASK2(:) .AND. ND(:,IMODE) > NUM_EPS(IMODE)
!
          IF(IMODE == 1) THEN
!
           IF((ICP == CP_SU).AND.(NMASCONDSUNUCSOL > 0)) THEN
            WHERE(MASK3(:))
!
             DELTAMS(:)=DELTAGC(:,JV)*NC(:,IMODE)/SUMNC(:)
!
             BUD_AER_MAS(:,NMASCONDSUNUCSOL)=                           &
             BUD_AER_MAS(:,NMASCONDSUNUCSOL)+DELTAMS(:)
!
            ENDWHERE
           ENDIF
!
           IF((ICP == CP_OC).AND.(NMASCONDOCNUCSOL > 0)) THEN
            WHERE(MASK3(:))
!
             DELTAMS(:)=DELTAGC(:,JV)*NC(:,IMODE)/SUMNC(:)
!
             BUD_AER_MAS(:,NMASCONDOCNUCSOL)=                           &
             BUD_AER_MAS(:,NMASCONDOCNUCSOL)+DELTAMS(:)
!
            ENDWHERE
           ENDIF
!
           IF((ICP == CP_SO).AND.(NMASCONDSONUCSOL > 0)) THEN
            WHERE(MASK3(:))
!
             DELTAMS(:)=DELTAGC(:,JV)*NC(:,IMODE)/SUMNC(:)
!
             BUD_AER_MAS(:,NMASCONDSONUCSOL)=                           &
             BUD_AER_MAS(:,NMASCONDSONUCSOL)+DELTAMS(:)
!
            ENDWHERE
           ENDIF
!
          ENDIF ! if IMODE=1
!
          IF(IMODE == 2) THEN
!
           IF((ICP == CP_SU).AND.(NMASCONDSUAITSOL > 0)) THEN
            WHERE(MASK3(:))
!
             DELTAMS(:)=DELTAGC(:,JV)*NC(:,IMODE)/SUMNC(:)
!
             BUD_AER_MAS(:,NMASCONDSUAITSOL)=                           &
             BUD_AER_MAS(:,NMASCONDSUAITSOL)+DELTAMS(:)
!
            ENDWHERE
           ENDIF
!
           IF((ICP == CP_SU).AND.(NMASCONDSUAITINS > 0)) THEN
            WHERE(MASK3(:))
!
             DELTAMI(:)=DELTAGC(:,JV)*NC(:,IMODE+3)/SUMNC(:)
!
             BUD_AER_MAS(:,NMASCONDSUAITINS)=                           &
             BUD_AER_MAS(:,NMASCONDSUAITINS)+DELTAMI(:)
!
             AGETERM1(:,IMODE-1,JV)=DELTAMI(:)
!
            ENDWHERE
           ENDIF
!
           IF((ICP == CP_OC).AND.(NMASCONDOCAITSOL > 0)) THEN
            WHERE(MASK3(:))
!
             DELTAMS(:)=DELTAGC(:,JV)*NC(:,IMODE)/SUMNC(:)
!
             BUD_AER_MAS(:,NMASCONDOCAITSOL)=                           &
             BUD_AER_MAS(:,NMASCONDOCAITSOL)+DELTAMS(:)
!
            ENDWHERE
           ENDIF
!
           IF((ICP == CP_OC).AND.(NMASCONDOCAITINS > 0)) THEN
            WHERE(MASK3(:))
!
             DELTAMI(:)=DELTAGC(:,JV)*NC(:,IMODE+3)/SUMNC(:)
!
             BUD_AER_MAS(:,NMASCONDOCAITINS)=                           &
             BUD_AER_MAS(:,NMASCONDOCAITINS)+DELTAMI(:)
!
             AGETERM1(:,IMODE-1,JV)=DELTAMI(:)
!
            ENDWHERE
           ENDIF
!
           IF((ICP == CP_SO).AND.(NMASCONDSOAITSOL > 0)) THEN
            WHERE(MASK3(:))
!
             DELTAMS(:)=DELTAGC(:,JV)*NC(:,IMODE)/SUMNC(:)
!
             BUD_AER_MAS(:,NMASCONDSOAITSOL)=                           &
             BUD_AER_MAS(:,NMASCONDSOAITSOL)+DELTAMS(:)
!
            ENDWHERE
           ENDIF
!
           IF((ICP == CP_SO).AND.(NMASCONDSOAITINS > 0)) THEN
            WHERE(MASK3(:))
!
             DELTAMI(:)=DELTAGC(:,JV)*NC(:,IMODE+3)/SUMNC(:)
!
             BUD_AER_MAS(:,NMASCONDSOAITINS)=                           &
             BUD_AER_MAS(:,NMASCONDSOAITINS)+DELTAMI(:)
!
             AGETERM1(:,IMODE-1,JV)=DELTAMI(:)
!
            ENDWHERE
           ENDIF
!
          ENDIF ! if IMODE=2
!
          IF(IMODE == 3) THEN
!
           IF((ICP == CP_SU).AND.(NMASCONDSUACCSOL > 0)) THEN
            WHERE(MASK3(:))
!
             DELTAMS(:)=DELTAGC(:,JV)*NC(:,IMODE)/SUMNC(:)
!
             BUD_AER_MAS(:,NMASCONDSUACCSOL)=                           &
             BUD_AER_MAS(:,NMASCONDSUACCSOL)+DELTAMS(:)
!
            ENDWHERE
           ENDIF
!
!!           IF((ICP == CP_SU).AND.(NMASCONDSUACCINS > 0)) THEN
!!            WHERE(MASK3(:))
!!!
!!             DELTAMI(:)=DELTAGC(:,JV)*NC(:,IMODE+3)/SUMNC(:)
!!!
!!             BUD_AER_MAS(:,NMASCONDSUACCINS)=                          &
!!     &       BUD_AER_MAS(:,NMASCONDSUACCINS)+DELTAMI(:)
!!!
!!             AGETERM1(:,IMODE-1,JV)=DELTAMI(:)
!!!
!!            ENDWHERE
!!           ENDIF
!
           IF((ICP == CP_OC).AND.(NMASCONDOCACCSOL > 0)) THEN
            WHERE(MASK3(:))
!
             DELTAMS(:)=DELTAGC(:,JV)*NC(:,IMODE)/SUMNC(:)
!
             BUD_AER_MAS(:,NMASCONDOCACCSOL)=                           &
             BUD_AER_MAS(:,NMASCONDOCACCSOL)+DELTAMS(:)
!
            ENDWHERE
           ENDIF
!
!!           IF((ICP == CP_OC).AND.(NMASCONDOCACCINS > 0)) THEN
!!            WHERE(MASK3(:))
!!!
!!             DELTAMI(:)=DELTAGC(:,JV)*NC(:,IMODE+3)/SUMNC(:)
!!!
!!             BUD_AER_MAS(:,NMASCONDOCACCINS)=                         &
!!     &       BUD_AER_MAS(:,NMASCONDOCACCINS)+DELTAMI(:)
!!!
!!             AGETERM1(:,IMODE-1,JV)=DELTAMI(:)
!!!
!!            ENDWHERE
!!           ENDIF
!
           IF((ICP == CP_SO).AND.(NMASCONDSOACCSOL > 0)) THEN
            WHERE(MASK3(:))
!
             DELTAMS(:)=DELTAGC(:,JV)*NC(:,IMODE)/SUMNC(:)
!
             BUD_AER_MAS(:,NMASCONDSOACCSOL)=                           &
             BUD_AER_MAS(:,NMASCONDSOACCSOL)+DELTAMS(:)
!
            ENDWHERE
           ENDIF
!
!!           IF((ICP == CP_SO).AND.(NMASCONDSOACCINS > 0)) THEN
!!            WHERE(MASK3(:))
!!!
!!             DELTAMI(:)=DELTAGC(:,JV)*NC(:,IMODE+3)/SUMNC(:)
!!!
!!             BUD_AER_MAS(:,NMASCONDSOACCINS)=                         &
!!     &       BUD_AER_MAS(:,NMASCONDSOACCINS)+DELTAMI(:)
!!!
!!             AGETERM1(:,IMODE-1,JV)=DELTAMI(:)
!!!
!!            ENDWHERE
!!           ENDIF
!
          ENDIF ! if IMODE=3
!
          IF(IMODE == 4) THEN
!
           IF((ICP == CP_SU).AND.(NMASCONDSUCORSOL > 0)) THEN
            WHERE(MASK3(:))
!
             DELTAMS(:)=DELTAGC(:,JV)*NC(:,IMODE)/SUMNC(:)
!
             BUD_AER_MAS(:,NMASCONDSUCORSOL)=                           &
             BUD_AER_MAS(:,NMASCONDSUCORSOL)+DELTAMS(:)
!
            ENDWHERE
           ENDIF
!
!!           IF((ICP == CP_SU).AND.(NMASCONDSUCORINS > 0)) THEN
!!            WHERE(MASK3(:))
!!!
!!             DELTAMI(:)=DELTAGC(:,JV)*NC(:,IMODE+3)/SUMNC(:)
!!!
!!             BUD_AER_MAS(:,NMASCONDSUCORINS)=                         &
!!     &       BUD_AER_MAS(:,NMASCONDSUCORINS)+DELTAMI(:)
!!!
!!             AGETERM1(:,IMODE-1,JV)=DELTAMI(:)
!!!
!!            ENDWHERE
!!           ENDIF
!
           IF((ICP == CP_OC).AND.(NMASCONDOCCORSOL > 0)) THEN
            WHERE(MASK3(:))
!
             DELTAMS(:)=DELTAGC(:,JV)*NC(:,IMODE)/SUMNC(:)
!
             BUD_AER_MAS(:,NMASCONDOCCORSOL)=                           &
             BUD_AER_MAS(:,NMASCONDOCCORSOL)+DELTAMS(:)
!
            ENDWHERE
           ENDIF
!
!!           IF((ICP == CP_OC).AND.(NMASCONDOCCORINS > 0)) THEN
!!            WHERE(MASK3(:))
!!!
!!             DELTAMI(:)=DELTAGC(:,JV)*NC(:,IMODE+3)/SUMNC(:)
!!!
!!             BUD_AER_MAS(:,NMASCONDOCCORINS)=                         &
!!     &       BUD_AER_MAS(:,NMASCONDOCCORINS)+DELTAMI(:)
!!!
!!             AGETERM1(:,IMODE-1,JV)=DELTAMI(:)
!!!
!!            ENDWHERE
!!           ENDIF
!
           IF((ICP == CP_SO).AND.(NMASCONDSOCORSOL > 0)) THEN
            WHERE(MASK3(:))
!
             DELTAMS(:)=DELTAGC(:,JV)*NC(:,IMODE)/SUMNC(:)
!
             BUD_AER_MAS(:,NMASCONDSOCORSOL)=                           &
             BUD_AER_MAS(:,NMASCONDSOCORSOL)+DELTAMS(:)
!
            ENDWHERE
           ENDIF
!
!!           IF((ICP == CP_SO).AND.(NMASCONDSOCORINS > 0)) THEN
!!            WHERE(MASK3(:))
!!!
!!             DELTAMI(:)=DELTAGC(:,JV)*NC(:,IMODE+3)/SUMNC(:)
!!!
!!             BUD_AER_MAS(:,NMASCONDSOCORINS)=                         &
!!     &       BUD_AER_MAS(:,NMASCONDSOCORINS)+DELTAMI(:)
!!!
!!             AGETERM1(:,IMODE-1,JV)=DELTAMI(:)
!!!
!!            ENDWHERE
!!           ENDIF
!
          ENDIF ! if IMODE=4
!
!         All mass condensed onto sol. & ins. goes to sol. mode
          DELTAM(:)=DELTAMS(:)+DELTAMI(:)
!
          WHERE(MASK3(:))
           MD(:,IMODE,ICP)=                                             &
            (MD(:,IMODE,ICP)*ND(:,IMODE)+DELTAM(:))/ND(:,IMODE)
           MDT(:,IMODE)=                                                &
            (MDT(:,IMODE)*ND(:,IMODE)+DELTAM(:))/ND(:,IMODE)
          ENDWHERE
!
         ENDIF ! if mode is present
!
        ENDDO ! IMODE=1,4 (soluble modes)

       ENDIF ! IF CONDENSABLE(JV)
!
      ENDDO ! DO JV=1,NCHEMG
!
      RETURN
      END SUBROUTINE UKCA_CONDEN
#endif
