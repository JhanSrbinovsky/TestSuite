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
!   Carries out ageing of particles in insoluble mode.
!   Calculates the number of particles which will be coated
!   by the total soluble material and transfers number and
!   mass accordingly.
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
      SUBROUTINE UKCA_AGEING(NBOX,ND,MD,MDT,                            &
       AGETERM1,AGETERM2,WETDP,VERBOSE,BUD_AER_MAS)
!----------------------------------------------------------
!
! Purpose
! ------
! Carries out ageing of particles in insoluble mode.
! Calculates the number of particles which will be coated
! by the total soluble material and transfers number and
! mass accordingly.
! Number of ptcls given by assuming 1 particle ages with
! a 1 molecule layer thickness of particle of geometric
! mean size of mode.  i.e. 1 ptcl requires r_g^2/r_molec^2
! = 10^4,10^6,10^8 molecules for Aitken, accum, coarse modes
!
! Inputs
! ------
! NBOX      : Number of grid boxes
! ND        : Aerosol ptcl no. concentration (ptcls per cc)
! MD        : Component median aerosol mass (molecules per ptcl)
! MDT       : Total median aerosol mass (molecules per ptcl)
! AGETERM1  : Depletion rate of each component (molecules cpt/cc/DTZ)
!             from condensation onto the 3 insoluble modes.
! AGETERM2  : Rate of accomodation of material to each insoluble mode
!             as a result of coagulation with smaller soluble modes
!             (in molecules cpt /cm3/DTZ)
! WETDP     : Wet diameter corresponding to mean sized particle
! VERBOSE   : Switch for whether to do various test print statements
!
! Outputs
! -------
! ND        : Updated number concentration [each mode] (ptcls/cc)
! MD        : Updated avg cpt   mass conc. [each mode] (molecules/ptcl)
! MDT       : Updated avg total mass conc. [each mode] (molecules/ptcl)
! BUD_AER_MAS : Aerosol mass budgets
!
! Local variables
! ---------------
! AGE1PTCL  : Mass of soluble material needed to age 1 ptcl (molecules)
! NAGED     : # of insoluble ptcls aged to soluble mode [total] (/cc)
! NAGED_JV  : # of insoluble ptcls aged to soluble mode [by gas](/cc)
! TOTAGE_JV : Ageing flux to cpt by particular gas (molecules gas/cc)
! TOTAGE1   : Total ageing flux to cpt by conden.  (molecules gas/cc)
! TOTAGE2   : Total ageing flux to cpt by coaguln. (molecules gas/cc)
! TOTAGE    : Total ageing flux to cpt (coag+cond) (molecules gas/cc)
! NDINSNEW  : # in ins mode after reduction due to ageing (/cc)
! NDSOLNEW  : # in corresp. sol mode after reduction due to ageing (/cc)
! F_MM      : Ratio of molar masses of condensable gas to aerosol cpt
! CP_COAG_ADDED : Switch for whether added on ageing flux by
!                 coagulation to that cpt already
!                 (loop over jv --- need to make sure only count once)
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! CONC_EPS  : Likewise, threshold for soluble material (molecules per cc)
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES    : Number of possible aerosol modes
! NCP       : Number of possible aerosol components
! DIMEN     : Molecular diameters of condensable components
! MODE      : Which modes are being carried
! SOLUBLE   : Logical variable defining which cpts are soluble
! COMPONENT : Which components are in each of modes
! MMID      : Mid-point masses for initial radius grid
! MFRAC_0   : Initial mass fraction to set when no particles.
! MM        : Molar masses of components (kg per mole)
! NUM_EPS   : Value of NEWN below which do not recalculate MD (per cc)
!             or carry out process
! CP_SU     : Component where sulfate is stored
! CP_BC     : Component in which black carbon is stored
! CP_OC     : Component in which primary organic carbon is stored
! CP_SO     : Component where secondary organic carbon is stored
!
! Inputted by module UKCA_SETUP_INDICES
! -------------------------------------
! NCHEMG    : Number of gas phase tracers for gas phase chemistry scheme
! MM_GAS    : Array of molar masses for gas phase species (kg/mol)
! MSEC_ORG  : Index of MM_GAS, WTRATC and S0G for SEC_ORG
! Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES
!
      IMPLICIT NONE
!
! .. Arguments
      INTEGER :: NBOX
      INTEGER :: VERBOSE
      REAL    :: ND(NBOX,NMODES)
      REAL    :: MD(NBOX,NMODES,NCP)
      REAL    :: MDT(NBOX,NMODES)
      REAL    :: AGETERM1(NBOX,3,NCHEMG)
      REAL    :: AGETERM2(NBOX,4,3,NCP)
      REAL    :: WETDP(NBOX,NMODES)
      REAL    :: BUD_AER_MAS(NBOX,NBUDAER)
!
! .. Local variables
      INTEGER :: JL,JV,IMODE,JMODE,ICP
      INTEGER :: CP_COAG_ADDED(NCP)
      REAL    :: TOTAGE(NCP)
      REAL    :: TOTAGE_JV
      REAL    :: TOTAGE1(NCP)
      REAL    :: TOTAGE2(NCP)
      REAL    :: AGE1PTCL
      REAL    :: NAGED
      REAL    :: NAGED_JV(NCHEMG)
      REAL    :: NDINSNEW
      REAL    :: NDSOLNEW
      REAL    :: F_MM
!
      DO IMODE=5,7 ! loop over insoluble modes
       IF(MODE(IMODE)) THEN
        DO JL=1,NBOX
         IF(ND(JL,IMODE) > NUM_EPS(IMODE)) THEN
          CP_COAG_ADDED(:)=0
          TOTAGE(:)=0.0
          TOTAGE1(:)=0.0 ! from condensation
          TOTAGE2(:)=0.0 ! from coagulation
          NAGED=0.0
          DO JV=1,NCHEMG
           TOTAGE_JV=0.0
           IF(CONDENSABLE(JV)) THEN

! .. Below add on amount of soluble material taken up by insoluble modes
! .. as result of condensation of this gas onto insoluble modes
            ICP=CONDENSABLE_CHOICE(JV)
            F_MM=MM_GAS(JV)/MM(ICP)
            TOTAGE_JV=AGETERM1(JL,IMODE-4,JV)/F_MM
            TOTAGE (ICP)=TOTAGE (ICP)+TOTAGE_JV
            TOTAGE1(ICP)=TOTAGE1(ICP)+TOTAGE_JV
! .. 100% of condensable material taken up by aerosol assumed soluble.
! .. AGETERM1 is from condensation onto insoluble modes
! .. Divide by F_MM because AGETERM1 is in "molecules" of cpt, whereas
! .. TOTAGE needs to be in molecules of condensible gas (H2SO4/SEC_ORG)
! .. since AGE1PTCL is in these units.

! .. Below add on amount of soluble material taken up by insoluble modes
! .. as result of coagulation with soluble modes
!
! .. IF(CP_COAG_ADDED(ICP) == 0) checks if already added on this
! .. aerosol component's coagulated material (make sure don't dblecount)
!
            IF(CP_COAG_ADDED(ICP) == 0) THEN
             DO JMODE=1,4
              TOTAGE_JV=TOTAGE_JV+                                      &
                 AGETERM2(JL,JMODE,IMODE-4,ICP)/F_MM
! .. add on amount coagulated
              TOTAGE (ICP)=TOTAGE (ICP)+                                &
                 AGETERM2(JL,JMODE,IMODE-4,ICP)/F_MM
              TOTAGE2(ICP)=TOTAGE2(ICP)+                                &
                 AGETERM2(JL,JMODE,IMODE-4,ICP)/F_MM
             ENDDO
             CP_COAG_ADDED(ICP)=1
            ENDIF
! .. 100% of condensable material taken up by aerosol assumed soluble.
! .. AGETERM2 is from coagn of soluble modes with larger insoluble modes
! .. Divide by F_MM because AGETERM2 is in "molecules" of cpt, whereas
! .. TOTAGE needs to be in molecules of condensible gas (H2SO4/SEC_ORG)
! .. since AGE1PTCL is in these units.
!
            AGE1PTCL=WETDP(JL,IMODE)*WETDP(JL,IMODE)/DIMEN(JV)/DIMEN(JV)
! above is number of molecules of condensable gas to age 1 particle
            NAGED_JV(JV)=TOTAGE_JV/AGE1PTCL
            NAGED=NAGED+NAGED_JV(JV)

           ENDIF ! if gas phase species is condensable
          ENDDO ! loop over gas phase species
!
          IF(NAGED > NUM_EPS(IMODE)) THEN
!
           IF(NAGED > ND(JL,IMODE)) THEN
            NAGED=ND(JL,IMODE) ! limit so no -ves
            TOTAGE(:)=TOTAGE(:)*ND(JL,IMODE)/NAGED
! above reduces ageing if limited by insoluble particles
           ENDIF
!
           NDINSNEW=ND(JL,IMODE)-NAGED
! set new insoluble mode no. (but don't update ND yet)
           NDSOLNEW=ND(JL,IMODE-3)+NAGED
! set new   soluble mode no. (but don't update ND yet)
!

           IF(SUM(TOTAGE) > 0.0) THEN
            DO ICP=1,NCP
! below transfers aged cpt masses from ins modes (doesn't include SU)
             IF(COMPONENT(IMODE-3,ICP)) THEN
! above if statement checks whether cpt is in corresponding soluble mode
              IF(IMODE == 5) THEN
               IF((ICP == CP_SU).AND.(NMASAGEDSUINTR52 > 0))            &
             BUD_AER_MAS(JL,NMASAGEDSUINTR52)=                          &
             BUD_AER_MAS(JL,NMASAGEDSUINTR52)+TOTAGE(ICP)*F_MM
               IF((ICP == CP_OC).AND.(NMASAGEDOCINTR52 > 0))            &
             BUD_AER_MAS(JL,NMASAGEDOCINTR52)=                          &
             BUD_AER_MAS(JL,NMASAGEDOCINTR52)+TOTAGE(ICP)*F_MM
               IF((ICP == CP_SO).AND.(NMASAGEDSOINTR52 > 0))            &
             BUD_AER_MAS(JL,NMASAGEDSOINTR52)=                          &
             BUD_AER_MAS(JL,NMASAGEDSOINTR52)+TOTAGE(ICP)*F_MM
!! above accounts for transfer of mass of condensed/coagulated material
!! multiply above by F_MM 'cos TOTAGE is in molecules of gas phase
!! species (H2SO4/SEC_ORG) whereas needs to be in "molecules" of
!! aerosol component (CP_SU/CP_OC/CP_SU)
               IF((ICP == CP_BC).AND.(NMASAGEDBCINTR52 > 0))            &
             BUD_AER_MAS(JL,NMASAGEDBCINTR52)=                          &
             BUD_AER_MAS(JL,NMASAGEDBCINTR52)+NAGED*MD(JL,IMODE,ICP)
               IF((ICP == CP_OC).AND.(NMASAGEDOCINTR52 > 0))            &
             BUD_AER_MAS(JL,NMASAGEDOCINTR52)=                          &
             BUD_AER_MAS(JL,NMASAGEDOCINTR52)+NAGED*MD(JL,IMODE,ICP)
!! .. above 2 lines account for transfer of mass due to aged aerosol
              ENDIF ! if mode is Aitken-insoluble

!! .. below calculates new cpt total masses in soluble modes due to
!! .. transfer of mass of particles which were in Ait-ins
!! .. (n.b. insoluble avg. masses unchanged)

              IF(COMPONENT(IMODE,ICP)) THEN ! if in insoluble mode
! if sol. mode cpt is in insoluble mode then update corresponding
! soluble mode cpt mass due to transfer from ins. mode & coag of sol.
               MD(JL,IMODE-3,ICP)=(ND(JL,IMODE-3)*MD(JL,IMODE-3,ICP)    &
                 +NAGED*MD(JL,IMODE,ICP)+TOTAGE2(ICP)*F_MM)/NDSOLNEW
!! .. only update from coag (TOTAGE2) 'cos MD already updated by conden
!! .. in UKCA_CONDEN (but used TOTAGE to age ptcls in insoluble modes
!! .. and used TOTAGE in budget terms)
              ELSE ! if not in insoluble mode
! if sol. mode component is not in insoluble mode then just update
! corresponding soluble mode cpt mass due to coag (& change in number)
! then still need to update avg mass to reflect reduction in # of ptcls
               MD(JL,IMODE-3,ICP)=(ND(JL,IMODE-3)*MD(JL,IMODE-3,ICP)    &
                                        +TOTAGE2(ICP)*F_MM)/NDSOLNEW
!! .. only update from coag (TOTAGE2) 'cos MD already updated by conden
!! .. in UKCA_CONDEN (but used TOTAGE to age ptcls in insoluble modes
!! .. and used TOTAGE in budget terms)
              ENDIF
             ENDIF ! if component is in corresponding soluble mode
            ENDDO ! loop over components
           ENDIF ! if total amount of accomodated material > 0
!
!! update ND and MDT for ins mode
           ND(JL,IMODE  )=NDINSNEW
           MDT(JL,IMODE)=0.0
           DO ICP=1,NCP
            IF(COMPONENT(IMODE,ICP)) THEN
             MDT(JL,IMODE)=MDT(JL,IMODE)+MD(JL,IMODE,ICP)
            ENDIF
           ENDDO ! end loop over cpts
!
!! update ND and MDT for sol mode
           ND(JL,IMODE-3)=NDSOLNEW
           MDT(JL,IMODE-3)=0.0
           DO ICP=1,NCP
            IF(COMPONENT(IMODE-3,ICP)) THEN
             MDT(JL,IMODE-3)=MDT(JL,IMODE-3)+MD(JL,IMODE-3,ICP)
            ENDIF
           ENDDO ! end loop over cpts
!
          ENDIF ! if number of aged particles > NUM_EPS(IMODE)
         ENDIF ! if some particles in insoluble modes (ND>epsilon)
        ENDDO ! end loop over boxes
       ENDIF ! if insoluble mode is present (IMODE=5,7)
      ENDDO ! Loop IMODE=5,7
!
! .. Note to accompany subroutine:
! .. ----------------------------
! ..
! .. Because the soluble mass was never transferred from the soluble
! .. modes to the insoluble modes in the coagulation routine, it
! .. needs to be transferred here to the appropriate soluble mode.
! .. i.e. need to move transferred mass from original soluble modes to
! .. soluble mode which corresponds to larger insoluble mode
! ..
! .. Since all condensable gas phase species mass taken up onto the
! .. insoluble modes is automatically added to the soluble modes in
! .. the condensation routine there is no need to do any transfer of
! .. condensed soluble mass here.
! ..
! .. However, if there were no particles initially in the soluble
! .. mode, then the mass will not have been added at the
! .. condensation stage. So in that case, after the number has
! .. been added, the mass from condensation does need to be added.
!
      RETURN
      END SUBROUTINE UKCA_AGEING
#endif
