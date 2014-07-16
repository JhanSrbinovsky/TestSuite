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
!    Calculates aerosol dry deposition and sedimentation.
!    Based on the parameterisation of Zhang et al (2001) which
!    uses the method in the model of Slinn (1982).
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
      SUBROUTINE UKCA_DDEPAER_INCL_SEDI(NBOX,ND,MD,MDT,RHOPAR,ZNOT,     &
       DTC,WETDP,USTR,PMID,PUPPER,PLOWER,T,SURTP,SEAICE,                &
       RHOA,MFPA,DVISC,BUD_AER_MAS,JLABOVE,SEDI_ON,SM)

!----------------------------------------------------------------------
!
! Purpose
! -------
! Calculates aerosol dry deposition and sedimentation.
! Based on the parameterisation of Zhang et al (2001) which
! uses the method in the model of Slinn (1982).
!
! Sedimentation is done using a simple explicit discretization
! which should be adequate for this process-split method.
!
! Evaluate deposition velocity in lowest level as:
!
! V_dep = V_g + 1/(A_r + S_r)
!
! where V_dep is the deposition velocity
!       V_g   is the gravitational velocity = rho_p*Dp^2*g*CF/(18*DVISC)
!       CF    is the Cunningham slip correction
!       DVISC is the dynamic viscosity
!       A_r   is the aerodynamic resitance
!       S_r   is the surface resistance
!       Dp    is the particle diameter
!       rho_p is the particle density
!       g     is the gravitational acceleration
!
! Evaluate S_r=1/{ 3 * ustar * (EB + EIM + EIN) }
!
! following parameterization by Zhang et al (2001) where
!
! EB,EIM,EIN are collection efficiencies for Brownian diffusion,
! impaction and interception respectively.
!
! EB = Sc^-YR where Sc is the particle Schmidt number = nu/D
!                                where nu = kinematic viscosity of air
!                                      D =  particle diffusion coeff.
!
!  and YR is surface-dependent constant, values as in Table 3 (Zhang01)
!         0.50 over water          (Land use category 13-14)
!         0.56 over forest         (Land use category  1- 5)
!         0.54 over grass and ice  (Land use category  6,12)
!
! EIM = { St/(ALPHA+St) }^2
!
!    where St is the Stokes number = V_g * ustar^2 / DVISC  (z0<1mm)
!                                  = V_g * ustar   / (g*CR) (z0>1mm)
!
!                                    [smooth & rough flow regimes]
!
!      and ALPHA,CR are surface-dependent constant, values as Table 3:
!         ALPHA=100.0, CR=0.0 over water [only divide by CR for veg]
!         ALPHA= 50.0, CR=0.0 over ice   [only divide by CR for veg]
!         ALPHA=  1.0, CR=0.005 over grass
!         ALPHA=  1.2, CR=0.002 over forest
!
! EIN = 0.5*Dp/CR
!
! Evaluates drydep & sedimentation for number & mass using 0th & 3rd
! order moment specific coefficients for modal aerosol as in Appendix 4
! Binkowski & Shankar (1995) JGR, vol 100, no D12, pp. 26,191--26,209.
!
! Note --- in this routine, sedimentation is included at all levels.
!
! Calls functions UKCA_DCOFF_PAR_AV_K and UKCA_VGRAV_AV_K
!
! Inputs :
! ------
! NBOX      : Number of grid boxes
! ND        : Initial no. concentration of aerosol mode (ptcls/cc)
! MD        : Avg cpt mass of aerosol ptcl in size mode (particle^-1)
! MDT       : Avg tot mass of aerosol ptcl in size mode (particle^-1)
! RHOPAR    : Particle density [incl. H2O & insoluble cpts] (kgm^-3)
! ZNOT      : Roughness length (m)
! DTC       : Chemical timestep (s)
! WETDP     : Wet diameter for ptcl with dry diameter DRYDP (m)
! USTR      : Friction velocity(ms-1)
! PMID      : Centre level pressure (Pa)
! PUPPER    : Pressure at box upper interface (Pa)
! PLOWER    : Pressure at box lower interface (Pa)
! T         : Centre level temperature (K)
! SURTP     : Surface type [0=sea-surf,1=land-surf,2=above-surf]
! SEAICE    : Fraction of horizontal gridbox area containing seaice
! RHOA      : Air density (kg/m3)
! MFPA      : Mean free path of air (m)
! DVISC     : Dynamic viscosity of air (kg m-1 s-1)
! JLABOVE   : Index of box directly above this grid box
! SEDI_ON   : Switch for whether aerosol sedimentation is on/off
! SM        : Grid box mass of air (kg)
!
! Calls subroutine GETROUGH to read roughness length
!
! Outputs
! -------
! Updated particle number density ND (/cm3)
! Updated particle avg mass MD (molecules/particle)
! BUD_AER_MAS: Aerosol mass budgets (mlcls/cc/tstep)
!
! Local Variables
! ---------------
! PS_AV_0    : 0th moment avg particle Schmidt Number
! PS_AV_3    : 3rd moment avg particle Schmidt Number
! KVISC      : Kinematic viscosity of air (m2 s-1)
! VGRAV_AV_0 : 0th moment avg grav. settling vel. (m s^-1)
! VGRAV_AV_3 : 3rd moment avg grav. settling vel. (m s^-1)
! VDEP_AV_0  : 0th moment avg deposition velocity (m s^-1)
! VDEP_AV_3  : 3rd moment avg deposition velocity (m s^-1)
! DCOEF_AV_0 : 0th moment avg particle diffusion coefficient(m2 s-1)
! DCOEF_AV_3 : 3rd moment avg particle diffusion coefficient(m2 s-1)
! SN_AV_0    : 0th moment avg Stokes number
! SN_AV_3    : 3rd moment avg Stokes number
! SR_AV_0    : 0th moment avg surface resistance
! SR_AV_3    : 3rd moment avg surface resistance
! EB_AV_0    : 0th moment avg collection eff. for Brownian diffusion
! EB_AV_3    : 3rd moment avg collection eff. for Brownian diffusion
! EIM_AV_0   : 0th moment avg collection eff. for impaction
! EIM_AV_3   : 3rd moment avg collection eff. for impaction
! EIN        : Collection eff. for interception
! AR         : Aerodynamic resistance
! MTOT       : Total aerosol mass conc [all cpts] (molecules/cm3)
! MCPTOT     : Total aersool mass conc [1 cpt] (molecules/cm3)
! NEWN       : Updated number concentration (/cm3)
! DZ         : Ht difference between box vertical interfaces (m)
! DZMID      : Ht difference between box lower interface & mid-level (m)
! SSIGMA     : Geometric standard deviation of mode
! CR,Y,ALPHA: aerosol deposition coefficients
!     [vary with land category & input via DATA statements]
! CR        : Characteristic radius of collectors (m)
! Y         : Parameter for calculating Brownian diffusion
! ALPHA     : Parameter for calculating EIM
! MASK      : Logical to define regions of domain to work on.
!
! Functions:
! ---------
! UKCA_DCOFF_PAR_AV_K
! UKCA_VGRAV_AV_K
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! GG        : Gravitational acceleration = 9.80665 ms^-2
! VKARMN    : Von Karman's constant = 0.4
! PPI       : 3.1415927
! RA        : Dry air gas constant = 287.05 Jkg^-1 K^-1
! ZBOLTZ    : Boltzman Constant (kg m2 s-2 K-1 molec-1)
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES    : Number of possible aerosol modes
! NCP       : Number of possible aerosol components
! MODE      : Defines which modes are set
! COMPONENT : Defines which cpts are allowed in each mode
! SIGMAG    : Geometric standard deviation of mode
! MM        : Molar masses of components (kg/mole)
! RHOCOMP   : Densities (dry) of each component (kg/m^3)
! NUM_EPS   : Value of ND_0 below which do not recalculate MD (per cc)
!                                              or carry out process
! CP_SU     : Index of component in which H2SO4 is stored
! CP_BC     : Index of component in which BC is stored
! CP_OC     : Index of component in which 1st OC cpt is stored
! CP_CL     : Index of component in which NaCl is stored
! CP_SO     : Index of component in which 2nd OC cpt is stored
!
! Inputted by module UKCA_SETUP_INDICES
! -------------------------------------
! Various indices for budget terms in BUD_AER_MAS
!
! References
! ----------
! Slinn, Atmos. En., 1982, 16, 1785-1794
! Zhang et al, Atmos. En., 2001, 35, 549-560
!
!----------------------------------------------------------------------
      USE UKCA_CONSTANTS,   ONLY: GG, VKARMN, PPI, RA, ZBOLTZ
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES
      IMPLICIT NONE
!
!     Subroutine interface
      INTEGER :: NBOX
      INTEGER :: SEDI_ON
      INTEGER :: JLABOVE(NBOX)
      REAL    :: ND(NBOX,NMODES)
      REAL    :: MD(NBOX,NMODES,NCP)
      REAL    :: MDT(NBOX,NMODES)
      REAL    :: RHOPAR(NBOX,NMODES)
      REAL    :: ZNOT(NBOX)
      REAL    :: DTC
      REAL    :: WETDP(NBOX,NMODES)
      REAL    :: USTR(NBOX)
      REAL    :: PMID(NBOX)
      REAL    :: PUPPER(NBOX)
      REAL    :: PLOWER(NBOX)
      REAL    :: T(NBOX)
      REAL    :: SURTP(NBOX)
      REAL    :: SEAICE(NBOX)
      REAL    :: RHOA(NBOX)
      REAL    :: MFPA(NBOX)
      REAL    :: DVISC(NBOX)
      REAL    :: BUD_AER_MAS(NBOX,NBUDAER)
      REAL    :: SM(NBOX)
!
!     Local Variables
      INTEGER :: ICAT(NBOX)
      INTEGER :: JL
      INTEGER :: IMODE
      INTEGER :: ICP
      INTEGER :: JL_UP
      LOGICAL :: MASK(NBOX)
      LOGICAL :: LOGIC1
      LOGICAL :: LOGIC2
      REAL    :: PS_AV_0
      REAL    :: PS_AV_3
      REAL    :: KVISC
      REAL    :: VGRAV_AV_0
      REAL    :: VGRAV_AV_3
      REAL    :: VGRAV_AV_0_UP
      REAL    :: VGRAV_AV_3_UP
      REAL    :: DCOEF_AV_0
      REAL    :: DCOEF_AV_3
      REAL    :: EB_AV_0
      REAL    :: EB_AV_3
      REAL    :: EIM_AV_0
      REAL    :: EIM_AV_3
      REAL    :: EIN
      REAL    :: SN_AV_0
      REAL    :: SN_AV_3
      REAL    :: AR
      REAL    :: SR_AV_0
      REAL    :: SR_AV_3
      REAL    :: VDEP_AV_0
      REAL    :: VDEP_AV_3
      REAL    :: UKCA_DCOFF_PAR_AV_K
      REAL    :: UKCA_VGRAV_AV_K
      REAL    :: DZMID
      REAL    :: DZ
      REAL    :: DZ_UP
      REAL    :: ND0(NBOX,NMODES)
      REAL    :: MD0(NBOX,NMODES,NCP)
      REAL    :: NDNEW
      REAL    :: TERMIN_1
      REAL    :: TERMIN_2
      REAL    :: TERMIN_N
      REAL    :: TERMOUT_N
      REAL    :: TERMIN_M(NCP)
      REAL    :: TERMOUT_M(NCP)
      REAL    :: DELNSEDI
      REAL    :: DELMSEDI(NCP)
      REAL    :: DELMDDEP
      REAL    :: DELNDDEPTOT
      REAL    :: DELMDDEPTOT
      REAL    :: VGRAV_LIM
      REAL    :: VGRAV_LIM_UP
      REAL, PARAMETER :: YR(5) = (/  0.5, 0.56, 0.54,0.0,0.54/)
      REAL, PARAMETER :: CR(5) = (/ 0.00,5.E-3,2.E-3,0.0,0.00/)
      REAL, PARAMETER :: ALPHA(5) = (/100.0,  1.0,  1.2,0.0,50.0/)
!
! Find out what category (water,forest,grass,desert)
! based on roughness length znot (desert not used at present)
! This should be updated in later version to read land type
! category directly. Desert not represented here.
!
! water/sea - z0<0.001m
      MASK(:)=(ZNOT(:) < 1.0E-3) ! water/sea
      WHERE(MASK(:)) ICAT(:)=1

! forests - z0>0.1m
      MASK(:)=(ZNOT(:) > 1.0E-1) ! forest
      WHERE(MASK(:)) ICAT(:)=2

! all other lands, grass 0.001<z0<0.1m
      MASK(:)=((ZNOT(:) >= 1.0E-3).AND.(ZNOT(:) <= 1.0E-1)) ! grass
      WHERE(MASK(:)) ICAT(:)=3

! If sea ice covers > 50% of sea surface, treat as sea ice
      MASK(:)=(SEAICE(:) > 0.5) ! seaice
      WHERE(MASK(:)) ICAT(:)=5
!
!     First copy original values of ND,MD to ND0,MD0
      ND0(:,:)=ND(:,:)
      MD0(:,:,:)=MD(:,:,:)
!
!     Loop over grid boxes
!
      DELNDDEPTOT=0.0 ! total number deposited to surface
      DELMDDEPTOT=0.0 ! total NaCl mass deposited to surface
!
      DO JL=1,NBOX
       DZ=(RA*T(JL)/GG)*LOG(PLOWER(JL)/PUPPER(JL))
       JL_UP=JLABOVE(JL)
!      Calculate kinematic viscosity of air
       KVISC=DVISC(JL)/RHOA(JL)
!
       IF(JL_UP > 0) THEN
        DZ_UP=(RA*T(JL_UP)/GG)*LOG(PLOWER(JL_UP)/PUPPER(JL_UP))
       ENDIF
!
       DO IMODE=1,NMODES
        IF(MODE(IMODE)) THEN
         LOGIC1=(ND0(JL,IMODE) > NUM_EPS(IMODE))
         IF(JL_UP > 0) LOGIC2=(ND0(JL_UP,IMODE) > NUM_EPS(IMODE))
         IF(JL_UP <= 0) LOGIC2=.FALSE.
         IF (LOGIC1.OR.LOGIC2) THEN ! if significant ptcls in box
!                                                       or in box above
!
! .. Calculate 0th & 3rd moment avg. grav. settling velocities
! DEPENDS ON: ukca_vgrav_av_k
          VGRAV_AV_0=UKCA_VGRAV_AV_K(0,WETDP(JL,IMODE),SIGMAG(IMODE),   &
                      DVISC(JL),MFPA(JL),RHOPAR(JL,IMODE))
! DEPENDS ON: ukca_vgrav_av_k
          VGRAV_AV_3=UKCA_VGRAV_AV_K(3,WETDP(JL,IMODE),SIGMAG(IMODE),   &
                      DVISC(JL),MFPA(JL),RHOPAR(JL,IMODE))
!
          IF(JL_UP > 0) THEN ! if not top model level
!
! .. Calc. 0th & 3rd moment avg. grav. settling velocities (box above)
! DEPENDS ON: ukca_vgrav_av_k
           VGRAV_AV_0_UP=UKCA_VGRAV_AV_K(0,WETDP(JL_UP,IMODE),          &
                            SIGMAG(IMODE),DVISC(JL_UP),                 &
                            MFPA(JL_UP),RHOPAR(JL_UP,IMODE))
! DEPENDS ON: ukca_vgrav_av_k
           VGRAV_AV_3_UP=UKCA_VGRAV_AV_K(3,WETDP(JL_UP,IMODE),          &
                            SIGMAG(IMODE),DVISC(JL_UP),                 &
                            MFPA(JL_UP),RHOPAR(JL_UP,IMODE))
          ENDIF
!
          IF(SURTP(JL) < 2.0) THEN ! if at surf. set VGRAV to ddep vel
!
!          Calculate 0th & 3rd moment avg particle diffusion coeffs
! DEPENDS ON: ukca_dcoff_par_av_k
           DCOEF_AV_0=UKCA_DCOFF_PAR_AV_K(0,WETDP(JL,IMODE),            &
                      SIGMAG(IMODE),T(JL),DVISC(JL),MFPA(JL))
! DEPENDS ON: ukca_dcoff_par_av_k
           DCOEF_AV_3=UKCA_DCOFF_PAR_AV_K(3,WETDP(JL,IMODE),            &
                      SIGMAG(IMODE),T(JL),DVISC(JL),MFPA(JL))
!          Calculate 0th and 3rd moment avg. particle Schmidt number
           PS_AV_0=KVISC/DCOEF_AV_0
           PS_AV_3=KVISC/DCOEF_AV_3
!          Calculate particle collection efficiencies
!          -- For Brownian Diffusion
           EB_AV_0=PS_AV_0**(-YR(ICAT(JL)))
           EB_AV_3=PS_AV_3**(-YR(ICAT(JL)))
!          -- For Impaction
           IF (ICAT(JL) == 1.OR.ICAT(JL) == 5) THEN
!            Calculate stokes number for smooth surfaces
             SN_AV_0=VGRAV_AV_0*USTR(JL)*USTR(JL)/DVISC(JL)
             SN_AV_3=VGRAV_AV_3*USTR(JL)*USTR(JL)/DVISC(JL)
           ELSEIF (ICAT(JL) == 2.OR.ICAT(JL) == 3) THEN
!            Calculate stokes number for vegetated surfcaes
             SN_AV_0=VGRAV_AV_0*USTR(JL)/(GG*CR(ICAT(JL)))
             SN_AV_3=VGRAV_AV_3*USTR(JL)/(GG*CR(ICAT(JL)))
           ENDIF
           EIM_AV_0=(SN_AV_0/(ALPHA(ICAT(JL))+SN_AV_0))**2
           EIM_AV_3=(SN_AV_3/(ALPHA(ICAT(JL))+SN_AV_3))**2
!          -- For Interception
           IF (ICAT(JL) == 1.OR.ICAT(JL) == 5) THEN
             EIN=0.0
           ELSEIF (ICAT(JL) == 2.OR.ICAT(JL) == 3) THEN
             EIN=0.5*(WETDP(JL,IMODE)*WETDP(JL,IMODE)                   &
                     /CR(ICAT(JL))/CR(ICAT(JL)))
           ENDIF
!
           DZMID=(RA*T(JL)/GG)*LOG(PLOWER(JL)/PMID(JL))
!          Calculate aerodynamic resistance
           AR=LOG(DZMID/ZNOT(JL))/(VKARMN*USTR(JL))
!          Calculate surface resistance
           SR_AV_0=1.0/(3.0*USTR(JL)*(EB_AV_0+EIM_AV_0+EIN))
           SR_AV_3=1.0/(3.0*USTR(JL)*(EB_AV_3+EIM_AV_3+EIN))
!          Calculate deposition velocity
           VDEP_AV_0=VGRAV_AV_0+1.0/(AR+SR_AV_0)
           VDEP_AV_3=VGRAV_AV_3+1.0/(AR+SR_AV_3)
!
!          Set gravitational velocity to deposition velocity
           VGRAV_AV_0=VDEP_AV_0
           VGRAV_AV_3=VDEP_AV_3
          ELSE
           IF(SEDI_ON == 0) THEN
            VGRAV_AV_0=0.0
            VGRAV_AV_3=0.0
            VGRAV_AV_0_UP=0.0
            VGRAV_AV_3_UP=0.0
           ENDIF
          ENDIF
!
! .. below limits grav. settling so only falls 0.5 box max [numerical]
          VGRAV_LIM=0.5*DZ/DTC
          IF(VGRAV_AV_0 > VGRAV_LIM) VGRAV_AV_0=VGRAV_LIM
          IF(VGRAV_AV_3 > VGRAV_LIM) VGRAV_AV_3=VGRAV_LIM
!
! .. below limits grav. settling so only falls 0.5 box max (box above)
          VGRAV_LIM_UP=0.5*DZ_UP/DTC
          IF(VGRAV_AV_0_UP > VGRAV_LIM_UP) VGRAV_AV_0_UP=VGRAV_LIM_UP
          IF(VGRAV_AV_3_UP > VGRAV_LIM_UP) VGRAV_AV_3_UP=VGRAV_LIM_UP
!
! .. Sediment number and mass for mode using 0th & 3rd moment VGRAV_AVs
          IF(JL_UP > 0) THEN ! if not top model level
           TERMIN_1=ND0(JL_UP,IMODE)/DZ_UP
           TERMIN_2=SM(JL_UP)*RHOA(JL)/SM(JL)/RHOA(JL_UP)
           TERMIN_N=TERMIN_1*TERMIN_2*DTC*VGRAV_AV_0_UP
           DO ICP=1,NCP
            IF(COMPONENT(IMODE,ICP)) THEN
             TERMIN_1=ND0(JL_UP,IMODE)*MD0(JL_UP,IMODE,ICP)/DZ_UP
             TERMIN_M(ICP)=TERMIN_1*TERMIN_2*DTC*VGRAV_AV_3_UP
            ENDIF
           ENDDO
          ELSE ! if top model level
           TERMIN_N=0.0
           TERMIN_M(:)=0.0
          ENDIF
          TERMOUT_N=(ND0(JL,IMODE)/DZ)*DTC*VGRAV_AV_0
          DO ICP=1,NCP
           IF(COMPONENT(IMODE,ICP)) THEN
            TERMOUT_M(ICP)=(ND0(JL,IMODE)*MD0(JL,IMODE,ICP)/DZ)         &
                          *DTC*VGRAV_AV_3
           ENDIF
          ENDDO
!
! .. below calculates net change in number and mass concentration
          DELNSEDI=-(TERMIN_N-TERMOUT_N)
          DO ICP=1,NCP
           IF(COMPONENT(IMODE,ICP)) THEN
            DELMSEDI(ICP)=-(TERMIN_M(ICP)-TERMOUT_M(ICP))
           ENDIF
          ENDDO
!
! .. below re-calculates number and component mass for each mode
          IF(ABS(DELNSEDI) > 0.0) THEN
           NDNEW=ND0(JL,IMODE)-DELNSEDI
           MDT(JL,IMODE)=0.0
           DO ICP=1,NCP
            IF(COMPONENT(IMODE,ICP)) THEN
             MD(JL,IMODE,ICP)=                                          &
       (MD0(JL,IMODE,ICP)*ND0(JL,IMODE)-DELMSEDI(ICP))/NDNEW
             MDT(JL,IMODE)=MDT(JL,IMODE)+MD(JL,IMODE,ICP)
            ENDIF ! IF COMPONENT(ICP)
           ENDDO ! loop over cpts
           ND(JL,IMODE)=NDNEW
!
           IF(SURTP(JL) < 2.0) THEN ! if gridbox at surface
            DELNDDEPTOT=DELNDDEPTOT+TERMOUT_N*SM(JL)/RHOA(JL)
            DELMDDEPTOT=DELMDDEPTOT+TERMOUT_M(CP_CL)*SM(JL)/RHOA(JL)
            DO ICP=1,NCP
             IF(COMPONENT(IMODE,ICP)) THEN
              DELMDDEP=TERMOUT_M(ICP)
              IF(ICP == CP_SU) THEN
               IF((IMODE == 1).AND.(NMASDDEPSUNUCSOL > 0))              &
                          BUD_AER_MAS(JL,NMASDDEPSUNUCSOL)=             &
                          BUD_AER_MAS(JL,NMASDDEPSUNUCSOL)+DELMDDEP
               IF((IMODE == 2).AND.(NMASDDEPSUAITSOL > 0))              &
                          BUD_AER_MAS(JL,NMASDDEPSUAITSOL)=             &
                          BUD_AER_MAS(JL,NMASDDEPSUAITSOL)+DELMDDEP
               IF((IMODE == 3).AND.(NMASDDEPSUACCSOL > 0))              &
                          BUD_AER_MAS(JL,NMASDDEPSUACCSOL)=             &
                          BUD_AER_MAS(JL,NMASDDEPSUACCSOL)+DELMDDEP
               IF((IMODE == 4).AND.(NMASDDEPSUCORSOL > 0))              &
                          BUD_AER_MAS(JL,NMASDDEPSUCORSOL)=             &
                          BUD_AER_MAS(JL,NMASDDEPSUCORSOL)+DELMDDEP
              ENDIF
              IF(ICP == CP_BC) THEN
               IF((IMODE == 2).AND.(NMASDDEPBCAITSOL > 0))              &
                          BUD_AER_MAS(JL,NMASDDEPBCAITSOL)=             &
                          BUD_AER_MAS(JL,NMASDDEPBCAITSOL)+DELMDDEP
               IF((IMODE == 3).AND.(NMASDDEPBCACCSOL > 0))              &
                          BUD_AER_MAS(JL,NMASDDEPBCACCSOL)=             &
                          BUD_AER_MAS(JL,NMASDDEPBCACCSOL)+DELMDDEP
               IF((IMODE == 4).AND.(NMASDDEPBCCORSOL > 0))              &
                          BUD_AER_MAS(JL,NMASDDEPBCCORSOL)=             &
                          BUD_AER_MAS(JL,NMASDDEPBCCORSOL)+DELMDDEP
               IF((IMODE == 5).AND.(NMASDDEPBCAITINS > 0))              &
                          BUD_AER_MAS(JL,NMASDDEPBCAITINS)=             &
                          BUD_AER_MAS(JL,NMASDDEPBCAITINS)+DELMDDEP
              ENDIF
              IF(ICP == CP_OC) THEN
               IF((IMODE == 1).AND.(NMASDDEPOCNUCSOL > 0))              &
                          BUD_AER_MAS(JL,NMASDDEPOCNUCSOL)=             &
                          BUD_AER_MAS(JL,NMASDDEPOCNUCSOL)+DELMDDEP
               IF((IMODE == 2).AND.(NMASDDEPOCAITSOL > 0))              &
                          BUD_AER_MAS(JL,NMASDDEPOCAITSOL)=             &
                          BUD_AER_MAS(JL,NMASDDEPOCAITSOL)+DELMDDEP
               IF((IMODE == 3).AND.(NMASDDEPOCACCSOL > 0))              &
                          BUD_AER_MAS(JL,NMASDDEPOCACCSOL)=             &
                          BUD_AER_MAS(JL,NMASDDEPOCACCSOL)+DELMDDEP
               IF((IMODE == 4).AND.(NMASDDEPOCCORSOL > 0))              &
                          BUD_AER_MAS(JL,NMASDDEPOCCORSOL)=             &
                          BUD_AER_MAS(JL,NMASDDEPOCCORSOL)+DELMDDEP
               IF((IMODE == 5).AND.(NMASDDEPOCAITINS > 0))              &
                          BUD_AER_MAS(JL,NMASDDEPOCAITINS)=             &
                          BUD_AER_MAS(JL,NMASDDEPOCAITINS)+DELMDDEP
              ENDIF
              IF(ICP == CP_CL) THEN
               IF((IMODE == 3).AND.(NMASDDEPSSACCSOL > 0))              &
                          BUD_AER_MAS(JL,NMASDDEPSSACCSOL)=             &
                          BUD_AER_MAS(JL,NMASDDEPSSACCSOL)+DELMDDEP
               IF((IMODE == 4).AND.(NMASDDEPSSCORSOL > 0))              &
                          BUD_AER_MAS(JL,NMASDDEPSSCORSOL)=             &
                          BUD_AER_MAS(JL,NMASDDEPSSCORSOL)+DELMDDEP
              ENDIF
              IF(ICP == CP_SO) THEN
               IF((IMODE == 1).AND.(NMASDDEPSONUCSOL > 0))              &
                          BUD_AER_MAS(JL,NMASDDEPSONUCSOL)=             &
                          BUD_AER_MAS(JL,NMASDDEPSONUCSOL)+DELMDDEP
               IF((IMODE == 2).AND.(NMASDDEPSOAITSOL > 0))              &
                          BUD_AER_MAS(JL,NMASDDEPSOAITSOL)=             &
                          BUD_AER_MAS(JL,NMASDDEPSOAITSOL)+DELMDDEP
               IF((IMODE == 3).AND.(NMASDDEPSOACCSOL > 0))              &
                          BUD_AER_MAS(JL,NMASDDEPSOACCSOL)=             &
                          BUD_AER_MAS(JL,NMASDDEPSOACCSOL)+DELMDDEP
               IF((IMODE == 4).AND.(NMASDDEPSOCORSOL > 0))              &
                          BUD_AER_MAS(JL,NMASDDEPSOCORSOL)=             &
                          BUD_AER_MAS(JL,NMASDDEPSOCORSOL)+DELMDDEP
              ENDIF
             ENDIF
            ENDDO
           ENDIF ! if gridbox at surface

          ENDIF ! if signficant change in number to be made

         ENDIF ! if ND>0.0 in this or above box
        ENDIF ! IF MODE(IMODE)
       ENDDO ! loop over modes
!
      ENDDO ! loop over grid boxes

      RETURN
      END SUBROUTINE UKCA_DDEPAER_INCL_SEDI
#endif
