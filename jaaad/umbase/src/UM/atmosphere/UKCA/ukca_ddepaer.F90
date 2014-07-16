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
!    Calculates aerosol dry deposition (and sedimentation).
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
      SUBROUTINE UKCA_DDEPAER(NBOX,ND,MD,MDT,                           &
          RHOPAR,ZNOT,SEAICE,                                           &
          DTC,WETDP,USTR,PMID,PUPPER,PLOWER,T,SURTP,                    &
          RHOA,MFPA,DVISC,BUD_AER_MAS)
!
! Purpose
! -------
! Calculates aerosol dry deposition (and sedimentation).
! Based on the parameterisation of Zhang et al (2001) which
! uses the method in the model of Slinn (1982).
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
! Note --- only evaluates sedimentation at lowest gridbox if this
!          routine is used --- sedimentation at higher levels neglected.
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
! SEAICE    : Fraction of horizontal gridbox area containing seaice
! DTC       : Chemical timestep (s)
! WETDP     : Wet diameter for ptcl with dry diameter DRYDP (m)
! USTR      : Friction velocity(ms-1)
! PMID      : Centre level pressure (Pa)
! PUPPER    : Pressure at box upper interface (Pa)
! PLOWER    : Pressure at box lower interface (Pa)
! T         : Centre level temperature (K)
! SURTP     : Surface type [0=sea-surf,1=land-surf,2=above-surf]
! RHOA      : Air density (kg/m3)
! MFPA      : Mean free path of air (m)
! DVISC     : Dynamic viscosity of air (kg m-1 s-1)
!
! Outputs
! -------
! Updated particle number density ND (/cm3)
! Updated particle avg mass MD (molecules/particle)
! BUD_AER_MAS : Aerosol mass budgets (mlcls/cc/tstep)
!
! Local Variables
! ---------------
! PS_AV_0    : 0th moment avg particle Schmidt Number
! PS_AV_3    : 3rd moment avg particle Schmidt Number
! KVISC      : Kinematic viscosity of air (m2 s-1)
! VGRAV_AV_0 : 0th moment avg grav. settling vel. (m/s)
! VGRAV_AV_3 : 3rd moment avg grav. settling vel. (m/s)
! VDEP_AV_0  : 0th moment avg deposition velocity (m/s)
! VDEP_AV_3  : 3rd moment avg deposition velocity (m/s)
! DCOEF_AV_0 : 0th moment avg particle diffusion coefficient(m2/s)
! DCOEF_AV_3 : 3rd moment avg particle diffusion coefficient(m2/s)
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
! NUM_EPS   : Value of NEWN below which do not recalculate MD (per cc)
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
      USE UKCA_CONSTANTS,  ONLY: GG, VKARMN, PPI, RA, ZBOLTZ
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES
!
      IMPLICIT NONE
!
!     Subroutine interface
      INTEGER :: NBOX
      REAL    :: ND(NBOX,NMODES)
      REAL    :: MD(NBOX,NMODES,NCP)
      REAL    :: MDT(NBOX,NMODES)
      REAL    :: RHOPAR(NBOX,NMODES)
      REAL    :: ZNOT(NBOX)
      REAL    :: SEAICE(NBOX)
      REAL    :: DTC
      REAL    :: WETDP(NBOX,NMODES)
      REAL    :: USTR(NBOX)
      REAL    :: PMID(NBOX)
      REAL    :: PUPPER(NBOX)
      REAL    :: PLOWER(NBOX)
      REAL    :: T(NBOX)
      REAL    :: SURTP(NBOX)
      REAL    :: RHOA(NBOX)
      REAL    :: MFPA(NBOX)
      REAL    :: DVISC(NBOX)
      REAL    :: BUD_AER_MAS(NBOX,NBUDAER)
      REAL    :: DELNDDEP
      REAL    :: DELMDDEP
!
!     Local Variables
      INTEGER :: ICAT(NBOX)
      INTEGER :: JL
      INTEGER :: IMODE
      INTEGER :: ICP
      LOGICAL :: MASK(NBOX)
      REAL    :: PS_AV_0
      REAL    :: PS_AV_3
      REAL    :: KVISC
      REAL    :: VGRAV_AV_0
      REAL    :: VGRAV_AV_3
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
      REAL    :: MTOT
      REAL    :: MCPTOT
      REAL    :: NEWN
      REAL    :: UKCA_DCOFF_PAR_AV_K
      REAL    :: UKCA_VGRAV_AV_K
      REAL    :: DZMID
      REAL    :: DZ
      REAL    :: SSIGMA
      REAL, PARAMETER :: YR(5) = (/  0.5, 0.56, 0.54,0.0,0.54/)
      REAL, PARAMETER :: CR(5) = (/ 0.00,5.E-3,2.E-3,0.0,0.00/)
      REAL, PARAMETER :: ALPHA(5) = (/100.0,  1.0,  1.2,0.0,50.0/)
!
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
!     Loop over grid boxes
      DO IMODE=1,NMODES
       IF(MODE(IMODE)) THEN
        SSIGMA=SIGMAG(IMODE)
        DO JL=1,NBOX
         IF(SURTP(JL) < 1.5) THEN  ! if at surface
!
          IF (ND(JL,IMODE) > NUM_EPS(IMODE)) THEN
!
!          Calculate 0th & 3rd moment avg particle diffusion coeffs
! DEPENDS ON: ukca_dcoff_par_av_k
           DCOEF_AV_0=UKCA_DCOFF_PAR_AV_K(0,WETDP(JL,IMODE),SSIGMA,     &
               T(JL),DVISC(JL),MFPA(JL))
! DEPENDS ON: ukca_dcoff_par_av_k
           DCOEF_AV_3=UKCA_DCOFF_PAR_AV_K(3,WETDP(JL,IMODE),SSIGMA,     &
               T(JL),DVISC(JL),MFPA(JL))
!          Calculate kinematic viscosity of air
           KVISC=DVISC(JL)/RHOA(JL)
!          Calculate 0th and 3rd moment avg. particle Schmidt number
           PS_AV_0=KVISC/DCOEF_AV_0
           PS_AV_3=KVISC/DCOEF_AV_3
!          Calculate 0th & 3rd moment avg. grav. settling velocities
! DEPENDS ON: ukca_vgrav_av_k
           VGRAV_AV_0=UKCA_VGRAV_AV_K(0,WETDP(JL,IMODE),SSIGMA,         &
               DVISC(JL),MFPA(JL),RHOPAR(JL,IMODE))
! DEPENDS ON: ukca_vgrav_av_k
           VGRAV_AV_3=UKCA_VGRAV_AV_K(3,WETDP(JL,IMODE),SSIGMA,         &
               DVISC(JL),MFPA(JL),RHOPAR(JL,IMODE))
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

!          Calculate aerodynamic resistance
           DZMID=(RA*T(JL)/GG)*LOG(PLOWER(JL)/PMID(JL))
           DZ=(RA*T(JL)/GG)*LOG(PLOWER(JL)/PUPPER(JL))
           AR=LOG(DZMID/ZNOT(JL))/(VKARMN*USTR(JL))
!          Calculate surface resistance
           SR_AV_0=1.0/(3.0*USTR(JL)*(EB_AV_0+EIM_AV_0+EIN))
           SR_AV_3=1.0/(3.0*USTR(JL)*(EB_AV_3+EIM_AV_3+EIN))
!          Calculate deposition velocity
           VDEP_AV_0=VGRAV_AV_0+1.0/(AR+SR_AV_0)
           VDEP_AV_3=VGRAV_AV_3+1.0/(AR+SR_AV_3)
           DELNDDEP=ND(JL,IMODE)*(1.0-EXP(-VDEP_AV_0*DTC/DZ))
           NEWN=ND(JL,IMODE)-DELNDDEP
           IF(NEWN > NUM_EPS(IMODE)) THEN
            MTOT=ND(JL,IMODE)*MDT(JL,IMODE)
            MDT(JL,IMODE)=MTOT*EXP(-VDEP_AV_3*DTC/DZ)/NEWN
            DO ICP=1,NCP
             IF(COMPONENT(IMODE,ICP)) THEN
              MCPTOT=ND(JL,IMODE)*MD(JL,IMODE,ICP)
              MD(JL,IMODE,ICP)=MCPTOT*EXP(-VDEP_AV_3*DTC/DZ)/NEWN
              DELMDDEP=MCPTOT*(1.0-EXP(-VDEP_AV_3*DTC/DZ))
              IF(ICP == CP_SU) THEN
               IF((IMODE == 1).AND.(NMASDDEPSUNUCSOL > 0))              &
       BUD_AER_MAS(JL,NMASDDEPSUNUCSOL)=                                &
       BUD_AER_MAS(JL,NMASDDEPSUNUCSOL)+DELMDDEP
               IF((IMODE == 2).AND.(NMASDDEPSUAITSOL > 0))              &
       BUD_AER_MAS(JL,NMASDDEPSUAITSOL)=                                &
       BUD_AER_MAS(JL,NMASDDEPSUAITSOL)+DELMDDEP
               IF((IMODE == 3).AND.(NMASDDEPSUACCSOL > 0))              &
       BUD_AER_MAS(JL,NMASDDEPSUACCSOL)=                                &
       BUD_AER_MAS(JL,NMASDDEPSUACCSOL)+DELMDDEP
               IF((IMODE == 4).AND.(NMASDDEPSUCORSOL > 0))              &
       BUD_AER_MAS(JL,NMASDDEPSUCORSOL)=                                &
       BUD_AER_MAS(JL,NMASDDEPSUCORSOL)+DELMDDEP
              ENDIF
              IF(ICP == CP_BC) THEN
               IF((IMODE == 2).AND.(NMASDDEPBCAITSOL > 0))              &
       BUD_AER_MAS(JL,NMASDDEPBCAITSOL)=                                &
       BUD_AER_MAS(JL,NMASDDEPBCAITSOL)+DELMDDEP
               IF((IMODE == 3).AND.(NMASDDEPBCACCSOL > 0))              &
       BUD_AER_MAS(JL,NMASDDEPBCACCSOL)=                                &
       BUD_AER_MAS(JL,NMASDDEPBCACCSOL)+DELMDDEP
               IF((IMODE == 4).AND.(NMASDDEPBCCORSOL > 0))              &
       BUD_AER_MAS(JL,NMASDDEPBCCORSOL)=                                &
       BUD_AER_MAS(JL,NMASDDEPBCCORSOL)+DELMDDEP
               IF((IMODE == 5).AND.(NMASDDEPBCAITINS > 0))              &
       BUD_AER_MAS(JL,NMASDDEPBCAITINS)=                                &
       BUD_AER_MAS(JL,NMASDDEPBCAITINS)+DELMDDEP
              ENDIF
              IF(ICP == CP_OC) THEN
               IF((IMODE == 1).AND.(NMASDDEPOCNUCSOL > 0))              &
       BUD_AER_MAS(JL,NMASDDEPOCNUCSOL)=                                &
       BUD_AER_MAS(JL,NMASDDEPOCNUCSOL)+DELMDDEP
               IF((IMODE == 2).AND.(NMASDDEPOCAITSOL > 0))              &
       BUD_AER_MAS(JL,NMASDDEPOCAITSOL)=                                &
       BUD_AER_MAS(JL,NMASDDEPOCAITSOL)+DELMDDEP
               IF((IMODE == 3).AND.(NMASDDEPOCACCSOL > 0))              &
       BUD_AER_MAS(JL,NMASDDEPOCACCSOL)=                                &
       BUD_AER_MAS(JL,NMASDDEPOCACCSOL)+DELMDDEP
               IF((IMODE == 4).AND.(NMASDDEPOCCORSOL > 0))              &
       BUD_AER_MAS(JL,NMASDDEPOCCORSOL)=                                &
       BUD_AER_MAS(JL,NMASDDEPOCCORSOL)+DELMDDEP
               IF((IMODE == 5).AND.(NMASDDEPOCAITINS > 0))              &
       BUD_AER_MAS(JL,NMASDDEPOCAITINS)=                                &
       BUD_AER_MAS(JL,NMASDDEPOCAITINS)+DELMDDEP
              ENDIF
              IF(ICP == CP_CL) THEN
               IF((IMODE == 3).AND.(NMASDDEPSSACCSOL > 0))              &
       BUD_AER_MAS(JL,NMASDDEPSSACCSOL)=                                &
       BUD_AER_MAS(JL,NMASDDEPSSACCSOL)+DELMDDEP
               IF((IMODE == 4).AND.(NMASDDEPSSCORSOL > 0))              &
       BUD_AER_MAS(JL,NMASDDEPSSCORSOL)=                                &
       BUD_AER_MAS(JL,NMASDDEPSSCORSOL)+DELMDDEP
              ENDIF
              IF(ICP == CP_SO) THEN
               IF((IMODE == 1).AND.(NMASDDEPSONUCSOL > 0))              &
       BUD_AER_MAS(JL,NMASDDEPSONUCSOL)=                                &
       BUD_AER_MAS(JL,NMASDDEPSONUCSOL)+DELMDDEP
               IF((IMODE == 2).AND.(NMASDDEPSOAITSOL > 0))              &
       BUD_AER_MAS(JL,NMASDDEPSOAITSOL)=                                &
       BUD_AER_MAS(JL,NMASDDEPSOAITSOL)+DELMDDEP
               IF((IMODE == 3).AND.(NMASDDEPSOACCSOL > 0))              &
       BUD_AER_MAS(JL,NMASDDEPSOACCSOL)=                                &
       BUD_AER_MAS(JL,NMASDDEPSOACCSOL)+DELMDDEP
               IF((IMODE == 4).AND.(NMASDDEPSOCORSOL > 0))              &
       BUD_AER_MAS(JL,NMASDDEPSOCORSOL)=                                &
       BUD_AER_MAS(JL,NMASDDEPSOCORSOL)+DELMDDEP
              ENDIF
             ENDIF
            ENDDO
            ND(JL,IMODE)=NEWN
           ENDIF ! if NEWN > NUM_EPS
          ENDIF ! if ND > NUM_EPS
         ENDIF ! if SURTP < 1.5 (at surface)
        ENDDO ! loop over grid boxes
       ENDIF ! if MODE(IMODE)
      ENDDO ! loop over modes

      RETURN
      END SUBROUTINE UKCA_DDEPAER
#endif
