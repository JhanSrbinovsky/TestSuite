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
!    Calculate the H2SO4-H2O binary nucleation rate of sulfate particles.
!    Parameterization based on
!    Method 1) Pandis et al (1994), JGR, vol.99, no. D8, pp 16,945-16,957.
!           2) Kulmala et al (1998), JGR, vol.103, no. D7, pp 8,301-8,307.
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
      SUBROUTINE UKCA_CALCNUCRATE(NBOX,DTZ,T,S,RH,AIRD,                 &
              H2SO4,DELTAH2SO4,JRATE)
!---------------------------------------------------------------
!
! Purpose
! -------
! Calculate the H2SO4-H2O binary nucleation rate of sulfate particles.
! Parameterization based on
! Method 1) Pandis et al (1994), JGR, vol.99, no. D8, pp 16,945-16,957.
!        2) Kulmala et al (1998), JGR, vol.103, no. D7, pp 8,301-8,307.
!
! Parameters
! ----------
! None
!
! Inputs:
! ------
! NBOX    : Number of grid boxes
! DTZ     : Timestep for nucl/cond competition (s)
! T       : Mid-level temperature (K)
! S       : Specific humidity (kg/kg)
! RH      : Relative humidity (dimensionless 0-1)
! AIRD    : Number density of air (per cc)
! H2SO4   : H2SO4 number density (cm-3)
!
! Outputs:
! -------
! H2SO4      : H2SO4 number density (cm-3)
! DELTAH2SO4 : Change in H2SO4 due to nucleation (molecules/DTZ)
! JRATE      : H2SO4 depletion rate by nucleation (/cc/s)
!
! Local Variables
! ---------------
! METHOD    : Switch : 1=use Pandis94 method, 2=use Kulmala98 method
! CORT      : Air temperature corrected to lie in range (K)
! CORRH     : Relative humidity corr. to lie in range (0-1)
! JPAN      : Pandis H2SO4/H2O nucleation rate (cm^-3 s^-1)
! JKUL      : Kulmala H2SO4/H2O nucleation rate (cm^-3 s^-1)
! H2SO4OLD2 : Old number density of H2SO4(g) (cm^-3)
! ACONS     : Constant in calculation of JKUL and JPAN
! BCONS     : Constant in calculation of JKUL and JPAN
! DELTA     : Parameter in calculation of JKUL
! A2        : Parameter in calculation of JKUL
! B2        : Parameter in calculation of JKUL
! C2        : Parameter in calculation of JKUL
! D2        : Parameter in calculation of JKUL
! E2        : Parameter in calculation of JKUL
! NCRIT     : H2SO4 no. dens. giving JKUL=1 cm^-3s^-1 (cm^-3)
! NSULF     : ln(H2SO4/NCRIT)
! PP        : partial pressure of H2SO4 (Pa)
! SVP       : saturation vapour pressure of H2SO4 (Pa)
! RAC       : relative acidity = PP/SVP
! NWP       : water vapour concentration (cm-3)
! XAL       : H2SO4 mole fraction
! THETA     : =chi in Kulmala et al (98) eq 20 = ln(JKUL)
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! ZBOLTZ  : Boltzman Constant (kg m2 s-2 K-1 molec-1)
! NMOL    : Number of molecules per particle at nucleation
! CONC_EPS: Threshold for concentration (molecules per cc)
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NCP     : Number of possible aerosol components
! CP_SU   : Index of component in which sulfate is stored
!
! References
! ----------
! * Pandis et al (1994) "The relationship between
!   DMS flux and CCN concentration in remote marine
!   regions", JGR, vol. 99, no. D8, pp. 16,945-16,957.
!
! * Kulmala et al (1998) "Parameterizations for
!   sulfuric acid/water nucleation rates",
!   JGR, vol. 103, no. D7, pp. 8301-8307.
!
!
!----------------------------------------------------------------------
      USE UKCA_CONSTANTS
      USE UKCA_MODE_SETUP

      IMPLICIT NONE

!     Arguments
      INTEGER :: NBOX
      REAL    :: H2SO4(NBOX)
      REAL    :: DTZ
      REAL    :: T(NBOX)
      REAL    :: S(NBOX)
      REAL    :: RH(NBOX)
      REAL    :: AIRD(NBOX)
      REAL    :: JRATE(NBOX)
      REAL    :: DELTAH2SO4(NBOX)

!     Local variables
      INTEGER :: JL
      INTEGER :: METHOD
      REAL :: CORT(NBOX)
      REAL :: CORRH(NBOX)
      REAL :: JPAN
      REAL :: H2SO4OLD2
      REAL :: ACONS
      REAL :: BCONS
      REAL :: DELTA
      REAL :: A2
      REAL :: B2
      REAL :: C2
      REAL :: D2
      REAL :: E2
      REAL :: NCRIT
      REAL :: NSULF
      REAL :: PP
      REAL :: SVP
      REAL :: NWP
      REAL :: XAL
      REAL :: RAC
      REAL :: THETA
      REAL :: JKUL
!
!     Set nucleation method
!     1=Pandis(94)
!     2=Kulmala(98)
!
      METHOD=2
!
      DO JL=1,NBOX
!      Correct relative humidity to lie within parametrisation range
       CORRH(JL)=RH(JL)
       IF (CORRH(JL) < 0.1) CORRH(JL)=0.1
       IF (CORRH(JL) > 0.9) CORRH(JL)=0.9
!      Correct temperature to lie within parametrisation range
       CORT(JL)=T(JL)
       IF (CORT(JL) < 233.0) CORT(JL)=233.0
       IF (CORT(JL) > 298.0) CORT(JL)=298.0
      ENDDO
!
      IF (METHOD == 1) THEN
       DO JL=1,NBOX
!       Calculate H2SO4-H2O binary nucleation rate using Pandis(94)
        ACONS=10.0**(7.0 - (64.24 + 4.7*CORRH(JL)))
        BCONS=6.13+1.95*CORRH(JL)
        JPAN=ACONS*(H2SO4(JL)**BCONS)
        IF (JPAN > 1.0E-5) THEN
!        Calculate new H2S04 concentration
         H2SO4OLD2=H2SO4(JL)
         H2SO4(JL)=(H2SO4(JL)**(1.0 - BCONS) +                          &
               (BCONS - 1.0)*NMOL*ACONS*DTZ )                           &
                ** ((1.0/(1.0 - BCONS) ) )
         IF (H2SO4(JL) > H2SO4OLD2) H2SO4(JL) = H2SO4OLD2
         IF (H2SO4(JL) < 0.0) H2SO4(JL) = 0.0
!        Calculate change in sulfuric acid
         DELTAH2SO4(JL)=H2SO4OLD2-H2SO4(JL)
        ELSE
         DELTAH2SO4(JL)=0.0
        ENDIF
       ENDDO
      ENDIF

      IF (METHOD == 2) THEN
       DO JL=1,NBOX
!       Calculate H2SO4-H2O binary nucleation rate using Kulmala(98)
        IF (H2SO4(JL) > CONC_EPS) THEN
!        Calculate parametrised variables
         DELTA=1.0+(CORT(JL)-273.15)/273.15
         A2=25.1289-(4890.8/CORT(JL))-2.2479*DELTA*CORRH(JL)
         B2=(7643.4/CORT(JL))-1.9712*DELTA/CORRH(JL)
         C2=-1743.3/CORT(JL)
!        Calculate H2SO4 needed to get J=1cm-3 s-1
         NCRIT=EXP (-14.5125 + 0.1335*CORT(JL) -                        &
          10.5462*CORRH(JL) + 1958.4*CORRH(JL)/CORT(JL))
         NSULF=LOG(H2SO4(JL)/NCRIT)
!        Calculate saturation vapour pressure of H2SO4 (Pa)
         SVP=EXP( 27.78492066 - 10156.0/T(JL))
!        Calculate water vapour concentration (cm-3)
         NWP=AIRD(JL)*1.609*S(JL)
!        Calculate partial pressure of H2SO4 (Pa)
         PP=H2SO4(JL)*1.E6*ZBOLTZ*T(JL)
!        Calculate relative acidity
         RAC=PP/SVP
         IF (NWP > 0.0) THEN
          D2=1.2233 - (0.0154*RAC/(RAC+CORRH(JL)))-                     &
                  0.0415*LOG(NWP) + 0.0016*CORT(JL)
         ELSE
          D2=1.2233 - (0.0154*RAC/(RAC+CORRH(JL)))                      &
                  + 0.0016*CORT(JL)
         ENDIF
         E2 = 0.0102
         XAL=D2+E2*LOG(H2SO4(JL))
! .. NSULF is ln(H2SO4/NCRIT)
         THETA=A2*NSULF+B2*XAL+C2
! .. JKUL is nucleation rate in particles/cc/s
         JKUL=EXP(THETA)
!         write(6,*) 'JKUL (must be > 1.0e-3 to nucleate)=',JKUL
!        Calculate new H2SO4
         IF (JKUL > 1.E-3) THEN
!          write(6,*) 'JKUL>1.0E-3,here1'
          ACONS=(EXP(B2*D2+C2))*(1.0/NCRIT)**A2
          BCONS=A2 + B2*E2
!         Calculate new H2S04 concentrations
          H2SO4OLD2 = H2SO4(JL)
          H2SO4(JL)=( H2SO4(JL)**(1.0-BCONS) +                          &
                  (BCONS - 1.0)*NMOL*ACONS*DTZ)                         &
                  **((1.0 / (1.0 - BCONS)))
!         Calculate change in sulfuric acid
          IF (H2SO4(JL) > H2SO4OLD2) H2SO4(JL)=H2SO4OLD2
          IF (H2SO4(JL) < 0.0) H2SO4(JL)=0.0
          DELTAH2SO4(JL)=H2SO4OLD2-H2SO4(JL)
!         Calculate depletion rate of H2SO4 (molecules/cc/s)
          JRATE(JL)=NMOL*ACONS*H2SO4OLD2**BCONS
         ELSE
          DELTAH2SO4(JL)=0.0
          JRATE(JL)=0.0
         ENDIF
        ELSE
         DELTAH2SO4(JL)=0.0
         JRATE(JL)=0.0
        ENDIF
       ENDDO
      ENDIF

      IF((METHOD /= 1).AND.(METHOD /= 2)) THEN
! DEPENDS ON: ereport
       CALL EREPORT('UKCA_CALCNUCRATE',1,'METHOD not 1 or 2')
      ENDIF

      RETURN
      END SUBROUTINE UKCA_CALCNUCRATE
#endif
