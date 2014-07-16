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
!    Subroutine which calculates the condensation coefficients
!    for condensable cpt vapour condensing on a particle with
!    radius RP. Includes a switch to either use Fuchs (1964)
!    or Fuchs-Sutugin (1971).
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
      SUBROUTINE UKCA_COND_COFF_V(NV,MASK,RP,TSQRT,AIRDM3,RHOA,         &
                                  MMCG,MMAIR,SE,DMOL,IFUCHS,CC)
!----------------------------------------------------------------
!
! Purpose
! -------
! Subroutine which calculates the condensation coefficients
! for condensable cpt vapour condensing on a particle with
! radius RP. Includes a switch to either use Fuchs (1964)
! or Fuchs-Sutugin (1971).
!
! Inputs
! ------
! NV     : total number of values
! MASK   : Mask where to calculate values
! RP     : Radius of aerosol particle
! TSQRT  : Square-root of mid-level air temperature (K)
! AIRDM3 : Number density of air (m3)
! RHOA   : Air density (kg/m3)
! MMCG   : Molar mass of condensing gas (kg per mole)
! MMAIR  : Molar mass of air (kg mol-1)
! SE     : Sticking efficiency [taken as 0.3 from Raes et al 1992]
! DMOL   : Molecular diameter of condensable (m)
! IFUCHS : Switch for Fuchs (1964) or Fuchs-Sutugin (1971)
!
! Outputs
! -------
! CC : Condensation coeff for cpt onto ptcl (m^3s^-1)
!
! Local Variables
! ---------------
! VEL_CP   : Thermal velocity of condensable gas (ms-1)
! MFP_CP   : Mean free path of condensable gas (m)
! DCOFF_CP : Diffusion coefficient of condensable gas in air (m2s-1)
! KN       : Knudsen number
! FKN      : Correction factor for molecular effects
! AKN      : Corr fac for limitations in interfacial mass transport
! DENOM    : Denominator in Fuchs (1964) condensation coeff expression
! ZZ       : Ratio of molar masses of condensable gas and air
! TERM1,TERM2,.. : Terms in evaluation of condensation coefficient
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! ZBOLTZ : Boltzmann's constant (kg m2 s-2 K-1 molec-1)
! RA     : Dry air gas constant = 287.05 Jkg^-1 K^-1
! RR     : Universal gas constant = 8.314 K mol^-1 K^-1
! PPI    : 3.1415927...
! AVC    : Avogadros constant (mol-1)
!
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS
!
      IMPLICIT NONE
!
! Subroutine interface:
      INTEGER :: NV
      INTEGER :: IFUCHS
      LOGICAL :: MASK(NV)
      REAL    :: RP(NV)
      REAL    :: TSQRT(NV)
      REAL    :: AIRDM3(NV)
      REAL    :: RHOA(NV)
      REAL    :: CC(NV)
      REAL    :: MMCG
      REAL    :: MMAIR
      REAL    :: SE
      REAL    :: DMOL
!
! .. Local variables
      REAL    :: VEL_CP(NV)
      REAL    :: MFP_CP(NV)
      REAL    :: DCOFF_CP(NV)
      REAL    :: KN(NV)
      REAL    :: FKN(NV)
      REAL    :: AKN(NV)
      REAL    :: DENOM(NV)
      REAL    :: ZZ
      REAL    :: TERM1
      REAL    :: TERM2
      REAL    :: TERM3
      REAL    :: TERM4
      REAL    :: TERM5
      REAL    :: TERM6

      TERM1=SQRT(8.0*RR/(PPI*MMCG))
! .. used in calcn of thermal velocity of condensable gas

      ZZ=MMCG/MMAIR
      TERM2=1.0/(PPI*SQRT(1.0+ZZ)*DMOL*DMOL)
! .. used in calcn of MFP of condensable gas (S&P, pg 457, eq 8.11)

      TERM3=(3.0/(8.0*AVC*DMOL*DMOL))
      TERM4=SQRT((RA*MMAIR*MMAIR/(2.0*PPI))*((MMCG+MMAIR)/MMCG))
      TERM5=TERM3*TERM4
! .. used in calcn of diffusion coefficient of condensable gas
!
      TERM6=4.0e6*PPI
! .. used in calculation of condensation coefficient
!
      CC(:)=0.0
      WHERE(MASK(:))
!       Calc. thermal velocity of condensable gas
        VEL_CP(:)=TERM1*TSQRT(:)
!       Calc. mean free path of condensable gas
        MFP_CP(:)=TERM2/AIRDM3(:)
!       Calculate diffusion coefficient of condensable gas
        DCOFF_CP(:)=TERM5*TSQRT(:)/RHOA(:)
      ENDWHERE
!
      IF(IFUCHS == 1) THEN
!      If IFUCHS=1 use basic Fuchs (1964) expression
       WHERE(MASK(:))
        DENOM(:)=                                                       &
       (4.0*DCOFF_CP(:)/(SE*VEL_CP(:)*RP(:)))+(RP(:)/(RP(:)+MFP_CP(:)))
!       Calc. condensation coefficient
        CC(:)=TERM6*DCOFF_CP(:)*RP(:)/DENOM(:)
       ENDWHERE
      ENDIF

      IF(IFUCHS == 2) THEN
!      If IFUCHS=2 use basic Fuchs-Sutugin (1971) expression
       WHERE(MASK(:))
!       Calculate Knudsen number of condensable gas w.r.t. particle
        KN(:)=MFP_CP(:)/RP(:)
!       Calc. corr. factor for molecular effects
        FKN(:)=(1.0+KN(:))/(1.0+1.71*KN(:)+1.33*KN(:)*KN(:))
!       Calc. corr. factor for limitations in interfacial mass transport
        AKN(:)=1.0/(1.0+1.33*KN(:)*FKN(:)*(1.0/SE-1.0))
!       Calc. condensation coefficient
        CC(:)=TERM6*DCOFF_CP(:)*RP(:)*FKN(:)*AKN(:)
       ENDWHERE
      ENDIF

      END SUBROUTINE UKCA_COND_COFF_V
#endif
