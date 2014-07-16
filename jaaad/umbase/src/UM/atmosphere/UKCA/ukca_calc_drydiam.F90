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
!   Calculates geometric mean dry diameter for multi-component
!   aerosol population which is lognormally distributed with
!   number concentration ND in mode, component mass concentration
!   MD in mode, component density RHOCOMP and component molecular mass MM
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
      SUBROUTINE UKCA_CALC_DRYDIAM(NBOX,ND,MD,MDT,DRYDP,DVOL,VERBOSE)

!-------------------------------------------------------------
!
!     Calculates geometric mean dry diameter for multi-component
!     aerosol population which is lognormally distributed with
!     number concentration ND in mode,
!     cpt mass concentration MD in mode,
!     cpt density RHOCOMP and cpt molecular mass MM
!
!     Calculate dry volume per particle using composition info as:
!
!     DVOL = sum_cp { MD(ICP)*MM(ICP)/(AVC*RHOCOMP(ICP)) }
!
!     Where AVC is Avogadro's constant. Then, from Jacobsen,
!     "Fundamentals of Atmospheric Modeling", pp. 412, have
!
!     dry volume conc. = ND*(PI/6)*(Dp^3)*exp{9/2 log^2(sigma_g)}
!
!     i.e. DVOL  = (PI/6)*(Dp^3)*exp{9/2 log^2(sigma_g)}
!
!     where Dp is the number mean dry diameter,
!     and sigma_g is the geometric standard deviation.
!
!     Then calcaulte Dp as:
!
!     Dp=CBRT( DVOL*(6/PI)/EXP{9.2 log^2(sigma_g)} )
!
!     Inputs
!     ------
!     NBOX      : Number of grid boxes
!     ND        : Aerosol ptcl no. concentration (ptcls per cc)
!     MD        : Component median aerosol mass (molecules per ptcl)
!     MDT       : Total median aerosol mass (molecules per ptcl)
!     VERBOSE   : Switch for level of verbosity
!
!     Outputs
!     -------
!     DRYDP     : Median particle dry diameter for each mode (m)
!     DVOL      : Median particle dry volume for each mode (m^3)
!
!     Local Variables
!     ---------------
!     None
!
!     Inputted by module UKCA_CONSTANTS
!     ---------------------------------
!     PPI       : 3.1415927
!     AVC       : Avogadro's constant (per mole)
!     MMSUL     : Molar mass of a pure H2SO4 aerosol (kg per mole)
!     RHOSUL    : Mass density of a pure H2SO4 aerosol (kg per m^3)
!
!     Inputted by module UKCA_MODE_SETUP
!     ----------------------------------
!     NMODES    : Number of aerosol modes
!     NCP       : Number of aerosol components
!     MM        : Molecular mass of each component
!     RHOCOMP   : Densities (dry) of each component (kg/m^3)
!     MODE      : Which modes are being carried
!     COMPONENT : Which components are in each of modes
!     DDPLIM0   : Lower limits for dry diameter for each mode (m)
!     DDPLIM1   : Upper limits for dry diameter for each mode (m)
!     MFRAC_0   : Initial mass fraction to set when no particles.
!     X         : EXP((9/2)*LOG^2(SIGMA_G))
!     NUM_EPS   : Value of NEWN below which do not carry out process
!
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS
      USE UKCA_MODE_SETUP

      IMPLICIT NONE

! Interface
      INTEGER :: NBOX
      INTEGER :: VERBOSE
      REAL    :: ND(NBOX,NMODES)
      REAL    :: MD(NBOX,NMODES,NCP)
      REAL    :: MDT(NBOX,NMODES)
      REAL    :: DRYDP(NBOX,NMODES)
      REAL    :: DVOL(NBOX,NMODES)

! Local variables
      INTEGER :: JL,IMODE,ICP
      REAL    :: CBRT
      REAL    :: DDPCUB(NBOX)
      REAL    :: SIXOVRPIX(NMODES)
      LOGICAL :: MASK(NBOX)
      LOGICAL :: TOOHI,TOOLO
      CHARACTER*20 :: STRERR
      REAL    :: RATIO1(NCP)
      REAL    :: RATIO2(NMODES)
!
! Below is over NCP
      RATIO1(:)=MM(:)/AVC/RHOCOMP(:)
!
! Below is over NMODES
      SIXOVRPIX(:)=6.0/(PPI*X(:))
      RATIO2(:)=MMSUL*MMID(:)/(AVC*RHOSUL)
!
      DO IMODE=1,NMODES
       IF(MODE(IMODE)) THEN
        MASK(:) = ND(:,IMODE) > NUM_EPS(IMODE)
        WHERE(MASK(:))
         DVOL(:,IMODE)=0.0
        ELSEWHERE
         DVOL(:,IMODE)=MMID(IMODE)*MMSUL/(AVC*RHOSUL)
        ENDWHERE
!     calculate particle dry volume using composition info
        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
          WHERE(MASK(:))                                                &
           DVOL(:,IMODE)=DVOL(:,IMODE)+RATIO1(ICP)*MD(:,IMODE,ICP)
         ENDIF
        ENDDO
!     DVOL calculates particle dry volume assuming pure H2SO4
       ELSE
        DVOL(:,IMODE)=RATIO2(IMODE)
       ENDIF
      ENDDO
!
      DO IMODE=1,NMODES
       DDPCUB(:)=SIXOVRPIX(IMODE)*DVOL(:,IMODE)
       DRYDP(:,IMODE)=CBRT(DDPCUB(:))
!!       DRYDP(:,IMODE)=DDPCUB(:)**(1.0/3.0)
      ENDDO
!
      RETURN
      END SUBROUTINE UKCA_CALC_DRYDIAM
#endif
