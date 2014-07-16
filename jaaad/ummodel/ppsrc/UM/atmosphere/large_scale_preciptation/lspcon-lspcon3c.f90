
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculates constants used in large-scale precipitation scheme.
           SUBROUTINE LSPCON(CX,CONSTP,x1i,x1ic,x1r,x2r,x4r,m_ci,   &
     &                       l_psd,ai,bi,aic,bic,                   &
     &                       lsp_ei,lsp_fi,lsp_eic,lsp_fic)

  Use cv_run_mod, Only:                                             &
    l_rp2

             IMPLICIT NONE
!
!  Description:
!     Calculates constants used within the LSP_ICE routine.
!
!  Method:
!     Calculate powers, gamma functions and constants from combinations
!     of physical parameters selected by the comdecks.
!
!  Current Code Owner:  Jonathan Wilkinson
!
!  Code description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables:







! C_LSPDRP start

      ! Microphysics parameters

      ! Drop size distribution for rain: N(D) =  N0 D^m exp(-lambda D)
      ! where N0 = X1R lambda^X2R  and m=X4R
!     REAL, PARAMETER :: X1R is set in the UMUI
!     REAL, PARAMETER :: X2R is set in the UMUI
!     REAL, PARAMETER :: X4R is set in the UMUI

      ! Drop size distribution for graupel: N(D) =  N0 D^m exp(-lambda D)
      ! where N0 = X1G lambda^X2G  and m=X4G
      REAL, PARAMETER :: X1G=5.E25
      REAL, PARAMETER :: X2G=-4.0
      REAL, PARAMETER :: X4G=2.5

      ! Particle size distribution for ice: N(D) = N0 D^m exp(-lambda D)
      ! where N0 = X1I TCG lambda^X2I, m=X4I and TCG=exp(- X3I T[deg C])
!     REAL, PARAMETER :: X1I is set in the UMUI
      REAL, PARAMETER :: X2I=0.0
      REAL, PARAMETER :: X3I=0.1222
      REAL, PARAMETER :: X4I=0.0
!     REAL, PARAMETER :: X1IC is set in the UMUI
      REAL, PARAMETER :: X2IC=0.0
      REAL, PARAMETER :: X3IC=0.1222
      REAL, PARAMETER :: X4IC=0.0

      ! Mass diameter relationship for graupel:  m(D) = AG D^BG
      REAL, PARAMETER :: AG=261.8
      REAL, PARAMETER :: BG=3.0

      ! Mass diameter relationship for ice:  m(D) = AI D^BI
      ! These are set in the UMUI (at 6.6). Values at 6.5 are listed below.
      ! Recommended values for the generic particle size distribution are
      ! from Brown and Francis and are
      ! AI=AIC=1.85E-2, BI=BIC=1.9. If l_calcfall is changed from .true.
      ! then the generic psd values should be set below to ci=cic=8.203
      ! and di=dic=0.2888
      ! REAL, PARAMETER :: AI=0.0444
      ! REAL, PARAMETER :: BI=2.1
      ! REAL, PARAMETER :: AIC=0.587
      ! REAL, PARAMETER :: BIC=2.45

      ! The area diameter relationships are only used if
      ! L_CALCFALL=.TRUE.
      ! Area diameter relationship for ice:  Area(D) = RI D^SI
      REAL, PARAMETER :: RI=0.131
      REAL, PARAMETER :: SI=1.88
      REAL, PARAMETER :: RIC=0.131
      REAL, PARAMETER :: SIC=1.88

      ! The Best/Reynolds relationships are only used if
      ! L_CALCFALL=.TRUE.
      ! Relationship between Best number and Reynolds number:
! Re(D) =LSP_EI(C) Be^LSP_FI(C)
      ! These values are set in the UMUI, but the default values are
      ! listed below. N.B. these have been renamed for VN7.3 with 
      ! 'LSP_' added from the previous versions to avoid conflicts 
      ! later in the code. 
      ! REAL, PARAMETER :: LSP_EI=0.2072  Set in the UMUI
      ! REAL, PARAMETER :: LSP_FI=0.638   Set in the UMUI
      ! REAL, PARAMETER :: LSP_EIC=0.2072 Set in the UMUI
      ! REAL, PARAMETER :: LSP_FIC=0.638  Set in the UMUI


      ! The fall speeds of ice particles are only used if
      ! L_CALCFALL=.FALSE.
      ! Fall speed diameter relationships for ice:
      ! vt(D) = CI D^DI
      REAL, PARAMETER :: CI0=14.3
      REAL, PARAMETER :: DI0=0.416
      REAL, PARAMETER :: CIC0=74.5
      REAL, PARAMETER :: DIC0=0.640

      ! Axial ratio (c-axis divided by a-axis) ESTIMATES. These are not
      ! consistent with those from the area diameter relationships.
      REAL, PARAMETER :: AR=1.0
      REAL, PARAMETER :: ARC=1.0

      ! Fall speed diameter relationship for rain: vt(D) = CR D^DR
      REAL, PARAMETER :: CR=386.8
      REAL, PARAMETER :: DR=0.67

      ! Fall speed diameter relationship for graupel: vt(D) = CG D^DG
      REAL, PARAMETER :: CG=253.0
      REAL, PARAMETER :: DG=0.734

      ! Do we wish to calculate the ice fall velocities?
      ! TRUE if calculate speeds, FALSE if specify speeds
      LOGICAL, PARAMETER :: L_CALCFALL=.TRUE.
!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

!*----------------------------------------------------------------------
      REAL RHO_WATER              ! DENSITY OF WATER (KG/M3)
!
      PARAMETER (RHO_WATER = 1000.0)
!
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------
!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

!*----------------------------------------------------------------------
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
! --------------------------COMDECK C_LSPDIF---------------------------
! input variables
!
! Values of reference variables
      REAL,PARAMETER:: AIR_DENSITY0=1.0             ! kg m-3
      REAL,PARAMETER:: AIR_VISCOSITY0=1.717E-5      ! kg m-1 s-1
      REAL,PARAMETER:: AIR_CONDUCTIVITY0=2.40E-2    ! J m-1 s-1 K-1
      REAL,PARAMETER:: AIR_DIFFUSIVITY0=2.21E-5     ! m2 s-1
      REAL,PARAMETER:: AIR_PRESSURE0=1.0E5          ! Pa

! Values of diffusional growth parameters
      ! Terms in deposition and sublimation
      REAL,PARAMETER:: APB1=(LC+LF)**2 * EPSILON /(R*AIR_CONDUCTIVITY0)
      REAL,PARAMETER:: APB2=(LC+LF) / AIR_CONDUCTIVITY0
      REAL,PARAMETER:: APB3=R/(EPSILON*AIR_PRESSURE0*AIR_DIFFUSIVITY0)
      ! Terms in evap of melting snow and rain
      REAL,PARAMETER:: APB4=LC**2*EPSILON/(R*AIR_CONDUCTIVITY0)
      REAL,PARAMETER:: APB5=LC /AIR_CONDUCTIVITY0
      REAL,PARAMETER:: APB6=R/(EPSILON*AIR_PRESSURE0*AIR_DIFFUSIVITY0)

! Values of numerical approximation to wet bulb temperature
      ! Numerical fit to wet bulb temperature
      REAL,PARAMETER:: TW1=1329.31
      REAL,PARAMETER:: TW2=0.0074615
      REAL,PARAMETER:: TW3=0.85E5
      ! Numerical fit to wet bulb temperature
      REAL,PARAMETER:: TW4=40.637
      REAL,PARAMETER:: TW5=275.0

! Ventilation parameters
      REAL,PARAMETER:: SC=0.6
! f(v)  =  VENT_ICE1 + VENT_ICE2  Sc**(1/3) * Re**(1/2)
      REAL,PARAMETER:: VENT_ICE1=0.65
      REAL,PARAMETER:: VENT_ICE2=0.44
! f(v)  =  VENT_RAIN1 + VENT_RAIN2  Sc**(1/3) * Re**(1/2)
      REAL,PARAMETER:: VENT_RAIN1=0.78
      REAL,PARAMETER:: VENT_RAIN2=0.31
!
! Subroutine arguments
!  Obtain the size of CONSTP and CX
! Start C_LSPSIZ
! Description: Include file containing idealised forcing options
! Author:      R. Forbes
!
! History:
! Version  Date      Comment
! -------  ----      -------
!   6.1    01/08/04  Increase dimension for rain/graupel.  R.Forbes
!   6.2    22/08/05  Include the step size between ice categories.
!                                                   Damian Wilson

! Sets up the size of arrays for CX and CONSTP
      REAL CX(100),CONSTP(100)
      INTEGER,PARAMETER:: ice_type_offset=20

! End C_LSPSIZ
!
      REAL, INTENT(IN)::                                                &
     &  x1i                                                             &
                  ! Intercept of aggregate size distribution
     &, x1ic                                                            &
                  ! Intercept of crystal size distribution
     &, x1r                                                             &
                  ! Intercept of raindrop size distribution
     &, x2r                                                             &
                  ! Scaling parameter of raindrop size distribution
     &, x4r                                                             &
                  ! Shape parameter of raindrop size distribution
     &, m_ci                                                            &
                  ! Used to modify ice fall speed
     &, ai, bi                                                          &
                  ! Aggregate mass-size information m(D) = ai D^bi
     &, aic, bic                                                        &
                  ! Crystal mass-size information m(D) = aic D^bic
     &, lsp_ei,  lsp_fi                                                 &
                  ! Ice aggregate Best number and Reynolds number 
                  ! relationship: Re(D) =LSP_EI Be^LSP_FI
     &, lsp_eic, lsp_fic                                                        
                  ! Ice crystal Best number and Reynolds number 
                  ! relationship: Re(D) =LSP_EIC Be^LSP_FIC


      LOGICAL, INTENT(IN):: l_psd   ! use generic particle size distrib.

!  Local scalars:
             INTEGER I
! Counter to print out the contents of CX and CONSTP
             REAL TEMP                                                  &
! Forms input to the GAMMAF routine which calculates gamma functions.
     &,           G1, G2, G3                                            &
     &,           GB1, GB2, GB3                                         &
     &,           GBC1, GBD1, GBDC1                                     &
     &,           GC1, GC2, GC3                                         &
     &,           GD3, GDC3, GD52, GDC52                                &
     &,           GDR3, GDR4, GDR52                                     &
     &,           GR2, GR4, GR5, GR6                                    &
     &,           GG1, GGB1, GGBD1                                      &
     &,           GG2, GDG52, GDG3, GG3                                 &
! Represents the gamma function of BI+DI+1 etc.
! Fall speed of ice particles parameters
     &,           CI,DI,CIC,DIC
!
! Function and subroutine calls:
             EXTERNAL GAMMAF
!
!- End of header
!
! Do we need to calculate fall speeds?
             IF (L_CALCFALL) THEN
! Define fall speeds
               CI=LSP_EI*AIR_VISCOSITY0**(1.0-2.0*LSP_FI)                   &
     &            *AIR_DENSITY0**(LSP_FI-1.0)                               &
     &            *(2.0*G)**LSP_FI*(AI/RI)**LSP_FI
               DI=LSP_FI*(BI+2.0-SI)-1.0
               CIC=LSP_EIC*AIR_VISCOSITY0**(1.0-2.0*LSP_FIC)                &
     &            *AIR_DENSITY0**(LSP_FIC-1.0)                              &
     &            *(2.0*G)**LSP_FIC*(AIC/RIC)**LSP_FIC
               DIC=LSP_FIC*(BIC+2.0-SIC)-1.0
! Modify fallspeeds for random parameters 2           
               IF (L_RP2) THEN 
                   CI=CI*M_CI
                   CIC=CIC*M_CI
               ENDIF       
             ELSE
! Use preset parameters
               CI=CI0
               DI=DI0
               CIC=CIC0
               DIC=DIC0
             ENDIF  ! Calculation of fall speeds
!
! CX values. 1-20 are for the crystal population. 21-40 are for the
! aggregate population. 41-60 are for rain. 61-80 are for graupel.
! 81-99 are for the generic ice size distribution.
!
! Crystals
             CX(1)=(BIC+1.0+X4IC-X2IC)/BIC
             CX(2)=-(X4IC+1.0-X2IC)/BIC
             CX(3)=DIC/(BIC+1.0+X4IC-X2IC)
             CX(4)=(2.0+X4IC-X2IC)/(BIC+1.0+X4IC-X2IC)
             CX(5)=(5.0+DIC+2.0*X4IC-2.0*X2IC)*0.5/(BIC+1.0+X4IC-X2IC)
             CX(6)=(3.0+DIC+X4IC-X2IC)/(BIC+1.0+X4IC-X2IC)
             CX(7)=1.0/(X2IC-X4IC-1.0-BIC)
             CX(8)=1.0+X4IC
             CX(9)=2.0+X4IC
             CX(10)=3.0+X4IC
             CX(11)=X2IC
             CX(12)=X3IC
             CX(13)=1.0+X4IC+BIC
             CX(14)=BIC

! Aggregates
             CX(23)=DI/(BI+1.0+X4I-X2I)
             CX(24)=(2.0+X4I-X2I)/(BI+1.0+X4I-X2I)
             CX(25)=(5.0+DI+2.0*X4I-2.0*X2I)*0.5/(BI+1.0+X4I-X2I)
             CX(26)=(3.0+DI+X4I-X2I)/(BI+1.0+X4I-X2I)
             CX(27)=1.0/(X2I-X4I-1.0-BI)
             CX(28)=1.0+X4I
             CX(29)=2.0+X4I
             CX(30)=3.0+X4I
             CX(31)=X2I
             CX(32)=X3I
             CX(33)=3.0+X4I+BI
             CX(34)=2.0+X4I+BI
             CX(35)=1.0+X4I+BI
! Rain
             CX(41)=DR/(4.0+DR-X2R+X4R)
             CX(42)=1.0/(4.0+DR-X2R+X4R)
             CX(43)=6.0+X4R
             CX(44)=5.0+X4R
             CX(45)=4.0+X4R
             CX(46)=X2R
             CX(47)=2.0+X4R-X2R
             CX(48)=1.0/(4.0+X4R+DR-X2R)
             CX(49)=(DR+5.0)*0.5-X2R+X4R
             CX(50)=(3.0+DR-X2R+X4R)/(4.0+DR-X2R+X4R)
! Rain mixing ratio
             CX(51)=DR/(4.0-X2R+X4R)
             CX(52)=1.0/(4.0-X2R+X4R)
             CX(53)=(3.0+DR-X2R+X4R)/(4.0-X2R+X4R)
!
! Graupel

             CX(63)=DG/(BG+1.0+X4G-X2G)
             CX(64)=(2.0+X4G-X2G)/(BG+1.0+X4G-X2G)
             CX(65)=(5.0+DG+2.0*X4G-2.0*X2G)*0.5/(BG+1.0+X4G-X2G)
             CX(66)=(3.0+DG+X4G-X2G)/(BG+1.0+X4G-X2G)
             CX(67)=1.0/(X2G-X4G-1.0-BG)
             CX(68)=1.0+X4G
             CX(69)=2.0+X4G
             CX(70)=3.0+X4G
             CX(71)=X2G

             CX(73)=3.0+X4G+BG
             CX(74)=2.0+X4G+BG
             CX(75)=1.0+X4G+BG

!
! Generic ice particle size distribution
             CX(81)=2.0+di 
             CX(82)=di+bi 
             CX(84)=1.0 
             CX(85)=1.0+0.5*(di+1.0)

! Define gamma values
! Crystals
             TEMP=1.0+X4IC
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GC1)
             TEMP=BIC+1.0+X4IC
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GBC1)
             TEMP=BIC+DIC+1.0+X4IC
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GBDC1)
             TEMP=2.0+X4IC
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GC2)
             TEMP=(DIC+5.0+2.0*X4IC)*0.5
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GDC52)
             TEMP=DIC+3.0+X4IC
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GDC3)
             TEMP=3.0+X4IC
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GC3)
! Aggregates
             TEMP=1.0+X4I
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,G1)
             TEMP=BI+1.0+X4I
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GB1)
             TEMP=BI+2.0+X4I
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GB2)
             TEMP=BI+3.0+X4I
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GB3)
             TEMP=BI+DI+1.0+X4I
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GBD1)
             TEMP=2.0+X4I
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,G2)
             TEMP=(DI+5.0+2.0*X4I)*0.5
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GD52)
             TEMP=DI+3.0+X4I
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GD3)
             TEMP=3.0+X4I
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,G3)
! Rain
             TEMP=DR+4.0+X4R
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GDR4)
             TEMP=4.0+X4R
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GR4)
             TEMP=6.0+X4R
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GR6)
             TEMP=5.0+X4R
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GR5)
             TEMP=2.0+X4R
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GR2)
             TEMP=(DR+5.0+2.0*X4R)*0.5
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GDR52)
             TEMP=DR+3.0+X4R
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GDR3)
! Graupel
             TEMP=1.0+X4G
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GG1)
             TEMP=BG+1.0+X4G
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GGB1)
             TEMP=BG+DG+1.0+X4G
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GGBD1)
             TEMP=2.0+X4G
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GG2)
             TEMP=(DG+5.0+2.0*X4G)*0.5
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GDG52)
             TEMP=DG+3.0+X4G
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GDG3)
             TEMP=3.0+X4G
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GG3)

!
! CONSTP values. 1-20 are for the crystal population. 21-40 are for the
! aggregate population. 41-60 are for rain. 61-80 are for graupel. 
! 81-99 are for the generic ice size distribution
!
! Crystals
             CONSTP(1)=X1IC
             CONSTP(2)=1.0/GC1
             CONSTP(3)=1.0/(AIC*GBC1)
             CONSTP(4)=CIC*GBDC1/GBC1
             CONSTP(5)=1.0/(AIC*X1IC*GBC1)
             CONSTP(6)=2.0*PI*X1IC
             CONSTP(7)=VENT_ICE1*GC2
             CONSTP(8)=VENT_ICE2*SC**(1.0/3.0)/AIR_VISCOSITY0**0.5      &
     &                 *SQRT(CIC)*GDC52
             CONSTP(9)=PI*0.25*X1IC*CIC*GDC3
             CONSTP(10)=5.0*GC1
             CONSTP(11)=2.0*GC2
             CONSTP(12)=0.25*GC3
             CONSTP(13)=PI**2*RHO_WATER*X1IC*X1R
             CONSTP(14)=2.0*PI*AIR_CONDUCTIVITY0/LF*X1IC
! Capacitance relative to spheres of same maximum dimension
! Formula depends on value of axial ratio
             CONSTP(15)=ARC
             IF (ARC  >   1.0) THEN
! Prolate
               CONSTP(15)=(1.0-(1.0/CONSTP(15))**2)**0.5 /              &
     &                    ALOG( CONSTP(15) + (CONSTP(15)**2-1.0)**0.5 )
             ELSEIF (ARC  ==  1.0) THEN
! Spherical
               CONSTP(15)=1.0
             ELSE
! Oblate
               CONSTP(15)=(1.0-CONSTP(15)**2)**0.5                      &
     &                    /ASIN((1.0-CONSTP(15)**2)**0.5)
             ENDIF
! Now adjust diffusional growth constants for capacitance
             CONSTP(6)=CONSTP(6)*CONSTP(15)
             CONSTP(14)=CONSTP(14)*CONSTP(15)
!            constp(20) is reserved for lsp_collection
!
! Values 16 to 23 are unused
! Aggregates
             CONSTP(24)=CI*GBD1/GB1
             CONSTP(25)=1.0/(AI*X1I*GB1)
             CONSTP(26)=2.0*PI*X1I
             CONSTP(27)=VENT_ICE1*G2
             CONSTP(28)=VENT_ICE2*SC**(1.0/3.0)/AIR_VISCOSITY0**0.5     &
     &                  *SQRT(CI)*GD52
             CONSTP(29)=PI*0.25*X1I*CI*GD3
             CONSTP(30)=5.0*G1
             CONSTP(31)=2.0*G2
             CONSTP(32)=0.25*G3
             CONSTP(33)=PI**2*RHO_WATER*X1I*X1R
             CONSTP(34)=2.0*PI*AIR_CONDUCTIVITY0/LF*X1I
! Capacitance relative to spheres of same maximum dimension
! Formula depends on value of axial ratio
             CONSTP(35)=AR
             IF (AR  >   1.0) THEN
! Prolate
               CONSTP(35)=(1.0-(1.0/CONSTP(35))**2)**0.5 /              &
     &                    ALOG( CONSTP(35) + (CONSTP(35)**2-1.0)**0.5 )
             ELSEIF (AR  ==  1.0) THEN
! Spherical
               CONSTP(35)=1.0
             ELSE
! Oblate
               CONSTP(35)=(1.0-CONSTP(35)**2)**0.5                      &
     &                    / ASIN((1.0-CONSTP(35)**2)**0.5)
             ENDIF
! Now adjust diffusional growth constants for capacitance
             CONSTP(26)=CONSTP(26)*CONSTP(35)
             CONSTP(34)=CONSTP(34)*CONSTP(35)
!
             CONSTP(36) = GC1*GB3
             CONSTP(37) = 2.0*GC2*GB2
             CONSTP(38) = GC3*GB1
             CONSTP(39) = AI*X1IC*X1I*PI/4.0
!            constp(40) is reserved for lsp_collection
! Rain
             CONSTP(41)=6.0*CR*GDR4/GR4
             CONSTP(42)=PI*RHO_WATER/6.0*X1R*GDR4*CR
             CONSTP(43)=1.0/120.0*GR6
             CONSTP(44)=1.0/24.0*GR5
             CONSTP(45)=1.0/6.0*GR4
             CONSTP(46)=2.0*PI*X1R
             CONSTP(47)=VENT_RAIN1*GR2
             CONSTP(48)=VENT_RAIN2*SC**(1.0/3.0)/AIR_VISCOSITY0**0.5    &
     &                  *GDR52*SQRT(CR)
             CONSTP(49)=PI*0.25*X1R*CR*GDR3
             CONSTP(50)=1.0/(PI*RHO_WATER*X1R*GR4/6.0)
! Values 51 to 60 are unused
!
! Graupel
             CONSTP(64) = CG*GGBD1/GGB1
             CONSTP(65) = 1.0/(AG*X1G*GGB1)
             CONSTP(67) = VENT_RAIN1*GG2
             CONSTP(68) = VENT_RAIN2*SC**(1.0/3.0)/AIR_VISCOSITY0**0.5  &
     &                    *SQRT(CG)*GDG52
             CONSTP(69) = PI*0.25*X1G*CG*GDG3
             CONSTP(74) = 2.0*PI*AIR_CONDUCTIVITY0*X1G/LF
             CONSTP(76) = GG1*GB3
             CONSTP(77) = 2.0*GG2*GB2
             CONSTP(78) = GG3*GB1
             CONSTP(79) = AI*X1G*X1I*PI/4.0
             CONSTP(80) = pi*pi/24.0*x1g*ag
!
! Generic ice particle size distribution 
             CONSTP(81)=pi*0.25*ci 
             CONSTP(82)=ci*ai 
             CONSTP(83)=2.0*pi*constp(35)
             CONSTP(84)=vent_ice1 
             CONSTP(85)=vent_ice2*sc**(1.0/3.0)/air_viscosity0**0.5     & 
     &                  *sqrt(ci)
             CONSTP(86)=pi*pi/24.0*rho_water*x1r
             CONSTP(87)=gr6
             CONSTP(88)=2.0*gr5
             CONSTP(89)=gr4
             CONSTP(90)=2.0*pi*constp(35)*air_conductivity0/Lf
             CONSTP(91)=gb3
             CONSTP(92)=2.0*gb2
             CONSTP(93)=gb1
! End the subroutine
             RETURN
           END SUBROUTINE LSPCON
!
!  SUBROUTINE GAMMAF--------------------------------------------------
!   PURPOSE: CALCULATES COMPLETE GAMMAF FUNCTION BY
!   A POLYNOMIAL APPROXIMATION
! --------------------------------------------------------------------
!
