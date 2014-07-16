#if defined(A04_3B)
! ------------------COMDECK C_LSPDRP-----------------------------------
!
! COX/GOLDING values - 1.5 x Heymsfield fallspeed
!
! Drop size distribution for rain: N(D) = N0 exp(-lambda D)
! where N0 = X1R lambda^X2R
!     REAL,PARAMETER:: X1R is set in the UMUI
!     REAL,PARAMETER:: X2R is set in the UMUI

! Particle size distribution for ice: N(D) = N0 D^m exp(-lambda D)
! where N0 = X1I exp( -X3I T[deg C]/8.18 ) lambda^X2I   and   m=X4I
!     Ice distn coeff
!     REAL,PARAMETER:: X1I is set in the UMUI
      REAL,PARAMETER:: X2I=0.0
!     Whether there is a temperature dependence to N0snow
      REAL,PARAMETER:: X3I=1.0
!     m parameter in gamma distribution
      REAL,PARAMETER:: X4I=0.0

! Mass diameter relationship for ice:  m(D) = AI D^BI
! These are set in the UMUI (at 6.6). Values at 6.5 are listed below.
!      REAL,PARAMETER:: AI=6.9E-2
!      REAL,PARAMETER:: BI=2.0

! Fall speed diameter relationship for ice: vt(D) = CI D^DI
! As CI is in common now initialised via blockdata in the 
! file c_lspdrpdt.h
       REAL :: CI
       COMMON/STPH_RP_2/CI
      REAL,PARAMETER:: DI=.527

! Fall speed diameter relationship for rain: vt(D) = CR D^DR
      REAL,PARAMETER:: CR=386.8
      REAL,PARAMETER:: DR=0.67

! C_LSPDRP end
#endif
#if defined(A04_3C) || defined(A04_3D)
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
#endif
