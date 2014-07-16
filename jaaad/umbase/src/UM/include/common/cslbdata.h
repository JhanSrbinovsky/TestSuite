! CSLBDATA start
! Description:
!   Constants passed via input namelist SLABLIST
!
! Current Code Owner: K.Williams
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.3   23/10/01   Added EVP variables
!   6.0   26/08/03   Added McPhee variables  (M. Crucifix)
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
      REAL :: DZ1       ! Slab ocean depth
      REAL :: AICEMIN   ! minimum ice concentration
      REAL :: AMXNORTH  ! maximum ice concentration (n hemisphere)
      REAL :: AMXSOUTH  ! maximum ice concentration (s hemisphere)
      REAL :: HICEMIN   ! minimum gridbox average ice depth
      REAL :: H0         ! minimum local ice depth

      LOGICAL :: CALIB      ! switch for calibration or corrected run

      ! Diffusion coeff. for turbulent ocean-ice heat flux.
      REAL :: eddydiff

      REAL :: epsilon    ! Minimum grid box mean ice depth (metres).
      REAL :: Ah_ice     ! Diffusion coeff. for diffusion of ice depth.
      REAL :: HCLIMIT    ! Limit for redistributing heat convergence.
      REAL :: AICEMIZFRY  ! Conc below which O2I is fixed wrt ice
      REAL :: MCPHEE_COEFF     !  Ocean-ice heat transfer coeff
      REAL :: MCPHEE_MINUSTAR  !  Min friction vel for McPhee scheme
      REAL :: DTE            ! Subcycling timestep (EVP)
      REAL :: CW             ! Cosine of turn. angle in water (EVP)
      REAL :: SW             ! Sine of turn. angle in water (EVP)
      REAL :: ECC2           ! Squared eccentricity of yield curve (EVP)
      REAL :: EYC            ! Coeff. for calc of parameter E (EVP)
      REAL :: CSTAR          ! c*, parameter in the defn of prss (EVP)
      REAL :: PSTAR          ! p*, parameter in the defn of prss (EVP)
      REAL :: QUAD_ICE_DRAG  ! Coefft of quadratic ice-ocean drag (EVP)


      COMMON /SLBLST/                                                   &
     &  DZ1,AICEMIN,AMXNORTH,AMXSOUTH,HICEMIN,H0,CALIB,                 &
     &  eddydiff,Ah_ice,HCLIMIT,AICEMIZFRY,                             &
     &  MCPHEE_COEFF,MCPHEE_MINUSTAR,                                   &
     &  DTE,CW,SW,ECC2,EYC,CSTAR,PSTAR,QUAD_ICE_DRAG
      NAMELIST/SLABLIST/DZ1,AICEMIN,AMXNORTH,                           &
     &  AMXSOUTH,HICEMIN,H0,CALIB,                                      &
     &  eddydiff,Ah_ice,HCLIMIT,AICEMIZFRY,                             &
     &  MCPHEE_COEFF,MCPHEE_MINUSTAR,                                   &
     &  DTE,CW,SW,ECC2,EYC,CSTAR,PSTAR,QUAD_ICE_DRAG
! CSLBDATA end
