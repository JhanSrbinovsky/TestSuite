
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE EX_COEF------------------------------------------------
!!!
!!!  Purpose: To calculate exchange coefficients for boundary layer
!!!           subroutine KMKH.
!!!
!!!  Programming standard: Unified Model Documentation Paper No 3
!!!
!!!  Documentation: UMDP No.24
!!!
!!!---------------------------------------------------------------------

!!  Arguments :-
      SUBROUTINE EX_COEF (                                              &
     & row_length,rows,off_x,off_y,BL_LEVELS                            &
     &,K_LOG_LAYR,LAMBDA_MIN,L_LAMBDAM2,L_FULL_LAMBDAS,lq_mix_bl        &
     &,LOCAL_FA,Prandtl,SIGMA_H,CCB,CCT,NTML_LOCAL,ISHEAR_BL            &
     &,L_SBLeq,L_SBLco,Muw_SBL,Mwt_SBL                                  &
     &,SBL_OP,FLANDG                                                    &
     &,CCA,DBDZ,DVDZM                                                   &
     &,RHO,ZH_LOCAL,Z_UV,Z_TQ,Z0M,H_BLEND                               &
     &,CUMULUS,NTPAR,NTML_NL                                            &
     &,U_P,V_P,V_S,FB_SURF,QW,TL                                        &
     &,FM_3D, FH_3D, L_subfilter_vert,L_subfilter_horiz                 &
     &,RHOKM,RHOKH                                                      &
! Variables for STPH_RP
     &,G0_RP,par_mezcla                                                 &
     &,LTIMER                                                           &
     &,BL_diag                                                          &
     &,nSCMDpkgs,L_SCMDiags                                             &
     & )
!
      USE BL_OPTION_MOD, ONLY : WeightLouisToLong

      Use cv_run_mod, Only:                                             &
          l_mom, l_rp,l_rp2

      Use bl_diags_mod, Only :                                          &
          strnewbldiag

      IMPLICIT NONE

      LOGICAL LTIMER
      LOGICAL                                                           &
     & L_SBLeq                                                          &
                   ! IN Switch for Equilibrium SBL model
     &,L_SBLco                                                          &
                   ! IN Switch for coupled gradient method in
!                  !    Equilibrium SBL model
     &,L_LAMBDAM2                                                       &
                   ! IN LambdaM=2*LambdaH (operational setting).
!                  !    Could use scaling factor for flexibility.
     &,L_FULL_LAMBDAS                                                   &
!                  ! IN Lambdas NOT reduced above NTML_LOCAL+1
     &,lq_mix_bl   ! IN True if mixing ratios used

      INTEGER                                                           &
     & row_length,rows,off_x,off_y                                      &
     &,BL_LEVELS                                                        &
                   ! IN maximum number of boundary layer levels
     &,ISHEAR_BL                                                        &
                   ! IN Switch for shear-dominated bl code
     &,K_LOG_LAYR                                                       &
                   ! IN num of levs requiring log-profile correction
     &,LOCAL_FA                                                         &
                   ! switch for free atmospheric mixing options
     &,Prandtl
                   ! switch for Prandtl number options

      INTEGER                                                           &
     & CCB(row_length,rows)                                             &
                                ! IN  Convective Cloud Base.
     &,CCT(row_length,rows)     ! IN  Convective Cloud Top.

      INTEGER                                                           &
     & NTML_NL(row_length,rows)                                         &
                                ! IN Number of model layers in the
!                                    turbulently mixed layer as
!                                    determined from the non-local
!                                    scheme.
     &,NTPAR(row_length,rows)   ! IN Top level of parcel ascent
!     Declaration of new BL diagnostics.
      Type (Strnewbldiag) :: BL_diag


! Start blopt8a

! Description:
!   Permissible settings for BL options.
!
! Current Code Owner: J. M. Edwards
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 6.2      27/01/06 Original code.  J. M. Edwards
!
      INTEGER, PARAMETER :: Off = 0  ! Switch disabled
      INTEGER, PARAMETER :: On  = 1  ! Switch enabled
!
      INTEGER, PARAMETER :: BrownGrant97 = 1
      INTEGER, PARAMETER :: BrownGrant97_limited = 2
!     Options for non-gradient stress following
!     Brown and Grant (1997), version 2 including a limit on its size
!
!     Options for flux gradient formulation
      INTEGER, PARAMETER :: Locketal2000   = 0
!       Flux gradients as in Lock et al. (2000)
      INTEGER, PARAMETER :: HoltBov1993 = 1
!       Flux gradients as in Lock et al (2000) but using
!       coefficients from Holtslag and Boville (1993)
      INTEGER, PARAMETER :: LockWhelan2006 = 2
!       Flux gradients as in Lock and Whelan (2006)
!
!     Options for form drag
      INTEGER, PARAMETER :: No_drag         = 0
      INTEGER, PARAMETER :: Effective_z0    = 1
      INTEGER, PARAMETER :: Explicit_stress = 2
!
!     Options for marine boundary layers
      INTEGER, PARAMETER :: Fixed_Z0T = 0
!       Stanard flixed value of thermal roughness length over sea
      INTEGER, PARAMETER :: SurfDivZ0T = 1
!       Thermal roughness length over sea defined from surface
!       divergence theory
      INTEGER, PARAMETER :: DynDiag_ZL = 1
!       The ratio of the height of the inversion to the surface
!       Obukhov length is used a dynamic criterion in the
!       diagnosis of BL types
!
      INTEGER, PARAMETER :: Use_Correct_Ustar = 2
!       Option under the COR_MO_ITER switch for the dust scheme
!       to use the correct ustar
!
!     Options for stable boundary layers
      INTEGER, PARAMETER ::  Long_tails           = 0
      INTEGER, PARAMETER ::  Sharpest             = 1
      INTEGER, PARAMETER ::  Sharp_sea_long_land  = 2
      INTEGER, PARAMETER ::  Mes_tails            = 3
      INTEGER, PARAMETER ::  Louis_tails          = 4
      INTEGER, PARAMETER ::  Depth_based          = 5
      INTEGER, PARAMETER ::  Sharp_sea_mes_land   = 6
      INTEGER, PARAMETER ::  LEM_stability        = 7
      INTEGER, PARAMETER ::  Sharp_sea_Louis_land = 8

!     Options for Prandtl number (in local Ri scheme)
      INTEGER, PARAMETER ::  Constant_SBL = 0
      INTEGER, PARAMETER ::  LockMailhot2004 = 1

! End blopt8a
      REAL                                                              &
     & LAMBDA_MIN                                                       &
                                ! IN Min value of length scale LAMBDA.
     &,CCA(row_length,rows)                                             &
                                ! IN Convective Cloud Amount.
     &,SIGMA_H(row_length,rows)                                         &
                                ! IN Standard deviation of subgrid 
                                !    orography (m) [= 2root2*ho2r2_orog]
     &,RHO(row_length,rows,BL_LEVELS)                                   &
                                ! IN density on theta levels;
!                               !    used in RHOKM so wet density
     &,DBDZ(row_length,rows,2:BL_LEVELS)                                &
                                ! IN Buoyancy gradient across lower
!                                     interface of layer.
     &,U_P(row_length,rows,BL_LEVELS)                                   &
                                      ! IN Westerly wind component
!                                     horizontally interpolated to
!                                     P-grid. (m/s)
     &,V_P(row_length,rows,BL_LEVELS)                                   &
                                      ! IN Southerly wind component
!                                     horizontally interpolated to
!                                     P-grid. (m/s)
     &,QW(row_length,rows,BL_LEVELS)                                    &
                                     ! IN Total water content
!                                      (kg per kg air).
     &,TL(row_length,rows,BL_LEVELS)                                    &
                                     ! IN Liquid/frozen water
!                                      temperature (K).
     &,Z_UV(row_length,rows,BL_LEVELS)                                  &
                                ! IN Z_UV(K) is height of u level k
     &,Z_TQ(row_length,rows,BL_LEVELS)                                  &
                                ! IN Z_TQ(K) is height of T,Q level k
!                               !    NOTE: RI(K) is held at Z_TQ(K-1)
     &,Z0M(row_length,rows)                                             &
                                ! IN Roughness length for momentum (m).
     &,H_BLEND(row_length,rows)                                         &
                                ! IN Blending height for effective
!                                     roughness length scheme
     &,DVDZM(row_length,rows,2:BL_LEVELS)                               &
                                         ! IN Modulus of wind shear.
     &,Muw_SBL,Mwt_SBL          ! IN Powers to use in prescription of
!                               !    equilibrium profiles of stress and
!                               !    buoyancy flux in Equilib. SBL model
      INTEGER                                                           &
     & SBL_OP                   ! IN stable boundary layer option

      REAL                                                              &
     & V_S(row_length,rows)                                             &
                                         ! IN Surface friction velocity
!                                          (m/s)
     &,FB_SURF(row_length,rows)                                         &
                                         ! IN Surface buoyancy flux over
!                                          density (m^2/s^3).
     &,FUNC(row_length,rows)                                            &
                                         ! 2D variable for SBL stabiliy
!                                        ! function options
     &,SHARP(row_length,rows)                                           &
                                         ! 2D variable for SHARP 
!                                        ! stabiliy function
     &, FM_3D(row_length,rows,BL_LEVELS)                                &
!             ! stability function for momentum transport.
!             ! level 1 value is dummy for use in diagnostics
     &, FH_3D(row_length,rows,BL_LEVELS)
!             ! stability function for heat and moisture.
!             ! level 1 value is dummy for use in diagnostics
      LOGICAL                                                           &
     & CUMULUS(row_length,rows)                                         &
                                ! INOUT Flag for boundary layer cumulus.
!                                     Can only be changed if ISHEAR_BL=1
     &, L_subfilter_vert                                                &
                            ! subgrid turbulence scheme in vertical     
     &, L_subfilter_horiz   
                            ! subgrid turbulence scheme in horizontal

      REAL                                                              &
     & RHOKM(1-off_x:row_length+off_x,1-off_y:rows+off_y,BL_LEVELS)     &
                                ! INOUT Layer K-1 - to - layer K
!                                        exchange coefficient for
!                                        momentum, on UV-grid with first
!                                        and last rows set to "missing
!                                        data"
     &,RHOKH(row_length,rows,BL_LEVELS)                                 &
                                ! INOUT Layer K-1 - to - layer K
!                               !       exchange coefficient for FTL.
!                               ! On OUT: still to be multiplied by rho 
!                               !    (if lq_mix_bl) and, for Ri-based 
!                               !    scheme, interpolated to rho
!                               !    levels in BDY_EXPL2
     &,ZH_LOCAL(row_length,rows)! INOUT Mixing layer height (m).

! Definition of variables for STPH_RP: 'Stochastic physics RP'
      REAL,INTENT(InOut) :: par_mezcla ! Used to modify the neutral
                                       ! mixing length
      REAL,INTENT(InOut) :: G0_RP ! Used to vary stability functions
      INTEGER                                                           &
     & NTML_LOCAL(row_length,rows)! OUT Number of model layers in the
!                                    turbulently mixed layer as
!                                    determined from the local
!                                    Richardson number profile.
!
! Additional variables for SCM diagnostics which are dummy in full UM
      INTEGER                                                           &
     &  nSCMDpkgs             ! No of SCM diagnostics packages

      LOGICAL                                                           &
     &  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages
!
!-----------------------------------------------------------------------
      REAL                                                              &
     &  FLANDG(ROW_LENGTH,ROWS)       !  IN Land fraction on all tiles.


      EXTERNAL SBLequil
      EXTERNAL TIMER

!-----------------------------------------------------------------------

!    Local and other symbolic constants :-
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
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
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

!*----------------------------------------------------------------------
! C_VKMAN start
      REAL,PARAMETER:: VKMAN=0.4 ! Von Karman's constant
! C_VKMAN end

      Character*(*), Parameter ::  RoutineName = 'ex_coef'

      REAL GRCP
      REAL EH,EM,G0,DH,DM,A_LAMBDA
      PARAMETER (                                                       &
     & EH=25.0                                                          &
                                ! Used in calc of stability function FH.
     &,EM=4.0                                                           &
                                ! Used in calc of stability function FM.
     &,A_LAMBDA=2.0                                                     &
                                ! used in calc of LAMBDA_EFF
     &,GRCP=G/CP                                                        &
                                ! Adiabatic lapse rate.
     &)

!  Equilibrium SBL model constants
      REAL    RtestMin
      INTEGER GN,NGz,kMINh
      PARAMETER (                                                       &
     & RtestMin=0.0                                                     &
                                ! Threshold for Rtest
     &,GN=19                                                            &
                                ! Size of "G"-tables (No. HonL values)
     &,NGz=90                                                           &
                                ! No. z/h steps in each "G" integration
     &,kMINh=2                                                          &
                                ! Level of minimum SBL height (>=2)
     &)

!  Define local storage.

!  Arrays.

      REAL                                                              &
     & RI(row_length,rows,2:BL_LEVELS)                                  &
                                           ! Local Richardson number.
     &,invMOsurf(row_length,rows)                                       &
                                           ! Inverse of sfce M-O length
!                    ! Note: Inverse is used so that neutral conditions
!                    !       can be handled (M-O length --> infinity)
     &,ZH_ESBL(row_length,rows)                                         &
                                           ! Ht of equilib SBL (sub-grid
     &,HLtab(GN)                                                        &
                                          ! Lookup tables (Gx calcs
     &,GHsav(GN,NGz),GMsav(GN,NGz)                                      &
                                          ! in equilib SBL scheme)
     &,THv_TQ(row_length,rows,BL_LEVELS)                                &
                                          ! Virtual potential
!                                         ! temperature on theta levels
     &,THv(row_length,rows,BL_LEVELS)                                   &
                                          ! THv_TQ interpolated to
!                                         ! U,V levels
     &,RIOUT(row_length,rows,BL_LEVELS)                                 &
                                          ! RI on theta levels for
                                          ! SCM diagnostic output
     &,PRANDTL_NUMBER(row_length,rows)                                  &
!                                         ! = KM/KH
     &,BL_weight(row_length,rows)
!                                         ! Fractional weight applied to
!                                         ! BL function, vs free atmos

      INTEGER                                                           &
     & NTML_ESBL(row_length,rows)         ! No. UV-levels inside
!                                         ! equilibrium SBL

      LOGICAL                                                           &
     & TOPBL(row_length,rows)             ! Flag for having reached
!                                    the top of the turbulently mixed
!                                    layer.
! Variables for stability function tails
      REAL                                                              &
     & FM_LOUIS                                                         &
!                                         ! FM calculated using Louis
     &,FM_SHARPEST                                                      &
!                                         ! FM calculated using SHARPEST
     &,Z_SCALE                                                          &
!                                         ! Scale height for interpoln
!                                         ! between stability functions
     &,ZPR
!                                         ! z/sigma_h

! Variables for boundary layer depth based formulation
      REAL                                                              &
     &  h_tkeb(row_length,rows)                                         &
                                     ! TKE budget based BL depth
     &, MOsurf(row_length,rows)                                         &
                                     ! surface Obukhov length
     &, diff_min(row_length,rows)

      REAL                                                              &
     & h_est                                                            &
     &,rifb                                                             &
                    ! Bulk flux Richardson number
     &,pr_n                                                             &
                    ! Neutral Prandtl number
     &,m_tau,m_buoy                                                     &
                    ! Indices for implied stress and buoyancy flux profs
     &,ind,diff

      REAL                                                              &
     & PRN, SUBA, SUBB, SUBC, SUBG, SUBH, RIC, RICINV, RIFAC
                    !Constants for LEM stability functions

!  Scalars.

      REAL                                                              &
     & ELH                                                              &
                 ! Mixing length for heat & moisture at lower layer bdy.
     &,ELM                                                              &
                 ! Mixing length for momentum at lower layer boundary.
     &,F_LOG                                                            &
                 ! Temporary in calculation of logarithmic correction
     &,FH                                                               &
                 ! (Value of) stability function for heat & moisture.
     &,FM                                                               &
                 ! (Value of) stability function for momentum transport.
     &,RTMRI                                                            &
                 ! Temporary in stability function calculation.
     &,VKZ                                                              &
                 ! Temporary in calculation of ELH.
     &,LAMBDAM                                                          &
                 ! Asymptotic mixing length for turbulent transport
!                  of momentum.
     &,LAMBDAH                                                          &
                 ! Asymptotic mixing length for turbulent transport
!                  of heat/moisture.
     &,LAMBDA_EFF! Effective mixing length used with effective
!                  roughness length scheme.

!     !Equilibrium SBL model temporary real scalar variables
      REAL                                                              &
     & Ztop,Zbot,zz,Zprev,DZ,HH                                         &
                                       ! Height variables
     &,U1,U2,Ubot,Ujmp                                                  &
                                       ! Velocity variables
     &,T1,T2,Tbot,Tjmp                                                  &
                                       ! Temperature variables
     &,Rib,Rilim,Rtest                                                  &
                                       ! Ri variables
     &,UWsfce, WTsfce, invMOsfce                                        &
                                       ! Surface flux variables
     &,USequil,WTequil,invMOequil                                       &
                                       ! Prescribed flux profiles
     &, KM, KH, KUU, KTT, KUT, KTU                                      &
                                       ! Full/partial eddy diffusivities
     &,PKM,PKH,PKUU,PKTT,PKUT,PKTU                                      &
                                       ! Scaled arguments of SBLequil
     &,Zhat,lhat,PHIe,PHIw,RiSBL                                        &
                                       ! Scaled arguments of SBLequil
     &,GM,GH,Gs,Gdz,Gz,G1,G2,H1,H2                                      &
                                       ! Temporaries for Gx calcs
     &,GHtab1,GHtab2,GMtab1,GMtab2                                      &
                                       ! Temporaries for Gx calcs
     &,CE                                                               &
                                       ! Constant in dissipation param.
!                                      ! (output by SBLequil)
     &,rpow,CB,CN                                                       &
                                       ! Constants in mixing length eqn
!                                      ! (output by SBLequil)
     &,tmp,HonL,DBDU,GonTv,ll,Mlamb                                     &
                                       ! Miscellaneous temporaries
     &,slope,DZ100,Zupper,fG0          ! Miscellaneous temporaries

      INTEGER                                                           &
     & I,j                                                              &
                   ! Loop counter (horizontal field index).
     &,K                                                                &
                 ! Loop counter (vertical level index).
     &,KM1                                                              &
                 ! K-1.
     &,MBL       ! Maximum number of model layers below mixed layer top.

!     !Equilibrium SBL model temporary integer scalar variables
      INTEGER                                                           &
     & kZtop,kZbot,GK,kG0                                               &
                          ! Temporary loop counters
     &,iERRSBL        ! SBLequil error status

!    !Equilibrium SBL model logical variables
      LOGICAL                                                           &
     & GcalcDone                                                        &
                      ! Calculation of Gx values has been performed
     &,PrevMin                                                          &
                      ! Previous Utest minimium has been detected
     &,subgrid        ! Will perform subgrid SBL depth calculation

!    !Switch to enable subgrid SBL depth diagnosis
      LOGICAL    sg_enabled
      PARAMETER (sg_enabled=.true.)

!    !Equilibrium SBL model SAVED variables
      SAVE HLtab,GHsav,GMsav,GcalcDone

!    !Equilibrium SBL model DATA statements
      DATA HLtab /0.0001,0.001,0.002,0.005,0.01,0.02,0.05,              &
     &            0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0,               &
     &            100.0,200.0,500.0/
      DATA GcalcDone /.FALSE./

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('EX_COEF ',103)
      ENDIF
!-----------------------------------------------------------------------
!   IF stochastic physics random parameters is used set the parameter
!   used to vary the stability function to a perturbed value, if not
!   use the standard setting.
!-----------------------------------------------------------------------
      IF (L_RP .or. L_RP2) THEN
       G0=G0_RP
      ELSE
       G0=10.0
      ENDIF

      DH=G0/EH                 ! Used in calc of stability function FH.
      DM=G0/EM                 ! Used in calc of stability function FM.


! Settings for LEM stability functions
      PRN = 0.7
      SUBA = 1.0/PRN
      SUBB = 40.0
      SUBC = 16.0
      SUBG = 1.2
      SUBH = 0.0
      RIC = 0.25
      RICINV = 1./RIC

!  Set MBL, "maximum number of boundary layer levels" for the purposes
!  of mixed layer height calculation.

      MBL = BL_LEVELS - 1

!-----------------------------------------------------------------------
!! 0. Initialise flag for having reached top of turbulently mixed layer
!!    and also the number of turbulently mixed layers and weight array
!-----------------------------------------------------------------------

      do j=1,rows
      DO I=1,row_length
        RIOUT(i,j,BL_LEVELS) = 0.0
        TOPBL(I,j)     = .FALSE.
        BL_weight(I,J) = 1.0
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
! Initialise 3D stability functions
!-----------------------------------------------------------------------
      If (L_subfilter_vert .or. L_subfilter_horiz) then
        Do k = 1, bl_levels
          Do j = 1, rows
            Do i = 1, row_length
              FM_3D(i,j,k) = 0.0
              FH_3D(i,j,k) = 0.0
            End Do
          End Do
        End Do
      End If
!-----------------------------------------------------------------------
!  1.1 Loop over levels calculating Richardson numbers.
!-----------------------------------------------------------------------
!
      DO K=2,BL_LEVELS
        KM1 = K-1
        do j=1,rows
        DO I=1,row_length

          RI(I,j,K) = DBDZ(I,j,K) / ( DVDZM(I,j,K)*DVDZM(I,j,K) )

! New boundary layer diagnostics
        IF (BL_diag%L_gradrich) BL_diag%gradrich(i,j,k)=ri(i,j,k)
        IF (BL_diag%L_dbdz)     BL_diag%dbdz(i,j,k)=dbdz(i,j,k)
        IF (BL_diag%L_dvdzm)    BL_diag%dvdzm(i,j,k)=dvdzm(i,j,k)

          RIOUT(i,j,k-1) = RI(i,j,k)  ! so RIOUT(K) is on theta-level K

!-----------------------------------------------------------------------
!! 1.2 If either a stable layer (Ri > 1) or the maximum boundary layer
!!     height has been reached, set boundary layer height (ZH_LOCAL) to
!!     the height of the lower boundary of the current layer
!-----------------------------------------------------------------------

          IF ( .NOT.TOPBL(I,j) .AND.                                    &
     &         (RI(I,j,K) >  1.0 .OR. K >  MBL) )  THEN
            TOPBL(I,j) = .TRUE.
            ZH_LOCAL(I,j) = Z_UV(I,j,K)
            NTML_LOCAL(I,j) = K-1
          ENDIF
        ENDDO  ! Loop over points
        ENDDO  ! Loop over points
      ENDDO  ! Loop over levels
!-----------------------------------------------------------------------
!! If NTML_LOCAL is greater than the top of the parcel ascent (NTPAR)
!! for a cumulus-capped layer, shear driven mixing is allowed to
!! dominate (if ISHEAR_BL=1 selected)
!-----------------------------------------------------------------------
      do j=1,rows
      DO I=1,row_length
        IF ( ISHEAR_BL  ==  1 .AND.                                     &
     &       NTML_LOCAL(I,j)  >   NTPAR(I,j) ) THEN
          CUMULUS(I,j) = .FALSE.
        ENDIF
      ENDDO
      ENDDO
!-----------------------------------------------------------------------
!! In CUMULUS layers the local scheme is capped at the LCL (given in
!! this case by NTML_NL).
!-----------------------------------------------------------------------
      do j=1,rows
      DO I=1,row_length
        IF ( CUMULUS(I,j) ) THEN
          NTML_LOCAL(I,j) = NTML_NL(I,j)
          ZH_LOCAL(I,j) = Z_UV(I,j,NTML_LOCAL(I,j)+1)

        ENDIF
      ENDDO
      ENDDO

!-----------------------------------------------------------------------
!! 2.  Richardson Number based local mixing scheme
!-----------------------------------------------------------------------
!! 2.0 Loop round "boundary" levels; calculate the stability-
!!     dependent turbulent mixing coefficients.
!-----------------------------------------------------------------------
! TKE budget based depth diagnosis
!
! Starting with the definition of the flux Richardson
! number, assuming similarity profiles for
! stress and buoyancy flux, and vertically integrating
! gives an expression for the stable boundary layer
! depth which is based just on surface fluxes and
! the wind speed change across the boundary layer.
!
! ----------------------------------------------------

      IF (SBL_OP  ==  Depth_based) THEN

! Index for assumed buoyancy profile
        m_buoy=1.0

! Index for assumed stress profile
        m_tau=1.0

! Effective bulk flux Richardson number
        rifb=0.3

        ind=m_buoy-m_tau+1.0

        DO J=1, ROWS
          DO I=1,ROW_LENGTH

! Set diff_min to a large initial value
            diff_min(i,j)=1000.0

! Surface Obukhov length
            MOsurf(I,j)= -V_S(I,j)*V_S(I,j)*V_S(I,j)                    &
     &                  /(VKMAN*FB_SURF(I,j))
          ENDDO
        ENDDO

        DO K=2,BL_LEVELS
          DO J=1, ROWS
            DO I=1,ROW_LENGTH

              ! The wind speed change from level k to the surface
              U1=sqrt(U_P(I,j,K)*U_P(I,j,K)+V_P(I,j,K)*V_P(I,j,K))

              ! h_est is the estimate of the stable boundary layer
              ! depth using the TKE based formula
              h_est=VKMAN*MOsurf(i,j)*ind*rifb*U1/V_S(i,j)

              ! Absolute difference between height and estimate
              diff=ABS(z_uv(i,j,k)-h_est)

              ! If h_est is closer than the previous closest value
              ! (diff_min) reset the h_tkeb to h_est

              IF (diff  <   diff_min(i,j)) THEN
                diff_min(i,j)=diff
                h_tkeb(i,j)=h_est
              ENDIF

            ENDDO
          ENDDO
        ENDDO

      ENDIF   ! SBL_OP = Depth_based

      DO K=2,BL_LEVELS

! ----------------------------------------------------------------
! Load up 2D array FUNC with selected stability function
!
!  SBL_OP                 Option
!
!  Long_tails             Long tails
!  Sharpest               SHARPEST function
!  Sharp_sea_long_land    SHARPEST over sea ; Long tails over land
!  Mes_tails              MESOSCALE model: Louis/SHARPEST blend
!  Louis_tails            Louis function
!  Depth_based            Boundary layer depth based formulation
!  Sharp_sea_mes_land     SHARP over sea; Mes over land
!  Sharp_sea_Louis_land   SHARP over sea; Louis over land
! ----------------------------------------------------------------

        SELECT CASE (SBL_OP)

          !--------------------------------------------
          ! LONG TAILS
          !--------------------------------------------
          CASE(Long_tails)

          DO J=1, ROWS
            DO I=1,ROW_LENGTH
               FUNC(I,J)=1.0 / ( 1.0 + G0 * RI(I,J,K) )
            ENDDO
          ENDDO

          !--------------------------------------------
          ! SHARP TAILS
          !--------------------------------------------
          CASE(Sharpest)

          DO J=1, ROWS
            DO I=1,ROW_LENGTH
              IF (RI(I,J,K)  <   1.0/G0) THEN
                FUNC(I,J) = 1.0 - 0.5 * G0 * RI(I,j,K)
              ELSE
                FUNC(I,J) = 1.0 / ( 2.0 * G0 * RI(I,j,K) )
              END IF
              FUNC(I,J)=FUNC(I,J)*FUNC(I,J)
            ENDDO
          ENDDO

          !--------------------------------------------
          ! SHARP over sea; long tails over land
          !--------------------------------------------
          CASE(Sharp_sea_long_land)

          DO J=1, ROWS
            DO I=1, ROW_LENGTH
              IF (FLANDG(i,j) < 0.5) THEN
                ! SHARPEST over sea
                IF (RI(I,J,K)  <   1.0/G0) THEN
                  FUNC(I,J) = 1.0 - 0.5 * G0 * RI(I,j,K)
                ELSE
                  FUNC(I,J) = 1.0 / ( 2.0 * G0 * RI(I,j,K) )
                END IF
                FUNC(I,J)=FUNC(I,J)*FUNC(I,J)
              ELSE
                ! Long tails over land
                FUNC(I,J)= 1.0 / ( 1.0 + G0 * RI(I,j,K) )
              ENDIF
            ENDDO
          ENDDO


          !--------------------------------------------
          ! MESOSCALE MODEL TAILS
          !--------------------------------------------
          CASE(Mes_tails)

          DO J=1, ROWS
            DO I=1,ROW_LENGTH
              ! Louis function
              FM = 1.0 / ( 1.0+5.0*RI(I,j,K) )
              FM_LOUIS = FM * FM
              ! code for SHARPEST
              IF (RI(I,j,K)  <   1.0/G0) THEN
                FM = 1.0 - 0.5 * G0 * RI(I,j,K)
              ELSE
                FM = 1.0 / ( 2.0 * G0 * RI(I,j,K) )
              END IF
              FM_SHARPEST = FM * FM
              ! Linear weighting function giving Louis
              ! at z=0, SHARPEST above Z_SCALE
              Z_SCALE = 200.0
              IF ( Z_TQ(I,j,K-1)  >=  Z_SCALE ) THEN
                FUNC(I,J) = FM_SHARPEST
              ELSE
                FUNC(I,J)= FM_LOUIS *( 1.0 - Z_TQ(I,j,K-1)/Z_SCALE )    &
     &                     + FM_SHARPEST * Z_TQ(I,j,K-1)/Z_SCALE
              ENDIF
            ENDDO
          ENDDO

          !--------------------------------------------
          ! LOUIS TAILS
          !--------------------------------------------
          CASE(Louis_tails)

          ! LOUIS FUNCTION
          DO J=1, ROWS
            DO I=1,ROW_LENGTH
            FUNC(I,J)=1.0 / ( 1.0+5.0*RI(I,j,K) )
            FUNC(I,J)=FUNC(I,J)*FUNC(I,J)
            ENDDO
          ENDDO

          !--------------------------------------------
          ! LONG TAILS FOR USE WITH DEPTH BASED SCHEME
          !--------------------------------------------
          CASE(Depth_based)
          ! LONG TAILS
          DO J=1, ROWS
            DO I=1,ROW_LENGTH
              FUNC(I,J)=1.0 / ( 1.0 + G0 * RI(I,J,K) )
            ENDDO
          ENDDO

          !--------------------------------------------
          ! SHARP TAILS OVER SEA; MES TAILS OVER LAND
          !--------------------------------------------
          CASE(Sharp_sea_mes_land)
          ! SHARP sea; MES land
          DO J=1, ROWS
            DO I=1,ROW_LENGTH
              ! Louis function
              FM = 1.0 / ( 1.0+5.0*RI(I,j,K) )
              FM_LOUIS = FM * FM
              ! code for SHARPEST
              IF (RI(I,j,K)  <   1.0/G0) THEN
                FM = 1.0 - 0.5 * G0 * RI(I,j,K)
              ELSE
                FM = 1.0 / ( 2.0 * G0 * RI(I,j,K) )
              END IF
              FM_SHARPEST = FM * FM
              ! Linear weighting function giving Louis at z=0,
              ! SHARPEST above Z_SCALE
              Z_SCALE = 200.0

              IF (FLANDG(i,j) < 0.5) THEN
                ! SHARP sea
                FUNC(i,j)=FM_SHARPEST
              ELSE
                ! MES land
                IF ( Z_TQ(I,j,K-1)  >=  Z_SCALE ) THEN
                  FUNC(I,J) = FM_SHARPEST
                ELSE
                  FUNC(I,J)= FM_LOUIS *( 1.0 - Z_TQ(I,j,K-1)/Z_SCALE )  &
     &                       + FM_SHARPEST * Z_TQ(I,j,K-1)/Z_SCALE
                ENDIF

              ENDIF  ! FLANDG(i,j) < 0.5


            ENDDO ! loop over row_length
          ENDDO ! loop over rows

          !--------------------------------------------
          ! SHARP TAILS OVER SEA; LOUIS TAILS OVER LAND
          !--------------------------------------------
          CASE(Sharp_sea_Louis_land)
          ! SHARP sea; Louis land
            DO J=1, ROWS
            DO I=1,ROW_LENGTH

              IF (FLANDG(i,j) < 0.5) THEN
                ! SHARP sea
                IF (RI(I,j,K)  <   1.0/G0) THEN
                  FM = 1.0 - 0.5 * G0 * RI(I,j,K)
                ELSE
                  FM = 1.0 / ( 2.0 * G0 * RI(I,j,K) )
                END IF
                FUNC(i,j)=FM * FM
              ELSE
                ! Louis land
                FM_LOUIS = 1.0 / ( 1.0+5.0*RI(I,j,K) )
                FUNC(I,J)= (1.0 - WeightLouisToLong) *                  &
     &                     FM_LOUIS * FM_LOUIS +                        &
     &                     WeightLouisToLong *                          &
     &                     1.0 / ( 1.0 + G0 * RI(I,j,K) )
              ENDIF  ! FLANDG(i,j) < 0.5

            ENDDO ! loop over row_length
            ENDDO ! loop over rows

        END SELECT ! SBL_OP

!------------------------------------------------------------------
! Additional code to allow the local Ri scheme to be used 
! sensibly in the free atmosphere, ie above the BL top,
! regardless of the tail option selected above.
!------------------------------------------------------------------

        IF (LOCAL_FA >= 1) THEN
!         !----------------------------
!         ! Calculate SHARPEST function
!         !----------------------------
          DO J=1, ROWS
          DO I=1,ROW_LENGTH
            IF (RI(I,J,K)  <   1.0/G0) THEN
              SHARP(I,J) = 1.0 - 0.5 * G0 * RI(I,j,K)
            ELSE
              SHARP(I,J) = 1.0 / ( 2.0 * G0 * RI(I,j,K) )
            END IF
            SHARP(I,J)=SHARP(I,J)*SHARP(I,J)
          ENDDO
          ENDDO

          do j=1,rows
          DO I=1,row_length
!           !---------------------------------------------------------
!           ! Optionally weight between standard tail (in FUNC) near 
!           ! the ground and SHARP tail at a height away from the 
!           ! surface
!           !---------------------------------------------------------
            IF (LOCAL_FA == 1) THEN
              Z_SCALE = 1000.0
!             !---------------------------------------------------------
!             ! Set to 1km to mimic old value of BL_LEVELS, 
!             ! gives BL_weight~0 by 2km, ~0.95 at 500m
!             !---------------------------------------------------------
              BL_weight(I,J) = 0.5*( 1.0 -                              &
     &                        TANH(3.*((Z_TQ(I,j,K-1)/Z_SCALE)-1.0) ) )
              FUNC(I,J) = FUNC(I,J) * BL_weight(I,J)                    &
     &             + SHARP(I,J)*( 1.0 - BL_weight(I,J) )
            ELSE IF (LOCAL_FA == 2 .OR. LOCAL_FA == 3) THEN
!             !-------------------------------------------------------
!             ! Weight dependent on subgrid orography
!             !-------------------------------------------------------
              FM = SHARP(I,J)
              IF ( SIGMA_H(I,J) .GT. 0.0 ) THEN
                ZPR = Z_TQ(I,j,K-1)/SIGMA_H(I,j)
                IF ( ZPR .LT. 2.0 ) THEN
                  BL_weight(I,J) = (1.0 - 0.5*ZPR) / (0.25 + 10.*ZPR)

                  FM = FUNC(I,J) * BL_weight(I,J)                       &
     &                 + SHARP(I,J)*( 1.0 - BL_weight(I,J) )
                END IF
              END IF
              FUNC(I,J) = FM
            END IF
          ENDDO
          ENDDO

        END IF

!---------------------------------------------------------------
! Set stable and neutral Prandtl number (=KM/KH) 
!---------------------------------------------------------------
        IF (PRANDTL == Constant_SBL) THEN

          PR_N = 1.0

          DO J=1,rows
          DO I=1,row_length
            PRANDTL_NUMBER(I,J) = 1.0
          ENDDO
          ENDDO

        ELSE IF (PRANDTL == LockMailhot2004) THEN

          PR_N=0.7

          DO J=1,rows
          DO I=1,row_length
            PRANDTL_NUMBER(I,J) = PR_N*( 1.0 + 2.0*RI(I,J,K) )
          ENDDO
          ENDDO

        END IF

        IF (SBL_OP  ==  Depth_based) PR_N=0.7

        KM1 = K-1
        do j=1,rows
        DO I=1,row_length
!-----------------------------------------------------------------------
!! 2.1 Calculate asymptotic mixing lengths LAMBDAM and LAMBDAH
!!     (may be equal or LambdaM=2*LambdaH (operational setting)).
!-----------------------------------------------------------------------
          
          IF (L_LAMBDAM2) THEN
            IF (L_RP .or. L_RP2) THEN
             LAMBDAM = MAX (LAMBDA_MIN , 2*par_mezcla*ZH_LOCAL(I,j))
             LAMBDAH = MAX (LAMBDA_MIN , par_mezcla*ZH_LOCAL(I,j))
            ELSE
             LAMBDAM = MAX (LAMBDA_MIN , 0.30*ZH_LOCAL(I,j))
             LAMBDAH = MAX (LAMBDA_MIN , 0.15*ZH_LOCAL(I,j))
            END IF
          ELSE
            IF (L_RP .or. L_RP2) THEN
             LAMBDAM = MAX ( LAMBDA_MIN , par_mezcla*ZH_LOCAL(I,j) )
            ELSE
             LAMBDAM = MAX ( LAMBDA_MIN , 0.15*ZH_LOCAL(I,j) )
            ENDIF
            LAMBDAH = LAMBDAM
          ENDIF
          
          LAMBDA_EFF = MAX (LAMBDAM, A_LAMBDA*H_BLEND(I,j) )
!-----------------------------------------------------------------------
!         ! Optionally reduce mixing length above local BL top
!-----------------------------------------------------------------------
          IF (K >= NTML_LOCAL(I,j)+2 .and. .NOT.L_FULL_LAMBDAS) THEN
            LAMBDAM = LAMBDA_MIN
            LAMBDAH = LAMBDA_MIN
            IF (Z_TQ(I,j,K-1) > A_LAMBDA*H_BLEND(I,j))                  &
     &                                        LAMBDA_EFF=LAMBDA_MIN
          END IF
          IF ( K >= NTML_LOCAL(I,j)+2 .and. L_FULL_LAMBDAS .and.        &
     &         LOCAL_FA == 1) THEN
!             ! Weight lambda to lambda_min with height
!             ! Assuming only LOCAL_FA = 1 will have L_FULL_LAMBDAS
!             ! If other LOCAL_FA options are coded here then 
!             ! changes must be included in section 5.3 of bdy_expl2

              LAMBDA_EFF = LAMBDA_EFF * BL_weight(I,j)                  &
     &                     + LAMBDA_MIN*( 1.0 - BL_weight(I,j) )
              LAMBDAH    = LAMBDAH * BL_weight(I,j)                     &
     &                     + LAMBDA_MIN*( 1.0 - BL_weight(I,j) )
          END IF
!-----------------------------------------------------------------------
!! 2.2 Calculate mixing lengths ELH, ELM coincident with RI(K) and so
!!     at Z_TQ(K-1)
!-----------------------------------------------------------------------

!  Incorporate log profile corrections to the vertical finite
!  differences into the definitions of ELM and ELH.
!  Note that ELH is calculated here for the unstable stability
!  functions - for RHOKH it will be calculated after interpolation
!  to uv-levels in BDYLYR.
!  To save computing logarithms for all K, the values of ELM and ELH
!  are unchanged for K > K_LOG_LAYR.

          IF (K  <=  K_LOG_LAYR) THEN
            VKZ   = VKMAN * ( Z_UV(I,j,K) - Z_UV(I,j,K-1) )
            F_LOG = LOG( ( Z_UV(I,j,K) + Z0M(I,j)   ) /                 &
     &                   ( Z_UV(I,j,K-1) + Z0M(I,j) ) )
            ELM = VKZ / ( F_LOG + VKZ/LAMBDA_EFF )
            ELH = VKZ / ( F_LOG + VKZ/LAMBDAH )
          ELSE
            VKZ = VKMAN * ( Z_TQ(I,j,K-1) + Z0M(I,j) )
            ELM = VKZ / (1.0 + VKZ/LAMBDA_EFF )
            ELH = VKZ / (1.0 + VKZ/LAMBDAH )
          ENDIF

!-----------------------------------------------------------------------
!! 2.4 Calculate (values of) stability functions FH, FM.
!-----------------------------------------------------------------------

          IF (RI(I,j,K)  >=  0.0) THEN
            RTMRI = 0.0

!           !---------------------------------------------------------
!           ! Set FM to the standard requested stability function and 
!           ! scale FH by the neutral Prandtl number, 
!           ! rather than keep FH the same and scale FM.
!           !---------------------------------------------------------
            FM = FUNC(I,J)
            FH = FM / PR_N

!           !-----------------------------------------------------------
!           ! Then, if requested, rescale FM by the full Prandtl number.
!           ! The reason for coding in this way is so that FM=1 at 
!           ! neutral and SHARP+Prandtl gives FH~1/Ri^2 and FM~1/Ri
!           !-----------------------------------------------------------
            FM = FH * PRANDTL_NUMBER(I,J)
!           !-----------------------------------------------------------
!           ! If convective cloud exists in layer K allow neutral mixing
!           ! of momentum between layers K-1 and K. This is to ensure
!           ! that a reasonable amount of momentum is mixed in the
!           ! presence of convection; it is not be required when
!           ! momentum transport is included in the convection scheme.
!           !-----------------------------------------------------------

            IF ( .NOT.L_MOM .AND. (CCA(I,j)  >   0.0) .AND.             &
     &           (K  >=  CCB(I,j)) .AND. (K  <   CCT(I,j)) )            &
     &         FM = 1.0
          ELSE
            RTMRI = (ELM/ELH) * SQRT ( -RI(I,j,K) )
            FM = 1.0 - ( G0*RI(I,j,K) / ( 1.0 + DM*RTMRI ) )
            FH = ( 1.0 - ( G0*RI(I,j,K) / ( 1.0 + DH*RTMRI ) ) ) / PR_N
          ENDIF

! If LEM_stability
          
          IF (SBL_OP == LEM_stability) THEN

            IF ((RI(I,j,K) >= 0.0) .AND. (RI(I,j,K)< Ric))  THEN
              RTMRI = 0.0
              RIFAC=(1.0-Ri(i,j,k)*RICINV)**4
              FM = RIFAC*(1.0-SUBH*Ri(i,j,k))
              FH = RIFAC*SUBA*(1.0-SUBG*Ri(i,j,k))
            ELSE IF (Ri(i,j,k) >= Ric) THEN
              FM=0.0
              FH=0.0
            ELSE
              FM = SQRT(1.0-SUBC*Ri(i,j,k))
              FH = SUBA*SQRT(1.0-SUBB*Ri(i,j,k)) 
            END IF           

          END IF ! LEM_stability 


!-----------------------------------------------------------------------
!! 2.5 Calculate exchange coefficients RHO*KM(K), RHO*KH(K)
!!     both on TH-level K-1 at this stage (RHOKH will be interpolated
!!     onto uv-levels and then be multiplied by ELH in BDYLYR)
!-----------------------------------------------------------------------

          If (L_subfilter_vert .or. L_subfilter_horiz) then
            FM_3D(I,j,K)=FM
            FH_3D(I,j,K)=FH
          End If
          RHOKM(I,j,K) = RHO(I,j,K-1) * ELM*ELM * DVDZM(I,j,K) * FM
          If (Lq_mix_bl) Then
!           ! Note "RHO" here is always wet density (RHO_TQ) so 
!           ! save multiplication of RHOKH to after interpolation
            RHOKH(I,j,K) =                  ELM * DVDZM(I,j,K) * FH
          Else
            RHOKH(I,j,K) = RHO(I,j,K-1)   * ELM * DVDZM(I,j,K) * FH
          End if
! -------------------------------------------
! Boundary layer depth based formulation
! -------------------------------------------

          IF (SBL_OP  ==  Depth_based .AND.                             &
     &        FB_SURF(i,j)  <=  0.0) THEN
            IF (z_tq(i,j,k-1) < h_tkeb(i,j)) THEN

              ! Formula for diffusion coefficient
              ! see Beare et al 2006, Boundary layer Met.

              KM = V_S(i,j) * VKMAN * z_tq(i,j,k-1) *                   &
     &                      ( (1.0-z_tq(i,j,k-1)/h_tkeb(i,j))**(1.5) )  &
     &                      /  (1.0 + 4.7*z_tq(i,j,k-1)/MOsurf(I,j))
              RHOKM(i,j,k)= RHO(I,j,K-1) * KM
              If (Lq_mix_bl) Then
!               ! Note "RHO" here is always wet density (RHO_TQ) so 
!               ! save multiplication of RHOKH to after interpolation
                RHOKH(i,j,k)= KM /(ELM*pr_n)
              Else
                RHOKH(I,j,K) = RHO(I,j,K-1)*KM /(ELM*pr_n)
              End if
            ELSE
              RHOKM(i,j,k)=0.0
              RHOKH(i,j,k)=0.0
            ENDIF

          ENDIF   !SBL_OP  ==  Depth_based


        ENDDO !I
        ENDDO !j
      ENDDO ! bl_levels

!-----------------------------------------------------------------------
!! 3.  Equilibrium Stable Boundary Layer (SBL) model.
!-----------------------------------------------------------------------
      IF (L_SBLeq) THEN

!-----------------------------------------------------------------------
!! 3.1 On first timestep only, calculate "G"-tables.
!!     Note that GcalcDone, HLtab, GMsav and GHsav are protected
!!     by a SAVE statement, so their values are retained
!!     in subsequent calls to this subprogram. Gx values
!!     are calculated, for each value of HonL, as contributions
!!     to integrals of z/h from 0.1 to 1.0 in NGz steps of size Gdz
!-----------------------------------------------------------------------
      IF (.NOT.GcalcDone) THEN
        Mlamb=1.5*Muw_SBL-Mwt_SBL
        Gdz=0.9/REAL(NGz)
        DO GK=1,GN
          HonL=HLtab(GK)
          Gz =0.1-0.5*Gdz
          DO K=1,NGz
            Gz  =Gz+Gdz
            Zhat=Gz*HonL/((1.0-Gz)**(Mlamb))
! DEPENDS ON: sblequil
            CALL SBLequil(Zhat,lhat,PKM,PKH,PKUU,PKTT,PKUT,PKTU,        &
     &                    PHIe,PHIw,RiSBL,CE,rpow,CB,CN,iERRSBL,LTIMER)
            tmp=PKH*lhat*VKMAN*Zhat
            GH=Gdz*((1.0-Gz)**(2.0*(Mwt_SBL-Muw_SBL)))/tmp
            tmp=PKM*lhat*VKMAN*Zhat
            GM=Gdz*((1.0-Gz)**(     Mwt_SBL-Muw_SBL ))/tmp
            GHsav(GK,K)=GH
            GMsav(GK,K)=GM
          ENDDO
        ENDDO
        GcalcDone=.TRUE.
      ENDIF

!-----------------------------------------------------------------------
!! 3.2 Diagnose depth of equilibrium SBL (ZH_ESBL) on UV-levels (Z_UV).
!!     Estimate for ZH_ESBL obtained by sub-grid interpolation.
!!     NTML_ESBL is then the first UV-level below ZH_ESBL
!-----------------------------------------------------------------------
      do j=1,rows
      DO I=1,row_length
        !Calculate theta_v and interpolate to UV grid
        DO K=1,BL_LEVELS
          THv_TQ(I,j,K)=(TL(I,j,K)+GRCP*Z_TQ(I,j,K))                    &
     &                    *(1.0+C_VIRTUAL*QW(I,j,K))
        ENDDO
        K=1
        T1=THv_TQ(I,j,K)
        T2=Thv_TQ(I,j,K+1)
        DZ=Z_TQ(I,j,K+1)-Z_TQ(I,j,K)
        slope=(T2-T1)/DZ
        Thv(I,j,K+1)=T1+slope*(Z_UV(I,j,K+1)-Z_TQ(I,j,K)) !Lin interp
        Thv(I,j,1)=T1-slope*(Z_TQ(I,j,1)-Z_UV(I,j,1))
        DO K=2,BL_LEVELS-1
          T1=THv_TQ(I,j,K)
          T2=Thv_TQ(I,j,K+1)
          DZ=Z_TQ(I,j,K+1)-Z_TQ(I,j,K)
          slope=(T2-T1)/DZ
          Thv(I,j,K+1)=T1+slope*(Z_UV(I,j,K+1)-Z_TQ(I,j,K)) !Lin interp
        ENDDO
        invMOsurf(I,j)=-VKMAN*FB_SURF(I,j)                              &
     &               /(V_S(I,j)*V_S(I,j)*V_S(I,j))
        invMOsfce=invMOsurf(I,j)

        TOPBL(I,j)    =.FALSE.
        ZH_ESBL(I,j)  =Z_UV(I,j,kMINh)+1.0
        NTML_ESBL(I,j)=kMINh

        IF (FB_SURF(I,j) <= 0.0) THEN !only for stable/neutral cases
          subgrid =.false.
          kZtop   =kMINh-1

          DO WHILE ((.NOT.TOPBL(I,j)).AND.(kZtop <= BL_LEVELS-1))
            Zprev=Z_UV(I,j,kZtop)
            kZtop=kZtop+1
            Ztop =Z_UV(I,j,kZtop)
            Zbot =0.1*Ztop !Height of surface layer top
            kZbot=0
            zz   =0.0
            DO WHILE (zz <  Zbot) !find kZbot (UV-level above Zbot)
              kZbot=kZbot+1
              zz=Z_UV(I,j,kZbot)
            ENDDO !find kZbot
            IF (kZbot >  1) THEN
              !Interpolation of U and T to Zbot
              K=kZbot-1
              U1=sqrt(U_P(I,j,K)*U_P(I,j,K)+V_P(I,j,K)*V_P(I,j,K))
              U2=sqrt(U_P(I,j,K+1)*U_P(I,j,K+1)                         &
     &               +V_P(I,j,K+1)*V_P(I,j,K+1))
              T1=THv(I,j,K)
              T2=Thv(I,j,K+1)
              DZ=Z_UV(I,j,K+1)-Z_UV(I,j,K)
              Ubot=U1+(U2-U1)*(Zbot-Z_UV(I,j,K))/DZ !Linear interp
              Tbot=T1+(T2-T1)*(Zbot-Z_UV(I,j,K))/DZ !Linear interp
              kG0=1
            ELSE !kZbot=1
              !Start integration at Z_UV(1) and truncate G-values
              Ubot=sqrt(U_P(I,j,1)*U_P(I,j,1)                           &
     &                 +V_P(I,j,1)*V_P(I,j,1))
              Tbot=THv(I,j,1)
              tmp=(real(NGz)/0.9)*((Z_UV(I,j,1)/Ztop)-0.1)
              kG0=max(int(tmp)+1,1)
              fG0=max(REAL(kG0)-tmp,0.0)
            ENDIF

            !Estimate GX values from lookup table
            HonL=Ztop*invMOsfce
            IF (HonL <  HLtab(1)) HonL=HLtab(1)
            GK=2
            DO WHILE ((HonL >  HLtab(GK)).AND.(GK <  GN))
              GK=GK+1
            ENDDO
            GMtab1=0.0
            GMtab2=0.0
            GHtab1=0.0
            GHtab2=0.0
            IF ((kG0 >  1).AND.(fG0 >  0.0)) THEN
              GMtab1=GMtab1+fG0*GMsav(GK-1,kG0-1)
              GMtab2=GMtab2+fG0*GMsav(GK,  kG0-1)
              GHtab1=GHtab1+fG0*GHsav(GK-1,kG0-1)
              GHtab2=GHtab2+fG0*GHsav(GK,  kG0-1)
            ENDIF
            DO K=kG0,NGz
              GMtab1=GMtab1+GMsav(GK-1,K)
              GMtab2=GMtab2+GMsav(GK,K)
              GHtab1=GHtab1+GHsav(GK-1,K)
              GHtab2=GHtab2+GHsav(GK,K)
            ENDDO
            H1=HLtab(GK-1)
            H2=HLtab(GK)
            G1=GMtab1
            G2=GMtab2
            Gs=(G2-G1)/(H2-H1)
            GM=Gs*(HonL-H2)+G2
            G1=GHtab1
            G2=GHtab2
            Gs=(G2-G1)/(H2-H1)
            GH=Gs*(HonL-H2)+G2
            Rilim=GH/(VKMAN*GM*GM)

            !Calculate Rtest
            U2=sqrt(U_P(I,j,kZtop)*U_P(I,j,kZtop)                       &
     &             +V_P(I,j,kZtop)*V_P(I,j,kZtop))
            T2=Thv(I,j,kZtop)
            Ujmp=U2-Ubot
            Tjmp=T2-Tbot
            Rtest=RtestMin-0.1e-10
            IF (abs(Ujmp) >= 0.1e-10) THEN
              Rib=Ztop*(G/THv(I,j,2))*Tjmp/(Ujmp*Ujmp)
              Rtest=Rilim-Rib
            ENDIF
            IF (Rtest <= RtestMin) THEN
              !Ztop is at or above the actual SBL top
              TOPBL(I,j)=.TRUE.
              ZH_ESBL(I,j)=Ztop+1.0
              NTML_ESBL(I,j)=kZtop
              !If we are above level kMINh,
              !estimate ZH_ESBL using subgrid scheme (below)
              IF (kZtop >  kMINh) subgrid=.true.
            ENDIF
          ENDDO

          IF (subgrid .AND. sg_enabled) THEN !perform subgrid diagnosis
            DZ100=(Ztop-Zprev)/100.0
            Rtest=RtestMin+1.0 !Ensures initial entry to WHILE loop
            Ztop=Zprev
            Zupper=ZH_ESBL(I,j)
            DO WHILE( (Rtest >  RtestMin).AND.                          &
     &               (Ztop <= Zupper-1.0-DZ100) )
              Zprev=Ztop
              Ztop=Ztop+DZ100
              Zbot=0.1*Ztop !Height of surface layer top
              kZbot=0
              zz=0.0
              DO WHILE (zz <  Zbot) !find kZbot (UV-level above Zbot)
                kZbot=kZbot+1
                zz=Z_UV(I,j,kZbot)
              ENDDO !find kZbot
              IF (kZbot >  1) THEN
                !Interpolation of U to Zbot
                K=kZbot-1
                U1=sqrt(U_P(I,j,K)*U_P(I,j,K)+V_P(I,j,K)*V_P(I,j,K))
                U2=sqrt(U_P(I,j,K+1)*U_P(I,j,K+1)                       &
     &                 +V_P(I,j,K+1)*V_P(I,j,K+1))
                T1=THv(I,j,K)
                T2=Thv(I,j,K+1)
                DZ=Z_UV(I,j,K+1)-Z_UV(I,j,K)
                Ubot=U1+(U2-U1)*(Zbot-Z_UV(I,j,K))/DZ !Linear interp
                Tbot=T1+(T2-T1)*(Zbot-Z_UV(I,j,K))/DZ !Linear interp
                kG0=1
              ELSE !kZbot=1
                !Start integration at Z_UV(1) and truncate G-values
                Ubot=sqrt(U_P(I,j,1)*U_P(I,j,1)                         &
     &                   +V_P(I,j,1)*V_P(I,j,1))
                Tbot=THv(I,j,1)
                tmp=(real(NGz)/0.9)*((Z_UV(I,j,1)/Ztop)-0.1)
                kG0=max(int(tmp)+1,1)
                fG0=max(REAL(kG0)-tmp,0.0)
              ENDIF

              !Estimate GX values from lookup table
              HonL=Ztop*invMOsfce
              IF (HonL <  HLtab(1)) HonL=HLtab(1)
              GK=2
              DO WHILE ((HonL >  HLtab(GK)).AND.(GK <  GN))
                GK=GK+1
              ENDDO
              GMtab1=0.0
              GMtab2=0.0
              GHtab1=0.0
              GHtab2=0.0
              IF ((kG0 >  1).AND.(fG0 >  0.0)) THEN
                GMtab1=GMtab1+fG0*GMsav(GK-1,kG0-1)
                GMtab2=GMtab2+fG0*GMsav(GK,  kG0-1)
                GHtab1=GHtab1+fG0*GHsav(GK-1,kG0-1)
                GHtab2=GHtab2+fG0*GHsav(GK,  kG0-1)
              ENDIF
              DO K=kG0,NGz
                GMtab1=GMtab1+GMsav(GK-1,K)
                GMtab2=GMtab2+GMsav(GK,K)
                GHtab1=GHtab1+GHsav(GK-1,K)
                GHtab2=GHtab2+GHsav(GK,K)
              ENDDO
              H1=HLtab(GK-1)
              H2=HLtab(GK)
              G1=GMtab1
              G2=GMtab2
              Gs=(G2-G1)/(H2-H1)
              GM=Gs*(HonL-H2)+G2
              G1=GHtab1
              G2=GHtab2
              Gs=(G2-G1)/(H2-H1)
              GH=Gs*(HonL-H2)+G2
              Rilim=GH/(VKMAN*GM*GM)

              !Calculate Rtest
              K=kZtop-1
              U1=sqrt(U_P(I,j,K)*U_P(I,j,K)                             &
     &               +V_P(I,j,K)*V_P(I,j,K))
              U2=sqrt(U_P(I,j,K+1)*U_P(I,j,K+1)                         &
     &               +V_P(I,j,K+1)*V_P(I,j,K+1))
              T1=Thv(I,j,K)
              T2=Thv(I,j,K+1)
              DZ=Z_UV(I,j,K+1)-Z_UV(I,j,K)
              slope=(T2-T1)/DZ
              T2=T1+slope*(Ztop-Z_UV(I,j,K))
              slope=(U2-U1)/DZ
              U2=U1+slope*(Ztop-Z_UV(I,j,K))
              Ujmp=U2-Ubot
              Tjmp=T2-Tbot
              Rtest=RtestMin-0.1e-10
              IF (abs(Ujmp) >= 0.1e-10) THEN
                Rib=Ztop*(G/THv(I,j,2))*Tjmp/(Ujmp*Ujmp)
                Rtest=Rilim-Rib
              ENDIF
              IF (Rtest <= RtestMin) ZH_ESBL(I,j)=Ztop
            ENDDO !DO WHILE
            IF (Rtest <= RtestMin) NTML_ESBL(I,j)=kZtop-1
          ENDIF !end subgrid ZH diagnosis

        ENDIF !only for stable/neutral cases
      ENDDO !I
      ENDDO !j

!-----------------------------------------------------------------------
!! 3.3 Prescribe local equilibrium fluxes and M-O length;
!!     calculate Zhat and invoke equilibrium model for SBL;
!!     convert output to diffusivities*density: RHO*Kx(K).
!-----------------------------------------------------------------------
      do j=1,rows
      DO I=1,row_length
        IF (FB_SURF(I,j) <= 0.0) THEN !only for stable/neutral cases

          DO K=2,BL_LEVELS !Set Kx back to zero
            RHOKM(I,j,K) = 0.0
            RHOKH(I,j,K) = 0.0
          ENDDO
          UWsfce=V_S(I,j)*V_S(I,j)
          WTsfce=FB_SURF(I,j) !Actually <w*THv>*(g/THv) at sfce
          invMOsfce=invMOsurf(I,j)
          HH=ZH_ESBL(I,j)
          DO K=2,NTML_ESBL(I,j)+1
            zz=Z_TQ(I,j,K-1)
            IF (zz <  HH) THEN !below ZH

              tmp=MAX(1.0-zz/HH,1.0E-3)
              USequil=sqrt(UWsfce*(tmp**Muw_SBL))
              WTequil=WTsfce*(tmp**Mwt_SBL) !Actually <w*THv>*(g/THv)
              invMOequil=invMOsfce/(tmp**(1.5*Muw_SBL-Mwt_SBL))
              Zhat   =zz*invMOequil
! DEPENDS ON: sblequil
              CALL SBLequil(Zhat,lhat,PKM,PKH,PKUU,PKTT,PKUT,PKTU,      &
     &                PHIe,PHIw,RiSBL,CE,rpow,CB,CN,iERRSBL,LTIMER)
              ll =lhat*VKMAN*zz
                KM =PKM *ll*USequil
                KH =PKH *ll*USequil
!               ! Note "RHO" here is always wet density (RHO_TQ)
                RHOKM(I,j,K) = RHO(I,j,K-1) * KM
                If (Lq_mix_bl) Then
!               ! Note "RHO" here is always wet density (RHO_TQ) so 
!               ! save multiplication of RHOKH to after interpolation
                  RHOKH(I,j,K) = KH
                Else
                  RHOKH(I,j,K) = RHO(I,j,K-1) * KH
                End If
              DBDU=-1.0
              IF ((DVDZM(I,j,K) >  1.e-10).AND.(DBDZ(I,j,K) >  1.e-10)) &
     &           DBDU=DBDZ(I,j,K)/DVDZM(I,j,K)
              IF (L_SBLco.AND.(DBDU >  0.0)) THEN
                KUU=PKUU*ll*USequil
                KTT=PKTT*ll*USequil
                KUT=PKUT*ll*ll !Actually KUT/(g/THv)
                tmp=WTequil/(USequil*USequil)
                GonTv=G/(  ( TL(I,j,K-1)+ GRCP * zz      )              &
     &                    *( 1.0 + C_VIRTUAL*QW(I,j,K-1) ) )
                KTU=-PKTU*ll*ll*tmp*tmp*GonTv*GonTv!Actually KTU*(g/THv)
                !Factor of g/THv in DBDU cancels in equations below,
                !due to factors of 1/(g/THv) in KUT, and g/THv in KTU
!               ! Note "RHO" here is always wet density (RHO_TQ)
                RHOKM(I,j,K) = RHO(I,j,K-1)*(KUU+KUT*DBDU)
                If (Lq_mix_bl) Then
!               ! Note "RHO" here is always wet density (RHO_TQ) so 
!               ! save multiplication of RHOKH to after interpolation
                  RHOKH(I,j,K) = (KTU/DBDU+KTT)
                Else
                  RHOKH(I,j,K) = RHO(I,j,K-1)*(KTU/DBDU+KTT)
                End If
              ENDIF

            ENDIF !below ZH
          ENDDO !K

        ENDIF !only for stable/neutral cases
      ENDDO !I
      ENDDO !j

!-----------------------------------------------------------------------
!! 3.4 Set values of ZH_LOCAL and NTML_LOCAL.
!-----------------------------------------------------------------------
      do j=1,rows
      DO I=1,row_length
        IF (FB_SURF(I,j) <= 0.0) THEN !only for stable/neutral cases
          ZH_LOCAL(I,j)=ZH_ESBL(I,j)
          K=NTML_ESBL(I,j)
          NTML_LOCAL(I,j)=K
          IF (Z_TQ(I,j,K) >  ZH_ESBL(I,j)) NTML_LOCAL(I,j)=K-1
        ENDIF !only for stable/neutral cases
      ENDDO !I
      ENDDO !j

      ENDIF !Equilibrium SBL model

!-----------------------------------------------------------------------
!! Finish up
!-----------------------------------------------------------------------

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('EX_COEF ',104)
      ENDIF

      RETURN
      END SUBROUTINE EX_COEF
