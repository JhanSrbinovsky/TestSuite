
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!!!  SUBROUTINE SF_OROG------------------------------------------------
!!!
!!!  Purpose: Calculate roughness lengths, blending height and wind
!!!           profile factor
!!!
! Modification History:
! Version Date     Change
!  5.2   15/11/00   New Deck         M. Best
!  6.1  08/12/03  Add !CDIR NODEP to force vectorisation. R Barnes
!  6.2  18/11/05  Replace l_z0_orog with FORMDRAG switch   S.B. Vosper
!!!
!!!--------------------------------------------------------------------
      SUBROUTINE SF_OROG(                                               &
     & ROW_LENGTH,ROWS,LAND_PTS,TILE_PTS,LAND_INDEX,TILE_INDEX,         &
     & FORMDRAG,FD_stab_dep,OROG_DRAG_PARAM,LTIMER,                     &
     & HO2R2_OROG,RIB,SIL_OROG,Z0M,Z1,                                  &
     & WIND_PROFILE_FACTOR,Z0M_EFF                                      &
     & )

      IMPLICIT NONE

      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                            ! IN Number of X points?
     &,ROWS                                                             &
                            ! IN Number of Y points?
     &,LAND_PTS                                                         &
                            ! IN Number of land points.
     &,TILE_PTS                                                         &
                            ! IN Number of tile points.
     &,LAND_INDEX(LAND_PTS)                                             &
                            ! IN Index of land points.
     &,TILE_INDEX(LAND_PTS)                                             &
                            ! IN Index of tile points
     &,FORMDRAG                                                         &
                            ! IN Switch for orographic drag
     &,FD_stab_dep          ! IN Switch to implement stability
!                           !    dependence of orographic form drag

      LOGICAL                                                           &
     & LTIMER               ! IN .TRUE. for timer diagnostics

      REAL                                                              &
     & HO2R2_OROG(LAND_PTS)                                             &
                            !IN Peak to trough height of unresolved
!                           !   orography divided by 2SQRT(2) (m).
     &,RIB(LAND_PTS)                                                    &
                            ! IN Bulk Richardson number for lowest layer
     &,SIL_OROG(LAND_PTS)                                               &
                            ! IN Silhouette area of unresolved orography
!                           !    per unit horizontal area
     &,OROG_DRAG_PARAM                                                  &
!                           ! IN Drag coefficient for orographic
!                           !    form drag
     &,Z0M(LAND_PTS)                                                    &
                            ! IN Roughness length for momentum (m).
     &,Z1(ROW_LENGTH,ROWS)  ! IN Height of lowest atmospheric level (m).

!  Output variables.

      REAL                                                              &
     & WIND_PROFILE_FACTOR(LAND_PTS)                                    &
!                           ! OUT For transforming effective surface
!                           !     transfer coefficients to those
!                           !     excluding form drag.
     &,Z0M_EFF(LAND_PTS)    ! OUT Effective roughness length for
!                                 momentum (m)

!  Work Variables

      INTEGER                                                           &
     & I,J                                                              &
                    ! Horizontal field index
     &,K                                                                &
                    ! Tile field index
     &,L            ! Land field index

      REAL                                                              &
     & H_BLEND_OROG                                                     &
                    ! Blending height
     &,RIB_FN                                                           &
                    ! Interpolation coefficient for 0 < RIB < RI_CRIT
     &,ZETA1                                                            &
                    ! Work space
     &,ZETA2                                                            &
                    ! More work space
     &,ZETA3        ! Even more work space

!   Common parameters

! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------
! C_VKMAN start
      REAL,PARAMETER:: VKMAN=0.4 ! Von Karman's constant
! C_VKMAN end
      ! Critical Richardson number, where Z0M_EFF=Z0M.
      ! Linear interpolation between RIB=0 and RI_CRIT
      REAL,PARAMETER:: RI_CRIT=0.5

      ! Tunable parameters in calculation of explicit orographic stresss
      REAL,PARAMETER:: ALPHA    = 12.0                                  &
                                             ! Tunable parameter for form
     &,                BETA     = 1.0                                   &
                                             ! drag calculation.
     &,                FD_DECAY = 0.3333                                &
                                             ! Height scale factors for
     &,                MAX_HT_SCALE  = 300.0                            &
                                             ! stress profiles
     &,                MIN_HT_SCALE  = 100.0
      LOGICAL,PARAMETER:: L_LOWHILL=.FALSE.  ! Set to .TRUE. for Wood &
!                                            ! Mason (1993) drag formula
!                                            ! Set to .FALSE. for steep
!                                            ! hill expression

!*----------------------------------------------------------------------
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
!   Local parameters

      REAL H_BLEND_MIN,H_BLEND_MAX
      PARAMETER (                                                       &
     & H_BLEND_MIN=0.0                                                  &
                             ! Minimun value of blending height
     &,H_BLEND_MAX=1000.0                                               &
                             ! Maximum value of blending height
     & )

!  External routines called :-
      EXTERNAL TIMER

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SF_OROG ',3)
      ENDIF

! Set blending height, effective roughness length and
! wind profile factor at land points.
!CDIR NODEP
      DO K=1,TILE_PTS
        L = TILE_INDEX(K)
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH

        WIND_PROFILE_FACTOR(L) = 1.0
        Z0M_EFF(L) = Z0M(L)

        IF (FORMDRAG ==  Effective_z0) THEN

          ZETA1 = 0.5 * OROG_DRAG_PARAM * SIL_OROG(L)
          ZETA2 = LOG ( 1.0 + Z1(I,J)/max(Z0M(L),1.e-6) )
          ZETA3 = 1.0 / SQRT ( ZETA1/(VKMAN*VKMAN) + 1.0/(ZETA2*ZETA2) )
          ZETA2 = 1.0 / EXP(ZETA3)
          H_BLEND_OROG = MAX ( Z1(I,J) / (1.0 - ZETA2) ,                &
     &                       HO2R2_OROG(L) * SQRT(2.0) )
          H_BLEND_OROG = MIN ( H_BLEND_MAX, H_BLEND_OROG )

! Apply simple stability correction to form drag if RIB is stable

          IF (SIL_OROG(L)  ==  RMDI) THEN
            ZETA1 = 0.0
          ELSE
            IF (FD_stab_dep == ON) THEN
              RIB_FN =  ( 1.0 - RIB(L) / RI_CRIT )
              IF (RIB_FN >  1.0) RIB_FN = 1.0
              IF (RIB_FN <  0.0) RIB_FN = 0.0
              ZETA1 = 0.5 * OROG_DRAG_PARAM * SIL_OROG(L) * RIB_FN
            ELSE 
              ZETA1 = 0.5 * OROG_DRAG_PARAM * SIL_OROG(L)
            END IF
          END IF

          Z0M_EFF(L) = H_BLEND_OROG / EXP ( VKMAN / SQRT ( ZETA1 +      &
     &                 (VKMAN / LOG ( H_BLEND_OROG / max(Z0M(L),1.e-6)) ) *       &
     &                 (VKMAN / LOG ( H_BLEND_OROG / max(Z0M(L),1.e-6)) ) ) )
          IF ( RIB(L) >  RI_CRIT .AND.                                  &
     &         FD_stab_dep == ON )  Z0M_EFF(L)=Z0M(L)

          WIND_PROFILE_FACTOR(L) = LOG( H_BLEND_OROG / Z0M_EFF(L) ) /   &
     &                             LOG( H_BLEND_OROG / max(Z0M(L),1.e-6) )

        END IF

      END DO

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SF_OROG ',4)
      END IF

      RETURN
      END SUBROUTINE SF_OROG