
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  Purpose: Calculate explicit flux of momentum in u or v direction
!!!
!!!
!!!  Model           Modification history
!!! version  Date
!!!  5.2   15/11/00   New Deck         M. Best
!!!  5.5   11/02/03   Add code for non-gradient stress parametrization.
!!!                               A.P.Lock
!!!  6.2   18/11/05   Add code for explicit orographic stress profile.
!!!                               S.B. Vosper
!!!  6.2   23/01/06   Improvements to Single Column Model Diagnostics
!!!                   System                          A. Kerr-Munslow
!!!
!!!  Programming standard: UMDP 3
!!!
!!!  Documentation: UM Documentation Paper No 24.
!!!
!!!---------------------------------------------------------------------

! SUBROUTINE EX_FLUX_UV

      SUBROUTINE EX_FLUX_UV (                                           &
     &  ROW_LENGTH                                                      &
     &, rows,off_x,off_y                                                &
     &, BL_LEVELS                                                       &
     &, U_V                                                             &
     &, ZH                                                              &
     &, RDZ_U_V                                                         &
     &, RHOKM_U_V                                                       &
     &, NG_STRESS, F_NGSTRESS_UV                                        &
     &, TAU_X_Y                                                         &
     &, TAU_XY_FD_UV, FORMDRAG                                          &
     &, TAU_GRAD,TAU_NON_GRAD                                           &
     &, LTIMER                                                          &
     &  )

      IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.

      LOGICAL                                                           &
     &  LTIMER                 ! IN Flag for TIMER diagnostics

      INTEGER                                                           &
     &  rows,off_x,off_y                                                &
     &, ROW_LENGTH                                                      &
                               ! IN No. of points in latitude row.
     &, BL_LEVELS                                                       &
                               ! IN No. of atmospheric levels for
!                                  which boundary layer fluxes are
!                                  calculated.
     &, FORMDRAG                                                        &
                               ! IN Switch for orographic form drag
     &, NG_STRESS
                               ! IN Switch for non-gradient stress


      REAL                                                              &
     &  RDZ_U_V (row_length,rows, 2:BL_LEVELS)                          &
!                                IN Reciprocal of the vertical
!                                   distance from level K-1 to
!                                   level K. (K > 1) on wind levels
     &, RHOKM_U_V (row_length,rows, BL_LEVELS)                          &
!                                IN Exchange coefficients for
!                                   momentum, on UV-grid with
!                                   first and last rows ignored.
!                                   for K>=2 (from KMKH).
     &, F_NGSTRESS_UV(row_length,rows,2:BL_LEVELS)                      &
!                                IN dimensionless function for
!                              !    non-gradient wind stress,
!                              !    either U or V depending on call
     &, U_V (1-off_x:row_length+off_x,1-off_y:rows+off_y, BL_LEVELS)    &
!                                IN Westerly_Southerly component of
!                                   wind.
     &, TAU_XY_FD_UV(row_length,rows, BL_LEVELS)                        &
!                              ! IN X/Y-component of form-drag stress
!                              !    at a UV point
     &, ZH(row_length,rows)    ! IN non-local BL depth


! ARGUMENTS WITH INTENT OUT. IE: INPUT VARIABLES CHANGED ON OUTPUT.

      REAL                                                              &
     &  TAU_X_Y (row_length,rows, BL_LEVELS)                            &
!                                OUT explicit x_y-component of
!                                    turbulent stress at levels
!                                    k-1/2; eg. TAUX(,1) is surface
!                                    stress. UV-grid, 1st and last rows
!                                    set to "missing data". (N/sq m)
     &, TAU_GRAD(row_length,rows,BL_LEVELS)                             &
!                                OUT k*du/dz grad stress (kg/m/s2)
     &, TAU_NON_GRAD(row_length,rows,BL_LEVELS)
!                                OUT Non-grad stress (kg/m/s2)

! LOCAL VARIABLES.
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

      REAL MAX_STRESS_GRAD
      PARAMETER ( MAX_STRESS_GRAD = 0.05 )
                               ! Maximum implied stress gradient
                               ! across the boundary layer used in the 
                               ! scaling of the non-gradient stresses
                               ! (m/s2)

      INTEGER                                                           &
     &  I                                                               &
     &, J                                                               &
     &, K                                                               &
     &, ERROR

      REAL                                                              &
     &  TAU_SURF(row_length,rows)                                       &
                               ! Explicit surface stress
     &, SIGN_TAU                                                        &
                               ! Sign of surface stress 
     &, BL_stress_grad
                               ! Stress gradient across BL

! External routines
      EXTERNAL TIMER


      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('EXFLUXUV',3)
      ENDIF

!-----------------------------------------------------------------------
!!  0.  Check that the scalars input to define the grid are consistent.
!-----------------------------------------------------------------------

      ERROR = 0

!-----------------------------------------------------------------------
!!  1.  Calculate "explicit" surface fluxes of momentum
!-----------------------------------------------------------------------

!     ! Set diagnostics to zero in level 1
      K=1
      DO J= 1, rows
      DO I= 1, row_length
        TAU_GRAD(I,j,k) = 0.0
        TAU_NON_GRAD(I,j,k) = 0.0
        TAU_SURF(I,J) = TAU_X_Y(I,j,1)
      END DO
      END DO

      IF (NG_STRESS == BrownGrant97_limited) THEN
!       ! Limit the explicitly calculated surface stress used to scale 
!       ! the non-gradient stress parametrization such that the implied 
!       ! stress gradient across the BL is less than MAX_STRESS_GRAD.
!       ! This has been found by experimentation to be sufficient to 
!       ! stop large T increments being generated in the dynamics 
!       ! solver, via large increments to W
        DO J= 1, rows
        DO I= 1, row_length
          SIGN_TAU = SIGN(1.0, TAU_X_Y(I,j,1) )
          BL_stress_grad = ABS( TAU_X_Y(I,j,1) )/ZH(I,j)
          BL_stress_grad = MIN( MAX_STRESS_GRAD, BL_stress_grad )
          TAU_SURF(I,J) = SIGN_TAU * ZH(I,j) * BL_stress_grad
        END DO
        END DO
      END IF

      DO K = 2,BL_LEVELS
        do j=1,rows
        do i=1,row_length

          TAU_GRAD(i,j,k) = RHOKM_U_V(i,j,k) *                          &
     &                    ( U_V(I,j,K) - U_V(I,j,K-1) ) *RDZ_U_V(I,j,K)
          TAU_NON_GRAD(i,j,k) = F_NGSTRESS_UV(i,j,k) * TAU_SURF(I,J)
          TAU_X_Y(i,j,k) = TAU_GRAD(i,j,k) + TAU_NON_GRAD(i,j,k)

! Add explicit orographic stress, noting that the surface stress
! is to be added later
           IF (FORMDRAG  ==  Explicit_stress) THEN
             TAU_X_Y(I,J,K) = TAU_X_Y(I,J,K) + TAU_XY_FD_UV(I,J,K)
           ENDIF

        END DO
        END DO
      END DO

    6 CONTINUE ! exit error point

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('EXFLUXUV',4)
      ENDIF
      RETURN

      END SUBROUTINE EX_FLUX_UV
