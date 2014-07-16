
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  Purpose: Calculate form drag profiles for distibuted
!!!           drag parametrization
!!!
!!!           Based on work by Wood, Brown and Hewer (2001),
!!!           Quart. J. Roy. Met. Soc., 127, 759--777.
!!!
!!!  Model              Modification history
!!! version  Date
!!!  6.2     23/11/05   New subroutine written by S.B. Vosper
!!!
!!!  Programming standard: UMDP 3
!!!
!!!---------------------------------------------------------------------

      SUBROUTINE FM_DRAG (                                              &
     &  ROW_LENGTH, ROWS, LAND_PTS, LAND_INDEX, FD_stab_dep             &
     &, BL_LEVELS                                                       &
     &, U_P, V_P, RHO_TQ, Z_UV, Z_TQ, Z0M                               &
     &, TAU_FD_X, TAU_FD_Y, ZH, RIB, SIL_OROG_LAND, OROG_DRAG_PARAM     &
     &, LTIMER, BL_diag )

      Use bl_diags_mod, Only :                                          &
          strnewbldiag

      IMPLICIT NONE
!
! Description:
! Turbulent form drag due to sub-grid orography
!
! Method:
! Calculates the form drag and stress profiles due to sub-grid scale
! orography. The stress is later added as an additional explicit
! stress to the boundary-layer equations in EXFXUV8A. The orographic
! stress is added to the surface stress in BDYLYR8A.
!
! Current Code Owner: Simon Vosper

! History:
! Version   Date     Comment
! -------   ----     -------
! 1.7       23/11/05 Original code, based on
!                    development work by Tom Allen.  Simon Vosper
!
! Code Description:
! Language: FORTRAN 77 + common extensions.
! This code is written to UMDP3 v6 programming standards.

      LOGICAL, Intent(IN):: LTIMER           ! Flag for TIMER diags

      INTEGER, Intent(IN):: FD_stab_dep      ! Switch to implement 
                                             ! stability dependence

      INTEGER, Intent(IN):: ROW_LENGTH, ROWS                            &
                                             ! Variables defining grid
     &,                     LAND_PTS                                    &
                                             ! Number of land points
     &,                     BL_LEVELS        ! No. of levels for which
!                                            ! boundary layer fluxes
!                                            ! are calculated

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
      Type (Strnewbldiag) :: BL_diag

      INTEGER, Intent(IN):: LAND_INDEX(land_pts) ! Index for compressed
!                                                ! land point array; ith
!                                                ! element holds
!                                                ! position in the FULL
!                                                ! field of the ith land
!                                                ! point to be processed

      REAL, Intent(IN)::                                                &
     & U_P(row_length,rows,BL_LEVELS)                                   &
                                      ! Wind component in x direction
!                                     ! horizontally interpolated to
!                                     ! P-grid (m/s)
     &,V_P(row_length,rows,BL_LEVELS)                                   &
                                      ! Wind component in y direction
!                                     ! horizontally interpolated to
!                                     ! P-grid (m/s)
     &,RHO_TQ(row_length,rows,BL_LEVELS)                                &
!                                     ! For a vertically staggered grid
!                                     ! with a u,v-level first above the
!                                     ! surface, RHO_TQ(*,K) is the
!                                     ! density of the k-th T,q-level
!                                     ! above the surface;
!                                     ! for an unstaggered grid the
!                                     ! densities at the layer interface
!                                     ! (half-levels) 1.5 to BL_LEVELS+0
!                                     ! should be input to elements 1 to
!                                     ! BL_LEVELS.
!                                     ! (Value for BL_LEVELS not used
!                                     ! in either case.)
     &,Z_UV(row_length,rows,BL_LEVELS)                                  &
!                                     ! For a vertically staggered grid
!                                     ! with a u,v-level first above the
!                                     ! surface, Z_UV(*,K) is the height
!                                     ! of the k-th u,v-level(half level
!                                     ! k-1/2) above the surface;
!                                     ! for an unstaggered grid the
!                                     ! heights of the half-levels
!                                     ! 0.5 to BL_LEVELS-0.5 should be
!                                     ! input to elements 1 to BL_LEVELS
!                                     ! (1st value not used in either
!                                     !  case.)
     &,Z_TQ(row_length,rows,BL_LEVELS)                                  &
!                                     ! For a vertically staggered grid
!                                     ! with a u,v-level first above the
!                                     ! surface, Z_TQ(*,K) is the height
!                                     ! of the k-th T,q-level (full
!                                     ! level k) above the surface;
!                                     ! for an unstaggered grid the
!                                     ! heights of the half levels
!                                     ! 1.5 to BL_LEVELS+0.5 should be
!                                     ! input to elements 1 to BL_LEVELS
!                                     ! (Value for BL_LEVELS not used
!                                     ! in either case.)
     &,Z0M(row_length,rows)                                             &
                                      ! Roughness length for momentum (m
     &,ZH(row_length,rows)                                              &
                                      ! Boundary layer height (actually
!                                     ! ZH_PREV)
     &,RIB(row_length,rows)                                             &
                                      ! Bulk Richardson number for
!                                     ! lowest layer
     &,SIL_OROG_LAND(land_pts)                                          &
!                                     ! Silhouette area of unresolved
!                                     ! orography per unit hoz. area
     &,OROG_DRAG_PARAM                ! Drag coefficient for form drag


      REAL, Intent(OUT)::                                               &
     & TAU_FD_X(row_length,rows,BL_LEVELS)                              &
                                           ! X-comp of orographic stress
     &,TAU_FD_Y(row_length,rows,BL_LEVELS) ! Y-comp of orographic stress
!                                          ! (N/m^2)

!-----------------------------------------------------------------------
!     External references:

      EXTERNAL TIMER

!-----------------------------------------------------------------------

!     Local and other symbolic constants:

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

!     Define local storage

!     Local arrays:

      REAL::                                                            &
     &  H_M(land_pts)                                                   &
                                     ! height at which form drag is
!                                    ! calculated (m)
     &, U_HM(land_pts)                                                  &
                                     ! X-component of wind at height h_m
     &, V_HM(land_pts)                                                  &
                                     ! Y-component of wind at height h_m
     &, FP_X(land_pts)                                                  &
                                     ! X-component of pressure force
     &, FP_Y(land_pts)               ! Y-component of pressure force

!     Local scalars:

      REAL::                                                            &
     &  ZETA                                                            &
                      ! Log(h_m/z0m)
     &, HEIGHT_FAC                                                      &
                      ! Height dependency factor for form drag
     &, WTA, WTB                                                        &
                      ! weights for interpolation between levels
     &, TAUSX,TAUSY                                                     &
                      ! Surface stress
     &, RIB_FN        ! Richardson number function for stability
!                     ! correction to drag

      INTEGER::                                                         &
     &  I                                                               &
                      ! Loop counter (horizontal field index)
     &, J                                                               &
                      ! Loop counter (offset within I loop)
     &, K                                                               &
                      ! Loop counter (vertical level index)
     &, L             ! Loop counter (horizontal land field index)

      IF (LTIMER) THEN
! DEPENDS ON: timer
       CALL TIMER('FM_DRAG ',103)
      ENDIF

!----------------------------------------------------------------
! 1. Calculate the height scale h_m and interpolate the wind and
!    density to this height.
!----------------------------------------------------------------

      DO l = 1, land_pts
        j = (land_index(l)-1)/row_length + 1
        i = land_index(l) - (j-1)*row_length

        h_m(l)=min(max_ht_scale,fd_decay*zh(i,j))
        h_m(l)=max(min_ht_scale,h_m(l))
      END DO

!     Interpolate to get U and V at z=h_m

      DO k = 2, BL_LEVELS
       DO l = 1, land_pts

        j = (land_index(l)-1)/row_length + 1
        i = land_index(l) - (j-1)*row_length

        IF(h_m(l) <= z_uv(i,j,k).and.h_m(l) >= z_uv(i,j,k-1)) THEN

          wta = ( h_m(l) - z_uv(i,j,k-1) )                              &
     &             /( z_uv(i,j,k) - z_uv(i,j,k-1) )
          wtb = ( z_uv(i,j,k) - h_m(l) )                                &
     &             /( z_uv(i,j,k) - z_uv(i,j,k-1) )
          u_hm(l)   = wta*u_p(i,j,k) + wtb*u_p(i,j,k-1)
          v_hm(l)   = wta*v_p(i,j,k) + wtb*v_p(i,j,k-1)

        ENDIF

       END DO
      END DO

!-----------------------------------------------------------------------
! 2. Calculate the pressure force from the wind components and density
!    at height h_m, the frontal silhouette area and surface roughness
!    length for momentum, z_0m.
!-----------------------------------------------------------------------
      DO l = 1, land_pts

         j = (land_index(l)-1)/row_length + 1
         i = land_index(l) - (j-1)*row_length

         IF ( FD_stab_dep == ON) THEN
           rib_fn=1.-rib(i,j)/ri_crit
           if(rib_fn >  1.)rib_fn=1.
           if(rib_fn <  0.)rib_fn=0.
         ELSE
           rib_fn=1.
         END IF

         zeta = log( h_m(l)/z0m(i,j) )

         IF(l_lowhill)THEN

!          Compute Wood and Mason (1993) low-hill drag expression

           tausx=(vkman/zeta)*(vkman/zeta)*u_hm(l)*                     &
     &           sqrt(u_hm(l)*u_hm(l)+v_hm(l)*v_hm(l))
           tausy=(vkman/zeta)*(vkman/zeta)*v_hm(l)*                     &
     &           sqrt(u_hm(l)*u_hm(l)+v_hm(l)*v_hm(l))
           fp_x(l) = rho_tq(i,j,1)*alpha*beta*pi*pi                     &
     &             *sil_orog_land(l)*sil_orog_land(l)                   &
     &             *rib_fn*tausx
           fp_y(l) = rho_tq(i,j,1)*alpha*beta*pi*pi                     &
     &             *sil_orog_land(l)*sil_orog_land(l)                   &
     &             *rib_fn*tausy
         ELSE

!          Compute steep hill drag expression

           tausx=u_hm(l)*sqrt(u_hm(l)*u_hm(l)+v_hm(l)*v_hm(l))
           tausy=v_hm(l)*sqrt(u_hm(l)*u_hm(l)+v_hm(l)*v_hm(l))
           fp_x(l)=0.5*rho_tq(i,j,1)*orog_drag_param*                   &
     &             sil_orog_land(l)*rib_fn*tausx
           fp_y(l)=0.5*rho_tq(i,j,1)*orog_drag_param*                   &
     &             sil_orog_land(l)*rib_fn*tausy

         ENDIF

         tau_fd_x(i,j,1) = fp_x(l)
         tau_fd_y(i,j,1) = fp_y(l)

      END DO

!-----------------------------------------------------------------------
! 3. Calculate the vertical profiles of the explicit orographic stress
!-----------------------------------------------------------------------
      DO k = 2,BL_LEVELS
       DO l = 1, land_pts

         j = (land_index(l)-1)/row_length + 1
         i = land_index(l) - (j-1)*row_length

         height_fac = exp(z_tq(i,j,k-1)/h_m(l))
         tau_fd_x(i,j,k) = tau_fd_x(i,j,1)/height_fac
         tau_fd_y(i,j,k) = tau_fd_y(i,j,1)/height_fac

       END DO
      END DO

      IF(LTIMER)THEN
! DEPENDS ON: timer
        CALL TIMER('FM_DRAG ',104)
      ENDIF

      RETURN
      END SUBROUTINE FM_DRAG
