
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE BDY_EXPL1----------------------------------------------
!!!
!!!  Purpose: Calculate explicit scalling parameters and additional
!!!           boundary layer information required by the surface
!!!           exchange scheme
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!   5.2   15/11/00  New deck.   M. Best
!!!   5.3   04/04/01  Grid-box mean buoyancy parameters passed out
!!!                   and CF passed in.     A.P.Lock
!!!   5.4   30/05/02  Remove more rapid mixing code and change comments
!!!                   to reflect spherical geometry for flux divergence
!  6.1  17/05/04  Changes to enable substepping.        M. Diamantakis.
!!!                               A.P.Lock
!!!
!!! Programming standard : unified model documentation paper No 3
!!!
!!!  Documentation: UMDP 24.
!!!
!!!---------------------------------------------------------------------

!    Arguments :-
      SUBROUTINE BDY_EXPL1 (                                            &

! IN mpp variables
     & halo_i, halo_j, off_x, off_y, row_length, rows, n_rows,          &
     & global_row_length, proc_row_group, at_extremity,                 &

!IN  Substepping information
     & Substep_Number,                                                  &

! IN values defining vertical grid of model atmosphere :
     & MODEL_DOMAIN,                                                    &
     & BL_LEVELS,                                                       &
     & r_rho_levels,r_theta_levels,                                     &
     & P,P_theta_levels, rho_rsq, rho_wet, rho_dry,                     &
     & SIN_THETA_LONGITUDE,COS_THETA_LONGITUDE,                         &

! IN U and V momentum fields.
     & U, V,                                                            &

! IN cloud data :
     & CF,Q,QCF,QCL,T,                                                  &

! IN everything not covered so far :
     & PSTAR,TIMESTEP,LQ_MIX_BL,                                        &

! OUT
     & DTRDZ_CHARNEY_GRID,RDZ_CHARNEY_GRID,DTRDZ_U,DTRDZ_V,             &
     & RDZ_U,RDZ_V,RHO_UV,RHO_TQ,RHO_DRY_TQ,DZL_charney,RDZ,            &
     & Z1_TQ,Z1_UV,Z_FULL,Z_HALF,Z_UV,Z_TQ,                             &
     & P_HALF,DELTAP,U_P,V_P,QW,TL,                                     &
     & BT,BQ,BT_CLD,BQ_CLD,BT_GB,BQ_GB,A_QS,A_DQSDT,DQSDT,              &

     & LTIMER                                                           &
     & )

      IMPLICIT NONE

!  Inputs :-

! (a) Defining horizontal grid and subset thereof to be processed.

      INTEGER                                                           &
! cjj additions - mpp variables.
     &  row_length                                                      &
                   ! Local number of points on a row
     &, rows                                                            &
                   ! Local number of rows in a theta field
     &, n_rows                                                          &
                   ! Local number of rows in a v field
     &, halo_i                                                          &
                   ! Size of halo in i direction.
     &, halo_j                                                          &
                   ! Size of halo in j direction.
     &, off_x                                                           &
                   ! Size of small halo in i
     &, off_y                                                           &
                   ! Size of small halo in j.
     &, global_row_length                                               &
                           ! number of points on a row
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &,Substep_Number

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         !  south,east or west of the processor grid

! Parameters
      Integer                                                           &
     &   PNorth,                                                        &
                      ! North processor address in the neighbor array
     &   PEast,                                                         &
                      ! East processor address in the neighbor array
     &   PSouth,                                                        &
                      ! South processor address in the neighbor array
     &   PWest,                                                         &
                      ! West processor address in the neighbor array
     &   NoDomain     ! Value in neighbor array if the domain has
                      !  no neighbor in this direction. Otherwise
                      !  the value will be the tid of the neighbor
      Parameter (                                                       &
     &   PNorth   = 1,                                                  &
     &   PEast    = 2,                                                  &
     &   PSouth   = 3,                                                  &
     &   PWest    = 4,                                                  &
     &   NoDomain = -1)

! (b) Defining vertical grid of model atmosphere.

      INTEGER                                                           &
     & BL_LEVELS                                                        &
                                   ! IN Max. no. of "boundary" levels
!                                     allowed.Assumed <= 30 for dim-
!                                     sioning of GAMMA in common deck
!                                     C_GAMMA used in SF_EXCH and KMKH
     &,MODEL_DOMAIN

      REAL                                                              &
     &  r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j,bl_levels+1)                  &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:bl_levels+1)              &
     &, COS_THETA_LONGITUDE (row_length,rows)                           &
     &, SIN_THETA_LONGITUDE (row_length,rows)                           &
     &, P(1-off_x:row_length+off_x,                                     &
     &    1-off_y:rows+off_y, BL_LEVELS+1)                              &
     &, P_theta_levels(row_length, rows, BL_LEVELS+1)                   &
     &, RHO_RSQ(1-off_x:row_length+off_x,                               &
     &          1-off_y:rows+off_y,BL_LEVELS+1)                         &
                                               ! IN Density * R**2
     &, rho_wet(row_length, rows, bl_levels+1)                          &
!                       ! IN wet density on rho levels (kg/m3)
     &, rho_dry(row_length, rows, bl_levels+1)                          &
!                       ! IN dry density on rho levels (kg/m3)
     &, U(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      BL_LEVELS)                                                  &
     &, V(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      BL_LEVELS)

! (e) Cloud data.

      REAL                                                              &
     & CF(row_length, rows, BL_LEVELS)                                  &
                                          ! IN Cloud fraction (decimal).
     &,QCF(row_length,rows,BL_LEVELS)                                   &
                                     ! IN Cloud ice (kg per kg air)
     &,QCL(row_length,rows,BL_LEVELS)                                   &
                                     ! IN Cloud liquid water
     &,Q(row_length,rows,BL_LEVELS)                                     &
                                     ! IN specific humidity
     &,T(row_length,rows,BL_LEVELS)       ! IN temperature
! Latest estimates to time level n+1 values

! (f) Atmospheric + any other data not covered so far, incl control.

      REAL                                                              &
     & PSTAR(row_length,rows)                                           &
                                    ! IN Surface pressure (Pascals).
     &,TIMESTEP                     ! IN Timestep (seconds).

      LOGICAL LQ_MIX_BL         ! IN switch for using mixing ratios

      LOGICAL LTIMER                ! Logical switch for TIMER diags


      REAL                                                              &
     & DTRDZ_CHARNEY_GRID(row_length,rows,BL_LEVELS)                    &
!                               ! OUT dt/(rho*r*r*dz) for scalar
!                               !     flux divergence
     &,RDZ_CHARNEY_GRID(row_length,rows,BL_LEVELS)                      &
!                               ! OUT RDZ(,1) is the reciprocal of the
!                                 height of level 1, i.e. of the
!                                 middle of layer 1.  For K > 1,
!                                 RDZ(,K) is the reciprocal
!                                 of the vertical distance
!                                 from level K-1 to level K.
     &,DTRDZ_U(row_length,rows,BL_LEVELS)                               &
                                           ! OUT dt/(rho*r*r*dz) for
     &,DTRDZ_V(row_length,n_rows,BL_LEVELS)                             &
                                           ! OUT U,V flux divergence
     &,RDZ_U(row_length,rows,2:BL_LEVELS)                               &
                                ! OUT 1/(Z_U(K)-Z_U(K-1)) for K > 1
     &,RDZ_V(row_length,n_rows,2:BL_LEVELS)                             &
                                ! OUT 1/(Z_V(K)-Z_V(K-1)) for K > 1
     &,RHO_UV(row_length,rows,BL_LEVELS+1)                              &
!                               ! OUT density on UV (ie. rho) levels;
!                               !    used in RHOKH so dry density if
!                               !    Lq_mix_bl is true
     &,RHO_TQ(row_length,rows,BL_LEVELS)                                &
!                               ! OUT density on TQ (ie. theta) levels;
!                               !    used in RHOKM so wet density
     &,RHO_DRY_TQ(row_length,rows,BL_LEVELS)                            &
!                               ! OUT density on TQ (ie. theta) levels;
!                               !    used in non-turb flux integration
!                               !    so dry density if Lq_mix_bl is true
     &,DZL_charney(row_length,rows,BL_LEVELS)                           &
                                ! OUT DZL(,K) is depth in m of theta
!                                 level K, i.e. distance from boundary
!                                 K-1/2 to boundary K+1/2.
     &,RDZ(row_length,rows,BL_LEVELS)                                   &
                                ! OUT RDZ(,1) is the reciprocal of the
!                                 height of level 1, i.e. of the
!                                 middle of layer 1.  For K > 1,
!                                 RDZ(,K) is the reciprocal
!                                 of the vertical distance
!                                 from level K-1 to level K.
     &,Z1_TQ(row_length,rows)                                           &
                                ! OUT Height of lowest theta level.
     &,Z1_UV(row_length,rows)                                           &
                                ! OUT Height of lowest u,v level.
     &,Z_FULL(row_length,rows,BL_LEVELS)                                &
                                ! OUT Z_FULL(*,K) is height of full
!                               ! level k.
     &,Z_HALF(row_length,rows,BL_LEVELS)                                &
                                ! OUT Z_HALF(*,K) is height of half
!                               ! level k-1/2.
     &,Z_UV(row_length,rows,BL_LEVELS)                                  &
                                ! OUT Z_UV(*,K) is height of half
!                               ! level k-1/2.
     &,Z_TQ(row_length,rows,BL_LEVELS)                                  &
                                ! OUT Z_TQ(*,K) is height of half
!                               ! level k+1/2.
!
     &,P_HALF(row_length,rows,BL_LEVELS)                                &
                                ! OUT P_HALF(*,K) is pressure at half
!                               ! level k-1/2.
     &,DELTAP(row_length,rows,BL_LEVELS)                                &
                                ! OUT Difference in pressure between
!                               ! levels
     &,U_P(row_length,rows,BL_LEVELS)                                   &
!                               ! OUT U on P-grid.
     &,V_P(row_length,rows,BL_LEVELS)                                   &
!                               ! OUT V on P-grid.
     &,QW(row_length, rows, BL_LEVELS)                                  &
!                               ! OUT Total water content
     &,TL(row_length, rows, BL_LEVELS)                                  &
!                               ! OUT Ice/liquid water temperature
     &,BT(row_length,rows,BL_LEVELS)                                    &
                                ! OUT A buoyancy parameter for clear air
!                               ! on p,T,q-levels (full levels).
     &,BQ(row_length,rows,BL_LEVELS)                                    &
                                ! OUT A buoyancy parameter for clear air
!                               ! on p,T,q-levels (full levels).
     &,BT_CLD(row_length,rows,BL_LEVELS)                                &
!                               ! OUT A buoyancy parameter for cloudy
!                               ! air on p,T,q-levels (full levels).
     &,BQ_CLD(row_length,rows,BL_LEVELS)                                &
                                ! OUT A buoyancy parameter for cloudy
!                               ! air on p,T,q-levels (full levels).
     &,BT_GB(row_length,rows,BL_LEVELS)                                 &
                                ! OUT A grid-box mean buoyancy parameter
!                               ! on p,T,q-levels (full levels).
     &,BQ_GB(row_length,rows,BL_LEVELS)                                 &
                                ! OUT A grid-box mean buoyancy parameter
!                               ! on p,T,q-levels (full levels).
     &,A_QS(row_length,rows,BL_LEVELS)                                  &
                                ! OUT Saturated lapse rate factor
!                               ! on p,T,q-levels (full levels).
     &,A_DQSDT(row_length,rows,BL_LEVELS)                               &
!                               ! OUT Saturated lapse rate factor
!                               ! on p,T,q-levels (full levels).
     &,DQSDT(row_length,rows,BL_LEVELS)
                                ! OUT Derivative of q_SAT w.r.t. T


!---------------------------------------------------------------------
!  External routines called :-

      EXTERNAL VERTICAL_DIFFS,U_TO_P,V_TO_P,Polar_vector_wind_n,        &
     & BOUY_TQ
      EXTERNAL TIMER

!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-

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
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
!---Soil layer thicknesses (m)
!---6 layers => CABLE else revert to MOSES
   REAL,PARAMETER:: DZSOIL(6) =(/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/)

! C_VKMAN start
      REAL,PARAMETER:: VKMAN=0.4 ! Von Karman's constant
! C_VKMAN end
! FLDTYPE definitions for the different field types recognised on the
! decomposition
      INTEGER,PARAMETER:: Nfld_max=7 ! maximum number of field types
      INTEGER,PARAMETER:: fld_type_p=1       ! grid on P points
      INTEGER,PARAMETER:: fld_type_u=2       ! grid on U points
      INTEGER,PARAMETER:: fld_type_v=3       ! grid on V points
      INTEGER,PARAMETER:: fld_type_comp_wave  = 4
                              ! Compressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_full_wave  = 5
                              ! Uncompressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_rim_wave   = 6
                              ! Boundary data for WAM Wave Field
      INTEGER,PARAMETER:: fld_type_r=7       ! grid on river points
      INTEGER,PARAMETER:: fld_type_unknown=-1! non-standard grid
! FLDTYPE end
! DOMTYP contains different model domain types
!
! Author : P.Burton
! History:
! Version  Date      Comment.
! 5.0      15/04/99  New comdeck
! 5.2      15/11/00  add bi_cyclic_lam domain   A. Malcolm

      INTEGER,PARAMETER:: mt_global        = 1
      INTEGER,PARAMETER:: mt_lam           = 2
      INTEGER,PARAMETER:: mt_cyclic_lam    = 3
      INTEGER,PARAMETER:: mt_bi_cyclic_lam = 4
      INTEGER,PARAMETER:: mt_single_column = 5
! DOMTYP end

! Derived local parameters.

      REAL LCRCP,LS,LSRCP

      PARAMETER (                                                       &
     & LCRCP=LC/CP                                                      &
                             ! Evaporation-to-dT conversion factor.
     &,LS=LF+LC                                                         &
                             ! Latent heat of sublimation.
     &,LSRCP=LS/CP                                                      &
                             ! Sublimation-to-dT conversion factor.
     &  )

!-----------------------------------------------------------------------

!  Workspace :-

!  Local scalars :-

      REAL                                                              &
     &  MAG_VECTOR_NP (bl_levels)                                       &
     &, DIR_VECTOR_NP (bl_levels)                                       &
     &, MAG_VECTOR_SP (bl_levels)                                       &
     &, DIR_VECTOR_SP (bl_levels)

      INTEGER                                                           &
     & I,J                                                              &
                  ! LOCAL Loop counter (horizontal field index).
     &,K          ! LOCAL Loop counter (vertical level index).

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BDY_EXPL1 ',3)
      ENDIF

!-----------------------------------------------------------------------
!! 1.  Perform calculations in what the documentation describes as
!!     subroutine Z_DZ.  In fact, a separate subroutine isn't used.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!! 1.2 Calculate layer depths and heights, and construct wind fields on
!!     P-grid.
!-----------------------------------------------------------------------

      If  ( Substep_Number  ==  1 ) Then
! DEPENDS ON: vertical_diffs
      CALL VERTICAL_DIFFS(                                              &
! IN Field dimensioning/pointers.
     &  rows, row_length, n_rows, halo_i, halo_j, off_x, off_y          &
     &, BL_LEVELS                                                       &
     &, TIMESTEP, LQ_MIX_BL                                             &
! IN Vertical coordinate information.
     &, r_rho_levels                                                    &
     &, r_theta_levels                                                  &
! IN Fields.
     &, RHO_RSQ, rho_wet, rho_dry                                       &
! OUT Vertical differences required by physics.
     &, RHO_uv, RHO_tq, RHO_DRY_tq                                      &
     &, DZL_charney, RDZ                                                &
     &, Z1_UV, Z1_TQ                                                    &
     &, RDZ_CHARNEY_GRID                                                &
     &, DTRDZ_CHARNEY_GRID                                              &
     &, dtrdz_u, dtrdz_v, rdz_u, rdz_v                                  &
     &  )

      DO K=1,BL_LEVELS
        DO J=1, rows
          DO I=1, row_length
            Z_FULL(I,j,K) = r_theta_levels(I,j,K)                       &
     &                      - r_theta_levels(I,j,0)
            Z_HALF(I,j,K) = r_rho_levels(I,j,K)                         &
     &                      - r_theta_levels(I,j,0)
            Z_UV(I,j,K) = z_half(i,j,k)
            Z_TQ(I,j,K) = z_full(i,j,k)
          ENDDO
        ENDDO
      ENDDO

! set pressure array.
      Do j = 1, rows
        DO I=1,row_length
          P_HALF(i,j,1) = pstar(i,j)
          DELTAP(i,j,1) = p(i,j,2) - pstar(i,j)
        ENDDO
      ENDDO
      DO K=2,BL_LEVELS
        Do j = 1, rows
          DO I=1,row_length
            P_HALF(i,j,k) = p(i,j,k)
            DELTAP(i,j,k) = p(i,j,k+1) - p(i,j,k)
          ENDDO
        ENDDO
      ENDDO  ! end of loop over bl_levels

      Endif ! Substep_Number  ==  1
! interpolate u and v to p grid.
! DEPENDS ON: u_to_p
      Call U_TO_P(u,row_length,rows,bl_levels,                          &
     &            off_x, off_y, model_domain,                           &
     &            at_extremity, u_p)

! DEPENDS ON: v_to_p
      Call V_TO_P(v,row_length,rows,n_rows,bl_levels,                   &
     &            off_x, off_y, model_domain,                           &
     &            at_extremity, v_p)

      IF(MODEL_DOMAIN  ==  mt_global) THEN

! Overwrite values of U_P, V_P at the poles with the magnitude of
! the vector wind.

! DEPENDS ON: polar_vector_wind_n
        Call Polar_vector_wind_n(                                       &
     &                       v,                                         &
     &                       sin_theta_longitude,                       &
     &                       cos_theta_longitude, row_length,           &
     &                       n_rows, bl_levels, mag_vector_np,          &
     &                       dir_vector_np, mag_vector_sp,              &
     &                       dir_vector_sp,                             &
     &                       off_x, off_y, global_row_length,           &
     &                       proc_row_group, at_extremity)
        If (at_extremity(PSouth) ) Then
          DO K=1,BL_LEVELS
            DO I=1,ROW_LENGTH
              V_P(I,1,K) = MAG_VECTOR_SP(k)
              U_P(I,1,K) = 0.0
            END DO
          End Do
        End If
        If (at_extremity(PNorth) ) Then
          DO K=1,BL_LEVELS
            DO I=1,ROW_LENGTH
              V_P(I,rows,K) = MAG_VECTOR_NP(k)
              U_P(I,rows,K) = 0.0
            END DO
          End Do
        End If

      ENDIF


!-----------------------------------------------------------------------
!! Calculate total water content, QW and Liquid water temperature, TL
!-----------------------------------------------------------------------

      DO K=1,BL_LEVELS
        Do j = 1, rows
          DO I=1,row_length
            QW(I,J,K) = Q(I,J,K) + QCL(I,J,K) + QCF(I,J,K)
                                ! P243.10
            TL(I,J,K) = T(I,J,K) - LCRCP*QCL(I,J,K) - LSRCP*QCF(I,J,K)
                                ! P243.9
          ENDDO
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!! 5.  Turbulent exchange coefficients and "explicit" fluxes between
!!     model layers in the boundary layer (P243b, routine KMKH).
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!! 5.1  Calculate bouyancy parameters BT and BQ.
!-----------------------------------------------------------------------
! rml 22/2/11 suggested fix from MRD for floating overflow in bouy_tq
      DQSDT=0.

! DEPENDS ON: bouy_tq
      CALL BOUY_TQ (                                                    &
     & row_length, rows, halo_i,halo_j                                  &
     &,BL_LEVELS, LQ_MIX_BL                                             &
     &,P_theta_levels,T,Q,QCF,QCL,CF                                    &
     &,BT,BQ,BT_CLD,BQ_CLD,BT_GB,BQ_GB,A_QS,A_DQSDT,DQSDT               &
     &,LTIMER                                                           &
     &  )



      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BDY_EXPL1 ',4)
      ENDIF

      RETURN
      END SUBROUTINE BDY_EXPL1
