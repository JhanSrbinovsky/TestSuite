#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!!!
!!! Subroutine BDY_IMPL3
!!!
!!!
!!!  Purpose: Calculate downward sweep of matrix for increments to
!!!           U, V, T and Q in the boundary layer for the 
!!!           unconditionally stable and non-oscillatory numerical solver
!!!
!!!  Model           Modification history
!!! version  Date
!!!  6.4   10/01/07   New Deck         M. Diamantakis
!!!
!!!  Programming standard: UMDP4
!!!
!!!  Documentation: 
!!!          http://www-nwp/~frmd/DR/Reports/new_BLsolver_guide.ps
!!!
!---------------------------------------------------------------------
!!  Arguments :-
      SUBROUTINE BDY_IMPL3 (                                            &
     & HALO_I,HALO_J,OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,BL_LEVELS       &
     &,L_CORRECT,R_RHO_LEVELS,R_THETA_LEVELS                            &
     &,Q,QCL,QCF,Q_LATEST,QCL_LATEST,QCF_LATEST,T,T_LATEST              &
     &,DTRDZ_CHARNEY_GRID,DTRDZ_U,DTRDZ_V,RHOKH,RHOKM_U,RHOKM_V         &
     &,RDZ_CHARNEY_GRID,RDZ_U,RDZ_V,GAMMA1,GAMMA2,GAMMA                 &
     &,DU_NT,DV_NT,DQW_NT,DTL_NT,FQW,FTL,TAU_X,TAU_Y                    &
     &,QW,TL,DQW1,DTL1,CT_CTQ,CTCTQ1,DQW,DTL,CQ_CM_U,CQ_CM_V,DU,DV      &
     &,LTIMER                                                           &
     &)

      IMPLICIT NONE

      LOGICAL                                                           &
     & L_correct                                                        &
     &,LTIMER

      INTEGER                                                           &
     &  ROW_LENGTH                                                      &
                                 ! IN Number of points on a row
     &, ROWS                                                            &
                                 ! IN Number of rows in a theta field
     &, N_ROWS                                                          &
                                 ! IN Number of rows in a v field
     &, HALO_I                                                          &
                                 ! IN Size of halo in i direction.
     &, HALO_J                                                          &
                                 ! IN Size of halo in j direction.
     &, OFF_X                                                           &
                                 ! IN Size of small halo in i
     &, OFF_Y                                                           &
                                 ! IN Size of small halo in j.
     &, BL_LEVELS                ! IN No. of atmospheric levels for
!                                !    which boundary layer fluxes are
!                                !    calculated.

      Real                                                              &
                                 ! IN vertical co-ordinates
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:bl_levels)                &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j, bl_levels)

      REAL                                                              &
     & GAMMA1(ROW_LENGTH,ROWS)                                          &
     &,GAMMA2(ROW_LENGTH,ROWS)                                          &
                                 ! IN new scheme weights.
     &,GAMMA(bl_levels)          ! IN standard implicit scheme weights.

      REAL                                                              &
     & Q(row_length,rows,BL_LEVELS)                                     &
!                                  ! IN specific humidity
     &,QCL(row_length,rows,BL_LEVELS)                                   &
!                                  ! IN Cloud liquid water
     &,QCF(row_length,rows,BL_LEVELS)                                   &
!                                  ! IN Cloud ice (kg per kg air)
     &,Q_latest(row_length,rows,BL_LEVELS)                              &
!                                  ! IN specific humidity
     &,QCL_latest(row_length,rows,BL_LEVELS)                            &
!                                  ! IN Cloud liquid water
     &,QCF_latest(row_length,rows,BL_LEVELS)                            &
                                   ! IN Cloud ice (kg per kg air)
     &,T(row_length,rows,BL_LEVELS)                                     &
                                   ! IN temperature
                                   !    Latest estimates to time
!                                  !    level n+1 values
     &,T_latest(row_length,rows,BL_LEVELS)                              &
!                                  ! IN temperature
     &,DTRDZ_CHARNEY_GRID(ROW_LENGTH,ROWS,BL_LEVELS)                    &
!                                  ! IN dz for bottom BL_LEVELS
     &,DTRDZ_U(ROW_LENGTH,ROWS,BL_LEVELS)                               &
!                                  ! IN -g.dt/dp for model wind layers
     &,DTRDZ_V(ROW_LENGTH,N_ROWS,BL_LEVELS)                             &
!                                  ! IN -g.dt/dp for model wind layers
     &,RHOKH(ROW_LENGTH,ROWS,2:BL_LEVELS)                               &
!                                  ! IN Exchange coeff for FTL above
!                                  !    surface.
     &,RHOKM_U(ROW_LENGTH,ROWS,2:BL_LEVELS)                             &
!                                  ! IN Exchange coefficients for
!                                  !    momentum, on U-grid with
!                                  !    first and last rows ignored.
!                                  !    for K>=2 (from KMKH).
     &,RHOKM_V(ROW_LENGTH,N_ROWS,2:BL_LEVELS)                           &
!                                  ! IN Exchange coefficients for
!                                  !    momentum, on V-grid with
!                                  !    first and last rows ignored.
!                                  !    for K>=2 (from KMKH).
     &,RDZ_CHARNEY_GRID(ROW_LENGTH,ROWS,BL_LEVELS)                      &
!                                  ! IN RDZ(,1) is the reciprocal of the
!                                  ! height of level 1, i.e. of the
!                                  ! middle of layer 1.  For K > 1,
!                                  ! RDZ(,K) is the reciprocal
!                                  ! of the vertical distance
!                                  ! from level K-1 to level K.
     &,RDZ_U(ROW_LENGTH,ROWS,2:BL_LEVELS)                               &
!                                  ! IN Reciprocal of the vertical
!                                  !    distance from level K-1 to
!                                  !    level K. (K > 1) on wind levels
     &,RDZ_V(ROW_LENGTH,N_ROWS,2:BL_LEVELS)                             &
!                                  ! IN Reciprocal of the vertical
!                                  !    distance from level K-1 to
!                                  !    level K. (K > 1) on wind levels
     &,DU_NT(1-OFF_X:ROW_LENGTH+OFF_X,                                  &
     &    1-OFF_Y:ROWS+OFF_Y,BL_LEVELS)                                 &
!                                  ! IN u non-turbulent increments.
     &,DV_NT(1-OFF_X:ROW_LENGTH+OFF_X,                                  &
     &    1-OFF_Y:N_ROWS+OFF_Y,BL_LEVELS)                               &
!                                  ! IN v non-turbulent increments.
     &,FQW(ROW_LENGTH,ROWS,BL_LEVELS)                                   &
!                                  ! IN Flux of QW (ie., for surface,
!                                  !    total evaporation). Kg/sq m/s
     &,FTL(ROW_LENGTH,ROWS,BL_LEVELS)                                   &
!                                  ! IN Flux of TL (ie., for surface,
!                                  !    H/Cp where H is sensible heat
!                                  !    in W per sq m).
     &,TAU_X(ROW_LENGTH,ROWS,BL_LEVELS)                                 &
!                                  ! IN x-component of turbulent
!                                  !    stress at levels k-1/2;
!                                  !    eg. TAUX(,1) is surface stress.
!                                  !    U-grid, 1st and last rows set
!                                  !    to "missing data". (N/sq m)
!                                  !    IN as "explicit" fluxes from
!                                  !    ex_flux_uv, OUT as "implicit
     &,TAU_Y(ROW_LENGTH,N_ROWS,BL_LEVELS)
!                                  ! IN y-component of turbulent
!                                  !    stress at levels k-1/2;
!                                  !    eg. TAUX(,1) is surface stress.
!                                  !    V-grid, 1st and last rows set
!                                  !    to "missing data". (N/sq m)
!                                  !    IN as "explicit" fluxes from
!                                  !    ex_flux_uv, OUT as "implicit


      REAL                                                              &
     & QW(row_length, rows, BL_LEVELS)                                  &
!                                  ! OUT total water
     &,TL(row_length, rows, BL_LEVELS)                                  &
                                   ! OUT liquid water temperature
     &,CT_CTQ(ROW_LENGTH,ROWS,BL_LEVELS)                                &
!                                  ! OUT Coefficient in T and q
!                                  !     tri-diagonal implicit matrix
     &,CTCTQ1(ROW_LENGTH,ROWS,BL_LEVELS)                                &
     &,CQ_CM_U(ROW_LENGTH,ROWS,BL_LEVELS)                               &
!                                  ! OUT Coefficient in U and V
!                                  !     tri-diagonal implicit matrix
     &,CQ_CM_V(ROW_LENGTH,N_ROWS,BL_LEVELS)                             &
!                                  ! OUT Coefficient in U and V
!                                  !     tri-diagonal implicit matrix
     &,DQW(ROW_LENGTH,ROWS,BL_LEVELS)                                   &
!                                  ! OUT BL increment to q field
     &,DTL(ROW_LENGTH,ROWS,BL_LEVELS)                                   &
!                                  ! OUT BL increment to T field
     &,DQW_NT(ROW_LENGTH,ROWS,BL_LEVELS)                                &
                                         ! NT incr to q field
     &,DTL_NT(ROW_LENGTH,ROWS,BL_LEVELS)                                &
                                         ! NT incr to T field
     &,DQW1(ROW_LENGTH,ROWS,bl_levels)                                  &
!                                  ! OUT 1 LEV BL increment to q field
     &,DTL1(ROW_LENGTH,ROWS,BL_LEVELS)                                  &
!                                  ! OUT 1 LEV BL increment to T field
     &,DU(1-OFF_X:ROW_LENGTH+OFF_X,                                     &
     &    1-OFF_Y:ROWS+OFF_Y,BL_LEVELS)                                 &
!                                  ! INOUT BL increment to u wind field
     &,DV(1-OFF_X:ROW_LENGTH+OFF_X,                                     &
     &    1-OFF_Y:N_ROWS+OFF_Y,BL_LEVELS)
!                                  ! INOUT BL increment to v wind field

!  External references :-
      EXTERNAL TIMER

      REAL                                                              &
     & r_theta_U(row_length,rows,0:bl_levels)                           &
                                               ! Vertical grids for U
     &,r_theta_V(row_length,n_rows,0:BL_LEVELS)! and V flux levels

!  Local scalars :-
      REAL                                                              &
     & AT                                                               &
                ! Matrix element in "T" row in eqn P244.79.
     &,RBT                                                              & 
                ! Reciprocal of BT' (eqns P244.107, 110, 113).
     &,AM                                                               &
                ! Matrix element in eqn P244.80.
     &,RBM                                                              & 
                ! Reciprocal of BM(') (eqns P244.81, 85, 89).
     &,R_SQ                                                             &
                ! square of height variables
     &,RR_SQ    ! 1/square of height variables

      INTEGER                                                           &
     & BLM1                                                             &
                ! BL_LEVELS minus 1.
     &,I,J                                                              &
                ! Loop counter (horizontal field index).
     &,K        ! Loop counter (vertical index).

#include "c_r_cp.h"
#include "c_lheat.h"

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


      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BDY_IMPL3 ',3)
      ENDIF

      IF ( L_CORRECT ) THEN

        DO K = 1,BL_LEVELS
          DO J= 1, rows
            DO I= 1, row_length
! Don't use QW, TL here as these are no longer at time level n
              DQW_NT(I,J,K) = Q_latest(I,j,K) + QCL_latest(I,j,K)       &
     &                      + QCF_latest(I,j,K)                         &
     &                      - Q(I,J,K) - QCL(I,J,K) - QCF(I,J,K)
              DTL_NT(I,J,K) = T_latest(I,J,K)                           &
     &             - LCRCP * QCL_latest(I,j,K)                          &
     &             - LSRCP * QCF_latest(I,j,K)                          &
     &             - ( T(I,J,K) - LCRCP*QCL(I,J,K) - LSRCP*QCF(I,J,K) )
            ENDDO
          ENDDO
        ENDDO
!
! Update explicit fluxes using predictor X* value as needed by the
! 2nd stage of the scheme. Note that: DTL=TL*-TL, DQW=QW*-QW etc
!
        Do K=2,BL_LEVELS
          Do J=1,ROWS
            Do I=1,ROW_LENGTH
              FTL(I,J,K) = FTL(I,J,K) - RHOKH(I,J,K) *                  &
     &          ( DTL(I,J,K) - DTL(I,J,K-1) ) * RDZ_CHARNEY_GRID(I,J,K)
              FQW(I,J,K) = FQW(I,J,K) - RHOKH(I,J,K) *                  &
     &          ( DQW(I,J,K) - DQW(I,J,K-1) ) * RDZ_CHARNEY_GRID(I,J,K)
              TAU_X(I,j,K) = TAU_X(I,J,K) + RHOKM_U(I,j,K) *            &
     &                    ( DU(I,j,K) - DU(I,j,K-1) ) *RDZ_U(I,j,K)
            ENDDO
          ENDDO
        ENDDO

        DO K=2, BL_LEVELS
          DO J=1,N_ROWS
            DO I=1,ROW_LENGTH
              TAU_Y(I,j,K) = TAU_Y(I,J,K) + RHOKM_V(I,j,K) *            &
     &                    ( DV(I,j,K) - DV(I,j,K-1) ) *RDZ_V(I,j,K)
            ENDDO
          ENDDO
        ENDDO

      ELSE

        DO K = 1,BL_LEVELS
          DO J= 1,ROWS 
            DO I= 1,ROW_LENGTH 
              QW(I,J,K) = Q(I,J,K) + QCL(I,J,K) + QCF(I,J,K)
              TL(I,J,K) = T(I,J,K) - LCRCP*QCL(I,J,K) - LSRCP*QCF(I,J,K)
              DQW_NT(I,J,K) = Q_latest(I,j,K) + QCL_latest(I,j,K)       &
     &                        + QCF_latest(I,j,K) - QW(I,J,K)
              DTL_NT(I,J,K) = T_latest(I,J,K)                           &
     &                        - LCRCP * QCL_latest(I,j,K)               &
     &                        - LSRCP * QCF_latest(I,j,K)               &
     &                        - TL(I,J,K)
            ENDDO
          ENDDO
        ENDDO

      ENDIF

      BLM1 = BL_LEVELS-1

!-----------------------------------------------------------------------
!!  1.0 Interpolate r_theta_levels to U,V columns
!-----------------------------------------------------------------------
#if !defined(SCMA)
! DEPENDS ON: p_to_u
      CALL P_TO_U(r_theta_levels,ROW_LENGTH,ROWS,BL_LEVELS+1,           &
     &            halo_i, halo_j, r_theta_U)

! DEPENDS ON: p_to_v
      CALL P_TO_V(r_theta_levels,ROW_LENGTH,ROWS,n_rows,                &
     &            BL_LEVELS+1, halo_i, halo_j, r_theta_V)
#else
      DO K=0,BL_LEVELS
        DO J= 1, rows
        DO I= 1, row_length
          r_theta_U(I,J,K) = r_theta_levels(I,J,K)
        ENDDO
        ENDDO
        DO J= 1, n_rows
        DO I= 1, row_length
          r_theta_V(I,J,K) = r_theta_levels(I,J,K)
        ENDDO
        ENDDO
      ENDDO
#endif
!-----------------------------------------------------------------------
!! 2.0 For simulations on a sphere we must use spherical geometry for
!!     vertical flux-divergences.   Thus, leaving out rho for
!!     simplicity, the standard cartesian flux-divergence:
!!          dQ(K)/dt = -(FQ(K+1)-FQ(K))/DZ
!!     becomes:
!!          dQ(K)/dt = -(r_flux(K+1)^2*FQ(K+1)-r_flux(K)^2*FQ(K))
!!                      / (r_full(K)^2 * DZ)
!!     In the code below it would be clearer to include the r^2
!!     explicitly where they occur in the equations but it is
!!     computationally cheaper to multiply the fluxes by r^2 in one go
!!     at the start and then divide out the r^2 at the end.
!-----------------------------------------------------------------------
      DO K=2,BL_LEVELS
       do j = 1,rows
       DO I = 1,row_length
         r_sq = r_rho_levels(i,j,k)*r_rho_levels(i,j,k)
         RHOKH(I,J,K) = r_sq * RHOKH(I,J,K)
         FQW(I,j,K)   = r_sq * FQW(I,j,K)
         FTL(I,j,K)   = r_sq * FTL(I,j,K)
       ENDDO
       ENDDO
      ENDDO

      DO K=2,BL_LEVELS
       do j = 1,rows
       DO I = 1,row_length
         r_sq           = r_theta_U(i,j,k-1)*r_theta_U(i,j,k-1)
         RHOKM_U(I,J,K) = r_sq * RHOKM_U(I,J,K)
         TAU_X(I,j,K)   = r_sq * TAU_X(I,j,K)
       ENDDO
       ENDDO
       do j = 1,n_rows
       DO I = 1,row_length
         r_sq           = r_theta_V(i,j,k-1)*r_theta_V(i,j,k-1)
         RHOKM_V(I,J,K) = r_sq * RHOKM_V(I,J,K)
         TAU_Y(I,j,K)   = r_sq * TAU_Y(I,j,K)
       ENDDO
       ENDDO
      ENDDO

! 
! scale surface scalar fluxes as done for levels 2,3, will be
! needed for the 2nd stage of the scheme.
!

      IF ( L_CORRECT ) THEN 
        DO J = 1, ROWS
          DO I = 1, ROW_LENGTH
            r_sq = r_rho_levels(i,j,1)*r_rho_levels(i,j,1)
            FQW(I,J,1)   = r_sq * FQW(I,J,1)
            FTL(I,J,1)   = r_sq * FTL(I,J,1)
          ENDDO
        ENDDO
      ENDIF 
!-----------------------------------------------------------------------
!!  3.0 Calculate matrix elements
!-----------------------------------------------------------------------

      DO J=1,ROWS
        DO I=1,ROW_LENGTH
! Include non-turbulent increments.
        DQW(I,J,BL_LEVELS) = ( DTRDZ_CHARNEY_GRID(I,J,BL_LEVELS) *      &
     &                       FQW(I,J,BL_LEVELS) + DQW_NT(I,J,BL_LEVELS) &
     &                       ) * GAMMA2(I,J)
        DTL(I,J,BL_LEVELS) = ( DTRDZ_CHARNEY_GRID(I,J,BL_LEVELS) *      &
     &                       FTL(I,J,BL_LEVELS) + DTL_NT(I,J,BL_LEVELS) &
     &                       ) * GAMMA2(I,J)
        CT_CTQ(I,J,BL_LEVELS) = -DTRDZ_CHARNEY_GRID(I,J,BL_LEVELS) *    &
     &         GAMMA1(I,J)*RHOKH(I,J,BL_LEVELS)*                        &
     &          RDZ_CHARNEY_GRID(I,J,BL_LEVELS)
        RBT = 1.0 / ( 1.0 - CT_CTQ(I,J,BL_LEVELS) )
        DQW(I,J,BL_LEVELS) = RBT * DQW(I,J,BL_LEVELS)
        DTL(I,J,BL_LEVELS) = RBT * DTL(I,J,BL_LEVELS)
        CT_CTQ(I,J,BL_LEVELS) = RBT * CT_CTQ(I,J,BL_LEVELS)    
        ENDDO
      ENDDO

      DO K=BLM1,2,-1
        DO J=1,ROWS
         DO I=1,ROW_LENGTH
            DQW(I,J,K) = ( -DTRDZ_CHARNEY_GRID(I,J,K)*                  &
     &            (FQW(I,J,K+1)-FQW(I,J,K))+DQW_NT(I,J,K) )*GAMMA2(I,J)
            DTL(I,J,K) = ( -DTRDZ_CHARNEY_GRID(I,J,K)*                  &
     &            (FTL(I,J,K+1)-FTL(I,J,K))+DTL_NT(I,J,K) )*GAMMA2(I,J)
            AT = -DTRDZ_CHARNEY_GRID(I,J,K) *                           &
     &            GAMMA1(I,J)*RHOKH(I,J,K+1)*RDZ_CHARNEY_GRID(I,J,K+1)
            CT_CTQ(I,J,K) = -DTRDZ_CHARNEY_GRID(I,J,K) *                &
     &            GAMMA1(I,J)*RHOKH(I,J,K)*RDZ_CHARNEY_GRID(I,J,K)
            RBT = 1.0 / ( 1.0 - CT_CTQ(I,J,K) -                         &
     &                             AT*( 1.0 + CT_CTQ(I,J,K+1) ) )
            DQW(I,J,K) = RBT * (DQW(I,J,K) - AT*DQW(I,J,K+1) )
            DTL(I,J,K) = RBT * (DTL(I,J,K) - AT*DTL(I,J,K+1) )
            CT_CTQ(I,J,K) = RBT * CT_CTQ(I,J,K)               
         ENDDO
        ENDDO
      ENDDO !blm1,2,-1
!
!-----------------------------------------------------------------------
!  Bottom model layer QW row of matrix equation.
!-----------------------------------------------------------------------
!
      IF ( .NOT. L_CORRECT ) THEN
!
!-----------------------------------------------------------------------
! The following calculations are only done on the 1st stage (predictor).
! Their purpose is to compute the surface scalar (T, Q) increments which 
! are needed by the surface scheme to compute the implicit scalar fluxes.
! The same implicit scalar fluxes are used as a boundary condition for
! discrete equations of the 2nd stage. 
! Due to the dependency to the surface scalar fluxes, the 1st stage 
! downward sweep remains incomplete. It is completed at the 
! beginning of bdy_impl4, since there the surface scalar fluxes
! have been fully updated. This is done by sf_impl2().
! NOTE: The standard scheme solver is used for this calculation. 
!       Incorporation of the new scheme for the scalar surface variables 
!       would have been a more preferable choice but currently not 
!       available.
!-----------------------------------------------------------------------
!
        DO J=1,ROWS
          DO I=1,ROW_LENGTH
! Include non-turbulent increments.
            DQW1(I,J,BL_LEVELS) = DTRDZ_CHARNEY_GRID(I,J,BL_LEVELS)*    &
     &                       FQW(I,J,BL_LEVELS) + DQW_NT(I,J,BL_LEVELS)
            DTL1(I,J,BL_LEVELS) = DTRDZ_CHARNEY_GRID(I,J,BL_LEVELS)*    &
     &                       FTL(I,J,BL_LEVELS) + DTL_NT(I,J,BL_LEVELS)
            CTCTQ1(I,J,BL_LEVELS) = -DTRDZ_CHARNEY_GRID(I,J,BL_LEVELS)* &
     &         GAMMA(BL_LEVELS)*RHOKH(I,J,BL_LEVELS)*                   &
     &         RDZ_CHARNEY_GRID(I,J,BL_LEVELS)
            RBT = 1.0 / ( 1.0 - CTCTQ1(I,J,BL_LEVELS) )
            DQW1(I,J,BL_LEVELS) = RBT * DQW(I,J,BL_LEVELS)
            DTL1(I,J,BL_LEVELS) = RBT * DTL(I,J,BL_LEVELS)
            CTCTQ1(I,J,BL_LEVELS) = RBT * CTCTQ1(I,J,BL_LEVELS)   
          ENDDO
        ENDDO

        DO K=BLM1,2,-1
          DO J=1,ROWS
            DO I=1,ROW_LENGTH
              DQW1(I,J,K) = -DTRDZ_CHARNEY_GRID(I,J,K) *                &
     &                  ( FQW(I,J,K+1) - FQW(I,J,K) ) + DQW_NT(I,J,K)
              DTL1(I,J,K) = -DTRDZ_CHARNEY_GRID(I,J,K) *                &
     &                  ( FTL(I,J,K+1) - FTL(I,J,K) ) + DTL_NT(I,J,K)
              AT = -DTRDZ_CHARNEY_GRID(I,J,K) *                         &
     &             GAMMA(K+1)*RHOKH(I,J,K+1)*RDZ_CHARNEY_GRID(I,J,K+1)
              CTCTQ1(I,J,K) = -DTRDZ_CHARNEY_GRID(I,J,K) *              &
     &                   GAMMA(K)*RHOKH(I,J,K)*RDZ_CHARNEY_GRID(I,J,K)
              RBT = 1.0 / ( 1.0 - CTCTQ1(I,J,K) -                       &
     &                             AT*( 1.0 + CTCTQ1(I,J,K+1) ) )
              DQW1(I,J,K) = RBT * (DQW1(I,J,K) - AT*DQW1(I,J,K+1) )
              DTL1(I,J,K) = RBT * (DTL1(I,J,K) - AT*DTL1(I,J,K+1) )
              CTCTQ1(I,J,K) = RBT * CTCTQ1(I,J,K)               
            ENDDO
          ENDDO
        ENDDO !blm1,2,-1

        DO J=1,ROWS
          DO I=1,ROW_LENGTH
            DQW1(I,J,1) = -DTRDZ_CHARNEY_GRID(I,J,1) * FQW(I,J,2) +     &
     &                    DQW_NT(I,J,1)
            DTL1(I,J,1) = -DTRDZ_CHARNEY_GRID(I,J,1) * FTL(I,J,2) +     &
     &                    DTL_NT(I,J,1)
            AT = -DTRDZ_CHARNEY_GRID(I,J,1) *                           &
     &                 GAMMA(2)*RHOKH(I,J,2)*RDZ_CHARNEY_GRID(I,J,2)
            RBT = 1.0 / ( 1.0 - AT*( 1.0 + CTCTQ1(I,J,2) ) )
            DQW1(I,J,1) = RBT * (DQW1(I,J,1) - AT*DQW1(I,J,2) )
            DTL1(I,J,1) = RBT * (DTL1(I,J,1) - AT*DTL1(I,J,2) )
!
! Now set CT_CTQ(1) to be r^2 * BETA
            r_sq = r_theta_levels(i,j,0)*r_theta_levels(i,j,0)
            CTCTQ1(I,J,1) = - r_sq * DTRDZ_CHARNEY_GRID(I,J,1) * RBT
          ENDDO
        ENDDO

      ELSE
!-----------------------------------------------------------------------
!
! The following calculations complete the downward sweep for  the surface 
! scalar variables. They apply on the 2nd stage of the scheme.
! The equivalent calculations for the 1st stage are done at the 
! beginning of bdy_impl4, since at this stage (bdy_impl3) the surface 
! scalar fluxes are not fully updated. This will be done by sf_impl2().
!
!-----------------------------------------------------------------------
!
        DO J=1,ROWS
          DO I=1,ROW_LENGTH
            DQW(I,J,1) = GAMMA2(I,J) * ( -DTRDZ_CHARNEY_GRID(I,J,1) *   &
     &          ( FQW(I,J,2) - FQW(I,J,1) ) + DQW_NT(I,J,1) )
            DTL(I,J,1) = GAMMA2(I,J) * ( -DTRDZ_CHARNEY_GRID(I,J,1) *   &
     &          ( FTL(I,J,2) - FTL(I,J,1) ) + DTL_NT(I,J,1) )
            AT = -DTRDZ_CHARNEY_GRID(I,J,1) *                           &
     &            GAMMA1(I,J)*RHOKH(I,J,2)*RDZ_CHARNEY_GRID(I,J,2)
            RBT = 1.0 / ( 1.0 - AT*( 1.0 + CT_CTQ(I,J,2) ) )
            DQW(I,J,1) = RBT * (DQW(I,J,1) - AT*DQW(I,J,2) )
            DTL(I,J,1) = RBT * (DTL(I,J,1) - AT*DTL(I,J,2) )
!
! Now set CT_CTQ(1) to be r^2 * BETA
            r_sq = r_theta_levels(i,j,0)*r_theta_levels(i,j,0)
            CT_CTQ(I,J,1) = - r_sq * DTRDZ_CHARNEY_GRID(I,J,1) * RBT
          ENDDO
        ENDDO

      ENDIF

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        DU(I,J,BL_LEVELS) = -DTRDZ_U(I,J,BL_LEVELS)*TAU_X(I,J,BL_LEVELS)
!
! addition of non-turbulent increments
!
        DU(I,J,BL_LEVELS) = GAMMA2(I,J) * ( DU(I,J,BL_LEVELS)           &
     &                       + DU_NT(I,J,BL_LEVELS) )
        CQ_CM_U(I,J,BL_LEVELS) = -DTRDZ_U(I,J,BL_LEVELS) *              &
     &        GAMMA1(I,J)*RHOKM_U(I,J,BL_LEVELS)*RDZ_U(I,J,BL_LEVELS)
        RBM = 1.0 / ( 1.0 - CQ_CM_U(I,J,BL_LEVELS) )
        DU(I,J,BL_LEVELS) = RBM * DU(I,J,BL_LEVELS)
        CQ_CM_U(I,J,BL_LEVELS) = RBM * CQ_CM_U(I,J,BL_LEVELS)
       ENDDO
      ENDDO

      DO K=BLM1,2,-1
        DO J=1,ROWS
         DO I=1,ROW_LENGTH
          DU(I,J,K) = DTRDZ_U(I,J,K) *                                  &
     &                   ( TAU_X(I,J,K+1) - TAU_X(I,J,K) )
! addition of non-turbulent increments
          DU(I,J,K) = GAMMA2(I,J) * (DU(I,J,K) + DU_NT(I,J,K))
          AM = -DTRDZ_U(I,J,K) * GAMMA1(I,J)*RHOKM_U(I,J,K+1)*          &
     &                    RDZ_U(I,J,K+1)
          CQ_CM_U(I,J,K) = -DTRDZ_U(I,J,K)*GAMMA1(I,J)*RHOKM_U(I,J,K)*  &
     &                      RDZ_U(I,J,K)
          RBM = 1.0 / ( 1.0 - CQ_CM_U(I,J,K) -                          &
     &                    AM*( 1.0 + CQ_CM_U(I,J,K+1) ) )
          DU(I,J,K) = RBM * ( DU(I,J,K) - AM*DU(I,J,K+1) )
          CQ_CM_U(I,J,K) = RBM * CQ_CM_U(I,J,K)
         ENDDO
        ENDDO
      ENDDO ! loop over 2,BLM1

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        DU(I,J,1) = DTRDZ_U(I,J,1) * TAU_X(I,J,2)
!
! addition of non-turbulent increments
!
        DU(I,J,1) = GAMMA2(I,J)*(DU(I,J,1) + DU_NT(I,J,1))
        AM = -DTRDZ_U(I,J,1) * GAMMA1(I,J)*RHOKM_U(I,J,2)               &
     &               *RDZ_U(I,J,2)
        RBM = 1.0 / ( 1.0 - AM *( 1.0 + CQ_CM_U(I,J,2) ) )
        DU(I,J,1) = RBM * ( DU(I,J,1) - AM*DU(I,J,2) )
!
! Now set CQ_CM_U(1) to be r^2 * BETA
!
        r_sq = r_theta_U(i,j,0)*r_theta_U(i,j,0)
        CQ_CM_U(I,J,1) = r_sq * DTRDZ_U(I,J,1) * RBM
       ENDDO
      ENDDO


      DO J=1,N_ROWS
       DO I=1,ROW_LENGTH
        DV(I,J,BL_LEVELS) = -DTRDZ_V(I,J,BL_LEVELS) *                   &
     &                     TAU_Y(I,J,BL_LEVELS)
! addition of non-turbulent increments
        DV(I,J,BL_LEVELS) = GAMMA2(I,J) * ( DV(I,J,BL_LEVELS)           &
     &                       + DV_NT(I,J,BL_LEVELS) )
        CQ_CM_V(I,J,BL_LEVELS) = -DTRDZ_V(I,J,BL_LEVELS) *              &
     &        GAMMA1(I,J)*RHOKM_V(I,J,BL_LEVELS)*RDZ_V(I,J,BL_LEVELS)
        RBM = 1.0 / ( 1.0 - CQ_CM_V(I,J,BL_LEVELS) )
        DV(I,J,BL_LEVELS) = RBM * DV(I,J,BL_LEVELS)
        CQ_CM_V(I,J,BL_LEVELS) = RBM * CQ_CM_V(I,J,BL_LEVELS)
       ENDDO
      ENDDO


      DO K=BLM1,2,-1
        DO J=1,N_ROWS
         DO I=1,ROW_LENGTH
          DV(I,J,K) = DTRDZ_V(I,J,K) *                                  &
     &                   ( TAU_Y(I,J,K+1) - TAU_Y(I,J,K) )
! addition of non-turbulent increments
          DV(I,J,K) = GAMMA2(I,J) * (DV(I,J,K) + DV_NT(I,J,K))
          AM = -DTRDZ_V(I,J,K) * GAMMA1(I,J)*RHOKM_V(I,J,K+1)*          &
     &                    RDZ_V(I,J,K+1)
          CQ_CM_V(I,J,K) = -DTRDZ_V(I,J,K)*GAMMA1(I,J)*RHOKM_V(I,J,K)*  &
     &          RDZ_V(I,J,K)
          RBM = 1.0 / ( 1.0 - CQ_CM_V(I,J,K) -                          &
     &                    AM*( 1.0 + CQ_CM_V(I,J,K+1) ) )
          DV(I,J,K) = RBM * ( DV(I,J,K) - AM*DV(I,J,K+1) )
          CQ_CM_V(I,J,K) = RBM * CQ_CM_V(I,J,K)
         ENDDO
        ENDDO
      ENDDO ! loop over 2,BLM1


      DO J=1,N_ROWS
       DO I=1,ROW_LENGTH
        DV(I,J,1) = DTRDZ_V(I,J,1) * TAU_Y(I,J,2)
! addition of non-turbulent increments
        DV(I,J,1) = GAMMA2(I,J) * (DV(I,J,1) + DV_NT(I,J,1))
        AM = -DTRDZ_V(I,J,1) * GAMMA1(I,J)*RHOKM_V(I,J,2)               &
     &               *RDZ_V(I,J,2)
        RBM = 1.0 / ( 1.0 - AM *( 1.0 + CQ_CM_V(I,J,2) ) )
        DV(I,J,1) = RBM * ( DV(I,J,1) - AM*DV(I,J,2) )
! Now set CQ_CM_V(1) to be r^2 * BETA
        r_sq = r_theta_V(i,j,0)*r_theta_V(i,j,0)
        CQ_CM_V(I,J,1) = r_sq * DTRDZ_V(I,J,1) * RBM
       ENDDO
      ENDDO

!-----------------------------------------------------------------------
!! 4.0 Return fluxes back to their true value by dividing by r*r
!-----------------------------------------------------------------------
      DO K=2,BL_LEVELS
       do j = 1,rows
       DO I = 1,row_length
         rr_sq = 1.0 / ( r_rho_levels(i,j,k)*r_rho_levels(i,j,k) )
         RHOKH(I,J,K) = RHOKH(I,J,K) * rr_sq
         FQW(I,j,K)   = FQW(I,j,K) * rr_sq
         FTL(I,j,K)   = FTL(I,j,K) * rr_sq
       ENDDO
       ENDDO
      ENDDO

      DO K=2,BL_LEVELS
        DO J = 1,ROWS
          DO I = 1,ROW_LENGTH
           rr_sq          = 1.0/(r_theta_U(i,j,k-1)*r_theta_U(i,j,k-1))
           RHOKM_U(I,J,K) = rr_sq * RHOKM_U(I,J,K)
           TAU_X(I,j,K)   = rr_sq * TAU_X(I,j,K)
          ENDDO
        ENDDO

        DO J = 1,N_ROWS
          DO I = 1,ROW_LENGTH
           rr_sq          = 1.0/(r_theta_V(i,j,k-1)*r_theta_V(i,j,k-1))
           RHOKM_V(I,J,K) = rr_sq * RHOKM_V(I,J,K)
           TAU_Y(I,j,K)   = rr_sq * TAU_Y(I,j,K)
          ENDDO
        ENDDO
      ENDDO

      IF ( L_CORRECT ) THEN
        DO J = 1, rows
          DO I = 1, row_length
            rr_sq = 1.0 / ( r_rho_levels(i,j,1)*r_rho_levels(i,j,1) )
            FQW(I,j,1) = FQW(I,j,1) * rr_sq
            FTL(I,j,1) = FTL(I,j,1) * rr_sq
          ENDDO
        ENDDO
      ENDIF

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BDY_IMPL3 ',4)
      ENDIF

      RETURN
      END SUBROUTINE BDY_IMPL3
#endif
