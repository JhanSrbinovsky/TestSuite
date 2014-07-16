
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE IM_BL_PT1 ----------------------------------------------
!!!
!!!  Purpose: Calculate downward sweep of matrix for increments to
!!!           T and Q in the boundary layer, using an
!!!           implicit numerical scheme.  The tridiagonal matrices are
!!!           inverted using simple Gaussian elimination.
!!!
!!!
!!!  Model           Modification history
!!! version  Date
!!!  5.2   15/11/00   New Deck         M. Best
!!!  5.4   29/08/02   Include spherical geometry.  Adrian Lock
!!!  5.5   10/03/03   Correct index of DTRDZ_CHARNEY_GRID from K to 1.
!!!                                                Adrian Lock
!!!
!!!
!!!  Programming standard: UM Documentation Paper No 4, Version 2,
!!!                        dated 18/1/90
!!!
!!!  System component covered: P244
!!!
!!!  Project task: P24
!!!
!!!  Documentation: UM Documentation Paper No 24.
!!!
!!!---------------------------------------------------------------------
!!  Arguments :-
      SUBROUTINE IM_BL_PT1 (                                            &
     & OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,BL_LEVELS                     &
     &,halo_i, halo_j,r_rho_levels,r_theta_levels                       &
     &,DTRDZ_CHARNEY_GRID,DTRDZ_U,DTRDZ_V                               &
     &,RHOKH,RHOKM_U,RHOKM_V                                            &
     &,RDZ_CHARNEY_GRID,RDZ_U,RDZ_V,GAMMA                               &
     &,DQW_NT,DTL_NT,DU_NT,DV_NT                                        &
     &,FQW,FTL,TAU_X,TAU_Y                                              &
     &,CT_CTQ,DQW,DTL,CQ_CM_U,CQ_CM_V,DU,DV                             &
     &,LTIMER                                                           &
     &)

      IMPLICIT NONE

      LOGICAL LTIMER

      INTEGER                                                           &
     &  ROW_LENGTH                                                      &
                                 ! IN Number of points on a row
     &, ROWS                                                            &
                                 ! IN Number of rows in a theta field
     &, N_ROWS                                                          &
                                 ! IN Number of rows in a v field
     &, OFF_X                                                           &
                                 ! IN Size of small halo in i
     &, OFF_Y                                                           &
                                 ! IN Size of small halo in j.
     &, halo_i                                                          &
                                 ! IN Size of halo in i direction
     &, halo_j                                                          &
                                 ! IN Size of halo in j direction
     &,BL_LEVELS                 ! IN No. of atmospheric levels for
!                                !    which boundary layer fluxes are
!                                !    calculated.

      Real                                                              &
                                 ! IN vertical co-ordinates
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:bl_levels)                &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j, bl_levels)

      REAL                                                              &
     & DTRDZ_CHARNEY_GRID(ROW_LENGTH,ROWS,BL_LEVELS)                    &
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
     &,GAMMA(BL_LEVELS)                                                 &
                                   ! IN Implicit weighting.
     &,DQW_NT(ROW_LENGTH,ROWS,BL_LEVELS)                                &
!                                  ! IN Non-turbulent increment for QW.
     &,DTL_NT(ROW_LENGTH,ROWS,BL_LEVELS)                                &
!                                  ! IN Non-turbulent increment for TL.
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
     & CT_CTQ(ROW_LENGTH,ROWS,BL_LEVELS)                                &
!                                  ! OUT Coefficient in T and q
!                                  !     tri-diagonal implicit matrix
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
     &,AT_w(ROW_LENGTH,ROWS)                                            &
     &,RBT                                                              & 
     &,RBT_w(ROW_LENGTH,ROWS)                                           & 
                ! Reciprocal of BT' (eqns P244.107, 110, 113).
     &,DQW_w(ROW_LENGTH,ROWS)                                           &
     &,DTL_w(ROW_LENGTH,ROWS)                                           &
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
!
! EAK
      integer ktau
      save ktau
      LOGICAL, save :: LFIRST = .true.

      IF ( LFIRST ) THEN
       ktau = 0
       LFIRST = .false.
      ENDIF
      ktau = ktau  + 1

!-----------------------------------------------------------------------
!!  0.  Check that the scalars input to define the grid are consistent.
!       See comments to routine SF_EXCH for details.
!-----------------------------------------------------------------------

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('IM_BL_PT1 ',3)
      ENDIF

      BLM1 = BL_LEVELS-1

!-----------------------------------------------------------------------
!!  1.0 Interpolate r_theta_levels to U,V columns
!-----------------------------------------------------------------------

! DEPENDS ON: p_to_u
      CALL P_TO_U(r_theta_levels,ROW_LENGTH,ROWS,BL_LEVELS+1,           &
     &            halo_i, halo_j, r_theta_U)

! DEPENDS ON: p_to_v
      CALL P_TO_V(r_theta_levels,ROW_LENGTH,ROWS,n_rows,                &
     &            BL_LEVELS+1, halo_i, halo_j, r_theta_V)
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
!-----------------------------------------------------------------------
!!  3.0 Calculate matrix elements
!-----------------------------------------------------------------------

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
! Include non-turbulent increments.
        DQW(I,J,BL_LEVELS) = DTRDZ_CHARNEY_GRID(I,J,BL_LEVELS) *        &
     &                       FQW(I,J,BL_LEVELS) + DQW_NT(I,J,BL_LEVELS)
        DTL(I,J,BL_LEVELS) = DTRDZ_CHARNEY_GRID(I,J,BL_LEVELS) *        &
     &                       FTL(I,J,BL_LEVELS) + DTL_NT(I,J,BL_LEVELS)

        CT_CTQ(I,J,BL_LEVELS) = -DTRDZ_CHARNEY_GRID(I,J,BL_LEVELS) *    &
     &         GAMMA(BL_LEVELS)*RHOKH(I,J,BL_LEVELS)*                   &
     &          RDZ_CHARNEY_GRID(I,J,BL_LEVELS)

        RBT = 1.0 / ( 1.0 - CT_CTQ(I,J,BL_LEVELS) )

        DQW(I,J,BL_LEVELS) = RBT * DQW(I,J,BL_LEVELS)
        DTL(I,J,BL_LEVELS) = RBT * DTL(I,J,BL_LEVELS)

        CT_CTQ(I,J,BL_LEVELS) = RBT * CT_CTQ(I,J,BL_LEVELS)    ! P244.1
       ENDDO
      ENDDO


      DO K=BLM1,2,-1
        DO J=1,ROWS
         DO I=1,ROW_LENGTH

            DQW(I,J,K) = -DTRDZ_CHARNEY_GRID(I,J,K) *                   &
     &                  ( FQW(I,J,K+1) - FQW(I,J,K) ) + DQW_NT(I,J,K)
            DTL(I,J,K) = -DTRDZ_CHARNEY_GRID(I,J,K) *                   &
     &                  ( FTL(I,J,K+1) - FTL(I,J,K) ) + DTL_NT(I,J,K)

            AT = -DTRDZ_CHARNEY_GRID(I,J,K) *                           &
     &             GAMMA(K+1)*RHOKH(I,J,K+1)*RDZ_CHARNEY_GRID(I,J,K+1)

            CT_CTQ(I,J,K) = -DTRDZ_CHARNEY_GRID(I,J,K) *                &
     &                   GAMMA(K)*RHOKH(I,J,K)*RDZ_CHARNEY_GRID(I,J,K)

            RBT = 1.0 / ( 1.0 - CT_CTQ(I,J,K) -                         &
     &                             AT*( 1.0 + CT_CTQ(I,J,K+1) ) )

            DQW(I,J,K) = RBT * (DQW(I,J,K) - AT*DQW(I,J,K+1) )
            DTL(I,J,K) = RBT * (DTL(I,J,K) - AT*DTL(I,J,K+1) )

            CT_CTQ(I,J,K) = RBT * CT_CTQ(I,J,K)               ! P244.1
         ENDDO
        ENDDO
      ENDDO !blm1,2,-1

!-----------------------------------------------------------------------
!! 3.3 Bottom model layer QW row of matrix equation.
!-----------------------------------------------------------------------

      DO J=1,ROWS
       DO I=1,ROW_LENGTH

            DQW(I,J,1) = -DTRDZ_CHARNEY_GRID(I,J,1) * FQW(I,J,2) +      &
     &                    DQW_NT(I,J,1)
            DTL(I,J,1) = -DTRDZ_CHARNEY_GRID(I,J,1) * FTL(I,J,2) +      &
     &                    DTL_NT(I,J,1)
            DQW_w(I,J) = DQW(I,J,1)
            DTL_w(I,J) = DTL(I,J,1)
            AT = -DTRDZ_CHARNEY_GRID(I,J,1) *                           &
     &                 GAMMA(2)*RHOKH(I,J,2)*RDZ_CHARNEY_GRID(I,J,2)
            AT_w(I,J) = AT
            RBT = 1.0 / ( 1.0 - AT*( 1.0 + CT_CTQ(I,J,2) ) )
            RBT_w(I,J) = RBT
            DQW(I,J,1) = RBT * (DQW(I,J,1) - AT*DQW(I,J,2) )
            DTL(I,J,1) = RBT * (DTL(I,J,1) - AT*DTL(I,J,2) )
!
! Now set CT_CTQ(1) to be r^2 * BETA
            r_sq = r_theta_levels(i,j,0)*r_theta_levels(i,j,0)
            CT_CTQ(I,J,1) = - r_sq * DTRDZ_CHARNEY_GRID(I,J,1) * RBT

       ENDDO
      ENDDO

!      if( ktau .gt. 275) then
!      print 287,ftl(4,2,1),ftl(4,2,2),fqw(4,2,1),fqw(4,2,2),     &
!      DTRDZ_CHARNEY_GRID(4,2,1),RDZ_CHARNEY_GRID(4,2,2),         &
!      DQW_NT(4,2,1),DTL_NT(4,2,1),GAMMA(2),AT_w(4,2),RBT_w(4,2), &
!      dtl_w(4,2),dqw_w(4,2),                                     & 
!      dtl(4,2,1),dqw(4,2,1),CT_CTQ(4,2,1),CT_CTQ(4,2,2)
!      print 288,ftl(28,57,1),ftl(28,57,2),fqw(28,57,1),fqw(28,57,2), &
!      DTRDZ_CHARNEY_GRID(28,57,1),RDZ_CHARNEY_GRID(28,57,2),         &
!      DQW_NT(28,57,1),DTL_NT(28,57,1),GAMMA(2),AT_w(28,57),RBT_w(28,57), &
!      dtl_w(28,57),dqw_w(28,57),                                     &
!      dtl(28,57,1),dqw(28,57,1),CT_CTQ(28,57,1),CT_CTQ(28,57,2)
!      print 289,ftl(35,56,1),ftl(35,56,2),fqw(35,56,1),fqw(35,56,2), &
!      DTRDZ_CHARNEY_GRID(35,56,1),RDZ_CHARNEY_GRID(35,56,2),         &
!      DQW_NT(35,56,1),DTL_NT(35,56,1),GAMMA(2),AT_w(35,56),RBT_w(35,56), &
!      dtl_w(35,56),dqw_w(35,56),                                     &
!      dtl(35,56,1),dqw(35,56,1),CT_CTQ(35,56,1),CT_CTQ(35,56,2)
!287   format('Imblpt1',2(f9.6,f11.6),x,f9.6,f6.3,x,f9.6,f6.3,f4.1,2f5.1, &
!              x,2(f5.2,f9.6),2f7.1)
!288   format('Imblpt2',2(f9.6,f11.6),x,f9.6,f6.3,x,f9.6,f6.3,f4.1,2f5.1, &
!              x,2(f5.2,f9.6),2f7.1)
!289   format('Imblpt3',2(f9.6,f11.6),x,f9.6,f6.3,x,f9.6,f6.3,f4.1,2f5.1, &
!              x,2(f5.2,f9.6),2f7.1)
!      endif


      DO J=1,ROWS
       DO I=1,ROW_LENGTH

        DU(I,J,BL_LEVELS) = -DTRDZ_U(I,J,BL_LEVELS) *                   &
     &                     TAU_X(I,J,BL_LEVELS)

! addition of non-turbulent increments
        DU(I,J,BL_LEVELS) = DU(I,J,BL_LEVELS)                           &
     &                       + DU_NT(I,J,BL_LEVELS)

        CQ_CM_U(I,J,BL_LEVELS) = -DTRDZ_U(I,J,BL_LEVELS) *              &
     &        GAMMA(BL_LEVELS)*                                         &
     &        RHOKM_U(I,J,BL_LEVELS)*RDZ_U(I,J,BL_LEVELS)

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
          DU(I,J,K) = DU(I,J,K) + DU_NT(I,J,K)

          AM = -DTRDZ_U(I,J,K) * GAMMA(K+1)*RHOKM_U(I,J,K+1)*           &
     &                    RDZ_U(I,J,K+1)

          CQ_CM_U(I,J,K) = -DTRDZ_U(I,J,K) * GAMMA(K)*RHOKM_U(I,J,K)*   &
     &          RDZ_U(I,J,K)

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

! addition of non-turbulent increments
        DU(I,J,1) = DU(I,J,1) + DU_NT(I,J,1)

        AM = -DTRDZ_U(I,J,1) * GAMMA(2)*RHOKM_U(I,J,2)                  &
     &               *RDZ_U(I,J,2)

        RBM = 1.0 / ( 1.0 - AM *( 1.0 + CQ_CM_U(I,J,2) ) )

        DU(I,J,1) = RBM * ( DU(I,J,1) - AM*DU(I,J,2) )

! Now set CQ_CM_U(1) to be r^2 * BETA
        r_sq = r_theta_U(i,j,0)*r_theta_U(i,j,0)
        CQ_CM_U(I,J,1) = r_sq * DTRDZ_U(I,J,1) * RBM
       ENDDO
      ENDDO




      DO J=1,N_ROWS
       DO I=1,ROW_LENGTH

        DV(I,J,BL_LEVELS) = -DTRDZ_V(I,J,BL_LEVELS) *                   &
     &                     TAU_Y(I,J,BL_LEVELS)

! addition of non-turbulent increments
        DV(I,J,BL_LEVELS) = DV(I,J,BL_LEVELS)                           &
     &                       + DV_NT(I,J,BL_LEVELS)

        CQ_CM_V(I,J,BL_LEVELS) = -DTRDZ_V(I,J,BL_LEVELS) *              &
     &        GAMMA(BL_LEVELS)*                                         &
     &        RHOKM_V(I,J,BL_LEVELS)*RDZ_V(I,J,BL_LEVELS)

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
          DV(I,J,K) = DV(I,J,K) + DV_NT(I,J,K)

          AM = -DTRDZ_V(I,J,K) * GAMMA(K+1)*RHOKM_V(I,J,K+1)*           &
     &                    RDZ_V(I,J,K+1)

          CQ_CM_V(I,J,K) = -DTRDZ_V(I,J,K) * GAMMA(K)*RHOKM_V(I,J,K)*   &
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
        DV(I,J,1) = DV(I,J,1) + DV_NT(I,J,1)

        AM = -DTRDZ_V(I,J,1) * GAMMA(2)*RHOKM_V(I,J,2)                  &
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
       do j = 1,rows
       DO I = 1,row_length
         rr_sq          = 1.0 / (r_theta_U(i,j,k-1)*r_theta_U(i,j,k-1))
         RHOKM_U(I,J,K) = rr_sq * RHOKM_U(I,J,K)
         TAU_X(I,j,K)   = rr_sq * TAU_X(I,j,K)
       ENDDO
       ENDDO
       do j = 1,n_rows
       DO I = 1,row_length
         rr_sq          = 1.0 / (r_theta_V(i,j,k-1)*r_theta_V(i,j,k-1))
         RHOKM_V(I,J,K) = rr_sq * RHOKM_V(I,J,K)
         TAU_Y(I,j,K)   = rr_sq * TAU_Y(I,j,K)
       ENDDO
       ENDDO
      ENDDO

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('IM_BL_PT1 ',4)
      ENDIF

      RETURN
      END SUBROUTINE IM_BL_PT1
