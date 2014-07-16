#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE BDY_IMPL1 ----------------------------------------------
!!!
!!!  Purpose: Calculate downward sweep of matrix for increments to
!!!           T and Q in the boundary layer, by calling subroutine
!!!           IM_BL_PT1.
!!!
!!!
!!!  Model           Modification history
!!! version  Date
!LL  5.2    15/11/00  New Deck     M. Best
!!!  5.4    15/05/02  Pass through heights arrays for implicit
!!!                   solver on a sphere   Adrian Lock
!!!  5.5    19/02/03  Remove redundent L_BL_LSPICE.  Adrian Lock
!  6.1  17/05/04  Changes to the BL solver to enable phys2 substepping.
!                                                       M. Diamantakis
!!!
!!!  JJ  <- Programmers of some or all of previous code or changes
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
      SUBROUTINE BDY_IMPL1 (                                            &
     & HALO_I,HALO_J,OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,BL_LEVELS       &
     &,L_phys2_substep                                                  &
     &,r_rho_levels,r_theta_levels                                      &
     &,Q,QCL,QCF,Q_LATEST,QCL_LATEST,QCF_LATEST                         &
     &,T,T_LATEST                                                       &
     &,DTRDZ_CHARNEY_GRID,DTRDZ_U,DTRDZ_V                               &
     &,RHOKH,RHOKM_U,RHOKM_V                                            &
     &,RDZ_CHARNEY_GRID,RDZ_U,RDZ_V,GAMMA                               &
     &,DU_NT,DV_NT                                                      &
     &,FQW,FTL,TAU_X,TAU_Y                                              &
     &,QW,TL                                                            &
     &,CT_CTQ,DQW,DTL,CQ_CM_U,CQ_CM_V,DU,DV                             &
     &,LTIMER                                                           &
     &)

      IMPLICIT NONE

      LOGICAL                                                           &
     & L_phys2_substep                                                  &
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
     &,RHOKH(ROW_LENGTH,ROWS,BL_LEVELS)                                 &
!                                  ! IN Exchange coeff for FTL above
!                                  !    surface.
     &,RHOKM_U(ROW_LENGTH,ROWS,BL_LEVELS)                               &
!                                  ! IN Exchange coefficients for
!                                  !    momentum, on U-grid with
!                                  !    first and last rows ignored.
!                                  !    for K>=2 (from KMKH).
     &,RHOKM_V(ROW_LENGTH,N_ROWS,BL_LEVELS)                             &
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

! Workspace
      Real                                                              &
     & DQW_NT(ROW_LENGTH,ROWS,BL_LEVELS)                                &
     &,DTL_NT(ROW_LENGTH,ROWS,BL_LEVELS)

!  External references :-
      EXTERNAL IM_BL_PT1
      EXTERNAL TIMER

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


!  Local scalars :-

      INTEGER                                                           &
     & I,J                                                              &
                ! Loop counter (horizontal field index).
     &,K        ! Loop counter (vertical index).


      integer ktau
      save ktau
      LOGICAL, save :: LFIRST = .true.


      IF ( LFIRST ) THEN
       ktau = 0
       LFIRST = .false.
      ENDIF
      ktau = ktau  + 1

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BDY_IMPL1 ',3)
      ENDIF

        DO K = 1,BL_LEVELS
          DO J= 1, rows
            DO I= 1, row_length
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

!     &,Q,QCL,QCF,Q_LATEST,QCL_LATEST,QCF_LATEST                         &
!     &,T,T_LATEST                                                       &
!     &,DTRDZ_CHARNEY_GRID,DTRDZ_U,DTRDZ_V                               &
!     &,RHOKH,RHOKM_U,RHOKM_V                                            &
!     &,RDZ_CHARNEY_GRID,RDZ_U,RDZ_V,GAMMA                               &
!     &,DU_NT,DV_NT                                                      &
!     &,FQW,FTL,TAU_X,TAU_Y                                              &
!     &,QW,TL                                                            &
!     &,CT_CTQ,DQW,DTL,CQ_CM_U,CQ_CM_V,DU,DV      

!      if( ktau .gt. 275) then
!      print 287,ftl(4,2,1),fqw(4,2,1),TAU_X(4,2,1),TAU_Y(4,2,1),   &
!      Q(4,2,1),Q_latest(4,2,1),T(4,2,1),T_latest(4,2,1),TL(4,2,1), &
!      QW(4,2,1),QCL_latest(4,2,1),QCF_latest(4,2,1),               &
!      RHOKH(4,2,1),RHOKM_U(4,2,1),                                 &
!      dtl(4,2,1),dqw(4,2,1),dtl_nt(4,2,1),dqw_nt(4,2,1) 
!      print 288,ftl(28,57,1),fqw(28,57,1),TAU_X(28,57,1),TAU_Y(28,57,1), &
!      Q(28,57,1),Q_latest(28,57,1),T(28,57,1),T_latest(28,57,1),  &
!      TL(28,57,1),                                                &
!      QW(28,57,1),QCL_latest(28,57,1),QCF_latest(28,57,1),        &
!      RHOKH(28,57,1),RHOKM_U(28,57,1),                            &
!      dtl(28,57,1),dqw(28,57,1),dtl_nt(28,57,1),dqw_nt(28,57,1)
!      print 289,ftl(35,56,1),fqw(35,56,1),TAU_X(35,56,1),TAU_Y(35,56,1), &
!      Q(35,56,1),Q_latest(35,56,1),T(35,56,1),T_latest(35,56,1), &
!      TL(35,56,1),                                               &
!      QW(35,56,1),QCL_latest(35,56,1),QCF_latest(35,56,1),       &
!      RHOKH(35,56,1),RHOKM_U(35,56,1),                           &
!      dtl(35,56,1),dqw(35,56,1),dtl_nt(35,56,1),dqw_nt(35,56,1)
!287   format('bdyipt1',f9.6,f10.7,2f7.3,2f9.6,x,3f6.1,3f9.6,x,2f6.3, &
!      2(f6.3,f10.7))
!288   format('bdyipt2',f9.6,f10.7,2f7.3,2f9.6,x,3f6.1,3f9.6,x,2f6.3, &
!      2(f6.3,f10.7))
!289   format('bdyipt3',f9.6,f10.7,2f7.3,2f9.6,x,3f6.1,3f9.6,x,2f6.3, &
!      2(f6.3,f10.7))
!      endif
!      if( ktau .gt. 1580) print *,'QCF_latest',ktau,QCF_latest(:,:,1)
!      if( ktau .gt. 1580) print *,'QCF',ktau,QCF(:,:,1)
!      if( ktau .gt. 1580) print *,'QCL_latest',ktau,QCL_latest(:,:,1)
!      if( ktau .gt. 1580) print *,'QCL',ktau,QCL(:,:,1)
!      if( ktau .gt. 1580) print *,'Q_latest',ktau,Q_latest(:,:,1)
!      if( ktau .gt. 1580) print *,'Q',ktau,Q(:,:,1)
!      if( ktau .gt. 1580) print *,'T_latest',ktau,T_latest(:,:,1)
!      if( ktau .gt. 1580) print *,'T',ktau,T(:,:,1)
!      if( ktau .gt. 1580) print *,'DU',ktau,DU(:,:,1)
!      if( ktau .gt. 1580) print *,'DU_NT',ktau,DU_NT(:,:,1)
!      if( ktau .gt. 1580) print *,'DV',ktau,DV(:,:,1)
!      if( ktau .gt. 1580) print *,'DV_NT',ktau,DV_NT(:,:,1)


! DEPENDS ON: im_bl_pt1
      CALL IM_BL_PT1 (                                                  &
     & OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,BL_LEVELS                     &
     &,halo_i,halo_j,r_rho_levels,r_theta_levels                        &
     &,DTRDZ_CHARNEY_GRID,DTRDZ_U,DTRDZ_V                               &
     &,RHOKH(1,1,2),RHOKM_U(1,1,2),RHOKM_V(1,1,2)                       &
     &,RDZ_CHARNEY_GRID,RDZ_U,RDZ_V,GAMMA                               &
     &,DQW_NT,DTL_NT,DU_NT,DV_NT                                        &
     &,FQW,FTL,TAU_X,TAU_Y                                              &
     &,CT_CTQ,DQW,DTL,CQ_CM_U,CQ_CM_V,DU,DV                             &
     &,LTIMER                                                           &
     &)


      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BDY_IMPL1 ',4)
      ENDIF

      RETURN
      END SUBROUTINE BDY_IMPL1
#endif
