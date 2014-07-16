#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!!! Subroutine BDY_IMPL4
!!!
!!!  Purpose: Calculate implicit correction to boundary layer fluxes of
!!!           heat, moisture and momentum for the unconditionally
!!!           stable and non-oscillatory numerical solver. 
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!  6.4 10/01/07     New deck introduced.   M. Diamantakis
!!!
!!!  Programming standard: UMDP4
!!!
!!!  Documentation: 
!!!          http://www-nwp/~frmd/DR/Reports/new_BLsolver_guide.ps
!!!
!!!---------------------------------------------------------------------
!!  Arguments :-
      SUBROUTINE BDY_IMPL4 (                                            &

! IN values defining field dimensions and subset to be processed :
     & HALO_I,HALO_J,OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,                &

! IN values defining vertical grid of model atmosphere :
     & BL_LEVELS,                                                       &

! IN Substepping information
     & L_CORRECT, L_FTL, L_FQW, L_TAUX, L_TAUY,                         &

! IN data :
     & GAMMA1,GAMMA2,R_RHO_LEVELS,RHOKH,RHOKM_U,RHOKM_V,                &
     & RDZ_CHARNEY_GRID,DTRDZ_CHARNEY_GRID,RDZ_U,RDZ_V,                 &

! INOUT data :
     & QW,TL,FQW,FTL,TAU_X,TAU_Y,                                       &
     & FQW_star,FTL_star,TAUX_star,TAUY_star,DU,DV,DU_STAR,DV_STAR,     &
     & CT_CTQ,CTCTQ1,DQW,DTL,DQW_NT,DTL_NT,CQ_CM_U,CQ_CM_V,             &

! OUT data:
     & T_LATEST,Q_LATEST,RHOKH_MIX,                                     &

! LOGICAL LTIMER
     & LTIMER                                                           &
     & )

      IMPLICIT NONE

!  Inputs :-

! (a) Defining horizontal grid and subset thereof to be processed.

      INTEGER                                                           &
     &  ROW_LENGTH                                                      &
                   ! Local number of points on a row
     &, ROWS                                                            &
                   ! Local number of rows in a theta field
     &, N_ROWS                                                          &
                   ! Local number of rows in a v field
     &, OFF_X                                                           &
                   ! Size of small halo in i
     &, OFF_Y                                                           &
                   ! Size of small halo in j
     &, HALO_I                                                          &
                   ! Size of halo in i
     &, HALO_J     
                   ! Size of halo in j

! (b) Defining vertical grid of model atmosphere.

      INTEGER                                                           &
     & BL_LEVELS                   ! IN Max. no. of "boundary" levels
!                                     allowed.
!  In :-

      Logical                                                           &
     & L_FTL                                                            &
     &,L_FQW                                                            &
     &,L_TAUX                                                           &
     &,L_TAUY                                                           &
     &,L_CORRECT


      REAL                                                              &
     & RHOKH(ROW_LENGTH,ROWS,2:BL_LEVELS)                               &
!                                  ! IN Exchange coeffs for moisture.
     &,RHOKM_U(ROW_LENGTH,ROWS,2:BL_LEVELS)                             &
!                                  ! IN Exchange coefficients for U
     &,RHOKM_V(ROW_LENGTH,N_ROWS,2:BL_LEVELS)                           &
!                                  ! IN Exchange coefficients for V
     &,RDZ_CHARNEY_GRID(ROW_LENGTH,ROWS,BL_LEVELS)                      &
!                                  ! IN RDZ(,1) is the reciprocal of the
!                                  ! height of level 1, i.e. of the
!                                  ! middle of layer 1.  For K > 1,
!                                  ! RDZ(,K) is the reciprocal
!                                  ! of the vertical distance
!                                  ! from level K-1 to level K.
     &,DTRDZ_CHARNEY_GRID(ROW_LENGTH,ROWS,BL_LEVELS)                    &
     &,R_RHO_LEVELS(1-halo_i:row_length+halo_i,                         &
     &               1-halo_j:rows+halo_j,BL_LEVELS)                    &
     &,RDZ_U(ROW_LENGTH,ROWS,2:BL_LEVELS)                               &
!                                  ! IN  RDZ (K > 1) on U-grid.
     &,RDZ_V(ROW_LENGTH,N_ROWS,2:BL_LEVELS)                             &
!                                  ! IN  RDZ (K > 1) on V-grid.
     &,GAMMA1(ROW_LENGTH,ROWS)                                          &
     &,GAMMA2(ROW_LENGTH,ROWS)                                          
                                   ! IN new scheme weights.

      LOGICAL LTIMER               ! Logical switch for TIMER diags

!  In/outs :-

      REAL                                                              &
     & QW(ROW_LENGTH,ROWS,BL_LEVELS)                                    &
!                                  ! INOUT Total water content, but
!                                  !       replaced by specific
!                                  !       humidity in LS_CLD.
     &,TL(ROW_LENGTH,ROWS,BL_LEVELS)                                    &
!                                  ! INOUT Ice/liquid water temperature,
!                                  !       but replaced by T in LS_CLD.
     &,FQW(ROW_LENGTH,ROWS,BL_LEVELS)                                   &
!                                  ! INOUT Moisture flux between layers
!                                  !       (kg per square metre per sec)
!                                  !       FQW(,1) is total water flux
!                                  !       from surface, 'E'.
     &,FTL(ROW_LENGTH,ROWS,BL_LEVELS)                                   &
!                                  ! INOUT FTL(,K) contains net
!                                  !       turbulent sensible heat flux
!                                  !       into layer K from below; so
!                                  !       FTL(,1) is the surface
!                                  !       sensible heat, H. (W/m2)
     &,TAU_X(ROW_LENGTH,ROWS,BL_LEVELS)                                 &
!                                  ! INOUT W'ly component of surface
!                                  !       wind stress (N/sq m).(On
!                                  !       UV-grid with first and last
!                                  !       rows undefined or at present,
!                                  !       set to  missing data
     &,TAU_Y(ROW_LENGTH,N_ROWS,BL_LEVELS)                               &
!                                  ! INOUT S'ly component of surface
!                                  !       wind stress (N/sq m).  On
!                                  !       UV-grid; comments as per TAUX
!                                  4 arrays below:
                                   ! INOUT temp arrays for diags
     &,FQW_star(row_length,rows,BL_LEVELS)                              &
     &,FTL_star(row_length,rows,BL_LEVELS)                              &
     &,TAUX_star(row_length,rows,BL_LEVELS)                             &
     &,TAUY_star(row_length,n_rows,BL_LEVELS)                           &
     &,CT_CTQ(ROW_LENGTH,ROWS,BL_LEVELS)                                &
!                                  ! INOUT Coefficient in T and q
!                                  !       tri-diagonal implicit matrix
     &,CTCTQ1(ROW_LENGTH,ROWS,BL_LEVELS)                                &
     &,CQ_CM_U(ROW_LENGTH,ROWS,BL_LEVELS)                               &
!                                  ! INOUT Coefficient in U tri-diagonal
!                                  !       implicit matrix
     &,CQ_CM_V(ROW_LENGTH,N_ROWS,BL_LEVELS)                             &
!                                  ! INOUT Coefficient in V tri-diagonal
!                                  !       implicit matrix
     &,DQW(ROW_LENGTH,ROWS,BL_LEVELS)                                   &
!                                  ! INOUT BL increment to q field
     &,DTL(ROW_LENGTH,ROWS,BL_LEVELS)                                   &
!                                  ! INOUT BL increment to T field
     &,DU(1-OFF_X:ROW_LENGTH+OFF_X,                                     &
     &    1-OFF_Y:ROWS+OFF_Y,BL_LEVELS)                                 &
!                                  ! INOUT BL increment to u wind field
     &,DV(1-OFF_X:ROW_LENGTH+OFF_X,                                     &
     &    1-OFF_Y:N_ROWS+OFF_Y,BL_LEVELS)                               &
!                                  ! INOUT BL increment to v wind field
     &,DU_STAR(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         BL_LEVELS)                                               &
                                        ! OUT BL incr to u wind field
     &,DV_STAR(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,          &
     &         BL_LEVELS)                                               &
     &,DQW_NT(ROW_LENGTH,ROWS,BL_LEVELS)                                &
                                        ! IN NT incr to qw
     &,DTL_NT(ROW_LENGTH,ROWS,BL_LEVELS)                                &
                                        ! IN NT incr to TL
     &,Q_latest(row_length,rows,BL_LEVELS)                              &
                                        ! OUT specific humidity
     &,T_latest(row_length,rows,BL_LEVELS)                              &
                                        ! OUT temperature
     &,RHOKH_mix(row_length,rows,BL_LEVELS)                             
                                        ! OUT Exch coeffs for moisture

!  External references :-
      EXTERNAL TIMER

!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-

#include "c_r_cp.h"

!  Local scalars :-

      REAL :: R_SQ, RBT, AT

      INTEGER                                                           &
     & I,J                                                              &
                  ! LOCAL Loop counter (horizontal field index).
     &,K          ! LOCAL Loop counter (vertical level index).

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BDY_IMPL4 ',3)
      ENDIF

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        DU(I,J,1) = DU(I,J,1) - CQ_CM_U(I,J,1)*TAU_X(I,J,1)
       ENDDO
      ENDDO

      DO J=1,N_ROWS
       DO I=1,ROW_LENGTH
        DV(I,J,1) = DV(I,J,1) - CQ_CM_V(I,J,1)*TAU_Y(I,J,1)
       ENDDO
      ENDDO

      DO K=2,BL_LEVELS
        DO J=1,ROWS
         DO I=1,ROW_LENGTH
          DU(I,J,K) = DU(I,J,K) - CQ_CM_U(I,J,K)*DU(I,J,K-1)
         ENDDO
        ENDDO
        DO J=1,N_ROWS
         DO I=1,ROW_LENGTH
          DV(I,J,K) = DV(I,J,K) - CQ_CM_V(I,J,K)*DV(I,J,K-1)
         ENDDO
        ENDDO
      ENDDO

      IF ( .NOT. L_CORRECT ) THEN 
!  1st stage: predictor
!  Keep a copy of computed taux_1.
! 
        DO J=1,ROWS
          DO I=1,ROW_LENGTH
            TAUX_star(I,J,1) = TAU_X(I,J,1)
          ENDDO
        ENDDO
        DO J=1,N_ROWS
          DO I=1, ROW_LENGTH
            TAUY_star(I,J,1) = TAU_Y(I,J,1)
          ENDDO
        ENDDO
!
! Save increment from 1st stage. Will be needed for the explicit flux 
! at the next stage.
!
        DO K=1,BL_LEVELS
          DO J=1,ROWS
            DO I=1,ROW_LENGTH
              DU_STAR(I,J,K) = DU(I,J,K)
            ENDDO
          ENDDO
          DO J=1,N_ROWS
            DO I=1,ROW_LENGTH
              DV_STAR(I,J,K) = DV(I,J,K)
            ENDDO
          ENDDO
        ENDDO
!---------------------------------------------------------------------
!
! Complete downward sweep of matrix for increments to TL and QW in the 
! boundary layer. It needs to be done here since the scalar fluxes 
! (FQW(,,1),FTL(,,1)) computed previously by the surface scheme are not
! available when the main downward sweep subroutine bdy_impl3 is first 
! invoked (at the 1st stage of the new solver). 
!
!---------------------------------------------------------------------
        DO J=1,ROWS
          DO I=1,ROW_LENGTH
            R_SQ = R_RHO_LEVELS(i,j,1)*R_RHO_LEVELS(i,j,1)
            RHOKH(I,J,2) = R_SQ * RHOKH(I,J,2)
            FQW(I,j,1)   = R_SQ * FQW(I,j,1)
            FTL(I,j,1)   = R_SQ * FTL(I,j,1)
            FQW(I,j,2)   = R_SQ * FQW(I,j,2)
            FTL(I,j,2)   = R_SQ * FTL(I,j,2)
            DQW(I,J,1) = GAMMA2(I,J) * ( -DTRDZ_CHARNEY_GRID(I,J,1) *   &
     &          ( FQW(I,J,2) - FQW(I,J,1) ) + DQW_NT(I,J,1) )
            DTL(I,J,1) = GAMMA2(I,J) * ( -DTRDZ_CHARNEY_GRID(I,J,1) *   &
     &          ( FTL(I,J,2) - FTL(I,J,1) ) + DTL_NT(I,J,1) )
            AT = -DTRDZ_CHARNEY_GRID(I,J,1) *                           &
     &               GAMMA1(I,J)*RHOKH(I,J,2)*RDZ_CHARNEY_GRID(I,J,2)
            RBT = 1.0 / ( 1.0 - AT*( 1.0 + CT_CTQ(I,J,2) ) )
            DQW(I,J,1) = RBT*(DQW(I,J,1) - AT*DQW(I,J,2) )
            DTL(I,J,1) = RBT*(DTL(I,J,1) - AT*DTL(I,J,2) )
            RHOKH(I,J,2) = RHOKH(I,J,2)/R_SQ
            FQW(I,j,1) = FQW(I,j,1)/R_SQ
            FTL(I,j,1) = FTL(I,j,1)/R_SQ
            FQW(I,j,2) = FQW(I,j,2)/R_SQ
            FTL(I,j,2) = FTL(I,j,2)/R_SQ
          ENDDO
        ENDDO

      ELSE ! L_CORRECT == TRUE: 2nd stage of the scheme
!
! Compute 2nd stage correction (total: from tn to tn+1)
! 
        DO K=1,BL_LEVELS
          DO J=1,ROWS
            DO I=1,ROW_LENGTH
              DU(I,J,K) = DU(I,J,K) + DU_STAR(I,J,K)
            END DO
          END DO
          DO J=1,N_ROWS
            DO I=1,ROW_LENGTH
              DV(I,J,K) = DV(I,J,K) + DV_STAR(I,J,K)
            END DO
          END DO
        END DO
!
! Compute total surface stress from tn to tn+1 (diagnostic):
! TAUX_star: flux from tn to t*.
! TAU_X: flux from t* to tn+1.
!
        DO J=1,ROWS
          DO I=1,ROW_LENGTH
            TAU_X(I,J,1) = TAU_X(I,J,1)+TAUX_star(I,J,1)
          END DO
        END DO
        DO J=1, N_ROWS
          DO I=1,ROW_LENGTH
            TAU_Y(I,J,1) = TAU_Y(I,J,1)+TAUY_star(I,J,1)
          END DO
        END DO

      ENDIF
!
! Update TL, QW and their increments
!
      DO J=1,ROWS
        DO I=1,ROW_LENGTH
          TL(I,J,1) = TL(I,J,1) + DTL(I,J,1)
          QW(I,J,1) = QW(I,J,1) + DQW(I,J,1)
        ENDDO
      ENDDO

      DO K=2,BL_LEVELS
        DO J=1,ROWS
          DO I=1,ROW_LENGTH
            DTL(I,J,K) = DTL(I,J,K) - CT_CTQ(I,J,K)*DTL(I,J,K-1)
            TL(I,J,K) = TL(I,J,K) + DTL(I,J,K)
            DQW(I,J,K) = DQW(I,J,K) - CT_CTQ(I,J,K)*DQW(I,J,K-1)
            QW(I,J,K) = QW(I,J,K) + DQW(I,J,K)
          ENDDO
        ENDDO
      ENDDO !bl_levels
!
! Calculate stress and flux diagnostics using the new scheme equations. 
! The fluxes are calculated only when requested.
!
      IF ( L_TAUX ) THEN
        IF ( .NOT. L_CORRECT ) THEN
          DO K=2,BL_LEVELS
            DO J=1,ROWS
              DO I=1,ROW_LENGTH
                TAUX_star(I,J,K) = GAMMA2(I,J)*TAU_X(I,J,K)             &
     &                       + GAMMA1(I,J)*RHOKM_U(I,J,K)*RDZ_U(I,J,K)  &
     &                       * ( DU(I,J,K)-DU(I,J,K-1) )
              ENDDO
            ENDDO
          ENDDO ! bl_levels
        ELSE
          DO K=2,BL_LEVELS
            DO J=1,ROWS
              DO I=1,ROW_LENGTH
                TAU_X(I,J,K)= TAUX_star(I,J,K)+GAMMA2(I,J)*TAU_X(I,J,K) &
     &                      + GAMMA1(I,J)*RHOKM_U(I,J,K)*RDZ_U(I,J,K)   &
     &                      * ( DU(I,J,K)-DU(I,J,K-1) )
              ENDDO
            ENDDO
          ENDDO ! bl_levels
        ENDIF ! L_correct
      ENDIF ! l_taux

      IF ( L_TAUY ) THEN
        IF ( .NOT. L_CORRECT ) THEN
          DO K=2,BL_LEVELS
            DO J=1,N_ROWS
              DO I=1,ROW_LENGTH
                TAUY_star(I,J,K) = GAMMA2(I,J)*TAU_Y(I,J,K)             &
     &                       + GAMMA1(I,J)*RHOKM_V(I,J,K)*RDZ_V(I,J,K)  &
     &                       * (DV(I,J,K)-DV(I,J,K-1))
              ENDDO
            ENDDO
          ENDDO ! bl_levels
        ELSE
          DO K=2,BL_LEVELS
            DO J=1,N_ROWS
              DO I=1,ROW_LENGTH
                TAU_Y(I,J,K)= TAUY_star(I,J,K)+GAMMA2(I,J)*TAU_Y(I,J,K) &
     &                       + GAMMA1(I,J)*RHOKM_V(I,J,K)*RDZ_V(I,J,K)  &
     &                       * (DV(I,J,K)-DV(I,J,K-1))
              ENDDO
            ENDDO
          ENDDO ! bl_levels
        ENDIF ! L_correct
      ENDIF ! l_tauy

      IF ( L_FTL ) THEN
        IF ( .NOT. L_CORRECT ) THEN
          DO K=2,BL_LEVELS
            DO J=1,ROWS
              DO I=1,ROW_LENGTH
                FTL_star(I,J,K) = GAMMA2(I,J)*FTL(I,J,K)                &
     &               - GAMMA1(I,J)*RHOKH(I,J,K)*RDZ_CHARNEY_GRID(I,J,K) &
     &                            * (DTL(I,J,K)-DTL(I,J,K-1))
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO K=2,BL_LEVELS
            DO J=1,ROWS
              DO I=1,ROW_LENGTH
                FTL(I,J,K) = FTL_star(I,J,K)+GAMMA2(I,J)*FTL(I,J,K)     &
     &               - GAMMA1(I,J)*RHOKH(I,J,K)*RDZ_CHARNEY_GRID(I,J,K) &
     &                            * (DTL(I,J,K)-DTL(I,J,K-1))
              ENDDO
            ENDDO
          ENDDO
        ENDIF ! L_correct
      ENDIF ! L_ftl

      IF ( L_FQW ) THEN
        IF ( .NOT. L_correct ) THEN
          DO K=2,BL_LEVELS
            DO J=1,ROWS
              DO I=1,ROW_LENGTH
                FQW_star(I,J,K) = GAMMA2(I,J)*FQW(I,J,K)                &
     &               - GAMMA1(I,J)*RHOKH(I,J,K)*RDZ_CHARNEY_GRID(I,J,K) &
     &               * (DQW(I,J,K)-DQW(I,J,K-1))
              ENDDO
            ENDDO
          ENDDO ! bl_levels
        ELSE
          DO K=2,BL_LEVELS
            DO J=1,ROWS
              DO I=1,ROW_LENGTH
                FQW(I,J,K) = FQW_star(I,J,K)+GAMMA2(I,J)*FQW(I,J,K)     &
     &               - GAMMA1(I,J)*RHOKH(I,J,K)*RDZ_CHARNEY_GRID(I,J,K) &
     &                 * (DQW(I,J,K)-DQW(I,J,K-1))
              ENDDO
            ENDDO
          ENDDO ! bl_levels
        ENDIF
      ENDIF


      IF ( L_CORRECT ) THEN 
!
!-----------------------------------------------------------------------
!!     Convert FTL to sensible heat flux in Watts per square metre.
!      Also, IMPL_CAL only updates FTL_TILE(*,1) and FQW_TILE(*,1)
!      over sea points, so copy this to remaining tiles
!-----------------------------------------------------------------------
        DO K=2,BL_LEVELS
!fpp$ Select(CONCUR)
          DO J=1,ROWS
            DO I=1,ROW_LENGTH
              FTL(I,J,K) = FTL(I,J,K)*CP
            ENDDO
          ENDDO
        ENDDO

!Copy T and Q from workspace to INOUT space.

        DO K=1,BL_LEVELS
!fpp$  Select(CONCUR)
          DO J=1,ROWS
            DO I=1,ROW_LENGTH
              T_LATEST(I,J,K)=TL(I,J,K)
              Q_LATEST(I,J,K)=QW(I,J,K)
            ENDDO
          ENDDO
        ENDDO

!-----------------------------------------------------------------------
!    Set RHOKH, the coefficients required for tracer mixing.
!    Required 5B and after due to change in contents of RHOKH in rest
!    of routine.
!-----------------------------------------------------------------------
        DO K = 2,BL_LEVELS
          DO J=1,ROWS
            DO I=1,ROW_LENGTH
              RHOKH_MIX(I,J,K) = RHOKH(I,J,K)*                          &
     &                       RDZ_CHARNEY_GRID(I,J,K)
            ENDDO
          ENDDO
        ENDDO

      ENDIF ! L_CORRECT 

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BDY_IMPL4 ',4)
      ENDIF

      RETURN
      END SUBROUTINE BDY_IMPL4
#endif
