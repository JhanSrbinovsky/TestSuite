#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE BDY_IMPL2----------------------------------------------
!!!
!!!  Purpose: Calculate implicit correction to boundary layer fluxes of
!!!           heat, moisture and momentum.
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!   5.2 15/11/00     New deck. M. Best
!  6.1  17/05/04  Pass switches for diagnostic calculations.
!                                                       M. Diamantakis
!!!  6.2  21/03/05  Pass through implicit scheme weights.
!!!                                       M. Diamantakis
!!!
!!!  6.2  21/03/05  Correct RHOKH_MIX used in TRMIX2C.dk
!!!                                        M. Diamantakis
!!!  Programming standard: Unified Model Documentation Paper No 4,
!!!                        Version ?, dated ?.
!!!
!!!  System component covered: P24.
!!!
!!!  Project task:
!!!
!!!  Documentation: UMDP 24.
!!!
!!!---------------------------------------------------------------------

!    Arguments :-
      SUBROUTINE BDY_IMPL2 (                                            &

! IN values defining field dimensions and subset to be processed :
     & OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,                              &

! IN values defining vertical grid of model atmosphere :
     & BL_LEVELS,                                                       &

! IN Substepping information
     & l_ftl, l_fqw, l_taux, l_tauy,                                    &

! IN data :
     & GAMMA,                                                           &
     & RHOKH,RHOKM_U,RHOKM_V,RDZ_CHARNEY_GRID,RDZ_U,RDZ_V,              &

! INOUT data :
     & QW,TL,FQW,FTL,TAUX,TAUY,                                         &
     & DU,DV,CT_CTQ,DQW,DTL,CQ_CM_U,CQ_CM_V,                            &

! OUT data :
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
     &, OFF_Y      ! Size of small halo in j.

! (b) Defining vertical grid of model atmosphere.

      INTEGER                                                           &
     & BL_LEVELS                   ! IN Max. no. of "boundary" levels
!                                     allowed.Assumed <= 30 for dim-
!                                     sioning of GAMMA in common deck
!                                     C_GAMMA used in SF_EXCH and KMKH

!  In :-

      Logical                                                           &
     & l_ftl                                                            &
     &,l_fqw                                                            &
     &,l_taux                                                           &
     &,l_tauy


      REAL                                                              &
     & RHOKH(ROW_LENGTH,ROWS,BL_LEVELS)                                 &
!                                  ! IN Exchange coeffs for moisture.
     &,RHOKM_U(ROW_LENGTH,ROWS,BL_LEVELS)                               &
!                                  ! IN Exchange coefficients for U
     &,RHOKM_V(ROW_LENGTH,N_ROWS,BL_LEVELS)                             &
!                                  ! IN Exchange coefficients for V
     &,RDZ_CHARNEY_GRID(ROW_LENGTH,ROWS,BL_LEVELS)                      &
!                                  ! IN RDZ(,1) is the reciprocal of the
!                                  ! height of level 1, i.e. of the
!                                  ! middle of layer 1.  For K > 1,
!                                  ! RDZ(,K) is the reciprocal
!                                  ! of the vertical distance
!                                  ! from level K-1 to level K.
     &,RDZ_U(ROW_LENGTH,ROWS,2:BL_LEVELS)                               &
!                                  ! IN  RDZ (K > 1) on U-grid.
     &,RDZ_V(ROW_LENGTH,N_ROWS,2:BL_LEVELS)                             &
!                                  ! IN  RDZ (K > 1) on V-grid.
     &,GAMMA(BL_LEVELS)         ! IN implicit weights

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
     &,TAUX(ROW_LENGTH,ROWS,BL_LEVELS)                                  &
!                                  ! INOUT W'ly component of surface
!                                  !       wind stress (N/sq m).(On
!                                  !       UV-grid with first and last
!                                  !       rows undefined or at present,
!                                  !       set to  missing data
     &,TAUY(ROW_LENGTH,N_ROWS,BL_LEVELS)                                &
!                                  ! INOUT S'ly component of surface
!                                  !       wind stress (N/sq m).  On
!                                  !       UV-grid; comments as per TAUX
     &,CT_CTQ(ROW_LENGTH,ROWS,BL_LEVELS)                                &
!                                  ! INOUT Coefficient in T and q
!                                  !       tri-diagonal implicit matrix
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
     &    1-OFF_Y:N_ROWS+OFF_Y,BL_LEVELS)
!                                  ! INOUT BL increment to v wind field

      REAL                                                              &
     & T_LATEST(ROW_LENGTH,ROWS,BL_LEVELS)                              &
!                                  ! OUT Temperature (K)
     &,Q_LATEST(ROW_LENGTH,ROWS,BL_LEVELS)                              &
!                                  ! OUT Specific Humidity (Kg/Kg)
     &,RHOKH_MIX(ROW_LENGTH,ROWS,BL_LEVELS)
                                   ! OUT Exchange coeffs for moisture.


!  External references :-
      EXTERNAL IM_BL_PT2
      EXTERNAL TIMER

!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-

#include "c_r_cp.h"

!  Local scalars :-

      INTEGER                                                           &
     & I,J                                                              &
                  ! LOCAL Loop counter (horizontal field index).
     &,K          ! LOCAL Loop counter (vertical level index).

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BDY_IMPL2 ',3)
      ENDIF

! DEPENDS ON: im_bl_pt2
      CALL IM_BL_PT2 (                                                  &
     & OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,BL_LEVELS                     &
     &,l_ftl,l_fqw,l_taux,l_tauy                                        &
     &,RHOKH(1,1,2),RHOKM_U(1,1,2),RHOKM_V(1,1,2)                       &
     &,RDZ_CHARNEY_GRID,RDZ_U,RDZ_V,GAMMA                               &
     &,CT_CTQ,DQW,DTL,CQ_CM_U,CQ_CM_V,DU,DV                             &
     &,FQW,FTL,TAUX,TAUY,QW,TL                                          &
     &,LTIMER                                                           &
     &)


!-----------------------------------------------------------------------
!! 6.1 Convert FTL to sensible heat flux in Watts per square metre.
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

!7.1 Copy T and Q from workspace to INOUT space.

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
! 10 Set RHOKH, the coefficients required for tracer mixing.
!    Required 5B and after due to change in contents of RHOKH in rest
!    of routine.
!-----------------------------------------------------------------------

      DO K = 2,BL_LEVELS
       DO J=1,ROWS
        DO I=1,ROW_LENGTH
          RHOKH_MIX(I,J,K) = RHOKH(I,J,K)*                              &
     &                       RDZ_CHARNEY_GRID(I,J,K)
        ENDDO
       ENDDO
      ENDDO

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BDY_IMPL2 ',4)
      ENDIF

      RETURN
      END SUBROUTINE BDY_IMPL2
#endif
