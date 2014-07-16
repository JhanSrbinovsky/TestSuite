#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!!!  SUBROUTINES STDEV1_SEA and STDEV1_LAND ----------------------------
!!!
!!!  Purpose: Calculate the standard deviations of layer 1 turbulent
!!!           fluctuations of temperature and humidity using approximate
!!!           formulae from first order closure.
!!!
!!!    Model            Modification history
!!!   version  date
!!!    5.2   15/11/00   New Deck         M. Best
!    5.3  25/04/01  Add coastal tiling. Nic Gedney
!!!
!!!    Programming standard:
!!!
!!!  -------------------------------------------------------------------
!

!!!  SUBROUTINE STDEV1_SEA ---------------------------------------------
!!!  Layer 1 standard deviations for sea and sea-ice
!!!  -------------------------------------------------------------------
      SUBROUTINE STDEV1_SEA (                                           &
     & ROW_LENGTH,ROWS,OFF_X,OFF_Y,FLANDG,                              &
     & BQ_1,BT_1,FQW_1,FTL_1,ICE_FRACT,RHOKM_1,RHOSTAR,VSHR,            &
     & Z0MSEA,Z0_ICE,Z1_TQ,                                             &
     & Q1_SD,T1_SD,LTIMER                                               &
     & )

      IMPLICIT NONE

      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                             ! IN Number of X points?
     &,ROWS                                                             &
                             ! IN Number of Y points?
     &,OFF_X                                                            &
                             ! Size of small halo in i.
     &,OFF_Y                 ! Size of small halo in j.

      LOGICAL                                                           &
     & LTIMER                     ! IN logical for TIMER

      REAL                                                              &
     & FLANDG(ROW_LENGTH,ROWS)                                          &
                                  ! IN Land fraction.
     &,BQ_1(ROW_LENGTH,ROWS)                                            &
                                  ! IN Buoyancy parameter.
     &,BT_1(ROW_LENGTH,ROWS)                                            &
                                  ! IN Buoyancy parameter.
     &,FQW_1(ROW_LENGTH,ROWS)                                           &
                                  ! IN Surface flux of QW.
     &,FTL_1(ROW_LENGTH,ROWS)                                           &
                                  ! IN Surface flux of TL.
     &,ICE_FRACT(ROW_LENGTH,ROWS)                                       &
                                  ! IN Fraction of gridbox which is
!                                 !    sea-ice.
     &,RHOKM_1(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)             &
!                                 ! IN Surface momentum exchange
!                                 !    coefficient.
     &,RHOSTAR(ROW_LENGTH,ROWS)                                         &
                                  ! IN Surface air density.
     &,VSHR(ROW_LENGTH,ROWS)                                            &
                                  ! IN Magnitude of surface-
!                                 !    to-lowest-level wind shear.
     &,Z0MSEA(ROW_LENGTH,ROWS)                                          &
                                  ! IN Sea roughness length.
     &,Z0_ICE(ROW_LENGTH,ROWS)                                          &
                                  ! IN Sea-ice roughness length.
     &,Z1_TQ(ROW_LENGTH,ROWS)     ! IN Height of lowest tq level.

      REAL                                                              &
     & Q1_SD(ROW_LENGTH,ROWS)                                           &
                                  ! OUT Standard deviation of
!                                 !     turbulent fluctuations of
!                                 !     surface layer specific
!                                 !     humidity (kg/kg).
     &,T1_SD(ROW_LENGTH,ROWS)     ! OUT Standard deviation of
!                                 !     turbulent fluctuations of
!                                 !     surface layer temperature (K).


!  External routines called :-
      EXTERNAL TIMER


#include "c_g.h"

!  Workspace --------------------------------------------------------
      INTEGER                                                           &
     & I,J                   ! Loop counter (horizontal field index).
      REAL                                                              &
     & VS                                                               &
                             ! Surface layer friction velocity
     &,VSF1_CUBED                                                       &
                             ! Cube of surface layer free convective
!                            ! scaling velocity
     &,WS1                                                              &
                             ! Turbulent velocity scale for surface
!                            ! layer
     &,Z0                    ! Roughness length

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('STDEV1  ',3)
      ENDIF

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        IF ( FLANDG(I,J) <  1.0 ) THEN

          Z0 = Z0MSEA(I,J)
          IF ( ICE_FRACT(I,J)  >   0. ) Z0 = Z0_ICE(I,J)
          VS = SQRT ( RHOKM_1(I,J)/RHOSTAR(I,J) * VSHR(I,J) )
          VSF1_CUBED = 1.25*G*(Z1_TQ(I,J) + Z0) *                       &
     &                ( BT_1(I,J)*FTL_1(I,J) +                          &
     &                   BQ_1(I,J)*FQW_1(I,J) )/RHOSTAR(I,J)
          IF ( VSF1_CUBED  >   0.0 ) THEN
            WS1 = ( VSF1_CUBED + VS*VS*VS ) ** (1.0/3.0)
            T1_SD(I,J) = MAX ( 0.0 ,                                    &
     &          (1.-FLANDG(I,J))*1.93*FTL_1(I,J) / (RHOSTAR(I,J)*WS1) )
            Q1_SD(I,J) = MAX ( 0.0 ,                                    &
     &          (1.-FLANDG(I,J))*1.93*FQW_1(I,J) / (RHOSTAR(I,J)*WS1) )

          ENDIF

        ENDIF
       ENDDO
      ENDDO

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('STDEV1  ',4)
      ENDIF

      RETURN
      END SUBROUTINE STDEV1_SEA

!!!  SUBROUTINE STDEV1_LAND --------------------------------------------
!!!  Layer 1 standard deviations for land tiles
!!!  -------------------------------------------------------------------
#endif
