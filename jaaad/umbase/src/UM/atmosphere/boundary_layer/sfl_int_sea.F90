#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINES SFL_INT_SEA AND SFL_INT_LAND--------------------------
!!!
!!!  Purpose: To calculate interpolation coefficients for 10m winds
!!!           and 1.5m temperature/specific humidity diagnostics.
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!   5.2  15/11/00   New Deck         M. Best
!    5.3  25/04/01  Add coastal tiling. Nic Gedney
!!!
!!!  Programming standard:
!!!
!!!  Logical component covered: Part of P243.
!!!
!!!  System Task:
!!!
!!!  External Documentation: UMDP No.24
!!!
!!!---------------------------------------------------------------------
!*L  Arguments :-
      SUBROUTINE SFL_INT_SEA (                                          &
     & ROW_LENGTH,ROWS,OFF_X,OFF_Y,FLANDG                               &
     &,VSHR,CD,CH,Z0M,Z0H,I_SCRN_T_DIAG                                 &
     &,RECIP_L_MO,V_S                                                   &
     &,ICE_FRACT,NSICE,SICE_INDEX,Z1_UV,Z1_TQ,DB_ICE                    &
     &,SU10,SV10,ST1P5,SQ1P5                                            &
     &,CDR10M,CHR1P5M,LTIMER                                            &
     &)
      IMPLICIT NONE

      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                            ! IN Number of X points?
     &,ROWS                                                             &
                            ! IN Number of Y points?
     &,OFF_X                                                            &
                            ! Size of small halo in i.
     &,OFF_Y                                                            &
                            ! Size of small halo in j.
     &,SICE_INDEX(ROW_LENGTH*ROWS,2)                                    &
                            ! IN Index of sea-ice points.
     &,NSICE
                            ! Number of sea-ice points.
!
      INTEGER, Intent(IN) :: I_SCRN_T_DIAG
!                           ! Method of diagnosing screen
!                           !  temperature



      REAL                                                              &
     & FLANDG(ROW_LENGTH,ROWS)                                          &
!                           ! IN Land fraction.
     &,Z0M(ROW_LENGTH,ROWS)                                             &
                            ! IN Roughness length for momentum (m).
     &,Z0H(ROW_LENGTH,ROWS)                                             &
                            ! IN Roughness length for heat and
!                           !    moisture (m).
     &,VSHR(ROW_LENGTH,ROWS)                                            &
                            ! IN Wind speed difference between the
!                           !    surface and the lowest wind level in
!                           !    the atmosphere (m/s).
     &,CD(ROW_LENGTH,ROWS)                                              &
                            ! IN Surface drag coefficient.
     &,CH(ROW_LENGTH,ROWS)                                              &
                            ! IN Surface transfer coefficient for heat
!                           !    and moisture.
     &,RECIP_L_MO(ROW_LENGTH,ROWS)                                      &
!                           ! IN Reciprocal of the Monin-Obukhov
!                           ! length (m)
     &,V_S(ROW_LENGTH,ROWS)                                             &
!                           ! IN Surface layer scaling velocity
!                           !    including orographic form drag (m/s).
     &,ICE_FRACT(ROW_LENGTH,ROWS)                                       &
!                           ! IN Fraction of gridbox which is
!                           !    sea-ice.
     &,Z1_TQ(ROW_LENGTH,ROWS)                                           &
!                           ! IN Height of lowest TQ level (m).
     &,Z1_UV(ROW_LENGTH,ROWS)                                           &
!                           ! IN Height of lowest UV level (m).
     &,DB_ICE(ROW_LENGTH,ROWS)
!                           ! IN Buoyancy difference between
!                           !    surface and lowest atmospheric
!                           !    level over sea-ice.

      LOGICAL                                                           &
     & SU10                                                             &
                                 ! IN 10m U-wind diagnostic flag
     &,SV10                                                             &
                                 ! IN 10m V-wind diagnostic flag
     &,ST1P5                                                            &
                                 ! IN screen temp diagnostic flag
     &,SQ1P5                                                            &
                                 ! IN screen specific humidity
!                                !    diagnostic flag
     &,LTIMER                    ! IN TIMER diagnostics flag
! Output variables
!
      REAL                                                              &
     & CDR10M(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y)              &
!                                ! OUT interpolation coefficicent for
!                                ! 10m wind
     &,CHR1P5M(ROW_LENGTH,ROWS)  ! OUT Interpolation coefficient for
!                                !     1.5m temperature

!*
!*L---------------------------------------------------------------------


!  External routines called :-
      EXTERNAL PHI_M_H_SEA
      EXTERNAL TIMER

!*
!*L---------------------------------------------------------------------
!    Local and other symbolic constants :-
#include "c_vkman.h"
#include "screendiag.h"
!
!  Define local storage.
!
!  (a) Local work arrays.
!
      REAL                                                              &
     & Z_WIND(ROW_LENGTH,ROWS)                                          &
                                 ! Height of wind observations.
     &,Z_TEMP(ROW_LENGTH,ROWS)                                          &
                                 ! Height of temperature and humidity
!                                ! observations.
     &,PHI_M_OBS(ROW_LENGTH,ROWS)                                       &
                                 ! Monin-Obukhov stability function for
!                                ! momentum integrated to the wind
!                                ! observation height.
     &,PHI_H_OBS(ROW_LENGTH,ROWS)! Monin-Obukhov stability function for
!                                ! scalars integrated to their
!                                ! observation height.
!
!  (b) Scalars.
!
      INTEGER                                                           &
     & I,J     ! Loop counter (horizontal field index).
      INTEGER                                                           &
     &  SI
      REAL                                                              &
     &  RIBI                     ! Bulk Richardson number of bottom
!                                ! layer, used only over ice.
!*
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SFL_INT   ',3)
      ENDIF
!
!-----------------------------------------------------------------------
!! 1. If diagnostics required calculate M-O stability functions at
!!    observation heights.
!-----------------------------------------------------------------------

      IF (SU10 .OR. SV10 .OR. ST1P5 .OR. SQ1P5) THEN
        DO J=1,ROWS
         DO I=1,ROW_LENGTH
          IF ( FLANDG(I,J) <  1.0 ) THEN
          Z_WIND(I,J) = Z_OBS_WIND
          Z_TEMP(I,J) = Z_OBS_TQ + Z0H(I,J) - Z0M(I,J)
          END IF
         ENDDO
        ENDDO
! DEPENDS ON: phi_m_h_sea
        CALL PHI_M_H_SEA (ROW_LENGTH,ROWS,FLANDG,                       &
     &                    RECIP_L_MO,Z_WIND,Z_TEMP,Z0M,Z0H,             &
     &                    PHI_M_OBS,PHI_H_OBS,LTIMER)
      ENDIF

!-----------------------------------------------------------------------
!! 2. If diagnostics required calculate interpolation coefficient
!!    for 1.5m screen temperature and specific humidity.
!-----------------------------------------------------------------------
!
      IF (ST1P5 .OR. SQ1P5) THEN
        DO J=1,ROWS
         DO I=1,ROW_LENGTH
          IF ( FLANDG(I,J) <  1.0 ) THEN
            CHR1P5M(I,J) = CH(I,J) * VSHR(I,J) *                        &
     &                          PHI_H_OBS(I,J)/(VKMAN*V_S(I,J))
          ENDIF
         ENDDO
        ENDDO
!
!       Modify the interpolation coefficient if decoupling is allowed.
        IF (I_SCRN_T_DIAG == IP_scrn_decpl_1) THEN
! 
!         Recalculate at decoupled points. Over grid-boxes with no 
!         sea-ice we do not expect decoupling, so the diagnostic 
!         defaults to the standard form. Once the grid-box contains 
!         any ice we follow the convention of the model that the 
!         diagnostic refers to the iced portion of the grid-box. 
! 
          DO SI=1,NSICE 
! 
            I = SICE_INDEX(SI,1) 
            J = SICE_INDEX(SI,2) 
! 
!           Test the bulk Richardson number of the lowest layer 
!           to see if decoupling is expected. This test uses DB_ICE 
!           from SF_RIB_SEA, but the Richardson number is 
!           recalculated to take proper account of the staggering 
!           of the bottom grid-levels. 
            IF ( ( FLANDG(I,J) < 1.0 ) .AND.                            & 
     &           ( ICE_FRACT(I,J) > 0.0 ) ) THEN 
! 
              RIBI = Z1_UV(I,J) * Z1_UV(I,J) * DB_ICE(I,J) /            & 
     &               ( VSHR(I,J) * VSHR(I,J) * Z1_TQ(I,J) ) 
              IF (RIBI > 0.25) THEN 
!               Note: This value is set for a screen level of 1.5m 
!               and has been fitted for the bottom level lying between 
!               1.5 and 20m. It should be recalibrated for coarser 
!               resolutions. 
                CHR1P5M(I,J) = 0.335+1.78/Z1_TQ(I,J)-1.19/Z1_TQ(I,J)**2 
              ENDIF 
            ENDIF 
! 
          ENDDO 
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
!! 3. If diagnostics required calculate interpolation coefficient
!!    for 10m winds.
!-----------------------------------------------------------------------
!
      IF ( SU10 .OR. SV10 ) THEN
        DO J=1,ROWS
         DO I=1,ROW_LENGTH
          IF ( FLANDG(I,J) <  1.0 ) THEN
            CDR10M(I,J) = (1.0-FLANDG(I,J)) * CD(I,J) * VSHR(I,J) *     &
     &                       PHI_M_OBS(I,J)/(VKMAN*V_S(I,J))
          ENDIF
         ENDDO
        ENDDO
      ENDIF
!
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SFL_INT ',4)
      ENDIF
      RETURN
      END SUBROUTINE SFL_INT_SEA

!!!
!!!---------------------------------------------------------------------
!*L  Arguments :-
#endif
