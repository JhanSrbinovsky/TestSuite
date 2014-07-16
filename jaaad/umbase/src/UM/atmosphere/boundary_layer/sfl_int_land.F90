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

!!!
!!!---------------------------------------------------------------------
!*L  Arguments :-
      SUBROUTINE SFL_INT_LAND (                                         &
     & ROW_LENGTH,ROWS,OFF_X,OFF_Y,LAND_PTS,TILE_PTS,                   &
     & TILE_INDEX,LAND_INDEX,FLANDG                                     &
     &,VSHR,CD_STD,CD,CH,TILE_FRAC                                      &
     &,Z0M,Z0M_STD,Z0H,I_SCRN_T_DIAG                                    &
     &,RECIP_L_MO,V_S,V_S_STD                                           &
     &,Z1_UV,Z1_TQ,DB                                                   &
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
     &,LAND_PTS                                                         &
                            ! IN Number of land points.
     &,TILE_PTS                                                         &
                            ! IN Number of tile points.
     &,TILE_INDEX(LAND_PTS)                                             &
                            ! IN Index of tile points.
     &,LAND_INDEX(LAND_PTS) ! IN Index of land points.
!
      INTEGER, Intent(IN) :: I_SCRN_T_DIAG
!                           ! Method of diagnosing screen
!                           !  temperature

      REAL                                                              &
     & FLANDG(ROW_LENGTH,ROWS)                                          &
!                           ! IN Land fraction
     &,Z0M(LAND_PTS)                                                    &
                            ! IN Roughness length for momentum (m).
     &,Z0H(LAND_PTS)                                                    &
                            ! IN Roughness length for heat and
!                         !    moisture (m).
     &,Z0M_STD(LAND_PTS  )                                              &
                            ! IN Roughness length for momentum without
!                         !    orographic component (m).
     &,VSHR(ROW_LENGTH,ROWS)                                            &
                            ! IN Wind speed difference between the
!                           !    surface and the lowest wind level in
!                           !    the atmosphere (m/s).
     &,CD(LAND_PTS)                                                     &
                          ! IN Surface drag coefficient.
     &,CH(LAND_PTS)                                                     &
                          ! IN Surface transfer coefficient for heat and
!                         !    moisture.
     &,CD_STD(LAND_PTS)                                                 &
                          ! IN Surface drag coefficient excluding
!                         !    orographic from drag.
     &,TILE_FRAC(LAND_PTS)                                              &
!                         ! IN Tile fraction.
     &,RECIP_L_MO(LAND_PTS)                                             &
!                        ! IN Reciprocal of the Monin-Obukhov length (m)
     &,V_S(LAND_PTS)                                                    &
                          ! IN Surface layer scaling velocity including
!                         !    orographic form drag (m/s).
     &,V_S_STD(LAND_PTS)                                                &
!                         ! IN Surface layer scaling velocity excluding
!                         !    orographic form drag (m/s).
     &,Z1_TQ(ROW_LENGTH,ROWS)                                           &
!                         ! IN Height of lowest TQ level (m).
     &,Z1_UV(ROW_LENGTH,ROWS)                                           &
!                         ! IN Height of lowest UV level (m).
     &,DB(LAND_PTS)       ! IN Buoyancy difference between
!                         !    surface and lowest atmospheric
!                         !    level

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
!                        ! OUT interpolation coefficicent for 10m wind
     &,CHR1P5M(LAND_PTS) ! OUT Interpolation coefficient for 1.5m
!                        !     temperature

!*
!*L---------------------------------------------------------------------


!  External routines called :-
      EXTERNAL PHI_M_H_LAND
      EXTERNAL TIMER

!*
!*L---------------------------------------------------------------------
!    Local and other symbolic constants :-
#include "c_vkman.h"
#include "screendiag.h"
      LOGICAL EFF_INT
      PARAMETER (EFF_INT = .FALSE.)
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
!                              ! observations.
     &,PHI_M_OBS(LAND_PTS)                                              &
                            ! Monin-Obukhov stability function for
!                           ! momentum integrated to the wind observatio
!                           ! height.
     &,PHI_H_OBS(LAND_PTS)  ! Monin-Obukhov stability function for
!                           ! scalars integrated to their observation
!                           ! height.
!
!  (b) Scalars.
!
      INTEGER                                                           &
     & I,J,K,L       ! Loop counter (horizontal field index).
      REAL                                                              &
     &   RIB             ! Bulk Richardson number of lowest layer

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
        DO K=1,TILE_PTS
          L = TILE_INDEX(K)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          Z_WIND(I,J) = Z_OBS_WIND
          Z_TEMP(I,J) = Z_OBS_TQ + Z0H(L) - Z0M(L)
        ENDDO
! DEPENDS ON: phi_m_h_land
        CALL PHI_M_H_LAND (ROW_LENGTH,ROWS,LAND_PTS,TILE_PTS,           &
     &                     TILE_INDEX,LAND_INDEX,                       &
     &                     RECIP_L_MO,Z_WIND,Z_TEMP,Z0M,Z0H,            &
     &                     PHI_M_OBS,PHI_H_OBS,LTIMER)
      ENDIF

!-----------------------------------------------------------------------
!! 2. If diagnostics required calculate interpolation coefficient
!!    for 1.5m screen temperature and specific humidity.
!-----------------------------------------------------------------------
!
      IF (ST1P5 .OR. SQ1P5) THEN
!
!       Calculate the screen temperature allowing for decoupling or
!       using pure surface similarity theory as the default. Seperate
!       blocks of code are used for efficiency.
!       
        IF (I_SCRN_T_DIAG == IP_scrn_decpl_1) THEN
!
          DO K=1,TILE_PTS
            L = TILE_INDEX(K)
            J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
            I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
            RIB = ( Z1_UV(I,J) * Z1_UV(I,J) * DB(L) ) /                 & 
     &            ( Z1_TQ(I,J) * VSHR(I,J) * VSHR(I,J) ) 
            IF (RIB> 0.25) THEN 
!             Allow for decoupling in very stable conditions 
!             based on the quasi-equilibrium radiative solution.
!             Note: This value is set for a screen level of 1.5m 
!             and has been fitted for the bottomlevel lying between 
!             1.5 and 20m. It should be recalibrated for coarser 
!             resolutions. 
              CHR1P5M(L) = 0.335+1.78/Z1_TQ(I,J)-1.19/Z1_TQ(I,J)**2 
            ELSE
!             Use pure surface similarity theory
              CHR1P5M(L) = CH(L) * VSHR(I,J) *                          &
     &                      PHI_H_OBS(L)/(VKMAN*V_S_STD(L))
            ENDIF
          ENDDO
!
        ELSE
!
          DO K=1,TILE_PTS
            L = TILE_INDEX(K)
            J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
            I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
            CHR1P5M(L) = CH(L) * VSHR(I,J) *                            &
     &                    PHI_H_OBS(L)/(VKMAN*V_S_STD(L))

          ENDDO
!
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
!! 3. If diagnostics required calculate interpolation coefficient
!!    for 10m winds.
!-----------------------------------------------------------------------
!
! EAK sea contribution to CDR10M already calculated in sfl_int_sea  
      IF ( (SU10 .OR. SV10) .AND. EFF_INT ) THEN
        DO K=1,TILE_PTS
          L = TILE_INDEX(K)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          CDR10M(I,J) = CDR10M(I,J) + FLANDG(I,J) * TILE_FRAC(L) *      &
     &              CD(L) * VSHR(I,J) * PHI_M_OBS(L)/(VKMAN*V_S(L))
        ENDDO
      ELSEIF ( (SU10 .OR. SV10) .AND. .NOT.EFF_INT ) THEN
        DO K=1,TILE_PTS
          L = TILE_INDEX(K)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          Z_TEMP(I,J) = Z_OBS_TQ + Z0H(L) - Z0M_STD(L)
        ENDDO
! DEPENDS ON: phi_m_h_land
        CALL PHI_M_H_LAND (ROW_LENGTH,ROWS,LAND_PTS,TILE_PTS,           &
     &                     TILE_INDEX,LAND_INDEX,                       &
     &                     RECIP_L_MO,Z_WIND,Z_TEMP,Z0M_STD,Z0H,        &
     &                     PHI_M_OBS,PHI_H_OBS,LTIMER)
        DO K=1,TILE_PTS
          L = TILE_INDEX(K)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          CDR10M(I,J) = CDR10M(I,J) + FLANDG(I,J) * TILE_FRAC(L) *      &
     &                CD_STD(L) * VSHR(I,J) * PHI_M_OBS(L)/             &
     &                     (VKMAN*V_S_STD(L))
        ENDDO
      ENDIF
!
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SFL_INT ',4)
      ENDIF
      RETURN
      END SUBROUTINE SFL_INT_LAND
#endif
