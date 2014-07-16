#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!   SUBROUTINES PHI_M_H_SEA AND PHI_M_H_LAND ------------------------
!!!
!!!  Purpose: Calculate the integrated froms of the Monin-Obukhov
!!!           stability functions for surface exchanges.
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!   5.2  15/11/00   New Deck         M. Best
!    5.3  25/04/01  Add coastal tiling. Nic Gedney
!  6.0  19/08/03  NEC SX-6 optimisation - force vectorisation of loop
!                 in phi_m_h_land.  R Barnes & J-C Rioual.
!!!
!!!  Programming standard:
!!!
!!!  System component covered: Part of P243.
!!!
!!!  Documentation: UM Documentation Paper No 24.
!!!
!*L  Arguments:---------------------------------------------------------

!!!
!*L  Arguments:---------------------------------------------------------
      SUBROUTINE PHI_M_H_LAND(                                          &
     & ROW_LENGTH,ROWS,LAND_PTS,TILE_PTS,TILE_INDEX,LAND_INDEX,         &
     & RECIP_L_MO,Z_UV,Z_TQ,Z0M,Z0H,PHI_M,PHI_H,LTIMER                  &
     &)
      IMPLICIT NONE

      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                            ! IN Number of X points?
     &,ROWS                                                             &
                            ! IN Number of Y points?
     &,LAND_PTS                                                         &
                            ! IN Number of land points.
     &,TILE_PTS                                                         &
                            ! IN Number of tile points.
     &,TILE_INDEX(LAND_PTS)                                             &
                            ! IN Index of tile points.
     &,LAND_INDEX(LAND_PTS) ! IN Index of land points.

      LOGICAL                                                           &
     & LTIMER


      REAL                                                              &
     & RECIP_L_MO(LAND_PTS)                                             &
!                       ! IN Reciprocal of the Monin-Obukhov length
!                       !     (m^-1).
     &,Z_UV(ROW_LENGTH,ROWS)                                            &
!                       ! IN Height of wind level above roughness
!                       !    height (m)
     &,Z_TQ(ROW_LENGTH,ROWS)                                            &
!                       ! IN Height of temperature, moisture and scalar
!                       !    lev above the roughness height (m).
     &,Z0M(LAND_PTS)                                                    &
                        ! IN Roughness length for momentum (m).
     &,Z0H(LAND_PTS)    ! IN Roughness length for heat/moisture/scalars
!                       !    (m)
!
      REAL                                                              &
     & PHI_M(LAND_PTS)                                                  &
                        ! OUT Stability function for momentum.
     &,PHI_H(LAND_PTS)  ! OUT Stability function for
!                       !     heat/moisture/scalars.
!
!*L  Workspace usage----------------------------------------------------
!    No work areas are required.
!
!*----------------------------------------------------------------------
!*L  External subprograms called:

      EXTERNAL TIMER

!*----------------------------------------------------------------------
!  Common and local physical constants.
!
      REAL A,B,C,D,C_OVER_D
      PARAMETER (                                                       &
     & A=1.0                                                            &
                      !
     &,B=2.0/3.0                                                        &
                      ! Constants used in the Beljaars and
     &,C=5.0                                                            &
                      ! Holtslag stable stability functions
     &,D=0.35                                                           &
                      !
     &,C_OVER_D=C/D                                                     &
                      !
     &)
!
!  Define local variables.
!
      INTEGER I,J,K,L   ! Loop counter; horizontal field index.
!
      REAL                                                              &
     & PHI_MN                                                           &
                      ! Neutral value of stability function for momentum
     &,PHI_HN                                                           &
                      ! Neutral value of stability function for scalars.
     &,ZETA_UV                                                          &
                      ! Temporary in calculation of PHI_M.
     &,ZETA_0M                                                          &
                      ! Temporary in calculation of PHI_M.
     &,ZETA_TQ                                                          &
                      ! Temporary in calculation of PHI_H.
     &,ZETA_0H                                                          &
                      ! Temporary in calculation of PHI_H.
     &,X_UV_SQ                                                          &
                      ! Temporary in calculation of PHI_M.
     &,X_0M_SQ                                                          &
                      ! Temporary in calculation of PHI_M.
     &,X_UV                                                             &
                      ! Temporary in calculation of PHI_M.
     &,X_0M                                                             &
                      ! Temporary in calculation of PHI_M.
     &,Y_TQ                                                             &
                      ! Temporary in calculation of PHI_H.
     &,Y_0H                                                             &
                      ! Temporary in calculation of PHI_H.
     &,PHI_H_FZ1                                                        &
                      ! Temporary in calculation of PHI_H.
     &,PHI_H_FZ0      ! Temporary in calculation of PHI_H.

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('PHI_M_H ',3)
      ENDIF

!CDIR NODEP
      DO K=1,TILE_PTS
        L = TILE_INDEX(K)
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
!
!-----------------------------------------------------------------------
!! 1. Calculate neutral values of PHI_M and PHI_H.
!-----------------------------------------------------------------------
!
        PHI_MN = LOG( (Z_UV(I,J) + Z0M(L)) / Z0M(L) )
        PHI_HN = LOG( (Z_TQ(I,J) + Z0M(L)) / Z0H(L) )
!
!-----------------------------------------------------------------------
!! 2. Calculate stability parameters.
!-----------------------------------------------------------------------
!
        ZETA_UV = (Z_UV(I,J) + Z0M(L)) * RECIP_L_MO(L)
        ZETA_TQ = (Z_TQ(I,J) + Z0M(L)) * RECIP_L_MO(L)
        ZETA_0M = Z0M(L) * RECIP_L_MO(L)
        ZETA_0H = Z0H(L) * RECIP_L_MO(L)
!
!-----------------------------------------------------------------------
!! 3. Calculate PHI_M and PHI_H for neutral and stable conditions.
!!    Formulation of Beljaars and Holtslag (1991).
!-----------------------------------------------------------------------
!
        IF (RECIP_L_MO(L)  >=  0.0) THEN
          PHI_M(L) = PHI_MN                                             &
     &               + A * (ZETA_UV - ZETA_0M)                          &
     &               + B * ( (ZETA_UV - C_OVER_D) * EXP(-D*ZETA_UV)     &
     &                      -(ZETA_0M - C_OVER_D) * EXP(-D*ZETA_0M) )
          PHI_H_FZ1 = SQRT(1.0 + (2.0/3.0)*A*ZETA_TQ)
          PHI_H_FZ0 = SQRT(1.0 + (2.0/3.0)*A*ZETA_0H)
          PHI_H(L) = PHI_HN +                                           &
     &                 PHI_H_FZ1*PHI_H_FZ1*PHI_H_FZ1                    &
     &               - PHI_H_FZ0*PHI_H_FZ0*PHI_H_FZ0                    &
     &               + B * ( (ZETA_TQ - C_OVER_D) * EXP(-D*ZETA_TQ)     &
     &                      -(ZETA_0H - C_OVER_D) * EXP(-D*ZETA_0H) )
!
!-----------------------------------------------------------------------
!! 4. Calculate PHI_M and PHI_H for unstable conditions.
!-----------------------------------------------------------------------
!
        ELSE

          X_UV_SQ = SQRT(1.0 - 16.0*ZETA_UV)
          X_0M_SQ = SQRT(1.0 - 16.0*ZETA_0M)
          X_UV = SQRT(X_UV_SQ)
          X_0M = SQRT(X_0M_SQ)
          PHI_M(L) = PHI_MN - 2.0*LOG( (1.0+X_UV) / (1.0+X_0M) )        &
     &                    - LOG( (1.0+X_UV_SQ) / (1.0+X_0M_SQ) )        &
     &                    + 2.0*( ATAN(X_UV) - ATAN(X_0M) )

          Y_TQ = SQRT(1.0 - 16.0*ZETA_TQ)
          Y_0H = SQRT(1.0 - 16.0*ZETA_0H)
          PHI_H(L) = PHI_HN - 2.0*LOG( (1.0+Y_TQ) / (1.0+Y_0H) )

        ENDIF

      ENDDO

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('PHI_M_H ',4)
      ENDIF

      RETURN
      END SUBROUTINE PHI_M_H_LAND
#endif
