
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!!!  SUBROUTINES SF_RIB_LAND and SF_RIB_SEA ---------------------------
!!!
!!!  Purpose: Calculate bulk Richardson number for surface layer
!!!
!!!
!!!    Model            Modification history
!!!   version  date
!!!    5.2   15/11/00   New Deck         M. Best
!    5.3  25/04/01  Add coastal tiling. Nic Gedney
!!!
!!!    Programming standard:
!!!
!!!  ------------------------------------------------------------------

!    SUBROUTINE SF_RIB_LAND--------------------------------------------
!
!    Calculate RIB for land tiles
!
!    ------------------------------------------------------------------

!    SUBROUTINE SF_RIB_SEA---------------------------------------------
!
!    Calculate RIB for sea, sea-ice and sea-ice leads
!
!    ------------------------------------------------------------------
      SUBROUTINE SF_RIB_SEA (                                           &
     & ROW_LENGTH,ROWS,FLANDG,NSICE,SICE_INDEX,                         &
     & BQ_1,BT_1,ICE_FRACT,QSTAR_ICE,QSTAR_SEA,QW_1,TL_1,TSTAR_SICE,    &
     & TSTAR_SEA,VSHR,Z0H_ICE,Z0H_SEA,Z0M_ICE,Z0M_SEA,Z1_TQ,Z1_UV,      &
     & RIB_SEA,RIB_ICE,DB_SEA,DB_ICE,LTIMER                             &
     & )

      IMPLICIT NONE

      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                                     ! IN Number of X points?
     &,ROWS                                                             &
                                     ! IN Number of Y points?
     &,NSICE                                                            &
                                     ! IN Number of sea-ice points.
     &,SICE_INDEX(ROW_LENGTH*ROWS,2) ! IN Index of sea-ice points.

      LOGICAL                                                           &
     & LTIMER                     ! IN logical for TIMER

      REAL                                                              &
     & FLANDG(ROW_LENGTH,ROWS)                                          &
                                  ! IN Land fraction on all pts.
     &,BQ_1(ROW_LENGTH,ROWS)                                            &
                                  ! IN A buoyancy parameter for lowest
!                                 !    atm level. ("beta-q twiddle").
     &,BT_1(ROW_LENGTH,ROWS)                                            &
                                  ! IN A buoyancy parameter for lowest
!                                 !    atm level. ("beta-T twiddle").
     &,ICE_FRACT(ROW_LENGTH,ROWS)                                       &
                                  ! IN Fraction of gridbox which is
!                                 !    sea-ice.
     &,QSTAR_ICE(ROW_LENGTH,ROWS)                                       &
                                  ! IN Surface saturated sp humidity
!                                 !    over sea-ice.
     &,QSTAR_SEA(ROW_LENGTH,ROWS)                                       &
                                  ! IN Surface saturated sp humidity
!                                 !    over sea and sea-ice leads.
     &,QW_1(ROW_LENGTH,ROWS)                                            &
                                  ! IN Total water content of lowest
!                                 !    atmospheric layer (kg per kg air)
     &,TL_1(ROW_LENGTH,ROWS)                                            &
                                  ! IN Liquid/frozen water temperature
!                                 !    for lowest atmospheric layer (K).
     &,TSTAR_SICE(ROW_LENGTH,ROWS)                                      &
                                   ! IN Surface temperature of
!                                 !    sea-ice (K).
     &,TSTAR_SEA(ROW_LENGTH,ROWS)                                       &
                                  ! IN Surface temperature of sea and
!                                 !    sea-ice leads (K).
     &,VSHR(ROW_LENGTH,ROWS)                                            &
                                  ! IN Magnitude of surface-
!                                 !    to-lowest-level wind shear.
     &,Z0H_ICE(ROW_LENGTH,ROWS)                                         &
                                  ! IN Roughness length for heat and
!                                 !    moisture transport over
!                                 !    sea-ice (m).
     &,Z0H_SEA(ROW_LENGTH,ROWS)                                         &
                                  ! IN Roughness length for heat and
!                                 !    moisture transport over sea or
!                                 !    sea-ice leads (m).
     &,Z0M_ICE(ROW_LENGTH,ROWS)                                         &
                                  ! IN Roughness length for momentum
!                                 !    over sea-ice (m).
     &,Z0M_SEA(ROW_LENGTH,ROWS)                                         &
                                  ! IN Roughness length for momentum
!                                 !    over sea or sea-ice leads (m).
     &,Z1_TQ(ROW_LENGTH,ROWS)                                           &
                                  ! IN Height of lowest TQ level (m).
     &,Z1_UV(ROW_LENGTH,ROWS)     ! IN Height of lowest UV level (m).

      REAL                                                              &
     & RIB_SEA(ROW_LENGTH,ROWS)                                         &
                                 ! OUT Bulk Richardson number for lowest
!                                !     layer over sea or sea-ice leads.
     &,RIB_ICE(ROW_LENGTH,ROWS)                                         &
                                 ! OUT Bulk Richardson number for lowest
!                                !     layer over sea-ice.
     &,DB_SEA(ROW_LENGTH,ROWS)                                          &
                                 ! OUT Buoyancy difference between
!                                !     surface and lowest atmospheric
!                                !     level over sea or sea-ice leads.
     &,DB_ICE(ROW_LENGTH,ROWS)   ! OUT Buoyancy difference between
!                                !     surface and lowest atmospheric
!                                !     level over sea-ice.


!  External routines called :-
      EXTERNAL TIMER


!  Symbolic constants -----------------------------------------------

!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
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

!  Workspace --------------------------------------------------------
      INTEGER                                                           &
     & I,J                                                              &
                           ! Horizontal field index.
     &,SI                  !Sea-ice field index.
      REAL                                                              &
     & DQ                                                               &
                           ! Sp humidity difference between surface
!                          ! and lowest atmospheric level (Q1 - Q*).
     &,DTEMP               ! Modified temperature difference between
!                          ! surface and lowest atmospheric level.

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SF_RIB  ',3)
      ENDIF

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        IF ( FLANDG(I,J) <  1.0 ) THEN
! Sea and sea-ice leads
          DTEMP = TL_1(I,J) - TSTAR_SEA(I,J)                            &
                                                           ! P243.118
     &            + (G/CP)*(Z1_TQ(I,J) + Z0M_SEA(I,J) - Z0H_SEA(I,J))
          DQ = QW_1(I,J) - QSTAR_SEA(I,J)                  ! P243.119
          DB_SEA(I,J) = G*( BT_1(I,J)*DTEMP + BQ_1(I,J)*DQ )
          RIB_SEA(I,J) = Z1_UV(I,J)*DB_SEA(I,J) /                       &
     &                     ( VSHR(I,J)*VSHR(I,J) )
        ENDIF
       ENDDO
      ENDDO

      DO SI=1,NSICE
        I = SICE_INDEX(SI,1)
        J = SICE_INDEX(SI,2)
! Sea-ice
        DTEMP = TL_1(I,J) - TSTAR_SICE(I,J)                             &
     &           + (G/CP)*(Z1_TQ(I,J) + Z0M_ICE(I,J) - Z0H_ICE(I,J))
        DQ = QW_1(I,J) - QSTAR_ICE(I,J)
        DB_ICE(I,J) = G*( BT_1(I,J)*DTEMP + BQ_1(I,J)*DQ )
        RIB_ICE(I,J) = Z1_UV(I,J)*DB_ICE(I,J) /                         &
     &                    ( VSHR(I,J) * VSHR(I,J) )
      ENDDO

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SF_RIB  ',4)
      ENDIF

      RETURN
      END SUBROUTINE SF_RIB_SEA
