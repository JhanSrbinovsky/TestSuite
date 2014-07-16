
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!   SUBROUTINE SCREEN_TQ----------------------------------------------
!!!
!!!  Purpose: Diagnose temperature and/or specific humidity at screen
!!!           height (1.5 metres), as requested via the STASH flags.
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!   5.2  15/11/00   New Deck         M. Best
!    5.3  25/04/01  Add coastal tiling. Nic Gedney
! 6.2      21/02/06    Changes for mixing ratios   A.P.Lock
!!!
!!!---------------------------------------------------------------------
      SUBROUTINE SCREEN_TQ (                                            &
     & ROW_LENGTH,ROWS,LAND_PTS,NTILES,                                 &
     & LAND_INDEX,TILE_INDEX,TILE_PTS,FLANDG,                           &
     & SQ1P5,ST1P5,CHR1P5M,CHR1P5M_SICE,PSTAR,QW_1,RESFT,               &
     & TILE_FRAC,TL_1,TSTAR_SSI,TSTAR_TILE,                             &
     & Z0HSSI,Z0H_TILE,Z0MSSI,Z0M_TILE,Z1,                              &
     & Q1P5M,Q1P5M_TILE,T1P5M,T1P5M_TILE,                               &
     & lq_mix_bl, l_cable                                               &
     & )

      IMPLICIT NONE

      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                            ! IN Number of X points?
     &,ROWS                                                             &
                            ! IN Number of Y points?
     &,LAND_PTS                                                         &
                            ! IN Number of land points to be processed.
     &,NTILES                                                           &
                            ! IN Number of tiles per land point.
     &,LAND_INDEX(LAND_PTS)                                             &
                            ! IN Index of land points.
     &,TILE_INDEX(LAND_PTS,NTILES)                                      &
!                           ! IN Index of tile points.
     &,TILE_PTS(NTILES)     ! IN Number of tile points.

      LOGICAL                                                           &
     & SQ1P5                                                            &
                            ! IN STASH flag for 1.5-metre sp humidity.
     &,ST1P5                ! IN STASH flag for 1.5-metre temperature.

      LOGICAL                                                           &
     & lq_mix_bl              ! TRUE if mixing ratios used in
!                             ! boundary layer code

      LOGICAL :: l_cable
      REAL                                                              &
     & FLANDG(ROW_LENGTH,ROWS)                                          &
!                           ! IN Fraction of gridbox which is land.
     &,CHR1P5M(LAND_PTS,NTILES)                                         &
!                           ! IN Ratio of coefficients for  calculation
!                           !    of 1.5 m T.
     &,CHR1P5M_SICE(ROW_LENGTH,ROWS)                                    &
!                           ! IN Ratio of coefficients for  calculation
!                           !    of 1.5 m T.
     &,PSTAR(ROW_LENGTH,ROWS)                                           &
                             ! IN Surface pressure (Pa).
     &,QW_1(ROW_LENGTH,ROWS)                                            &
                             ! IN Total water content of lowest
!                                 atmospheric layer (kg per kg air).
     &,RESFT(LAND_PTS,NTILES)                                           &
!                           ! IN Surface resistance factor.
     &,TILE_FRAC(LAND_PTS,NTILES)                                       &
!                           ! IN Tile fractions.
     &,TL_1(ROW_LENGTH,ROWS)                                            &
                            ! IN Liquid/frozen water temperature for
!                                lowest atmospheric layer (K).
     &,TSTAR_SSI(ROW_LENGTH,ROWS)                                       &
!                           ! IN Sea/sea-ice mean sfc temperature (K).
     &,TSTAR_TILE(LAND_PTS,NTILES)                                      &
!                           ! IN Tile surface temperatures (K).
     &,Z0HSSI(ROW_LENGTH,ROWS)                                          &
                            ! IN Roughness length for heat and
!                           !    moisture (m).
     &,Z0H_TILE(LAND_PTS,NTILES)                                        &
!                           ! IN Tile roughness lengths for heat and
!                           !    moisture (m).
     &,Z0MSSI(ROW_LENGTH,ROWS)                                          &
                            ! IN Roughness length for momentum (m).
     &,Z0M_TILE(LAND_PTS,NTILES)                                        &
!                           ! IN Tile roughness lengths for momentum (m)
     &,Z1(ROW_LENGTH,ROWS)  ! IN Height of lowest atmospheric level (m).

      REAL                                                              &
     & Q1P5M(ROW_LENGTH,ROWS)                                           &
                             ! OUT Specific humidity at screen height
!                           !      of 1.5 metres (kg water per kg air).
     &,Q1P5M_TILE(LAND_PTS,NTILES)                                      &
!                           ! OUT Q1P5M over land tiles.
     &,T1P5M(ROW_LENGTH,ROWS)                                           &
                             ! OUT Temperature at screen height of
!                           !     1.5 metres (K).
     &,T1P5M_TILE(LAND_PTS,NTILES)
!                           ! OUT T1P5M over land tiles.


!  External routines called :-
      EXTERNAL QSAT_mix


      REAL                                                              &
     & CER1P5M                                                          &
                            ! Ratio of coefficients reqd for
!                           ! calculation of 1.5 m Q.
     &,PSTAR_LAND(LAND_PTS)                                             &
                            ! Surface pressure for land points.
     &,QS(ROW_LENGTH,ROWS)                                              &
                            ! Surface saturated sp humidity.
     &,QS_TILE(LAND_PTS)    ! Surface saturated sp humidity.

       INTEGER                                                          &
     & I,J                                                              &
                   ! Loop counter (horizontal field index).
     &,K                                                                &
                   ! Loop counter (tile point index).
     &,L                                                                &
                   ! Loop counter (land point field index).
     &,N           ! Loop counter (tile index).

! Local and other symbolic constants used :-

!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
! C_HT_M constants for subroutine SF_EXCH

      ! height of 10m level for diagnostic calculations (m).
      REAL,PARAMETER:: Z10M  = 10.0

      ! height of 1.5m level for diagnostic calculations (m).
!      REAL,PARAMETER:: Z1P5M = 1.5
      REAL,PARAMETER:: Z1P5M = 2.0

! C_HT_M end
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
      REAL GRCP
      PARAMETER ( GRCP = G / CP )

!-----------------------------------------------------------------------
! Diagnose local and GBM temperatures at 1.5 m if requested via ST1P5
!----------------------------------------------------------------------
      IF (ST1P5) THEN

        DO J=1,ROWS
         DO I=1,ROW_LENGTH
          T1P5M(I,J) = 0.
          IF (FLANDG(I,J) <  1.0 ) THEN
            T1P5M(I,J) = (1.-FLANDG(I,J))*                              &
     &        (TSTAR_SSI(I,J) - GRCP*Z1P5M +                            &
     &        CHR1P5M_SICE(I,J) *  (TL_1(I,J) - TSTAR_SSI(I,J) +        &
     &          GRCP*(Z1(I,J)+Z0MSSI(I,J)-Z0HSSI(I,J))))
          ENDIF
         ENDDO
        ENDDO

        DO N=1,NTILES
          if ( .not. l_cable ) then
             DO L=1,LAND_PTS
                T1P5M_TILE(L,N) = 0.
             ENDDO
          end if
          DO K=1,TILE_PTS(N)
            L = TILE_INDEX(K,N)
            J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
            I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
            if ( .not. l_cable ) then
            T1P5M_TILE(L,N) = TSTAR_TILE(L,N) - GRCP*Z1P5M +            &
     &                    CHR1P5M(L,N)*( TL_1(I,J) - TSTAR_TILE(L,N) +  &
     &                    GRCP*(Z1(I,J)+Z0M_TILE(L,N)-Z0H_TILE(L,N)) )
            end if
            T1P5M(I,J) = T1P5M(I,J)                                     &
     &        + FLANDG(I,J)*TILE_FRAC(L,N)*T1P5M_TILE(L,N)
          ENDDO
        ENDDO

      ENDIF

!-----------------------------------------------------------------------
! Diagnose local and GBM humidities at 1.5 m if requested via SQ1P5
!-----------------------------------------------------------------------
      IF (SQ1P5) THEN

! DEPENDS ON: qsat_mix
        CALL QSAT_mix(QS,TSTAR_SSI,PSTAR,ROW_LENGTH*ROWS,lq_mix_bl)
        DO J=1,ROWS
         DO I=1,ROW_LENGTH
          Q1P5M(I,J) = 0.
          IF (FLANDG(I,J) <  1.0 ) THEN
            CER1P5M = CHR1P5M_SICE(I,J) - 1.
            Q1P5M(I,J) = (1.-FLANDG(I,J))*                              &
     &        (QW_1(I,J) + CER1P5M*( QW_1(I,J) - QS(I,J) ))
          ENDIF
         ENDDO
        ENDDO

        DO L=1,LAND_PTS
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          PSTAR_LAND(L) = PSTAR(I,J)
        ENDDO

        DO N=1,NTILES
          if ( .not. l_cable ) then
             DO L=1,LAND_PTS
                Q1P5M_TILE(L,N) = 0.
             ENDDO
          end if
! DEPENDS ON: qsat_mix
          CALL QSAT_mix(QS_TILE,TSTAR_TILE(1,N),PSTAR_LAND,LAND_PTS     &
     &    ,lq_mix_bl)
          DO K=1,TILE_PTS(N)
            L = TILE_INDEX(K,N)
            J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
            I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
            CER1P5M = RESFT(L,N)*(CHR1P5M(L,N) - 1.)
            if ( .not. l_cable ) then
               Q1P5M_TILE(L,N) = QW_1(I,J) +                               &
     &                        CER1P5M*( QW_1(I,J) - QS_TILE(L) )
            end if
            Q1P5M(I,J) = Q1P5M(I,J)                                     &
     &        + FLANDG(I,J)*TILE_FRAC(L,N)*Q1P5M_TILE(L,N)
          ENDDO
        ENDDO

      ENDIF

      RETURN
      END SUBROUTINE SCREEN_TQ
