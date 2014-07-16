
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!**********************************************************************
! Routine to calculate the aerodynamic resistance
!
!
!    Model            Modification history
!   version  date
!    5.2   15/11/00   New Deck         M. Best
!
!    Programming standard:
!
!**********************************************************************
      SUBROUTINE RAERO (ROW_LENGTH,ROWS,LAND_PTS,LAND_INDEX             &
     &,                 VEG_PTS,VEG_INDEX                               &
     &,                 RIB,WIND,Z0H,Z0M,Z1,RA)

      IMPLICIT NONE

      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                                  ! IN Number of points on a row
     &,ROWS                                                             &
                                  ! IN Number of rows in a theta field
     &,LAND_PTS                                                         &
                                  ! IN Total number of land points.
     &,LAND_INDEX(LAND_PTS)                                             &
                                  ! IN Index of land points on the
!                                 !    P-grid.
     &,VEG_PTS                                                          &
                                  ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_PTS)                                              &
                                  ! IN Index of vegetated points
!                                 !    on the land grid.
     &,I,J,K,L                    ! WORK Loop counters.

      REAL                                                              &
     & RIB(ROW_LENGTH,ROWS)                                             &
                                  ! IN Bulk Richardson number.
     &,WIND(ROW_LENGTH,ROWS)                                            &
                                  ! IN Windspeed (m/s).
     &,Z0H(LAND_PTS)                                                    &
                                  ! IN Roughness length for heat (m).
     &,Z0M(LAND_PTS)                                                    &
                                  ! IN Roughness length for momentum (m)
     &,Z1(ROW_LENGTH,ROWS)                                              &
                                  ! IN Reference level (m).
     &,RA(LAND_PTS)                                                     &
                                  ! OUT Aerodynamic resistance (s/m).
     &,BH                                                               &
                                  ! WORK Stability coefficient.
     &,CHN(LAND_PTS)                                                    &
                                  ! WORK Neutral drag coefficient.
     &,FH(LAND_PTS)                                                     &
                                  ! WORK Stability factor.
     &,ZETAH,ZETAM                ! WORK Tempories in calculation of
!                                 !      CHN.
!-----------------------------------------------------------------------
! Local parameters
!-----------------------------------------------------------------------
      REAL                                                              &
     & AH,CZ                                                            &
                                  ! Stability parameters.
     &,KARMAN                     ! Von Karman's constant.
      PARAMETER (AH = 10.0, CZ = 4.0, KARMAN = 0.4)

!-----------------------------------------------------------------------
! Calculate the neutral bulk tranfer coefficient.
!-----------------------------------------------------------------------
      DO K=1,VEG_PTS
        L = VEG_INDEX(K)
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        ZETAM = LOG((Z1(I,J) + Z0M(L)) / Z0M(L))
        ZETAH = LOG((Z1(I,J) + Z0M(L)) / Z0H(L))
        CHN(L) = (KARMAN * KARMAN) / (ZETAH * ZETAM)
      ENDDO

!-----------------------------------------------------------------------
! Calculate the stability factor.
!-----------------------------------------------------------------------
      DO K=1,VEG_PTS
        L = VEG_INDEX(K)
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        BH = AH * CHN(L) * CZ * SQRT (Z1(I,J) / Z0H(L))
        IF (RIB(I,J)  >=  0.0) THEN
          FH(L) = 1.0 / (1 + AH * RIB(I,J))
        ELSE
          FH(L) = 1 - AH * RIB(I,J) / (1 + BH * SQRT(-RIB(I,J)))
        ENDIF
      ENDDO

!-----------------------------------------------------------------------
! Calculate the aerodynamic resistance.
!-----------------------------------------------------------------------
      DO K=1,VEG_PTS
        L = VEG_INDEX(K)
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        RA(L) = 1.0 / (FH(L) * CHN(L) * WIND(I,J))
      ENDDO

      RETURN
      END SUBROUTINE RAERO
