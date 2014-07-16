
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!**********************************************************************
! Subroutine to adjust canopy conductance and soil moisture extraction
! for soil evaporation beneath vegetation.
!
!    Model            Modification history
!   version  date
!    5.2   15/11/00   New Deck         M. Best
!
!    Programming standard:
!
!**********************************************************************

      SUBROUTINE SOIL_EVAP (NPNTS,NSHYD,TILE_PTS,TILE_INDEX,            &
     &                      GSOIL,LAI,GS,WT_EXT)

      IMPLICIT NONE

      INTEGER                                                           &
     & NPNTS                                                            &
                            ! IN Number of gridpoints.
     &,NSHYD                                                            &
                            ! IN Number of soil moisture layers.
     &,TILE_PTS                                                         &
                            ! IN Number of points containing the
!                           !    given surface type.
     &,TILE_INDEX(NPNTS)    ! IN Indices on the land grid of the
!                           !    points containing the given
!                           !    surface type.

       REAL                                                             &
     & GSOIL(NPNTS)                                                     &
                            ! IN Soil surface conductance (m/s).
     &,LAI(NPNTS)           ! IN Leaf area index.

      REAL                                                              &
     & GS(NPNTS)                                                        &
                            ! INOUT Surface conductance (m/s).
     &,WT_EXT(NPNTS,NSHYD)  ! INOUT Fraction of evapotranspiration
!                           !       extracted from each soil layer.

       REAL                                                             &
     & FSOIL(NPNTS)         ! Fraction of ground below canopy
!                           ! contributing to evaporation.

       INTEGER                                                          &
     & J,K,L                ! Loop indices

!CDIR NODEP - order independent so can vectorise
        DO J=1,TILE_PTS
          L=TILE_INDEX(J)
          FSOIL(L) = EXP(-0.5*LAI(L))
        ENDDO

        DO K=2,NSHYD
          DO J=1,TILE_PTS
            L=TILE_INDEX(J)
            WT_EXT(L,K) = GS(L)*WT_EXT(L,K)/(GS(L) + FSOIL(L)*GSOIL(L))
          ENDDO
        ENDDO

!CDIR NODEP
        DO J=1,TILE_PTS
          L=TILE_INDEX(J)
          WT_EXT(L,1) = (GS(L)*WT_EXT(L,1) + FSOIL(L)*GSOIL(L))         &
     &                   / (GS(L) + FSOIL(L)*GSOIL(L))
          GS(L) = GS(L) + FSOIL(L)*GSOIL(L)
        ENDDO

      RETURN
      END SUBROUTINE SOIL_EVAP
