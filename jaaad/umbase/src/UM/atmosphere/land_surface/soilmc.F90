#if defined(A08_7A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!**********************************************************************
!
! Description:
!     Diagnoses the soil moisture in a layer at the surface
!
!    Model            Modification history
!   version  date
!    5.2   15/11/00   New Deck         M. Best
!
!    Programming standard:
!
!**********************************************************************

!    SUBROUTINE SOILMC-------------------------------------------------

      SUBROUTINE SOILMC ( NPNTS,NSHYD,SOIL_PTS,SOIL_INDEX,              &
     &                    DZ,STHU,V_SAT,V_WILT,SMC )

      IMPLICIT NONE

      INTEGER                                                           &
     & NPNTS                                                            &
                            ! IN Number of gridpoints.
     &,NSHYD                                                            &
                            ! IN Number of soil moisture levels.
     &,SOIL_PTS                                                         &
                            ! IN Number of soil points.
     &,SOIL_INDEX(NPNTS)    ! IN Array of soil points.

      REAL                                                              &
     & DZ(NSHYD)                                                        &
                            ! IN Thicknesses of the soil layers (m).
     &,STHU(NPNTS,NSHYD)                                                &
                            ! IN Unfrozen soil moisture content of
!                           !    each layer as a frac. of saturation.
     &,V_SAT(NPNTS)                                                     &
                            ! IN Volumetric soil moisture conc. at
!                           !    saturation (m3 H2O/m3 soil).
     &,V_WILT(NPNTS)        ! IN Volumetric soil moisture conc. below
!                           !    which stomata close (m3 H2O/m3 soil).

      REAL                                                              &
     & SMC(NPNTS)           ! OUT Soil moisture (kg/m2).

      REAL                                                              &
     & Z1,Z2                                                            &
                            ! WORK Depth of the top and bottom of the
!                           !      soil layers (m).
     &,ZSMC                 ! WORK Depth of layer for soil moisture
!                           !      diagnostic (m).
      PARAMETER ( ZSMC = 1. )

      INTEGER                                                           &
     & I,J,N                ! WORK Loop counters

#include "c_densty.h"

      DO I=1,NPNTS
        SMC(I) = 0.
      ENDDO

      Z2 = 0.
      DO N=1,NSHYD
        Z1 = Z2
        Z2 = Z2 + DZ(N)
        IF ( Z2 <  ZSMC ) THEN
!CDIR NODEP
          DO J=1,SOIL_PTS
            I = SOIL_INDEX(J)
            SMC(I) = SMC(I) + RHO_WATER * DZ(N) *                       &
     &                               ( STHU(I,N)*V_SAT(I) - V_WILT(I) )
          ENDDO
        ELSEIF ( Z2 >= ZSMC .AND. Z1 <  ZSMC ) THEN
!CDIR NODEP
          DO J=1,SOIL_PTS
            I = SOIL_INDEX(J)
            SMC(I) = SMC(I) + RHO_WATER * ( ZSMC - Z1 ) *               &
     &                               ( STHU(I,N)*V_SAT(I) - V_WILT(I) )
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE SOILMC
#endif
