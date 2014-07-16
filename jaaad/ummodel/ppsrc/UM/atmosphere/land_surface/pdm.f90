
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE PDM----------------------------------------------------
!
! Subroutine Interface:

      SUBROUTINE PDM(                                                   &
     &               NPNTS,SOIL_PTS,SOIL_INDEX,NSHYD,                   &
     &               TOT_TFALL,SNOW_MELT,INFILTRO,TIMESTEP,             &
     &               V_SAT,SATEXRO,STHU,STHF)

!
! Description:
!     Calculates the saturation excess runoff using PDM.
!     See Moore, 1985
!
! Documentation :
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.5    03/02/03   New Deck         Nic Gedney
!  6.1  17/08/04  Add SSFM code                           Ian Pearman
!  6.2  02/02/06  Remove SSFM.        P.Selwood
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!

      IMPLICIT NONE

! Global variables:
!---Soil layer thicknesses (m)
!---6 layers => CABLE else revert to MOSES
   REAL,PARAMETER:: DZSOIL(6) =(/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/)

      REAL RHO_WATER              ! DENSITY OF WATER (KG/M3)
!
      PARAMETER (RHO_WATER = 1000.0)
!
! C_PDM start
! 5.5 17/02/03    Required for large-scale hydrology L_PDM code.
!
! Soil layer thickness for PDM (m):
      REAL,PARAMETER :: DZ_PDM = 1.0
! Shape factor for PDM:
       REAL,PARAMETER :: B_PDM = 1.0
!
! C_PDM end

! Subroutine arguments
!   Scalar arguments with intent(IN) :

      INTEGER                                                           &
     & NPNTS                                                            &
                           ! IN Number of gridpoints.
     &,NSHYD                                                            &
                           ! IN Number of soil moisture levels.
     &,SOIL_PTS            ! IN Number of soil points.

      REAL                                                              &
     & TIMESTEP            ! IN Model timestep (s).

!   Array arguments with intent(IN) :

      INTEGER                                                           &
     & SOIL_INDEX(NPNTS)   ! IN Array of soil points.

      REAL                                                              &
     & TOT_TFALL(NPNTS)                                                 &
                            ! IN Cumulative canopy throughfall
                            !     (kg/m2/s).
     &,SNOW_MELT(NPNTS)                                                 &
                            ! IN GBM snow melt (kg/m2/s).
     &,INFILTRO(NPNTS)                                                  &
                            ! IN Surface runoff by infiltration
                            !    excess (kg/m2/s).
     &,V_SAT(NPNTS)                                                     &
                            ! IN Volumetric soil moisture
!                           !    concentration at saturation
!                           !    (m^3 water/m^3 soil).
     &,STHU(NPNTS,NSHYD)                                                &
                            ! IN Unfrozen soil moisture content of
!                           !    each layer as fraction of saturation.
     &,STHF(NPNTS,NSHYD)                                                &
                            ! IN Frozen soil moisture content of
!                           !    each layer as fraction of saturation.

!   Array arguments with intent(OUT) :
     &,SATEXRO(NPNTS)       ! OUT Sat. excess runoff (kg/m2/s).

! Local scalars:
      INTEGER                                                           &
     & I,J,N                ! WORK Loop counters.

      REAL                                                              &
     & P_IN                                                             &
                       ! WORK Water reaching soil surface (m^3/m^3).
     &,THCAP_MAX                                                        &
                       ! WORK Used in equation for TH_STAR (m^3/m^3).
     &,TH_STAR                                                          &
                       ! WORK Critical moisture storage (m^3/m^3).
     &,TH_STAR_P                                                        &
                       ! WORK Crit. TH_STAR after rain input (m^3/m^3).
     &,STH_PDM                                                          &
                       ! WORK Soil moisture for PDM / sat. value
     &,DZ_PDM_T                                                         &
                       ! WORK Soil depth temporary
     &,DZ_PDM_B        ! WORK Soil depth temporary
!

      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)

        P_IN = (TOT_TFALL(I)+SNOW_MELT(I)-INFILTRO(I))                  &
     &    *TIMESTEP/(DZ_PDM*RHO_WATER)

        IF (P_IN  >   0.0) THEN
!
! DZ_PDM_T is the residual soil depth for the PDM below the top
! of the current layer.
!
          DZ_PDM_T = DZ_PDM
          STH_PDM = 0.0

          DO N=1,NSHYD

            DZ_PDM_B = DZ_PDM_T - DZSOIL(N)
!
! DZ_PDM_B is the residual soil depth for the PDM below the
! bottom of the current layer
!
            IF (DZ_PDM_B  >=  0.0) THEN
              STH_PDM = STH_PDM + (STHU(I,N)+STHF(I,N))*DZSOIL(N)
            ELSE
              STH_PDM = STH_PDM+(STHU(I,N)+STHF(I,N))*MAX(0.,DZ_PDM_T)
            ENDIF
!
            DZ_PDM_T = DZ_PDM_B
!
          ENDDO ! Loop over soil layers
!
          STH_PDM = STH_PDM/DZ_PDM
!
          THCAP_MAX = (B_PDM+1.) * V_SAT(I)
          IF(STH_PDM  >   1.0) STH_PDM = 1.0
          TH_STAR = THCAP_MAX * (1.-(1.-STH_PDM)**(1/(B_PDM+1.)))
          TH_STAR_P = TH_STAR + P_IN

          IF(TH_STAR_P  <   THCAP_MAX) THEN
            SATEXRO(I) = P_IN - V_SAT(I)*                               &
     &        (1. - (1. - TH_STAR_P/THCAP_MAX)**(B_PDM+1.) - STH_PDM)
          ELSE
            SATEXRO(I) = P_IN - V_SAT(I)*(1. - STH_PDM)
          ENDIF
!
          IF (SATEXRO(I)  <   0.) SATEXRO(I) = 0.
!
          SATEXRO(I) = SATEXRO(I) / TIMESTEP * (DZ_PDM*RHO_WATER)

        ENDIF

      ENDDO


      RETURN
      END SUBROUTINE PDM

