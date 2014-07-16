#if defined(A08_7A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE SOIL_HTC------------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE SOIL_HTC (                                             &
     & NPNTS,NSHYD,NTILES,SOIL_PTS,SOIL_INDEX,TILE_PTS,TILE_INDEX       &
     &,BEXP,DZ,FRAC,HCAP,HCON,SATHH,SNOW_TILE,SURF_HT_FLUX,TIMESTEP     &
     &,V_SAT,W_FLUX,SMCL,STHU,STHF,TSOIL                                &
     &,LTIMER                                                           &
     &)

      IMPLICIT NONE
!
! Description:
!     Updates deep soil temperatures, frozen and unfrozen
!     frozen soil water content.
!
! Documentation : UM Documentation Paper 25
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.2    15/11/00   New Deck         M. Best
!  6.1  08/12/03  Add !CDIR NODEP to force vectorisation. R Barnes
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! System component covered: P25
! System Task: P25
!

! Global variables:
#include "c_soilh.h"

! Subroutine arguments
!   Scalar arguments with intent(IN) :
      INTEGER                                                           &
     & NPNTS                                                            &
                            ! IN Number of gridpoints.
     &,NSHYD                                                            &
                            ! IN Number of soil moisture levels.
     &,NTILES                                                           &
                            ! IN Number of tiles.
     &,SOIL_PTS             ! IN Number of soil points.

      REAL                                                              &
     & TIMESTEP             ! IN Model timestep (s).


!   Array arguments with intent(IN) :
      INTEGER                                                           &
     & SOIL_INDEX(NPNTS)                                                &
                            ! IN Array of soil points.
     &,TILE_PTS(NTILES)                                                 &
                            ! IN Number of tile points.
     &,TILE_INDEX(NPNTS,NTILES)
!                           ! IN Index of tile points.

      REAL                                                              &
     & BEXP(NPNTS)                                                      &
                            ! IN Clapp-Hornberger exponent.
     &,DZ(NSHYD)                                                        &
                            ! IN Thicknesses of the soil layers (m).
     &,FRAC(NPNTS,NTILES)                                               &
                            ! IN Tile fractions.
     &,HCAP(NPNTS)                                                      &
                            ! IN Soil heat capacity (J/K/m3).
     &,HCON(NPNTS)                                                      &
                            ! IN Soil thermal conductivity (W/m/K).
     &,SATHH(NPNTS)                                                     &
                            ! IN Saturated soil water pressure (m).
     &,SMCL(NPNTS,NSHYD)                                                &
                            ! IN Soil moisture content of each
!                           !    layer (kg/m2).
     &,SNOW_TILE(NPNTS,NTILES)                                          &
!                           ! IN Lying snow on tiles (kg/m2).
     &,SURF_HT_FLUX(NPNTS)                                              &
                            ! IN Net downward surface heat flux (W/m2).
     &,V_SAT(NPNTS)                                                     &
                            ! IN Volumetric soil moisture
!                           !    concentration at saturation
!                           !    (m3 H2O/m3 soil).
     &,W_FLUX(NPNTS,0:NSHYD)! IN The fluxes of water between layers
!                           !    (kg/m2/s).
!
      LOGICAL LTIMER        ! Logical switch for TIMER diags

!   Array arguments with intent(OUT) :

!   Array arguments with intent(INOUT) :
      REAL                                                              &
     & STHF(NPNTS,NSHYD)                                                &
                            ! INOUT Frozen soil moisture content of
!                           !       each layer as a fraction of
!                           !       saturation.
     &,STHU(NPNTS,NSHYD)                                                &
                            ! INOUT Unfrozen soil moisture content of
!                           !       each layer as a fraction of
!                           !       saturation.
     &,TSOIL(NPNTS,NSHYD)   ! INOUT Sub-surface temperatures (K).

!---------------------------------------------------------------------
!  External routines called :-

      EXTERNAL HEAT_CON,GAUSS
      EXTERNAL TIMER

!-----------------------------------------------------------------------
! Local scalars:
      INTEGER                                                           &
     & I,J,M,N                                                          &
                            ! WORK Loop counters.
     &,ITER_PTS                                                         &
                            ! WORK Number of soil points which require
!                           !      iteration.
     &,MMAX                                                             &
                            ! WORK Maximum number of iterations on
!                           !      temperature.
     &,FIRST_SOIL_PT        !      First soil point

      REAL                                                              &
     & FACUR                                                            &
                            ! WORK Required flux conservation accuracy
!                           !      (W/m2).
     &,SI_TILE                                                          &
                            ! WORK Tile snow insulation factor.
     &,SNOW_DEPTH                                                       &
                            ! WORK Depth of lying snow (m).
     &,TACUR                ! WORK Required accuracy of temperature
!                           !      calculation (Celsius).
      PARAMETER( MMAX=3, FACUR=0.01, TACUR=0.00000 )

! Local arrays:
      INTEGER                                                           &
     & ITER_INDEX(NPNTS)    ! WORK Array of soil points which require
!                           !      iteration.

      REAL                                                              &
     & CEACUR(NPNTS)                                                    &
                            ! WORK Flux conservation accuracy of the
!                           !      calculation (W/m2)
     &,DHSL0(NPNTS,NSHYD)                                               &
                            ! WORK Total heat increment to the layer
!                           !      (J/m2/timestep)
     &,DHSL(NPNTS,NSHYD)                                                &
                            ! WORK The heat available to update the
!                           !      layer temperature (J/m2/timestep)
     &,DSMCLF(NPNTS,NSHYD)                                              &
                            ! WORK The increment to the layer frozen
!                           !      soil moisture (kg/m2/timestep).
     &,DTHU(NPNTS,NSHYD)                                                &
                            ! WORK Rate of change of volumetric unfrozen
!                           !      soil moisture concentration with
!                           !      temperature (m3 liquid H2O/m3 soil/K)
     &,DTSL(NPNTS,NSHYD)                                                &
                            ! WORK The increment to the layer temperatur
!                           !      (K/timestep).
     &,DTSLMAX(NPNTS,NSHYD)                                             &
                            ! WORK Maximum value of DTSL (K/timestep).
     &,DTSLMIN(NPNTS,NSHYD)                                             &
                            ! WORK Minimum value of DTSL (K/timestep).
     &,HCAPT(NPNTS,NSHYD)                                               &
                           ! WORK The total volumetric heat capacity
!                           !      (soil+water) of the layer (J/m3/K).
     &,HC(NPNTS,NSHYD)                                                  &
                            ! WORK The thermal conductivity of each
!                           !      layer (W/m/K).
     &,HCONS(NPNTS)                                                     &
                            ! WORK The thermal conductivity between
!                           !      adjacent soil layers (W/m/K).
     &,H_FLUX(NPNTS,0:NSHYD)                                            &
                             !WORK The fluxes of heat between layers
!                           !      (W/m2).
     &,HADV(NPNTS,NSHYD)                                                &
                            ! WORK Heat flux due to moisture advection
!                           !      (W/m2).
     &,SIFACT(NPNTS)                                                    &
                            ! WORK Snow insulation factor.
     &,SMCLF(NPNTS,NSHYD)                                               &
                            ! WORK Frozen moisture content of each
!                           !      soil layer (kg/m2).
     &,SMCLF0(NPNTS,NSHYD)                                              &
                            ! WORK Previous value of SMCLF (kg/m2).
     &,SMCLSAT(NPNTS,NSHYD)                                             &
                            ! WORK The saturation moisture content of
!                           !      each layer (kg/m2).
     &,SMCLU(NPNTS,NSHYD)                                               &
                            ! WORK Unfrozen moisture content of each
!                           !      soil layer (kg/m2).
     &,SMCLU0(NPNTS,NSHYD)                                              &
                            ! WORK Previous value of SMCLU (kg/m2).
     &,SMCFU(NPNTS)                                                     &
                            ! WORK Fractional saturation (unfrozen water
!                           !      at layer boundaries.
     &,SMCFF(NPNTS)                                                     &
                            ! WORK Fractional saturation (frozen water)
!                           !      at layer boundaries.
     &,TMAX(NPNTS,NSHYD)                                                &
                            ! WORK Temperature above which all water is
!                           !      unfrozen (Celsius)
     &,TSL(NPNTS,0:NSHYD)                                               &
                            ! WORK Soil layer temperatures (Celsius)
!                           !      TSL(0) temperature of incoming water.
!                           !      TSL(1:NSHYD) sub-surface soil
!                           !      temperatures .
     &,TSL0(NPNTS,0:NSHYD)  ! WORK Previous value of TSL (Celsius).
      LOGICAL                                                           &
     & ITER(NPNTS)          ! WORK .T. on points requiring iterations.

!-----------------------------------------------------------------------
! Variables required for the implicit calculation.
!-----------------------------------------------------------------------
      REAL                                                              &
     & DHFLUX_DTSL1(NPNTS,0:NSHYD),DHFLUX_DTSL2(NPNTS,0:NSHYD)          &
     &,DHADV_DTSL0(NPNTS,NSHYD),DHADV_DTSL1(NPNTS,NSHYD)                &
     &,DHADV_DTSL2(NPNTS,NSHYD)                                         &
                                ! WORK Rate of change of the explicit
!                           ! fluxes with the layer temperatures
!                           ! (W/m2/K).
     &,A(NPNTS,NSHYD),B(NPNTS,NSHYD),C(NPNTS,NSHYD),D(NPNTS,NSHYD)      &
!                           ! WORK Matrix elements.
     &,GAMCON               ! WORK Forward timestep weighting constant.
! Local parameters:
      REAL                                                              &
     & GAMMA                ! Forward timestep weighting.
      PARAMETER (GAMMA=1.0)


#include "c_densty.h"
#include "c_lheat.h"
#include "c_perma.h"
#include "c_0_dg_c.h"

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SOILHTC ',103)
      ENDIF

!-----------------------------------------------------------------------
! Initialisations
!-----------------------------------------------------------------------
      FIRST_SOIL_PT=SOIL_INDEX(1)

      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        TSL(I,0)=TSOIL(I,1)-ZERODEGC
      ENDDO

      DO N=1,NSHYD
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
!-----------------------------------------------------------------------
! Define soil layer temperatures TSL (in celsius).
!-----------------------------------------------------------------------
          TSL(I,N)=TSOIL(I,N)-ZERODEGC
        ENDDO
      ENDDO

      DO N=1,NSHYD
! CDIR$ IVDEP here would force vectorization but changes results!
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
!-----------------------------------------------------------------------
! Diagnose the frozen and unfrozen water.
!-----------------------------------------------------------------------

          SMCLSAT(I,N)=RHO_WATER*DZ(N)*V_SAT(I)
          SMCLF(I,N)=SMCLSAT(I,N)*STHF(I,N)
          SMCLU(I,N)=SMCL(I,N)-SMCLF(I,N)
        ENDDO                                  !  J=1,SOIL_PTS
      ENDDO                                    ! N=1,NSHYD
!      print 11,SMCL,SMCLU,STHU,SMCLF,STHF,SMCLSAT
!11    format(1x,'chekingsthu',25f7.2)


!-----------------------------------------------------------------------
! Initialise the array of points for which calculations are required.
!-----------------------------------------------------------------------
      DO I=1,NPNTS
        ITER(I)=.FALSE.
      ENDDO

      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        ITER(I)=.TRUE.
      ENDDO

!----------------------------------------------------------------------
! Calculate the heat conductivity between adjacent layers
!----------------------------------------------------------------------

      DO N=1,NSHYD-1
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          SMCFU(I)=(DZ(N+1)*STHU(I,N)+DZ(N)*STHU(I,N+1))                &
     &              /(DZ(N+1)+DZ(N))
          SMCFF(I)=(DZ(N+1)*STHF(I,N)+DZ(N)*STHF(I,N+1))                &
     &              /(DZ(N+1)+DZ(N))
        ENDDO
! DEPENDS ON: heat_con
        CALL HEAT_CON(NPNTS-FIRST_SOIL_PT+1,HCON(FIRST_SOIL_PT)         &
     &,              SMCFU(FIRST_SOIL_PT),SMCFF(FIRST_SOIL_PT)          &
     &,              V_SAT(FIRST_SOIL_PT),HCONS(FIRST_SOIL_PT),LTIMER)
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          HC(I,N)=HCONS(I)
        ENDDO
      ENDDO

!--------------------------------------------------------------------
! Calculate the snow insulation factor
!--------------------------------------------------------------------
      DO I=1,NPNTS
        SIFACT(I) = 0.
      ENDDO
      DO N=1,NTILES
        DO J=1,TILE_PTS(N)
          I = TILE_INDEX(J,N)
          SI_TILE = 1.
          IF (SNOW_TILE(I,N) >  0. .AND. V_SAT(I) /= 0.) THEN
            SNOW_DEPTH = SNOW_TILE(I,N) / RHO_SNOW
            IF (SNOW_DEPTH  <=  (0.5*DZ(1))) THEN
              SI_TILE = 1. / ( 1. + 2*SNOW_DEPTH/(DZ(1) + DZ(2)) )
            ELSE
              SI_TILE =(DZ(1) + DZ(2)) /                                &
     &                 ( HC(I,1)*(2*SNOW_DEPTH - DZ(1))/SNOW_HCON       &
     &                   + 2*DZ(1) + DZ(2) )
            ENDIF
          ENDIF
          SIFACT(I) = SIFACT(I) + FRAC(I,N)*SI_TILE
        ENDDO
      ENDDO

!--------------------------------------------------------------------
! Calculate heat fluxes across layer boundaries
!--------------------------------------------------------------------
      DO N=1,NSHYD-1
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          H_FLUX(I,N)=-HC(I,N)*2.0*(TSL(I,N+1)-TSL(I,N))                &
     &                            /(DZ(N+1)+DZ(N))
          DHFLUX_DTSL1(I,N)=HC(I,N)*2.0/(DZ(N+1)+DZ(N))
          DHFLUX_DTSL2(I,N)=-HC(I,N)*2.0/(DZ(N+1)+DZ(N))
        ENDDO
      ENDDO
!DIR$ IVDEP
!CDIR NODEP
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        H_FLUX(I,NSHYD)=0.0
        H_FLUX(I,0) = SURF_HT_FLUX(I)
        H_FLUX(I,1) = SIFACT(I)*H_FLUX(I,1)
        DHFLUX_DTSL1(I,NSHYD)=0.0
        DHFLUX_DTSL2(I,NSHYD)=0.0
        DHFLUX_DTSL1(I,0)=0.0
        DHFLUX_DTSL2(I,0)=0.0
        DHFLUX_DTSL1(I,1)=SIFACT(I)*DHFLUX_DTSL1(I,1)
        DHFLUX_DTSL2(I,1)=SIFACT(I)*DHFLUX_DTSL2(I,1)
      ENDDO

!--------------------------------------------------------------------
! Calculate the advection of heat by moisture fluxes
!--------------------------------------------------------------------
      DO N=2,NSHYD-1
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          HADV(I,N)=HCAPW*DZ(N)*                                        &
     &    (W_FLUX(I,N-1)*(TSL(I,N-1)-TSL(I,N))/(DZ(N)+DZ(N-1))          &
     &     +W_FLUX(I,N)*(TSL(I,N)-TSL(I,N+1))/(DZ(N)+DZ(N+1)))
          DHADV_DTSL0(I,N)=HCAPW*DZ(N)*W_FLUX(I,N-1)/(DZ(N)+DZ(N-1))
          DHADV_DTSL1(I,N)=HCAPW*DZ(N)*                                 &
     &                     (-W_FLUX(I,N-1)/(DZ(N)+DZ(N-1))              &
     &                      +W_FLUX(I,N)/(DZ(N)+DZ(N+1)))
          DHADV_DTSL2(I,N)=-HCAPW*DZ(N)*W_FLUX(I,N)/(DZ(N)+DZ(N+1))
        ENDDO
      ENDDO

!DIR$ IVDEP
!CDIR NODEP
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        HADV(I,1)=HCAPW*DZ(1)*                                          &
     &   (W_FLUX(I,0)*(TSL(I,0)-TSL(I,1))/DZ(1)                         &
     &   +W_FLUX(I,1)*(TSL(I,1)-TSL(I,2))/(DZ(1)+DZ(2)))
        DHADV_DTSL0(I,1)=0.0
        DHADV_DTSL1(I,1)=HCAPW*DZ(1)*                                   &
     &    (-W_FLUX(I,0)/DZ(1)+W_FLUX(I,1)/(DZ(1)+DZ(2)))
        DHADV_DTSL2(I,1)=-HCAPW*DZ(1)*W_FLUX(I,1)/(DZ(1)+DZ(2))
        HADV(I,NSHYD)=HCAPW*DZ(NSHYD)*                                  &
     &               W_FLUX(I,NSHYD-1)*(TSL(I,NSHYD-1)-TSL(I,NSHYD))    &
     &               /(DZ(NSHYD)+DZ(NSHYD-1))
        DHADV_DTSL0(I,NSHYD)=HCAPW*DZ(NSHYD)*W_FLUX(I,NSHYD-1)          &
     &                       /(DZ(NSHYD)+DZ(NSHYD-1))
        DHADV_DTSL1(I,NSHYD)=-HCAPW*DZ(NSHYD)*W_FLUX(I,NSHYD-1)         &
     &                       /(DZ(NSHYD)+DZ(NSHYD-1))
        DHADV_DTSL2(I,NSHYD)=0.0
      ENDDO

      DO N=1,NSHYD
!DIR$ IVDEP
!CDIR NODEP
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
!-----------------------------------------------------------------------
! Calculate TMAX, the temperature above which all soil water is
! unfrozen. Check that (SMCLSAT/SMCL)**BEXP will not overflow when SMCL
! is very small. The function EPSILON  gives the number of type (real)
! of SMCL that is negligeable compared to 1.
!-----------------------------------------------------------------------
#if defined(SCMA) && !defined(T3E)
          IF (SMCL(I,N)  >   SMCLSAT(I,N)*1E-20) THEN
#else
          IF (SMCL(I,N)  >   SMCLSAT(I,N)*(EPSILON(SMCL(I,N)) ) )THEN
#endif

            TMAX(I,N)=-SATHH(I)/DPSIDT                                  &
     &             *(SMCLSAT(I,N)/SMCL(I,N))**(BEXP(I))
          ELSE
            TMAX(I,N)=-ZERODEGC
          ENDIF
          TMAX(I,N)=MAX(TMAX(I,N),-ZERODEGC)

          DHSL0(I,N)=TIMESTEP*(H_FLUX(I,N-1)-H_FLUX(I,N)+HADV(I,N))

          DHSL(I,N)=DHSL0(I,N)

        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Iteration loop
!-----------------------------------------------------------------------
      DO M=1,MMAX


!-----------------------------------------------------------------------
! Define the array of points which fail to meet the flux criterion.
!-----------------------------------------------------------------------
      ITER_PTS=0
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)

        IF (ITER(I)) THEN
          ITER_PTS=ITER_PTS+1
          ITER_INDEX(ITER_PTS)=I
        ENDIF
        ITER(I)=.FALSE.

      ENDDO

!-----------------------------------------------------------------------
! Update calculations at these points.
!-----------------------------------------------------------------------
      DO N=1,NSHYD

!CDIR NODEP
        DO J=1,ITER_PTS
          I=ITER_INDEX(J)

          TSL0(I,N)=TSL(I,N)
          SMCLF0(I,N)=SMCLF(I,N)
          SMCLU0(I,N)=SMCLU(I,N)

          DTSLMAX(I,N)=1.0E4-TSL(I,N)
          DTSLMIN(I,N)=-ZERODEGC-TSL(I,N)

          IF (TSL(I,N) >  TMAX(I,N)) THEN         ! All water unfrozen
            DTHU(I,N)=0.0
            DTSLMIN(I,N)=TMAX(I,N)-TSL(I,N)
          ELSEIF (TSL(I,N) == TMAX(I,N).AND.                            &
                                                  ! Remains unfrozen
     &            DHSL(I,N) >= 0.0) THEN
            DTHU(I,N)=0.0
          ELSE                                    ! Phase changes
            DTHU(I,N)=DPSIDT*SMCLSAT(I,N)                               &
     &               /(BEXP(I)*SATHH(I)*RHO_WATER*DZ(N))                &
     &            *(-DPSIDT*TSL(I,N)/SATHH(I))**(-1.0/BEXP(I)-1.0)
            DTSLMAX(I,N)=TMAX(I,N)-TSL(I,N)
          ENDIF

          HCAPT(I,N)=HCAP(I)+(HCAPW-HCAPI)*SMCLU(I,N)/DZ(N)             &
     &            +HCAPI*SMCL(I,N)/DZ(N)                                &
     &            +RHO_WATER*DTHU(I,N)*((HCAPW-HCAPI)*TSL(I,N)+LF)


!-----------------------------------------------------------------------
! Calculate the matrix elements required for the implicit update.
!-----------------------------------------------------------------------
          GAMCON=GAMMA*TIMESTEP/(HCAPT(I,N)*DZ(N))
          A(I,N)=-GAMCON*(DHFLUX_DTSL1(I,N-1)+DHADV_DTSL0(I,N))
          B(I,N)=1.0-GAMCON*(DHFLUX_DTSL2(I,N-1)-DHFLUX_DTSL1(I,N)      &
     &                                          +DHADV_DTSL1(I,N))
          C(I,N)=GAMCON*(DHFLUX_DTSL2(I,N)+DHADV_DTSL2(I,N))
          D(I,N)=1.0/(HCAPT(I,N)*DZ(N))*DHSL(I,N)

        ENDDO
      ENDDO


!-----------------------------------------------------------------------
! Solve the triadiagonal matrix equation.
!-----------------------------------------------------------------------
! DEPENDS ON: gauss
      CALL GAUSS(NSHYD,NPNTS,ITER_PTS,ITER_INDEX,A,B,C,D                &
     &,          DTSLMIN,DTSLMAX,DTSL)

!-----------------------------------------------------------------------
! Diagnose the implicit DHSL
!-----------------------------------------------------------------------
      DO N=2,NSHYD-1
        DO J=1,ITER_PTS
          I=ITER_INDEX(J)
          DHSL(I,N)=DHSL(I,N)-DZ(N)*HCAPT(I,N)*(A(I,N)*DTSL(I,N-1)      &
     &                    +(B(I,N)-1)*DTSL(I,N)+C(I,N)*DTSL(I,N+1))
        ENDDO
      ENDDO

      DO J=1,ITER_PTS
        I=ITER_INDEX(J)
        DHSL(I,1)=DHSL(I,1)-DZ(1)*HCAPT(I,1)*(                          &
     &                  +(B(I,1)-1)*DTSL(I,1)+C(I,1)*DTSL(I,2))
        DHSL(I,NSHYD)=DHSL(I,NSHYD)-DZ(NSHYD)*HCAPT(I,NSHYD)*           &
     &  (A(I,NSHYD)*DTSL(I,NSHYD-1)+(B(I,NSHYD)-1)*DTSL(I,NSHYD))
      ENDDO

!-----------------------------------------------------------------------
! Update the layer temperatures
!-----------------------------------------------------------------------
      DO N=1,NSHYD
!CDIR NODEP
        DO J=1,ITER_PTS
          I=ITER_INDEX(J)

            TSL(I,N)=TSL(I,N)+DTSL(I,N)

!-----------------------------------------------------------------------
! If the temperature increment is small and frozen water exists
! assume that the excess energy goes into phase change
!-----------------------------------------------------------------------
            IF (ABS(DTSL(I,N)) <  TACUR.AND.                            &
     &          TSL(I,N) <= TMAX(I,N)) THEN
              DSMCLF(I,N)=-DHSL(I,N)/((HCAPW-HCAPI)*TSL(I,N)+LF)
              DSMCLF(I,N)=MAX(DSMCLF(I,N),-SMCLF(I,N))
              DSMCLF(I,N)=MIN(DSMCLF(I,N),SMCLU(I,N))
              SMCLU(I,N)=SMCLU(I,N)-DSMCLF(I,N)
              SMCLF(I,N)=SMCLF(I,N)+DSMCLF(I,N)
            ENDIF

!-----------------------------------------------------------------------
! Diagnose unfrozen and frozen water contents
!-----------------------------------------------------------------------
            IF (TSL(I,N) >= TMAX(I,N)) THEN
              SMCLU(I,N)=SMCL(I,N)
              SMCLF(I,N)=0.0
            ELSE
              SMCLU(I,N)=SMCLSAT(I,N)                                   &
     &                *(-DPSIDT*TSL(I,N)/SATHH(I))**(-1.0/BEXP(I))
              SMCLF(I,N)=SMCL(I,N)-SMCLU(I,N)
            ENDIF

!-----------------------------------------------------------------------
! Calculate the error in heat conservation
!-----------------------------------------------------------------------
            DSMCLF(I,N)=SMCLF(I,N)-SMCLF0(I,N)
            DHSL(I,N)=DHSL(I,N)-(HCAP(I)*DZ(N)+HCAPW*SMCLU0(I,N)        &
     &                          +HCAPI*SMCLF0(I,N))*DTSL(I,N)           &
     &                     -DSMCLF(I,N)*((HCAPI-HCAPW)*TSL0(I,N)-LF)

!-----------------------------------------------------------------------
! Calculate the error in flux conservation
!-----------------------------------------------------------------------
          CEACUR(I)=ABS(DHSL(I,N))/TIMESTEP

          IF (CEACUR(I)  >   FACUR) THEN
            ITER(I)=.TRUE.
          ENDIF

        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! End of iteration loop
!-----------------------------------------------------------------------
      ENDDO

!-----------------------------------------------------------------------
! Diagnose soil temperatures (K) and fractional values of unfrozen and
! frozen water.
!-----------------------------------------------------------------------

      DO N=1,NSHYD
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          TSOIL(I,N)=TSL(I,N)+ZERODEGC
          STHU(I,N)=SMCLU(I,N)/SMCLSAT(I,N)
          STHF(I,N)=SMCLF(I,N)/SMCLSAT(I,N)
        ENDDO
      ENDDO
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SOILHTC ',104)
      ENDIF

      RETURN
      END SUBROUTINE SOIL_HTC
#endif
