#if defined(A19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Version 2A of vegetation section: models leaf phenology and vegetation
! competition
!
! Subroutine Interface:
      SUBROUTINE VEG(LAND_PTS,LAND_INDEX,NTILES,CAN_MODEL               &
     &,              row_length, rows                                   &
     &,              A_STEP,ASTEPS_SINCE_TRIFFID                        &
     &,              PHENOL_PERIOD,TRIFFID_PERIOD                       &
     &,              L_PHENOL,L_TRIFFID,L_TRIF_EQ                       &
     &,              ATIMESTEP,FRACA,SATCON                             &
     &,              G_LEAF_AC,G_LEAF_PHEN_AC,NPP_AC                    &
     &,              RESP_S_AC,RESP_W_AC                                &
     &,              CS,FRAC,LAI,CLAY_FRAC,HT                           &
     &,              CATCH_S,CATCH_T,INFIL_T,Z0_T                       &
     &,              C_VEG,CV,LIT_C,LIT_C_MN,G_LEAF_DAY,G_LEAF_PHEN     &
     &,              LAI_PHEN,G_LEAF_DR_OUT,NPP_DR_OUT,RESP_W_DR_OUT    &
     &,              RESP_S_DR_OUT                                      &
     &               )


      IMPLICIT NONE
!
! Description:
!   Updates Leaf Area Index for Plant Functional Types (PFTs) and uses
!   this to derive new vegetation parameters for PFTs along with gridbox
!   mean values where appropriate.
!
! Method:
!   Calls PHENOL which models phenolgy and updates Leaf Area Index
!   (LAI), then calls TRIFFID to update vegetation and soil fractions,
!   LAI, canopy height, veg and soil carbon and carbon fluxes.  Passes
!   fractions, LAI and canopy height to SPARM which derives the
!   vegetation parameters for each PFT and also the gridbox means where
!   this is required.
!
! Current Code Owner:  Richard Betts
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   4.4    8/10/97   Original code.  Richard Betts
!   4.5   12/05/98   Find total fraction of gridbox covered by
!                    vegetation or soil, use this to derive indices of
!                    land points on which TRIFFID may operate, and pass
!                    both to TRIFFID.  Initialise top and bottom rows
!                    for all variables.  Richard Betts
!   4.5   30/06/98   Add second call to TILEPTS to update TILE_INDEX
!                    after TRIFFID.  Richard Betts
!   4.5    6/08/98   Call SWAPB_LAND to update halo regions of input
!                    fields.   Richard Betts
!   4.5   23/11/98   Output G_LEAF_DAY, G_LEAF_PHEN, LAI_PHEN,
!                    G_LEAF_DR_OUT, NPP_DR_OUT, RESP_W_DR_OUT and
!                    RESP_S_DR_OUT as diagnostics.  Richard Betts
!   5.2   15/11/00   Re-Written for New Dynamics      M. Best
!   5.4   11/04/02   Canopy snow model for needleleaf tree tile.
!                    R. Essery
!  6.2  01/03/06  replace single pool soil carbon model with RothC.
!                                                           C.D. Jones
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.


      INTEGER                                                           &
     & LAND_PTS                                                         &
                             ! IN Number of land points.
     &,NTILES                                                           &
                             ! IN Number of land-surface tiles.
     &,CAN_MODEL                                                        &
                                    ! IN Swith for thermal vegetation
     &,row_length                                                       &
                        ! number of points on a row
     &,rows                                                             &
                        ! number of rows in a theta field
     &,A_STEP                                                           &
                             ! IN Atmospheric timestep number.
     &,ASTEPS_SINCE_TRIFFID                                             &
                             ! INOUT Number of atmosphere
!                                    timesteps since last call
!                                           to TRIFFID.
     &,PHENOL_PERIOD                                                    &
                             ! IN Phenology period (days).
     &,TRIFFID_PERIOD        ! IN TRIFFID period (days).

#include "nstypes.h"

      INTEGER                                                           &
     & LAND_INDEX(LAND_PTS)         ! IN I=LAND_INDEX(L) => the Ith
!                                   !    P-point is the Lth land
!                                   !    point.

      INTEGER                                                           &
     & I,J,K,L,N                                                        &
                                    ! WORK loop counters.
     &,KITER                        ! WORK Number of TRIFFID iterations.

      LOGICAL                                                           &
     & L_PHENOL                                                         &
                                    ! IN .T. for interactive leaf
!                                   !    phenology.
     &,L_TRIFFID                                                        &
                                    ! IN .T. for interactive vegetation.
     &,L_TRIF_EQ                    ! IN .T. for vegetation equilibrium.

      REAL                                                              &
     & ATIMESTEP                                                        &
                                    ! IN Atmospheric timestep (s).
     &,FRACA(LAND_PTS)                                                  &
                                    ! IN Fraction of agriculture.
     &,SATCON(LAND_PTS)                                                 &
                                    ! IN Saturated hydraulic
!                                   !    conductivity (kg/m2/s).
     &,G_LEAF_AC(LAND_PTS,NPFT)                                         &
                                    ! INOUT Accumulated leaf turnover
!                                   !       rate.
     &,G_LEAF_PHEN_AC(LAND_PTS,NPFT)                                    &
                                    ! INOUT Accumulated leaf turnover
!                                   !       rate including phenology.
     &,NPP_AC(LAND_PTS,NPFT)                                            &
                                    ! INOUT Accumulated NPP (kg C/m2).
     &,RESP_W_AC(LAND_PTS,NPFT)                                         &
                                    ! INOUT Accumulated wood respiration
!                                   !       (kg C/m2).
     &,RESP_S_AC(LAND_PTS,4)                                            &
                                 ! INOUT Accumulated soil respiration
!                                       !       (kg C/m2).
     &,CS(LAND_PTS,4)                                                   &
                                 ! INOUT Soil carbon content
     &,CS_TOT(LAND_PTS)                                                 &
                                 ! INOUT Soil carbon content
!                                   !       (kg C/m2).
     &,FRAC(LAND_PTS,NTYPE)                                             &
                                    ! INOUT Fractions of surface types.
     &,LAI(LAND_PTS,NPFT)                                               &
                                    ! INOUT LAI of plant functional
!                                   !       types.
     &,HT(LAND_PTS,NPFT)                                                &
                                    ! INOUT Height of plant functional
!                                   !       types (m).
     &,CATCH_S(LAND_PTS)                                                &
                                    ! OUT Snow capacity for NLT tile
!                                   !     (kg/m2).
     &,CATCH_T(LAND_PTS,NTILES)                                         &
                                    ! OUT Canopy capacity for tiles
!                                   !     (kg/m2).
     &,INFIL_T(LAND_PTS,NTILES)                                         &
                                    ! OUT Maximum surface infiltration
!                                   !     rate for tiles (kg/m2/s).
     &,Z0_T(LAND_PTS,NTILES)                                            &
                                    ! OUT Roughness length for tiles (m)
     &,C_VEG(LAND_PTS,NPFT)                                             &
                                    ! OUT Total carbon content of
!                                   !     the vegetation (kg C/m2).
     &,CV(LAND_PTS)                                                     &
                                    ! OUT Gridbox mean vegetation
!                                   !     carbon (kg C/m2).
     &,G_LEAF_DAY(LAND_PTS,NPFT)                                        &
                                    ! OUT Mean leaf turnover rate for
!                                   !      input to PHENOL (/360days).
     &,G_LEAF_DR_OUT(LAND_PTS,NPFT)                                     &
                                    ! OUT Mean leaf turnover rate for
!                                   !       driving TRIFFID (/360days).
     &,LAI_PHEN(LAND_PTS,NPFT)                                          &
                                    ! OUT LAI of PFTs after phenology.
     &,LIT_C(LAND_PTS,NPFT)                                             &
                                    ! OUT Carbon Litter
!                                   !     (kg C/m2/360days).
     &,LIT_C_MN(LAND_PTS)                                               &
                                    ! OUT Gridbox mean carbon litter
!                                   !     (kg C/m2/360days).
     &,NPP_DR_OUT(LAND_PTS,NPFT)                                        &
                                    ! OUT Mean NPP for driving TRIFFID
!                                   !     (kg C/m2/360days).
     &,RESP_W_DR_OUT(LAND_PTS,NPFT)                                     &
                                    ! OUT Mean wood respiration for
!                                   !       driving TRIFFID
!                                   !       (kg C/m2/360days).
     &,RESP_S_DR_OUT(LAND_PTS,5)                                        &
                                 ! OUT Mean soil resp for
     &,CLAY_FRAC(ROW_LENGTH, ROWS)                                      &
                                        ! IN  Clay fraction of soil
     &,CLAY_LAND(LAND_PTS)                                              &
                                        ! IN  clay on land points
     &,RESP_FRAC(LAND_PTS)            ! respired fraction of RESP_S
!                                   !     driving TRIFFID
!                                   !     (kg C/m2/360days).

      INTEGER                                                           &
     & NSTEP_PHEN                                                       &
                                    ! WORK Number of atmospheric
!                                   !      timesteps between calls to
!                                   !      PHENOL.
     &,NSTEP_TRIF                                                       &
                                    ! WORK Number of atmospheric
!                                   !      timesteps between calls to
!                                   !      TRIFFID.
     &,TILE_PTS(NTYPE)                                                  &
                                    ! WORK Number of land points which
!                                   !      include the nth surface type.
     &,TILE_INDEX(LAND_PTS,NTYPE)                                       &
                                    ! WORK Indices of land points which
!                                   !      include the nth surface type.
     &,TRIF_PTS                                                         &
                                    ! WORK Number of points on which
!                                   !      TRIFFID may operate
     &,TRIF_INDEX(LAND_PTS)         ! WORK Indices of land points on
!                                   !      which TRIFFID may operate

      REAL                                                              &
     & DTIME_PHEN                                                       &
                                    ! WORK The phenology timestep (yr).
     &,FORW                                                             &
                                    ! WORK Forward timestep weighting
!                                   !      for TRIFFID.
     &,GAMMA                                                            &
                                    ! WORK Inverse TRIFFID timestep
!                                   !      (/360days).
     &,GAM_TRIF                                                         &
                                    ! WORK Inverse TRIFFID coupling
!                                   !      timestep (/360days).
     &,RATIO                                                            &
                                    ! WORK Ratio of fractional
!                                   !      coverage before to that
!                                   !      after TRIFFID.
     &,DCS(LAND_PTS)                                                    &
                                    ! WORK Change in soil carbon
!                                   !      (kg C/m2).
     &,FRAC_AGRIC(LAND_PTS)                                             &
                                    ! WORK Fraction of agriculture as
!                                   !      seen by TRIFFID.
     &,FRAC_OLD(LAND_PTS,NTYPE)                                         &
                                    ! WORK Fractions of surface types
!                                   !      before the call to TRIFFID.
     &,FRAC_VS(LAND_PTS)                                                &
                                    ! WORK Total fraction of gridbox
!                                   !      covered by veg or soil.
     &,G_LEAF_PHEN(LAND_PTS,NPFT)                                       &
                                    ! WORK Mean leaf turnover rate over
!                                   !      phenology period (/360days).
     &,G_LEAF_DR(LAND_PTS,NPFT)                                         &
                                    ! WORK Mean leaf turnover rate
!                                   !      for driving TRIFFID
!                                   !      (/360days).
     &,NPP_DR(LAND_PTS,NPFT)                                            &
                                    ! WORK Mean NPP for driving
!                                   !      TRIFFID (kg C/m2/360days).
     &,RESP_W_DR(LAND_PTS,NPFT)                                         &
                                    ! WORK Mean wood respiration for
!                                   !      driving TRIFFID
!                                   !      (kg C/m2/360days).
     &,RESP_S_DR(LAND_PTS,5)     ! WORK Mean soil resp for
!                                   !      driving TRIFFID
!                                   !      (kg C/m2/360days).
!-----------------------------------------------------------------------
! Local parameters
!-----------------------------------------------------------------------
#include "descent.h"
#include "seed.h"

      LOGICAL                                                           &
     & AGRIC                        ! .T. for TRIFFID to see agriculture
!                                   ! .F. for natural vegetation.
      PARAMETER (AGRIC=.TRUE.)

!-----------------------------------------------------------------------
! Initialisations
!-----------------------------------------------------------------------
      DO N=1,NPFT
        DO L=1,LAND_PTS
          G_LEAF_PHEN(L,N)=0.0
          G_LEAF_DAY(L,N)=0.0
          G_LEAF_DR(L,N)=0.0
          NPP_DR(L,N)=0.0
          RESP_W_DR(L,N)=0.0
          C_VEG(L,N)=0.0
          LIT_C(L,N)=0.0
        ENDDO
      ENDDO

      DO N=1,NTILES
        DO L=1,LAND_PTS
          CATCH_T(L,N)=0.0
          INFIL_T(L,N)=0.0
          Z0_T(L,N)=0.0
        ENDDO
      ENDDO
      
      IF (CAN_MODEL == 4) THEN
        DO L=1,LAND_PTS
          CATCH_S(L)=0.0
        ENDDO
      ENDIF
      
      DO L=1,LAND_PTS
        CV(L)=0.0
        LIT_C_MN(L)=0.0
        FRAC_VS(L) = 0.0
        FRAC_AGRIC(L) = 0.0
      ENDDO
      DO N=1,5
        DO L=1,LAND_PTS
          RESP_S_DR(L,N)=0.0
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Calculate the number of atmospheric timesteps between calls to PHENOL
! and TRIFFID.
!-----------------------------------------------------------------------
      NSTEP_PHEN=INT(86400.0*PHENOL_PERIOD/ATIMESTEP)
      NSTEP_TRIF=INT(86400.0*TRIFFID_PERIOD/ATIMESTEP)


!-----------------------------------------------------------------------
! Find total fraction of gridbox covered by vegetation and soil, and use
! this to set indices of land points on which TRIFFID may operate
!-----------------------------------------------------------------------
      TRIF_PTS = 0
      DO L=1,LAND_PTS
        DO N=1,NPFT
          FRAC_VS(L) = FRAC_VS(L) + FRAC(L,N)
        ENDDO
        N=SOIL
        FRAC_VS(L) = FRAC_VS(L) + FRAC(L,N)
        IF (FRAC_VS(L) >= (NPFT*FRAC_MIN)) THEN
          TRIF_PTS = TRIF_PTS + 1
          TRIF_INDEX(TRIF_PTS) = L
        ENDIF
      ENDDO

!-----------------------------------------------------------------------
! Create the TILE_INDEX array of land points with each surface type
!-----------------------------------------------------------------------
! DEPENDS ON: tilepts
      CALL TILEPTS(LAND_PTS,FRAC,TILE_PTS,TILE_INDEX)

      IF (L_PHENOL .AND. MOD(A_STEP,NSTEP_PHEN) == 0) THEN

!-----------------------------------------------------------------------
! Calculate the phenology timestep in years.
!-----------------------------------------------------------------------
        DTIME_PHEN=FLOAT(PHENOL_PERIOD)/360.0

        DO N=1,NPFT

!-----------------------------------------------------------------------
! Calculate the mean turnover rate and update the leaf phenological
! state, and take copy of updated LAI field for output as diagnostic.
!-----------------------------------------------------------------------
          DO J=1,TILE_PTS(N)
            L=TILE_INDEX(J,N)
            G_LEAF_DAY(L,N)=G_LEAF_AC(L,N)/DTIME_PHEN
          ENDDO

          WRITE(6,*) 'Calling phenology'

! DEPENDS ON: phenol
          CALL PHENOL (LAND_PTS,TILE_PTS(N),TILE_INDEX(1,N),N,          &
     &                 G_LEAF_DAY(1,N),HT(1,N),DTIME_PHEN,              &
     &                 G_LEAF_PHEN(1,N),LAI(1,N))

          WRITE(6,*) 'Phenology completed normally'

          DO L=1,LAND_PTS
            LAI_PHEN(L,N)=LAI(L,N)
          ENDDO

!-----------------------------------------------------------------------
! Increment the leaf turnover rate for driving TRIFFID and reset the
! accumulation over atmospheric model timesteps to zero.
!-----------------------------------------------------------------------
          DO J=1,TILE_PTS(N)
            L=TILE_INDEX(J,N)
            G_LEAF_PHEN_AC(L,N)=G_LEAF_PHEN_AC(L,N)                     &
     &                    +G_LEAF_PHEN(L,N)*DTIME_PHEN
          ENDDO

          DO L=1,LAND_PTS
            G_LEAF_AC(L,N)=0.0
          ENDDO

        ENDDO
      ENDIF

!-----------------------------------------------------------------------
! Call TRIFFID vegetation model to update vegetation and terrestrial
! carbon storage.
!-----------------------------------------------------------------------
      IF (L_TRIFFID .AND.                                               &
     &   (ASTEPS_SINCE_TRIFFID == NSTEP_TRIF)) THEN

!-----------------------------------------------------------------------
! Calculate the TRIFFID inverse coupling timestep.
!-----------------------------------------------------------------------
        GAM_TRIF=360.0/FLOAT(TRIFFID_PERIOD)

!-----------------------------------------------------------------------
! Diagnose the mean fluxes over the coupling period.
!-----------------------------------------------------------------------
        DO L=1,LAND_PTS
          RESP_S_DR(L,1)=RESP_S_AC(L,1)*GAM_TRIF
          RESP_S_DR(L,2)=RESP_S_AC(L,2)*GAM_TRIF
          RESP_S_DR(L,3)=RESP_S_AC(L,3)*GAM_TRIF
          RESP_S_DR(L,4)=RESP_S_AC(L,4)*GAM_TRIF
        ENDDO

        DO N=1,NPFT
          DO J=1,TILE_PTS(N)
            L=TILE_INDEX(J,N)
            G_LEAF_DR(L,N)=G_LEAF_PHEN_AC(L,N)*GAM_TRIF
            NPP_DR(L,N)=NPP_AC(L,N)*GAM_TRIF
            RESP_W_DR(L,N)=RESP_W_AC(L,N)*GAM_TRIF
          ENDDO
        ENDDO

!-----------------------------------------------------------------------
! Diagnose the mean leaf turnover rates over the coupling period.
!-----------------------------------------------------------------------
        IF (L_PHENOL) THEN
          DO N=1,NPFT
            DO J=1,TILE_PTS(N)
              L=TILE_INDEX(J,N)
              G_LEAF_DR(L,N)=G_LEAF_PHEN_AC(L,N)*GAM_TRIF
            ENDDO
          ENDDO
        ELSE
          DO N=1,NPFT
            DO J=1,TILE_PTS(N)
              L=TILE_INDEX(J,N)
              G_LEAF_DR(L,N)=G_LEAF_AC(L,N)*GAM_TRIF
            ENDDO
          ENDDO
        ENDIF

!-----------------------------------------------------------------------
! Define the agricultural regions.
!-----------------------------------------------------------------------
        IF (AGRIC) THEN
        DO L=1,LAND_PTS
            FRAC_AGRIC(L)=FRACA(L)
        ENDDO
        ENDIF

!-----------------------------------------------------------------------
! Take copies of TRIFFID input variables for output as diagnostics.
!-----------------------------------------------------------------------
        DO N=1,NPFT
          DO L=1,LAND_PTS
            G_LEAF_DR_OUT(L,N)=G_LEAF_DR(L,N)
            NPP_DR_OUT(L,N)=NPP_DR(L,N)
            RESP_W_DR_OUT(L,N)=RESP_W_DR(L,N)
            FRAC_OLD(L,N)=FRAC(L,N)
          ENDDO
        ENDDO
        DO N=1,5
          DO L=1,LAND_PTS
            RESP_S_DR_OUT(L,N)=RESP_S_DR(L,N)
          ENDDO
        ENDDO
        DO L=1,LAND_PTS
          DCS(L)=CS(L,1)+CS(L,2)+CS(L,3)+CS(L,4)
!          RESP_S_DR_OUT(L,5) = RESP_S_DR_OUT(L,1) + RESP_S_DR_OUT(L,2)
!     &                       + RESP_S_DR_OUT(L,3) + RESP_S_DR_OUT(L,4)
        ENDDO

!-----------------------------------------------------------------------
! Select timestep and forward timestep weighting parameters for
! equilibrium or dynamic vegetation and call TRIFFID.
!-----------------------------------------------------------------------

! calculate RESP_FRAC.
! NOTE - here we want to use (1-X) - i.e. the fraction of
!        respiration NOT released to the atmos, so RESP_FRAC
!        here is 1-RESP_FRAC as calc'd in BL_CTL.

      DO L=1,LAND_PTS
        j=(land_index(l)-1)/row_length + 1
        i=land_index(l) - (j-1)*row_length
        clay_land(l) = clay_frac(i,j)
! NOTE!! NASTY HACK FOR NOW UNTIL ANCILLARY FILE IS SORTED OUT...
        IF (CLAY_LAND(L)  <   -1e9) CLAY_LAND(L)=0.23
        RESP_FRAC(L) = 1.0 / (4.0895 + 2.672*                           &
     &                   exp(-0.0786 * 100.0*clay_land(l)))
      ENDDO

        IF (L_TRIF_EQ) THEN
          FORW=1.0
          GAMMA=GAMMA_EQ
          KITER=ITER_EQ
        ELSE
          FORW=0.0
          GAMMA=GAM_TRIF
          KITER=1
        ENDIF

        DO K=1,KITER

          WRITE(6,*) 'Calling TRIFFID'

! DEPENDS ON: triffid
          CALL TRIFFID (LAND_PTS,TRIF_PTS,TRIF_INDEX,FORW,GAMMA         &
     &,                 FRAC_VS,FRAC_AGRIC,G_LEAF_DR,NPP_DR,RESP_S_DR   &
     &,                 RESP_W_DR,CS,FRAC,HT,LAI,RESP_FRAC              &
     &,                 C_VEG,CV,LIT_C,LIT_C_MN)

          WRITE(6,*) 'TRIFFID completed normally'

        ENDDO

!-----------------------------------------------------------------------
! Update TILE_INDEX for new surface type fractions.
!-----------------------------------------------------------------------
! DEPENDS ON: tilepts
        CALL TILEPTS(LAND_PTS,FRAC,TILE_PTS,TILE_INDEX)


!-----------------------------------------------------------------------
! Reset the accumulation fluxes to zero in equilibrium mode.
!-----------------------------------------------------------------------
        IF (L_TRIF_EQ) THEN

        DO N=1,4
          DO L=1,LAND_PTS
            RESP_S_AC(L,N)=0.0
          ENDDO
        ENDDO

        DO N=1,NPFT
          DO L=1,LAND_PTS
            NPP_AC(L,N)=0.0
              RESP_W_AC(L,N)=0.0
          ENDDO
        ENDDO

!-----------------------------------------------------------------------
! Reset the accumulation fluxes to the TRIFFID-diagnosed corrections
! in dynamic mode. Such corrections will be typically associated with
! the total depletion of a carbon reservoir or with the maintenance
! of the plant seed fraction. Any residual for the 4 component pools is
! assumed to be in proportion to the pool size itself.
!
! (in the case of zero total soil carbon, use 1e-10 to prevent zero
!  divide. This will have no impact because it is multiplied by the
!  constituent stores which are also zero...)
!-----------------------------------------------------------------------
        ELSE

          DO L=1,LAND_PTS
            CS_TOT(L)=MAX(1e-10,CS(L,1)+CS(L,2)+CS(L,3)+CS(L,4))
            DCS(L)=CS_TOT(L)-DCS(L)
            RESP_S_DR(L,1)=(1-RESP_FRAC(L)) * (RESP_S_DR_OUT(L,5))
            RESP_S_AC(L,1)=(LIT_C_MN(L)-GAMMA*DCS(L)-RESP_S_DR(L,1)) *  &
     &                       CS(L,1) / (GAMMA*CS_TOT(L))
            RESP_S_AC(L,2)=(LIT_C_MN(L)-GAMMA*DCS(L)-RESP_S_DR(L,1)) *  &
     &                       CS(L,2) / (GAMMA*CS_TOT(L))
            RESP_S_AC(L,3)=(LIT_C_MN(L)-GAMMA*DCS(L)-RESP_S_DR(L,1)) *  &
     &                       CS(L,3) / (GAMMA*CS_TOT(L))
            RESP_S_AC(L,4)=(LIT_C_MN(L)-GAMMA*DCS(L)-RESP_S_DR(L,1)) *  &
     &                       CS(L,4) / (GAMMA*CS_TOT(L))
          ENDDO

          DO N=1,NPFT
            DO J=1,TILE_PTS(N)
              L=TILE_INDEX(J,N)
              RATIO=FRAC_OLD(L,N)/FRAC(L,N)
              NPP_AC(L,N)=RATIO*(NPP_DR(L,N)-NPP_DR_OUT(L,N))/GAM_TRIF
              RESP_W_AC(L,N)=RATIO*(RESP_W_DR(L,N)-RESP_W_DR_OUT(L,N))  &
     &                             /GAM_TRIF
            ENDDO
          ENDDO

        ENDIF

!-----------------------------------------------------------------------
! Reset the accumulated leaf turnover rates to zero.
!-----------------------------------------------------------------------
        IF (L_PHENOL) THEN
          DO N=1,NPFT
            DO L=1,LAND_PTS
              G_LEAF_PHEN_AC(L,N)=0.0
            ENDDO
          ENDDO
        ELSE
          DO N=1,NPFT
            DO L=1,LAND_PTS
              G_LEAF_AC(L,N)=0.0
            ENDDO
          ENDDO
        ENDIF

        ASTEPS_SINCE_TRIFFID=0

      ENDIF

!-----------------------------------------------------------------------
! Calculate gridbox mean vegetation parameters from fractions of
! surface functional types
!-----------------------------------------------------------------------
! DEPENDS ON: sparm
      CALL SPARM (LAND_PTS,NTILES,CAN_MODEL,TILE_PTS,TILE_INDEX         &
     &,           FRAC,HT,LAI,SATCON,CATCH_S,CATCH_T,INFIL_T,Z0_T)

      RETURN
      END SUBROUTINE VEG
#endif
