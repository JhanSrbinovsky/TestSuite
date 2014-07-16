#if defined(A19_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Version 1A of vegetation section: models leaf phenology
!
! Subroutine Interface:
      SUBROUTINE VEG(LAND_PTS,LAND_INDEX,NTILES,CAN_MODEL               &
     &,              A_STEP,PHENOL_PERIOD,L_PHENOL                      &
     &,              ATIMESTEP,SATCON                                   &
     &,              G_LEAF_AC,FRAC,LAI,HT                              &
     &,              CATCH_S,CATCH_T,INFIL_T,Z0_T                       &
     &,              G_LEAF_DAY,G_LEAF_PHEN,LAI_PHEN                    &
     &               )


      IMPLICIT NONE
!
! Description:
!   Updates Leaf Area Index for Plant Functional Types (PFTs) and uses
!   this to derive new vegetation parameters for PFTs along with gridbox
!   mean values where appropriate.
!
! Method:
!   Calls PHENOL which models phenology and updates Leaf Area Index
!   (LAI), then passes new LAI into SPARM along with canopy height
!   and fractional cover of Plant Functional Types.  SPARM uses this to
!   derive the vegetation parameters for each PFT, and also derives
!   gridbox means where this is required.
!
! Current Code Owner:  Richard Betts
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   4.4    8/10/97   Original code.  Richard Betts
!   4.5   16/09/98   Call SWAPB_LAND to update halo regions of input
!                    fields.   Richard Betts
!   4.5   23/11/98   Output G_LEAF_DAY, G_LEAF_PHEN and LAI_PHEN as
!                    diagnostics.  Richard Betts
!   5.2   15/11/00   Re-Written for New Dynamics      M. Best
!   5.4   11/04/02   Canopy snow model for needleleaf tree tile.
!                    R. Essery
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.


      INTEGER                                                           &
     & LAND_PTS                                                         &
                             ! IN Number of land points to be processed.
     &,NTILES                                                           &
                             ! IN Number of land-surface tiles.
     &,CAN_MODEL                                                        &
                                    ! IN Swith for thermal vegetation
     &,A_STEP                                                           &
                             ! IN Atmospheric timestep number.
     &,PHENOL_PERIOD         ! IN Phenology period (days).

#include "nstypes.h"

      INTEGER                                                           &
     & LAND_INDEX(LAND_PTS)         ! IN Index of land points

      INTEGER                                                           &
     & I,J,L,N                      ! WORK loop counters.

      LOGICAL                                                           &
     & L_PHENOL                     ! IN .T. for interactive leaf
!                                   !    phenology.
      REAL                                                              &
     & ATIMESTEP                                                        &
                                    ! IN Atmospheric timestep (s).
     &,SATCON(LAND_PTS)                                                 &
                                    ! IN Saturated hydraulic
!                                   !    conductivity (kg/m2/s).
     &,G_LEAF_AC(LAND_PTS,NPFT)                                         &
                                    ! INOUT Accumulated leaf turnover
!                                   !       rate.
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
     &,LAI_PHEN(LAND_PTS,NPFT)                                          &
                                    ! OUT LAI of PFTs after phenology.
!                                   !     Required as separate variable
!                                   !     for top-level argument list
!                                   !     matching with VEG_IC2A.
     &,Z0_T(LAND_PTS,NTILES)        ! OUT Roughness length for tiles (m)

      INTEGER                                                           &
     & NSTEP_PHEN                                                       &
                                    ! WORK Number of atmospheric
!                                   !      timesteps between calls to
!                                   !      PHENOL.
     &,TILE_PTS(NTYPE)                                                  &
                                    ! WORK Number of land points which
!                                   !      include the nth surface type.
     &,TILE_INDEX(LAND_PTS,NTYPE)   ! WORK Indices of land points which
!                                   !      include the nth surface type.

      REAL                                                              &
     & DTIME_PHEN                                                       &
                                    ! WORK The phenology timestep (yr).
     &,G_LEAF_DAY(LAND_PTS,NPFT)                                        &
                                    ! WORK Mean leaf turnover rate for
!                                   !      input to PHENOL (/360days).
     &,G_LEAF_PHEN(LAND_PTS,NPFT)   ! WORK Mean leaf turnover rate over
!                                   !      phenology period (/360days).

!-----------------------------------------------------------------------
! Initialisations
!-----------------------------------------------------------------------
      DO N=1,NTILES
        DO L=1,LAND_PTS
          CATCH_S(L)=0.0
          CATCH_T(L,N)=0.0
          INFIL_T(L,N)=0.0
          Z0_T(L,N)=0.0
        ENDDO
      ENDDO

      DO N=1,NPFT
        DO L=1,LAND_PTS
          G_LEAF_PHEN(L,N)=0.0
          G_LEAF_DAY(L,N)=0.0
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Calculate the number of atmospheric timesteps between calls to PHENOL
! and TRIFFID.
!-----------------------------------------------------------------------
      NSTEP_PHEN=INT(86400.0*PHENOL_PERIOD/ATIMESTEP)


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
! state.
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
! Reset the accumulation over atmospheric model timesteps to zero.
!-----------------------------------------------------------------------
          DO L=1,LAND_PTS
            G_LEAF_AC(L,N)=0.0
          ENDDO
        ENDDO
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
