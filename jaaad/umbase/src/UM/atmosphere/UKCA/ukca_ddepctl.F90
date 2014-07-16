#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT*******************************
!
!  Purpose: To assign dry deposition rates in s-1 to array dpd in 
!           ASAD module asad_mod. This is passed for use in UKCA chemistry.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!           Called from UKCA_CHEMISTRY_CTRL.
!
! Current code owner: Olaf Morgenstern/Fiona O'Connor
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
      SUBROUTINE UKCA_DRYDEP_CTL(row_length, rows, bl_levels,         &
        land_points, land_index, tile_pts, tile_index, timestep,      &
        sinlat, tile_frac, t_surf, p_surf, dzl, zbl, surf_hf, u_s, q, &
        stcon, soilmc_lp, fland, seaice_frac, laift_lp, canhtft_lp,   &
        z0tile_lp, t0tile_lp, canwctile_lp,                           &
        nlev_in_bl, zdryrt)

      IMPLICIT NONE

#include "parvars.h"
#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "c_0_dg_c.h"

      INTEGER, INTENT(IN) :: row_length, rows
      INTEGER, INTENT(IN) :: bl_levels
      INTEGER, INTENT(IN) :: land_points
      INTEGER, INTENT(IN) :: land_index(land_points)

      INTEGER, DIMENSION(ntype), INTENT(IN)             :: tile_pts
      INTEGER, DIMENSION(land_points,ntype), INTENT(IN) :: tile_index

      REAL, INTENT(IN)                                  :: timestep

      REAL, DIMENSION(row_length,rows), INTENT(IN)      :: sinlat
      REAL, DIMENSION(land_points,ntype), INTENT(IN)    :: tile_frac
      REAL, DIMENSION(row_length,rows), INTENT(IN)      :: t_surf
      REAL, DIMENSION(row_length,rows), INTENT(IN)      :: p_surf
      REAL, DIMENSION(row_length,rows,bl_levels), INTENT(IN) :: dzl
      REAL, DIMENSION(row_length,rows), INTENT(IN)      :: zbl
      REAL, DIMENSION(row_length,rows), INTENT(IN)      :: surf_hf
      REAL, DIMENSION(row_length,rows), INTENT(IN)      :: u_s
      REAL, DIMENSION(row_length,rows), INTENT(IN)      :: q
      REAL, DIMENSION(row_length,rows), INTENT(IN)      :: seaice_frac
      REAL, DIMENSION(row_length,rows,npft), INTENT(IN) :: stcon
      REAL, DIMENSION(land_points), INTENT(IN)       :: soilmc_lp
      REAL, DIMENSION(land_points), INTENT(IN)       :: fland
      REAL, DIMENSION(land_points,npft), INTENT(IN)  :: laift_lp
      REAL, DIMENSION(land_points,npft), INTENT(IN)  :: canhtft_lp
      REAL, DIMENSION(land_points,ntype), INTENT(IN) :: z0tile_lp
      REAL, DIMENSION(land_points,ntype), INTENT(IN) :: t0tile_lp
      REAL, DIMENSION(land_points,ntype), INTENT(IN) :: canwctile_lp

      INTEGER, DIMENSION(row_length,rows), INTENT(OUT) :: nlev_in_bl

      REAL, DIMENSION(row_length,rows,jpdd), INTENT(OUT) :: zdryrt

!     Local variables

      INTEGER :: i, j, k, l, m, n
      INTEGER :: nwater = 7
      INTEGER, DIMENSION(ntype-1) :: land_tiles =                &
        (/ 1, 2, 3, 4, 5, 6, 8, 9 /)

      REAL :: seafrac ! Fraction of sea in grid square
      REAL :: lftotal ! Sum of land tile fractions

      REAL, DIMENSION(row_length,rows)       :: soilmc
      REAL, DIMENSION(row_length,rows)       :: land_fraction
      REAL, DIMENSION(row_length,rows,npft)  :: lai_ft
      REAL, DIMENSION(row_length,rows,npft)  :: canht_ft
      REAL, DIMENSION(row_length,rows,npft)  :: canwc_ft
      REAL, DIMENSION(row_length,rows,ntype) :: z0tile
      REAL, DIMENSION(row_length,rows,ntype) :: t0tile
      REAL, DIMENSION(row_length,rows)       :: o3_stom_frac

! Aerodynamic resistance (s m-1)
      REAL, DIMENSION(row_length,rows,ntype) :: ra

! Quasi-laminar resistance (s m-1)
      REAL, DIMENSION(row_length,rows,jpdd) :: rb

! Surface resistance (s m-1)
      REAL, DIMENSION(row_length,rows,ntype,jpdd) :: rc

! Global surface fraction array
      REAL, DIMENSION(row_length,rows,ntype) :: gsf
!
! Set up global surface fraction array, gsf
! Tile fractions from STASH (in array tile_frac) are only defined over land.
! Set water fractions over the sea. Add sea ice and sea fractions to gsf.
! Adjust tile fractions so they add up to 1. If gsf has the water tile
! fraction = 1.0, set the land fraction to 0.0
!
      land_fraction = 0.0
      soilmc        = 0.0
      lai_ft        = 0.0
      canht_ft      = 0.0
      canwc_ft      = 0.0
      z0tile        = 0.0
      t0tile        = 0.0
      gsf           = 0.0
!
! Expand arrays from land points to lon-lat grid
!
      DO l = 1, land_points
        j = (land_index(l)-1)/row_length + 1
        i = land_index(l) - (j-1)*row_length
        soilmc(i,j) = soilmc_lp(l)
        land_fraction(i,j) = fland(l)
      END DO
!
      DO m = 1, npft
        DO n = 1, tile_pts(m)
          l = tile_index(n,m)
          j = (land_index(l)-1)/row_length + 1
          i = land_index(l) - (j-1)*row_length
          lai_ft(i,j,m) = laift_lp(l,m)
          canht_ft(i,j,m) = canhtft_lp(l,m)
          canwc_ft(i,j,m) = canwctile_lp(l,m)
        END DO
      END DO
!
      DO m = 1, ntype
        DO n = 1, tile_pts(m)
          l = tile_index(n,m)
          j = (land_index(l)-1)/row_length + 1
          i = land_index(l) - (j-1)*row_length
          z0tile(i,j,m) = z0tile_lp(l,m)
          t0tile(i,j,m) = t0tile_lp(l,m)
          gsf(i,j,m) = tile_frac(l,m)
        END DO
      END DO
!
! Add in sea ice data and coastal land fractions to gsf array.
!
      DO k = 1, rows
        DO i = 1, row_length
          IF (gsf(i,k,nwater) == 1.0) land_fraction(i,k) = 0.0
          seafrac = 1.0 - land_fraction(i,k)
          IF (land_fraction(i,k) < 1.0 .AND. gsf(i,k,nwater) < 1.0) THEN
            gsf(i,k,nwater) = seafrac
            IF (land_fraction(i,k) > 0.0) THEN
              lftotal = 0.0
!CDIR EXPAND=8
              DO j = 1, ntype-1
                n = land_tiles(j)
                lftotal = lftotal + gsf(i,k,n)
              END DO
!CDIR EXPAND=8
              DO j = 1, ntype-1
                n = land_tiles(j)
                gsf(i,k,n) = gsf(i,k,n) * land_fraction(i,k) / lftotal
              END DO
            END IF
          END IF
          IF (seaice_frac(i,k) > 0.0) THEN
            gsf(i,k,nwater) = (1.0 - seaice_frac(i,k)) * seafrac
            gsf(i,k,ntype) = gsf(i,k,ntype) + seaice_frac(i,k) * seafrac
          END IF
        END DO
      END DO
!
! Set all tile temperatures to t0 where undefined.
!
      DO n = 1, ntype
        DO k = 1, rows
          DO i = 1, row_length
            IF (t0tile(i,k,n) < 100.0) t0tile(i,k,n) = t_surf(i,k)
          END DO
        END DO
      END DO
!
! Set up tile temperatures where a mixture of sea and sea ice is
! present. Set sea to freezing temperature (tfs) and ice to sea ice
! temperature.
!
      DO k = 1, rows
        DO i = 1, row_length
          IF (seaice_frac(i,k) > 0.0) THEN
!           t0tile(i,k,ntype) = sicetemp(i,k)
            t0tile(i,k,nwater) = tfs
          END IF
        END DO
      END DO
!
! Calculate aerodynamic and quasi-laminar resistances (Ra, Rb).
!
! DEPENDS ON: ukca_aerod
      CALL UKCA_AEROD(row_length, rows, nwater, t_surf, p_surf,  &
        surf_hf, u_s, canht_ft, gsf, z0tile, ra, rb)
!
! Calculate surface resistance term Rc
!
! DEPENDS ON: ukca_surf_dd_res
      CALL UKCA_SURF_DD_RES(row_length, rows, sinlat,            &
        t_surf, p_surf, q, soilmc, gsf, stcon, t0tile,           &
        lai_ft, canwc_ft, rc, o3_stom_frac)
!
! Combine resistance terms to calculate overall dry deposition
! velocity, and hence first-order rate constant describing
! dry deposition rate
!
! DEPENDS ON: ukca_ddcalc
      CALL UKCA_DDCALC(row_length, rows, bl_levels, timestep,    &
        dzl, zbl, gsf, ra, rb, rc,                               &
        nlev_in_bl, zdryrt)
!
      END SUBROUTINE UKCA_DRYDEP_CTL
#endif
