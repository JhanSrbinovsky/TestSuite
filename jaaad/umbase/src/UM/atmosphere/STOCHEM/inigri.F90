#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE INIGRI(eta_theta_levels,eta_rho_levels,total_mass)
! ---------------------------------------------------------------------
!-
!-   Purpose and Methods : SET UP GRIDS
!-
!-   Inputs  : None
!-   Outputs : LONG,LAT,LONGM,LONGM_HALF,LATM,LATM_HALF, and
!-             Stochem_Grid_Area
!-             are passed to module IN_STOCHEM_GRD
!-   Controls:
!-
!- History:
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  3.5    09/12/93  Created.  W.J. Collins
!  5.5    31/05/01  Altered for C grid staggering and passed outputs
!                   via module. C.E. Johnson
!  5.5    13/02/04  Calls HEIGHT_INI to initialise height arrays.
!                   K. Ketelsen.
!  6.1    20/10/04  Minor tidying of code. M.G. Sanderson
!  6.2    28/03/06  Minor changes for vn6.2  M.G. Sanderson
!
!-
!VVV  V5.0  INIGRI 30/V/01
! ---------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE HEIGHT_MOD                 !kk

      IMPLICIT NONE
! ---------------------------------------------------------------------
      ! eta on theta levels
      REAL, INTENT(IN) :: eta_theta_levels(0:model_levels)
      ! eta on rho levels
      REAL, INTENT(IN) :: eta_rho_levels(model_levels)
      ! mass of atmosphere up to model levels
      REAL, INTENT(IN) :: total_mass 

#include "typsize.h"
#include "typcona.h"

      INTEGER :: i ! loop counters
      INTEGER :: j !

      REAL    :: h

! Set up Eulerian grid boundaries:
!         Longitude: 0-360, W-E
!         Latitude:  0-180, S-N            ! Reverse direction in ND

! Met grid parameters - use v-grid.
      nmetlong= glsize(1,fld_type_u)
      nmetlat = glsize(2,fld_type_u)
      nmetlev = glsize(3,fld_type_u)

! The UM lat and long arrays are set in Inigri
      ALLOCATE(longm(1-halo_i:nmetlong+halo_i))       ! U
      ALLOCATE(longm_half(1-halo_i:nmetlong+halo_i))  ! V,W,T & Q
      ALLOCATE(latm(1-halo_j:nmetlat+halo_j))         ! V
      ALLOCATE(latm_half(1-halo_j:nmetlat+halo_j))    ! U,W,T & Q
      ALLOCATE(eta_theta(0:nmetlev))
      ALLOCATE(eta_rho(nmetlev))

      dlongm = delta_lambda / pi_over_180
      dlatm = delta_phi / pi_over_180

! STOCHEM grid
      long=(/(dlong*i,i=0,nlong-1)/)
      lat=(/(dlat*j,j=0,mnlat)/)

! Met (i.e. UM) grids
      longm_half = (/(dlongm*(i-1),i=1-halo_i,nmetlong+halo_i)/)
      longm = longm_half + dlongm/2.0

      latm_half = (/(dlatm*(j-1),j=1-halo_j,nmetlat+halo_j)/)
      latm = (/(dlatm/2.0+dlatm*(j-1),j=1-halo_j,nmetlat+halo_j)/)

      eta_theta = eta_theta_levels
      eta_rho = eta_rho_levels

!kk
      CALL HEIGHT_INI        ! Initialise Height for table lookup

      lmolec = (na*total_mass/mair) / ncell
! Calculate area of grid squares in each latitude band of STOCHEM
      h = 2.0*pi*(Earth_Radius**2) / REAL(nlong)

      stochem_grid_area = h * (SIN((90.0-lat(0:mnlat-1))*pi_over_180)   &
     &                      - SIN((90.0-lat(1:mnlat))*pi_over_180))

      END SUBROUTINE INIGRI
#endif
