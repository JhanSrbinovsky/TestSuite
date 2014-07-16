#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      MODULE IN_STOCHEM_GRD
! Description:
! Specify physical constants and grid parameters
!
! Current Code Owner: C.E. Johnson
!
! History:
! Version  Date     Comment
! -------  ----     -------
!
!  5.0   24/05/01  New dynamics module for STOCHEM. C.E. Johnson
!  6.0   27/04/05  Fortran 90 format. M.G. Sanderson
!  6.0   17/09/03   Correct non-standard F90 continuation lines.
!                   Introduce standard UM modification history.
!                                                         P.Dando
!  6.1   20/08/04   Changes to allow vectorisation on SX6.
!                                                   M. Sanderson.
!  6.2   15/11/05  Extra variables added. M.G. Sanderson
! ----------------------------------------------------------------------
      IMPLICIT NONE
! ----------------------------------------------------------------------
#include "gccom.h"
#include "parvars.h"
#include "c_a.h"
#include "c_g.h"
#include "c_pi.h"
#include "c_rmol.h"
#include "cmaxsize.h"
#include "c_0_dg_c.h"
#include "nstypes.h"

! STOCHEM processor arrangement

      INTEGER, SAVE :: ncol      ! Number of longitude columns of pes
      INTEGER, SAVE :: numrows   ! Number of latitude rows of pes

! STOCHEM grid parameters
      INTEGER, PARAMETER :: nlong=96 ! No. of points in longitude
      INTEGER, PARAMETER :: mnlat=72 ! No. of points in latitude
      INTEGER, PARAMETER :: nlev=20  ! No. of levels
      INTEGER, SAVE      :: nlpe     ! No of latitude points per PE
      INTEGER, SAVE      :: nlnpe    ! No of longitude points per PE
      INTEGER, SAVE      :: max_stochem_points ! max no gridpts per PE

      INTEGER, SAVE      :: nclprc          ! Number of cells per pe
      INTEGER, PARAMETER :: mxcoll=65536    ! max size of arrays
                                            ! used in GCOM routines
      INTEGER            :: nstep=0         ! No. of steps
      INTEGER, PARAMETER :: ncell=100000    ! No. Lagrangian cells
      INTEGER, PARAMETER :: nhr=1           ! No. Hours in Met data

! Met grid parameters
      INTEGER, SAVE      :: nmetlong ! No. of longitude points
      INTEGER, SAVE      :: nmetlat  ! No. of latitude points
      INTEGER, SAVE      :: nmetlev  ! No. of levels
      INTEGER, PARAMETER :: nwater = 7 ! Water is tile 7.

! Number of Met. latitude points per PE:
      INTEGER, SAVE      :: nlatpe

! Number of Met. longitude points per PE:
      INTEGER, SAVE      :: nlonpe

! Number of Met. latitude points in use on this PE.
! Dimension nlatpe is set to maximum over all PEs, and hence arrays have
! an extra row on all PEs except PE 'nproc'.
      INTEGER, SAVE      :: rowspe

! PE responsible for grid sq
      INTEGER, DIMENSION(nlong,mnlat), SAVE :: procmap
      INTEGER, SAVE :: lnbound ! lowest Met longitude index on this pe
      INTEGER, SAVE :: lobound ! lowest Met latitude index on this pe
      INTEGER, SAVE :: lndat   ! lowest Stochem longitude index on pe
      INTEGER, SAVE :: ltdat   ! lowest Stochem latitude index on pe

! ----------------------------------------------------------------------
! Constants
      INTEGER, DIMENSION(12), SAVE :: daym                ! No. days
      REAL, SAVE                   :: daym_all            ! Sum of
      REAL, SAVE                   :: ysec                ! Seconds in
      REAL, PARAMETER              :: hoursec=3600.0      ! Seconds in
      REAL, PARAMETER              :: daysec=24.0*hoursec ! Seconds in

      REAL            :: rgc=rmol          ! Ideal gas constant
      REAL, PARAMETER :: na=6.022e23       ! Avogadro no. molecules/mol
      REAL, PARAMETER :: pstar=101325.     ! Std. Surface Press
      REAL, PARAMETER :: h_planck=6.6218E-34 ! Plancks constant
      REAL, PARAMETER :: c_light=2.997925E+08 ! speed of light
      REAL, PARAMETER :: dob=2.69e20       ! molecules O3/m2
      REAL, PARAMETER :: xx_co2=360.0e-6   ! [CO2]=360 ppmv
      REAL, PARAMETER :: rho_h2o=1.0e-3    ! density H2O kg/cm3
      REAL, PARAMETER :: etamed=0.07       ! Medium cloud boundary
      REAL, PARAMETER :: etahigh=0.18      ! High cloud boundary
      REAL, PARAMETER :: vdiff_bl=1.0e-01  ! Diffusion coeff in BL
      REAL, PARAMETER :: vdiff_trop=1.0E-02 ! Diffusion coeff in trop
      REAL, PARAMETER :: hdiff_bl=5300.0   ! Diffusion coeff in BL
      REAL, PARAMETER :: hdiff_trop=5300.0/4.0 ! Diffusion coeff in trop
      REAL, PARAMETER :: z_blplus=10.0     ! Extra b.l. mixing height m
      REAL, PARAMETER :: strat_tau=10*daysec ! Relax time for O3_conc

! Diffusion coefficients:
      REAL, PARAMETER :: mix_fact1=1.0e-2
      REAL, PARAMETER :: mix_fact2=1.0e-6

! Lightning profile parameters
      INTEGER, PARAMETER :: nlprofs=3
      INTEGER, PARAMETER :: nlkms=16

! STOCHEM time step parameters
      REAL, PARAMETER :: stochem_advection_step=1800.0
      REAL, PARAMETER :: stochem_chemistry_step=5.0

! STOCHEM
      REAL, PARAMETER :: dlat=180.0/mnlat  ! Latitude grid spacing
      REAL, PARAMETER :: dlong=360.0/nlong ! Longitude grid spacing

! UM
      REAL, SAVE      :: dlatm             ! Latitude grid spacing
      REAL, SAVE      :: dlongm            ! Longitude grid spacing

! Mean tropopause heights from -90 to 90 at 10 deg resolution
      REAL, PARAMETER, DIMENSION(20) :: mean_tropz = (/                 &
     &   7.50E03,    7.67E03,    7.78E03,    8.17E03,    8.89E03,       &
     &  10.56E03,   13.61E03,   15.39E03,   15.83E03,   16.11E03,       &
     &  15.83E03,   15.39E03,   13.61E03,   10.56E03,    8.89E03,       &
     &   8.17E03,    7.78E03,    7.67E03,    7.50E03,    7.50E03/)

! ----------------------------------------------------------------------
! UM Eta levels are set in INIGRI
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: eta_theta
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: eta_rho

! STOCHEM eta levels:
      REAL, DIMENSION(0:nlev) :: eta_stochem=(/                         &
     &         .000000  ,.010955  ,.022404  ,.034372                    &
     &        ,.046940  ,.060105  ,.073992  ,.088777                    &
     &        ,.104523  ,.121367  ,.139473  ,.159074                    &
     &        ,.180493  ,.204140  ,.231056  ,.261682                    &
     &        ,.297937  ,.344160  ,.410264  ,.522665                    &
     &        ,1.000000/)

! These constants are filled in Ini_Lev_Const
      REAL, DIMENSION(0:nlev), SAVE    :: zsea_stochem
      REAL, DIMENSION(0:nlev), SAVE    :: ck_stochem
      REAL, DIMENSION(nlev), SAVE      :: zsea_stochem_half
      REAL, DIMENSION(nlev), SAVE      :: ck_stochem_half

! The STOCHEM lat and long arrays are set in Inigri
      REAL, DIMENSION(nlong), SAVE     :: long        ! Stochem
      REAL, DIMENSION(0:mnlat), SAVE   :: lat         ! Stochem
      REAL, DIMENSION(mnlat), SAVE     :: stochem_grid_area

! The UM lat and long arrays are set in Inigri
      REAL, DIMENSION(:),ALLOCATABLE, SAVE  :: longm       ! U
      REAL, DIMENSION(:),ALLOCATABLE, SAVE  :: longm_half  ! V,W,T & Q
      REAL, DIMENSION(:),ALLOCATABLE, SAVE  :: latm        ! V
      REAL, DIMENSION(:),ALLOCATABLE, SAVE  :: latm_half   ! U,W,T & Q

! No. of molecules per Lagrangian cell:
      REAL, SAVE  :: lmolec

! Random number generation
      INTEGER, PARAMETER :: ransize = 128  ! FRANV needs 128 bytes


      CONTAINS

        SUBROUTINE SET_DAYM_ALL
          IMPLICIT NONE

          daym_all = SUM(daym)

          RETURN
        END SUBROUTINE SET_DAYM_ALL

      END MODULE IN_STOCHEM_GRD
#endif
