#if defined(OASIS3)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************

MODULE oasis3_atm_data_mod
  ! Description:
  ! Useful data and arrays for use with OASIS coupling.
  !
  ! Author: R. Hill
  !
  ! Current Code Owner : R. Hill
  !
  !=====================================================================
  use mod_prism

  IMPLICIT NONE

  INTEGER(kind=ip_intwp_p) :: OASIS_COMP_ID !INTEGER :: OASIS_COMP_ID
  INTEGER :: OASIS_COMMUNICATOR
  INTEGER :: OASIS_GRID_ID


  INTEGER :: OASIS_PE  ! MYPE in OASIS 
  ! context operations
  INTEGER :: OASIS_CNTLPE
  INTEGER :: OASIS_NPROC
  INTEGER :: OASIS_FLD_TYPE

  ! The UM atmos doesn't have land sea masks available
  ! to us as a matter of course. So we have to go
  ! through the agony of calculating these ourselves
  ! on each grid point type. Here, we define the arrays
  ! which hold the masks.

  ! Mask arrays, local
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: OASIS_TMASK
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: OASIS_UMASK
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: OASIS_VMASK

  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: UM_TMASK
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: UM_UMASK
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: UM_VMASK


  ! Mask arrays, global.
  LOGICAL, DIMENSION(:,:), ALLOCATABLE :: GMASKU
  LOGICAL, DIMENSION(:,:), ALLOCATABLE :: GMASKV

  LOGICAL :: L_MASTER

  ! The current prism timestep.
  REAL :: PRISM_TIMESTEP
  INTEGER :: PRISM_NSEC  ! number of model seconds

  ! Indices used to point to various coupling fields.
  INTEGER :: vind_ocn_freezen(5)
  INTEGER :: vind_ocn_snowthickn(5)
  INTEGER :: vind_ocn_hicen(5)
  INTEGER :: vind_topmeltn(5)
  INTEGER :: vind_botmeltn(5)

!outgoing fields
  INTEGER , PARAMETER :: vind_heatflux = 1
  INTEGER , PARAMETER :: vind_taux = 2
  INTEGER , PARAMETER :: vind_tauy = 3
  INTEGER , PARAMETER :: vind_wme = 4
  INTEGER , PARAMETER :: vind_pen_solar = 5
  INTEGER , PARAMETER :: vind_pme = 6
  INTEGER , PARAMETER :: vind_runoff = 7
  INTEGER , PARAMETER :: vind_tsnow = 8
  INTEGER , PARAMETER :: vind_sublim = 9
  INTEGER , PARAMETER :: vind_topmelt = 10
  INTEGER , PARAMETER :: vind_botmelt = 11
  INTEGER , PARAMETER :: vind_co2 = 12
  INTEGER , PARAMETER :: vind_train = 13
  INTEGER , PARAMETER :: vind_lhflx = 14
  INTEGER , PARAMETER :: vind_evap2d = 15
  ! auscom coupling fields : outgoing
  INTEGER , PARAMETER :: vind_swflx = 16
  INTEGER , PARAMETER :: vind_lwflx = 17
  INTEGER , PARAMETER :: vind_shflx = 18
  INTEGER , PARAMETER :: vind_press = 19

  INTEGER , PARAMETER :: vind_wnd10 = 20

!incoming fields
  INTEGER , PARAMETER :: vind_ocn_sst = 100
  INTEGER , PARAMETER :: vind_ocn_u = 101
  INTEGER , PARAMETER :: vind_ocn_v = 102
  INTEGER , PARAMETER :: vind_ocn_freeze = 103
  ! auscom coupling fields : incoming
  INTEGER , PARAMETER :: vind_sal = 105

  INTEGER , PARAMETER :: vind_ocn_co2 = 106
  INTEGER , PARAMETER :: vind_ocn_co2fx = 107

  ! The use of DATA statements is not ideal in the following
  ! but we want to set up arrays and PARAMETER
  ! is no good for that.
  DATA vind_ocn_freezen    /131,132,133,134,135/
  DATA vind_ocn_snowthickn /136,137,138,139,140/
  DATA vind_ocn_hicen      /141,142,143,144,145/

  DATA vind_topmeltn       /31,32,33,34,35/
  DATA vind_botmeltn       /36,37,38,39,40/

  INTEGER  :: var_ind

  ! Arrays containing the "vind" numbers of the 
  ! transients we want to use. These are used to 
  ! define our transients using a simple loop.

  INTEGER , PARAMETER :: vind_max=200

  INTEGER :: vind_out(vind_max)
  CHARACTER (len=1) :: vind_out_type(vind_max)
  CHARACTER (len=8) :: vind_out_name(vind_max)
  INTEGER :: vind_in(vind_max)
  CHARACTER (len=1) :: vind_in_type(vind_max)
  CHARACTER (len=8) :: vind_in_name(vind_max)


  ! Dimensions for coupling transient arrays.
  ! These will simply hold copies of lasize, but make 
  ! code less long winded.
  INTEGER :: oasis_jmt
  INTEGER :: oasis_jmt_u
  INTEGER :: oasis_jmt_v
  INTEGER :: oasis_imt

  ! Outgoing coupling arrays
  REAL,DIMENSION(:,:),ALLOCATABLE :: taux
  REAL,DIMENSION(:,:),ALLOCATABLE :: tauy

  REAL,DIMENSION(:,:),ALLOCATABLE :: solar2d
  REAL,DIMENSION(:,:),ALLOCATABLE :: blue2d
  REAL,DIMENSION(:,:),ALLOCATABLE :: evap2d
  REAL,DIMENSION(:,:),ALLOCATABLE :: longwave2d
  REAL,DIMENSION(:,:),ALLOCATABLE :: sensible2d
  REAL,DIMENSION(:,:),ALLOCATABLE :: sublim
  REAL,DIMENSION(:,:),ALLOCATABLE :: heatflux
  REAL,DIMENSION(:,:),ALLOCATABLE :: latentflux

  REAL,DIMENSION(:,:),ALLOCATABLE :: rainls
  REAL,DIMENSION(:,:),ALLOCATABLE :: snowls
  REAL,DIMENSION(:,:),ALLOCATABLE :: rainconv
  REAL,DIMENSION(:,:),ALLOCATABLE :: snowconv
  REAL,DIMENSION(:,:),ALLOCATABLE :: totalrain
  REAL,DIMENSION(:,:),ALLOCATABLE :: totalsnow
  REAL,DIMENSION(:,:),ALLOCATABLE :: riverout
  REAL,DIMENSION(:,:),ALLOCATABLE :: wme
  REAL,DIMENSION(:,:),ALLOCATABLE :: um_co2
  REAL,DIMENSION(:,:),ALLOCATABLE :: um_wnd10

  REAL,DIMENSION(:,:,:),ALLOCATABLE :: topmeltn
  REAL,DIMENSION(:,:,:),ALLOCATABLE :: botmeltn

  ! Incoming coupling transients
  REAL,DIMENSION(:,:),ALLOCATABLE :: ocn_sst , ocn_sst_orig
  REAL,DIMENSION(:,:),ALLOCATABLE :: ocn_freeze
  REAL,DIMENSION(:,:,:),ALLOCATABLE :: ocn_freezen
  REAL,DIMENSION(:,:),ALLOCATABLE :: ocn_hice
  REAL,DIMENSION(:,:,:),ALLOCATABLE :: ocn_hicen
  REAL,DIMENSION(:,:),ALLOCATABLE :: ocn_snowthick
  REAL,DIMENSION(:,:,:),ALLOCATABLE :: ocn_snowthickn
  REAL,DIMENSION(:,:),ALLOCATABLE :: ocn_u
  REAL,DIMENSION(:,:),ALLOCATABLE :: ocn_v
  REAL,DIMENSION(:,:),ALLOCATABLE :: ocn_co2
  REAL,DIMENSION(:,:),ALLOCATABLE :: ocn_co2fx
  REAL,DIMENSION(:,:),ALLOCATABLE :: tstar_local
  REAL,DIMENSION(:,:),ALLOCATABLE :: tstar_ssi
  REAL,DIMENSION(:,:),ALLOCATABLE :: tstar_sice_local
  REAL,DIMENSION(:,:),ALLOCATABLE :: tstar_land_local

  ! Local land-fraction      
  REAL,DIMENSION(:,:),ALLOCATABLE :: fland_loc

END MODULE oasis3_atm_data_mod
#endif
