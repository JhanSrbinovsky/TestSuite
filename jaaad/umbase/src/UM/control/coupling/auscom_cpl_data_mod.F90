#if defined(OASIS3) && defined(ACCESS)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************

MODULE auscom_cpl_data_mod
  ! Description:
  ! Useful data and arrays for use with AusCOM/OASIS coupling.
  !
  !=====================================================================

  IMPLICIT NONE

  REAL(KIND=8) :: access_tfs         ! sea water freezing temperature (Celsius)
  !dhb599 20110601: fra298's new cloud fraction increment scaling
  real(kind=8) :: xfactor = 1.0            
  !dhb599 20110812: define Solar Constant here and make it read in from namelist
  real(kind=8) :: SC=1365.0   !0.0
  !*set default SC=0 to quickly verify if this namelist parameter takes effect!
  !dhb599 20110815: define CO2 level here
  real(kind=8) :: co2_init=370.0     !PD level
  real(kind=8) :: VOLCTS_val=97.33   !average volcanic forcing for control run

  LOGICAL(KIND=4) :: ocn_sss         ! switch to receive salinity field from ice
  LOGICAL(KIND=4) :: sdump_enable    ! outgoing cpl fields netcdf dump toggle
  LOGICAL(KIND=4) :: rdump_enable    ! incoming cpl fields netcdf dump toggle

  namelist/coupling/access_tfs,ocn_sss,sdump_enable,rdump_enable,xfactor,SC,co2_init,VOLCTS_val

!moved to oasis3_atm_data_mod.F90 to avoid duplicated vind when adding new fields
!  ! auscom coupling fields : outgoing
!  INTEGER , PARAMETER :: vind_swflx = 16
!  INTEGER , PARAMETER :: vind_lwflx = 17
!  INTEGER , PARAMETER :: vind_shflx = 18
!  INTEGER , PARAMETER :: vind_press = 19
!
!  ! auscom coupling fields : incoming
!  INTEGER , PARAMETER :: vind_sal = 25

  ! auscom coupling data buffers : outgoing
  REAL,DIMENSION(:,:),ALLOCATABLE :: auscom_swflx
  REAL,DIMENSION(:,:),ALLOCATABLE :: auscom_lwflx
  REAL,DIMENSION(:,:),ALLOCATABLE :: auscom_shflx
  REAL,DIMENSION(:,:),ALLOCATABLE :: auscom_pressure

  ! auscom coupling data buffers: incoming
  REAL,DIMENSION(:,:),ALLOCATABLE :: auscom_salinity

  ! send dump variable ids for coupling fields (outgoing)
  integer(kind=8) :: dump_hflx, dump_solflx, dump_runoff, dump_wme
  integer(kind=8) :: dump_train, dump_tsnow, dump_evap, dump_lhflx
  integer(kind=8) :: dump_top(5), dump_bot(5), dump_taux, dump_tauy
  integer(kind=8) :: dump_co2, dump_wnd10

  ! additional fields for auscom coupling (ACCESS)
  integer(kind=8) :: dump_swflx, dump_lwflx, dump_shflx, dump_press

  ! recv dump variable ids for coupling fields (incoming)
  integer(kind=8) :: dump_sst, dump_frzn(5), dump_snwtn(5)
  integer(kind=8) :: dump_hicn(5), dump_suno, dump_svno
  integer(kind=8) :: dump_ocn_co2, dump_ocn_co2fx

  ! additional fields for auscom coupling (ACCESS)
  integer(kind=8) :: dump_sss

  ! override l_oasis
  logical(kind=8) :: l_auscom

END MODULE auscom_cpl_data_mod
#endif
