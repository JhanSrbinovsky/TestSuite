#if defined(OASIS3)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
MODULE OASIS3_ATMOS_INIT
  !
  ! Description:    
  !     Define modules from the OASIS3 system
  !     containing variables and routines which we want to use.
  !
  ! Author: R. Hill
  !
  ! Current Code Owner : R. Hill
  !
  !=====================================================================

  USE mod_prism

  USE oasis3_atm_data_mod
#if defined(ACCESS)
  USE auscom_cpl_data_mod
  USE dump_sent
  USE dump_received
#endif

  IMPLICIT NONE

#include "c_kinds.h"

  ! Allocate space for three different grid types
  ! (T, U, V)
  INTEGER (kind=integer32), PARAMETER :: n_grids=6
  INTEGER (kind=integer32)  :: grid_id(n_grids)
  INTEGER (kind=integer32), PARAMETER :: t_grid = 1
  INTEGER (kind=integer32), PARAMETER :: u_grid = 2
  INTEGER (kind=integer32), PARAMETER :: v_grid = 3

  INTEGER (kind=integer32)  :: n_invar
  INTEGER (kind=integer32)  :: n_outvar
  INTEGER (kind=integer32)  :: um_start_date(6)

  INTEGER (kind=integer32), PARAMETER :: n_tranvar=vind_max !150 ! max number of
  ! trans vars

  INTEGER (kind=integer32)  :: var_id(n_tranvar)  ! In/outgoing
  ! transients
  LOGICAL :: no_neighbour_N

CONTAINS
  SUBROUTINE OASIS3_grid32(MODEL_BASIS_TIME,ierr64    &
     &, valid_sh64,actual_sh64,nGridDims64                &
     &,glsize_x,datastart,Ndim_max,MYPE)


    IMPLICIT NONE
#include <netcdf.inc>

    INTEGER (kind=integer64) :: ierr64
    INTEGER (kind=integer64) :: MYPE
    INTEGER (kind=integer32) :: ierror
    INTEGER (kind=integer32) :: info

    INTEGER (kind=integer64) :: MODEL_BASIS_TIME(6)

    INTEGER (kind=integer64) :: Ndim_max
    INTEGER (kind=integer64) :: datastart(Ndim_max)
    INTEGER (kind=integer64) :: glsize_x

    ! This module expects to define an Arakawa C grid:
    !   - this is an instance of PRISM_reglonlatvrt
    !   - 3D grid
    !   - points are created for the first time here
    !   - each grid cell is a cube w/ 8 vertices
    !   - C grid has t, u, v points

    INTEGER (kind=integer64) :: nGridDims64
    INTEGER (kind=integer64) :: t_points_64,u_points_64,v_points_64
    INTEGER (kind=integer32) :: t_points_id,u_points_id,v_points_id

    INTEGER (kind=integer64) :: valid_sh64(2,nGridDims64)
    INTEGER (kind=integer32) :: valid_shape(2,nGridDims64)
    INTEGER (kind=integer64) :: actual_sh64(2,nGridDims64)
    INTEGER (kind=integer32) :: actual_shape(2,nGridDims64)

    ! V grid has 1 fewer row
    INTEGER (kind=integer32) :: valid_shape_v(2,nGridDims64)
    INTEGER (kind=integer32) :: actual_shape_v(2,nGridDims64)

    INTEGER (kind=integer32) :: var_nodims(2)
    INTEGER (kind=integer32) :: var_shape(4)

    ! local variables...
    INTEGER (kind=integer32) :: i,j,k,ip1, istatus
    INTEGER (kind=integer32) :: max_i

    ! Define arrays for t, u and v masks. In OASIS3 terms, 1
    ! indicates land and 0 indicates sea.

    REAL (kind=real64), DIMENSION(actual_sh64(2,1),actual_sh64(2,2))  &
       & :: rla_lon, rla_lat, r_area

    REAL (kind=real64), DIMENSION(actual_sh64(2,1),actual_sh64(2,2),4)&
       & :: rla_clon, rla_clat

    REAL (kind=real64), DIMENSION(actual_sh64(2,1),actual_sh64(2,2)-1)&
       & :: rla_lon_v, rla_lat_v, r_area_v

    REAL (kind=real64),                     &
       & DIMENSION(actual_sh64(2,1),actual_sh64(2,2)-1,4)  &
       & :: rla_clon_v, rla_clat_v

    REAL (kind=real64), DIMENSION(actual_sh64(2,1),actual_sh64(2,2))  &
       & :: t_lon, t_lat, u_lon, u_lat, t_area, u_area

    INTEGER (kind=integer32), &
       &         DIMENSION(actual_sh64(2,1),actual_sh64(2,2)) ::    &
       &         t_mask, u_mask

    REAL (kind=real64), DIMENSION(actual_sh64(2,1),actual_sh64(2,2),4)&
       & :: tc_lon, tc_lat, uc_lon, uc_lat

    REAL (kind=real64), &
       & DIMENSION(actual_sh64(2,1),actual_sh64(2,2)-1)  &
       & :: v_lon, v_lat, v_area

    INTEGER (kind=integer32), &
       & DIMENSION(actual_sh64(2,1),actual_sh64(2,2)-1)  &
       & :: v_mask

    REAL (kind=real64), &
       & DIMENSION(actual_sh64(2,1),actual_sh64(2,2)-1,4)  &
       & :: vc_lon, vc_lat

    INTEGER :: nc_flag  ! Flag used in OASIS3 prism calls
    INTEGER :: nc_file_id
    INTEGER :: nc_var_id
    INTEGER :: grid_file_flag
    INTEGER :: ist(4), icnt(4)

    INTEGER (kind=integer32), DIMENSION(:), ALLOCATABLE :: il_paral 
    ! Decomposition for each proc
    INTEGER :: ig_nsegments  ! Number of segments of process
    !decomposition
    INTEGER :: ig_parsize    ! Size of array decomposition
    INTEGER :: id_nbcplproc  ! Number of processes involved in the
    !coupling
    INTEGER (kind=integer32):: partition_id    ! Local partition ID
    INTEGER (kind=integer32):: partition_id_v  ! Local partition ID 
    ! for v points
    INTEGER (kind=integer32):: partition
    INTEGER :: id_length     ! Size of partial field for each process
    INTEGER :: id_rank       ! Rank of process
    INTEGER :: ld_mparout    ! Unit of log file
    LOGICAL :: ld_comparal
    INTEGER :: ib

    REAL (kind=real32), &
       & DIMENSION(actual_sh64(2,1),actual_sh64(2,2)) &
       & :: data_2d

    REAL (kind=real32), &
       & DIMENSION(actual_sh64(2,1),actual_sh64(2,2),1) &
       & :: data_3d



    REAL (kind=real64), &
       & DIMENSION(actual_sh64(2,1),actual_sh64(2,2)) &
       & :: data_2d64

    REAL (kind=real64), &
       & DIMENSION(actual_sh64(2,1),actual_sh64(2,2),1) &
       & :: data_3d64

    um_start_date = MODEL_BASIS_TIME

    ! Make local 32 bit copies of valid and actual shape arrays
    ! which can then be used directly in PRISM OASIS3 calls copies.
    valid_shape(:,:) = valid_sh64(:,:)
    actual_shape(:,:) = actual_sh64(:,:)

    ! Set valid shape for v grid
    valid_shape_v(:,:) = valid_shape(:,:)
    actual_shape_v(:,:) = actual_shape(:,:)

    ! If this PE is at the northern domain limit
    ! subtract 1 row from the v grid dimension.
    ! What this effectively means is that we're saying
    ! if a PE owns a T point, then it must also
    ! own the v point to the north of the t point
    ! (and also the u point to the east).

    IF (L_MASTER) THEN
      ! If we're coupling through a master pe
      ! then the master pe needs to reduce the size
      ! of v grid fields by 1 row. In fact
      ! this will be perfomed on all PEs
      ! so that all PEs area aware of the
      ! the global size required which could be useful
      ! for later operations
      valid_shape_v(2,2) = valid_shape_v(2,2)-1
      actual_shape_v(2,2) = actual_shape_v(2,2)-1
    ELSE
      ! If we're coupling through each PE in parallel
      ! only the northern-most PE(s) trim a row
      ! from V grid fields
      IF (no_neighbour_N) THEN
        valid_shape_v(2,2) = valid_shape_v(2,2)-1
        actual_shape_v(2,2) = actual_shape_v(2,2)-1
      END IF
    END IF

    IF (MYPE.EQ.OASIS_CNTLPE) THEN
      ! If I'm the master PE then check whether the OASIS
      ! grid-definition netcdf files are present in the run directory
      ! and if not create them.

      !Oasis3-mct has not supported prism_start_grids_writing() which always returns TRUE
      !CALL prism_start_grids_writing(grid_file_flag)
      grid_file_flag = 0

      IF (grid_file_flag .EQ. 1) THEN

        istatus=NF_OPEN('grids_in.nc', NF_NOWRITE, &
           &                       nc_file_id)
        IF (istatus/=0) WRITE(6,*) "Open grids file status=",istatus

        ist(1:4) = 1
        icnt(1)  = valid_shape(2,1)
        icnt(2)  = valid_shape(2,2)
        icnt(3)  = 4
        icnt(4)  = 1

        istatus=NF_INQ_VARID(nc_file_id, 'atm3.lon' , nc_var_id)
        istatus=NF_GET_VARA_DOUBLE (nc_file_id, nc_var_id, ist(1:2),&
           icnt(1:2), t_lon(:,:))

        istatus=NF_INQ_VARID(nc_file_id, 'atm3.lat' , nc_var_id)
        istatus=NF_GET_VARA_DOUBLE (nc_file_id, nc_var_id, ist(1:2),&
           icnt(1:2), t_lat(:,:))

        ! Get corners for T points
        istatus=NF_INQ_VARID(nc_file_id, 'atm3.clo' , nc_var_id)
        istatus=NF_GET_VARA_DOUBLE (nc_file_id, nc_var_id, ist(1:3),&
           icnt(1:3), tc_lon(:,:,:))

        istatus=NF_INQ_VARID(nc_file_id, 'atm3.cla' , nc_var_id)
        istatus=NF_GET_VARA_DOUBLE (nc_file_id, nc_var_id, ist(1:3),&
           icnt(1:3), tc_lat(:,:,:))


        ! Write everything we've just read, using the OASIS3 calls
        rla_lon(:,:) = t_lon(:,:)
        rla_lat(:,:) = t_lat(:,:)
        rla_clon(:,:,:) = tc_lon(:,:,:)
        rla_clat(:,:,:) = tc_lat(:,:,:)
        CALL prism_write_grid('atm3', icnt(1), icnt(2), rla_lon, &
           rla_lat)
        CALL prism_write_corner('atm3', icnt(1), icnt(2), 4, &
           rla_clon, rla_clat)



        istatus=NF_INQ_VARID(nc_file_id, 'aum3.lon' , nc_var_id)
        istatus=NF_GET_VARA_DOUBLE (nc_file_id, nc_var_id, ist(1:2),&
           icnt(1:2), u_lon(:,:))

        istatus=NF_INQ_VARID(nc_file_id, 'aum3.lat' , nc_var_id)
        istatus=NF_GET_VARA_DOUBLE (nc_file_id, nc_var_id, ist(1:2),&
           icnt(1:2), u_lat(:,:))

        ! Get corners for U points
        istatus=NF_INQ_VARID(nc_file_id, 'aum3.clo' , nc_var_id)
        istatus=NF_GET_VARA_DOUBLE (nc_file_id, nc_var_id, ist(1:3),&
           icnt(1:3), uc_lon(:,:,:))

        istatus=NF_INQ_VARID(nc_file_id, 'aum3.cla' , nc_var_id)
        istatus=NF_GET_VARA_DOUBLE (nc_file_id, nc_var_id, ist(1:3),&
           icnt(1:3), uc_lat(:,:,:))

        rla_lon(:,:) = u_lon(:,:)
        rla_lat(:,:) = u_lat(:,:)
        rla_clon(:,:,:) = uc_lon(:,:,:)
        rla_clat(:,:,:) = uc_lat(:,:,:)
        CALL prism_write_grid('aum3', icnt(1), icnt(2), rla_lon, &
           rla_lat)
        CALL prism_write_corner('aum3', icnt(1), icnt(2), 4, &
           rla_clon, rla_clat)


        ! Adjust our Y dimension to cater for the fact that
        ! there is 1 row fewer for V points.
        icnt(2) = valid_shape_v(2,2)

        istatus=NF_INQ_VARID(nc_file_id, 'avm3.lon' , nc_var_id)
        istatus=NF_GET_VARA_DOUBLE (nc_file_id, nc_var_id, ist(1:2),&
           icnt(1:2), v_lon(:,:))

        istatus=NF_INQ_VARID(nc_file_id, 'avm3.lat' , nc_var_id)
        istatus=NF_GET_VARA_DOUBLE (nc_file_id, nc_var_id, ist(1:2),&
           icnt(1:2), v_lat(:,:))

        ! Get corners for V points
        istatus=NF_INQ_VARID(nc_file_id, 'avm3.clo' , nc_var_id)
        istatus=NF_GET_VARA_DOUBLE (nc_file_id, nc_var_id, ist(1:3),&
           icnt(1:3), vc_lon(:,:,:))

        istatus=NF_INQ_VARID(nc_file_id, 'avm3.cla' , nc_var_id)
        istatus=NF_GET_VARA_DOUBLE (nc_file_id, nc_var_id, ist(1:3),&
           icnt(1:3), vc_lat(:,:,:))

        rla_lon_v(:,:) = v_lon(:,:)
        rla_lat_v(:,:) = v_lat(:,:)
        rla_clon_v(:,:,:) = vc_lon(:,:,:)
        rla_clat_v(:,:,:) = vc_lat(:,:,:)
        CALL prism_write_grid('avm3', icnt(1), icnt(2), rla_lon_v, &
           rla_lat_v)
        CALL prism_write_corner('avm3', icnt(1), icnt(2), 4, &
           rla_clon_v, rla_clat_v)

        ! Close our input NC file.
        istatus=NF_CLOSE(nc_file_id)

        ! Now define the land-sea masks for the T, U and V points
        istatus=NF_OPEN('masks_in.nc', NF_NOWRITE, &
           nc_file_id)

        ist(:) = 1
        icnt(1)  = valid_shape(2,1)
        icnt(2)  = valid_shape(2,2)
        icnt(3:4) = 1


        istatus=NF_INQ_VARID(nc_file_id, 'atm3.msk' , nc_var_id)
        istatus=NF_GET_VARA_INT (nc_file_id, nc_var_id, ist(1:2), &
           icnt(1:2), &
           t_mask(:,:))

        ! Define L-S mask on T points for OASIS3
        CALL prism_write_mask('atm3', icnt(1), icnt(2),  t_mask)


        istatus=NF_INQ_VARID(nc_file_id, 'aum3.msk' , nc_var_id)
        istatus=NF_GET_VARA_INT (nc_file_id, nc_var_id, ist(1:2), &
           icnt(1:2), &
           u_mask(:,:))

        ! Define L-S mask on U points for OASIS3
        CALL prism_write_mask('aum3', icnt(1), icnt(2),  u_mask)

        ! Get the V mask (1 row less)
        icnt(2) = valid_shape_v(2,2)
        istatus=NF_INQ_VARID(nc_file_id, 'avm3.msk' , nc_var_id)
        istatus=NF_GET_VARA_INT (nc_file_id, nc_var_id, ist(1:2), &
           icnt(1:2), &
           v_mask(:,:))

        ! Define L-S mask on V points for OASIS3
        CALL prism_write_mask('avm3', icnt(1), icnt(2),v_mask)

        istatus=NF_CLOSE(nc_file_id)

        ! Read and define the areas for each set of grid points.
        ! Now define the areas for the T, U and V points
        istatus=NF_OPEN('areas_in.nc', NF_NOWRITE, &
           nc_file_id)

        ist(:) = 1
        icnt(1)  = valid_shape(2,1)
        icnt(2)  = valid_shape(2,2)
        icnt(3:4) = 1

        ! Read in and write T box areas
        istatus=NF_INQ_VARID(nc_file_id, 'atm3.srf' , nc_var_id)
        istatus=NF_GET_VARA_DOUBLE (nc_file_id, nc_var_id, ist(1:2),&
           icnt(1:2), t_area(:,:))

        r_area(:,:) = t_area(:,:)
        CALL prism_write_area('atm3', icnt(1), icnt(2), r_area)


        ! Read in and write U box areas
        istatus=NF_INQ_VARID(nc_file_id, 'aum3.srf' , nc_var_id)
        istatus=NF_GET_VARA_DOUBLE (nc_file_id, nc_var_id, ist(1:2),&
           icnt(1:2), u_area(:,:))

        r_area(:,:) = u_area(:,:)
        CALL prism_write_area('aum3', icnt(1), icnt(2), r_area)


        ! Read in and write V box areas - adjust dimensions
        icnt(2)  = valid_shape_v(2,2)
        istatus=NF_INQ_VARID(nc_file_id, 'avm3.srf' , nc_var_id)
        istatus=NF_GET_VARA_DOUBLE (nc_file_id, nc_var_id, ist(1:2),&
           icnt(1:2), v_area(:,:))

        r_area_v(:,:) = v_area(:,:)
        CALL prism_write_area('avm3', icnt(1), icnt(2), r_area_v)


        istatus=NF_CLOSE(nc_file_id)


      END IF ! IF grids.nc file needed writing

      ! End the grids writing
      ! This call is absolutely critical. Without it the code
      ! will deadlock in prism_enddef_proto and give you
      ! a traceback to an MPI_WAIT and no clues about
      ! what the real problem is.
      CALL prism_terminate_grids_writing()


    END IF

    ! So now we should have written all necessary grid related data
    ! to the standard OASIS3 NC files if they did not already exist.

    ierror = 0  ! Initialise error flag to zero because some OASIS
    ! routines do not do this explicitily which
    ! can cause us problems. 

    IF ((MYPE.EQ.OASIS_CNTLPE).OR.(.NOT.L_MASTER)) THEN
      ! If we're not coupling through the master, then each PE
      ! couples in parallel, hence they all need to call this.
      ! Otherwise only the master calls this code.

      ig_nsegments = 1
      ig_parsize = 5
      ALLOCATE(il_paral(ig_parsize))

      IF (L_MASTER) THEN

#if defined(ACCESS)
        ! Oasis3 handling of single-processor 2D partition
        ! is very inefficient so we use serial instead.
        il_paral(1) = 0
        il_paral(2) = 0
        il_paral(3) = valid_shape(2,1)*valid_shape(2,2)

#else
        ! Set up partitions in TWO dimensions NOT 1 as in the
        ! standard OASIS3 toy models which makes visualisation
        ! awkward.
        ! We couple through master - no decomposition necessary
        il_paral(1) = 2
        il_paral(2) = 0
        il_paral(3) = valid_shape(2,1)
        il_paral(4) = valid_shape(2,2)
        il_paral(5) = valid_shape(2,1)
#endif

        ! We have to define two separate partitions -
        ! The first one for the T and U points, the second one for
        ! the V points which has one row fewer.
        CALL prism_def_partition_proto (partition_id, il_paral, &
           ierror)

#if defined(ACCESS)
        il_paral(3) = valid_shape_v(2,1)*valid_shape_v(2,2)
#else
        il_paral(4) = valid_shape_v(2,2)
#endif

        CALL prism_def_partition_proto (partition_id_v, il_paral, &
           ierror)

      ELSE
        ! Set up partitions in TWO dimensions NOT 1 as in the
        ! standard OASIS3 toy models which makes visualisation
        ! awkward.

        ! See OASIS user guide for the meanings of the following
        ! array components and possible alternatives.
        ! 1 = partition type (we use box)
        ! 2 = data offset
        ! 3 = local x dimension
        ! 4 = local y dimension
        ! 5 = global x dimension

        il_paral(1) = 2
        il_paral(2) = (datastart(2)-1)*glsize_x+datastart(1)-1 
        il_paral(3) = valid_shape(2,1)
        il_paral(4) = valid_shape(2,2)
        il_paral(5) = glsize_x         

        ! We have to define two separate partitions -
        ! The first one for the T and U points, the second one for
        ! the V points which has one row fewer.
        CALL prism_def_partition_proto (partition_id, il_paral, &
           ierror)

        il_paral(4) = valid_shape_v(2,2)

        CALL prism_def_partition_proto (partition_id_v, il_paral, &
           ierror)
      END IF

      DEALLOCATE(il_paral)
      !=================================================================
      ! NOTE: OASIS3 will issue warnings about transients with
      ! names longer  than 8 characters and actually ignores them
      ! because it seems to truncate everything to 8 characters
      ! which means the names then don't correspond to anything
      ! in the namcouple file! 
      ! If your name here doesn't correspond with the name in namcouple
      ! it will not necessarily tell you, or crash. Quite often it will run
      ! apparently OK even to the extent that a receiving component will
      ! get data (apparently zeros) from somewhere (who knows where?)
      ! and proceed OK. That's a nightmare to spot so be VERY
      ! careful about spelling transient field names.  
      ! Here we define the transient variables which are actually required
      ! in an attempt to make the code slightly more flexible.
      !=================================================================
      vind_out(1) = vind_heatflux
      vind_out_name(1) = "heatflux"

      vind_out(2) = vind_pen_solar
      vind_out_name(2) = "pen_sol"

      vind_out(3) = vind_runoff
      vind_out_name(3) = "runoff"

      vind_out(4) = vind_wme
      vind_out_name(4) = "wme"

      vind_out(5) = vind_train
      vind_out_name(5) = "train"

      vind_out(6) = vind_tsnow
      vind_out_name(6) = "tsnow"

      vind_out(7) = vind_evap2d
      vind_out_name(7) = "evap2d"

      vind_out(8) = vind_lhflx
      vind_out_name(8) = "lhflx"

      vind_out(9:13)  = vind_topmeltn(1:5)
      vind_out_name(9) = "tmlt01"
      vind_out_name(10) = "tmlt02"
      vind_out_name(11) = "tmlt03"
      vind_out_name(12) = "tmlt04"
      vind_out_name(13) = "tmlt05"

      vind_out(14:18) = vind_botmeltn(1:5)
      vind_out_name(14) = "bmlt01"
      vind_out_name(15) = "bmlt02"
      vind_out_name(16) = "bmlt03"
      vind_out_name(17) = "bmlt04"
      vind_out_name(18) = "bmlt05"

      vind_out(19) = vind_taux
      vind_out_name(19) = "taux"

      vind_out(20) = vind_tauy
      vind_out_name(20) = "tauy"

      vind_out_type(1:18) = "T"
      vind_out_type(19) = "U"
      vind_out_type(20) = "V"

#if defined(ACCESS)
! gol124: auscom coupling
      vind_out(21) = vind_swflx
      vind_out_name(21) = "swflx"

      vind_out(22) = vind_lwflx
      vind_out_name(22) = "lwflx"

      vind_out(23) = vind_shflx
      vind_out_name(23) = "shflx"

      vind_out(24) = vind_press
      vind_out_name(24) = "press"

      vind_out_type(21:24) = "T"

      vind_out(25) = vind_co2
      vind_out_name(25) = "co2_a"
      vind_out_type(25) = "T"
      vind_out(26) = vind_wnd10
      vind_out_name(26) = "wnd_a"
      vind_out_type(26) = "V"

#endif

      var_nodims(1) = 2
      var_nodims(2) = 1

      var_id(:) = 0 ! Initialise - OASIS should assign these

      var_shape(1) = 1
      var_shape(2) = valid_shape(2,1)
      var_shape(3) = 1

#if defined(ACCESS)
      DO I = 1, 26
#else
      DO I = 1, 20
#endif

        IF (vind_out_type(I) == "V") THEN
          var_shape(4) = valid_shape_v(2,2)
          partition = partition_id_v
        ELSE
          var_shape(4) = valid_shape(2,2)
          partition = partition_id
        END IF

        ! Define output field.
        CALL prism_def_var_proto(var_id(vind_out(I)),         &
           &        vind_out_name(I), partition,                        &
           &        var_nodims, PRISM_Out,  var_shape,                  &
           &        PRISM_Real, ierror )

      ENDDO

#if defined(ACCESS)
      IF (sdump_enable) THEN
          CALL add_svar (vind_out_name(1),  dump_hflx, .false.)
          CALL add_svar (vind_out_name(2),  dump_solflx, .false.)
          CALL add_svar (vind_out_name(3),  dump_runoff, .false.)
          CALL add_svar (vind_out_name(4),  dump_wme, .false.)
          CALL add_svar (vind_out_name(5),  dump_train, .false.)
          CALL add_svar (vind_out_name(6),  dump_tsnow, .false.)
          CALL add_svar (vind_out_name(7),  dump_evap, .false.)
          CALL add_svar (vind_out_name(8),  dump_lhflx, .false.)
          CALL add_svar (vind_out_name(9),  dump_top(1), .false.)
          CALL add_svar (vind_out_name(10), dump_top(2), .false.)
          CALL add_svar (vind_out_name(11), dump_top(3), .false.)
          CALL add_svar (vind_out_name(12), dump_top(4), .false.)
          CALL add_svar (vind_out_name(13), dump_top(5), .false.)
          CALL add_svar (vind_out_name(14), dump_bot(1), .false.)
          CALL add_svar (vind_out_name(15), dump_bot(2), .false.)
          CALL add_svar (vind_out_name(16), dump_bot(3), .false.)
          CALL add_svar (vind_out_name(17), dump_bot(4), .false.)
          CALL add_svar (vind_out_name(18), dump_bot(5), .false.)
          CALL add_svar (vind_out_name(19), dump_taux, .false.)
          CALL add_svar (vind_out_name(20), dump_tauy, .true.)

          ! gol124: auscom coupling fields
          CALL add_svar (vind_out_name(21), dump_swflx, .false.)
          CALL add_svar (vind_out_name(22), dump_lwflx, .false.)
          CALL add_svar (vind_out_name(23), dump_shflx, .false.)
          CALL add_svar (vind_out_name(24), dump_press, .false.)
          CALL add_svar (vind_out_name(25), dump_co2, .false.)
          CALL add_svar (vind_out_name(26), dump_wnd10, .true.)

          CALL end_sdef
      END IF
#endif

      ! ===============================================
      ! Define incoming fields
      ! ===============================================
      vind_in(1) = vind_ocn_sst
      vind_in_name(1) = "ocn_sst"

      vind_in(2:6) = vind_ocn_freezen(1:5)
      vind_in_name(2) = "ofrzn01"
      vind_in_name(3) = "ofrzn02"
      vind_in_name(4) = "ofrzn03"
      vind_in_name(5) = "ofrzn04"
      vind_in_name(6) = "ofrzn05"

      vind_in(7:11) = vind_ocn_snowthickn(1:5)
      vind_in_name(7) = "osnwtn01"
      vind_in_name(8) = "osnwtn02"
      vind_in_name(9) = "osnwtn03"
      vind_in_name(10) = "osnwtn04"
      vind_in_name(11) = "osnwtn05"

      vind_in(12:16) = vind_ocn_hicen(1:5)
      vind_in_name(12) = "ohicn01"
      vind_in_name(13) = "ohicn02"
      vind_in_name(14) = "ohicn03"
      vind_in_name(15) = "ohicn04"
      vind_in_name(16) = "ohicn05"

      vind_in(17) = vind_ocn_u
      vind_in_name(17) = "sunocean"
      vind_in(18) = vind_ocn_v
      vind_in_name(18) = "svnocean"

      vind_in_type(1:16) = "T"
      vind_in_type(17) = "U"
      vind_in_type(18) = "V"

#if defined(ACCESS)
      vind_in(19) = vind_ocn_co2
      vind_in_name(19) = "co2_ia"
      vind_in(20) = vind_ocn_co2fx
      vind_in_name(20) = "co2fx_ia"
      vind_in_type(19:20) = "T"
#endif

#if defined(ACCESS)
! gol124: auscom coupling
      IF (ocn_sss) THEN
          vind_in(21) = vind_sal
          vind_in_name(21) = "ocn_sss"
          vind_in_type(21) = "T"
          max_i = 21
      ELSE
          max_i = 20
      END IF
#endif

#if defined(ACCESS)
      DO I = 1, max_i
#else
      DO I = 1 , 18
#endif
        IF (vind_in_type(I) == "V") THEN
          var_shape(4) = valid_shape_v(2,2)
          partition = partition_id_v
        ELSE
          var_shape(4) = valid_shape(2,2)
          partition = partition_id
        END IF

        CALL prism_def_var_proto (var_id(vind_in(I)),    &
           &        vind_in_name(I), partition,                    &
           &        var_nodims, PRISM_In, var_shape,               &
           &        PRISM_Real, ierror )

      ENDDO

#if defined(ACCESS)
      IF (rdump_enable) THEN
          CALL add_rvar (vind_in_name(1),  dump_sst, .false.)
          CALL add_rvar (vind_in_name(2),  dump_frzn(1), .false.)
          CALL add_rvar (vind_in_name(3),  dump_frzn(2), .false.)
          CALL add_rvar (vind_in_name(4),  dump_frzn(3), .false.)
          CALL add_rvar (vind_in_name(5),  dump_frzn(4), .false.)
          CALL add_rvar (vind_in_name(6),  dump_frzn(5), .false.)
          CALL add_rvar (vind_in_name(7),  dump_snwtn(1), .false.)
          CALL add_rvar (vind_in_name(8),  dump_snwtn(2), .false.)
          CALL add_rvar (vind_in_name(9),  dump_snwtn(3), .false.)
          CALL add_rvar (vind_in_name(10),  dump_snwtn(4), .false.)
          CALL add_rvar (vind_in_name(11),  dump_snwtn(5), .false.)
          CALL add_rvar (vind_in_name(12),  dump_hicn(1), .false.)
          CALL add_rvar (vind_in_name(13),  dump_hicn(2), .false.)
          CALL add_rvar (vind_in_name(14),  dump_hicn(3), .false.)
          CALL add_rvar (vind_in_name(15),  dump_hicn(4), .false.)
          CALL add_rvar (vind_in_name(16),  dump_hicn(5), .false.)
          CALL add_rvar (vind_in_name(17),  dump_suno, .false.)
          CALL add_rvar (vind_in_name(18),  dump_svno, .true.)
          CALL add_rvar (vind_in_name(19),  dump_ocn_co2, .false.)
          CALL add_rvar (vind_in_name(20),  dump_ocn_co2fx, .false.)

          ! gol124: auscom coupling fields
          IF (ocn_sss) THEN
              CALL add_rvar (vind_in_name(21), dump_sss, .false.)
          END IF

          CALL end_rdef
      END IF
#endif

      ! Finish the prism definition phase and perform inter component
      ! integrity checking.

      CALL prism_enddef_proto(ierror)
    END IF  ! coupling through master only or
    ! in parallel without explicit gather/scatter.


    ! Set the start time. For OASIS3 this is simply based on the
    ! number of seconds this run, so at the start of the
    ! run it will be zero regardless of the true model date
    ! It may be that we need to do some cross checking with other
    ! components in a coupled model system, in which case this
    ! is one suitable point to emply such a test.  
    prism_nsec = 0

  END SUBROUTINE OASIS3_grid32

  !======================================================================
  SUBROUTINE oasis3_advance_date32(ierr64)
    !
    ! Description: This routine advances the prism date and calculates
    !              therebuild_ompi133.o402718 appropriate time for use in  put/get operations.
    !
    !======================================================================


    IMPLICIT NONE
    INTEGER (kind=integer64) :: ierr64
    INTEGER (kind=integer32) :: ierror

    ! Advance the date/time as used by prism

    prism_nsec = prism_nsec + NINT(prism_timestep)

    RETURN
  END SUBROUTINE oasis3_advance_date32
  !=======================================================================
  SUBROUTINE OASIS3_PUT32(data_64,rowl_64,rows_64,ierr64)
    !
    ! Description: PUT data to OASIS in 64 bit form (prism double
    !              precision).
    !=======================================================================
    !USE mod_prism_put_proto

    IMPLICIT NONE

    INTEGER (kind=integer64) :: I, J, rowl_64
    INTEGER (kind=integer64) :: rows_64

    INTEGER (kind=integer64) :: ierr64
    INTEGER (kind=integer32) :: ierror, var_id32, prism_nsec32

    REAL (kind=real64),DIMENSION(1:rowl_64,1:rows_64,1) :: data_64
    REAL (kind=real64),DIMENSION(1:rowl_64,1:rows_64) :: data_2d

    ierror=0

    DO J = 1, rows_64
      DO I = 1, rowl_64
        data_2d(I,J) = data_64(I,J,1)
      ENDDO
    ENDDO


    var_id32 = var_id(var_ind)
    prism_nsec32 = prism_nsec

    write(6,*) "OASIS3_PUT32 at sec:", prism_nsec32
    ! PUT the outgoing coupling field
    CALL prism_put_proto(var_id32,prism_nsec32,data_2d,ierror)

    IF (ierror .lt. 0) THEN
      WRITE(6,*) "return code from prism_put_proto(data_64)",ierror
      WRITE(6,*) "field: ", var_id(var_ind),prism_nsec
    END IF

    ierr64 = ierror

    RETURN

  END SUBROUTINE OASIS3_PUT32
  !=======================================================================
  SUBROUTINE OASIS3_GET32(data_recv,rowl_64,rows_64,ierr64)

    !
    ! Description: Get date from OASIS in 64 bit form (double
    !              precision as far as OASIS is concerned, but
    !              standard precision for the UM).
    !=======================================================================
    !USE mod_prism_get_proto

    IMPLICIT NONE

    INTEGER (kind=integer64) :: I,J,rowl_64
    INTEGER (kind=integer64) :: rows_64

    INTEGER (kind=integer64) :: ierr64
    INTEGER (kind=integer32) :: ierror,var_id32, prism_nsec32

    REAL (kind=real64),DIMENSION(1:rowl_64,1:rows_64)::data_recv
    REAL (kind=real32),DIMENSION(1:rowl_64,1:rows_64)::data_2d

    INTEGER (kind=integer64) :: top

#include "cmaxsize.h"
#include "csubmodl.h"
#include "chsunits.h"
#include "chistory.h"
#include "cntlall.h"
#include "cntlatm.h"

    ierror=0
    var_id32=var_id(var_ind)
    prism_nsec32=prism_nsec

    IF (LTIMER) THEN
       top=3
! DEPENDS ON: timer
       CALL TIMER('GET_PROTO',top)
    END IF

    write(6,*), "OASIS3_GET32 at sec: ", prism_nsec32
    ! GET the incoming coupling field
    CALL prism_get_proto(var_id32,prism_nsec32,data_recv,ierror)

    IF (LTIMER) THEN
       top=4
! DEPENDS ON: timer
       CALL TIMER('GET_PROTO',top)
    END IF

    IF (ierror .lt. 0) THEN
      WRITE(6,*) "return code from prism_get_proto(data_64)",ierror
      WRITE(6,*) "field: ", var_id(var_ind),prism_nsec
    END IF

    ierr64 = ierror
    RETURN

  END SUBROUTINE OASIS3_get32

  SUBROUTINE OASIS3_um_init(comm64,ierror,cmessage)

    IMPLICIT NONE
    !======================================================================
    ! Description: Initialisation of OASIS3 for UM model components
    !======================================================================


    INTEGER (kind=integer64), INTENT(out) :: ierror
    INTEGER (kind=integer32)              :: ierr32
    INTEGER (kind=integer32)              :: comm32
    INTEGER (kind=integer32)              :: comp_id
    LOGICAL                               :: pinit
    INTEGER (kind=integer64), INTENT(out) :: comm64
    CHARACTER(len=*),INTENT(out) :: cmessage

    !
    ! Init Component
    !
    ! Character string  descriptors for model (and component)
    ! name. Note OASIS3 is extraordinarily restrictive about this
    ! You mustnt exceed 6 characters in the model name which makes
    ! meaningful names rather difficult.
    CHARACTER(LEN=6) :: model_name , comp_name
    !

!    model_name="toyatm"
    model_name="um7.3x"

    comp_name="atmos"  ! Component name is only required if running
    ! with OASIS4

    ierr32=0
!    CALL prism_init_comp_proto(OASIS_comp_id,'toyatm', ierr32)
    CALL prism_init_comp_proto(OASIS_comp_id,model_name, ierr32)

    IF (ierror /= 0) THEN
      cmessage = "error in prism_init_comp_proto"
    ELSE
      cmessage = "PRISM_init_comp_proto called from oasis3_um_init"
    END IF

    ! OASIS3 does not require a call to initialise the component. 
    ! cf. OASIS4.
    !      write(6,*)'init_oasis - b4 prism_init_comp'
    !      call PRISM_Init_comp (comp_id, comp_name, ierr32 )
    !      ierror = ierr32
    !      write(6,*)'init_oasis - af prism_init_comp',comp_id,ierror
    !      if (ierror /= 0) then
    !        cmessage = "error in PRISM_init_comp from init_oasis4"
    !      else
    !        cmessage = "PRISM_init_comp called from init_oasis4"
    !      end if

    ! Get the communicator allocated to this model component
    CALL PRISM_get_localcomm_proto(comm32,ierr32 )
    comm64=comm32

    RETURN

  END SUBROUTINE OASIS3_um_init

END MODULE OASIS3_ATMOS_INIT
#endif
