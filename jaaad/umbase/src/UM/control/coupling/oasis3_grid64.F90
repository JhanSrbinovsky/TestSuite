#if defined(OASIS3)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE OASIS3_GRID64(comp_id,MODEL_BASIS_TIME,ierror)

  ! Description: Set up grid shape/sizes for defining grids using
  !              OASIS. Computes all the values needed for OASIS3 
  !              grid definition.  However the OASIS3 software
  !              needs 32bit integers, so the oasis calls are 
  !              wrapped inside subroutine grid32_oasis4, compiled 
  !              -dw with explicit kinds, where integers are 
  !              converted from 64-32bit.
  !
  ! Author: R. Hill
  !
  ! Current Code Owner : R. Hill
  !
  !=====================================================================

  USE oasis3_atmos_init
  USE oasis3_atm_data_mod
#if defined(ACCESS)
  USE auscom_cpl_data_mod
  USE dump_sent
  USE dump_received
  USE mpl, Only: MPL_REAL,MPL_COMM_WORLD
#endif

  IMPLICIT NONE

  !
  INTEGER         :: comp_id
  INTEGER         :: ierror
  INTEGER         :: MODEL_BASIS_TIME(6)

#include "csubmodl.h"
#include "parvars.h"
#include "typsize.h"
#include "version.h"
#include "model.h"
#include "cntlatm.h"

  INTEGER            :: grid_type
  INTEGER, PARAMETER :: nGridDims = 3

  INTEGER, SAVE :: valid_shape(2,nGridDims)
  INTEGER, SAVE :: actual_shape(2,nGridDims)

  ! local variables...
  INTEGER :: i,j,k

#if defined(ACCESS)

#include "c_0_dg_c.h"
! DEPENDS ON: OASIS3_UMVARS_INIT
  call oasis3_umvars_init() !broadcast data in "input_atm.nml" 
  access_tfs = TFS
#endif

  ! valid_shape = actual_shape since I'm not using halos.
  valid_shape(1,1) = 1
  valid_shape(1,2) = 1
  valid_shape(1,3) = 1

  ! Valid shapes are the global dimsions for OASIS3
  ! Since, if we write a grids.nc, areas.nc or masks.nc
  ! file we'll do it via PE 0 only.

  IF (L_COUPLE_MASTER) THEN
    write(6,*) 'grid64: setup for master coupling, pe = ',mype
    valid_shape(2,1) = glsize(1,fld_type_p)
    valid_shape(2,2) = glsize(2,fld_type_p)
    valid_shape(2,3) = 1
  ELSE
    write(6,*) 'grid64: setup for parallel coupling, pe = ',mype
    valid_shape(2,1) = lasize(1,fld_type_p,halo_type_no_halo)
    valid_shape(2,2) = lasize(2,fld_type_p,halo_type_no_halo)
    valid_shape(2,3) = 1
  END IF
  write(6,*) 'grid64: partition shape, pe,x,y = ',mype, &
 & valid_shape(2,1),valid_shape(2,2)

  ! Now set up the actual shape. These can contain halos
  ! if required and therefore must be >= valid shape equivalents
  ! but for now we keep things the same.
  actual_shape(1,1) = 1
  actual_shape(1,2) = 1
  actual_shape(1,3) = 1
  actual_shape(2,1) = valid_shape(2,1)     ! I (E-W)
  actual_shape(2,2) = valid_shape(2,2)     ! J (S-N)
  actual_shape(2,3) = valid_shape(2,3)     ! K (Vertical)

  ! Detect whether this PE is at the northern limit of our
  ! global domain.
  no_neighbour_N = at_extremity(Pnorth)

  L_MASTER = L_COUPLE_MASTER

#if defined(ACCESS)
  IF (sdump_enable) THEN
      CALL start_sdump (glsize(1,fld_type_p), glsize(2,fld_type_p), mype)
  END IF
  IF (rdump_enable) THEN
      CALL start_rdump (glsize(1,fld_type_p), glsize(2,fld_type_p), mype)
  END IF
#endif

  ! Call the routine where our grid will actually be defined
  CALL oasis3_grid32(MODEL_BASIS_TIME,ierror                        &
     &, valid_shape,actual_shape,nGridDims                              &
     &, glsize(1,fld_type_p),datastart,Ndim_max, MYPE)

  RETURN
END SUBROUTINE OASIS3_GRID64
#endif
