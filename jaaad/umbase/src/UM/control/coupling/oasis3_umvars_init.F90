#if defined(OASIS3) && defined(ACCESS)


SUBROUTINE oasis3_umvars_init()

  USE oasis3_atm_data_mod
  USE auscom_cpl_data_mod

IMPLICIT NONE

#include "parvars.h"
#include "c_0_dg_c.h"


  INTEGER :: ierr, my_comm
  logical :: ll_exist
  integer, dimension(10) :: logical_vals
  real, dimension(10) :: real_vals

! gol124: reset AUSCOM coupling options to default values
  access_tfs = TFS
  xfactor = 1.0       !dhb599: actually already defined (default)
  ocn_sss=.false.
  sdump_enable = .false.
  rdump_enable = .false.

! gol124: read AUSCOM coupling options
  IF (MYPE.EQ.OASIS_CNTLPE) THEN
      INQUIRE (file="input_atm.nml", exist=ll_exist)

      IF (ll_exist) THEN
          open(unit=99,file="input_atm.nml",form="formatted",status="old")
          read (99, coupling)
          close(unit=99)
      END IF 
! convert to Kelvin as expected
      IF (access_tfs < 0.0) THEN
          access_tfs = ZeroDegC + access_tfs
      END IF

  END IF

  logical_vals(1) = ocn_sss
  logical_vals(2) = sdump_enable
  logical_vals(3) = rdump_enable
  real_vals(1) = access_tfs
  real_vals(2) = xfactor
  real_vals(3) = SC
  real_vals(4) = co2_init
  real_vals(5) = VOLCTS_val

  call gc_ibcast(1, 3, 0, nproc, ierr, logical_vals)
  call gc_rbcast(2, 5, 0, nproc, ierr, real_vals)

  ocn_sss = logical_vals(1)
  sdump_enable = logical_vals(2)
  rdump_enable = logical_vals(3)
  access_tfs = real_vals(1)
  xfactor = real_vals(2)
  SC = real_vals(3)
  co2_init = real_vals(4)
  VOLCTS_val = real_vals(5)

  write (6,*) "logical values from bcast: ocn_sss, isdump_enable, rdump_enable"
  write (6,*)  ocn_sss, sdump_enable, rdump_enable
  write (6,*) "real values from bcast: access_tfs, xfactor, SC, co2_init,VOLCTS_val"
  write (6,*)  access_tfs, xfactor, SC, co2_init,VOLCTS_val

END 

#endif

