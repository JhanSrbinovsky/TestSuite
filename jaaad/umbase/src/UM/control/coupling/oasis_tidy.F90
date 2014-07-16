#if defined(OASIS3) || defined(OASIS4)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE oasis_tidy(                                        &
#include "argd1.h"
#include "argsts.h"
#include "argduma.h"
#include "argptra.h"
  & icode,cmessage)

  USE oasis3_atm_data_mod
#if defined(ACCESS)
  USE auscom_cpl_data_mod
  USE dump_sent
  USE dump_received
#endif

  IMPLICIT NONE

  !
  ! Description: Tidy up at the end of an OASIS-based coupled run.
  !
  ! Author: R. Hill
  ! Current Code Owner : R. Hill 
  !
  !--------------------------------------------------------------------
#include "csubmodl.h"
#include "cntlatm.h"
#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
#include "typsts.h"
#include "typduma.h"
#include "typptra.h"

  ! Dummy arguments - only used with stub version of this routine
  ! and hence theoretically never used.

  INTEGER       icode       ! OUT - Error return code
  CHARACTER*(*) cmessage    ! OUT - Error return message



  ! Deallocate arrays used in the processing of
  ! outgoing coupling data
  DEALLOCATE(taux)
  DEALLOCATE(tauy)
  DEALLOCATE(solar2d)
  DEALLOCATE(blue2d)
  DEALLOCATE(evap2d)
  DEALLOCATE(longwave2d)
  DEALLOCATE(sensible2d)
  DEALLOCATE(heatflux)
  DEALLOCATE(sublim)
  DEALLOCATE(latentflux)
  DEALLOCATE(rainls)
  DEALLOCATE(snowls)
  DEALLOCATE(rainconv)
  DEALLOCATE(snowconv)
  DEALLOCATE(totalrain)
  DEALLOCATE(totalsnow)
  DEALLOCATE(riverout)
  DEALLOCATE(wme)
  DEALLOCATE(um_co2)
  DEALLOCATE(um_wnd10)
  DEALLOCATE(topmeltn)
  DEALLOCATE(botmeltn)

#if defined(ACCESS)
! gol124: auscom coupling (outgoing)
  DEALLOCATE(auscom_swflx)
  DEALLOCATE(auscom_lwflx)
  DEALLOCATE(auscom_shflx)
  DEALLOCATE(auscom_pressure)
#endif

  ! Deallocate arrays used in the processing of
  ! incoming coupling data
  DEALLOCATE(ocn_hice)
  DEALLOCATE(ocn_hicen)
  DEALLOCATE(ocn_freeze)
  DEALLOCATE(ocn_freezen)
  DEALLOCATE(ocn_snowthick)
  DEALLOCATE(ocn_snowthickn)
  DEALLOCATE(ocn_sst)
  DEALLOCATE(ocn_sst_orig)
  DEALLOCATE(ocn_u)
  DEALLOCATE(ocn_v)
  DEALLOCATE(ocn_co2)
  DEALLOCATE(ocn_co2fx)
  DEALLOCATE(tstar_local)
  DEALLOCATE(tstar_ssi)
  DEALLOCATE(tstar_sice_local)
  DEALLOCATE(tstar_land_local)
  DEALLOCATE(fland_loc)

#if defined(ACCESS)

! gol124: auscom coupling (incoming)
  DEALLOCATE(auscom_salinity)

  IF (sdump_enable) THEN
      call end_sdump
  END IF
  IF (rdump_enable) THEN
      call end_rdump
  END IF

#endif

END SUBROUTINE oasis_tidy
#endif
