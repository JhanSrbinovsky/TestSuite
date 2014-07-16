! TYPOCAC declares variables used by ocean assimilation
      REAL :: O_OBS(JO_MAX_OBS_VAL)   ! main obs array
      REAL :: AICE_OBS(IMT,JMT,JO_NMAX_OBS_ICE)
      REAL :: ICE_OBS_DYS(JO_NMAX_OBS_ICE) ! obs time - whole days
      REAL :: ICE_OBS_SCS(JO_NMAX_OBS_ICE) ! obs time - rem secs
      REAL :: O_COV(JO_LEN_COV)       ! climate and covariance flds
      REAL :: O_LON_C(JO_MAX_COLS_C)  ! climate grid longitudes
      REAL :: O_LAT_C(JO_MAX_ROWS_C)  ! climate grid latitudes
      REAL :: O_DEP_LEVS_C(JO_MAX_LEVS_C) ! climate grid depths
      REAL :: O_LON_M(IMT)     ! model grid longitudes
      REAL :: O_LAT_M(JMT)     ! model grid latitudes
      REAL :: O_DEP_LEVS_M(KM) ! model tracer grid depths
#include "comocac.h"
! TYPOCAC end
