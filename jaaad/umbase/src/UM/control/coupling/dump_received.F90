MODULE dump_received

   USE netcdf

   IMPLICIT NONE

   character (len = *), parameter :: FILE_NAME = "fields_i2a_in_atm.nc"

   integer(kind=4), parameter :: NDIMS = 3
   integer(kind=4), parameter :: NRECS = 12
!   integer(kind=4), parameter :: NLONS = 96, NLATS = 73

   character (len = *), parameter :: LAT_NAME = "latitude"
   character (len = *), parameter :: LATV_NAME = "latitude_v"
   character (len = *), parameter :: LON_NAME = "longitude"
   character (len = *), parameter :: REC_NAME = "time"
   integer(kind=4) :: lon_dimid, lat_dimid, latv_dimid, rec_dimid

   integer(kind=4) :: start(NDIMS), count(NDIMS)

!   real(kind=8) :: lats(NLATS), lons(NLONS)
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dr_vlats, dr_lats, dr_lons
   integer(kind=4) :: lon_varid, lat_varid, latv_varid

   character (len = *), parameter :: TEMP_NAME="temp"
   character (len = *), parameter :: RAW_NAME ="raw"

   integer(kind=4) :: temp_varid
   integer(kind=4) :: temp_varid_raw
   integer(kind=4) :: dimids(NDIMS)
   integer(kind=4) :: vdimids(NDIMS)

   character (len = *), parameter :: UNITS = "units"
   character (len = *), parameter :: RANGE = "valid_range"
   character (len = *), parameter :: FILLV = "_FillValue"
   character (len = *), parameter :: MISSING = "missing_value"
   character (len = *), parameter :: TEMP_UNITS = "kelvin"
   character (len = *), parameter :: LAT_UNITS = "degrees_north"
   character (len = *), parameter :: LON_UNITS = "degrees_east"

   real(kind=8), parameter :: START_LAT = -90.0, START_LON = 0.0

   integer(kind=4) :: ncid

   integer(kind=4) :: last_rec

   integer(kind=4) :: lmype

   integer(kind=8) :: dr_glx, dr_gly
CONTAINS

SUBROUTINE rcheck(status)

   USE oasis3_atm_data_mod

   IMPLICIT NONE

#include "cntlatm.h"

   integer, intent ( in) :: status
     
   IF (LMYPE .NE. OASIS_CNTLPE) RETURN

   if(status /= nf90_noerr) then
      print *, 'sdump: ',trim(nf90_strerror(status))
      stop 2
   end if
END SUBROUTINE rcheck

SUBROUTINE start_rdump (alons, alats, mype)

   USE oasis3_atm_data_mod

   IMPLICIT NONE

#include "cntlatm.h"

   integer(kind=8)  :: alons, alats

   integer (kind=8) :: nlons, nlats, mype

   real(kind=8) :: lat_delta, lon_delta
   real(kind=8) :: vrange(2), nodata
   real(kind=4) :: fill
   integer(kind=4) :: lat, lon

   dr_glx = alons
   dr_gly = alats
   nlons = alons
   nlats = alats
   lmype = mype

   IF (LMYPE .NE. OASIS_CNTLPE) RETURN

!   nodata = -1.0737e+09
   nodata = -1.0e+34
!   fill = -1.0e+34
   fill = NF90_FILL_REAL
   lat_delta = 180.0 / NLATS
   lon_delta = 360.0 / NLONS

   ALLOCATE (dr_lats(nlats))
   do lat = 1, NLATS
      dr_lats(lat) = START_LAT + (lat - 1) * lat_delta
   end do

   ALLOCATE (dr_vlats(nlats-1))
   do lat = 1, NLATS-1
      dr_vlats(lat) = dr_lats(lat+1)
   end do

   ALLOCATE (dr_lons(nlons))
   do lon = 1, NLONS
      dr_lons(lon) = START_LON + (lon - 1) * lon_delta
   end do

   lat = NLATS
   lon = NLONS

   call rcheck( nf90_create(FILE_NAME, NF90_CLOBBER, ncid) )
   call rcheck( nf90_def_dim(ncid, LAT_NAME, lat, lat_dimid) )
   ! v-type grid
   lat = NLATS - 1
   call rcheck( nf90_def_dim(ncid, LATV_NAME, lat, latv_dimid) )
   call rcheck( nf90_def_dim(ncid, LON_NAME, lon, lon_dimid) )
   call rcheck( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )

   call rcheck( nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, lat_varid) )
   call rcheck( nf90_def_var(ncid, LATV_NAME, NF90_REAL, latv_dimid, latv_varid) )
   call rcheck( nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lon_varid) )

   call rcheck( nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS) )
   call rcheck( nf90_put_att(ncid, latv_varid, UNITS, LAT_UNITS) )
   call rcheck( nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS) )

   dimids = (/ lon_dimid, lat_dimid, rec_dimid /)
   vdimids = (/ lon_dimid, latv_dimid, rec_dimid /)

   count = (/ lon, lat, 1 /)
   last_rec = 1
   start = (/ 1, 1, last_rec /)

END SUBROUTINE start_rdump

SUBROUTINE add_rvar (vname, varid, vgrid)

   USE oasis3_atm_data_mod

   IMPLICIT NONE

#include "cntlatm.h"

   character(len=*) :: vname
   integer(kind=8) :: varid
   logical :: vgrid

   integer(kind=4) :: lvarid

   IF (LMYPE .NE. OASIS_CNTLPE) RETURN

   if (vgrid) then
       call rcheck( nf90_def_var(ncid, vname, NF90_REAL, vdimids, lvarid))
   else
       call rcheck( nf90_def_var(ncid, vname, NF90_REAL, dimids, lvarid))
   end if
   varid = lvarid

END SUBROUTINE add_rvar

SUBROUTINE end_rdef

   USE oasis3_atm_data_mod

   IMPLICIT NONE

#include "cntlatm.h"

   IF (LMYPE .NE. OASIS_CNTLPE) RETURN

   call rcheck( nf90_enddef(ncid) )

   call rcheck( nf90_put_var(ncid, lat_varid, dr_lats) )
   call rcheck( nf90_put_var(ncid, latv_varid, dr_vlats) )
   call rcheck( nf90_put_var(ncid, lon_varid, dr_lons) )

END SUBROUTINE end_rdef

SUBROUTINE inc_rrec

   USE oasis3_atm_data_mod

   IMPLICIT NONE

#include "cntlatm.h"

   IF (LMYPE .NE. OASIS_CNTLPE) RETURN

   last_rec = last_rec + 1

END SUBROUTINE inc_rrec

! Output received data to netcdf file
SUBROUTINE write_rdata ( varid, ft, data_64, l_cols, l_rows, vgrid)

   USE oasis3_atm_data_mod
 
   IMPLICIT NONE

#include "parvars.h"
#include "decomptp.h"
#include "decompdb.h"

#include "cntlatm.h"

   integer(kind=8) :: varid ! netcdf variable id
   integer(kind=8) :: ft ! field type
   integer(kind=8) :: l_cols, l_rows ! local buffer
   logical :: vgrid
!   integer(kind=8) :: s_cols, s_rows ! send buffer
   real(kind=8) :: data_64(l_cols,l_rows,1)

   INTEGER(kind=8) :: halo
   INTEGER :: info
   CHARACTER*(80) :: cmessage
   REAL(kind=8) ,ALLOCATABLE ,DIMENSION(:,:)  :: data_send

   integer(kind=4) :: lvarid

   integer(kind=8) :: gly

   if (vgrid) then
           gly = dr_gly - 1
   else
           gly = dr_gly
   end if
   count = (/ dr_glx, gly, 1 /)

   ALLOCATE (DATA_SEND(dr_glx,gly))

   ! gol124: this module is compiled with default 32 bit integers
   ! for netcdf compatibility
   ! but some subs/constants are written assuming integers by
   ! default are 64 bit so we need this line below to enforce it!
   halo = halo_type_no_halo
! DEPENDS ON: gather_field
   CALL GATHER_FIELD(data_64(1,1,1),data_send,                    &
                  l_cols,                                         &
                  l_rows,                                         &
                  dr_glx,                                         &
                  gly,                                            &
                  ft,halo,                                        &
                  OASIS_CNTLPE,GC_ALL_PROC_GROUP,info,cmessage)

   IF (LMYPE == OASIS_CNTLPE) THEN
      lvarid = varid
      start = (/ 1, 1, last_rec /)
      call rcheck( nf90_put_var(ncid, lvarid, data_send, &
              start = start, &
              count = count) )
      call rcheck( nf90_sync(ncid))
   END IF

   deallocate(data_send)

   RETURN
END SUBROUTINE write_rdata

SUBROUTINE end_rdump

   USE oasis3_atm_data_mod

   IMPLICIT NONE

#include "cntlatm.h"

   IF (LMYPE .NE. OASIS_CNTLPE) RETURN

   call rcheck( nf90_close(ncid) )
   deallocate (dr_lats)
   deallocate (dr_vlats)
   deallocate (dr_lons)

END SUBROUTINE end_rdump

END MODULE dump_received
