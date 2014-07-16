MODULE dump_sent

   USE netcdf

   IMPLICIT NONE

   character (len = *), parameter :: FILE_NAME = "fields_a2i_in_atm.nc"

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
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ds_vlats, ds_lats, ds_lons
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

   integer(kind=8) :: ds_glx, ds_gly

CONTAINS

SUBROUTINE scheck(status)

   USE oasis3_atm_data_mod

   IMPLICIT NONE

#include "cntlatm.h"

   integer, intent ( in) :: status
     
   IF (LMYPE .NE. OASIS_CNTLPE) RETURN

   if(status /= nf90_noerr) then
      print *, 'sdump: ',trim(nf90_strerror(status))
      stop 2
   end if
END SUBROUTINE scheck

SUBROUTINE start_sdump (alons, alats, mype)

   USE oasis3_atm_data_mod

   IMPLICIT NONE

#include "cntlatm.h"

   integer(kind=8) :: alons, alats

   integer (kind=8) :: nlons, nlats, mype

   real(kind=8) :: lat_delta, lon_delta
   real(kind=8) :: vrange(2), nodata
   real(kind=4) :: fill
   integer(kind=4) :: lat, lon

   ds_glx = alons
   ds_gly = alats
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

   ALLOCATE (ds_lats(nlats))
   do lat = 1, NLATS
      ds_lats(lat) = START_LAT + (lat - 1) * lat_delta
   end do

   ALLOCATE (ds_vlats(nlats-1))
   do lat = 1, NLATS-1
      ds_vlats(lat) = ds_lats(lat+1)
   end do

   ALLOCATE (ds_lons(nlons))
   do lon = 1, NLONS
      ds_lons(lon) = START_LON + (lon - 1) * lon_delta
   end do

   lat = NLATS
   lon = NLONS

   call scheck( nf90_create(FILE_NAME, NF90_CLOBBER, ncid) )
   call scheck( nf90_def_dim(ncid, LAT_NAME, lat, lat_dimid) )
   ! v-type grid
   lat = NLATS - 1
   call scheck( nf90_def_dim(ncid, LATV_NAME, lat, latv_dimid) )
   call scheck( nf90_def_dim(ncid, LON_NAME, lon, lon_dimid) )
   call scheck( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )

   call scheck( nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, lat_varid) )
   call scheck( nf90_def_var(ncid, LATV_NAME, NF90_REAL, latv_dimid, latv_varid) )
   call scheck( nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lon_varid) )

   call scheck( nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS) )
   call scheck( nf90_put_att(ncid, latv_varid, UNITS, LAT_UNITS) )
   call scheck( nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS) )

   dimids = (/ lon_dimid, lat_dimid, rec_dimid /)
   vdimids = (/ lon_dimid, latv_dimid, rec_dimid /)

   count = (/ lon, lat, 1 /)
   last_rec = 1
   start = (/ 1, 1, last_rec /)

END SUBROUTINE start_sdump

SUBROUTINE add_svar (vname, varid, vgrid)

   USE oasis3_atm_data_mod

   IMPLICIT NONE

#include "cntlatm.h"

   character(len=*) :: vname
   integer(kind=8) :: varid
   logical :: vgrid

   integer(kind=4) :: lvarid

   IF (LMYPE .NE. OASIS_CNTLPE) RETURN

   if (vgrid) then
       call scheck( nf90_def_var(ncid, vname, NF90_REAL, vdimids, lvarid))
   else
       call scheck( nf90_def_var(ncid, vname, NF90_REAL, dimids, lvarid))
   end if
   varid = lvarid

END SUBROUTINE add_svar

SUBROUTINE end_sdef

   USE oasis3_atm_data_mod

   IMPLICIT NONE

#include "cntlatm.h"

   IF (LMYPE .NE. OASIS_CNTLPE) RETURN

   call scheck( nf90_enddef(ncid) )

   call scheck( nf90_put_var(ncid, lat_varid, ds_lats) )
   call scheck( nf90_put_var(ncid, latv_varid, ds_vlats) )
   call scheck( nf90_put_var(ncid, lon_varid, ds_lons) )

END SUBROUTINE end_sdef

SUBROUTINE inc_srec

   USE oasis3_atm_data_mod

   IMPLICIT NONE

#include "cntlatm.h"

   IF (LMYPE .NE. OASIS_CNTLPE) RETURN

   last_rec = last_rec + 1

END SUBROUTINE inc_srec

! Output received data to netcdf file
SUBROUTINE write_sdata ( varid, ft, data_64, l_cols, l_rows, vgrid)

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

   INTEGER :: info
   CHARACTER*(80) :: cmessage
   REAL(kind=8) ,ALLOCATABLE ,DIMENSION(:,:)  :: data_send

   integer(kind=8) :: halo
   integer(kind=4) :: lvarid

   integer(kind=8) :: gly

   if (vgrid) then
           gly = ds_gly - 1
   else
           gly = ds_gly
   end if
   count = (/ ds_glx, gly, 1 /)

   ! gather all data on single cpu that will write to the disk

   ALLOCATE (DATA_SEND(ds_glx, gly))

   ! gol124: this module is compiled with default 32 bit integers
   ! for netcdf compatibility
   ! but some subs/constants are written assuming integers by
   ! default are 64 bit so we need this line below to enforce it!
   halo = halo_type_no_halo
! DEPENDS ON: gather_field
   CALL GATHER_FIELD(data_64(1,1,1),data_send,                    &
                  l_cols,                                         &
                  l_rows,                                         &
                  ds_glx,                                         &
                  gly,                                            &
                  ft,halo,                                        &
                  OASIS_CNTLPE,GC_ALL_PROC_GROUP,info,cmessage)

   IF (LMYPE == OASIS_CNTLPE) THEN
      lvarid = varid
      start = (/ 1, 1, last_rec /)
      call scheck( nf90_put_var(ncid, lvarid, data_send, &
              start = start, &
              count = count) )
      call scheck( nf90_sync(ncid))
   END IF

   deallocate(data_send)

   RETURN
END SUBROUTINE write_sdata

SUBROUTINE end_sdump

   USE oasis3_atm_data_mod

   IMPLICIT NONE

#include "cntlatm.h"

   IF (LMYPE .NE. OASIS_CNTLPE) RETURN

   call scheck( nf90_close(ncid) )
   deallocate (ds_lats)
   deallocate (ds_vlats)
   deallocate (ds_lons)

END SUBROUTINE end_sdump

END MODULE dump_sent
