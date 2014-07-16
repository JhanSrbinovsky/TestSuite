#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(UTILIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! + Scatters any type of field from one processor to many processors
!
! Subroutine interface:

      SUBROUTINE GENERAL_SCATTER_FIELD (                                &
     &  LOCAL_FIELD , GLOBAL_FIELD ,                                    &
     &  LOCAL_SIZE , GLOBAL_SIZE ,                                      &
     &  LEVELS ,                                                        &
     &  ADDR_INFO ,                                                     &
     &  SCATTER_PE ,                                                    &
     &  ICODE , CMESSAGE)

      IMPLICIT NONE

! Description:
! Takes a general field on a single processor and decomposes (scatters)
! it over many processors
!
! Current code owner : P.Burton
!
! Subroutine arguments:

      INTEGER                                                           &
     &  GLOBAL_SIZE                                                     &
                        ! IN:  size of GLOBAL FIELD
     &, LEVELS                                                          &
                           ! IN:  How many levels of data to do
     &, SCATTER_PE(LEVELS)                                              &
                           ! IN: Which PE to scatter each level from
     &, LOCAL_SIZE                                                      &
                        ! OUT: size of LOCAL_FIELD
     &, ICODE           ! OUT: return code, 0=OK

! Required for dimensioning ADDR_INFO
#include "d1_addr.h"

      INTEGER                                                           &
     &  ADDR_INFO(D1_LIST_LEN)   ! IN: addressing info about field

      REAL                                                              &
     &  GLOBAL_FIELD(GLOBAL_SIZE,*)                                     &
                                    ! IN:  field to scatter
     &, LOCAL_FIELD(*)            ! OUT: my local part of fld

      CHARACTER*(80)                                                    &
     &  CMESSAGE        ! OUT : Error message

! Parameters and common blocks

#include "cppxref.h"
#include "stparam.h"

#include "gccom.h"

#include "parvars.h"
#include "typsize.h"
#include "atm_lsm.h"

! Local variables

      INTEGER                                                           &
     &  grid_type                                                       &
                  ! grid type of field being scattered
     &, grid_code                                                       &
                  ! ppx grid code of field being scattered
     &, halo_type                                                       &
                  ! halo type of field being scattered
     &, info                                                            &
             ! return code from GCOM routines
     &, dummy                                                           &
              ! dummy variables - ignored return values
     &, north,south,east,west                                           &
                               ! domain limits for STASH output
     &, mean_type                                                       &
                  ! spatial meaning type on diagnostic
     &, k                                                               &
           ! loop over levels
     &, my_k                                                            &
             ! value of k for GLOBAL_FIELD on this PE
     &, my_k_temp                                                       &
                  ! Temp. my_k for safe references.
     &, rim_type                                                        &
                  ! RIMWIDTH type for LBC field
     &, loc_len_rim                                                     &
                    ! length of local rimdata for single level
     &, iproc                                                           &
                  ! Loop variable over processors
     &, i         ! Loop variable for debugging

      INTEGER                                                           &
     &  GET_FLD_TYPE  ! function for finding field type

      REAL                                                              &
     &  buf_expand(Max2DFieldSize)                                      &
     &, buf_expand_local(Max2DFieldSize)

!===================================================================

      grid_code=ADDR_INFO(d1_grid_type)
! DEPENDS ON: get_fld_type
      grid_type=GET_FLD_TYPE(grid_code)
      halo_type=ADDR_INFO(d1_halo_type)

!-------------------------------------------------------------------

! Timeseries data

      IF ((ADDR_INFO(d1_object_type)  ==  diagnostic) .AND.             &
     &   ((ADDR_INFO(d1_proc_no_code)  ==  st_time_series_code).OR.     &
     &   (ADDR_INFO(d1_proc_no_code)  ==  st_time_series_mean))) THEN
! Copy the data from GLOBAL_FIELD to PE 0

! Multi-level fields not supported here
      IF (LEVELS  >   1) THEN
        WRITE(6,*) 'GENERAL_SCATTER_FIELD : Cannot have more than ',    &
     &             '1 level for scattering timeseries data.'
        ICODE=1
        CMESSAGE='GENERAL_SCATTER_FIELD : Multi-level timeseries field'
        GOTO 9999
      ENDIF

         IF (mype  ==  SCATTER_PE(1)) THEN

          info=GC_NONE
          CALL GC_RSEND(99,GLOBAL_SIZE,0,info,LOCAL_FIELD,              &
     &                  GLOBAL_FIELD)

        ENDIF

        CALL GC_SSYNC(nproc,info)

        IF (mype  ==  0) THEN

          info=GC_NONE
          CALL GC_RRECV(99,GLOBAL_SIZE,SCATTER_PE,info,                 &
     &                  LOCAL_FIELD,GLOBAL_FIELD)

        ENDIF

        LOCAL_SIZE=GLOBAL_SIZE

!-------------------------------------------------------------------

! Surface (land points only) fields

      ELSEIF                                                            &
     &  (grid_code  ==  ppx_atm_compressed) THEN

       my_k=0
       DO k=1,LEVELS

       IF (mype  ==  SCATTER_PE(k)) THEN

         my_k=my_k+1

! DEPENDS ON: from_land_points
         CALL FROM_LAND_POINTS(buf_expand,GLOBAL_FIELD(1,my_k),         &
     &                         atmos_landmask,                          &
     &                          glsize(1,grid_type)*glsize(2,grid_type),&
     &                          dummy)

! SCATTER_PE now contains the expanded version of the full field

        ENDIF

! Now scatter this to all the other processors, putting the local
! part of the field into the array buf_expand_local

! DEPENDS ON: scatter_field
        CALL SCATTER_FIELD(buf_expand_local , buf_expand,               &
     &                     lasize(1,grid_type,halo_type),               &
     &                     lasize(2,grid_type,halo_type),               &
     &                     glsize(1,grid_type),glsize(2,grid_type),     &
     &                     grid_type,halo_type,                         &
     &                     SCATTER_PE(k),GC_ALL_PROC_GROUP,             &
     &                     ICODE,CMESSAGE)

        IF (ICODE  /=  0) THEN
          WRITE(6,*)                                                    &
     &      'GENERAL_SCATTER_FIELD : Error detected in call to ',       &
     &      'SCATTER_FIELD '
          WRITE(6,*) 'Return code : ',ICODE
          WRITE(6,*) 'Error message : ',CMESSAGE

          ICODE=1
          CMESSAGE='GENERAL_SCATTER_FIELD : Error scattering field'
          GOTO 9999
        ENDIF

! Fill the halo area of the local field

! DEPENDS ON: swap_bounds
        CALL SWAP_BOUNDS(buf_expand_local,                              &
     &                  lasize(1,grid_type,halo_type),                  &
     &                  lasize(2,grid_type,halo_type),                  &
     &                  1,                                              &
     &                  halosize(1,halo_type),                          &
     &                  halosize(2,halo_type),                          &
     &                  grid_type,.FALSE.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(buf_expand_local,                      &
     &                     lasize(1,grid_type,halo_type),               &
     &                     lasize(2,grid_type,halo_type),               &
     &                     1,halosize(1,halo_type),                     &
     &                     halosize(2,halo_type))


! Pack the local field down to local land points and put
! the packed field into LOCAL_FIELD

#if defined (CUMF) || defined (PUMF) || defined (MERGE)
        local_land_field=dummy
#endif
! DEPENDS ON: to_land_points
       CALL TO_LAND_POINTS(buf_expand_local,                            &
     &                     LOCAL_FIELD(1+(k-1)*local_land_field),       &
     &                      atmos_landmask_local,                       &
     &                      lasize(1,grid_type,halo_type)*              &
     &                      lasize(2,grid_type,halo_type),              &
     &                      dummy)
        ENDDO ! k : loop over levels

        local_size = local_land_field

!-------------------------------------------------------------------

! Atmosphere Lateral boundary fields

      ELSEIF                                                            &
     &  ((grid_code  ==  ppx_atm_lbc_theta) .OR.                        &
     &   (grid_code  ==  ppx_atm_lbc_u) .OR.                            &
     &   (grid_code  ==  ppx_atm_lbc_v) .OR.                            &
     &   (grid_code  ==  ppx_atm_lbc_orog)) THEN

        IF (grid_code  ==  ppx_atm_lbc_orog) THEN
          rim_type=rima_type_orog
        ELSE
          rim_type=rima_type_norm
        ENDIF

! DEPENDS ON: scatter_atmos_lbcs
        CALL SCATTER_ATMOS_LBCS(GLOBAL_FIELD,LOCAL_FIELD,               &
     &                          global_LENRIMA(grid_type,halo_type,     &
     &                                         rim_type),               &
     &                          LENRIMA(grid_type,halo_type,rim_type),  &
     &                          LEVELS,LEVELS,                          &
     &                          grid_type,halo_type,rim_type,           &
     &                          SCATTER_PE,                             &
     &                          ICODE,CMESSAGE)

        IF (ICODE  /=  0) THEN
          WRITE(6,*)                                                    &
     &      'GENERAL_SCATTER_FIELD : Error detected in call to ',       &
     &      'SCATTER_ATMOS_LBCS '
          WRITE(6,*) 'Return code : ',ICODE
          WRITE(6,*) 'Error message : ',CMESSAGE

          ICODE=2
          CMESSAGE='GENERAL_SCATTER_FIELD : Error scattering LBCs'
          GOTO 9999
        ENDIF


        LOCAL_SIZE=LENRIMDATA_A

!-------------------------------------------------------------------

      ELSEIF                                                            &
!     atmosphere grids
     &  ((grid_code  ==  ppx_atm_tall)  .OR.                            &
     &   (grid_code  ==  ppx_atm_tland) .OR.                            &
     &   (grid_code  ==  ppx_atm_tsea)  .OR.                            &
     &   (grid_code  ==  ppx_atm_uall) .OR.                             &
     &   (grid_code  ==  ppx_atm_cuall) .OR.                            &
     &   (grid_code  ==  ppx_atm_cvall) .OR.                            &
     &   (grid_code  ==  ppx_atm_tzonal) .OR.                           &
     &   (grid_code  ==  ppx_atm_tmerid) .OR.                           &
     &   (grid_code  ==  ppx_atm_compressed) .OR.                       &
     &   (grid_code  ==  ppx_atm_ozone) .OR.                            &
     &   (grid_code  ==  ppx_atm_river))                                &
     &   THEN

        LOCAL_SIZE=ADDR_INFO(d1_length)/ADDR_INFO(d1_no_levels)


        IF (ADDR_INFO(d1_object_type)  ==  diagnostic) THEN
          north=ADDR_INFO(d1_north_code)
          south=ADDR_INFO(d1_south_code)
          east=ADDR_INFO(d1_east_code)
          west=ADDR_INFO(d1_west_code)

          mean_type=ADDR_INFO(d1_gridpoint_code)/10
          IF (mean_type  ==  2) THEN ! zonal mean
            east=west
          ELSEIF (mean_type  ==  3) THEN ! meridional mean
            south=north
          ELSEIF (mean_type  >=  4) THEN ! field/global mean
            east=west
            south=north
          ENDIF

        ELSE
          north=glsize(2,grid_type)
          west=1
          east=glsize(1,grid_type)
          south=1
        ENDIF

        my_k=0
        DO k=1,LEVELS

        IF (mype  ==  SCATTER_PE(k)) THEN
           my_k=my_k+1
           my_k_temp = my_k
        ELSE
           my_k_temp=1
        END IF

!       STASH_SCATTER_FIELD can distribute whole fields, or subarea
!       fields

! DEPENDS ON: stash_scatter_field
        CALL STASH_SCATTER_FIELD(                                       &
     &    LOCAL_FIELD(1+(k-1)*LOCAL_SIZE),GLOBAL_FIELD(1,my_k_temp),    &
     &    LOCAL_SIZE,GLOBAL_SIZE,1,                                     &
     &    north,east,south,west,                                        &
     &    grid_code,halo_type,SCATTER_PE(k),ICODE,CMESSAGE)

        IF (ICODE  /=  0) THEN
          WRITE(6,*)                                                    &
     &      'GENERAL_SCATTER_FIELD : Error detected in call to ',       &
     &      'STASH_SCATTER_FIELD '
          WRITE(6,*) 'Return code : ',ICODE
          WRITE(6,*) 'Error message : ',CMESSAGE

          ICODE=4
          CMESSAGE='GENERAL_SCATTER_FIELD : Error scattering field'
          GOTO 9999
        ENDIF
      ENDDO ! k : loop over levels

!-------------------------------------------------------------------
! Any other type of field
      ELSE

        ICODE=10
        CMESSAGE='GENERAL_SCATTER_FIELD : Field type not recognized'

      ENDIF

 9999 CONTINUE

      RETURN

      END SUBROUTINE GENERAL_SCATTER_FIELD
#endif
