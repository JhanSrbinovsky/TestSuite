#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(UTILIO) || defined(FLDIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! + Gathers any type of field from many processors to one processor
!
! Subroutine interface:

      SUBROUTINE GENERAL_GATHER_FIELD (                                 &
     &  LOCAL_FIELD , GLOBAL_FIELD ,                                    &
     &  LOCAL_SIZE , GLOBAL_SIZE ,                                      &
     &  LEVELS ,                                                        &
     &  ADDR_INFO ,                                                     &
     &  GATHER_PE ,                                                     &
     &  ICODE , CMESSAGE)

      IMPLICIT NONE

! Description:
! Takes a general decomposed field on many processors and gathers it
! to a single processor.
!
! Current code owner : P.Burton
!
! Subroutine arguments:

      INTEGER                                                           &
     &  GLOBAL_SIZE                                                     &
                        ! IN:  size of GLOBAL FIELD
     &, LEVELS                                                          &
                        ! IN:  How many levels of data to do
     &, GATHER_PE(LEVELS)                                               &
                          ! IN: Which PE to gather each level to
     &, LOCAL_SIZE                                                      &
                        ! OUT: size of LOCAL_FIELD
     &, ICODE           ! OUT: return code, 0=OK

! Required for dimensioning ADDR_INFO
#include "d1_addr.h"

      INTEGER                                                           &
     &  ADDR_INFO(D1_LIST_LEN)   ! IN: addressing info about field

      REAL                                                              &
     &  LOCAL_FIELD(*)                                                  &
                                  ! IN: my local part of field
     &, GLOBAL_FIELD(GLOBAL_SIZE,*)
                        ! OUT:  array to gather field to


      CHARACTER*(80)                                                    &
     &  CMESSAGE        ! OUT : Error message

! Parameters and common blocks

#include "cppxref.h"
#include "stparam.h"

#include "gccom.h"

#include "parvars.h"
#include "typsize.h"
#include "atm_lsm.h"
#include "i_stgfld.h"

! Local variables

      INTEGER                                                           &
     &  grid_type                                                       &
                  ! grid type of field being gathered
     &, grid_code                                                       &
                  ! ppx grid code of field being gathered
     &, halo_type                                                       &
                  ! halo type of field being gathered
     &, k                                                               &
           ! loop over levels
     &, my_k                                                            &
             ! value of k for GLOBAL_FIELD on this PE
     &, my_k_temp                                                       &
             ! Temp. my_k for safe references
     &, info                                                            &
             ! return code from GCOM routines
     &, dummy                                                           &
              ! dummy variables - ignored return values
     &, north,south,east,west                                           &
                               ! domain limits for STASH output
     &, mean_type                                                       &
                  ! spatial meaning type on diagnostic
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

! Multi-level fields not supported here
      IF (LEVELS  >   1) THEN
        WRITE(6,*) 'GENERAL_GATHER_FIELD : Cannot have more than ',     &
     &             '1 level for gathering timeseries data.'
        ICODE=1
        CMESSAGE='GENERAL_GATHER_FIELD : Multi-level timeseries field'
        GOTO 9999
      ENDIF

! Copy the data to GLOBAL_FIELD from PE 0

      CALL GC_SSYNC(nproc,info)

      IF (mype  ==  0) THEN

        info=GC_NONE
        CALL GC_RSEND(99,GLOBAL_SIZE,GATHER_PE,info,GLOBAL_FIELD,       &
     &                LOCAL_FIELD)

      ENDIF

      IF (mype  ==  GATHER_PE(1)) THEN

        info=GC_NONE
        CALL GC_RRECV(99,GLOBAL_SIZE,0,info,GLOBAL_FIELD,LOCAL_FIELD)

      ENDIF

      LOCAL_SIZE=GLOBAL_SIZE

!-------------------------------------------------------------------

! Surface (land points only) fields

      ELSEIF                                                            &
     &  (grid_code  ==  ppx_atm_compressed) THEN

        my_k=0
        DO k=1,LEVELS
! Unpack the local field out to full (local) field size and
! put this into the array buf_expand_local

! DEPENDS ON: from_land_points
        CALL from_land_points(buf_expand_local,                         &
     &                        LOCAL_FIELD(1+(k-1)*local_land_field),    &
     &                        atmos_landmask_local,                     &
     &                        lasize(1,grid_type,halo_type)*            &
     &                        lasize(2,grid_type,halo_type),            &
     &                        dummy)

! Now gather in all the processors local fields into the global
! field (array buf_expand)

! DEPENDS ON: gather_field
        CALL GATHER_FIELD(buf_expand_local,buf_expand,                  &
     &                    lasize(1,grid_type,halo_type),                &
     &                    lasize(2,grid_type,halo_type),                &
     &                    glsize(1,grid_type),glsize(2,grid_type),      &
     &                    grid_type,halo_type,                          &
     &                    GATHER_PE(k),GC_ALL_PROC_GROUP,               &
     &                    ICODE,CMESSAGE)

        IF (ICODE  /=  0) THEN
          WRITE(6,*)                                                    &
     &       'GENERAL_GATHER_FIELD : Error detected in call to ',       &
     &       'GATHER_FIELD (land point field)'
          WRITE(6,*) 'Return code : ',ICODE
          WRITE(6,*) 'Error message : ',CMESSAGE

          ICODE=1
          CMESSAGE='GENERAL_GATHER_FIELD : Error gathering field'
          GOTO 9999
        ENDIF

! And now pack the global field (buf_expand) back to land points
! and put into the array GLOBAL_FIELD.

        IF (mype  ==  GATHER_PE(k)) THEN
          my_k=my_k+1
! DEPENDS ON: to_land_points
          CALL to_land_points(buf_expand,GLOBAL_FIELD(1,my_k),          &
     &                        atmos_landmask,                           &
     &                        glsize(1,grid_type)*glsize(2,grid_type),  &
     &                        dummy)
        ENDIF
        ENDDO ! k : loop over levels

        local_size = local_land_field

!-------------------------------------------------------------------

! Atmosphere Lateral boundary fields

      ELSEIF                                                            &
     &  (grid_type  ==  ppx_atm_rim) THEN

! DEPENDS ON: gather_atmos_lbcs
        CALL GATHER_ATMOS_LBCS(GLOBAL_FIELD,global_LENRIMDATA_A,        &
     &                         LOCAL_FIELD,LENRIMDATA_A,                &
     &                         GATHER_PE,                               &
     &                         ICODE,CMESSAGE)

        IF (ICODE  /=  0) THEN
          WRITE(6,*)                                                    &
     &      'GENERAL_GATHER_FIELD : Error detected in call to ',        &
     &      'GATHER_ATMOS_LBCS'
          WRITE(6,*) 'Return code : ',ICODE
          WRITE(6,*) 'Error message : ',CMESSAGE

          ICODE=2
          CMESSAGE='GENERAL_GATHER_FIELD : Error gathering Atmos LBCs'
          GOTO 9999
        ENDIF

        LOCAL_SIZE=LENRIMDATA_A

!-------------------------------------------------------------------

! Ocean Lateral boundary fields

      ELSEIF                                                            &
     &  (grid_type  ==  ppx_ocn_rim) THEN
      WRITE(6,*)'GENERAL_GATHER_FIELD : Error attempting to ',          &
     & 'gather ocean boundary data. Data not stored in dump ',          &
     & 'from version 5.3'
      ICODE=3
      CMESSAGE='GENERAL_GATHER_FIELD : Error gathering Ocean LBCs'
      GO TO 9999

!-------------------------------------------------------------------
! "Normal" fields

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
     &   (grid_code  ==  ppx_atm_river) .OR.                            &
!     ocean grids
     &   (grid_code  ==  ppx_ocn_tall) .OR.                             &
     &   (grid_code  ==  ppx_ocn_tfield) .OR.                           &
     &   (grid_code  ==  ppx_ocn_tzonal) .OR.                           &
     &   (grid_code  ==  ppx_ocn_uzonal) .OR.                           &
     &   (grid_code  ==  ppx_ocn_tmerid) .OR.                           &
     &   (grid_code  ==  ppx_ocn_uall) .OR.                             &
     &   (grid_code  ==  ppx_ocn_cuall) .OR.                            &
     &   (grid_code  ==  ppx_ocn_ufield) .OR.                           &
     &   (grid_code  ==  ppx_ocn_umerid) .OR.                           &
     &   (grid_code  ==  ppx_ocn_cvall) .OR.                            &
!     Full WAM Wave Model Grid - Land Sea Mask only
     &   (grid_code  ==  ppx_wam_all))                                  &
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
            north=south
          ELSEIF (mean_type  >=  4) THEN ! field/global mean
            east=west
            north=south
          ENDIF

        ELSE
          north=glsize(2,grid_type)
          west=1
          east=glsize(1,grid_type)
          south=1
        ENDIF

        my_k=0
        DO k=1,LEVELS

        IF (mype  ==  GATHER_PE(k)) THEN
          my_k      = my_k+1
          my_k_temp = my_k
        ELSE
          my_k_temp = 1
        END IF

!       STASH_GATHER_FIELD can distribute whole fields, or subarea
!       fields

! DEPENDS ON: stash_gather_field
        CALL STASH_GATHER_FIELD(                                        &
     &    LOCAL_FIELD(1+(k-1)*LOCAL_SIZE),GLOBAL_FIELD(1,my_k_temp),    &
     &    LOCAL_SIZE,GLOBAL_SIZE,1,                                     &
     &    north,east,south,west,                                        &
     &    grid_code,halo_type,GATHER_PE(k),.TRUE.,                      &
     &    ICODE=ICODE, CMESSAGE=CMESSAGE)

        IF (ICODE  /=  0) THEN
          WRITE(6,*)                                                    &
     &      'GENERAL_GATHER_FIELD : Error detected in call to ',        &
     &      'STASH_GATHER_FIELD'
          WRITE(6,*) 'Return code : ',ICODE
          WRITE(6,*) 'Error message : ',CMESSAGE

          ICODE=4
          CMESSAGE='GENERAL_GATHER_FIELD : Error gathering field'
          GOTO 9999
        ENDIF
        ENDDO ! k : loop over levels

!-------------------------------------------------------------------
! Any other type of field
      ELSE

        ICODE=10
        CMESSAGE='GENERAL_GATHER_FIELD : Field type not recognized'

      ENDIF

 9999 CONTINUE

      RETURN

      END SUBROUTINE GENERAL_GATHER_FIELD

#endif
