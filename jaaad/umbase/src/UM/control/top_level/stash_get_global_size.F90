#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C)
#if defined(MPP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculates the "global" size of STASHed data.
!
! Subroutine Interface:
      SUBROUTINE STASH_GET_GLOBAL_SIZE(                                 &
     &  GLOBAL_NORTH_IN , GLOBAL_EAST_IN ,                              &
     &  GLOBAL_SOUTH_IN , GLOBAL_WEST_IN ,                              &
     &  LEVELS_IN ,                                                     &
     &  GRIDPOINT_CODE , PROCESSING_CODE ,                              &
     &  GLOBAL_FIELD_SIZE ,                                             &
     &  ICODE , CMESSAGE)

      IMPLICIT NONE

! Description:
! Calculates the global (ie. total size on disk) size of a
! STASH request.
!
! Method:
! Using the PROCESSING_CODE to indicate the type of STASH request,
! the GRIDPOINT_CODE to indicate the grid type of the data,
! and the subdomain limits, the total size of the data is calculated.
! Current code owner : P.Burton
!
! History:
!  Model    Date      Modification history from model version 4.2
!  version
!   4.2     28/10/96  New DECK created for MPP version of STASH
!                     P.Burton
!   4.3     14/3/97   Various fixes                       P.Burton
!   5.0     20/09/99  Correct for South->North ordering      P.Burton
!   6.0    16/10/03  Add def for new NEC opt section C96_1C. R Barnes
!
! Subroutine arguments:

      INTEGER                                                           &
     &  GLOBAL_NORTH_IN                                                 &
                          ! IN: specification of subdomain boundaries
     &, GLOBAL_EAST_IN                                                  &
                          ! IN: ""
     &, GLOBAL_SOUTH_IN                                                 &
                          ! IN: ""
     &, GLOBAL_WEST_IN                                                  &
                          ! IN: ""
     &, LEVELS_IN                                                       &
                          ! IN: number of levels
     &, GRIDPOINT_CODE                                                  &
                          ! IN: indicates the output grid type
     &, PROCESSING_CODE   ! IN: indicates the type of STASH processing

      INTEGER                                                           &
     &  GLOBAL_FIELD_SIZE ! OUT: size of STASH data on disk

      INTEGER                                                           &
     &  ICODE             ! OUT: Return code (0=OK)

      CHARACTER*80                                                      &
     &  CMESSAGE          ! OUT: Error message

! Parameters and common blocks
#include "stparam.h"
#include "parvars.h"


! Local variables

      INTEGER                                                           &
! copies of input arguments, which get modified according the
! type of output grid
     &  global_north,global_east,global_south,global_west               &
     &, levels

! ------------------------------------------------------------------

      global_north = GLOBAL_NORTH_IN
      global_east  = GLOBAL_EAST_IN
      global_south = GLOBAL_SOUTH_IN
      global_west  = GLOBAL_WEST_IN
      levels       = LEVELS_IN

! Fix wrap-arounds s.t. east > west

      IF (global_west  >   global_east)                                 &
     &  global_east=global_east+glsize(1,fld_type_p)

! Full field or subdomain output:

      IF ((PROCESSING_CODE  ==  st_replace_code) .OR.                   &
     &    (PROCESSING_CODE  ==  st_accum_code) .OR.                     &
     &    (PROCESSING_CODE  ==  st_time_mean_code) .OR.                 &
     &    (PROCESSING_CODE  ==  st_max_code) .OR.                       &
     &    (PROCESSING_CODE  ==  st_min_code)) THEN

        IF ((GRIDPOINT_CODE  >=  vert_mean_base) .AND.                  &
                                                         ! vertical
     &      (GRIDPOINT_CODE  <   zonal_mean_base)) THEN ! mean
          levels=1

        ELSEIF ((GRIDPOINT_CODE  >=  zonal_mean_base) .AND.             &
                                                             ! zonal
     &          (GRIDPOINT_CODE  <   merid_mean_base)) THEN  ! mean
          global_east=global_west

        ELSEIF ((GRIDPOINT_CODE  >=  merid_mean_base) .AND.             &
                                                            ! merid.
     &          (GRIDPOINT_CODE  <   field_mean_base)) THEN ! mean
          global_south=global_north

        ELSEIF ((GRIDPOINT_CODE  >=  field_mean_base)  .AND.            &
                                                             ! field
     &          (GRIDPOINT_CODE  <   global_mean_base)) THEN ! fmean
          global_east=global_west
          global_south=global_north

        ELSEIF (GRIDPOINT_CODE  >=  global_mean_base) THEN
          global_east=global_west
          global_south=global_north
          levels=1

        ELSEIF (GRIDPOINT_CODE  >   global_mean_top) THEN
          ICODE=1
          WRITE(6,*) 'Grid type ',GRIDPOINT_CODE,                       &
     &               ' not yet supported by MPP STASH'
          CMESSAGE='Unsupported grid type'
          GOTO 9999
        ENDIF

        GLOBAL_FIELD_SIZE=(global_east-global_west+1)*                  &
     &                    (global_north-global_south+1)*                &
     &                    levels

      ELSE
        ICODE=2
        WRITE(6,*) 'Processing code ',PROCESSING_CODE,                  &
     &             ' not yet supported by MPP STASH'
        CMESSAGE='Unsupported processing code'
        GOTO 9999
      ENDIF

 9999 CONTINUE

      RETURN

      END SUBROUTINE STASH_GET_GLOBAL_SIZE

#endif
#endif
