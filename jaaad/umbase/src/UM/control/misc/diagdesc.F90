#if defined(C70_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine DIAGDESC -----------------------------------------------
!LL
!LL  Purpose: Prints a formatted diagnostic description using the name
!LL           of a diagnostic plus it's PPXREF and STASH record.  Gives
!LL           a hardcopy record of the diagnostics included in a run.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 6.1
!LL
!LL  Author:   T.Johns            Date:           14 January 1992
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.1  05/02/93  Correct minor bug in printout for climate mean tag.
!LL                  Print out pseudo-level information.
!LL  3.1   3/02/93 : added comdeck CHSUNITS to define NUNITS for i/o.
!LL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
!LL                  1.Removes the limit on primary STASH item numbers.
!LL                  2.Removes the assumption that (section,item)
!LL                    defines the sub-model.
!LL                  3.Thus allows for user-prognostics.
!LL   3.4  13/01/94  Replace hardwired gridcodes by ppx_ parameters, and
!LL                  cover all options.   T. Johns
!     4.4  02/12/96 Add daily mean timeseries R. A. Stratton.
!LL   5.2  22/08/00 Change format for extra pp units. R Rawlins
!     5.2  06/11/00 Change row listing from N->S to S->N. R Rawlins
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: C401
!LL
!LL  Project task: C4
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage handling and diagnostic
!LL                                 system (STASH)
!LLEND --------------------------------------------------------------
!
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE DIAGDESC(seqno,name,stlist,ppxref,                     &
     &           stash_levels,num_stash_levels,num_level_lists,         &
     &           stash_pseudo_levels,num_stash_pseudo,num_pseudo_lists, &
     &           sttabl,nsttims,nsttabl,                                &
     &           stash_series,stash_series_rec_len,stash_series_len,    &
     &           stash_series_index,stash_ser_index_size)
!
      IMPLICIT NONE
!
      CHARACTER*36                                                      &
     &    name                                  ! IN  diagnostic name
      INTEGER                                                           &
     &    seqno,                                                        &
                                                ! IN  sequence number
     &    stlist(*),                                                    &
                                                ! IN  STASHlist record
     &    ppxref(*)                             ! IN  PPXREF record
!
! STASH levels list information
      INTEGER                                                           &
     &       num_stash_levels                                           &
                                              ! IN Max levels in a list
     &,      num_level_lists                                            &
                                              ! IN Number of lists
     &,      stash_levels(num_stash_levels+1,num_level_lists)
! STASH pseudo-levels list information
      INTEGER                                                           &
     &       num_stash_pseudo                                           &
                                              ! IN Max ps-levs in a list
     &,      num_pseudo_lists                                           &
                                              ! IN No of ps-lev lists
     &,      stash_pseudo_levels(num_stash_pseudo+1,num_pseudo_lists)
! STASH time list information
      INTEGER                                                           &
     &       nsttims                                                    &
                                              ! IN Max times in a list
     &,      nsttabl                                                    &
                                              ! IN Number of lists
     &,      sttabl(nsttims,nsttabl)
! STASH timeseries information
      INTEGER                                                           &
     &       stash_series_len                                           &
                                              ! IN Total no of records
     &,      stash_series_rec_len                                       &
                                              ! IN Length of each record
     &,      stash_series(stash_series_rec_len,stash_series_len)        &
!                                             ! IN array of records
     &,      stash_ser_index_size                                       &
                                              ! IN No of index records
     &,      stash_series_index(2,stash_ser_index_size)
!
#include "stparam.h"
#include "cppxref.h"
#include "csubmodl.h"
#include "chsunits.h"
#include "ccontrol.h"
!
! Local variables
!
      CHARACTER*132 line,line1,line2 ! Encoded line of information
      CHARACTER*80  ch            ! Working character string variable
      INTEGER i1,i2,k             ! Array indices
      INTEGER j                   ! Code value
      INTEGER time_list,lev_list                                        &
                                  ! pointers to time and levels lists
     &,       plev_list                                                 &
                                  ! pointer  to pseudo-level list
     &,       tser_list           ! pointer  to time series record list
      INTEGER ntimes              ! no of times in a time list
      INTEGER packing_profile     ! packing profile for output PPfield
!
!L----------------------------------------------------------------------
!L 0. Write header if sequence no indicates first item
!L
      IF (seqno == 1) THEN
        WRITE(6,1000)
 1000   FORMAT(                                                         &
     &  '   ********************************************************'/  &
     &  '   ********************************************************'/  &
     &  '   **                                                    **'/  &
     &  '   **    LIST OF USER-DEFINED DIAGNOSTICS IN THIS RUN    **'/  &
     &  '   **                                                    **'/  &
     &  '   ********************************************************'/  &
     &  '   ********************************************************'/  &
     &  '   **                                                    **'/  &
     &  '   ** NOTES:                                             **'/  &
     &  '   **   Time processing details are in timesteps, where  **'/  &
     &  '   **     ... represents "for ever".                     **'/  &
     &  '   **   Spatial processing domain is in gridpoints.      **'/  &
     &  '   **                                                    **'/  &
     &  '   ********************************************************'/  &
     &  '   ********************************************************'// &
     &'=================================================================&
     &==========================================================')
      ENDIF
!L----------------------------------------------------------------------
!L 1. For each diagnostic processing request in the STASHlist,
!L    print the diagnostic name followed by a summary of the processing
!L    information on 3 lines.
!L
!L 1.0 If diagnostic is not required for output, exit routine
!L
      IF (stlist(st_proc_no_code) == 0) GOTO 999
!L
!L 1.1 Line 1.
!L
      line=' '
! #No
      i1=2
      i2=4
      write(ch,'(i3)') seqno
      line(i1:i2)=ch(1:1+i2-i1)
! Name
      i1=i2+2
      i2=i1+36-1
      line(i1:i2)=name
! Submodel
      i1=i2+2
      i2=i1+8-1
      j=stlist(st_sect_no_code)
      IF      (ppxref(ppx_model_number) == ocean_im) THEN
        ch=' OCEAN  '
      ELSE IF (ppxref(ppx_model_number) ==  slab_im) THEN
        ch=' SLAB   '
      ELSE IF (ppxref(ppx_model_number) == atmos_im) THEN
        ch=' ATMOS  '
      ELSE IF (ppxref(ppx_model_number) == wave_im) THEN
        ch=' WAVE   '
      ELSE
        WRITE(6,*)' Error in DIAGDES. Unknown model'
        ch=' UNKNOWN'
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
! Item
      i1=i2+2
      i2=i1+4-1
      j=stlist(st_item_code)
      write(ch,'(i4)') j
      line(i1:i2)=ch(1:1+i2-i1)
! Section
      i1=i2+2
      i2=i1+7-1
      j=stlist(st_sect_no_code)
      write(ch,'(i7)') j
      line(i1:i2)=ch(1:1+i2-i1)
! PPfcode
      i1=i2+2
      i2=i1+7-1
      j=ppxref(ppx_field_code)
      write(ch,'(i7)') j
      line(i1:i2)=ch(1:1+i2-i1)
! Datatype
      i1=i2+2
      i2=i1+8-1
      j=ppxref(ppx_data_type)
      IF (j == 1.OR.j == 4) THEN
        ch='  REAL  '
      ELSEIF (j == 2.OR.j == 5) THEN
        ch='INTEGER '
      ELSEIF (j == 3) THEN
        ch='LOGICAL '
      ELSE
        ch='UNKNOWN '
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
! Gridtype
      i1=i2+2
      i2=i1+8-1
      j=ppxref(ppx_grid_type)
      IF (j == ppx_atm_nonstd.OR.j == ppx_ocn_nonstd) THEN
        ch=' NONSTD '
      ELSEIF ((j >  ppx_atm_nonstd.AND.j <= ppx_atm_tsea) .OR.          &
     &         j == ppx_atm_compressed.OR.j == ppx_atm_ozone) THEN
        ch=' P-GRID '
      ELSEIF (j >= ppx_atm_uall.AND.j <= ppx_atm_usea) THEN
        ch=' UV-GRID'
      ELSEIF (j == ppx_atm_cuall.OR.j == ppx_ocn_cuall) THEN
        ch=' CU-GRID'
      ELSEIF (j == ppx_atm_cvall.OR.j == ppx_ocn_cvall) THEN
        ch=' CV-GRID'
      ELSEIF (j == ppx_atm_tzonal) THEN
        ch=' PZ-GRID'
      ELSEIF (j == ppx_atm_uzonal) THEN
        ch=' UZ-GRID'
      ELSEIF (j == ppx_atm_tmerid) THEN
        ch=' PM-GRID'
      ELSEIF (j == ppx_atm_umerid) THEN
        ch=' UM-GRID'
      ELSEIF (j == ppx_atm_rim.OR.j == ppx_ocn_rim) THEN
        ch='   RIM  '
      ELSEIF (j == ppx_ocn_tcomp.OR.j == ppx_ocn_tall.OR.               &
     &        j == ppx_ocn_tfield) THEN
        ch=' T-GRID '
      ELSEIF (j == ppx_ocn_tzonal) THEN
        ch=' TZ-GRID'
      ELSEIF (j == ppx_ocn_uzonal) THEN
        ch=' UZ-GRID'
      ELSEIF (j == ppx_ocn_tmerid) THEN
        ch=' TM-GRID'
      ELSEIF (j == ppx_ocn_umerid) THEN
        ch=' UM-GRID'
      ELSEIF (j == ppx_ocn_ucomp.OR.j == ppx_ocn_uall.OR.               &
     &        j == ppx_ocn_ufield) THEN
        ch=' UV-GRID'
      ELSEIF (j == ppx_atm_scalar.OR.j == ppx_ocn_scalar) THEN
        ch=' SCALAR '
      ELSEIF (j == ppx_wam_all.OR.j == ppx_wam_sea) THEN
        ch=' WAVE   '
      ELSEIF (j == ppx_wam_rim) THEN
        ch=' RIM    '
      ELSE
        ch=' UNKNOWN'
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
! Leveltype
      i1=i2+2
      i2=i1+9-1
      j=ppxref(ppx_lv_code)
      IF (j == ppx_full_level) THEN
        ch='FULLLEVEL'
      ELSEIF (j == ppx_half_level) THEN
        ch='HALFLEVEL'
      ELSE
        ch='STD-LEVEL'
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
! Meto8LV
      i1=i2+2
      i2=i1+7-1
      j=ppxref(ppx_meto8_levelcode)
      write(ch,'(i7)') j
      line(i1:i2)=ch(1:1+i2-i1)
! Meto8FC
      i1=i2+2
      i2=i1+7-1
      j=ppxref(ppx_meto8_fieldcode)
      write(ch,'(i7)') j
      line(i1:i2)=ch(1:1+i2-i1)
! PackAcc
      i1=i2+2
      i2=i1+7-1
      j=stlist(st_output_code)
      IF (j == 1) THEN
        IF (stlist(st_macrotag) >= 1000) THEN
          packing_profile=pp_pack_code(27)
        ELSE
          packing_profile=0
        ENDIF
      ELSEIF(j == 2) THEN
        packing_profile=0
      ELSEIF(j <  0) THEN
        packing_profile=pp_pack_code(-j)
      ELSE
        packing_profile=0
      ENDIF
      IF (packing_profile == 0) THEN
        ch='       '
      ELSE
        j=ppxref(ppx_packing_acc+packing_profile-1)
        write(ch,'(i7)') j
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
!
      line1=line
!L
!L 1.2 Line 2.
!L
      line=' '
! Time-processing
      i1=2
      i2=16
      j=stlist(st_proc_no_code)
      tser_list=0
      IF (j == st_replace_code) THEN
        ch='   EXTRACT     '
      ELSEIF (j == st_accum_code) THEN
        ch=' ACCUMULATION  '
      ELSEIF (j == st_time_mean_code) THEN
        ch='  TIME MEAN    '
      ELSEIF (j == st_time_series_code) THEN
        write(ch,'(''  TIME SERIES  '')')
        tser_list=stlist(st_series_ptr)
      ELSEIF (j == st_max_code) THEN
        ch='MAX OVER PERIOD'
      ELSEIF (j == st_min_code) THEN
        ch='MIN OVER PERIOD'
      ELSEIF (j == st_append_traj_code) THEN
        ch='  TRAJECTORY   '
      ELSEIF (j == st_variance_code) THEN
        ch=' TIME VARIANCE '
      ELSEIF (j == st_time_series_mean) THEN
        ch='MEAN TIMESERIES'
      ELSE
        ch='  UNKNOWN      '
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
! -From-
      i1=i2+2
      i2=i1+6-1
      IF (stlist(st_freq_code) <  0) THEN
        ch='      '
      ELSE
        j=stlist(st_start_time_code)
        write(ch,'(i6)') j
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
! --To--
      i1=i2+2
      i2=i1+6-1
      IF (stlist(st_freq_code) <  0) THEN
        ch='      '
      ELSE
        j=stlist(st_end_time_code)
        IF (j == st_infinite_time) THEN
          ch='  ... '
        ELSE
          write(ch,'(i6)') j
        ENDIF
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
! Frequency
      i1=i2+2
      i2=i1+9-1
      j=stlist(st_freq_code)
      ! Extra flag for dump frequency
      IF (j <  0 .and. j /= -9999) THEN
        j=-j
        write(ch,'(''TIME LIST'')')
        time_list=j
      ELSE
        write(ch,'(i9)') j
        time_list=0
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
! Period
      i1=i2+2
      i2=i1+6-1
      IF (stlist(st_freq_code) <  0) THEN
        ch='      '
      ELSE
        j=stlist(st_period_code)
        IF (stlist(st_proc_no_code) == st_replace_code) THEN
          ch='      '
        ELSEIF (j == st_infinite_time) THEN
          ch='  ... '
        ELSE
          write(ch,'(i6)') j
        ENDIF
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
! __Source__
      i1=i2+2
      i2=i1+10-1
      j=stlist(st_input_code)
      IF (j == 0) THEN
        ch='PROGNOSTIC'
      ELSEIF(j == 1) THEN
        ch='  STWORK  '
      ELSEIF(j <  0) THEN
        j=-j
        write(ch,'(''DUMP #'',i4)') j
      ELSE
        ch=' UNKNOWN  '
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
! ___Destination___
      i1=i2+2
      i2=i1+17-1
      j=stlist(st_output_code)
      IF (j == 1) THEN
        IF (stlist(st_macrotag) >= 1000) THEN
          ch='MEAN PP VIA DUMP'
        ELSEIF (stlist(st_macrotag) >  0) THEN
          write(ch,'(''DUMP WITH TAG '',i3)') stlist(st_macrotag)
        ELSE
          ch='      DUMP       '
        ENDIF
      ELSEIF(j == 2) THEN
        ch='   SECONDARY     '
      ELSEIF(j <  0) THEN
        j=-j
        IF (j == 27) THEN
          ch='MEAN PP (DIRECT) '
        ELSE
          write(ch,'(''   PP UNIT'',i3)') j
        ENDIF
      ELSE
        ch='  UNKNOWN  '
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
!
      line2=line
!L
!L 1.3 Line 3.
!L
      line=' '
! Spatial-Processing
      i1=2
      i2=19
      j=stlist(st_gridpoint_code)
      IF (j >= extract_base.AND.j <  extract_top) THEN
        ch='    FULL FIELD    '
      ELSEIF (j >= vert_mean_base.AND.j <  vert_mean_top) THEN
        ch='  VERTICAL MEAN   '
      ELSEIF (j >= zonal_mean_base.AND.j <  zonal_mean_top) THEN
        ch='   ZONAL MEAN     '
      ELSEIF (j >= merid_mean_base.AND.j <  merid_mean_top) THEN
        ch=' MERIDIONAL MEAN  '
      ELSEIF (j >= field_mean_base.AND.j <  field_mean_top) THEN
        ch=' FIELD MEAN - 2D  '
      ELSEIF (j >= global_mean_base.AND.j <  global_mean_top) THEN
        ch=' GLOBAL MEAN - 3D '
      ELSE
        ch='  ** UNKNOWN **   '
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
! Levels-domain
      i1=i2+2
      i2=i1+13-1
      j=stlist(st_output_bottom)
      lev_list=0
      IF (j == st_special_code) THEN
        ch='STANDARD LEV '
      ELSEIF (j >  0) THEN
        write(ch,'(''LEVELS '',i2,''-'',i2)') j,stlist(st_output_top)
      ELSEIF (j <  0) THEN
        j=-j
        write(ch,'('' LEVELS LIST '')')
        lev_list=j
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
! Pseudo-levels
      i1=i2+2
      i2=i1+15-1
      j=stlist(st_pseudo_out)
      plev_list=0
      IF (j >  0) THEN
        ch='PSEUDO-LEV LIST'
        plev_list=j
      ELSE
        ch='     NONE      '
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
! Horizontal-domain.....
      i1=i2+2
      i2=i1+23-1
      write(ch,'(''ROW:'',i3,''-'',i3,'' COL:'',i3,''-'',i3)')          &
     &  stlist(st_south_code),stlist(st_north_code),                    &
     &  stlist(st_west_code),stlist(st_east_code)
      line(i1:i2)=ch(1:1+i2-i1)
! Weighting
      i1=i2+2
      i2=i1+9-1
      j=stlist(st_weight_code)
      IF (j == stash_weight_null_code) THEN
        ch='  NONE   '
      ELSEIF (j == stash_weight_area_code) THEN
        ch='  AREA   '
      ELSEIF (j == stash_weight_volume_code) THEN
        ch=' VOLUME  '
      ELSEIF (j == stash_weight_mass_code) THEN
        ch='  MASS   '
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
! Masking
      i1=i2+2
      i2=i1+7-1
      j=mod(stlist(st_gridpoint_code),block_size)
      IF (j == stash_null_mask_code) THEN
        ch=' NONE  '
      ELSEIF (j == stash_land_mask_code) THEN
        ch=' LAND  '
      ELSEIF (j == stash_sea_mask_code) THEN
        ch='  SEA  '
      ELSE
        ch='UNKNOWN'
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
!L
!L 1.4 Print the main part of the summary
!L
      WRITE(6,1010) line1,line2,line
 1010 FORMAT(' #No ',                                                   &
     &      'Diagnostic Description-------------- Submodel Item Section &
     &PPfcode Datatype Gridtype Leveltype MetO8lv MetO8fc Packacc'/     &
     &  a124/                                                           &
     &' Time-processing -From- --To-- Frequency Period --Source-- ---Des&
     &tination---                                               '/      &
     &  a124/                                                           &
     &' Spatial-processing Levels-domain -Pseudo-levels- ---Horizontal-d&
     &omain--- Weighting Masking                                '/      &
     &  a124)
!L
!L 1.5 Print associated time and levels lists if appropriate
!L
!L 1.5.1 Time list
!L
      IF (time_list /= 0) THEN
        DO j=1,nsttims
          IF (sttabl(j,time_list) == st_end_of_list) THEN
            ntimes=j-1
            GOTO 210
          ENDIF
        ENDDO
  210   CONTINUE
        WRITE(6,'('' ***** TIME LIST ***** '',i3,                       &
     &            '' times are as follows:-'')') ntimes
        i1=1
        i2=8
        DO j=1,ntimes
          IF (i1 == 1) line=' '
          WRITE(ch,'(1x,i7)') sttabl(j,time_list)
          line(i1:i2)=ch(1:8)
          i1=i1+8
          i2=i2+8
          IF (i2 >  80) THEN
            i1=1
            i2=8
            WRITE(6,'(a80)') line
          ENDIF
        ENDDO
        IF (i2 <= 80) WRITE(6,'(a80)') line
      ENDIF
!L
!L 1.5.2 Levels list
!L
      IF (lev_list /= 0) THEN
        write(6,'('' ***** LEVELS LIST ***** '',i3,                     &
     &          '' levels are as follows:-'')') stash_levels(1,lev_list)
        i1=1
        i2=8
        DO j=2,1+stash_levels(1,lev_list)
          IF (i1 == 1) line=' '
          write(ch,'(1x,i7)') stash_levels(j,lev_list)
          line(i1:i2)=ch(1:8)
          i1=i1+8
          i2=i2+8
          IF (i2 >  80) THEN
            i1=1
            i2=8
            write(6,'(a80)') line
          ENDIF
        ENDDO
        IF (i2 <= 80) write(6,'(a80)') line
      ENDIF
!L
!L 1.5.3 Pseudo-levels list
!L
      IF (plev_list /= 0) THEN
        write(6,'('' ***** PSEUDO-LEVELS LIST ***** '',i3,              &
     &          '' pseudo-levels are as follows:-'')')                  &
     &          stash_pseudo_levels(1,plev_list)
        i1=1
        i2=8
        DO j=2,1+stash_pseudo_levels(1,plev_list)
          IF (i1 == 1) line=' '
          write(ch,'(1x,i7)') stash_pseudo_levels(j,plev_list)
          line(i1:i2)=ch(1:8)
          i1=i1+8
          i2=i2+8
          IF (i2 >  80) THEN
            i1=1
            i2=8
            write(6,'(a80)') line
          ENDIF
        ENDDO
        IF (i2 <= 80) write(6,'(a80)') line
      ENDIF
!L
!L 1.5.4 Time series subdomain record list
!L
      IF (tser_list /= 0) THEN
        i1=stash_series_index(1,tser_list)
        i2=stash_series_index(2,tser_list)
        WRITE(6,'('' ***** TIME SERIES ***** '',i3,                     &
     & '' subdomain records are as follows:-''/                         &
     & '' Record      North/South       West/ East     Bottom/  Top'')')&
     &    i2
        DO j=1,i2
          WRITE(6,'(3x,i4,1x,3(5x,i5,1x,i5,1x))')                       &
     &        j,(stash_series(3+k,i1+j-1),k=1,6)
        ENDDO
      ENDIF
!L
!L 1.5.5 Print final ruler line
!L
      WRITE(6,1020)
 1020 FORMAT(                                                           &
     &'=================================================================&
     &==========================================================')
!
 999  CONTINUE
      RETURN
      END SUBROUTINE DIAGDESC
#endif
