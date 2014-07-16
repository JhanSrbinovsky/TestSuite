
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: TEMPORAL -------------------------------------------------
!LL
!LL  Purpose: Control routine to handle temporal processing options
!LL           within STASH.  Its input and output arguments look like
!LL           1D arrays (ie. all the data should be in contiguous areas
!LL           of memory).  Lower level service routines are called to
!LL           perform the individual processing options.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Author:   S.Tett
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.1  24/02/93  Change name of variable 'end' to 'last_ts' (ST).
!     4.4  25/11/96  Add processing code option 8 - daily mean
!                    timeseries. R A Stratton.
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: D72
!LL
!LL  Project task: D7
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage handling and diagnostic
!LL                                 system (STASH)
!LL
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE TEMPORAL(variable,result,size,extra_size,              &
     &  control,control_size,ocean,                                     &
     &  timestep,error,errmssg,start,amdi)
!
      IMPLICIT NONE
!
      INTEGER size                  ! IN  size of arrays
      REAL variable(size)           ! IN  data array
      REAL result(size)             ! OUT output array
      INTEGER extra_size            ! IN size of extra data
      INTEGER control_size          ! IN  size of control
      INTEGER control(control_size) ! IN  control
      INTEGER timestep              ! IN  present value of timestep
      INTEGER error                 ! OUT error code
      CHARACTER*(*) errmssg         ! OUT error message
      REAL amdi                     ! IN  missing data indicator
      LOGICAL ocean                 ! IN  true if ocean diagnostic
      LOGICAL start                 ! OUT true if start timestep
!*----------------------------------------------------------------------
!LL  Comdeck: STERR ----------------------------------------------------
!LL
!LL  Purpose: PARAMETER names for STASH processing error codes;
!LL           fatal errors have positive codes, warnings negative.
!LL
!LL  Author:   S.Tett
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.3  16/09/93  Add st_illegal_weight error code.
!LL                   Added st_no_data for mpp code
!LL                   (means a processor does not contain any data
!LL                    for a given subdomain)                 P.Burton
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: D70
!LL
!LL  Project task: D7
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage handling and diagnostic
!LL                                 system (STASH)
!
! Warning codes
!
         integer st_upper_less_lower ! warning code for bad domain
         parameter(st_upper_less_lower=-1)

         integer st_not_supported ! warning code for unsupported routine
         parameter(st_not_supported=-2)
         integer st_no_data,st_nd ! indicates no data on a processor
         parameter(st_no_data=-3,st_nd=-3)
!
! Error codes
!
         integer st_bad_array_param ! error code for dodgy array params
         parameter(st_bad_array_param=1)

         integer st_bad_address     ! error code for address violation
         parameter(st_bad_address=2)

         integer st_unknown ! error code for unknown option
         parameter(st_unknown=3)

         integer st_bad_wraparound ! error code for illegal wraparound
         parameter(st_bad_wraparound=4)

         integer st_illegal_weight ! error code for illegal weighting
         parameter(st_illegal_weight=9)

         integer unknown_weight ! error code for an unknown weight
         parameter(unknown_weight=10)

         integer unknown_mask ! error code for an unknown mask
         parameter(unknown_mask=11)

         integer unknown_processing ! error code for unknown processing
         parameter(unknown_processing=12)

         integer nonsense ! error code for general nonsense request
         parameter(nonsense=13)
! STPARAM
!
!  Purpose: Meaningful PARAMETER names for STASH processing routines.
!           Both a long name and short name have been declared, to
!           reduce the existence of "magic" numbers in STASH.
!           Format is that first the address of the item is declare in
!           both long and short form. example is;
!             integer st_item_code,s_item  !Item number (declaration)
!             parameter(st_item_code=3,s_item=3)
!
!  Author:   S.Tett             Date:           22 January 1991
!
!  Model            Modification history from model version 3.0:
! version  Date
!   3.5    Mar. 95  Sub-models project.
!                   st_model_code=28 added to STLIST addresses
!                                   S.J.Swarbrick
!   4.2    27/11/96 mpp code: Added new stlist "magic numbers" :
!                   st_dump_output_length, st_dump_output_addr
!                                                       P.Burton
!   4.4    23/09/97 Add st_offset_code to the STASH list
!                   S.D. Mullerworth
!    4.4  02/12/96 Time mean timeseries added R A Stratton.
!    4.5  23/01/98 Added new stlist magic number
!                  st_dump_level_output_length
!    4.5  23/01/98 A
!    5.5  28/02/03 Original Modifications for WAM. M.Holt
!         06/08/00 Modification for parallelisation of WAM.
!                          Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!    6.0  08/09/03 Add st_riv_grid 23. C.Bunton
!
!  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
!  Logical components covered: D70
!
!  Project task: D7
!
!  External documentation:
!    Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                 system (STASH)
!--------------------------------------------------------------

      ! Internal model number address
      INTEGER,PARAMETER:: st_model_code = 28
      INTEGER,PARAMETER:: s_modl        = 28

      ! Section Number address
      INTEGER,PARAMETER:: st_sect_no_code = 2
      INTEGER,PARAMETER:: s_sect          = 2
      INTEGER,PARAMETER:: st_sect_code    = 2

      INTEGER,PARAMETER:: st_item_code=1,s_item=1 ! Item number address

      ! Processing Code address
      INTEGER,PARAMETER:: st_proc_no_code=3,s_proc=3

      ! subsidiary codes for st_proc_no_code now
      INTEGER,PARAMETER:: st_replace_code=1
      INTEGER,PARAMETER:: st_accum_code=2
      INTEGER,PARAMETER:: st_time_mean_code=3
      INTEGER,PARAMETER:: st_time_series_code=4
      INTEGER,PARAMETER:: st_max_code=5
      INTEGER,PARAMETER:: st_min_code=6
      INTEGER,PARAMETER:: st_append_traj_code=7
      INTEGER,PARAMETER:: st_time_series_mean=8
      INTEGER,PARAMETER:: st_variance_code=9

      ! Frequency (Input & output) addres
      INTEGER,PARAMETER:: st_freq_code=4,s_freq=4

      ! Offset for sampling
      INTEGER,PARAMETER:: st_offset_code=30,s_offs=30

      ! start timestep address
      INTEGER,PARAMETER:: st_start_time_code=5,s_times=5

      ! end timestep address
      INTEGER,PARAMETER:: st_end_time_code=6,s_timee=6

      ! period in timesteps address
      INTEGER,PARAMETER:: st_period_code=7,s_period=7

      ! infinite end/period value
      INTEGER,PARAMETER:: st_infinite_time=-1

      INTEGER,PARAMETER:: st_end_of_list=-1 !end-of-list marker in times

      ! grid point stuff
      ! gridpoint info address
      INTEGER,PARAMETER:: st_gridpoint_code=8,s_grid=8

      ! now subsid grid point stuff
      ! no masking done
      INTEGER,PARAMETER:: stash_null_mask_code=1,s_nomask=1

      ! land mask conds
      INTEGER,PARAMETER:: stash_land_mask_code=2,s_lndms=2

      ! sea mask code
      INTEGER,PARAMETER:: stash_sea_mask_code=3,s_seams =3

      ! processing options

      ! size of block for gridpoint code
      INTEGER,PARAMETER:: block_size=10

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: extract_base=block_size*0

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: extract_top=block_size*1

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_base=block_size*1

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_top=block_size*2

      ! max code for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_base=block_size*2

      ! base codes for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_top=block_size*3

      ! max code for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_base=block_size*3

      ! base codes for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_top=block_size*4

      ! max code for field mean subroutine
      INTEGER,PARAMETER:: field_mean_base=block_size*4

      ! base codes for field mean subroutine
      INTEGER,PARAMETER:: field_mean_top=block_size*5

      ! max code for global mean subroutine
      INTEGER,PARAMETER:: global_mean_base=block_size*5

      ! base codes for global mean subroutine
      INTEGER,PARAMETER:: global_mean_top=block_size*6

      ! Weighting

      ! weighting info address
      INTEGER,PARAMETER:: st_weight_code=9,s_weight=9

      INTEGER,PARAMETER:: stash_weight_null_code  =0,s_noweight  =0
      INTEGER,PARAMETER:: stash_weight_area_code  =1,s_areaweight=1
      INTEGER,PARAMETER:: stash_weight_volume_code=2,s_volweight =2
      INTEGER,PARAMETER:: stash_weight_mass_code  =3,s_massweight=3

      ! Domain definition

      ! row addresses
      INTEGER,PARAMETER:: st_north_code=12,s_north=12
      INTEGER,PARAMETER:: st_south_code=13,s_south=13
      INTEGER,PARAMETER:: st_west_code =14,s_west =14
      INTEGER,PARAMETER:: st_east_code =15,s_east =15

      ! Levels

      ! input bottom level address
      INTEGER,PARAMETER:: st_input_bottom=10,s_bottom =10

      ! special code
      INTEGER,PARAMETER:: st_special_code=100,s_special=100

      ! input top level address
      INTEGER,PARAMETER:: st_input_top=11,s_top=11

      ! output bottom level address
      INTEGER,PARAMETER:: st_output_bottom=21,s_outbot=21

      ! output top level address
      INTEGER,PARAMETER:: st_output_top=22,s_outtop=22

      INTEGER,PARAMETER:: st_model_level_code=1,s_model=1

      ! code for pressure leve
      INTEGER,PARAMETER:: st_pressure_level_code=2,s_press=2

      ! code for height levels
      INTEGER,PARAMETER:: st_height_level_code=3,s_height=3

      ! input code addres
      INTEGER,PARAMETER:: st_input_code=16,s_input=16

      ! input length of diagnostic address
      INTEGER,PARAMETER:: st_input_length=17,s_length=17

      ! output code address
      INTEGER,PARAMETER:: st_output_code=18,s_output=18

      ! Pointer to D1 addressing information
      ! Pos of item in D1 for relevant submodel
      INTEGER,PARAMETER:: st_position_in_d1=29,st_d1pos=29

      ! Output destination options

      INTEGER,PARAMETER:: st_dump=1
      INTEGER,PARAMETER:: st_secondary=2

      ! output length of diagnostic address
      INTEGER,PARAMETER:: st_output_length=19,s_outlen=19
         integer st_dump_output_length,s_doutlen ! output length on
         parameter(st_dump_output_length=32,s_doutlen=32)  ! dump
         integer st_dump_level_output_length,s_dlevoutlen
         parameter(st_dump_level_output_length=33,s_dlevoutlen=33)
! output length of a single level on dump

         integer st_output_addr,s_outadd ! start locn of diag after stas
         parameter(st_output_addr=20,s_outadd=20)       ! output address
         integer st_dump_output_addr,s_doutadd ! output address on
         parameter(st_dump_output_addr=31,s_doutadd=31)  ! dump

      ! ptr to dump lookup header address
      INTEGER,PARAMETER:: st_lookup_ptr=23

      ! ptr into stash_series where control data address
      INTEGER,PARAMETER:: st_series_ptr=24

      ! subsid stuff for time series
      INTEGER,PARAMETER:: series_grid_type=1
      INTEGER,PARAMETER:: series_grid_code=0
      INTEGER,PARAMETER:: series_long_code=1
      INTEGER,PARAMETER:: series_size=2
      INTEGER,PARAMETER:: series_proc_code=3
      INTEGER,PARAMETER:: series_north=4
      INTEGER,PARAMETER:: series_south=5
      INTEGER,PARAMETER:: series_west=6
      INTEGER,PARAMETER:: series_east=7
      INTEGER,PARAMETER:: series_list_start=8
      INTEGER,PARAMETER:: series_list_end=9
      INTEGER,PARAMETER:: record_size=9

      ! Miscellaneous parameters

      ! system/user tag field in stlist address
      INTEGER,PARAMETER:: st_macrotag=25

      ! Pseudo-level list pointers

      ! pseudo-levels input list address
      INTEGER,PARAMETER:: st_pseudo_in=26

      ! pseudo-levels output list address
      INTEGER,PARAMETER:: st_pseudo_out=27

      ! Internal horizontal gridtype codes common to all diagnostics

      INTEGER,PARAMETER:: st_tp_grid =1 ! T-p grid
      INTEGER,PARAMETER:: st_uv_grid =2 ! u-v grid
      INTEGER,PARAMETER:: st_cu_grid =3 ! C-grid u point
      INTEGER,PARAMETER:: st_cv_grid =4 ! C-grid v point
      INTEGER,PARAMETER:: st_zt_grid =5 ! Zonal T-grid
      INTEGER,PARAMETER:: st_zu_grid =6 ! Zonal u-grid
      INTEGER,PARAMETER:: st_mt_grid =7 ! Meridional T-grid
      INTEGER,PARAMETER:: st_mu_grid =8 ! Meridional u-grid
      INTEGER,PARAMETER:: st_riv_grid= 23    ! river_routing grid
      INTEGER,PARAMETER:: st_scalar  =9 ! Scalar (ie. single value)
      INTEGER,PARAMETER:: st_wam_all= 60    ! Wam Field on Full Grid
      INTEGER,PARAMETER:: st_wam_sea= 62    ! Wam Field on Sea Points

! STPARAM end
!
! Subroutines called
!
      EXTERNAL staccum,stmax,stmin
!
! Local variables
!
      LOGICAL masking        ! indicator for masking (ie. missing data)
      INTEGER proc_code      ! value of processing code
      INTEGER mask_code      ! value of masking code
      REAL divisor           ! divisor for the time mean (1/period)
      INTEGER mod_period     ! timesteps since start modulo period.
      INTEGER start_time     ! value of start time
      INTEGER i              ! loop counter
      LOGICAL last_ts        ! true if end timestep
      INTEGER proc_size      ! size of data to be processed
!L---------------------------------------------------------------------
!L 1. Set processing option code and select appropriate service routine
!L
      proc_size=size-extra_size
      proc_code=control(st_proc_no_code)
!
!  Replace (null processing)
!
      IF (proc_code == st_replace_code) THEN
        DO i=1,size
          result(i)=variable(i)
        ENDDO
        start=(control(st_start_time_code) == timestep)
!
!  Mean/accumulation
!
      ELSEIF (proc_code == st_accum_code.or.                            &
     &        proc_code == st_time_mean_code) THEN
        start_time=control(st_start_time_code)
        IF (control(st_period_code) == st_infinite_time) THEN
          start=(timestep == start_time)
          last_ts=.FALSE.
        ELSE
          mod_period=mod(timestep-start_time,control(st_period_code))
          start=(mod_period == 0)
          last_ts=(mod_period == (control(st_period_code)-              &
     &                        control(st_freq_code)))
        ENDIF
        mask_code=control(st_gridpoint_code)
        mask_code=mod(mask_code,block_size)
        masking=(mask_code /= stash_null_mask_code).or.ocean
        IF (start) THEN      ! first timestep.
          DO i=1,size
            result(i)=variable(i)
          ENDDO
        ELSE
! DEPENDS ON: staccum
          CALL STACCUM(variable,result,proc_size,masking,amdi)
          DO i=proc_size+1,size
            result(i)=variable(i) ! copy over the extra data (if any)
          ENDDO
        ENDIF
!  Normalise at end of mean period
        IF (last_ts.and.proc_code == st_time_mean_code) THEN
          divisor=(float(control(st_freq_code))/                        &
     &             float(control(st_period_code)))
! If field is masked test for MDI, otherwise don't
          IF (masking) THEN
            DO i=1,proc_size
              IF (result(i) /= amdi) THEN
                result(i)=result(i)*divisor
              ENDIF
            ENDDO
          ELSE
            DO i=1,proc_size
              result(i)=result(i)*divisor
            ENDDO
          ENDIF
        ENDIF
!
!  Maximum
!
      ELSEIF (proc_code == st_max_code) THEN
        start_time=control(st_start_time_code)
        mod_period=mod(timestep-start_time,control(st_period_code))
        start=(mod_period == 0)
        IF (start) THEN
          DO i=1,size
            result(i)=variable(i)
          ENDDO
        ELSE
          mask_code=control(st_gridpoint_code)
          mask_code=mod(mask_code,block_size)
          masking=(mask_code /= stash_null_mask_code).or.ocean
! DEPENDS ON: stmax
          CALL STMAX(variable,result,proc_size,masking,amdi)
          DO i=proc_size+1,size
            result(i)=variable(i) ! copy over the extra data (if any)
          ENDDO
        ENDIF
!
!  Minimum
!
      ELSEIF (proc_code == st_min_code) THEN
        start_time=control(st_start_time_code)
        mod_period=mod(timestep-start_time,control(st_period_code))
        start=(mod_period == 0)
        IF (start) THEN
          DO i=1,size
            result(i)=variable(i)
          ENDDO
        ELSE
          mask_code=control(st_gridpoint_code)
          mask_code=mod(mask_code,block_size)
          masking=(mask_code /= stash_null_mask_code).or.ocean
! DEPENDS ON: stmin
          CALL STMIN(variable,result,proc_size,masking,amdi)
          DO i=proc_size+1,size
            result(i)=variable(i) ! copy over the extra data (if any)
          ENDDO
        ENDIF
!
!  Timeseries (append)
!
      ELSEIF (proc_code == st_time_series_code) THEN
        DO i=1,size
! Note that on start timestep this will include the extra data
          result(i)=variable(i)
        ENDDO
        start_time=control(st_start_time_code)
        mod_period=mod(timestep-start_time,control(st_period_code))
        start=(mod_period == 0)
        last_ts=(mod_period == (control(st_period_code)-                &
     &                      control(st_freq_code)))
!
!  Append trajectories
!
      ELSEIF (proc_code == st_append_traj_code) THEN
        start_time=control(st_start_time_code)
        mod_period=mod(timestep-start_time,control(st_period_code))
        start=(mod_period == 0)
        last_ts=(mod_period == (control(st_period_code)-                &
     &                      control(st_freq_code)))
        error=st_not_supported
        write(errmssg,100)' do not support append trajects'
        goto 999
!
!  Timeseries (append) - option 8 daily mean
!
      ELSEIF (proc_code == st_time_series_mean) THEN

        DO i=1,size
! Note that on start timestep this will include the extra data
          result(i)=variable(i)
        ENDDO
        start_time=control(st_start_time_code)
        mod_period=mod(timestep-start_time,control(st_period_code))
        start=(mod_period == 0)
        last_ts=(mod_period == (control(st_period_code)-                &
     &                      control(st_freq_code)))
!
!  Error condition
!
      ELSE
        error=unknown_processing
        write(errmssg,101)' unknown processing code',proc_code
        goto 999
      ENDIF
!
999   CONTINUE   ! jump for errors
!
100   FORMAT('TEMPORAL : >>> FATAL ERROR <<<',a30)
101   FORMAT('TEMPORAL : >>> FATAL ERROR <<<',a30,i5)
!
      RETURN
      END SUBROUTINE TEMPORAL
