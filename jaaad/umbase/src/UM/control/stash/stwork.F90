#if defined(C84_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine STWORK -------------------------------------------------
!LL
!LL  Purpose: Processes all the STASHlist entries for a particular
!LL item and section number after a timestep.  Raw input diagnostics to
!LL STWORK come either from D1, the main data array, or STASH_WORK, a
!LL work array dimensioned by the control level and passed in.  Each
!LL field is spatially processed, temporally processed and then output.
!LL The output destination can either be to an address in D1 or to a PP
!LL fieldsfile on a given unit.  In either case a PP-type LOOKUP header
!LL is generated to describe the contents of the output field.  Now
!LL handles atmosphere or ocean diagnostics according to arguments
!LL passed by calling routine.

! Method:
! Input information describing all STASH processing requests for a
! single model variable (STASH item,STASH section,internal model).
! Data input is either D1 array (primary variables) or STWORK array
! (diagnostic variables). Note that extraction of ocean/wave model
! fields requires 3-D de-compression using STOCGT/STWVGT routines.
!
! 0. Initialisation. Generate flags describing the grid location
!    for this STASH variable.
! Start loop over STASH requests for this STASH variable......
! 1. Determine input and output lengths and levels, and relative
!    positions, dependent on type of processing:
!    simple extraction, spatial, time-series, MOS.
!    Determine which type of processing required.
! 2. If spatial processing:
! 2.1  Timeseries. Extracts data and processes via MULTI_SPATIAL
!      routine.
! 2.2  Standard spatial processing: loop over output levels (unless
!      multi-level processing required) and process via SPATIAL.
! 3. If not spatial processing:
! 3.1  MOS output (output on pre-determined mask of locations).
! 3.2  Direct extraction/copying.
! 4. Output.
! 4.1  Determine output unit, PP lookup location and output lengths.
! 4.2  Pack output data if requested, via:
!      PP_FILE - COEX (WGDOS packing), or
!      GRIB_FILE - GRIB packing.
!      Set up PP lookup values in PP_HEAD
!      Write out data to PP file.
! 4.3  Write PP lookup field to PP file.
! 4.4  Output to D1 array (either to model dump or secondary space).
!      Time-processing (eg time means, accumulations) through TEMPORAL
!      routine if requested.
! End of loop over STASH requests..................
! 9. Error exits.
!
!LL
!LL  Author:   T.Johns/S.Tett
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered : C3, C4, C8
!LL
!LL  Project task: C4
!LL
!LL  External documentation : UMDP no C4
!LL
!*L  Interface and arguments: ------------------------------------------
      SUBROUTINE STWORK (                                               &
#include "argppx.h"
     &  D1,LEN_TOT,STASH_WORK,STASH_WORK_LEN,LENOUT,                    &
     &  global_LENOUT,                                                  &
     &  len_ocwork,                                                     &
     &  IS,IM,ILSTART,ILEND,STEP,steps_per_period,secs_per_period,      &
     &  previous_time,                                                  &
     &  STLIST,LEN_STLIST,TOTITEMS,SI,NSECTS,NITEMS,                    &
     &  STASH_LEVELS,NUM_STASH_LEVELS,NUM_LEVEL_LISTS,                  &
     &  STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,NUM_PSEUDO_LISTS,          &
     &  MAX_STASH_LEVS,STTABL,NSTTIMS,                                  &
     &  NSTTABL,STASH_SERIES,STASH_SERIES_LEN,                          &
     &  stash_series_rec_len,stash_series_index,stash_ser_index_size,   &
     &  MOS_MASK,MOS_MASK_LEN,MOS_OUTPUT_LENGTH,                        &
     &  PP_PACK_CODE,MODEL_FT_UNIT,FT_STEPS,FT_FIRSTSTEP,               &
     &  FIXHD,INTHD,REALHD,LEN_FIXHD,LEN_INTHD,LEN_REALHD,              &
     &  LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC,                              &
     &  LOOKUP,RLOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,                         &
     &  PP_LEN2_LOOKUP,NUNITS,PP_LEN2_LOOK,                             &
     &  LCYCLIC,LRLE,                                                   &
     &  T_rows,U_ROWS,ROW_LENGTH,T_field,U_FIELD,T_levels,              &
     &  RIVER_ROWS, RIVER_ROW_LENGTH,                                   &
     &  FCST_PRD,RUN_INDIC_OP,ELF,FT_LASTFIELD,                         &
     &  sm_ident,im_ident,dump_pack,                                    &
     &  stsuparrlen, stsuparr, istsuparr, sa_idx, sa_idxlen,            &
     &  ldump,                                                          &
     &  ICODE,CMESSAGE)

      USE FIELD_BUFF_MOD, ONLY :                                        &
     &    ATTACH_IPPLOOK,                                               &
     &    ATTACH_FXH

      IMPLICIT NONE
!
#include "csubmodl.h"
      INTEGER                                                           &
     &  TOTITEMS                                                        &
                           !IN    MAX NO OF ITEMS IN STASHLIST
     &, NSECTS                                                          &
                           !IN    MAX NO OF SECTIONS
     &, NITEMS                                                          &
                           !IN    MAX NO OF ITEMS IN A SECTION
     &, LEN_TOT                                                         &
                           !IN    LENGTH OF REAL DATA ARRAY D1
     &, len_ocwork                                                      &
                           !IN    Length of ocean work array OCWORK
     &, LEN_FIXHD                                                       &
                           !IN    LENGTH OF FIXED CONSTANTS
     &, LEN_INTHD                                                       &
                           !IN    LENGTH OF INTEGER CONSTANTS
     &, LEN_REALHD                                                      &
                           !IN   LENGTH OF REAL CONSTANTS
     &, LEN1_LEVDEPC                                                    &
                           !IN First dimension of LEVDEPC
     &, LEN2_LEVDEPC                                                    &
                           !IN Second dimension of LEVDEPC
     &, NUM_STASH_LEVELS                                                &
                           !IN   DIMENSION OF STASH_LEVELS
     &, NUM_LEVEL_LISTS                                                 &
                           !IN   DIMENSION OF STASH_LEVELS
     &, NUM_STASH_PSEUDO                                                &
                           !IN   Maximum no of pseudo-levels in a list
     &, NUM_PSEUDO_LISTS                                                &
                           !IN   Number of pseudo-level lists
     &, MAX_STASH_LEVS                                                  &
                           !IN   Max no of output levels for any diag
     &, LEN1_LOOKUP                                                     &
                           !IN   First dimension of LOOKUP/IPPLOOK
     &, LEN2_LOOKUP                                                     &
                           !IN   Second dimension of LOOKUP
     &, PP_LEN2_LOOKUP                                                  &
                           !IN   Largest poss. value in PP_LEN2_LOOK
     &, NUNITS                                                          &
                           !IN   Max i/o FT unit no
     &, PP_LEN2_LOOK(20:NUNITS)                                         &
                               !IN   Individual PP_LEN2_LOOKs per unit
     &, PP_PACK_CODE(20:NUNITS)                                         &
                               !IN   Packing code per unit
     &, FT_LASTFIELD(20:NUNITS)                                         &
                               !IN   Current write posn in each PP file
     &, FT_STEPS(20:NUNITS)                                             &
                               !IN   File reinitialisation freq per unit
     &, FT_FIRSTSTEP(20:NUNITS)                                         &
                               !IN   First step file initialised
     &, NSTTIMS                                                         &
                           !IN   Number of times against to test
     &, NSTTABL                                                         &
                           !IN   Number of STASH timetables
     &, NUM_WORDS                                                       &
                           !IN    Number of 64 Bit words to hold DATA
     &, sm_ident                                                        &
                           !IN    Submodel identifier
     &, im_ident                                                        &
                           !IN    Internal model identifier
     &, dump_pack                                                       &
                           !IN    Packing Indicator for Dump
     &, sa_idxlen                                                       &
                           !IN    Superarray index length
     &, sa_idx(sa_idxlen)                                               &
                           !IN    Superarray index
     &, stsuparrlen                                                     &
                           !IN    Superarray index length
     &, istsuparr(stsuparrlen)!IN Integer superarray
      CHARACTER*80                                                      &
     &  MODEL_FT_UNIT(NUNITS)  !IN   Current table of file associations

      INTEGER                                                           &
     &  FIXHD(LEN_FIXHD)                                                &
                           !IN    ARRAY OF FIXED CONSTANTS
     &, INTHD(LEN_INTHD)                                                &
                           !IN    ARRAY OF integer CONSTANTS
     &, ILSTART                                                         &
                           !IN    START OF LOOP OVER ENTRIES
     &, ILEND                                                           &
                           !IN    END OF LOOP OVER ENTRIES
     &, IS                                                              &
                           !IN    SECTION NUMBERS
     &, IM                                                              &
                           !IN    ITEM NUMBER
     &, STEP                                                            &
                           !IN    MODEL STEP NUMBER
     &, steps_per_period                                                &
                           !IN    No of steps in defining period
     &, secs_per_period                                                 &
                           !IN    No of secs in period (define timestep)
     &, previous_time(7)                                                &
                           !IN    Time at start of current step
     &, LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)                                 &
                                        ! Integer LOOKUP headers
     &, RLOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)                                &
                                         ! Real version of LOOKUP
     &, ICODE              !OUT   RETURN CODE FROM ROUTINE
      LOGICAL, intent(in) :: LDUMP
      INTEGER, POINTER :: PP_FIXHD(:)
!
      INTEGER                                                           &
     &  LENOUT                                                          &
                            !IN     Length of largest workfield needed
     &, global_LENOUT                                                   &
                            !IN     Output length of largest field
     &, T_field                                                         &
                            !IN     NO OF TEMP/PRESS POINTS
     &, U_FIELD                                                         &
                            !IN     NO OF U,V POINTS
     &, ROW_LENGTH                                                      &
                            !IN     No of points per row
     &, U_ROWS                                                          &
                            !IN     No of U,V, rows
     &, T_rows                                                          &
                            !IN     No of PRESS/TEMP rows
     &, T_levels                                                        &
                            !IN     No of model Press/Temp levels
     &, STASH_WORK_LEN                                                  &
                            !IN     LENGTH of STASH_WORK
     &, MOS_MASK_LEN                                                    &
                            !IN     Size of MOS_MASK array
     &, MOS_OUTPUT_LENGTH                                               &
                            !IN     No of MOS data pts extracted
     &, LEN_STLIST                                                      &
                            !IN     No of entries in STASHlist
     &, STLIST(LEN_STLIST,TOTITEMS)                                     &
                                    !IN STASHLIST
     &, SI(NITEMS,0:NSECTS,N_INTERNAL_MODEL)                            &
                                             !IN     STASH IN ADDRESS
     &, STTABL(NSTTIMS,NSTTABL)                                         &
                               !IN  STASH TIME TABLES
     &, STASH_LEVELS(NUM_STASH_LEVELS+1,NUM_LEVEL_LISTS)                &
     &, STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,NUM_PSEUDO_LISTS)        &
     &, MOS_MASK(MOS_MASK_LEN)                                          &
                               ! IN mask used to output data on MOS grid
     &, RUN_INDIC_OP                                                    &
                            !IN     Operational Run indicator (ITAB)
! STASH timeseries information
     &, STASH_SERIES_LEN                                                &
                             ! IN no of STASH SERIES records
     &, stash_series_rec_len                                            &
                             ! IN length of each record
     &, STASH_SERIES(stash_series_rec_len,STASH_SERIES_LEN)             &
!                            ! IN individual sample records
     &, stash_ser_index_size                                            &
                             ! IN no. of index blocks
     &, stash_series_index(2,stash_ser_index_size)                      &
!                            ! IN index block (1=start, 2=no of records)
     &, EXPPXI                                                          &
                             ! Function to extract ppxref info
     &, im_index                                                        &
                             ! Internal model index number
     &, N1                   ! Packing Indicator for Lookup(21)

! UM6.5 MODEL_ANALYSIS_HRS changed to REAL -
!             requires FCST_PRD changed to REAL also 
      REAL FCST_PRD  !IN     Forecast period

      INTEGER, INTENT(IN)  :: RIVER_ROWS         ! River routeing
      INTEGER, INTENT(IN)  :: RIVER_ROW_LENGTH   ! dimensions
!
      CHARACTER*(80) CMESSAGE ! OUT MESSAGE FROM ROUTINE
!
!
      LOGICAL                                                           &
     &  LCYCLIC                                                         &
                      !IN TRUE if cyclic EW BCs
     &, ELF                                                             &
                      !IN True if the input grid is rotated Equatorial
     &, LRLE                                                            &
                       !IN Run Length Encoding
     &, packing_hold
!

      REAL                                                              &
     &  D1(LEN_TOT)                                                     &
                                    !IN  REAL DATA ARRAY
     &, REALHD(LEN_REALHD)                                              &
                                    !IN  ARRAY OF REAL CONSTANTS
     &, LEVDEPC(LEN1_LEVDEPC*LEN2_LEVDEPC+1)                            &
                                                !IN level dep constants
     &, STASH_WORK(STASH_WORK_LEN)                                      &
                                    !IN    INPUT work array to STASH
     &, stsuparr(stsuparrlen)       !IN  Real superarray

!
!*----------------------------------------------------------------------
!
! Common blocks and PARAMETERs
!
#include "stparam.h"
#include "sterr.h"
#include "cppxref.h"
! Contains *CALL VERSION
#include "ppxlook.h"
#include "clookadd.h"
#include "c_mdi.h"
#include "parvars.h"
#include "decomptp.h"
#include "i_stgfld.h"
#include "comvgrid.h"

      INTEGER EXTRA_VAR_DATA ! Size of extra grid data

!
! External function:
      INTEGER GET_FLD_TYPE
!
! Local variables
!
      REAL                                                              &
     & PPFIELD(LENOUT)                                                  &
                           ! Main internal work array
     &, LEVEL(max_stash_levs)                                           &
                               ! The levels of the data as REAL nos
     &, sample_prd                                                      &
                           ! Sampling period in hours for means, etc
     &, A_IO               ! The output code from the unit command.
!

! I/O buffer workspace arrays - used for holding the data to be
! written to disk.

      REAL                                                              &
     &  buf(global_LENOUT)

!DIR$ CACHE_ALIGN buf

      INTEGER, POINTER :: IPPLOOK(:,:)
      INTEGER                                                           &
     & N_ROWS_OUT                                                       &
                          ! No of rows used for a diagnostic
     &, N_COLS_OUT                                                      &
                           ! No of cols used PPHORIZ=N_ROWS*N_COLS_
     &, SROW_IN,SROW_OUT                                                &
                           ! North, South, East and West
     &, NROW_IN,NROW_OUT                                                &
                           ! subdomain limits in the horizontal sense
     &, WCOL_IN,WCOL_OUT                                                &
                           ! corresponding to the subarea before being
     &, ECOL_IN,ECOL_OUT                                                &
                           ! processed (IN) and after processing (OUT)
     &, LEV_IN_ADDR                                                     &
                           ! The num of pts skipped in the Input
     &, GR                                                              &
                           ! Grid point code
     &, LBPROC_COMP(14)                                                 &
                           ! Array of 1/0 denoting LBPROC components
     &, UNITPP                                                          &
                           ! PPinit number (also used in PP_FILE)
     &, LENBUF                                                          &
                           ! PPHORIZ_OUT rnd to 512 words (used PP)
     &, COMP_ACCRCY                                                     &
                           ! PACKING ACCURACY IN POWER OF 2
     &, PPHORIZ_OUT                                                     &
                           ! No of points in the output field
     &, PPHORIZ_IN                                                      &
                           ! No of points in the input field
     &, IWA                                                             &
                           ! Record no used inSETPOS
     &, LEN_BUF_WORDS                                                   &
                           ! Number of 64 Bit words (rounded to 512)
     &, NUM_LEVS_OUT                                                    &
                           ! Number of output levels
     &, NUM_LEVS_IN                                                     &
                           ! Number of input levels
     &, this_index_lev                                                  &
                           ! index of level in output array for multi-
                           ! level processing in SPATIAL. Note:
                           ! loop over output levels:-
                           !  this_index_level=1,num_levs_out
                           !  level_list(this_ ) is model level or
                           !                        pressure level
                           !  index_lev(this_ ) is index to level in
                           !  input array
     &, INDEX_LEV(MAX_STASH_LEVS)                                       &
                                     ! Index used to relate input and
!                                    ! output levels
     &, level_list(MAX_STASH_LEVS)                                      &
                                     ! model level for each output level
     &, pseudo_level(MAX_STASH_LEVS)                                    &
                                     ! pseudo-level at each output level
     &, lv                                                              &
                                     ! LV code for field from STASHmster
     &, samples                                                         &
                                     ! no of samples (timeseries/trajec)
     &, icurrll_dump_ptr                                                &
                                     ! pointer to mother record LOOKUP
     &, start_time(7)                                                   &
                                     ! start time for PPheader in
!                                    ! timeseries/time mean/accumulation
     &, no_records                                                      &
                             ! no of records processed by multi_spatial
     &, record_start                                                    &
                             ! the start record for multi_spatial
     &, Jp_1                                                            &
                             ! Jpointer to first level of P in D1
     &, Jpstar                                                          &
                             ! Jpointer to PSTAR in D1
     &, len_field                                                       &
                          ! Holds the amount of data available
                          ! to PP_FILE
     &, packing_type_hold                                               &
                          ! Holds the packing type aftr packing by
                          ! (STASH_)GATHER_FIELD
     &, num_out           ! Number of 32-bit IBM words from COEX
      REAL                                                              &
     & Dummy3D(1,1,1), Dummy2D(1,1)
!
!dir$ cache_align ipplook
#include "cntl_io.h"

      INTEGER                                                           &
     &        JL,II,IL,JJ,IT,                                           &
     &        ntab,                                                     &
                           ! Number of the STASH TIMES table
     &        IKOUNT,                                                   &
                           ! Local  counter
     &        POINTS,                                                   &
                           ! No of points in a field
     &        KL,                                                       &
                           ! Local  counter
     &        ML,                                                       &
                           ! Local  counter
     &        LEN_IO,                                                   &
                           ! The length of data transferred.
     &        ILPREV,                                                   &
                           !The counter of the first of a pair of STLIST
     &        ILCURR,                                                   &
                           ! The current value of IL
     &        IWL,                                                      &
                           ! The word address of the LOOKUP table
     &        LBVCL,                                                    &
                           ! Vertical coordinate code in LOOKUP table
     &        ICURRLL,                                                  &
                           ! Current position in the PP Lookup table
     &        I                                                         &
                             ! loop variable
     &       ,J                                                         &
                           ! Level indicator used in call to GRIB_FILE
     &       ,PACKING_TYPE                                              &
                           ! 0 No packing, 1 for WGDOS, 3 for GRIB
     &       ,GRIB_PACKING ! ppxref profile number used to determine
                           ! grib packing method
      INTEGER VX,VY,VZ                                                  &
                           ! SIZES OF ARRAYS.
     &,       st_grid                                                   &
                           ! Horizontal grid type, (eg. T-p or u-v)
     &,LEN_PPNAME
      INTEGER INPUT_CODE   ! VALUE OF INPUT_CODE IN STASHLIST
      INTEGER ADDR         ! ADDRESS OF STASH VARIABLE IN EITHER DUMP OR
      INTEGER ADDR_OUT     ! ADDRESS OF SPATIALLY PROCESSED FIELD
      INTEGER ELAP_TIME    ! NO OF TIMESTEPS ELAPSED IN PERIOD.
      INTEGER SERIES_PTR   ! THE ADDRESS IN STASH_SERIES WHERE DOMIN INF
      INTEGER INDEX_SIZE   ! THE NUMBER OF LEVELS IN THE INDEX.
      INTEGER BASE_LEVEL   ! Base model level needed for mass weighting
      INTEGER top_level    ! Top model level for 3D ocean decompress
      INTEGER base_level0,top_level0 ! Ref base/top levels in levs loop
      INTEGER what_proc    ! What kind of processing will be done
      INTEGER what_mean    ! What kind of meaning will be done
      INTEGER output_code  ! Output destination code from STLIST
      INTEGER expected_len         ! expected length of output field
      INTEGER ocnlev_bottom  ! first ocean level diagnosed
      LOGICAL                                                           &
     &        S_F                                                       &
                             ! TRUE for items requiring processing
     &,       OCEAN                                                     &
                             ! TRUE if processing an ocean diagnostic
     &,       LLPROC                                                    &
                             ! TRUE if spatial processing required
     &,       lnullproc                                                 &
                             ! TRUE if null processing indicated
     &,       lfullfield                                                &
                             ! TRUE if output field on full horiz domain
     &,       lmasswt                                                   &
                             ! TRUE if level-by-level mass-weights exist
     &,       start_step                                                &
                             ! TRUE at start of time period (TEMPORAL)
     &,       MOS                                                       &
                             ! TRUE if MOS output is required
     &,       PACKING                                                   &
                             ! TRUE if packing required
     &,       GRIB_OUT                                                  &
                           ! TRUE if output to be in GRIB code.
     &,       ROTATE         ! TRUE if input data to be rotated
      CHARACTER*14                                                      &
     &        PPNAME           ! PPfile name decoded from MODEL_FT_UNIT
      CHARACTER*80                                                      &
     &        STRING         ! PPfile name decoded from MODEL_FT_UNIT

      integer expected_extra ! expected length of extra data
      INTEGER extraw        ! number of extra words this timestep
      INTEGER extraw_hdr    ! number of extra words for the header
      INTEGER data_type_code ! ppx_data_type code copied from PPX file
      INTEGER rotatecode   ! code for rotated grid
      INTEGER NT_DIM         ! Number of tracers
      INTEGER pointer_dummy  ! dummy pointer variable for ocean

      INTEGER                                                           &
! local versions of the subdomain limits
     &  local_NROW_OUT,local_SROW_OUT,local_ECOL_OUT,local_WCOL_OUT     &
     &, local_NROW_IN,local_SROW_IN,local_ECOL_IN,local_WCOL_IN         &
! global versions of the total size of output
     &, global_N_ROWS_OUT,global_N_COLS_OUT, global_PPHORIZ_OUT         &
! MOS variables
     &, global_NROWS                                                    &
     &, info ! return variable from GCOM
      INTEGER, PARAMETER  :: current_io_pe = 0

      INTEGER                                                           &
     & grid_type_code                                                   &
                           ! grid type of field being processed
     &,fld_type                                                         &
                           ! field type: u-,v- or p- location on C grid
     &,no_of_levels_masswt                                              &
                           ! no. levels for mass weights array
     &, halo_type  ! halo type of the field being processed

      INTEGER                                                           &
     &  open_flag                                                       &
                           ! is PP file open already?
     &, istat              ! GCOM return flag

!L----------------------------------------------------------------------
!L 0. Initialise: set constants relating to input grid type and size
!L
!L 0.1  Set up internal grid type st_grid and input field size
!L      according to the master GR code for the diagnostic
!L

! Get internal model index
      im_index = internal_model_index(im_ident)
      NT_DIM = (sa_idx(2) - sa_idx(1))/2

      OCEAN = .false.

! [Care needed! stsuparr passed as REAL but INTEGER value here.]
      Jp_1   = stsuparr(sa_idx(7))
      Jpstar = stsuparr(sa_idx(8))


! Get STASHmaster gridtype code
! DEPENDS ON: exppxi
      GR = EXPPXI( im_ident, IS, IM, ppx_grid_type,                     &
#include "argppx.h"
     &             icode, cmessage)
      grid_type_code=GR

! Get field type, ie u or v or p location in Arakawa C grid staggering
! DEPENDS ON: get_fld_type
      fld_type=GET_FLD_TYPE(grid_type_code) ! Function

! Find out what halo_type (and hence halo width) this field has

! DEPENDS ON: exppxi
      halo_type=EXPPXI( im_ident, IS, IM, ppx_halo_type,                &
#include "argppx.h"
     &                  icode, cmessage)

! Determine the input length pphoriz_in and stash grid staggering for
! this field.

      IF (GR == ppx_atm_tall.OR.GR == ppx_atm_tland.OR.                 &
     &    GR == ppx_atm_tsea) THEN
! Atmosphere data on T-grid
        st_grid=st_tp_grid
        PPHORIZ_IN=T_field
      ELSEIF (GR == ppx_atm_uall.OR.GR == ppx_atm_uland.OR.             &
     &        GR == ppx_atm_usea) THEN
! Atmosphere data on U-grid
        st_grid=st_uv_grid
        PPHORIZ_IN=U_FIELD
      ELSEIF(GR == ppx_atm_compressed) THEN
! Atmosphere data on T-grid (compressed)
        st_grid=st_tp_grid
        PPHORIZ_IN=T_field
      ELSEIF(GR == ppx_atm_cuall) THEN
! Atmosphere data on C-grid (u-points)
        st_grid=st_cu_grid
        PPHORIZ_IN=T_field
      ELSEIF(GR == ppx_atm_cvall) THEN
! Atmosphere data on C-grid (v-points)
        st_grid=st_cv_grid
        PPHORIZ_IN=U_field
      ELSEIF(GR == ppx_atm_tzonal) THEN
! Atmosphere zonal data on T-grid
        st_grid=st_zt_grid
        PPHORIZ_IN=T_rows
      ELSEIF(GR == ppx_atm_uzonal) THEN
! Atmosphere zonal data on u-grid
        st_grid=st_zu_grid
        PPHORIZ_IN=u_rows
      ELSEIF(GR == ppx_atm_tmerid) THEN
! Atmosphere meridional data on T-grid
        st_grid=st_mt_grid
        PPHORIZ_IN=row_length
      ELSEIF(GR == ppx_atm_umerid) THEN
! Atmosphere meridional data on u-grid
        st_grid=st_mu_grid
        PPHORIZ_IN=row_length
      ELSEIF(GR == ppx_atm_scalar) THEN
! Atmosphere scalar
        st_grid=st_scalar
        PPHORIZ_IN=1
      ELSEIF (GR == ppx_atm_river) THEN
! Atmosphere river routing
        st_grid=st_riv_grid
        PPHORIZ_IN=river_rows * river_row_length
      ELSE
! Unknown grid type
        ICODE=1
        CMESSAGE='STWORK   : Unknown grid type found in STASHmaster'
        GOTO 999
      ENDIF

! The input length pphoriz_in has been calculated explicitly depending
! on grid_type where there is an implicit assumption that diagnostic
! fields have no halos.

! For input fields in section 0, ie primary fields held in D1, this is
! not necessarily the case and input lengths are extracted using an
! addressing service routine dependent on grid and halo codes:
      IF(IS == 0 .OR. IS == 33 .OR. IS == 34) THEN ! section 0,33 or 34
! DEPENDS ON: addrln
         CALL ADDRLN(grid_type_code,halo_type,PPHORIZ_IN,local_data)
      ENDIF

!L
!L 0.2 Set up ROTATE to flag fields which are rotated (eg. ELF winds)
!L     (this is used to set alternative fieldcodes in PPHEAD)
!L

! DEPENDS ON: exppxi
      rotatecode = EXPPXI( im_ident, IS, IM, ppx_rotate_code,           &
#include "argppx.h"
     &             icode, cmessage)

      IF (rotatecode  ==  ppx_elf_rotated .AND. ELF)                    &
     &  THEN
        ROTATE=.TRUE.
      ELSE
        ROTATE=.FALSE.
      ENDIF


!L----------------------------------------------------------------------
!L 1. Loop over entries with this section/item number
!L
      DO 200 IL=ILSTART,ILEND  !loop over num entries for each item/sec
      extraw=0 ! no extra data by default
! Set MOS flag if output PP unit indicates output to MOS file
      MOS=(STLIST(st_output_code,IL) == -89)
!L
!L 1.1 Set up S_F which has to be set for each STASHLIST entry. The
!L     STASHFLAG is be set for a particular ITEM/SECTION the S_F for
!L     the STASHLIST entry.
!L
      S_F=.FALSE.
      IF (STLIST(st_freq_code,IL) == 1.AND.                             &
     &    STEP >= STLIST(st_start_time_code,IL).AND.                    &
     &   (STEP <= STLIST(st_end_time_code,IL).OR.                       &
     &    STLIST(st_end_time_code,IL) == st_infinite_time)) THEN
!       ... if required every step between start and end
        S_F=.TRUE.
      ELSEIF(STLIST(st_freq_code,IL) == -9999 ) THEN
         ! Note that this test overlaps with the next one so don't
         ! use STLIST(st_freq_code,IL) == -9999 .and. ldump
         S_F = ldump
      ELSEIF(STLIST(st_freq_code,IL) <  0) THEN
!       ... if required at specified times and this is one of them
        NTAB=-STLIST(st_freq_code,IL)
        DO 220 IT=1,NSTTIMS
          IF (STTABL(IT,NTAB) == st_end_of_list) GOTO 230
          IF (STEP == STTABL(IT,NTAB)) S_F=.TRUE.
  220   CONTINUE
      ELSEIF (STLIST(st_freq_code,IL) >  0) THEN
        IF   (MOD((STEP-STLIST(st_start_time_code,IL)),                 &
     &             STLIST(st_freq_code,IL)) == 0.AND.                   &
     &        STEP >= STLIST(st_start_time_code,IL).AND.                &
     &       (STEP <= STLIST(st_end_time_code,IL).OR.                   &
     &        STLIST(st_end_time_code,IL) == st_infinite_time))         &
!       ... if required every N timesteps and this is one of them
     &  S_F=.TRUE.
      ENDIF
  230 CONTINUE
!
!  S_F now set - Start of IF (S_F) block .......
!
      IF(S_F) THEN

!L
!L 1.2 Find number of input and output levels and relative positions
!L     and set up levels and pseudo-levels arrays for PPheaders.
!L     Set indicator lmasswt if level-by-level mass weighting possible
!L     - only currently available with atmosphere model full levels.
!L
! special case of mean timeseries leave ilcurr pointing to il

        ilcurr=il   ! The current STASHlist entry IL
        IF (STLIST(st_input_code,IL) <  0.and.                          &
     &      STLIST(st_proc_no_code,IL) /= st_time_series_mean) THEN
          ilcurr=-STLIST(st_input_code,il) ! points to prev entry
        ENDIF
!

! Get STASHmaster lbvc code
! DEPENDS ON: exppxi
      lbvcl = EXPPXI( im_ident, IS, IM, ppx_lbvc_code,                  &
#include "argppx.h"
     &             icode, cmessage)

! Get STASHmaster lv code
! DEPENDS ON: exppxi
      lv = EXPPXI( im_ident, IS, IM, ppx_lv_code,                       &
#include "argppx.h"
     &             icode, cmessage)

! DEPENDS ON: stlevels
        CALL STLEVELS(stlist(1,ilcurr),len_stlist,                      &
     &      stash_levels,num_stash_levels,num_level_lists,              &
     &      stash_pseudo_levels,num_stash_pseudo,num_pseudo_lists,      &
     &      max_stash_levs,num_levs_in,num_levs_out,index_size,         &
     &      index_lev,level_list,lbvcl,level,pseudo_level,              &
     &      icode,cmessage)
        IF (icode >  0) goto 999
        vz=num_levs_in

! Set switch to indicate whether level by level mass weights are
! available ( whether or not mass weighting requested)
        no_of_levels_masswt=1 ! default no. of levels for pre-
                              ! calculating mass weights array
        IF (lv == ppx_full_level) THEN
          lmasswt=.TRUE.
! Level by level mass weights available and mass weighting required:
! Extra space only needed if mass weighting requested.
          IF(stlist(st_weight_code,il) == stash_weight_mass_code) THEN
             no_of_levels_masswt = index_size  ! req. no. model levels
          ENDIF
        ELSE
          lmasswt=.FALSE.
        ENDIF

!L
!L 1.3 Find the horizontal dimensions for the output grid
!L     and the input field subdomain limits in which processing happens.
!L
        WCOL_IN= STLIST(st_west_code,ILCURR)    ! Input subdomain limits
        ECOL_IN= STLIST(st_east_code,ILCURR)
        NROW_IN= STLIST(st_north_code,ILCURR)
        SROW_IN= STLIST(st_south_code,ILCURR)
!
!L
!L 1.3.2 Other output types need to calculate lengths in detail
!L       (ie. number of rows, columns, and horizontal field size)
!L       according to processing options
!L

! Calculate local versions of the subdomain limits and area

! DEPENDS ON: global_to_local_subdomain
          CALL GLOBAL_TO_LOCAL_SUBDOMAIN(                               &
     &      .TRUE.,.TRUE.,                                              &
     &      grid_type_code,halo_type,mype,                              &
     &      SROW_IN,ECOL_IN,NROW_IN,WCOL_IN,                            &
     &      local_SROW_IN,local_ECOL_IN,                                &
     &      local_NROW_IN,local_WCOL_IN)

          what_proc=stlist(st_proc_no_code,ilcurr)
          what_mean=(stlist(st_gridpoint_code,ilcurr)/block_size)*      &
     &       block_size
          samples=0         ! Initialise value for non-timeseries
!L
!L 1.3.2.1 Time series or trajectory processing
!L
          IF (what_proc == st_time_series_code.or.                      &
     &        what_proc == st_append_traj_code.or.                      &
     &        what_proc == st_time_series_mean) THEN
!L
!L 1.3.2.2 Compute number of samples in period for timeseries for
!L         input to PP_HEAD.
!L         No of output rows and cols are set to the no of points
!L         in each time sample and number of time samples in the
!L         period spanned by the output field, respectively.
!L
            samples=stlist(st_period_code,ilcurr)/                      &
     &             stlist(st_freq_code,ilcurr)
            wcol_out=1
            ecol_out=samples
            nrow_out=stlist(st_output_length,ilcurr)/samples
            srow_out=1

            local_WCOL_OUT=1
            local_ECOL_OUT=samples
            local_NROW_OUT=stlist(st_output_length,ilcurr)/samples
            local_SROW_OUT=1

!L 1.3.2.3 Multi spatial processing of some other type (not supported)
          ELSEIF (stlist(st_series_ptr,ilcurr) /= 0) THEN
            ICODE=1323
            CMESSAGE='STWORK   : Illegal timeseries processing selected'
            GOTO 999
!L 1.3.2.4 Primary record requesting an extract
          ELSEIF (what_mean == extract_base) THEN
            WCOL_OUT= WCOL_IN
            ECOL_OUT= ECOL_IN
            NROW_OUT= NROW_IN
            SROW_OUT= SROW_IN

            local_WCOL_OUT= local_WCOL_IN
            local_ECOL_OUT= local_ECOL_IN
            local_NROW_OUT= local_NROW_IN
            local_SROW_OUT= local_SROW_IN

!L 1.3.2.5 Primary record requesting a vertical mean
          ELSEIF (what_mean == vert_mean_base) THEN
            WCOL_OUT= WCOL_IN
            ECOL_OUT= ECOL_IN
            NROW_OUT= NROW_IN
            SROW_OUT= SROW_IN

            local_WCOL_OUT= local_WCOL_IN
            local_ECOL_OUT= local_ECOL_IN
            local_NROW_OUT= local_NROW_IN
            local_SROW_OUT= local_SROW_IN

!L 1.3.2.6 Primary record requesting a zonal mean
          ELSEIF (what_mean == zonal_mean_base) THEN
            WCOL_OUT= 1
            ECOL_OUT= 1
            NROW_OUT= NROW_IN
            SROW_OUT= SROW_IN

            local_WCOL_OUT= 1
            local_ECOL_OUT= 1
            local_NROW_OUT= local_NROW_IN
            local_SROW_OUT= local_SROW_IN

!L 1.3.2.7 Primary record requesting a meridional mean
          ELSEIF (what_mean == merid_mean_base) THEN
            WCOL_OUT= WCOL_IN
            ECOL_OUT= ECOL_IN
            NROW_OUT= 1
            SROW_OUT= 1

            local_WCOL_OUT= local_WCOL_IN
            local_ECOL_OUT= local_ECOL_IN
            local_NROW_OUT= 1
            local_SROW_OUT= 1

!L 1.3.2.8 Primary record requesting a global mean
          ELSEIF (what_mean == global_mean_base) THEN
            WCOL_OUT= 1
            ECOL_OUT= 1
            NROW_OUT= 1
            SROW_OUT= 1

            local_WCOL_OUT= 1
            local_ECOL_OUT= 1
            local_NROW_OUT= 1
            local_SROW_OUT= 1

!L 1.3.2.9 Primary record requesting a field mean
          ELSEIF (what_mean == field_mean_base) THEN
            WCOL_OUT= 1
            ECOL_OUT= 1
            NROW_OUT= 1
            SROW_OUT= 1

            local_WCOL_OUT= 1
            local_ECOL_OUT= 1
            local_NROW_OUT= 1
            local_SROW_OUT= 1

!L 1.3.2.10 Error trap for unknown request
          ELSE         ! Invalid option
            icode=st_unknown
            write(cmessage,87)'unknown option in setup',what_mean
            goto 999 ! jump to error return
          ENDIF
!L
!L 1.3.3 Compute expected length. This differs from total output length
!L       when data is appended from multiple timesteps into the same
!L       field, being output_length/number_of_appends in this case.
!L
          IF (stlist(st_output_code,il) >= 0.and.                       &
     &        (what_proc == st_time_series_code.or.                     &
     &        what_proc == st_append_traj_code.or.                      &
     &        what_proc == st_time_series_mean)) THEN
          series_ptr=stlist(st_series_ptr,il) !set up ptr to stashseries
            expected_extra=(stash_series_index(2,series_ptr)+1)*6
            extraw_hdr=expected_extra
            expected_len=((stlist(st_output_length,ilcurr)              &
     &        -expected_extra)*                                         &
     &      stlist(st_freq_code,ilcurr))/stlist(st_period_code,ilcurr)
          ELSE
            expected_len=stlist(st_output_length,ilcurr)
            expected_extra=0 ! no extra data for non timeseries stuff
            extraw_hdr=0
          ENDIF
!L
!L 1.3.6 Compute number of rows and columns and field size for output
!L       - first adjust easternmost column if field wraps EW
!L
          IF (WCOL_IN  >  ECOL_IN .AND.LCYCLIC)                         &
     &      ECOL_IN =ECOL_IN + glsize(1,fld_type)

          IF (WCOL_OUT >  ECOL_OUT.AND.LCYCLIC)                         &
     &      ECOL_OUT=ECOL_OUT + glsize(1,fld_type)

!

          IF (local_WCOL_OUT  >   local_ECOL_OUT)                       &
     &      local_ECOL_OUT=local_ECOL_OUT+blsize(1,fld_type)

          N_ROWS_OUT = local_NROW_OUT - local_SROW_OUT + 1
          N_COLS_OUT = local_ECOL_OUT - local_WCOL_OUT + 1
          global_N_ROWS_OUT =  NROW_OUT - SROW_OUT + 1
          global_N_COLS_OUT = ECOL_OUT - WCOL_OUT + 1

          PPHORIZ_OUT= N_ROWS_OUT*N_COLS_OUT
          global_PPHORIZ_OUT=global_N_ROWS_OUT*global_N_COLS_OUT

!L
!L 1.4 Check to see if any processing is required.
!L     Set flag LLPROC if some SPATIAL processing indicated.
!L     NB: If input and output bottom levels differ (or the input and
!L         output pseudo-levels lists differ), level-by-level
!L         processing in the SPATIAL loop IS required.
!L         MULTI-SPATIAL processing is always required for timeseries.
!L
        lfullfield=((st_grid == st_tp_grid .OR. st_grid == st_cu_grid   &
     &              .OR. st_grid == st_riv_grid)                        &
     &       .AND.  stlist(st_west_code,il) == 1.and.                   &
     &              stlist(st_east_code,il) == glsize(1,fld_type).and.  &
     &              stlist(st_south_code,il) == 1.and.                  &
     &              stlist(st_north_code,il) == glsize(2,fld_type)).or. &
     &             ((st_grid == st_uv_grid .OR. st_grid == st_cv_grid)  &
     &       .AND.  stlist(st_west_code,il) == 1.and.                   &
     &              stlist(st_east_code,il) == glsize(1,fld_type).and.  &
     &              stlist(st_north_code,il) == glsize(2,fld_type).and. &
     &              stlist(st_south_code,il) == 1) .OR.                 &
     &             ((st_grid == st_zt_grid)                             &
     &       .AND.  stlist(st_north_code,il) == glsize(2,fld_type).and. &
     &              stlist(st_south_code,il) == 1) .OR.                 &
     &             ((st_grid == st_zu_grid)                             &
     &       .AND.  stlist(st_north_code,il) == glsize(2,fld_type).and. &
     &              stlist(st_south_code,il) == 1) .OR.                 &
     &             ((st_grid == st_mt_grid .OR. st_grid == st_mu_grid)  &
     &       .AND.  stlist(st_west_code,il) == 1.and.                   &
     &              stlist(st_east_code,il) == glsize(1,fld_type)) .OR. &
     &             (st_grid == st_scalar)

        lnullproc= lfullfield .AND.                                     &
     &             (stlist(st_input_bottom,il) ==                       &
     &                stlist(st_output_bottom,il)) .and.                &
     &             (stlist(st_pseudo_in,il) ==                          &
     &                stlist(st_pseudo_out,il)) .and.                   &
     &             (stlist(st_gridpoint_code,il) ==                     &
     &               (extract_base+stash_null_mask_code) .and.          &
     &              stlist(st_weight_code,IL) ==                        &
     &                stash_weight_null_code )
        IF (STLIST(st_series_ptr,IL) >  0) THEN
          lnullproc=.FALSE.     ! Timeseries always requires processing
        ENDIF
!  LLPROC must be false for MOS output, output from a prev STLIST
!  or simple extraction of full field with no weighting
        IF ((STLIST(st_input_code,IL) < 0 .and.                         &
     &      STLIST(st_proc_no_code,IL) /= st_time_series_mean) .or.     &
     &      lnullproc) THEN
          LLPROC=.FALSE.
        ELSE
          LLPROC=.TRUE.
        ENDIF
!L
!L 1.5 Check that no spatial processing is requested if the input field
!L     is of integer or logical type -- these types of fields can be
!L     passed directly through STASH, for example for coupling purposes,
!L     but no arithmetic is allowed at present.
!L
! DEPENDS ON: exppxi
        data_type_code=EXPPXI(im_ident,IS,IM,ppx_data_type,             &
#include "argppx.h"
     &                        icode,cmessage)
!
        IF ((       data_type_code  == ppx_type_int    .OR.             &
     &              data_type_code  == ppx_type_log)   .AND.            &
     &       .NOT. lnullproc)  THEN
          ICODE=st_not_supported
          CMESSAGE='STWORK  : Spatial processing of INT/LOGICAL illegal'
          GOTO 999
        ENDIF
!L----------------------------------------------------------------------
!L 2. Perform spatial processing (loop over output levels)
      ! If this is a variable horizontal grid set flags
      ! to ensure that BDX and BDY are set properly
      ! in the pp header. Also work out the size of
      ! extra data which will be generated. This must be done
      ! here because the size is needed for the LOOKUP headers
      ! which are created long before the actual ouput data
      ! is set up.
      VAR_GRID_TYPE = 0
      EXTRA_VAR_DATA = 0

      IF (X_VAR_GRID .OR. Y_VAR_GRID) THEN

         ! Each dimesion has 3 lots of extra grid data -
         ! grid coordinates, lower and upper box boundaries.
         ! May need to incorporate test to exclude certain
         ! grid types (eg LBCs) if they're not appropriate.
         IF (fld_type == fld_type_p) THEN
            VAR_GRID_TYPE = 1 ! It's a T field
            IF (X_VAR_GRID) EXTRA_VAR_DATA =                            &
     &                      EXTRA_VAR_DATA + (3*(global_N_COLS_OUT+1))
            IF (Y_VAR_GRID) EXTRA_VAR_DATA =                            &
     &                      EXTRA_VAR_DATA + (3*(global_N_ROWS_OUT+1))
         ELSEIF (fld_type == fld_type_u.OR.                             &
     &           fld_type == fld_type_v) THEN
            VAR_GRID_TYPE = 2 ! It's a U,V field
            IF (X_VAR_GRID) EXTRA_VAR_DATA =                            &
     &                      EXTRA_VAR_DATA + (3*(global_N_COLS_OUT+1))
            IF (Y_VAR_GRID) EXTRA_VAR_DATA =                            &
     &                      EXTRA_VAR_DATA + (3*(global_N_ROWS_OUT+1))
         ENDIF

         extraw_hdr = extraw_hdr + EXTRA_VAR_DATA

      ENDIF
!L
        IF (LLPROC) THEN   ! Processing is required
          input_code=stlist(st_input_code,il)
!L make sure no volume processing asked for as not supported
!L this will need adding at some point
          IF (stlist(st_weight_code,il) == stash_weight_volume_code)    &
     &      THEN
            icode=st_not_supported
            cmessage='STWORK  : volume processing not supported'
            goto 999
          ENDIF
! Work out vx,vy (depends on kind of grid)
          vx=lasize(1,fld_type,halo_type)
          vy=lasize(2,fld_type,halo_type)

          IF ((st_grid == st_zt_grid) .OR.                              &
     &        (st_grid == st_zu_grid)) THEN
            vx=1
          ELSEIF ((st_grid == st_mt_grid) .OR.                          &
     &            (st_grid == st_mu_grid)) THEN
            vy=1
          ELSEIF (st_grid == st_scalar) THEN
            vx=1
            vy=1
          ENDIF

!L Work out if this is the first timestep in a timeseries.
!L This is required so that the extra data can be generated
!
          series_ptr=stlist(st_series_ptr,il)
          IF (series_ptr >  0) THEN ! multi spatial processing reqd.
!L recompute expected sizes
            elap_time=step-stlist(st_start_time_code,il)
            elap_time=mod(elap_time,stlist(st_period_code,il))
            start_step=(elap_time == 0)
            IF (start_step) THEN
              expected_len=stlist(st_output_length,ilcurr)
            ELSE
              expected_extra=0  ! reset to zero as no extra data
              expected_len=((stlist(st_output_length,ilcurr)            &
     &          -((stash_series_index(2,series_ptr)+1)*6))*             &
     &        stlist(st_freq_code,ilcurr))/stlist(st_period_code,ilcurr)
            ENDIF
!L
!L 2.1 Timeseries extraction section follows
!L
            no_records=stash_series_index(2,series_ptr)
            record_start=stash_series_index(1,series_ptr)
!L
!L 2.1.1 Process a primary field from D1 (timeseries)
!L
            IF (input_code == 0) THEN
                 addr=si(im,is,im_index)
! DEPENDS ON: multi_spatial
                CALL MULTI_SPATIAL(d1(addr),                            &
     &            vx,vy,vz,grid_type_code,st_grid,fld_type,halo_type,   &
     &            halosize(1,halo_type),halosize(2,halo_type),          &
     &            lcyclic,lmasswt,                                      &
     &            pphoriz_out,num_levs_out,                             &
     &                   no_of_levels_masswt,                           &
!     Extra arguments for atmos sub-model
     &            D1(Jp_1),D1(Jpstar),                                  &
                                              ! pressure,pstar
     &            stsuparr(sa_idx(9)),                                  &
                                              ! cos_v_latitude
     &            stsuparr(sa_idx(10)),                                 &
                                              ! cos_theta_latitude
     &            stsuparr(sa_idx(11)),                                 &
                                              ! land mask
     &            stsuparr(sa_idx(12)),                                 &
                                              ! sea mask
!
     &            row_length,T_rows,u_rows,T_levels,                    &
     &            ppfield,lenout,                                       &
     &            rmdi,stlist(1,il),len_stlist,                         &
     &            stash_series(1,record_start),                         &
     &            stash_series_rec_len,no_records,                      &
     &            index_size,index_lev,level_list,                      &
     &            start_step,extraw,n_rows_out,n_cols_out,              &
     &            realhd,len_realhd,inthd,len_inthd,ocean,              &
     &            icode,cmessage)
!L
!L 2.1.2 Process a field from STASHWORK (timeseries)
!L
            ELSEIF (input_code == 1) THEN
                 addr=si(im,is,im_index)
! DEPENDS ON: multi_spatial
                CALL MULTI_SPATIAL(stash_work(addr),                    &
     &            vx,vy,vz,grid_type_code,st_grid,fld_type,halo_type,   &
     &            halosize(1,halo_type),halosize(2,halo_type),          &
     &            lcyclic,lmasswt,                                      &
     &            pphoriz_out,num_levs_out,                             &
     &                   no_of_levels_masswt,                           &
!     Extra arguments for atmos sub-model
     &            D1(Jp_1),D1(Jpstar),                                  &
                                              ! pressure,pstar
     &            stsuparr(sa_idx(9)),                                  &
                                              ! cos_v_latitude
     &            stsuparr(sa_idx(10)),                                 &
                                              ! cos_theta_latitude
     &            stsuparr(sa_idx(11)),                                 &
                                              ! land mask
     &            stsuparr(sa_idx(12)),                                 &
                                              ! sea mask
!
     &            row_length,T_rows,u_rows,T_levels,                    &
     &            ppfield,lenout,                                       &
     &            rmdi,stlist(1,il),len_stlist,                         &
     &            stash_series(1,record_start),                         &
     &            stash_series_rec_len,no_records,                      &
     &            index_size,index_lev,level_list,                      &
     &            start_step,extraw,n_rows_out,n_cols_out,              &
     &            realhd,len_realhd,inthd,len_inthd,ocean,              &
     &            icode,cmessage)
            ELSEIF (input_code <  0) then
!L
!L 2.1.3 Process a field from previously STASHed position in D1
!L     (currently unsupported since diagnostic-of-diagnostic)
!L
             IF (what_proc == st_time_series_mean) THEN
! special case of mean timeseries
!   Mother record
                 ILPREV=-stlist(st_input_code,IL)
! address of mother record in D1
                 addr=stlist(20,ILPREV)
! DEPENDS ON: multi_spatial
                CALL MULTI_SPATIAL(D1(addr),                            &
     &            vx,vy,vz,grid_type_code,st_grid,fld_type,halo_type,   &
     &            halosize(1,halo_type),halosize(2,halo_type),          &
     &            lcyclic,lmasswt,                                      &
     &            pphoriz_out,num_levs_out,                             &
     &                   no_of_levels_masswt,                           &
!     Extra arguments for atmos sub-model
     &            D1(Jp_1),D1(Jpstar),                                  &
                                              ! pressure,pstar
     &            stsuparr(sa_idx(9)),                                  &
                                              ! cos_v_latitude
     &            stsuparr(sa_idx(10)),                                 &
                                              ! cos_theta_latitude
     &            stsuparr(sa_idx(11)),                                 &
                                              ! land mask
     &            stsuparr(sa_idx(12)),                                 &
                                              ! sea mask
!
     &            row_length,T_rows,u_rows,T_levels,                    &
     &            ppfield,lenout,                                       &
     &            rmdi,stlist(1,il),len_stlist,                         &
     &            stash_series(1,record_start),                         &
     &            stash_series_rec_len,no_records,                      &
     &            index_size,index_lev,level_list,                      &
     &            start_step,extraw,n_rows_out,n_cols_out,              &
     &            realhd,len_realhd,inthd,len_inthd,ocean,              &
     &            icode,cmessage)
             ELSE
              icode=st_not_supported
              cmessage='STWORK1  : diag-of-diagnostic unsupported'
              goto 999 ! jump to error return
             ENDIF
            ELSE
              icode=st_unknown
              write(cmessage,87)'unknown input option',input_code
              goto 999
            ENDIF
            if (icode /= 0) goto 999 ! error exit
            IF (start_step.AND.(extraw /= expected_extra)) THEN
              icode=st_bad_array_param
              write(cmessage,89) extraw,expected_extra
 89           FORMAT('STWORK : Inconsistent length for extra data ',    &
     &          i8,1x,i8)
              goto 999
            endif
!
!L         pphoriz_out has been computed by multi_spatial
!L         it is the size of the output field
          ELSE  ! do "normal" spatial processing
!L
!L 2.2 Standard spatial processing section follows
!L
!L In ocean case, a possibly 3D decompression is required (to top_level)
!L
!L If multi-level processing (ie. vertical, global mean) is performed by
!L SPATIAL, the input field is passed in with the original start address
!L but if single-level processing is done by SPATIAL, the field is
!L passed in with an address pointing to the single level required.
!L
            base_level0=stlist(st_input_bottom,il)
            what_proc=stlist(st_gridpoint_code,il)
!
            addr_out=1                   ! Initialise output address
!
            DO kl=1,num_levs_out         ! --- Start of levels loop ---
! Work out model level if model level range, otherwise set to 1
              IF (base_level0 <  0.OR.base_level0 == st_special_code)   &
     &        THEN
                base_level=1
              ELSE
                base_level=base_level0+index_lev(kl)-1
              ENDIF
! Pass level index into SPATIAL (instead of base_level)
            this_index_lev = kl
!L
!L 2.2.1 Process a primary field from D1 
!L
              IF (input_code == 0) THEN
                  IF ((what_proc <  vert_mean_top .AND.                 &
     &                 what_proc >  vert_mean_base) .OR.                &
     &                (what_proc <  global_mean_top .AND.               &
     &                 what_proc >  global_mean_base)) THEN
                    addr=si(im,is,im_index)
                  ELSE
                    addr=si(im,is,im_index)+                            &
     &                   (index_lev(kl)-1)*pphoriz_in
                  ENDIF
                  IF (addr <  1.or.addr >  len_tot) THEN
                    icode=st_bad_address
                    cmessage='STWORK  : D1 address out of bounds'
                    goto 999
                  ENDIF
! DEPENDS ON: spatial
                  CALL SPATIAL(d1(addr),vx,vy,vz,                       &
     &                grid_type_code,st_grid,                           &
     &                fld_type,halo_type,                               &
     &                halosize(1,halo_type),halosize(2,halo_type),      &
     &                lcyclic,lmasswt,                                  &
     &                n_cols_out,n_rows_out,this_index_lev,             &
     &                level_list,index_lev,index_size,                  &
     &                   no_of_levels_masswt,                           &
!     Extra arguments for atmos sub-model
     &            D1(Jp_1),D1(Jpstar),                                  &
                                              ! pressure,pstar
     &            stsuparr(sa_idx(9)),                                  &
                                              ! cos_v_latitude
     &            stsuparr(sa_idx(10)),                                 &
                                              ! cos_theta_latitude
     &            stsuparr(sa_idx(11)),                                 &
                                              ! land mask
     &            stsuparr(sa_idx(12)),                                 &
                                              ! sea mask
!
     &                row_length,T_rows,u_rows,                         &
     &                blsize(2,fld_type),T_levels,                      &
     &                ppfield(addr_out),pphoriz_out,                    &
     &                stlist(1,il),len_stlist,rmdi,                     &
     &                icode,cmessage)
!L
!L 2.2.2 Process a field from STASHWORK (or ocwork if OCEAN)
!L
              ELSEIF (input_code == 1) THEN
                  IF ((what_proc <  vert_mean_top .AND.                 &
     &                 what_proc >  vert_mean_base) .OR.                &
     &                (what_proc <  global_mean_top .AND.               &
     &                 what_proc >  global_mean_base)) THEN
                    addr=si(im,is,im_index)
                  ELSE
                    addr=si(im,is,im_index)+                            &
     &                   (index_lev(kl)-1)*pphoriz_in
                  ENDIF
                  IF (addr <  1.or.addr >  stash_work_len) THEN
                    icode=st_bad_address
                    cmessage='STWORK  : STASHWORK addr out of bounds'
                    goto 999
                  ENDIF
! DEPENDS ON: spatial
                  CALL SPATIAL(stash_work(addr),vx,vy,vz,               &
     &                grid_type_code,st_grid,                           &
     &                fld_type,halo_type,                               &
     &                halosize(1,halo_type),halosize(2,halo_type),      &
     &                lcyclic,lmasswt,                                  &
     &                n_cols_out,n_rows_out,this_index_lev,             &
     &              level_list,index_lev,index_size,                    &
     &                   no_of_levels_masswt,                           &
!     Extra arguments for atmos sub-model
     &            D1(Jp_1),D1(Jpstar),                                  &
                                              ! pressure,pstar
     &            stsuparr(sa_idx(9)),                                  &
                                              ! cos_v_latitude
     &            stsuparr(sa_idx(10)),                                 &
                                              ! cos_theta_latitude
     &            stsuparr(sa_idx(11)),                                 &
                                              ! land mask
     &            stsuparr(sa_idx(12)),                                 &
                                              ! sea mask
     &              row_length,T_rows,u_rows,                           &
     &              blsize(2,fld_type),T_Levels,                        &
     &              ppfield(addr_out),pphoriz_out,                      &
     &              stlist(1,il),len_stlist,rmdi,                       &
     &              icode,cmessage)
              ELSEIF (input_code <  0) THEN
!L
!L 2.2.3 Process a field from previously STASHed position in D1
!L     (currently unsupported since diagnostic-of-diagnostic)
!L
                icode=st_not_supported
                cmessage='STWORK1  : diag-of-diagnostic unsupported'
                goto 999
              ELSE
                icode=st_unknown
                write(cmessage,87)'unknown input option',input_code
87              format('STWORK1 : >>FATAL ERROR <<',a,1x,i5)
                goto 999
              ENDIF
!
              IF (icode >  0) goto 999 ! Trap error
!
!L compute pphoriz_out
!L pphoriz_out is the size of the output vector
!L we should not be doing timeseries processing here.
!L
!L NOTE: n_cols_out and n_rows_out should agree with values calculated
!L       before, but are not checked for consistency.
!
              pphoriz_out=n_cols_out*n_rows_out
              addr_out=addr_out+pphoriz_out ! increment output address
            ENDDO                      ! --- End of levels loop ---
!
          ENDIF         ! End of multi-spatial/spatial IF block
!
          IF (icode >  0) goto 999     ! Trap processing error
!L
!L 2.3 Set length of output field and check against expected length
!L
!L check that extrawords are the same as the expected number of extrawor

! Calculate size of global pphoriz_out - the size on disk

          IF (what_proc  ==  st_time_series_code .OR.                   &
     &        what_proc  ==  st_time_series_mean) THEN
            global_pphoriz_out=pphoriz_out
            global_n_rows_out=n_rows_out
            global_n_cols_out=n_cols_out
          ELSE
! DEPENDS ON: stash_get_global_size
            CALL STASH_GET_GLOBAL_SIZE(                                 &
     &       stlist(st_north_code,il) , stlist(st_east_code,il),        &
     &       stlist(st_south_code,il) , stlist(st_west_code,il),        &
     &       1,                                                         &
     &       STLIST(st_gridpoint_code,il) , STLIST(st_proc_no_code,il), &
     &       global_pphoriz_out,                                        &
     &       ICODE, CMESSAGE)

            IF (icode  /=  0) goto 999

          ENDIF


          IF (pphoriz_out*num_levs_out /= expected_len) THEN
            icode=st_bad_array_param
            write(cmessage,88) pphoriz_out*num_levs_out,expected_len
 88         FORMAT('STWORK   : Inconsistent length for output field ',  &
     &             i8,1x,i8)
            goto 999
          ENDIF
        ELSE
!L----------------------------------------------------------------------
!L 3. No SPATIAL processing - extract output field by direct copy
!L
!L 3.1 MOS output. (No longer supported)
!L
!L 3.2 Other output - determine input source
!L
            input_code=STLIST(st_input_code,IL)
!L
!L 3.2.1 Other fields are simply copied
!L
            IF (input_code == 0) THEN
! Simple extraction with no weighting from primary field in D1
! except for those needing special extraction on funny grids
! (Ocean, Wave)
              addr=si(im,is,im_index)
              DO JL=1,STLIST(st_output_length,IL)
                PPFIELD(JL)=D1(addr+JL-1)
              ENDDO
            ELSEIF (input_code == 1) THEN
! Simple extraction with no weighting from STASH_WORK
! except for those needing special extraction on funny grids
! (Ocean, Wave)
              addr=si(im,is,im_index)
              DO JL=1,STLIST(st_output_length,IL)
                PPFIELD(JL)=STASH_WORK(addr+JL-1)
              ENDDO
            ELSEIF (input_code <  0) THEN
! Previously STASHed entry in D1
! for all sub-models as diagnostic D1 is always on a proper grid.

              addr=STLIST(st_output_addr,-input_code)
              DO JL=1,STLIST(st_output_length,IL)
                PPFIELD(JL)=D1(addr+JL-1)
              ENDDO
            ELSE
! Illegal input code
              ICODE=st_unknown
              CMESSAGE='STWORK   : Unknown input code encountered'
              GOTO 999
            ENDIF
        ENDIF      ! End of LLPROC IF BLOCK    ************************

!L-----------------------------------------------------------------
!L 4. OUTPUT section.
!L
!L    The data is in PPFIELD with a length LENOUT.
!L    The horizontal field size PPHORIZ_OUT and number of output levels
!L    NUM_LEVS_OUT were calculated in section 1.
!L    Output option depends on the STLIST code.
!L
!L 4.0 Find mother STASHlist record if necessary.
!L
!   Packing_type not set from PP_FILE in the ELSE IF part of block.
            PACKING_TYPE = 0  ! Default is unpacked.
        IF(STLIST(st_input_code,IL) <  0) THEN ! Second of two STLIST
          ILPREV=-STLIST(st_input_code,IL)
        ELSE
          ilprev=IL ! no daughter record
        ENDIF
!L
!L 4.0.1 Set up LBPROC sub-components based on STASH processing info.
!L
        DO JJ=1,14
          LBPROC_COMP(JJ)=0
        ENDDO
!
        IF(STLIST(st_gridpoint_code,ilprev) >= zonal_mean_base)         &
     &    LBPROC_COMP(7)=1
!
        IF((STLIST(st_gridpoint_code,ilprev) >= vert_mean_base .AND.    &
     &      STLIST(st_gridpoint_code,ilprev) <  vert_mean_top) .OR.     &
     &    (STLIST(st_gridpoint_code,ilprev) >= global_mean_base .AND.   &
     &     STLIST(st_gridpoint_code,ilprev) <  global_mean_top))        &
     &    LBPROC_COMP(12)=1
!
        IF((STLIST(st_proc_no_code,ilprev) == st_accum_code) .OR.       &
     &     (STLIST(st_proc_no_code,ilprev) == st_time_mean_code).OR.    &
     &     (STLIST(st_proc_no_code,ilprev) == st_time_series_mean))     &
     &    LBPROC_COMP(8)=1
!
        IF(STLIST(st_proc_no_code,ilprev) == st_min_code)               &
     &    LBPROC_COMP(13)=1
!
        IF(STLIST(st_proc_no_code,ilprev) == st_max_code)               &
     &    LBPROC_COMP(14)=1
!
        output_code=STLIST(st_output_code,IL)
!L
!L 4.1 OUTPUT to PPfile
!L
        IF (output_code <  0) THEN                  ! PP Output
! Find appropriate dump header if a daughter record
          IF (il /= ilcurr) then
            icurrll_dump_ptr=stlist(st_lookup_ptr,ilcurr)
          ENDIF
!L
!L 4.1.0 Determine output PP unit and associated filename; OPEN file
!L
! If preattached files are used the file is left open by PPCTL following
! the initial OPEN; if reinitialised files are used the unit must be
! OPENed and CLOSEd explicitly every time it is used.
!
          UNITPP=-output_code

! attach the lookup table
            IF (mype == current_io_pe) THEN
               CALL ATTACH_IPPLOOK(IPPLOOK, UNITPP)
            Else
              IWA = IMDI
              IWL = IMDI
            END IF

          IF (FT_STEPS(UNITPP) /= 0) THEN ! Filename generated by model
! Check if re-initialised file stream has been opened yet
            IF(STEP <  FT_FIRSTSTEP(UNITPP)) THEN ! File stream opened?
               ICODE=1
               CMESSAGE='STWORK  : Re-initialised file not yet created'
               write(6,*)                                               &
     &          'STWORK  : FATAL ERROR. Attempt to write to ',          &
     &          're-initialised file stream before file first opened:'
               write(6,*)                                               &
     &          '        : Check that output on unit ',unitpp,' is not',&
     &          ' requested before first initialisation of output file:'
               write(6,*)                                               &
     &          '        :  See UMUI window (Initialisation of PP file',&
     &          's) accessed from (Post Processing) from (Submodel ',   &
     &          'independent).'
               GO TO 999
             ENDIF                                ! File stream opened?
            STRING=MODEL_FT_UNIT(UNITPP)
            DO JJ=80,1,-1
              IF (STRING(JJ:JJ) == '/') GOTO 411
            ENDDO
            ICODE=1
            CMESSAGE='STWORK  : Illegal output PPfile name'
            GOTO 999
 411        CONTINUE
            IF (JJ >  66) THEN
              PPNAME=STRING(JJ+1:80)
            ELSE
              PPNAME=STRING(JJ+1:JJ+14)
            ENDIF
            LEN_PPNAME=LEN(PPNAME)
! To avoid unnecessary open/close/read/write's, only
! open this file once more, after the initialisation.
!
            CALL IS_UNIT_OPEN(UNITPP, OPEN_FLAG)
            CALL GC_IBCAST(1008, 1, current_io_pe, nproc, istat,        &
     &                     open_flag)

            IF (OPEN_FLAG == 1) THEN   ! file is closed - open it
! DEPENDS ON: file_open
              CALL FILE_OPEN(UNITPP,PPNAME,LEN_PPNAME,1,1,ICODE)
            END IF
          ENDIF
!L
!L 4.1.1 Read in the pp fixed-length header
!L
          CALL ATTACH_FXH(PP_FIXHD,UNITPP)

!L
!L 4.1.2 Find the first available pp lookup record.
!L
          ICURRLL=FT_LASTFIELD(UNITPP) ! Position of the last field
          ICURRLL=ICURRLL+1            ! Position of the next field
!L
!L 4.1.3 Find the first available position for the next data record(s)
!L       by reading last pp lookup record.
!L
          IWL= PP_FIXHD(150)-1  ! NOTE for BUFFIN I/O the start address
!                               ! is zero for word 1. This is pointer
!                               ! to start of lookups.
          IF(ICURRLL == 1) THEN      ! First record
            IWA=PP_FIXHD(160)-1   ! Pointer to start of data

          ELSE
! Pointer to next available data location in output file
           IF (mype == current_io_pe) THEN
              IWA = IPPLOOK(LBEGIN, ICURRLL-1)+                         &
     &              IPPLOOK(LBNREC, ICURRLL-1)
           ELSE
              IWA = IMDI
           END IF


          ENDIF                     ! Test on first record
!L
!L 4.1.4 If a daughter record is being processed then recover
!L         size information from dump LOOKUP header referenced by
!L         mother record (unless MOS output)
!L
          IF (il /= ilcurr) THEN
            extraw_hdr=lookup(lbext,icurrll_dump_ptr)

            global_pphoriz_out=lookup(lblrec,icurrll_dump_ptr)
            global_n_rows_out=lookup(lbrow,icurrll_dump_ptr)
            global_n_cols_out=lookup(lbnpt,icurrll_dump_ptr)
            IF (what_proc == st_time_series_mean) then
! As work is done on only PE 0 and copy to buf 3 uses pphoriz_out
! this must be reset.
              pphoriz_out=global_pphoriz_out
              n_rows_out=global_n_rows_out
              n_cols_out=global_n_cols_out
            endif
            IF (what_proc  ==  st_time_series_code) THEN
              pphoriz_out=global_pphoriz_out
              n_rows_out=global_n_rows_out
              n_cols_out=global_n_cols_out
            ENDIF

          ENDIF
!L
!L 4.1.5 Check PP_PACK_CODE for GRIB output. Set GRIB flag and reset
!L       PP_PACK_CODE to give packing profile.
       IF(PP_PACK_CODE(UNITPP) >= 100)THEN
         GRIB_OUT=.TRUE.
         PP_PACK_CODE(UNITPP)=PP_PACK_CODE(UNITPP)-100
         GRIB_PACKING=PP_PACK_CODE(UNITPP)
       ELSE
         GRIB_OUT=.FALSE.
       ENDIF
!L
!L 4.1.6 Set packing accuracy for output data field and buffer length
!L       Multiple packing profiles are held in STASHmaster and chosen on
!L       a per-unit basis through PP_PACK_CODE.  Profile 0 means
!L       unpacked. If the field has any extra data switch off packing.
!L
          IF (PP_PACK_CODE(UNITPP) == 0.OR.extraw_hdr /= 0) THEN
            PACKING=.FALSE.
            COMP_ACCRCY=-99
          ELSE
            PACKING=.TRUE.
! DEPENDS ON: exppxi
      comp_accrcy= EXPPXI( im_ident, IS, IM,                            &
     &                     ppx_packing_acc+PP_PACK_CODE(UNITPP)-1,      &
#include "argppx.h"
     &             icode, cmessage)
          ENDIF
          IF(GRIB_OUT) THEN  ! reset packing code
            PP_PACK_CODE(UNITPP)=PP_PACK_CODE(UNITPP)+100
          ENDIF

          LENBUF=((global_PPHORIZ_OUT+um_sector_size-1)/um_sector_size)*&
     &     um_sector_size
!        ! Output length before pack

!L
!L 4.2 Select routine to output data using logical GRIB.
!L     If data to be output in grib code then call GRIB_FILE
!L     If data to be output in PP   code then call PP_FILE
!L
        DO  II=1,NUM_LEVS_OUT           ! --- Start of levels loop ---

! Gather together distributed field to PE 0
! Distributed data is in PP_FIELD, gathered data will be
! in the buf array

          IF ( (what_proc == st_time_series_code) .OR.                  &
     &    (what_proc == st_time_series_mean)  ) THEN
! if it's either MOS or timeseries output - just copy on PE 0

            IF (mype  ==  0) THEN
              DO I=1,PPHORIZ_OUT
                buf(I)=PPFIELD(I+(II-1)*PPHORIZ_OUT)
              ENDDO
            ENDIF

          ELSE ! not timeseries output

            packing_type=0
            num_words=0
            num_out=0
! DEPENDS ON: stash_gather_field
            CALL STASH_GATHER_FIELD (                                   &
     &        PPFIELD(1+(II-1)*PPHORIZ_OUT) , buf,                      &
     &        PPHORIZ_OUT , global_PPHORIZ_OUT , 1,                     &
     &        stlist(st_south_code,il)+global_N_ROWS_OUT-1 ,            &
     &        stlist(st_west_code,il)+global_N_COLS_OUT-1,              &
     &        stlist(st_south_code,il),                                 &
     &        stlist(st_west_code,il),                                  &
     &        grid_type_code,halo_type,current_io_pe,.TRUE.,            &
     &        PACKING, IM_IDENT, LRLE, PACKING_TYPE,                    &
     &        NUM_OUT,                                                  &
     &        COMP_ACCRCY, RMDI,                                        &
     &        ICODE=ICODE, CMESSAGE=CMESSAGE)

            IF (ICODE  >   0) THEN
              WRITE(6,*) 'Error occured in STASH while gathering ',     &
     &                   'data for output.'
              GOTO 999
            ENDIF

          ENDIF 

!  Reset index for PP_HEAD if there is a levels list of hybrid levels
          IF (stlist(st_output_bottom,il) <  0.and.                     &
     &         lbvcl  ==  ppx_lbvc_hybrid) THEN
            JJ = level_list(II)
!  or a range of model levels, as they may not be consecutive,
          ELSE IF (stlist(st_output_bottom,il) >= 1.and.                &
     &             stlist(st_output_top,il) <= T_levels) THEN
            JJ = level_list(II)
!  otherwise use of level index is alright
          ELSE
            JJ = II
          END IF
!
!    Check that PP output file has sufficient headers pre-allocated
       IF(ICURRLL >  PP_FIXHD(152)) THEN
          ICODE=4
          WRITE(6,*) 'ERROR detected in routine STWORK: stop model'
          WRITE(6,*) ': No. of output fields (=',ICURRLL,')',           &
     &          ' exceeds no. of reserved PP headers for unit ',UNITPP
          CMESSAGE='STWORK   : NO. OF FIELDS EXCEEDS RESERVED HEADERS'
          GOTO 999
       ENDIF     ! end  no. of pp fields check
!
       IF(GRIB_OUT) THEN
!
! NOTE cannot pack data into grib before pphead correctly setup
!
         NUM_WORDS = -99   ! ie unset before call to pp_head
         PACKING_TYPE=3
       ELSE
!L Pack data into PP code.


!
! Set the normal default lengths for no packing
!
          packing_hold=packing
          packing_type_hold=packing_type
          len_buf_words=lenbuf
          len_field=global_PPHORIZ_OUT
!
! Check if we have already packed this data, and if so
! record the packing flag, and turn off packing
!
          if(packing.and.packing_type == 1) then
!
            packing=.false.
            packing_type=0
!
            num_words=(num_out+1)/2 ! Round up to the nearest 64 Bit
                                    ! CRAY Words
            len_buf_words=((num_words+um_sector_size-1)/um_sector_size)*&
     &       um_sector_size
            len_field=num_words
!
          endif ! packing.and.packing_type == 1
! DEPENDS ON: pp_file
          CALL PP_FILE(buf,                                             &
     &      len_buf_words, num_words, rmdi, comp_accrcy,                &
     &      len_field, unitpp, iwa,                                     &
     &      global_N_COLS_OUT,global_N_ROWS_OUT,                        &
     &     PACKING,im_ident,LRLE,                                       &
     &     PACKING_TYPE,current_io_pe,EXTRA_VAR_DATA                    &
     &    ,SROW_OUT,WCOL_OUT,ICODE,CMESSAGE)

!
! Restore the packing flags
!
        if(packing_type == 0)then
          packing=packing_hold
          packing_type=packing_type_hold
        endif
! Make sure all processors get the return code

          CALL GC_IBCAST(101,1,0,nproc,info,ICODE)


! Num_words is the no of 64 bit words required
          IF(ICODE >  0) THEN
            GOTO 999
          ENDIF
       ENDIF
      ! Add any EXTRA data concerning variable grids to the extra data
      ! for lookup header

          LEN_BUF_WORDS=((NUM_WORDS+um_sector_size-1)/um_sector_size)*  &
     &     um_sector_size ! No of words output
!L
!L 4.2.1 Set STASH processing codes and sampling period for PPheader
!L
          GR=STLIST(st_gridpoint_code,ILPREV)! Grid point code
! Any time-processed field has a (non-zero) sample_prd set -
! this will be translated by PP_HEAD into an LBTIM subcode IB of 2
          sample_prd=0.0
          IF(STLIST(st_proc_no_code,ILPREV) >  st_replace_code) THEN
            sample_prd=REAL(STLIST(st_freq_code,ILPREV)*secs_per_period)&
     &                /REAL(steps_per_period*3600)
          ENDIF
!L
!L 4.2.2 Verification time comes from fixhd(28), current time fixhd(21)
!L       2 cases that require consideration here:
!L
!L      (1) this record is not a daughter record.
!L          in which case, set start_step=.true., verif time from fixhd
!L          present time also from fixhd
!L
!L      (2) this record IS a daughter record.
!L          in which case, will need to retreive info on start_time
!L          from dump
!L
          start_step=.true.
! We need to protect the call
          IF (mype == current_io_pe) THEN
          IF (il == ilcurr) THEN ! not daughter record
! DEPENDS ON: pp_head
            CALL PP_HEAD(                                               &
#include "argppx.h"
     &        im_ident,FIXHD,INTHD,REALHD,                              &
     &        LEN_INTHD,LEN_REALHD,                                     &
     &        IM,IS,GR,lfullfield,                                      &
     &        level(II),pseudo_level(II),                               &
     &        samples,start_step,fixhd(28),fixhd(21),LEN1_LOOKUP,       &
     &        extraw_hdr,IPPLOOK(1,ICURRLL),IPPLOOK(1,ICURRLL),         &
     &        global_N_COLS_OUT,NUM_WORDS,LEN_BUF_WORDS,                &
     &        global_N_ROWS_OUT,NROW_IN,SROW_IN,WCOL_IN,ECOL_IN,        &
     &        lbproc_comp,sample_prd,                                   &
     &        FCST_PRD,COMP_ACCRCY,PACKING_TYPE,                        &
     &        st_grid,IWA,                                              &
! superarray 1-4: zseak_rho,Ck_rho,zseak_theta,Ck_theta
     &        stsuparr(sa_idx(1)),stsuparr(sa_idx(2)),                  &
     &        stsuparr(sa_idx(3)),stsuparr(sa_idx(4)),                  &
     &        T_levels,JJ,ROTATE,ELF,                                   &
     &        OCEAN,LEVDEPC,LEN1_LEVDEPC,                               &
     &        ICODE,CMESSAGE)
          ELSE ! daughter record so start time is in dump
! set up start_time from data in LOOKUP(lbyr,icurrll_dump_ptr)
            start_time(1)=LOOKUP(lbyr,icurrll_dump_ptr)
            start_time(2)=LOOKUP(lbmon,icurrll_dump_ptr)
            start_time(3)=LOOKUP(lbdat,icurrll_dump_ptr)
            start_time(4)=LOOKUP(lbhr,icurrll_dump_ptr)
            start_time(5)=LOOKUP(lbmin,icurrll_dump_ptr)
            start_time(7)=LOOKUP(lbday,icurrll_dump_ptr)
! DEPENDS ON: pp_head
            CALL PP_HEAD(                                               &
#include "argppx.h"
     &        im_ident,FIXHD,INTHD,REALHD,                              &
     &        LEN_INTHD,LEN_REALHD,                                     &
     &        IM,IS,GR,lfullfield,                                      &
     &        level(II),pseudo_level(II),                               &
     &        samples,start_step,start_time,                            &
     &        fixhd(28),LEN1_LOOKUP,                                    &
     &        extraw_hdr,IPPLOOK(1,ICURRLL),IPPLOOK(1,ICURRLL),         &
     &        global_N_COLS_OUT,NUM_WORDS,LEN_BUF_WORDS,                &
     &        global_N_ROWS_OUT,NROW_IN,SROW_IN,WCOL_IN,ECOL_IN,        &
     &        lbproc_comp,sample_prd,                                   &
     &        FCST_PRD,COMP_ACCRCY,PACKING_TYPE,                        &
     &        st_grid,IWA,                                              &
! superarray 1-4: zseak_rho,Ck_rho,zseak_theta,Ck_theta
     &        stsuparr(sa_idx(1)),stsuparr(sa_idx(2)),                  &
     &        stsuparr(sa_idx(3)),stsuparr(sa_idx(4)),                  &
     &        T_levels,JJ,                                              &
     &        ROTATE,ELF,OCEAN,LEVDEPC,LEN1_LEVDEPC,                    &
     &        ICODE,CMESSAGE)
           ENDIF ! end if block over daughter/mother record

           IF(ICODE >  0) GOTO 999 ! An error has occured
          END IF ! mype == current_io_pe


           IF(GRIB_OUT) THEN
!  Now safe to call grib coder as pphead correctly set apart from
!  length of data
!L Pack data into grib code

             IF (mype  ==  0) THEN
! DEPENDS ON: grib_file
               CALL GRIB_FILE(LEN1_LOOKUP,PP_LEN2_LOOKUP,               &
     &                  IPPLOOK,IPPLOOK,ICURRLL,                        &
     &                  buf,global_PPHORIZ_OUT,                         &
     &                  LENBUF,NUM_WORDS,UNITPP,IWA,GRIB_PACKING,       &
     &                  ICODE,CMESSAGE)
             ENDIF  ! (IF mype == 0)

! Make sure all processors get the return code

          CALL GC_IBCAST(101,1,0,nproc,info,ICODE)

             IF(ICODE >  0)THEN
               CMESSAGE='STWORK  : Error in GRIB_FILE'
               GOTO 990
             ENDIF

           ENDIF       ! end of grib_out

           IF (mype == current_io_pe) THEN
             IWA = IPPLOOK(LBEGIN,ICURRLL) + IPPLOOK(LBNREC,ICURRLL)
           ELSE
             IWA = IMDI
           END IF
           ICURRLL=ICURRLL+1  ! Update the counter
           icurrll_dump_ptr=icurrll_dump_ptr+1
! strictly only needs doing if a daughter record
!
         ENDDO                           ! --- End of levels loop ---
!
         FT_LASTFIELD(UNITPP)=ICURRLL-1  ! Position of the last field
!
        ELSEIF (output_code == st_dump.OR.output_code == st_secondary)  &
     &  THEN
!L
!L 4.4 OUTPUT to dump or secondary D1 space - this implies some
!L     TEMPORAL processing possibly.  If destination is secondary D1
!L     space, there will be no associated LOOKUP header to update.
!L
! Length is calculated from STASHlist
! NB: Full field length must be coded here, even for partial timeseries
!
          NUM_WORDS=STLIST(st_dump_output_length,IL)/NUM_LEVS_OUT
          ICURRLL=STLIST(st_lookup_ptr,IL) ! Location of DUMP header
!
          DO  II=1,NUM_LEVS_OUT          ! --- Start of levels loop ---
!  Reset index for PP_HEAD if there is a levels list of hybrid levels
                IF (stlist(st_output_bottom,il)  <   0 .AND.            &
     &              lbvcl == ppx_lbvc_hybrid) THEN
              JJ = level_list(II)
!  or a range of model levels, as they may not be consecutive,
            ELSE IF (stlist(st_output_bottom,il) >= 1.and.              &
     &               stlist(st_output_top,il) <= T_levels) THEN
              JJ = level_list(II)
!  otherwise use of level index is alright
            ELSE
              JJ = II
            END IF
            addr=stlist(st_output_addr,il) ! start address
            IF (what_proc == st_time_series_code.or.                    &
     &          what_proc == st_time_series_mean) THEN

!L
!L 4.4.1 Timeseries addresses are incremented according to timestep
!L
              IF (stlist(st_freq_code,il) <  1) THEN
                icode=st_not_supported
                cmessage=                                               &
     &              'STWORK  : STASHtime for timeseries not supported'
                goto 999 ! got an error so jump to return
              ENDIF
              elap_time=step-stlist(st_start_time_code,il)
              elap_time=(mod(elap_time,stlist(st_period_code,il)))/     &
     &          stlist(st_freq_code,il)
              addr=addr+(elap_time*pphoriz_out)
!L on the first time step of a timeseries processing
!L pphoriz_out is the length of the entire output vector
!L including extra data -- on other timesteps it is
!L the length of a single record (data for just one timestep)
            ENDIF
!L
!L 4.4.2 TEMPORAL processing from ppfield array to D1
!L
! DEPENDS ON: temporal
            CALL TEMPORAL(ppfield(1+(ii-1)*pphoriz_out),                &
     &        d1(addr+(ii-1)*pphoriz_out),pphoriz_out,extraw,           &
     &        stlist(1,il),len_stlist,OCEAN,step,                       &
     &        icode,cmessage,start_step,rmdi)

            IF (icode >  0) goto 999
!L
!L 4.4.3 Set up LOOKUP header if destination is main part of D1
!L
            IF (output_code == st_dump) THEN
!L
!L 4.4.3 Set other information for input to PPHEAD
!L
              GR=STLIST(st_gridpoint_code,ILPREV)! Grid point code
! Any time-processed field has a (non-zero) sample_prd set -
! this will be translated by PP_HEAD into an LBTIM subcode IB of 2
              sample_prd=0.0
              IF (STLIST(st_proc_no_code,ILPREV) >  st_replace_code)    &
     &        THEN
                sample_prd=REAL(STLIST(st_freq_code,ILPREV)*            &
     &                     secs_per_period)/REAL(steps_per_period*3600)
              ENDIF

! Address of whole field is calculated from STASHlist
              IWA=STLIST(st_dump_output_addr,IL)+(II-1)*NUM_WORDS

!L
!L 4.4.4 Call PPHEAD to set LOOKUP header for field STASHed to D1.
!L       Here pass previous_time as well as start_step from temporal
!L       if start_step is true then start time will be updated.
!L       Value of end time is unimportant as that is handled properly
!L       when data is written out to pp file.
!L       Note that LBNREC is hardwired to 0 and so too is BACC.
!L
         

! DEPENDS ON: pp_head
            CALL PP_HEAD(                                               &
#include "argppx.h"
     &        im_ident,FIXHD,INTHD,REALHD,                              &
     &          LEN_INTHD,LEN_REALHD,                                   &
     &          IM,IS,GR,lfullfield,                                    &
     &        level(II),pseudo_level(II),                               &
     &          samples,start_step,previous_time,fixhd(28),LEN1_LOOKUP, &
     &          extraw_hdr,LOOKUP(1,ICURRLL),RLOOKUP(1,ICURRLL),        &
     &          global_N_COLS_OUT,NUM_WORDS,0,                          &
     &          global_N_ROWS_OUT,NROW_IN,SROW_IN,WCOL_IN,ECOL_IN,      &
     &          lbproc_comp,sample_prd,                                 &
     &          FCST_PRD,0,PACKING_TYPE,                                &
     &        st_grid,IWA,                                              &
! superarray 1-4: zseak_rho,Ck_rho,zseak_theta,Ck_theta
     &        stsuparr(sa_idx(1)),stsuparr(sa_idx(2)),                  &
     &        stsuparr(sa_idx(3)),stsuparr(sa_idx(4)),                  &
     &        T_levels,JJ,                                              &
     &          ROTATE,ELF,OCEAN,LEVDEPC,LEN1_LEVDEPC,                  &
     &          ICODE,CMESSAGE)
!
         
              IF(ICODE >  0) goto 999 ! An error has occured

! Only (optionally) pack fields if no extra words of data
              IF (extraw_hdr  ==  0) THEN
                LOOKUP(LBPACK,ICURRLL) =                                &
! DEPENDS ON: exppxi
     &            EXPPXI( im_ident, IS, IM, ppx_dump_packing,           &
#include "argppx.h"
     &             icode, cmessage)
                IF (DUMP_PACK == 3 ) THEN
!                 Override packing indicator from PPXREF
                  N1 = 0   !   No packing
                  LOOKUP(LBPACK,ICURRLL) =                              &
     &            (LOOKUP(LBPACK,ICURRLL)/10)*10 + N1
                ENDIF
              ELSE
                LOOKUP(LBPACK,ICURRLL)=0
              ENDIF
! Set data type (REAL/INTEGER) from STASHmaster (-ve for timeseries)
              IF (STLIST(st_series_ptr,ilprev) >  0.OR.                 &
     &   STLIST(st_proc_no_code,ilprev) == st_time_series_mean) THEN
                LOOKUP(DATA_TYPE,ICURRLL) =                             &
! DEPENDS ON: exppxi
     &           -EXPPXI( im_ident, IS, IM, ppx_data_type,              &
#include "argppx.h"
     &             icode, cmessage)
              ELSE
                LOOKUP(DATA_TYPE,ICURRLL) =                             &
! DEPENDS ON: exppxi
     &            EXPPXI( im_ident, IS, IM, ppx_data_type,              &
#include "argppx.h"
     &             icode, cmessage)
              ENDIF
              ICURRLL=ICURRLL+1 ! Update the counter for the next field
!
            ENDIF               ! End of IF output_code=dump
!
          ENDDO                           ! --- End of levels loop ---
!
        ELSE
          ICODE=9
          CMESSAGE='STWORK  : Illegal output destination in STLIST'
          GOTO 999
        ENDIF    ! End of STLIST output destination IF block
!
      ENDIF      ! END OF S_F IF Block ---------------------------------
!L
!L 5. End of loop over STASHlist entries - RETURN to calling routine
!L
      NULLIFY(PP_FIXHD)
  200 CONTINUE
!
  999 CONTINUE
      RETURN
!L
!L 9. IO error exits
!L
  990 WRITE(CMESSAGE,'("STWORK  : Error opening output PP file on unit "&
     &                ,I2)') UNITPP
      RETURN
      END SUBROUTINE STWORK
#endif
