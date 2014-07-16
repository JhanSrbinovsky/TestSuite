

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL SUBROUTINE TIMER ------------------------------------------------
!LL
!LL                    Purpose:
!LL Allows the recording of time spent in any section of the program
!LL Two types of timings are supported:
!LL non-inclusive : if a timed section of code (1) contains another
!LL                 timed section of code (2), then the timing for
!LL                 section (1) will not include the time spent in
!LL                 section (2). This is the normal use for the timer
!LL                 routine in the UM up to vn3.4
!LL inclusive     : allows the user to measure the time taken between
!LL                 any two points in the code, irrespective of any
!LL                 other calls to the timer routine within the timed
!LL                 section
!LL
!LL NB: Non-inclusive timers DO INCLUDE any inclusive timer sections
!LL     contained within them. If this section of code should not be
!LL     included, then also time it with a non-inclusive timer
!LL
!LL Timer now also records the time spent in itself
!LL Parameters:
!LL section_name - 20 byte character string containing name of
!LL                timer reference
!LL
!LL action:
!LL  1 -> first call to timer (timer initialisation)
!LL  2 -> last call to timer (prints out the collected data)
!LL  3 -> non-inclusive start timer
!LL  4 -> non-inclusive end timer
!LL  5 -> inclusive start timer
!LL  6 -> inclusive end timer
!LL
!LL Timer should be called with action=1 before the first executable
!LL statement, and with action=2 after the last executable statement.
!LL
!LL   Model               Modification History
!LL  version    Date
!LL   4.1       16/08/94  Based on old UM timer routine
!LL   4.2       08/10/96  Corrected intermediate timer error in
!LL                       elapsed times.
!LL                       Corrected size of message arg. in
!LL                       TIMER_OUTPUT.
!LL                       P.Burton
!LL   4.3       23/01/97  Added better overview for mpp runs
!LL                       Corrected T3E wallclock time calculation
!LL                       P.Burton
!LL   4.5       17/04/98  Added barrier to allow imbalance to be
!LL                       included in the correct timer section
!LL                                                     P.Burton
!LL   4.5       09/07/98  Replaced missing array index for
!LL                       last_ni_wallclock_time_elapsed.
!LL                                                     P. Burton
!LL   5.1       10/04/00  New reconfiguration support. P.Selwood.
!LL   5.2       06/04/99  Moved the initial synchronisation in
!LL                       TIMER_OUTPUT to make sure that the
!LL                       data is not overwritten, if PE 0 is getting
!LL                       from the other PE's.
!LL                                           Author: Bob Carruthers
!LL   5.3       08/06/01  Move get_cpu_time and get_wallclock_time
!LL                       functions into timefn2a.  A van der Wal
!     5.3       26/07/01  Take an internal copy of "action" variable
!                         so as not to modify it (intent in). P.Selwood
!     6.0       12/09/03  Add FLDCALC def. D Robinson
!     6.1       19/10/04  Add MAKEBC def, R. Sempers
!LL
!LL  Author : Paul Burton
!LL


!*********************************************************************


      SUBROUTINE TIMER_OUTPUT(                                          &
     &  in_number_of_timers, ni_number_of_timers,                       &
     &  in_cpu_time_elapsed, ni_cpu_time_elapsed,                       &
     &  in_wallclock_time_elapsed, ni_wallclock_time_elapsed,           &
     &  in_number_of_times_timed, ni_number_of_times_timed,             &
     &  in_timer_name, ni_timer_name,                                   &
     &  action,message)

      IMPLICIT NONE

! Arguments

      INTEGER                                                           &
     &  in_number_of_timers                                             &
                             ! IN number of inclusive timers
     &, ni_number_of_timers                                             &
                             ! IN number of non-inclusive timers
     &, in_number_of_times_timed(in_number_of_timers)                   &
!                            ! IN number of times timed - inclusive
     &, ni_number_of_times_timed(ni_number_of_timers)                   &
!                            ! IN number of times timed - non-incl.
     &, action  ! final output or intermediate

      REAL                                                              &
     &  in_cpu_time_elapsed(in_number_of_timers)                        &
!                            ! IN elapsed inclusive CPU time
     &, ni_cpu_time_elapsed(ni_number_of_timers)                        &
!                            ! IN elapsed non-inclusive CPU time
     &, in_wallclock_time_elapsed(in_number_of_timers)                  &
!                            ! IN elapsed inclusive wallclock time
     &, ni_wallclock_time_elapsed(ni_number_of_timers)

      CHARACTER*20                                                      &
     &  in_timer_name(in_number_of_timers)                              &
!                            ! IN name of timed section - inclusive
     &, ni_timer_name(ni_number_of_timers)
!                            ! IN name of timed section - non-incl.

      CHARACTER*(*)                                                     &
     &   message              ! IN message to print


! Local variables

      INTEGER max_timers
      PARAMETER(max_timers=300)

      INTEGER last_call_to_timer,intermediate_output
      PARAMETER(last_call_to_timer=2,                                   &
     &          intermediate_output=7)

      INTEGER                                                           &
     &  number_of_timers                                                &
     &, local_number_of_times_timed(max_timers)

      REAL                                                              &
     &  local_cpu_time_elapsed(max_timers)                              &
     &, local_wallclock_time_elapsed(max_timers)

      CHARACTER*20                                                      &
     &  local_timer_name(max_timers)

! Variables required for using intermediate timers
! They record the values on the last call to this routine
      INTEGER                                                           &
     &  last_in_number_of_times_timed(max_timers)                       &
     &, last_ni_number_of_times_timed(max_timers)

      REAL                                                              &
     &  last_in_cpu_time_elapsed(max_timers)                            &
     &, last_ni_cpu_time_elapsed(max_timers)                            &
     &, last_in_wallclock_time_elapsed(max_timers)                      &
     &, last_ni_wallclock_time_elapsed(max_timers)

      LOGICAL                                                           &
     &  first_intermediate_timer_call

      DATA first_intermediate_timer_call /.TRUE./
      DATA last_in_number_of_times_timed /max_timers*0/
      DATA last_ni_number_of_times_timed /max_timers*0/
      DATA last_in_cpu_time_elapsed /max_timers*0.0/
      DATA last_ni_cpu_time_elapsed /max_timers*0.0/
      DATA last_in_wallclock_time_elapsed /max_timers*0.0/
      DATA last_ni_wallclock_time_elapsed /max_timers*0.0/

      SAVE                                                              &
     &  last_in_number_of_times_timed, last_ni_number_of_times_timed,   &
     &  last_in_cpu_time_elapsed, last_ni_cpu_time_elapsed,             &
     &  last_in_wallclock_time_elapsed, last_ni_wallclock_time_elapsed, &
     &  first_intermediate_timer_call



      INTEGER sortwork_int    ! work variable for sort
      REAL    sortwork_real   ! work variable for sort
      CHARACTER*20 sortwork_char ! work variable for sort

      REAL total_cpu_time,                                              &
                                   ! total cpu time spent in program
     &     total_wallclock_time,                                        &
                                   ! total wallclock time spent in
!                                  ! program
     &     average_cpu_elapsed,                                         &
                                   ! average cpu elapsed time
     &     average_wallclock_elapsed,                                   &
                                      ! average wallclock elapsed time
     &     percent_of_cpu_total,                                        &
                                   ! % of cpu time spent in a section
     &     percent_of_wallclock_total,                                  &
                                  ! % of wallclock time spent in a
!                                 ! section
     &     speed_up               ! speed_up=cpu/wallclock

! ------------------------ Comdeck PARVARS -------------------------
! Parameters and common blocks required by the mpp-UM
! Added new comdeck AMAXSIZE required for new arrays in PARCOMM
! Add non-mpp option
!                                                      P.Burton
!========================== COMDECK PARPARM ====================
!   Description:
!
!   This COMDECK contains PARAMETERs for the mpp-UM
!
!   Two sets of parameters are set up -
!     i)  for the mpp-UM itself.
!     ii) for the interface to the Message Passing Software.
!
      !=================================================================
      ! Parameters needed for the mpp-UM
      !=================================================================
      ! maximum number of spatial dimensions
      INTEGER,PARAMETER:: Ndim_max = 3 ! 3d data

      ! number of different halo types
      INTEGER,PARAMETER:: NHalo_max = 3 ! for N.D. atmos. model

      INTEGER,PARAMETER:: halo_type_single   = 1
      INTEGER,PARAMETER:: halo_type_extended = 2
      INTEGER,PARAMETER:: halo_type_no_halo  = 3

! FLDTYPE definitions for the different field types recognised on the
! decomposition
      INTEGER,PARAMETER:: Nfld_max=7 ! maximum number of field types
      INTEGER,PARAMETER:: fld_type_p=1       ! grid on P points
      INTEGER,PARAMETER:: fld_type_u=2       ! grid on U points
      INTEGER,PARAMETER:: fld_type_v=3       ! grid on V points
      INTEGER,PARAMETER:: fld_type_comp_wave  = 4
                              ! Compressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_full_wave  = 5
                              ! Uncompressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_rim_wave   = 6
                              ! Boundary data for WAM Wave Field
      INTEGER,PARAMETER:: fld_type_r=7       ! grid on river points
      INTEGER,PARAMETER:: fld_type_unknown=-1! non-standard grid
! FLDTYPE end

      ! Used in addressing to indicate if calculation is for a local or
      ! global (ie. disk dump) size

      INTEGER,PARAMETER:: local_data=1
      INTEGER,PARAMETER:: global_dump_data=2

      ! maximum permitted size of a halo
      INTEGER,PARAMETER:: Max_Halo_Size=10

      !=================================================================
      ! Parameters needed for the Message Passing Software
      !=================================================================
      INTEGER,PARAMETER:: Maxproc = 512 ! Max number of processors

      ! Processor addresses in the neighbour array
      INTEGER,PARAMETER:: PNorth   = 1
      INTEGER,PARAMETER:: PEast    = 2
      INTEGER,PARAMETER:: PSouth   = 3
      INTEGER,PARAMETER:: PWest    = 4

      ! Value in neighbour array if the domain has  no neighbour in this
      ! direction. Otherwise the value will be the tid of the neighbor
      INTEGER,PARAMETER:: NoDomain = -1

      INTEGER,PARAMETER:: BC_STATIC   = 1 ! Static boundary conditions
      INTEGER,PARAMETER:: BC_CYCLIC   = 2 ! Cyclic boundary conditions
      INTEGER,PARAMETER:: BC_OVERPOLE = 3 ! Transfer over pole
! PARPARM end
!====================== COMDECK AMAXSIZE ========================
! Description
!   This comdeck provides parameters giving the maximum likely sizes
!   of key UM resolution variables, useful for sizing static arrays.
!
!   History:
!   Model    Date     Modification history
!  version
!   4.2      18/11/96 New comdeck created.  P.Burton
!   4.3      24/01/97 Define MaxFieldSize to be a quarter of the
!                     SHMEM common block size.         P.Burton
!   4.4      3/7/97   Add MaxFieldSizeMes. Deborah Salmond
!   4.5     12/01/98  Added new variables, and changed sizes to
!                     correspond to global hi-res forecast - current
!                     largest configuration.                P.Burton
!                     Changed MAX_SHMEM_COMMON_SIZE to 3000000
!                     required for operational data assimilation.
!                                                           P.Burton
!   5.0     29/04/99  Changed variable names:
!                       P_ROWS_MAX -> ROWS_MAX
!                       P_LEVELS_MAX -> MODEL_LEVELS_MAX
!                       Q_LEVELS_MAX -> WET_LEVELS_MAX
!                       MaxHaloSize -> MaxHaloArea
!                     Removed variable:
!                       HALO_MAX (use PARPARM Max_Halo_Size instead)
!    5.0   29/04/99  Remove mpp #define
!    5.3   05/12/01  Remove MaxFieldSize, MaxFieldSizeMes and
!                    Max3DFieldSize.  S.Cusack
!    5.5   22/01/03  Increase ROW_LENGTH_MAX and HORIZ_DIM_MAX
!                    from 432 to 548. D Robinson.
!    6.1   31/08/04  Allow up to 100 levels.  R.Barnes
!    6.2   13/02/06  Increase max values of row_length and
!                    rows to cope with FOAM high res, as well
!                    as Global N320 and NAE.  M Martin.
!    6.2   24/11/05  Use max function for horiz_dim_max. R Barnes
!    6.2     11/01/06 Remove max_shmem_common_size here and
!                     in rdobs2.   Camilla Mathison/R Barnes
!

! Maximum sector size for I/O
      INTEGER,PARAMETER:: IO_SECTOR_SIZE_MAX=4096
      INTEGER,PARAMETER:: ROW_LENGTH_MAX   = 840 ! Maximum row length
      INTEGER,PARAMETER:: ROWS_MAX         = 600 ! Max no of rows

      ! MAX(ROW_LENGTH_MAX,ROWS_MAX)
      INTEGER,PARAMETER:: HORIZ_DIM_MAX=MAX(ROW_LENGTH_MAX,ROWS_MAX)

      INTEGER,PARAMETER:: MODEL_LEVELS_MAX = 100 ! Max no of total levels
      INTEGER,PARAMETER:: WET_LEVELS_MAX   = 100 ! Max no of wet levels
      INTEGER, PARAMETER :: Max2DFieldSize = ROW_LENGTH_MAX*ROWS_MAX +  &
     &  IO_SECTOR_SIZE_MAX
      INTEGER, PARAMETER :: MaxHaloArea    = HORIZ_DIM_MAX*Max_Halo_Size
!========================== COMDECK PARCOMM ====================
!
! *** NOTE : This comdeck requires comdeck PARPARM to be *CALLed
!            first.
!
!   Description:
!
!   This COMDECK contains COMMON blocks for the mpp-UM
!
!
!   Two COMMON blocks are defined:
!     i)  UM_PARVAR holds information required by the
!         Parallel Unified Model itself
!     ii) MP_PARVAR holds information required by the interface to
!         the Message Passing Software used by the PUM
!
!   Key concepts used in the inline documentation are:
!     o global data - the entire data domain processed by the UM
!     o LOCAL data - the fragment of the global data which is
!       stored by this particular process
!     o PERSONAL data - the fragment of the LOCAL data which is
!       updated by this particular process
!     o HALO data - a halo around the PERSONAL data which forms
!       the LOCAL data
!
!     Acronyms used:
!     LPG - Logical Process Grid, this is the grid of logical
!           processors; each logical processor handles one of the
!           decomposed parts of the global data. It does not
!           necessarily represent a physical grid of processors.
!
!   History:
!
!   4.1      27/1/96  New comdeck based on second section of
!                     old PARVARS.   P.Burton
!   4.2     19/08/96  Removed some unused variables, and added
!                     current_decomp_type variable to allow use
!                     of flexible decompositions.
!                     Added nproc_max to indicate the max. number
!                     of processors used for mpp-UM
!                                                      P.Burton
!   5.0     12/04/99  - Added halosize array to common block
!                     - Added halo_i and halo_j to common block
!                     - Added fld_type dimension to glsize
!                     - Added halo type dimension to lasize
!                     - Added fld_type dimension to lasize
!                     - Replaced blsizep/blsizeu by blsize with
!                       extra fld_type dimension
!                     - Replace attop etc. with at_extremity
!                     - Added g_pe_index to common block
!                                                      P.Burton
!   5.1     22/05/00  Removed DATA statement and put in BLKDATA
!                                                      P.Burton
!   5.1     26/01/00  - Renamed g_pe_index -> g_pe_index_EW
!                     - Added g_pe_index_NS
!                                                     P.Burton
!   5.2     02/08/00  Added g_at_extremity        P.Burton
!   5.3     14/09/01  Added sb_model_domain   P.Burton
!   5.5     06/08/00  Modification for parallelisation of WAM.
!                     Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!   5.5     30/01/03  Generalised datastart   P.Selwood.
!   5.5     07/02/03  SX now uses PARCOMM instead of SXCOMM    E.Leung
!   6.0     18/09/03  F90-fy continuation lines.               P.Dando
!   6.2     23/11/05  Removed all references to the wavemodel.
!                     T.Edwards
!
! -------------------- COMMON BLOCKS --------------------
!
! =======================================================
! Common block for the Parallel Unified Model
! =======================================================

      INTEGER :: first_comp_pe       ! top left pe in LPG
      INTEGER :: last_comp_pe        ! bottom right pe in LPG
      INTEGER :: current_decomp_type ! current decomposition type
      INTEGER :: Offx                ! standard halo size in East-West
      INTEGER :: Offy                ! standard halo size in North-South
      INTEGER :: halo_i              ! extended halo size in East-West
      INTEGER :: halo_j              ! extended halo size in North-South
      INTEGER :: halosize(Ndim_max,NHalo_max) ! available halo sizes
      INTEGER :: glsize(Ndim_max,Nfld_max) ! global data size
      INTEGER :: lasize(Ndim_max,Nfld_max,NHalo_max) ! local data size
      INTEGER :: blsize(Ndim_max,Nfld_max) ! personal data size

      ! position of personal data in global data (in terms of standard
      ! Fortran array notation
      INTEGER :: datastart(Ndim_max)

      ! Generalised version of datastart for *all* fieldtypes
      INTEGER :: datastart_f(Ndim_max,Nfld_max)

      INTEGER :: gridsize(Ndim_max)  ! size of the LPG in each dimension

      ! position of this process in the LPG 0,1,2,...,nproc_x-1 etc.
      INTEGER :: gridpos(Ndim_max)

      INTEGER :: sb_model_domain


      ! logicals indicating if a processor is at the edge of the LPG
      LOGICAL :: at_extremity(4)

      COMMON /UM_PARVAR/                                                &
     &  first_comp_pe, last_comp_pe, current_decomp_type, Offx, Offy,   &
     &  halo_i, halo_j, halosize, glsize, lasize, blsize, datastart,    &
     &  datastart_f, gridsize, gridpos                                  &
     &,                 at_extremity,sb_model_domain

      ! Common block for the Message Passing Software

      ! type of boundary (cyclic or static) in each direction
      INTEGER :: bound(Ndim_max)

      ! global copy of local data size
      INTEGER :: g_lasize(Ndim_max,Nfld_max,NHalo_max,0:maxproc)

      ! global copy of personal data area
      INTEGER :: g_blsize(Ndim_max,Nfld_max,0:maxproc)

      ! global copy of datastart
      INTEGER :: g_datastart(Ndim_max,0:maxproc)

      ! global copy of datastart_f
      INTEGER :: g_datastart_f(Ndim_max,Nfld_max,0:maxproc)

      INTEGER :: g_gridpos(Ndim_max,0:maxproc) ! global copy of gridpos

      ! Which processor column a given point is in: 0 -> nproc_x-1
      INTEGER :: g_pe_index_EW(1-Max_Halo_Size:                         &
     &  ROW_LENGTH_MAX+Max_Halo_Size)

      ! Which processor row a given point is in: 0 -> nproc_y-1
      INTEGER :: g_pe_index_NS(1-Max_Halo_Size:ROWS_MAX+Max_Halo_Size)

      INTEGER :: nproc      ! number of processors in current decomp
      INTEGER :: mype      ! number of this processor (starting from 0)
      INTEGER :: nproc_max  ! maximum number of processors
      INTEGER :: nproc_x    ! number of processors in x-direction
      INTEGER :: nproc_y    ! number of processors in y-direction

      ! array with the tids of the four neighbours in the horizontal
      ! plane
      INTEGER :: neighbour(4)

      INTEGER :: gc_proc_row_group  ! GID for procs along a proc row
      INTEGER :: gc_proc_col_group  ! GID for procs along a proc col
      INTEGER :: gc_all_proc_group  ! GID for all procs

      ! at_extremity for each processor
      LOGICAL :: g_at_extremity(4,0:maxproc)

      COMMON /MP_PARVAR/                                                &
     &  bound,                                                          &
     &  g_lasize,g_blsize,g_datastart, g_datastart_f, g_gridpos,        &
     &  g_pe_index_EW,g_pe_index_NS,                                    &
     &  nproc_max,nproc_x,nproc_y,                                      &
     &  neighbour,gc_proc_row_group,                                    &
     &  gc_proc_col_group, gc_all_proc_group                            &
     &  ,nproc,mype                                                     &
     &, g_at_extremity

! PARCOMM end

! These are the declarations for mpp timer

      INTEGER info,                                                     &
     &  wallclock_max_pe(max_timers),wallclock_min_pe(max_timers),      &
     &  cpu_max_pe(max_timers),cpu_min_pe(max_timers)

      REAL wallclock_mean(max_timers),cpu_mean(max_timers),             &
     &     wallclock_median(max_timers),cpu_median(max_timers),         &
     &     wallclock_sd(max_timers),cpu_sd(max_timers),                 &
     &     wallclock_max(max_timers),wallclock_min(max_timers),         &
     &     cpu_max(max_timers),cpu_min(max_timers),                     &
     &     cpu_total(max_timers),speedup(max_timers),                   &
     &     efficiency(max_timers)

      INTEGER                                                           &
     &  summ_n_timers                                                   &
                         ! number of routines for ni summary
     &, routine_id          ! routine id on this processor

      REAL                                                              &
     &  wallclock_times(0:MAXPROC)                                      &
                                    ! wallclock time from each proc
     &, cpu_times(0:MAXPROC)                                            &
                                    ! cpu time from each proc
     &, total_cpu,max_wall       ! total cpu, maxumum wallclock times

      CHARACTER*20 summ_section(max_timers)  ! names of sections
      COMMON /MPP_TIMER/ summ_n_timers,total_cpu,max_wall,              &
     &                   wallclock_times,cpu_times,                     &
     &                   summ_section

! Variables for loops etc.
      INTEGER I,J,K,timer_kind

! Check to see if this is an intermediate output, and the first
! time it has been called
      IF ((action  ==  intermediate_output) .AND.                       &
     &    (first_intermediate_timer_call  ) ) THEN
! Copy the arguments into the last_* arrays
        first_intermediate_timer_call=.FALSE.

        DO I=1,in_number_of_timers
          last_in_number_of_times_timed(I)=in_number_of_times_timed(I)
          last_in_cpu_time_elapsed(I)=in_cpu_time_elapsed(I)
          last_in_wallclock_time_elapsed(I)=                            &
     &      in_wallclock_time_elapsed(I)
        ENDDO

        DO I=1,ni_number_of_timers
          last_ni_number_of_times_timed(I)=ni_number_of_times_timed(I)
          last_ni_cpu_time_elapsed(I)=ni_cpu_time_elapsed(I)
          last_ni_wallclock_time_elapsed(I)=                            &
     &      ni_wallclock_time_elapsed(I)
        ENDDO

        GOTO 9999  ! jump to end - no output on first call
      ENDIF

      WRITE(6,*)
      WRITE(6,*) '******************************************'
      WRITE(6,*)

      DO timer_kind=1,2  ! 1 is non-inclusive and 2 is inclusive
! Copy arguments into local arrays
        IF (action  ==  last_call_to_timer) THEN
          WRITE(6,*) 'END OF RUN - TIMER OUTPUT'
          WRITE(6,*) 'Timer information is for whole run'
          IF (timer_kind  ==  1) THEN  ! non-inclusive timer
            number_of_timers=ni_number_of_timers
            DO I=1,number_of_timers
              local_timer_name(I)=ni_timer_name(I)
              local_cpu_time_elapsed(I)=ni_cpu_time_elapsed(I)
              local_wallclock_time_elapsed(I)=                          &
     &          ni_wallclock_time_elapsed(I)
              local_number_of_times_timed(I)=                           &
     &          ni_number_of_times_timed(I)
            ENDDO
          ELSE ! timer_kind  ==  2 - inclusive timer
            number_of_timers=in_number_of_timers
            DO I=1,number_of_timers
              local_timer_name(I)=in_timer_name(I)
              local_cpu_time_elapsed(I)=in_cpu_time_elapsed(I)
              local_wallclock_time_elapsed(I)=                          &
     &          in_wallclock_time_elapsed(I)
              local_number_of_times_timed(I)=                           &
     &          in_number_of_times_timed(I)
            ENDDO
          ENDIF ! which timer kind this was
        ELSE  ! this is an intermediate output call
          WRITE(6,*) 'INTERMEDIATE TIMER OUTPUT :',message
          WRITE(6,*) 'Timer information is only for code executed ',    &
     &               'since last intermediate timer output.'
          IF (timer_kind  ==  1) THEN  ! non-inclusive timer
            number_of_timers=ni_number_of_timers
            DO I=1,number_of_timers
              local_timer_name(I)=ni_timer_name(I)
              local_cpu_time_elapsed(I)=ni_cpu_time_elapsed(I)-         &
     &                                  last_ni_cpu_time_elapsed(I)
              local_wallclock_time_elapsed(I)=                          &
     &          ni_wallclock_time_elapsed(I)-                           &
     &          last_ni_wallclock_time_elapsed(I)
              local_number_of_times_timed(I)=                           &
     &          ni_number_of_times_timed(I)-                            &
     &          last_ni_number_of_times_timed(I)
            ENDDO
          ELSE ! timer kind  ==  2 - inclusive timer
            number_of_timers=in_number_of_timers
            DO I=1,number_of_timers
              local_timer_name(I)=in_timer_name(I)
              local_cpu_time_elapsed(I)=in_cpu_time_elapsed(I)-         &
     &                                  last_in_cpu_time_elapsed(I)
              local_wallclock_time_elapsed(I)=                          &
     &          in_wallclock_time_elapsed(I)-                           &
     &          last_in_wallclock_time_elapsed(I)
              local_number_of_times_timed(I)=                           &
     &          in_number_of_times_timed(I)-                            &
     &          last_in_number_of_times_timed(I)
            ENDDO
          ENDIF  ! what timer type
        ENDIF  ! what action to perform

! Do work for non-inclusive timers

! Calculate the total time in the program (based on non-inclusive
! timers)
        IF (timer_kind  ==  1) THEN
          total_cpu_time = 0.0
          total_wallclock_time = 0.0
          DO I=1,number_of_timers
          total_cpu_time = total_cpu_time + local_cpu_time_elapsed(I)
          total_wallclock_time =                                        &
     &      total_wallclock_time + local_wallclock_time_elapsed(I)
          ENDDO

          WRITE(6,*) 'PE ',mype,' Elapsed CPU Time: ',                  &
     &               total_cpu_time
          WRITE(6,*) 'PE ',mype,'  Elapsed Wallclock Time: ',           &
     &                total_wallclock_time

! Calculate the total cpu time over all processors and the
! maximum elapsed time - so allowing a speedup to be caclulated

          total_cpu=total_cpu_time
          max_wall=total_wallclock_time

          CALL GC_RSUM(1,nproc,info,total_cpu)
          CALL GC_RMAX(1,nproc,info,max_wall)

          max_wall=MAX(max_wall,0.000001)
          WRITE(6,*)
          WRITE(6,*) 'Total Elapsed CPU Time: ',                        &
     &               total_cpu
          WRITE(6,*) 'Maximum Elapsed Wallclock Time: ',                &
     &               max_wall
          WRITE(6,*) 'Speedup: ',total_cpu/max_wall
          WRITE(6,*) '--------------------------------------------'

        ENDIF

! Sort subroutines into time order (based on wallclock time)

        DO I=1,number_of_timers-1
          DO J=(I+1),number_of_timers
            IF (local_wallclock_time_elapsed(J)  >                      &
     &          local_wallclock_time_elapsed(I)) THEN

!             Swap the two entries
              sortwork_real = local_cpu_time_elapsed(I)
              local_cpu_time_elapsed(I) = local_cpu_time_elapsed(J)
              local_cpu_time_elapsed(J) = sortwork_real

              sortwork_real = local_wallclock_time_elapsed(I)
              local_wallclock_time_elapsed(I) =                         &
     &          local_wallclock_time_elapsed(J)
              local_wallclock_time_elapsed(J) = sortwork_real

              sortwork_int = local_number_of_times_timed(I)
              local_number_of_times_timed(I) =                          &
     &          local_number_of_times_timed(J)
              local_number_of_times_timed(J) = sortwork_int

              sortwork_char = local_timer_name(I)
              local_timer_name(I) = local_timer_name(J)
              local_timer_name(J) = sortwork_char

            ENDIF
          ENDDO
        ENDDO

 20     FORMAT(20X, A45,I4)
 21     FORMAT(3X,'ROUTINE',14X,'CALLS',2X,'TOT CPU',4X,                &
     &       'AVERAGE',3X,'TOT WALL',2X,'AVERAGE',2X,                   &
     &       '% CPU',4X,'% WALL',4X,'SPEED-UP')
 22     FORMAT(3X,'ROUTINE',14X,'CALLS',2X,'TOT CPU',4X,                &
     &         'AVERAGE',3X,'TOT WALL',2X,'AVERAGE',2X,'SPEED-UP')
 23     FORMAT(I3,1X,A20,1X,I4,4(2X,F8.2),2(2X,F6.2),4X,F6.2)
        IF (timer_kind  ==  1) THEN
          WRITE(6,20) 'Non-Inclusive Timer Summary for PE ',mype
          WRITE(6,21)
        ELSE
          WRITE(6,20) 'Inclusive Timer Summary for PE ',mype
          WRITE(6,22)
        ENDIF

        DO I=1,number_of_timers
          IF (local_number_of_times_timed(I)  /=  0) THEN
            average_cpu_elapsed =  local_cpu_time_elapsed(I)/           &
     &                             local_number_of_times_timed(I)
            average_wallclock_elapsed = local_wallclock_time_elapsed(I)/&
     &                                  local_number_of_times_timed(I)
          ELSE
             average_cpu_elapsed = 0.0
             average_wallclock_elapsed = 0.0
          ENDIF

          IF (local_wallclock_time_elapsed(I)  >   0) THEN
            speed_up=local_cpu_time_elapsed(I)/                         &
     &               local_wallclock_time_elapsed(I)
          ELSE
            speed_up=1.0
          ENDIF

          IF (timer_kind  ==  1) THEN  ! non-inclusive timer has some
!                                      ! extra output

            percent_of_cpu_total = 100.0*local_cpu_time_elapsed(I)/     &
     &                             total_cpu_time
            percent_of_wallclock_total =                                &
     &        100.0*local_wallclock_time_elapsed(I)/                    &
     &        total_wallclock_time


            WRITE(6,23) I,local_timer_name(I),                          &
     &                  local_number_of_times_timed(I),                 &
     &                  local_cpu_time_elapsed(I),average_cpu_elapsed,  &
     &                  local_wallclock_time_elapsed(I),                &
     &                  average_wallclock_elapsed,                      &
     &                  percent_of_cpu_total,                           &
     &                  percent_of_wallclock_total,speed_up

          ELSE ! inclusive timer has slightly less to output

            WRITE(6,23) I,local_timer_name(I),                          &
     &                  local_number_of_times_timed(I),                 &
     &                  local_cpu_time_elapsed(I),average_cpu_elapsed,  &
     &                  local_wallclock_time_elapsed(I),                &
     &                  average_wallclock_elapsed,speed_up

          ENDIF

        ENDDO



! And now to assemble an overall timing assesment on PE0
! Each PE sends it total wallclock and cpu time spent in each routine
! to PE0, which calculates the average, s.d., max and min, and
! sorts on the basis of the average wallclock time
!
!
! We'll use the list of routines that PE0 already has as the master
! list.

        IF (mype  ==  0) THEN
          WRITE(6,*)
          WRITE(6,*) 'MPP Timing information : '
          WRITE(6,*)  nproc,' processors in configuration ',nproc_x,    &
     &                ' x ',nproc_y

          summ_n_timers=number_of_timers
          DO I=1,summ_n_timers
            summ_section(I)=local_timer_name(I)
          ENDDO
        ENDIF

! tell everyone else how many routines to do summary on - and which
! routines they are
        CALL GC_IBCAST(3213,1,0,nproc,info,summ_n_timers)
        CALL GC_CBCAST(3214,20*summ_n_timers,0,nproc,info,              &
     &                 summ_section)


        DO I=1,summ_n_timers
!
! Ensure that we do not change the output data area on the
! T3E before PE 0 has had a chance to get the data, if that is the
! protocol being used
!
          CALL GC_GSYNC (nproc,info)

! which section_ref is this for me?

          routine_id=0
          DO J=1,number_of_timers
            IF (local_timer_name(J)  ==  summ_section(I))               &
     &        routine_id=J
          ENDDO

          IF (routine_id  >   0) THEN
            wallclock_times(mype)=                                      &
     &        local_wallclock_time_elapsed(routine_id)
            cpu_times(mype)=local_cpu_time_elapsed(routine_id)
          ELSE
            wallclock_times(mype)=0.0
            cpu_times(mype)=0.0
          ENDIF

! send my information to PE 0.
          CALL GC_RSEND(1000+mype,1,0,info,wallclock_times(mype),       &
     &                wallclock_times(mype))
          CALL GC_GSYNC(nproc,info)

          IF (mype  ==  0) THEN
            DO J=0,nproc-1
              CALL GC_RRECV(1000+J,1,J,info,wallclock_times(J),         &
     &                      wallclock_times(J))
            ENDDO
          ENDIF
          CALL GC_GSYNC(nproc,info)

          CALL GC_RSEND(10000+mype,1,0,info,cpu_times(mype),            &
     &                cpu_times(mype))
          CALL GC_GSYNC(nproc,info)

          IF (mype  ==  0) THEN
            DO J=0,nproc-1
              CALL GC_RRECV(10000+J,1,J,info,cpu_times(J),              &
     &                      cpu_times(J))
            ENDDO
          ENDIF

          IF (mype  ==  0) THEN
! collect all the information - and start calculating the statistics
            wallclock_mean(I)=0.0
            cpu_total(I)=0.0
            wallclock_max(I)=-1.0E30
            wallclock_min(I)=1.0E30
            cpu_max(I)=-1.0E30
            cpu_min(I)=1.0E30

            DO J=0,nproc-1

              wallclock_mean(I)=wallclock_mean(I)+wallclock_times(J)
              cpu_total(I)=cpu_total(I)+cpu_times(J)

              IF (wallclock_times(J) >  wallclock_max(I)) THEN
                wallclock_max(I)=wallclock_times(J)
                wallclock_max_pe(I)=J
              ENDIF
              IF (wallclock_times(J) <  wallclock_min(I)) THEN
                wallclock_min(I)=wallclock_times(J)
                wallclock_min_pe(I)=J
              ENDIF
              IF (cpu_times(J) >  cpu_max(I)) THEN
                cpu_max(I)=cpu_times(J)
                cpu_max_pe(I)=J
              ENDIF
              IF (cpu_times(J) <  cpu_min(I)) THEN
                cpu_min(I)=cpu_times(J)
                cpu_min_pe(I)=J
              ENDIF

            ENDDO ! loop over processors

            IF (wallclock_max(I)  >   0.0) THEN
              speedup(I)=cpu_total(I)/wallclock_max(I)
            ELSE
              speedup(I)=1.0
            ENDIF
            efficiency(I)=speedup(I)/nproc

! and calculate the statistics
! first calculate the means
            wallclock_mean(I)=wallclock_mean(I)/nproc
            cpu_mean(I)=cpu_total(I)/nproc
! To stop a divide by zero later:
            IF (wallclock_mean(I)  ==  0.0) wallclock_mean(I)=1.0E-20
            IF (cpu_mean(I)  ==  0.0) cpu_mean(I)=1.0E-20
! and now the standard deviation
            wallclock_sd(I)=0.0
            cpu_sd(I)=0.0
            DO J=0,nproc-1
              wallclock_sd(I)=wallclock_sd(I)+                          &
     &          (wallclock_times(J)-wallclock_mean(I))*                 &
     &          (wallclock_times(J)-wallclock_mean(I))
              cpu_sd(I)=cpu_sd(I)+(cpu_times(J)-cpu_mean(I))*           &
     &                      (cpu_times(J)-cpu_mean(I))
            ENDDO
            wallclock_sd(I)=SQRT(wallclock_sd(I)/nproc)
            cpu_sd(I)=SQRT(cpu_sd(I)/nproc)

! Calculate the median
            DO J=0,nproc-2
              DO K=J+1,nproc-1
                IF (wallclock_times(K)  >   wallclock_times(J)) THEN
                  sortwork_real=wallclock_times(J)
                  wallclock_times(J)=wallclock_times(K)
                  wallclock_times(K)=sortwork_real
                ENDIF
                IF (cpu_times(K)  >   cpu_times(J)) THEN
                  sortwork_real=cpu_times(J)
                  cpu_times(J)=cpu_times(K)
                  cpu_times(K)=sortwork_real
                ENDIF
              ENDDO
            ENDDO

            IF (MOD(nproc,2)  ==  0) THEN
              wallclock_median(I)=(wallclock_times((nproc/2)-1)+        &
     &                             wallclock_times(nproc/2))*0.5
              cpu_median(I)=(cpu_times((nproc/2)-1)+                    &
     &                       cpu_times(nproc/2))*0.5
            ELSE
              wallclock_median(I)=wallclock_times(nproc/2)
              cpu_median(I)=cpu_times(nproc/2)
            ENDIF

          ENDIF ! am I PE 0?

        ENDDO ! loop over sections

! Sort and output the information on PE 0

        IF (mype  ==  0) THEN

          DO I=1,summ_n_timers-1
            DO J=(I+1),summ_n_timers
              IF (wallclock_max(J)  >   wallclock_max(I)) THEN

! Swap the entries I and J

              sortwork_char=summ_section(I)
              summ_section(I)=summ_section(J)
              summ_section(J)=sortwork_char

              sortwork_real=wallclock_mean(I)
              wallclock_mean(I)=wallclock_mean(J)
              wallclock_mean(J)=sortwork_real

              sortwork_real=wallclock_median(I)
              wallclock_median(I)=wallclock_median(J)
              wallclock_median(J)=sortwork_real

              sortwork_real=wallclock_sd(I)
              wallclock_sd(I)=wallclock_sd(J)
              wallclock_sd(J)=sortwork_real

              sortwork_real=wallclock_max(I)
              wallclock_max(I)=wallclock_max(J)
              wallclock_max(J)=sortwork_real

              sortwork_real=wallclock_min(I)
              wallclock_min(I)=wallclock_min(J)
              wallclock_min(J)=sortwork_real

              sortwork_int=wallclock_min_pe(I)
              wallclock_min_pe(I)=wallclock_min_pe(J)
              wallclock_min_pe(J)=sortwork_int

              sortwork_int=wallclock_max_pe(I)
              wallclock_max_pe(I)=wallclock_max_pe(J)
              wallclock_max_pe(J)=sortwork_int

              sortwork_real=cpu_mean(I)
              cpu_mean(I)=cpu_mean(J)
              cpu_mean(J)=sortwork_real

              sortwork_real=cpu_median(I)
              cpu_median(I)=cpu_median(J)
              cpu_median(J)=sortwork_real

              sortwork_real=cpu_sd(I)
              cpu_sd(I)=cpu_sd(J)
              cpu_sd(J)=sortwork_real

              sortwork_real=cpu_max(I)
              cpu_max(I)=cpu_max(J)
              cpu_max(J)=sortwork_real

              sortwork_real=cpu_min(I)
              cpu_min(I)=cpu_min(J)
              cpu_min(J)=sortwork_real

              sortwork_real=cpu_total(I)
              cpu_total(I)=cpu_total(J)
              cpu_total(J)=sortwork_real

              sortwork_real=speedup(I)
              speedup(I)=speedup(J)
              speedup(J)=sortwork_real

              sortwork_real=efficiency(I)
              efficiency(I)=efficiency(J)
              efficiency(J)=sortwork_real

              sortwork_int=cpu_min_pe(I)
              cpu_min_pe(I)=cpu_min_pe(J)
              cpu_min_pe(J)=sortwork_int

              sortwork_int=cpu_max_pe(I)
              cpu_max_pe(I)=cpu_max_pe(J)
              cpu_max_pe(J)=sortwork_int

              ENDIF
            ENDDO
          ENDDO

! and write out the information
          WRITE(6,*)
          IF (timer_kind  ==  1) THEN
            WRITE(6,*) 'MPP : None Inclusive timer summary'
          ELSE
            WRITE(6,*) 'MPP : Inclusive timer summary'
          ENDIF

          WRITE(6,*)
          WRITE(6,*)  'WALLCLOCK  TIMES'
          WRITE(6,40)
          DO I=1,summ_n_timers

            WRITE(6,41) I,summ_section(I),                              &
     &                 wallclock_mean(I),wallclock_median(I),           &
     &                 wallclock_sd(I),                                 &
     &                 (wallclock_sd(I)/wallclock_mean(I))*100.0,       &
     &                 wallclock_max(I),wallclock_max_pe(I),            &
     &                 wallclock_min(I),wallclock_min_pe(I)
          ENDDO

          WRITE(6,*)
          WRITE(6,*)  'CPU TIMES (sorted by wallclock times)'
          WRITE(6,40)
          DO I=1,summ_n_timers
            WRITE(6,41) I,summ_section(I),                              &
     &                 cpu_mean(I),cpu_median(I),                       &
     &                 cpu_sd(I),                                       &
     &                 (cpu_sd(I)/cpu_mean(I))*100.0,                   &
     &                 cpu_max(I),cpu_max_pe(I),                        &
     &                 cpu_min(I),cpu_min_pe(I)
          ENDDO

          WRITE(6,*)
          WRITE(6,*) 'PARALLEL SPEEDUP SUMMARY ',                       &
     &             '(sorted by wallclock times)'
          WRITE(6,50)
          DO I=1,summ_n_timers
            WRITE(6,51) I,summ_section(I),cpu_total(I),                 &
     &                  wallclock_max(I),speedup(I),                    &
     &                  efficiency(I)
          ENDDO


 40       FORMAT(4X,'ROUTINE',19X,'MEAN',3X,'MEDIAN',7X,                &
     &           'SD',3X,'% of mean',6X,'MAX',3X,                       &
     &           '(PE)',6X,'MIN',3X,'(PE)')

 41       FORMAT(I3,1X,A20,1X,                                          &
     &           3(1X,F8.2),5X,F6.2,'%',                                &
     &           2(1X,F8.2,1X,'(',I4,')'))

 50       FORMAT(4X,'ROUTINE',14X,'CPU TOTAL',3X,                       &
     &           'WALLCLOCK MAX',3X,'SPEEDUP',3X,                       &
     &           'PARALLEL EFFICIENCY')

 51       FORMAT(I3,1X,A20,1X,                                          &
     &           1X,F8.2,8X,F8.2,2X,F8.2,14X,F8.2)

        ENDIF
        WRITE(6,*)


      ENDDO ! loop over timer kind

! Finally copy the timer info into the last_* arrays so that the
! intermediate timer can calculate the timings since this point

      DO I=1,in_number_of_timers
        last_in_number_of_times_timed(I)=in_number_of_times_timed(I)
        last_in_cpu_time_elapsed(I)=in_cpu_time_elapsed(I)
        last_in_wallclock_time_elapsed(I)=                              &
     &    in_wallclock_time_elapsed(I)
      ENDDO

      DO I=1,ni_number_of_timers
        last_ni_number_of_times_timed(I)=ni_number_of_times_timed(I)
        last_ni_cpu_time_elapsed(I)=ni_cpu_time_elapsed(I)
        last_ni_wallclock_time_elapsed(I)=                              &
     &    ni_wallclock_time_elapsed(I)
      ENDDO


 9999 CONTINUE

      RETURN

      END SUBROUTINE TIMER_OUTPUT
