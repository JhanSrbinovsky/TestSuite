#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Initialise the SCM output diagnostic system

      subroutine setup_diags(row_length,rows,model_levels,              &
     &     wet_model_levels,boundary_layer_levels,sm_levels,            &
     &     st_levels,land_points,ntiles,n_vis_thresh,                   &
     &     cloud_levels,                                                &
     &     total_nsteps,timestep,full_daysteps,ntrad,                   &
     &     a_sw_radstep_prog,a_sw_radstep_diag,ntrad1,                  &
     &     daycount,stepcount,num_substeps,main_diag_switch,            &
     &     netcdf_chunksize,nSCMDpkgs,L_SCMDiags,                       &
     &     SCMop)

      implicit none

! Description:
      ! Perform all pre-timestepping initialisation for the output
      ! diagnostic system, i.e. initialise SCMop. In particular, the
      ! DIAGS namelist is read in, and the output streams and domain
      ! profiles are set up. The output files are -not- initialised
      ! here since this requires information which is not available
      ! until the end of at least one timestep (see dump_streams_init).
      ! All arguments are INTENT IN, apart from main_diag_switch,
      ! the L_SCMDiags logicals and SCMop which are INTENT OUT.
! Method:
      ! At the time of writing documentation is available from
      ! http://www-nwp/~frlj
! Owner: Luke Jones
! History:
! Version  Date      Comment
! =======  ====      =======
! 5.5      06/02/03  Original code (Luke Jones)
! 6.0      05/09/03  Modify calls to SCMoutput. Luke Jones.
! 6.0      05/09/03  Adapt to changes in contents of include files.
!                    Luke Jones.
! 6.0      02/09/03  Added new functionality to diagnostic system,
!                    controlled by new namelist: DIAGS. Luke Jones.
! 6.2      29/11/05  Add cloud_levels domain profile.  Include
!                    L_SCMDiags logicals to DIAGS namelist for
!                    packages system.          A. Kerr-Munslow
! 6.2      19/01/06  Enabled outputting of SCM diagnostics on
!                    multiple sub-steps. Luke Jones.
! 6.2      08/02/08  Outputting variables a_sw_radstep_prog
!                    and a_sw_radstep_diag for new versions
!                    3C and 3Z of the radiation code.
!                                            (J.-C. Thelen)
!
! Code Description:
      ! Language: Fortran77 with some bits of Fortran90

! SCMop_type is defined in here...
#include "s_scmoptype_defn.h"

      type(SCMop_type) :: SCMop ! OUT The derived-type structure
                                ! containing all the diagnostic
                                ! information
      integer row_length,rows,model_levels,wet_model_levels,            &
     &     boundary_layer_levels,sm_levels,st_levels,land_points,       &
     &     ntiles,n_vis_thresh,cloud_levels,total_nsteps,full_daysteps, &
     &     ntrad, a_sw_radstep_prog, a_sw_radstep_diag,ntrad1
      real timestep
      integer :: num_substeps

      ! SCMop has pointers which will point to stepcount and daycount
      ! so that it can know which step it's on
      integer, target :: daycount,stepcount

! Parameters and stuff used around and about the internal
! workings of the diagnostic system
#include "s_scmop_internal.h"

      integer, parameter :: listlength=5000 ! The maximum length of the
                                ! accept and reject lists.
      integer i,j,            & ! Counters.
     &     pos,               & ! Current position in string.
     &     last_comma,        & ! Position of last comma.
     &     n_words              ! No. of comma separated words found.
      character (len=5) suffix  ! Holds frmt-dependent filename suffix
      character (len=listlength) list ! Temporarily holds either an
                                      ! accept or reject list
      character (len=lsname) name ! Temporarily holds a (short) name
                                  ! of a diagnostic

      Character*(*), Parameter ::  RoutineName = 'setup_diags'

      integer istatus             ! Error code

      ! Some constants which will be used to scale diagnostics
      real, parameter :: rsec_day=86400.
      real ntimestepsperday,oneKsecday

      integer nSCMDpkgs             ! No of SCM diagnostics packages
                                    ! (set in S_MAIN)
      logical L_SCMDiags(nSCMDpkgs) ! Logicals array for SCM
                                    ! diagnostics packages
!
!     Note that s_scmop.h is included in s_scmop_internal.h
!
!----------------------------------------------------------------------
      ! The DIAGS namelist variable declarations...
!----------------------------------------------------------------------
      integer main_diag_switch  ! If set to zero the diagnostic system
                         ! is  off: no diagnostics will be calculated
                         ! or output. Default=1.
      integer strm_switch(maxnstreams) ! If set to zero the respective
                         ! stream is off. Default=1 for stream 1,
                         ! 0 for all others.
      integer strm_unit(maxnstreams) ! The output unit to which the
                         ! stream will be written. Default=36+n.
      character (len=lfilename) strm_filename(maxnstreams) ! The name
                         ! of the file to which the stream
                         ! will be written. Default=compiler
                         ! dependent.
      integer strm_format(maxnstreams) ! Format of the output data:
                         ! 0=can be read into PV-Wave with scmread2.pro
                         ! 1=can be easily read with human eye
                         ! 2=SSFM format, compatible with FSSI d'base
                         ! Default=0.
      integer strm_dumpstep(maxnstreams) ! Dumping period of the
                         ! stream, i.e. number of timesteps between
                         ! dumps. If negative it's treated as a number
                         ! of seconds and is converted to the closest
                         ! possible non-zero number of timesteps. If
                         ! the number of dumps per day is `non-integer',
                         ! a warning will be printed. Default=1.
      integer strm_heed_hardwired(maxnstreams) ! If zero then
                         ! diagnostics sent to this stream by the
                         ! hard-wired stream list provided in the call
                         ! to SCMoutput are ignored. Default=1.
      integer strm_heed_acceptlist(maxnstreams) ! If zero then the
                         ! contents of strm_acceptlist are ignored.
      integer strm_heed_rejectlist(maxnstreams) ! If zero then the
                         ! contents of strm_rejectlist are ignored.
      character (len=listlength) strm_acceptlist(maxnstreams) ! Comma-
                         ! separated list of diagnostic (short) names
                         ! (as specifed in corresponding calls to
                         ! SCMoutput) which are to be sucked into the
                         ! respective stream. This provides a way of
                         ! sending a diagnostic to a stream even if
                         ! that stream is not specified as a
                         ! destination in the hard-wired list of the
                         ! respective SCMoutput call. Default=''
      character (len=listlength) strm_rejectlist(maxnstreams) ! Like
                         ! strm_acceptlist but this is a list of
                         ! diagnostics which are to be prevented from
                         ! being sent to the respective stream. If
                         ! there is any conflict between rejectlist
                         ! and either acceptlist and/or the
                         ! "hard-wired" preferences, rejectlist wins.
                         ! Default=''

      integer netcdf_chunksize ! The ChunkSize input to NF90_Create(),
                         ! the routine that creates a NetCDF file.
                         ! Controls a space versus time trade-off:
                         ! memory allocated in the netcdf library
                         !versus number of system calls.

                         ! Individual logicals for SCM diagnostics
                         ! packages:
      Logical L_SCMDiag_gen,     & ! General diagnostics  - default true
     &        L_SCMDiag_rad,     & ! Radiation            - default false
     &        L_SCMDiag_bl,      & ! Boundary layer       - default false
     &        L_SCMDiag_surf,    & ! Surface              - default false
     &        L_SCMDiag_land,    & ! Land points only     - default false
                                   ! reset in S_MAIN if not land point
     &        L_SCMDiag_sea,     & ! Sea points only      - default false
                                   ! reset in S_MAIN if not sea point
     &        L_SCMDiag_lsp,     & ! Large scale precip   - default false
     &        L_SCMDiag_conv,    & ! Convection           - default false
     &        L_SCMDiag_lscld,   & ! Large scale cloud    - default false
     &        L_SCMDiag_PC2,     & ! PC2                  - default false
                                   ! set in S_MAIN to false if PC2 off
     &        L_SCMDiag_forc,    & ! Forcing              - default false
     &        L_SCMDiag_incs       ! Increments           - default false

!----------------------------------------------------------------------
      ! The DIAGS namelist variable data statements...
!----------------------------------------------------------------------
      ! main_diag_switch, netcdf_chunksize and L_SCMDiags are
      ! initialised with a assignments below. They cannot be initialised
      ! with data statements because they are outputs from this routine.

      data strm_switch          /maxnstreams*0/ ! strm_switch(1) is
                                         ! set non-zero below.
      data strm_unit            /maxnstreams*0/ ! Default stream units
                                         ! are set in a do-loop below.
      data strm_filename        /maxnstreams*Default/
      data strm_format          /maxnstreams*3/
      data strm_dumpstep        /maxnstreams*1/
      data strm_heed_hardwired  /maxnstreams*1/
      data strm_heed_acceptlist /maxnstreams*1/
      data strm_heed_rejectlist /maxnstreams*1/
      data strm_acceptlist      /maxnstreams*''/
      data strm_rejectlist      /maxnstreams*''/

!----------------------------------------------------------------------
      ! The DIAGS namelist...
!----------------------------------------------------------------------
      Namelist/DIAGS/ main_diag_switch,                                 &
     &     strm_switch,                                                 &
     &     strm_unit,                                                   &
     &     strm_filename,                                               &
     &     strm_format,                                                 &
     &     strm_dumpstep,                                               &
     &     strm_heed_hardwired,                                         &
     &     strm_heed_acceptlist,                                        &
     &     strm_heed_rejectlist,                                        &
     &     strm_acceptlist,                                             &
     &     strm_rejectlist,                                             &
     &     netcdf_chunksize,                                            &
     &     L_SCMDiag_gen,                                               &
     &     L_SCMDiag_rad,                                               &
     &     L_SCMDiag_bl,                                                &
     &     L_SCMDiag_surf,                                              &
     &     L_SCMDiag_land,                                              &
     &     L_SCMDiag_sea,                                               &
     &     L_SCMDiag_lsp,                                               &
     &     L_SCMDiag_conv,                                              &
     &     L_SCMDiag_lscld,                                             &
     &     L_SCMDiag_PC2,                                               &
     &     L_SCMDiag_forc,                                              &
     &     L_SCMDiag_incs

      SCMop%first_pass=.true.  ! Diagnostic list not yet finalised
      SCMop%on=.false.         ! Timestepping not yet started
      SCMop%daycount=>daycount
      SCMop%stepcount=>stepcount
      SCMop%full_daysteps=full_daysteps
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      SCMop%a_sw_radstep_prog=a_sw_radstep_prog
      SCMop%a_sw_radstep_diag=a_sw_radstep_diag
#else
      SCMop%ntrad=ntrad
#endif
      SCMop%ntrad1=ntrad1
      SCMop%maxnentries=0   ! No memory allocated for diagnostics yet
      SCMop%nentries=0      ! No diagnostics requested to be stored yet
      SCMop%n_output=0      ! No diagnostics requested to be output yet
      SCMop%nSCMoutput=0    ! No calls to SCMoutput yet
      SCMop%num_substeps=num_substeps ! Expected number of sub-steps
      SCMop%substep_number=0          ! No substepping started yet

      ! Initialise one element of the domain profiles to flag
      ! undefined profiles
      do i=1,maxndomprof
         SCMop%d_lev1(i)=-1
      enddo

      ! Allocate memory for streams
      allocate(SCMop%strm(maxnstreams))

      ! Allocate some token space to the diag_mem array (ensures the
      ! SIZE(diag_mem) operation in SCMoutput doesn't fall over the
      ! first time)
      allocate(SCMop%diag_mem(1,2))


      ! Modify default values of some namelist variables...

      ! Diagnostic system on by default (cannot be set by data
      ! statement because it is intent out).
      main_diag_switch=1

      ! Default NetCDF file chunk size
#if defined(NEC)
      netcdf_chunksize=32000000
#else
      netcdf_chunksize=8192
#endif

      ! L_SCMDiags logicals (intent out) are also set to default values
      L_SCMDiag_gen     = .TRUE.
      L_SCMDiag_rad     = .FALSE.
      L_SCMDiag_bl      = .FALSE.
      L_SCMDiag_surf    = .FALSE.
      L_SCMDiag_land    = .FALSE.
      L_SCMDiag_sea     = .FALSE.
      L_SCMDiag_lsp     = .FALSE.
      L_SCMDiag_conv    = .FALSE.
      L_SCMDiag_lscld   = .FALSE.
      L_SCMDiag_PC2     = .FALSE.
      L_SCMDiag_forc    = .FALSE.
      L_SCMDiag_incs    = .FALSE.

      ! Stream 1 will be on by default
      strm_switch(1)=1
      ! The default output units for each stream (37 upwards)...
      do i=1,maxnstreams
         strm_unit(i)=36+i
      enddo

      ! Read the DIAGS namelist from namelist.scm
      open(10,file='namelist.scm',iostat=istatus)
      if (istatus == 0) then
         read(10,DIAGS,iostat=istatus)
         close(10)
         if (istatus /= 0) then
            write(6,*)'------------------------------------------------'
            write(6,*)'SETUP_DIAGS WARNING: DIAGS namelist not found in'
            write(6,*)'namelist file. Diagnostic system will operate on'
            write(6,*)'defaults.'
            write(6,*)'------------------------------------------------'
         endif
      else
         write(6,*) ' S_SETUP_DIAGS: ERROR OPENING FILE ON UNIT 10'
         write(6,*) ' FILENAME = namelist.scm'
         write(6,*) ' IOSTAT =',ISTATUS
         write(6,*) ' Continuing with default DIAGS namelist settings'
      endif

      ! No point doing any more if the diagnostic system has been
      ! switched off
      if (main_diag_switch == 0) RETURN

      ! NetCDF I/O is done using C routines which don't utilise Fortran
      ! unit numbers. Set to -1 so that this is reflected in the stream
      ! summary print-out.
      do i=1,maxnstreams
         if (strm_format(i) == 4) then
            strm_unit(i)=-1
         endif
      enddo

      ! Set up the output streams...
      do i=1,maxnstreams
         SCMop%strm(i)%switch=strm_switch(i)      ! Stream on/off
         SCMop%strm(i)%op_unit=strm_unit(i)       ! Output unit
         SCMop%strm(i)%format=strm_format(i)      ! Output format
         ! Output filename...
         if (strm_filename(i) /= Default) then
            SCMop%strm(i)%filename=strm_filename(i)
         else
            ! The filename suffix will be format dependent
            suffix='.out'      ! Unrecognised format
            if (SCMop%strm(i)%format == 0) then
               suffix='.dat'   ! Old Wave format
            elseif (SCMop%strm(i)%format == 1) then
               suffix='.txt'   ! Easy to read format
            elseif (SCMop%strm(i)%format == 2) then
               suffix='.fssi'  ! Format for FSSI
            elseif (SCMop%strm(i)%format == 3) then
               suffix='.dat'   ! New Wave format
            elseif (SCMop%strm(i)%format == 4) then
               suffix='.nc'    ! NetCDF format
            endif
            write(SCMop%strm(i)%filename,'(A,I2.2,A)')'stream',i,suffix
         endif

         ! If the dumping period is positive, take this as a number
         ! of timesteps. If it is negative take it as a number of
         ! seconds. If it is zero, replace it with the number of
         ! timesteps in the whole run (i.e. one dump on the
         ! last step).
         if (strm_dumpstep(i) == 0) strm_dumpstep(i)=total_nsteps
         if (strm_dumpstep(i) >  0) then
            ! The dumping period has been specified in timesteps
            SCMop%strm(i)%dump_step=strm_dumpstep(i)
         elseif (strm_dumpstep(i) <  0) then
            ! The dumping period has been specified in seconds, convert
            ! to nearest non-zero number of timesteps
            SCMop%strm(i)%dump_step=                                    &
     &           max(nint(-strm_dumpstep(i)/timestep),1)
         endif

         ! Accept diagnostics sent to stream by respective SCMoutput
         ! parameter? (0=no)
         SCMop%strm(i)%heed_hardwired=strm_heed_hardwired(i)

         ! Accept diagnostics sent to stream by namelist? (0=no)
         SCMop%strm(i)%heed_acceptlist=strm_heed_acceptlist(i)

         ! Reject diagnostics sent to stream by any method? (0=no)
         SCMop%strm(i)%heed_rejectlist=strm_heed_rejectlist(i)

         ! No. of diagnostics sent to this stream (none until at least
         ! first call to SCMoutput)
         SCMop%strm(i)%n_output=0

         ! We now want to extract the individual diagnostic names from
         ! the comma-separated lists in strm_acceptlist(i) and
         ! strm_rejectlist(i). We will scan through both strings
         ! twice, first to count the number of names so the correct
         ! space can be allocated (j=1 and j=3) and then again to put
         ! the names into the allocated arrays (j=2 and j=4).
         do j=1,4

            IF (J <= 2) THEN
               ! First two loops work on the acceptlist
               list=strm_acceptlist(i)
            ELSE
               ! Next two are indepedent of the first and work
               ! on the reject list
               list=strm_rejectlist(i)
            ENDIF

            ! We're going to scan through the list of
            ! comma-separated words one character at a time
            ! starting from the beginning.
            ! Make sure there's a trailing comma on the list
            list=trim(list)//','
            pos=1               ! Current position
            last_comma=0        ! Pos'n of last comma found
            n_words=0           ! No. of words found so far
            ! Loop through the characters of the string
            do while (pos <= len_trim(list))

#if defined(LINUX) || defined(NEC) || defined(IBM)
               if (list(pos:pos) == '\\') then
#else
               if (list(pos:pos) == '\') then
#endif
                  ! This character is a backslash - it escapes the
                  ! next character causing us to jump forward one.
                  pos=pos+1
               elseif (list(pos:pos) == ' '.and.                        &
     &                 last_comma == pos-1) then
                  ! This is a space trailing a comma - ignore by
                  ! "moving" the last comma forward one
                  last_comma=pos
               elseif (list(pos:pos) == ',') then
                  ! We've found a (non-escaped) comma. Extract the name
                  ! enclosed by this comma and the last.
                  name=list(last_comma+1:pos-1)
                  ! Update the position of the last comma found
                  last_comma=pos
                  ! Count the number of names we've found
                  n_words=n_words+1
                  ! On the second pass for each list, fill the allocated
                  ! arrays with the names.
                  IF (J == 2) SCMop%strm(i)%accept_list(n_words)=name
                  IF (J == 4) SCMop%strm(i)%reject_list(n_words)=name
               endif
               ! Update the position
               pos=pos+1

            enddo ! over characters in a list

            ! After the first pass for each list, allocate the arrays
            ! into which the names will be put on the second pass.
            IF (J == 1) allocate(SCMop%strm(i)%accept_list(n_words))
            IF (J == 3) allocate(SCMop%strm(i)%reject_list(n_words))

         enddo ! j=1,4

         ! FSSI-format output files require a dump for the first
         ! timestep as well as at the end of every dumping period. The
         ! way this has been coded in routine dump_streams means that
         ! if the dumping period is equal to one then this additional
         ! dump will not be produced and the file will be one line
         ! shorter than expected. I am assured that the SSFM will never
         ! be run with a dumping period of one, but it's worth printing
         ! a warning in case it ever is.
         if (SCMop%strm(i)%format == 2.and.                             &
     &        SCMop%strm(i)%dump_step == 1) then
            write(6,*)'---'
            write(6,*)'setup_diags WARNING: a FSSI-format output file'
            write(6,*)'has been requested with a dumping period of'
            write(6,*)'one. This file will not have an additional dump'
            write(6,*)'for the first timestep and so may be one line'
            write(6,*)'shorter than expected.'
            write(6,*)'---'
         endif

      enddo                     ! i=1,maxnstreams



      ! Set up the domain profiles
! DEPENDS ON: define_domprof
      call define_domprof(d_all,'all_levels', &
     &     1,row_length,                      & ! Every column
     &     1,rows,                            & ! Every row
     &     1,model_levels,SCMop)                ! Every level

! DEPENDS ON: define_domprof
      call define_domprof(d_sl,'single_level',&
     &     1,row_length,                      & ! Every column
     &     1,rows,                            & ! Every row
     &     1,1,SCMop)                           ! Only the lowest level

! DEPENDS ON: define_domprof
      call define_domprof(d_wet,'wet_levels',                           &
     &     1,row_length,                                                &
     &     1,rows,                                                      &
     &     1,wet_model_levels,SCMop)

! DEPENDS ON: define_domprof
      call define_domprof(d_bl,'bl_levels',                             &
     &     1,row_length,                                                &
     &     1,rows,                                                      &
     &     1,boundary_layer_levels,SCMop)

! DEPENDS ON: define_domprof
      call define_domprof(d_soilm,'soil_moist_levels',                  &
     &     1,row_length,                                                &
     &     1,rows,                                                      &
     &     1,sm_levels,SCMop)

! DEPENDS ON: define_domprof
      call define_domprof(d_soilt,'soil_temp_levels',                   &
     &     1,row_length,                                                &
     &     1,rows,                                                      &
     &     1,st_levels,SCMop)

! DEPENDS ON: define_domprof
      call define_domprof(d_land,'land_points',                         &
     &     1,land_points,       & ! This is a fudge !
     &     1,1,                                                         &
     &     1,1,SCMop)

! DEPENDS ON: define_domprof
      call define_domprof(d_allxtra,'all_levs_plus1',                   &
     &     1,row_length,                                                &
     &     1,rows,                                                      &
     &     1,model_levels+1,SCMop)

! DEPENDS ON: define_domprof
      call define_domprof(d_tile,'tile_types',                          &
     &     1,row_length,                                                &
     &     1,rows,                                                      &
     &     1,ntiles,SCMop)

! DEPENDS ON: define_domprof
      call define_domprof(d_point,'single_point',1,1,1,1,1,1,SCMop)

! DEPENDS ON: define_domprof
      call define_domprof(d_vis,'vis_thresholds',                       &
     &     1,row_length,                                                &
     &     1,rows,                                                      &
     &     1,n_vis_thresh,SCMop)

! DEPENDS ON: define_domprof
      call define_domprof(d_cloud,'cloud_levels',                       &
     &     1,row_length,                                                &
     &     1,rows,                                                      &
     &     1,cloud_levels,SCMop)

      ! Some diagnostics are just constants which will be combined
      ! with other, variable diagnostics. May as well just call
      ! SCMoutput for them once then, rather than every timestep.

! DEPENDS ON: scmoutput
      call SCMoutput(rsec_day,                                          &
           'sec_day','No. of seconds per day','-',                      &
           t_const,d_point,DoNotWrite(default_streams),'',RoutineName)

      oneKsecday=1000.*rsec_day
! DEPENDS ON: scmoutput
      call SCMoutput(oneKsecday,                                        &
           'oneKsecday','1000.*rsec_day','-',                           &
            t_const,d_point,DoNotWrite(default_streams),'',RoutineName)

      ntimestepsperday=rsec_day/timestep
! DEPENDS ON: scmoutput
      call SCMoutput(ntimestepsperday,                                  &
           'ntspday','No. of timesteps per day','-',                    &
           t_const,d_point,DoNotWrite(default_streams),'',RoutineName)

      ! Place details from read in diagnostics packages logicals
      ! into logicals array
      ! Note: general is default true, the rest default false
      L_SCMDiags(SCMDiag_gen)  = L_SCMDiag_gen
      L_SCMDiags(SCMDiag_rad)  = L_SCMDiag_rad
      L_SCMDiags(SCMDiag_bl)   = L_SCMDiag_bl
      L_SCMDiags(SCMDiag_surf) = L_SCMDiag_surf
      L_SCMDiags(SCMDiag_land) = L_SCMDiag_land
      L_SCMDiags(SCMDiag_sea)  = L_SCMDiag_sea
      L_SCMDiags(SCMDiag_lsp)  = L_SCMDiag_lsp
      L_SCMDiags(SCMDiag_conv) = L_SCMDiag_conv
      L_SCMDiags(SCMDiag_lscld)= L_SCMDiag_lscld
      L_SCMDiags(SCMDiag_PC2)  = L_SCMDiag_PC2
      L_SCMDiags(SCMDiag_forc) = L_SCMDiag_forc
      L_SCMDiags(SCMDiag_incs) = L_SCMDiag_incs

      return

      END SUBROUTINE setup_diags
#endif
