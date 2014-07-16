#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Write SCM diagnostic data to output files

      subroutine dump_streams(SCMop,day,time,                           &
     &           row_length,rows,model_levels,dayno_init,               &
     &           steptime,site,lat,lon,                                 &
     &           timeinit,year,lcal360)

      USE NetCDF
      implicit none

! Description:
      ! Write the diagnostic output data to file for all streams
      ! with a dumping period ending now.
! Method:
      ! At the time of writing documentation is available from
      ! http://www-nwp/~frlj
! Owner: Luke Jones
! History:
! Version  Date      Comment
! =======  ====      =======
! 5.5      06/02/03  Original code (Luke Jones)
! 6.0      05/09/03  Adapt to changes in contents of include files.
!                    Luke Jones.
! 6.0      02/09/03  Added new functionality to SCM output diagnostic
!                    system. System now controlled by new namelist,
!                    DIAGS. Luke Jones.
! 6.1      06/08/04  One line found to be non-portable, and one found
!                    to have a bug when used with row_length*rows>1.
!                    Both fixed. L. Jones.
! 6.2      19/01/06  Enabled outputting of SCM diagnostics on multiple
!                    sub-steps. Luke Jones.
!
! Code Description:
      ! Language: Fortran77 with some bits of Fortran90

! SCMop_type is defined in here...
#include "s_scmoptype_defn.h"

      ! NOTE: on the NEC this routine may be compiled in 32 bits
      ! because requests to create 32 bit integers (necessary to
      ! communicate with the NetCDF library routines - see incdf) are
      ! ignored by the compiler when compiling at 64 bits. Thus all
      ! input and output variables are explicitly declared as 64 bit
      ! using the KIND types i64, r64 and l64 delcared in
      ! s_scmoptype_defn.h

      type(SCMop_type) :: SCMop ! IN The derived-type structure
                                ! containing all the diagnostic
                                ! information
      integer(i64) model_levels ! IN The number of model levels
      integer(i64) day          ! IN Day
      real(r64) time            ! IN Time in seconds

      integer unit

! Parameters and stuff used around and about the internal
! workings of the diagnostic system
#include "s_scmop_internal.h"

      integer N,M,i,j,k,l       ! Counters
      integer(i64) N64          ! 64-bit copy of N
      integer dump              ! The dump number (1:ndumps)
      integer ncols,nrows,nlevs ! No. of columns,rows,levels of a
                                ! diagnostic
      integer nlines            ! Number of lines written to file
      integer dump_step         ! A dumping period

      ! Local variables for NetCDF-format output. Integers that are
      ! passed into NetCDF library routines should have the same
      ! precision as those used internally within the library (which
      ! may have been compiled at 32 bit), hence the use of incdf.
      integer(incdf) :: Status,shape(4)

      ! Local variables for FSSI-format output
      type fssi_dump_type
         integer :: ndiags
         integer, pointer :: diags(:)
         integer :: dumpsize
      end type fssi_dump_type
      type(fssi_dump_type) :: fssi_dump(maxnstreams)
      save fssi_dump
      integer, allocatable :: diags(:)
      integer ndiags

      integer(i64) row_length,rows
      integer(i64) dayno_init        ! Initial value of DAYNUMBER
      integer(i64) steptime          ! Elapsed time (=stepcount*timestep)
      integer(i64) site(row_length,rows) ! SSFM site WMO number
      integer(i64) timeinit          ! Start time of integration
      integer(i64) year              ! Year (4 digits!!!)

      real(r64) lat(row_length,rows) ! SSFM site latitude
      real(r64) lon(row_length,rows) ! SSFM site longitude

      logical(l64) lcal360           ! Flag for 360 day calendar

      INTEGER(i64) ElapsedDays       ! Whole days elapsed from start
                                     ! of run
      INTEGER(i64) ElapsedSecs       ! Whole seconds elapsed from start
                                     ! of run
      INTEGER(i64) BasisTimeDays     ! Whole days to basis time from
                                     ! start of calendar
      INTEGER(i64) BasisTimeSecs     ! Seconds in day at basis time
      INTEGER RecLen                 ! File record length
      INTEGER(i64) CurrentYear       ! Year at current time in model run
      INTEGER(i64) Month             ! Month at current time in
                                     ! model run
      INTEGER(i64) Day2              ! Day at current time in model run
      INTEGER(i64) DayNumber         ! Day number at current time in
                                     ! model run
      INTEGER(i64) Hour              ! Hour at current time in model run
      INTEGER(i64) Minute            ! Minute at current time in
                                     ! model run
      INTEGER(i64) Second            ! Second at current time in
                                     ! model run
      INTEGER isec                   ! Second at current time in
                                     ! model run
      INTEGER Forecastime            ! Forecast time in `hundred'
                                     ! hours
      INTEGER DumpSize               ! Size of dump containing
                                     ! diagnostics required for output
                                     ! to file
      INTEGER HdrOffset              ! Offset for heading at start of
                                     ! each output line, containing
                                     ! site, date, time, etc.
      PARAMETER (HdrOffset = 11*12)

      character*200                                                     &
     & fmt605                                                           &
     &,fmt620                                                           &
     &,fmt650

      character (len=100) fmt,fmt2 ! Hold format specifiers
      character (len=3) clsname ! Holds the short-name length as string

      ! Character function that incorporates the substep number into
      ! the short name of a diagnostic. Somewhat longer than a normal
      ! short name.
      character (len=lsname+10) :: add_substep_to_sname,short_name

      ! Reset no. of calls to SCMoutput (this should be done at the
      ! end of every timestep since it acts as the number of calls
      ! to SCMoutput within the current timestep)
      SCMop%nSCMoutput=0

      ! Loop over the streams
      do M=1,maxnstreams

         ! If this stream is not switched on then skip to the next
         if (SCMop%strm(M)%switch == 0) then
            CYCLE
         endif
         ! If there are no diagnostics in this stream then skip to
         ! the next
         if (SCMop%strm(M)%n_output <= 0) then
            CYCLE
         endif
         ! If we are not at the end of the dumping period of this stream
         ! (or, in the case of FSSI-format streams, it's not the first
         ! timestep) then skip to the next.
         dump_step=SCMop%strm(M)%dump_step
         if (mod(stepnmbr(SCMop)-1,dump_step)+1 /= dump_step            &
     &        .and..not.                                                &
     &        (SCMop%strm(M)%format == 2.and.stepnmbr(SCMop) == 1)) then
            CYCLE
         endif

         ! The output unit for this stream
         unit=SCMop%strm(M)%op_unit

         ! What we write to the unit depends on the stream format

!-------------------------------------------------------------------
         ! Formats 0 and 3: designed for reading by PV-wave routine
         ! scmread2.pro
!-------------------------------------------------------------------
         if (SCMop%strm(M)%format == 0.or.SCMop%strm(M)%format == 3)    &
     &        then

            ! Write the day and time
            write(unit,'(I9,1PE18.10E3)') day,time

            ! Write the diagnostic data
            do N=1,SCMop%nentries
               ! Is this entry to be sent to this stream, and
               ! not been flagged as one not to output?
               if (StreamIsOn(SCMop%streams(N),M).and..not.             &
     &              NotWritten(SCMop%streams(N))) then
                  ! Yes, write it.
                  write(unit,'(5(1PE18.10E3,","))')                     &
     &                 (SCMop%diag(N)%dump(i),                          &
     &                 i=1,SCMop%nelements(N))
               endif
            enddo

!-------------------------------------------------------------------
         ! Format=1 : designed for easy visual perusal.
!-------------------------------------------------------------------
         elseif (SCMop%strm(M)%format == 1) then

            ! Which dump is this?
            dump=stepnmbr(SCMop)/SCMop%strm(M)%dump_step

            ! Write the dump and step numbers
            write(unit,*) 'Dump= ',dump,' ; Timestep=',                 &
     &           stepnmbr(SCMop)

            ! Write the day and time
            write(unit,*) 'Day= ',day,'; Time= ',time

            ! We will keep tabs on the number of lines we've
            ! written in this dump
            nlines=0

            ! Write the diagnostic data
            do N=1,SCMop%nentries
               ! Is this entry to be written to this stream's
               ! output file?
               if (StreamIsOn(SCMop%streams(N),M).and..not.             &
     &              NotWritten(SCMop%streams(N))) then
                  ! Yes.

                  ! Get the number of columns, rows and levels
                  ! this diagnostic is defined on
                  ncols=SCMop%ncols(N)
                  nrows=SCMop%nrows(N)
                  nlevs=SCMop%nlevs(N)
                  ! Ensure that the sub-step is appended to the short
                  ! name of each diagnostic entry, if required.
! DEPENDS ON: add_substep_to_sname
                  N64=N
                  short_name=add_substep_to_sname(SCMop,N64)
                  ! Create a format to get all the numbers
                  ! corresponding to one point on one line
                  ! (i.e. all levels on one line)
                  write(fmt ,'(A,I3,A)')                                &
     &                 '(I3,1X,A,',nlevs,                               &
     &                 '(1PE18.10E3,","))'
                  write(clsname,'(I3)')lsname-6
                  write(fmt2,'(A,A,A,I3,A)')                            &
     &                 '(',clsname,'X,A,I4,A,',model_levels,            &
     &                 '(" Level=",I3,8X," "))'
                  do i=1,ncols
                     do j=1,nrows
                        if (mod(nlines,15) == 0) then
                            write(unit,fmt2)'Dump=',dump,':',           &
     &                          (k,k=1,model_levels)
                        endif
                        nlines=nlines+1
                        write(unit,fmt)SCmop%sname_id(N),               &
     &                       short_name,(SCMop%diag(N)%dump             &
     &                       ((k-1)*nrows*ncols+(j-1)*ncols+i),         &
     &                       k=1,nlevs)
                     enddo
                  enddo
               endif
            enddo

!-------------------------------------------------------------------
         ! Format=2 : format required for FSSI database
!-------------------------------------------------------------------
         elseif (SCMop%strm(M)%format == 2) then

            ! Calculate times and dates required for output
            ElapsedDays = dayno_init + SCMop%daycount-1                 &
     &           + int((TimeInit + StepTime) / 86400 )
            ElapsedSecs = mod(TimeInit + StepTime, int(86400,i64))

            ! Calculate Model Basis Time (i.e. begining of current year)
            ! relative to calender zero
! DEPENDS ON: time2sec
            call TIME2SEC(Year,1,0,0,0,0,0,0,                           &
     &           BasisTimeDays,BasisTimeSecs,LCal360)

            ! Calculate full date at current time in model integration
! DEPENDS ON: sec2time
            CALL SEC2TIME(ElapsedDays,ElapsedSecs,BasisTimeDays,        &
     &           BasisTimeSecs,CurrentYear,Month,Day2,Hour,             &
     &           Minute,Second,DayNumber,LCal360)

            If (SCMop%stepcount == 1 .and. SCMop%daycount == 1) Then
               ISec = 0
            EndIf

            ! Calculate forecast time in `hundred' hours
            Forecastime=((SCMop%daycount-1)*24+int(StepTime/3600))*100  &
     &           +mod(int(StepTime/60.0),60)

            ! If this is the first timestep then make a list of the
            ! diagnostics for output to this stream, diags(1:ndiags),
            ! together with the total size of one dump for one site,
            ! DumpSize. Save this info in the structure fssi_dump(M)
            ! so we don't have to re-calculate it next time.
            if (stepnmbr(SCMop) == 1) then
               ndiags=0
               DumpSize=0
               allocate(diags(SCMop%nentries))

               ! Unless heed_acceptlist=2 and the acceptlist has
               ! non-zero size then output all diagnostics in the same
               ! order as the calls to SCMoutput.
               if (SCMop%strm(M)%heed_hardwired /= 0.or.                &
     &              SCMop%strm(M)%heed_acceptlist == 0) then
                  do N=1,SCMop%nentries
                     ! Is this entry to be written to this stream's
                     ! output file?
                     if (StreamIsOn(SCMop%streams(N),M).and..not.       &
     &                    NotWritten(SCMop%streams(N))) then
                        ! Yes.
                        ndiags=ndiags+1
                        diags(ndiags)=N
                        ! This diagnostic will contribute a number of
                        ! elements to the dump for one site equal to
                        ! the number of levels it's defined on.
                        DumpSize=DumpSize+SCMop%nlevs(N)
                     endif
                  enddo
               else
                  ! Output only those diagnostics specified by the
                  ! accept_list and do it in the order given in the
                  ! accept_list.
                  do i=1,size(SCMop%strm(M)%accept_list)
                     do N=1,SCMop%nentries
                        ! Does the entry have the right short name?
                        if (trim(SCMop%sname(N)) ==                     &
     &                       trim(SCMop%strm(M)%accept_list(i))) then
                           ! Is this entry to be written to this
                           ! stream's output file?
                           if (StreamIsOn(SCMop%streams(N),M).and..not. &
     &                          NotWritten(SCMop%streams(N))) then
                              ! Yes.
                              ndiags=ndiags+1
                              diags(ndiags)=N
                              ! This diagnostic will contribute a
                              ! number of elements to the dump for
                              ! one site equal to the number of levels
                              ! it's defined on.
                              DumpSize=DumpSize+SCMop%nlevs(N)
                           endif
                        endif
                     enddo
                  enddo
               endif

               ! Store ndiags, diags(1:ndiags) and DumpSize in
               ! fssi_dump(M) so we don't have to execute this section
               ! of code again.
               fssi_dump(M)%ndiags=ndiags
               allocate(fssi_dump(M)%diags(ndiags))
               fssi_dump(M)%diags=diags(1:ndiags)
               fssi_dump(M)%dumpsize=DumpSize
               deallocate(diags)
            endif               ! if (stepnmbr(SCMop) == 1)

            ! Format statements for output
            fmt605=                                                     &
     &           '(2a12,2(i2.2,''/''),i4.4,1x,                          &
     &           i2.2,'':'',i2.2,''.'',i2.2,'' Diags'',i4,              &
     &           '' RL'',i8)'
            Write(fmt620,                                               &
     &           fmt='("(2x,''Namelist'',7x,''Day'',8x,''Month'',"      &
     &           "7x,''Year'',8x,''Hour'',"                             &
     &           "7x,''Minute'',6x,''Second'',9x,''FT'',8x,             &
     &           ''Site'',7x,''Lat'',8x,''Long'',4x,",                  &
     &           I4,"(5x,i3,'':'',i3))")') fssi_dump(M)%dumpsize
            Write(fmt650,                                               &
     &           fmt='("(a12,8(i12),2F12.5,",I4.4,                      &
     &           "(g12.5))")') fssi_dump(M)%dumpsize

            ! Write output to file, one line for each site at each
            ! time, all diagnostics (single and multiple level) on
            ! the same line
            Do i=1,row_length
               Do j=1,rows

                  ! Perform initialisation of the file if this is the
                  ! first (i,j) cycle on the first timestep
                  If (SCMop%stepcount == 1.and.SCMop%daycount == 1      &
     &                 .and.i == 1.and.j == 1) Then

                     RecLen = fssi_dump(M)%dumpsize*12 + HdrOffset
                     Open(Unit=unit,Recl=RecLen,Access='Sequential',    &
     &                    File=SCMop%strm(M)%filename)

                     ! Write out initial headings for outputs to file
                     ! (a dummy string is substituted for DiagRefFile at
                     ! this stage)
                     Write(unit,fmt=fmt605)'SSM         ',              &
     &                    'XXXXXXXXXXXX',Day2,Month,CurrentYear,        &
     &                    Hour,Minute,ISec,SCMop%strm(M)%n_output,RecLen
                     Write(unit,fmt=fmt620)                             &
     &                    ((SCMop%sname_id(fssi_dump(M)%diags(l)),k,    &
     &                    k=1,SCMop%nlevs(fssi_dump(M)%diags(l))),      &
     &                    l=1,fssi_dump(M)%ndiags)
                  EndIf

                  ! Write out diagnostics data to file
                  ! (a dummy string is substituted for DiagRefFile at
                  ! this stage).
                  Write(unit,fmt=fmt650)'XXXXXXXXXXXX',Day2,Month,      &
     &                 CurrentYear,Hour,Minute,Second,Forecastime,      &
     &                 Site(i,j),Lat(i,j),Lon(i,j),                     &
                  ! What follows is a little bit complicated: two
                  ! nested implied do-loops. The inner index is k,
                  ! representing vertical level, and the outer is l,
                  ! representing a given output diagnostic. For a given
                  ! k and l, what is written out is the element of the
                  ! dump array of lth output diagnostic corresponding
                  ! to point(i,j) and level k. Get it?
     &                 (                                                &
     &                   (                                              &
     &                     SCMop%diag(fssi_dump(M)%diags(l))%dump       &
     &                       (                                          &
     &                        (k-1)*SCMop%nrows(fssi_dump(M)%diags(l))* &
     &                        SCMop%ncols(fssi_dump(M)%diags(l)) +      &
     &                        (j-1)*SCMop%ncols(fssi_dump(M)%diags(l))+i&
     &                       ),                                         &
     &                     k=1,SCMop%nlevs(fssi_dump(M)%diags(l))       &
     &                   ),                                             &
     &                 l=1,fssi_dump(M)%ndiags                          &
     &                 )

               EndDo            ! j=1,rows
            EndDo               ! i=1,row_length

!-------------------------------------------------------------------
         ! Format=4 : NetCDF
!-------------------------------------------------------------------
         elseif (SCMop%strm(M)%format == 4) then

            ! Which dump is this?
            dump=stepnmbr(SCMop)/SCMop%strm(M)%dump_step

            ! Write out each diagnostic
            do N=1,SCMop%nentries
               ! Is this entry to be sent to this stream, and
               ! not been flagged as one not to output?
               if (StreamIsOn(SCMop%streams(N),M).and..not.             &
     &              NotWritten(SCMop%streams(N))) then
                  ! Yes. Get the proper shape of the variable (at the
                  ! moment it's being held as a 1D array).
                  shape=int((/SCMop%ncols(N),                           &
     &                        SCMop%nrows(N),                           &
     &                        SCMop%nlevs(N),                           &
     &                        int(1,i64)/),incdf)
                  Status = Nf90_Put_Var(                                &
     &                 NcID = int(SCMop%strm(M)%op_unit,incdf),         &
     &                 VarID = SCMop%netcdf_id(N),                      &
     &                 start = int((/1,1,1,dump/),incdf),               &
     &                 count = shape,                                   &
     &                 Values = reshape(SCMop%diag(N)%dump,shape))
               endif
            enddo

         else
            print*,'dump_streams ERROR: unknown format '//              &
     &           'for stream ',M

         endif                  ! if (SCMop%strm(M) == ...)

      enddo ! M=1,maxnstreams

      return
      END SUBROUTINE dump_streams
#endif
