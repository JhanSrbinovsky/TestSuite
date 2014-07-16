#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Initialise SCM diagnostic output files

      subroutine dump_streams_init(SCMop,                               &
     &     row_length,rows,model_levels,wet_levels,                     &
     &     bl_levels,cloud_levels,ozone_levels,st_levels,               &
     &     sm_levels,ntiles,                                            &
     &     year_init,month_init,day_init,hour_init,min_init,sec_init,   &
     &     timestep,ndayin,nminin,nsecin,sec_day,tot_nsteps,            &
     &     ntrad,a_sw_radstep_prog,a_sw_radstep_diag,                   &
     &     z_top_of_model,first_constant_r_rho_level,                   &
     &     eta_theta,eta_rho,orog,r_theta_levels,r_rho_levels,          &
     &     netcdf_chunksize)

      USE NetCDF
      implicit none

! Description:
      ! Perform all steps necessary for subsequent writing of
      ! output stream data : write the headers of the data files
      ! and any auxilliary files that may accompany them
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
! 6.1      06/08/04  Output model start time to format-3 dumps.
!                    L. Jones.
! 6.1      06/08/04  Altered three format statements to allow for
!                    row_length*rows > 999. L. Jones.
! 6.2      19/01/06  Enabled outputting of SCM diagnostics on multiple
!                    sub-steps. Luke Jones.
! 6.2      08/02/06  Functionality for improved time stepping, radiative
!                    forcing and radiance code  in version 3C and 3Z
!                    added. Only timestepping and radiances are 
!                    activated in the SCM. (J.-C. Thelen)
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

      type(SCMop_type) :: SCMop ! INOUT The derived-type structure
                                ! containing all the diagnostic
                                ! information
      ! IN Horizontal and vertical model domain...
      integer(i64)                                                      &
     &     row_length,rows,model_levels,wet_levels,bl_levels,           &
     &     cloud_levels,ozone_levels,st_levels,sm_levels,ntiles
      integer(i64)         &
     &   year_init         & ! IN Initial year
     &  ,month_init        & ! IN Initial month
     &  ,day_init          & ! IN Initial day
     &  ,hour_init         & ! IN Initial hour
     &  ,min_init          & ! IN Initial minute
     &  ,sec_init            ! IN Initial second
      real(r64) timestep     ! IN Timestep
      integer(i64)         &
     &   ndayin            & ! IN No. of days requested in run
     &  ,nminin            & ! IN No. of minutes requested in run
     &  ,nsecin            & ! IN No. of seconds requested in run
     &  ,sec_day           & ! IN No. of seconds in a day (Why this
                             ! isn't in an include file I don't know)
     &  ,tot_nsteps        & ! IN Total no. of steps to be done
     &  ,ntrad             & ! IN No. of steps between calls to 
                             ! radiation (3A)
     &  ,a_sw_radstep_prog & ! IN No. of timesteps between prognostic 
                             ! calls to radiation (3C/3Z)
     &  ,a_sw_radstep_diag   ! IN No. of timesteps between diagnostic 
                             ! calls to radiation (3C/3Z)
      real(r64)            &
     &   z_top_of_model      ! IN The height of the "top of the
                             ! atmosphere"
      integer(i64)         & ! IN Lowest rho level that has constant
     &   first_constant_r_rho_level                         ! height
      real(r64)                     &
     &   eta_theta(model_levels+1)  & ! IN The etas of the theta and
     &  ,eta_rho(model_levels)      & ! rho levels
     &  ,orog(row_length,rows)      & ! IN Orography height
     ! IN The physical heights of the theta and rho levels...
     &  ,r_theta_levels(row_length,rows,0:model_levels) &
     &  ,r_rho_levels(row_length,rows,model_levels)
      ! IN The ChunkSize input to NF90_Create(), the routine that
      ! creates a NetCDF file. Controls a space versus time trade-off:
      ! memory allocated in the netcdf library versus number of system
      ! calls.
      integer(i64) netcdf_chunksize

! Parameters and stuff used around and about the internal
! workings of the diagnostic system
#include "s_scmop_internal.h"

      integer GET_UM_VERSION    ! Function subroutine which returns
                                ! the UM version as defined by
                                ! environment variable $VN.
      integer um_version_env    ! Stores UM version from $VN.
      character (len=3) um_version_c ! The UM version as a string.
      character (len=100) fmt   ! Holds format specifiers
      integer i,N,M,          & ! Counters
     &     nhourin              ! No. of hours requested in run
      integer(i64) N64          ! 64-bit counter for passing to other 
                                ! routines compiled in 64 bit
      integer unit
      integer(i64) diags(SCMop%nentries),ndiags
      real ndump

      ! Copy of netcdf_chunksize at the correct precision for passing
      ! into NetCCDF routines. May be altered by NF90_Create.
      integer(incdf) netcdf_chunksize_io

      ! NetCDF identifiers for each dimension of each domain
      ! profile. Integers that are passed into NetCDF library routines
      ! should have the same precision as those used internally within the
      ! library (which may have been compiled at 32 bit), hence the use of
      ! incdf.
      integer (incdf) :: netcdf_dimension_time,                         &
     &                   netcdf_dimension(3,maxndomprof)

      ! Other NetCDF-format related variables
      integer (incdf) :: Status,NcID
      character (len=5) :: c_unit

      ! Character function that incorporates the substep number into
      ! the short name of a diagnostic. Somewhat longer than a normal
      ! short name.
      Character (len=lsname+10) :: ADD_SUBSTEP_TO_SNAME

      ! Temporary variables to hold output from ADD_SUBSTEP_TO_SNAME
      Character (len=lsname+10) :: short_name
      Character (len=lsname+10), allocatable :: short_names(:)

      ! Get the UM version from the environment variable $VN
! DEPENDS ON: get_um_version
      um_version_env=get_um_version()

      ! Check for a sensible version number
      If (um_version_env < 0) Then
         write(*,'(A)')'---'
         write(6,'(A)')'dump_streams_init: WARNING, the '           //  &
     &                 'environment variable $VN does not seem'
         write(6,'(A)')'to be set to a valid UM version number. I ' //  &
     &                 'will assume this is UM VN6.2,'
         write(6,'(A)')'but if incorrect it could lead to '         //  &
     &                 'unreadable output files.'
         write(*,'(A)')'---'
         um_version_env=602
      endif

      ! Convert the UM version number to character form
      write(um_version_c,'(I1,A1,I1)')int(um_version_env/100.),'.',     &
     &     mod(int(um_version_env),10)

      ! Write out the UM version number
      write(*,'(A)')'=============='
      write(*,'(A)')'UM Version='//trim(um_version_c)
      write(*,'(A)')'=============='
      write(*,'(A)')' '

      ! Give a summary of the output streams
      if (SCMop%n_output /= 0) then
         write(*,'(65("-"))')
         write(*,'(A)')'Summary of open streams:'
         write(*,'(65("-"))')
         write(*,'(A)')'No.| Unit | Time between dumps | No. dgnstcs '//&
     &        '| Format | Filename '
         write(*,'(A)')'   |      | (steps) |(minutes) |             '//&
     &        '|        |          '
         write(*,'(A)')'---+------+---------+----------+-------------'//&
     &        '+--------+----------'
         do i=1,maxnstreams
            if (SCMop%strm(i)%switch /= 0) then
               ! Write the stream output unit into a character variable,
               ! unless this is a NetCDF stream, in which case the unit
               ! is not used.
               if (SCMop%strm(i)%format /= 4) then
                  write(c_unit,'(I5)') SCMop%strm(i)%op_unit
               else
                  c_unit='n/a'
                  c_unit=adjustr(c_unit)
               endif
               write(*,'(I2," |",A5," |",I8," | ",G8.3," |",I12,        &
     &              " |",I7," | ",A)')                                  &
     &              i,c_unit,SCMop%strm(i)%dump_step,                   &
     &              SCMop%strm(i)%dump_step*timestep/60.,               &
     &              SCMop%strm(i)%n_output,SCMop%strm(i)%format,        &
     &              trim(SCMop%strm(i)%filename)
            endif
         enddo
         write(*,'(65("-"))')
         write(*,'(A)')' '

      else
         print*,'dump_streams_init: No diagnostics to be written.'
         ! Nothing to do if there's no diagnostics to be written
         goto 999
      endif

      ! write the scumlist file
! DEPENDS ON: write_scumlist
      N64=0
      call write_scumlist(SCMop,N64)
      ! write the no. of levels, etc into a file
! DEPENDS ON: write_domain_data
      call write_domain_data(row_length,rows,model_levels,              &
     &     z_top_of_model,first_constant_r_rho_level,                   &
     &     orog,eta_theta,eta_rho)

      ! Calculate nhourin from nminin and nsecin
      nhourin = int(nminin/60.0+nsecin/3600.0)

      ! Write the headers of the output data files
      do N=1,maxnstreams
         if (SCMop%strm(N)%switch /= 0.and.                             &
     &        SCMop%strm(N)%n_output >  0) then
            ! This stream is open and there are diagnostics
            ! to be written out from it
            unit=SCMop%strm(N)%op_unit
            ! Check this unit is free (not applicable for NetCDF
            ! files)
            do M=1,N-1
               if (SCMop%strm(M)%switch /= 0.and.                       &
     &              SCMop%strm(M)%format /= 4.and.                      &
     &              SCMop%strm(M)%n_output >  0.and.                    &
     &              SCMop%strm(M)%op_unit == unit) then
                  print*,'dump_streams_init ERROR: two output files '// &
     &                 'with the same unit numbers: ',N,M,unit
                  print*,'Output redirected to unit 103!'
                  unit=103
                  SCMop%strm(N)%op_unit=unit
               endif
            enddo

            ! The header depends on the format

!-------------------------------------------------------------------
            ! Format=0: designed for reading by PV-wave routine
            ! scmread2.pro
!-------------------------------------------------------------------
            if (SCMop%strm(N)%format == 0) then

               ! Open the stream's output file
               open(unit=unit,file=SCMop%strm(N)%filename)

               ! Write the UM version and the stream format
               write(unit,'(A,I1)')'vn'//um_version_c//',old_format'
               ! Write the total no. of days, hours, dumps per day and
               ! steps, plus the timestep, the dumping period of this
               ! stream and the no. of diagnostics per dump
               ndump=sec_day/(SCMop%strm(N)%dump_step*timestep)
               write(unit,                                              &
     &              '(I4,1X,I4,1X,F7.3,1X,I7,1X,F9.2,1X,I9,1X,I4)')     &
     &              ndayin,nhourin,ndump,tot_nsteps,timestep,           &
     &              SCMop%strm(N)%dump_step,SCMop%strm(N)%n_output
               ! Make a list the diagnostics we're going to write
               ndiags=0
               do i=1,SCMop%nentries
                  if (StreamIsOn(SCMop%streams(i),N).and..not.          &
     &                 NotWritten(SCMop%streams(i))) then
                     ndiags=ndiags+1
                     diags(ndiags)=i
                  endif
               enddo
               if (ndiags /= SCMop%strm(N)%n_output) then
                  print*,'dump_streams_init ERROR: an inconsistency '// &
     &                 'has ocurred !',ndiags,SCMop%strm(N)%n_output,   &
     &                 N,SCMop%nentries
                  ! Switch the diagnostic system off so dodgy data is
                  ! not mistakenly used in good faith.
                  print*,'Switching diagnostic system off!'
                  SCMop%on=.false.
               endif
               ! Write that list to the file
               write(fmt,'(A,I5,A)') '(',ndiags,'(I5,1X))'
               write(unit,fmt)(SCMop%sname_id(diags(i)),i=1,ndiags)

!-------------------------------------------------------------------
            ! Format=1 : designed for easy visual perusal.
!-------------------------------------------------------------------
            elseif (SCMop%strm(N)%format == 1) then

               ! Open the stream's output file
               open(unit=unit,file=SCMop%strm(N)%filename)

               write(unit,'(A,A)')'UM version=',um_version_c
               ! Write the total no. of days, hours, dumps per day and
               ! steps, plus the timestep, the dumping period of this
               ! stream and the no. of diagnostics per dump
               ndump=sec_day/(SCMop%strm(N)%dump_step*timestep)
               write(unit,'(A)')'General run information...'
               write(unit,'(A,I4,1X,A,I4,1X,A,F7.3,1X,A,I7)')           &
     &              'No. of days=',ndayin,'No. of hrs=',nhourin,        &
     &              'No. of dumps=',ndump,'Total number of steps=',     &
     &              tot_nsteps
               write(unit,'(A,F9.2,1X,A,I9,1X,A,I4)')                   &
     &              'Timestep=',timestep,                               &
     &              'No. of steps between dumps=',                      &
     &              SCMop%strm(N)%dump_step,                            &
     &              'No. of diags outputting to this stream=',          &
     &              SCMop%strm(N)%n_output
               ! list the diagnostics we're going to write
               ndiags=0
               do i=1,SCMop%nentries
                  if (StreamIsOn(SCMop%streams(i),N).and..not.          &
     &                 NotWritten(SCMop%streams(i))) then
                     ndiags=ndiags+1
                     diags(ndiags)=i
                  endif
               enddo
               if (ndiags /= SCMop%strm(N)%n_output) then
                  print*,'dump_streams_init ERROR: an inconsistency '// &
     &                 'has ocurred !',ndiags,SCMop%strm(N)%n_output,   &
     &                 N,SCMop%nentries
                  ! Switch the diagnostic system off so dodgy data is
                  ! not mistakenly used in good faith.
                  print*,'Switching diagnostic system off!'
                  SCMop%on=.false.
               endif
               write(unit,'(A)')'Diagnostics in this file... '//        &
     &              '(subset of entries in scumlist file)'

               ! Ensure that the sub-step is appended to the short
               ! name of each diagnostic entry, if required.

               Allocate(short_names(ndiags))

               Do i=1,ndiags
! DEPENDS ON: add_substep_to_sname
                 short_names(i)=ADD_SUBSTEP_TO_SNAME(SCMop,diags(i))
               End Do ! i

               ! Write the information
               write(unit,'(I3,1X,A,1X,A)')                             &
     &              (SCMop%sname_id(diags(i)),short_names(i),           &
     &              SCMop%lname(diags(i)),i=1,ndiags)

               Deallocate(short_names)

!-------------------------------------------------------------------
            ! Format=2 : format required for FSSI database
!-------------------------------------------------------------------
            elseif (SCMop%strm(N)%format == 2) then

               ! Initialisation done in DUMP_STREAMS for the
               ! time being.

!-------------------------------------------------------------------
            ! Format=3 : new format designed for reading by PV-wave
            ! routine scmread2.pro
!-------------------------------------------------------------------
            elseif (SCMop%strm(N)%format == 3) then

               ! Open the stream's output file
               open(unit=unit,file=SCMop%strm(N)%filename)

               ! Write the UM version, plus any extra strings necessary
               ! for PV-Wave to know what to expect in the file.
               write(unit,'(A,I1)')'vn'//um_version_c//                 &
     &              ',start_time_present'

               ! Write information about the size and shape of
               ! the domain
               write(unit,'(I4,I4,8I3)') row_length,rows,               &
     &              model_levels,wet_levels,bl_levels,cloud_levels,     &
     &              ozone_levels,st_levels,sm_levels,ntiles

               ! Write the starting time of the run as specified in
               ! the INDATA namelist
               write(unit,'(I6,1x,I4,1x,I5,1x,I4,1x,I4,1x,I4,1x)')      &
     &              year_init,month_init,day_init,hour_init,min_init,   &
     &              sec_init

               ! Write information pertaining to the length of
               ! the model run and time in general
!#if defined(A01_3C) || defined(A02_3C) \
! || defined(A01_3Z) || defined(A02_3Z)
! When 3C/3Z radiation is implemented in the SCM, it will require a 
! change to the formatting of this output file as below. But since the
! code doesn't yet exist in SCMoutput to handle 3C/3Z
! radiation-timestep arrays, we'll leave this commented out for now.
!               Write(unit,'(I9,1x,F9.2,1x,I9,1x,I4,1x,I4                &
!     &                                             ,1xI5,1x,I5,1x,I5)') &
!     &              tot_nsteps,timestep,SCMop%strm(N)%dump_step,        &
!     &              a_sw_radstep_prog,a_sw_radstep_diag,                &
!     &              ndayin,nminin,nsecin
!#else
               Write(unit,'(I9,1x,F9.2,1x,I9,1x,I4,1x,I5,1x,I5,1x,I5)') &
     &              tot_nsteps,timestep,SCMop%strm(N)%dump_step,ntrad,  &
     &              ndayin,nminin,nsecin
!#endif

               ! Write information about the levels and orography...
               write(fmt,'(A,I3,A)')                                    &
     &            '(',model_levels+1 ,'(1PE18.10E3," "))'
               write(unit,fmt) eta_theta

               write(fmt,'(A,I3,A)')                                    &
     &            '(',model_levels   ,'(1PE18.10E3," "))'
               write(unit,fmt) eta_rho

               write(fmt,'(A,I6,A)')                                    &
     &            '(',row_length*rows,'(1PE18.10E3," "))'
               write(unit,fmt) orog

               write(fmt,'(A,I6,A)')'(',row_length*rows*                &
     &              (model_levels+1),'(1PE18.10E3," "))'
               write(unit,fmt) r_theta_levels

               write(fmt,'(A,I6,A)')'(',row_length*rows*                &
     &              (model_levels  ),'(1PE18.10E3," "))'
               write(unit,fmt) r_rho_levels

               ! Write data into the file pertaining to the
               ! individual diagnostics that it will contain.
! DEPENDS ON: write_scumlist
               N64=N
               call WRITE_SCUMLIST(SCMop,N64)

!-------------------------------------------------------------------
            ! Format=4 : NetCDF
!-------------------------------------------------------------------
            elseif (SCMop%strm(N)%format == 4) then

               ! Create the NetCDF file. Chunksize is a parameter,
               ! the size of which may greatly affect the I/O speed 
               ! of the NetCDF file writing.
               netcdf_chunksize_io=netcdf_chunksize
               Status = Nf90_Create(SCMop%strm(N)%filename,             &
     &              Nf90_Clobber,NcID,Chunksize=netcdf_chunksize_io)

               ! Declare the three spatial dimensions of each domain
               ! profile
               do i=1,maxndomprof
                 ! Check this domain has been defined.
                 if (SCMop%d_lev1(i) /= -1) then
                   Status = Nf90_Def_Dim(                               &
     &              NcID = NcID,                                        &
     &              Name = trim(adjustl(SCMop%d_name(i)))//'_i',        &
     &              Len =int(SCMop%d_rowa2(i)-SCMop%d_rowa1(i)+1,incdf),&
     &              DimID = netcdf_dimension(1,i))                      
                   Status = Nf90_Def_Dim(                               &
     &              NcID = NcID,                                        &
     &              Name = trim(adjustl(SCMop%d_name(i)))//'_j',        &
     &              Len =int(SCMop%d_rowb2(i)-SCMop%d_rowb1(i)+1,incdf),&
     &              DimID = netcdf_dimension(2,i))
                   Status = Nf90_Def_Dim(                               &
     &              NcID = NcID,                                        &
     &              Name = trim(adjustl(SCMop%d_name(i)))//'_k',        &
     &              Len =int(SCMop%d_lev2(i)-SCMop%d_lev1(i)+1,incdf),  &
     &              DimID = netcdf_dimension(3,i))
                 endif
               enddo

               ! Declare the time dimension
               Status = Nf90_Def_Dim(NcID = NcID, Name = "time",        &
     &              Len = int(tot_nsteps/SCMop%strm(N)%dump_step,incdf),&
     &              DimID = netcdf_dimension_time)

               ! Define each diagnostic that will go to this stream
               do i=1,SCMop%nentries
                  ! Is this entry to be sent to this stream, and
                  ! not been flagged as one not to output?
                  if (StreamIsOn(SCMop%streams(i),N).and..not.          &
     &                 NotWritten(SCMop%streams(i))) then
                     ! Get the short name with the substep number
                     ! appended if necessary.
                     N64=i
                     short_name=ADD_SUBSTEP_TO_SNAME(SCMop,N64)
                     ! Define this variable
                     Status = Nf90_Def_Var(NcID = NcID,                 &
     &                    Name = short_name,                            &
     &                    XType = Nf90_Float,                           &
     &                    DimIDs = (/                                   &
     &                      netcdf_dimension(1,SCMop%domprof(i)),       &
     &                      netcdf_dimension(2,SCMop%domprof(i)),       &
     &                      netcdf_dimension(3,SCMop%domprof(i)),       &
     &                      netcdf_dimension_time/),                    &
     &                    VarID = SCMop%netcdf_id(i))
                     if (Status /= nf90_NoErr) then
                        print*,'ERROR DEFINING NETCDF VARIABLE: ',      &
     &                       short_name
                     endif
                     ! Define its attributes
                     Status = Nf90_Put_Att(NcID = NcID,                 &
     &                    VarID = SCMop%netcdf_id(i), Name = "units",   &
     &                    Values = SCMop%units(i))
                     Status = Nf90_Put_Att(NcID = NcID,                 &
     &                    VarID = SCMop%netcdf_id(i),                   &
     &                    Name = "long_name", Values = SCMop%lname(i))
                  endif
               enddo

               ! End of file definition
               Status = Nf90_EndDef(NcID = NcID)

               ! Record the file ID so it can be used later for output.
               SCMop%strm(N)%op_unit=NcID

            else
               print*,'dump_streams_init ERROR: unknown format '//      &
     &              'for stream ',N
            endif

         endif                  ! SCMop%strm(N)%switch /= 0
      enddo                     ! do N=1,maxnstreams

 999  continue
      return
      END SUBROUTINE dump_streams_init
#endif
