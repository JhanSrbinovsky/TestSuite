#if defined(MAKEBC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! *********************************************************************
!+ Program MAKEBC : Top-level program to create boundary dataset
!                   from model analyses/dumps.
!


!+ Subroutine MAKEBC : Creates a boundary dataset from model dumps
!
! Subroutine Interface :

!+ Subroutine LOOP_OVER_DUMPS : Loop over dumps to get boundary data
!
! Subroutine Interface :
       subroutine loop_over_dumps (n_dumps,nhours,unit_no_bc,um_versn,  &
#include "argppx.h"
#include "arginfa.h"
     &            no_lams)

! Introduce type for times
      USE makebc_time_mod, ONLY: &
        time
 
      IMPLICIT NONE
!
! Description : Loop over the dumps and get the boundary conditions
!
! Method : For each dump, GET_BC is called to read in the data from
!          the dump and generate the boundary conditions.
!
! Current Code Owner : Dave Robinson, NWP
!
! History :
! Version    Date    Comment
! -------    ----    -------
!   4.5    18/02/98  Subroutine MAKEBC in 4.4 split into MAKEBC and
!                    LOOP_OVER_DUMPS. D. Robinson.
!   6.0     05/09/03   Upgrade to vn 6.0. R. Sempers
!   6.1     13/08/04 Code for reading Fields Files and changes
!                    to allow LBCs for multiple LAMs in single run.
!                    R. Sempers
!   6.1     22/10/04 Fix potential data/validity time problems with
!                    fields files
!                    R. Sempers
!   6.2    20/01/06  Tidy code, increase use of ereport move from
!                    icode to errorstatus
!   6.2     23/11/05 Allow reading from reinitialised fields
!                    files
!   6.2    14/02/05  Upgrades to allow generation of murk and
!                    PC2 lbcs
!                    R. Sempers

!   6.2     22/02/05 Obtain grid data from lookups instead of
!                    integer and real headers to allow use of sub-
!                    area STASH
!                    R. Sempers
!
! Code Description :
! Language : FORTRAN 77 + common extensions
! This code is written to UMDP3 v6 programming standards.
!
! Declarations :
!
! Global Variables :
!
#include "cmaxsize.h"
#include "parvars.h"
#include "typsize.h"
#include "chsunits.h"
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "ctime.h"
#include "cntlall.h"
#include "decomptp.h"
      Data TOT_LEVELS /-99/

#include "typinfa.h"
#include "cprintst.h"
#include "typptra.h"
#include "cintfa.h"

! Subroutine arguments
!   Scalar arguments with intent(in) :

      Integer n_dumps     ! No of model dumps
      Integer nhours      ! No of hours between dumps
      Integer unit_no_bc  ! Unit No for Boundary dataset
      Integer um_versn    ! UM Version Boundary dataset for
      integer :: no_lams  ! No of LAMS for which LBCs are to
                          ! be generated in this run
! Add cntlatm which contains l_pc2, microphysics and l_murk logicals
! to switch on pc2 and/or murk lbc generation
#include "cntlatm.h"

!   Array arguments with intent(in) :



      Integer errorstatus          ! Error code
      Character*80 cmessage  ! Error Message

!   Local scalars

      Integer unit_no        !  Unit no for input dump
      integer :: read_only=0  !  Input dumps - read only
      integer :: read_write=1 !  Output Boundary File - read & write
      integer :: len_env=6    !  Length of env. variable
      integer :: env_var=0    !  Indicator that filename is in env var
      Character*6 env        !  Env Variable for input dump filename


      Integer j,jdump        !  Loop indices
      integer :: inthd(46)
      Integer yy,mm,dd,hr,mn,ss,day_no  !  Time/date for first dump
      Integer elapsed_days   !  No of days elapsed
      Integer elapsed_secs   !  No of secs elapsed
      Integer len_actual     !  Length of data read in BUFFIN
      Integer HALOX,HALOY
      Integer p_rows

      Real a                 !  Return code from BUFFIN

!   Local dynamic arrays

      Integer fixhd(len_fixhd)  !  Fixed header from dump

! Declarations for vars that have been removed from TYPSIZE.
      integer :: p_field
      integer :: u_field

! Declarations for fieldsfile reading
      integer :: test(64)
      integer :: dim2_lookup
      integer :: len_io
      integer :: returncode
      integer :: open_loop, close_loop
      real :: a_io
      character(len=80),parameter :: routinename='loop_over_dumps'
      character(len=6)  :: bc_filename
! Max_progs gives the size required for indexing
! arrays in uniform manner
      integer :: max_progs

! Target_time is the lbc time required to select lookups
!         pointing to the correct data
! End_time holds the expected value of target_time at end
!         of LBC generation
      type(time) :: target_time
      type(time) :: end_time
! Flag to indicate all required data is present
      logical :: alldata
      integer :: iloop
      integer :: jcount
      integer, allocatable :: valid_starts(:)
      integer :: start_hour
      integer :: end_hour
!-  End of Header
! Set max_progs
      max_progs=20


! Initialise all items in target_time to 0
      target_time%year =0
      target_time%month =0
      target_time%day =0
      target_time%hour =0
      target_time%min =0
      target_time%sec =0
      target_time%day_no =0

!     Open Boundary Dataset
      write (6,*) ' '

!     Loop over number of LAMs, opening boundary file for each
      do open_loop=1,no_lams
        write(bc_filename,'(a,i1)')'BCFIL',open_loop
        if(printstatus >= prstatus_diag)then
          write(6,*)' bc_filename= ',bc_filename
        endif
! DEPENDS ON: file_open
        call file_open                                                  &
     &     (unit_no_bc+(open_loop-1),bc_filename,len_env,read_write,    &
     &      env_var,errorstatus)
        if (errorstatus /= 0) then
          write (6,*) 'Error in opening Boundary Dataset on unit no ',  &
     &    unit_no_bc
! DEPENDS ON: ereport
          call ereport(routinename,errorstatus,cmessage)
        endif
      enddo

!     Loop over model dumps
      do jdump=1,n_dumps
      write (6,*) ' '
      write (6,*) ' Processing dump no ',jdump

!     Unit number for this dump
      unit_no = jdump+30

!     Open the dump
      env = 'FILE  '
      write (env(5:6),'(I2)') unit_no
      write (6,*) ' '
! DEPENDS ON: file_open
      call file_open (unit_no,env,len_env,read_only,env_var,errorstatus)
      if (errorstatus /= 0) then
        write (6,*) 'Error in opening dump on unit no ',unit_no
! DEPENDS ON: ereport
        call ereport(routinename,errorstatus,cmessage)
      endif

!     Read in fixed header from this dump
! DEPENDS ON: setpos
      call setpos (unit_no,0,errorstatus)
      if (errorstatus >  0) then
        write (6,*) 'Error in SETPOS for Fixed Header.'
! DEPENDS ON: ereport
        call ereport(routinename,errorstatus,cmessage)
      endif

! DEPENDS ON: read_flh
      call read_flh (unit_no,fixhd,len_fixhd,errorstatus,cmessage)
      if (errorstatus >  0) then
        write (6,*) 'Error in READ_FLH for dump ',jdump
! DEPENDS ON: ereport
        call ereport(routinename,errorstatus,cmessage)
      endif

      if (fixhd(5) == 1)then
        write(6,*)'this is a Dump'
      elseif (fixhd(5) == 3)then
        write(6,*)'this is a FieldsFile'
      else
        cmessage= 'Invalid input file type'
        returncode=1
! DEPENDS ON: ereport
        call Ereport( RoutineName, ReturnCode, Cmessage )
      endif
!     For first dump only
!     Set model basis time to be date/time in first dump and use
!     to initialise variables required for time processing
      if (jdump == 1) then
        if(fixhd(150) > 0)then
! DEPENDS ON: setpos
          call setpos(unit_no,fixhd(150)-1,errorstatus)
          if (errorstatus >  0) then
            write (6,*) 'Error in SETPOS for lookups.'
! DEPENDS ON: ereport
            call ereport(routinename,errorstatus,cmessage)
          endif
! DEPENDS ON: buffin
          call buffin(unit_no, test, fixhd( 151 ),                      &
     &                        Len_IO, a_IO )
          if (a_IO /= -1.0) then
            write (6,*) 'Problem reading first lookup.'
            write (6,*) 'Return code from buffin = ',a_IO
            write (6,*) 'Length of data read by buffin = ',len_actual
            errorstatus=nint(a_IO)
! DEPENDS ON: ereport
            call ereport(routinename,errorstatus,cmessage)
          endif

          model_basis_time(1)=test(1)
          model_basis_time(2)=test(2)
          model_basis_time(3)=test(3)
          model_basis_time(4)=test(4)
          model_basis_time(5)=0
          model_basis_time(6)=0
        else
          write(6,*)'ERROR: No Lookups to get basis time from'
        endif
        model_basis_time(5)=0
        model_basis_time(6)=0

! Set the target_time from the model basis time
        target_time%year = model_basis_time(1)
        target_time%month= model_basis_time(2)
        target_time%day  = model_basis_time(3)
        target_time%hour = model_basis_time(4)
        target_time%min  = model_basis_time(5)
        target_time%sec  = model_basis_time(6)
! Set day_no directly from lookup
        target_time%day_no=test(6)

! Make sure earliest start time and tatest end time are chosen
        allocate(valid_starts(n_intf_a))
        valid_starts(:)=-1
        jcount=1
        do iloop=1,8
          if(a_intf_end_hr(iloop) /= 0)then
            valid_starts(jcount)=a_intf_start_hr(iloop)
            jcount=jcount+1
          endif
        enddo
        if(printstatus >= prstatus_diag)then
          write(6,*)'valid_starts=',valid_starts
        endif
        end_hour=maxval(a_intf_end_hr)
        start_hour=minval(valid_starts)
        deallocate(valid_starts)

! Call calc_new_time to generate validity time of first LBC
!         time to search for in the lookups
! DEPENDS ON: calc_new_time
          call calc_new_time(start_hour,                                &
     &                    0,0,                                          &
     &                    target_time%year,                             &
     &                    target_time%month,                            &
     &                    target_time%day,                              &
     &                    target_time%hour,                             &
     &                    target_time%min,                              &
     &                    target_time%sec,                              &
     &                    target_time%day_no,                           &
     &                    lcal360)


! Set the initial values for end_time from the model
! basis time
        end_time%year = model_basis_time(1)
        end_time%month= model_basis_time(2)
        end_time%day  = model_basis_time(3)
        end_time%hour = model_basis_time(4)
        end_time%min  = model_basis_time(5)
        end_time%sec  = model_basis_time(6)
! Set day_no directly from lookup
        end_time%day_no=test(6)

! Call calc_new_time to add correct ammount of time to make
! end_time correspond value of target_time expected after all
! LBCs have been created.
! DEPENDS ON: calc_new_time
          call calc_new_time(end_hour+a_intf_freq_hr,                   &
     &                    0+a_intf_freq_mn,0+a_intf_freq_sc,            &
     &                    end_time%year,                                &
     &                    end_time%month,                               &
     &                    end_time%day,                                 &
     &                    end_time%hour,                                &
     &                    end_time%min,                                 &
     &                    end_time%sec,                                 &
     &                    end_time%day_no,                              &
     &                    lcal360)

        i_year   = model_basis_time(1)
        i_month  = model_basis_time(2)
        i_day    = model_basis_time(3)
        i_hour   = model_basis_time(4)
        i_minute = model_basis_time(5)
        i_second = model_basis_time(6)

        basis_time_days = 0
        basis_time_secs = 0
! DEPENDS ON: time2sec
        call time2sec(i_year,i_month,i_day,i_hour,i_minute,i_second,    &
     &     basis_time_days,basis_time_secs,elapsed_days,elapsed_secs,   &
     &     lcal360)

        basis_time_days = elapsed_days
        basis_time_secs = elapsed_secs
        elapsed_days = 0
        elapsed_secs = 0
! DEPENDS ON: sec2time
        call sec2time(elapsed_days,elapsed_secs,                        &
     &                basis_time_days,basis_time_secs,                  &
     &                yy,mm,dd,hr,mn,ss,day_no,lcal360)

      endif

!     Remove negative dimensions, if any.
      do j=100,256
      if (fixhd(j) <  0) fixhd(j)=0
      enddo

!     Get header dimensions from Fixed Header
      a_len_inthd    = fixhd(101)
      a_len_realhd   = fixhd(106)
      a_len1_levdepc = fixhd(111)
      a_len2_levdepc = fixhd(112)
      a_len1_rowdepc = fixhd(116)
      a_len2_rowdepc = fixhd(117)
      a_len1_coldepc = fixhd(121)
      a_len2_coldepc = fixhd(122)
      a_len1_flddepc = fixhd(126)
      a_len2_flddepc = fixhd(127)
      a_len_extcnst  = fixhd(131)
      a_len_cfi1     = fixhd(141)
      a_len_cfi2     = fixhd(143)
      a_len_cfi3     = fixhd(145)
      a_len2_lookup  = fixhd(152)
      a_len_data     = fixhd(161)

!     Get length of data in this dump
      len_tot  = a_len_data

!     Read in Integer Constants for this dump
! DEPENDS ON: setpos
      call setpos(unit_no,fixhd(100)-1,errorstatus)
      if (errorstatus >  0) then
        write (6,*) 'Error in SETPOS for Integer Constants array.'
! DEPENDS ON: ereport
        call ereport(routinename,errorstatus,cmessage)
      endif

! DEPENDS ON: buffin
      call buffin (unit_no,inthd(1),a_len_inthd,len_actual,a)
      if (a /= -1.0) then
        write (6,*) 'Problem with reading Integer Constants array.'
        write (6,*) 'Return code from buffin = ',a
        write (6,*) 'Length of data read by buffin = ',len_actual
        errorstatus=nint(a)
! DEPENDS ON: ereport
        call ereport(routinename,errorstatus,cmessage)
      endif

!     Get model grid for this dump
      row_length = test(19)   ! from lookup
      rows       = test(18)   ! from lookup
      n_rows     = rows-1
      p_field    = row_length * rows
      u_field    = row_length * n_rows
      u_field_intfa = u_field
      write (6,*) ' u_field_intfa set to ',u_field_intfa

!     Get model levels for this dump
      model_levels = inthd(8)
      wet_levels   = inthd(9)
      tr_levels  = inthd(12)

      if(printstatus >= prstatus_normal)then
        write (6,*) ' '
        write (6,*) ' Model Grid/Levels in this dump. '
        write (6,*) ' row_length   = ',row_length
        write (6,*) ' rows         = ',rows
        write (6,*) ' n_rows       = ',n_rows
        write (6,*) ' model_levels = ',model_levels
        write (6,*) ' wet_levels   = ',wet_levels
        write (6,*) ' tr_levels  = ',tr_levels
        write (6,*) ' p_field    = ',p_field
        write (6,*) ' u_field    = ',u_field
      endif


! Set pointers into level dependent constants array
      jetatheta=1
      jetarho=jetatheta+model_levels+1

!     Ensure TR_LEVELS > 0 to prevent zero dynamic allocation.
      if (tr_levels <= 0) then
        tr_levels = 1
      endif

! Get decompostion information
! ----------------------------
! TOT_LEVELS not used in SX
! DEPENDS ON: decompose_smexe
      CALL DECOMPOSE_SMEXE( row_length, rows,                           &
     &                      0, 0, tot_levels)

!     Proceed to read model data from dump and get boundary conditions
! DEPENDS ON: get_bc
      call get_bc (jdump,unit_no,unit_no_bc,um_versn,                   &
     &  a_len2_lookup,                                                  &
     &  Fixhd_intfa,Inthd_intfa,Lookup_intfa,                           &
     &  Realhd_intfa,Levdepc_intfa,                                     &
     &  Rowdepc_intfa,Coldepc_intfa,                                    &
     &  Intf_akh,Intf_bkh,Intf_ak,Intf_bk,                              &
#include "argppx.h"
! Pass target_time and end time
     &  target_time,                                                    &
     &  end_time,                                                       &
! Pass indicator alldata
     &  alldata,                                                        &
     &  start_hour,                                                     &
! Pass nhours in
     &  nhours,                                                         &
     &  lbc_eta_theta,lbc_eta_rho,                                      &
     &  max_progs)


!     Close the dump
! DEPENDS ON: file_close
      call file_close (unit_no,env,len_env,env_var,0,errorstatus)
      if (errorstatus /= 0) then
        write (6,*) 'Error in closing dump on unit no ',unit_no
! DEPENDS ON: ereport
        call ereport(routinename,errorstatus,cmessage)
      endif

      enddo   !    End of loop over dumps

! Generate meaningful error message if target_time hasn't
! reached the expected value. Traps cases where data for an lbc
! time is missing. Call ereport.
      if(.not. alldata)then
        if(end_time%year  /= target_time%year .or.                      &
     &     end_time%month /= target_time%month .or.                     &
     &     end_time%day   /= target_time%day .or.                       &
     &     end_time%hour  /= target_time%hour .or.                      &
     &     end_time%min   /= target_time%min .or.                       &
     &     end_time%sec   /= target_time%sec .or.                       &
     &     end_time%day_no/= target_time%day_no                         &
     &    )then
          write(6,*)'Error - Data missing for time'
          write(6,*)target_time%year,' ',target_time%month,' ',         &
     &     target_time%day,' ',target_time%hour,':',                    &
     &     target_time%min,':',target_time%sec
          cmessage='ERROR data missing for one or more LBC time'
          errorstatus=1
! DEPENDS ON: ereport
          call ereport(routinename,errorstatus,cmessage)
        endif
      endif

      return
      END SUBROUTINE loop_over_dumps

!+ Subroutine GET_BC : Get boundary conditions from model dump
!
! Subroutine Interface :
#endif
