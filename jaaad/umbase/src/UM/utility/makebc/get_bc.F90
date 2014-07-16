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

!+ Subroutine GET_BC : Get boundary conditions from model dump
!
! Subroutine Interface :
      subroutine get_bc (jdump,unit_no,unit_no_bc,um_versn,             &
     &  len2_lookup,                                                    &
     &  Fixhd_intfa,Inthd_intfa,Lookup_intfa,                           &
     &  Realhd_intfa,Levdepc_intfa,                                     &
     &  Rowdepc_intfa,Coldepc_intfa,                                    &
     &  Intf_akh,Intf_bkh,Intf_ak,Intf_bk,                              &
#include "argppx.h"
! Pass target_time, end_time and alldata
     &  target_time,                                                    &
     &  end_time,                                                       &
     &  alldata,                                                        &
     &  start_hour,                                                     &
! Pass nhours in
     &  nhours,                                                         &
     &  lbc_eta_theta,lbc_eta_rho,                                      &
     &  max_progs)

! Introduce type for times
      USE makebc_time_mod, ONLY: &
        time

      IMPLICIT NONE
!
! Description : Get boundary conditions from model dump
!
! Method : For each dump, read in data through UM_READDUMP and
!          generate boundary conditions through GEN_INTF. Also
!          calls IN_INTF to initialise boundary dataset.
!
! Current Code Owner : Dave Robinson, NWP
!
! Code Description :
! Language : FORTRAN 77 + common extensions
! This code is written to UMDP3 v6 programming standards.
!
! Declarations :
!
! Global Variables :
!
#include "csubmodl.h"
#include "cmaxsize.h"
#include "cintfa.h"
#include "parvars.h"
#include "typsize.h"
! Replace use of typd1.h, allowing use of allocatable arrays
! for D1, ID1, and LD1
#include "caltsubm.h"
      REAL, allocatable    ::  D1(:) ! D1
      LOGICAL, allocatable :: LD1(:) ! Logical D1
      INTEGER, allocatable :: ID1(:) ! Integer D1
#include "d1_addr.h"
      ! D1 addressing array and number of objects in each submodel
      INTEGER :: D1_ADDR(D1_LIST_LEN,N_OBJ_D1_MAX,                      &
     &  ALT_N_SUBMODEL_PARTITION)

      INTEGER :: NO_OBJ_D1(ALT_N_SUBMODEL_PARTITION_MAX)
      COMMON/common_D1_ADDRESS/ NO_OBJ_D1

! Variable for calculating correct length for d1 type arrays
      integer :: calc_length
#include "typduma.h"
#include "typinfa.h"
#include "typsts.h"
#include "typptra.h"
#include "chsunits.h"
#include "ccontrol.h"
#include "ppxlook.h"
#include "cprintst.h"
#include "clookadd.h"

#include "ctime.h"

! Subroutine arguments
!   Scalar arguments with intent(in) :

      Integer jdump            !  No of dump being processed.
      Integer unit_no          !  Unit No for input dump
      Integer unit_no_bc       !  Unit No for boundary dataset
      Integer um_versn         !  UM Version Boundary Dataset for
      Integer len2_lookup
! Define max_progs
      integer :: max_progs

      Integer :: errorstatus
      Character*80 cmessage    !  Error Message
      Character*80 routinename

!   Array arguments with intent(out) :

!   Local parameters :

!   Local scalars :

      Integer submodel_id      !  Sub model identifier
      Integer internal_model   !  Internal model identifier
      integer, parameter :: idummy=1
      integer num_levels(NITEMS,N_INTERNAL_MODEL)
      integer num_levels_pointers(max_progs)
      Logical readhdr          !  T : Read headers from dump
      integer :: nftout, jintf
      integer n_flds,fld,disused,start_block
      Integer mpp_lookup(mpp_len1_lookup,len2_lookup)
                               !  dummy lookup


      integer :: u_field

      integer ipt , jrow, i, j, ipt0, no_rows, no_levels
      integer col, row, lev
      integer halo_x, halo_y
      real, allocatable :: data_in (:,:,:)
      real, allocatable :: data    (:,:,:)

! Variables and equivalence to source grid data from lookups
      integer :: ndx
      real :: r_db(6)

      integer d1_pointers (max_progs)
      integer :: halo_size_out(max_progs,3)
! Variables for reading wgdos
      integer :: wgdos_expand

! Declare target_time - validity time of required data
      type(time) :: target_time
      type(time) :: end_time

! Target_time and end_time in days and seconds for comparison
      integer :: target_days, target_secs
      integer :: end_days, end_secs
      integer :: start_hour
! Declare next_dump which allows program to move on to
! next dump file
      logical :: next_dump

! Declare flag for all required data being present
      logical :: alldata

! Maximum number of times to loop over in dump. Prevents infinite loop.
      integer :: times_loop

! Declare nhours and advance
      integer :: nhours
      integer, save :: advance

!-  End of Header
! Set routine name for error reporting
      routinename='get_bc'
! Initialise advance to 0
      if(jdump == 1)then
        advance = 0
      endif

      submodel_id    = 7    ! indicates SX
      internal_model = 1

! Inialise next_dump to .false. on entering get_bc
      next_dump=.false.
! Initialise flag for all required data being present to .false.
      alldata=.false.
! Initialise calc_length to 0
      calc_length=0
      
!     Go to start of dump
! DEPENDS ON: setpos
      call setpos (unit_no,0,errorstatus)
      if (errorstatus >  0) then
        write (6,*) 'Error in SETPOS for Model Dump.'
! DEPENDS ON: ereport
        call ereport(routinename,errorstatus,cmessage)
      endif

!Read in headers from the dump

! DEPENDS ON: timer
        call timer ( 'ReadHead', 3 )


! DEPENDS ON: readhead
        CALL READHEAD(unit_no,a_FIXHD,LEN_FIXHD,                        &
     &                a_INTHD,a_LEN_INTHD,                              &
     &                a_REALHD,a_LEN_REALHD,                            &
     &                a_LEVDEPC,a_LEN1_LEVDEPC,a_LEN2_LEVDEPC,          &
     &                a_ROWDEPC,a_LEN1_ROWDEPC,a_LEN2_ROWDEPC,          &
     &                a_COLDEPC,a_LEN1_COLDEPC,a_LEN2_COLDEPC,          &
     &                a_FLDDEPC,a_LEN1_FLDDEPC,a_LEN2_FLDDEPC,          &
     &                a_EXTCNST,a_LEN_EXTCNST,                          &
     &                a_DUMPHIST,LEN_DUMPHIST,                          &
     &                a_CFI1,a_LEN_CFI1,                                &
     &                a_CFI2,a_LEN_CFI2,                                &
     &                a_CFI3,a_LEN_CFI3,                                &
     &                a_LOOKUP,LEN1_LOOKUP,a_LEN2_LOOKUP,               &
     &                a_LEN_DATA,                                       &
#include "argppx.h"
     &                START_BLOCK,ERRORSTATUS,CMESSAGE)

! DEPENDS ON: timer
        call timer ( 'ReadHead', 4 )

        IF (ERRORSTATUS  /=  0) THEN
          WRITE(6,*) 'GET_BC : Error reading dump header ',             &
     &               'on unit ',unit_no
          WRITE(6,*) 'Return code from READHEAD was ',ERRORSTATUS,      &
     &               ' and error message was ',CMESSAGE
          ERRORSTATUS=3
          CMESSAGE='Error reading dump header'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,Errorstatus,Cmessage)
        ENDIF

! Initiate loop to loop over times in dump (max 900 loops)
      do times_loop=1,900
! DEPENDS ON: time2sec
           call time2sec(                                               &
     &                    target_time%year,                             &
     &                    target_time%month,                            &
     &                    target_time%day,                              &
     &                    target_time%hour,                             &
     &                    target_time%min,                              &
     &                    target_time%sec,                              &
     &      0,0,target_days,target_secs,                                &
     &      lcal360)
! DEPENDS ON: time2sec
           call time2sec(                                               &
     &                    end_time%year,                                &
     &                    end_time%month,                               &
     &                    end_time%day,                                 &
     &                    end_time%hour,                                &
     &                    end_time%min,                                 &
     &                    end_time%sec,                                 &
     &      0,0,end_days,end_secs,                                      &
     &      lcal360)
      if(target_days > end_days .or.                                    &
     &    (target_days == end_days .and. target_secs >= end_secs))then
        if(printstatus >= prstatus_diag)then
          write(6,*)'done last time'
          write(6,*)'end_time'
          write(6,*) end_time%year,                                     &
     &               end_time%month,                                    &
     &               end_time%day,                                      &
     &               end_time%hour,                                     &
     &               end_time%min,                                      &
     &               end_time%sec,                                      &
     &               end_time%day_no
          write(6,*)'target_time'
          write(6,*) target_time%year,                                  &
     &               target_time%month,                                 &
     &               target_time%day,                                   &
     &               target_time%hour,                                  &
     &               target_time%min,                                   &
     &               target_time%sec,                                   &
     &               target_time%day_no
          write(6,*)'end_days=',end_days,'end_secs=',end_secs
          write(6,*)'target_days=',target_days,                         &
     &                           'target_secs=',target_secs

        endif
        alldata=.true.
        exit
      endif
! DEPENDS ON: set_ppindex
       call set_ppindex(jorog,ju,jv,jw,                                 &
     &           jrho,jtheta,jq,jqcl,jqcf,                              &
     &           jexner_rho_levels,                                     &
     &           ju_adv,jv_adv,jw_adv,                                  &
     &           jqcf2,jqrain,jqgraup,                                  &
     &           jmurk,                                                 &
     &           jcf_bulk,jcf_liquid,jcf_frozen,                        &
     &           jtracer,                                               &
     &           nitems,ppindex,                                        &
     &           len1_lookup,a_len2_lookup,a_lookup,                    &
     &           num_levels,num_levels_pointers,                        &
     &           d1_pointers,                                           &
#include "argppx.h"
! Next_dump logical and time index passed to set_ppindex
     &           next_dump,target_time,                                 &
! Add new calculated_length for d1 allocation
     &           calc_length,                                           &
     &           halo_size_out,                                         &
     &           l_pc2,l_murk,max_progs,                                &
     &           l_mcr_qcf2,l_mcr_qrain,l_mcr_qgraup)



! If next_dump is true, reset and exit the loop
      if(next_dump)then
        next_dump=.false.
        exit
      endif

! Set len_tot to the value calculated in set_ppindex
! and reset calculated length to zero ready for next time
      len_tot = calc_length
      calc_length = 0

! Allocate the D1 and equivalent arrays
      allocate(d1(len_tot))
      allocate(id1(len_tot))
      allocate(ld1(len_tot))

! Calculate stepim
!     Set STEPim(1) for this dump (controlled by nhours)
      STEPim(1) = start_hour                                            &
     &              + advance*nhours

      advance=advance+1
      if(printstatus >= prstatus_diag)then
        write(6,*) ' '
        write(6,*) ' pointers from set_ppindex'
        write(6,*) ' jorog =  ',jorog
        write(6,*) ' ju=      ',ju
        write(6,*) ' jv=      ',jv
        write(6,*) ' jw=      ',jw
        write(6,*) ' jrho=    ',jrho
        write(6,*) ' jtheta=  ',jtheta
        write(6,*) ' jq=      ',jq
        write(6,*) ' jqcl=    ',jqcl
        write(6,*) ' jqcf=    ',jqcf
        write(6,*) ' jexner_rho_levels=  ',jexner_rho_levels
        write(6,*) ' ju_adv=  ',ju_adv
        write(6,*) ' jv_adv=  ',jv_adv
        write(6,*) ' jw_adv=  ',jw_adv
        write(6,*) ' jqcf2=   ',jqcf2
        write(6,*) ' jqrain=  ',jqrain
        write(6,*) ' jqgraup= ',jqgraup
        write(6,*) ' jmurk=   ',jmurk
        write(6,*) ' jcf_bulk=  ',jcf_bulk
        write(6,*) ' jcf_liquid=',jcf_liquid
        write(6,*) ' jcf_frozen=',jcf_frozen

        write(6,*) ' jtracer= ',jtracer(1,1)

        write(6,*)' '
        write(6,*)' ppindex from set_ppindex'
        write(6,*)' ppindex(2) = ',ppindex(2,1)
        write(6,*)' ppindex(3) = ',ppindex(3,1)
        write(6,*)' ppindex(4) = ',ppindex(4,1)
        write(6,*)' ppindex(10) = ',ppindex(10,1)
        write(6,*)' ppindex(12) = ',ppindex(12,1)
        write(6,*)' ppindex(33) = ',ppindex(33,1)
        write(6,*)' ppindex(150) =',ppindex(150,1)
        write(6,*)' ppindex(253) =',ppindex(253,1)
        write(6,*)' ppindex(254) =',ppindex(254,1)
        write(6,*)' ppindex(255) =',ppindex(255,1)
        write(6,*)' ppindex(256) =',ppindex(256,1)
        write(6,*)' ppindex(257) =',ppindex(257,1)
        write(6,*)' ppindex(258) =',ppindex(258,1)
        write(6,*)' ppindex(271) =',ppindex(271,1)
        write(6,*)' ppindex(272) =',ppindex(272,1)
        write(6,*)' ppindex(273) =',ppindex(273,1)
        write(6,*)' ppindex(90) =',ppindex(90,1)
        write(6,*)' ppindex(266) =',ppindex(266,1)
        write(6,*)' ppindex(267) =',ppindex(267,1)
        write(6,*)' ppindex(268) =',ppindex(268,1)
        write(6,*)' '
        write(6,*)' num_levels from set_ppindex'
        write(6,*)' num_levels(2) = ',num_levels(2,1)
        write(6,*)' num_levels(3) = ',num_levels(3,1)
        write(6,*)' num_levels(4) = ',num_levels(4,1)
        write(6,*)' num_levels(10) = ',num_levels(10,1)
        write(6,*)' num_levels(12) = ',num_levels(12,1)
        write(6,*)' num_levels(33) = ',num_levels(33,1)
        write(6,*)' num_levels(150) =',num_levels(150,1)
        write(6,*)' num_levels(253) =',num_levels(253,1)
        write(6,*)' num_levels(254) =',num_levels(254,1)
        write(6,*)' num_levels(255) =',num_levels(255,1)
        write(6,*)' num_levels(256) =',num_levels(256,1)
        write(6,*)' num_levels(257) =',num_levels(257,1)
        write(6,*)' num_levels(258) =',num_levels(258,1)
        write(6,*)' num_levels(271) =',num_levels(271,1)
        write(6,*)' num_levels(272) =',num_levels(272,1)
        write(6,*)' num_levels(273) =',num_levels(273,1)
        write(6,*)' num_levels(90) =',num_levels(90,1)
        write(6,*)' num_levels(266) =',num_levels(266,1)
        write(6,*)' num_levels(267) =',num_levels(267,1)
        write(6,*)' num_levels(268) =',num_levels(268,1)
      endif

      n_flds=a_FIXHD(152)

! For reading in WGDOS files, set wgdos_expand
      wgdos_expand=1

      do fld=1,max_progs
        if(num_levels_pointers(fld) /= -1)then
          no_rows   = a_lookup(18,ppindex(num_levels_pointers(fld),1))
          no_levels = num_levels(num_levels_pointers(fld),1)
          Allocate ( data_in (1:row_length, 1:no_rows, no_levels) )

! DEPENDS ON: timer
          call timer ( 'ReadFlds', 3 )

! DEPENDS ON: readflds
          call READFLDS(unit_no,                                        &
     &                  no_levels,                                      &
     &                  ppindex(num_levels_pointers(fld),1),            &
     &                  a_lookup,                                       &
     &                  len1_lookup,                                    &
     &                  data_in,                                        &
     &                  disused,a_fixhd,                                &
#include "argppx.h"
! Pass in wgdos_expand (called expand in readflds) for reading
! in WGDOS files
     &                  wgdos_expand,                                   &
     &                  ERRORSTATUS,CMESSAGE)

! DEPENDS ON: timer
          call timer ( 'ReadFlds', 4 )
          if(errorstatus /= 0) then
            write(6,*)'Problem in readflds called from get_bc'
! DEPENDS ON: ereport
            call ereport(routinename,errorstatus,cmessage)
          endif

!     Allocate grid with haloes


          halo_x =  halo_size_out(fld,1) !1
          halo_y =  halo_size_out(fld,2) !1
          Allocate ( data (1-halo_x : row_length+halo_x ,               &
     &                     1-halo_y : no_rows+halo_y,                   &
     &                     no_levels ) )

          data(:,:,:) = 0.0

!     Now copy into grid with haloes.

          Do lev = 1, no_levels
            Do row = 1, no_rows
              Do col = 1, row_length
                data(col,row,lev) = data_in(col,row,lev)
              End Do
            End Do
          End Do

!     Now copy into D1

          ipt0 = d1_pointers(fld) - 1

          ipt  = 0
          Do lev = 1, no_levels
            Do row = 1-halo_y, no_rows+halo_y
              Do col = 1-halo_x, row_length+halo_x
                ipt = ipt + 1
                d1(ipt0+ipt) = data(col,row,lev)
              End Do
            End Do
          End Do

        ipt=0


          deallocate ( data_in)
          deallocate ( data)


        endif
      enddo

!     Set up the Headers in the boundary dataset.

!     Ensure that A_REALHD defines the right grid, even for fieldfiles
!     with sub-areas. Use values from orography lookup.
      ndx = ppindex(33,1)
      r_db(1) = transfer(a_lookup(BDX, ndx),r_db(1))
      r_db(2) = transfer(a_lookup(BDY, ndx),r_db(2))
      r_db(3) = transfer(a_lookup(BZY, ndx),r_db(3))
      r_db(4) = transfer(a_lookup(BZX, ndx),r_db(4))
      r_db(5) = transfer(a_lookup(BPLAT, ndx),r_db(5))
      r_db(6) = transfer(a_lookup(BPLON, ndx),r_db(6))

      A_REALHD(rh_deltaEW) = r_db(1)
      A_REALHD(rh_deltaNS) = r_db(2)
      A_REALHD(rh_baselat) = r_db(3) + r_db(2)
      A_REALHD(rh_baselong)= r_db(4) + r_db(1)
      A_REALHD(rh_rotlat)  = r_db(5)
      A_REALHD(rh_rotlong) = r_db(6)

      if(printstatus >= prstatus_diag)then
        write(6,*) 'Redefined A_REALHD:'
        write(6,*)'rh_baselat=',rh_baselat
        write(6,*)'rh_baselong=',rh_baselong
        write(6,*) 'A_REALHD(rh_deltaEW)=',A_REALHD(rh_deltaEW)
        write(6,*) 'A_REALHD(rh_deltaNS)=',A_REALHD(rh_deltaNS)
        write(6,*) 'A_REALHD(rh_baselat)=',A_REALHD(rh_baselat)
        write(6,*) 'A_REALHD(rh_baselong)=',A_REALHD(rh_baselong)
        write(6,*) 'A_REALHD(rh_rotlat)=',A_REALHD(rh_rotlat)
        write(6,*) 'A_REALHD(rh_rotlong)=',A_REALHD(rh_rotlong)
      endif

!    Check if vertical interpolation required.
! DEPENDS ON: lbc_chk_vert_interp
          Call LBC_Chk_Vert_Interp (                                    &
#include "arginfa.h"
#include "argduma.h"
#include "argptra.h"
     &    idummy )

        write (6,*) ' '
        write (6,*) ' Dump No ',jdump,' : Calling IN_INTF.'

        do jintf = 1,n_intf_a

! DEPENDS ON: intf_unit
          call intf_unit (1,jintf,nftout)
          write (6,*) ' nftout from intf_unit ',nftout



          if(jdump == 1)then
! DEPENDS ON: in_intf
            call in_intf (                                              &
#include "argduma.h"
#include "arginfa.h"
     &                  unit_no_bc+jintf-1,errorstatus,cmessage)


            if (errorstatus >  0) then
              write (6,*) 'Error in IN_INTF.'
! DEPENDS ON: ereport
              call ereport(routinename,errorstatus,cmessage)
            endif
          endif

!       Set UM Version for Boundary dataset
          fixhd_intfa(12,jintf) = um_versn

        enddo  !jintf = 1,n_intf_a



      write (6,*) ' '
      write (6,*) ' Dump No ',jdump,' : Calling GEN_INTF.'

! DEPENDS ON: timer
      call timer ( 'Gen_Intf', 3 )

!     Call GEN_INTF to generate boundary conditions for this dump
! DEPENDS ON: gen_intf
      call gen_intf (                                                   &
#include "argd1.h"
#include "argduma.h"
#include "argsts.h"
#include "argptra.h"
#include "arginfa.h"
#include "argppx.h"
     &              internal_model,errorstatus,cmessage)

! DEPENDS ON: timer
      call timer ( 'Gen_Intf', 4 )

      if (errorstatus >  0) then
        write (6,*) 'Error in GEN_INTF.'
! DEPENDS ON: ereport
        call ereport(routinename,errorstatus,cmessage)

      endif

 9999 continue
! Deallocate the D1 and equivalent arrays
      deallocate(d1)
      deallocate(id1)
      deallocate(ld1)
! Call calc_new_time to generate the next validity time
!         search for in the lookups
! DEPENDS ON: calc_new_time
          call calc_new_time(a_intf_freq_hr,                            &
     &                    a_intf_freq_mn,a_intf_freq_sc,                &
     &                    target_time%year,                             &
     &                    target_time%month,                            &
     &                    target_time%day,                              &
     &                    target_time%hour,                             &
     &                    target_time%min,                              &
     &                    target_time%sec,                              &
     &                    target_time%day_no,                           &
     &                    lcal360)
! End of infinite loop over times in dump
      enddo
      return
      END SUBROUTINE get_bc
#endif
