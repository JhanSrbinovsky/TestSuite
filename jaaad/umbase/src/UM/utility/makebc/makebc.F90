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
      subroutine makebc(ppxRecs,n_dumps,nhours,                         &
     &                  um_versn,sub_hr_int,no_lams)


      IMPLICIT NONE
!
! Description : Control routine for MAKEBC utility.
!
! Method : Read in DUMP2BOUND namelist. Call INTF_CTL which reads
!          in INTFCNST namelist. Loop over dumps to generate
!          boundary conditions.
!          No longer reads in DUMP2BOUND in this routine - frpz
!
! Current Code Owner : Dave Robinson, NWP
!
! History :
! Version    Date    Comment
! -------    ----    -------
!   4.4    10/10/97  Original Code
!   6.0    05/09/03  Upgraded code. R. Sempers
!   6.1    13/08/04  Changes made to allow multi LAM
!                    LBC generation. Read of DUMP2BOUND moved
!                    to MAIN_MAKEBC. TIMER call(s) added. R. Sempers.
!   6.2    14/02/05  Upgrades to allow generation of murk and
!                    PC2 lbcs
!                    R. Sempers
!   6.2    01/12/05  Add write statements after call to intf_ctl
!                    to write out intf_ExtHalo_NS and EW.
!                    R. Sempers
!    6.2  20/01/06  Tidy up code. Increase use of ereport
!                   Use errorstatus in preference to icode
!
! Code Description :
! Language : FORTRAN 77 + common extensions
! This code is written to UMDP3 v6 programming standards.
!
! Declarations :
!
! Global Variables :
!
#include "parparm.h"
#include "typsize.h"
#include "csubmodl.h"
#include "cmaxsize.h"
#include "chsunits.h"
#include "chistory.h"
#include "ctime.h"
#include "clookadd.h"
#include "c_mdi.h"
#include "ccontrol.h"
#include "cintfa.h"
#include "typinfa.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "typduma.h"
#include "typptra.h"
#include "cprintst.h"

! Subroutine arguments
!   Scalar arguments with intent(in) :

!   Array arguments with intent(in) :

!   Scalar arguments with intent(inout) :

!   Array arguments with intent(inout) :

!   Scalar arguments with intent(out) :

      character(len=80), parameter :: routinename='makebc'
      Integer errorstatus            !  Error code
      Character*80 cmessage    !  Error Message

!   Local scalars :

      Integer j,jintf          !  Loop indices
      Integer irow_number      !  Row number, required for GETPPX
      Integer internal_model   !  Internal Model Identifier
      Integer ME_GC,NPROC_GC

      Character*100 PAREXE_ENV  ! hold name of the // exec script


!     Required for I/O
      Integer unit_no_bc  !  Unit No for boundary dataset



!     DUMP2BOUND namelist for MAKEBC Program
      Integer  n_dumps      !  No of model dumps
      Integer  nhours       !  No of hours between dumps
      Integer  um_versn     !  UM Version Boundary Dataset for

      integer :: sub_hr_int   ! No of lbcs desired per hour
      integer :: no_lams      ! No of lams to generate lbcs for
                              ! in this run

!     lcal360  defined in CNTLALL/CNTLATM

!-  End of Header

      PAREXE_ENV=' '

      CALL GC_INIT(PAREXE_ENV,ME_GC,NPROC_GC)

      internal_model = 1

!     Initialise LLBOUTim in CNTLGEN
      LLBOUTim(internal_model)=.true.

!     Initialise variables to nullify A_STEPS_PER_HR in INTF_CTL
      if(sub_hr_int == 0)then
        STEPS_PER_PERIODim(internal_model) = 1
      else
        STEPS_PER_PERIODim(internal_model) = sub_hr_int
      endif
      write(6,*)'steps per hour=',  STEPS_PER_PERIODim(internal_model)
      SECS_PER_PERIODim(internal_model)  = 3600

!     Initialise STEPIM to correspond to first dump
      STEPim(internal_model)=0

!     Use UM Unit No 140-147 for Atmos Boundary Datasets 1-8
      unit_no_bc = 140

!     Initialise variables in CNTLALL for this unit no.
!     Reinitialising of boundary dataset not supported yet.
      TYPE_LETTER_1(140:147) = 'b'
      FT_STEPS     (140:147) = 0
      FT_FIRSTSTEP (140:147) = 0
      FT_OUTPUT    (140:147) = 'N'

! DEPENDS ON: timer
      call timer ( 'Intf_Ctl', 3 )

!     Get model grid for which boundary conditions are required
! DEPENDS ON: intf_ctl
      call intf_ctl (                                                   &
#include "arginfa.h"
#include "argduma.h"
#include "argptra.h"
     &               errorstatus,cmessage)

! DEPENDS ON: timer
      call timer ( 'Intf_Ctl', 4 )

!     Print out INTFCNTL namelist variables (Read in intf_ctl)
      if(printstatus >= prstatus_diag)then
        write (6,*) ' '
        write (6,*) ' Namelist INTFCNSTA read in'
        do jintf=1,n_intf_a
          write (6,*)' '
          write (6,*) ' For area ',jintf
          write (6,*) ' a_intf_start_hr  ',A_INTF_START_HR(JINTF)
          write (6,*) ' a_intf_freq_hr   ',A_INTF_FREQ_HR(JINTF)
          write (6,*) ' a_intf_end_hr    ',A_INTF_END_HR(JINTF)
          write (6,*) ' intf_p_rows      ',INTF_P_ROWS(JINTF)
          write (6,*) ' intf_row_length  ',INTF_ROW_LENGTH(JINTF)
          write (6,*) ' intf_p_levels    ',INTF_P_LEVELS(JINTF)
          write (6,*) ' intf_q_levels    ',INTF_Q_LEVELS(JINTF)
          write (6,*) ' intf_tr_levels   ',INTF_TR_LEVELS(JINTF)
          write (6,*) ' intf_firstlat    ',INTF_FIRSTLAT(JINTF)
          write (6,*) ' intf_firstlong   ',INTF_FIRSTLONG(JINTF)
          write (6,*) ' intf_nsspace     ',INTF_NSSPACE(JINTF)
          write (6,*) ' intf_ewspace     ',INTF_EWSPACE(JINTF)
          write (6,*) ' intf_polelat     ',INTF_POLELAT(JINTF)
          write (6,*) ' intf_polelong    ',INTF_POLELONG(JINTF)
          write (6,*) ' intf_pack        ',INTF_PACK(JINTF)
          write (6,*) ' intfwidtha       ',INTFWIDTHA(JINTF)
          write (6,*) ' intf_ExtHalo_NS  ',INTF_EXTHALO_NS(JINTF)
          write (6,*) ' intf_ExtHalo_EW  ',INTF_EXTHALO_EW(JINTF)
          write (6,*) ' intf_vert_interp ',INTF_VERT_INTERP(JINTF)
        enddo
        write(6,*)' '
      endif


!     Read StashMaster file
      irow_number=0
! DEPENDS ON: getppx
      call getppx (22,2,'STASHmaster_A',irow_number,                    &
#include "argppx.h"
     &  errorstatus,cmessage)

      if (errorstatus >  0) then
        write (6,*) 'Error in GETPPX.'
! DEPENDS ON: ereport
        call ereport(routinename,errorstatus,cmessage)
      endif

! DEPENDS ON: loop_over_dumps
       call loop_over_dumps (n_dumps,nhours,unit_no_bc,um_versn,        &
#include "argppx.h"
#include "arginfa.h"
     &      no_lams)

 9999  continue

       return
       END SUBROUTINE makebc

!+ Subroutine LOOP_OVER_DUMPS : Loop over dumps to get boundary data
!
! Subroutine Interface :

!+ Subroutine GET_BC : Get boundary conditions from model dump
!
! Subroutine Interface :
#endif
