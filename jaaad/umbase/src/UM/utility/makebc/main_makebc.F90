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
      Program MAIN_MAKEBC

      IMPLICIT NONE
!
! Description : Create a boundary dataset from UM model analyses
!               or dumps.
!
! Method : For each dump, boundary conditions are generated through
!          GEN_INTF for the area specified in the INFTCNST namelist.
!          This routine initialises various variables in TYPSIZE
!          before it can be used in the lower routines.
!
! Current Code Owner : Dave Robinson, NWP
!
! History :
! Version    Date    Comment
! -------    ----    -------
!   4.4    10/10/97  Original Code
!   4.5    07/08/98  Call new subroutine LOOP_OVER_DUMPS from MAKEBC.
!                    Adapt to new 4.5 changes. Use new unit number 140.
!                    Call new routine DERV_INTF_A. Rename CINTF to
!                    CINTFA. Read in env var UM_SECTOR_SIZE.
!                    D. Robinson.
!LL   5.0  3/6/99    Remove ICODE and CMESSAGE from UM_READDUMP
!LL                                                        P.Burton
!LL  5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!    6.0  05/09/03  Upgrade the code to version 6.0   R.Sempers
!LL  6.0  02/07/03  Make MPP as the only option for makebc
!LL                                              E.Leung
!    6.1  13/08/04  Upgrade to version 6.1.
!                   DUMP2BOUND now read in MAIN_MAKEBC and
!   6.2    20/01/06  Tidy code, increase use of errorstatus, move
!                    from icode to errorstatus
!                   calls to INIT_PRINTSTATUS and TIMER added.
!                   R. Sempers
!    6.2  14/02/05  Upgrades to allow generation of murk,
!                   microphysics and PC2 lbcs
!                   R. Sempers
!    6.2  25/05/05  Move include SX_DT to BLKDATA and remove duplicate
!                   declarations of TOT_LEVELS and ROW_LENGTH. P.Dando
!
!
! Code Description :
! Language : FORTRAN 77 + common extensions
! This code is written to UMDP3 v6 programming standards.
!
! Declarations :
!
      Integer internal_model   !  Internal Model Identifier
      Integer ppxRecs          !  No of stashmaster records
      Integer icode            !  Error code
      integer :: errorstatus

      Character*80 cmessage    !  Error Message
      character(len=80),parameter :: routinename='main_makebc'
      Character*8  c_um_sector_size  ! Char variable to read env var

#include "parparm.h"
#include "typsize.h"
#include "csubmodl.h"
#include "cntl_io.h"
#include "typstsz.h"
#include "cprintst.h"

      integer :: u_field

      CHARACTER(LEN=*), PARAMETER :: ProgName = "MakeBC"

#include "chsunits.h"
#include "cntlall.h"

! new bits shifting dump2bound read to this routine
!     DUMP2BOUND namelist for MAKEBC Program
      Integer  n_dumps      !  No of model dumps
      Integer  nhours       !  No of hours between dumps
      Integer  um_versn     !  UM Version Boundary Dataset for

      integer :: sub_hr_int   ! No of lbcs desired per hour
      integer :: no_lams      ! No of lams to generate lbcs for
                              ! in this run
! Add cntlatm which contains l_pc2 and l_murk logicals to switch
! on pc2 and/or murk lbc generation
#include "cntlatm.h"

!     lcal360  defined in CNTLALL/
      NAMELIST /DUMP2BOUND/ n_dumps,nhours,um_versn,lcal360,            &
     &no_lams,sub_hr_int                                                &
     &,l_pc2,l_murk,l_mcr_qcf2,l_mcr_qrain,l_mcr_qgraup

! Number of optional lbcs (ie pc2, murk) required
      integer :: num_optional_lbcs_out

!-  End of header

! DEPENDS ON: timer
      Call Timer ( ProgName, 1 )
! DEPENDS ON: initprintstatus
      Call InitPrintStatus


      write (6,*) ' ##########################################'
      write (6,*) ' Running MAKEBC Utility to create a'
      write (6,*) ' Boundary Dataset from Model Output'
      write (6,*) ' ##########################################'
      write (6,*) ' '

      icode = 0

!     Only Atmosphere Model catered for
      n_internal_model = 1
      internal_model = 1
      internal_model_index(internal_model) = 1

!     Determine no of Atmos records in STASHmaster file
      ppxRecs=1
! DEPENDS ON: hdppxrf
      call hdppxrf(22,'STASHmaster_A',ppxRecs,icode,cmessage)
      if (icode >  0) then
        write (6,*) 'Error in HDPPXRF for STASHmaster_A.'
! DEPENDS ON: ereport
        call ereport(routinename,icode,cmessage)
      endif

!     Get the current sector size for disk I/O
      CALL FORT_GET_ENV('UM_SECTOR_SIZE',14,c_um_sector_size,8,icode)
      IF (icode  /=  0) THEN
        WRITE(6,*) ' Warning : Environment variable UM_SECTOR_SIZE',    &
     &             ' has not been set.'
        WRITE(6,*) ' Setting UM_SECTOR_SIZE to 2048'
        um_sector_size=2048
      ELSE
        READ(c_um_sector_size,'(I4)') um_sector_size
        write (6,*) ' '
        write (6,*) ' UM_SECTOR_SIZE is set to ',um_sector_size
      ENDIF

!     Initialise variables in TYPSIZE
      nsects = 20
      nitems = 512
      n_req_items = 20
      n_ppxrecs = 20
      totitems = 20
      nsttims = 20
      nsttabl = 20
      num_stash_pseudo = 1
      num_pseudo_lists = 1
      nstash_series_records = 1
      nstash_series_block = 1
      mos_mask_len = 1

!     Dimensions of Headers in Boundary dataset
!     Integer/Real Constants
      pp_len_inthd = 46
      pp_len_realhd = 38
      a_len_inthd=pp_len_inthd
      a_len_realhd=pp_len_realhd

!     Level Dependent Constants array (Second dimension)
      intf_len2_levdepc = 4
!     Row/Col Dependent Constants array (Second dimension)
      intf_len2_rowdepc = 2
      intf_len2_coldepc = 2
      
! Read in the DUMP2BOUND namelist
!     Defaults for DUMP2BOUND namelist
      n_dumps  = 0
      nhours   = 0
      um_versn = 601
      lcal360  = .false.
      sub_hr_int = 0
      no_lams = 1
      l_pc2 = .false.
      l_murk = .false.
      l_mcr_qcf2 = .false.
      l_mcr_qrain = .false.
      l_mcr_qgraup= .false.

!     Read in DUMP2BOUND namelist and print
      rewind 5
      read  (5,DUMP2BOUND)
      write (6,*) ' '
      write (6,*) 'Namelist DUMPBOUND read in '
      write (6,DUMP2BOUND)

      if (no_lams > 8)then
        cmessage = 'Too many LAMs (>8) requested'
        errorstatus = 1
! DEPENDS ON: ereport
        call ereport(routinename, errorstatus, cmessage)
      endif

!     Check namelist
      if (n_dumps == 0 .or. nhours == 0) then
        write (6,*) ' Error in setting DUMP2BOUND namelist'
        write (6,*) ' Both N_DUMPS and NHOURS must be set'
        write (6,*) ' N_DUMPS ',N_DUMPS,' NHOURS ',NHOURS
        go to 9999  !  Return
      endif

!     No of areas requiring boundary conditions.
!     Allow for multiple LAMs with number from dump2bound
      n_intf_a = no_lams

!     Derive data lengths.
!     U_FIELD is not known yet : Set to 1 for DERV_INTF_A
      U_FIELD = 1

! DEPENDS ON: derv_intf_a
      CALL DERV_INTF_A (TOT_LEN_INTFA_P,TOT_LEN_INTFA_U,                &
     &     MAX_INTF_MODEL_LEVELS,MAX_LBCROW_LENGTH,MAX_LBCROWS,         &
     &     N_INTF_A,U_FIELD,U_FIELD_INTFA)

!     Length of super arrays.
      len_a_spsts =1
      len_a_ixsts =1

!     No of data types for which boundary conditions required
!     Assume no tracer variables and set default number of
!     optional lbcs to zero
      tr_vars = 0
      num_optional_lbcs_out=0

! Add the number of optional lbc in this run depending on l_pc2
! and l_murk
      if(l_mcr_qcf2)then
        num_optional_lbcs_out=num_optional_lbcs_out+1
      endif
      if(l_mcr_qrain)then
        num_optional_lbcs_out=num_optional_lbcs_out+1
      endif
      if(l_mcr_qgraup)then
        num_optional_lbcs_out=num_optional_lbcs_out+1
      endif
      if(l_pc2)then
        num_optional_lbcs_out=num_optional_lbcs_out+3
      endif
      if(l_murk)then
        num_optional_lbcs_out=num_optional_lbcs_out+1
      endif

      intf_lookupsa=13+num_optional_lbcs_out+tr_vars

! DEPENDS ON: makebc
      call makebc(ppxRecs,n_dumps,nhours,                               &
     &            um_versn,sub_hr_int,no_lams)

 9999 continue

! DEPENDS ON: timer
      Call Timer ( ProgName, 2 )


      write (6,*) ' '
      write (6,*) ' ##################################'
      write (6,*) ' MAKEBC program completed normally.'
      write (6,*) ' ##################################'

      stop
      END PROGRAM MAIN_MAKEBC


!+ Subroutine MAKEBC : Creates a boundary dataset from model dumps
!
! Subroutine Interface :

!+ Subroutine LOOP_OVER_DUMPS : Loop over dumps to get boundary data
!
! Subroutine Interface :

!+ Subroutine GET_BC : Get boundary conditions from model dump
!
! Subroutine Interface :
#endif
