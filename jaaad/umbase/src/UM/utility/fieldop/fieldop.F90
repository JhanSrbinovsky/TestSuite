#if defined(FLDOP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
       PROGRAM FIELDOP
       IMPLICIT NONE
!
! Routine: fieldop -------------------------------------------------
!
! Description:
! To read two model dumps or direct access fieldsfiles with unpacked
! or packed (wgdos,grib,cray 32 bits) data and write out to a new file
! the difference, sum or product of the data values. Alternatively
! if a single dataset is read the data may be divided by an integer.
!
! Method:
!
! Current Code Owner: I Edmond
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 4.3       20/2/97  Code added to enable fieldop (1) to perform
!  arithmetic operations on individual levels of specified fields
!  (indicated by lblev and lbuser3), (2) specify which data/
!  validity times (file1 or file2) are written out to output file.
!  (See fieldop script for details).
!  As the dimensions of a wgdos packed field depends on the size of
!  the data values; having operated on fieldsfiles it is
!  necessary to write data out to new addresses. Corrections were
!  made to code to do this.   Ian Edmond.
!    4.4 17/7/97 Fix to subroutine READFF to read wfio dumpfiles.
!        Initialise icode =0 now that it is not done in hdppxrf  IE
!LL  4.4   Oct. 1997 Changed error handling from routine HDPPXRF
!LL                  so only fatal (+ve) errors are handled.
!LL                                             Shaun de Witt
!    4.5  14/07/98  Initialised UM_SECTOR_SIZE from BLKDATA
!                   (A Van der Wal)
!    4.5  23/11/98  Use UM_SECTOR_SIZE instead of hardwired 512 in
!                   pp_file. Automatically skip mathematical
!                   operations on land-sea mask, to prevent failure
!                   with fieldsfiles.
!                   Author D.M. Goddard
!LL  5.2  7/11/99  Enable run length encoding of ocean fieldsfiles to
!LL                compress the sequences of mdi values that represent
!LL                the land points. Ian Edmond
!LL  5.0  06/05/99  Pass RoutineName to EREPORT. M Gallani
!    5.3  14/02/02  Enable MPP as the only option for
!                   small executables         E.Leung
!    5.4  30/08/02  Pass dummy variables to PP_FILE in place
!                   of fields used when processing variable horizontal
!                   grid "extra" data. The "extra data" vectors
!                   relating to grid spacings are not set in dump
!                   fields at present, so for now, we ensure that
!                   existing functionality is not damaged. PP files
!                   created via this process will thus not contain the
!                   "extra" grid data which pp files generated via a
!                   main UM run would contain.                R. Hill
!    5.4  03/09/02  Arguments are added to DECOMPOSE_SMEXE
!                                                  E.Leung
!    5.5  25/04/03  Grib data format not supported on non-CRAY
!                   platform                         E.Leung
!    6.0  11/09/03  Replaced call to ABORT with EREPORT call.  P.Dando
!    6.0  02/02/04  Surround call to degrib with NECSX6 DEF as degrib
!                   only currently available on the SX6. W Roseblade.
!    6.1  18/08/04  Re-name PP_FILE subroutine as FLDOP_PP_FILE to avoid
!                   confusion with routine provided in deck PPFILE1A.
!                   Remove dummy variables passed in call to PP_FILE.
!                   Re-enable Grib data on non-CRAY platforms.   P.Dando
!    6.2  10/04/05 Removed calls to ABORT.  J. Gill
!    6.2  04/02/05  Fixed Cray 32 bit unpacking routine by replacing
!                   with a call to buffin32 in this case. T. Edwards
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: <appropriate code>
! System Task:              <appropriate code>
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Routine arguments
!   Scalar arguments
      INTEGER                                                           &
     & i,                                                               &
                         ! Counter.
     & len2_lookup,                                                     &
                      ! Size of the lookup on the file
     & len2_lookup2,                                                    &
                      ! Size of the lookup on the file
     & max_len2_lookup,                                                 &
                       ! Size of the lookup on the file
     & LEN_INTHD,                                                       &
     & LEN_REALHD,                                                      &
     & LEN1_LEVDPC,                                                     &
     & LEN2_LEVDPC,                                                     &
     & pp_unit_out,                                                     &
                         ! Unit no of output file; value varies
                         ! - depends on 1 or 2 i/p files.
     & icode,                                                           &
                         ! Return code
     & data_add1,                                                       &
                         ! The word address of the data.
     & data_add2,                                                       &
                         ! The word address of the data.
     & iwa,                                                             &
                         ! Word address in call setpos
     & iwa2,                                                            &
                         ! Word address in call setpos
     & len_io,                                                          &
                         ! Length of IO done
     & l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,                          &
     & l13,l14,l15,l16,l17,l18,l19,l20,                                 &
     & stash1,stash2,stash3,stash4,stash5,                              &
                                                ! Stash codes of fields
     & stash6,stash7,stash8,stash9,stash10,                             &
                                                ! which are not operated
     & stash11,stash12,stash13,stash14,stash15,                         &
                                                ! upon.
     & stash16,stash17,stash18,stash19,stash20,                         &
                                                !
     & divisor,                                                         &
                        ! Integer divisor for data in file 1 if required
     & err,                                                             &
                        ! Error code.
     & OpenStatus                                                       &
     &,ustash

      REAL                                                              &
     &     a_io            ! status returned by buffin

      CHARACTER                                                         &
     &     cmessage*80                                                  &
                          ! Error message from lower routines
     &    ,op*8            ! Operation type +,-,*

      CHARACTER NOMLIST*80
      CHARACTER*(*) RoutineName
      PARAMETER (RoutineName='FIELDOP')

      LOGICAL                                                           &
     & nfields                                                          &
     &,tfields                                                          &
     &,llev                                                             &
     &,Tcopy

! Parameters:
      INTEGER len_fixhd       ! Length of fixed length header
        PARAMETER(len_fixhd=256)

      INTEGER len1_lookup     ! First dim. of the lookup of 1st dump
        PARAMETER(len1_lookup=64)

      INTEGER len1_lookup2    ! First dim. of the lookup of 2nd dump
        PARAMETER(len1_lookup2=64)

      INTEGER pp_unit1        ! Unit number of input dump/fieldsfile.
        PARAMETER(pp_unit1=20)

      INTEGER pp_unit2        ! Unit number of 2nd i/p dump/fieldsfile.
        PARAMETER(pp_unit2=21)

! Array arguments:
      INTEGER                                                           &
     & pp_fixhd(len_fixhd),                                             &
                                !  Fixed length header of 1st file.
     & pp_fixhd2(len_fixhd)     !  Fixed length header of 2nd file.

      CHARACTER*100 PAREXE_ENV  ! hold name of the // exec script
      INTEGER ME_GC,NPROC_GC
! Function & Subroutine calls:
      External readff,setpos,ioerror,fieldop_main,gc_init

#include "cntl_io.h"

      DATA l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,                      &
     &     l13,l14,l15,l16,l17,l18,l19,l20  /                           &
     &     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                                &
     &     0, 0, 0, 0, 0, 0, 0, 0, 0, 0  /

      DATA stash1,stash2,stash3,stash4,stash5,                          &
     &     stash6,stash7,stash8,stash9,stash10,                         &
     &     stash11,stash12,stash13,stash14,stash15,                     &
     &     stash16,stash17,stash18,stash19,stash20  /                   &
     &     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                                &
     &     0, 0, 0, 0, 0, 0, 0, 0, 0, 0  /

!- End of header

!------------------------------------------------------------------
! Read in the Fixed Length Headers of files 1 and 2.
!------------------------------------------------------------------

      namelist/CONTROL/op,divisor,nfields,tfields,llev,Tcopy
      namelist/STASHES/stash1,stash2,stash3,stash4,stash5,              &
     &                 stash6,stash7,stash8,stash9,stash10,             &
     &                 stash11,stash12,stash13,stash14,stash15,         &
     &                 stash16,stash17,stash18,stash19,stash20
      namelist/LEVELS/l1,l2,l3,l4,l5,                                   &
     &                 l6,l7,l8,l9,l10,                                 &
     &                 l11,l12,l13,l14,l15,                             &
     &                 l16,l17,l18,l19,l20
      namelist/USTSFILE/ustash

      ! Initialise error code.
      icode = 0


      PAREXE_ENV=' '
      CALL GC_INIT(PAREXE_ENV,ME_GC,NPROC_GC)

      Call GET_FILE(5,NOMLIST,80,ICODE)
      OPEN(UNIT=5,FILE=NOMLIST)
      read(5,CONTROL)
      read(5,STASHES)
      read(5,LEVELS)
      read(5,USTSFILE)

      len2_lookup  = 0
      len2_lookup2 = 0
      ! Open 1st dump or fieldsfile.
! DEPENDS ON: file_open
      call file_open(pp_unit1,'FILE1',5,0,0,err)

      If (op  /=  'idiv    ') then

       ! Open 2nd dump or fieldsfile and output file.
! DEPENDS ON: file_open
        call file_open(pp_unit2,'FILE2',5,0,0,err)

        pp_unit_out =22
! DEPENDS ON: file_open
        call file_open(pp_unit_out,'FILE3',5,1,0,err)
      Else

       ! Only one dump or fieldsfile needed, open output file.
        pp_unit_out=21
! DEPENDS ON: file_open
        call file_open(pp_unit_out,'FILE2',5,1,0,err)

      End if

      ! Read fixed header of first file.
! DEPENDS ON: buffin
      call buffin(pp_unit1,pp_fixhd,len_fixhd,len_io,a_io)

      ! Error check.
      If (a_io  /=  -1.0 .or. len_io  /=  len_fixhd) then

! DEPENDS ON: ioerror
        call ioerror('Buffer in fixed length header',a_io,len_io,       &
     &                len_fixhd)
        cmessage ='FIELDOP : I/O error reading Fixed Length Header'
        icode =2
        write(*,*)' I/O error reading Fixed Length Header'
! DEPENDS ON: abort_io
        call abort_io('FIELDOP1',cmessage,icode,pp_unit1)

      End if

      data_add1      = pp_fixhd(160) -1 ! Start address for the data.
      iwa            = pp_fixhd(150) -1 ! Start address of lookup table.
      len2_lookup = pp_fixhd(152)    ! 2nd dim of lookup of file1.

      write(*,*)' dump type=',pp_fixhd(5),                              &
     &       ' 3=fieldsfile,1=dump,2=time mean dump,4=ancil,5=bound'

      If (op  /=  'idiv    ') then

       ! Read fixed header of second file.
! DEPENDS ON: buffin
        call buffin(pp_unit2,pp_fixhd2,len_fixhd,len_io,a_io)

       ! Error check.
        If(a_io  /=  -1.0 .or. len_io  /=  len_fixhd) then

! DEPENDS ON: ioerror
          call ioerror('Buffer in fixed length header2',a_io,len_io,    &
     &                  len_fixhd)
          cmessage='FIELDOP : I/O error reading Fixed Length Header'
          icode=2
          write(*,*)' I/O error reading Fixed Length Header'
! DEPENDS ON: abort_io
        call abort_io('FIELDOP1',cmessage,icode,pp_unit2)

        End if

        data_add2       = pp_fixhd2(160)-1 ! Start address for the data.
        iwa2            = pp_fixhd2(150)-1 ! Start address of lookup.
        len2_lookup2 = pp_fixhd2(152)   ! 2nd dim of lookup of file2.

        ! Compare fixed length headers
        write(6,*)' '
        write(6,*)'Fixed Length Header:'

        Do i =1,len_fixhd

          If (pp_fixhd(i)  /=  pp_fixhd2(i)) then

            write(6,*)'Item=',I,pp_fixhd(I),pp_fixhd2(I)
            ! Abort if files 1 and 2 have different indicators for
            ! dataset type.
            If (i == 5) then
              CMessage = 'ERROR: Different dataset types'
              ICODE = 1
! DEPENDS ON: ereport
              call EReport(RoutineName,ICODE,CMessage)
            End if

          End if

        End do ! i

      End if

      max_len2_lookup=MAX(len2_lookup,len2_lookup2)

! DEPENDS ON: fieldop_main
      call fieldop_main(len2_lookup,                                    &
                                     !IN 2nd dim of lookup of file1.
     &                  max_len2_lookup,                                &
                                        !IN 1st dim of lookup (file1).
     &                  len1_lookup,                                    &
                                        !IN 1st dim of lookup (file1).
     &                  data_add1,                                      &
                                        !IN Start address for data
                                        !   in 1st file.
     &                  pp_fixhd,                                       &
                                        !IN Fixed header of 1st file.
     &                  pp_fixhd2,                                      &
                                        !IN Fixed header of 2nd file.
     &                  len_fixhd,                                      &
                                        !IN Fixed header length
     &                  pp_unit2,                                       &
                                        !IN Unit no. of 2nd i/p dataset.
     &                  op,                                             &
                                        !IN Operation type +,-,* (char)
     &                  iwa,                                            &
                                        !IN Start address for the lookup
                                        !   table of 1st file.
     &                  pp_unit1,                                       &
                                        !IN Unit no. of 1st i/p dataset.
     &                  pp_unit_out,                                    &
                                        !IN Unit number of o/p file.
     &                  divisor,                                        &
                                        !IN Integer divisor for data
                                        !   in file 1 if required.
     &                  len2_lookup2,                                   &
                                     !IN 2nd dim of lookup of file2.
     &                  len1_lookup2,                                   &
                                        !IN 1st dim of lookup (file2).
     &                  data_add2,                                      &
                                        !IN Start address for data
                                        !   in 2nd file.
     &                  iwa2,                                           &
                                        !IN Start address for the lookup
                                        !   table of 2nd file.
     &                  nfields,                                        &
     &                  tfields,                                        &
     &                  llev,                                           &
     &                  Tcopy,                                          &
     &                  l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,         &
     &                  l13,l14,l15,l16,l17,l18,l19,l20,                &
     &                  stash1,stash2,stash3,stash4,stash5,             &
     &                  stash6,stash7,stash8,stash9,stash10,            &
     &                  stash11,stash12,stash13,stash14,stash15,        &
     &                  stash16,stash17,stash18,stash19,stash20,        &
     &                  ustash,                                         &
                                        !IN =1 if user STASHmaster file
     &                  icode,cmessage) !   =0 if no user STASHmaster

      If (icode /= 0) then

! DEPENDS ON: ereport
        call ereport(RoutineName,icode,cmessage)

      End if

      STOP
      END PROGRAM FIELDOP
!
! Subroutine interface:

!
! Subroutine interface:

!
! Subroutine interface:

!
! Subroutine interface:

!
! Subroutine interface:
!
! Subroutine interface:
!
!Subroutine interface:
!
!Subroutine interface:
#endif
