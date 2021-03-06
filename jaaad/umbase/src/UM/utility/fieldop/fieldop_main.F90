#if defined(FLDOP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine interface:
      SUBROUTINE fieldop_main(pp_len2_lookup,                           &
     &                        max_len2_lookup,                          &
     &                        len1_lookup,                              &
     &                        data_add1,                                &
     &                        pp_fixhd,                                 &
     &                        pp_fixhd2,                                &
     &                        len_fixhd,                                &
     &                        pp_unit2,                                 &
     &                        op,                                       &
     &                        iwa,                                      &
     &                        pp_unit1,                                 &
     &                        pp_unit_out,                              &
     &                        divisor,                                  &
     &                        pp_len2_lookup2,                          &
     &                        len1_lookup2,                             &
     &                        data_add2,                                &
     &                        iwa2,                                     &
     &                        nfields,                                  &
     &                        tfields,                                  &
     &                        llev,                                     &
     &                        Tcopy,                                    &
     &                        l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,   &
     &                        l13,l14,l15,l16,l17,l18,l19,l20,          &
     &                        stash1,stash2,stash3,stash4,stash5,       &
     &                        stash6,stash7,stash8,stash9,stash10,      &
     &                        stash11,stash12,stash13,stash14,stash15,  &
     &                        stash16,stash17,stash18,stash19,stash20,  &
     &                        ustash,                                   &
     &                        icode,cmessage)
      IMPLICIT NONE
!
! Description: Read in the lookup tables for files 1 and 2, find
!              lengths of current record, checking that the number of
!              values in field agree in both files. Obtain also, the
!              PPXREF codes for each field.
! Method:
!
! Current Code Owner: I Edmond
!
! History:
! Version   Date     Comment
! -------   ----     -------
! <version> <date>   Original code. <Your name>
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
! 1.0 Global variables (*CALLed COMDECKs etc...):
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "cstash.h"
#include "clookadd.h"
#include "c_mdi.h"

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER                                                           &
     & len1_lookup,                                                     &
                         ! 1st dimension of lookup
     & len1_lookup2,                                                    &
                         ! 1st dimension of lookup
     & pp_len2_lookup,                                                  &
                         ! 2nd dimension of lookup
     & pp_len2_lookup2,                                                 &
                         ! 2nd dimension of lookup
     & max_len2_lookup,                                                 &
                         ! 2nd dimension of lookup
     & len_fixhd,                                                       &
     & pp_unit1,                                                        &
                         ! unit no of required fieldsfile
     & pp_unit2,                                                        &
                         ! unit no of required fieldsfile
     & pp_unit_out,                                                     &
                         ! unit no of output file
     & data_add1,                                                       &
                         ! Start address of data in 1st file; Word
                         ! address of the data.
     & data_add2,                                                       &
                         ! Start address of data in 2nd file; Word
                         ! address of the data.
     & iwa,                                                             &
                         ! Start address of lookup in 1st file; Word
                         ! address in call setpos
     & iwa2,                                                            &
                         ! Start address of lookup in 2nd file; Word
                         ! address in call setpos
     & stash1,stash2,stash3,stash4,stash5,                              &
                                                 ! Stash codes of fields
     & stash6,stash7,stash8,stash9,stash10,                             &
                                                 ! not operated upon.
     & stash11,stash12,stash13,stash14,stash15,                         &
                                                 !
     & stash16,stash17,stash18,stash19,stash20,                         &
                                                 !
     & l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,                          &
     & l13,l14,l15,l16,l17,l18,l19,l20,                                 &
     & divisor           ! Integer divisor (file1 only)

      CHARACTER                                                         &
     & op*8

!   Array  arguments with intent(in):
      INTEGER                                                           &
     & pp_fixhd(len_fixhd),                                             &
                                         ! Dump/fieldsfile fixed header
     & pp_fixhd2(len_fixhd),                                            &
                                         ! Dump/fieldsfile fixed header
     & lookup(len1_lookup,pp_len2_lookup),                              &
                                             ! Integer lookup of file1.
     & lookup2(len1_lookup2,pp_len2_lookup2) ! Integer lookup of file2.

!   Array  arguments with intent(out):

!   ErrorStatus
      INTEGER                                                           &
     & icode

      CHARACTER                                                         &
     & cmessage*80

! Local parameters:
      INTEGER nft1,nft2
        PARAMETER(nft1=22, nft2=2)

      CHARACTER (Len=*), Parameter  :: RoutineName='FIELDOP'

! Local scalars:
       INTEGER                                                          &
     & i,j,k,                                                           &
                            ! local counters
     & nent,                                                            &
                            ! No of entries in the printfile
     & len1,                                                            &
                            ! Number of fields in File1
     & len2,                                                            &
                            ! Number of fields in File2
     & err,                                                             &
                            ! error code.
     & num_values,                                                      &
                            ! No. of points in data field
     & idim,                                                            &
                            ! num_values rounded to an even no
     & max_len,                                                         &
                            ! used to dimension the data array
     & len_i,                                                           &
                            ! No of data points in a fieldsfile field
                            ! used to find max_len.
     & entry_no,                                                        &
                            ! lookup entry no of the field.
     & entry_no2,                                                       &
                            ! lookup entry no of the field.
     & dummy,                                                           &
     & len_io,                                                          &
                            ! actual no of words transferred by IO.
     & len_io_expected,                                                 &
                            ! expected no of words transferred by IO
     & exppxi,                                                          &
     & rownumber                                                        &
     & ,ustash

      REAL                                                              &
     & a_io                 ! status returned by buffin / buffout

      CHARACTER                                                         &
     & exppxc*(36)

      LOGICAL                                                           &
     & model_flag,                                                      &
                            ! flag - set to true if model dump
     & lmore,                                                           &
     & l_copy,                                                          &
     & nfields,                                                         &
     & tfields,                                                         &
     & llev,                                                            &
     & ignore,                                                          &
     & Tcopy

! Local dynamic arrays:
      INTEGER                                                           &
     & pos1(max_len2_lookup),                                           &
                                ! Array of field positions in lookup1.
     & pos2(max_len2_lookup)    ! Equivalent field positions in lookup2.

! Function & Subroutine calls:
      External readff,read_write,writeff,ioerror,readstm

!- End of header

! Alter data and validity times in lookup tables so that lookup(1)
! -> lookup(14) taken from second file rather than first if
! logical Tcopy is TRUE.
        If (Tcopy) then

         ! Copy fixed header from pp_fixhd2 into pp_fixhd.
         Do j=21,41
          pp_fixhd(j) = pp_fixhd2(j)
         End do
! DEPENDS ON: setpos
         call setpos(pp_unit_out,0,icode)
! DEPENDS ON: buffout
         call buffout(pp_unit_out,pp_fixhd(1),len_fixhd,len_io,a_io)

         !  Check for I/O errors
         If (a_io  /=  -1.0 .or. len_io  /=  len_fixhd) then

! DEPENDS ON: ioerror
           call ioerror('buffer out of fixed header',a_io,len_io,       &
     &                   len_fixhd)
           cmessage='FIELDOP: I/O error'
           icode=25
           return
         End if
        End if

      ! Read in the lookup table of file1 if first time through
! DEPENDS ON: setpos
      call setpos(pp_unit1,iwa,icode)

      len_io_expected=pp_len2_lookup*len1_lookup
! DEPENDS ON: buffin
      call buffin(pp_unit1,lookup,len_io_expected,len_io,a_io)

      If (a_io /= -1.0 .or. len_io  /=  len_io_expected) then

! DEPENDS ON: ioerror
        call ioerror('Buffer in lookup table   ',a_io,len_io,           &
     &                len_io_expected )
        cmessage='fieldop_main: I/O error reading lookup table  '
        icode=3
        write(*,*)' I/O error reading lookup table'
        return

      End if

! Find which internal models are present and read in information from
! STASHmaster and user-STASHmaster files required by writflds.

      ! Find which internal models are present.
      internal_model_index(1) = 0
      internal_model_index(2) = 0
      internal_model_index(3) = 0
      internal_model_index(4) = 0
      n_internal_model = 1

      if (pp_fixhd(12) <  400)then

        do i =1, pp_len2_lookup

        ! Section 0: Prognostic fields.
          if(lookup(42,i) <= 100.or.                                    &
     &      (lookup(42,i) >= 200.and.lookup(42,i) <= 205))then
            lookup(45,i)=1

          else if((lookup(42,i) >  100.and.lookup(42,i) <= 176).or.     &
     &            (lookup(42,i) >= 180.and.lookup(42,i) <  200))then
            lookup(45,i)=2

          else if((lookup(42,i) >= 177.and.lookup(42,i) <= 179).or.     &
     &            (lookup(42,i) >= 210.and.lookup(42,i) <= 212))then
            lookup(45,i)=3

          ! Sections 1 - 99: Diagnostic fields
          else if(lookup(42,i) >= 1000.and.lookup(42,i) <= 29999)then
            if((lookup(42,i) >= 21177.and.lookup(42,i) <= 21179).or.    &
     &         (lookup(42,i) >= 21225.and.lookup(42,i) <= 21227).or.    &
     &         (lookup(42,i) >= 22177.and.lookup(42,i) <= 22179).or.    &
     &         (lookup(42,i) >= 22225.and.lookup(42,i) <= 22227).or.    &
     &         (lookup(42,i) >= 23177.and.lookup(42,i) <= 23179).or.    &
     &         (lookup(42,i) >= 23225.and.lookup(42,i) <= 23227).or.    &
     &         (lookup(42,i) >= 24177.and.lookup(42,i) <= 24179).or.    &
     &         (lookup(42,i) >= 24225.and.lookup(42,i) <= 24227))then
              lookup(45,i)=3        !Slab diagnostic

            else
              lookup(45,i)=1        !Atmosphere diagnostic

            end if

          else if(lookup(42,i) >= 30000.and.lookup(42,i) <= 99999)then
            if(lookup(42,i) >= 40000.and.lookup(42,i) <= 40999)then
              lookup(45,i)=3        !Slab diagnostic

            else
              lookup(45,i)=2        !Ocean diagnostic

            end if

          else
            write(6,*) 'WARNING: User defined field found - ',          &
     &                 'STASH code : ', lookup(42,i)
            write(6,*) ' Internal model number can not be defined.'
            write(6,*) ' Setting internal model number to atmosphere.'
            lookup(45,i)=1

          end if

        end do

      end if

      do i =1, pp_len2_lookup
        l_copy=.true.
        do j =1,n_internal_model
          if (lookup(45,i) == internal_model_index(j)) then
            l_copy =.false.
          end if
        end do
        if (l_copy) then
          internal_model_index(n_internal_model) = lookup(45,i)
          n_internal_model = n_internal_model +1
        end if
      end do
      n_internal_model = n_internal_model -1

      ! Read in STASHmaster file
      ppxRecs=1
      RowNumber=0
      do i=1,n_internal_model

        if(internal_model_index(i) == 1)then
! DEPENDS ON: hdppxrf
          call hdppxrf(nft1,'STASHmaster_A',ppxRecs,icode,cmessage)
        else if(internal_model_index(i) == 2)then
! DEPENDS ON: hdppxrf
          call hdppxrf(nft1,'STASHmaster_O',ppxRecs,icode,cmessage)
        else if(internal_model_index(i) == 3)then
! DEPENDS ON: hdppxrf
          call hdppxrf(nft1,'STASHmaster_S',ppxRecs,icode,cmessage)
        end if

        if(ICODE >  0)THEN
! DEPENDS ON: ereport
          CALL EReport(RoutineName,ICODE,CMessage)
        end if

        if(internal_model_index(i) == 1)then
! DEPENDS ON: getppx
          call getppx(nft1,nft2,'STASHmaster_A',RowNumber,              &
#include "argppx.h"
     &                icode,cmessage)
        else if(internal_model_index(I) == 2)then
! DEPENDS ON: getppx
          call getppx(nft1,nft2,'STASHmaster_O',RowNumber,              &
#include "argppx.h"
     &                icode,cmessage)
        else if(internal_model_index(I) == 3)then
! DEPENDS ON: getppx
          call getppx(nft1,nft2,'STASHmaster_S',RowNumber,              &
#include "argppx.h"
     &                icode,cmessage)
        end if

        if(icode /= 0)then
          write(6,*) cmessage
! DEPENDS ON: ereport
          CALL EREPORT('FIELDOP_MAIN', icode,                           &
     &     cmessage)

        end if

      end do


      !User STASHmaster
      if (ustash /= 0) then

! DEPENDS ON: hdppxrf
        call hdppxrf(0,' ',ppxRecs,icode,cmessage)

        if(icode /= 0)then
          write(6,*) cmessage
! DEPENDS ON: ereport
          CALL EREPORT('FIELDOP_MAIN', icode,                           &
     &     cmessage)

        end if

! DEPENDS ON: getppx
        call getppx(0,nft2,' ',RowNumber,                               &
#include "argppx.h"
     &              icode,cmessage)

        if(icode /= 0)then
          write(6,*) cmessage
! DEPENDS ON: ereport
          CALL EREPORT('FIELDOP_MAIN', icode,                           &
     &     cmessage)

        end if

      end if


      ! Read in the lookup table of file2 if first time through.
      If (op  /=  'idiv    ') then

! DEPENDS ON: setpos
        call setpos(pp_unit2,iwa2,icode)

        len_io_expected=pp_len2_lookup2*len1_lookup2
! DEPENDS ON: buffin
        call buffin(pp_unit2,lookup2,len_io_expected,len_io,a_io)

        If(a_io  /=  -1.0 .or. len_io  /=  len_io_expected) then

! DEPENDS ON: ioerror
          call ioerror('Buffer in lookup table2   ',a_io,len_io,        &
     &                 len_io_expected )
          cmessage='fieldop_main: I/O error reading lookup table  '
          icode=3
          write(*,*)' I/O error reading lookup table'
          return

        End if
      End if

      ! Calculate the number of fields in File1
      len1=0
      Do i =1,pp_len2_lookup

        pos1(i)=-1
        If (lookup(lbrow,i)  /=  -99) then

          len1 =len1 +1
        Else

          goto 2
        End if

      End do
    2 continue

      If (OP  /=  'idiv    ') then

        ! Calculate the number of fields in File2
        len2 =0
        Do i=1,pp_len2_lookup2

          pos2(i)=-1
          If (lookup2(lbrow,i)  /=  -99) then

            len2 =len2 +1
          Else

            goto 3
          End if

        End do ! i
    3   continue

        ! Find positions of corresponding fields in files 1 and 2;
        ! Store field positions in pos1 and equivalent file2 field
        ! positions in pos2.
        Do k =1,len1

          Do i =1,len2


            If ((lookup(42,k)  ==  lookup2(42,i)) .and.                 &
     &          (lookup(18,k)  ==  lookup2(18,i)) .and.                 &
     &          (lookup(19,k)  ==  lookup2(19,i)) .and.                 &
     &          (lookup(32,k)  ==  lookup2(32,i)) .and.                 &
     &          (lookup(33,k)  ==  lookup2(33,i))) then

               Do j=1,k
                If (i  ==  pos2(j)) goto 4
               End do
               pos1(k) =k
               pos2(k) =i
                goto 5

            End if

    4       continue
          End do
    5     continue
        End do

      Else
        Do k=1,len1
          pos1(k) =k
        End do

      End if

! Note lbrow=18,lbnpt=19
! For a DUMP lblrec will hold original no of data points.
! LBNREC will be set to zero.
!
! For a PP_file lblrec will hold the no of CRAY words needed to hold
! the data. The original field size will be rows*columns.
! If the data is not packed then lblrec=lbrow*lbnpt+lbext, where
! lbext will be greater than 0 for timeseries (which are never packed).
!  !! WARNING lbext - may be -32768 MISSING VALUE !!

      ! Set model_flag and reset UNPACK if DUMP
      If(pp_fixhd(5) /= 3) then

        model_flag=.true.       ! Model dump
        write(*,*)'Model dump - UNPACK set true '
      Else

        model_flag=.false.      ! Fieldsfile
      End if

      ! Find maximum field length to dimension array 'field' holding
      ! the data for each field - prevents writing outside bounds
      ! of array.
      max_len=0
      Do i =1,len1

        If (model_flag) then

          If (lookup(lblrec,i) >  max_len) then
            max_len =lookup(lblrec,i)
          End if

          ! For datafiles 1 and 2, check that number of values in
          ! field agrees.
        Else

          len_i =lookup(lbrow,i) *lookup(lbnpt,i)+lookup(lbext,i)
          If (len_i >  max_len) then
            max_len =len_i
          End if

        End if

      End do

      ! Set the length of data record (see above).
      ! Loop thro all the entries within the field
      Do i=1,len1

        ignore=.false.
        If (pos1(i)  ==  -1) ignore=.true.

        If (model_flag) then
          num_values =lookup(lblrec,i)

        Else
          num_values =lookup(lbrow,i) *lookup(lbnpt,i)                  &
     &                +lookup(lbext,i)
          End if

        If (lookup(lbext,i)  >   0) then

         ! got some extra data; check to see that we don't have
         ! packing if we have extra data....
          If(lookup(lbrow,i)*lookup(lbnpt,i)+                           &

     &      lookup(lbext,i)  /=  lookup(lblrec,i)) then
            cmessage='fieldop_main: Packing of extra data not supported'
            icode=1
            return

          End if
        End if

        idim =((num_values+1)/2) *2
        entry_no  =i
         entry_no2 =pos2(i)

! Alter data and validity times in lookup tables so that lookup(1)
! -> lookup(14) taken from second file rather than first if
! logical Tcopy is TRUE.
        If (Tcopy) then

         If (.not.model_flag) then
          ! Copy lookup tables from lookup2 to lookup.
          If (entry_no2  /=  -1) then
           Do j=1,14
            lookup(j,entry_no) = lookup2(j,entry_no2)
           End do
          Else
           cmessage='ERROR with -T option as fields dont match'
! DEPENDS ON: ereport
           CALL EREPORT('FIELDOP_MAIN', 1003,                           &
     &      cmessage)

          End if
         Else
          Do j=1,14
           lookup(j,entry_no) = lookup2(j,1)
          End do
         End if

! DEPENDS ON: setpos
         call setpos(pp_unit_out,pp_fixhd(150)+(i-1)                    &
     &               *len1_lookup-1,icode)
! DEPENDS ON: buffout
         call buffout(pp_unit_out,lookup(1,i),pp_fixhd(151)             &
     &                ,len_io,a_io)

      !  Check for I/O errors
         If (a_io  /=  -1.0 .or. len_io  /=  pp_fixhd(151)) then

! DEPENDS ON: ioerror
           call ioerror('buffer out of lookup table',a_io,len_io,       &
     &                  pp_fixhd(151))
           cmessage='FIELDOP: I/O error'
           icode=25
           return
         End if
        End if

! DEPENDS ON: read_write
        call read_write(num_values,                                     &
                                        ! IN Length of packed field
     &                  pp_unit1,                                       &
                                        ! IN Unit no of 1st I/P file.
     &                  pp_unit2,                                       &
                                        ! IN Unit no of 2nd I/P file.
     &                  len1_lookup,                                    &
                                        ! IN 1st dim of lookup of file1.
     &                  pp_len2_lookup,                                 &
                                        ! IN 2nd dim of lookup of file1.
     &                  len_fixhd,                                      &
     &                  pp_fixhd,                                       &
                                        ! IN Fixed header of file1.
     &                  lookup,                                         &
                                        ! IN Lookup table file1.
                                        !    (integer part used).
     &                  pp_fixhd2,                                      &
                                        ! IN Fixed header of file2.
     &                  lookup2,                                        &
                                        ! IN Lookup table file2.
                                        !    (Integer part used).
     &                  lookup2,                                        &
                                        ! IN (Real part used).
     &                  len1_lookup2,                                   &
                                        ! IN 1st dim of lookup of file2.
     &                  pp_len2_lookup2,                                &
                                        ! IN 2nd dim of lookup of file1.
     &                  op,                                             &
                                        ! IN Operation type.
     &                  lookup,                                         &
                                        ! IN Lookup table file1.
                                        !    (real part used).
     &                  entry_no,                                       &
                                        ! IN Posn of field in lookup1.
     &                  entry_no2,                                      &
                                        ! IN Posn of field in lookup2.
     &                  data_add2,                                      &
                                     ! IN Start address of data file2.
     &                  data_add1,                                      &
                                     ! IN Start address of data file1.
     &                  model_flag,                                     &
                                     ! IN TRUE (dump).FALSE (fieldsfile)
     &                  nfields,                                        &
     &                  tfields,                                        &
     &                  llev,                                           &
     &                  ignore,                                         &
     &                  pp_unit_out,                                    &
                                     ! IN Unit no. of O/P file.
     &                  max_len,                                        &
     &                  divisor,                                        &
                                     ! IN Integer divisor if specified.
     &                  lookup(63,entry_no),                            &
     &                  l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,         &
     &                  l13,l14,l15,l16,l17,l18,l19,l20,                &
     &                  stash1,stash2,stash3,stash4,stash5,             &
     &                  stash6,stash7,stash8,stash9,stash10,            &
     &                  stash11,stash12,stash13,stash14,stash15,        &
     &                  stash16,stash17,stash18,stash19,stash20,        &
#include "argppx.h"
     &                  icode,cmessage)     ! Error code/message.

        If(icode /= 0) then
          return
        End if

      End do ! i

      RETURN
      END SUBROUTINE fieldop_main

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
