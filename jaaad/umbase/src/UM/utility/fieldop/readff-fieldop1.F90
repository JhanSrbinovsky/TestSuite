#if defined(FLDOP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine interface:

!
! Subroutine interface:

!
! Subroutine interface:

!
! Subroutine interface:
      subroutine readff(pp_unit1,                                       &
     &                  field,                                          &
     &                  idim,                                           &
     &                  entry_no,                                       &
     &                  ilabel,                                         &
     &                  rlabel,                                         &
     &                  pp_len2_lookup,                                 &
     &                  len1_lookup,                                    &
     &                  len_fixhd,                                      &
     &                  pp_fixhd,                                       &
     &                  lookup,                                         &
     &                  rookup,                                         &
     &                  data_add1,                                      &
     &                  model_flag,                                     &
     &                  max_len_ilabel,                                 &
     &                  max_len_rlabel,                                 &
     &                  max_len,                                        &
     &                  pppak,                                          &
     &                  len_ilabel,                                     &
     &                  len_rlabel,                                     &
     &                  iwa,                                            &
     &                  icode,cmessage)
      IMPLICIT NONE
!
!
! Description: To read a direct access PP file.
!
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
#include "clookadd.h"
#include "c_mdi.h"

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER                                                           &
     & len1_lookup,                                                     &
                              ! first dimension of the lookup
     & pp_len2_lookup,                                                  &
                              ! secnd dimension of the lookup
     & pp_unit1,                                                        &
                              ! unit no of required fieldsfile
     & idim,                                                            &
                              ! Size of data field (rounded)
     & max_len_rlabel,                                                  &
                              ! max sixe of rlabel
     & max_len_ilabel,                                                  &
                              ! max sixe of ilabel
     & max_len,                                                         &
     & data_add1,                                                       &
                              ! The word address of the data.
     & entry_no,                                                        &
                              ! Lookup entry no of the Field.
     & len_fixhd,                                                       &
     & lookup(len1_lookup,pp_len2_lookup) ! integer lookup

      REAL                                                              &
     & rookup(len1_lookup,pp_len2_lookup) ! real lookup

      LOGICAL                                                           &
     & model_flag             ! True => Dump False =>Fieldsfile

!   Array  arguments with intent(in):
      INTEGER                                                           &
     & pp_fixhd(len_fixhd)    ! fixed header

!   Scalar arguments with intent(out):
      INTEGER                                                           &
     & len_rlabel,                                                      &
                              ! actual size of rlabel
     & len_ilabel,                                                      &
                              ! actual size of ilabel
     & pppak

!   Array  arguments with intent(out):
      INTEGER                                                           &
     & ilabel(max_len_ilabel) ! integer part of lookup

      REAL                                                              &
     & field(idim),                                                     &
                           ! array holding final output data.
     & rlabel(max_len_rlabel) ! real part of lookup

!   ErrorStatus
      INTEGER                                                           &
     & icode                  ! error code

      CHARACTER                                                         &
     & cmessage*80          ! error message

! Local scalars:
       INTEGER                                                          &
     & i,j,                                                             &
                              ! Local counters
     & pack_type,                                                       &
                              ! packing type N1 of LBPACK
     & num_cray_words,                                                  &
                              ! number of words for field
     & nvals,                                                           &
                              ! number of points in a data field
     & iwa,                                                             &
                              ! Word address in call setpos
     & length_of_data,                                                  &
                              ! Length of a particular field
     & addr,                                                            &
                              ! Address of a field in the data store
     & pos_rlabel,                                                      &
                              ! position of first REAL in PPhdr
     & pack_type_i            ! packing type N1 of LBPACK

      REAL                                                              &
     & amdi                   ! Missing data indicator for lookup

! Function & Subroutine calls:
      External setpos,read_rec,ioerror,coex,integer_to_real

!- End of header

      amdi=rookup(bmdi,entry_no)
      If (amdi /= rmdi) write(*,*)' NON STANDARD MISSING DATA USED'

      pack_type = MOD(lookup(lbpack,entry_no),10)

      ! Reading a model type dump
      ! A model dump has no direct addressing only relative.

      If(model_flag) then

! Old Format dumpfiles
        if((lookup(lbnrec,entry_no) == 0) .or.                          &
! Prog lookups in dump before vn3.2:
     &    ((lookup(lbnrec,entry_no) == imdi) .and.                      &
     &                             (pp_fixhd(12) <= 301))) then

        If(pack_type == 2) then            ! 32 bit packing.

          num_cray_words = (lookup(lblrec,entry_no)+1)/2
        Else if (pack_type >  0) then

          num_cray_words = lookup(lblrec,entry_no)/2
        Else

          num_cray_words = lookup(lblrec,entry_no)
        End if

        nvals = lookup(lblrec,entry_no) ! No of data points
        addr=data_add1

        If (entry_no >  1) then

          Do i =1,entry_no-1

            pack_type_i = MOD(lookup(LBPACK,I),10)
            If (pack_type_i  ==  2) then ! 32 Bit packed

              length_of_data = (lookup(lblrec,I)+1)/2
            Else

              length_of_data = lookup(lblrec,I)
            End if

            addr = addr + length_of_data

          End do ! i
        Else       !  If the first entry.

          addr = data_add1
          If (pack_type  ==  2) then ! 32 Bit packed

            length_of_data = (lookup(lblrec,1)+1)/2
          Else

            length_of_data=lookup(lblrec,1)
          End if

          write(*,*)'  length_of_data  ',length_of_data

        End if

        iwa=addr  ! Not -1 as this is already done in dump

      Else
! New format Dumpfiles (vn4.4 onwards)

        If(pack_type == 2) then            ! 32 bit packing.
          num_cray_words=(lookup(lblrec,entry_no)+1)/2
        Elseif(pack_type >  0) then
          num_cray_words=lookup(lblrec,entry_no)/2
        Else
          num_cray_words=lookup(lblrec,entry_no)
        Endif
        iwa = lookup(LBEGIN,entry_no)
        nvals = lookup(lbrow,entry_no) * lookup(lbnpt,entry_no)
      Endif
      Else ! Reading a PP type file.

        num_cray_words = lookup(lblrec,entry_no) ! PP type file
        iwa = lookup(LBEGIN,entry_no)
        nvals = lookup(lbrow,entry_no) * lookup(lbnpt,entry_no)         &
     &                                 + lookup(lbext,entry_no)

      End if

  107 FORMAT(' ENTRY NO=',I5,'num_cray_words= ',I6,'nvals=',I6)

        If (idim  <   num_cray_words) then

          icode = num_cray_words
          cmessage ='readff  Idim to small icode holds correct value'
          goto 9999

        End if

      icode=0
! DEPENDS ON: read_rec
      call read_rec(field,                                              &
                                     ! OUT array holding data
     &              num_cray_words,                                     &
                                     ! IN No of CRAY words holding data
     &              iwa,                                                &
                                  ! IN WORD address of field to be read
     &              pp_unit1,                                           &
                                  ! IN unit no of the file
     &              max_len,                                            &
     &              icode,                                              &
                                  ! IN/OUT
     &              pack_type)

 2212 FORMAT('  FIELDS FILE NUMBER ',I2,'  ON UNIT',I2,2X,'BEING read')

        If (icode == 0) then

          pos_rlabel = MOD(lookup(lbrel,entry_no),100)

          ! Treat lookup(45) (submodel identifier) as an integer.
          POS_RLABEL=46


          len_rlabel=1+len1_lookup-pos_rlabel
          len_ilabel=len1_lookup-len_rlabel

          Do i=1,len_ilabel
            ilabel(i)=lookup(i,entry_no)
          End do

!         check for valid release number
          if (ilabel(lbrel) <  1) then

            write(*,*)' resetting LBREL from',ilabel(lbrel),' to 2'
            ilabel(lbrel)=2

          endif

          Do i=1,len_rlabel
            rlabel(i)=rookup(i+pos_rlabel-1,entry_no)
          End do

        End if

        ! At this point field holds the data either packed or un-packed
        ! Is the packing indicator set and is un-packing required?
        ! If so then the data is temp un-packed into a work array of
        ! length idim
        If (pack_type >  0) then       ! Is the field packed.

! DEPENDS ON: un_pack
            call un_pack(pack_type,                                     &
                                       ! IN packing type N1 of LBPACK
     &                   idim,                                          &
                                       ! IN length of unpacked pp buffer
     &                   field,                                         &
                                       ! IN/OUT I/P contains packed data
                                       ! Output contains un-packed data.
     &                   num_cray_words,                                &
                                         ! IN length of input field
     &                   ilabel,                                        &
                                       ! IN holds integer part of lookup
     &                   len_ilabel,                                    &
                                       ! IN length of ilabel array
     &                   amdi,                                          &
                                       ! IN Missing data indicator.
     &                   pp_fixhd,                                      &
                                       ! IN PPfile fixed length header
     &                   len_fixhd,                                     &
     &                   pppak,                                         &
     &                   icode,cmessage)  ! IN/OUT

        Else if(lookup(data_type,entry_no) == 2) then !Fld is integer

! DEPENDS ON: integer_to_real
          call integer_to_real(idim,                                    &
                                       ! IN full unpacked size of field
     &                         field,                                   &
                                       ! IN contains integer data.
     &                         field,                                   &
                                       ! OUT contains Real data.
     &                         nvals,                                   &
                                       ! IN no of values in field
     &                         max_len,                                 &
     &                         ilabel,                                  &
                                       ! IN/OUT integer part of lookup
     &                         icode)  ! IN/OUT error code

        End if

 9999 continue
  100 FORMAT(//,32X,'   ARRAY        ',//,32(16F5.0/))
  101 FORMAT(//,32X,'   lookup       ',//,32(16I5/))
  103 FORMAT('   LENIN  ',I12)

      RETURN
      END SUBROUTINE readff

!
! Subroutine interface:
!
! Subroutine interface:
!
!Subroutine interface:
!
!Subroutine interface:
#endif
