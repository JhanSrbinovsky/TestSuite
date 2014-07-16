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
      subroutine writeff(pp_unit_out,                                   &
     &               field,                                             &
     &               idim,                                              &
     &               entry_no,                                          &
     &               data_add,                                          &
     &               lookup,                                            &
     &               len_fixhd,                                         &
     &               fixhd,                                             &
     &               len2_lookup,                                       &
     &               len1_lookup,                                       &
     &               n_rows_out,                                        &
     &               n_cols_out,                                        &
     &               packed,                                            &
     &               max_len,                                           &
     &               comp_accry,                                        &
     &               op,                                                &
#include "argppx.h"
     &               icode,cmessage)

      IMPLICIT NONE
!
!
! Description: To ouput a field to a UM dump or fieldsfile, with the
!              data written in packed (wgdos,grib or cray 32 bits) or
!              unpacked form.
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
#include "clookadd.h"
#include "c_mdi.h"
#include "parparm.h"
#include "amaxsize.h"
#include "comvgrid.h"

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER                                                           &
     & n_rows_out,                                                      &
     & n_cols_out,                                                      &
     & len_fixhd,                                                       &
     & pp_unit_out,                                                     &
     & len2_lookup,                                                     &
     & len1_lookup,                                                     &
     & entry_no,                                                        &
     & grib_packing,                                                    &
     & max_len,                                                         &
     & comp_accry,                                                      &
     & idim,                                                            &
     & data_add,                                                        &
     & exppxi

      CHARACTER                                                         &
     & op*(8)

      CHARACTER                                                         &
     & exppxc*(36)

      LOGICAL                                                           &
     & packed

!   Array  arguments with intent(in):
      INTEGER                                                           &
     & lookup(len1_lookup,len2_lookup),                                 &
     & fixhd(len_fixhd),                                                &
     & ifield(max_len)

      REAL                                                              &
     & field(max_len)

!   ErrorStatus
      INTEGER                                                           &
     & icode

      REAL                                                              &
     & A

      CHARACTER                                                         &
     & cmessage*80

! Local scalars:
      INTEGER                                                           &
     & comp_accrcy,                                                     &
                          !
     & num_words,                                                       &
                          !
     & pack_type,                                                       &
                          ! Packing type N1 of LBPACK
     & len_io,                                                          &
                          !
     & isize,                                                           &
     & tot_levels,                                                      &
                          ! total number of levels
     & i                  !

! Function & Subroutine calls:
      External writflds

!- End of header

#include "cntl_io.h"
      icode       = 0
      num_words   = -99
      pack_type   = MOD(lookup(lbpack,entry_no),10)
      If (pack_type >  0) packed =.true.

      ! Method of GRIB packing - use width method, with simple packing
      ! to be similar to the ECMWF MARS archive.
      grib_packing=6

      If((lookup(44,entry_no) <  0) .or.                                &
     &                     (lookup(44,entry_no) >  100)) then

        lookup(44,entry_no) = 0
        lookup(44,entry_no) = lookup(44,entry_no) + 1
      Else

        lookup(44,entry_no) = lookup(44,entry_no) + 1
      End if

! DEPENDS ON: setpos
      call setpos(pp_unit_out,fixhd(150)+(entry_no-1)                   &
     &            *len1_lookup-1,icode)
! DEPENDS ON: buffout
      call buffout(pp_unit_out,lookup(1,entry_no),fixhd(151)            &
     &             ,len_io,A)

      ! Check for I/O errors
      If (A  /=  -1.0 .or. len_io  /=  fixhd(151)) then

! DEPENDS ON: ioerror
        call ioerror('buffer out of lookup table',A,len_io,             &
     &               fixhd(151))
        cmessage='FIELDOP: I/O error'
        icode=25
        return
      End if

      If (pack_type  == 1 .or.pack_type  ==  4 )then
        isize=(((((idim+um_sector_size-1)/um_sector_size)*              &
     &    um_sector_size) +1)/2) * 2
      ! Deactivate varible grid extra data processing.
      VAR_GRID_TYPE = 0

       ! Data packed using WGDOS method and written to O/P file.
! DEPENDS ON: fldop_pp_file
        call fldop_pp_file(field,                                       &
                                    ! IN Array to store expanded data
     &               isize,                                             &
                                    ! IN length of pp buffer
     &               num_words,                                         &
                                    ! IN No of 64bit words of data
     &               rmdi,                                              &
                                    ! IN Missing data
     &               comp_accry,                                        &
                                   ! IN PPXREF accuracy code.
     &               idim,                                              &
                                   ! IN length of pp buffer
     &               pp_unit_out,                                       &
                                    ! IN Unit no of O/P field.
     &               data_add,                                          &
                                   ! IN Word address of data (file1)
     &               n_cols_out,                                        &
                                    ! IN
     &               n_rows_out,                                        &
                                   ! IN
     &               packed,                                            &
                                    ! IN TRUE - packing required.
     &               pack_type,                                         &
                                    ! IN WGDOS packed data.
     &               lookup,                                            &
                                   ! IN lookup headers of file1.
     &               len1_lookup,                                       &
                                   ! IN
     &               len2_lookup,                                       &
                                   ! IN
     &               entry_no,                                          &
                                    ! IN
     &               icode,cmessage)

        If (icode >  0) then
          cmessage='FIELDOP : Error in FLDOP_PP_FILE'
! DEPENDS ON: ereport
          CALL EREPORT('WRITEFF', icode,                                &
     &     cmessage)

        End if

      Else if (pack_type == 3) then

       ! Data compressed using the GRIB method, written to O/P file.
! DEPENDS ON: grib_file
        call grib_file(len1_lookup,                                     &
                                       ! IN
     &                 len2_lookup,                                     &
                                       ! IN
     &                 lookup,                                          &
                                       ! IN
     &                 lookup,                                          &
                                       ! IN
     &                 entry_no,                                        &
                                       ! IN Posn of field in lookup.
     &                 field,                                           &
                                       ! IN Unpacked output data.
     &                 max_len,                                         &
                                       ! IN Length of pp buffer
     &                 max_len,                                         &
                                       ! IN
     &                 num_words,                                       &
                                       ! IN No of 64bit words of data
     &                 pp_unit_out,                                     &
                                       ! IN Unit no of O/P field.
     &                 pp_unit_out,                                     &
                                       ! IN Word address of record.
     &                 grib_packing,                                    &
                                       !
     &                 icode,cmessage) ! IN

        If (icode >  0) then
          cmessage='FIELDOP : Error in GRIB_FILE'
! DEPENDS ON: ereport
          CALL EREPORT('WRITEFF', icode,                                &
     &     cmessage)

        End if

      Else if ((pack_type == 0).or.(pack_type == 2)) then
        ! Update lookup header data lengths and addressing for
        ! unpacked data in fieldsfile.

        If ((fixhd(5) == 3).and.(pack_type == 0)) then
          If (entry_no  ==  1) then
            lookup(29,entry_no) = data_add
          Else
            lookup(29,entry_no) = lookup(29,entry_no-1)                 &
     &                            + lookup(30,entry_no-1)
          End If
          lookup(40,entry_no) = lookup(29,entry_no)
        End If

! DEPENDS ON: decompose_smexe
        call decompose_smexe(n_cols_out, n_rows_out,0,0,1)
        If (lookup(data_type,entry_no)  ==  2) then

        do i=1,idim
          ifield(i)=field(i)
        enddo

! DEPENDS ON: writflds
        call writflds(pp_unit_out,1,entry_no,lookup,len1_lookup,        &
     &                  ifield,lookup(lblrec,entry_no),fixhd,           &
#include "argppx.h"
     &                icode,cmessage)
        Else

       ! Data unpacked.
! DEPENDS ON: writflds
        call writflds(pp_unit_out,1,entry_no,lookup,len1_lookup,        &
     &                  field,lookup(lblrec,entry_no),fixhd,            &
#include "argppx.h"
     &                icode,cmessage)
        End If

        If (icode >  0) then
          cmessage='FIELDOP : Error in MODEL DUMP'
! DEPENDS ON: ereport
          CALL EREPORT('WRITEFF', icode,                                &
     &     cmessage)

        End if

      Else
        cmessage='FIELDOP : Pack type not supported'
! DEPENDS ON: ereport
        CALL EREPORT('WRITEFF', 1003,                                   &
     &   cmessage)

      End if

! DEPENDS ON: setpos
      call setpos(pp_unit_out,fixhd(150)+(entry_no-1)                   &
     &            *len1_lookup-1,icode)
! DEPENDS ON: buffout
      call buffout(pp_unit_out,lookup(1,entry_no),fixhd(151)            &
     &             ,len_io,A)

      ! Check for I/O errors
      If (A  /=  -1.0 .or. len_io  /=  fixhd(151)) then

! DEPENDS ON: ioerror
        call ioerror('buffer out of lookup table',A,len_io,             &
     &               fixhd(151))
        cmessage='FIELDOP: I/O error'
        icode=25
        return
      End if

      If (op  ==  'add     ') then

        fixhd(15) = 100
      Else if (op  ==  'subtract') then

        fixhd(15) = 200
      Else if (op  ==  'multiply') then

        fixhd(15) = 300
      Else if (op  ==  'idiv') then

        fixhd(15) = 400
      End if

! DEPENDS ON: setpos
      call setpos(pp_unit_out,0,icode)
! DEPENDS ON: buffout
      call buffout(pp_unit_out,fixhd(1),len_fixhd,len_io,A)

      ! Check for I/O errors
      If(A  /=  -1.0 .or. len_io  /=  len_fixhd) then
! DEPENDS ON: ioerror
        call ioerror('buffer out of fixed length header',A,len_io       &
     &               ,len_fixhd)
        cmessage='FIELDOP: I/O error'
        icode=1
        return
      End if
      RETURN
      END SUBROUTINE writeff

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
