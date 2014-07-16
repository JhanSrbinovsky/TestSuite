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

!
! Subroutine interface:
!
! Subroutine interface:
      subroutine un_pack(pack_type,                                     &
     &                   npoints,                                       &
     &                   pdata,                                         &
     &                   num_cray_words,                                &
     &                   ilabel,                                        &
     &                   len_ilabel,                                    &
     &                   amdi,                                          &
     &                   pp_fixhd,                                      &
     &                   len_fixhd,                                     &
     &                   pppak,                                         &
     &                   icode,cmessage)
      IMPLICIT NONE
!
! Description: To unpack data from the input array pdata and return
!              the data in pdata.
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

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER                                                           &
     & npoints,                                                         &
                               ! full unpacked size of a pdata
     & max_len,                                                         &
     & num_cray_words,                                                  &
                            ! length of input pdata
     & len_fixhd,                                                       &
     & len_ilabel           ! length of ilabel array

      REAL                                                              &
     & amdi                 ! Missing data indicator.

!   Scalar arguments with intent(in):
      INTEGER                                                           &
     & pp_fixhd(len_fixhd)  ! PPfile fixed length header

!   Scalar arguments with intent(in/out):
      INTEGER                                                           &
     & pack_type            ! Type of packing used

!   Array  arguments with intent(in/out):
      INTEGER                                                           &
     & ilabel(len_ilabel)


!   ErrorStatus
      INTEGER                                                           &
     & icode

      CHARACTER                                                         &
     & cmessage*80

! Local scalars:
      INTEGER                                                           &
     & num_unpack_values,                                               &
                              ! Number of numbers originally packed
     & i,                                                               &
                              ! loop counter
     & ixx,                                                             &
                              ! Returned X dimension from COEX
     & iyy,                                                             &
                              ! Returned Y dimension from COEX
     & idum,                                                            &
                              ! Dummy variable
     & pppak                  ! Packing acc

! Local parameters:
      INTEGER len_full_word   ! The length of a FULL_WORD
        PARAMETER(len_full_word=64)

! Local arrays:
      REAL                                                              &
     & field(npoints)                                                   &
                            !WORK array used for un_packing
     &,pdata(npoints)       ! Input contains packed data.

! Function & Subroutine calls:

!- End of header

      If (pack_type == 1) then     ! WGDOS packing

! DEPENDS ON: coex
        call coex(field,                                                &
                                   ! OUT
     &            npoints,                                              &
                                   ! IN
     &            pdata,                                                &
                                   ! IN
     &            npoints,                                              &
                                   ! IN
     &            ixx,iyy,                                              &
                                   ! OUT
     &            idum,                                                 &
     &            pppak,                                                &
                                   ! OUT
     &            .false.,                                              &
                                   ! IN
     &            amdi,                                                 &
                                   ! IN
     &            len_full_word,                                        &
                                   ! IN
     &            icode,                                                &
                                   ! OUT
     &            cmessage)        ! OUT

        num_unpack_values = ixx * iyy
        ilabel(lblrec) = ilabel(lbrow) * ilabel(lbnpt) + ilabel(lbext)
      Else if (pack_type  ==  2) then !  32 Bit CRAY packing


        num_cray_words = num_cray_words*2
! DEPENDS ON: expand21
        call EXPAND21(num_cray_words,                                   &
                                             ! IN
                      pdata,                                            &
                                             ! IN
                      field)                 ! OUT

        num_unpack_values = num_cray_words


      Else if (pack_type  ==  3) then !  GRIB packing

#if defined(NECSX6)
! DEPENDS ON: degrib
        call degrib(pdata,                                              &
                                        ! IN
     &              field,                                              &
                                   ! OUT
     &              npoints,                                            &
                                           ! IN
     &              num_cray_words,                                     &
                                        ! IN
     &              ilabel,                                             &
                                        ! IN
     &              amdi,                                               &
                                        ! IN
     &              num_unpack_values,                                  &
                                        ! IN
     &              len_full_word)      ! IN
#else
        WRITE(6,*) 'Grib unpacking only supported on NEC SX6.'
! DEPENDS ON: ereport
        CALL EREPORT('UN_PACK', 1000,                                   &
     &   cmessage)

#endif

      Else if (pack_type  ==  4) then ! Run length encoded

! DEPENDS ON: runlen_decode
        call runlen_decode(field,                                       &
                                   ! OUT
     &            npoints,                                              &
                                   ! IN
     &            pdata,                                                &
                                   ! IN
     &            ilabel(lblrec),                                       &
                                          ! IN
     &            amdi,                                                 &
                                   ! IN
     &            icode,                                                &
                                   ! OUT
     &            cmessage)        ! OUT

        num_unpack_values =  npoints
        ilabel(lblrec) = npoints
      Else

        icode=6
        cmessage=' UNPACK - packing type not yet supported'
      End if

      ! Write unpacked data back into array pdata.
      Do i =1,num_unpack_values
       pdata(i) = field(i)
      End do

      ilabel(data_type) =1                ! data must now be real
      ilabel(LBPACK)    =ilabel(LBPACK)-pack_type ! data no
      pack_type         =0                        ! longer packed

      RETURN
      END SUBROUTINE un_pack
!
!Subroutine interface:
!
!Subroutine interface:
#endif
