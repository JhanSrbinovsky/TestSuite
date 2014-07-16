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
      subroutine read_rec(field,                                        &
     &                    num_cray_words,                               &
     &                    iwa,                                          &
     &                    pp_unit1,                                     &
     &                    max_len,                                      &
     &                    icode,                                        &
     &                    pack_type)

      IMPLICIT NONE
!
!
! Description: To read a data record from a  pp file/dump.
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
! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER                                                           &
     & num_cray_words,                                                  &
                              !IN  No of CRAY words holding the data
     & max_len,                                                         &
     & pp_unit1,                                                        &
                              !IN  unit no of the PP FILE
     & iwa,                                                             &
                              !IN  WORD address of field to be read
     & pack_type              !IN  packing type of the field

!   Scalar arguments with intent(out):
      INTEGER                                                           &
     & icode                  !OUT error code

!   Array arguments with intent(out):
      REAL                                                              &
     & field(max_len)         !OUT array holding data

! Local scalars.
      INTEGER                                                           &
     & i,j,                                                             &
                                ! local counter
     & len_io                   ! length of data read by buffin

      REAL                                                              &
     & a_io                     ! return code from buffin

! Function & Subroutine calls:
      External setpos,buffin

!- End of header

! DEPENDS ON: setpos
      call setpos(pp_unit1,iwa,icode)
      ! We have to read in the fields differently on little
      ! endian machines if there are 32 bits packed into a
      ! 64bit field, (results in swapping of numbers.
      ! Check if we are using this type of cray packing and read
      ! in twice as many 32 bit numbers as 64 bit
      if(pack_type  ==  2) then
! DEPENDS ON: buffin32
        call buffin32(pp_unit1,field,num_cray_words*2,len_io,a_io)
      else
! DEPENDS ON: buffin
        call buffin(pp_unit1,field,num_cray_words,len_io,a_io)
      end if
      RETURN
      END SUBROUTINE read_rec
!
! Subroutine interface:
!
!Subroutine interface:
!
!Subroutine interface:
#endif
