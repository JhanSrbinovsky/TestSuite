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
!
!Subroutine interface:
      subroutine integer_to_real(idim,integer_field,field,nvals,        &
     &                           max_len,ilabel,icode)
      IMPLICIT NONE
!
! Description: Converts integer data into real.
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
     & idim,                                                            &
                            !IN full unpacked size of a field
     & max_len,                                                         &
     & nvals                !IN no of values in an input field

!   Array  arguments with intent(in):
      INTEGER                                                           &
     & integer_field(max_len)  ! contains integer data.

!   Scalar arguments with intent(out):
      INTEGER                                                           &
     & icode                !OUT error code

!   Array arguments with intent(out):
      INTEGER                                                           &
     & ilabel(44)           !OUT integer part of lookup

      REAL                                                              &
     & field(max_len)          !OUT contains Real data.

! Local scalars:
      INTEGER                                                           &
     & i                    ! loop counter

!- End of header

      Do  i =1,nvals
        field(i) = integer_field(i)
      End do

      ilabel(data_type) =1       ! The data type must now be real
      icode=0

      RETURN
      END SUBROUTINE integer_to_real
!
!Subroutine interface:
#endif
