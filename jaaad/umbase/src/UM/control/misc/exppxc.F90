#if defined(C70_1A) || defined(UTILIO) || defined(FLDOP)               \
 || defined(VAROPSVER)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Integer Function to extract data from lookup array PPXI
!
! Function Interface:

!---------------------------------------------------------------------
!+ Character Function to extract names from lookup array PPXC
!
! Function Interface:
      CHARACTER*36 FUNCTION EXPPXC(Im_ident,section,item,               &
#include "argppx.h"
     &                                   ErrorStatus ,CMESSAGE)
      IMPLICIT NONE
!
! Description:
!   Extracts a diagnostic name from ppxref lookup array PPXC.
!
! Method:
!   The required name is identified by the function arguments
!   Im_ident, section, item. The appropriate row in PPXC is found
!   from the 3-d pointer array PPXPTR as PPXPTR(m,s,i).
!
! Current code owner:  S.J.Swarbrick
!
! Code description:
!   FORTRAN 77 + common Fortran 90 extensions.
!   Written to UM programming standards version 7.
!
! System component covered:
! System task:               Sub-Models Project
!
! Global Variables:
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"

! Function arguments:
!   Scalar arguments with intent(in):
      INTEGER Im_ident    ! Internal model identifier (absolute)
      INTEGER section     ! STASH section no.
      INTEGER item        ! STASH item no.

!   Scalar arguments with intent(out):
      CHARACTER*80 CMESSAGE

! Local scalars
      INTEGER row         ! Row no. in PPXC array
      INTEGER I           ! Loop counter
      INTEGER Im_index    ! Internal model index

! Error status:
      INTEGER ErrorStatus !+ve = fatal error

!- End of Header ---------------------------------------------------

      IF (Im_ident <= 0 .OR. section <  0 .OR. item <= 0) THEN
        IF (Im_ident <= 0) WRITE(6,*) 'EXPPXC: INVALID Im_ident'
        IF (section  <  0) WRITE(6,*) 'EXPPXC: INVALID SECTION NO.'
        IF (item     <= 0) WRITE(6,*) 'EXPPXC: INVALID ITEM NO.'
        WRITE(6,*)                                                      &
     & 'Im_ident ',Im_ident,' section ',section,' item ',item
        ErrorStatus=1
        CMESSAGE='ERROR EXPPXC: INVALID STASH RECORD ID'
      ELSE

! Obtain row no. in PPXC array
#if defined(CONVIEEE) || defined(CONVPP)           \
 || defined(CUMF) || defined(PUMF) || defined(MERGE)
      row = PPXPTR(Im_ident,section,item)
      IF (ROW /= 0) THEN
#else
      Im_index = INTERNAL_MODEL_INDEX(Im_ident)
      row = PPXPTR(Im_index,section,item)
      IF (row <= 0) THEN
        WRITE(6,*) 'ERROR in EXPPXC: INVALID row VALUE: ',row
        WRITE(6,*) 'Model,Sec,Item: ',Im_ident,section,item
!        ErrorStatus = 1
!        CMESSAGE='ERROR EXPPXC: INVALID row VALUE'
      END IF
#endif

! Obtain required name

      IF (row >  0) THEN
        DO I = 1,PPXREF_CHARLEN
          EXPPXC(I:I) = PPXC(row,I)
        END DO
      END IF
#if defined(CONVIEEE) || defined(CONVPP)           \
 || defined(CUMF) || defined(PUMF) || defined(MERGE)
      ELSE IF ((ITEM >= 177.AND.ITEM <= 179).OR.                        &
     &         (ITEM >= 301.AND.ITEM <= 324)) THEN
        ErrorStatus=100
      ELSE
        ErrorStatus=101
      END IF
#endif
      END IF
      RETURN
      END FUNCTION EXPPXC

!---------------------------------------------------------------------
#endif
