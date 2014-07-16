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
      INTEGER FUNCTION EXPPXI(Im_ident,section,item,element,            &
#include "argppx.h"
     &                              ErrorStatus ,CMESSAGE)
      IMPLICIT NONE
!
! Description:
!   Extracts an individual data value from ppxref lookup array PPXI.
!
! Method:
!   The required data element is identified by the function arguments
!   Im_ident, section, item, element. The appropriate row in PPXI is
!   found from the 3-d pointer array PPXPTR as PPXPTR(m,s,i). The
!   address of the required element in PPXI is then given by
!   (row, element).
!
!
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
      INTEGER element     ! Position of required value in PPXI row

!   Scalar arguments with intent(out):
      CHARACTER*80 CMESSAGE

! Error status:
      INTEGER ErrorStatus !+ve = fatal error

! Local scalars
      INTEGER row         ! Row no. in PPXI array
      INTEGER Im_index    ! Internal model index

!- End of Header ---------------------------------------------------

      ErrorStatus = 0
      IF (Im_ident <= 0 .OR. section <  0 .OR. item <= 0) THEN
        IF (Im_ident <= 0) WRITE(6,*) 'EXPPXI: INVALID Im_ident'
        IF (section  <  0) WRITE(6,*) 'EXPPXI: INVALID SECTION NO.'
        IF (item     <= 0) WRITE(6,*) 'EXPPXI: INVALID ITEM NO.'
        WRITE(6,*)                                                      &
     & 'Im_ident ',Im_ident,' section ',section,' item ',item
        ErrorStatus=1
        CMESSAGE='ERROR EXPPXI: INVALID STASH RECORD ID'
      ELSE

! Obtain row no. in PPXI array
#if defined(CONVIEEE) || defined(CONVPP)           \
 || defined(CUMF) || defined(PUMF) || defined(MERGE)
      row = PPXPTR(Im_ident,section,item)
#else
      Im_index = INTERNAL_MODEL_INDEX(Im_ident)
      row = PPXPTR(Im_index,section,item)
      IF (row <= 0) THEN
        WRITE(6,*) 'ERROR EXPPXI: INVALID row VALUE: ',row
        WRITE(6,*) 'Im_ident,Sec,Item: ',Im_ident,section,item
!        ErrorStatus = 1
!        CMESSAGE='ERROR EXPPXI: INVALID row VALUE'
      END IF
#endif

! Obtain required data value
      IF (row >  0) THEN
        EXPPXI = PPXI(row,element)
      END IF
      END IF
      RETURN
      END FUNCTION EXPPXI

!---------------------------------------------------------------------
!+ Character Function to extract names from lookup array PPXC
!
! Function Interface:

!---------------------------------------------------------------------
#endif
