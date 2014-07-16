#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Bubble sort of STASH list by model/section/item/input code
!
! Subroutine Interface:

      SUBROUTINE ORDER(NRECS)
      IMPLICIT NONE

!  Description:
!    Bubble sort of STASH list by model/section/item/input code
!
!  Method:
!  Loops through records in LIST_S (STASH list) and compares values of
!  model number, then section, then item, then input code. If the values
!  are not in descending order for any one of these categories,
!   the records are interchanged. A logical LSWAP is used to
!   determine when the sort is complete.
!    Called by DUPLIC, STPROC
!
! Current code owner:  S.J.Swarbrick
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Mar. 95    Original code.  S.J.Swarbrick
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
!  Global variables
!

#include "csubmodl.h"
#include "version.h"
#include "cstash.h"
#include "stextend.h"
#include "stparam.h"

!  Subroutine arguments
!    Scalar argument with intent(in):
      INTEGER NRECS       !  No. of records in STASH list

!  Local scalars:

      INTEGER LIST_T !  Used for swapping values in consecutive records
      LOGICAL LSWAP  !  If TRUE, order not complete
      INTEGER I
      INTEGER IREC
      INTEGER IREC1

!- End of Header ------------------------------------------------------

      LSWAP=.TRUE.

  100 CONTINUE            ! return point for order loop

        LSWAP=.FALSE.

!  Loop over records in LIST_S

        DO IREC=NRECS-1,1,-1
          IREC1=IREC+1

!  Sort by internal model

          IF(LIST_S(st_model_code,IREC) >                               &
     &       LIST_S(st_model_code,IREC1)) THEN
             LSWAP=.TRUE.
            DO I=1,NELEMP+1
              LIST_T=LIST_S(I,IREC)
              LIST_S(I,IREC)=LIST_S(I,IREC1)
              LIST_S(I,IREC1)=LIST_T
            END DO
          END IF

!  Sort by section for same internal model

          IF( (LIST_S(st_model_code,IREC) ==                            &
     &         LIST_S(st_model_code,IREC1)).AND.                        &
     &        (LIST_S(st_sect_no_code,IREC) >                           &
     &         LIST_S(st_sect_no_code,IREC1)) ) THEN
            LSWAP=.TRUE.
            DO I=1,NELEMP+1
              LIST_T=LIST_S(I,IREC)
              LIST_S(I,IREC)=LIST_S(I,IREC1)
              LIST_S(I,IREC1)=LIST_T
            END DO
          END IF

!  Sort by item for same model, section

          IF( (LIST_S(st_model_code,IREC) ==                            &
     &         LIST_S(st_model_code,IREC1)).AND.                        &
     &        (LIST_S(st_sect_no_code,IREC) ==                          &
     &         LIST_S(st_sect_no_code,IREC1)).AND.                      &
     &        (LIST_S(st_item_code,IREC) >                              &
     &         LIST_S(st_item_code,IREC1)) ) THEN
            LSWAP=.TRUE.
            DO I=1,NELEMP+1
              LIST_T=LIST_S(I,IREC)
              LIST_S(I,IREC)=LIST_S(I,IREC1)
              LIST_S(I,IREC1)=LIST_T
            END DO
          END IF

!  Sort by input code when model, section, item are correct

          IF( (LIST_S(st_input_code,IREC) <                             &
     &         LIST_S(st_input_code,IREC1)).AND.                        &
     &        (LIST_S(st_model_code,IREC) ==                            &
     &         LIST_S(st_model_code,IREC1)).AND.                        &
     &        (LIST_S(st_item_code,IREC) ==                             &
     &         LIST_S(st_item_code,IREC1)).AND.                         &
     &        (LIST_S(st_sect_no_code,IREC) ==                          &
     &         LIST_S(st_sect_no_code,IREC1)) ) THEN
            LSWAP=.TRUE.
            DO I=1,NELEMP+1
              LIST_T=LIST_S(I,IREC)
              LIST_S(I,IREC)=LIST_S(I,IREC1)
              LIST_S(I,IREC1)=LIST_T
            END DO
          END IF

        END DO     !  Loop over records

      IF (LSWAP) GOTO 100

      RETURN
      END SUBROUTINE ORDER
#endif
