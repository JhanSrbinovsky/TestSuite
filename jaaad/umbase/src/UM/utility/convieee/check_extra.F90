#if defined(PUMF) || defined(FLDC) || defined(IEEE)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine interface:



!LL  Routine: CHECK_EXTRA ----------------------------------------------
!LL
!LL  Purpose: To check that code is correct for vector
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model               Modification history:
!LL version  Date
!LL  6.0     02/07/03    Move subroutine CHECK_EXTRA and function
!LL                      INT_FROM_REAL in deck FIELDCOS to deck PREXTRA
!LL                      (minimize exec_xref decks on NEC)     E.Leung
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered:
!LL
!LL  Project task:
!LL
!LL  External documentation:
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
      SUBROUTINE CHECK_EXTRA(CODE,DATA_VALUES,ICODE,CMESSAGE)
      IMPLICIT NONE
#include "comextra.h"
      LOGICAL VALID_TYPE ! Flag to indicate valid vector type
      INTEGER I ! LOOP counter
!     arguments
      CHARACTER                                                         &
     &     CMESSAGE*(*)         !OUT error message
      INTEGER                                                           &
     &     CODE                                                         &
                                !IN Code to be checked
     &    ,DATA_VALUES                                                  &
                                !IN Number of data values in vector
     &    ,ICODE                !OUT error code
!     Local variables
      INTEGER                                                           &
     &     TYPE

      DATA_VALUES=CODE/1000
      TYPE=CODE-DATA_VALUES*1000
      VALID_TYPE = .FALSE.
      DO I = 1, NO_EXTRA_VECTORS
         IF (TYPE == EXTRA_VECTOR(I)) VALID_TYPE = .TRUE.
      ENDDO

      IF (.NOT.VALID_TYPE) THEN
         ICODE = 1
         CMESSAGE='CHECK_DATA: Unsupported extra data vector code'
         WRITE(6,*) "Unsupported Extra Data code", TYPE
      ENDIF

      RETURN
      END SUBROUTINE CHECK_EXTRA



#endif
