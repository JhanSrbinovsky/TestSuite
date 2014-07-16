#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to perform a shell sort.
!
! Method:
!       The standard shell sorting algorithm is implemented.
!
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SHELL_SORT(N, POINTER, KEY)
!
!
!
      IMPLICIT NONE
!
! Include header file
#include "c_kinds.h"
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N
!           Number of elements
      INTEGER, INTENT(INOUT) ::                                         &
     &    POINTER(N)
!           Pointer to succeeding elements
      REAL  (Real64), INTENT(IN) ::                                     &
     &    KEY(N)
!           Key for sorting
!
!     Local variables.
      INTEGER                                                           &
     &    GAP                                                           &
!           Searching interval
     &  , POINTER_TEMP                                                  &
!           Temporary value of pointer
     &  , J                                                             &
!           Loop variable
     &  , K
!           Loop variable
!
!
      IF (N == 1) THEN
        POINTER(1)=1
        RETURN
      ENDIF
!
      GAP=N
      DO WHILE(GAP >= 2)
        GAP=GAP/2
        DO J=GAP, N-1
          DO K=J-GAP+1, 1, -GAP
            IF (KEY(POINTER(K)) >  KEY(POINTER(K+GAP))) THEN
                    POINTER_TEMP=POINTER(K)
              POINTER(K)=POINTER(K+GAP)
              POINTER(K+GAP)=POINTER_TEMP
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SHELL_SORT
#endif
