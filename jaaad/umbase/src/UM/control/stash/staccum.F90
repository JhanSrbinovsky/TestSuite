#if defined(C84_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: STACCUM  -------------------------------------------------
!LL
!LL  Purpose: Accumulates fields within STASH (temporal service routine)
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Author:   S.Tett/T.Johns
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: D721
!LL
!LL  Project task: D7
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage handling and diagnostic
!LL                                 system (STASH)
!LL
!*L  Interface and arguments: ------------------------------------------
      SUBROUTINE STACCUM(fieldin,result,size,masking,amdi)
!
      IMPLICIT NONE
!
      INTEGER size              ! IN size of fieldin and result.
      REAL fieldin(size)        ! IN  input field to be processed
      REAL result(size)         ! OUT where accum is done.
      LOGICAL masking           ! IN true if masked (ie. MDI possible)
      REAL amdi                 ! IN missing data indicator
!* ---------------------------------------------------------------------
!
!  Local variables
!
      INTEGER i                 ! loop count
!L----------------------------------------------------------------------
!L 1.1 loop over array size, if either result or fieldin is amdi, set
!L    result to amdi else accumulate it in result.
!L
      IF (masking) THEN
        DO i=1,size
          IF ((result(i) /= amdi).and.(fieldin(i) /= amdi)) THEN
            result(i)=result(i)+fieldin(i)
          ELSE
            result(i)=amdi
          ENDIF
        ENDDO
      ELSE
!L
!L 1.2 loop over array size, accumulate result without checking for
!L     missing data
!L
        DO i=1,size
          result(i)=result(i)+fieldin(i)
        ENDDO
      ENDIF
!
      RETURN
      END SUBROUTINE STACCUM
#endif
