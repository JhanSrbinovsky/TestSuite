#if defined(C84_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: STMAX ----------------------------------------------------
!LL
!LL  Purpose: Computes the point-by-point maximum in time of a field
!LL           by comparing the field at the current time with the
!LL           maximum so far (STASH TEMPORAL service routine)
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Author:   S.Tett/T.Johns
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: D722
!LL
!LL  Project task: D7
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage handling and diagnostic
!LL                                 system (STASH)
!LL
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE STMAX(fieldin,result,size,masking,amdi)
!
      IMPLICIT NONE
!
      INTEGER size             ! IN size of fieldin and result.
      REAL fieldin(size)       ! IN input field
      REAL result(size)        ! OUT output field (maximum)
      LOGICAL masking          ! IN true if masked (ie. MDI possible)
      REAL amdi                ! IN missing data indicator
!*----------------------------------------------------------------------
!
! Local variables
!
      INTEGER i ! loop count
!-----------------------------------------------------------------------
!L 1.1 loop over array size, if either result or fieldin is amdi set
!L     result to amdi, else set result to maximum of fieldin and result
!L
      IF (masking) THEN
        DO i=1,size
          IF ((result(i) /= amdi).and.(fieldin(i) /= amdi)) THEN
            result(i)=max(result(i),fieldin(i))
          ELSE
            result(i)=amdi
          ENDIF
        ENDDO
      ELSE
!L
!L 1.2 loop over array size, set result to maximum of fieldin and result
!L     without checking for missing data
!L
        DO i=1,size
          result(i)=max(result(i),fieldin(i))
        ENDDO
      ENDIF
!
      RETURN
      END SUBROUTINE STMAX
#endif
