#if defined(C70_1A) || defined(FLDIO) || defined(UTILHIST) \
 || defined(UTILIO) || defined(FLUXPROC)
!  Routine: EREPORT --------------------------------------------------
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!  Purpose: Reports error exit code and message at end of model run.
!
!
!  Programming standard: UM Doc Paper 1, version 1 (15/1/90)
!
!  Logical components covered: C0
!
!  Project task: C0
!
!  External documentation: On-line UM document C0 - The top-level
!                          control system
!
!  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE EREPORT (RoutineName,ICODE,CMESSAGE)
      IMPLICIT NONE
      INTEGER ICODE            ! In - Error code from model
      CHARACTER *(*)                                                    &
     &  RoutineName                                                     &
                                ! IN : Name of routine with error
     &, CMESSAGE                ! IN : Error message
!
#include "parvars.h"
!*----------------------------------------------------------------------
!  Local variables
      INTEGER                                                           &
     &  i,j                                                             &
                    ! loop bounds
     &, len_message                                                     &
                    ! length of text within CMESSAGE
     &, flush_code  ! return code from FLUSH

!
!L----------------------------------------------------------------------
!L 1. Write informative message summarising completion state of model
!L

      WRITE(6,1000)

      IF (ICODE  <   0) THEN ! Warning message

        WRITE(6,*) 'UM WARNING :'
        WRITE(6,*) 'Routine generating warning: ',RoutineName
        WRITE(6,*) 'Warning code: ',ICODE
        WRITE(6,*) 'Warning message: '

      ELSEIF (ICODE  >   0) THEN ! Errors also go to Unit0 (stderr)

        WRITE(6,*) 'UM ERROR (Model aborting) :'
        WRITE(6,*) 'Routine generating error: ',RoutineName
        WRITE(6,*) 'Error code: ',ICODE
        WRITE(6,*) 'Error message: '

        IF (mype  ==  0) THEN
          WRITE(0,1000)
          WRITE(0,*) 'UM ERROR (Model aborting) :'
          WRITE(0,*) 'Routine generating error: ',RoutineName
          WRITE(0,*) 'Error code: ',ICODE
          WRITE(0,*) 'Error message: '
        ENDIF

      ENDIF

      len_message=len_trim(CMESSAGE)

      DO i=1,len_message,80
        j=MIN(len_message,i+79)

        WRITE(6,'(a)') CMESSAGE(i:j)
        IF (ICODE  >   0) THEN ! Error
          IF (mype  ==  0) THEN
            WRITE(0,'(a)') CMESSAGE(i:j)
          ENDIF
        ENDIF

      ENDDO

      WRITE(6,1000)
      IF (mype  ==  0) THEN
        IF (ICODE  >   0) WRITE(0,1000)
      ENDIF

! DEPENDS ON: um_fort_flush
      CALL UM_FORT_FLUSH(6,flush_code)

      IF (ICODE  >   0) THEN ! Fatal error
#if defined(T3E)
        CALL ABORT('T3E Hard Abort')
#else
        CALL GC_ABORT(mype, nproc, cmessage)
#endif
      ENDIF


      IF (ICODE  <   0) ICODE=0   ! reset error code if warning

 1000 FORMAT(" ****************************************",               &
     &       "*****************************************")
      RETURN
!L----------------------------------------------------------------------
      END SUBROUTINE EREPORT
#endif
