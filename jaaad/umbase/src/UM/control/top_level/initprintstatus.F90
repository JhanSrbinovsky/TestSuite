#if defined(CONTROL) || defined(UTILIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initialise print status for standard output
!
! Subroutine Interface:
      Subroutine InitPrintStatus

      IMPLICIT NONE

!
! Description:
!   Initialises value of PrintStatus variable to control subsequent
!   printing of standard output to unit 6.
! Method:
! Use environment variable available from UMUI supplied to top level
! SCRIPT file to set switch PrintStatus, accessed locally by both
! fortran and C code.
! Current Code Owner: R Rawlins
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.0   29/04/99   Original code based on operational modset used at
!                   vn4.4. R Rawlins
!  5.1   11/05/00   Added defined(UTILIO). D.P.Matthews
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Global variables (*CALLed COMDECKs etc...):
! Print status information in CPRINTST:
#include "cprintst.h"
! Subroutine arguments
! Local parameters:
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='InitPrintStatus')
! Local scalars:
      INTEGER ErrorStatus       ! Error return code
      CHARACTER*80 PRINT_STATUS   ! Text content of env variable
! Local dynamic arrays:
! Function & Subroutine calls:
      External Fort_Get_Env                                             &
     &        ,Set_Printstatus
!- End of header

      CALL Fort_Get_Env('PRINT_STATUS',12,PRINT_STATUS,80,ErrorStatus)

      IF(ErrorStatus /= 0) THEN
        WRITE(6,*) RoutineName,': Warning: problem in reading ',        &
     &  'environment variable PRINT_STATUS, ErrorStatus=',ErrorStatus
      ENDIF

      IF    (PRINT_STATUS(1:13) == 'PrStatus_Diag') THEN
        PrintStatus = PrStatus_Diag
      ELSEIF(PRINT_STATUS(1:15) == 'PrStatus_Normal') THEN
        PrintStatus = PrStatus_Normal
      ELSEIF(PRINT_STATUS(1:13) == 'PrStatus_Oper') THEN
        PrintStatus = PrStatus_Oper
      ELSE
        PrintStatus = PrStatus_Min
      ENDIF

! Set PrintStatus in C code for control of I/O messages
      CALL Set_Printstatus(PrintStatus)

      IF(PrintStatus >= PrStatus_Oper) THEN
        write(6,*) RoutineName,': PrintStatus initialised=',PrintStatus
      ENDIF

      RETURN
      END SUBROUTINE InitPrintStatus
#endif
