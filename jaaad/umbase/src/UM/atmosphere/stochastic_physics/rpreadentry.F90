#if defined(CONTROL) && defined(ATMOS) && defined(A35_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine RPReadEntry(value)
! ------------------------------
! This routine reads in the random seed which is the passed
! up to the calling routine as a subroutine argument.

SUBROUTINE RPReadEntry(seed,seedsize)
IMPLICIT NONE

#include "cntlatm.h"
#include "parvars.h"
#include "cprintst.h"

! Arguments
 INTEGER, INTENT(IN)     :: seedsize       ! The size of the random seed
 INTEGER(8), INTENT(OUT) :: seed(seedsize) ! The random seed

! Local Variables
 INTEGER           :: i               ! loop variable
 INTEGER           :: icode           ! 0 if successful read
 LOGICAL           :: opened, named
 INTEGER           :: rlength
 CHARACTER(LEN=80) :: fname 
 CHARACTER(LEN=9)  :: dir_acc, rw_sts, r_sts
 
! If the unit is open the write the value to the file stream
 IF  (L_RPSEED_READ) THEN 
  IF (mype == 0) THEN
    DO i=1,seedsize
      INQUIRE(unit=149,NAMED=named, OPENED=opened, NAME=fname, &
              READ=r_sts,READWRITE=rw_sts,RECL=rlength,DIRECT=dir_acc)
      IF (I == 1) THEN
       IF (PrintStatus  >  PrStatus_Normal) THEN              
        WRITE(6,*) 'Readwrite', rw_sts
        WRITE(6,*) 'Rec length', rlength
        WRITE(6,*) 'read', r_sts
        WRITE(6,*) 'Direct_acces', dir_acc
        WRITE(6,*) 'opened', opened
        WRITE(6,*) 'named', named
        WRITE(6,*) 'fname',fname              
       ENDIF  ! Print
      ENDIF! i=1
      READ(149,'(256i20)',advance='yes') seed(i)
    ENDDO    ! seedsize
   ENDIF     ! mype 
 ENDIF       ! L_RPSEED_READ
 END SUBROUTINE RPReadEntry
#endif
