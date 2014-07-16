#if defined(ATMOS) && defined(A35_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine RPOpenInput
! -----------------------
! Opens the file stream constructing a filename from
! the the RPSEED variables that has been set in a
! namelist. If this is blank or just spaces, then
! construct a name from the $DATAW and $RUNID
! environment variables.
!
! Random seed will only be read in if L_RPSEED_READ is true
! this implies that L_RPSEED_WRITE must be false to prevent 
! overwriting the seed read in
!
Subroutine RPOpenInput()
IMPLICIT NONE

#include "chsunits.h"
#include "clfhist.h"
#include "cntlatm.h"
#include "parvars.h"

! Local variables

INTEGER           :: icode       ! Error return strings
CHARACTER(len=64) :: datawString ! To hold DATAW env var
CHARACTER(len=5)  :: runidString ! To hold RUNID env var
CHARACTER(len=80) :: fileName    ! Name of filename
LOGICAL           :: fileExists  ! Check if file already exists

! Checks if output is required, and if unit is already open
 IF (mype == 0) THEN
  IF (L_RPSEED_READ) THEN
    L_RPSEED_WRITE=.false.
    ! -----------------------
    ! Determine file name 
    ! -----------------------
    Filename=''
    CALL fort_get_env('RPSEED',6,Filename,64,icode)
    Filename=trim(Filename)
    
    IF (icode == 0) THEN
        WRITE(6, *) 'Successfully retrieved filename', Filename
    ELSE 
        WRITE(6, *) 'Problem retrieving filename from RPSEED'
        WRITE(6, *) 'RPSEED =', RPSEED
        WRITE(6, *) 'Filename =', Filename 

        CALL fort_get_env('DATAW',5,datawString,80,icode)
        IF (icode /= 0) WRITE(6,*) 'RPSEED IN: Failed to get value of $DATAW'
    
        CALL fort_get_env('RUNID',5,runidString,5,icode)
        IF (icode /= 0) WRITE(6,*) 'RPSEED IN: Failed to get value of $RUNID'

        fileName = trim(datawString) // "/" // trim(runidString) // ".seed"
        WRITE(6, *) 'Will use Filename=', Filename
    ENDIF
    
    IF (len_trim(fileName) == 0) THEN
      CALL fort_get_env('DATAW',5,datawString,64,icode)
      IF (icode /= 0) THEN
        WRITE(6,*) 'RPSEED IN: Failed to get value of $DATAW'
      ENDIF
      CALL fort_get_env('RUNID',5,runidString,5,icode)
      IF (icode /= 0) then
        WRITE(6,*) 'RPSEED IN: Failed to get value of $RUNID'
      ENDIF
      fileName = trim(datawString) // "/" // trim(runidString) // ".seed"
       WRITE(6,*) 'RPSEED IN FILE :',filename
    ENDIF
    
    ! -----------------------
    !  Check that the file exists 
    ! -----------------------
    INQUIRE(file=fileName,exist=fileExists)
    IF (fileExists) THEN
       CLOSE(unit=149)
       OPEN(unit=149,file=fileName,status='old',action='read',     &
            recl=20,iostat=icode)

       IF (icode == 0)  WRITE(6,*) 'RPSEED IN: File has been opened'
       IF (icode /= 0)  WRITE(6,*) 'RPSEED IN: Unable to open ..', Filename
     ELSE  
       WRITE(6,*) "RPSEED IN: "//filename//                        &
                  "does not exist cannot read"//                   &
                  " old seed ! Set L_RPSEED_READ=FALSE in umui   "
    ENDIF                        
  ENDIF  
 ELSE
   ! If we are not PE0 then ensure L_RPSEED_READ is set to false
   L_RPSEED_READ = .false.
 ENDIF
END SUBROUTINE RPOpenInput
#endif
