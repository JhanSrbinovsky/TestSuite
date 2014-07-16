#if defined(ATMOS) && defined(A35_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine RPOpenOutput
! -----------------------
! Opens the file stream constructing a filename from
! the the RPSEED variables that has been set in a
! namelist. If this is blank or just spaces, then
! construct a name from the $DATAW and $RUNID
! environment variables.
!
! Random seed will always be written out if L_RPSEED_WRITE is true
! however should only be written out once per run ,i.e, the
! first time STPH_RP is called .
!------------------------------------------------------------------- 
Subroutine RPOpenOutput()
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
  IF (L_RPSEED_WRITE) THEN

    Filename=''
    CALL fort_get_env('RPSEED',6,Filename,80,icode)
    Filename=trim(Filename)
    
    IF (icode == 0) THEN
        WRITE(6, *) 'Successfully retrieved filename', Filename
    ELSE 
        WRITE(6, *) 'Problem retrieving filename from RPSEED'
        WRITE(6, *) 'RPSEED =', RPSEED
        WRITE(6, *) 'Filename =', Filename 

        CALL fort_get_env('DATAW',5,datawString,64,icode)
        IF (icode /= 0) WRITE(6,*) 'RPSEED OUT: Failed to get value of $DATAW'
    
        CALL fort_get_env('RUNID',5,runidString,5,icode)
        IF (icode /= 0) WRITE(6,*) 'RPSEED OUT: Failed to get value of $RUNID'

        fileName = trim(datawString) // "/" // trim(runidString) // ".seed"
        WRITE(6, *) 'Will use Filename=', Filename
    ENDIF
    
    
    IF (len_trim(fileName) == 0) THEN
      WRITE(6,*) 'length of filename is equal to zero'
      CALL fort_get_env('DATAW',5,datawString,64,icode)
      IF (icode /= 0) WRITE(6,*) 'RPSEED OUT: Failed to get value of $DATAW'
      WRITE(6,*) 'DATAW',datawString
    
      CALL fort_get_env('RUNID',5,runidString,5,icode)
      IF (icode /= 0) WRITE(6,*) 'RPSEED OUT: Failed to get value of $RUNID'
      WRITE(6,*) 'RUNID',runidstring

      fileName = trim(datawString) // "/" // trim(runidString) // ".seed"
      WRITE(6, *) 'Filename', Filename
    ENDIF
    
    !-------------------------------------------------------------------------
    ! Note if CRUN we assume a new random seed is written out at the beginning
    ! of the CRUN to a new file. So we treat no differently to a standard run.
    !------------------------------------------------------------------------
    OPEN(unit=149,file=fileName,status='replace',action='write', iostat=icode)
    IF(icode > 0) WRITE(6,*) "RPSEED OUT: Failed to open ", Filename
    IF(icode == 0) WRITE(6,*) "RPSEED OUT: Successfully opened file ", Filename
    
  ENDIF  ! L_RPSEED_WRITE
 ELSE
! If we are not PE0 then ensure L_SEEDWRITE is set to false
   L_RPSEED_WRITE = .false.
 ENDIF
END SUBROUTINE RPOpenOutput

#endif
