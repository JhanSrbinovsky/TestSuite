! ----------------------- Comdeck: NATFORCE  ----------------------------
! Description: COMDECK containing a common block for natural forcing filenames
!
! Author : C.D.Jones
!
! History:
! Version  Date      Comment.
!  6.6.2  11/06/09  New code. C.D.Jones

      CHARACTER*(120) file_scvary
      CHARACTER*(120) file_volcts

      COMMON /FILENATFORCE/                                                 &
     &  file_scvary, file_volcts
      NAMELIST /FILENATFORCE/                                               &
     &  file_scvary, file_volcts
