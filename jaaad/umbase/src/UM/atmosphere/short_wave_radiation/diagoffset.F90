#if defined(A70_1C) || defined(A70_1Z) \
 || defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:
!   Set the offset off the item numbers for the forcing (and radiance)
!   diagnostics in STASH.
!
!
! Method:
!
! Current Code Owner: Jean-Claude Thelen
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 6.2       13/02/06  Original code.  Jean-Claude Thelen
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.
!
! Declarations:
!
! Global variables (#include statements etc):
MODULE DIAGOFFSET
!
!     Offsets for diagnostics from successive calls to radiation code
      INTEGER, Parameter :: diagnostic_offset = 200
!
END MODULE DIAGOFFSET
#endif
