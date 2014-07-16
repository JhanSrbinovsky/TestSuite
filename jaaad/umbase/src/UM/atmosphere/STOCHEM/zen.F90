#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      REAL FUNCTION ZEN(time,xlat,xlong)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Calculates local solar zenith angle for
!-                         particular time
!-
!-   Inputs  : TIME,xlat,xlong
!-   Outputs : none
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  3.4    13/01/95  Created.  W.J. Collins
!  5.2    04/09/00  Automatic 360/365 day selection.  C.E. Johnson
!  5.5    16/02/04  Function type declared directly.  K. Ketelsen
!  6.1    21/10/04  Reformatted code. M.G. Sanderson
!
!-
!VVV  V2.6  ZEN 4/IX/00  360/365 selection
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!----------------------------------------------------------------------
      REAL, INTENT(IN) :: TIME
      REAL, INTENT(IN) :: XLAT,XLONG  ! lat & long are arrays in GRD

      REAL :: xlha, decl, cosine

! Local hour angle
      xlha = (1.0+time/4.32e4+xlong/180.0)*pi

! Declination from DAYM sum automatically selects 360/365 day calendar
      decl = -0.4142*COS(pi+(2.0*pi*time/(daym_all*86400.0)))
      cosine = COS(xlha) * COS(xlat*pi/180.0) * COS(decl)               &
     &  + SIN(xlat*pi/180.0) * SIN(decl)
      ZEN = ACOS(cosine)

      END FUNCTION ZEN
#endif
