! CGRIBTAB start
!
! Description:
!  Holds a lookup table for converting between Unified Model
! stash codes and some other grib table of codes
! Current Code Owner: R.A.Stratton
!
! History:
! Version  Date     Comment
! -------  ----     -------
!  4.0    31/03/95 :Original code. R.A.Stratton
!
! Declarations:
      INTEGER, PARAMETER:: MAX_SECT_GRBTAB=16
      INTEGER, PARAMETER:: MAX_ITEM_GRBTAB=300

      INTEGER :: GRIB_TABLE(0:MAX_SECT_GRBTAB,MAX_ITEM_GRBTAB)

      COMMON /GRIBTAB/ GRIB_TABLE

! CGRIBTAB end
