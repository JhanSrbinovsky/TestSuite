#if defined(A03_8A)||defined(A03_8B)||defined(A03_8C) \
  ||defined(A19_1A)||defined(A19_2A)
! Description:
!   This deck sets up the parameter LB
!
! Current Code Owner: Z Gardner
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 5.3      16/08/01 Added in header and changed definitions. Z Gardner
! 5.5      17/04/03 Remove reference to obsolete section
!                   A03_7A. T.White
!
! Declarations:
! Start blend_h
! Description:
!   This file sets the value of the variable LB
!
! Current Code Owner:
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   5.3   25/09/01  Portability changes.  Z. Gardner

      REAL,PARAMETER:: LB = 20.0 ! Blending height (m).
#endif
! End Blend_h
