#if defined(PUMF) || defined(FLDC) || defined(IEEE)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine interface:



!LL  Routine: CHECK_EXTRA ----------------------------------------------
!LL
!LL  Purpose: To check that code is correct for vector
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model               Modification history:
!LL version  Date
!LL  6.0     02/07/03    Move subroutine CHECK_EXTRA and function
!LL                      INT_FROM_REAL in deck FIELDCOS to deck PREXTRA
!LL                      (minimize exec_xref decks on NEC)     E.Leung
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered:
!LL
!LL  Project task:
!LL
!LL  External documentation:
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------


      INTEGER FUNCTION INT_FROM_REAL(number)
! function to return the integer EQUIVALENCE of a real number
! This function has been moved from deck FIELDCOS to PREXTRA
      integer number
      int_from_real=number
      RETURN
      END FUNCTION INT_FROM_REAL

#endif
