#if defined(C97_4A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.3   03/10/01   Add Implicit none, headers and history. R.Sharp
!   6.2   15/08/05   Remove RECON def. P.Selwood
      SUBROUTINE TIMER(SUB,I)
! Dummy routine called when no timing info is required
      Implicit None
!
      CHARACTER*(*) SUB
      INTEGER I
      RETURN
      END SUBROUTINE TIMER
#endif
