
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! A collection of routines to flush buffers of STASH data
!
! Current Code Owner: Paul Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.0   15/01/04   Original version. JC Rioual, P. Selwood
!   6.1   16/08/04   Updated for C IO with buffering of lookups.
!                                      JC Rioual, P. Selwood
!   6.2   17/05/06   Fix out of bounds for lookups when closing
!                    read-only files. P.Selwood.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!

!----------------------------------------------------------------------
! Subroutine flush_lookup
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Subroutine test_and_flush_pp
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Subroutine flush_all_pp
!----------------------------------------------------------------------

      Subroutine flush_all_pp()

! Description:
!    This routine runs thorough all possible STASH units and ensures
!    all buffered data is flushed to disk

      Use Field_Buff_Mod, Only :                                        &
     &    BUFF_MIN,                                                     &
     &    BUFF_MAX

      Implicit None

! Local variables
      Integer :: i

      Do i = BUFF_MIN, BUFF_MAX
! DEPENDS ON: test_and_flush_pp
         Call test_and_flush_pp(i)
      End Do

      End Subroutine flush_all_pp

!----------------------------------------------------------------------
! Subroutine init_buffers_pp
!----------------------------------------------------------------------


