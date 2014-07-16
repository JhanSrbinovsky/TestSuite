
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
       Subroutine test_and_flush_pp(unit)

! Description:
!   This routine checks if the unit is being used, and if it is,
!   flushes it to disk

       Implicit None

! Subroutine Arguments
       Integer, Intent(In) :: unit

! Local variables
       Integer :: ierr
       Integer :: status       ! is file open

       Call is_unit_open(unit,status)
       If (status == 0) then
         ! flush the lookup table to the disk if the file is
         ! still open but keep the info in memory.
! DEPENDS ON: flush_lookup
         Call flush_lookup(unit)

         ! Setpos to 0 will force flushing of data
! DEPENDS ON: setpos
         Call setpos(unit,0,ierr)
       End If

       End Subroutine test_and_flush_pp

!----------------------------------------------------------------------
! Subroutine flush_all_pp
!----------------------------------------------------------------------


!----------------------------------------------------------------------
! Subroutine init_buffers_pp
!----------------------------------------------------------------------


