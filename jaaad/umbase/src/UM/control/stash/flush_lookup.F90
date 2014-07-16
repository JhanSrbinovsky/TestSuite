#if defined(C84_1A)
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
       Subroutine flush_lookup(unit)
! Description:
!   This routine flushes bufffered lookups to disk

          Use Field_Buff_Mod, Only :                                    &
     &        pp_buff

          Use flumerun  ! FLUME-STASH

          Implicit None

          Integer :: unit

          Integer :: iwa, icode, nwords
          Integer :: len_io
          Real :: iostat

          ! Only do the flushing if we actually are storing
          ! a lookup table.
          If (associated(pp_buff(unit) % pp_lookup_table)) Then

            iwa=pp_buff(unit) % pp_lookup_address
            nwords=pp_buff(unit) % len1_pp_lookup                       &
     &            *pp_buff(unit) % len2_pp_lookup

! DEPENDS ON: setpos
            call setpos(unit,iwa,icode)
            call buffout_single(unit,pp_buff(unit) % pp_lookup_table,   &
     &                          nwords,len_io,iostat)
            IF (Flume_run) Call setpos(unit,0,icode) ! FLUME-STASH
          End If


       End Subroutine flush_lookup

!----------------------------------------------------------------------
! Subroutine test_and_flush_pp
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Subroutine flush_all_pp
!----------------------------------------------------------------------


!----------------------------------------------------------------------
! Subroutine init_buffers_pp
!----------------------------------------------------------------------


#endif
