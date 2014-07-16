

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: PR_F_C1A ------------------------------------------------
!LL
!LL  Purpose: To print messages generated by the I/O routines in C
!LL
!LL  Author:  Bob Carruthers, Cray Research.   Date: 31 June 1997
!LL
!LL Version Date      Modification history
!LL  4.5    09/10/98  Added def UTILHIST to top def line to allow
!LL                   history executables to be built. K Rogers
!    6.0    16/09/03  Use Fortran intrinsic len_trim instead of
!                     get_char_len. D.Robinson
!    6.2    23/11/05  Removed DIAG92 diagnostic define T.Edwards
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
      SUBROUTINE PRINT_FROM_C(unit, message)
      implicit none
!
      integer unit              ! the nuit number that the message
                                ! applies to, or -1 if there is no
                                ! unit number.
!
      character*(*) message     ! the text of the message
!
      character*256 old_message ! the text of the last message
!
      integer                                                           &
     &  old_unit                                                        &
                                ! the unit number for the last message
     & ,i                                                               &
                                ! loop index
     & ,message_count                                                   &
                                ! count of the number of times the last
                                ! message appeared, after the first time
     & ,l                                                               &
                                ! length of the message from
                                ! get_char_len
     & ,old_l                   ! length of the old message
!
      save old_unit, old_message, message_count
      data old_unit/-2/, message_count/0/, old_l/1/,                    &
     & old_message/'? ? ? ? ? ? ? ? ?'/
!
      l=len_trim(message)
      if(old_unit == unit .and. old_message(1:l) == message(1:l)) then
        message_count=message_count+1
      else
        if(message_count >  0) then
          write(6,'(a,'' - Repeated '',i3,                              &
     &     '' Time(s)'')') old_message(1:old_l), message_count
        endif
        if(unit /= old_unit) write(6,'('' '')')
        write(6,'(a)') message(1:l)
        message_count=0
        old_message=message
        old_l=l
        old_unit=unit
      endif
      return
      END SUBROUTINE PRINT_FROM_C
