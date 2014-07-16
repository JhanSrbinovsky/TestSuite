#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Signal end of sub-stepping to SCM output diagnostic system

      SUBROUTINE scm_substepping_end
      USE global_SCMop
      IMPLICIT NONE

! Description:
      ! In order to be able to output diagnostics that are defined
      ! within sub-stepped parts of the model, the diagnostic system
      ! needs to know which parts of the model are sub-stepped and
      ! which are not. The routine scm_substep_start is used to let
      ! the system know that a particular sub-step has started. This
      ! routine is called to let the system know that sub-stepping
      ! has finished, and is not intended to be called from inside the 
      ! sub-stepping loop like scm_substep_start, but from outside. 
      ! This pair of routines is used for every section of the code 
      ! that is sub-stepped.

! Method:
! Owner: Luke Jones
! History:
! Version  Date      Comment
! =======  ====      =======
! 6.2      19/01/06  Original code (Luke Jones)
! Code Description:
      ! Language: Fortran77 with some bits of Fortran90

      ! A sub-step number of zero will indicate we are no longer in a
      ! sub-stepping part of the code.

      SCMop%substep_number=0

      RETURN
      END SUBROUTINE scm_substepping_end
#endif
