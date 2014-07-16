#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Append sub-step number to diagnostic's short name

      function add_substep_to_sname(SCMop,I)
      implicit none

! Description:
      ! Diagnostics may be defined on multiple sub-steps, in which
      ! case the different entries for the different sub-steps need to
      ! have the sub-step number appended to their short name before
      ! being output so that they can be distinguished.  This routine
      ! returns the short name of diagnostic entry with, if it's
      ! defined within a sub-stepped part of the code, its sub-step
      ! number appended to the end.
! Method:
! Owner: Luke Jones
! History:
! Version  Date      Comment
! =======  ====      =======
! 6.2      19/01/06  Original code (Luke Jones)
! Code Description:
      ! Language: Fortran77 with some bits of Fortran90

! SCMop_type is defined in here...
#include "s_scmoptype_defn.h"

      ! INOUT: The derived-type structure containing all the diagnostic
      ! information
      type(SCMop_type) :: SCMop
      ! IN: The entry in SCMop we want the (possibly altered) short
      ! name of.
      Integer :: I

      ! The output string should be somewhat longer than a normal sname
      ! to account for the appended substep number.
      Character (len=lsname+10) :: add_substep_to_sname
      Character (len=8) :: c_substep

      ! Add the sub-step number onto the diagnostic entry's short name
      ! if its sub-step number is greater than one, or it is equal to
      ! one but it is expected there will be multiple sub-steps
      If (SCMop%substep(I) > 1 .OR.                                     &
     &     (SCMop%substep(I) == 1 .AND. SCMop%num_substeps > 1)) Then
         write(c_substep,'(I8)') SCMop%substep(I)
         write(add_substep_to_sname,'(A)')                              &
     &      trim(SCMop%sname(I))//'_'//trim(adjustl(c_substep))
      Else
         ! This diagnostic is defined on only one substep, or is outside
         ! of the sub-stepped part(s) of the code. We'll leave its name
         ! as it is then.
         add_substep_to_sname=SCMop%sname(I)
      End If

      RETURN
      END FUNCTION add_substep_to_sname
#endif
