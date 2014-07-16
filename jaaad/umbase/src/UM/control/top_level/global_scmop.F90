#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ SCM module for making SCMop accessible to SCMoutput from SCM_main.

      module global_SCMop
      implicit none
!
! Description:
!   All the SCM diagnostic information is stored in a single
!   structure, SCMop, first declared in SCM_main. This module acts
!   like a common block to allow both SCM_main and SCMoutput access
!   to this structure without it having to be passed down the call tree.
!   A module has to be used because a derived-type structure cannot be
!   put in a common.
! Method:
!   At the time of writing documentation about the whole SCM diagnostic
!   system is available from http://www-nwp/~frlj
! Owner: Luke Jones
! History:
! Version  Date      Comment
! =======  ====      =======
! 6.0      05/09/03  Original code. Luke Jones.
! Code Description:
      ! Language: Fixed form Fortran 90.

! The definition of SCMop_type...
#include "s_scmoptype_defn.h"

      ! Declare SCMop
      type(SCMop_type) :: SCMop

      end module global_SCMop
#endif
