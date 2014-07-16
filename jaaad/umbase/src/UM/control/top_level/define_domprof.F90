#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Set up a domain profile for SCM diagnostics

      subroutine define_domprof(tindex,tname,                           &
     &     rowa1,rowa2,                                                 &
     &     rowb1,rowb2,                                                 &
     &     lev1,lev2,SCMop)
      implicit none

! Description:
      ! Set up a domain profile, i.e. a spatial region over which
      ! diagnostics may be defined. This will be stored in SCMop
      ! and in future can be referred to using tindex.
! Method:
      ! At the time of writing documentation is available from
      ! http://www-nwp/~frlj
! Owner: Luke Jones
! History:
! Version  Date      Comment
! =======  ====      =======
! 5.5      06/02/03  Original code (Luke Jones)
! 6.0      05/09/03  Adapt to changes in contents of include files.
!                    Luke Jones.
! 6.0      29/10/03  Added check that tindex does not exceed
!                    maxndomprof. Luke Jones.
! Code Description:
      ! Language: Fortran77 with some bits of Fortran90

      integer tindex            ! IN An index for the profile
      character (len=*) tname   ! IN A name for the profile
      integer rowa1,rowa2,    & ! IN The horizontal area over which
     &     rowb1,rowb2          ! diagnostics of this type are defined
      integer lev1,lev2         ! IN The vertical range over which " " "

! SCMop_type is defined in here...
#include "s_scmoptype_defn.h"

      type(SCMop_type) :: SCMop ! INOUT The derived-type structure
                                ! containing all the diagnostic
                                ! information

! Parameters and stuff used around and about the internal
! workings of the diagnostic system
#include "s_scmop_internal.h"

      ! Check that tindex is not greater than maxndomprof (i.e.
      ! avoid array out-of-bounds errors)
      if (tindex >  maxndomprof) then
         print*,'define_domprof ERROR: trying to define a domain'
         print*,'profile with an index above the maximum. Increase'
         print*,'maxndomprof or use a lower index. Domain profile'
         print*,'not defined:',trim(tname)
         goto 999
      endif

      ! For safety, enforce a rule that no diagnostic may be defined
      ! unless *every* domain has been defined. This is to ensure
      ! SCMop%nelements(x) can never become incorrect by a domain being
      ! re-defined.
      if (SCMop%nentries >  0) then
         print*,'define_domprof ERROR: nentries is non-zero, setting ', &
     &        'to zero now. Some diagnostics may be discarded',         &
     &        SCMop%nentries
         SCMop%nentries=0
      endif

      ! Fill the relevant section of the SCMop structure
      SCMop%d_name(tindex)=tname
      SCMop%d_rowa1(tindex)=rowa1
      SCMop%d_rowa2(tindex)=rowa2
      SCMop%d_rowb1(tindex)=rowb1
      SCMop%d_rowb2(tindex)=rowb2
      SCMop%d_lev1(tindex)=lev1
      SCMop%d_lev2(tindex)=lev2

 999  continue
      return
      END SUBROUTINE define_domprof
#endif
