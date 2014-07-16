#if defined(SCMA)
! Start of include file: s_scmop_internal.h
! Description:
!  Defines/declares some parameters and statement functions used
!  around and about the internal workings of the SCM diagnostic system.
! Current Code Owner: Luke Jones
!
! History:
! Version  Date      Comment
! =======  ====      =======
! 6.0      29/08/03  Original code (actually, mostly ripped out of
!                    s_scmop.h). Luke Jones

! Include the parameters that are used in the calls
! to SCMoutput...
#include "s_scmop.h"

      ! Statement function to translate daycount and stepcount into
      ! an integer representing the number of timesteps since the
      ! start of the run
      integer stepnmbr
      stepnmbr(SCMop)=(SCMop%daycount-1)*SCMop%full_daysteps+           &
     &     SCMop%stepcount

!-------- Stream inquiry/manipulation statement functions -----------

      ! These need to be declared for the statement functions below
      integer(i64) streamlist

      ! Statement function to modify an encoded stream list
      ! to flag that the streams should not be written to file
      integer DoNotWrite
      DoNotWrite(streamlist)=streamlist+2**(inot_written-1)

      ! Statement function to test whether a particular stream is
      ! switched on in an encoded stream list
      logical StreamIsOn
      StreamIsOn(streamlist,strm)=                                      &
     &     (streamlist-int(streamlist/Stream(strm+1))*Stream(strm+1))/  &
     &     Stream(strm).ge.1

      ! Statement function to test whether the "do not write" flag
      ! has been switched on in an encoded stream list
      logical NotWritten
      NotWritten(streamlist)=StreamIsOn(streamlist,inot_written)

      ! Statement function to test whether any stream is switched
      ! on in an encoded stream list
      logical AnyStreamOn
      AnyStreamOn(streamlist)=                                          &
     &     (.not.NotWritten(streamlist).and.streamlist.ne.0)            &
     &     .or.                                                         &
     &     (NotWritten(streamlist).and.streamlist.ne.DoNotWrite(0))

! End of include file: s_scmop_internal.h
#endif
