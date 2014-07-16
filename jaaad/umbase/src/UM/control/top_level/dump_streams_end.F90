#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Close SCM output files

      subroutine DUMP_STREAMS_END(SCMop)

      USE NetCDF
      implicit none

! Description:
      ! Close all open SCM output files
! Method:
      ! At the time of writing documentation is available from
      ! http://www-nwp/~frlj
! Owner: Luke Jones
! History:
! Version  Date      Comment
! =======  ====      =======
! 6.3      24/08/06  Original code (Luke Jones)
! Code Description:
      ! Language: Fortran77 with some bits of Fortran90

! SCMop_type is defined in here...
#include "s_scmoptype_defn.h"

      type(SCMop_type) :: SCMop ! INOUT The derived-type structure
                                ! containing all the diagnostic
                                ! information

! Parameters and stuff used around and about the internal
! workings of the diagnostic system.
#include "s_scmop_internal.h"

      integer :: N,Status

      ! Loop over each output stream
      do N=1,maxnstreams

         ! Is this stream switched on and are any diagnostics being sent
         ! to it?
         if (SCMop%strm(N)%switch /= 0.and.                             &
     &        SCMop%strm(N)%n_output > 0) then

            ! This stream was opened. Close it now.
            if (SCMop%strm(N)%format >= 0.and.                          &
     &           SCMop%strm(N)%format <= 3) then
               close(SCMop%strm(N)%op_unit)

            elseif (SCMop%strm(N)%format == 4) then
               ! Close the NetCDF file
               Status = Nf90_Close(int(SCMop%strm(N)%op_unit,incdf))

            else
               print*,'dump_streams_end ERROR: unknown format '//       &
     &              'for stream ',N
            endif

         endif
      enddo

      return
      end
#endif
