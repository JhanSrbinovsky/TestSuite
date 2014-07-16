#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Obtain indices of diagnostic entries in SCMop given SCMoutput inputs

      subroutine getdistinctdiags(sname_i,lname,units,timprof,domprof,  &
     &     streams,sname2_i,                                            &
     &     ndistinct_diags,distinct_diags,SCMop)
      implicit none

! Description:
      ! On the first call to this routine with a given sname_i, create
      ! the corresponding diagnostic entries in SCMop and return their
      ! indices in distinct_diags(1:ndistinct_diags). On subsequent
      ! calls, simply look up the previously created entries and return
      ! in distinct_diags(1:ndistinct_diags).
! Method:
      ! For all streams to which the diagnostic is to be sent, look
      ! for a corresponding entry in SCMop with the correct dump_step.
      ! If one does not exist create it by calling NEWDIAG. Record its
      ! index if it's different to the index found in any previous
      ! cycle of the loop over streams. Thus build up a list of the
      ! entries in SCMop to which this diagnostic is associated, and
      ! return.
      ! At the time of writing documentation is available from
      ! http://www-nwp/~frlj
! Owner: Luke Jones
! History:
! Version  Date      Comment
! =======  ====      =======
! 5.5      06/02/03  Original code (Luke Jones)
! 6.0      05/09/03  Adapt to changes in contents of include files.
!                    Luke Jones.
! 6.0      02/09/03  Added new functionality controlled by namelist
!                    DIAGS. Luke Jones.
! 6.2      19/01/06  Enabled outputting of SCM diagnostics on multiple
!                    sub-steps.                            Luke Jones.
!
! Code Description:
      ! Language: Fortran77 with some bits of Fortran90

! SCMop_type is defined in here...
#include "s_scmoptype_defn.h"

      ! See SCMoutput for a description of the input variables
      character(len=*) :: sname_i,lname,units,sname2_i ! IN
      integer timprof,domprof,streams                  ! IN
      integer ndistinct_diags,distinct_diags(maxnstreams) ! OUT
      type(SCMop_type) :: SCMop ! INOUT The derived-type structure
                                ! containing all the diagnostic
                                ! information

! Parameters and stuff used around and about the internal
! workings of the diagnostic system...
#include "s_scmop_internal.h"

      ! d will equal this if the entry could not be created
      integer, parameter :: unset=0

      integer i,j,d,initial_nentries
      character(len=lsname) :: sname,sname2

      logical distinct

      ndistinct_diags=0

      ! Elsewhere we are going to rely on the "sname" for any
      ! diagnostic being exactly lsname characters long, so here we
      ! make them the right length.
      sname=sname_i
      sname2=sname2_i
      if (len(sname_i) >  lsname.and.SCMop%first_pass) then
         ! Warn the user that we have done this
         print*,'getdistinctdiags WARNING: diagnostic name '//          &
     &        'truncated: ',sname_i(1:len(sname_i)),' -> ',             &
     &        sname
      endif

      ! If this is the first timestep, check this "sname" has not been
      ! used in a previous call to SCMoutput
      if (SCMop%first_pass) then
         do i=1,SCMop%nentries
            If (SCMop%sname(i) == sname .AND.                           &
     &           SCMop%substep(i) == SCMop%substep_number) Then
               print*,'getdistinctdiags ERROR: same sname used in >1 '//&
     &              'call to SCMoutput on same sub-step:',sname
               goto 999
            endif
         enddo
      endif

      ! Get list of streams to which the diagnostic is to be
      ! sent. Nominally this is determined by the variable, streams
      ! (the so-called "hard-wired" list of streams). But namelist
      ! information can alter this by choosing to ignore the
      ! hard-wired list and/or to add/remove diagnostics to/from
      ! the list. The two lists of extra diagnostics to be added to,
      ! and rejected from, stream X are in SCMop%strm(X)%accept_list and
      ! SCMop%strm(X)%reject_list respectively.

      streamlist=0 ! This will be the final list of streams encoded into
                   ! one integer. For now it is set to "no streams".

      ! Loop over all streams.
      do j=1,maxnstreams

         ! If this stream is closed do nothing.
         if (SCMop%strm(j)%switch == 0) then
            CYCLE ! Go to next value of j.
         endif

         ! Shall we pay attention to the list of streams as
         ! specified in the call to SCMoutput?
         if (SCMop%strm(j)%heed_hardwired /= 0) then
            ! Yes. Add stream j if requested in the call.
            if (StreamIsOn(streams,j).and.                              &
     &           .not.StreamIsOn(streamlist,j)) then
               streamlist=streamlist+Stream(j)
            endif
            ! If it has been requested in the call that this
            ! diagnostic should not be written out, make sure this
            ! is so.
            if (NotWritten(streams).and..not.NotWritten(streamlist))    &
     &           then
               streamlist=DoNotWrite(streamlist)
            endif
         endif

         ! Shall we pay attention to the list of diagnostics
         ! requested by namelist to be sent to this stream?
         if (SCMop%strm(j)%heed_acceptlist /= 0) then
           ! Yes.

           ! If the diagnostic is not already being sent to this stream
            if (.not.StreamIsOn(streamlist,j)) then
               ! then check if it's in the list of extra diagnostics
               ! to suck into this stream.

               ! Loop over the the names of diagnostics to be accepted
               do i=1,size(SCMop%strm(j)%accept_list)
                  ! Do we have a name match?
                  if (trim(SCMop%strm(j)%accept_list(i)) ==             &
     &                 trim(sname)) then
                     ! The name was in the list for stream j. Send
                     ! this diagnostic to stream j.
                     streamlist=streamlist+Stream(j)
                     EXIT   ! this inner loop
                  endif
               enddo

            endif

         endif ! SCMop%strm(j)%heed_acceptlist /= 0

         ! Shall we pay attention to list of diagnostics requested
         ! by namelist to be prevented from going to this stream?
         if (SCMop%strm(j)%heed_rejectlist /= 0) then
            ! Yes.

            ! If the diagnostic is being sent to this stream
            if (StreamIsOn(streamlist,j)) then
               ! then check if it's in the list of diagnostics
               ! which are not to be sent to this stream

               do i=1,size(SCMop%strm(j)%reject_list)
                  if (trim(SCMop%strm(j)%reject_list(i)) ==             &
     &                 trim(sname)) then
                     ! This diagnostic is not to be sent to this stream.
                     streamlist=streamlist-Stream(j)
                     EXIT
                  endif
               enddo

            endif ! StreamIsOn(streamlist,j)

         endif    ! SCMop%strm(j)%heed_rejectlist /= 0

      enddo       ! j=1,maxnstreams

      ! From now on we will use streamlist instead of streams as the
      ! integer representing the list of streams to which the
      ! diagnostic is to be sent.

      initial_nentries=SCMop%nentries
      ! Loop over all output streams
      do j=1,maxnstreams
         ! j represents a stream, is the diagnostic to go
         ! to this stream?
         if (StreamIsOn(streamlist,j).and.                              &
     &        SCMop%strm(j)%switch /= 0) then
            ! Yes, look to see if the entry already exists in SCMop
            d=unset ! (label d as unset for now)
            do i=1,SCMop%nentries
               if (SCMop%sname(i) == sname.and.                         &
     &              SCMop%substep(i) == SCMop%substep_number .AND.      &
     &              SCMop%dump_step(i) == SCMop%strm(j)%dump_step) then
                    ! This diagnostic exists.
                  d=i
                  exit
               endif
            enddo
            if (d == unset) then
               if (SCMop%first_pass) then
                  ! We didn't find an existing diagnostic fitting the
                  ! inputs, create a new diagnostic entry in SCMop
! DEPENDS ON: newdiag
                  call newdiag(sname,lname,units,timprof,domprof,j,     &
     &                 NotWritten(streamlist),                          &
     &                 sname2(1:lsname*min(len(sname2_i),1)),d,SCMop)
               else
                  ! Should not be having new diagnostics beyond
                  ! the first timestep
                  print*,'getdistinctdiags ERROR: new diag after '//    &
     &                 'first timestep:',                               &
     &                 sname(1:len(sname)),                             &
     &                 SCMop%stepcount,SCMop%daycount,d
                  goto 999
               endif
            else
               if (SCMop%first_pass) then
                  ! d is the index of a diagnostic with the same sname
                  ! as given in the input parameter list and the same
                  ! dump_step as stream j. Send diagnostic d to stream
                  ! j as well as it's existing streams then.
                  if (.not.StreamIsOn(SCMop%streams(d),j)) then
                     SCMop%streams(d)=SCMop%streams(d)+Stream(j)
                     ! Increment n_output to count the number we're
                     ! going to output to this stream
                     if (.not.NotWritten(streamlist)) then
                        SCMop%strm(j)%n_output=                         &
     &                       SCMop%strm(j)%n_output+1
                     endif
                  else
                     ! Diagnostic d is already going to this stream,
                     ! this should not ocurr.
                     print*,'getdistinctdiags ERROR: Same diag. sent '//&
     &                    'to same stream twice',d,sname(1:len(sname)), &
     &                    j,SCMop%streams(d)
                  endif
               else
                  ! Check that diagnostic d is set up to go to
                  ! this stream
                  if (.not.StreamIsOn(SCMop%streams(d),j)) then
                     print*,'getdistinctdiags ERROR: the requested '//  &
     &                    'streams for this diagnostic have changed',   &
     &                    d,sname(1:len(sname)),j,SCMop%streams(d)
                  endif
               endif
            endif

            ! We should now have a value for d, but it may still be
            ! unset if an error ocurred in newdiag
            if (d /= unset) then
               ! Is this diagnostic distinct? i.e. is the value of d at
               ! this point different than for any previous value of j
               ! in this loop?
               distinct=.true.
               do i=1,ndistinct_diags
                  if (distinct_diags(i) == d) then
                     distinct=.false.
                     exit
                  endif
               enddo
               if (distinct) then
                  ndistinct_diags=ndistinct_diags+1
                  distinct_diags(ndistinct_diags)=d
                  ! Make a requirement that diagnostic entries must be
                  ! created in order (see use of diag_mem in SCMoutput)
                  if (ndistinct_diags >  1) then
                     if (d /= distinct_diags(ndistinct_diags-1)+1) then
                        print*,'getdistinctdiags ERROR: '               &
     &                       //'non-consecutive diags',                 &
     &                       distinct_diags(1:ndistinct_diags)
                     endif
                  endif
               endif
            endif

         endif                  ! (StreamIsOn(streams,j))
      enddo                     ! j=1,maxnstreams

 999  continue
      return

      END SUBROUTINE getdistinctdiags
#endif
