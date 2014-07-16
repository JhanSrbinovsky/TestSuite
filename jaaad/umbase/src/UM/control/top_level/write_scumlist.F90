#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Write out info about the output diagnostics for one or all streams

      subroutine write_scumlist(SCMop,istrm)
      implicit none

! Description:
      ! For a given stream, writes information about each domain
      ! profile and diagnostic that is active in that stream to the
      ! stream's output file as part of its header. PV-wave can then
      ! read this info to make sense of the rest of the data file. If
      ! the requested stream number is zero, a new file will be opened
      ! and all the same information will be written to it, but now
      ! pertaining to all streams. Such a file will list all domain
      ! profiles and diagnostics, and will have additional comments.
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
! 6.0      02/09/03  Previously, routine only wrote data to a new file
!                    called scumlist.dat. Now has option of writing to
!                    a new file -or- to the existing output file of an
!                    output stream. Luke Jones.
! 6.2      19/01/06  Enabled outputting of SCM diagnostics on multiple
!                    sub-steps.                            Luke Jones.
! Code Description:
      ! Language: Fortran77 with some bits of Fortran90

! SCMop_type is defined in here...
#include "s_scmoptype_defn.h"

      type(SCMop_type) :: SCMop ! IN The derived-type structure
                                ! containing all the diagnostic
                                ! information
      integer istrm  ! IN The stream in question. Zero=all streams.

      integer ndiags ! The number of diagnostic entries found in
                     ! SCMop which we will write about.
      integer diags(SCMop%nentries) ! Contains in 1:niags the indices
                     ! of the entries in SCMop
      integer domprof_sparse(maxndomprof) ! A sparse array indicating
                     ! the domain profiles possessed by the ndiags
                     ! entries.
      integer i,j    ! Counters

      ! Character function that incorporates the substep number into
      ! the short name of a diagnostic. Somewhat longer than a normal
      ! short name.
      Character (len=lsname+10) :: add_substep_to_sname,short_name

      ! The longest possible sname after the substep has been appended
      Integer :: sname_length

! Parameters and stuff used around and about the internal
! workings of the diagnostic system
#include "s_scmop_internal.h"

      integer unit ! The unit to which we'll write
      character (len=100) fmt
      character (len=lsname) sname
      Integer substep

      ! Make a list the diagnostics we're going to write and
      ! which domain profiles they use. ndiags will hold the total
      ! number of diagnostic entries in SCMop fitting the bill, and
      ! domprof_sparse will be a sparse array indicating which domain
      ! profiles are used by these diagnostics.
      ndiags=0
      domprof_sparse=0
      ! If a specific stream has been specified...
      if (istrm /= 0) then
         ! ... select only those diagnostics going to that stream
         do i=1,SCMop%nentries
            if (StreamIsOn(SCMop%streams(i),istrm).and..not.            &
     &           NotWritten(SCMop%streams(i))) then
               ndiags=ndiags+1
               diags(ndiags)=i
               domprof_sparse(SCMop%domprof(i))=1
            endif
         enddo
         ! Make a check just for the sake of it
         if (ndiags /= SCMop%strm(istrm)%n_output) then
            print*,'dump_streams_init ERROR: an inconsistency '//       &
     &           'has ocurred !',ndiags,SCMop%strm(istrm)%n_output,     &
     &           istrm,SCMop%nentries
            ! Switch the diagnostic system off so dodgy data is
            ! not mistakenly used in good faith.
            print*,' -> switching diagnostic system off!'
            SCMop%on=.false.
         endif
      else
         ! If strm=0 then use all diagnostics which are going to
         ! at least one stream. But since the same diagnostic can
         ! go to several streams, and thus have several entries in
         ! SCMop, then make a condition that if the i'th entry has
         ! the same short name as the (i-1)'th, then skip the i'th
         ! to avoid multiple identical lines. We -do- want multiple
         ! lines for different sub-steps though, so allow identical
         ! consecutive snames if their sub-step numbers are different.

         sname=''
         substep=-999

         do i=1,SCMop%nentries
            if (AnyStreamOn(SCMop%streams(i)).and.                      &
     &           (SCMop%sname(i) /= sname .OR.                          &
     &            SCMop%substep(i) /= substep)) Then

               ndiags=ndiags+1
               diags(ndiags)=i
               domprof_sparse(SCMop%domprof(i))=1
               sname=SCMop%sname(i)
               substep=SCMop%substep(i)

            endif
         enddo
      endif

      ! Determine the unit number to write to
      if (istrm /= 0) then
         ! Write to the unit of the specified stream. The file
         ! should already be open.
         unit=SCMop%strm(istrm)%op_unit
      else
         ! Write to a new file we'll attach to unit 10
         unit=10
         open (unit=unit,file='scumlist.dat')
      endif

!-------------------------------------------------------------
!     Write about the domain profiles
!-------------------------------------------------------------

      ! Write the number of domain profiles used by
      ! diagnostics in this stream
      if (istrm == 0) write(unit,*)'No. of domain profiles:'
                      write(unit,'(I3)') sum(domprof_sparse)

      ! Write the format with which we will write some of
      ! the upcoming lines
      fmt='(I3,1X,A15,I3,1X,I3,1X,I3,1X,I3,1X,I3,1X,I3)'
      if (istrm == 0) write(unit,'(A)')'Line format:'
                      write(unit,'(A)') trim(fmt)

      ! Write info about each domain profile
      if (istrm == 0) write(unit,*)'List of domain profiles:'
      do i=1,maxndomprof
         if (domprof_sparse(i) == 1) then
            write(unit,fmt)                                             &
     &           i,SCMop%d_name(i),                                     &
     &           SCMop%d_rowa1(i),SCMop%d_rowa2(i),                     &
     &           SCMop%d_rowb1(i),SCMop%d_rowb2(i),                     &
     &           SCMop%d_lev1(i),SCMop%d_lev2(i)
         endif
      enddo

!-------------------------------------------------------------
!     Write about the diagnostics themselves
!-------------------------------------------------------------

      ! Write the number
      if (istrm == 0) write(unit,*)'No. of diagnostics:'
                      write(unit,'(I3)') ndiags

      ! If there are multiple sub-steps then the short-name will
      ! have "_#N" tagged onto it, where N is the sub-step number.
      ! In this case we need to increase the amount of space made
      ! available for the short name.

      If (SCMop%num_substeps <= 1) Then
        sname_length=lsname
      Else
        sname_length=lsname+2+(int(log10(float(SCMop%num_substeps)))+1)
      End If

      ! Compose the format of the line that will describe
      ! each diagnostic

      write(fmt,'(A,I2, A,I2, A,I2, A)')'(I3' //                        &
     &     ',1X,A',sname_length,                                        &
     &     ',1X,A',llname,                                              &
     &     ',1X,A',lunits,                                              &
     &     ',I2)'
      ! Write the format to the file
      if (istrm == 0) write(unit,'(A)')'Line format:'
                      write(unit,'(A)') trim(fmt)

      ! Write a 1-line description of each diagnostic
      ! consisting of a unique integer ID, its short name,
      ! its long name, its units and its domain profile.
      if (istrm == 0) write(unit,*)'List of diagnostics '//             &
     &                '(i,sname,lname,units,domprof):'
      do i=1,ndiags
         j=diags(i)
! DEPENDS ON: add_substep_to_sname
         short_name=add_substep_to_sname(SCMop,j)
         write(unit,fmt)                                                &
     &        SCMop%sname_id(j),short_name,                             &
     &        SCMop%lname(j),SCMop%units(j),SCMop%domprof(j)
      enddo

      if (istrm == 0) then
         ! We need to close the new file we opened above
         close(unit)
      endif

      return
      END SUBROUTINE write_scumlist
#endif
