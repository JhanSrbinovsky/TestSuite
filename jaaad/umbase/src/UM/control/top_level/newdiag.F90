#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Create a new SCM diagnostic entry in SCMop.

      subroutine newdiag(sname,lname,units,timprof,domprof,             &
     &     istrm,lnot_written,sname2,d,SCMop)
      implicit none

! Description:
      ! Create a new diagnostic entry in SCMop. Returns index of
      ! newly created diagnostic, which is unchanged from its input 
      ! value if an error ocurrs and the entry could not be created.
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
! 6.0      02/09/03  Set SCMop%ncols, %nrows and %nlevs. Luke Jones.
! 6.2      19/01/06  Enabled outputting of SCM diagnostics on multiple
!                    sub-steps.                            Luke Jones.
! Code Description:
      ! Language: Fortran77 with some bits of Fortran90

! SCMop_type is defined in here...
#include "s_scmoptype_defn.h"

      type(SCMop_type) :: SCMop ! INOUT The derived-type structure
                                ! containing all the diagnostic
                                ! information
      character (len=lsname) sname ! IN Short name for the diagnostic
      character (len=*) lname      ! IN Long name for the diagnostic
      character (len=*) units      ! IN Units of the diagnostic
      character (len=*) sname2     ! IN Short name of a previously
      ! defined diagnostic used in the construction of this one
      integer timprof,   & ! IN Time profile for the diagnostic
     &     domprof,      & ! IN Domain profile for the diagnostic
     &     istrm           ! IN Stream to which diagnostic is to be sent
      logical lnot_written ! IN If true, diagnostic will not be
                           ! written out
      integer d            ! IN/OUT Index of the newly created
                           ! diagnostic entry in SCMop. Unchanged from
                           ! input value if entry could not be created.

! Parameters and stuff used around and about the internal
! workings of the diagnostic system
#include "s_scmop_internal.h"

      integer d_input,i,d1,d2,d3
      Logical                                                           &
     &  sname2_found                                                    &
     &, right_dumping_period

      d_input=d ! record input value of d

      if (SCMop%nentries == SCMop%maxnentries) then
         ! No room at the inn, make the inn bigger.
! DEPENDS ON: expand_scmop
         call expand_SCMop(SCMop)
      endif
      ! Increment the recorded number of entries
      SCMop%nentries=SCMop%nentries+1
      ! d is the index of this new entry
      d=SCMop%nentries

      ! Record the details of this entry in the respective arrays...
      SCMop%sname(d)=sname
      SCMop%lname(d)=lname
      SCMop%units(d)=units
      SCMop%domprof(d)=domprof

      if (timprof <  only_radsteps) then
         ! This is a normal diagnostic - calculated on every timestep
         SCMop%timprof(d)=timprof
         SCMop%only_radsteps(d)=.false.
      else
         ! This diagnostic is based on an array only valid on
         ! radiation timesteps
         SCMop%timprof(d)=timprof-only_radsteps
         SCMop%only_radsteps(d)=.true.
      endif

      SCMop%streams(d)=Stream(istrm)
      if (lnot_written) then
         SCMop%streams(d)=DoNotWrite(SCMop%streams(d))
      endif

      ! The dumping period of the diagnostic is that of the stream
      ! it's being sent to
      SCMop%dump_step(d)=SCMop%strm(istrm)%dump_step
      SCMop%nadd2dump(d)=SCMop%strm(istrm)%dump_step

      ! The dimensions of the diagnostic array can be obtained from
      ! the domain profile
      d1=SCMop%d_rowa2(domprof)-SCMop%d_rowa1(domprof)+1
      d2=SCMop%d_rowb2(domprof)-SCMop%d_rowb1(domprof)+1
      d3=SCMop%d_lev2 (domprof)-SCMop%d_lev1 (domprof)+1
      SCMop%ncols(d)=d1
      SCMop%nrows(d)=d2
      SCMop%nlevs(d)=d3
      SCMop%nelements(d)=d1*d2*d3

      ! Allocate the space for the dump array
      allocate(SCMop%diag(d)%dump(d1*d2*d3))
      ! Initialise it for initialisation's sake
      SCMop%diag(d)%dump=-999.

      ! lastencounter will be set in SCMoutput
      SCMop%lastencounter(d)=-1

      ! Set the substep number that this entry is being created for
      SCMop%substep(d)=SCMop%substep_number

      ! wd will be the index of a diagnostic entry upon which this
      ! diagnostic depends (set to zero if sname2 is a null string)
      SCMop%wd(d)=0
      if (len(sname2) >  0) then
         ! Find the index of the weighting diagnostic
         sname2_found=.false. ! (flags if at least found the right name)
         right_dumping_period=.FALSE. ! (" " " " right dumping period)

         do i=1,d-1
            if (SCMop%sname(i) == sname2) then
               ! We have found a diagnostic of the correct name
               sname2_found=.true.
               ! But it must have the same dumping period or be constant
               if (SCMop%dump_step(i) == SCMop%strm(istrm)%dump_step    &
     &              .or.                                                &
     &              SCMop%timprof(i) == t_const) then

                  right_dumping_period=.TRUE.

                  ! But it must also be defined on the same substep as
                  ! the current entry, or have been defined outside of
                  ! a sub-stepped part of the model.

                  If (SCMop%substep(i) == SCMop%substep_number .OR.     &
     &                SCMop%substep(i) == 0) Then

                     ! We have found an entry with the correct name &
                     ! substep
                     SCMop%wd(d)=i
                     exit

                  End If
               endif
            endif
         enddo
         if (SCMop%wd(d) == 0) then
            write(*,*)' '
            write(*,*)'************ERROR IN ROUTINE NEWDIAG************'
            write(*,*)'* The following error stems from a call to '//   &
     &           'SCMoutput...'
            write(*,*)'* You have requested that diagnostic "',         &
     &           trim(sname),'" be dependent '
            write(*,*)'* on diagnostic "',trim(sname2),                 &
     &           '" (which is non-constant),'
            if (.not.sname2_found) then
               write(*,*)'* but the latter diagnostic has not yet '//   &
     &              'been defined.'
               write(*,*)'* Diagnostic "',trim(sname),'" will '//       &
     &              'therefore not be calculated or output.'
               write(*,*)'* Note: this message will be '//              &
     &              'repeated for every stream you have '
               write(*,*)'* requested for this diagnostic'
            Else If (.NOT. right_dumping_period) Then
               write(*,*)'* but you have requested that the former '//  &
     &              'be sent to a stream with a'
               write(*,*)'* dumping period of ',                        &
     &              SCMop%strm(istrm)%dump_step,                        &
     &              ', while the latter is not. This is not'
               write(*,*)'* permitted: if diagnostic "A" is to '//      &
     &              'depend on diagnostic "B", and'
               write(*,*)'* diagnostic "B" is not constant, then '//    &
     &              'diagnostic "B" must be calculated'
               write(*,*)'* for every dumping period that is '//        &
     &              'requested for diagnostic "A". This'
               write(*,*)'* can be guaranteed by sending "B" to all '// &
     &              'the streams which you have'
               write(*,*)'* requested for "A".'
               write(*,*)'* Diagnostic "',trim(sname),'" will '//       &
     &              'therefore not be calculated or'
               write(*,*)'* output with a dumping period of ',          &
     &              SCMop%strm(istrm)%dump_step
            Else If (SCMop%substep_number /= 0) Then
               write(6,*)'* but while that diagnostic appears to be ' //&
     &                     'defined within a sub-stepped'
               write(6,*)'* part of the code, no entry could be '     //&
     &                     'found for it on the current substep '
               write(6,*)'* of ',SCMop%substep_number
               write(*,*)'**************************************'//     &
     &              '**********'
               write(*,*)' '
            Else
               write(6,*)'* but the latter is defined within a '      //&
     &                     'sub-stepped part of the '
               write(6,*)'* code while the former is not. This does ' //&
     &                     'not make sense.'
            endif
            write(6,*)'**************************************'        //&
     &                '**********'
            write(6,*)' '
            d=d_input
            SCMop%nentries=SCMop%nentries-1
            goto 999
         endif
      else
         if (timprof == t_div.or.timprof == t_mult.or.                  &
     &        timprof == t_acc_div.or.                                  &
     &        timprof == t_acc_mult) then
            write(*,'(A,1X,A,1X,A,1X,I3)')                              &
     &           'newdiag ERROR: you have requested this '//            &
     &           'diag to be dependent on another but have not '//      &
     &           'specified which.',trim(sname),timprof
            d=d_input
            SCMop%nentries=SCMop%nentries-1
            goto 999
         endif
      endif

      ! Assign this diagnostic an integer unique to its sname
      ! (do this by searching for a previously defined diagnostic
      ! with the same sname, if you find it give it that number,
      ! if you don't give it the highest number you came across
      ! plus one)
      SCMop%sname_id(d)=0
      do i=1,SCMop%nentries-1
         if (SCMop%sname(i) == sname) then
            ! i has the same sname as d, give it the same sname_id
            SCMop%sname_id(d)=SCMop%sname_id(i)
            exit
         else
            ! record the largest sname_id, but store as -ve to
            ! indicate that we have not found a match for sname
            SCMop%sname_id(d)=                                          &
     &           min(SCMop%sname_id(d),-SCMop%sname_id(i))
         endif
      enddo
      if (SCMop%sname_id(d) <= 0) then
         ! A diagnostic with the same sname was not found: assign a
         ! value one larger than the current largest value of sname_id
         SCMop%sname_id(d)=-SCMop%sname_id(d)+1
         ! Record the number of diagnostics that will be output (with
         ! no double counting from the same diagnostic being
         ! calculated with different dumping periods)
         if (.not.lnot_written) SCMop%n_output=SCMop%n_output+1
      endif
      ! Record the number of diagnostics that will be sent to
      ! this stream
      if (.not.lnot_written) SCMop%strm(istrm)%n_output=                &
     &     SCMop%strm(istrm)%n_output+1

 999  continue

      return
      END SUBROUTINE newdiag
#endif
