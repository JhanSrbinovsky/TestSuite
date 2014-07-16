#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Create a SCM output diagnostic

      subroutine SCMoutput(x, sname, lname, units, timprof, domprof     &
               , streams, sname2, calling_routine)

      USE global_SCMop
      USE scm_utils, Only: scm_trap_nan

      implicit none

! Description:
      ! Create output diagnostic based on given inputs. This routine
      ! cannot be called more than once per timestep with the same
      ! "sname" unless the call is inside a sub-stepped part of the 
      ! model and the sub-stepping has been delimited by calls to 
      ! scm_substep_start and scm_substepping_end. The order of calls 
      ! to SCMoutput should not change between timesteps.
! Method:
      ! At the time of writing documentation is available from
      ! http://www-nwp/~frlj

! Original Author: Luke Jones
! Current Owner  : SCM Code owner

! Code Description:
      ! Language: Fortran77 with some bits of Fortran90

      integer domprof           ! IN Domain profile for the diagnostic
      real x  &  ! IN Variable from which diagnostic will be constructed
     &    ((SCMop%d_rowa2(domprof)-SCMop%d_rowa1(domprof)+1)*           &
     &     (SCMop%d_rowb2(domprof)-SCMop%d_rowb1(domprof)+1)*           &
     &     (SCMop%d_lev2 (domprof)-SCMop%d_lev1 (domprof)+1))

      character(len=*) :: sname,& ! IN Short name for the diagnostic,
                                  ! this should be unique
     &     lname,               & ! IN Long name for the diagnostic
     &     units,               & ! IN Units of the diagnostic
     &     sname2                 ! IN Short name of another, previously
                                  ! defined diagnostic which will be used
                                  ! in the construction of this one
                                  ! according to the time profile

      character(len=*) :: calling_routine
                                  ! IN Routine that has called scmoutput

      integer timprof,          & ! IN The time profile for the diagnostic
     &     streams                ! IN An encoded integer specifying
                                  ! which output streams the diagnostic
                                  ! is to go to
! Parameters and stuff used around and about the internal
! workings of the diagnostic system
#include "s_scmop_internal.h"


      integer d,i,j               ! General use
      character(len=30) :: sdum0

      logical startperiod,endperiod ! Will be used to flag if we are at
                                ! the start or end of a dumping period

      integer ndistinct_diags,distinct_diags(maxnstreams),              &
     &     nThroughPeriod,ntrad,ntrad1
      logical call_add2dump,order_changing

      ! Will hold the contents of SCMop%diag_mem while being
      ! re-allocated
      integer, allocatable :: itemp(:,:)

      ! Perform no action if SCMop is not turned on, except in the
      ! case of a constant diagnostic
      if (.not.SCMop%on.and.timprof /= t_const) goto 999

      ! Enforce a rule that sname must be at least one character long
      if (len(sname) == 0) then
         if (SCMop%first_pass) then
            print*,'SCMoutput ERROR: sname is a null string, this '//   &
     &           'is not allowed, diagnostic ignored: ',                &
     &           lname(1:len(lname))
         endif
         goto 999
      endif

#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      ! If the array is only valid on radiation time-steps then we
      ! need to know which time-steps those are. That code has not yet
      ! been implemented for the 3C and 3Z radiation schemes.
      if (timprof > only_radsteps)
         write(*,*)' '
         write(*,*)'**********WARNING FROM ROUTINE SCMOUTPUT**********'
         write(*,*)'At present the code does not exist inside the '
         write(*,*)'diagnostic system to output variables on radiation '
         write(*,*)'time-steps when using the 3C or 3Z radiation '
         write(*,*)'schemes. This diagnostic will be ignored: ',sname
         write(*,*)'************************************************'
         goto 999
      endif
#endif

      ! Increment recorded no. of calls to this routine this timestep
      if (SCMop%on) SCMop%nSCMoutput=SCMop%nSCMoutput+1

      ! If the requested list of streams span a range of dumping
      ! periods then the given inputs will correspond to more than one
      ! diagnostic entry in SCMop. Thus we need to get
      ! distinct_diags(1:ndistinct_diags) - the indices of those
      ! entries. The routine getdistinctdiags can do this but, if this
      ! is not the first timestep, and assuming the order of the calls
      ! to this routine doesn't change between timesteps, we can just
      ! use our memory of a previous timestep instead.

      order_changing=.false.
      if (.not.SCMop%first_pass) then
         ! Use memory of previous timesteps to know which diagnostic
         ! entries this call pertains to.
         ndistinct_diags=SCMop%diag_mem(SCMop%nSCMoutput,1)
         do i=1,ndistinct_diags
            distinct_diags(i)=SCMop%diag_mem(SCMop%nSCMoutput,2)+i-1

            ! Check we're right...
            if (trim(sname(1:min(lsname,len(sname)))) /=                &
     &           trim(SCMop%sname(distinct_diags(i))))                  &
     &           then
               write(6,*)'*******************************************'
               write(*,*) 'SCMoutput warning: the order of the calls'// &
     &              ' to SCMoutput seems to be changing:'
               write(*,'(A,I4,A)') '  On step ',stepnmbr(SCMop),        &
     &            ' expected '//                                        &
     &            trim(SCMop%sname(distinct_diags(i)))//                &
     &            ', but encountered '//sname(1:min(lsname,len(sname)))
               write(6,*)'*******************************************'
               order_changing=.TRUE.

            Else If (SCMop%substep_number /=                            &
     &               SCMop%substep(distinct_diags(i))) Then
               write(6,*)'*******************************************'
               write(6,'(A)') ' SCMoutput error: the order of the '  // &
     &                        'sub-steps seems to be changing.'
               write(6,'(A,I5,A,I3,A,I3)')                              &
     &              ' On step ',stepnmbr(SCMop),' diagnostic '       // &
     &                sname(1:min(lsname,len(sname)))                // &
     &              ' expected sub-step ',                              &
     &                SCMop%substep(distinct_diags(i)),                 &
     &              ' but got sub-step ',SCMop%substep_number
               write(6,'(A)') ' This does not make sense and is a '  // &
     &                        'sign of a potentially serious problem.'
               write(6,*)'*******************************************'
               order_changing=.true.
            endif
         enddo
      endif

      if (SCMop%first_pass.or.order_changing) then
         ! Either this is the first timestep or the order of the calls
         ! to SCMoutput is changing
! DEPENDS ON: getdistinctdiags
         call getdistinctdiags(sname,lname,units,timprof,domprof, & ! IN
     &        streams,sname2,                                     & ! IN
     &        ndistinct_diags,distinct_diags,SCMop)      ! OUT,OUT,INOUT
         ! Don't want to do this next bit for constant diagnostics
         ! (which should be declared when the system is off)
         if (SCMop%on) then
            ! Store the values of ndistinct_diags and distinct_diags(1)
            ! for future reference (to avoid unnecessary calls to
            ! getdistinctdiags). But is there enough space in the
            ! SCMop%diag_mem array?
            if (size(SCMop%diag_mem,1) == SCMop%nSCMoutput-1) then
               ! No. Make the array bigger...
               ! Allocate a temporary array with the same size and
               ! shape as diag_mem
               allocate(itemp(SCMop%nSCMoutput-1,2))
               ! Copy the contents of diag_mem into it
               itemp=SCMop%diag_mem
               ! Re-allocate diag_mem with a larger size
               deallocate(SCMop%diag_mem)
               allocate(SCMop%diag_mem(SCMop%nSCMoutput+49,2))
               ! Copy the original contents back in
               SCMop%diag_mem(1:SCMop%nSCMoutput-1,1:2)=itemp
            endif
            SCMop%diag_mem(SCMop%nSCMoutput,1)=ndistinct_diags
            SCMop%diag_mem(SCMop%nSCMoutput,2)=distinct_diags(1)
         endif
      endif

      ! From here on in, none of the input parameters are
      ! referred to at all, they have been distilled to
      ! distinct_diags(1:ndistinct_diags) and information
      ! in SCMop

      ntrad=SCMop%ntrad        ! No. of timesteps between calls to rad'n
      ntrad1=SCMop%ntrad1      ! Timestep containing 1st call to rad'n

      do i=1,ndistinct_diags
         d=distinct_diags(i)

         ! If this is not the first time we've seen this diagnostic,
         ! check the last time was the previous timestep
         if (SCMop%lastencounter(d) >= 0.and.                           &
     &        SCMop%lastencounter(d) /= stepnmbr(SCMop)-1) then
            print*,'SCMoutput ERROR: last encounter with this '//       &
     &           'diagnostic was not last timestep: ',sname,            &
     &           SCMop%lastencounter(d),stepnmbr(SCMop),                &
     &           (SCMop%daycount-1),SCMop%full_daysteps,                &
     &           SCMop%stepcount
         endif

         ! Record the fact that this diagnostic was seen by
         ! this routine on this timestep.
         SCMop%lastencounter(d)=stepnmbr(SCMop)

         ! Calculate how many timesteps we are through the current
         ! dumping period. 1=first time step, dump_step=last timestep
         ! (and so a dump will occur this timestep).
         nThroughPeriod=mod(stepnmbr(SCMop)-1,SCMop%dump_step(d))+1

         ! Decide whether we are at the start of a dumping period, at
         ! the end of of a dumping period, and whether we need to call
         ! add2dump (using nThroughPeriod this is trivial for most
         ! diagnostics, but has added complications in the case of
         ! diagnostics only calculated on radiation timesteps)
         startperiod=.false.
         endperiod=.false.
         call_add2dump=.true.
         if (.not.SCMop%only_radsteps(d)) then
            ! Diagnostic d is a normal diagnostic, valid at every
            ! timestep
            if (nThroughPeriod == 1) startperiod=.true.
            if (nThroughPeriod == SCMop%dump_step(d)) endperiod=.true.

         else
            ! Diagnostic d is based on a variable which only has
            ! valid values on radiation timesteps.

            if (mod(SCMop%stepcount-ntrad1,ntrad) /= 0)                 &
     &           then
               ! This is not a radiation timestep, assume input
               ! array x contains nonsense information - do not
               ! call add2dump.
               call_add2dump=.false.
            else
               ! The criteria for startperiod and endperiod are now
               ! altered slightly, since startperiod must be true if
               ! this is the first radiation time step during this
               ! dumping period, and endperiod must be true if this is
               ! the last radiation time step during this dumping
               ! period.
               if (nThroughPeriod-1 <  ntrad) startperiod=.true.
               if (SCMop%dump_step(d)-nThroughPeriod <  ntrad)          &
     &                                          endperiod=.true.

               ! This is a radiation timestep and so we can call
               ! add2dump, but the no. of timesteps by which to divide
               ! in order to calculate the average (or whatever) is
               ! not dump_step, but the no. of times this part of the
               ! code has been reached during this dumping period,
               ! stored in SCMop%nadd2dump(d) (which for normal
               ! diagnostics is set to dump_step in newdiag)
               if (startperiod) then
                  SCMop%nadd2dump(d)=1
               else
                  SCMop%nadd2dump(d)=SCMop%nadd2dump(d)+1
               endif

            endif            ! (mod(SCMop%stepcount-ntrad1,ntrad) /= 0)
         endif             ! (.not.SCMop%only_radsteps(d))

         if (call_add2dump) then
! DEPENDS ON: add2dump
            call add2dump(x,SCMop%nelements(d),d,SCMop,                 &
     &           startperiod,endperiod)


            Do j=1, SCMop%nelements(d)

              write(sdum0,*) x(j)

              If (index(sdum0, 'NaN') /= 0) Then
                call scm_trap_nan(sname, calling_routine) 
                Exit
              End If

              If (index(sdum0, 'nan') /= 0) Then
                call scm_trap_nan(sname, calling_routine) 
                Exit
              End If

              If (index(sdum0, 'NAN') /= 0) Then
                call scm_trap_nan(sname, calling_routine) 
                Exit
              End If

            End do

         endif

      enddo

 999  continue
      return
      END SUBROUTINE SCMoutput
#endif
