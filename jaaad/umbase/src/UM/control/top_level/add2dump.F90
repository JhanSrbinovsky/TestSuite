#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Update the dump array of a SCM diagnostic entry in SCMop

      subroutine add2dump(x,nelements,d,SCMop,startperiod,endperiod)
      implicit none

! Description:
      ! Take input array x and add to the dump array in the manner
      ! specified by the time profile. This performs this timestep's
      ! role in constructing the diagnostic.
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
! Code Description:
      ! Language: Fortran77 with some bits of Fortran90

! SCMop_type is defined in here...
#include "s_scmoptype_defn.h"

      type(SCMop_type) :: SCMop ! INOUT The derived-type structure
                                ! containing all the diagnostic
                                ! information
      integer nelements         ! The no. of elements in array x
      real x(nelements)         ! IN The array from which the diagnostic
                                ! is constructed
      integer d                 ! IN The diagnostic index in SCMop
      logical startperiod,endperiod ! IN These flag the start and end
                                ! of the dumping period, when special
                                ! actions may need to be taken

! Parameters and stuff used around and about the internal
! workings of the diagnostic system
#include "s_scmop_internal.h"

      real, pointer :: dump(:)  ! The dump array in which the
                                ! diagnostic is constructed
      integer timprof           ! The time profile of the diagnostic
      integer nadd2dump         ! The number by which to divide
                                ! accumulated diag's in order to get
                                ! the mean
      integer wnelements        ! The no. of elements in array wdump
      real, pointer :: wdump(:) ! A dump array of another
                                ! diagnostic which may be used in the
                                ! construction of this diagnostic
      integer wtimprof          ! The timeprofile of the other diag.
      integer wnadd2dump        ! As nadd2dump, but for the other
                                ! diagnostic
      integer i                 ! Counter

      ! We don't want lots of '%'s around, make some abbreviations...
      dump=>SCMop%diag(d)%dump
      nelements=SCMop%nelements(d)
      timprof=SCMop%timprof(d)
      nadd2dump=SCMop%nadd2dump(d)
      if (SCMop%wd(d) >  0) then
         wdump=>SCMop%diag(SCMop%wd(d))%dump
         wnelements=SCMop%nelements(SCMop%wd(d))
         wtimprof=SCMop%timprof(SCMop%wd(d))
         wnadd2dump=SCMop%nadd2dump(SCMop%wd(d))
      else
         nullify(wdump)
         wnelements=0
         wtimprof=0
         wnadd2dump=0
      endif

!-T_CONST--A constant value---------------------------------------------
      if (timprof == t_const) then
         do i=1,nelements
            dump(i)=x(i)
         enddo

!-T_INST--Instaneous value----------------------------------------------
      elseif (timprof == t_inst) then
         if (startperiod) then
            ! Set the dump to a default value. Shouldn't really be
            ! necessary, but handy for checking everything's working OK.
            do i=1,nelements
               dump(i)=-999.
            enddo
         endif
         if (endperiod) then
            ! The dump will ocurr at this timestep, so we want the
            ! value now.
            do i=1,nelements
               dump(i)=x(i)
            enddo
         endif

!-T_ACC--Accumulated value----------------------------------------------
      elseif (timprof == t_acc) then
         if (startperiod) then
            ! Take this timestep's values
            do i=1,nelements
               dump(i)=x(i)
            enddo
         else
            ! Add to the current values
            do i=1,nelements
               dump(i)=dump(i)+x(i)
            enddo
         endif

!-T_AVG--Value averaged over a period----------------------------------
      elseif (timprof == t_avg) then
         if (startperiod) then
            ! Take this timestep's values
            do i=1,nelements
               dump(i)=x(i)
            enddo
         else
            ! Add to the current values
            do i=1,nelements
               dump(i)=dump(i)+x(i)
            enddo
         endif
         ! And form average if it's the end of the period
         if (endperiod) then
            do i=1,nelements
               dump(i)=dump(i)/nadd2dump
            enddo
         endif

!-T_MIN--The minimum value over a time period---------------------------
      elseif (timprof == t_min) then
         if (startperiod) then
            ! Take this timestep's values
            do i=1,nelements
               dump(i)=x(i)
            enddo
         else
            ! Overwrite the current values if smaller
            do i=1,nelements
               if (x(i) <  dump(i)) then
                  dump(i)=x(i)
               endif
            enddo
         endif

!-T_MAX--The maximum value over a time period---------------------------
      elseif (timprof == t_max) then
         if (startperiod) then
            ! Take this timestep's values
            do i=1,nelements
               dump(i)=x(i)
            enddo
         else
            ! Overwrite the current value if bigger
            do i=1,nelements
               if (x(i) >  dump(i)) then
                  dump(i)=x(i)
               endif
            enddo
         endif

!-T_DIV-T_MULT--The accumulated value of one diagnostic multiplied or---
!---------------divided by the value of another-------------------------
      elseif (timprof == t_div    .or.timprof == t_mult .or.            &
     &        timprof == t_acc_div.or.timprof == t_acc_mult) then
         if (startperiod) then
            do i=1,nelements
               dump(i)=x(i)
            enddo
         else
            ! Add to the current value
            do i=1,nelements
               dump(i)=dump(i)+x(i)
            enddo
         endif
         if (endperiod) then
            if (timprof == t_mult.or.timprof == t_acc_mult) then
               ! Multiply the two diagnostics
               do i=1,nelements
                  dump(i)=dump(i)*wdump(mod(i-1,wnelements)+1)
               enddo
            elseif (timprof == t_div.or.timprof == t_acc_div) then
               ! Divide the two diagnostics
               do i=1,nelements
                  if (wdump(mod(i-1,wnelements)+1) /= 0) then
                     dump(i)=dump(i)/                                   &
     &                    wdump(mod(i-1,wnelements)+1)
                  elseif (dump(i) /= 0) then
                     write(*,'(A,1X,A,1X,I3,1X,A,1X,I6)')               &
     &                    'add2dump WARNING:divide by zero avoided:'
                  endif
               enddo
            endif
            ! Divide by dumping period?
            if (timprof == t_div.or.timprof == t_mult) then
               do i=1,nelements
                  dump(i)=dump(i)/nadd2dump
               enddo
            endif
         endif

!-----------------------------------------------------------------------
      else
         print*,'add2dump ERROR: I do not know about this ',            &
     &        'temporal type:',timprof
      endif
!-----------------------------------------------------------------------

      return
      END SUBROUTINE add2dump
#endif
