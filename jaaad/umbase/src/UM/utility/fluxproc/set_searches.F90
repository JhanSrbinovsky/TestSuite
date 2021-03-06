#if defined(FLUXPROC) || defined(FLXPLIN)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
! History:
! version  date         change
! 4.5      03/09/98     New code
!
! Author:     M. J. Bell
!----------------------------------------------------------------------
! contains routines: set_searches
!
! Purpose: Flux processing routine.
!          To determine the number of times spiral_s should be called
!          for a particular grid
!----------------------------------------------------------------------
      subroutine set_searches ( ncols, nrows, LCyclic, sea_land,        &
     &           lsm, n_pts_unres, index_unres, max_no_searches,        &
     &           n_calls_spiral, n_pts_spiral, icode)

      implicit none

! declaration of argument list

! intent IN
      integer ncols    ! IN # of columns in target field
      integer nrows    ! IN # of rows in target field
      logical LCyclic  ! IN T => target grid is cyclic
      integer sea_land ! IN =0 for sea field  =1/-1 for land field
      integer lsm (ncols * nrows) ! IN land sea mask of target grid
      integer n_pts_unres  ! IN # of unresolved points on target grid
      integer index_unres(n_pts_unres) ! IN indices of unresolved points
      integer max_no_searches ! IN maximum number of searches allowed

! intent OUT or IN/OUT
      integer n_calls_spiral  ! OUT number of searches required
      integer n_pts_spiral(max_no_searches) ! OUT nsearch to use
                                            !     on each search
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected


! declaration of parameters
#include "c_mdi.h"

! declaration of globals used
#include "cunitnos.h"
#include "cmess.h"

! declaration of local arrays
       real field ( ncols * nrows )    ! local field (on target grid)
       integer index_work ( n_pts_unres ) ! working set of
                        ! unresolved pts (updated by call to spiral_s)
       integer n_pts_unresolved (0:max_no_searches)  !  number of pts
                        ! unresolved after each call

! declaration of local scalars
       integer i, iloop   ! loop indices
       integer n_pts_work ! working copy of # of unresolved points
                          ! (updated by each call to spiral_s)
       integer nsearch    ! number of points in search "radius"

       external spiral_s
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'set_searches'  ! subroutine name for error messages

! 1. Set up a data field with rmdi values at unresolved points
!    and zero values elsewhere.

      do i = 1, ncols * nrows
        field ( i ) = 0.0
      end do

      do i = 1, n_pts_unres
        field ( index_unres ( i ) ) = rmdi
      end do

! 2. Make working copies of index_unres and n_pts_unres
      n_pts_work = n_pts_unres
      n_pts_unresolved(0) = n_pts_unres

      do i = 1, n_pts_unres
        index_work ( i )  = index_unres ( i )
      end do

! 3. Start loop which increments the number of points in search radius;
!    double search radius at the start of each loop

      do iloop = 1, max_no_searches
        if (iloop  ==  1) then
          nsearch = 1
        else
          nsearch = 2 * nsearch
        end if

! 4. call spiral_s

! DEPENDS ON: spiral_s
        call spiral_s(lsm,index_work,n_pts_work,                        &
     &       nrows,ncols,field,nsearch,sea_land,LCyclic)

! save value of nsearch and number of unresolved points
        n_pts_spiral (iloop) = nsearch
        n_pts_unresolved(iloop) = n_pts_work

! 5. if # of unresolved points is zero, exit the routine
        if ( n_pts_work  ==  0) then
           n_calls_spiral = iloop
           write(UnStd,*)CStd,CSub, ' step 5.: ', n_calls_spiral,       &
     &  ' calls to spiral_s resolve all sea points on this grid. ',     &
     &  ' Points per call are: ',                                       &
     &  (n_pts_unresolved(i), i=0,iloop-1)
           go to 9999
        end if

! 6. end loop incrementing search radius and return with an error

      end do !  iloop

      icode = 33
      write(UnErr,*)CErr,CSub, ' step 6.: ', n_pts_work,                &
     &  ' unresolved points still exist after max_no_searches. ',       &
     &  ' Points per call are: ',                                       &
     &  (n_pts_unresolved(i), i=0,max_no_searches)

9999  continue
      return
      END SUBROUTINE set_searches
!----------------------------------------------------------------------
#endif
