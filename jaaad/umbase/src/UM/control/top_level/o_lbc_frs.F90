#if defined(OCEAN)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Set open boundaries using flow relaxation scheme
!
! Subroutine Interface
      subroutine o_lbc_frs(ncols,nrows,n_cols_bdy,n_rows_bdy,           &
     &     levn, d1_field, f_type,                                      &
     &     levels_dataset,rimwidth,rimweights,                          &
     &     l_obdy_north,l_obdy_east,l_obdy_south,l_obdy_west,           &
     &     bdy_n, bdy_e, bdy_s, bdy_w,                                  &
     &     icode, cmessage)

      implicit none
!---------------------------------------------------------------------
! Program: o_lbc_frs
!
! Purpose: This routine creates open north/east/south/west
!          boundaries using the Flow Relaxation Scheme (FRS).
!          This is the scheme which is presently used in the
!          atmosphere LAM's and is based on the subroutine
!          mergefld called from boundval.
!
! Author: Mike Bell
!
! History
! Model   Date      Modification history from model version 4.5
! version
! 4.5   07/07/98    New subroutine/deck.  M.Bell,S.Ineson
! 5.2   31/07/00    Remove l_apply_lsm and land_val from argument
!                   list, adjust definitions of fld_t etc and adjust
!                   checks for data at land points. M J Bell
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!   6.0   01/09/03   Optimisations for SX6 - moving writes
!                    outside loops to allow vectorisation. R. Hill

!
!--------------------------------------------------------------------
! Subroutine Arguments:

      integer ncols           ! IN number of columns (on tracer grid)
      integer nrows           ! IN number of rows
      integer n_cols_bdy      ! IN number of cols in boundary data
      integer n_rows_bdy      ! IN number of rows in boundary data
      integer levn            ! IN present level number

      real d1_field(ncols*nrows) ! IN/OUT section of the D1 array to
                                 ! be updated

      integer f_type          ! IN: field type indicator
      real levels_dataset(ncols,nrows) ! IN levels dataset for this
                                       ! field type

      integer rimwidth               ! IN width of rim to update(max 4)
      real rimweights(rimwidth)      ! IN rim weights for each row

      logical l_obdy_north           ! IN T=> update north lateral bdy
      logical l_obdy_east            ! IN T=> update east lateral bdy
      logical l_obdy_south           ! IN T=> update south lateral bdy
      logical l_obdy_west            ! IN T=> update west lateral bdy

      real bdy_n(n_cols_bdy*rimwidth) ! IN north part of boundary data
      real bdy_e(n_rows_bdy*rimwidth) ! IN east part of boundary data
      real bdy_s(n_cols_bdy*rimwidth) ! IN south part of boundary data
      real bdy_w(n_rows_bdy*rimwidth) ! IN west part of boundary data


      integer icode             ! OUT Return code
      character*(80) cmessage   ! OUT Error message
!L--------------------------------------------------------------------

!L Global variables
#include "parvars.h"
#if !defined(MPP)
      integer Offx, Offy
      logical atright, atleft, atbase, attop
#endif
#include "c_mdi.h"
#include "coctol.h"

! Local parameters
      INTEGER                                                           &
     &   fld_t                                                          &
                         ! indicates a tracer
     &,  fld_u                                                          &
                         ! indicates a velocity
     &,  fld_sf          ! indicates a streamfn
      PARAMETER (fld_t=1, fld_u=2, fld_sf=3)

#include "decomptp.h"
#include "decompdb.h"
! Local scalars

      integer i             ! loop over columns
      integer j             ! loop over rows
      integer irim          ! pointer to position in bdy data arrays
      integer ifld          ! pointer to position in field to update
      real rwt              ! rimweight for present boundary pt

      integer row_start  ! start row (not counting offsets) for E & W
      integer row_end    ! end row (not counting offsets) for E & W bdys
      integer irim_start ! initial value of irim for E and W bdys -
                         ! depends on value of l_obdy_south.

!----------------------------------------------------------------------
!L 0.0 Set grid off-sets and indicators showing if domain lies next to
!L     each model boundary; non-mpp settings for offsets are zero and
!L     for boundary indicators are true. This allows the same code
!L     to be used for MPP and non-MPP cases.
!L

#if !defined(MPP)
      Offx = 0
      Offy = 0
      atright = .true.
      atleft  = .true.
      atbase  = .true.
      attop   = .true.
#endif

!L 0.1 Set row_start and row_end and irim_start

      if (at_extremity(peast).and. l_obdy_east .or.                     &
     &     at_extremity(pwest) .and. l_obdy_west) then

        row_start = 1
        irim_start=1
        if (at_extremity(psouth) .and. l_obdy_south) then
          row_start = rimwidth + 1
          irim_start = rimwidth * rimwidth + 1
        endif

        row_end = n_rows_bdy   !  first guess !
        if (at_extremity(pnorth) .and. l_obdy_north) then
          row_end = n_rows_bdy - rimwidth
        else if (at_extremity(pnorth) .and.                             &
     &                       .not. l_obdy_north .and.                   &
     &            f_type  ==  fld_u ) then
          row_end = n_rows_bdy - 1
        end if

      end if  ! if east or west boundaries to update


!L 1.1 Northern boundary: set all points including those in corners at
!L     east or west boundaries

      if ( l_obdy_north .and. at_extremity(pnorth)) then

       irim = 1

       do j = n_rows_bdy-rimwidth+1,n_rows_bdy
         do i = 1, n_cols_bdy

! note: ncols IS the number of columns of data stored for all cases
           ifld = i + Offx + (j+Offy-1)*ncols

           rwt = rimweights(n_rows_bdy+1-j) ! first guess

! rim weights are reduced in the corners
           if (at_extremity(pwest) .and.                                &
     &                       i  <   n_rows_bdy+1-j) then
             rwt = rimweights(i)
           else if ( at_extremity(pwest) .and.                          &
     &               i  >   n_cols_bdy - n_rows_bdy + j ) then
             rwt = rimweights(n_cols_bdy + 1 - i)
           end if

           if ( levels_dataset(i+Offx,j+Offy)  >=  levn) then
             d1_field(ifld) = bdy_n(irim)*rwt +                         &
     &                        d1_field(ifld)*(1.0-rwt)
             if ( abs(bdy_n(irim)-rmdi)  <   abs(rmdi*TOL_SMALL) )then
               if ( f_type  /=  fld_sf) then
                 icode = 1
               end if
             end if
           end if

           irim=irim+1

         enddo   ! i
        enddo  ! j

        if (icode == 1)  write(6,*) 'o_frs_lbcs: FATAL ERROR north:'
       endif  ! northern boundary

!L 1.2 Southern boundary: as for northern boundary


      if (l_obdy_south .and. at_extremity(psouth) ) then

      irim=1

      do j = 1,rimwidth
         do i = 1, n_cols_bdy

           ifld = i + Offx + (j+Offy-1)*ncols

           rwt = rimweights(j)

           if (at_extremity(pwest) .and. i  <   j) then
             rwt = rimweights(i)
           else if (at_extremity(peast) .and.                           &
     &               i  >   n_cols_bdy + 1 - j ) then
             rwt = rimweights(n_cols_bdy + 1 - i)
           endif

           if ( levels_dataset(i+Offx,j+Offy)  >=  levn) then
             d1_field(ifld) = bdy_s(irim)*rwt +                         &
     &                        d1_field(ifld)*(1.0-rwt)
             if ( abs(bdy_s(irim)-rmdi)  <   abs(rmdi*TOL_SMALL) )then
               if ( f_type  /=  fld_sf) then
                 icode = 1
               end if
             end if
           end if

           irim=irim+1

         enddo  ! i
        enddo  ! j
        if (icode == 1)  write(6,*) 'o_frs_lbcs: FATAL ERROR south:'
       endif  ! southern boundary

!L 1.3 Eastern boundary: does not calculate values for corner points
!L     if they have already been set.

           if (at_extremity(peast) .and. l_obdy_east) then

        irim = irim_start

        do j = row_start,row_end

          do i = n_cols_bdy - rimwidth + 1, n_cols_bdy

            ifld = i + Offx + (j+Offy-1)*ncols
            rwt = rimweights(n_cols_bdy+1-i)

           if ( levels_dataset(i+Offx,j+Offy)  >=  levn) then
             d1_field(ifld) = bdy_e(irim)*rwt +                         &
     &                        d1_field(ifld)*(1.0-rwt)
             if ( abs(bdy_e(irim)-rmdi)  <   abs(rmdi*TOL_SMALL) )then
               if ( f_type  /=  fld_sf) then
                 icode = 1
               end if
             end if
            end if

            irim=irim+1

          enddo
        enddo

        if (icode == 1)  write(6,*) 'o_frs_lbcs: FATAL ERROR east:'
       endif

! 1.4 Western boundary: same as for eastern boundary

       if (at_extremity(pwest) .and.l_obdy_west) then

        irim =  irim_start

        do j = row_start,row_end

          do i = 1,rimwidth

            ifld = i + Offx + (j+Offy-1)*ncols
            rwt = rimweights(i)

            if ( levels_dataset(i+Offx,j+Offy)  >=  levn) then
             d1_field(ifld) = bdy_w(irim)*rwt +                         &
     &                        d1_field(ifld)*(1.0-rwt)
             if ( abs(bdy_w(irim)-rmdi)  <   abs(rmdi*TOL_SMALL) )then
               if ( f_type  /=  fld_sf) then
                 icode = 1
               end if
             end if
            end if

            irim=irim+1

          enddo   ! i
        enddo   ! j
        if (icode == 1)  write(6,*) 'o_frs_lbcs: FATAL ERROR west:'

       endif   ! western boundary

       return

       END SUBROUTINE o_lbc_frs
#endif
