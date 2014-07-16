#if defined(FLUXPROC) || defined(FLXPLPR)
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
! contains routines: interp_time
!
! Purpose: Flux processing routine.
!          Sets fieldout = field1 * weight1 + field2 * weight2
!          and changes date in lookup table to that of
!          validity time input
!----------------------------------------------------------------------
      subroutine interp_time(Int_Head, ncols, nrows, rmdi,              &
#include "avaltim.h"
     &         weight1, weight2, Field1, Field2, FieldOut)

      implicit none

! declaration of parameters
#include "clookadd.h"
#include "plookups.h"

! declaration of argument list
      integer Int_Head(Len_IntHd) ! IN/OUT   date is changed
      integer ncols               ! IN  number of columns
      integer nrows               ! IN  number of rows
      real    rmdi                ! IN  real missing data indicator
! validity time to insert in Lookup table: intent IN
#include "cvaltim.h"
      real weight1   ! IN weight to give to 1st climate field
      real weight2   ! IN weight to give to 2nd climate field
      real Field1(ncols,nrows)   ! IN  1st field
      real Field2(ncols,nrows)   ! IN  2nd field
      real FieldOut(ncols,nrows) ! OUT interpolated field

! declaration of local scalars
      integer jrow  ! row number
      integer icol  ! column number
!----------------------------------------------------------------------

! 1. do time interpolation

      do jrow = 1, nrows
        do icol = 1, ncols
          if ( Field1 (icol, jrow)  /=  rmdi .and.                      &
     &         Field2 (icol, jrow)  /=  rmdi       ) then

            FieldOut (icol, jrow) = weight1 * Field1 (icol, jrow)       &
     &                            + weight2 * Field2 (icol, jrow)

          else
            FieldOut (icol, jrow) = rmdi

          end if
        end do    ! icol
      end do      ! jrow

! 2. set validity time in integer lookup table

      Int_Head(LBYR)  = ValidYear
      Int_Head(LBMON) = ValidMonth
      Int_Head(LBDAT) = ValidDay
      Int_Head(LBHR)  = ValidHour
      Int_Head(LBMIN) = ValidMin

      return
      END SUBROUTINE interp_time
!----------------------------------------------------------------------
#endif
