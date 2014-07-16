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
! 5.3      10/10/00     Changes to enable flux selection and grid
!                       interpolation to be split between 2 programs.
!                       M. J. Bell and A. Hines.
! 5.3      14/08/01     Allow atmosphere fields on B or C grids MJB
!
! Author:     S. A. Spall
!----------------------------------------------------------------------
! contains routines: windspeed
!
! Purpose: Flux processing routine.
!          To produce a pp field containing:
!            wind speed (x-direction)
!            wind speed (y-direction)
!          for all the fields required
!----------------------------------------------------------------------
      subroutine windspd(                                               &
#include "afields.h"
#include "argppx.h"
     &                 icode )

      implicit none

! declaration of argument list

! array dimensions, lsms, interpolation coeffs etc. : all intent IN
#include "cfields.h"

      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! declaration of parameters
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "cfdcodes.h"
#include "plookups.h"

! declaration of globals used
#include "cunitnos.h"
#include "cmess.h"
#include "cvaloff.h"
#include "cdebug.h"
#include "cselct.h"
#include "clookadd.h"
#include "clookups.h"


! declaration of local arrays
      integer Int_Head_wspx(Len_IntHd)  ! integer part of lookup table
      integer Int_Head_wspy(Len_IntHd)  ! integer part of lookup table
      real Real_Head_wspx(Len_RealHd)   ! real part of lookup table
      real Real_Head_wspy(Len_RealHd)   ! real part of lookup table
! windspeeds on true atmosphere grids
      real wndsp_u(ncols,nrowsu)
      real wndsp_v(ncols,nrowsv)
! windspeeds with components calculated at coincident points
      real windspeedu(ncols, nrowsuv)   ! wind speed fields
      real windspeedv(ncols, nrowsuv)   ! wind speed fields
#if !defined(FLXPLPR)
      real wndu_tmp(ncols, nrowsuv)   ! wind spd on reg. lat-long grid
      real wndv_tmp(ncols, nrowsuv)   ! wind spd on reg. lat-long grid
#endif

! declaration of local scalars

      integer ivt           ! loop index over validity times
      integer IVTOffHr      ! offset of validity time from reference
      integer IOutUnit      ! output unit

      integer i       ! loop index for columns
      integer j       ! loop index for rows

      logical ldebug        ! T => output debugging info (set in 0.)

! declaration of externals
      external read_fields, write_one_field

!----------------------------------------------------------------------
! 0. Preliminaries
!----------------------------------------------------------------------
      CSub = 'windspeed'  ! subroutine name for error messages

      ldebug = l_windspd_dbg     ! set by debug input control file

!----------------------------------------------------------------------
! 1. start loop over validity times
!----------------------------------------------------------------------
      do ivt = 1, NoValidTimes

        IVTOffHr = IValidOffHr(ivt)
        IOutUnit = IOutUnitOff(ivt) + UnitWindspdOut

!----------------------------------------------------------------------
! 2.1 read in windspeeds and convert to B grid values if necessary
!----------------------------------------------------------------------

! DEPENDS ON: read_fields
        call read_fields(StCWindSpeedU, IVTOffHr,                       &
     &                 ldebug, Int_Head_wspx, Real_Head_wspx,           &
     &                 ncols, nrowsu,                                   &
     &                 wndsp_u,                                         &
#include "argppx.h"
     &                 icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2.1 unable to read u wind speed'
          icode = 1006
          go to 9999
        end if

! DEPENDS ON: read_fields
        call read_fields(StCWindSpeedV, IVTOffHr,                       &
     &                 ldebug, Int_Head_wspy, Real_Head_wspy,           &
     &                 ncols, nrowsv,                                   &
     &                 wndsp_v,                                         &
#include "argppx.h"
     &                 icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2.1 unable to read v wind speed'
          icode = 1006
          go to 9999
        end if

        if ( l_b_grid ) then
          do j=1, nrowsuv
            do i = 1, ncols
              windspeedu(i,j) = wndsp_u(i,j)
              windspeedv(i,j) = wndsp_v(i,j)
            end do
          end do

        else   ! l_b_grid
          do j=1, nrowsuv
            do i = 1, ncols-1
              windspeedu(i,j) = 0.5*(wndsp_u(i,j)+wndsp_u(i,j+1))
              windspeedv(i,j) = 0.5*(wndsp_v(i,j)+wndsp_v(i+1,j))
            end do
            if (LCyclic) then
              windspeedu(ncols,j)=windspeedu(1,j)
              windspeedv(ncols,j)=windspeedv(1,j)
            else
              windspeedu(ncols,j)=0.0 ! missing data
              windspeedv(ncols,j)=0.0 ! missing data
            end if
          end do

! Update Lookup tables to reflect changes to grid
          Real_Head_wspx(BZY-Len_IntHd)=                                &
     &          Real_Head_wspy(BZY-Len_IntHd)
          Real_Head_wspy(BZX-Len_IntHd)=                                &
     &          Real_Head_wspx(BZX-Len_IntHd)

        end if  ! l_b_grid


#if !defined(FLXPLPR)
!----------------------------------------------------------------------
! 2.2 Rotate wind vectors if rotated grids are used
!----------------------------------------------------------------------

      if (rotg) then
! DEPENDS ON: w_eqtoll
        call w_eqtoll(coef_angle1, coef_angle2, windspeedu,             &
     &           windspeedv, wndu_tmp, wndv_tmp, ncols*nrowsuv,         &
     &           ncols*nrowsuv)
      else
        do j = 1, nrowsuv
          do i = 1, ncols
            wndu_tmp(i,j)=windspeedu(i,j)
            wndv_tmp(i,j)=windspeedv(i,j)
          enddo
        enddo
      endif

      if (rotgO) then
! DEPENDS ON: w_lltoeq
        call w_lltoeq(coef_angle3, coef_angle4, wndu_tmp,               &
     &           wndv_tmp, windspeedu, windspeedv, ncols*nrowsuv,       &
     &           ncols*nrowsuv)
      else
        do j = 1, nrowsuv
          do i = 1, ncols
            windspeedu(i,j)=wndu_tmp(i,j)
            windspeedv(i,j)=wndv_tmp(i,j)
          enddo
        enddo
      endif
#endif

!----------------------------------------------------------------------
! 2.3 write out U + V component of wind speed
!----------------------------------------------------------------------

! DEPENDS ON: write_one_field
        call write_one_field (                                          &
#include "afields.h"
     &       OutStCWSPX, FFWSPX, PPWSPX, IVTOffHr,                      &
     &          IUGrid, nrowsuv,                                        &
     &          Int_Head_wspx, Real_Head_wspx, IOutUnit, ldebug,        &
     &       windspeedu, icode)
        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2.2 unable to write U component of wind speed'
          icode = 1103
          go to 9999
        end if
! DEPENDS ON: write_one_field
        call write_one_field (                                          &
#include "afields.h"
     &       OutStCWSPY, FFWSPY, PPWSPY, IVTOffHr,                      &
     &          IUGrid, nrowsuv,                                        &
     &          Int_Head_wspy, Real_Head_wspy, IOutUnit, ldebug,        &
     &       windspeedv, icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2.2 unable to write V component of windstresses'
          icode = 1104
          go to 9999
        end if

!----------------------------------------------------------------------
! Last. end loop over validity times
!----------------------------------------------------------------------
        enddo    !  ivt

9999  continue
      return
      END SUBROUTINE windspd
!----------------------------------------------------------------------
#endif
