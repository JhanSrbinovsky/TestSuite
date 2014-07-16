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
!                       interpolation to be split between two programs.
!                       M. J. Bell and A. Hines.
! 5.3      14/08/01     Allow atmosphere fields on B or C grids MJB
!
! Author:     L. Gregorious
!----------------------------------------------------------------------
! contains routines: winds
!
! Purpose: To produce a pp field containing:
!            wind stress (x-direction)
!            wind stress (y-direction)
!            wind mixing energy
!          for all the fields required.
!          Addition of rotated grid (S. Spall)
!----------------------------------------------------------------------
      subroutine winds(                                                 &
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
      integer Int_Head_taux(Len_IntHd)  ! integer part of lookup table
      integer Int_Head_tauy(Len_IntHd)  ! integer part of lookup table
      integer Int_Head_WME(Len_IntHd)  ! integer part of lookup table
      real Real_Head_taux(Len_RealHd)   ! real part of lookup table
      real Real_Head_tauy(Len_RealHd)   ! real part of lookup table
      real Real_Head_WME(Len_RealHd)   ! real part of lookup table
! windstresses on true atmosphere grids
      real winds_u(ncols,nrowsu)
      real winds_v(ncols,nrowsv)
! windstresses with components calculated at coincident points
      real windstressu(ncols, nrowsuv)   ! wind stress fields (taux)
      real windstressv(ncols, nrowsuv)   ! wind stress fields (tauy)
      real wind_mixing_energy(ncols,nrowst) ! wind mixing energy
#if !defined(FLXPLPR)
      real wndu_tmp(ncols, nrowsuv)   ! wind str. on reg. lat-long gri
      real wndv_tmp(ncols, nrowsuv)   ! wind str. on reg. lat-long gri
#endif

! declaration of local scalars

      integer ivt           ! loop index over validity times
      integer IVTOffHr      ! offset of validity time from reference
      integer IOutUnit      ! output unit

      integer i       ! loop index for columns
      integer j       ! loop index for rows

      logical ldebug        ! T => output debugging info (set in 0.)
      logical l_leads       ! T => using minleadsfrac
                            ! F => using minicefrac

! declaration of externals
      external read_fields, write_one_field, read_leads_flds

!----------------------------------------------------------------------
! 0. Preliminaries
!----------------------------------------------------------------------
      CSub = 'winds'  ! subroutine name for error messages

      ldebug = l_winds_dbg     ! set by debug input control file

!----------------------------------------------------------------------
! 1. start loop over validity times
!----------------------------------------------------------------------
      do ivt = 1, NoValidTimes

        IVTOffHr = IValidOffHr(ivt)
        IOutUnit = IOutUnitOff(ivt) + UnitWindsOut

!----------------------------------------------------------------------
! 2.1 read in windstresses and convert to B grid values if necessary
!----------------------------------------------------------------------

! DEPENDS ON: read_fields
        call read_fields(StCWindStressU, IVTOffHr,                      &
     &                 ldebug, Int_Head_taux, Real_Head_taux,           &
     &                 ncols, nrowsu,                                   &
     &                 winds_u,                                         &
#include "argppx.h"
     &                 icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2.1 unable to read u windstresses'
          icode = 1006
          go to 9999
        end if

! DEPENDS ON: read_fields
        call read_fields(StCWindStressV, IVTOffHr,                      &
     &                 ldebug, Int_Head_tauy, Real_Head_tauy,           &
     &                 ncols, nrowsv,                                   &
     &                 winds_v,                                         &
#include "argppx.h"
     &                 icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2.1 unable to read v windstresses'
          icode = 1006
          go to 9999
        end if

        if ( l_b_grid ) then
          do j=1, nrowsuv
            do i = 1, ncols
              windstressu(i,j) = winds_u(i,j)
              windstressv(i,j) = winds_v(i,j)
            end do
          end do

        else   ! l_b_grid
          do j=1, nrowsuv
            do i = 1, ncols-1
              windstressu(i,j) = 0.5*(winds_u(i,j)+winds_u(i,j+1))
              windstressv(i,j) = 0.5*(winds_v(i,j)+winds_v(i+1,j))
            end do
            if (LCyclic) then
              windstressu(ncols,j)=windstressu(1,j)
              windstressv(ncols,j)=windstressv(1,j)
            else
              windstressu(ncols,j)=0.0 ! missing data
              windstressv(ncols,j)=0.0 ! missing data
            end if
          end do

! Update Lookup tables to reflect changes to grid
          Real_Head_taux(BZY-Len_IntHd)=                                &
     &          Real_Head_tauy(BZY-Len_IntHd)
          Real_Head_tauy(BZX-Len_IntHd)=                                &
     &          Real_Head_taux(BZX-Len_IntHd)

        end if  ! l_b_grid


#if !defined(FLXPLPR)
!----------------------------------------------------------------------
! 2.2 Rotate wind vectors if rotated grids are used
!----------------------------------------------------------------------

      if (rotg) then
! DEPENDS ON: w_eqtoll
        call w_eqtoll(coef_angle1, coef_angle2, windstressu,            &
     &           windstressv, wndu_tmp, wndv_tmp, ncols*nrowsuv,        &
     &           ncols*nrowsuv)
      else
        do j = 1, nrowsuv
          do i = 1, ncols
            wndu_tmp(i,j)=windstressu(i,j)
            wndv_tmp(i,j)=windstressv(i,j)
          enddo
        enddo
      endif

      if (rotgO) then
! DEPENDS ON: w_lltoeq
        call w_lltoeq(coef_angle3, coef_angle4, wndu_tmp,               &
     &           wndv_tmp, windstressu, windstressv, ncols*nrowsuv,     &
     &           ncols*nrowsuv)
      else
        do j = 1, nrowsuv
          do i = 1, ncols
            windstressu(i,j)=wndu_tmp(i,j)
            windstressv(i,j)=wndv_tmp(i,j)
          enddo
        enddo
      endif
#endif

!----------------------------------------------------------------------
! 2.3 write out U + V component of windstress
!----------------------------------------------------------------------

! DEPENDS ON: write_one_field
        call write_one_field (                                          &
#include "afields.h"
     &       OutStCTAUX, FFTAUX, PPTAUX, IVTOffHr,                      &
     &          IUGrid, nrowsuv,                                        &
     &          Int_Head_taux, Real_Head_taux, IOutUnit, ldebug,        &
     &       windstressu, icode)
        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2.2 unable to write U component of windstresses'
          icode = 1103
          go to 9999
        end if
! DEPENDS ON: write_one_field
        call write_one_field (                                          &
#include "afields.h"
     &       OutStCTAUY, FFTAUY, PPTAUY, IVTOffHr,                      &
     &          IUGrid, nrowsuv,                                        &
     &          Int_Head_tauy, Real_Head_tauy, IOutUnit, ldebug,        &
     &       windstressv, icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2.2 unable to write V component of windstresses'
          icode = 1104
          go to 9999
        end if

!----------------------------------------------------------------------
! 3. Read in wind mixing energy
!----------------------------------------------------------------------
        l_leads = .true.                ! set to use minleadsfrac
! DEPENDS ON: read_leads_flds
        call read_leads_flds(StCWindMixEng, StCAICE,                    &
     &                    IVTOffHr, ldebug,                             &
     &                    l_leads,Int_Head_WME,                         &
     &                    Real_Head_WME, ncols, nrowst,                 &
     &                    wind_mixing_energy,                           &
#include "argppx.h"
     &                    icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3. unable to read Wind Mixing Energy'
          icode = 1007
          go to 9999
        end if

!----------------------------------------------------------------------
! 3.1 Write out wind mixing energy
!----------------------------------------------------------------------
! DEPENDS ON: write_one_field
        call write_one_field (                                          &
#include "afields.h"
     &       OutStCWME, FFWME, PPWME, IVTOffHr,                         &
     &          ITGrid, nrowst,                                         &
     &          Int_Head_WME, Real_Head_WME, IOutUnit, ldebug,          &
     &       wind_mixing_energy, icode)
        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3.1. unable to write wind mixing energy'
          icode = 1105
          go to 9999
        end if

!----------------------------------------------------------------------
! Last. end loop over validity times
!----------------------------------------------------------------------
        enddo    !  ivt

9999  continue
      return
      END SUBROUTINE winds
!----------------------------------------------------------------------
#endif
