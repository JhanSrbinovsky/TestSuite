#if defined(FLUXPROC) || defined(FLXPLPR) || defined(FLXPLIN)
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
! 6.1      30/07/04     Correct call to spiral_s to use LCyclicO.
!                       A. Hines.
!
! Author:     M. J. Bell
!----------------------------------------------------------------------
! contains routines: write_one_field,write_pp
!
! Purpose: Flux processing routine.
!          write_one_field:
!          This routine writes out one pp-field to the required output
!          pp file. The pp-header is atered for the correct stashcode
!          and the field written to the file using the routine write_pp.
!          The routine also sets land points to missing data and checks
!          for consistency in the field dimensions.
!          write_pp: Writes out a pp header then a pp field
!----------------------------------------------------------------------
      subroutine write_one_field (                                      &
#include "afields.h"
     &       StCode, FFCode, PPCode, IVTOffHr,                          &
     &       IGridtype, nrows,                                          &
     &       Int_Head, Real_Head, IOutUnit, ldebug,                     &
     &       field_atm, icode )

      implicit none

! declaration of parameters used in argument list
#include "plookups.h"

! declaration of argument list

! array dimensions, lsms, interpolation coeffs etc. : all intent IN
#include "cfields.h"

! field codes to insert in integer header that is output
      integer StCode   ! IN stash code
      integer FFCode   ! IN Met O 8 field code
      integer PPCode   ! IN PP package code
      integer IVTOffHr ! IN offset of validity time from reference

! other input
      integer IGridtype  ! IN  grid type (0 = tracer, 1 = velocity)
      integer nrows      ! IN  number of rows in input field

      integer Int_Head(Len_IntHd)   ! IN integer part of lookup table
      real Real_Head(Len_RealHd)    ! IN real part of lookup table
      integer IOutUnit   ! IN  output unit
      logical ldebug     ! IN  T => output debugging info

      real field_atm( ncols, nrows ) ! IN  field on NWP grid
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected


! declaration of parameters
#include "c_mdi.h"

! declaration of globals used
#include "cunitnos.h"
#include "cmess.h"
#include "clookups.h"
#include "cvaloff.h"

! declaration of local arrays
#if !defined(FLXPLPR) && !defined(FLXPLIN)
      real field_ocean( ncolsO, nrowstO ) ! used for both t and u cases
      integer index_unres ( ncolsO * nrowstO ) ! indices to unresolved
                                               ! points on ocean grid
#endif
!----------------------------------------------------------------------

! declaration of local scalars

      integer ipts          ! loop index over unresolved points
      integer isearch       ! loop index over calls to spiral_s
      integer nsearch       ! # of pts in search "radius"
      integer n_pts_unres   ! local counter of # of unresolved points
      integer ncolsOut      ! # of columns in output field
      integer nrowsOut      ! # of rows in output field

#if defined(FLXPLPR) || defined(FLXPLIN)
      external write_pp,amend_lookup_flux
#else
      external lsm_set,h_int_lsm,spiral_s,write_pp
#endif
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'write_one_field'  ! subroutine name for error messages

! 0.1 check that nrows and IGridtype are consistent
      if (  IGridtype  ==  0 ) then

        if ( nrows  /=  nrowst ) then
          icode = 44
          write(UnErr,*)CErr,CSub,                                      &
     &       ' 0.1.1 nrows and IGridtype inconsistent: ',               &
     &       ' nrows, nrowst, IGridtype =', nrows, nrowst, IGridtype
          go to 9999
        end if

      else if ( IGridtype  ==  1 ) then

        if ( nrows  /=  nrowsuv ) then
          icode = 45
          write(UnErr,*)CErr,CSub,                                      &
     &       ' 0.1.2 nrows and IGridtype inconsistent: ',               &
     &       ' nrows, nrowsuv, IGridtype =', nrows, nrowsuv, IGridtype
          go to 9999
        end if

      else

        icode = 46
        write(UnErr,*)CErr,CSub,                                        &
     &       ' 0.1.3 not coded for IGridtype =', IGridtype
          go to 9999

      end if ! IGridtype

#if defined(FLXPLPR) || defined(FLXPLIN)
! 1. Amend grid information in lookup table
! DEPENDS ON: amend_lookup_flux
      call amend_lookup_flux(                                           &
     &                       Int_Head, Real_Head,                       &
     &                       StCode, FFCode, PPCode, IVTOffHr )


! 2. write out pp field on atmosphere grid
! DEPENDS ON: write_pp
      call write_pp(IOutUnit, Int_Head, Real_Head,                      &
     &              ncols, nrows, field_atm, icode)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 2. error writing out a pp header and field  '
        go to 9999
      end if
#else
! 1. Set land points to missing data (use atmosphere grids)

! for tracer grid
      if ( IGridtype  ==  0 ) then
! DEPENDS ON: lsm_set
        call lsm_set( ncols, nrows, lsmt, ILandPt,                      &
     &       rmdi, ldebug, field_atm )

      else if ( IGridtype  ==  1 ) then
! DEPENDS ON: lsm_set
        call lsm_set( ncols, nrows, lsmu, ILandPt,                      &
     &       rmdi, ldebug, field_atm )

      end if

! 2. Interpolate to ocean grid

      if ( IGridtype  ==  0) then

        ncolsOut = ncolsO
        nrowsOut = nrowstO
! DEPENDS ON: h_int_lsm
        call h_int_lsm(nrowst,ncols,ncolsOut*nrowsOut, rmdi,            &
     &     index_bl_t,index_br_t, field_atm,                            &
     &     weight_bl_t,weight_br_t,weight_tl_t,weight_tr_t,             &
     &     lsmtO,                                                       &
     &     field_ocean)

      else if ( IGridtype  ==  1) then

        ncolsOut = ncolsO
        nrowsOut = nrowsuO
! DEPENDS ON: h_int_lsm
        call h_int_lsm(nrowsu,ncols,ncolsOut*nrowsOut, rmdi,            &
     &     index_bl_u,index_br_u, field_atm,                            &
     &     weight_bl_u,weight_br_u,weight_tl_u,weight_tr_u,             &
     &     lsmuO,                                                       &
     &     field_ocean)

      end if

! 3. fill in coastal values

! 3.1 for a tracer grid
      if ( IGridtype  ==  0) then

! 3.1.1 copy unresolved points into a local array (which is
!       updated by each call to spiral_s)

        n_pts_unres = n_pts_unres_t
        do ipts = 1, n_pts_unres
          index_unres(ipts) = index_unres_t(ipts)
        end do

! 3.1.2 do spiral searches

        do isearch = 1, n_calls_spiral_t

          nsearch = n_pts_spiral_t(isearch)

! DEPENDS ON: spiral_s
          call spiral_s(lsmtO,index_unres,n_pts_unres,                  &
     &      nrowsOut,ncolsOut,field_ocean,nsearch,ISeaPt,LCyclicO)

        end do ! isearch

! 3.2 for a velocity grid
      else if ( IGridtype  ==  1) then

! 3.2.1 copy unresolved points into a local array (which is
!       updated by each call to spiral_s)

        n_pts_unres = n_pts_unres_u
        do ipts = 1, n_pts_unres
          index_unres(ipts) = index_unres_u(ipts)
        end do

! 3.2.2 do spiral searches

        do isearch = 1, n_calls_spiral_u

          nsearch = n_pts_spiral_u(isearch)

! DEPENDS ON: spiral_s
          call spiral_s(lsmuO,index_unres,n_pts_unres,                  &
     &      nrowsOut,ncolsOut,field_ocean,nsearch,ISeaPt,LCyclicO)

        end do ! isearch


      end if  ! IGridtype

! 4. Reset missing data values at land points if user has
!    chosen to do so

      if ( output_land_value   /=  rmdi ) then
        if ( IGridtype  ==  0) then
! DEPENDS ON: lsm_set
          call lsm_set( ncolsOut, nrowsOut, lsmtO, ILandPt,             &
     &                  output_land_value, ldebug, field_ocean )
        else if ( IGridtype  ==  1 ) then
! DEPENDS ON: lsm_set
          call lsm_set( ncolsOut, nrowsOut, lsmuO, ILandPt,             &
     &                  output_land_value, ldebug, field_ocean )
        end if
      end if

! 5. Amend grid information in lookup table
      if ( IGridtype  ==  0) then
! DEPENDS ON: amend_lookup
        call amend_lookup (  LookuplsmtO, Int_Head, Real_Head,          &
     &                       output_land_value,                         &
     &                       StCode, FFCode, PPCode, IVTOffHr )

      else if ( IGridtype  ==  1 ) then
! DEPENDS ON: amend_lookup
        call amend_lookup (  LookuplsmuO, Int_Head, Real_Head,          &
     &                       output_land_value,                         &
     &                       StCode, FFCode, PPCode, IVTOffHr )

      end if

! 6. write out filled pp field on ocean grid
! DEPENDS ON: write_pp
      call write_pp(IOutUnit, Int_Head, Real_Head,                      &
     &              ncolsOut, nrowsOut, field_ocean, icode)
      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 5. error writing out a pp header and field  '
        go to 9999
      end if
#endif

9999  continue
      return
      END SUBROUTINE write_one_field
!----------------------------------------------------------------------
!----------------------------------------------------------------------
#endif
