#if defined(FLXPLIN)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! contains routines: field_interpolate_write
!
! Purpose: Flux processing routine.
!          field_interpolate_write
!          This routine writes out one pp-field to the required output
!          pp file; it interpolates the field to a new grid, fills in
!          missing data points, checks for consistency in the field
!          dimensions and updates the Real and Integer Headers so that
!          they describe the new grid.
!
!    Model            Modification history:
!   version  Date
!    5.3  15/10/01  New deck. A. Hines
!
!    Programming standard :
!
!    Logical components covered :
!
!    System task:
!
!    External documentation:
!----------------------------------------------------------------------
      subroutine field_interpolate_write(                               &
#include "afields.h"
     &       Int_Head, Real_Head, ldebug, IGridtype, nrows,             &
     &       IOutUnit, field_atm, icode )

      implicit none

! declaration of parameters used in argument list
#include "plookups.h"

! declaration of argument list

#include "cfields.h"
      integer Int_Head(Len_IntHd)   ! IN integer part of lookup table
      real Real_Head(Len_RealHd)    ! IN real part of lookup table

      logical ldebug     ! IN  T => output debugging info

      integer IGridtype  ! IN  grid type (0 = tracer, 1 = velocity)
      integer nrows      ! IN  number of rows in input field

      integer IOutUnit   ! IN  output unit
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
      real field_ocean( ncolsO, nrowstO ) ! used for both t and u cases
      integer index_unres ( ncolsO * nrowstO ) ! indices to unresolved
                                               ! points on ocean grid
!----------------------------------------------------------------------

! declaration of local scalars

      integer ipts          ! loop index over unresolved points
      integer isearch       ! loop index over calls to spiral_s
      integer nsearch       ! no of pts in search "radius"
      integer n_pts_unres   ! local counter of no of unresolved points
      integer ncolsOut      ! no of columns in output field
      integer nrowsOut      ! no of rows in output field

      external lsm_set, h_int_lsm, spiral_s, write_pp
!              amend_lookup
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'write_one_field'  ! subroutine name for error messages

      output_land_value = rmdi

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
        call h_int_lsm(nrowsuv,ncols,ncolsOut*nrowsOut, rmdi,           &
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
     &      nrowsOut,ncolsOut,field_ocean,nsearch,ISeaPt,LCyclic)

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
     &      nrowsOut,ncolsOut,field_ocean,nsearch,ISeaPt,LCyclic)

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
! DEPENDS ON: amend_lookup_grid
        call amend_lookup_grid (  LookuplsmtO, Int_Head, Real_Head,     &
     &                       output_land_value)

      else if ( IGridtype  ==  1 ) then
! DEPENDS ON: amend_lookup_grid
        call amend_lookup_grid (  LookuplsmuO, Int_Head, Real_Head,     &
     &                       output_land_value)

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

9999  continue
      return
      END SUBROUTINE field_interpolate_write
!----------------------------------------------------------------------
#endif
