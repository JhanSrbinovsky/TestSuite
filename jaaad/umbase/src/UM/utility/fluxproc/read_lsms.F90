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
! Author:     M. J. Bell
! 5.3      14/08/01     Allow atmosphere fields on B or C grids MJ
! 6.1      30/07/04     Correct call to set_searches to use LCyclicO.
!                       A. Hines.
!----------------------------------------------------------------------
! deck: RDLSMS
!
! contains routines: read_lsms
!
! Purpose: reads land / sea masks, calculates coefficients for
!          interpolation between grids and indices for "unresolved"
!          seapoints on ocean grid which are not surrounded by seapoints
!          on the atmosphere grid.
!          Addition to handle rotated grids (S. Spall)
!----------------------------------------------------------------------
      subroutine read_lsms(                                             &
#include "afields.h"
#include "argppx.h"
     &           icode)

      implicit none

! declaration of arguments

#include "cfields.h"

      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! parameters
#include "plookups.h"
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"

! Globals
#include "cunitnos.h"
#include "cmess.h"
#include "clookups.h"
#include "clookadd.h"

! declaration of local arrays

! arrays required as input by h_int_co
       real phi_full(ncolsO*nrowstO)      ! { allocated with largest
       real lambda_full(ncolsO*nrowstO)   ! { row dimension

! temporary arrays of positions on atmosphere grid required for w_coeff
       real lambda_eqA(ncols*nrowsuv) ! Longitude on atmos  ==  grid
       real phi_eqA(ncols*nrowsuv)    ! Latitude on atmos  ==  grid
       real lambda_tmp1(ncols*nrowsuv) ! Long. on reg. lat-long grid
       real phi_tmp1(ncols*nrowsuv)    ! Latitude on reg. lat-long grid
       real lambda_tmp2(ncols*nrowsuv) ! Longitude on ocean eq. grid
       real phi_tmp2(ncols*nrowsuv)    ! Latitude on ocean eq. grid

! arrays required in the conversion of lat-long for ocean points
       real phi_eq(ncolsO*nrowstO)      ! {
       real lambda_eq(ncolsO*nrowstO)   ! { allocated with largest
       real phi_ll(ncolsO*nrowstO)      ! { row dimension
       real lambda_ll(ncolsO*nrowstO)   ! {

! arrays output by coast_aj which are not subsequently used
! they are allocated with largest ocean grid dimensions (i.e. tracer)
       integer index_targ(ncolsO*nrowstO)
       integer index_srce(ncols*nrowst)
!       integer coastal_points(ncolsO*nrowstO)
       integer coastal_points
       integer index_land_unres(ncolsO*nrowstO)

! declaration of local scalars
       logical mask    ! T => land sea mask is provided
       integer i       ! loop index for columns
       integer j       ! loop index for rows
       integer ij      ! loop index for points in 2D field

! scalar output by coast_aj which is not subsequently used
       integer n_pts_unres_land   ! number of unresolved land points

       external read_lsm_anc, set_lsmu, h_int_co,                       &
     &         coast_aj, set_searches, copy_to_real, eqtoll, lltoeq

!----------------------------------------------------------------------

! 0. Preliminaries
      CSub = 'read_lsms'  ! subroutine name for error messages

!----------------------------------------------------------------------
! 1. Read land / sea masks
!----------------------------------------------------------------------

! 1.1 read atmosphere tracer land / sea mask and calculate
!    grid coordinates
! DEPENDS ON: read_lsm_anc
      call read_lsm_anc(UnitNWPlsmt, Len_FixHd, Len1_Lookup, FixHdlsmt, &
     &       Lookuplsmt, ncols, nrowst, lsmt, lambda_t, phi_t,          &
#include "argppx.h"
     &       icode)

! check icode
      if (icode  >   0)then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 1.1 failed reading NWP tracer land / sea mask'
        go to 9999
      end if

! 1.2 set atmosphere velocity land / sea mask from tracer land / sea
!     mask and calculate grid coordinates

! DEPENDS ON: set_lsmu
      call set_lsmu (  ncols, nrowst, nrowsuv, LCyclic,                 &
     &                 lambda_t, phi_t, lsmt, ILandPt,                  &
     &                 lsmu, lambda_u, phi_u )



! 1.3 read ocean tracer land / sea mask and calculate
!    grid coordinates
! DEPENDS ON: read_lsm_anc
      call read_lsm_anc(UnitFOAMlsmt, Len_FixHd, Len1_Lookup,           &
     &       FixHdlsmtO, LookuplsmtO, ncolsO, nrowstO, lsmtO,           &
     &       lambda_tO, phi_tO,                                         &
#include "argppx.h"
     &       icode)

! check icode
      if (icode  >   0)then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 1.3 failed reading ocean tracer land / sea mask'
        icode = icode + 2000
        go to 9999
      end if

! 1.4 read ocean velocity land / sea mask and calculate
!    grid coordinates
! DEPENDS ON: read_lsm_anc
      call read_lsm_anc(UnitFOAMlsmu, Len_FixHd, Len1_Lookup,           &
     &       FixHdlsmuO, LookuplsmuO, ncolsO, nrowsuO, lsmuO,           &
     &       lambda_uO, phi_uO,                                         &
#include "argppx.h"
     &       icode)

! check icode
      if (icode  >   0)then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 1.4 failed reading ocean velocity land / sea mask'
        icode = icode + 2000
        go to 9999
      end if

!----------------------------------------------------------------------
! 2. Find if the atmosphere grid is rotated
!---------------------------------------------------------------------

! 2.1 Get the position of the poles from the
!     atmosphere lsm header

! DEPENDS ON: copy_to_real
      call copy_to_real ( Lookuplsmt(BPLAT), pole_lat )
! DEPENDS ON: copy_to_real
      call copy_to_real ( Lookuplsmt(BPLON), pole_lon )

! 2.2 Find if a rotated grid is being used

      rotg=.true.
      if ( pole_lat  >   89.99 .and. pole_lat  <   90.01) then
        rotg=.false.
      end if

      if (pole_lat  <   -1.0e5) then
        rotg=.false.
      end if

! 2.3 Do error checking on the positions of the poles

      if (rotg) then

        if ( pole_lat  >   90.0 .or. pole_lat  <   -90.0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2.3 incorrect latitude of pole in atmos lsm header'
          icode = icode + 2000
          go to 9999
        end if

        if ( pole_lon  >   360.0 .or. pole_lon  <   -360.0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2.3 incorrect longitude of pole in atmos lsm header'
          icode = icode + 2000
          go to 9999
        end if

      end if

! 2.4 Write out details of the type of grid

      if (rotg) then
        write(UnStd,*)CStd//CSub//'Atmosphere on a rotated grid;',      &
     &    ' BPLAT = ', pole_lat, '; BPLON = ', pole_lon
      else
        write(UnStd,*)CStd//CSub//'Atmosphere on a non-rotated grid'
      end if

!----------------------------------------------------------------------
! 3. Find if the ocean grid is rotated
!---------------------------------------------------------------------

! 3.1 Get the position of the poles from the
!     ocean tracer grid lsm header

! DEPENDS ON: copy_to_real
      call copy_to_real ( LookuplsmtO(BPLAT), poleO_lat )
! DEPENDS ON: copy_to_real
      call copy_to_real ( LookuplsmtO(BPLON), poleO_lon )

! 3.2 Find if a rotated grid is being used

      rotgO=.true.
      if ( poleO_lat  >   89.99 .and. poleO_lat  <   90.01) then
        rotgO=.false.
      end if

      if (poleO_lat  <   -1.0e5) then
        rotgO=.false.
      end if

! 3.3 Do error checking on the positions of the poles

      if (rotgO) then

        if ( poleO_lat  >   90.0 .or. poleO_lat  <   -90.0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3.3 incorrect latitude of pole in ocean lsm header'
          icode = icode + 2000
          go to 9999
        end if

        if ( poleO_lon  >   360.0 .or. poleO_lon  <   -360.0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3.3 incorrect longitude of pole in ocean lsm header'
          icode = icode + 2000
          go to 9999
        end if

      end if

! 3.4 Write out details of the type of grid

      if (rotgO) then
        write(UnStd,*)CStd//CSub//'Ocean on a rotated grid;',           &
     &    ' BPLAT = ', poleO_lat, '; BPLON = ', poleO_lon
      else
        write(UnStd,*)CStd//CSub//'Ocean on a non-rotated grid'
      end if

!----------------------------------------------------------------------
! 4. Calculate interpolation coefficients for interpolation
!    from atmosphere to ocean tracer grids
!----------------------------------------------------------------------

! 4.1 prepare target grid coordinates required
!     by h_int_co for tracer grid

      do j = 1, nrowstO
        do i = 1, ncolsO
          ij = i + (j-1) * ncolsO
          phi_eq( ij ) = phi_tO ( j )
          lambda_eq( ij ) =  lambda_tO ( i )
       end do
      end do

! 4.2 If the ocean uses a rotated grid, convert the ocean lat-long
!     vector to a standard lat-long grid

      if (rotgO) then
! DEPENDS ON: eqtoll
        call eqtoll(phi_eq, lambda_eq, phi_ll, lambda_ll,               &
     &                      poleO_lat, poleO_lon, ncolsO*nrowstO)
      else
        do i = 1, ncolsO*nrowstO
          phi_ll( i ) = phi_eq( i )
          lambda_ll( i ) =  lambda_eq ( i )
        end do
      end if

! 4.3 If the atmosphere uses a rotated grid, convert the ocean standard
!     lat-long to the atmosphere rotated grid

      if (rotg) then
! DEPENDS ON: lltoeq
        call lltoeq(phi_ll, lambda_ll, phi_full, lambda_full,           &
     &                      pole_lat, pole_lon, ncolsO*nrowstO)
      else
        do i = 1, ncolsO*nrowstO
          phi_full( i ) = phi_ll( i )
          lambda_full( i ) =  lambda_ll ( i )
        end do
      end if

! 4.4 Convert target longitude to correct range

      do i = 1, ncolsO*nrowstO
        lambda_full(i)=mod(lambda_full(i)-lambda_t(1)+720.,360.)        &
     &                         +lambda_t(1)
      end do

! 4.5 Calculate interpolation coefficients for tracer grids

! DEPENDS ON: h_int_co
      call h_int_co(index_bl_t,index_br_t,                              &
     & weight_tr_t,weight_br_t,weight_tl_t,weight_bl_t,                 &
     & lambda_t, phi_t, lambda_full, phi_full,                          &
     & ncols,nrowst,ncolsO*nrowstO,LCyclic)

!----------------------------------------------------------------------
! 5. Calculate interpolation coefficients for interpolation
!    from atmosphere to ocean velocity grids. Atmosphere velocity
!    components must be defined at coincident points.
!----------------------------------------------------------------------

! 5.1 prepare target grid coordinates required
!     by h_int_co for velocity grid

      do j = 1, nrowsuO
        do i = 1, ncolsO
          ij = i + (j-1) * ncolsO
          phi_eq( ij ) = phi_uO ( j )
          lambda_eq (ij ) =  lambda_uO ( i )
       end do
      end do

! 5.2 If the ocean uses a rotated grid, convert the ocean lat-long
!     vector to a standard lat-long grid

      if (rotgO) then
! DEPENDS ON: eqtoll
        call eqtoll(phi_eq, lambda_eq, phi_ll, lambda_ll,               &
     &                      poleO_lat, poleO_lon, ncolsO*nrowsuO)
      else
        do i = 1, ncolsO*nrowsuO
          phi_ll( i ) = phi_eq( i )
          lambda_ll( i ) =  lambda_eq ( i )
        end do
      end if

! 5.3 If the atmosphere uses a rotated grid, convert the ocean standard
!     lat-long to the atmosphere rotated grid

      if (rotg) then
! DEPENDS ON: lltoeq
        call lltoeq(phi_ll, lambda_ll, phi_full, lambda_full,           &
     &                      pole_lat, pole_lon, ncolsO*nrowsuO)
      else
        do i = 1, ncolsO*nrowsuO
          phi_full( i ) = phi_ll( i )
          lambda_full( i ) =  lambda_ll ( i )
        end do
      end if

! 5.4 Convert target longitude to correct range

      do i = 1, ncolsO*nrowsuO
        lambda_full(i)=mod(lambda_full(i)-lambda_u(1)+720.,360.)        &
     &                         +lambda_u(1)
      end do

! 5.5 Calculate interpolation coefficients for velocity grids
! DEPENDS ON: h_int_co
      call h_int_co(index_bl_u,index_br_u,                              &
     & weight_tr_u,weight_br_u,weight_tl_u,weight_bl_u,                 &
     & lambda_u, phi_u, lambda_full, phi_full,                          &
     & ncols,nrowsuv,ncolsO*nrowsuO,LCyclic)

!----------------------------------------------------------------------
! 6. Calculate the coefficients for rotating wind vectors to align
!    with the ocean grid
!----------------------------------------------------------------------

! 6.1 Set up the coefficients for atmosphere to reg. lat-long

      do j = 1, nrowsuv
        do i = 1, ncols
          ij = i + (j-1) * ncols
          phi_eqA( ij ) = phi_u ( j )
          lambda_eqA ( ij ) =  lambda_u ( i )
       end do
      end do

      if (rotg) then
! DEPENDS ON: eqtoll
        call eqtoll(phi_eqA, lambda_eqA, phi_tmp1, lambda_tmp1,         &
     &                      pole_lat, pole_lon, ncols*nrowsuv)
! DEPENDS ON: w_coeff
        call w_coeff(coef_angle1, coef_angle2, lambda_tmp1,             &
     &             lambda_eqA, pole_lat, pole_lon, ncols*nrowsuv)
      else
        do ij = 1, ncols*nrowsuv
          lambda_tmp1(ij)=lambda_eqA(ij)
          phi_tmp1(ij)=phi_eqA(ij)
        enddo
      endif

! 6.2 Set up the coefficients for reg. lat-long to ocean

      if (rotgO) then
! DEPENDS ON: lltoeq
        call lltoeq(phi_tmp1, lambda_tmp1, phi_tmp2, lambda_tmp2,       &
     &                      poleO_lat, poleO_lon, ncols*nrowsuv)
! DEPENDS ON: w_coeff
        call w_coeff(coef_angle3, coef_angle4, lambda_tmp1,             &
     &             lambda_tmp2, poleO_lat, poleO_lon, ncols*nrowsuv)
      endif

!----------------------------------------------------------------------
! 7. Calculate indices for unresolved points i.e. seapoints on ocean
!    grid which are not surrounded by seapoints on the atmosphere grid
!----------------------------------------------------------------------

! 7.1 Calculate indices for unresolved points for tracer grids

      mask = .true.  ! land / sea mask for target grid is to be used

! DEPENDS ON: coast_aj
      call coast_aj (index_bl_t,index_br_t,                             &
     & weight_tr_t,weight_br_t,weight_tl_t,weight_bl_t,                 &
     & ncols,nrowst,ncolsO*nrowstO,                                     &
     & lsmt,lsmtO,                                                      &
     & index_targ,index_srce,coastal_points,mask,                       &
     & index_unres_t,n_pts_unres_t,                                     &
     & index_land_unres,n_pts_unres_land)

! 7.2 Calculate indices for unresolved points for velocity grids

! DEPENDS ON: coast_aj
      call coast_aj (index_bl_u,index_br_u,                             &
     & weight_tr_u,weight_br_u,weight_tl_u,weight_bl_u,                 &
     & ncols,nrowsuv,ncolsO*nrowsuO,                                    &
     & lsmu,lsmuO,                                                      &
     & index_targ,index_srce,coastal_points,mask,                       &
     & index_unres_u,n_pts_unres_u,                                     &
     & index_land_unres,n_pts_unres_land)

!----------------------------------------------------------------------
! 8. Determine number of searchs needed to fill in unresolved points
!----------------------------------------------------------------------

! 8.1 on tracer grid
! DEPENDS ON: set_searches
      call set_searches ( ncolsO, nrowstO, LCyclicO, ISeaPt,            &
     &     lsmtO, n_pts_unres_t, index_unres_t, max_no_searches,        &
     &     n_calls_spiral_t, n_pts_spiral_t, icode)

      if (icode  >   0)then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 4.1 unable to fill in all unresolved points '
        icode = icode + 2000
        go to 9999
      end if


! 8.2 on velocity grid
! DEPENDS ON: set_searches
      call set_searches ( ncolsO, nrowsuO, LCyclicO, ISeaPt,            &
     &     lsmuO, n_pts_unres_u, index_unres_u, max_no_searches,        &
     &     n_calls_spiral_u, n_pts_spiral_u, icode)

      if (icode  >   0)then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 4.2 unable to fill in all unresolved points '
        icode = icode + 2000
        go to 9999
      end if

9999  continue
      return
      END SUBROUTINE read_lsms
!----------------------------------------------------------------------
#endif
