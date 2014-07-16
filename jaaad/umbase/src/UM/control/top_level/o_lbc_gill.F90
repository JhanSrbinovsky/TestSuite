#if defined(OCEAN)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculate Gill open boundary condition
!
      subroutine o_lbc_gill(ncols,nrows,n_cols_bdy,n_rows_bdy,          &
     &     levn,field,                                                  &
     &     levels_dataset,xgrid,ygrid,f_type,                           &
     &     l_obdy_north,l_obdy_east,l_obdy_south,l_obdy_west,           &
     &     bdy_n, bdy_e, bdy_s, bdy_w,                                  &
     &     icode,cmessage)

      implicit none
!
!
! Purpose: This routine creates open north and/or south
!          boundaries using the implicit Gill scheme.
!          The scheme was devised for use in tropical ocean
!          models. See UMDP No. 48 for details on the Gill
!          scheme.
!
! Author : Catherine Jones
!
! History
! Model    Date     Modification history from model version 4.5
! version
!  5.2   31/07/00   Adjust definitions of fld_t etc. for consistency
!                   with o_lbc_frs
!  4.5   17/6/98    New subroutine/deck. C.Jones,S.Ineson
!
!
! Subroutine Arguments:
      integer ncols           ! IN number of columns
      integer nrows           ! IN number of rows
      integer n_cols_bdy      ! IN number of cols in boundary data
      integer n_rows_bdy      ! IN number of rows in boundary data
      integer levn            ! IN present level number

      real field(ncols,nrows)       ! IN/OUT field to be updated
      real levels_dataset(ncols,nrows) ! IN levels dataset for this
                                       ! field type

      real xgrid(ncols)              ! IN grid spacing along row
      real ygrid(nrows)              ! IN grid spacing along columns

      integer f_type                 ! IN field type indicator

      logical l_obdy_north           ! IN T=> update north lateral bdy
      logical l_obdy_east            ! IN T=> update east lateral bdy
      logical l_obdy_south           ! IN T=> update south lateral bdy
      logical l_obdy_west            ! IN T=> update west lateral bdy

      real bdy_n(n_cols_bdy)  ! IN north part of boundary data
      real bdy_e(n_rows_bdy)  ! IN east part of boundary data
      real bdy_s(n_cols_bdy)  ! IN south part of boundary data
      real bdy_w(n_rows_bdy)  ! IN west part of boundary data

      integer icode             ! OUT Return code
      character*(80) cmessage   ! OUT Error message

!  Global variables
! for obdy_gill_mu and obdy_gill_lamda
#include "umscalar.h"
#include "parvars.h"
#include "cocnindx.h"

#include "decomptp.h"
#include "decompdb.h"

#if !defined(MPP)
      integer Offx, Offy
      logical atright, atleft, atbase, attop
#endif

! Local parameters
      INTEGER                                                           &
     &   fld_t                                                          &
                         ! indicates a tracer
     &,  fld_u                                                          &
                         ! indicates a velocity
     &,  fld_sf          ! indicates a streamfn
      PARAMETER (fld_t=1, fld_u=2, fld_sf=3)

! Local scalars

      integer i             ! local loop counter
      integer j             ! local loop counter

      real x                ! local constant
      real y                ! local constant
      real z                ! local constant

!-------------------------------------------------------------------
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

! 1.0 For tracer data update boundaries using the Gill scheme

      if (f_type  ==  fld_t .or. f_type  ==  fld_sf ) then

! 1.1 Northern boundary

        if (l_obdy_north .and.at_extremity(PNorth)) then

         x=1.0 - (1.0 - obdy_gill_mu) * obdy_gill_lamda * ygrid(j_jmt)
         y=1.0 + (1.0 - obdy_gill_mu) * obdy_gill_lamda * ygrid(j_jmt)
         z=2.0 * obdy_gill_lamda * ygrid(j_jmt)

           do i = 1, ncols

            field(i,j_jmt) =                                            &
     &       (1./x)*(field(i,j_jmt-2)*y +                               &
     &       z*(obdy_gill_mu*field(i,j_jmt-1) -                         &
     &       bdy_n(i)))


            if (levels_dataset(i,j_jmt-1)  <   levn) then

              field(i,j_jmt) = 0.0          ! this is a land point

            endif

           enddo

         endif          ! Northern boundary

! 1.2 Southern boundary


         if (l_obdy_south .and. at_extremity(PSouth)) then

           x = 1.0 - ( 1.0-obdy_gill_mu )*obdy_gill_lamda*ygrid(j_1)
           y = 1.0 + ( 1.0-obdy_gill_mu )*obdy_gill_lamda*ygrid(j_1)
           z = 2.0 * obdy_gill_lamda * ygrid(j_1)

           do i = 1, ncols

            field(i,j_1) = (1./x)*(field(i,j_1+2)*y +                   &
     &                   z*(obdy_gill_mu*field(i,j_1+1) - bdy_s(i)))


            if (levels_dataset(i,j_1+1)  <   levn) then

              field(i,j_1) = 0.0          ! this is a land point

            endif

           enddo

         endif          ! southern boundary

         if (l_obdy_east .or. l_obdy_east) then
           icode=1
           write(6,*)' OLBCGIL1: E/W code not implemented for GILL'
           cmessage=' OLBCGIL1:  E/W code not implemented for GILL'
         endif

      endif                  ! updating tracer fields

      return

      END SUBROUTINE o_lbc_gill
#endif
