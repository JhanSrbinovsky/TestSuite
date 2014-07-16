#if defined(OCEAN) || defined(UTILIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ code to dimension arrays to read interface data
!
      subroutine o_lbc_sizes( km ,nt, rimwidtho,                        &
     &       numside_rowso, numside_colso,                              &
     &       RIM_LOOKUPSO, n_obdy_t_grd, n_obdy_u_grd )

      implicit none
!
! Description:
! Sets some variables used to dimension ocean LAM arrays.
! Called from DERVSIZE and ADDRLN
!
!L Input Arguments  - all INTENT OUT
      integer km        ! # of vertical levels for tracers, u and v
      integer nt        ! # of tracers (including T & S)
      integer rimwidtho ! rim width of boundary data to be read

!L Output arguments  - all INTENT IN
      integer numside_rowso ! # of active boundary rows
      integer numside_colso ! # of active boundary columns
      integer RIM_LOOKUPSO  ! # of ocean bdy lookup tables
      integer n_obdy_t_grd  ! # of ocean tracer grid bdy fields
      integer n_obdy_u_grd  ! # of ocean velocity grid bdy fields

!L Global variables (only used for input)
#if !defined(RECON)
#include "cntlocn.h"
#endif

!--------------------------------------------------------
! 1. calculate numside_rowso and numside_colso
      numside_rowso = 0
      if ( l_obdy_north ) then
        numside_rowso = numside_rowso + rimwidtho
      end if
      if ( l_obdy_south ) then
        numside_rowso = numside_rowso + rimwidtho
      end if
      numside_colso = 0
      if ( l_obdy_east) then
        numside_colso = numside_colso + rimwidtho
      end if
      if ( l_obdy_west) then
        numside_colso = numside_colso + rimwidtho
      end if

! 2. find the number of tracer grid and velocity grid boundary
!    fields and the number of lookup tables

      RIM_LOOKUPSO = 0
      n_obdy_t_grd = 0
      n_obdy_u_grd = 0

      if ( l_obdy_tracer ) then
        n_obdy_t_grd = n_obdy_t_grd + nt * km
        RIM_LOOKUPSO = RIM_LOOKUPSO + nt
      end if

      if ( l_obdy_uv ) then
        n_obdy_u_grd = n_obdy_u_grd + 2 * km
        RIM_LOOKUPSO = RIM_LOOKUPSO + 2
      end if

      if ( l_obdy_STREAM ) then
        n_obdy_t_grd = n_obdy_t_grd + 2
        RIM_LOOKUPSO = RIM_LOOKUPSO + 2
      end if

      if ( l_obdy_ice ) then
        n_obdy_t_grd = n_obdy_t_grd + 3
        RIM_LOOKUPSO = RIM_LOOKUPSO + 3
      end if

      return
      END SUBROUTINE o_lbc_sizes
#endif
