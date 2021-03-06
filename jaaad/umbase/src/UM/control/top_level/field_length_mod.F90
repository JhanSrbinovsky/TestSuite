! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module field_length_mod

  implicit none

!  interface
!    function field_length(grid_type,halo_type,levels)
!      integer :: field_length ! function output
!      integer, intent(in) :: grid_type
!      integer, intent(in) :: halo_type
!      integer, intent(in) :: levels
!    end function field_length
!  end interface

      ! the set of grid types that can be specified in a stashmater file
      integer, parameter :: theta_points = 1 !1
      integer, parameter :: theta_points_land_only = 2 !2
      integer, parameter :: theta_points_sea_only = 3 !3
      integer, parameter :: zonal_theta_points = 4 !4
      integer, parameter :: merid_theta_points = 5 !5
      integer, parameter :: uv_points = 6 !11
      integer, parameter :: uv_points_land_only = 7 !12
      integer, parameter :: uv_points_sea_only = 8 !13
      integer, parameter :: zonal_uv_points = 9 !14
      integer, parameter :: merid_uv_points = 10 !15
      integer, parameter :: scalar_points = 11 !17
      integer, parameter :: u_points = 12 !18
      integer, parameter :: v_points = 13 !19
      integer, parameter :: land_points = 14 !21
      integer, parameter :: ozone_points = 15 !22
      integer, parameter :: river_points = 16 !23
      integer, parameter :: lbc_points = 17 !25
      integer, parameter :: lbc_theta_points = 18 !26
      integer, parameter :: lbc_u_points = 19 !27
      integer, parameter :: lbc_v_points = 20 !28
      integer, parameter :: lbc_orography_points = 21 !29
      integer, parameter :: non_standard_ocn_points = 22 !30
      integer, parameter :: comp_ocn_mass_points = 23 !31
      integer, parameter :: comp_ocn_velocity_points = 24 !32
      integer, parameter :: cyclic_ocn_mass_points = 25 !36
      integer, parameter :: cyclic_ocn_velocity_points = 26 !37
      integer, parameter :: cyclic_ocn_u_points = 27 !38
      integer, parameter :: cyclic_ocn_v_points = 28 !39
      integer, parameter :: ocn_mass_points = 29 !41
      integer, parameter :: ocn_velocity_points = 30 !42
      integer, parameter :: zonal_ocn_mass_points = 31 !43
      integer, parameter :: zonal_ocn_velocity_points = 32 !44
      integer, parameter :: merid_ocn_mass_points = 33 !45
      integer, parameter :: merid_ocn_velocity_points = 34 !46
      integer, parameter :: scalar_ocn_points = 35 !47
      integer, parameter :: lbc_ocn_points = 36 !51
      integer, parameter :: lbc_ocn_t_points = 37 !52
      integer, parameter :: lbc_ocn_uv_points = 38 !53
      integer, parameter :: wave_points = 39 !60
      integer, parameter :: wave_points_sea_only = 40 !62
      integer, parameter :: lbc_wave_points = 41 !65

      ! the set of halo types that can be specified in a stashmater file
      integer, parameter :: single_halo = 1
      integer, parameter :: extended_halo = 2
      integer, parameter :: no_halo = 3

#if defined(ATMOS)

contains

   integer function field_length(grid_type,halo_type,levels)

      implicit none

#include <parvars.h>
#include <typsize.h>
#include <typptra.h>

! Function arguments:
      integer, intent(in) :: grid_type
      integer, intent(in) :: halo_type
      integer, intent(in) :: levels

! Local variables:
      integer :: length 
      integer :: x, y, halo_x, halo_y
 
      select case (grid_type)
        case (theta_points                                              &
     &        ,theta_points_land_only ,theta_points_sea_only            &
     &        ,uv_points ,uv_points_land_only ,uv_points_sea_only       &
     &        ,ozone_points                                             &
     &       )
          x = row_length
          y = rows
        case (u_points)
          x = row_length
          y = rows
        case (v_points)
          x = row_length
          y = n_rows
        case (land_points)
          ! this is a special case (assume there are no halos)
          if (halo_type .ne. no_halo) then
            write(6,*) "WARNING MSG FROM FIELD_LENGTH()"
            write(6,*) "POSSIBLE ERROR DETERMINING FIELD LENGTH"
            write(6,*) "UNHANDLED HALO_TYPE (FOR LAND_FIELDS):",halo_type
            write(6,*) "ASSUMING 'no_halo'"
            write(6,*) "END WARNING MSG"       
          endif
          length = land_field * levels
          field_length = length
          return          
        case (river_points)
          x = river_row_length
          y = river_rows
        !! EVENTUALLY NEED TO ADD OCEAN, WAVE FIELDS
        case default !! unspecified/unhandled grid_type 
          write(6,*) "WARNING MSG FROM FIELD_LENGTH()"
          write(6,*) "POSSIBLE ERROR DETERMINING FIELD LENGTH"
          write(6,*) "UNHANDLED GRID_TYPE: ",grid_type
          write(6,*) "ASSUMING (ROW_LENGTH*ROWS)" 
          write(6,*) "END WARNING MSG"
          x = row_length
          y = rows
      end select

      select case (halo_type)
        case (extended_halo)
          halo_x = halo_i
          halo_y = halo_j
        case (single_halo)
          halo_x = offx
          halo_y = offy
        case default ! (no_halo)
          halo_x = 0
          halo_y = 0
      end select

      length =  (((2 * halo_x) + x) * ((2 * halo_y) + y)) * levels
      field_length = length

      end function field_length

!- End of Function code

#endif

end module field_length_mod
