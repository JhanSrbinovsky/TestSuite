#if defined(A32_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Unrotates the model winds if on a rotated grid.
!
! Subroutine Interface:

      Subroutine LBC_Unrotate_Model_Winds (                             &
     &     u                                                            &
     &,    v                                                            &
     &,    u_size                                                       &
     &,    v_size                                                       &
     &,    row_length                                                   &
     &,    rows                                                         &
     &,    rows_v                                                       &
     &,    model_levels                                                 &
     &,    src_halo_type                                                &
     &,    src_pole_lat                                                 &
     &,    src_pole_long                                                &
     &,    src_first_lat                                                &
     &,    src_first_long                                               &
     &,    src_delta_lat                                                &
     &,    src_delta_long                                               &
     & )

      Implicit NONE
!
! Description:
!   Unrotates the model winds if on a rotated grid.
!
! Method:
!   1. Linear Interpolation of u-comp from u-grid to p-grid.
!   2. Linear Interpolation of v-comp from v-grid to p-grid.
!      (u and v now on same grid points)
!   3. Call W_EqToLL to unrotate the winds.
!   4. Linear Interpolation of u-comp from p-grid to u-grid.
!   5. Linear Interpolation of v-comp from p-grid to v-grid.
!      (u and v now back on u and v grids)
!
! Original Author : Dave Robinson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.5    16/10/02  Original code. Dave Robinson
!   6.1    07/12/04  Correct halosize used to calculate lambda_rot
!                    and phi_rot. Dave Robinson
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables :

#include "cmaxsize.h"
#include "parvars.h"

! Subroutine arguments

      Integer :: u_size
      Integer :: v_size
      Integer :: row_length
      Integer :: rows
      Integer :: rows_v
      Integer :: model_levels
      Integer :: src_halo_type

      Real    :: src_delta_lat
      Real    :: src_delta_long
      Real    :: src_first_lat
      Real    :: src_first_long
      Real    :: src_pole_lat
      Real    :: src_pole_long

      Real :: u (u_size, model_levels)
      Real :: v (v_size, model_levels)

! Local parameters:

      Character (Len=*), Parameter ::                                   &
     &                   RoutineName = 'LBC_Unrotate_Model_Winds'

! Local scalars:

      Integer :: ipt, i, j   ! Loop indices
      Integer :: ij          ! Grid point number
      Integer :: points      ! No of grid points
      Integer :: halo_x      ! Halo size in EW direction
      Integer :: halo_y      ! Halo size in NS direction
      Integer :: level       ! Loop index
      Integer :: ErrorStatus ! Error Code

      Character (Len=80) :: CMessage

! Local dynamic arrays:

      Real, dimension (:), allocatable :: coeff3
      Real, dimension (:), allocatable :: coeff4
      Real, dimension (:), allocatable :: lambda_rot
      Real, dimension (:), allocatable :: phi_rot
      Real, dimension (:), allocatable :: lambda_true
      Real, dimension (:), allocatable :: phi_true
      Real, dimension (:), allocatable :: u_p
      Real, dimension (:), allocatable :: v_p
      Real, dimension (:), allocatable :: u_p_rot
      Real, dimension (:), allocatable :: v_p_rot
!
!- End of header

      ErrorStatus = 0
      CMessage = ' '

! --------------------
! Allocate work arrays
! --------------------

      points = g_lasize(1,fld_type_p,src_halo_type,mype) *              &
     &         g_lasize(2,fld_type_p,src_halo_type,mype)

      allocate ( coeff3      (points) )
      allocate ( coeff4      (points) )
      allocate ( lambda_rot  (points) )
      allocate ( phi_rot     (points) )
      allocate ( lambda_true (points) )
      allocate ( phi_true    (points) )

! -------------------
! Get halo dimensions
! -------------------

      halo_x = halosize (1, src_halo_type)
      halo_y = halosize (2, src_halo_type)

! -------------------------------------------
! Set up lat/longs for the rotated model grid
! -------------------------------------------

      ipt = 0
      Do j = 1,g_lasize(2,fld_type_p,src_halo_type,mype)
        Do i = 1,g_lasize(1,fld_type_p,src_halo_type,mype)

          ipt =ipt + 1

          Lambda_Rot(ipt) = Src_First_Long + Src_Delta_Long *           &
     &                     (I + g_datastart(1,mype) - halo_x - 2)

          Phi_Rot(ipt) = Src_First_Lat + Src_Delta_Lat *                &
     &                  (J + g_datastart(2,mype) - halo_y - 2)

        End Do
      End Do

! -------------------------------------
! Get true lat/longs for the model grid
! -------------------------------------

! DEPENDS ON: eqtoll
      Call EqToLL (Phi_Rot, Lambda_Rot, Phi_True, Lambda_True,          &
     &             Src_Pole_Lat, Src_Pole_Long, Points)

! ----------------------------------
! Compute wind rotation coefficients
! ----------------------------------

! DEPENDS ON: w_coeff
      Call W_Coeff (Coeff3, Coeff4, Lambda_True, Lambda_Rot,            &
     &              Src_Pole_Lat, Src_Pole_Long, Points)

      deallocate (lambda_rot)
      deallocate (phi_rot)
      deallocate (lambda_true)
      deallocate (phi_true)

! --------------------------------------
! Allocate work arrays for interpolation
! --------------------------------------

      allocate ( u_p     (points) )
      allocate ( v_p     (points) )
      allocate ( u_p_rot (points) )
      allocate ( v_p_rot (points) )

      Do level = 1, model_levels

! --------------------------------------------------
! Interpolate u from u-grid to p-grid ; u to u_p_rot
! --------------------------------------------------

! DEPENDS ON: lbc_u_to_p
        call lbc_u_to_p (u(1,level), u_p_rot,                           &
     &                   row_length, rows, halo_x, halo_y)

! --------------------------------------------------
! Interpolate v from v-grid to p-grid ; v to v_p_rot
! --------------------------------------------------

! DEPENDS ON: lbc_v_to_p
        call lbc_v_to_p (v(1,level), v_p_rot,                           &
     &                   row_length, rows, rows_v, halo_x, halo_y )

! -----------------------------------------------------
! Rotate winds to standard lat-lon ; u/v_p_rot to u/v_p
! -----------------------------------------------------

! DEPENDS ON: w_eqtoll
        call w_eqtoll(coeff3, coeff4, u_p_rot, v_p_rot, u_p, v_p,       &
     &                points, points)

! ---------------------------------------------------
! Interpolate u from p-grid back to u-grid ; u_p to u
! ---------------------------------------------------

! DEPENDS ON: lbc_p_to_u
        call lbc_p_to_u (u_p, u(1,level),                               &
     &                   row_length, rows, halo_x, halo_y)

! ---------------------------------------------------
! Interpolate v from p-grid back to v_grid ; v_p to v
! ---------------------------------------------------

! DEPENDS ON: lbc_p_to_v
        call lbc_p_to_v (v_p, v(1,level),                               &
     &                   row_length, rows, rows_v, halo_x, halo_y)

      End Do  !  Loop over levels

      deallocate ( u_p )
      deallocate ( v_p )
      deallocate ( u_p_rot )
      deallocate ( v_p_rot )
      deallocate ( coeff3 )
      deallocate ( coeff4 )

      Return
      END SUBROUTINE LBC_Unrotate_Model_Winds
#endif
