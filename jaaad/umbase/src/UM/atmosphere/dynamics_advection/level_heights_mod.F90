#if defined(A12_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Constants necessary for model level heights in advection and other schemes.

Module level_heights_Mod

! Description:
! This module is used to hold levels values set in atmos_init.
!
! Method:
! The height levels arrays are calculated in routine control/top_level/
! set_levels called from atmos_init, 
! and used in dynamics_advection and other routines.
!
! Current Code Owner: R.Barnes
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Implicit None

! Arguments

! Heights Arrays

Real, Allocatable  ::  eta_theta_levels(:)    ! eta values of theta levels
Real, Allocatable  ::  eta_rho_levels  (:)    ! eta values of rho levels
Real, Allocatable  ::  r_theta_levels (:,:,:) ! height of theta levels
Real, Allocatable  ::  r_rho_levels   (:,:,:) ! height of rho levels
Real, Allocatable  ::  r_at_u     (:,:,:) ! height at u points on rho levels
Real, Allocatable  ::  r_at_v     (:,:,:) ! height at v points on rho levels

End Module level_heights_Mod
#endif
