#if defined(A12_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Constants necessary for Coriolis terms in advection scheme.

Module trignometric_Mod

! Description:
! This module is used to hold trignometric values set in atmos_init.
!
! Method:
! The trignometric arrays are calculated in routine control/top_level/
! set_trig called from atmos_init, 
! and used in dynamics_advection and other routines.
!
! Current Code Owner: R.Barnes
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Implicit None

! Arguments

! Trignometric Arrays
! Most names are self-defining with 

Real, Allocatable  ::  cos_theta_latitude (:,:) 
Real, Allocatable  ::  sec_theta_latitude (:,:)
Real, Allocatable  ::  FV_cos_theta_latitude (:,:)
Real, Allocatable  ::  FV_sec_theta_latitude (:,:)
Real, Allocatable  ::  sin_theta_latitude (:,:)
Real, Allocatable  ::  tan_theta_latitude (:,:)
Real, Allocatable  ::  sin_v_latitude (:,:)
Real, Allocatable  ::  tan_v_latitude (:,:)
Real, Allocatable  ::  cos_v_latitude (:,:)
Real, Allocatable  ::  sec_v_latitude (:,:)
Real, Allocatable  ::  cos_theta_longitude (:,:)
Real, Allocatable  ::  sin_theta_longitude (:,:)
Real, Allocatable  ::  cos_u_longitude (:,:)
Real, Allocatable  ::  sin_u_longitude (:,:)
Real, Allocatable  ::  true_latitude (:,:)
Real, Allocatable  ::  true_longitude (:,:)

End Module trignometric_Mod
#endif
