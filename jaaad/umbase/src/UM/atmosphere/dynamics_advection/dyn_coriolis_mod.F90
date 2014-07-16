#if defined(A12_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Constants necessary for Coriolis terms in advection scheme.

Module dyn_coriolis_Mod

! Description:
! This module is used to hold Coriolis values set in atmos_init.
!
! Method:
! The f1/2/3 arrays are calculated in routine control/top_level/
! set_coriolis called from atmos_init, 
! and used in dynamics_advection routines.
!
! Current Code Owner: R.Barnes
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Implicit None

! Arguments

! Coriolis Arrays

Real, Allocatable  ::  f1_at_v (:,:)  ! components of coriolis force
Real, Allocatable  ::  f2_at_u (:,:)
Real, Allocatable  ::  f3_at_u (:,:)
Real, Allocatable  ::  f3_at_v (:,:)

End Module dyn_coriolis_Mod
#endif
