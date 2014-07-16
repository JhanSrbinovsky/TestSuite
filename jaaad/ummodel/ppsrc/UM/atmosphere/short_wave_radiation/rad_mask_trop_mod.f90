
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Constants necessary for radiation mask and trop levels in radiation scheme.

Module rad_mask_trop_Mod

! Description:
! This module is used to hold values set in atmos_init.
!
! Method:
! The arrays are calculated in routine control/top_level/
! set_radmask called from atmos_init, 
! and used in short_ and long_wave_radiation routines.
!
! Current Code Owner: R.Barnes
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Implicit None

! Arguments

! Radiation masking arrays and tropopause levels

Real,    Allocatable  ::  es_space_interp(:,:,:) 
Logical, Allocatable  ::  rad_mask (:,:)
Integer  ::  min_trop_level  ! first reference model level about 700h
Integer  ::  max_trop_level  ! first reference model level about 50hp

End Module rad_mask_trop_Mod
