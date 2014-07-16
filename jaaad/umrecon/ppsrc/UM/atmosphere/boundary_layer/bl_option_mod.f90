! *****************************COPYRIGHT********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT********************************
!
!+ Data module for switches/options concerned with the BL scheme.

MODULE bl_option_mod

! Description:
!   Module containing runtime options/data used by the boundary 
!   layer scheme.
!
! Method:
!   Switches and associated data values used by the boundary layer 
!   scheme are defined here and assigned default values. These may 
!   be overridden by namelist input.
!
!   Any routine wishing to use these options may do so with the 'USE'
!   statement.
!
! Current Code Owner: Adrian Lock
!
! Code Description:
!   Language: FORTRAN 90
!
! Declarations:

  IMPLICIT NONE

!======================================================================
! Real values set from RUN_BL 
!======================================================================
  REAL :: WeightLouisToLong = 0.0
!                                     ! Weighting of the Louis tail
!                                     ! towards long tails:
!                                     ! 0.0 = Louis tails
!                                     ! 1.0 = Long tails
!
END MODULE bl_option_mod
