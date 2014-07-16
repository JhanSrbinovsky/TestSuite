#if defined(FRAMES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module containing information on LBC variables for Frames

MODULE LBC_Frames_mod

! Description : This module stores the information on each LBC
!               field as required for the FRAMES program.
!
! Owner: Thomas Green
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 

  IMPLICIT NONE

! ------------------------------------------------------------------
! WARNING : Currently set for 13 variables - this will need to
! extended for any new LBC fields used beyond STASH codes 32/013. 
! ------------------------------------------------------------------

  INTEGER, Parameter :: Num_LBC_Vars = 13

  TYPE LBC_Fields

!    For each LBC variable :

     INTEGER :: Model               ! Model Identifier
     INTEGER :: Section             ! Section Number
     INTEGER :: StashCode           ! Stash code (in section 32) 
     INTEGER :: Sect0_StashCode     ! Corresponding Section 0 stash code
     INTEGER :: Fld_Type            ! Field type
     INTEGER :: Halo_Type           ! Halo type
     INTEGER :: N_Levels            ! Number of levels
     INTEGER :: Bot_lev             ! Bottom level number
     INTEGER :: Top_Lev             ! Top level number

  END TYPE LBC_Fields

  TYPE ( LBC_Fields ) :: LBC_Variables (Num_LBC_Vars)

END MODULE LBC_Frames_mod
#endif
