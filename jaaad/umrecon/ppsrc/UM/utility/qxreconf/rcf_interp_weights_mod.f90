
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Defines interpolation weight variables

Module Rcf_Interp_Weights_Mod

! Description:
! This module defines the interpolation data-types, containing
! weights etc for interpolations.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Implicit None

! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

!-------------------------------------------------------------------
! Horizontal interpolation information
!-------------------------------------------------------------------

! Some parameters used for array sizing for weights
Integer, Parameter        :: Idim = 4
Integer, Parameter        :: Icof = 2

! These are the horizontal interpolation methods that could be
! used by reconfiguration.

Integer, Parameter        :: bilinear          = 1
Integer, Parameter        :: area_weighted     = 2
Integer, Parameter        :: nearest_neighbour = 3
Integer, Parameter        :: unset             = imdi


Integer, Save            :: smcp_int_method = unset

Logical, Save            :: h_int_active  ! Is interpolation on or off?
Integer, Save            :: h_int_method = unset  ! Interpolation method

  ! Weights and indexes etc for bi-linear interp
Integer, Save, Allocatable, Target   :: bl_index_b_l( :, : )
Integer, Save, Allocatable, Target   :: bl_index_b_r( :, : )
Integer, Save, Allocatable, Target   :: bl_index_nearest( : )
Real, Save, Allocatable, Target      :: weight_t_r( :, : )
Real, Save, Allocatable, Target      :: weight_b_r( :, : )
Real, Save, Allocatable, Target      :: weight_t_l( :, : )
Real, Save, Allocatable, Target      :: weight_b_l( :, : )

  ! Weights and indexes etc for area-weighted interp
    ! Fixed size
Real, Save                   :: aw_area_box( idim )

    ! Allocatables
Integer, Save, Allocatable, Target   :: aw_index_targ_lhs( :, : )
Integer, Save, Allocatable, Target   :: aw_index_targ_top( :, : )
Real, Save, Allocatable, Target      :: aw_colat_t( :, : )
Real, Save, Allocatable, Target      :: aw_long_l( :, : )

  ! Rotation coefficients
Real, Save, Allocatable      :: Coeff1( : )
Real, Save, Allocatable      :: Coeff2( : )
Real, Save, Allocatable      :: Coeff3( : )
Real, Save, Allocatable      :: Coeff4( : )

End Module Rcf_Interp_Weights_Mod

