#if defined(ATMOS)
! History:
! Version  Date    Comment
!  3.4   18/05/94  Add new field in A_SPCON. J F Thomson.
!  5.0   02/06/99  Add new fields for C-P C-grid dynamics. M L Gallani.
!  5.3   01/10/01  Add fields for chequerboard radiation. S Cusack
! 6.1  04/08/04  Add extra arrays diff_coeff_u, diff_coeff_v
!                                             Terry Davies
! 6.2  05/01/06  Add true_latitude. Yongming Tang
! 6.2  14/12/06  Add separate arrays VarRes Array co-ordinates
!                                                       Terry Davies

! artcona.h contained constants for the atmosphere.
! As of vn6.6 these constants have moved to a set of modules:
! LEVEL_HEIGHTS_MOD, TRIGNOMETRIC_MOD, DYN_CORIOLIS_MOD, DYN_VAR_RES_MOD,
! DIFF_COEFF_MOD, AD_MASK_TROP_MOD, ROT_COEFF_MOD

#endif
