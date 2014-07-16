#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Deallocate the SW extinction diagnostics
!
! Current Code Owner: James Manners
!
! Description of Code:
!    FORTRAN 95
!
!-----------------------------------------------------------------------
!
Subroutine deallocate_swdiag_ext(j_sw)

  Use sw_diag_mod
  Implicit None

  Integer, Intent(in) :: j_sw              ! call to SW radiation

! Extinction diagnostics

  If (Associated(SW_diag(j_sw)%cloud_extinction)) &
    Deallocate(SW_diag(j_sw)%cloud_extinction)
  If (Associated(SW_diag(j_sw)%cloud_weight_extinction)) &
    Deallocate(SW_diag(j_sw)%cloud_weight_extinction)
  If (Associated(SW_diag(j_sw)%ls_cloud_extinction)) &
    Deallocate(SW_diag(j_sw)%ls_cloud_extinction)
  If (Associated(SW_diag(j_sw)%ls_cloud_weight_extinction)) &
    Deallocate(SW_diag(j_sw)%ls_cloud_weight_extinction)
  If (Associated(SW_diag(j_sw)%cnv_cloud_extinction)) &
    Deallocate(SW_diag(j_sw)%cnv_cloud_extinction)
  If (Associated(SW_diag(j_sw)%cnv_cloud_weight_extinction)) &
    Deallocate(SW_diag(j_sw)%cnv_cloud_weight_extinction)

End Subroutine deallocate_swdiag_ext
#endif
