#if defined(A01_3C) || defined(A01_3Z) \
 || defined(A02_3C) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
Module sw_diag_mod

  Use max_calls
  Use swrdiag_mod

  Implicit None

!   Declaration of SW diagnostics.

!   Those calculated within the radiation code and not used to
!   advance the integration in any way are contained within
!   a structure defined in swrdiag_mod.

    Type (StrSWDiag), Dimension(npd_swcall), save :: SW_diag  

End Module sw_diag_mod
#endif
