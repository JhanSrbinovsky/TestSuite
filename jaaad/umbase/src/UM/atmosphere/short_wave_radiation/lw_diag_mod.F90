#if defined(A01_3C) || defined(A01_3Z) \
 || defined(A02_3C) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
Module lw_diag_mod

  Use max_calls
  Use lwrdiag_mod

  Implicit None

!     Declaration of LW diagnostics.
!     Those quantities which are purely diagnostic are bundled into
!     a structure.

      Type (StrLWDiag), Dimension(npd_lwcall), save :: LW_diag

End Module lw_diag_mod
#endif
