#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Data module containing vertical interpolation control variables

Module Rcf_V_Int_Ctl_Mod

! Description:
!  Contains variables used for control of vertical interpolation.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.2   30/11/00   Allow non-extrapolating linear interpolation.
!                    P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Implicit None

Logical                 :: v_int_active
Logical                 :: v_int_active_soil   ! switch for soil interp
Integer                 :: v_int_order

! Orders for vertical interpolation
Integer, Parameter      :: Linear  = 1
Integer, Parameter      :: Linear_NoEx = 2  ! Linear w/out extrapolation
Integer, Parameter      :: Cubic   = 3
Integer, Parameter      :: Quintic = 5

End Module Rcf_V_Int_Ctl_Mod
#endif
