
! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE SHALLOW_GRAD_STRESS------------------------------------
!LL
!LL  PURPOSE:  CALCULATES THE GRADIENT COMPONENT OF THE STRESS
!LL            DUE TO SHALLOW CONVECTION
!LL
!LL  CALLED FROM: CONVECT
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL
!LL  MODEL            MODIFICATION HISTORY:
!LL VERSION  DATE
!LL
!LL   5.4  7/8/2002   New deck for revised convection scheme 4A
!LL
!LL                                     A.L.M. Grant
!     5.5  20/02/03   Replaced #ENDIF with #endif.      P.Dando
!     6.0   05/08/03   NEC optimisation - replace function F_W by MIN
!                      & rewrite UW,VW loop.  R Barnes.
!     6.2  02/09/05   Part of version 5A. R A Stratton
!
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL
!LL  DOCUMENTATION : http://hc0500/~hadag/cmt_param.ps.gz
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!

      REAL FUNCTION F_W(X)   ! ENSEMBLE VERTICAL VELOCITY PROFILE
      IMPLICIT NONE
!
! Description:
!   Calculates non-dimenisonal vertical velocity profile for shallow cum
!
! Method:
!   Derived fromt large-eddy simulations.
!
! Current Code Owner: A. Grant
!
! History:
! Version    Date     Comment
! -------    ----     -------
! Vn5.4    02/09/02   Original code. <A. Grant>
! Vn5.5    02/05/03   Remove declaration of F_W for nag compilation.
!                     (M.Hughes)
!
! Code Description:
!   Language: FORTRAN 77 + some CRAY extensions
!   This code is written to UMDP3 v6 programming standards.
!
       REAL                                                             &


     &     X             ! non-dimesional height in cumulus layer.
      IF(X <= 1.0) THEN
       F_W=6.0*X
      ELSE
       F_W=6.0
      ENDIF
      RETURN
      END FUNCTION F_W
