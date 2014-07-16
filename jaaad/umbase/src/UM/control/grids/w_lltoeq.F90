#if defined(C92_2A) || defined(MAKEBC) || defined(RECON)               \
 || defined(VAROPSVER)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine W_LLTOEQ----------------------------------------------
!LL
!LL  Purpose:  Calculates u and v components of wind on equatorial
!LL            (eq) latitude longitude grid by rotating wind
!LL            components on standard latitude-longitude (eq)
!LL            grid.
!LL
!LL  Written by A. Dickinson
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL  4.1   31/05/96     The number of v points to be processed on a
!LL                     C grid differs by row_length. u,v therefore
!LL                     treated separately.
!LL                     Author I.Edmond       Reviewer D. Goddard
!LL
!LL Programming standard :
!LL
!LL Logical components covered : S134
!LL
!LL Project task :
!LL
!LL External documentation:
!LL
!LL  Documentation: The transformation formulae are described in
!LL                 unified model on-line documentation paper S1.
!LL
!LLEND -----------------------------------------------------------------
!
!*L  Arguments:--------------------------------------------------------
      SUBROUTINE W_LLTOEQ(COEFF1,COEFF2,U,V,U_EQ,V_EQ,POINTS,POINTS2)

      IMPLICIT NONE

      INTEGER                                                           &
     & POINTS                                                           &
                         !IN  Number of points to be processed
     &,POINTS2    ! IN  Number of v points to be processed

      REAL                                                              &
     & COEFF1(POINTS)                                                   &
                         !IN  Coefficient of rotation no 1
     &,COEFF2(POINTS)                                                   &
                         !IN  Coefficient of rotation no 2
     &,U_EQ(POINTS)                                                     &
                         !OUT u component of wind on equatorial grid
     &,V_EQ(POINTS)                                                     &
                         !OUT v component of wind on equatorial grid
     &,U(POINTS)                                                        &
                         !IN  u component of wind on lat-lon grid
     &,V(POINTS)         !IN  v component of wind on lat-lon grid
! Workspace usage:-----------------------------------------------------
! None
#include "c_mdi.h"
!----------------------------------------------------------------------
! External subroutines called:-----------------------------------------
! None
!*---------------------------------------------------------------------
! Define local varables:-----------------------------------------------
      INTEGER I
!----------------------------------------------------------------------

!L 1. Transform wind components
!
! Formulae used are from eq (4.13)

      DO 100 I = 1,POINTS
      IF ( U(I)  ==  RMDI .OR. V(I)  ==  RMDI ) THEN
        U_EQ(I) = RMDI
      ELSE
       U_EQ(I)=COEFF1(I)*U(I)-COEFF2(I)*V(I)
      END IF
100   CONTINUE
      ! On a C grid number of u,v points processed differs by row_length
      DO I = 1,POINTS2
      IF ( U(I)  ==  RMDI .OR. V(I)  ==  RMDI ) THEN
        V_EQ(I) = RMDI
      ELSE
       V_EQ(I)=COEFF1(I)*V(I)+COEFF2(I)*U(I)
      END IF
      END DO

      RETURN
      END SUBROUTINE W_LLTOEQ
#endif
