#if defined(CONTROL) || defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Convert lat/long specification to row/column specification
! Subroutine Interface:


!- End of subroutine code ----------------------------------------


!+Use subroutine LLTOEQ to convert to equatorial lat/lon for LAM
! Subroutine Interface:

      SUBROUTINE LLTOLL(LNTHLL,LSTHLL,LESTLL,LWSTLL,                    &
     &                  PHI_POLE,LAMBDA_POLE)

! Description:
!
! Method:
!
! Current code owner:  S.J.Swarbrick
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Apr. 95    Original code.  S.J.Swarbrick
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project

!  Subroutine arguments:

!    Scalar arguments with intent(in):
      REAL    PHI_POLE    !  Latitude of pole in equatorial system
      REAL    LAMBDA_POLE !  Longitude do.

!    Scalar arguments with intent(inout):
      INTEGER LNTHLL
      INTEGER LSTHLL
      INTEGER LESTLL
      INTEGER LWSTLL

!  Local parameters:
      INTEGER   POINTS
      PARAMETER(POINTS=9)

!  Local arrays:
      REAL PHI      (POINTS)
      REAL LAMBDA   (POINTS)
      REAL LAMBDA_EQ(POINTS)
      REAL PHI_EQ   (POINTS)

!- End of Header -----------------------------------------------------


      PHI(1)=LNTHLL
      PHI(2)=LNTHLL
      PHI(3)=LNTHLL
      PHI(4)=(LNTHLL+LSTHLL)/2
      PHI(5)=(LNTHLL+LSTHLL)/2
      PHI(6)=LSTHLL
      PHI(7)=LSTHLL
      PHI(8)=LSTHLL
      PHI(9)=(LNTHLL+LSTHLL)/2
      LAMBDA(1)=LWSTLL
      IF(LWSTLL <  LESTLL) THEN
        LAMBDA(2)=(LWSTLL+LESTLL)/2
      ELSE
        LAMBDA(2)=(LWSTLL+LESTLL-360)/2
        IF(LAMBDA(2) <  0) LAMBDA(2)=LAMBDA(2)+360
      END IF
      LAMBDA(3)=LESTLL
      LAMBDA(4)=LWSTLL
      LAMBDA(5)=LESTLL
      LAMBDA(6)=LWSTLL
      LAMBDA(7)=LAMBDA(2)
      LAMBDA(8)=LESTLL
      LAMBDA(9)=LAMBDA(2)

! DEPENDS ON: lltoeq
      CALL LLTOEQ                                                       &
     &       (PHI,LAMBDA,PHI_EQ,LAMBDA_EQ,PHI_POLE,LAMBDA_POLE,POINTS)

      IF(LAMBDA_EQ(3) <  LAMBDA_EQ(2)) LAMBDA_EQ(2)=LAMBDA_EQ(2)-360.
      IF(LAMBDA_EQ(2) <  LAMBDA_EQ(1)) LAMBDA_EQ(1)=LAMBDA_EQ(1)-360.
      IF(LAMBDA_EQ(5) <  LAMBDA_EQ(9)) LAMBDA_EQ(9)=LAMBDA_EQ(9)-360.
      IF(LAMBDA_EQ(9) <  LAMBDA_EQ(4)) LAMBDA_EQ(4)=LAMBDA_EQ(4)-360.
      IF(LAMBDA_EQ(8) <  LAMBDA_EQ(7)) LAMBDA_EQ(7)=LAMBDA_EQ(7)-360.
      IF(LAMBDA_EQ(7) <  LAMBDA_EQ(6)) LAMBDA_EQ(6)=LAMBDA_EQ(6)-360.

      RNTHLL=PHI_EQ(1)
      RSTHLL=PHI_EQ(1)
      RESTLL=LAMBDA_EQ(1)
      RWSTLL=LAMBDA_EQ(1)

      DO I=2,8
        RNTHLL=MAX(RNTHLL,PHI_EQ(I))
        RSTHLL=MIN(RSTHLL,PHI_EQ(I))
        RWSTLL=MIN(RWSTLL,LAMBDA_EQ(I))
        RESTLL=MAX(RESTLL,LAMBDA_EQ(I))
      END DO

      IF(RWSTLL <  0)  RWSTLL=RWSTLL+360.
      IF(RESTLL <  0)  RESTLL=RESTLL+360.
      RNTHLL=RNTHLL+.999
      LNTHLL=RNTHLL
      LNTHLL=MIN(90,LNTHLL)
      LSTHLL=RSTHLL
      LWSTLL=RWSTLL
      RESTLL=RESTLL+.999
      LESTLL=RESTLL
      LESTLL=MIN(360,LESTLL)

      RETURN
      END SUBROUTINE LLTOLL

!- End of subroutine code ---------------------------------------------
#endif
