#if defined(C70_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+
! Subroutine Interface:
      SUBROUTINE LEVSRT(TYPE,NLEVS,IL,RL)
      IMPLICIT NONE
! Description:
!
! Method:
!
! Current code owner:  S.J.Swarbrick
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Mar. 95    Original code.  S.J.Swarbrick
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project

! Subroutine arguments:

!   Scalar arguments with intent(in):

      CHARACTER*1 TYPE
      INTEGER     NLEVS

!   Array arguments with intent(inout):

      REAL        RL(NLEVS)
      INTEGER     IL(NLEVS)

! Local variables:

      LOGICAL     LSWAP
      INTEGER     I
      INTEGER     J
      INTEGER     ILT
      REAL        RLT

!- End of Header ----------------------------------------------------

      DO 100 I=1,NLEVS
        LSWAP=.FALSE.
        DO 200 J=1,NLEVS-1
          IF(TYPE == 'I') THEN
            IF(IL(J) >  IL(J+1)) THEN
              LSWAP=.TRUE.
              ILT=IL(J)
              IL(J)=IL(J+1)
              IL(J+1)=ILT
            END IF
          ELSE
            IF(RL(J) <  RL(J+1)) THEN
              LSWAP=.TRUE.
              RLT=RL(J)
              RL(J)=RL(J+1)
              RL(J+1)=RLT
            END IF
          END IF
  200   CONTINUE
        IF(.NOT.LSWAP) RETURN
  100 CONTINUE
      RETURN
      END SUBROUTINE LEVSRT
#endif
