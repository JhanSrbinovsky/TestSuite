#if defined(A12_2A)
!+ Interface to random number routine. Returns same random numbers on
! all platforms. Based on the VAR routine of the same name.

SUBROUTINE Var_RandomNumber(IntRand, Out,gi,gj,gk)

! Description:
!   Generate random numbers that are consistent across all machine
!
! Method:
!   TBD
!
! Owner: Manager of Operational Data Assimilation
!
! History:
! Date     Ticket Comment
! -------- ------ -------
! 16/11/05    319 Imported to FCM - previous history available on the web
!             319 Stephen Oxley
! -------- ------ End History
!
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
!
! Code Description:
!   Language:           Fortran 95.
!   Software Standards: VTDP 3.
!
! Location: VarMod_PFInfo

IMPLICIT NONE

INTEGER :: m = 86436
INTEGER :: a = 1093
INTEGER :: c = 18257
Integer, INTENT(INOUT) :: IntRand
Integer, INTENT(IN) :: gi,gj,gk

! Subroutine arguments

REAL, INTENT(OUT) :: Out(gi,gj,gk)

! Local variables

INTEGER :: i
INTEGER :: j
INTEGER :: k
REAL    :: rm


rm=REAL(m)

DO i=1,gi
  DO j=1,gj
    DO k=1,gk
      IntRand=MOD(IntRand*a+c,m)
      Out(i,j,k)=REAL(IntRand)/rm
    END DO
  END DO
END DO

RETURN
END SUBROUTINE Var_RandomNumber

#endif
