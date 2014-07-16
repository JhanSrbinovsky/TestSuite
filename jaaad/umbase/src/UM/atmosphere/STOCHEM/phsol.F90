#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE PHSOL(kh,ke,ll,tc,m,cthno3,ctso2,ctnh3,so4,hp,l_crit)
!-----------------------------------------------------------------------
!
!   Purpose and Methods:
!     To solve H+ by bissection routine, the method is one of trial and
!     error search to find a minima in the modulus of the imposed minus
!     the returned value of [H+], then use bissection to obtain the
!     value to the supplied tolerance.
!
!   Inputs:    KH,KE,LL,M,TC,CTs
!
!   Outputs:   HP
!
!
! Current Owner of Code: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  4.4    18/03/97  Created.  C.E. Johnson
!  5.5    03/11/03  Partial vectorisation of code. M.G. Sanderson
!  5.5    12/03/04  Final code vectorisation. K. Ketelsen
!  6.1    20/10/04  No change.
!
!VVV  V2.2  PH3 15/III/00
!-----------------------------------------------------------------------
      USE IN_STOCHEM_CHM
      IMPLICIT NONE
!-----------------------------------------------------------------------
      REAL, INTENT(IN) :: l_crit
      REAL, DIMENSION(chunk), INTENT(IN) :: ll,tc,cthno3,ctso2,ctnh3,   &
     &  m,so4
      REAL, DIMENSION(6,chunk), INTENT(IN) :: kh  ! Henry's Law coeffs
      REAL, DIMENSION(6,chunk), INTENT(IN) :: ke  ! Equilibrium constant
      REAL, DIMENSION(chunk),  INTENT(OUT) :: hp  ! [H+]

      INTEGER :: i, j                             ! Loop counts
      INTEGER :: maxcount = 20                    ! Max number iteration
      REAL, DIMENSION(20,chunk) :: x, y
      REAL, DIMENSION(5,chunk) :: x2, y2
      REAL :: hp1
      REAL :: xtol = 0.001                        ! Fractional tolerance
      REAL :: half = 1.0 / 2.0                    ! One half
      REAL :: fifth = 1.0 / 5.0                   ! One fifth

      LOGICAL, DIMENSION(chunk) :: todo

! Set at largest expected [H+] *25:
      hp = 1.0

      todo = (ll > l_crit)

      DO i = 1, maxcount
        DO j = 1, chunk

! Find first minima.
          IF (todo(j)) THEN
            hp(j) = hp(j) * fifth ! decrease HP by a factor of 5 each ti
            x(i,j) = hp(j)
! DEPENDS ON: phcalc
            CALL PHCALC(x(i,j),hp1,cthno3(j),ctso2(j),ctnh3(j),so4(j),  &
     &        ll(j),m(j),tc(j),kh(:,j),ke(:,j))
            y(i,j) = ABS(x(i,j) - hp1)
            IF (i > 2) THEN
              IF (y(i,j) > y(i-1,j)) THEN
                x2(1:3,j) = x(i-2:i,j)
                y2(1:3,j) = y(i-2:i,j)
                todo(j) = .FALSE.
              END IF
            END IF
          END IF
        END DO
        IF (.NOT.ANY(todo)) EXIT
        IF (i > maxcount) THEN
          WRITE(6,*) '**** Error in PHSOL ****'
          STOP 'error'
        END IF
      END DO

      todo = (ll > l_crit)
      WHERE (y2(2,:) > y2(1,:) .OR. y2(2,:) > y2(3,:))
        hp = 1.0e-6    ! In case of trouble.
        todo = .FALSE.
      END WHERE

! Iterate bisection until the interval is within the tolerance.
      DO i = 1, maxcount
        IF (.NOT.ANY(todo)) EXIT
        DO j = 1, chunk
          IF (todo(j)) THEN
            x2(4,j) = (x2(1,j) + x2(2,j)) * half
! DEPENDS ON: phcalc
            CALL PHCALC(x2(4,j),hp1,cthno3(j),ctso2(j),ctnh3(j),        &
     &        so4(j),ll(j),m(j),tc(j),kh(:,j),ke(:,j))
            y2(4,j) = abs(x2(4,j) - hp1)
            x2(5,j) = (x2(3,j) + x2(2,j)) * half
! DEPENDS ON: phcalc
            CALL PHCALC(x2(5,j),hp1,cthno3(j),ctso2(j),ctnh3(j),        &
     &        so4(j),ll(j),m(j),tc(j),kh(:,j),ke(:,j))
            y2(5,j) = ABS(x2(5,j) - hp1)
            IF (y2(2,j) < y2(4,j) .AND. y2(2,j) < y2(5,j)) THEN
              x2(1,j) = x2(4,j)
              x2(3,j) = x2(5,j)
              y2(1,j) = y2(4,j)
              y2(3,j) = y2(5,j)
            ELSE IF (y2(4,j) < y2(5,j)) THEN
              x2(3,j) = x2(2,j)
              x2(2,j) = x2(4,j)
              y2(3,j) = y2(2,j)
              y2(2,j) = y2(4,j)
              IF ((x2(1,j) - x2(2,j)) / x2(2,j) < xtol) THEN
                todo(j) = .FALSE.
              END IF
            ELSE
              x2(1,j) = x2(2,j)
              x2(2,j) = x2(5,j)
              y2(1,j) = y2(2,j)
              y2(2,j) = y2(5,j)
              IF ((x2(1,j) - x2(2,j)) / x2(2,j) < xtol) THEN
                todo(j) = .FALSE.
              END IF
            END IF
            hp(j) = x2(2,j)
          END IF
        END DO
      END DO
!
      END SUBROUTINE PHSOL
#endif
