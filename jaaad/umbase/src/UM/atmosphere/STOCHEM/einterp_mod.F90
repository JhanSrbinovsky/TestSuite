#if defined(A25_1A)
      MODULE EINTERP_MOD

      PRIVATE

      INTERFACE EINTERP
        MODULE PROCEDURE EINTERP
      END INTERFACE

      PUBLIC EINTERP

      CONTAINS
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE EINTERP(nn,l,k,p,et,pp,x)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Returns single column temperature profile on
!-                         stochem vertical grid and surface temperature
!-                         and pressure at centre of Eulerian xy grid.
!-
!-   Inputs  : L,K,T0,T,P0,P_TH
!-   Outputs : ET,PP
!-   Controls:
!-
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.0  11/10/95   Created. C.E. Johnson
!  4.4  07/07/97   Converted to F90 and parallelised. W.J. Collins
!  5.5  23/02/04   Now defined as a module procedure. K. Ketelsen
!  6.1  24/09/04   Minor tidying of code. M.G. Sanderson
!
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_INTF
      USE INTERP_MOD
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: nn
      INTEGER, DIMENSION(nn),                   INTENT(IN) :: k
      INTEGER, DIMENSION(nn),                   INTENT(IN) :: l
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: p
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN)  :: x
      REAL, DIMENSION(nn,0:nlev),               INTENT(OUT) :: et
      REAL, DIMENSION(nn,0:nlev),               INTENT(OUT) :: pp

      INTEGER :: i,m

      REAL, DIMENSION(4,nn)                       :: rpos
      LOGICAL, DIMENSION(nn)                      :: todo
      REAL, DIMENSION(nn)                         :: et1, pp1

      todo      = .true.
      rpos(4,:) = 0.0

!kk   X-array set outside subroutine
! Calculate position in middle of Stochem Eulerian grid.
      DO i=1,nn
        IF (l(i) <nlong) THEN
          rpos(1,i) = (long(l(i))+long(l(i)+1)) / 2.0
        ELSE
          rpos(1,i) = (long(l(i))+long(MOD(l(i),nlong)+1)+360.) / 2.0
        END IF
        rpos(2,i) = (lat(k(i)-1)+lat(k(i))) / 2.0
      END DO

! Compose array from T and T0, then interpolate to get ET profile
!kk Use Vector Version of INTERP
      DO m = 0,nlev
        DO i=1,nn
          rpos(3,i) = eta_stochem(m)
        END DO
        CALL INTERP_V(nn,todo,rpos,x,et1,.TRUE.,.FALSE.)
        et(:,m) = et1
      ENDDO

! Compose array from P_TH and P0, then interpolate to get PP profile
!kk Use Vector Version of INTERP
      DO m = 0,nlev
        DO i=1,nn
          rpos(3,i) = eta_stochem(m)
        END DO
        CALL INTERP_V(nn,todo,rpos,p,pp1,.TRUE.,.FALSE.)
        pp(:,m) = pp1
      ENDDO

      END SUBROUTINE EINTERP

      END MODULE EINTERP_MOD
#endif
