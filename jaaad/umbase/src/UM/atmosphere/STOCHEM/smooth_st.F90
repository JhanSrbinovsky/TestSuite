#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE SMOOTH_ST(CONC,NNN)
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : SMOOTHS HOLES
!-
!-   Inputs  : CONC,NNN,LAT,LONG
!-   Outputs : CONC
!-   Controls:
!-
!-   Created   7-JUN-1995   Bill Collins
!-   Updated   6-AUG-1996   Bill Collins  Parameters now in INCLUDE
!-   Updated  16-JAN-1998   Bill Collins
!-                           Only smooths in longitude since there
!-                           are latitude divisions between PEs
!-   Updated  17-OCT-1997   Bill Collins  Converted to Fortran 90
!-   Updated  10-MAR-1998   Bill Collins
!-                           Reduced size of Eulerian arrays
!-                           no longer smooth in latitude.
!-
!VVV  V5.2.1  SMOOTH_ST 16/VIII/01
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_OUT
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER, DIMENSION(NLONG,NLPE,NLEV), INTENT(IN) :: NNN
      REAL, DIMENSION(NUMCHEM,NLONG,NLPE,NLEV), INTENT(INOUT) :: CONC

      REAL, DIMENSION(0:NLONG/2-1) :: WEIGHT
      REAL, DIMENSION(NUMCHEM) :: SUMCON
      REAL :: DIST,SUMW,LAT1,LAT2,R
      INTEGER :: I,J,K,L,II,JJ,ILON,JLAT,MLONG,N

! For grid squares with no cells, take average of surrounding squares
! that do contain cells

! Radius for smoothing (in degrees subtended at centre of earth)
      R=2.
      N=0
      DO J=1,NLPE
        LAT1=pi_over_180*(LAT(J+ltdat-1)+LAT(J-1+ltdat-1))/2.
        MLONG=INT(1.5/SIN(LAT1)+0.5)
        IF(MLONG>NLONG/2-1) MLONG=NLONG/2-1
        DO II=0,NLONG/2-1
          DIST=COS(LAT1)**2+SIN(LAT1)**2*COS(II*360.*pi_over_180/NLONG)
          IF(DIST>1.0) DIST=1.0
          DIST=ACOS(DIST)/pi_over_180
          WEIGHT(II)=EXP(-DIST**2/R**2)
          IF(WEIGHT(II)<1E-20) WEIGHT(II)=0.
        ENDDO
        DO K=1,NLEV
          DO I=1,NLONG
            SUMW=0.
            SUMCON=0.
! calculate coordinates of surrounding squares
            DO II=-MLONG,MLONG
              ILON=MOD(I+II+NLONG-1,NLONG)+1
              JLAT=J
! Add up the weights
              IF(NNN(ILON,JLAT,K)/=0) THEN
                SUMW=SUMW+WEIGHT(ABS(II))
! concentration is the average
                SUMCON=SUMCON+CONC(:,ILON,JLAT,K)*WEIGHT(ABS(II))
              ENDIF
            ENDDO
            IF(SUMW>0.) THEN
              CONC(:,I,J,K)=SUMCON/SUMW
            ELSE
!              WRITE(6,*) 'SUMW=0. I,J,K=',I,J,K
              N=N+1
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      WRITE(6,*) 'SMOOTH_ST: Grid elements with no cells: ',N

      END SUBROUTINE SMOOTH_ST
#endif
