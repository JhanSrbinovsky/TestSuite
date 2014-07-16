#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE STATION(xx,pos,cellno,bl,orog,stlon,stlat,nstation,    &
     &    day,month,year,clist,nchem,nfill)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : To write station data to file
!-
!-   Inputs  :
!-   Outputs :
!-   Controls:
!-
!-   Created  24-FEB-1997   Bill Collins
!-   Revised  24-NOV-2000 Colin Johnson  Bill's mods + FINDNAME
!-   Revised  12-SEP-2001 Colin Johnson  BL diagnosed with GETMETPOINT
!-
!VVV  V5.0  STATION 12/IX/01 - Bl diagnosed with GETMETPOINT
!     6.2   21/10/05   Replaced GSYNC with SSYNC. P.Selwood
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_OUT
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: nstation
      INTEGER, INTENT(IN) :: month
      INTEGER, INTENT(IN) :: year
      INTEGER, INTENT(IN) :: nchem
      INTEGER, INTENT(IN) :: nfill
      INTEGER, DIMENSION(nclprc), INTENT(IN)  :: cellno
      INTEGER, DIMENSION(numchem), INTENT(IN) :: clist

      REAL, INTENT(IN) :: day
      REAL, DIMENSION(nc,nclprc),     INTENT(IN) :: xx
      REAL, DIMENSION(4,nclprc),      INTENT(IN) :: pos
      REAL, DIMENSION(numstat),       INTENT(IN) :: stlon,stlat
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: bl
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: orog

      INTEGER :: n
      INTEGER :: i
      INTEGER :: j
      INTEGER :: k
      INTEGER :: blcode
      INTEGER :: info

      REAL :: cos_c
      REAL :: sin_c
      REAL :: cos_d
      REAL :: sin_d
      REAL :: cosang
      REAL :: cosmax
      REAL :: zorog
      REAL :: zbl
      REAL :: maxang=1.0 ! Maximum solid angle in degrees
      REAL, DIMENSION(numstat) :: cos_a
      REAL, DIMENSION(numstat) :: sin_a
      REAL, DIMENSION(numstat) :: cos_b
      REAL, DIMENSION(numstat) :: sin_b
      CHARACTER(LEN=14) :: filename

! DEPENDS ON: findname
      CALL FINDNAME('s','m','c',0,0,filename,month,year)

      cosmax=COS(pi_over_180*maxang)
      cos_a=COS(pi_over_180*(stlat+90.))
      sin_a=SIN(pi_over_180*(stlat+90.))
      cos_b=COS(pi_over_180*stlon)
      sin_b=SIN(pi_over_180*stlon)

      DO k=0,nproc-1
        info=GC_SHM_PUT
        CALL GC_SSYNC(nproc,info)
        OPEN(59,FILE=filename,STATUS='UNKNOWN',                         &
     &      POSITION='APPEND')
        IF(mype==k) THEN
          DO n=1,nfill
            cos_c=COS(pi_over_180*pos(2,n))
            sin_c=SIN(pi_over_180*pos(2,n))
            cos_d=COS(pi_over_180*pos(1,n))
            sin_d=SIN(pi_over_180*pos(1,n))
! DEPENDS ON: getmetpoint
            Zorog=GETMETPOINT(pos(:,n),orog,.TRUE.)
! DEPENDS ON: getmetpoint
            Zbl=GETMETPOINT(pos(:,n),bl,.TRUE.)
            DO i=1,nstation
              cosang=cos_c*cos_a(i)+                                    &
     &          sin_c*sin_a(i)*(cos_b(i)*cos_d+sin_b(i)*sin_d)
              IF(cosang>cosmax)THEN
! Indicies for met grids
                blcode=0
                IF (pos(4,n)-Earth_Radius-Zorog <= Zbl) blcode=1
                WRITE(59,'(I3,I6,I5,F8.3,F12.5,I3,I6)') i,year,month,   &
     &               day,pos(3,n),blcode,cellno(n)
                DO j=1,nchem
                  WRITE(59,*) xx(clist(j),n)
                ENDDO
              ENDIF       ! COSANG
            ENDDO        ! NSTATION
          ENDDO         ! NFILL
          CLOSE(59)
        ENDIF         ! mype
      ENDDO          ! K

      END SUBROUTINE STATION
#endif
