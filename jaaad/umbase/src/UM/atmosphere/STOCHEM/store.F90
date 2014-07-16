#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE STORE(estore,emiss,nnn,nbl,time,month)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : STORE EMISSIONS IF NO CELLS ARE PRESENT
!-
!-   Inputs  : EMISS,NNN
!-   Outputs : ESTORE
!-   Controls:
!-
!-   Created   9-DEC-1993   W.J. Collins
!-   Updated   6-dec 1994   David Stevenson to store empty top layer
!-                          o3 & hno3 emissions
!-   Updated   8-MAR-1995   Bill Collins
!-                           Now looks to see if there are any cells
!-                          in the boundary layer
!-   Updated  10-MAR-1995   Bill Collins
!-                           Don't store isoprene if it is night
!-   Replaced 25-JUL-1995   Colin Johnson
!-                           From Bill, updates to species list.
!-   Updated   6-AUG-1996   Bill Collins  Parameters now in INCLUDE
!-   Updated  17-OCT-1997   Bill Collins  Converted to Fortran 90
!-   Updated  10-MAR-1998   Bill Collins
!-                           Reduced size of Eulerian arrays
!-   Updated  22-OCT-1999   Colin Johnson Removed species magic Nos.
!-
!VVV  V2.3  STORE 22/X/99 - Removed species magic nos.
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: month
      INTEGER, DIMENSION(nlnpe,nlpe,nlev), INTENT(IN) :: nnn
      INTEGER, DIMENSION(nlnpe,nlpe),      INTENT(IN) :: nbl

      REAL, INTENT(IN) :: time
      REAL, DIMENSION(nc,nlnpe,nlpe),  INTENT(IN) :: emiss
      REAL, DIMENSION(nc,nlnpe,nlpe), INTENT(OUT) :: estore

      INTEGER :: i
      INTEGER :: j
      INTEGER :: k
      REAL :: thet
      REAL :: ti
      REAL :: rlong
      REAL :: rlat

!  FIND grid boxes IN boundary LAYER with no cells (not O3 or HNO3)
!+ FIND grid boxes IN top LAYER with no cells (O3 & HNO3)

      DO k=1,nc
        IF(k==i_o3.OR.k==i_hno3) THEN
! check for cells in top level for stratospheric emissions
          WHERE(nnn(:,:,nlev)==0)
            estore(k,:,:)=estore(k,:,:)+emiss(k,:,:)
          ELSE WHERE
            estore(k,:,:)=0.
          END WHERE
        ELSE IF(k/=i_c5h8) THEN
          WHERE(nbl==0)
            estore(k,:,:)=estore(k,:,:)+emiss(k,:,:)
          ELSE WHERE
            estore(k,:,:)=0.
          END WHERE
        ELSE
! Test to see if it is daytime for isoprene emissions
! DEPENDS ON: secs
          ti=SECS(15.0,month,1)+MOD(time,daysec)
          DO i=1,nlnpe
            DO j=1,nlpe
              rlat=(lat(j+ltdat-1-1)+lat(j+ltdat-1))/2.-90.
              IF(i+lndat-1<nlong) THEN
                rlong=(long(i+lndat-1)+long(i+lndat-1+1))/2.
              ELSE
                rlong=(long(i+lndat-1)+                                 &
     &                 long(MOD(i+lndat-1,nlong)+1)+360.)/2.
              ENDIF
! DEPENDS ON: zen
              thet=ZEN(ti,rlat,rlong)
              IF(nbl(i,j)==0) THEN
                estore(k,i,j)=estore(k,i,j)+emiss(k,i,j)*               &
     &            MAX(COS(thet),0.) ! i.e. set to 0. if COS(THET)<0.
              ELSE
                estore(k,i,j)=0.
              ENDIF
            END DO
          END DO
        END IF
      END DO
      END SUBROUTINE STORE
#endif
