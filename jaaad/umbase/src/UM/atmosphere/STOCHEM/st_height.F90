#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      INTEGER FUNCTION ST_HEIGHT(pos,eta_array)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Calculate height indicies on eta grids
!-
!-   Inputs  : POS,ETA_ARRAY
!-   Outputs : HEIGHT
!-   Controls:
!-
!
! Current Code Owner: Bill Collins
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.0  29/05/01   New Dynamics version.  C.E. Johnson
!  6.1  06/08/04   Renamed function as ST_HEIGHT to avoid clash with
!                  existing variable. M.G. Sanderson
!  6.2  31/01/06   Place names of interface blocks on end statements
!                  behind comments for portability. T. Edwards
!-
!VVV  V5.0  HEIGHT 29/V/01 - New Dynamics version
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!----------------------------------------------------------------------
      REAL, INTENT(IN) :: POS
      CHARACTER(*), INTENT(IN) :: ETA_ARRAY

      INTEGER :: N1
      CHARACTER(LEN=72) :: cmessage

      N1=1
      IF (ETA_ARRAY == 'Eta_theta') THEN
!        write(6,*) eta_array
        DO
          IF(POS<=Eta_theta(N1)) EXIT
          N1=N1+1
        ENDDO
        st_height=n1-1
      ELSE IF (ETA_ARRAY == 'Eta_rho') THEN
!        write(6,*) eta_array
        IF (POS > Eta_rho(NMETLEV)) THEN
          ST_HEIGHT=nmetlev
        ELSE
          DO
            IF(POS<=Eta_rho(N1)) EXIT
            N1=N1+1
          ENDDO
          ST_HEIGHT=n1-1
        ENDIF
      ELSE IF (ETA_ARRAY == 'Eta_stochem') THEN
!        write(6,*) eta_array
        DO
          IF(POS<=Eta_stochem(N1)) EXIT
          N1=N1+1
        ENDDO
        ST_HEIGHT=n1
      ELSE
        cmessage='Unknown eta array in st_height function'
        write(6,*) cmessage,eta_array
! DEPENDS ON: ereport
        CALL EREPORT('st_height',1,cmessage)
      ENDIF

      END FUNCTION ST_HEIGHT

#endif
