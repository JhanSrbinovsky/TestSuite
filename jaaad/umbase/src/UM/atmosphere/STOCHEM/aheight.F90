#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE AHEIGHT(asize,pos,eta_arr_name,level_num)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Calculate height indicies on eta grids
!-   This version expects an array of values
!-
!-   Inputs  : asize,pos,eta_arr_name
!-   Outputs : level_num
!-   Controls:
!-
!
! Current Owner of Code: M.G. Sanderson
!
! History:
! Version   Date                    Comment
!  5.5    08/01/04  Created.  M.G. Sanderson
!  5.5    22/01/04  Updated.  K. Ketelsen
!  6.1    20/10/04  No Change.
!-
!VVV  V1.0  AHEIGHT 08/I/04 - New Dynamics version
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: asize
      REAL, DIMENSION(asize), INTENT(IN) :: pos
      CHARACTER(*), INTENT(IN) :: eta_arr_name
      INTEGER, DIMENSION(asize), INTENT(OUT) :: level_num

      INTEGER :: j, nlevels
      REAL, DIMENSION(0:nmetlev) :: eta_array
      CHARACTER(LEN=72) :: cmessage

! Set up nlevels and eta_array with appropriate values
      eta_array = 1.0
      IF (eta_arr_name == 'Eta_theta') THEN
        nlevels = nmetlev
        eta_array = eta_theta
      ELSE IF (eta_arr_name == 'Eta_rho') THEN
        nlevels = nmetlev
        eta_array(1:nmetlev) = eta_rho
      ELSE IF (eta_arr_name == 'Eta_stochem') THEN
        nlevels = nlev
        eta_array(0:nlev) = eta_stochem
      ELSE
        cmessage='Unknown eta array in height function'
        write(6,*) cmessage,eta_arr_name
! DEPENDS ON: ereport
        CALL EREPORT('HEIGHT',1,cmessage)
      END IF

! Calculate vertical level number. Adjust values as necessary
      level_num = 1
      DO j = 1, nlevels
        WHERE (pos > eta_array(j)) level_num = level_num + 1
      END DO
      IF (eta_arr_name == 'Eta_theta' .OR. eta_arr_name == 'Eta_rho')   &
     &  level_num = level_num - 1
      IF (eta_arr_name == 'Eta_rho') THEN
        WHERE (pos > eta_rho(nmetlev)) level_num = nmetlev
      END IF

      END SUBROUTINE AHEIGHT
#endif
