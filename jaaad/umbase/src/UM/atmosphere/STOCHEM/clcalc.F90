#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE CLCALC(cloud,dc,cc)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Calculate 3d cloud field on data grid
!-
!-   Inputs  : DC,CC,CB,CT
!-   Outputs : CLOUD
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.2    08/08/96  Created. W.J. Collins
!  5.1    25/04/01  Split 3D cloud cover into layers. W.J. Collins
!  6.1    06/08/04  Changed function name to ST_HEIGHT. M.G. Sanderson
!
!VVV  V5.2  CLCALC 28/8/01 ND version
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTERFACE
! DEPENDS ON: st_height
        INTEGER FUNCTION ST_HEIGHT(pos,eta_array)
          CHARACTER(*), INTENT(IN) :: eta_array
          REAL, INTENT(IN) :: pos
        END FUNCTION ST_HEIGHT
      END INTERFACE

      REAL, DIMENSION(nlonpe,nlatpe,nmetlev), INTENT(IN) :: dc
      REAL, DIMENSION(nlonpe,nlatpe,nmetlev), INTENT(IN) :: cc
      REAL, DIMENSION(nlnpe,nlpe,nlev),      INTENT(OUT) :: cloud

      INTEGER :: i
      INTEGER :: j
      INTEGER :: kmed
      INTEGER :: khigh
      INTEGER, DIMENSION(1) :: k

      REAL, DIMENSION(nlnpe,nlpe,nlev) :: dcd
      REAL, DIMENSION(nlnpe,nlpe,nlev) :: ccd

! DEPENDS ON: st_height
      kmed = ST_HEIGHT(etamed,'Eta_stochem')
! DEPENDS ON: st_height
      khigh = ST_HEIGHT(etahigh,'Eta_stochem')
      cloud = 0.0

! Convert dynamic and convective cloud cover on to stochem grid
! DC (stash 265) is on theta levels
! DEPENDS ON: met2data
      CALL MET2DATA(dcd,dc,nmetlev,nlev,.FALSE.)

! CC (stash 5212) is on rho levels
! DEPENDS ON: met2data
      CALL MET2DATA(ccd,cc,nmetlev,nlev,.TRUE.)

! Only want one non-zero cloud layer from dynamic and convective cloud i
! each band (low,medium,high). So take the levels with the maximum cloud
      DO j = 1,nlpe
        DO i = 1,nlnpe
! Calculate maximum low cloud amount
          k = MAXLOC(dcd(i,j,:kmed))
          cloud(i,j,k) = cloud(i,j,k) + dcd(i,j,k)
          k = MAXLOC(ccd(i,j,:kmed))
          cloud(i,j,k) = cloud(i,j,k) + ccd(i,j,k)
! Calculate maximum medium cloud amount
          k = MAXLOC(dcd(i,j,kmed:khigh)) + kmed-1
          cloud(i,j,k) = cloud(i,j,k) + dcd(i,j,k)
          k = MAXLOC(ccd(i,j,kmed:khigh)) + kmed-1
          cloud(i,j,k) = cloud(i,j,k) + ccd(i,j,k)
! Calculate maximum hight cloud amount
          k = MAXLOC(dcd(i,j,khigh:)) + khigh-1
          cloud(i,j,k) = cloud(i,j,k) + dcd(i,j,k)
          k = MAXLOC(ccd(i,j,khigh:)) + khigh-1
          cloud(i,j,k) = cloud(i,j,k) + ccd(i,j,k)
        END DO
      END DO
      cloud = MIN(cloud,1.0)

      END SUBROUTINE CLCALC
#endif
