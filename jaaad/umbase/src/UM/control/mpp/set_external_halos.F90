#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Routine for setting the External Halos to a provided
! constant value.

      SUBROUTINE SET_EXTERNAL_HALOS(                                    &
     &  FIELD,ROW_LENGTH,ROWS,LEVELS,                                   &
     &  HALO_X, HALO_Y, VALUE)

      IMPLICIT NONE

! Purpose:
!  Set the external halos (those around the edge of model domain) to
!  a provided constant value.
!  This is useful for LAM models because SWAPBOUNDS doesn't fill these
!  halos, they are normally filled with LBC data. However, not all
!  fields have LBC data applied to them, and such fields need to have
!  sensible numbers put in these external halo regions.
!  This subroutine is similar in principle to FILL_EXTERNAL_HALOS, but
!  a provided value rather than the adjacent point's value is used.
!
! Author: Nicola Robertson
! Current code owner:  UM system development team.
!
! History
! Date       Version   Comment
! --------   -------   -------
! 26/04/02   5.4       New subroutine added. Nicola Robertson
!  6.0  17/09/03  Add def for new NEC opt section C96_1C. R Barnes
!
! Arguments:

      INTEGER, INTENT(IN) ::                                            &
     &  ROW_LENGTH                                                      &
                         ! number of points on a row
                         !        (not including halos)
     &, ROWS                                                            &
                         ! number of rows in a theta field
                         !        (not including halos)
     &, LEVELS                                                          &
                         ! number of model levels
     &, HALO_X                                                          &
                         ! size of halo in "i" direction
     &, HALO_Y           ! size of halo in "j" direction

      REAL, INTENT(IN) ::                                               &
     &  VALUE            ! value to set halos to.

      REAL, INTENT(INOUT) ::                                            &
     &  FIELD(1-HALO_X:ROW_LENGTH+HALO_X,                               &
     &        1-HALO_Y:ROWS+HALO_Y,                                     &
     &        LEVELS)    ! Field to have its halos updated


! Comdecks
#include "parvars.h"

! Local variables

      INTEGER ::                                                        &
     &  i,j,k            ! loop indicies

!----------------------------------------------------------------

! Western halo region

      IF (neighbour(PWest)  ==  NoDomain) THEN
        DO k=1,LEVELS
          DO j=1-HALO_Y,ROWS+HALO_Y
            DO i=1-HALO_X,0
              FIELD(i,j,k)=VALUE
            ENDDO ! i
          ENDDO ! j
        ENDDO ! k
      ENDIF ! IF (neighbour(PWest)  ==  NoDomain)

! Eastern halo region

      IF (neighbour(PEast)  ==  NoDomain) THEN
        DO k=1,LEVELS
          DO j=1-HALO_Y,ROWS+HALO_Y
            DO i=ROW_LENGTH+1,ROW_LENGTH+HALO_X
              FIELD(i,j,k)=VALUE
            ENDDO ! i
          ENDDO ! j
        ENDDO ! k
      ENDIF ! IF (neighbour(PEast)  ==  NoDomain)

! Northern halo region

      IF (neighbour(PNorth)  ==  NoDomain) THEN
        DO k=1,LEVELS
          DO j=ROWS+1,ROWS+HALO_Y
            DO i=1-HALO_X,ROW_LENGTH+HALO_X
              FIELD(i,j,k)=VALUE
            ENDDO ! i
          ENDDO ! j
        ENDDO ! k
      ENDIF ! IF (neighbour(PNorth)  ==  NoDomain)

! Southern halo region

      IF (neighbour(PSouth)  ==  NoDomain) THEN
        DO k=1,LEVELS
          DO j=1-HALO_Y,0
            DO i=1-HALO_X,ROW_LENGTH+HALO_X
              FIELD(i,j,k)=VALUE
            ENDDO ! i
          ENDDO ! j
        ENDDO ! k
      ENDIF ! IF (neighbour(PSouth)  ==  NoDomain)

      RETURN
      END SUBROUTINE SET_EXTERNAL_HALOS

#endif
