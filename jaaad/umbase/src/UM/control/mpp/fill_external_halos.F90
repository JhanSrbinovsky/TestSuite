#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(UTILIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Routine for filling in External Halos

      SUBROUTINE FILL_EXTERNAL_HALOS(                                   &

     &  FIELD,ROW_LENGTH,ROWS,LEVELS,                                   &
     &  HALO_X, HALO_Y)

      IMPLICIT NONE

! Purpose:
!  Fills external halos (those around the edge of model domain) with
!  sensible (copy of interior points) numbers.
!  This is useful for LAM models when SWAPBOUNDS doesn't fill these
!  halos as they are normally filled with LBC data. However, not all
!  fields have LBC data applied to them, and such fields need to have
!  sensible numbers put in these external halo regions.

! Author: Paul Burton
! Current code owner: Paul Burton
!
! History
! Date       Version   Comment
! --------   -------   -------
! 07/06/00   5.1       New deck added. Paul Burton
! 06/09/02   5.4       SX now uses this deck - add def UTILIO.
!                                                     E.Leung
!  6.0  17/09/03  Add def for new NEC opt section C96_1C. R Barnes
!
! Arguments:

      INTEGER                                                           &
     &  ROW_LENGTH                                                      &
                         ! IN: number of points on a row
                         !     (not including halos)
     &, ROWS                                                            &
                         ! IN: number of rows in a theta field
                         !     (not including halos)
     &, LEVELS                                                          &
                         ! IN: number of model levels
     &, HALO_X                                                          &
                         ! IN: size of halo in "i" direction
     &, HALO_Y           ! IN: size of halo in "j" direction

      REAL                                                              &
     &  FIELD(1-HALO_X:ROW_LENGTH+HALO_X,                               &
     &        1-HALO_Y:ROWS+HALO_Y,                                     &
     &        LEVELS)    ! IN/OUT : Field to have its halos updated


! Comdecks
#include "parvars.h"

! Local variables

      INTEGER                                                           &
     &  i,j,k            ! loop indicies

!----------------------------------------------------------------

! Western halo region

      IF (neighbour(PWest)  ==  NoDomain) THEN

        DO k=1,LEVELS
          DO j=1-HALO_Y,ROWS+HALO_Y
            DO i=1-HALO_X,0
              FIELD(i,j,k)=FIELD(i+HALO_X,j,k)
            ENDDO ! i
          ENDDO ! j
        ENDDO ! k

      ENDIF ! IF (neighbour(PWest)  ==  NoDomain)

! Eastern halo region

      IF (neighbour(PEast)  ==  NoDomain) THEN

        DO k=1,LEVELS
          DO j=1-HALO_Y,ROWS+HALO_Y
            DO i=ROW_LENGTH+1,ROW_LENGTH+HALO_X
              FIELD(i,j,k)=FIELD(i-HALO_X,j,k)
            ENDDO ! i
          ENDDO ! j
        ENDDO ! k

      ENDIF ! IF (neighbour(PEast)  ==  NoDomain)

! Northern halo region

      IF (neighbour(PNorth)  ==  NoDomain) THEN

        DO k=1,LEVELS
          DO j=ROWS+1,ROWS+HALO_Y
            DO i=1-HALO_X,ROW_LENGTH+HALO_X
              FIELD(i,j,k)=FIELD(i,j-HALO_Y,k)
            ENDDO ! i
          ENDDO ! j
        ENDDO ! k

      ENDIF ! IF (neighbour(PNorth)  ==  NoDomain)

! Southern halo region

      IF (neighbour(PSouth)  ==  NoDomain) THEN

        DO k=1,LEVELS
          DO j=1-HALO_Y,0
            DO i=1-HALO_X,ROW_LENGTH+HALO_X
              FIELD(i,j,k)=FIELD(i,j+HALO_Y,k)
            ENDDO ! i
          ENDDO ! j
        ENDDO ! k

      ENDIF ! IF (neighbour(PSouth)  ==  NoDomain)


      RETURN
      END SUBROUTINE FILL_EXTERNAL_HALOS
#endif
