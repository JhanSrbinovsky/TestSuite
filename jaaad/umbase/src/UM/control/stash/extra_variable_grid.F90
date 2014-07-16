#if defined(C84_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE EXTRA_VARIABLE_GRID(buf3,LEN_IN                        &
     &            ,SROW_OUT,N_ROWS_OUT,WCOL_OUT,N_COLS_OUT)

      IMPLICIT NONE
!
! Description: This routine controls the addition of
!              variable horizontal grid information in the
!              "extra data" part of a PP field.
!
! Author :     R Hill
!
! Modification History
! UM Version   Date     Description
! ----------  --------  --------------------------------
!   5.4       Mar 2002  Original Code.
!==========================================================
#include "parvars.h"
#include "comvgrid.h"

      INTEGER :: LEN_IN              ! IN size of buf3
      INTEGER :: SROW_OUT            ! IN 1st southern row
      INTEGER :: N_ROWS_OUT          ! IN No. of rows to output
      INTEGER :: WCOL_OUT            ! IN 1st western col
      INTEGER :: N_COLS_OUT          ! IN No. of cols to output

      REAL :: BUF3(LEN_IN)           ! IN/OUT data for output

      INTEGER :: EXTRA_START         ! Local start point for
                                     ! extra data
      INTEGER :: VECTOR_SIZE         ! Local size of vector


      ! The start point for any extra data is at the first point
      ! after the incoming data.
      EXTRA_START = 1

      ! If grid is variable in the X direction, add appropriate
      ! extra data.
      IF (X_VAR_GRID) THEN

         ! Coordinates first
! DEPENDS ON: extra_vgrid_info
         CALL EXTRA_VGRID_INFO(BUF3(EXTRA_START),N_COLS_OUT             &
     &              ,X_GRID(WCOL_OUT,VAR_GRID_TYPE),VECTOR_SIZE         &
     &              ,x_coord_vector)
         EXTRA_START = EXTRA_START + VECTOR_SIZE

         ! Lower boundaries
! DEPENDS ON: extra_vgrid_info
         CALL EXTRA_VGRID_INFO(BUF3(EXTRA_START),N_COLS_OUT             &
     &              ,X_BOUNDARY(WCOL_OUT,VAR_GRID_TYPE),VECTOR_SIZE     &
     &              ,x_lbnd_vector)
         EXTRA_START = EXTRA_START + VECTOR_SIZE

         ! Upper boundaries
! DEPENDS ON: extra_vgrid_info
         CALL EXTRA_VGRID_INFO(BUF3(EXTRA_START),N_COLS_OUT             &
     &              ,X_BOUNDARY(WCOL_OUT+1,VAR_GRID_TYPE),VECTOR_SIZE   &
     &              ,x_ubnd_vector)

         EXTRA_START = EXTRA_START + VECTOR_SIZE
      ENDIF

      IF (Y_VAR_GRID) THEN
         ! Coordinates first
! DEPENDS ON: extra_vgrid_info
         CALL EXTRA_VGRID_INFO(BUF3(EXTRA_START),N_ROWS_OUT             &
     &              ,Y_GRID(SROW_OUT,VAR_GRID_TYPE),VECTOR_SIZE         &
     &              ,y_coord_vector)
         EXTRA_START = EXTRA_START + VECTOR_SIZE

         ! Lower boundaries
! DEPENDS ON: extra_vgrid_info
         CALL EXTRA_VGRID_INFO(BUF3(EXTRA_START),N_ROWS_OUT             &
     &              ,Y_BOUNDARY(SROW_OUT,VAR_GRID_TYPE),VECTOR_SIZE     &
     &              ,y_lbnd_vector)
         EXTRA_START = EXTRA_START + VECTOR_SIZE

         ! Upper boundaries
! DEPENDS ON: extra_vgrid_info
         CALL EXTRA_VGRID_INFO(BUF3(EXTRA_START),N_ROWS_OUT             &
     &              ,Y_BOUNDARY(SROW_OUT+1,VAR_GRID_TYPE),VECTOR_SIZE   &
     &              ,y_ubnd_vector)
         EXTRA_START = EXTRA_START + VECTOR_SIZE
      ENDIF

      RETURN
      END SUBROUTINE EXTRA_VARIABLE_GRID

#endif
