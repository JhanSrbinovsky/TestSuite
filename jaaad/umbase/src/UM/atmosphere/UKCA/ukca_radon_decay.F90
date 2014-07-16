#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  To do radioactive decay of Rn-222 to Pb-210
!  Part of the UKCA model. UKCA is a community model supported
!  by The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Current Code Owner:       Colin Johnson/Olaf Morgenstern
!                           Fiona O'Connor
!
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
! Subroutine Interface:
      SUBROUTINE RADON_DECAY(                                          &
       off_x, row_length, off_y, rows, tr_model_levels,                &
       RadonBefore, RadonAfter, TimeStep)


       INTEGER, INTENT(IN) :: off_x
       INTEGER, INTENT(IN) :: row_length
       INTEGER, INTENT(IN) :: off_y
       INTEGER, INTENT(IN) :: rows
       INTEGER, INTENT(IN) :: tr_model_levels

       REAL, INTENT(IN)    :: TimeStep    

       REAL, INTENT(INOUT) :: RadonBefore(1-off_x:row_length+off_x,    &
                                          1-off_y:rows+off_y,          &
                                          tr_model_levels) !mmr radon
       REAL, INTENT(INOUT) :: RadonAfter  (1-off_x:row_length+off_x,   &
                                           1-off_y:rows+off_y,         &
                                           tr_model_levels)   !mmr lead
!      LOCAL VARIABLES.

       INTEGER :: i
       INTEGER :: j
       INTEGER :: k

       REAL, PARAMETER :: rate=2.1E-6
       REAL            :: max_before
       REAL            :: min_before
       REAL            :: max_after
       REAL            :: min_after
       REAL            :: Delta           ! Amount of radon decayed.

       DO k = 1, tr_model_levels
        max_before = RadonBefore(1,1,k)
        min_before = RadonBefore(1,1,k)
        max_after  = RadonAfter(1,1,k)
        min_after  = RadonAfter(1,1,k)
        DO j = 1, rows
          DO i = 1, row_length
            IF (max_before <  RadonBefore(i,j,k)) THEN
                max_before = RadonBefore(i,j,k)
            END IF
            IF (min_before >  RadonBefore(i,j,k)) THEN        
                min_before = RadonBefore(i,j,k)
            END IF
            IF (max_after <  RadonAfter(i,j,k)) THEN
                max_after = RadonAfter(i,j,k)
            END IF
            IF (min_after >  RadonAfter(i,j,k)) THEN
                min_after = RadonAfter(i,j,k)
            END IF
          ENDDO
        ENDDO
       ENDDO          ! end of level loop

! This loop cycles through all points on all levels,
! not avoiding the polar points.

       DO k = 1, tr_model_levels
        DO j = 1, rows
          DO i = 1, row_length
            Delta         = (rate * TimeStep) * RadonBefore(i,j,k)
            RadonBefore(i,j,k) = RadonBefore(i,j,k) - Delta
            RadonAfter(i,j,k)  = RadonAfter(i,j,k) + Delta
          ENDDO
        ENDDO
       ENDDO

       DO k = 1, tr_model_levels
        max_before = RadonBefore(1,1,k)
        min_before = RadonBefore(1,1,k)
        max_after  = RadonAfter(1,1,k)
        min_after  = RadonAfter(1,1,k)
        DO j = 1, rows
          DO i = 1, row_length
            IF (max_before <  RadonBefore(i,j,k)) THEN
                max_before = RadonBefore(i,j,k)
            END IF
            IF (min_before >  RadonBefore(i,j,k)) THEN
                min_before = RadonBefore(i,j,k)
            END IF
            IF (max_after <  RadonAfter(i,j,k)) THEN
                max_after = RadonAfter(i,j,k)
            END IF
            IF (min_after >  RadonAfter(i,j,k)) THEN
                min_after = RadonAfter(i,j,k)
            END IF
          ENDDO
        ENDDO
        WRITE(6,*)'After chem, Level, max/min before, max/min after=', &
          k,max_before,min_before,max_after,min_after
       ENDDO                ! end of level loop

       RETURN
       END SUBROUTINE RADON_DECAY
#endif
