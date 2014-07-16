#if defined(A10_2A) || defined(A10_2B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine solve_backward_y1
! Purpose : Portable solver for dynamics
!
! Author : Paul Burton (based on original code by Mark Mawson)
!
! Version    Date     Modification history from model version 5.2
! -------    ----     -------------------------------------------
! 5.2        15/09/00 Deck rewritten for portability    P.Burton
! 5.4        23/04/02 Introduced pipelining optimisation  P.Burton
! 6.2        10/03/06  Add section 10_2B           A. Malcolm
!
! ---------------------------------------------------------
      SUBROUTINE solve_backward_y1(                                     &
     &                              soln, a0_y, a1_y,                   &
     &                              row_length, rows, model_levels,     &
     &                              off_x,off_y)

      IMPLICIT NONE

! Arguments:

      INTEGER                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, off_x,off_y

      REAL                                                              &
     &  soln(1-off_x:row_length+off_x,                                  &
                                                 ! IN/OUT
     &       1-off_y:rows+off_y,                                        &
     &       model_levels)                                              &
     &, a0_y(row_length,rows,model_levels)                              &
                                                 ! IN
     &, a1_y(row_length,rows,model_levels)       ! IN

#include "parvars.h"

! Locals:
      INTEGER                                                           &
     &  i,j,k                                                           &
     &, j_start,j_end

#if defined(T3E)
!DIR$ SYMMETRIC receive_data
      REAL                                                              &
     &  receive_data(glsize(1,fld_type_p),MODEL_LEVELS)
                                  ! Data from neighbouring processor

!DIR$ SYMMETRIC level_sync
      INTEGER                                                           &
     &  level_sync(MODEL_LEVELS)                                        &
                                  ! Non-zero indicates there is data in
                                  ! halo ready to be used
     &, ready_val                 ! Used for setting level_sync
#else
      REAL                                                              &
     &  send_data(ROW_LENGTH,MODEL_LEVELS),                             &
                                               ! Data to send
     &  receive_data(ROW_LENGTH,MODEL_LEVELS)  ! Received data

      INTEGER                                                           &
     &  info ! GCOM return code
#endif

!------------------------------------------------------------
#if defined(T3E)
      level_sync=0

      CALL barrier()  ! Make sure everyone's ready to start

      ready_val=1
#endif


      j_end=1
      j_start=rows-1
      IF (at_extremity(PSouth)) j_end=2
      IF (at_extremity(PNorth)) j_start=rows-2

      DO k=1,MODEL_LEVELS

        IF (.NOT. at_extremity(PNorth)) THEN
          ! Everyone except Northernmost processor waits to receive
          ! data from their North
#if defined(T3E)
!DIR$ SUPPRESS level_sync
          DO WHILE (level_sync(k)  ==  0) ! wait for data to arrive
!DIR$ SUPPRESS level_sync
          ENDDO                           ! from neighbouring PE
#else
          CALL GC_RRECV(k,ROW_LENGTH,neighbour(PNorth),info,            &
     &                  receive_data(1,k),                              &
     &                  send_data(1,k))
#endif
        ENDIF

        IF (.NOT. at_extremity(PNorth)) THEN
          ! First row, use the data I've received from
          ! neighbouring processor

          j=ROWS
          DO i=1,ROW_LENGTH
            soln(i,j,k) = a0_y(i,j,k)*(soln(i,j,k)-                     &
     &                                 a1_y(i,j,k)*                     &
     &                                 receive_data(i,k))
          ENDDO
        ENDIF

        DO j=j_start,j_end,-1
          DO i=1,ROW_LENGTH
            soln(i,j,k)=a0_y(i,j,k)*(soln(i,j,k)-                       &
     &                               a1_y(i,j,k)*soln(i,j+1,k))
          ENDDO
        ENDDO

        IF (.NOT. at_extremity(PSouth)) THEN
           ! Everybody except Southernmost processor sends
           ! data to the processor to their South

#if defined(T3E)
          CALL shmem_put(                                               &
     &      receive_data(1,k),soln(1,1,k),                              &
     &      ROW_LENGTH,neighbour(PSouth))

          ! And set the remote level_sync(k) to indicate data
          ! is ready to be used

          CALL shmem_put(                                               &
     &      level_sync(k),ready_val,1,neighbour(PSouth))
#else
          DO i=1,ROW_LENGTH
            send_data(i,k)=soln(i,1,k)
          ENDDO

          CALL GC_RSEND(k,ROW_LENGTH,neighbour(PSouth),info,            &
     &                  receive_data(1,k),                              &
     &                  send_data(1,k))
#endif
        ENDIF ! IF (.NOT. at_extremity(PSouth))

      ENDDO ! k

      RETURN
      END SUBROUTINE solve_backward_y1
#endif
