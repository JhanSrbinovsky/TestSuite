#if defined(A10_2A) || defined(A10_2B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine solve_forward_y1
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

      SUBROUTINE solve_forward_y1(                                      &
     &                             a0_y, a1_y, a2_y,factor_y,           &
     &                             row_length, rows, model_levels)

      IMPLICIT NONE

! Arguments:

      INTEGER                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels

      REAL                                                              &
     &  a0_y(row_length,rows,model_levels)                              &
                                                 ! IN/OUT
     &, a1_y(row_length,rows,model_levels)                              &
                                                 ! IN
     &, a2_y(row_length,rows,model_levels)                              &
                                                 ! IN
     &, factor_y(row_length,rows,model_levels)   ! OUT  ! Include files

#include "parvars.h"

! Locals:
      INTEGER                                                           &
     &  i,j,k                                                           &
     &, j_start,j_end

#if defined(T3E)
!DIR$ SYMMETRIC receive_data
      REAL                                                              &
     &  receive_data(2,glsize(1,fld_type_p),MODEL_LEVELS)
                                  ! Data from neighbouring processor
                                  ! Index 1: a0_y
                                  ! Index 2: a1_y

!DIR$ SYMMETRIC level_sync
      INTEGER                                                           &
     &  level_sync(MODEL_LEVELS)                                        &
                                  ! Non-zero indicates there is data in
                                  ! halo ready to be used
     &, ready_val                 ! Used for setting level_sync
#else
      REAL                                                              &
     &  send_data(2,ROW_LENGTH,MODEL_LEVELS),                           &
                                                 ! Data to send
     &  receive_data(2,ROW_LENGTH,MODEL_LEVELS)  ! Received data
                                                 ! Index 1: a0_y
                                                 ! Index 2: a1_y
      INTEGER                                                           &
     &  info ! GCOM return code
#endif

!------------------------------------------------------------

#if defined(T3E)
      level_sync=0

      CALL barrier()  ! Make sure everyone's ready to start

      ready_val=1
#endif

      j_start=2
      j_end=rows
      IF (at_extremity(PSouth)) j_start=3
      IF (at_extremity(PNorth)) j_end=rows-1

      DO k=1,MODEL_LEVELS

        IF (.NOT. at_extremity(PSouth)) THEN
          ! Everyone except Southernmost processor waits to receive
          ! data from their South
#if defined(T3E)
!DIR$ SUPPRESS level_sync
          DO WHILE (level_sync(k)  ==  0) ! wait for data to arrive
!DIR$ SUPPRESS level_sync
          ENDDO                           ! from neighbouring PE
#else
          CALL GC_RRECV(k,2*ROW_LENGTH,neighbour(PSouth),info,          &
     &                  receive_data(1,1,k),                            &
     &                  send_data(1,1,k))
#endif
        ENDIF

        IF (.NOT. at_extremity(PSouth)) THEN
          ! First row, use the data I've received from
          ! neighbouring processor

          DO i=1,ROW_LENGTH
            factor_y(i,1,k) = a2_y(i,1,k)*receive_data(1,i,k)
            a0_y(i,1,k)     = 1.0/(a0_y(i,1,k)-factor_y(i,1,k)*         &
     &                             receive_data(2,i,k))
          ENDDO
        ENDIF

        DO j=j_start,j_end
          DO i=1,ROW_LENGTH
            factor_y(i,j,k) = a2_y(i,j,k)*a0_y(i,j-1,k)
            a0_y(i,j,k)     = 1.0/(a0_y(i,j,k)-factor_y(i,j,k)*         &
     &                             a1_y(i,j-1,k))
          ENDDO
        ENDDO

        IF (.NOT. at_extremity(PNorth)) THEN
           ! Everybody except Northernmost processor sends
           ! data to the processor to their North

#if defined(T3E)
          CALL shmem_iput(                                              &
     &      receive_data(1,1,k),a0_y(1,ROWS,k),                         &
     &      2,1,ROW_LENGTH,neighbour(PNorth))
          CALL shmem_iput(                                              &
     &      receive_data(2,1,k),a1_y(1,ROWS,k),                         &
     &      2,1,ROW_LENGTH,neighbour(PNorth))

          ! And set the remote level_sync(k) to indicate data
          ! is ready to be used

          CALL shmem_put(                                               &
     &      level_sync(k),ready_val,1,neighbour(PNorth))
#else
          DO i=1,ROW_LENGTH
            send_data(1,i,k)=a0_y(i,ROWS,k)
            send_data(2,i,k)=a1_y(i,ROWS,k)
          ENDDO

          CALL GC_RSEND(k,2*ROW_LENGTH,neighbour(PNorth),info,          &
     &                  receive_data(1,1,k),                            &
     &                  send_data(1,1,k))
#endif
        ENDIF ! IF (.NOT. at_extremity(PNorth))

      ENDDO ! k

      RETURN
      END SUBROUTINE solve_forward_y1
#endif
