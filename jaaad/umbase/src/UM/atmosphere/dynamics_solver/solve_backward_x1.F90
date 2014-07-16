#if defined(A10_2A) || defined(A10_2B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine solve_backward_x1
! Purpose : Portable solver for dynamics
!
! Author : Paul Burton (based on original code by Mark Mawson)
!
! Version    Date     Modification history from model version 5.2
! -------    ----     -------------------------------------------
! 5.2        15/09/00 Deck rewritten for portability    P.Burton
! 5.3        19/10/01 Optimisations for speed of execution.  S. Cusack
! 5.4        23/04/02 Introduced pipelining optimisation  P.Burton
!  6.0  18/08/03  NEC SX-6 optimisation - change dimension to avoid
!                 bank conflict.  R Barnes & J-C Rioual.
!                 change step size for t3e code      A. Malcolm
! 6.0        12/12/03  Vectorisation optimisation. P.Selwood.
! 6.2        10/03/06  Add section 10_2B           A. Malcolm
!
! ---------------------------------------------------------

      SUBROUTINE solve_backward_x1(                                     &
     &                              soln, a0_x, a1_x,                   &
     &                              row_length, rows, model_levels,     &
     &                              j_start, j_end, off_x, off_y)

      IMPLICIT NONE

! Arguments:

      INTEGER                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, j_start,j_end                                                   &
     &, off_x,off_y

      REAL                                                              &
     &  soln(1-off_x:row_length+off_x+1,                                &
                                                   ! IN/OUT
     &       1-off_y:rows+off_y,                                        &
     &       model_levels)                                              &
     &, a0_x(row_length+1,rows,model_levels)                            &
                                                   ! IN
     &, a1_x(row_length+1,rows,model_levels)       ! IN

! Include files
#include "parvars.h"

! Locals:
      INTEGER                                                           &
     &  i,j,k                                                           &
     &, i_start,i_end

#if defined(T3E)
!DIR$ SYMMETRIC receive_data
      REAL                                                              &
     &  receive_data(glsize(2,fld_type_p),MODEL_LEVELS)
                                  ! Data from neighbouring processor

!DIR$ SYMMETRIC level_sync
      INTEGER                                                           &
     &  level_sync(MODEL_LEVELS)                                        &
                                  ! Non-zero indicates there is data in
                                  ! halo ready to be used
     &, ready_val                 ! Used for setting level_sync
#else
      REAL                                                              &
     &  send_data(ROWS,MODEL_LEVELS),                                   &
                                         ! Data to send
     &  receive_data(ROWS,MODEL_LEVELS)  ! Received data

      INTEGER                                                           &
     &  info ! GCOM return code
#endif

!------------------------------------------------------------

#if defined(T3E)
      level_sync=0

      CALL barrier()  ! Make sure everyone's ready to start

      ready_val=1
#endif

      i_end=1
      IF (at_extremity(PEast)) THEN
        i_start=row_length-2
      ELSE
        i_start=row_length-1
      ENDIF

      DO k=1,MODEL_LEVELS

        IF (.NOT. at_extremity(PEast)) THEN
          ! Everyone except Easternmost processor waits to receive
          ! data from their East
#if defined(T3E)
!DIR$ SUPPRESS level_sync
          DO WHILE (level_sync(k)  ==  0) ! wait for data to arrive
!DIR$ SUPPRESS level_sync
          ENDDO                           ! from neighbouring PE
#else
          CALL GC_RRECV(k,j_end-j_start+1,neighbour(PEast),info,        &
     &                  receive_data(j_start,k),                        &
     &                  send_data(j_start,k))
#endif
        ENDIF

        IF (.NOT. at_extremity(PEast)) THEN
          DO j = j_start, j_end
            ! First point, use the data I've received from
            ! neighbouring processor
            i=i_start+1
            soln(i,j,k)=a0_x(i,j,k)*(soln(i,j,k)-                       &
     &                  a1_x(i,j,k)*                                    &
     &                  receive_data(j,k))
          END DO
        END IF

        DO j = j_start, j_end
          DO i=i_start,i_end,-1
            soln(i,j,k)=a0_x(i,j,k)*(soln(i,j,k)-                       &
     &                               a1_x(i,j,k)*                       &
     &                               soln(i+1,j,k))
          ENDDO
        ENDDO

        IF (.NOT. at_extremity(PWest)) THEN
           ! Everybody except Westernmost processor sends
           ! data to the processor to their West

#if defined(T3E)
          CALL shmem_iput(                                              &
     &      receive_data(j_start,k),soln(1,j_start,k),                  &
     &      1,ROW_LENGTH+2*Off_x+1,                                     &
     &      (j_end-j_start)+1,neighbour(PWest))

          ! And set the remote level_sync(k) to indicate data
          ! is ready to be used

          CALL shmem_put(                                               &
     &      level_sync(k),ready_val,1,neighbour(PWest))
#else
          DO j=j_start,j_end
            send_data(j,k)=soln(1,j,k)
          ENDDO

          CALL GC_RSEND(k,j_end-j_start+1,neighbour(PWest),info,        &
     &                  receive_data(j_start,k),                        &
     &                  send_data(j_start,k))
#endif
        ENDIF ! IF (.NOT. at_extremity(PWest))

      ENDDO ! k
      RETURN
      END SUBROUTINE solve_backward_x1
#endif
