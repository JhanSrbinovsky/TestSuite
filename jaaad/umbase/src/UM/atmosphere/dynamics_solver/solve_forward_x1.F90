#if defined(A10_2A) || defined(A10_2B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine solve_forward_x1
! Purpose : Portable solver for dynamics
!
! Author : Paul Burton (based on original code by Mark Mawson)
!
! Version    Date     Modification history from model version 5.2
! -------    ----     -------------------------------------------
! 5.2        15/09/00 Deck rewritten for portability    P.Burton
! 5.4        23/04/02 Introduced pipelining optimisation  P.Burton
!  6.0  18/08/03  NEC SX-6 optimisation - change dimension to avoid
!                 bank conflict.  R Barnes & J-C Rioual.
!                 change step size for t3e code      A. Malcolm
! 6.2        10/03/06  Add section 10_2B           A. Malcolm
!
! ---------------------------------------------------------
      SUBROUTINE solve_forward_x1(                                      &
     &                             a0_x, a1_x, a2_x,                    &
     &                             factor_x, F_vector_x,                &
     &                             row_length, rows, model_levels,      &
     &                             j_start, j_end)

      IMPLICIT NONE

! Arguments:

      INTEGER                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, j_start,j_end

      REAL                                                              &
     &  a0_x(row_length+1,rows,model_levels)                            &
                                                   ! IN/OUT
     &, a1_x(row_length+1,rows,model_levels)                            &
                                                   ! IN
     &, a2_x(row_length,rows,model_levels)                              &
                                                 ! IN
     &, factor_x(row_length+1,rows,model_levels)                        &
                                                   ! OUT
     &, F_vector_x(row_length+1,rows,model_levels) ! IN/OUT

! Include files
#include "parvars.h"

! Locals:
      INTEGER                                                           &
     &  i,j,k                                                           &
     &, i_start,i_end

#if defined(T3E)
!DIR$ SYMMETRIC receive_data
      REAL                                                              &
     &  receive_data(3,glsize(2,fld_type_p),MODEL_LEVELS)
                                  ! Data from neighbouring processor
                                  ! Index 1: a0_x
                                  ! Index 2: a1_x
                                  ! Index 3: F_vector_x

!DIR$ SYMMETRIC level_sync
      INTEGER                                                           &
     &  level_sync(MODEL_LEVELS)                                        &
                                  ! Non-zero indicates there is data in
                                  ! halo ready to be used
     &, ready_val                 ! Used for setting level_sync
#else
      REAL                                                              &
     &  send_data(3,ROWS,MODEL_LEVELS),                                 &
                                           ! Data to send
     &  receive_data(3,ROWS,MODEL_LEVELS)  ! Received data
                                           ! Index 1: a0_x
                                           ! Index 2: a1_x
                                           ! Index 3: F_vector_x

      INTEGER                                                           &
     &  info ! GCOM return code
#endif

!------------------------------------------------------------

#if defined(T3E)
      level_sync=0

      CALL barrier()  ! Make sure everyone's ready to start

      ready_val=1
#endif

      i_start=2
      IF (at_extremity(PEast)) THEN
        i_end=row_length-1
      ELSE
        i_end=row_length
      ENDIF

      DO k=1,MODEL_LEVELS

        IF (.NOT. at_extremity(PWest)) THEN
          ! Everyone except Westernmost processor waits to receive
          ! data from their West
#if defined(T3E)
!DIR$ SUPPRESS level_sync
          DO WHILE (level_sync(k)  ==  0) ! wait for data to arrive
!DIR$ SUPPRESS level_sync
          ENDDO                           ! from neighbouring PE
#else
          CALL GC_RRECV(k,3*(j_end-j_start+1),neighbour(PWest),info,    &
     &                  receive_data(1,j_start,k),                      &
     &                  send_data(1,j_start,k))
#endif
        ENDIF

        DO j=j_start,j_end

          IF (.NOT. at_extremity(PWest)) THEN
            ! First point, use the data I've received from
            ! neighbouring processor
            factor_x(1,j,k)   = a2_x(1,j,k)*receive_data(1,j,k)
            a0_x(1,j,k)       = 1.0/(a0_x(1,j,k)-factor_x(1,j,k)*       &
     &                                           receive_data(2,j,k))
            F_vector_x(1,j,k) = F_vector_x(1,j,k)-                      &
     &                          factor_x(1,j,k)*receive_data(3,j,k)
          ENDIF

          DO i=i_start,i_end
            factor_x(i,j,k)   = a2_x(i,j,k) * a0_x(i-1,j,k)
            a0_x(i,j,k)       = 1.0/(a0_x(i,j,k) - factor_x(i,j,k)*     &
     &                                             a1_x(i-1,j,k))
            F_vector_x(i,j,k) = F_vector_x(i,j,k) -                     &
     &                          factor_x(i,j,k)*F_vector_x(i-1,j,k)
          ENDDO

        ENDDO

        IF (.NOT. at_extremity(PEast)) THEN
          ! Everybody except Easternmost processor sends
          ! data to the processor to their East

#if defined(T3E)
          CALL shmem_iput(                                              &
     &      receive_data(1,j_start,k),a0_x(ROW_LENGTH,j_start,k),       &
     &      3,ROW_LENGTH+1,                                             &
     &      (j_end-j_start)+1,neighbour(PEast))
          CALL shmem_iput(                                              &
     &      receive_data(2,j_start,k),a1_x(ROW_LENGTH,j_start,k),       &
     &      3,ROW_LENGTH+1,                                             &
     &      (j_end-j_start)+1,neighbour(PEast))
          CALL shmem_iput(                                              &
     &      receive_data(3,j_start,k),F_vector_x(ROW_LENGTH,j_start,k), &
     &      3,ROW_LENGTH+1,                                             &
     &      (j_end-j_start)+1,neighbour(PEast))

          ! And set the remote level_sync(k) to indicate data
          ! is ready to be used

          CALL shmem_put(                                               &
     &      level_sync(k),ready_val,1,neighbour(PEast))
#else
          DO j=j_start,j_end
            send_data(1,j,k)=a0_x(ROW_LENGTH,j,k)
            send_data(2,j,k)=a1_x(ROW_LENGTH,j,k)
            send_data(3,j,k)=F_vector_x(ROW_LENGTH,j,k)
          ENDDO

          CALL GC_RSEND(k,3*(j_end-j_start+1),neighbour(PEast),info,    &
     &                  receive_data(1,j_start,k),                      &
     &                  send_data(1,j_start,k))
#endif
        ENDIF ! IF (.NOT. at_extremity(PEast))

      ENDDO ! k

      RETURN
      END SUBROUTINE solve_forward_x1
#endif
