#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Deletes duplicate diags & times; checks for overlap levs & times.
!
! Subroutine Interface:

      SUBROUTINE DUPLIC(NRECS,NTIMES,NLEVELS)
      IMPLICIT NONE
!
! Description:
!   Deletes duplicate diagnostic entries from STASH list; deletes
!   duplicate STASH times;
!   checks for overlap of levels and times, to some extent.
!   Called by STPROC.
!   Input : NRECS   No. of STASH list records
!           NTIMES  No. of STASH times
!           NLEVELS No. of STASH levels
!           LIST_S  STASH list array with prelim. pointer system
!           ITIM_S  STASH times array
!   Output: NRECS   Reduced no. of STASH list records
!           NTIMES  Reduced no. of STASH times
!           NLEVEL  Reduced no. of STASH levels
!           ITIM_S  Reduced STASH times array
!           LIST_S  Reduced STASH list with prelim. pointers,
!                           consistent with STASH times.
!
! Method:
!
!   (a) STASH times tables in ITIM_S.
!       The times at which STASH processing occurs for a diagnostic
!   IREC may be specified by the entries (start_time,end_time,period)
!   in LIST_S.
!       Alternatively, if LIST_S(st_freq_code,IREC) has value '-n',
!   then STASH processing times for this diagnostic are given by a
!   'times table' in ITIM_S.
!   (In such a case, the above 3 entries in LIST_S are ignored).
!   The times table is given by column 'n' of ITIM_S, i.e.,
!   ITIM_S(time,n).
!   In this routine, the logical array entry LTUSE(n) is set to
!   .TRUE. if col 'n' of ITIM_S contains a times table. Any column
!   of ITIM_S which does not contain a times table is filled with
!   -1's. The cols which contain times tables are then shuffled along,
!   so that they occupy the first NTIMES cols of ITIM_S. The pointers
!   in LIST_S(st_freq_code,IREC) are altered accordingly.
!
!   (b) STASH levels lists in LEVLST_S.
!       The levels on which STASH processing occurs for a diagnostic
!   IREC is specified by the entries (output_bottom, output_top) in
!   LIST_S.
!     If LIST_S(bot)=m, then output is on a range of model levels,
!   with level m as the bottom level, and LIST_S(top) points to the
!   top output model level.
!     If LIST_S(bot)=-n, then there is a levels list in col 'n' of
!   LEVLST_S, and LIST_S(top) contains a code value indicating the
!   type of levels (model, pressures, heights or theta). Each levels
!   list also has a corresponding entry in LLISTTY, indicating whether
!   the list is real or integer.
!     In this routine, the cols of LEVLST_S which contain levels lists
!   are shuffled along so that they occupy the first NLEVELS cols of
!   LEVLST_S. The pointers in LIST_S(output_bottom,IREC), and the
!   entries in LLISTTY, are altered accordingly.
!
! Current code owner:  S.J.Swarbrick
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Mar. 95    Original code.  S.J.Swarbrick
!   5.2     31/01/01   Split .OR. statement into two parts because new
!                      compiler otherwise evaluates second part before
!                      the first causing potential divide by zero error
!                      S.D.Mullerworth
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Global variables:

#include "csubmodl.h"
#include "cppxref.h"
#include "version.h"
#include "cstash.h"
#include "stextend.h"
#include "stparam.h"

! Subroutine Arguments:
!
!   Scalar arguments with intent(InOut):

      INTEGER NRECS      ! No. of STASH list records
      INTEGER NTIMES     ! No. of STASH times
      INTEGER NLEVELS    ! No. of STASH levels

! Local scalars:

      LOGICAL LTRPT
      LOGICAL TESTFLAG
      INTEGER I
      INTEGER I1
      INTEGER I2
      INTEGER IEND
      INTEGER IITM
      INTEGER IL
      INTEGER IREC
      INTEGER ISEC
      INTEGER ISTR
      INTEGER IT
      INTEGER ITAGS1
      INTEGER ITAGS2
      INTEGER ITAGU1
      INTEGER ITAGU2
      INTEGER ITEND1
      INTEGER ITEND2
      INTEGER MODL
      INTEGER NLEVSW
      INTEGER NRECSW
      INTEGER NTIMESW

! Local arrays:

      LOGICAL LTUSE(2*NPROFTP+2)  ! LTUSE(n) set to .T. if column n
                                  ! in ITIM_S contains a STASH times
                                  ! table.

! Function and subroutine calls:

      EXTERNAL SINDX,ORDER

!- End of Header -----------------------------------------------------


! Initialise LTUSE array

      DO I=1,NTIMES
        LTUSE(I)=.FALSE.
      END DO


! Blank out unused STASH times

      DO IREC=1,NRECS
        IF (LIST_S(st_freq_code,IREC) <  0) THEN   ! STASH times table
          LTUSE(-LIST_S(st_freq_code,IREC))=.TRUE. !  exists for IREC
        END IF
      END DO

      DO I=1,NTIMES
        IF(.NOT.LTUSE(I)) ITIM_S(1,I)=-1  ! Fill unused columns in
      END DO                              ! ITIM_S with -1 in each row.


! Delete blank STASH times

      NTIMESW=1

      DO IT=1,NTIMES

! If col 'IT' contains a times table, find
! corresponding record IREC in LIST_S, and replace entry '-IT'
! by '-NTIMESW'. In each case, NTIMESW <= IT.

        IF(ITIM_S(1,IT) /= -1) THEN
          DO IREC=1,NRECS

            IF (LIST_S(st_freq_code,IREC) == -IT) THEN
                LIST_S(st_freq_code,IREC)=-NTIMESW
            END IF

          END DO

          IF (IT /= NTIMESW) THEN
 ! Move times table in col 'IT' to col 'NTIMESW'. Hence array
 ! ITIM_S is compressed.
            DO I=1,NTIMEP
              ITIM_S(I,NTIMESW)=ITIM_S(I,IT)
            END DO
          END IF

          NTIMESW=NTIMESW+1
        END IF
      END DO

      NTIMES=NTIMESW-1  ! No. of STASH-times tables remaining, so far


! Delete blank STASH levels

      NLEVSW=1

      DO IL=1,NLEVELS

! If col 'IL' of LEVLST_S contains a levs list, then find corresponding
!  record IREC in LIST_S and replace entry '-IL' by '-NLEVSW'. In each
! case, NLEVSW <= IL.

        IF(LEVLST_S(1,IL) /= 0) THEN
          DO IREC=1,NRECS
            IF (LIST_S(st_output_bottom,IREC) == -IL) THEN
                LIST_S(st_output_bottom,IREC)=-NLEVSW
            END IF
          END DO
          IF(IL /= NLEVSW) THEN
 ! Move levels list in col 'IL' to col 'NLEVSW'. Hence array
 ! LEVLST_S is compressed.
            DO I=1,NLEVP_S
              LEVLST_S(I,NLEVSW)=LEVLST_S(I,IL)
            END DO
            LLISTTY(NLEVSW)=LLISTTY(IL) ! Move corresponding entry in
          END IF                        ! LLISTTY

          NLEVSW=NLEVSW+1
        END IF
      END DO

      NLEVELS=NLEVSW-1


! Check for duplication/overlap of STASH levels

      NRECSW=NRECS

      DO MODL  = 1,N_INTERNAL_MODEL_MAX
      DO ISEC  = 0,NSECTP
      DO IITM  = 1,NITEMP

        IF(INDX_S(2,MODL,ISEC,IITM) >= 2) THEN  !More than one STASH rec
                                                !  for (model,sec,item)
          ISTR=     INDX_S(1,MODL,ISEC,IITM)    !1st record with m,s,i
          IEND=ISTR+INDX_S(2,MODL,ISEC,IITM)-1  !Last record with m,s,i

          DO I1=ISTR,IEND-1

           ITAGS1=LIST_S(st_macrotag,I1)/1000              ! System tag
           ITAGU1=LIST_S(st_macrotag,I1)-1000*ITAGS1       ! User tag

           IF (LIST_S(st_model_code,I1) <= N_INTERNAL_MODEL_MAX) THEN
!              Not flagged redundant
            DO I2=I1+1,IEND

              ITAGS2=LIST_S(st_macrotag,I2)/1000           ! System tag
              ITAGU2=LIST_S(st_macrotag,I2)-1000*ITAGS2    ! User tag

              IF((LIST_S(st_proc_no_code,I1) ==                         &
     &            LIST_S(st_proc_no_code,I2)).AND.                      &
     &           (LIST_S(st_freq_code,I1) ==                            &
     &            LIST_S(st_freq_code,I2)).AND.                         &
     &           (LIST_S(st_period_code,I1) ==                          &
     &            LIST_S(st_period_code,I2)).AND.                       &
     &           (LIST_S(st_gridpoint_code,I1) ==                       &
     &            LIST_S(st_gridpoint_code,I2)).AND.                    &
     &           (LIST_S(st_weight_code,I1) ==                          &
     &            LIST_S(st_weight_code,I2)).AND.                       &
     &           (LIST_S(st_north_code,I1) ==                           &
     &            LIST_S(st_north_code,I2)).AND.                        &
     &           (LIST_S(st_south_code,I1) ==                           &
     &            LIST_S(st_south_code,I2)).AND.                        &
     &           (LIST_S(st_west_code,I1) ==                            &
     &            LIST_S(st_west_code,I2)).AND.                         &
     &           (LIST_S(st_east_code,I1) ==                            &
     &            LIST_S(st_east_code,I2)).AND.                         &
     &           (LIST_S(st_input_code,I1) ==                           &
     &            LIST_S(st_input_code,I2)).AND.                        &
     &           (LIST_S(st_output_code,I1) ==                          &
     &            LIST_S(st_output_code,I2)).AND.                       &
     &           (LIST_S(st_series_ptr,I1) ==                           &
     &            LIST_S(st_series_ptr,I2)).AND.                        &
     &           (LIST_S(st_pseudo_out,I1) ==                           &
     &            LIST_S(st_pseudo_out,I2)).AND.                        &
     &        ((ITAGS1 == ITAGS2).OR.(ITAGS1 == 0).OR.                  &
     &        (ITAGS2 == 0)).AND.                                       &
     &        ((ITAGU1 == ITAGU2).OR.(ITAGU1 == 0).OR.                  &
     &        (ITAGU2 == 0)).AND.                                       &
     &      (LIST_S(st_model_code,I2) <= N_INTERNAL_MODEL_MAX)) THEN
!            Not flagged redundant

! If they are the same in all but time and level

                ITEND1=LIST_S(st_end_time_code,I1)
                ITEND2=LIST_S(st_end_time_code,I2)

                IF(ITEND1 == -1) ITEND1=                                &
     &             LIST_S(st_start_time_code,I2)+1     ! Force overlap

                IF(ITEND2 == -1) ITEND2=                                &
     &             LIST_S(st_start_time_code,I1)+1     ! Force overlap

! Where period_code is zero we have to prevent second part from
! being evaluated so break this OR statement into two parts:
                TESTFLAG=.FALSE.
                IF((LIST_S(st_period_code,I1) == 0).OR.                 &
     &            (LIST_S(st_period_code,I1) == -1)) THEN
                  TESTFLAG=.TRUE.
                ELSEIF((MOD(LIST_S(st_start_time_code,I2)-              &
     &               LIST_S(st_start_time_code,I1),                     &
     &               LIST_S(st_period_code,I1)) == 0)) THEN
                  TESTFLAG=.TRUE.
                ENDIF

                IF((.NOT.((LIST_S(st_start_time_code,I1)                &
     &                                            >  ITEND2).OR.        &
     &            (ITEND1 <  LIST_S(st_start_time_code,I2))).OR.        &
     &              (LIST_S(st_output_code,I1) >  0)).AND.              &
     &          (MOD(LIST_S(st_start_time_code,I2)-                     &
     &               LIST_S(st_start_time_code,I1),                     &
     &               LIST_S(st_freq_code,I1)) == 0).AND.                &
     &               TESTFLAG.AND.                                      &
     &              (LIST_S(st_output_bottom,I2) ==                     &
     &               LIST_S(st_output_bottom,I1)).AND.                  &
     &              (LIST_S(st_output_top,I2) ==                        &
     &               LIST_S(st_output_top,I1))) THEN

! (Times overlap or in dump) and overlay in freq & period
!                            and levels the same

                  IF(ITAGU1 == 0) THEN
                    LIST_S(st_macrotag,I1)=ITAGU2
                    ITAGU1=ITAGU2
                  ELSE
                    LIST_S(st_macrotag,I1)=ITAGU1
                  END IF

                  IF(ITAGS1 == 0) THEN
                    LIST_S(st_macrotag,I1)=                             &
     &              ITAGS2*1000+LIST_S(st_macrotag,I1)
                    ITAGS1=ITAGS2
                  ELSE
                    LIST_S(st_macrotag,I1)=                             &
     &              ITAGS1*1000+LIST_S(st_macrotag,I1)
                  END IF

                  LIST_S(st_start_time_code,I1)=                        &
     &            MIN(LIST_S(st_start_time_code,I1),                    &
     &                LIST_S(st_start_time_code,I2))

                  IF((LIST_S(st_end_time_code,I1) == -1).OR.            &
     &               (LIST_S(st_end_time_code,I2) == -1)) THEN
                      LIST_S(st_end_time_code,I1)=-1
                  ELSE
                    LIST_S(st_end_time_code,I1)=                        &
     &              MAX(LIST_S(st_end_time_code,I1),                    &
     &                  LIST_S(st_end_time_code,I2))
                  END IF

                  LIST_S(st_model_code,I2)=N_INTERNAL_MODEL_MAX+1
! Sets model id to be greater than no of models,
! so that this diag is put at the end of any sorted list.

                  NRECSW=NRECSW-1

                  DO I=ISTR,IEND                      !Change pointers
                    IF(LIST_S(st_input_code,I )   ==                    &
     &                -LIST_S(NELEMP+1     ,I2)) THEN
                       LIST_S(st_input_code,I ) =                       &
     &                -LIST_S(NELEMP+1     ,I1)
                    END IF
                  END DO
                END IF

              END IF   ! I1,I2 comparison
            END DO     ! I2
           END IF      ! I1 Not flagged redundant
          END DO       ! I1

        END IF   ! More than one STASH record for m,s,i
      END DO     ! Items
      END DO     ! Sections
      END DO     ! Models

! Remove unwanted records (i.e., those flagged redundant)

! DEPENDS ON: order
      CALL ORDER(NRECS)
      NRECS=NRECSW
! DEPENDS ON: sindx
      CALL SINDX(NRECS)
!
      RETURN
      END SUBROUTINE DUPLIC
#endif
