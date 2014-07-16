#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE CPUSTATS(NAME,L)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : To time subroutines
!-
!-   Inputs  : NAME:  the name of the subroutine (12 chars.)
!-             L:     the stage to be implemented
!-
!-   Outputs : Prints out cpu and wallclock times for named routines.
!-   Controls: Initialises automatically
!-             L=1 Start timing
!-             L=2 End timing
!-             L=3 Print results
!-             Add 100 to these actions for synchronization of all PEs.
!-             Set DEL_MPP_OUTPUT in SCRIPT file to false for all pe o/p
!-
!
! Current Owner of Code: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  4.5    01/05/98  Created.  C.E. Johnson
!  5.0    02/11/00  Added wallclock times & simplified.  C.E. Johnson
!  5.5    09/11/03  Changed timing subroutines to F90 intrinsic and
!                   NEC routine. M.G. Sanderson
!  6.1    21/10/04  No change.
!-
!VVV  V2.6 CPUSTATS 2/XI/00  wallclock times added.
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER,          INTENT(IN)              :: l
      CHARACTER(LEN=12), INTENT(IN)             :: name

      INTEGER                                   :: i
      INTEGER                                   :: isub
      INTEGER, PARAMETER            :: nsubs=38
      INTEGER                                   :: k
      INTEGER                       :: m
      INTEGER                                   :: info
      INTEGER                       :: rtime
      INTEGER                       :: rtime_rate
      INTEGER,DIMENSION(nsubs),SAVE :: count

      REAL,   DIMENSION(nsubs),SAVE :: cpu
      REAL,   DIMENSION(nsubs),SAVE :: before
      REAL,   DIMENSION(nsubs),SAVE :: wall
      REAL,   DIMENSION(nsubs),SAVE :: wall_b
      REAL,   DIMENSION(nsubs)      :: after
      REAL,   DIMENSION(nsubs)      :: wall_a
      REAL                                      :: sum
      REAL                          :: dd1
      REAL                          :: dd2
      CHARACTER(LEN=12), DIMENSION(nsubs) :: subnames =                 &
     &(/'MAIN        ','INIT_STOCHEM','PROCESS_MET ',                   &
     &               'STEP        ','RKO4_SINGLE ','RKO4        ',      &
     &               'CLMIX2      ','SWAP        ','JVALS       ',      &
     &               'STRATCALC2  ','SWAP2       ','DERIV       ',      &
     &               'PLOT        ','POOL        ','SURFCALC    ',      &
     &               'DRYDEP      ','ADVFLUX     ','GETNNN      ',      &
     &               'MISS        ','STRATMASK   ','CLCALC      ',      &
     &               'TEMP        ','WATER       ','MOLDENSE    ',      &
     &               'EGET        ','DDCALC      ','DWCALC      ',      &
     &               'AINDEX      ','STORE       ','STORED      ',      &
     &               'STATS       ','MCALC       ','READO3      ',      &
     &               'O3COLUMN    ','LIGHTREAD   ','DUMP        ',      &
     &  'MONTHLY_OUTP','ZERO_STATS  '/)
      LOGICAL                       :: first=.true.
! Logical to turn on processor synchronisation for all subroutines
      LOGICAL                       :: lallsync=.false.

      IF(first) THEN               ! initialise
        cpu=0.0
        wall=0.0
        count=0
        first=.FALSE.
      END IF

      DO i=1,nsubs
        IF (subnames(i)==name) THEN
          isub=i
          EXIT
        END IF
        IF (i == nsubs) THEN
          WRITE(6,*) 'CPUSTATS: SUBROUTINE ',name,' NOT FOUND'
          RETURN
        END IF
      END DO

! Only synchronize if L > 99 or LAllSync is true
      IF (l > 99) THEN
        CALL GC_GSYNC(nproc,info)
        k = l - 100
      ELSE IF (lallsync) THEN
        CALL GC_GSYNC(nproc,info)
        k = l
      ELSE
        k = l
      END IF

      IF (k == 0) THEN            ! initialise/reset
        cpu = 0.0
        wall = 0.0
        count = 0
      ELSE IF (k==1) THEN         ! start of routine
        CALL CLOCK(dd1)           ! Get current CPU clock time (s)
        before(isub) = dd1
        CALL SYSTEM_CLOCK(COUNT=rtime) ! Real time (cycles)
        wall_b(isub) = REAL(rtime)
        count(isub) = count(isub) + 1
      ELSE IF (k==2) THEN         ! end of routine
        CALL CLOCK(dd1)           ! Get CPU clock time (s)
        after(isub) = dd1
        CALL SYSTEM_CLOCK(COUNT=rtime) ! Real time (cycles)
        wall_a(isub) = REAL(rtime)
        cpu(isub)=cpu(isub)+(after(isub)-before(isub))
        wall(isub)=wall(isub)+wall_a(isub)-wall_b(isub)
      ELSE IF (k==3) THEN
! Print out the results for all processers
! Need to retain all fort6 o/p to read this, do this by setting
! UM_DEL_MPP_OUTPUT in SCRIPT file to false.

        CALL SYSTEM_CLOCK(COUNT_RATE=rtime_rate) ! Get cycles per second
        wall = wall / REAL(rtime_rate)  ! Convert to s
        WRITE(6,*) ' '
        WRITE(6,*) 'CPU TIMINGS OF NAMED ROUTINES FROM CPUSTATS:'
        WRITE(6,*) 'PE:',mype
        WRITE(6,*) 'SUBROUTINE     CPU              %     WALL',        &
     &             '     No. of CALLS'
        WRITE(6,*) '------------------------------------------',        &
     &             '-----------------'
        sum = 0.0
        DO i=4, nsubs
          sum = sum + cpu(i)
        END DO
        DO i=1,nsubs
          IF (count(1)>0) THEN
              WRITE(6,'(A12,1PE12.3,0PF12.2,1PE12.3,I12)')              &
     &        subnames(i),cpu(i),(cpu(i)/cpu(1))*100.0,                 &
     &        wall(i),count(i)
            ELSE
              WRITE(6,'(A12,1PE12.3,A12,1PE12.3,I12)')                  &
     &        subnames(i),cpu(i),'   ******   ',wall(i),count(i)
          END IF
        END DO
        IF (cpu(1) > 0.0) THEN
          sum = cpu(1)-sum-cpu(2)-cpu(3)
            WRITE(6,'(A12,1PE12.3,0PF12.2)')                            &
     &       'RESIDUAL: ',sum,(sum/cpu(1))*100.0
        END IF
          WRITE(6,*) '++++++++++++++++++++++++++++++++++++++++++',      &
     &               '+++++++++++++++++'
          WRITE(6,*) ' '
      ELSE
          WRITE(6,*) ' *** ERROR IN CPUSTATS: WRONG VALUE OF K: ',K
      END IF

      RETURN
      END SUBROUTINE CPUSTATS
#endif
