#if defined(C97_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      SUBROUTINE TIMER(SUB,I_ARG)
!FPP$ NOCONCUR R
!   ..................................................................
!   SUBROUTINE TIMER
!   ----------------
!   A PROGRAM TO RECORD THE FLOW THROUGH AND TIME SPENT IN EACH
!   SUBROUTINE OF A MULTI SUBROUTINED PROGRAM.
!   CALLS TO TIMER MUST BE MADE BEFORE THE FIRST EXECUTABLE STATEMENT
!   OF EACH SUBROUTINE AND BEFORE EACH RETURN POINT.
!
!   PARAMETERS:
!   SUB - 8 BYTE CHARACTER STRING CONTAINING NAME OF CALLING
!         PROGRAM WHICH HAS BEEN LEFT ADJUSTED
!
!   I  -  I=1 - FIRST CALL FROM MAIN (OR THE HIGHEST LEVEL PROGRAM)
!         I=2 - LAST CALL FROM MAIN ( OR THE HIGHEST LEVEL PROGRAM)
!         I=3 - FIRST CALL FROM LOWER SUBROUTINE
!         I=4 CALL BEFORE RETURN STATEMENTS IN LOWER SUBROUTINE
!
!   ---NOTE :THE CRAY FACILITY PERFTRACE IS MORE APPROPRIATE FOR SINGLE
!   TASKED TIMING DIAGNOSTICS, BUT A BREAKDOWN OF ELAPSE TIME BY
!   SUBROUTINE OF MULTITASKED JOBS IS UNAVAILABLE WITHOUT THIS HOMEMADE
!   TIMER ROUTINE
!
!   REVISED TO CONFORM WITH UM SOFTWARE STANDARDS, THIS VERSION RUNS
!   ON BOTH THE CRAY & HP WORKSTATIONS. IF THERE ARE > 200 DIFFERENT
!   SUBROUTINE CALLS TO TIMER, FINAL TABLE REPLACED BY WARNING MESSAGE.
!   ..................................................................
!   ---AUTHOR OF THIS REVISION IAN EDMOND
!   ---DATE 12/04/94
!
!    Model            Modification history from model version 3.0:
!   version  Date
!LL   4.1       12/03/96  Changes to avoid clashes with TIMER2A
!LL                       P.Burton
!LL   4.5       17/04/98  Modified to cope with action (I) > 100
!LL                       intended for MPP timer.       P.Burton
!  4.5  12/06/98  Use CLOCK(...,0,2) for CPUtime on Fujitsu.
!                                        RBarnes@ecmwf.int
!  5.3  08/06/01  Replace calls to timef and second.  A van der Wal
!     5.3       26/07/01  Take an internal copy of "action" variable
!                         so as not to modify it (intent in). P.Selwood
!     6.2       15/08/05  Remove RECON def. P.Selwood
!     6.4       09/01/07  Remove MAKEBC def. P.Selwood
!
!   Programming standard :
!
!   Logical components covered : D66
!
!   Project task :
!
!   External documentation:
!
!  ---------------------------------------------------------------------
!

       IMPLICIT NONE

       INTEGER KLIMIT
       PARAMETER(KLIMIT=200)
       CHARACTER*8 SUB,SUBNAME(KLIMIT),RETURNAM(KLIMIT),SWORK
       INTEGER I_ARG    ! Argument to routine
       INTEGER NENTER(KLIMIT),I,L,J,K,NSUBS,ISTOP,NCALLS,IWORK
       REAL ELAPSE(KLIMIT),TOTE,ELPEND,ELPSTART,CPUSTART,TOT,AT,CPUEND
       REAL TOTELAP,AVELAP,PCENT,RWORK,TOTLAP,SPEEDUP,P,TIME(KLIMIT)
       REAL UWORK
       SAVE SUBNAME,RETURNAM,NENTER,ELAPSE,TIME
       SAVE K,J,NSUBS,ELPSTART,ISTOP,CPUSTART

      REAL, EXTERNAL :: Get_CPU_Time
      REAL, EXTERNAL :: Get_Wallclock_Time
!      -----------------------------------------------------------------
       I = I_ARG   ! allow changes to input argument

      IF (I  >   100) I=I-100
       IF (I == 1) THEN

!      First call to timer from the main program
!      Set up initial values of variables

         K         = 1
         J         = 0
         ISTOP     = 0
         NSUBS     = 1
         SUBNAME(1)= SUB

         DO L=1,KLIMIT
           ELAPSE(L) = 0.0
           TIME(L) = 0.0
           NENTER(L) = 0

         ENDDO

      NENTER(1)=1
! DEPENDS ON: get_wallclock_time
      ELPSTART = Get_Wallclock_Time()
! DEPENDS ON: get_cpu_time
      CPUSTART = Get_CPU_Time()

!      -----------------------------------------------------------------
       ELSEIF ((I == 2).AND.(ISTOP == 0)) THEN

!      Last call to timer from main program
!      Print out table of results

! DEPENDS ON: get_cpu_time
      CPUEND = Get_CPU_Time()
! DEPENDS ON: get_wallclock_time
      ELPEND = Get_Wallclock_Time()
         ELAPSE(1)=ELAPSE(1)+(ELPEND-ELPSTART)*1.E-3
         TIME(1)=TIME(1)+CPUEND-CPUSTART

!        Calculate total time in program
         TOTE=0.0
         TOT=0.0
         DO K=1,NSUBS

           TOTE=TOTE+ELAPSE(K)
           TOT=TOT+TIME(K)

         ENDDO

!        Sort subroutines into time order

         DO K=1,(NSUBS-1)

           DO J=(K+1),NSUBS

             IF (TIME(J) >  TIME(K)) THEN

!              Swap the values:
               RWORK=TIME(K)
               TIME(K)=TIME(J)
               TIME(J)=RWORK
               UWORK=ELAPSE(K)
               ELAPSE(K)=ELAPSE(J)
               ELAPSE(J)=UWORK
               IWORK=NENTER(K)
               NENTER(K)=NENTER(J)
               NENTER(J)=IWORK
               SWORK=SUBNAME(K)
               SUBNAME(K)=SUBNAME(J)
               SUBNAME(J)=SWORK

             ENDIF

           ENDDO

         ENDDO


!        Output timing information


      WRITE(*,'(''1'',//,20X,'' FLOW TRACE SUMMARY'',/)')
      WRITE(*,'(4X,''ROUTINE'',6X,''CPU'',6X,''%'',3X,                  &
     &  ''CALLS'',2X,''AVERAGE'',4X,''ELAPSE'',4X,''%''                 &
     &  ,6X,''AVERAGE'',1X,''CPU'')')
      WRITE(*,'(17X,''TIME'',4X,''CPU'',9X,''CPUTIME'',4X,              &
     &  ''TIME'',3X,''ELAPSE'',4X,''ELAPSE'',2X,''SPEEDUP'')')

         DO K=1,NSUBS

           SUB=SUBNAME(K)
           TOTLAP=TIME(K)
           TOTELAP=ELAPSE(K)
           NCALLS=NENTER(K)
           AVELAP=TOTELAP/NCALLS
           P=100.0*TOTLAP/TOT
           PCENT=100.0*TOTELAP/TOTE
           AT=TIME(K)/NENTER(K)
           IF (AVELAP == 0.0) THEN
             SPEEDUP=0.0
           ELSE
             SPEEDUP=AT/AVELAP
           ENDIF
           WRITE(*,'(/,T1,I3,T5,A8,T13,F10.4,T25,F5.2,T30,              &
     &       I5,T35,F10.4,T45,F10.4,T57,F5.2,T62,F10.4,T74,F5.2)')      &
     &       K,SUB,TOTLAP,P,NENTER(K),AT,TOTELAP,PCENT,AVELAP,SPEEDUP

         ENDDO

         SPEEDUP=TOT/TOTE
         WRITE(*,'(/,T3,''**TOTAL'',T12,F11.4,T44,F11.4,                &
     &     T74,F5.2)')TOT,TOTE,SPEEDUP

!      -----------------------------------------------------------------
       ELSEIF ((I == 3).AND.(ISTOP == 0)) THEN

!      First call in subroutine

!        Switch off timer
! DEPENDS ON: get_cpu_time
         CPUEND = Get_CPU_Time()
! DEPENDS ON: get_wallclock_time
         ELPEND = Get_Wallclock_Time()
      ELAPSE(K)=ELAPSE(K)+(ELPEND-ELPSTART)*1.E-3
         TIME(K)=TIME(K)+CPUEND-CPUSTART

!        Save name of calling subroutine
      J=J+1
      RETURNAM(J)=SUBNAME(K)

!        Check subroutine name
         DO K=1,NSUBS

           IF (SUBNAME(K) == SUB) GOTO 10

         ENDDO

!        New subroutine entered
      NSUBS=NSUBS+1
         IF (NSUBS  <=  KLIMIT) THEN

      SUBNAME(NSUBS)=SUB
      K=NSUBS

         ELSE

           WRITE(*,'(''WARNING: More than''I4                           &
     &       '' different subroutine calls to TIMER'')')KLIMIT
           ISTOP=1
           GOTO 9999

         ENDIF

 10      CONTINUE

!        Start timer for subroutine
      NENTER(K)=NENTER(K)+1
! DEPENDS ON: get_wallclock_time
      ELPSTART = Get_Wallclock_Time()
! DEPENDS ON: get_cpu_time
      CPUSTART = Get_CPU_Time()
!      -----------------------------------------------------------------
       ELSEIF ((I == 4).AND.(ISTOP == 0)) THEN

!      Return from subroutine

!        Stop timer
! DEPENDS ON: get_cpu_time
      CPUEND = Get_CPU_Time()
! DEPENDS ON: get_wallclock_time
      ELPEND = Get_Wallclock_Time()
      ELAPSE(K)=ELAPSE(K)+(ELPEND-ELPSTART)*1.E-3
      TIME(K)=TIME(K)+CPUEND-CPUSTART

!        Find name of calling program
         DO K=1,NSUBS

           IF (SUBNAME(K) == RETURNAM(J)) GOTO 11

         ENDDO

         WRITE(*,'(3X,''Calling prog:-'',1X,A8,1X,''not found           &
     &     ,now in'',1X,A8/3X,''TIMER being DISABLED for the rest       &
     &     of this run'')')RETURNAM(J),SUBNAME(J+1)
      ISTOP=1
         GOTO 9999
 11      CONTINUE

!        Start timer for calling program
! DEPENDS ON: get_wallclock_time
      ELPSTART = Get_Wallclock_Time()
! DEPENDS ON: get_cpu_time
      CPUSTART = Get_CPU_Time()
      J=J-1

!      -----------------------------------------------------------------
       ELSEIF ((I <  1).OR.(I >  6)) THEN

!      If I<1 or I>6then there is an error. If 4<I<=6 then this call
!      to TIMER is ignored. These values of I are recognised by the
!      TIMER3A version.

         WRITE(*,'(3X,                                                  &
     &     ''Illegal call to TIMER by subroutine'',1X,A8)')SUB

       ENDIF

 9999    CONTINUE

       RETURN
       END SUBROUTINE TIMER


#endif
