#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: SETTSCTL -------------------------------------------------
!LL
!LL  Purpose: Sets timestep loop control switches and STASHflags.
!LL           Note STEP on entry is the values at the N+1 (ie. updated)
!LL           timelevel; control switches are therefore set using
!LL           A_STEP/O_STEP, whereas physical switches are set using
!LL           A_STEP-1/O_STEP-1 to ensure correct synchronisation.
!LL           Note also that step on entry is N (i.e. not updated) when
!CL           called from INITIAL.
!CL
!LL  Author:   T.C.Johns
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (07/9/90)
!LL
!LL  Logical components covered: C0
!LL
!LL  Project task: C0
!LL
!LL  External documentation: On-line UM document C0 - The top-level
!LL                          control system
!LL
!LLEND------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE SETTSCTL (                                             &
#include "argduma.h"
#include "argsts.h"
#include "arginfa.h"
     &       internal_model,LINITIAL,MEANLEV,ICODE,CMESSAGE,            &
     &       isSTASH)  ! FLUME-STASH 
!
! FLUME-STASH 
USE flumerun
!
      IMPLICIT NONE
!
      LOGICAL isSTASH  ! FLUME-STASH true for Flume-STASH process
      INTEGER internal_model   ! IN  - internal model identifier
      LOGICAL LINITIAL         ! IN  - true in called from INITIAL
      INTEGER MEANLEV          ! OUT - Mean level indicator
#include "cmaxsize.h"
! Needed for N_INTERNAL_MODEL
#include "csubmodl.h"
#include "parparm.h"
#include "typsize.h"
#include "typduma.h"
! Contains *CALL CPPXREF
#include "typsts.h"
#include "typinfa.h"
      INTEGER ICODE            ! Out - Return code
      CHARACTER*(80) CMESSAGE  ! Out - Return error message
!
!*----------------------------------------------------------------------
!  Common blocks
!
#include "chsunits.h"
#include "chistory.h"
#include "ccontrol.h"
#include "ctime.h"
#include "stparam.h"
#include "clookadd.h"
#include "ctfilt.h"
!
!
!  Local variables
!
      INTEGER I,NFTUNIT        ! Loop counters
      INTEGER ITIME            ! Loop index over selected dump times
      INTEGER STEP             ! A_STEP or O_STEP for atmos/ocean
      INTEGER IS,II,IL,IM,NTAB,IT,IE  ! STASH variables
      INTEGER modl             ! Int model no, read from STASH list arra
      INTEGER MODJOBREL        ! modulus of JOBREL_STEP()
      INTEGER STEP_ASSIM       ! model step rel. to assim start
!
      INTEGER      JINTF            ! Interface area index
      INTEGER   ANCIL_REF_DAYS,ANCIL_REF_SECS
      INTEGER   ANCIL_OFFSET_STEPS ! offset of ref. from basis time
      INTEGER MINS_PER_STEP    ! minutes per timestep
      INTEGER MONTHS_IN        ! Number of months into forecast

! Temporary variables for development of 3.5 start
      INTEGER    SECS_PER_DAY
      PARAMETER (SECS_PER_DAY = 24*3600)
      INTEGER A_S_Steps       ! Atmosphere steps: slab steps

      INTEGER SECS_PER_PERIOD
      INTEGER STEPS_PER_PERIOD
      INTEGER DUMPFREQ
      INTEGER OFFSET_DUMPS
      INTEGER EXITFREQ
      INTEGER TARGET_END_STEP
      INTEGER ANCILLARY_STEPS
      INTEGER BOUNDARY_STEPS
      INTEGER BNDARY_OFFSET
      INTEGER DUMPTIMES(DUMPTIMES_LEN1)
      INTEGER MEANFREQ (MEANFREQ_LEN1)
      INTEGER PRINTFREQ(PRINTFREQ_LEN1)
      INTEGER JOBREL_STEP(JOBREL_LEN1)
! Temporary variables for development of 3.5 end

      LOGICAL IAU_ResetDT ! If .TRUE., reset data time.
      LOGICAL ST_TIME_RESET ! Flag for resetting stash times for dump period means
      INTEGER, save :: RESET_STEP = 0, NEW_RESET_STEP
      INTEGER, save :: month_len, i_month_plus

      ICODE=0
! Initialise internal_model (was done in SETGRCTL)
      internal_model = 1 ! Atmosphere

!L 1. Set timestep loop top-level control switches
!L
!L 1.0  Initialise control switches which are shared in coupled runs
!L
!
      STEP=                        STEPim(internal_model)
      SECS_PER_PERIOD=  SECS_PER_PERIODim(internal_model)
      STEPS_PER_PERIOD=STEPS_PER_PERIODim(internal_model)
      mins_per_step=SECS_PER_PERIOD/(60*STEPS_PER_PERIOD)
      DUMPFREQ=                DUMPFREQim(internal_model)
      OFFSET_DUMPS=        OFFSET_DUMPSim(internal_model)
      EXITFREQ=                EXITFREQim(internal_model)
      TARGET_END_STEP=  TARGET_END_STEPim(internal_model)
      ANCILLARY_STEPS=  ANCILLARY_STEPSim(internal_model)
      BOUNDARY_STEPS =   BOUNDARY_STEPSim(internal_model)
      BNDARY_OFFSET  =   BNDARY_OFFSETim(internal_model)
      DO I=1,DUMPTIMES_LEN1
        DUMPTIMES(I)=       DUMPTIMESim(I,internal_model)
      ENDDO ! I
      DO I=1,MEANFREQ_LEN1
        MEANFREQ(I)=         MEANFREQim(I,internal_model)
      ENDDO ! I
      DO I=1,PRINTFREQ_LEN1
        PRINTFREQ(I)=       PRINTFREQim(I,internal_model)
      ENDDO ! I
      DO I=1,JOBREL_LEN1
        JOBREL_STEP(I)=   JOBREL_STEPim(I,internal_model)
      ENDDO ! I
      MEANLEV=0
!
      LASSIMILATION=.FALSE.
      LDUMP=        .FALSE.
      LEXIT=        .FALSE.
      LMEAN=        .FALSE.
      LHISTORY=     .FALSE.
      LPRINT=       .FALSE.
      LANCILLARY=   .FALSE.
      LBOUNDARY=    .FALSE.
      LINTERFACE=   .FALSE.
      LJOBRELEASE = .FALSE.
!L
!L 1.1  Set up PPfile switches for the timestep
!L
      LPP =.FALSE.
! FLUME-STASH - in a Flume run, for LINITIAL, this loop is executed 
!               only if the parallel-STASH process is active.
!             - when not LINITIAL, it is called to allow reinitialisation
!               of lateral boundary condition output files by the model.
     IF (isSTASH .OR. .not.Flume_run .OR. .not.LINITIAL) THEN

      DO NFTUNIT=20,NUNITS
        IF (FT_INPUT(NFTUNIT) == 'Y'.OR.FT_OUTPUT(NFTUNIT) == 'Y') THEN

          ! Initialise:
          LPP_SELECT(NFTUNIT) = .FALSE.

          ! Allow mix of real-month & regular reinit. periods on
          ! different units:
          IF (FT_STEPS(NFTUNIT) <  0) THEN ! Real-month reinit.

            ! Select files to be reinitialised on this timestep:
            IF (STEP == 0 .AND. FT_FIRSTSTEP(NFTUNIT) == 0) THEN

              LPP_SELECT(NFTUNIT) =.TRUE.

            ELSEIF (.NOT.LCAL360                                        &
     &              .AND.(NFTUNIT >= 60 .AND. NFTUNIT <  70)            &
     &              .AND.(I_DAY     ==  1)                              &
     &              .AND.(I_HOUR    ==  (MINS_PER_STEP)/60)             &
     &              .AND.(I_MINUTE  ==  MOD(MINS_PER_STEP,60))          &
     &              .AND.(STEP /= 1)) THEN

              MONTHS_IN = I_MONTH - MODEL_BASIS_TIME(2) +               &
                 12 * (I_YEAR - MODEL_BASIS_TIME(1))

              if (MONTHS_IN < 0) MONTHS_IN = MONTHS_IN + 12

              IF (MONTHS_IN >= FT_FIRSTSTEP(NFTUNIT) .AND.              &
     &            (FT_LASTSTEP(NFTUNIT) <  0 .OR.                       &
     &             MONTHS_IN <= FT_LASTSTEP(NFTUNIT))) THEN

                IF (FT_STEPS(NFTUNIT) == -1) THEN
                  LPP_SELECT(NFTUNIT) = .TRUE. ! Months
                ELSEIF (FT_STEPS(NFTUNIT) == -3 .AND.                   &
     &                  MOD(MONTHS_IN -                                 &
     &                      FT_FIRSTSTEP(NFTUNIT),3)  == 0) THEN
                  LPP_SELECT(NFTUNIT) = .TRUE. ! Seasons
                ELSEIF (FT_STEPS(NFTUNIT) == -12 .AND.                  &
     &                  MOD(MONTHS_IN -                                 &
     &                      FT_FIRSTSTEP(NFTUNIT),12) == 0) THEN
                  LPP_SELECT(NFTUNIT) = .TRUE. ! Years
                ENDIF ! Of FT_STEPS(NFTUNIT) = reinit period

              ENDIF ! Of MONTHS_IN within reinitialisation limits

            ENDIF  ! of lcal360 and nftunit etc.

          ELSE IF (FT_STEPS(NFTUNIT) >  0) THEN ! Regular reinit.

            ! Select files to be reinitialised on this timestep:
            IF (STEP == 0 .AND. FT_FIRSTSTEP(NFTUNIT) == 0) THEN

              LPP_SELECT(NFTUNIT) = .TRUE.

            ELSEIF ((STEP-1) >= FT_FIRSTSTEP(NFTUNIT) .AND.             &
     &              (FT_LASTSTEP(NFTUNIT) <  0 .OR.                     &
     &               (STEP-1) <= FT_LASTSTEP(NFTUNIT))) THEN

              IF (MOD((STEP-1) -                                        &
     &                FT_FIRSTSTEP(NFTUNIT),FT_STEPS(NFTUNIT)) == 0)    &
     &          LPP_SELECT(NFTUNIT) = .TRUE.

              ! Do not reinitialise files on step 1 if they were
              ! initialised at step 0:
              IF (STEP == 1 .AND. FT_FIRSTSTEP(NFTUNIT) == 0)           &
     &          LPP_SELECT(NFTUNIT) = .FALSE.

              ! Sub model id must correspond with TYPE_LETTER_2:
              LPP_SELECT(NFTUNIT) =                                     &
     &          LPP_SELECT(NFTUNIT)                                     &
     &          .AND. ( (INTERNAL_MODEL == ATMOS_IM .AND.               &
     &                   TYPE_LETTER_2(NFTUNIT) == 'a')                 &
     &                .OR.                                              &
     &                  (INTERNAL_MODEL == OCEAN_IM .AND.               &
     &                   TYPE_LETTER_2(NFTUNIT) == 'o') )

            ENDIF
!L
!L 1.1.1  Set up PPfile switches for boundary output files - this needs
!L        to be protected on SLAB model timesteps
!L
            IF (TYPE_LETTER_1(NFTUNIT) == 'b') THEN ! Boundary File

!               Get interface area number
! DEPENDS ON: intf_area
                call intf_area( internal_model, NFTUNIT, JINTF)

                LPP_SELECT(NFTUNIT) = ( LPP_SELECT(NFTUNIT)             &
     &          .AND. .NOT.                                             &
!               Do not reinitialise first file on timestep ft_firststep
     &          (STEP-1-FT_FIRSTSTEP(NFTUNIT) == 0) )                   &
     &          .OR.                                                    &
!               Initialise first file if start of sequence is offset
!               from beginning of run
     &          (STEP-1+INTERFACE_STEPSim(JINTF,internal_model) -       &
     &           FT_FIRSTSTEP(NFTUNIT)  ==  0)                          &
     &          .OR.                                                    &
!               Select boundary file if incomplete on continuation run
     &          (LINITIAL                                               &
     &           .AND.                                                  &
     &           STEP+INTERFACE_STEPSim(JINTF,internal_model) >         &
     &                                         FT_FIRSTSTEP(NFTUNIT)    &
     &           .AND.                                                  &
     &           STEP <= INTERFACE_LSTEPim(JINTF,internal_model)        &
     &           .AND.                                                  &
     &           (STEP-FT_FIRSTSTEP(NFTUNIT) == 0 .OR.                  &
     &           MOD(STEP-FT_FIRSTSTEP(NFTUNIT),FT_STEPS(NFTUNIT)) /= 0)&
     &            )

            ENDIF

          ELSE  ! for files not reinitialised, ie. ft_steps(nftunit)=0

!           Initialise at step 0
            LPP_SELECT(NFTUNIT) = STEP == 0 .AND.                       &
     &      (FT_STEPS(NFTUNIT) == 0.OR.FT_FIRSTSTEP(NFTUNIT) == 0)

!           Select boundary file if incomplete on continuation run
            IF (LINITIAL .AND.                                          &
     &          TYPE_LETTER_1(NFTUNIT) == 'b') THEN ! Boundary File

!             Get interface area number
! DEPENDS ON: intf_area
              call intf_area( internal_model, NFTUNIT, JINTF)

              LPP_SELECT(NFTUNIT) = LPP_SELECT(NFTUNIT) .OR.            &
     &        (STEP >  0.AND.STEP <=                                    &
     &                       INTERFACE_LSTEPim(JINTF,internal_model))

            ENDIF

          END IF  ! of FT_STEPS(NFTUNIT) lt, gt or =0, ie. reinit. type
        ELSE ! of FT_INPUT(NFTUNIT) or FT_OUTPUT(NFTUNIT) =Y
          LPP_SELECT(NFTUNIT)=.FALSE.
        ENDIF
        LPP = LPP .OR. LPP_SELECT(NFTUNIT)
      END DO  ! of loop over nftunit from 20 to nunits

      END IF  ! FLUME-isSTASH  
!L
!L 1.2   Set switches for general internal models.
!L       For coupled models dump related switches can only be set when
!L       the last internal model in a submodel has completed its group
!L       of timesteps. For coupled models the only safe restart point
!L       is at the completion of all groups within a model timestep.

      IF(N_INTERNAL_MODEL == 1.OR.(                                     &
                                           ! if not coupled model, or
     &   LAST_IM_IN_SM(internal_model).AND.                             &
                                            ! last model in submodel
     &   MOD(STEP,STEPS_PER_PERIOD) == 0))                              &
                                            ! and last step in group
     &   THEN
!
!  LDUMP   : Write-up dump on this timestep
!
        IF (DUMPFREQ == -9999) THEN
           ! -9999 is used to signal dumping at the end of every calendar month
           LDUMP=   .NOT. LINITIAL                                      &
     &              .AND.(I_DAY     ==  1)                              &
     &              .AND.(I_HOUR    ==  0)                              &
     &              .AND.(I_MINUTE  ==  0)          
           LDUMP = LDUMP .or. STEP==0
       ELSE IF (DUMPFREQ >  0) THEN
          LDUMP=       (MOD(STEP,DUMPFREQ)    == 0)
        ELSE
          LDUMP=.FALSE.
          DO ITIME=1,DUMPTIMES_LEN1
            LDUMP=LDUMP.OR.(STEP == DUMPTIMES(ITIME))
          ENDDO
        ENDIF

!
!  LMEAN   : Perform climate-meaning from dumps on this timestep
!  LHISTORY: Write-up history file on this timestep
!
        IF (DUMPFREQ >  0.AND.MEANFREQ(1) >  0) THEN
          LMEAN=     (MOD(STEP,DUMPFREQ)       == 0)
          LHISTORY=  (MOD(STEP+OFFSET_DUMPS*DUMPFREQ,                   &
     &                    DUMPFREQ*MEANFREQ(1)) == 0)
        ELSE
          LMEAN=     .FALSE.
          LHISTORY=   LDUMP
        ENDIF
!             For coupled models, only allow history write-ups if model
!             timestep complete, ie last internal model.
      IF(LHISTORY.AND.N_INTERNAL_MODEL >  1.AND.                        &
                                                       ! Coupled model
     &   internal_model /= INTERNAL_MODEL_LIST(N_INTERNAL_MODEL)) THEN
         LHISTORY=   .FALSE.
      ENDIF               ! Test on last internal model in coupled
!
!  LPRINT  : Write a formatted print from dump on this timestep
!
        IF (LDUMP) THEN
          IF (DUMPFREQ >  0.AND.PRINTFREQ(1) >  0) THEN
            LPRINT=  (MOD(STEP,DUMPFREQ*PRINTFREQ(1)) == 0)
          ELSE
            LPRINT=  .FALSE.
          ENDIF
        ELSE
          LPRINT=    .FALSE.
        ENDIF
!
!  LEXIT   : Check for exit condition on this timestep
!
        IF (EXITFREQ >  0) THEN ! Implies climate meaning

          LEXIT=  ( (MOD(STEP,EXITFREQ)    == 0)  .OR.                  &
     &              (STEP  >=  TARGET_END_STEP) )

        ELSE                    ! No climate meaning

          LEXIT=    (STEP  >=  TARGET_END_STEP)

        ENDIF
!
!  LANCILLARY: Update ancillary fields on this timestep
!
!L  Convert ancillary reference time to days & secs
! DEPENDS ON: time2sec
        CALL TIME2SEC(ANCIL_REFTIME(1),ANCIL_REFTIME(2),                &
     &                ANCIL_REFTIME(3),ANCIL_REFTIME(4),                &
     &                ANCIL_REFTIME(5),ANCIL_REFTIME(6),                &
     &                    0,0,ANCIL_REF_DAYS,ANCIL_REF_SECS,LCAL360)

!L  Compute offset in timesteps of basis time from ancillary ref.time
! DEPENDS ON: tim2step
        CALL TIM2STEP(BASIS_TIME_DAYS-ANCIL_REF_DAYS,                   &
     &                BASIS_TIME_SECS-ANCIL_REF_SECS,                   &
     &       STEPS_PER_PERIOD,SECS_PER_PERIOD,ANCIL_OFFSET_STEPS)

        IF (ANCILLARY_STEPS >  0.AND.STEP >  0)                         &
     &LANCILLARY=(MOD(STEP+ANCIL_OFFSET_STEPS,ANCILLARY_STEPS) == 0)

      ENDIF     ! Test for non-coupled or coupled + last step in group

!
!  LBOUNDARY    : Update boundary fields on this timestep
!
      IF (BOUNDARY_STEPS >  0)                                          &
     &    LBOUNDARY= (MOD(STEP+BNDARY_OFFSET,BOUNDARY_STEPS)  == 0)
!
!  LJOBRELEASE  : Release user jobs on this timestep
!
        DO I=1,JOBREL_LEN1
         MODJOBREL=ABS(JOBREL_STEP(I))
         IF (MODJOBREL /= 0)THEN
          LJOBRELEASE = LJOBRELEASE .OR.                                &
     &      (STEP == JOBREL_STEP(I) .AND. JOBREL_STEP(I) >  0).OR.      &
     &    (MOD(STEP,MODJOBREL) == 0.AND. JOBREL_STEP(I) <  0)
         ENDIF
        ENDDO


!L
!L 1.2.1  Set switches for atmosphere timestep
!L
#if defined(ATMOS)
      IF (internal_model == atmos_im) THEN
!L 1.2.2 Set switches for all cases

! Energy correction switches
! Set on the assumption that the energy correction is evaluated at the
! end of a timestep
        IF (A_ENERGYSTEPS >  0) THEN
! true if this is the step on which to do energy calculation
          LENERGY =  (MOD(step,a_energysteps) == 0)
! True if this is the step after last energy calculation
          LFLUX_RESET = (MOD(step-1,a_energysteps) == 0)
        ENDIF

!
!  LASSIMILATION: Perform data assimilation on this timestep
!
        LASSIMILATION= (STEP >  ASSIM_FIRSTSTEPim(a_im) .AND.           &
     &   (STEP-ASSIM_FIRSTSTEPim(a_im)  <=                              &
     &         ASSIM_STEPSim(a_im)+ASSIM_EXTRASTEPSim(a_im)))
        LASSIMILATION=LASSIMILATION.AND.                                &
     &               ( MODEL_ASSIM_MODE == 'Atmosphere' .OR.            &
     &                 MODEL_ASSIM_MODE == 'Coupled   ')
!L
!L      Reset Data Time fields in dump header at new analysis time
!L      ( also in history block )
!L      ( This now includes resetting prognostic field LOOKUP headers )
!L      NB: done even if assimilation suppressed by RUN_ASSIM_MODE
!L      Set MEANLEV to -1  and LDUMP to TRUE to force dump of analysis
!L      - otherwise set MEANLEV to 0 (ie. instantaneous)
!L
        IAU_ResetDT = .FALSE.
        IF (L_IAU .OR. RUN_ASSIM_MODE  ==  "NoIAU     ") THEN
          IF (STEP  ==  IAU_DTResetStep) IAU_ResetDT = .TRUE.
        END IF
        IF ( (STEP == ASSIM_FIRSTSTEPim(a_im)+ASSIM_STEPSim(a_im)       &
     &        .AND.LASSIMILATION) .OR. IAU_ResetDT) THEN

          DO I=21,27
            A_FIXHD(I)=A_FIXHD(I+7)
          ENDDO
          DO I=1,6
            MODEL_DATA_TIME(I)=A_FIXHD(20+I)
          ENDDO
          DO I=1,A_PROG_LOOKUP
            A_LOOKUP(LBYRD ,I)=A_FIXHD(21)
            A_LOOKUP(LBMOND,I)=A_FIXHD(22)
            A_LOOKUP(LBDATD,I)=A_FIXHD(23)
            A_LOOKUP(LBHRD ,I)=A_FIXHD(24)
            A_LOOKUP(LBMIND,I)=A_FIXHD(25)
            A_LOOKUP(LBDAYD,I)=A_FIXHD(27)
          ENDDO
          MEANLEV=-1
        ELSE
          MEANLEV=0
        ENDIF
!L
!L      Suppress assimilation using RUN_ASSIM_MODE if necessary
        LASSIMILATION=LASSIMILATION.AND.                                &
     &               ((RUN_ASSIM_MODE == 'Atmosphere').OR.              &
     &                (RUN_ASSIM_MODE == 'Coupled   '))
!
!  L_SW_RADIATE : Perform SW-radiation on this timestep
!  L_LW_RADIATE : Perform LW-radiation on this timestep
!
        IF (A_SW_RADSTEP >  0)                                          &
     &    L_SW_RADIATE=(MOD(STEP-1,A_SW_RADSTEP)  == 0)
        IF (A_LW_RADSTEP >  0)                                          &
     &    L_LW_RADIATE=(MOD(STEP-1,A_LW_RADSTEP)  == 0)
        IF (A_SW_RADSTEP_DIAG >  0)                                     &
     &    L_SW_RADIATE_DIAG=(MOD(STEP-1,A_SW_RADSTEP_DIAG)  == 0)
        IF (A_LW_RADSTEP_DIAG >  0)                                     &
     &    L_LW_RADIATE_DIAG=(MOD(STEP-1,A_LW_RADSTEP_DIAG)  == 0)

        IF (A_SW_RADSTEP_PROG >  0)                                     &
     &    L_SW_RADIATE_PROG=(MOD(STEP-1,A_SW_RADSTEP_PROG)  == 0)
        IF (A_LW_RADSTEP_PROG >  0)                                     &
     &    L_LW_RADIATE_PROG=(MOD(STEP-1,A_LW_RADSTEP_PROG)  == 0)
!
!  LINTERFACE: Output interface fields on this timestep (>1 file poss.)
!
        DO JINTF=1,N_INTF_A
          IF (INTERFACE_STEPSim(JINTF,atmos_im) >  0) THEN
            LINTERFACE = LINTERFACE .OR.                                &
     &      (MOD(STEP-INTERFACE_FSTEPim(JINTF,atmos_im),                &
     &           INTERFACE_STEPSim(JINTF,atmos_im)) == 0                &
     &      .AND. STEP >= INTERFACE_FSTEPim(JINTF,atmos_im)             &
     &      .AND. STEP <= INTERFACE_LSTEPim(JINTF,atmos_im) )
! To allow output of interface files from model,
! even if Flume parallel-STASH output is being used for fieldsfiles.
            IF (LINITIAL) THEN
              LPP_SELECT(JINTF+139) = .TRUE.
            END IF
          ENDIF
        ENDDO
      ENDIF ! Test for atmos_im
#endif
!L
!L----------------------------------------------------------------------
!L 2. Set STASHflags to activate diagnostics this timestep
!L
!     IS is section number
!     IE is item number within section
!     IL is item number within STASHlist
!     II is counter within given section/item sublist for repeated items
!     IM is cumulative active item number within section
!
!L  Clear all STASHflags

      DO IS=0,NSECTS
      DO IM=0,NITEMS
        SF(IM,IS) = .FALSE.
      ENDDO
      ENDDO

!L  Loop over all items in STASHlist, enabling flags for diagnostics
!L  which are active this step --
!L    note that atmosphere/ocean diagnostics must be discriminated

      IF (DUMPFREQ == -9999                                             &
         ! -9999 is used to signal dumping at the end of every calendar month
         ! Test on iday==2 here so don't muck up the end of the previous month
         ! Assumes period is set larger than one day initially.
     &              .AND.(I_DAY     ==  1)                              &
     &              .AND.(I_HOUR    ==  0)                              &
     &              .AND.(I_MINUTE  ==  0) ) THEN
         ST_TIME_RESET = .TRUE.
         if ( linitial ) then
            ! So works for restarts as well
            reset_step = step
         end if
         ! Get length of this month in days. Use month+1 because of offset in
         ! setperlen
! DEPENDS ON: setperlen
         if (lcal360) then
            month_len = 30
         else
            i_month_plus = i_month+1
            if (i_month_plus == 13) then
               i_month_plus = 1
            end if
            call setperlen(1,i_month_plus,i_year,month_len)
         end if
         ! Convert length to time steps
         month_len = month_len * 1440 / mins_per_step
         new_reset_step = step + month_len
      ELSE
         ST_TIME_RESET = .FALSE.
      ENDIF



      DO IL =1,TOTITEMS

        if  ( DUMPFREQ == -9999 .and.  step > reset_step .and.           &
              STLIST(st_macrotag,il)==-999 ) then
           ! If using end of month dump and this is a dump freq variable
           ! change the period
           STLIST(st_period_code,IL) = month_len
           STLIST(st_start_time_code,IL) = step
        end if
        modl=STLIST(st_model_code  ,IL)
        IS  =STLIST(st_sect_no_code,IL)
        IM  =STLIST(st_item_code   ,IL)
! Skip diagnostics which are not active
        IF(STLIST(st_proc_no_code,IL) == 0) GOTO 200
! Skip diagnostics which don't correspond to the submodel this timestep
!       IF(ISUBMODL /= modl) GOTO 200
        IF(Internal_model /= modl) GOTO 200
!
!     STASHflag is off by default.
!     But reset ...
!
!       ... if required every timestep between start and end
!
        IF (STLIST(st_freq_code,IL) == 1) THEN
          IF (STEP >= STLIST(st_start_time_code,IL).AND.                &
     &       (STEP <= STLIST(st_end_time_code,IL).OR.                   &
     &        STLIST(st_end_time_code,IL) == st_infinite_time))         &
     &    THEN
            SF(IM,IS)=.TRUE.
            SF(0 ,IS)=.TRUE.
          ENDIF
!
!       ... if required at specified times and this is one of them
!
        ELSEIF(STLIST(st_freq_code,IL) <  0) THEN
           IF (STLIST(st_freq_code,IL)==-9999) THEN
              IF (.not. LINITIAL)THEN
                 ! End of month dump
                 SF(IM,IS)=SF(IM,IS).or.LDUMP
                 SF(0 ,IS)=SF(0 ,IS).or.LDUMP
              ENDIF
           ELSE
              NTAB=-STLIST(st_freq_code,IL)
              DO IT=1,NSTTIMS
                 IF(STTABL(IT,NTAB) == st_end_of_list) GOTO 200
                 IF(STEP == STTABL(IT,NTAB)) THEN
                    SF(IM,IS)=.TRUE.
                    SF(0 ,IS)=.TRUE.
                 ENDIF
              ENDDO
           ENDIF
!        
!       ... if required every N timesteps and this is one of them
!
        ELSEIF (STLIST(st_freq_code,IL) >  0) THEN
          IF  (MOD((STEP-STLIST(st_start_time_code,IL)),                &
     &              STLIST(st_freq_code,IL)) == 0.AND.                  &
     &         STEP >= STLIST(st_start_time_code,IL).AND.               &
     &        (STEP <= STLIST(st_end_time_code,IL).OR.                  &
     &         STLIST(st_end_time_code,IL) == st_infinite_time))        &
     &    THEN
            SF(IM,IS)=.TRUE.
            SF(0 ,IS)=.TRUE.
          ENDIF
        ENDIF
!
!     Next item
!
  200 CONTINUE
      ENDDO
      if  ( step > reset_step ) then
         reset_step = new_reset_step
      end if
!L----------------------------------------------------------------------
!L 3. If errors occurred in setting STASHflags, set error return code
      IF (ICODE == 1) CMESSAGE='SETTSCTL: STASH code error'
 9999 CONTINUE
      RETURN
!L----------------------------------------------------------------------
      END SUBROUTINE SETTSCTL
#endif
