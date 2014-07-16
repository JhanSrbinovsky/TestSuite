#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: DUMPCTL --------------------------------------------------
!LL
!LL  Purpose: Controls the production and naming of output dump files.
!LL           Also selectively adds dump files to the list of dumps
!LL           for processing by the external dump server process.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Author:   T.C.Johns
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: C2
!LL
!LL  Project task: C2
!LL
!LL  External documentation: On-line UM document C0 - The top-level
!LL                          control system; On-line document C2 -
!LL                          Dump Handling.
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE DUMPCTL (                                              &
#include "argd1.h"
#include "argduma.h"
#include "arglndm.h"
#include "argptra.h"
#include "argsts.h"
#include "argppx.h"
     &           I_AO,MEANLEV,lwritd1,tmpfname,writestep,ICODE,CMESSAGE)
!
      IMPLICIT NONE

#include "parvars.h"

!*L Arguments
#include "cmaxsize.h"
#include "csubmodl.h"
#include "typsize.h"
#include "nstypes.h"
#include "typd1.h"
#include "typduma.h"
#include "typlndm.h"
#include "typptra.h"
#include "typsts.h"
#include "ppxlook.h"
!
      INTEGER I_AO             ! IN  - Atmosphere/Ocean indicator
      INTEGER MEANLEV          ! IN  - Mean period level for dump
      INTEGER writestep        ! IN  - Timestep on which to write it
      CHARACTER*14 tmpfname    ! IN  - Name of file to be written
                               !       containing temporary copy of D1
      LOGICAL lwritd1          ! IN  - True if doing a WRITD1 temp write

      INTEGER ICODE            ! OUT - Error return code
      CHARACTER*80 CMESSAGE
!
!*----------------------------------------------------------------------
!  Common blocks
!
#include "chsunits.h"
#include "chistory.h"
#include "ccontrol.h"
#include "cenvir.h"
#include "clookadd.h"
#include "c_mdi.h"
#include "ctime.h"
#include "ctfilt.h"
! Including cprintst.h so we can use PrintStatus Variable
#include "cprintst.h"
!
!  Local variables
!
      LOGICAL LARCHIVE         ! WORK  - Switch for dump archiving
      LOGICAL LUNITTYPE        ! True if unit number can be superceded
!                              ! dump or ppfile
      LOGICAL LKEEPATM         ! True if atmos safe restart dump needs
!                              ! to be kept until the ocean safe
!                              ! restart dump is produced
      LOGICAL LDELATM          ! True if atmos safe restart dump can
!                              ! now be deleted ie. ocean safe restart
!                              ! dump has been produced.
      CHARACTER*80 FILENAME
      CHARACTER*1     C1
      INTEGER MYEAR,MMONTH,MDAY,MHOUR,MMIN,MSEC  ! Creation date
      INTEGER     I,J              ! Loop counters
!
      INTEGER       NFTOUT,                                             &
                                ! Output unit number
     &              BUFLEN                                              &
                                ! Length of i/o buffer for WRITDUMP
     &             ,TOGGLE                                              &
                                ! Dummy argument for GET_NAME
     &             ,REINIT_STEPS                                        &
                                ! Dummy argument for GET_NAME
     &,LEN_DUMPNAME                                                     &
                     !No of characters in file name
     &,ERROR                                                            &
                     !Error code returned by OPEN
     &,STEP                                                             &
                     !Step number
     &,archdump_monfreq                                                 &
                           ! dump archiving frequency (in months)
     &,archdump_monoffset                                               &
                           ! dump archiving offset (in months)
     &,disk_address                                                     &
                                       ! Current rounded disk address
     &,number_of_data_words_on_disk                                     &
                                       ! Number of data words on disk
     &,number_of_data_words_in_memory  ! Number of Data Words in memory
      INTEGER internal_model
      INTEGER im      ! temporary internal model id for ocean or slab


      CHARACTER*1   FILETYPE    ! Code letter for file type
      CHARACTER*1   LETTER_3    ! dummy argument for GET_NAME
      CHARACTER*14  DUMPNAME    ! Model generated dump name
      INTEGER STP1im(N_INTERNAL_MODEL_MAX)!NO OF STEPS SINCE im DUMP
      INTEGER STP2im(N_INTERNAL_MODEL_MAX)!NO OF STEPS BETWEEN im
                                          !PERIOD 1 MEANS
      INTEGER D1_ADDR_SUBMODEL_ID  ! submodel id in D1_ADDR array
!L
!L----------------------------------------------------------------------
! Get name of pipe
      CALL GET_FILE(8,FILENAME,80,ICODE)
!L
!L 1. (Temporary correction in preparation for removing temporary
!L    history copy from dump headers, since not used within model
!L    and history file size now (from vn3.1 on) exceeds reserved space).
!L    Write missing data indicators into dump copy of temporary history
!L    block to prevent earlier overwriting error.
!L
#if defined(ATMOS)
      IF (I_AO == 1) THEN
        DO I=1,LEN_DUMPHIST+1
          A_DUMPHIST(I) = RMDI
        ENDDO
      ENDIF
#endif
!L----------------------------------------------------------------------
!L 2. Set LOOKUP header data and validity times from FIXHD
!L
!L ** This section deleted due to conflict with STASH use of dump
!L ** LOOKUP headers to retain timestamp information across timesteps.
!L
!L----------------------------------------------------------------------
!L 3. Construct dump name from model information using defined
!L     naming convention.
!L
      IF (MEANLEV /= -1) THEN             ! Analyses already named
        IF (lwritd1) THEN
          STEP = STEPim(I_AO)
          IF (STEP  /=  writestep) GOTO 999
          WRITE(DUMPNAME,1011) STEP
 1011     FORMAT('..........',i4.4)
          DO I=1,10
            IF(TMPFNAME(I:I) /= ' ') DUMPNAME(I:I)=TMPFNAME(I:I)
          END DO
        ELSE
!
        FILETYPE='d'                      ! Indicates dump
        TOGGLE=1
        REINIT_STEPS=0                    ! default dummy value
        LETTER_3='a'                      ! default dummy value
      IF ( I_AO == 1 .AND. (MODEL_ASSIM_MODE == "Atmosphere" .OR.       &
     &                      MODEL_ASSIM_MODE == "Coupled   " .OR.       &
     &                      RUN_ASSIM_MODE   == "NoIAU     " .OR.       &
     &                      L_IAU ) ) THEN
        ! UM6.5 : MODEL_ANALYSIS_HRS replaced by MODEL_ANALYSIS_MINS
        MODEL_ANALYSIS_HRS = REAL(MODEL_ANALYSIS_MINS)/60.0
! DEPENDS ON: get_name
        CALL GET_NAME(EXPT_ID,JOB_ID,I_AO,MEANLEV,TOGGLE,               &
     &       REINIT_STEPS,FILETYPE,LETTER_3,MODEL_STATUS,               &
     &       TIME_CONVENTION,MODEL_ANALYSIS_HRS,DUMPNAME,ICODE,CMESSAGE,&
     &       LCAL360)
      ELSE
! DEPENDS ON: get_name
        CALL GET_NAME(EXPT_ID,JOB_ID,I_AO,MEANLEV,TOGGLE,               &
     &       REINIT_STEPS,FILETYPE,LETTER_3,MODEL_STATUS,               &
     &       TIME_CONVENTION,0,DUMPNAME,ICODE,CMESSAGE,LCAL360)
      ENDIF
      IF (ICODE >  0) GOTO 999
        ENDIF
      ELSE
! Initialise dumpname to prevent problems later.
        DUMPNAME='    '
!
      ENDIF
!L----------------------------------------------------------------------
!L 4. Assign dump name to appropriate IO unit and open for write.
!L    NB: Analysis dumps have preassigned names
!L
      IF (MEANLEV == 0 .OR. lwritd1) THEN
!      Cater for instantaneous dumps
        IF(I_AO == 1)THEN
          NFTOUT=22    ! Atmos
        ENDIF
      ELSEIF (MEANLEV == -1) THEN
!      Cater for analysis dumps
        IF(I_AO == 1)THEN
          NFTOUT=28    ! Atmos
        ENDIF
      ELSE
!      Cater for mean dump
        NFTOUT=27
      ENDIF

! Flush all pp files' buffers to ensure that data is consistent
! at restart points.
! DEPENDS ON: flush_all_pp
      Call flush_all_pp()

#if defined(ATMOS)
!--compute the new addresses and lengths
      if(i_ao == 1) then
! DEPENDS ON: set_dumpfile_address
        call set_dumpfile_address(                                      &
     &   a_fixhd, len_fixhd,                                            &
     &   a_lookup, len1_lookup, a_len2_lookup,                          &
     &   number_of_data_words_in_memory, number_of_data_words_on_disk,  &
     &   disk_address)
      endif
#endif

!--output the new length of this dumpfile
#if defined(MPP)
      if(mype == 0) then
#endif
      if(PrintStatus >= PrStatus_Diag) then
        write(6,9921) trim(dumpname), nftout, disk_address
        write(0,9921) trim(dumpname), nftout, disk_address
9921    format(/'Dumpfile Size for File ',a,' on Unit ',i4,             &
     &   ' to be set to ',i10,' Words')
      end if
        call set_dumpfile_length(nftout, disk_address)
#if defined(MPP)
      endif
#endif
!L
!L 4.1 Open unit for dump : different call required if an analysis
!L     since name pre-assigned through environment variable
!L
      IF (MEANLEV /= -1) THEN
      WRITE(6,*)"DUMPCTL: Opening new file ",DUMPNAME," on unit ",NFTOUT
      LEN_DUMPNAME=LEN(DUMPNAME)
! DEPENDS ON: file_open
      CALL FILE_OPEN(NFTOUT,DUMPNAME,LEN_DUMPNAME,1,1,ERROR)
      IF (ERROR /= 0) GOTO 900
        ICODE=0

      ELSE

! DEPENDS ON: file_open
        CALL FILE_OPEN(NFTOUT,FT_ENVIRON(NFTOUT),                       &
     &            LEN_FT_ENVIR(NFTOUT),1,0,ERROR)
        IF(ERROR /= 0) GOTO 900
        ICODE=0
      ENDIF
!L----------------------------------------------------------------------
!L 5. Write dump on appropriate unit putting timestamp in header
!L
#if defined(ATMOS)
      IF (I_AO == 1) THEN

! Creation date and time

        CALL DATE_TIME(A_FIXHD(35),A_FIXHD(36),A_FIXHD(37),             &
     &  A_FIXHD(38),A_FIXHD(39),A_FIXHD(40))

! Maximum length of field, required for IO buffer

        BUFLEN=A_LOOKUP(LBLREC,1)
        IF (A_LEN2_LOOKUP >  1) THEN
          DO I=2,A_LEN2_LOOKUP
            BUFLEN=MAX(BUFLEN,A_LOOKUP(LBLREC,I))
          ENDDO
        ENDIF

        IF (MEANLEV >  0) A_FIXHD(5)=2    ! Set FIXHD(5) for mean dump

      D1_ADDR_SUBMODEL_ID = SUBMODEL_FOR_SM(atmos_sm)

! DEPENDS ON: um_writdump
        CALL UM_WRITDUMP(NFTOUT,A_FIXHD,LEN_FIXHD,                      &
     &                A_INTHD,A_LEN_INTHD,                              &
     &                A_REALHD,A_LEN_REALHD,                            &
     &                A_LEVDEPC,A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,          &
     &                A_ROWDEPC,A_LEN1_ROWDEPC,A_LEN2_ROWDEPC,          &
     &                A_COLDEPC,A_LEN1_COLDEPC,A_LEN2_COLDEPC,          &
     &                A_FLDDEPC,A_LEN1_FLDDEPC,A_LEN2_FLDDEPC,          &
     &                A_EXTCNST,A_LEN_EXTCNST,                          &
     &                A_DUMPHIST,LEN_DUMPHIST,                          &
     &                A_CFI1,A_LEN_CFI1,                                &
     &                A_CFI2,A_LEN_CFI2,                                &
     &                A_CFI3,A_LEN_CFI3,                                &
     &                A_LOOKUP,LEN1_LOOKUP,A_LEN2_LOOKUP,               &
#if defined(MPP)
     &                A_MPP_LOOKUP,MPP_LEN1_LOOKUP,                     &
#endif
     &                BUFLEN,                                           &
#include "argppx.h"
     &                atmos_sm,                                         &
     &                NO_OBJ_D1(D1_ADDR_SUBMODEL_ID),                   &
     &                D1_ADDR(1,1,D1_ADDR_SUBMODEL_ID),                 &
     &                A_LEN_DATA,D1)

        A_FIXHD(5)=1    ! Set FIXHD(5) back to instantaneous dump

        IF (MEANLEV == 0) THEN
          ARESTART='ARESTART: $DATAM/'
          ARESTART(18:31)=DUMPNAME
        ELSEIF (MEANLEV == -1) THEN
          ARESTART='ARESTART: '
          ARESTART(11:80)=ATMANL(11:80)
        ELSEIF (NFTOUT == 27) THEN
! Special case for mean dump file
          AOMEAN = 'AOMEAN  : $DATAM/'
          AOMEAN(18:31) = DUMPNAME
        ENDIF

      ENDIF

#endif
!L
!L 5.1 Close unit
!L

      IF(MEANLEV /= -1) THEN    ! Not analysis

      LEN_DUMPNAME=LEN(DUMPNAME)
! DEPENDS ON: file_close
        CALL FILE_CLOSE(NFTOUT,DUMPNAME,LEN_DUMPNAME,1,0,ICODE)
      ELSE                      ! analysis

! DEPENDS ON: file_close
        CALL FILE_CLOSE(NFTOUT,FT_ENVIRON(NFTOUT),LEN_FT_ENVIR(NFTOUT), &
     &  0,0,ICODE)
      ENDIF

!--now set the current length back to zero after we have done
!  the o/p
#if defined(MPP)
      if(mype == 0) then
#endif
        call set_dumpfile_length(nftout, 0)
#if defined(MPP)
      endif
#endif

!LL 5.2 Exit here for WRITD1 temporary writes of dumps since
!LL     no requests need to be sent to the archive server

      IF (lwritd1) GOTO 999

!L----------------------------------------------------------------------
!L 6. Construct dump processing requests and send to slave task
!L
!L 6.1 Delete previous last-restart-dump from disk (slave request),
!L     and update last-restart-dump to be current dump
!L     (exception is: first dump not to be deleted if operational)
!L

      IF (NFTOUT == 22.OR.NFTOUT == 42.OR.(NFTOUT >  59.AND.            &
     &    NFTOUT <  68)) THEN
!       Instantaneous dump (atmos/ocean) or daily PP file? (add
!       extra unit numbers 68 and 69 if so but PP files not used here)
        LUNITTYPE=.TRUE.
      ELSE
!       All other files
        LUNITTYPE=.FALSE.
      ENDIF

      LKEEPATM = .FALSE.
      LDELATM  = .TRUE.


      IF (MEANLEV <= 0) THEN  ! Instantaneous dump or analysis

#if defined(ATMOS)
        IF (I_AO == 1) THEN   ! Atmos submodel

          IF (MEANLEV == 0.AND.FT_SELECT(NFTOUT) == "Y") THEN
!           Meaning not switched on and files on this unit number
!           to be deleted

            IF ((internal_model  == ocean_im)  .AND.                    &
     &          (steps_per_periodim(a_im)  /=  dumpfreqim(a_im))) THEN
!             Coupled model where coupling frequency is not the
!             same as the dump frequency.  Therefore always need to
!             keep the atmos dump until the ocean dump has been
!             written at the next dump time, for restartability.
              LKEEPATM=.true.
              CURRATMim(A_IM)=DUMPNAME
            ELSE
              LKEEPATM=.false.
            ENDIF


            IF (MEANFREQim(1,A_IM)  /=  0 ) THEN
!             Period 1 means are being calculated in this run

              STP1im(A_IM)=STEPim(A_IM)+(DUMPFREQim(A_IM)*              &
     &                             (OFFSET_DUMPSim(A_IM)-1))
              STP2im(A_IM)=MEANFREQim(1,A_IM)*DUMPFREQim(A_IM)

              IF (MOD(STP1im(A_IM),STP2im(A_IM))  ==  0) THEN
!               The last dump was at a period 1 mean point so a
!               new safe restart dump was created.  Therefore do
!               not delete this latest safe restart dump but the
!               previous safe restart dump in END_DUMPim

                IF ((internal_model  == ocean_im)  .AND.                &
     &            (steps_per_periodim(a_im)  /=  dumpfreqim(a_im)))THEN
                  LASTATMim(A_IM)=SAFEDMPim(A_IM)
                ELSE
                  END_DUMPim(A_IM)=SAFEDMPim(A_IM)
                ENDIF
              ENDIF ! (MOD(STP1im(A_IM),STP2im(A_IM))  ==  0)

            ENDIF ! (MEANFREQim(1,A_IM)  /=  0 )

#if defined(MPP)
            IF (mype == 0) THEN

#endif
            IF(MEANFREQim(1,A_IM)  /=  0 ) THEN ! Period1 means are on,
                                               ! so safe to use STP2im.
            IF (.NOT.(MOD((STP1im(A_IM)+DUMPFREQim(A_IM)),STP2im(A_IM)) &
     &            ==  0)) THEN
              IF (END_DUMPim(A_IM) /= "              " .AND.            &
     &           .NOT. LKEEPATM ) THEN
!               Filename to be deleted is not blank and is not to
!               be kept until the ocean dump for the current step
!               is written
                WRITE (8,610) END_DUMPim(A_IM) ! Delete request
#if defined(T3E)
! DEPENDS ON: um_fort_flush
                call um_fort_flush(8, icode)
#else
                CLOSE(8)
                OPEN(8,FILE=FILENAME)
#endif
              ENDIF
            ENDIF
            ELSE ! MEANFREQim=0 so climate means are off. No STP2im.
              IF (END_DUMPim(A_IM) /= "              " .AND.            &
     &           .NOT. LKEEPATM ) THEN
                WRITE (8,610) END_DUMPim(A_IM) ! Write delete request
#if defined(T3E)
! DEPENDS ON: um_fort_flush
                call um_fort_flush(8, icode)
#else
                CLOSE(8)
                OPEN(8,FILE=FILENAME)
#endif
              ENDIF   ! END_DUMPim  /=  blank
            ENDIF ! of the MEANFREQim(1,A_IM) test.
#if defined(MPP)

            ENDIF  !  (mype  ==  0)
#endif
          ENDIF  ! (MEANLEV == 0.AND.FT_SELECT(NFTOUT) == "Y")


          IF (H_STEPim(A_IM) /= DUMPTIMESim(1,A_IM).OR.                 &
     &        MODEL_STATUS /= "OPERATIONAL   ") THEN
!           This step is not the first one in the list of dumptimes
!           or the model is not operational

            IF (MEANFREQim(1,A_IM) /= 0.AND.FT_SELECT(NFTOUT) == "Y")   &
     &      THEN
!             Period 1 Means switched on and files on this unit
!             number to be deleted

              IF (MOD((STP1im(A_IM)+DUMPFREQim(A_IM)),STP2im(A_IM))     &
     &            ==  0) THEN
!               File is at a Period 1 mean point.  Update names of
!               old and new safe restart points.
                SAFEDMPim(A_IM)=NEWSAFEim(A_IM)
                NEWSAFEim(A_IM)=DUMPNAME
                LASTDMPim(A_IM) = END_DUMPim(A_IM)
              ELSE
                END_DUMPim(A_IM)=DUMPNAME
              ENDIF
      ELSE
      END_DUMPim(A_IM)=DUMPNAME
            ENDIF

          ENDIF ! (H_STEPim(A_IM) /= DUMPTIMESim(1,A_IM).OR
                ! MODEL_STATUS /= "OPERATIONAL   ")

        ENDIF  ! (I_AO == 1)
#endif
 610    FORMAT('%%% ',A14,' DELETE')

      ENDIF ! (MEANLEV <= 0)

!
!L 6.2 If current dump is to be archived construct archiving request
!L     followed by delete request if appropriate
!L
      LARCHIVE = .FALSE.

      IF (I_AO == 1 .or. I_AO == 2 .or. I_AO == 4) THEN

        IF (MEANLEV <= 0 .AND. DUMPFREQim(I_AO) >  0                    &
     &      .AND. ARCHDUMP_FREQim(I_AO) >  0) THEN
!         No meaning, regular dumping, regular archiving

          if (lclimrealyr) then ! get freq + offset in terms of real
                                ! months from 360d-style freq + offset
            if (i_day  ==  1 .and. i_hour  ==  0) then ! end of month
              archdump_monfreq=(DUMPFREQim(I_AO)*ARCHDUMP_FREQim(I_AO)* &
     &                         SECS_PER_PERIODim(I_AO))/                &
     &                         (30*86400*STEPS_PER_PERIODim(I_AO))
            archdump_monoffset=(DUMPFREQim(I_AO)*ARCHDUMP_OFFSETim(I_AO)&
     &                         *SECS_PER_PERIODim(I_AO))/               &
     &                         (30*86400*STEPS_PER_PERIODim(I_AO))
!  N.B. i_month is used below, not (i_month-1), because offset is from
!  _start_ of 1st month
              if(mod((i_month-(MODEL_BASIS_TIME(2)+archdump_monoffset)),&
     &               archdump_monfreq)  ==  0) then
                LARCHIVE = .true.
              endif
            endif
          else IF ( H_STEPim(I_AO)/DUMPFREQim(I_AO)  >=                 &
     &         ARCHDUMP_OFFSETim(I_AO) ) THEN
!           Have passed the timestep from which to start archiving
!           Calculate whether this timestep is an archive time

            LARCHIVE = (MOD((H_STEPim(I_AO)/                            &
     &                 DUMPFREQim(I_AO)-ARCHDUMP_OFFSETim(I_AO)),       &
     &                  ARCHDUMP_FREQim(I_AO)) == 0)
          ENDIF

        ELSEIF (MEANLEV <= 0 .AND. DUMPFREQim(I_AO) == 0) THEN
!         No meaning, no regular dumping.  Will archive if archiving
!         frequency is greater than 0.
          LARCHIVE= (ARCHDUMP_FREQim(I_AO) >  0)

        ELSEIF (MEANLEV >  0) THEN
!         Meaning on. Archive is mean archive frequncy > 0.
          LARCHIVE= (MEANARCHim(MEANLEV,I_AO) == 1)
        ENDIF
#if defined(MPP)
        IF (mype == 0) THEN
#endif
          IF (LARCHIVE) THEN
!           Archiving turned on
            WRITE(8,620) DUMPNAME  ! archive request
#if defined(T3E)
! DEPENDS ON: um_fort_flush
            call um_fort_flush(8, icode)
#else
      CLOSE(8)
            OPEN(8,FILE=FILENAME)
#endif
            IF (MEANLEV >  0) THEN
!             Meaning turned on
              IF (LUNITTYPE) THEN
!               Correct sort of unit to be deleted
                IF (FT_SELECT(NFTOUT) == "Y") THEN
!                 Files on this unit number is to be deleted
                  WRITE(8,610) DUMPNAME  ! Delete request
                ENDIF
              ELSE
                WRITE(8,610) DUMPNAME ! Delete request
              ENDIF
#if defined(T3E)
! DEPENDS ON: um_fort_flush
              call um_fort_flush(8, icode)
#else
              CLOSE(8)
              OPEN(8,FILE=FILENAME)
#endif

            ENDIF

          ELSE  ! LARCHIVE false

            IF (MEANLEV >  0) THEN
!             Meaning turned on
              IF (LUNITTYPE) THEN
!               Correct sort of unit to be deleted
                IF (FT_SELECT(NFTOUT) == "Y") THEN
!                 Files on this unit number is to be deleted
                  WRITE(8,610) DUMPNAME ! Delete request
                ENDIF
              ELSE
                WRITE(8,610) DUMPNAME ! Delete request
              ENDIF
#if defined(T3E)
! DEPENDS ON: um_fort_flush
              call um_fort_flush(8, icode)
#else
              CLOSE(8)
              OPEN(8,FILE=FILENAME)
#endif
            ENDIF

          ENDIF
#if defined(MPP)
        ENDIF
#endif
      ENDIF ! (I_AO == 1 .or. I_AO == 2 .or. I_AO == 4)

 620  FORMAT('%%% ',A14,' ARCHIVE DUMP')

!
!     Normal return
!
      RETURN
!
!     Error returns
!
 900  ICODE=1
      CMESSAGE="DUMPCTL : Fail to open output dump - may already exist"
 999  CONTINUE
      RETURN
!L----------------------------------------------------------------------
      END SUBROUTINE DUMPCTL
#endif
