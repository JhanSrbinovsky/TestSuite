#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL    Subroutine: MEANCTL-----------------------------------------
!LL
!LL    Purpose: To accumulate partial sums and create time-meaned data
!LL
!LL    Tested under compiler:           Tested under OS version:
!LL    cft77                            UNICOS 5.1
!LL
!LL   Programming standard: UM Doc Paper 3 vn2 (7/9/90)
!LL
!LL  Logical system components covered: C5
!LL
!LL    Project tasks: C5,C51,C52
!LL
!LL    External documentation: UMDP C5 - Control of means calculations
!LL
!LLEND ------------------------------------------------------------
!*L    Interface and arguments:
      SUBROUTINE MEANCTL (                                              &
#include "argd1.h"
#include "argduma.h"
#include "argptra.h"
#include "argsts.h"
#include "arglndm.h"
#include "arginfa.h"
#include "argppx.h"
     &           IND_IM,MEANLEV,ICODE,CMESSAGE)
!
      IMPLICIT NONE

!*L Arguments

#include "clookadd.h"
#include "parvars.h"
#include "cmaxsize.h"
#include "csubmodl.h"
#include "cntl_io.h"
#include "cintfa.h"
#include "typsize.h"
#include "typd1.h"
#include "typduma.h"
#include "typptra.h"
#include "typsts.h"
#include "typlndm.h"
#include "typinfa.h"
#include "ppxlook.h"
#include "decomptp.h"

      INTEGER                                                           &
     &       IND_IM,                                                    &
                               ! IN Internal model indicator
     &       MEANLEV,                                                   &
                               ! INOUT  Mean level indicator
     &       ICODE             ! OUT return code; successful=0, error> 0

!
      CHARACTER*80                                                      &
     &       CMESSAGE          ! OUT Error message if ICODE > 0
!
!      Common blocks
!
#include "ctime.h"
#include "chsunits.h"
#include "chistory.h"
#include "ccontrol.h"
#include "cmeanctl.h"
#include "cenvir.h"
!
!      Local variables and arrays
!

      INTEGER                                                           &
     &  HEAD_LEN                                                        &
                                ! Length of each record in header
     &  ,HEAD_SIZE              ! Size of header on disk

      PARAMETER(                                                        &
     &  HEAD_LEN=3                                                      &
     &  )

      INTEGER                                                           &
     &  HEAD_OUT(HEAD_LEN,TOTITEMS) ! Header info for output ps file

      INTEGER                                                           &
     &       IFIND,                                                     &
                                  ! Loop counter
     &       INDEXL,                                                    &
                                  ! Loop index
     &       STEP_DUMPS,                                                &
                                  ! Timestep (in multiples of restart
                                  !           dump frequency)
     &       NMEANS,                                                    &
                                  ! No. of means chosen (fixed, unless
                                  ! there is a mean offset)
     &       RESIDU,                                                    &
                                  ! Reference for partial sum data
     &       MEANS_TOTAL,                                               &
                                  ! Number of mean dumps this timestep
     &       NMVALS(4),                                                 &
                                  ! Absolute meaning periods (in
                                  ! multiples of restart dump frequency)
     &       PS_FLAG(4),                                                &
                                  ! Flag for partial sum updating
     &       FT_PS(4,2),                                                &
                                  ! Unit numbers for dumps of partial
                                  ! sum data (read/write)
     &       PERIODLEN,                                                 &
                                  ! Length in days of current period N
     &       PERIODLENDM,                                               &
                                  ! Length in dumps of period 1 or
                                  ! in days of current period N>1
     &       DUMPS_PER_DAY                                              &
                                  ! Number of restart dumps per day
#if defined(MPP)
     &       ,orig_decomp                                               &
                                  ! Used to check for change in
     &       ,new_decomp                                                &
                                  ! decomposition
#endif
     &       ,TAG                                                       &
                                  ! Stash tag
     &       ,PTD1                                                      &
                                  ! Pointer to D1_ADDR information
     &       ,MODNUM                                                    &
                                  ! Pointer to D1_ADDR submodel
     &       ,PTL                                                       &
                                  ! Pointer to LOOKUP information
     &       ,global_length                                             &
                                  ! Length of global field
     &       ,CITEMS              ! Counter for no of objects to mean

      INTEGER                                                           &
     &       INDEX_READ,                                                &
                                  ! Specific index no for reading
     &       INDEX_WRITE          ! Specific index no for writing
      INTEGER                                                           &
     &       FT_READ,                                                   &
                                  ! Unit number for partial sum read
     &       FT_WRITE,                                                  &
                                  ! Unit number for partial sum write
     &       FT_DELETE,                                                 &
                                  ! Unit number for partial sum delete
     &       FT_SSD                                                     &
                                  ! Unit no for transfer of inst. data
     &,LEN_PSNAME                                                       &
     &,disk_address                                                     &
                                       ! Current rounded disk address
     &,number_of_data_words_on_disk                                     &
                                       ! Number of data words on disk
     &,number_of_data_words_in_memory                                   &
                                       ! Number of Data Words in memory
     &,maximum_file_length             ! Maximum file length for FT_SSD
      INTEGER internal_model
      INTEGER im      ! temporary internal model id for ocean or slab
      INTEGER                                                           &
     &       REINIT_STEPS                                               &
                                  ! dummy input for GET_NAME
     &,      D1_ADDR_SUBMODEL_ID  ! Submodel number in D1_ADDR
      LOGICAL                                                           &
     &       LMEANINC             ! increment MEANS_TOTAL or not

#include "stparam.h"

      INTEGER                                                           &
     &  IE                                                              &
            ! loop counter over items
     &, maxsize                                                         &
                ! maximum dump output length
     &  ,lmaxsize                                                       &
                                ! Maximum output length per level
     &  ,totsize                                                        &
                                ! Size required for partial sum dumps
     &  ,levsize                ! Size of each level
      INTEGER dump_freq
      INTEGER mean_freq1,mean_freq2,mean_freq3,mean_freq4
      INTEGER STASHoutputTime,MEAN1FREQts

!
!      Character data
!
      CHARACTER*14 PSNAME_READ,PSNAME_WRITE,PSNAME_DELETE
      CHARACTER*1                                                       &
     &       BLANK                                                      &
     &,      LETTER_3             ! dummy input for GET_NAME
      CHARACTER*80 FILENAME       ! Name of pipe to server
!
      BLANK=' '
      CMESSAGE=' '
      LETTER_3='a'
      REINIT_STEPS=0

! Get name of pipe
      CALL GET_FILE(8,FILENAME,80,ICODE)
!
!      Set up unit number for instantaneous
!      dump transfer
!
      FT_SSD=17
!
!      Define mode of use for means program
!
      IF(IND_IM == 1)THEN
        WRITE(6,*)'MEANCTL: ***** Called in ATMOSPHERIC mode *****'
      ENDIF
!L
!L----------------------------------------------------------------------
!L     Find out which mean datasets need to be created
!L     on this timestep (if any) and set MEANS_TOTAL accordingly
!L----------------------------------------------------------------------
!L
!
      INDEXL=RUN_MEANCTL_RESTART         ! Zero in normal circumstances
      ICODE=0
      MEANS_TOTAL=0
!
!      Initially check validity of call to subroutine
!
        IF(DUMPFREQim(IND_IM) == 0)THEN ! Is mean dump production off?
          ICODE=1
          CMESSAGE='MEANCTL: Invalid call to subroutine'
      WRITE(6,*) 'MEANCTL: DUMPFREQ(',IND_IM,')= ',DUMPFREQim(IND_IM)
          GO TO 9999
        ELSEIF(MOD(STEPim(IND_IM),DUMPFREQim(IND_IM)) /= 0)THEN
          ICODE=2        ! This is not a dumping timestep
          CMESSAGE='MEANCTL: Incorrect timestep to call subroutine'
          WRITE(6,*) 'MEANCTL: STEP is not a multiple of DUMPFREQ'
          WRITE(6,*) '         STEP(',IND_IM,')= ',STEPim(IND_IM),      &
     &                   ' DUMPFREQ(',IND_IM,')= ', DUMPFREQim(IND_IM)
          GO TO 9999
        ELSE
          STEP_DUMPS = (STEPim(IND_IM)/DUMPFREQim(IND_IM))+             &
     &                   OFFSET_DUMPSim(IND_IM)
        END IF
!
!      Pick up number of means chosen from history file (MEAN_NUMBERim)
!      or number determined by the offset from the reference time whilst
!      the staggered start of means production unwinds  (MEAN_OFFSETim)
!
        DO IFIND=1,MEAN_NUMBERim(IND_IM)
          IF(IFIND == 1)NMVALS(IFIND) = MEANFREQim(IFIND,IND_IM)
          IF(IFIND >  1)NMVALS(IFIND) = MEANFREQim(IFIND,IND_IM)*       &
     &                                   NMVALS(IFIND-1)
        END DO
!
      if (lclimrealyr) then
        DUMPS_PER_DAY=(24*3600*STEPS_PER_PERIODim(IND_IM))/             &
     &                (DUMPFREQim(IND_IM)*SECS_PER_PERIODim(IND_IM))
      end if
!
        IF(MEAN_OFFSETim(IND_IM) == MEAN_NUMBERim(IND_IM))THEN
          NMEANS=MEAN_NUMBERim(IND_IM)
        ELSE  ! there is a non-zero offset so increment mean_offset
              ! when run is partway into the current mean period.
          DO IFIND=MEAN_OFFSETim(IND_IM)+1,MEAN_NUMBERim(IND_IM)
            if (lclimrealyr) then
              if ((ifind  ==  1 .and. i_day  ==  2 .and.                &
     &           (dumps_per_day  <=  1)) .or.                           &
                                                        ! 24h dumps
     &           (ifind  ==  1 .and. i_day  ==  1 .and.                 &
                                                        ! dumps lt 24h
     &           i_hour  ==  (24/dumps_per_day))) then
                mean_offsetim(IND_IM)=mean_offsetim(IND_IM)+1 ! months
              elseif (ifind  ==  2 .and.                                &
     &               (i_day  ==  1) .and. i_hour  ==  0 .and.           &
     &               mod(((i_month-1)-MEAN_REFTIMEim(2,IND_IM)),3) == 0)&
     &               then
                mean_offsetim(IND_IM)=mean_offsetim(IND_IM)+1 ! seasons
              elseif (ifind  ==  3 .and.                                &
     &               (i_day  ==  1) .and. (i_hour  ==  0) .and.         &
     &              mod(((i_month-3)-MEAN_REFTIMEim(2,IND_IM)),12) == 0)&
     &               then
                mean_offsetim(IND_IM)=mean_offsetim(IND_IM)+1 ! years
              endif  ! of test on ifind, i_day, etc.
            else  ! for 360d year, use array of no. of mean dumps
              IF(STEP_DUMPS <  0)THEN  ! mean ref time not reached yet
              RESIDU=1-NMVALS(IFIND)
              ELSE                     ! mean ref time has been passed
              RESIDU=1
            ENDIF
            IF(MOD(STEP_DUMPS,NMVALS(IFIND)) == RESIDU)THEN
              MEAN_OFFSETim(IND_IM)=MEAN_OFFSETim(IND_IM)+1
            ENDIF
            endif  ! end of checking whether to increment MEAN_OFFSETim
          ENDDO  ! end of loop over IFIND
          NMEANS=MEAN_OFFSETim(IND_IM)
          WRITE(6,*)' MEAN_OFFSET(',IND_IM,')=',MEAN_OFFSETim(IND_IM)
        ENDIF  ! end of setting NMEANS
!L
!L      If no processing is required (because of staggered
!L      start in means production) then skip to end of subroutine
!L
      IF(NMEANS == 0)THEN
        ICODE=-1
        CMESSAGE='MEANCTL: No accumulation/meaning done this step'
        WRITE(6,*)'MEANCTL: No accm/meaning due to staggered start'
        GO TO 9999
      ENDIF
!L
!L      Output message whilst stagger unwinds
!L
        DO IFIND=1,MEAN_NUMBERim(IND_IM)
          IF(NMEANS <  IFIND)THEN
            WRITE(6,*)'MEANCTL: Period_',IFIND,' mean not activated',   &
     &      BLANK,'because of staggered start in means production'
          ENDIF
        ENDDO
!
      DO IFIND=1,NMEANS
!
        PS_FLAG(IFIND)=0
        lmeaninc=.false. ! Initialise to avoid false positive results
!
! Find out if the end of any meaning period has been reached. If so,
! increment MEANS_TOTAL
        if (lclimrealyr .and. (ifind  ==  1) .and.                      &
     &     (i_day  ==  1) .and. (i_hour == 0)                           &
                                                ! but not 1st day of run
     &     .and. (STEPim(IND_IM)  >   STEPS_PER_PERIODim(IND_IM))) then
          lmeaninc=.true.                               ! monthly
        elseif (lclimrealyr .and. (ifind  ==  2) .and.                  &
     &         (i_day  ==  1) .and. (i_hour  ==  0) .and.               &
                                                          ! seasonal
     &         (STEPim(IND_IM)  >   STEPS_PER_PERIODim(IND_IM)) .and.   &
     &         mod((i_month-MEAN_REFTIMEim(2,IND_IM)),3) == 0) then
          lmeaninc=.true.
        elseif (lclimrealyr .and. (ifind  ==  3) .and.                  &
     &         (i_day  ==  1) .and. (i_hour  ==  0) .and.               &
                                                          ! annual
     &         (STEPim(IND_IM)  >   STEPS_PER_PERIODim(IND_IM)) .and.   &
     &         mod((i_month-MEAN_REFTIMEim(2,IND_IM)),12) == 0) then
          lmeaninc=.true.
        elseif (.not. lclimrealyr .and.                                 &
                                            ! 360-day calendar
     &         MOD(STEP_DUMPS,NMVALS(IFIND)) == 0) then
          lmeaninc=.true.
        endif  ! end of IF test on reaching end of meaning period
        if (lmeaninc) then
          MEANS_TOTAL=MEANS_TOTAL+1
          IF(RUN_MEANCTL_RESTART >  1.AND.                              &
     &       RUN_MEANCTL_RESTART >  IFIND)THEN
            WRITE(6,*) 'MEANCTL: Period_',IFIND,' mean dump',BLANK,     &
     &                        'already created this timestep'
          ELSE
            WRITE(6,*) 'MEANCTL: Period_',IFIND,' mean dump',BLANK,     &
     &                        'to be created this timestep'
          ENDIF
        endif  ! end of IF test on lmeaninc
!
! Find out if run is one period(N-1) into the period(N), in which case
! set PS_FLAG(IFIND)=1 because there will not already be a file for
! that partial sum on the disk, so ACUMPS must get it from D1.
!
        if (lclimrealyr) then
          if ((ifind  ==  1 .and. i_day  ==  2 .and.                    &
     &       (dumps_per_day  <=  1))                                    &
                                                         ! 24h dumps
     &       .or. (ifind  ==  1 .and. i_day  ==  1 .and.                &
                                                         ! dumps lt 24h
     &       i_hour  ==  (24/dumps_per_day))) then
            ps_flag(1)=1                                   ! months
          elseif (ifind  ==  2 .and.                                    &
     &           (i_day  ==  1) .and. (i_hour  ==  0) .and.             &
     &           mod(((i_month-1)-MEAN_REFTIMEim(2,IND_IM)),3) == 0)    &
     &           then
            ps_flag(2)=1                                   ! seasons
          elseif (ifind  ==  3 .and.                                    &
     &           (i_day  ==  1) .and. (i_hour  ==  0) .and.             &
     &           mod(((i_month-3)-MEAN_REFTIMEim(2,IND_IM)),12) == 0)   &
     &           then
            ps_flag(3)=1                                   ! years
          endif
        else  ! for 360d year, check using array of no. of mean dumps
        IF(STEP_DUMPS <  0)THEN
          IF(IFIND == 1)RESIDU=1-NMVALS(IFIND)
          IF(IFIND /= 1)RESIDU=NMVALS(IFIND-1)-NMVALS(IFIND)
        ELSE
          IF(IFIND == 1)RESIDU=1
          IF(IFIND /= 1)RESIDU=NMVALS(IFIND-1)
        ENDIF
        IF(MOD(STEP_DUMPS,NMVALS(IFIND)) == RESIDU)THEN
          PS_FLAG(IFIND)=1
        ENDIF
        endif  ! end of test on lclimrealyr and setting of PS_FLAG
!
!
      ENDDO       ! end of loop over IFIND from 1 to NMEANS
!
!      Set up unit numbers for partial sum dumps
!
      IF (IND_IM == 1) THEN
        FT_PS(1,1)=23
        FT_PS(1,2)=24
      ENDIF
      DO IFIND=2,4
        FT_PS(IFIND,1)=25
        FT_PS(IFIND,2)=26
      ENDDO
! Units must still be alternated as 3 may be open at a time
      FT_PS(3,1)=45
      FT_PS(3,2)=46

! check STASH items output time against restart dump period
      dump_freq=DUMPFREQim(ind_im)     ! Unit: timestep
      mean_freq1=MEANFREQim(1,ind_im)
      mean_freq2=MEANFREQim(2,ind_im)
      mean_freq3=MEANFREQim(3,ind_im)
      mean_freq4=MEANFREQim(4,ind_im)
! Calculating mean_freq1 in unit of timestep
      MEAN1FREQts=mean_freq1*dump_freq
! Create header for partial sum file
#if defined(ATMOS)
      IF (IND_IM == 1) THEN
! Set pointer to submodel info in D1_ADDR
        MODNUM=SUBMODEL_FOR_SM(IND_IM)
        CITEMS=0
        DO IE=1,TOTITEMS
          TAG=STLIST(st_macrotag,IE)/1000
! Get pointer to element in D1_ADDR array
          PTD1=STLIST(st_d1pos,IE)
          IF(tag /= 0.AND.                                              &
     &      STLIST(s_modl,IE) == D1_ADDR(d1_imodl,PTD1,MODNUM))THEN
          STASHoutputTime=STLIST(st_freq_code,ie)  ! Unit: timestep
          IF (STASHoutputTime >= MEAN1FREQts) THEN
            WRITE(6,*)'ERROR: STASH output time frequency greater'
            WRITE(6,*)'       than climate mean period frequency!'
            WRITE(6,*)'STASH output frequency=',STASHoutputTime,'(ts)'  &
     &               ,'  Period 1 period=',MEAN1FREQts,'(ts)'
            WRITE(6,*)'Please check using VERIFY DIAGNOSTICS in UMUI'
            ICODE=4
            CMESSAGE='MEANCTL: ERROR - STASHoutput > CMEAN period'
            GO TO 9999
          END IF
! Tagged for meaning and submodel information matches.
! Counter for no of variables that will be processed
            CITEMS=CITEMS+1
            PTL=D1_ADDR(d1_lookup_ptr,PTD1,MODNUM)
#if defined(MPP)
            global_length=STLIST(st_dump_level_output_length,IE)
#else
            global_length=STLIST(st_output_length,IE)/                  &
     &        D1_ADDR(d1_no_levels,PTD1,MODNUM)
#endif
            global_length=global_length+1
            HEAD_OUT(1,CITEMS)=GLOBAL_LENGTH
            HEAD_OUT(2,CITEMS)=1
            HEAD_OUT(3,CITEMS)=(((GLOBAL_LENGTH+um_sector_size)         &
     &        /um_sector_size)*um_sector_size)
            IF (MOD((A_LOOKUP(LBPACK,PTL)),10)  ==  2) THEN
              IF (A_LOOKUP(DATA_TYPE,PTL)  ==  1) THEN
! Data to be packed so reduce data length accordingly
                GLOBAL_LENGTH=(GLOBAL_LENGTH+1)/2
                HEAD_OUT(1,CITEMS)=GLOBAL_LENGTH
                HEAD_OUT(2,CITEMS)=2
                HEAD_OUT(3,CITEMS)=(((GLOBAL_LENGTH+um_sector_size)     &
     &            /um_sector_size)*um_sector_size)
              END IF
            END IF
          END IF
        END DO
      END IF
#endif

! Find the largest output field size to dimension I/O buffer array
! and the total size of storage required to set file size
        lmaxsize=1
        head_size=CITEMS*HEAD_LEN+UM_SECTOR_SIZE+2
        head_size=(head_size/UM_SECTOR_SIZE)*UM_SECTOR_SIZE
        totsize=head_size
        CITEMS=0

        DO IE=1,TOTITEMS
          tag=STLIST(st_macrotag,IE)/1000
! Get pointer to element in D1_ADDR array
          PTD1=STLIST(st_d1pos,IE)
          IF (tag /= 0.AND.                                             &
     &      STLIST(s_modl,IE) == D1_ADDR(d1_imodl,PTD1,MODNUM)) THEN
! Tagged for meaning and submodel information matches.
            CITEMS=CITEMS+1
#if defined(MPP)
            lmaxsize=max(lmaxsize,                                      &
     &        STLIST(st_dump_level_output_length,IE))
#else
            lmaxsize=max(lmaxsize,STLIST(st_output_length,IE)/          &
     &        D1_ADDR(d1_no_levels,PTD1,MODNUM))
#endif
            totsize=totsize+                                            &
     &        HEAD_OUT(3,CITEMS)*D1_ADDR(d1_no_levels,PTD1,MODNUM)
          END IF
        END DO
! Adding 1 extra space allowing for checksum
        lmaxsize=lmaxsize+1

        lmaxsize=((lmaxsize+UM_SECTOR_SIZE)/UM_SECTOR_SIZE)             &
     &    *UM_SECTOR_SIZE

!L
!L**********************************************************************
!L                  LOGICAL SUB-PROCESS C51
!L     Start of default process: updating period_1 partial sum data
!L**********************************************************************
!L
!
      IF (RUN_MEANCTL_RESTART == 0) THEN
!
        INDEX_READ=RUN_MEANCTL_INDICim(1,IND_IM)
        INDEX_WRITE=3-RUN_MEANCTL_INDICim(1,IND_IM)
        FT_READ=FT_PS(1,INDEX_READ)
        FT_WRITE=FT_PS(1,INDEX_WRITE)
!
!      Temporary check on unit numbers
!
        IF (PS_FLAG(1) /= 1) THEN
          WRITE(6,*) 'Period_1 data read from unit number ',FT_READ
        END IF
        WRITE(6,*) 'Period_1 data written to unit number ',FT_WRITE
!
!L
!L                           STEP 1
!L     Update or create period_1 partial sum data and write out
!L     to period_1 partial sum dump
!L
#if defined(MPP)

        orig_decomp=current_decomp_type
        new_decomp=orig_decomp

        IF (IND_IM==1 .AND. orig_decomp/=decomp_standard_atmos) THEN
          new_decomp = decomp_standard_atmos
        ELSE IF (IND_IM==2 .AND. orig_decomp/=decomp_nowrap_ocean) THEN
          new_decomp = decomp_nowrap_ocean
        END IF
        IF (new_decomp  /=  orig_decomp) THEN
! DEPENDS ON: change_decomposition
          CALL CHANGE_DECOMPOSITION(new_decomp,icode)
          IF (ICODE /= 0) THEN
            WRITE(6,*)'ERROR : MEANCTL'
            WRITE(6,*)'Failed to change decomposition to ',new_decomp
            CMESSAGE='MEANCTL1: Failed to change decomposition'
            GO TO 9999
          END IF
        END IF
#endif
        IF (MEANS_TOTAL >  0) THEN
!L
!L                 STEP 1
!L     Copy instantaneous dump to SSD
!L
!L Strictly, only climate mean tagged diagnostics need to be saved.

! Compute the maximum length for the FT_SSD File
          maximum_file_length=0
#if defined(ATMOS)
          maximum_file_length=max(maximum_file_length,                  &
     &      a_len_data+(model_levels+1)*THETA_FIELD_SIZE)
#endif
! Set the length of the file needed
          call set_dumpfile_length(ft_ssd, maximum_file_length)
!
#if defined(ATMOS)
          IF (IND_IM == 1) THEN
! DEPENDS ON: transout
            CALL TRANSOUT(                                              &
#include "argd1.h"
     &        A_LEN_DATA+(MODEL_LEVELS+1)*THETA_FIELD_SIZE,FT_SSD,IND_IM&
     &        ,ICODE,CMESSAGE)
          END IF

#endif

!      Check return code from TRANSOUT
          IF (ICODE /= 0) THEN
            RUN_MEANCTL_RESTART=1
            WRITE(6,*) 'MEANCTL: RESTART AT PERIOD_',RUN_MEANCTL_RESTART
            GO TO 9999
          END IF

! DEPENDS ON: file_close
          CALL FILE_CLOSE(FT_SSD,FT_ENVIRON(FT_SSD),                    &
     &      LEN_FT_ENVIR(FT_SSD),0,0,ICODE)

        END IF
!
!
!--preset the file lengths prior to the open
!
!--pass the new file length to the I/O routines
        call set_dumpfile_length(ft_read , totsize)
        call set_dumpfile_length(ft_write, totsize)
!
!      Open input and output partial sum files (preassigned names)
!
! DEPENDS ON: file_open
        CALL FILE_OPEN(FT_READ,FT_ENVIRON(FT_READ),                     &
     &                 LEN_FT_ENVIR(FT_READ),1,0,ICODE)
        IF (ICODE /= 0) GO TO 9999
! DEPENDS ON: file_open
        CALL FILE_OPEN(FT_WRITE,FT_ENVIRON(FT_WRITE),                   &
     &                 LEN_FT_ENVIR(FT_WRITE),1,0,ICODE)
        IF (ICODE /= 0) GO TO 9999
!--zero the file lengths
        call set_dumpfile_length(ft_read , 0)
        call set_dumpfile_length(ft_write, 0)

#if defined(ATMOS)
        IF (IND_IM == 1) THEN
          D1_ADDR_SUBMODEL_ID = SUBMODEL_FOR_SM(1)
! DEPENDS ON: acumps
          CALL ACUMPS(                                                  &
     &  NO_OBJ_D1(D1_ADDR_SUBMODEL_ID),                                 &
     &  D1_ADDR(1,1,D1_ADDR_SUBMODEL_ID),                               &
     &  A_LEN_DATA,D1,D1,D1,                                            &
     &  LMAXSIZE,MEANS_TOTAL,                                           &
     &  PS_FLAG(1),FT_READ,FT_WRITE,LCLIMREALYR,MEANLEV,                &
     &  I_MONTH,I_YEAR,                                                 &
     &  HEAD_OUT,HEAD_LEN,HEAD_SIZE,                                    &
     &  STEPim(IND_IM),CITEMS,A_FIXHD(12),                              &
#include "argsts.h"
     &  ICODE,CMESSAGE)
        END IF
#endif
!
!      Check return code from ACUMPS
!
        IF (ICODE /= 0) THEN
          WRITE(6,*) 'MEANCTL: RESTART AT PERIOD_',RUN_MEANCTL_RESTART
          GO TO 9999
        END IF
!
!      Close input and output partial sum files
!
! DEPENDS ON: file_close
        CALL FILE_CLOSE(FT_READ,FT_ENVIRON(FT_READ),                    &
     &                  LEN_FT_ENVIR(FT_READ),0,1,ICODE)
! DEPENDS ON: file_close
        CALL FILE_CLOSE(FT_WRITE,FT_ENVIRON(FT_WRITE),                  &
     &                  LEN_FT_ENVIR(FT_WRITE),0,0,ICODE)
!
        INDEXL=INDEXL+1
!
!      Update RUN_MEANCTL_INDICim for period_1 data
!
        RUN_MEANCTL_INDICim(1,IND_IM)=3-RUN_MEANCTL_INDICim(1,IND_IM)
!
      END IF
!L
!L**********************************************************************
!L     End of default process: updating period_1 partial sum data
!L**********************************************************************
!L
      IF (MEANS_TOTAL >  0) THEN
!L
!L**********************************************************************
!L                  LOGICAL SUB-PROCESS C52
!L     Start of means processing and updating of subsequent
!L     partial sum dump
!L**********************************************************************
!
        DO MEANLEV=INDEXL,MEANS_TOTAL
!
          INDEX_READ=RUN_MEANCTL_INDICim(MEANLEV,IND_IM)
          FT_READ=FT_PS(MEANLEV,INDEX_READ)
!
!      Temporary check on unit number
!
        WRITE(6,*) 'Period_',MEANLEV,' data read:unit number ',FT_READ
!L
!L                             STEP 2
!L     Generate period_N time-meaned data and store in main data block
!L
!       If real-period meaning selected, find length of current period
        IF (LCLIMREALYR) THEN
! DEPENDS ON: setperlen
          CALL SETPERLEN (MEANLEV,I_MONTH,I_YEAR,PERIODLEN)
          if (meanlev == 1) then ! divisor only needs to be in terms
                                 ! of restart dumps for Period_1
            PERIODLENDM=PERIODLEN*DUMPS_PER_DAY
          else
            PERIODLENDM=PERIODLEN
          end if
        END IF
#if defined(ATMOS)
          IF (IND_IM == 1) THEN
!
!      Open input partial sum file (preassigned or calculated name)
!
!--pass the new file length to the I/O routines
            call set_dumpfile_length(ft_read , totsize)
            IF (MEANLEV == 1) THEN
! DEPENDS ON: file_open
        CALL FILE_OPEN(FT_READ,FT_ENVIRON(FT_READ),                     &
     &            LEN_FT_ENVIR(FT_READ),1,0,ICODE)
        IF (ICODE /= 0) GO TO 9999
            ELSE
! DEPENDS ON: get_name
              CALL GET_NAME(EXPT_ID,JOB_ID,IND_IM,MEANLEV,INDEX_READ,   &
     &                      REINIT_STEPS,'s',LETTER_3,                  &
     &                      MODEL_STATUS,TIME_CONVENTION,               &
     &                       0,PSNAME_READ,ICODE,CMESSAGE,LCAL360)
              IF (ICODE >  0) GO TO 9999
              LEN_PSNAME=LEN(PSNAME_READ)
! DEPENDS ON: file_open
              CALL FILE_OPEN(FT_READ,PSNAME_READ,LEN_PSNAME,1,1,ICODE)
              IF (ICODE /= 0) GO TO 9999
            ENDIF
!--unset the file length in the I/O routines
            call set_dumpfile_length(ft_read , 0)
!
            IF (LCLIMREALYR) THEN  ! Real-period meaning selected
! DEPENDS ON: meanps
         CALL MEANPS(                                                   &
     &          NO_OBJ_D1(D1_ADDR_SUBMODEL_ID),                         &
     &          D1_ADDR(1,1,D1_ADDR_SUBMODEL_ID),                       &
     &          A_LEN_DATA,D1,D1,D1,                                    &
#include "argsts.h"
     &          PERIODLENDM                                             &
     &          )
            ELSE            ! 360d year meaning selected
! DEPENDS ON: meanps
              CALL MEANPS(                                              &
     &          NO_OBJ_D1(D1_ADDR_SUBMODEL_ID),                         &
     &          D1_ADDR(1,1,D1_ADDR_SUBMODEL_ID),                       &
     &          A_LEN_DATA,D1,D1,D1,                                    &
#include "argsts.h"
     &          MEANFREQim(MEANLEV,IND_IM)                              &
     &          )

            END IF            ! end of test on LCLIMREALYR
          END IF            ! end of test on IND_IM == 1
#endif
!
!      Check return code from MEANPS
!
          IF (ICODE /= 0) THEN
            RUN_MEANCTL_RESTART=MEANLEV
            WRITE(6,*) 'MEANCTL: RESTART AT PERIOD_',RUN_MEANCTL_RESTART
            GO TO 9999
          ENDIF
!
        PSNAME_DELETE=PSNAME_READ
! DEPENDS ON: setpos
          CALL SETPOS(FT_READ,0,ICODE)
          FT_DELETE=FT_READ
!L
!L                          STEP 3.1
!L     Calculate mean diagnostics and extract PPfields from mean data
!L
#if defined(ATMOS)
          IF (IND_IM == 1) THEN
! [p_exner calculated here from SETEXNER pre-vn5.0. Not now required.]

#if defined(MPP)
! Find the largest output field size to dimension I/O buffer array
            maxsize=1

            DO IE=1,TOTITEMS
              tag=STLIST(st_macrotag,IE)/1000
              IF (MOD(tag/(2**(MEANLEV-1)),2)  ==  1) THEN
                maxsize=MAX(maxsize,                                    &
     &                      STLIST(st_dump_level_output_length,IE))
              END IF
            END DO
#endif
! Extract mean diagnostics (both normal and sections 21-24/41-44)
! DEPENDS ON: meandiag
            CALL MEANDIAG (                                             &
#include "argd1.h"
#include "argduma.h"
#include "argsts.h"
#include "argptra.h"
#include "arglndm.h"
#include "arginfa.h"
#include "argppx.h"
     &       IND_IM,MEANLEV,PP_LEN2_MEANim(MEANLEV,IND_IM),STEP_DUMPS,  &
     &           NMVALS(MEANLEV),                                       &
#if defined(MPP)
     &           maxsize,                                               &
#endif
     &           ICODE,CMESSAGE)
          END IF
#endif
          IF (ICODE >  0) GO TO 9999
!L
!L                          STEP 4
!L     Check to see if period_N+1 partial sum data needs to updated
!L     or created.
!L     If so, proceed and write out to period_N+1 partial sum dump
!L
          IF (MEANLEV /= NMEANS) THEN
!
            INDEX_READ=RUN_MEANCTL_INDICim(MEANLEV+1,IND_IM)
            INDEX_WRITE=3-RUN_MEANCTL_INDICim(MEANLEV+1,IND_IM)
            FT_READ=FT_PS(MEANLEV+1,INDEX_READ)
            FT_WRITE=FT_PS(MEANLEV+1,INDEX_WRITE)
!
!      Temporary check on unit numbers
!
        IF (PS_FLAG(MEANLEV+1) /= 1) THEN
      WRITE(6,*) 'Period_',MEANLEV+1,' data read:unit number ',FT_READ
        END IF
        WRITE(6,*) 'Period_',MEANLEV+1,' data written:unit number ',    &
     &          FT_WRITE
!
#if defined(ATMOS)
            IF (IND_IM == 1) THEN
!
!      Open input and output partial sum files (calculated names)
!
! DEPENDS ON: get_name
              CALL GET_NAME(EXPT_ID,JOB_ID,IND_IM,MEANLEV+1,INDEX_READ, &
     &                      REINIT_STEPS,'s',LETTER_3,                  &
     &                      MODEL_STATUS,TIME_CONVENTION,               &
     &                       0,PSNAME_READ,ICODE,CMESSAGE,LCAL360)
              IF (ICODE >  0) GO TO 9999
!--pass the new file length to the I/O routines
              call set_dumpfile_length(ft_read , totsize)
              LEN_PSNAME=LEN(PSNAME_READ)
! DEPENDS ON: file_open
              CALL FILE_OPEN(FT_READ,PSNAME_READ,LEN_PSNAME,1,1,ICODE)
              IF (ICODE /= 0) GO TO 9999
!--unset the file length in the I/O routines
              call set_dumpfile_length(ft_read , 0)
!
! DEPENDS ON: get_name
              CALL GET_NAME(EXPT_ID,JOB_ID,IND_IM,MEANLEV+1,INDEX_WRITE,&
     &                      REINIT_STEPS,'s',LETTER_3,                  &
     &                      MODEL_STATUS,TIME_CONVENTION,               &
     &                       0,PSNAME_WRITE,ICODE,CMESSAGE,LCAL360)
              IF (ICODE >  0) GO TO 9999
!
!--pass the new file length to the I/O routines
              call set_dumpfile_length(ft_write, totsize)
              LEN_PSNAME=LEN(PSNAME_WRITE)
! DEPENDS ON: file_open
              CALL FILE_OPEN(FT_WRITE,PSNAME_WRITE,LEN_PSNAME,1,1,ICODE)
              IF (ICODE /= 0) GO TO 9999
!--unset the file length in the I/O routines
              call set_dumpfile_length(ft_write, 0)
!
      D1_ADDR_SUBMODEL_ID = SUBMODEL_FOR_SM(1)
! DEPENDS ON: acumps
          CALL ACUMPS(                                                  &
     &  NO_OBJ_D1(D1_ADDR_SUBMODEL_ID),                                 &
     &  D1_ADDR(1,1,D1_ADDR_SUBMODEL_ID),                               &
     &  A_LEN_DATA,D1,D1,D1,                                            &
     &  LMAXSIZE,MEANS_TOTAL,                                           &
     &  PS_FLAG(MEANLEV+1),FT_READ,FT_WRITE,                            &
     &  LCLIMREALYR,MEANLEV,I_MONTH,I_YEAR,                             &
     &  HEAD_OUT,HEAD_LEN,HEAD_SIZE,                                    &
     &  STEPim(IND_IM),CITEMS,A_FIXHD(12),                              &
#include "argsts.h"
     &  ICODE,CMESSAGE)
          END IF
#endif
!
!      Check return code from ACUMPS
!
        IF (ICODE /= 0) THEN
        WRITE(6,*) 'MEANCTL: RESTART AT PERIOD_',RUN_MEANCTL_RESTART
          GO TO 9999
        END IF

!
!      Update RUN_MEANCTL_INDICim for period_N+1 data
!
            RUN_MEANCTL_INDICim(MEANLEV+1,IND_IM)=                      &
     &        3-RUN_MEANCTL_INDICim(MEANLEV+1,IND_IM)
!
          END IF
!L
!
!      Decide disposition of period_N+1 partial sum dumps
!      NB: for restartability it is NOT safe to delete period 2+ sums
!
          IF (MEANLEV /= NMEANS) THEN
            LEN_PSNAME=LEN(PSNAME_READ)
! DEPENDS ON: setpos
            CALL SETPOS(FT_READ,0,ICODE)
! DEPENDS ON: file_close
            CALL FILE_CLOSE(FT_READ,PSNAME_READ,LEN_PSNAME,1,0,ICODE)
            LEN_PSNAME=LEN(PSNAME_WRITE)
! DEPENDS ON: setpos
            CALL SETPOS(FT_WRITE,0,ICODE)
! DEPENDS ON: file_close
            CALL FILE_CLOSE(FT_WRITE,PSNAME_WRITE,LEN_PSNAME,1,0,ICODE)
          END IF
!
!      Decide disposition of remaining period_N partial sum dump
!      NB: for restartability it is NOT safe to delete period 2 sums
!
          IF (MEANLEV >= 2) THEN
            LEN_PSNAME=LEN(PSNAME_DELETE)
! DEPENDS ON: file_close
            CALL FILE_CLOSE(FT_DELETE,PSNAME_DELETE,LEN_PSNAME,1,0,     &
     &      ICODE)
          ELSE
! DEPENDS ON: file_close
            CALL FILE_CLOSE(FT_DELETE,FT_ENVIRON(FT_DELETE),            &
     &      LEN_FT_ENVIR(FT_DELETE),0,0,ICODE)
          END IF
!
        END DO
!L
!L                    STEP 6
!L     Read back instantaneous dump from SSD
!L
#if defined(ATMOS)
        IF (IND_IM == 1) THEN
! DEPENDS ON: transin
          CALL TRANSIN(                                                 &
#include "argd1.h"
     &                  A_LEN_DATA+(MODEL_LEVELS+1)*THETA_FIELD_SIZE,   &
     &                  FT_SSD,IND_IM                                   &
     &    ,ICODE,CMESSAGE)
          A_FIXHD(5)=1
        END IF
#endif
!
!      Check return code from TRANSIN
!
        IF (ICODE /= 0) THEN
          RUN_MEANCTL_RESTART=0
          WRITE(6,*) 'MEANCTL: MEANS COMPLETE',BLANK,                   &
     &             '- RECOVERY OF INSTANTANEOUS DATA HAS FAILED'
          GO TO 9999
        END IF
!
!       CALL SETPOS(FT_SSD,0,ICODE)
! DEPENDS ON: file_close
        CALL FILE_CLOSE(FT_SSD,FT_ENVIRON(FT_SSD),LEN_FT_ENVIR(FT_SSD), &
     &  0,0,ICODE)
!
        IF (ICODE  ==  0) THEN
! Meaning has been successful so it is now safe to delete the restart
! dumps from the previous dump time.

      internal_model=atmos_im

#if defined(MPP)
          IF (mype == 0) THEN
#endif

            IF (IND_IM  ==  A_IM) THEN

              IF ((internal_model  ==  atmos_im).OR.                    &
     &             (internal_model  ==  slab_im))  THEN
                IF (LASTDMPim(A_IM) /= "              ") THEN
                  WRITE(8,890) LASTDMPim(A_IM)
                  CLOSE(8)
                END IF
                OPEN(8,FILE=FILENAME)
              END IF


            ELSE IF (IND_IM  ==  O_IM) THEN

              IF (LASTDMPim(O_IM) /= "              ") THEN
                WRITE(8,890) LASTDMPim(O_IM)
                CLOSE(8)
              END IF
              OPEN(8,FILE=FILENAME)

              IF (internal_model  ==  atmos_im)  THEN


!               There is > 1 internal model ie. coupled then delete the
!               last atmos dump too (should put stronger test here)
                IF (LASTDMPim(A_IM) /= "              ") THEN
                  WRITE(8,890) LASTDMPim(A_IM)
                  CLOSE(8)
                END IF
                OPEN(8,FILE=FILENAME)
              END IF


            ELSE IF (IND_IM  ==  W_IM) THEN

              IF (LASTDMPim(W_IM) /= "              ") THEN
                WRITE(8,890) LASTDMPim(W_IM)
                CLOSE(8)
              END IF
              OPEN(8,FILE=FILENAME)

            END IF

#if defined(MPP)
          END IF ! (mype == 0)
#endif

 890      FORMAT('%%% ',A14,' DELETE')

        END IF
      END IF
!L
!L**********************************************************************
!L     End of means processing and updating of subsequent
!L     partial sum dumps
!L**********************************************************************
!L
!      Reset RUN_MEANCTL_RESTART to zero
!
      RUN_MEANCTL_RESTART=0
!

#if defined(MPP)
      orig_decomp=current_decomp_type
      new_decomp=orig_decomp
      IF (IND_IM == 1.AND.orig_decomp /= decomp_standard_atmos)THEN
        new_decomp=decomp_standard_atmos
      ELSE IF (IND_IM==2.AND.orig_decomp/=decomp_standard_ocean) THEN
        new_decomp=decomp_standard_ocean
      END IF
      IF (new_decomp  /=  orig_decomp) THEN
! DEPENDS ON: change_decomposition
        CALL CHANGE_DECOMPOSITION(new_decomp,icode)
        IF (ICODE /= 0) THEN
          WRITE(6,*)'ERROR : MEANCTL'
          WRITE(6,*)'Failed to change decomposition to ',new_decomp
          CMESSAGE='MEANCTL1: Failed to change decomposition'
        END IF
      END IF
#endif

 9999 CONTINUE
!
!      Reset MEANLEV to zero
!
      MEANLEV=0
      RETURN
      END SUBROUTINE MEANCTL
#endif
