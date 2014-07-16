#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE NUM_OBS   NUM_OBS2   CHECK_OBS  ------------------------
!LL
!LL  3 Subroutines in deck : NUM_OBS , NUM_OBS2  and  CHECK_OBS
!LL
!LL  Purpose : Calculate no of observations and no of observation
!LL            values in AC observation files. These values are
!LL            used to dimension any observation data arrays in the
!LL            assimilation code.
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL   3.2   25/05/93  Call SETPOS to get to start of obs files. For
!LL                   new format files only. D. Robinson
!LL           8/7/93  Eliminate QA FORTRAN complaints    S Bell
!LL
!LL   3.3   23/11/93  Include ACOBS/MODEL level and grid consistency
!LL                   checks.  G. Bason.
!LL   4.0   21/2/95  Add control for CHECK_OBS,        G. Bason.
!LL                  and relax check of real constants between file
!LL                  header and model values.          R. Rawlins.
!LL
!LL   4.0   20/07/95  Remove reference to old format acobs GBason
!LL
!LL   4.0   02/11/95  Remove redundant code. (N.Farnon)
!     4.1   18/06/96  Changes to cope with changes in STASH addressing
!                     Author D.M. Goddard.
!    4.2 25/11/96: T3E mods Stuart Bell
!LL
!LL  4.1  4/09/96:  Port to CRAY T3E  Deborah Salmond
!    4.4  3/7/97:  Set PER_FILE_TNDVMAX.   Deborah Salmond
!    4.5  28/7/98: In MPP runs do on 1pe &en broadcast. Deborah Salmond
!    4.5  5/8/98:  increase OBS_FILE_NAME,OBS_DIR_NAME size. Stuart Bell
!    5.2 22/3/01:  amend glsize. Bruce Macpherson
!    5.3 25/9/01:  Portability changes.     Z. Gardner
!    5.3 03/10/01:  comment out call to check_obs
!                   until revised for ND. B Macpherson
!    6.0 10/10/03  : Replace SHMEM with GCOM for SX6 Clive Jones
!    6.1 18/08/04: Declare function LOC as an INTEGER.       P.Dando
!    6.2 21/10/05: Replace GSYNC with SSYNC. P.Selwood
!
!LL  Programming Standard : UM Doc Paper No 3 ; Version 4 ; 5/2/92
!LL
!LL  Project Task : P3
!LL
!LL  Documentation : UM Documentation Paper No ??
!LLEND------------------------------------------------------------------
      SUBROUTINE NUM_OBS (NFILES,NOBSMAX,TNDVMAX,P_LEVELS,Q_LEVELS,     &
     &                    P_ROWS,ROW_LENGTH,AK,BK,REALHD1,REALHD2,      &
     &                    REALHD3,REALHD4,REALHD5,REALHD6,              &
#include "argppx.h"
     &                    ICODE,CMESSAGE)
!L----------------------------------------------------------------------
!     Do not autotask this routine
!FPP$ NOCONCUR R
!
      IMPLICIT NONE

#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
!     AC Comdecks
#include "acparm.h"
#include "comobs.h"
#include "comacp.h"
#include "parvars.h"
#include "chsunits.h"
#include "cenvir.h"
#include "clfhist.h"
!-----------------------------------------------------------------------
!     ARGUMENTS
!-----------------------------------------------------------------------
      INTEGER NFILES        ! IN  No of obs files to be read
      INTEGER NOBSMAX       ! OUT Total no of obs in obs files
      INTEGER TNDVMAX       ! OUT Total no of data values in obs files
      INTEGER P_LEVELS      ! IN  No of model levels
      INTEGER Q_LEVELS      ! IN  No of model wet levels
      INTEGER P_ROWS        ! IN  No of model rows
      INTEGER ROW_LENGTH    ! IN  No of points on row
      REAL    AK(P_LEVELS)  ! IN  Vertical grid
      REAL    BK(P_LEVELS)  !
      REAL REALHD1,REALHD2  ! IN  Horizontal grid
      REAL REALHD3,REALHD4
      REAL REALHD5,REALHD6
      INTEGER ICODE          ! OUT Return code
      CHARACTER*256 CMESSAGE ! OUT Error message
!-----------------------------------------------------------------------
!     LOCAL VARIABLES
!-----------------------------------------------------------------------
      INTEGER IDUMMY
      INTEGER IO_STAT
      REAL DATALEVS(P_LEVELS+1,NOBTYPMX)
      INTEGER JF,JOBT,JLEV   !  Loop counters
      INTEGER UNIT_NO        !  Unit no for obs file
      INTEGER ISPARE(7)
      INTEGER ISTAT          ! GCOM status
!-----------------------------------------------------------------------
      INTEGER TNDV, INDV  (MAX_NUM_ACOB_FILES)
      INTEGER TNOBS,INOBS (MAX_NUM_ACOB_FILES)
!-----------------------------------------------------------------------
#include "dump_len.h"
      PARAMETER (LEN_FIXHD=256)
      INTEGER FIXHD(LEN_FIXHD)
      INTEGER LEN_DATA
!-----------------------------------------------------------------------
      INTEGER OBS_FILE_YY,OBS_FILE_MM,OBS_FILE_DD
      INTEGER OBS_FILE_HH,OBS_FILE_MIN,OBS_FILE_SEC
      INTEGER DIRNUM, FILENUM
      CHARACTER*256 OBS_FILE_NAME
      CHARACTER*256 OBS_DIR_NAME
      LOGICAL LFILE_EXISTS
      INTEGER ENVVAR
      INTEGER LEN_OBS_FILE_NAME
      INTEGER LEN_OBS_DIR_NAME
      INTEGER IFOUND
#if defined(MPP)
      integer ibcast_buf(4),iadd1,iadd2,len_get
#endif
#if !defined(MPP)
      INTEGER  mype
      PARAMETER (mype = 0 ) ! always zero in non-MPP code
#endif
!-----------------------------------------------------------------------
!     FUNCTION DECLARATIONS
!-----------------------------------------------------------------------
      INTEGER LOC  ! Function to get the value of an address.
                   ! Only available on CRAY and Met Office NEC SX-6.
!-----------------------------------------------------------------------
      EXTERNAL GET_DIM, READ_FLH, TIMER, NUM_OBS2, SETPOS               &
     &         ,FORT_GET_ENV
#if defined(MPP)
        EXTERNAL OPEN_SINGLE
#else
        EXTERNAL FILE_OPEN
#endif
!-----------------------------------------------------------------------
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER ('NUMOBS  ',3)
!-----------------------------------------------------------------------

      OBS_UNITNO=70
      if(mype == 0)then
      PRINT *, 'READ IN AC OBS FILES - Headers only'

      IF (OBS_FORMAT  ==  3) THEN
        NFILES = NFILES * NUM_OB_FILE_TYPES
      END IF

      NUM_USED_FILES = 0
      ENVVAR = 1

      DO 2000 JF=1,NFILES   !  Loop over observation files

        IF (OBS_FORMAT  ==  3) THEN
          DIRNUM=INT((JF-1)/(NUM_OB_FILE_TYPES))+1
          FILENUM=MOD(JF-1,NUM_OB_FILE_TYPES)+1
          CALL FORT_GET_ENV(FT_ENVIRON(OBS_UNITNO+DIRNUM-1),            &
     &                      LEN_FT_ENVIR(OBS_UNITNO+DIRNUM-1),          &
     &                      OBS_DIR_NAME,256,ICODE)
          LEN_OBS_DIR_NAME=(INDEX(OBS_DIR_NAME," ") - 1)
          OBS_FILE_NAME = OBS_DIR_NAME(1:LEN_OBS_DIR_NAME)//'/'//       &
     &    OB_FILE_TYPE(FILENUM)(1:INDEX(OB_FILE_TYPE(FILENUM)," ")-1)   &
     &    //'.acobs'
          LEN_OBS_FILE_NAME=(INDEX(OBS_FILE_NAME," ") - 1)
        ELSE
          CALL FORT_GET_ENV(FT_ENVIRON(OBS_UNITNO+JF-1),                &
     &                      LEN_FT_ENVIR(OBS_UNITNO+JF-1),              &
     &                      OBS_FILE_NAME,256,ICODE)
          LEN_OBS_FILE_NAME=(INDEX(OBS_FILE_NAME," ") - 1)
        END IF

        INQUIRE(FILE=OBS_FILE_NAME,IOSTAT=ICODE,EXIST=LFILE_EXISTS)
        IF (ICODE  /=  0) GOTO 9999

        IF (LFILE_EXISTS) THEN
        PRINT *,OBS_FILE_NAME,LEN_OBS_FILE_NAME

          NUM_USED_FILES = NUM_USED_FILES + 1
          USED_FILES(NUM_USED_FILES) = OBS_FILE_NAME
          FILENAME_LEN(NUM_USED_FILES) = LEN_OBS_FILE_NAME

#if defined(MPP)
          CALL OPEN_SINGLE(OBS_UNITNO,OBS_FILE_NAME,                    &
     &              FILENAME_LEN(NUM_USED_FILES),0,ENVVAR,ICODE)
! set flag so BUFFIN called from READ_FLH and READHEAD
! does not broadcast the header data
          call set_unit_bcast_flag(OBS_UNITNO)

#else
! DEPENDS ON: file_open
          CALL FILE_OPEN(OBS_UNITNO,OBS_FILE_NAME,                      &
     &              FILENAME_LEN(NUM_USED_FILES),0,ENVVAR,ICODE)
#endif

          INOBS(NUM_USED_FILES) = 0
          INDV (NUM_USED_FILES) = 0

!       Get unit number of observation file
          UNIT_NO = OBS_UNITNO

        PRINT '(A,I8)', ' AC OBS FILE - UNIT NO :',UNIT_NO

        NOBTYP = 0
!       Read File contents
!       If FILE is empty - proceed to read next file.

!         Go to start of obs file
! DEPENDS ON: setpos
          CALL SETPOS (UNIT_NO,0,ICODE)

!         Read in fixed length header (FLH)
! DEPENDS ON: read_flh
          CALL READ_FLH (UNIT_NO,FIXHD,LEN_FIXHD,                       &
     &                   ICODE,CMESSAGE)
          IF (ICODE >  0) GOTO 9999

!         Get dimensions of all data set components from FLH
! DEPENDS ON: get_dim
          CALL GET_DIM (FIXHD,                                          &
#include "dump_ar2.h"
     &                  LEN_DATA)

!         Get date/time of observation file
          OBS_FILE_YY     = FIXHD(21)
          OBS_FILE_MM     = FIXHD(22)
          OBS_FILE_DD     = FIXHD(23)
          OBS_FILE_HH     = FIXHD(24)
          OBS_FILE_MIN    = FIXHD(25)
          OBS_FILE_SEC    = FIXHD(26)


! DEPENDS ON: num_obs2
          CALL NUM_OBS2 (UNIT_NO,P_LEVELS,Q_LEVELS,P_ROWS,ROW_LENGTH,   &
     &                   AK,BK,REALHD1,REALHD2,REALHD3,REALHD4,         &
     &                   REALHD5,REALHD6,                               &
#include "dump_ar2.h"
     &                   LEN_DATA,NOBTYP,TNOBS,TNDV,                    &
#include "argppx.h"
     &                   ICODE,CMESSAGE)

          IF (ICODE >  0) GOTO 9999

!       Record values for this file
          INOBS(NUM_USED_FILES) = TNOBS
          INDV (NUM_USED_FILES) = TNDV

        PRINT *, ' '
        PRINT '(A,I3.2,I2.2,A,I4.2,A,I2.2,A,I4)',                       &
     &                     ' TIME and DATE         :',                  &
     &  OBS_FILE_HH,OBS_FILE_MIN,'Z',                                   &
     &  OBS_FILE_DD,'/',OBS_FILE_MM,'/',OBS_FILE_YY
        PRINT '(A,I8)', ' No of Obs Types       :',NOBTYP
        PRINT '(A,I8)', ' No of Observations    :',TNOBS
        PRINT '(A,I8)', ' No of Data Values     :',TNDV

#if defined(MPP)
! clear the BUFFIN broadcast flag
          call clear_unit_bcast_flag(OBS_UNITNO)
#endif
! DEPENDS ON: file_close
          CALL FILE_CLOSE(OBS_UNITNO,OBS_FILE_NAME,                     &
     &                    FILENAME_LEN(NUM_USED_FILES),ENVVAR,0,ICODE)
          END IF
      NFILES=NUM_USED_FILES

 2000 CONTINUE

      PRINT '(/,'' Obs file no       '',5I8)',                          &
     &      (JF,JF=1,NFILES)
      PRINT '(  '' No of obs         '',5I8)',                          &
     &      (INOBS(JF),JF=1,NFILES)
      PRINT '(  '' No of data values '',5I8)',                          &
     &      (INDV(JF),JF=1,NFILES)

!     Add up INOBS and INDV to get NOBSMAX and TNDVMAX

      NOBSMAX = 0
      TNDVMAX = 0

      DO JF=1,NFILES
        IF (INOBS(JF) >= 0) THEN
          NOBSMAX = NOBSMAX + INOBS(JF)
        ENDIF
        IF (INDV(JF) >= 0) THEN
          TNDVMAX = TNDVMAX + INDV (JF)
        ENDIF
        PER_FILE_TNDVMAX(JF)=INDV (JF)
      ENDDO

      IF (NOBSMAX == 0) THEN
        NOBSMAX = 1   !  Reset to 1 for dynamic allocation
        PRINT *, ' NOBSMAX = ',NOBSMAX
        PRINT *, ' Reset to 1 prevent allocation problems.'
      ENDIF

      IF (TNDVMAX == 0) THEN
        TNDVMAX = 1   !  Reset to 1 for dynamic allocation
        PRINT *, ' TNDVMAX = ',TNDVMAX
        PRINT *, ' Reset to 1 prevent allocation problems.'
      ENDIF

      PRINT '(/,'' Values calculated by NUM_OBS '',/,                   &
     &           '' Total no of obs (NOBSMAX)         = '',I9,/,        &
     &           '' Total no of data values (TNDVMAX) = '',I9/)',       &
     &            NOBSMAX,TNDVMAX

9999  CONTINUE
#if defined(MPP)
      ibcast_buf(1)=NFILES
      ibcast_buf(2)=NOBSMAX
      ibcast_buf(3)=TNDVMAX
      ibcast_buf(4)=ICODE

      endif ! if(mype == 0)

      CALL GC_SSYNC(NPROC,ISTAT)

! Broadcast the NUM_OBS output arguments to all pes
      CALL GC_IBCAST(1,4,0,NPROC,ISTAT,ibcast_buf)

      NUM_USED_FILES=ibcast_buf(1)
      iadd1=loc(NOBTYP)
      iadd2=loc(NUM_USED_FILES)
      len_get=(iadd2-iadd1)/8+1

! Broadcast COMOBS common block to all pes
      CALL GC_IBCAST(1,len_get,0,NPROC,ISTAT,NOBTYP)

      CALL GC_SSYNC(NPROC,ISTAT)

      NFILES=ibcast_buf(1)
      NOBSMAX=ibcast_buf(2)
      TNDVMAX=ibcast_buf(3)
      ICODE=ibcast_buf(4)

#else
      endif ! if(mype == 0)
#endif

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER ('NUMOBS  ',4)
      RETURN
      END SUBROUTINE NUM_OBS
!LL  SUBROUTINE NUM_OBS2 -----------------------------------------------
!LL
!LL  Written 14/7/92 by Dave Robinson.
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  Programming Standard : UM Doc Paper No 3 ; Version 4 ; 5/2/92
!LL
!LL  Project Task : P3
!LL
!LL  Purpose : Read in header section of observation file and get
!LL            number of observations and data values in this obs file.
!LL
!LL  Documentation : UM Documentation Paper No ??
!LLEND------------------------------------------------------------------
!-----------------------------------------------------------------------
#endif
