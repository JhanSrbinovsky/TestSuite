#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  4 Subroutines in deck : RDOBS, RDOBS2, RDOBS3 and DAYS -----
!LL
!LL  Purpose : Read from ACOBS Files,reformat and place OBS header
!LL            details in COMOBS. The bulk of the required OBS data
!LL            is put into dynamic work array OBS for transmission via
!LL            argument list to GETOBS. OBS is written out to a cache
!LL            file for subsequent reading at later timesteps.
!LL            Thus reread of ACOBS files only required intermittently
!LL            (The routine DAYS does a dd/mm/yy to dayno)
!LL
!LL  For use on Cray
!LL  For Cray - Global  ; Enable defs GLOBAL
!LL
!LL S.Bell      <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL   3.2   25/05/93  Call SETPOS to get to start of new format obs
!LL                   files. Call READ_FLH to get FIXHD(162). Correct
!LL                   argument list and dimensions in RDOBS3. Tidy up
!LL                   print out of observation levels. D Robinson.
!LL           8/7/93  Eliminate QA FORTRAN complaints   S Bell
!LL
!LL   3.2   10/07/93  Cater for type 306    S Bell
!LL   3.3   18/04/94  Modify arglist to READDUMP.    T Johns
!LL   3.3   25/11/93  Correct setting of headers in RDOBS2 Greg Bason
!LL   4.0   20/07/95  Remove references to old format acobs GBason
!LL
!LL   4.0   02/11/95  Remove redundant code. (N.Farnon)
!     4.1   18/06/96  Changes to cope with changes in STASH addressing
!                     Author D.M. Goddard.

!LL  4.1   31/05/96     The number of v points to be processed on a
!LL                     C grid differs from u by row_length. u,v
!LL                     dimensioned separately in call to WLLTOEQ.
!LL                     Requirement for VAR.
!LL                     Author I.Edmond       Reviewer D. Goddard
!     4.2   25/11/96  T3E mods + adjust time window  Stuart Bell

!LL
!LL  4.1  4/09/96:  Port to CRAY T3E  Deborah Salmond
!LL  4.2   5/9/96     Fix overwriting triggered by 4.1 mod above. S.Bell
!LL  4.3   19/03/97   Removed MPP args from READDUMP     P.Burton
!    4.3  31/1/97: Bugfix for LAM on T3E. Stuart Bell
!    4.4  2/6/97: Fix to use less memory and allow more obs. SB/DS
!    4.4  17/10/97: Fix to ensure obs on 1 PE only. Deborah Salmond
!    5.2  12/11/00: change attop,base,right,left to elements of
!                   at_extremity              B Macpherson
!    5.3  09/07/01: amend for S->N ND grid order
!                     Bruce Macpherson
!    5.3  05/12/01:  Remove reference to the shmcomm & iovarsac include
!                    files, use a local dynamic array rather than
!                    common block for array 'work' in rdobs2.  S.Cusack
!    6.0 11/09/03:   Removed double ? for IBM cpp.             P.Dando
!    6.0  10/10/03:  Replace SHMEM with GCOM for SX6. Clive Jones
!    6.0  30/12/03:  Make argument list to READACOBS consistent with
!                    subroutine. Dave Robinson.
!    6.1  17/08/04:  Amend boundaries for correct allocation of
!                    observations to processors. Adam Maycock,
!                    (J.Bornemann lodged).
!    6.2  21/10/05:  Replace GSYNC with SSYNC. P.Selwood
!    6.2  15/08/05   Free format fixes. P.Selwood
!    6.2  24/01/06:  Replace large hard-wired MAX_SHMEM_SIZE and
!                    obs(num)dim with dynamic allocation. R Barnes
!
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Logical components covered:
!LL
!LL  Project Task : P3
!LL
!LL  External documentation:
!LL
!LLEND------------------------------------------------------------------
      SUBROUTINE RDOBS2(NFILES,TIMESTEP_NO,OBS,OBS_FLAG,TNDV,           &
     &                  TNOBS,TNDVMAX,NOBSMAX,P_LEVELS,Q_LEVELS,TIMEREL,&
#include "argppx.h"
     &               ICODE,CMESSAGE)
!L----------------------------------------------------------------------
!L   INTENT IN:
!L      NFILES   : NO OF AC OBSERVATION FILES TO BE READ
!L      TIMESTEP_NO : TIMESTEP NUMBER
!L      TIMEREL     : relative time for this timestep
!L      TNDVMAX  : MAX SIZE OF OBS ARRAY
!L      NOBSMAX  : MAX NO OF OBS (FOR DIMENSIONING)
!L      P_LEVELS : NUMBER OF MODEL LEVELS
!L      Q_LEVELS : NUMBER OF WET MODEL LEVELS
!L   INTENT OUT:
!L      TNDV     : ACTUAL SIZE OF OBS ARRAY
!L      TNOBS    : ACTUAL NO OF OBS
!L      OBS      : OBS array
!L      OBS_FLAG : do not use flags
!L      ICODE/CMESSAGE: for error processing
!L----------------------------------------------------------------------
!L
!
!     Do not autotask this routine
!FPP$ NOCONCUR R
!
      IMPLICIT NONE

! EXTERNAL SUBROUTINE CALLS
        EXTERNAL DAYS,SETDAC,READ_FLH,RDOBS3,GET_DIM
#if !defined(GLOBAL)
        EXTERNAL EQTOLL,LLTOEQ,W_COEFF,W_LLTOEQ
#endif
#if defined(MPP)
        EXTERNAL OPEN_SINGLE
#else
        EXTERNAL FILE_OPEN
#endif
#include "acparm.h"
#include "comobs.h"
#include "comacp.h"
#include "commg.h"
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#if defined(MPP)
#include "parvars.h"
#include "mppac.h"
      integer common_length
#endif

!-----------------------------------------------------------------------
!     ARGUMENTS
!-----------------------------------------------------------------------
      INTEGER TNDV,TNOBS,MAX_NDV
      INTEGER P_LEVELS,Q_LEVELS
      INTEGER TNDVMAX,NOBSMAX
      INTEGER NFILES
      INTEGER TIMESTEP_NO
      INTEGER ICODE
      CHARACTER*256 CMESSAGE
!-----------------------------------------------------------------------
!     LOCAL VARIABLES
!-----------------------------------------------------------------------
      INTEGER KTYPE
      INTEGER IO_STAT
      INTEGER JF,JOBT,JOBT2,JDV,J,JOB,JTYP,JLEV,JACT,JMOT
      INTEGER ITOT0,ITOT1,ITOTAL
      INTEGER IRDAY,IFDAY
      INTEGER IPT,IPC,IFILE,I,IWORKSP,IOBT,IJF
      INTEGER IPT_THIS_FILE
      INTEGER MAX_NDATAV
      INTEGER INDVMAX
      REAL                                                              &
#if !defined(MPP)
     & ZZLATMX,ZZLATMN,                                                 &
#else
     & ZZLATMX,ZZLATMN,ZZLONGMX,ZZLONGMN,                               &
#endif
     & TIMEREL,TIMEADJ,TWSTART,TWEND,TGETOBB,TGETOBA
#if defined(MPP)
      COMMON/ZZ/ZZLATMX,ZZLATMN,ZZLONGMX,ZZLONGMN
      INTEGER iproc,istat,imsg
#endif
!-----------------------------------------------------------------------
#include "c_pi.h"
!-----------------------------------------------------------------------
      REAL DATALEVS (P_LEVELS+1,NOBTYPMX)
      INTEGER IREF     (NOBTYPMX,NDATAVMX,NFILES)
      INTEGER INOBS    (NOBTYPMX,NFILES)
      INTEGER IOBSTYP  (NOBTYPMX,NFILES)
      INTEGER INDATAV  (NOBTYPMX,NFILES)
      INTEGER INOBLEV  (NOBTYPMX,NFILES)
      INTEGER IOBLVTYP (NOBTYPMX,NFILES)
      INTEGER INDVHDR  (NFILES)
      INTEGER INOBTYP  (NFILES)
      INTEGER IMAXNLEV (NFILES)
      REAL PLEVELS     (P_LEVELS+1,NOBTYPMX,NFILES)
      INTEGER KOBSTYP  (NOBTYPMX)
      LOGICAL LEMPTY   (NFILES)
!-----------------------------------------------------------------------
      INTEGER OBS_FILE_YY,OBS_FILE_MM,OBS_FILE_DD
      INTEGER OBS_FILE_HH,OBS_FILE_MIN,OBS_FILE_SEC
      INTEGER IP_LAT,IP_LONG,IP_TIME,IP_TYPE,IP_MOT,IP_U,IP_V
      INTEGER IP_MOT1,IP_MOT2,IP_TYPE1,IP_TYPE2
      INTEGER IDUMMY
!-----------------------------------------------------------------------
#include "dump_len.h"
      PARAMETER (LEN_FIXHD=256)
      INTEGER FIXHD(LEN_FIXHD)
      INTEGER LEN_DATA
!-----------------------------------------------------------------------
!     DYNAMIC ALLOCATION
      INTEGER OBS_FLAG(NOBSMAX)
      INTEGER IWORK(NOBSMAX)
      REAL OBS(TNDVMAX)
      REAL WORK(TNDVMAX+2048)
#if !defined(GLOBAL)
      REAL U_WRK(NOBSMAX),V_WRK(NOBSMAX)
      REAL WRKLAT(NOBSMAX),WRKLON(NOBSMAX)
      REAL COEFF1(NOBSMAX),COEFF2(NOBSMAX)
#endif
#if !defined(MPP)
      INTEGER mype
      PARAMETER (mype = 0 ) ! always zero in non-MPP code
#endif
      INTEGER ISPARE(7)
      INTEGER ENVVAR
!-----------------------------------------------------------------------
!L           SECTION 1: COPY INPUT FILES TO WORK AREA, ETC.
!-----------------------------------------------------------------------

!       Read in the AC Observation files and merge the observations.

        IF (DIAG_RDOBS >= 1.AND.mype == 0) THEN
          PRINT *, ' '
          PRINT *, 'READ IN AC OBS FILES - TIMESTEP : ',TIMESTEP_NO
        ENDIF

!       Set up time for next read of observation files
        TIMENEXT = TIMEREL + TIMEINT

!       Initilaise to zero to allow for empty files.
        DO JF=1,NFILES
          LEMPTY  (JF) = .FALSE.
          INDVHDR (JF) = 0
          INOBTYP (JF) = 0
          IMAXNLEV(JF) = 0
          DO JOBT=1,NOBTYPMX
            IOBSTYP  (JOBT,JF) = 0
            INDATAV  (JOBT,JF) = 0
            INOBS    (JOBT,JF) = 0
            INOBLEV  (JOBT,JF) = 0
            IOBLVTYP (JOBT,JF) = 0
            DO JLEV=1,P_LEVELS+1
              PLEVELS(JLEV,JOBT,JF) = 0.0
            ENDDO
          ENDDO
        ENDDO

!       Initialise IREF
        DO JF = 1,NFILES
          DO JDV = 1,NDATAVMX
            DO JOBT = 1,NOBTYPMX
              IREF(JOBT,JDV,JF) = 0
            ENDDO
          ENDDO
        ENDDO

!       IPT IS THE AMOUNT OF WORK SPACE USED IN THE ARRAY WORK SO FAR
        IPT=0

!       JF LOOP: ONE CYCLE FOR EACH INPUT FILE.
        DO 2000 JF=1,NUM_USED_FILES

        IFILE = OBS_UNITNO

!       Get unit number of observation file

        ENVVAR = 1
        if(mype == 0)then
#if defined(MPP)
        CALL OPEN_SINGLE(IFILE,USED_FILES(JF),                          &
     &                   FILENAME_LEN(JF),0,ENVVAR,ICODE)
#else
! DEPENDS ON: file_open
        CALL FILE_OPEN(IFILE,USED_FILES(JF),                            &
     &                 FILENAME_LEN(JF),0,ENVVAR,ICODE)
#endif
        endif

!       READ FILE CONTENTS TO BUFFER.
!       IF FILE IS EMPTY - PROCEED TO READ NEXT FILE

        NOBTYP = 0
        TNDV   = 0
        TNOBS  = 0
        NDVHDR = 0
        DO JOBT=1,NOBTYPMX
          OBSTYP  (JOBT) = 0
          NDATAV  (JOBT) = 0
          NOBLEV  (JOBT) = 0
          OBLEVTYP(JOBT) = 0
          NOBS    (JOBT) = 0
        ENDDO

!-----------------------------------------------------------------------
        IF (OBS_FORMAT == 2 .OR. OBS_FORMAT == 3) THEN

!         Go to start of obs file
! DEPENDS ON: setpos
          CALL SETPOS (IFILE,0,ICODE)

!         Read in fixed length header
! DEPENDS ON: read_flh
          CALL READ_FLH (IFILE,FIXHD,LEN_FIXHD,ICODE,CMESSAGE)
          IF (ICODE >  0) GOTO 9999

!         Get dimensions of data set components in Fixed length header
! DEPENDS ON: get_dim
          CALL GET_DIM (FIXHD,                                          &
#include "dump_ar2.h"
     &                  LEN_DATA)

          OBS_FILE_YY  = FIXHD(21)   ! )
          OBS_FILE_MM  = FIXHD(22)   ! ) Date/time
          OBS_FILE_DD  = FIXHD(23)   ! ) of
          OBS_FILE_HH  = FIXHD(24)   ! ) Observation
          OBS_FILE_MIN = FIXHD(25)   ! ) File
          OBS_FILE_SEC = FIXHD(26)   ! )

          MAX_NDV      = FIXHD(162)  !  Max total no of data values

!         If MAX_NDV = 0, set to 1. (Obs file with no observations)
!         Prevents dynamic allocation with 0 in READACOBS
          IF (MAX_NDV == 0) MAX_NDV = 1

#if defined(MPP)
!
!  check dimension of WORK
!
      common_length = 2048+IPT+PER_FILE_TNDVMAX(JF)
! 2048 allows space for headers when no or very few obs; was 3000.
      if (common_length  >   TNDVMAX+2048) then
      write(6,*)'RDOBS2 - common_length,IPT,PER_FILE_TNDVMAX,'          &
     &,'TNDVMAX+2048 ',common_length,IPT,PER_FILE_TNDVMAX,              &
     & TNDVMAX+2048
!
! important failure messages to PEn,OUT and OUT2 output streams
      WRITE(0,*)' RDOBS2 : Insufficient space in WORK array'
      WRITE(0,*) ' Recode or rerun with fewer obs'
      WRITE(0,*)'TNDVMAX',TNDVMAX,' should be >= ',common_length-2048
      WRITE(6,*)' RDOBS2 : Insufficient space in WORK array'
      WRITE(6,*) ' Recode or rerun with fewer obs'
      WRITE(6,*)'TNDVMAX',TNDVMAX,' should be >= ',common_length-2048
      WRITE(7,*)' RDOBS2 : Insufficient space in WORK array'
      WRITE(7,*) ' Recode or rerun with fewer obs'
      ICODE = 1
      CMESSAGE =' RDOBS2 : Insufficient space in WORK array'
      GO TO 9999
      else
       IF (DIAG_RDOBS >= 1) THEN
! diagnostic message directed to every PEn output stream
       IPC=(common_length*100)/(TNDVMAX+2048)
       WRITE(6,*)'RDOBS2:% of WORK used so far for reading obs=',IPC
       ENDIF
      endif
#endif
! DEPENDS ON: rdobs3
          CALL RDOBS3(IFILE,NOBTYPMX,NOBTYP,OBSTYP,NOBS,NDATAV,         &
     &                NOBLEV,OBLEVTYP,DATALEVS,WORK(IPT+1),             &
#include "dump_ar2.h"
     &                LEN_DATA,MAX_NDV,TNOBS,MISSD,P_LEVELS,            &
#include "argppx.h"
     &                ICODE,CMESSAGE                                    &
#if defined(MPP)
     &                         ,IPT                                     &
#endif
     &                             )
          IF (ICODE >  0) GOTO 9999

          NDVHDR = 5                 ! No of data values in header
          MAXNLEV1 = LEN1_LEVDEPC-2  ! Max no of levels + 1

        ELSE
          ICODE=1
          CMESSAGE='RDOBS2: ILLEGAL OBS_FORMAT'
          GOTO 9999
        ENDIF

!       Convert any old obs type numbers to new obs type numbers.
!       For version 2.6 onwards.
        IF (NOBTYP >  0) THEN
          DO JOBT=1,NOBTYP
            IF (OBSTYP(JOBT) == 501) THEN
              OBSTYP(JOBT) = 302
             if(mype == 0)PRINT *, 'Type 501 in Obs file changed to 302'
            ENDIF
            IF (OBSTYP(JOBT) == 502) THEN
              OBSTYP(JOBT) = 305
             if(mype == 0)PRINT *, 'Type 502 in Obs file changed to 305'
            ENDIF
          ENDDO
        ENDIF

!       PRINT CONTENTS OF FILE HEADER.
        IF (DIAG_RDOBS >= 1.AND.mype == 0) THEN
          PRINT *, ' '
          PRINT '(A,I8)', ' AC OBS FILE - UNIT NO :',IFILE
          PRINT '(A,I8)', ' NO OF OBS TYPES       :',NOBTYP
          PRINT '(A,I8)', ' TOTAL NO OF OBS       :',TNOBS
          PRINT *, ' '
          PRINT '(A,I3.2,A,I2.2,A,I4)',                                 &
     &      ' DATE :',OBS_FILE_DD,'/',OBS_FILE_MM,'/',OBS_FILE_YY
          PRINT '(A,I3.2,I2.2,A)',                                      &
     &      ' TIME :',OBS_FILE_HH,OBS_FILE_MIN,'Z'
          PRINT *, ' '

          IF (NOBTYP >  0) THEN

            PRINT '(A,T30,18I6/(T30,18I6))',                            &
     &           ' AC OBS TYPES      :',(OBSTYP(I),I=1,NOBTYP)
            PRINT '(A,T30,18I6/(T30,18I6))',                            &
     &            ' NO OF LEVELS      :',(NOBLEV(I),I=1,NOBTYP)
            PRINT '(A,T30,18I6/(T30,18I6))',                            &
     &            ' NO OF DATA VALUES :',(NDATAV(I),I=1,NOBTYP)
            PRINT '(A,T30,18I6/(T30,18I6))',                            &
     &            ' OBS LEVEL TYPE    :',(OBLEVTYP(I),I=1,NOBTYP)
            PRINT '(A,T30,18I6/(T30,18I6))',                            &
     &           ' NO OF OBS         :',(NOBS(I),I=1,NOBTYP)

          ENDIF

        ENDIF

        IF (TIMESTEP_NO == 1 .AND. TNOBS == 0) THEN
          IF (OBSTYP(1) == 406) THEN
            WRITE (7,'(A,I3.2,I2.2,A,I3.2,A,I2.2,A,I4,A)')              &
     &        ' AC Observation File (MOPS) - ',                         &
     &        OBS_FILE_HH,OBS_FILE_MIN,'Z',                             &
     &        OBS_FILE_DD,'/',OBS_FILE_MM,'/',OBS_FILE_YY,              &
     &        ' - No MOPS ACOBS data'
          ELSE
          WRITE (7,'(A,I3.2,I2.2,A,I3.2,A,I2.2,A,I4,A)')                &
     &        ' AC Observation File - ',                                &
     &        OBS_FILE_HH,OBS_FILE_MIN,'Z',                             &
     &        OBS_FILE_DD,'/',OBS_FILE_MM,'/',OBS_FILE_YY,              &
     &        ' - No standard ACOBS data'
          ENDIF
        ENDIF

!       Check that no of observation types in file (NOBTYP) does not
!       exceed maximum allowed. (NOBTYPMX)
        IF (NOBTYP >  NOBTYPMX) THEN
          ICODE = 1
          CMESSAGE = ' RDOBS2 : TOO MANY OBSERVATION TYPES IN OBS FILE'
          if(mype == 0)then
          PRINT *,' RDOBS2 : TOO MANY OBSERVATION TYPES IN OBS FILE'
          PRINT '(A,I5)',' NO OF OBS TYPES = ',NOBTYP
          PRINT '(A,I5)',' MAXIMUM ALLOWED = ',NOBTYPMX
          endif
          GO TO 9999
        ENDIF

!       Check that no of data values for each obs type (NDATAV) does
!       not exceed maximum allowed. (NDATAVMX)
        IF (NOBTYP >  0) THEN
          DO JOBT=1,NOBTYP
          IF (NDATAV(JOBT) >  NDATAVMX) THEN
            ICODE = 1
            CMESSAGE = ' RDOBS2 : TOO MANY DATA VALUES IN OBS FILE'
            if(mype == 0)then
            PRINT *, ' RDOBS2 : Too many Data values in Obs file'
            PRINT *, ' Observation Types =  ',OBSTYP(JOBT)
            PRINT *, ' No of Data Values =  ',NDATAV(JOBT)
            PRINT *, ' Maximum allowed   =  ',NDATAVMX
            endif
            GO TO 9999
          ENDIF
          ENDDO
        ENDIF

!       Store information for this obs file
        INOBTYP(JF) = NOBTYP
        INDVHDR(JF) = NDVHDR
        IMAXNLEV(JF) = MAXNLEV1
        DO JOBT=1,NOBTYP
          IOBSTYP(JOBT,JF)  = OBSTYP(JOBT)
          INDATAV(JOBT,JF)  = NDATAV(JOBT)
          INOBLEV(JOBT,JF)  = NOBLEV(JOBT)
          IOBLVTYP(JOBT,JF) = OBLEVTYP(JOBT)
          DO JLEV=1,MAXNLEV1
            PLEVELS(JLEV,JOBT,JF) = DATALEVS(JLEV,JOBT)
          ENDDO
        ENDDO

!       EVALUATE POINTERS TO EACH SUB-VECTOR.
!       IPT=0 BEFORE READING FIRST FILE

        IPT_THIS_FILE = IPT

          DO JOBT=1,NOBTYP
            DO JDV=1,NDATAV(JOBT)
              IREF(JOBT,JDV,JF)=IPT
              IPT = IPT+NOBS(JOBT)
            ENDDO
          ENDDO

#if !defined(MPP)
!       TEST WORK SPACE USED SO FAR
        IF (IPT >  TNDVMAX) THEN
          ICODE = 1
          CMESSAGE = ' RDOBS2 : WORK SPACE FOR OBSERVATIONS EXCEEDED'
! required in every pe
          WRITE(6,*)' RDOBS2 : WORK SPACE FOR OBSERVATIONS EXCEEDED'
          WRITE(6,*)' INCREASE VALUE OF TNDVMAX'
          WRITE(6,'(A,I8)')' WORK SPACE USED      =',IPT
          WRITE(6,'(A,I8)')' WORK SPACE AVAILABLE =',TNDVMAX
          GO TO 9999
        ENDIF
#endif

!       CONVERT REFERENCE DATE & FILE DATE TO ELAPSED DAYS

! DEPENDS ON: days
        CALL DAYS (OBS_REF_DD, OBS_REF_MM, OBS_REF_YY, IRDAY)
! DEPENDS ON: days
        CALL DAYS (OBS_FILE_DD,OBS_FILE_MM,OBS_FILE_YY,IFDAY)

!       FIND RELATIVE TIME ADJUSTMENT.

        TIMEADJ = 1440*(IFDAY-IRDAY) +                                  &
     &   60*(OBS_FILE_HH-OBS_REF_HH) + OBS_FILE_MIN-OBS_REF_MIN

        IF (NOBTYP >  0) THEN

!         LOOP OVER OBSERVATION TYPES

          DO 1500 JOBT=1,NOBTYP

            TGETOBB = 0.0
            TGETOBA = 0.0
            DO JOBT2=1,NOBTYPMX
              IF (OBSTYP(JOBT) == MASTER_AC_TYPES(JOBT2)) THEN
                TGETOBB = DEF_TGETOBB(JOBT2)
                TGETOBA = DEF_TGETOBA(JOBT2)
              ENDIF
            ENDDO
            IF (TGETOBB == 0.0 .OR. TGETOBA == 0.0) THEN
              ICODE = 1
              CMESSAGE = 'RDOBS2 : TGETOBB and/or TGETOBA = 0.0 ?'
              GO TO 9999
            ENDIF

!         Set up time window to get observations required for assim.
          TWSTART  = TIMEREL - TGETOBA
          TWEND    = TIMEREL + TIMEINT + TGETOBB

!       FOR WIND OBS RESTRICT DOMAIN TO INSIDE WIND GRID
#if !defined(MPP)
        IF (OBSTYP(JOBT)/100 == 3) THEN
          ZZLATMX = OBS_LAT_N-DLAT*0.5
          ZZLATMN = OBS_LAT_S+DLAT*0.5
        ELSE
          ZZLATMX = OBS_LAT_N
          ZZLATMN = OBS_LAT_S
        ENDIF
#else
        IF (OBSTYP(JOBT)/100 == 3) THEN
           ZZLATMX = lat_n-dlat-dlat*0.5
           ZZLATMN = lat_s-dlat*0.5
         if(at_extremity(PNorth))ZZLATMX=OBS_LAT_N-DLAT*0.5
         if(at_extremity(PSouth))ZZLATMN=OBS_LAT_S+DLAT*0.5
#if defined(GLOBAL)
           ZZLONGMX=long_e+dlong*0.5
           ZZLONGMN=long_w+dlong+dlong*0.5
#else
           ZZLONGMX=long_e_model+dlong*0.5
           ZZLONGMN=long_w_model+dlong+dlong*0.5
         if(at_extremity(PEast))ZZLONGMX=OBS_LONG_E
         if(at_extremity(PWest))ZZLONGMN=OBS_LONG_W
#endif
        ELSE
           ZZLATMX = lat_n
           ZZLATMN = lat_s-0.5*dlat
         if(at_extremity(PNorth))ZZLATMX=OBS_LAT_N
         if(at_extremity(PSouth))ZZLATMN=OBS_LAT_S
#if defined(GLOBAL)
           ZZLONGMX=long_e
           ZZLONGMN=long_w+dlong
#else
           ZZLONGMX= long_e_model+0.5*dlong
           ZZLONGMN= long_w_model
         if(at_extremity(PEast))ZZLONGMX=OBS_LONG_E
         if(at_extremity(PWest))ZZLONGMN=OBS_LONG_W
#endif
        ENDIF
#if !defined(GLOBAL)
        if(ZZLONGMN >= 360.0)ZZLONGMN=ZZLONGMN-360.0
        if(ZZLONGMX >= 360.0)ZZLONGMX=ZZLONGMX-360.0
#endif
#endif

        IP_LAT  = 1
        IP_LONG = 2
        IP_TIME = 3

        IF (NOBS(JOBT) >  0) THEN

#if defined(MPP)
!
!  Ensure boundaries are numerically the same.

        CALL GC_SSYNC(NPROC,ISTAT)

! process PEs not along eastern boundary
        IF(.NOT.AT_EXTREMITY(PEAST))THEN

          DO IPROC=1,NPROC-1

            IMSG=JOBT*100+IPROC  ! message tag

            IF(mype == IPROC) THEN
              IF(.NOT.AT_EXTREMITY(PWEST)) THEN
! this PE receives if not along western boundary
                CALL GC_RRECV(IMSG,1,IPROC-1,ISTAT,ZZLONGMN,ZZLONGMX)
              ENDIF
            ELSEIF(mype == IPROC-1) THEN
! PEs not along eastern boundary send to adjacent PE
              CALL GC_RSEND(IMSG,1,IPROC,ISTAT,ZZLONGMN,ZZLONGMX)
            ENDIF

          ENDDO

        ELSEIF(.NOT.AT_EXTREMITY(PWEST)) THEN
! PEs along eastern boundary can receive only

          DO IPROC=1,NPROC-1
            IF(mype == IPROC) THEN        ! this PE receives
              IMSG=JOBT*100+IPROC
              CALL GC_RRECV(IMSG,1,IPROC-1,ISTAT,ZZLONGMN,ZZLONGMX)
            ENDIF
          ENDDO

        ENDIF

        CALL GC_SSYNC(NPROC,ISTAT)

! Now along southern/northern boundaries

! process PEs not along southern boundary
        IF(.NOT.AT_EXTREMITY(PSOUTH))THEN

          DO IPROC=NPROC-1,0,-1

            IMSG=JOBT*1000+IPROC  ! message tag

            IF(mype == IPROC-nproc_x) THEN
              IF(.NOT.AT_EXTREMITY(PNORTH)) THEN
! this PE receives if not along northern boundary
                CALL GC_RRECV(IMSG,1,IPROC,ISTAT,ZZLATMX,ZZLATMN)
              ENDIF
            ELSEIF(mype == IPROC) THEN
! PEs not along southern boundary send to adjacent PE
              CALL GC_RSEND(IMSG,1,IPROC-nproc_x,ISTAT,ZZLATMX,ZZLATMN)
            ENDIF

          ENDDO

        ELSE
! PEs along southern boundary can receive only

          DO IPROC=NPROC-1,0,-1
            IF(mype == IPROC-nproc_x) THEN        ! this PE receives
              IMSG=JOBT*1000+IPROC
              CALL GC_RRECV(IMSG,1,IPROC,ISTAT,ZZLATMX,ZZLATMN)
            ENDIF
          ENDDO

        ENDIF

        CALL GC_SSYNC(NPROC,ISTAT)
#endif

#if !defined(GLOBAL)
!       ROTATE REAL LAT/LON OF OBS TO ELF CO-ORDINATES
!       WRITE TRANSFORMED LAT,LON BACK TO SAME AREA OF WORK ARRAY
!       OUTPUT LONGITUDES FROM LLTOEQ ARE IN RANGE 0-->360 DEGREES.
! DEPENDS ON: lltoeq
        CALL LLTOEQ                                                     &
     &  (WORK(IREF(JOBT,IP_LAT,JF)+1),WORK(IREF(JOBT,IP_LONG,JF)+1),    &
     &   WORK(IREF(JOBT,IP_LAT,JF)+1),WORK(IREF(JOBT,IP_LONG,JF)+1),    &
     &   ELFPLAT,ELFPLON,NOBS(JOBT) )
#endif

        DO 1400 JOB=1,NOBS(JOBT)

!       Reset Observation time so times are relative to start of assm.

        WORK(IREF(JOBT,IP_TIME,JF)+JOB) =                               &
     &  WORK(IREF(JOBT,IP_TIME,JF)+JOB)+TIMEADJ

!       Test if observation in area and time window.
!   -0.5 on timewindow helps make results same on different machines

#if !defined(MPP)
        IF ( OBS_LONG_W  <   OBS_LONG_E ) THEN
#else
         IF ( ZZLONGMN  <   ZZLONGMX ) THEN
#endif

          IF ( WORK(IREF(JOBT,IP_LAT, JF)+JOB)  <   ZZLATMX .AND.       &
     &         WORK(IREF(JOBT,IP_LAT, JF)+JOB)  >=  ZZLATMN .AND.       &
#if !defined(MPP)
     &         WORK(IREF(JOBT,IP_LONG,JF)+JOB)  <   OBS_LONG_E  .AND.   &
     &         WORK(IREF(JOBT,IP_LONG,JF)+JOB)  >=  OBS_LONG_W  .AND.   &
#else
     &         WORK(IREF(JOBT,IP_LONG,JF)+JOB)  <   ZZLONGMX  .AND.     &
     &         WORK(IREF(JOBT,IP_LONG,JF)+JOB)  >=  ZZLONGMN  .AND.     &
#endif
     &         WORK(IREF(JOBT,IP_TIME,JF)+JOB)  <   TWEND-0.5 .AND.     &
     &         WORK(IREF(JOBT,IP_TIME,JF)+JOB)  >   TWSTART+0.5  ) THEN

!           Count no of obs in area and time window and
!           record those observations.

            INOBS(JOBT,JF)        = INOBS(JOBT,JF)+1
            IWORK(INOBS(JOBT,JF)) = JOB

          ENDIF

        ELSE

          IF ( WORK(IREF(JOBT,IP_LAT, JF)+JOB)  <   ZZLATMX  .AND.      &
     &         WORK(IREF(JOBT,IP_LAT, JF)+JOB)  >=  ZZLATMN  .AND.      &
#if !defined(MPP)
     &       ( WORK(IREF(JOBT,IP_LONG,JF)+JOB)  <   OBS_LONG_E .OR.     &
     &         WORK(IREF(JOBT,IP_LONG,JF)+JOB)  >=  OBS_LONG_W ) .AND.  &
#else
     &       ( WORK(IREF(JOBT,IP_LONG,JF)+JOB)  <   ZZLONGMX .OR.       &
     &         WORK(IREF(JOBT,IP_LONG,JF)+JOB)  >=  ZZLONGMN ) .AND.    &
#endif
     &         WORK(IREF(JOBT,IP_TIME,JF)+JOB)  <   TWEND-0.5 .AND.     &
     &         WORK(IREF(JOBT,IP_TIME,JF)+JOB)  >   TWSTART+0.5  ) THEN

!           Count no of obs in area and time window and
!           record those observations.

            INOBS(JOBT,JF)        = INOBS(JOBT,JF)+1
            IWORK(INOBS(JOBT,JF)) = JOB
          ENDIF
        ENDIF

 1400   CONTINUE

!       Compress out observations not required.
        DO JDV=1,NDATAV(JOBT)
        DO JOB=1,INOBS(JOBT,JF)
        WORK(IREF(JOBT,JDV,JF)+JOB)=WORK(IREF(JOBT,JDV,JF)+IWORK(JOB))
        ENDDO
        ENDDO

        ENDIF

        CALL GC_SSYNC(NPROC,ISTAT)

 1500   CONTINUE
        ENDIF

!       Print no of observations in time window for file JF.
#if defined(MPP)
        DO JOBT=1,NOBTYP
        CountA(JOBT)=INOBS(JOBT,JF)
        ENDDO

        CALL GC_SSYNC(NPROC,ISTAT)

        If(mype == 0)Then
          Do JOBT=1,NOBTYP
            CountC(JOBT)=0
          EndDo
        Endif

        Do iproc=1,nproc-1
          IF(mype == 0) THEN
            CALL GC_IRECV(IPROC,NOBTYP,IPROC,ISTAT,COUNTB,COUNTA)
            Do JOB=1,NOBTYP
              CountC(JOB)=CountC(JOB)+CountB(JOB)
            EndDo
          ELSEIF(mype == IPROC) THEN
            CALL GC_ISEND(IPROC,NOBTYP,0,ISTAT,COUNTB,COUNTA)
          ENDIF
          IF (DIAG_RDOBS >  1.AND.mype == 0) THEN
            PRINT '(A,I3,T30,18I6/(T30,18I6))',                         &
     &            ' NO OF OBS IN T.W: pe=',iproc,                       &
     &            (CountB(JOBT),JOBT=1,NOBTYP)
          ENDIF
          Do JOBT=1,NOBTYP
            CountC(JOBT)=CountC(JOBT)+CountB(JOBT)
          EndDo
        EndDo
        If(mype == 0) THEN
          PRINT *, ' '
          PRINT '(A,T30,18I6/(T30,18I6))',                              &
     &          ' NO OF OBS IN T.W: total:',                            &
     &          (CountC(JOBT),JOBT=1,NOBTYP)
        ENDIF

        CALL GC_SSYNC(NPROC,ISTAT)
#else
        PRINT *, ' '
        PRINT '(A,T30,18I6/(T30,18I6))',                                &
     &        ' NO OF OBS IN T.W  :',(INOBS(JOBT,JF),JOBT=1,NOBTYP)
#endif

!       Compress out unused work space
!       ==============================

!       Following commented line may be brought back in future

#if !defined(MPP)
!DR     Find maximum NDATAV in AC Obs file
!DR     MAX_NDATAV = INDATAV(1,JF)
!DR     DO JOBT=2,NOBTYP
!DR       MAX_NDATAV = MAX (MAX_NDATAV,INDATAV(JOBT,JF))
!DR     ENDDO
!DR
!DR     IPT = IPT_THIS_FILE
!DR     DO JDV=1,MAX_NDATAV
!DR       DO JOBT=1,NOBTYP
!DR       IF (JDV <= NDATAV(JOBT)) THEN
!DR         IF (INOBS(JOBT,JF) >  0) THEN
!DR           DO JOB=1,INOBS(JOBT,JF)
!DR           WORK(IPT+JOB) = WORK(IREF(JOBT,JDV,JF)+JOB)
!DR           ENDDO
!DR         ENDIF
!DR         IREF(JOBT,JDV,JF)=IPT
!DR         IPT = IPT+INOBS(JOBT,JF)
!DR       ENDIF
!DR       ENDDO
!DR     ENDDO
#else
        IPT = IPT_THIS_FILE
        DO JOBT=1,NOBTYP
          DO JDV=1,NDATAV(JOBT)
            IF (INOBS(JOBT,JF) >  0) THEN
              DO JOB=1,INOBS(JOBT,JF)
              WORK(IPT+JOB) = WORK(IREF(JOBT,JDV,JF)+JOB)
              ENDDO
            ENDIF
            IREF(JOBT,JDV,JF)=IPT
            IPT = IPT+INOBS(JOBT,JF)
          ENDDO
        ENDDO
#endif
        IF (DIAG_RDOBS >= 1) THEN
          IPC = (IPT*100)/TNDVMAX
! diagnostic message directed to every PEn output stream
          WRITE(6,*)'RDOBS2:% of OBS array required=',IPC
        ENDIF

          if(mype == 0)then
! DEPENDS ON: file_close
          CALL FILE_CLOSE(IFILE,USED_FILES(JF),                         &
     &                    FILENAME_LEN(JF),ENVVAR,0,ICODE)
          endif

 2000   CONTINUE
!-----------------------------------------------------------------------

!         DETERMINE LIST OF AC OBS TYPES FOR MERGED FILES
!         NOBTYP IS NOW THE NO OF AC OBS TYPES IN THE MERGED LIST

!         INITIAL LIST IS LIST FOR FIRST FILE
          NOBTYP=INOBTYP(1)
          IF (NOBTYP >  0) THEN
            DO JOBT=1,NOBTYP
              KOBSTYP(JOBT) = IOBSTYP(JOBT,1)
            ENDDO
          ENDIF

        IF (NFILES >  1) THEN

!         GET FULL LIST FROM OTHER FILES
          DO JF=2,NFILES
          IF (.NOT.LEMPTY(JF)) THEN
          DO 2030 JTYP=1,INOBTYP(JF)
          DO 2010 JOBT=1,NOBTYP
            IF (IOBSTYP(JTYP,JF) == KOBSTYP(JOBT)) GO TO 2020
 2010     CONTINUE
          NOBTYP = NOBTYP+1
          IF (NOBTYP <= NOBTYPMX) THEN
            KOBSTYP(NOBTYP) = IOBSTYP(JTYP,JF)
          ELSE
            ICODE = 1
            CMESSAGE =                                                  &
     &      ' RDOBS2 : MAX NO OF AC OBS TYPES REACHED FOR MERGED FILES'
            if(mype == 0)PRINT *,                                       &
     &      ' RDOBS : MAX NO OF AC OBS TYPES REACHED FOR MERGED FILES'
            GO TO 9999
          ENDIF
 2020     CONTINUE
 2030     CONTINUE
          ENDIF
          ENDDO

        ENDIF

!       NOBTYP is now set to the no of obs types to be assimilated
!       which is the same as NACT.
        NOBTYP = NACT

        DO JOBT=1,NOBTYP
          OBSTYP(JOBT)   = -1
          NOBS  (JOBT)   =  0
          NDATAV(JOBT)   = -1
          NOBLEV(JOBT)   = -1
          OBLEVTYP(JOBT) = -1
          DO JLEV=1,MAXNLEV1
            DATALEVS(JLEV,JOBT) = MISSD
          ENDDO
        ENDDO

        IF (NOBTYP >  0) THEN

!         OBSTYP is now the list of obs types to be assimilated
!         which is the same as LACT.
          DO JACT=1,NACT
            OBSTYP(JACT) = LACT(JACT)
          ENDDO

          DO JOBT=1,NOBTYP
            DO JF=1,NFILES
              IF (.NOT.LEMPTY(JF)) THEN
              DO JTYP=1,INOBTYP(JF)
                IF (IOBSTYP(JTYP,JF) == OBSTYP(JOBT)) THEN

                  NDATAV(JOBT)   = INDATAV(JTYP,JF)
                  NOBLEV(JOBT)   = INOBLEV(JTYP,JF)
                  OBLEVTYP(JOBT) = IOBLVTYP(JTYP,JF)

                  DO JLEV=1,MAXNLEV1
                    DATALEVS(JLEV,JOBT) = PLEVELS(JLEV,JTYP,JF)
                  ENDDO
                  GO TO 2110

                ENDIF
              ENDDO
              ENDIF
            ENDDO
 2110       CONTINUE
          ENDDO

        ENDIF

!       Check that data is consistent for obs types being assimilated.

        IF (NFILES >  1 .AND. NOBTYP >  0) THEN

          DO JF=1,NFILES

!         ChecK NDVHDR

          IF (.NOT.LEMPTY(JF) .AND. INDVHDR(JF) /= NDVHDR) THEN
            ICODE = 1
            CMESSAGE =' RDOBS2 : MIS-MATCH IN AC OBS FILES - NDVHDR'
            if(mype == 0)then
            PRINT *,' '
            PRINT *,' RDOBS2 : DIFFERENT VALUES OF NDVHDR ?'
            PRINT '(A,5I5)',' NDVHDR =',(INDVHDR(IJF),IJF=1,NFILES)
            endif
            GO TO 9999
          ENDIF

!         Check MAXNLEV1

          IF (.NOT.LEMPTY(JF) .AND. IMAXNLEV(JF) /= MAXNLEV1) THEN
            ICODE = 1
            CMESSAGE =' RDOBS2 : MIS-MATCH IN AC OBS FILES - MAXNLEV1'
            if(mype == 0)then
            PRINT *,' '
            PRINT *,' RDOBS2 : DIFFERENT VALUES OF MAXNLEV1'
            PRINT '(A,5I5)',' MAXNLEV1 =',(IMAXNLEV(IJF),IJF=1,NFILES)
            endif
            GO TO 9999
          ENDIF

          ENDDO

          DO JOBT=1,NOBTYP
            DO JF=1,NFILES
            IF (INOBTYP(JF) >  0) THEN
              DO JTYP=1,INOBTYP(JF)

                IF (IOBSTYP(JTYP,JF) == OBSTYP(JOBT)) THEN

!               Check no of data values (NDATAV)

                IF (INDATAV(JTYP,JF) /= NDATAV(JOBT)) THEN
                  ICODE = 1
                  CMESSAGE =' RDOBS2 : MIS-MATCH IN OBS FILES - NDATAV'
                  if(mype == 0)then
                  PRINT *, ' RDOBS2 : Different No of Data Values.'
                  PRINT *, '        : See Obs Type ',OBSTYP(JOBT)
                  endif
                  GO TO 9999
                ENDIF

!               Check no of observation levels (NOBLEV)

                IF (INOBLEV(JTYP,JF) /= NOBLEV(JOBT)) THEN
                  ICODE = 1
                  CMESSAGE =' RDOBS2 : MIS-MATCH IN OBS FILES - NOBLEV'
                  if(mype == 0)then
                  PRINT *, ' RDOBS2 : Different No of Obs levels.'
                  PRINT *, '        : See Obs Type ',OBSTYP(JOBT)
                  endif
                  GO TO 9999
                ENDIF

!               Check Observation level type (OBLEVTYP)

                IF (IOBLVTYP(JTYP,JF) /= OBLEVTYP(JOBT)) THEN
                  ICODE = 1
                  CMESSAGE =                                            &
     &            ' RDOBS2 : MIS-MATCH IN OBS FILES - OBLEVTYP'
                  if(mype == 0)then
                  PRINT *,                                              &
     &            ' RDOBS2 : Different Observation Level Type'
                  PRINT *, '        : See Obs Type ',OBSTYP(JOBT)
                  endif
                  GO TO 9999
                ENDIF

!               Check pressure levels of observations (DATALEVS)

                DO JLEV=1,MAXNLEV1
                 IF (PLEVELS(JLEV,JTYP,JF) /= DATALEVS(JLEV,JOBT)) THEN
                   ICODE = 1
                   CMESSAGE =                                           &
     &             ' RDOBS2 : MIS-MATCH IN OBS FILES - DATALEVS'
                  if(mype == 0)then
                   PRINT *, ' RDOBS2 : Different levels in Obs files'
                   PRINT *, '        : See Obs Type ',OBSTYP(JOBT)
                   PRINT '(A,I5,A,2F10.1)', ' LEVEL',JLEV,              &
     &             ' ; Pressures =',PLEVELS(JLEV,JTYP,JF),              &
     &             DATALEVS(JLEV,JOBT)
                   endif
                   GO TO 9999
                 ENDIF
                ENDDO

                ENDIF

              ENDDO  !  End of loop over JTYP
            ENDIF
            ENDDO  !  End of loop over JF
          ENDDO  !  End of loop over JOBT

        ENDIF

!-----------------------------------------------------------------------
!L           SECTION 2: MERGE INPUT DATA & SET UP HEADER.
!-----------------------------------------------------------------------

!       Merge observation data and put into output buffer OBS
        TNDV=0
        DO JOBT=1,NOBTYP
          DO JDV=1,NDATAV(JOBT)
            DO JF=1,NFILES
              DO JTYP=1,INOBTYP(JF)
                IF (IOBSTYP(JTYP,JF) == OBSTYP(JOBT) .AND.              &
     &              INOBS(JTYP,JF)   >  0) THEN

                  IF (TNDV+INOBS(JTYP,JF) <= TNDVMAX) THEN

                    DO JOB=1,INOBS(JTYP,JF)
                      OBS(TNDV+JOB) = WORK(IREF(JTYP,JDV,JF)+JOB)
                    ENDDO
                    TNDV = TNDV+INOBS(JTYP,JF)

                  ELSE

                    ICODE = 1
                    CMESSAGE =                                          &
     &              ' RDOBS2 : Insufficient space in OBS array'
! important failure messages to PEn,OUT and OUT2 output streams
                    WRITE (0,*)                                         &
     &              ' RDOBS2 : Insufficient space in OBS array'
                    WRITE (0,*) ' Recode or rerun with fewer obs'
                    WRITE (6,*)                                         &
     &              ' RDOBS2 : Insufficient space in OBS array'
                    WRITE (6,*) ' Recode or rerun with fewer obs'
                    WRITE (7,*)                                         &
     &              ' RDOBS2 : Insufficient space in OBS array'
                    WRITE (7,*) ' Recode or rerun with fewer obs'
                    GO TO 9999

                  ENDIF
                ENDIF
              ENDDO   !  Loop over JTYP
            ENDDO   !  Loop over JF
          ENDDO   !  Loop over JDV
        ENDDO   !  Loop over JOBT
      WRITE(6,*)'temp check OBS ',TNDV,TNDVMAX,(TNDV*TNDVMAX)*100.0,'%'

!       Count total no of obs for each obs type
        DO JOBT=1,NOBTYP
          NOBS(JOBT)=0
          DO JF=1,NFILES
            DO JTYP=1,INOBTYP(JF)
            IF (IOBSTYP(JTYP,JF) == OBSTYP(JOBT)) THEN
              NOBS(JOBT) = NOBS(JOBT)+INOBS(JTYP,JF)
            ENDIF
            ENDDO
          ENDDO
        ENDDO

!       Get total no of obs (TNOBS)
        TNOBS = 0
        DO JOBT=1,NOBTYP
          TNOBS = TNOBS + NOBS(JOBT)
        ENDDO

!       Set up pointers to start of data for each obs type
        MDISPOBT(1)=0
        DO JOBT=2,NOBTYP
          MDISPOBT(JOBT) = MDISPOBT(JOBT-1)+NOBS(JOBT-1)*NDATAV(JOBT-1)
        ENDDO

!       Offset to first obs for each type
        OBS_NO_ST(1)=0
        DO JOBT=2,NOBTYP
          OBS_NO_ST(JOBT) = OBS_NO_ST(JOBT-1)+NOBS(JOBT-1)
        ENDDO


!       Check no of observations against maximum allowed (NOBSMAX)
        IF (TNOBS >  NOBSMAX) THEN
          ICODE = 1
          CMESSAGE = ' RDOBS2 : Insufficient space in OBS_FLAG array'
! check in every pe?
          WRITE(0,*)' RDOBS2 : Insufficient space in OBS_FLAG array'
          WRITE(0,*)' Recode or rerun with fewer obs'
          WRITE(0,'(A,I8)')' NO OF OBS   = ',TNOBS
          WRITE(6,*)' RDOBS2 : Insufficient space in OBS_FLAG array'
          WRITE(6,*)' Recode or rerun with fewer obs'
          WRITE(6,'(A,I8)')' NO OF OBS   = ',TNOBS
          WRITE(7,*)' RDOBS2 : Insufficient space in OBS_FLAG array'
          WRITE(7,*)' Recode or rerun with fewer obs'
          WRITE(7,'(A,I8)')' NO OF OBS   = ',TNOBS
          WRITE(0,'(A,I8)')' MAX ALLOWED = ',NOBSMAX
          WRITE(6,'(A,I8)')' MAX ALLOWED = ',NOBSMAX
          WRITE(7,'(A,I8)')' MAX ALLOWED = ',NOBSMAX
          GO TO 9999
         ELSE
          IF (DIAG_RDOBS >= 1) THEN
           IPC = (TNOBS*100)/NOBSMAX
! diagnostic message directed to every PEn output stream
          WRITE(6,*)'RDOBS2:% of OBS_FLAG array required=',IPC
         ENDIF
        ENDIF

!       Print summary of output file.
        IF (DIAG_RDOBS >  0) THEN
#if defined(MPP)
        DO JOBT=1,NOBTYP
        CountA(JOBT)=NOBS(JOBT)
        ENDDO

        CALL GC_SSYNC(NPROC,ISTAT)

        If(mype == 0)Then
          PRINT *, ' '
          PRINT '(A,I3.2,A,I2.2,A,I4)',                                 &
     &    ' REF DATE :',OBS_REF_DD,'/',OBS_REF_MM,'/',OBS_REF_YY
          PRINT '(A,I3.2,I2.2,A)',                                      &
     &    ' REF TIME :',OBS_REF_HH,OBS_REF_MIN,'Z'
          PRINT *, ' '
          PRINT '(A,F8.1,A)', ' REL TIME      :',TIMEREL,'M'
          PRINT '(A,F8.1,A)', ' REL T.W START :',TWSTART,'M'
          PRINT '(A,F8.1,A)', ' REL T.W END   :',TWEND,'M'
          PRINT *, ' '
          PRINT '(A,T30,18I6/(T30,18I6))',                              &
     &          ' AC OBS TYPES :',(OBSTYP(JOBT),JOBT=1,NOBTYP)
          PRINT *, ' '
          Do JOBT=1,NOBTYP
            CountC(JOBT)=0
          EndDo
        ENDIF

        CALL GC_SSYNC(NPROC,ISTAT)

        Do iproc=1,nproc-1
          IF(mype == 0) THEN
            CALL GC_IRECV(IPROC,NOBTYP,IPROC,ISTAT,COUNTB,COUNTA)
            Do JOB=1,5
              CountC(JOB)=CountC(JOB)+CountB(JOB)
            EndDo
          ELSEIF(mype == IPROC) THEN
            CALL GC_ISEND(IPROC,NOBTYP,0,ISTAT,COUNTB,COUNTA)
          ENDIF

          IF (DIAG_RDOBS >  1) THEN
          PRINT '(A,I3,T30,18I6/(T30,18I6))',                           &
     &        ' NO OF OBS    : pe=',iproc,(CountB(JOBT),JOBT=1,NOBTYP)
          ENDIF
          Do JOBT=1,NOBTYP
            CountC(JOBT)=CountC(JOBT)+CountB(JOBT)
          EndDo
        EndDo
        if(mype == 0) PRINT '(A,T30,18I6/(T30,18I6))',                  &
     &         ' NO OF OBS    : total:',(CountC(JOBT),JOBT=1,NOBTYP)
!        Get total no of obs (all pe's)
        If(NOBTYP >= 2)Then
          Do JOBT=2,NOBTYP
            CountC(1) = CountC(1) + CountC(JOBT)
          Enddo
        Endif  !NOBTYP >= 2
        IF(mype == 0) THEN
          PRINT *, ' '
          PRINT '(A,I8,A)', ' TOTAL NO OF OBS IN T.W :',CountC(1)
!        Any observations to assimilate ?
         If (CountC(1) == 0) Then
          PRINT *, ' '
          PRINT *, 'Timestep ',TIMESTEP_NO
          PRINT *, 'There are no observations to be assimilated.'
          WRITE (7,*) 'Timestep ',TIMESTEP_NO
          WRITE (7,*) 'There are no observations to be assimilated.'
         Endif !TNOBS == 0
        ENDIF

        CALL GC_SSYNC(NPROC,ISTAT)
#else
          PRINT *, ' '
          PRINT '(A,I3.2,A,I2.2,A,I4)',                                 &
     &    ' REF DATE :',OBS_REF_DD,'/',OBS_REF_MM,'/',OBS_REF_YY
          PRINT '(A,I3.2,I2.2,A)',                                      &
     &    ' REF TIME :',OBS_REF_HH,OBS_REF_MIN,'Z'
          PRINT *, ' '
          PRINT '(A,F8.1,A)', ' REL TIME      :',TIMEREL,'M'
          PRINT '(A,F8.1,A)', ' REL T.W START :',TWSTART,'M'
          PRINT '(A,F8.1,A)', ' REL T.W END   :',TWEND,'M'
          PRINT *, ' '
          PRINT '(A,T30,18I6/(T30,18I6))',                              &
     &          ' AC OBS TYPES :',(OBSTYP(JOBT),JOBT=1,NOBTYP)
          PRINT *, ' '
          PRINT '(A,T30,18I6/(T30,18I6))',                              &
     &          ' NO OF OBS    :',(NOBS(JOBT),JOBT=1,NOBTYP)
          PRINT *, ' '
          PRINT '(A,I8,A)', ' TOTAL NO OF OBS IN T.W :',TNOBS
!        Any observations to assimilate ?
         If (TNOBS == 0) Then
          PRINT *, ' '
          PRINT *, 'Timestep ',TIMESTEP_NO
          PRINT *, 'There are no observations to be assimilated.'
          WRITE (7,*) 'Timestep ',TIMESTEP_NO
          WRITE (7,*) 'There are no observations to be assimilated.'
         Endif !TNOBS == 0
#endif
        ENDIF

#if !defined(GLOBAL)
        DO JOBT=1,NOBTYP
          KTYPE=OBSTYP(JOBT)
          IF (KTYPE/100 == 3 .AND. NOBS(JOBT) >  0) THEN
!           ROTATE WIND COMPONENTS TO RESOLVE COMPONENTS ALONG ELF GRID
!           LAT/LON VALUES IN OBS ARRAY ARE FOR ELF GRID.
!           FIRST GET REAL LAT/LON OF WIND OBS IN WRKLAT/LON.
!           (ALL LAT/LON VALUES ARE IN DEGREES AT THIS STAGE)
            IP_LAT  = MDISPOBT(JOBT)
            IP_LONG = MDISPOBT(JOBT) + NOBS(JOBT)
! DEPENDS ON: eqtoll
            CALL EQTOLL ( OBS(IP_LAT+1),OBS(IP_LONG+1),                 &
     &                    WRKLAT,WRKLON,                                &
     &                    ELFPLAT,ELFPLON,NOBS(JOBT) )
!           NEXT CALCULATE COEFFICIENTS OF ROTATION.
! DEPENDS ON: w_coeff
            CALL W_COEFF( COEFF1,COEFF2,WRKLON,OBS(IP_LONG+1),          &
     &                    ELFPLAT,ELFPLON,NOBS(JOBT) )
!           NOW ROTATE WIND COMPONENTS.
            IF (KTYPE == 301) THEN
              IP_U = MDISPOBT(JOBT) + NOBS(JOBT)*NDVHDR
              IP_V = IP_U           + NOBS(JOBT)*NOBLEV(JOBT)
              DO JLEV=1,NOBLEV(JOBT)
! DEPENDS ON: w_lltoeq
                CALL W_LLTOEQ(COEFF1,COEFF2,                            &
     &          OBS(IP_U+1), OBS(IP_V+1), U_WRK, V_WRK,                 &
     &          NOBS(JOBT), NOBS(JOBT) )
                DO JOB=1,NOBS(JOBT)
                 OBS(IP_U+JOB)=U_WRK(JOB)
                 OBS(IP_V+JOB)=V_WRK(JOB)
                ENDDO
                IP_U = IP_U + NOBS(JOBT)
                IP_V = IP_V + NOBS(JOBT)
              ENDDO
            ELSEIF (KTYPE == 303) THEN
              IP_U = NDVHDR+2
              IP_V = NDVHDR+3
              IP_U = MDISPOBT(JOBT) + NOBS(JOBT)*(IP_U-1)
              IP_V = MDISPOBT(JOBT) + NOBS(JOBT)*(IP_V-1)
! DEPENDS ON: w_lltoeq
              CALL W_LLTOEQ ( COEFF1,COEFF2,                            &
     &             OBS(IP_U+1), OBS(IP_V+1), U_WRK, V_WRK,              &
     &             NOBS(JOBT), NOBS(JOBT) )
              DO JOB=1,NOBS(JOBT)
               OBS(IP_U+JOB)=U_WRK(JOB)
               OBS(IP_V+JOB)=V_WRK(JOB)
              ENDDO
            ELSEIF (KTYPE == 311) THEN
              IP_U = MDISPOBT(JOBT) + NOBS(JOBT)*NDVHDR
              IP_V = IP_U           + NOBS(JOBT)*NOBLEV(JOBT)
              DO JLEV=1,NOBLEV(JOBT)
! DEPENDS ON: w_lltoeq
                CALL W_LLTOEQ(COEFF1,COEFF2,                            &
     &             OBS(IP_U+1), OBS(IP_V+1), U_WRK, V_WRK,              &
     &             NOBS(JOBT), NOBS(JOBT) )
                DO JOB=1,NOBS(JOBT)
                 OBS(IP_U+JOB)=U_WRK(JOB)
                 OBS(IP_V+JOB)=V_WRK(JOB)
                ENDDO
                IP_U = IP_U + NOBS(JOBT)
                IP_V = IP_V + NOBS(JOBT)
              ENDDO
            ELSEIF (KTYPE == 302 .OR. KTYPE == 304 .OR.                 &
     &              KTYPE == 305 .OR. KTYPE == 306 ) THEN
              IP_U = MDISPOBT(JOBT) + NOBS(JOBT)*NDVHDR
              IP_V = IP_U           + NOBS(JOBT)
! DEPENDS ON: w_lltoeq
              CALL W_LLTOEQ ( COEFF1,COEFF2,                            &
     &             OBS(IP_U+1), OBS(IP_V+1), U_WRK, V_WRK,              &
     &             NOBS(JOBT), NOBS(JOBT) )
              DO JOB=1,NOBS(JOBT)
               OBS(IP_U+JOB)=U_WRK(JOB)
               OBS(IP_V+JOB)=V_WRK(JOB)
              ENDDO
            ELSE
            ICODE=1
            CMESSAGE='RDOBS2: a wind type remains unrotated'
            GOTO 9999
            ENDIF
          ENDIF
        ENDDO
#endif

        IF (NOBTYP >  0) THEN

          DO JOBT=1,NOBTYP

            IF (NOBS(JOBT) >  0) THEN

!             Convert Obs Latitudes (degrees) to Co-latitudes (radians)
              IP_LAT = MDISPOBT(JOBT)
              DO JOB=1,NOBS(JOBT)
                OBS(IP_LAT+JOB) = (90.0-OBS(IP_LAT+JOB))*PI_OVER_180
              ENDDO

!             Convert Obs Longitudes to Radians in range 0 to 2*PI.
              IP_LONG = MDISPOBT(JOBT)+NOBS(JOBT)
              DO JOB=1,NOBS(JOBT)
                IF (OBS(IP_LONG+JOB)  <   0.0)                          &
     &          OBS(IP_LONG+JOB) = OBS(IP_LONG+JOB)+360.0
                OBS(IP_LONG+JOB) = OBS(IP_LONG+JOB)*PI_OVER_180
              ENDDO

            ENDIF

          ENDDO

!L        SET UP OBLEVELS AND OBLAYERB (LEVELS IN PASCALS)
!         ------------------------------------------------
          DO JOBT=1,NOBTYP

          IF (OBLEVTYP(JOBT) == 3 .OR. OBLEVTYP(JOBT) == 4) THEN

          IF (OBLEVTYP(JOBT) == 3) THEN

            DO JLEV=1,NOBLEV(JOBT)
              OBLEVELS(JLEV,JOBT) = DATALEVS(JLEV,JOBT)
            ENDDO

            DO JLEV=2,NOBLEV(JOBT)
              OBLAYERB(JLEV,JOBT) =                                     &
     &        SQRT ( OBLEVELS(JLEV-1,JOBT) * OBLEVELS(JLEV,JOBT) )
            ENDDO

            OBLAYERB(1,JOBT) =  OBLEVELS(1,JOBT) *                      &
     &       OBLEVELS(1,JOBT) / OBLAYERB(2,JOBT)

            OBLAYERB(NOBLEV(JOBT)+1,JOBT) =                             &
     &             OBLEVELS(NOBLEV(JOBT),JOBT) *                        &
     &             OBLEVELS(NOBLEV(JOBT),JOBT) /                        &
     &             OBLAYERB(NOBLEV(JOBT),JOBT)

          ENDIF

          IF (OBLEVTYP(JOBT) == 4) THEN

            DO JLEV=1,NOBLEV(JOBT)+1
              OBLAYERB(JLEV,JOBT) = DATALEVS(JLEV,JOBT)
            ENDDO

            DO JLEV=1,NOBLEV(JOBT)
              OBLEVELS(JLEV,JOBT) =                                     &
     &        SQRT ( OBLAYERB(JLEV,JOBT) * OBLAYERB(JLEV+1,JOBT) )
            ENDDO

          ENDIF

          if(mype == 0)then
          PRINT *, ' '
          PRINT *, ' Observation Type ',OBSTYP(JOBT)
          PRINT *, ' '
          PRINT '(A,/,(1X,5F8.1))', '  Levels (mb) =',                  &
     &           (OBLEVELS(JLEV,JOBT)*0.01,JLEV=1,NOBLEV(JOBT))
          PRINT '(A,/,(1X,5F8.1))', '  Layer boundaries (mb) =',        &
     &           (OBLAYERB(JLEV,JOBT)*0.01,JLEV=1,NOBLEV(JOBT)+1)

          endif
          ENDIF

          ENDDO

      ENDIF
!     ====================================================
!     SET UP ARRAY NERLEV1 WHICH POINTS THE FIRST DATA VALUE
!     CORRESPONDING TO THE FIRST LEVEL OF ERROR RATIO FOR THE OBS TYPE

      IF (NOBTYP >  0) THEN
        DO JOBT=1,NOBTYP
          NERLEV1(JOBT) = NDATAV(JOBT)-NOBLEV(JOBT)+1
        ENDDO
      ENDIF
!     ====================================================
      IF (NACT >  0) THEN

        if(mype == 0)then
        PRINT   10, (LACT(J),J=1,NACT)
 10     FORMAT('0TYPES TO BE PROCESSED : LACT    =',15I5)
        PRINT   12, (NOBLEV(J),J=1,NACT)
 12     FORMAT('          NO OF LEVELS : NOBLEV  =',15I5)
        PRINT   13, (NDATAV(J),J=1,NACT)
 13     FORMAT('  NO OF DATA VARIABLES : NDATAV  =',15I5)
        PRINT   14, (NERLEV1(J),J=1,NACT)
 14     FORMAT('       FIRST ERROR LEV : NERLEV1 =',15I5)
        endif

      ENDIF
!     ====================================================
!     CALL SETDAC TO SET UP FOR ANY DIAGNOSTICS REQUIRED
! DEPENDS ON: setdac
      CALL SETDAC (OBS,TNDV)
!     ====================================================

!     RE-USE IWORK IN THIS SECTION

      IF (TNOBS >  0) THEN

        DO JOB=1,TNOBS
          OBS_FLAG(JOB)=0
        ENDDO

!       Convert Analysis Types into INTEGER.

        IP_TYPE = 4
        DO JOBT=1,NOBTYP
          IF (NOBS(JOBT) >  0) THEN
            IP_MOT = MDISPOBT(JOBT)+(IP_TYPE-1)*NOBS(JOBT)
            DO JOB=1,NOBS(JOBT)
              IWORK(OBS_NO_ST(JOBT)+JOB) = NINT( OBS(IP_MOT+JOB) )
            ENDDO
          ENDIF
        ENDDO

        ITOT0=0
        ITOT1=0
        DO JMOT=1,NANALTYP
          IF (IOMITOBS(JMOT) >  0) THEN
            DO JOB=1,TNOBS
              IF (IWORK(JOB) == IOMITOBS(JMOT)) OBS_FLAG(JOB)=1
            ENDDO
            IF (DIAG_RDOBS >= 1) THEN
              ITOT1=0
              DO JOB=1,TNOBS
                IF (OBS_FLAG(JOB) == 1) ITOT1=ITOT1+1
              ENDDO
              ITOTAL=ITOT1-ITOT0
! in every pe?
              WRITE(6,'(A,I5,A,I6)')' ANALYSIS TYPE ',IOMITOBS(JMOT),   &
     &        ' : NO OF OBSERVATIONS FLAGGED (DO NOT USE)  -',ITOTAL
              ITOT0=ITOT1
            ENDIF
          ENDIF
        ENDDO


!       OPTIONAL PRINT OUT.
        IF (DIAG_RDOBS >= 1 .AND. ITOT1 /= 0) THEN
! in every pe?
          WRITE(6, '(A,I6)')                                            &
     &    ' TOTAL NO OF OBSERVATIONS FLAGGED (DO NOT USE)  -',ITOT1

        ENDIF

      ENDIF
!-----------------------------------------------------------------------

#if !defined(MPP)
      REWIND IUNITNO
      WRITE (IUNITNO) TNOBS,TNDV
      IF (TNOBS >  0 .AND. TNDV >  0) THEN

!       NO OBSERVATIONS TO BE SAVED

        WRITE (IUNITNO) (OBS_FLAG(J),J=1,TNOBS),                        &
     &                  (OBS(J),     J=1,TNDV)
        IF (DIAG_RDOBS >  1)                                            &
     &  WRITE (6,*) ' OBS and OBS_FLAG written to unit no ',IUNITNO

      ENDIF
#endif

 9999 CONTINUE
      RETURN
      END SUBROUTINE RDOBS2
#endif
