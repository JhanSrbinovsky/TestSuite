#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Read the basis file information
! Subroutine Interface:

      SUBROUTINE RDBASIS                                                &
     & (IU,CMESSAGE,ErrorStatus)
      IMPLICIT NONE

! Description:
!
! Method:
!
! Current code owner:  S.J.Swarbrick
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Mar. 95    Original code.  S.J.Swarbrick
!   4.0     Sept.95    Original code.  S.J.Swarbrick
!  4.0  18/10/95  Add ErrorStatus to GET_FILE call. RTHBarnes
!   4.1     June 96    DISCT_LEV function used to check for model
!                         levels       S.J.Swarbrick
!   4.4     Sep. 97    Add IOFF to namelist for offset to sampling
!                      frequency. S.D. Mullerworth
!   4.5     July 98    Remove call to INTRFACE (A Van der Wal)
!   4.5    30/10/97    Read stash data on PE 0 for the T3E
!                      and distribute it.
!                        Author: Bob Carruthers, Cray Research
!LL  5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!    6.0  11/09/03  Add file='empty' for IBM (provided by Zoe Chaplin)
!                                                           P.Dando
!    6.2  23/11/05  Removed DIAG190 define and diagnostic. T. Edwards
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Global variables:

#include "c_mdi.h"
#include "lenfil.h"
#include "csubmodl.h"
#include "version.h"
#include "cstash.h"
#include "parparm.h"
#include "typsize.h"
#include "model.h"
#include "stextend.h"
#include "cocnindx.h"

! Subroutine arguments

!   Scalar arguments with intent(in):
      INTEGER IU       ! Unit no. of stash basis file
      INTEGER IE
!   Scalar arguments with intent(out):
      CHARACTER*80 CMESSAGE      ! Error return message

!   Error status:
      INTEGER        ErrorStatus ! Error return code

! Local variables:
      LOGICAL MODEL_LEV  !TRUE for model levels
      INTEGER I,J,L
      INTEGER IDUM
      INTEGER IOSTAT
      INTEGER NtsRecs    !Counter for time series records
#if defined(MPP) && defined(T3E)
!
      integer shmem_n_pes, msg, info, nproc, shmem_my_pe, mype, k
#endif

!   Namelist STASHNUM
      INTEGER            NUM_REQ,NUM_DOM,NUM_TIM,NUM_USE
      NAMELIST/STASHNUM/ NUM_REQ,NUM_DOM,NUM_TIM,NUM_USE
#if defined(MPP) && defined(T3E)
      common/shmem_stashnum/ num_req,num_dom,num_tim,num_use
!dir$ cache_align /shmem_stashnum/
#endif

!   Namelist STREQ: STASH requests
      INTEGER         IMOD,ISEC,ITEM,IDOM,ITIM,IUSE
      NAMELIST/STREQ/ IMOD,ISEC,ITEM,IDOM,ITIM,IUSE
#if defined(MPP) && defined(T3E)
      common/shmem_streq/ imod,isec,item,idom,itim,iuse
!dir$ cache_align /shmem_streq/
#endif

!   Namelist TIME: Time profiles
      CHARACTER*8 NAME
      CHARACTER*2 UNT1,UNT2,UNT3
      INTEGER     ITYP,ISAM,INTV,IOPT
      INTEGER     ISTR,IEND,IFRE,IOFF,ITIMES,ISER(NTIMEP)
      NAMELIST/TIME/ITYP,ISAM,INTV,UNT1  ,UNT2,UNT3,IOPT                &
     &             ,ISTR,IEND,IFRE,IOFF,ITIMES,ISER,NAME

#if defined(MPP) && defined(T3E)
      common/shmem_time/ ityp,isam,intv,iopt                            &
     &             ,istr,iend,ifre,ioff,itimes,iser
!dir$ cache_align /shmem_time/
      common/shmem_time_c/ name, unt1, unt2, unt3
!dir$ cache_align /shmem_time_c/
#endif
!   Namelist DOMAIN: Domain profiles
      INTEGER     IOPL           !Level type code
      INTEGER     ILEVS          !Flag for range/selected model levels
      INTEGER     LEVB,LEVT      !Bottom/top levels for range
      INTEGER     PLT            !Pseudo level type code
      INTEGER     IOPA                !Horizontal domain type code
      INTEGER     INTH,ISTH,IEST,IWST !Horiz domain limits (IOPA=9,10)
      INTEGER     IMSK                !Grid type code
      INTEGER     IMN                 !Spatial meaning code
      INTEGER     IWT                 !Weighting code
      INTEGER     LEVLST (NLEVP)      !Levels lists array: integer
      REAL        RLEVLST(NLEVP)      !                    real
      INTEGER     PSLIST (NPSLEVP)    !                    pseudo
      CHARACTER*1 TS                    !Flag for time series domain
      INTEGER     TSNUM                 !No. of time ser doms in prof
      INTEGER     TBLIM (tsdp),TTLIM (tsdp) !TS dom limits (top/bottom)
      REAL        TBLIMR(tsdp),TTLIMR(tsdp) !ditto for real levels
      INTEGER     TNLIM (tsdp),TSLIM (tsdp) !TS dom limits (N/S)
      INTEGER     TWLIM (tsdp),TELIM (tsdp) !TS dom limits (E/W)

      NAMELIST/DOMAIN/IOPL ,ILEVS ,LEVB  ,LEVT  ,PLT    ,IOPA ,IMSK ,   &
     &                IMN  ,IWT   ,TS    ,LEVLST,RLEVLST,NAME ,         &
     &                INTH ,ISTH  ,IEST  ,IWST  ,PSLIST ,               &
     &                TSNUM,TBLIM ,TTLIM ,TNLIM ,TSLIM  ,TELIM,TWLIM,   &
     &                      TBLIMR,TTLIMR
#if defined(MPP) && defined(T3E)

      common/shmem_domain/ IOPL ,ILEVS ,LEVB  ,LEVT  ,PLT    ,IOPA ,    &
     &              IMSK ,                                              &
     &              IMN  ,IWT   ,LEVLST,RLEVLST,                        &
     &              INTH ,ISTH  ,IEST  ,IWST  ,PSLIST ,                 &
     &              TSNUM,TBLIM ,TTLIM ,TNLIM ,TSLIM  ,TELIM,TWLIM,     &
     &                    TBLIMR,TTLIMR
!dir$ cache_align /shmem_domain/
      common/shmem_domain_c/ ts
!dir$ cache_align /shmem_domain_c/
#endif

!   Namelist USE: Useage profiles
      INTEGER      LOCN,IUNT
      NAMELIST/USE/LOCN,IUNT,NAME
#if defined(MPP) && defined(T3E)
      common/shmem_use/locn, iunt
!dir$ cache_align /shmem_use/
#endif

! Function and subroutine calls:
      LOGICAL  DISCT_LEV
      EXTERNAL GET_FILE

!- End of Header ------------------------------------------------------

! Initialisation
      NDIAG   =0
      NTPROF  =0
      NDPROF  =0
      NUPROF  =0
      NUM_REQ =0
      NUM_TIM =0
      NUM_DOM =0
      NUM_USE =0
      NtsRecs =0

      DO I = 1,NDIAGPM
        MODL_B(I)=0
        ISEC_B(I)=0
        ITEM_B(I)=0
        ITIM_B(I)=0
        IDOM_B(I)=0
        IUSE_B(I)=0
      END DO

        ITIMES   =0
      DO I=1,NPROFTP
        TIMPRO(I)='        '
        ITYP_T(I)=0
        INTV_T(I)=0
        UNT1_T(I)='  '
        ISAM_T(I)=0
        UNT2_T(I)='  '
        IOPT_T(I)=0
        ISTR_T(I)=0
        IEND_T(I)=0
        IFRE_T(I)=0
        IOFF_T(I)=0
        UNT3_T(I)='  '
        ITIM_T(I)=0
        MODL_T(I)=0
        DO J = 1,NTIMEP
          ISER_T(J,I)=0
        END DO
      END DO

      DO I=1,NPROFDP
        DOMPRO  (I)='        '
        IOPL_D  (I)=0
        LEVB_D  (I)=0
        LEVT_D  (I)=0
        PLT_D   (I)=0
        IOPA_D  (I)=0
        INTH_D  (I)=0
        ISTH_D  (I)=0
        IEST_D  (I)=0
        IWST_D  (I)=0
        IMSK_D  (I)=0
        IMN_D   (I)=0
        IWT_D   (I)=0
        TS_D    (I)=' '
        PLLEN_D (I)=0
        PLPOS_D (I)=0
        ILEV_D  (I)=0
      END DO
      DO I = 1,NLEVP
      DO J = 1,NPROFDP
        LEVLST_D (I,J)=0
        RLEVLST_D(I,J)=0
      END DO
      END DO
      DO I = 1,NPSLEVP
      DO J = 1,NPSLISTP
        PSLIST_D(I,J)=0
      END DO
      END DO

      DO I = 1,NPROFUP
        USEPRO(I)='        '
        LOCN_U(I)=0
        IUNT_U(I)=0
      END DO

#if defined(MPP) && defined(T3E)
      mype=shmem_my_pe()
      nproc=shmem_n_pes()
      if(mype == 0) then
#endif
#if defined(IBM)
      file='empty '
#endif
      CALL GET_FILE(IU,FILE,80,ErrorStatus)   ! Get name for stash file

! Rewind stash control file
      REWIND(IU)
#if defined(MPP) && defined(T3E)
      endif
#endif

!STASH control file header namelist
#if defined(MPP) && defined(T3E)
  99  continue
      if(mype == 0) READ(IU,STASHNUM)
!
      msg=7010
      info=0
      call gc_ibcast(msg, 4, 0, nproc, info, num_req)
!
#else
  99  READ(IU,STASHNUM)
#endif
      IF (NUM_REQ == -1) GOTO 999
      NDIAG =NDIAG +NUM_REQ
      NTPROF=NTPROF+NUM_TIM
      NDPROF=NDPROF+NUM_DOM
      NUPROF=NUPROF+NUM_USE
      IF (NDIAG  >  NRECDP ) THEN
        WRITE(6,*) 'NUMBER OF DIAGNOSTIC REQUESTS EXCEEDS LIMIT OF ',   &
     &  NRECDP ,' SOME HAVE BEEN IGNORED'
        GO TO 999
      END IF
      IF (NDPROF >  NPROFDP) THEN
        WRITE(6,*) 'ERROR IN STASHC:'
        WRITE(6,*) 'NUMBER OF DOMAIN PROFILES EXCEEDS LIMIT OF ',       &
     &  NPROFDP
        WRITE(6,*) 'ARRAYS WILL BE OVERWRITTEN'
        ErrorStatus=1
        GO TO 9999
      END IF
      IF (NUPROF >  NPROFUP) THEN
        WRITE(6,*) 'ERROR IN STASHC:'
        WRITE(6,*) 'NUMBER OF USEAGE PROFILES EXCEEDS LIMIT OF ',       &
     &  NPROFUP
        WRITE(6,*) 'ARRAYS WILL BE OVERWRITTEN'
        ErrorStatus=1
        GO TO 9999
      END IF
      IF (NTPROF >  NPROFTP) THEN
        WRITE(6,*) 'ERROR IN STASHC:'
        WRITE(6,*) 'NUMBER OF TIME PROFILES EXCEEDS LIMIT OF ',         &
     &  NPROFTP
        WRITE(6,*) 'ARRAYS WILL BE OVERWRITTEN'
        ErrorStatus=1
        GO TO 9999
      END IF

!STASH request namelists
      IF (NUM_REQ >  0) THEN
      DO I = NDIAG-NUM_REQ+1,NDIAG
        IMOD=0
        ISEC=0
        ITEM=0
        ITIM=0
        IDOM=0
        IUSE=0
#if defined(MPP) && defined(T3E)
        if(mype == 0) READ(IU,STREQ)
!
        msg=7011
        info=0
        call gc_ibcast(msg, 6, 0, nproc, info, imod)
!
#else
        READ(IU,STREQ)
#endif
        MODL_B(I)=IMOD
        ISEC_B(I)=ISEC
        ITEM_B(I)=ITEM
        IF (ITIM /= 0) THEN
          ITIM     =ITIM+NTPROF-NUM_TIM
        END IF
        IDOM     =IDOM+NDPROF-NUM_DOM
        IUSE     =IUSE+NUPROF-NUM_USE
        ITIM_B(I)=ITIM
        IDOM_B(I)=IDOM
        IUSE_B(I)=IUSE
      END DO
      END IF

!Time profile namelists
      IF (NUM_TIM >  0) THEN
      DO I = NTPROF-NUM_TIM+1,NTPROF
        NAME='        '
        ISAM=0
        INTV=0
        IOPT=0
        ISTR=0
        IEND=0
        IFRE=0
        IOFF=0
        DO J = 1,NTIMEP
          ISER(J)=0
        END DO
!  Read namelist
#if defined(MPP) && defined(T3E)
        if(mype == 0) READ(IU,TIME)
!
        msg=7011
        info=0
        call gc_ibcast(msg, 9+ntimep, 0, nproc, info, ityp)
        msg=7012
        info=0
        call gc_cbcast(msg, 14, 0, nproc, info, name)
!
#else
        READ(IU,TIME)
#endif
        TIMPRO(I)=NAME
        ITYP_T(I)=ITYP
        IF (ITYP /= 1) THEN
!  Diagnostic output is time-processed
          ISAM_T(I)=ISAM  !Sampling frequency
          UNT2_T(I)=UNT2
          INTV_T(I)=INTV  !Processing interval
          UNT1_T(I)=UNT1
        END IF
!  Diag. output time option
        IOPT_T(I)=IOPT
        UNT3_T(I)=UNT3
        IF      (IOPT == 1) THEN
!    Regular output time interval
          ISTR_T(I)=ISTR
          IEND_T(I)=IEND
          IFRE_T(I)=IFRE
          IOFF_T(I)=IOFF
        ELSE IF (IOPT == 2) THEN
!    Specified list of output times
!     Length of times table
          ITIM_T(I)=ITIMES
!     Internal model label for times table
          IF (NDIAG >  0) THEN
            MODL_T(I)=MODL_B(NDIAG)
          END IF
!     Times table
          DO J = 1,ITIMES
            ISER_T(J,I)=ISER(J)
          END DO
        END IF
      END DO
      END IF

!Domain profile namelists
      IF (NUM_DOM >  0) THEN
      DO I = NDPROF-NUM_DOM+1,NDPROF
!Initialise
        NAME ='        '
        IOPL =0
        LEVB =0
        LEVT =0
        ILEVS=0
          LEVLST (1)= IMDI
         RLEVLST (1)= RMDI
        DO J = 2,NLEVP
          LEVLST (J)= 0
          RLEVLST(J)= 0.0
        END DO
        DO J = 1,NPSLEVP
          PSLIST(J)=0
        END DO
        DO J = 1,tsdp
          TBLIM (J)=0
          TTLIM (J)=0
          TBLIMR(J)=0.
          TTLIMR(J)=0.
          TNLIM (J)=0
          TSLIM (J)=0
          TELIM (J)=0
          TWLIM (J)=0
        END DO
!Read namelist
#if defined(MPP) && defined(T3E)
        if(mype == 0) READ(IU,DOMAIN)
!
        msg=7013
        info=0
        j=loc(ttlimr(tsdp))
        k=loc(iopl)
        j=(j-k+8)/8
        if(mype == 0) write(6,*) 'Length of DOMAIN is ',j,              &
     &   ' Words'
        call gc_ibcast(msg, j, 0, nproc, info, iopl)
        msg=7014
        info=0
        call gc_cbcast(msg, 1, 0, nproc, info, ts)
        msg=7015
        info=0
        call gc_cbcast(msg, 8, 0, nproc, info, name)
!
#else
        READ(IU,DOMAIN)
#endif
!Check for errors in levels lists
! DEPENDS ON: disct_lev
        MODEL_LEV=DISCT_LEV(IOPL,ErrorStatus,CMESSAGE)
        IF (MODEL_LEV) THEN
! Model levels
          IF (ILEVS == 2) THEN
            IF ( LEVLST(1) == IMDI ) THEN
              WRITE(6,*)                                                &
     &       'ERROR,RDBASIS: LEVELS LIST IN DOMAIN PROFILE '            &
     &       ,I,' HAS NO ENTRIES'
              CMESSAGE='ERROR,RDBASIS: LEVELS LIST HAS NO ENTRIES'
              ErrorStatus=1
              GO TO 9999
            END IF
          END IF
        ELSE IF (IOPL /= 5) THEN
          IF (RLEVLST(1) == RMDI) THEN
            WRITE(6,*)                                                  &
     &     'ERROR,RDBASIS: LEVELS LIST IN DOMAIN PROFILE '              &
     &     ,I,' HAS NO ENTRIES'
            CMESSAGE='ERROR,RDBASIS: LEVELS LIST HAS NO ENTRIES'
            ErrorStatus=1
            GO TO 9999
          END IF
        END IF
!Profile name
        DOMPRO(I)=NAME
!Store level type code in IOPL_D
        IOPL_D(I)=IOPL
! DEPENDS ON: disct_lev
        MODEL_LEV=DISCT_LEV(IOPL,ErrorStatus,CMESSAGE)
        IF (MODEL_LEV) THEN
!Integer levels
          ILEV_D(I)=ILEVS
          IF (ILEVS == 1) THEN
!  Range of model levels
            LEVB_D(I)=LEVB
            LEVT_D(I)=LEVT
          END IF
          IF (ILEVS == 2) THEN
!  List of selected model levels
            LEVB_D(I)=-1
            DO J=1,NLEVP
              LEVLST_D(J,I) = LEVLST(J)
              IF (LEVLST(J) >  0) THEN
!  Store no. of levels in LEVT_D(I)
                LEVT_D(I)=LEVT_D(I)+1
              END IF
            END DO
          END IF
        ELSE IF (IOPL /= 5) THEN
!Real levels
          LEVB_D(I)=-1
          DO J=1,NLEVP
            RLEVLST_D(J,I) = RLEVLST(J)
              IF (RLEVLST(J) >  0.0) THEN
!  Store no. of levels in LEVT_D(I)
                LEVT_D(I)=LEVT_D(I)+1
              END IF
          END DO
        END IF
!Store pseudo level type code in PLT_D
        PLT_D (I)=PLT
        IF (PLT >  0) THEN
!Domain profile 'I' has pseudo levels list
!    Count total no. of pseudo levels lists
        NPSLISTS = NPSLISTS + 1
!    Store list in column 'NPSLISTS' of PSLIST_D
          DO J=1,NPSLEVP
            PSLIST_D(J,NPSLISTS) = PSLIST (J)
!    PPLEN_D(I) stores length of ps.lev.list for domain 'I'
            IF (PSLIST(J) >  0) THEN
              PLLEN_D (I)        = PLLEN_D(I) + 1
            END IF
          END DO
!    PLPOS(I) stores the column no. in PSLIST_D for dom. prof. 'I'
          PLPOS_D(I) = NPSLISTS
        END IF
!Store horizontal domain type in IOPA_D
        IOPA_D(I)=IOPA
        IF (IOPA == 9.OR.IOPA == 10) THEN
!    Specified area
          INTH_D(I)=INTH
          ISTH_D(I)=ISTH
          IEST_D(I)=IEST
          IWST_D(I)=IWST
        END IF
        IMSK_D(I)=IMSK ! Gridpoint option
        IMN_D (I)=IMN  ! Meaning option
        IWT_D (I)=IWT  ! Weighting option
        TS_D  (I)=TS   ! Time series domain
        IF (TS_D(I)  ==  'Y') THEN
!This domain profile has a block of time series domains
!  Store time series data for current profile in _TS arrays
          NSERIES    = NSERIES+1        ! Time series block number:
          NPOS_TS(I) = NSERIES          !      used as a pointer
          NRECS_TS(NSERIES) = TSNUM     ! No. of records in ts block
          NSERREC_S  = NSERREC_S+TSNUM  ! Cumulative total ts records
          IF (NSERREC_S <= NTimSerP) THEN
            DO J = 1,TSNUM
              IF (J <= tsdp) THEN
                NtsRecs = NtsRecs+1
! DEPENDS ON: disct_lev
                MODEL_LEV=DISCT_LEV(IOPL,ErrorStatus,CMESSAGE)
                IF (MODEL_LEV) THEN
                  BLIM_TS (NtsRecs)= TBLIM (J)
                  TLIM_TS (NtsRecs)= TTLIM (J)
                ELSE IF (IOPL /= 5) THEN
                  BLIMR_TS(NtsRecs)= TBLIMR(J)
                  TLIMR_TS(NtsRecs)= TTLIMR(J)
                END IF
                NLIM_TS(NtsRecs) = TNLIM(J)
                SLIM_TS(NtsRecs) = TSLIM(J)
                ELIM_TS(NtsRecs) = TELIM(J)
                WLIM_TS(NtsRecs) = TWLIM(J)
              ELSE
                WRITE(6,*)                                              &
     &         'MESSAGE FROM ROUTINE RDBASIS: ',                        &
     &         'no. of time series in domain profile ',I,               &
     &         ' exceeds allowed limit of ',tsdp,' some will be',       &
     &         ' ignored'
              END IF
            END DO
          ELSE
            WRITE(6,*)                                                  &
     &     'TIMSER: total no. of time series requested exceeds ',       &
     &     'allowed limit of ',NTimSerP,'; some will be ignored.'
          END IF
        ELSE
          NPOS_TS  (I)=-1
        END IF
      END DO
      END IF
      NSERBLK_S = NSERIES

!Useage profile namelists
      IF (NUM_USE >  0) THEN
        NAME='        '
        LOCN=0
        IUNT=0
      DO I = NUPROF-NUM_USE+1,NUPROF
#if defined(MPP) && defined(T3E)
        if(mype == 0) READ(IU,USE)
!
        msg=7016
        info=0
        call gc_ibcast(msg, 2, 0, nproc, info, locn)
        msg=7017
        info=0
        call gc_cbcast(msg, 8, 0, nproc, info, name)
!
#else
        READ(IU,USE)
#endif
        USEPRO(I)=NAME
        LOCN_U(I)=LOCN
        IUNT_U(I)=IUNT
      END DO
      END IF
      GO TO 99

 999  CONTINUE

!Initialise model config. arrays before reading STSHCOMP
      DO I = 0,NSECTP
        ATMOS_SR(I)='  '
        OCEAN_SR(I)='  '
         SLAB_SR(I)='  '
         WAVE_SR(I)='  '
        INDEP_SR(I)='  '
      END DO
#if defined(MPP) && defined(T3E)
      if(mype == 0) then
#endif
      READ(5,STSHCOMP)
#if defined(MPP) && defined(T3E)
      endif
!
      msg=7018
      info=0
      j=loc(oasfldid(4))
      k=loc(run_target_end)
      j=(j-k+8)/8
      if(mype == 0) write(6,*) 'Length of STSHCOMM is ',j,              &
     & ' Words'
      call gc_ibcast(msg, j, 0, nproc, info, run_target_end)
      msg=7019
      info=0
      j=loc(wave_sr(nsectp))
      k=loc(bspmsl)
      j=j-k+1
      if(mype == 0) write(6,*) 'Length of STSHCHAR is ',j,              &
     & ' Bytes'
      call gc_cbcast(msg, j, 0, nproc, info, bspmsl)
!
#endif

#if defined(MPP) && defined(OCEAN)
        NROWSO = JFIN - JST + 1 + (2*O_NS_HALO)
#endif
#if defined(MPP) && defined(T3E)
      if(mype == 0) then
#endif
      CLOSE(UNIT=IU,STATUS='KEEP',IOSTAT=IOSTAT)
#if defined(MPP) && defined(T3E)
      endif
#endif

 9999 RETURN
      END SUBROUTINE RDBASIS

!- End of subroutine code ------------------------------------------
#endif
