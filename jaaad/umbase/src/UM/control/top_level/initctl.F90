#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE INITCTL-------------------------------------------------
!LL
!LL  PROGRAMMING STANDARD: UNIFIED MODEL DP NO. 3, VERSION 3
!LL
!LL  SYSTEM TASK: C4
!LL
!LL  SYSTEM COMPONENTS: C30, C40
!LL
!LL  PURPOSE:   Initialises STASH control arrays from STASH control
!LL            file.
!LL
!LL  EXTERNAL DOCUMENTATION: UMDP NO. C4 VERSION NO. 4
!LL
!LLEND-------------------------------------------------------------

      SUBROUTINE INITCTL(                                               &
#include "argsts.h"
#include "argppx.h"
#include "argd1.h"
     &                  ICODE,CMESSAGE)

      IMPLICIT NONE

! Defines N_INTERNAL_MODEL to dimension STASH arrays
#include "csubmodl.h"

!L Arguments

#include "parparm.h"
#include "typsize.h"
! Contains *CALL CPPXREF
#include "typsts.h"
! Contains *CALL VERSION
#include "ppxlook.h"


      INTEGER                                                           &
     &    ICODE                  ! OUT: Error return code
!
      CHARACTER*256                                                     &
     &    CMESSAGE               ! OUT: Error return message

#include "chsunits.h"
#include "ccontrol.h"
#include "clookadd.h"
#include "chistory.h"
#include "stparam.h"
#include "c_mdi.h"
#include "cstash.h"
! Declares arrays used in STASH_PROC code (LIST_S etc.);
#include "stextend.h"
                !   also contains common block STEXTEND
! For accessing D1 addressing array
#include "typd1.h"
! Print status information in CPRINTST:
#include "cprintst.h"

! External subroutines called

      INTEGER      EXPPXI
      CHARACTER*36 EXPPXC

!  Local arrays

!  STASH input lengths
      INTEGER SI_LEN(NITEMS,0:NSECTS,N_INTERNAL_MODEL_MAX)
      INTEGER ppxref_dat(PPXREF_CODELEN)

! Local variables

      CHARACTER*36 NAME
      REAL                                                              &
     &       REAL_LEVELS(NUM_STASH_LEVELS,NUM_LEVEL_LISTS)

      INTEGER                                                           &
     &        NUM_LISTS,                                                &
     &        NUM_LEVELS,                                               &
     &        NUM_PSEUDO_LEVELS,                                        &
     &        N_TABLES,                                                 &
     &        I,                                                        &
     &        IPK,                                                      &
     &        KK,                                                       &
     &        L,                                                        &
     &        IS,                                                       &
     &        IE,                                                       &
     &        II,                                                       &
     &        SM,                                                       &
     &        IOBJ,                                                     &
     &        ISEC,                                                     &
     &        ITM,                                                      &
     &        Im_ident,                                                 &
     &        Sm_ident,                                                 &
     &        IX,                                                       &
     &        ISTEP,                                                    &
     &        IL,                                                       &
     &        IM,                                                       &
     &        JJ,                                                       &
     &        INPUT_LENGTH

      INTEGER IR1               ! loop start
      INTEGER IR,J,K            ! loop count
      CHARACTER*1  VAR_TYPE
      CHARACTER*80 FILENAME
      LOGICAL INT_MOD_INCLUDED  ! Flag to indicate whether a particular
                                !   internal model is included

!-----------------------------------------------------------------------

!  1. Assign STASHlist and associated lists to appropriate UM arrays

!     Initialise STLIST to zero

      DO   II = 1,TOTITEMS
        DO IE = 1,LEN_STLIST
          STLIST(IE,II)=0
        END DO
      END DO

!     Assign STASH list to STLIST

      DO   I = 1,TOTITEMS
        DO J = 1,LEN_STLIST
          STLIST(J,I) = LIST_S(J,I)
        END DO
      END DO

!     Assign STASH times tables to STTABL

      DO   I = 1,NSTTABL
        DO J = 1,NSTTIMS
          STTABL(J,I) = ITIM_S(J,I)
        END DO
      END DO

!     Assign STASH levels lists to STASH_LEVELS

      DO   II=1,NUM_LEVEL_LISTS      ! Initialise STASH_LEVELS to -99
        DO JJ=1,NUM_STASH_LEVELS+1
          STASH_LEVELS(JJ,II)=-99
        END DO
      END DO

      DO   I = 1,NUM_LEVEL_LISTS
        NUM_LEVELS        = LEVLST_S(1,I)
        STASH_LEVELS(1,I) = NUM_LEVELS
        DO J = 1,NUM_LEVELS
          IF      (LLISTTY(I) == 'R') THEN
            REAL_LEVELS (J  ,I) = RLEVLST_S  (J+1,I)
            STASH_LEVELS(J+1,I) =(REAL_LEVELS(J  ,I)+0.0001)*1000.0
          ELSE IF (LLISTTY(I) == 'I') THEN
            STASH_LEVELS(J+1,I) = LEVLST_S   (J+1,I)
          END IF
        END DO
      END DO

!     Store STASH pseudo levels lists in STASH_PSEUDO_LEVELS

      DO   II=1,NUM_PSEUDO_LISTS       ! Initialise STASH_PSEUDO_LEVELS
        DO JJ=1,NUM_STASH_PSEUDO+1     !   to -99
          STASH_PSEUDO_LEVELS(JJ,II)=-99
        ENDDO
      ENDDO

      DO   I = 1,NUM_PSEUDO_LISTS
        NUM_LEVELS = LENPLST(I)
        STASH_PSEUDO_LEVELS(1,I) = NUM_LEVELS
        DO J = 1,NUM_LEVELS
          STASH_PSEUDO_LEVELS(J+1,I) = PSLIST_D(J,I)
        END DO
      END DO

!Transfer time series data to STASH_SERIES array
      IF (NSERIES >  0) THEN
!There are timeseries domains
!  Loop over STASHC domain profiles
        DO I=1,NDPROF
          IF (NPOS_TS(I)  >   0) THEN
!  This domain profile has a block of time series domains
!    J=time series block identifier (pointer):
            J=NPOS_TS(I)
!    STASH_SERIES_INDEX(1,J)=sequence no. of first record for
!                            ts block J in STASH_SERIES
            IF (J == 1) THEN
              STASH_SERIES_INDEX(1,J)=1
            ELSE
              STASH_SERIES_INDEX(1,J)= STASH_SERIES_INDEX(1,J-1)        &
     &                                +STASH_SERIES_INDEX(2,J-1)
            END IF
!    STASH_SERIES_INDEX(2,J)=no. of records in ts block J
            STASH_SERIES_INDEX(2,J)=NRECS_TS(J)
            IR1  =STASH_SERIES_INDEX(1,J)
            DO IR=IR1,IR1+NRECS_TS(J)-1
              STASH_SERIES(1,IR)=IG_TS
              STASH_SERIES(2,IR)=I1_TS
              STASH_SERIES(3,IR)=I51_TS
              STASH_SERIES(4,IR)=NLIM_TS(IR)
              STASH_SERIES(5,IR)=SLIM_TS(IR)
              STASH_SERIES(6,IR)=WLIM_TS(IR)
              STASH_SERIES(7,IR)=ELIM_TS(IR)
              STASH_SERIES(8,IR)=BLIM_TS(IR)
              STASH_SERIES(9,IR)=TLIM_TS(IR)
            END DO
          END IF
        END DO
      END IF

!     Initialise STINDEX and SI

      DO IE=1,NITEMS
      DO IS=0,NSECTS
      DO IM=1,N_INTERNAL_MODEL
        DO II=1,2
          STINDEX(II,IE,IS,IM)=0
        END DO
        SI    (IE,IS,IM)=1
        SI_LEN(IE,IS,IM)=0
        IF (IS == 0) THEN
          PPINDEX(IE,IM)=0
        END IF
      END DO
      END DO
      END DO

! 2. Read STASHindex and compute STASHWORK array lengths.
!    The Lth. row in STINDEX, SI, SI_LEN, PPINDEX, corresponds to the
!        Lth. internal model in INTERNAL_MODEL_LIST.
!    Output a formatted description of the selected diagnostics.

      ii  =0       ! Counter for checking no. of diags. printed
      L   =0       ! Counter for rows in STINDEX, SI, etc.
      DO K=1,N_INTERNAL_MODEL_MAX
      INT_MOD_INCLUDED=.FALSE.
! Find out whether int. model  'K' is included. If it is:
!  Set logical flag, increment row number
      DO KK=1,N_INTERNAL_MODEL_MAX
        IF (INTERNAL_MODEL_LIST(KK) == K) THEN
            INT_MOD_INCLUDED=.TRUE.
            L = L + 1
        END IF
      END DO
      IF (INT_MOD_INCLUDED) THEN
       DO J=0,NSECTS                        ! NSECTS=NSECTP (WSTLST)
       DO I=1,NITEMS                        ! NITEMS=NITEMP (WSTLST)
        IF (IN_S(1,K,J,I) >= 1) THEN        ! Entry in STASH list
          ii = ii + 1
          STINDEX(1,I,J,L) = INDX_S (1,K,J,I)  ! STASH index
          STINDEX(2,I,J,L) = INDX_S (2,K,J,I)
          SI     (  I,J,L) = IN_S   (1,K,J,I)  ! STASH lengths and
          SI_LEN (  I,J,L) = IN_S   (2,K,J,I)  !   addresses in D1
          IF (J == 0 .AND. PPIND_S(K,I) /= 0) THEN
            PPINDEX(I,  L) = PPIND_S(  K,  I)  ! Index for pp header
          END IF                               !   array
! Extract ppxref information to be passed into
!   diagnostic description routine
! DEPENDS ON: exppxc
          NAME = EXPPXC(K,J,I,                                          &
#include "argppx.h"
     &                         ICODE,CMESSAGE)
          ppxref_dat(ppx_model_number) = K
! DEPENDS ON: exppxi
          ppxref_dat(ppx_field_code  ) = EXPPXI(K,J,I,ppx_field_code  , &
#include "argppx.h"
     &                          ICODE,CMESSAGE)
! DEPENDS ON: exppxi
          ppxref_dat(ppx_data_type   ) = EXPPXI(K,J,I,ppx_data_type   , &
#include "argppx.h"
     &                          ICODE,CMESSAGE)
! DEPENDS ON: exppxi
          ppxref_dat(ppx_grid_type   ) = EXPPXI(K,J,I,ppx_grid_type   , &
#include "argppx.h"
     &                          ICODE,CMESSAGE)
! DEPENDS ON: exppxi
          ppxref_dat(ppx_lv_code     ) = EXPPXI(K,J,I,ppx_lv_code     , &
#include "argppx.h"
     &                          ICODE,CMESSAGE)
! DEPENDS ON: exppxi
          ppxref_dat(ppx_cf_levelcode) = EXPPXI(K,J,I,ppx_cf_levelcode, &
#include "argppx.h"
     &                          ICODE,CMESSAGE)
! DEPENDS ON: exppxi
          ppxref_dat(ppx_cf_fieldcode) = EXPPXI(K,J,I,ppx_cf_fieldcode, &
#include "argppx.h"
     &                          ICODE,CMESSAGE)
          DO IPK = 0,9
! DEPENDS ON: exppxi
          ppxref_dat(ppx_pack_acc+IPK) = EXPPXI(K,J,I,ppx_pack_acc+IPK, &
#include "argppx.h"
     &                          ICODE,CMESSAGE)
          END DO

          IF(PrintStatus >= PrStatus_Normal) THEN
!  Write a formatted description of the diagnostic to output file.
          DO IX=STINDEX(1,I,J,L),                                       &
                                                           ! Loop over
     &          STINDEX(1,I,J,L)+STINDEX(2,I,J,L)-1        !   entries
! DEPENDS ON: timer
                                   IF (LTIMER) CALL TIMER('DIAGDESC',3)
! DEPENDS ON: diagdesc
          CALL DIAGDESC(IX,NAME,STLIST(1,IX),ppxref_dat(1),             &
     &    stash_levels,num_stash_levels,num_level_lists,                &
     &    stash_pseudo_levels,num_stash_pseudo,num_pseudo_lists,        &
     &    sttabl,nsttims,nsttabl,                                       &
     &    stash_series,time_series_rec_len,nstash_series_records,       &
     &    stash_series_index,nstash_series_block)
! DEPENDS ON: timer
                                   IF (LTIMER) CALL TIMER('DIAGDESC',4)
          ENDDO
          ENDIF  ! PrintStatus test

        ELSE IF (IN_S(1,K,J,I) == -1) THEN
          ii = ii + 1
        END IF                       !  INDX_S(1,K,J,I) >= 1
        IF (ICODE >  0) GO TO 9999
       END DO                       !  Items
       END DO                       !  Sections

      END IF                       !  INT_MOD_INCLUDED
      END DO                       !  Models

      IF (ii /= N_PPXRECS) THEN
        WRITE(6,*) ' Error in INITCTL: N_PPXRECS not correct  ',II
        CMESSAGE='INITCTL  : N_PPXRECS not correct               '
        ICODE=1
        GO TO 9999
      END IF

!Assign values to PP_LEN2_LOOK, FT_OUTPUT
      DO I = 20,NUNITS
        ! If a bigger value than the computed one has been supplied
        ! by the user then we will use this.
        IF (PP_LEN2_LOOK(I) < PPlen2LkUp(I)) THEN
          PP_LEN2_LOOK(I) = PPlen2LkUp(I)
        END IF

        FT_OUTPUT(I)=FTOutUnit (I)
      END DO


!L 2.1 Find the max length in STASH_WORK and store in STASH_MAXLEN

      DO IM=1,N_INTERNAL_MODEL
      DO IS=1,NSECTS  !  Note not section Zero as the data is in D1
        STASH_MAXLEN(IS,IM)=1
        DO IE=1,NITEMS  ! Again only data not in D1
          IF(STINDEX(1,IE,IS,IM) /= 0) THEN
!         Item is in STASHlist ...
            IF (STLIST(st_input_code,STINDEX(1,IE,IS,IM)) == 1) THEN
!           ...  and input from STASHwork
!             ... input length not from ST_LIST as this is post STOCGT
              INPUT_LENGTH=SI_LEN(IE,IS,IM)
              STASH_MAXLEN(IS,IM)=STASH_MAXLEN(IS,IM)+INPUT_LENGTH
            ENDIF
          ENDIF
        END DO
      END DO
      END DO

!L
!L
!L----------------------------------------------------------------------
!L  3.   Set derived control variables for use in STASH/STWORK
!L
!L       Set PP_LEN2_LOOKUP to maximum PP_LEN2_LOOK value for any PP
!L       unit referenced in the STASHlist (minimum value possible is 8).
!L       Set MAX_STASH_LEVS to the maximum possible no of output levels
!L       for any diagnostic, allowing for possible pseudo-levels.
!L
      PP_LEN2_LOOKUP=8
      MAX_STASH_LEVS=1
      DO II=1,TOTITEMS
        IF (STLIST(st_output_code,II) <  0) THEN
! output is to PP file
          IF (PP_LEN2_LOOK(-STLIST(st_output_code,II))                  &
     &         >  PP_LEN2_LOOKUP)                                       &
     &        PP_LEN2_LOOKUP=PP_LEN2_LOOK(-STLIST(st_output_code,II))
        ENDIF
!
!       Input levels list/range is always longer than output
        IF (STLIST(st_input_bottom,II) == st_special_code) THEN
!          On special level
           NUM_LEVELS=1
        ELSE IF (STLIST(st_input_bottom,II) <  0) THEN
!          Using levels list, element 1 holds length.
           NUM_LEVELS=STASH_LEVELS(1,-STLIST(st_input_bottom,II))
        ELSE
!          Range
           NUM_LEVELS=                                                  &
     &        STLIST(st_input_top,II)-STLIST(st_input_bottom,II)+1
        END IF
        IF (STLIST(st_pseudo_in,II) /= 0) THEN
!          On pseudo levels
           NUM_PSEUDO_LEVELS=STASH_PSEUDO_LEVELS(1,                     &
     &                       STLIST(st_pseudo_in,II))
        ELSE
!          Not on pseudo levels
           NUM_PSEUDO_LEVELS=1
        END IF
        MAX_STASH_LEVS=MAX(MAX_STASH_LEVS,NUM_LEVELS*NUM_PSEUDO_LEVELS)
      END DO
! Round PP_LEN2_LOOKUP up to a multiple of 8
      PP_LEN2_LOOKUP=((PP_LEN2_LOOKUP+7)/8)*8


! DEPENDS ON: fill_d1_array
      CALL FILL_D1_ARRAY(                                               &
#include "argsts.h"
#include "argppx.h"
#include "argd1.h"
     &                  ICODE,CMESSAGE)
      IF (ICODE /= 0) GO TO 9999

! Initialise PEER files if required
! DEPENDS ON: peer_initialise
      CALL PEER_INITIALISE(ICODE,CMESSAGE)
      IF (ICODE /= 0) GO TO 9999

!----------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE INITCTL
#endif
