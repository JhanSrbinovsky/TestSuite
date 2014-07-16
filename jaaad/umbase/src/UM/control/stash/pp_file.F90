#if defined(C84_1A) || defined(FLDOP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE PP_FILE -----------------------------------------
!LL
!LL  Purpose:- To output a field to a PP_FILE
!LL
!LL  Tested under compiler CFT77
!LL  Tested under OS version 5.1
!LL
!LL TJ, RR      <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL   3.2  19/04/93  Code for new real missing data indicator (TCJ).
!LL  3.4   04/08/94  No packing indicator change from -26 to -99  PJS
!LL  4.1   22/11/96  Modify I/O calls for MPP use  P.Burton
!LL  4.3   30/04/97  Added code to use UM_SECTOR_SIZE to make transfers
!LL                  well-formed.
!LL                  B. Carruthers  Cray Research.
!LL  4.4   16/06/97  Add processing after the write, so
!LL                  that all the processors know the answer
!LL                    Author: Bob Carruthers, Cray Rsearch.
!LL  4.5   28/05/98  Code for parallel processing in COEX Packing
!LL                    Author: Paul Burton & Bob Carruthers
!LL  5.2   11/02/00 Code for Run Length Encoding  I. Edmond
!LL
!LL   5.1 May. 00  Coex now returns an error code set greater than 0
!LL                if invalid data is being stored. (i.e if the number
!LL                of bits per value used in the wgdos packed data is
!LL                >31). PPFILE checks the return code from coex and
!LL                throws an exception by returning the error code from
!LL                coex (summed over all pes). Ian Edmond
!    5.4    11/04/02 Set up extra data vectors for variable
!                    horizontal grids (only ocean catered for
!                    at this stage). Also prevent RLE when
!                    extra data vectors are present.  R. Hill
!    5.5    14/01/03 Ensure WFIO output buffer is extended
!                    to acommodate extra data. R. Hill
!    5.5    20/02/03 Allowed for inclusion of PARVARS comdeck
!                    on non-T3E MPP platforms to provide variables
!                    needed by COMVGRID comdeck.      P.Dando
!    6.0    15/01/03  Replace calls to SETPOS and BUFFOUT with
!                     buffered versions for PP files.
!                     JC Rioual, P. Selwood
!    6.1    13/07/04 Add packing in 2 stage gather.
!                    P.Selwood/B.Carruthers (Cray)
!    6.1    28/08/04 Buffered IO in C. JC Rioual P Selwood
!   6.2    19/01/05  Allow RLE to be active even if EXTRA data
!                    is present in a record (but not other
!                    forms of packing which could still result
!                    in corruption of the extra data) R. Hill

!LL  Programming standard: U M DOC  Paper NO. 4,
!LL
!LL  Logical components covered C4
!LL
!LL  Project TASK: C4
!LL
!LL  External documentation  C4
!LL
!LLEND-------------------------------------------------------------

!
!*L  ARGUMENTS:---------------------------------------------------
      SUBROUTINE PP_FILE(PPFIELD,LENBUF,NUM_WORDS,RMDI,COMP_ACCRCY,     &
     &PPHORIZ_OUT,UNITPP,IWA,N_COLS_OUT,N_ROWS_OUT,PACKING,             &
     & im_ident,LRLE,                                                   &
     &     PACKING_TYPE,current_io_pe,LEN_EXTRA                         &
     &    ,SROW_OUT,WCOL_OUT,ICODE,CMESSAGE)
      IMPLICIT NONE


      INTEGER LEN_EXTRA ! IN size of expected extra data
      INTEGER SROW_OUT  ! 1st southern row to output
      INTEGER WCOL_OUT  ! 1st western column to output
      CHARACTER*(80) CMESSAGE !OUT OUT MESSAGE FROM ROUTINE
!
      LOGICAL                                                           &
     &  PACKING                                                         &
                           !IN OVERALL Packing switch (T if pckng reqd)
     &,  LRLE

      INTEGER                                                           &
     &  ICODE                                                           &
                           !IN    RETURN CODE FROM ROUTINE
     &, LENBUF                                                          &
                           !IN     LENGTH OFF PP BUFFER
     &, UNITPP                                                          &
                           !IN     OUTPUT PP UNIT NUMBER
     &, LEN_IO                                                          &
                           !NOT USED, BUT NEEDED FOR BUFFOUT CALL
     &, im_ident                                                        &
                           !IN    Internal model identifier
     &, current_io_pe                                                   &
                           !IN     PE which will do the I/O
     &, ocode                                                           &
                        !Copy of the input value of ICODE
     &, ncode                                                           &
                        ! seperate item for each processor
     &, istat           ! Isum return status, GC_OK on success.
!
      INTEGER                                                           &
     &  N_ROWS_OUT                                                      &
                      !IN   PPHORIZ_OUT=N_ROWS_OUT*N_COLS_OUT
     &, N_COLS_OUT                                                      &
                      !IN    PPHORIZ_OUT=N_COLS_OUT*N_ROWS_OUT
     &, NUM_OUT                                                         &
                      !IN    NUMBER OF COMPRESSED (32 BIT) WORDS
     &, COMP_ACCRCY                                                     &
                      !IN    PACKING ACCURACY IN POWER OF 2
     &, U_ROWS                                                          &
                      !IN    NO OF U,V, ROWS
     &, P_ROWS                                                          &
                      !IN    PRESS/TEMP ROWS
     &, PPHORIZ_OUT                                                     &
                      !IN    SIZE OF OUTPUT FIELD
     &, NUM_WORDS                                                       &
                      !IN    NUMBER OF 64 BIT WORDS WORTH OF DATA
     &, PACKING_TYPE  ! OUT set to 1 if WGDOS packing else set to zero.
!
      REAL                                                              &
     &  PPFIELD(PPHORIZ_OUT)                                            &
                               !INOUT ARRAY TO STORE PPDATA
     &, BUFOUT(LENBUF+LEN_EXTRA)                                        &
                                  !OUTPUT PP BUFFER (ROUNDED UP)
     &, RMDI                   !IN     Missing data indicator
!
!
#include "parvars.h"

#include "comvgrid.h"

!
!dir$ cache_align bufout
#include "cntl_io.h"
#include "chsunits.h"
#include "csubmodl.h"
!*---------------------------------------------------------------------

!*L  WORKSPACE USAGE:-------------------------------------------------
!   DEFINE LOCAL WORKSPACE ARRAYS: 1 REAL ARRAY
!   AT FULL FIELD LENGTH
!
!*---------------------------------------------------------------------
!
!*L EXTERNAL SUBROUTINES CALLED---------------------------------------
      EXTERNAL SETPOS,COEX,BUFFOUT
!*------------------------------------------------------------------
!L  MAXIMUM VECTOR LENGTH ASSUMED IS (ROWS-1) * ROWLENGTH
!L---------------------------------------------------------------------
!----------------------------------------------------------------------
!    DEFINE LOCAL VARIABLES
      INTEGER                                                           &
     &  ML                                                              &
                      !     LONGITUDE COUNTER
     &, JL                                                              &
                      !     LATITUDE COUNTER
     &, IWA                                                             &
                      !     RECORD NUMBER
     &, II                                                              &
                      !     COUNTER
     &, LENGTH_FULLWRD                                                  &
                      !     LENGTH IN BITS OF FULLWORD VAR
     &, LEN_BUF_WORDS !     NUM_WORDS ROUNDED BY 512 AND ACTUALLY

      INTEGER                                                           &
     &  JJ            !     Local counter

      REAL                                                              &
     &  IX            !     RETURN VALUE FROM UNIT COMMAND

      LOGICAL                                                           &
     &  UV                 !

!    Check if an error has already been encountered, and get out
!    if it has.
      ocode = 0
      IF (icode  >   0) then
         goto 999
      ELSE IF (icode  <   0)then
         ocode = icode
         icode = 0
      END IF

!
!
!    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS IPPLOOK GOES
!    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
!    IWA=325119 then IWA=325120 to 388095 then 388096 etc
!
!======================================================================
      LENGTH_FULLWRD=64   !   LENGTH IN BITS OF FULLWORD VAR
!L   At this point packing,if required,will be done using the WGDOS
!L   method of packing.
      PACKING_TYPE=0
! Note the value of -26 corresponds to -15 (F) in ppxref.
! The packing acuracy is scaled to allow greater accuracy.
! Packing will only be attempted if there are at least 2 points per row
! in the PPfield.
!
!      Run length encoding applies only to unpacked ocean fields, but
!      most ocean fields remain unpacked even when packing profiles
!      5 or 6 are set. Hence selecting, for example, both packing
!      profile 5 and run length encoding makes sense.
      IF(PACKING .AND. N_COLS_OUT  >=  2) THEN
        ! 1. Climate wgdos packing has been selected for current file
        !    stream via UMUI
        IF(COMP_ACCRCY  <=  -99 .AND. LRLE .AND.                        &
     &     im_ident  ==  ocean_im )THEN
           ! 2. STASH packing profile for the field is -99
           ! 3. Run Length Encoding has been selected
           ! 4. Submodel is Ocean
           PACKING_TYPE = 4
        ELSE IF(COMP_ACCRCY  >   -99) THEN
           ! 2. STASH packing profile for the field is set.
           PACKING_TYPE = 1
        ENDIF
      ELSE
         ! 1. Packing may or may not have been selected for current
         !    file stream. This section of code is not executed when
         !    GRIB packing selected.
         IF (LRLE .AND.                                                 &
     &       im_ident  ==  ocean_im) THEN
           ! 2. Run Length Encoding has been selected
           ! 3. Submodel is Ocean
           ! This should be safe at least as far as output
           ! generation goes - the issue will be post
           ! processing functionality - does that have
           ! all the necessary intelligence to unpack
           ! RLE fields AND process extra data
           PACKING_TYPE = 4
         ENDIF
      END IF
!
      IF(PACKING_TYPE == 1)THEN
#if defined(PRINT84)
        WRITE(6,*)'*********  PPOUT PACKING REQD***********  '
#endif
#if defined(T3E) && defined(MPP)
! DEPENDS ON: mpp_coex
      CALL MPP_COEX(PPFIELD,PPHORIZ_OUT,BUFOUT,LENBUF,N_COLS_OUT,       &
     &              N_ROWS_OUT,NUM_OUT,COMP_ACCRCY,.TRUE.,RMDI,         &
     &              1,1,current_io_pe,1,icode,cmessage)
       call gc_isum(1,nproc,istat,icode)

       IF(icode >  0) GOTO 999
#else
! DEPENDS ON: coex
        CALL COEX(PPFIELD,PPHORIZ_OUT,BUFOUT,LENBUF,N_COLS_OUT,         &
     &  N_ROWS_OUT,NUM_OUT,COMP_ACCRCY,.TRUE.,RMDI,LENGTH_FULLWRD,      &
     &  icode,cmessage)
#endif

        NUM_WORDS=(NUM_OUT+1)/2 ! Round up to the nearest 64 Bit CRAY Wd
!  COEX returns the number of IBM words needed to hold the packed data
!                             ~~~
        LEN_BUF_WORDS=((NUM_WORDS+um_sector_size-1)/um_sector_size)*    &
     &   um_sector_size
#if defined(PRINT84)
        WRITE(6,*)'NUM_OUT',NUM_OUT
#endif
      ELSE IF(PACKING_TYPE == 4)THEN
        if (mype  ==  current_io_pe) then
! DEPENDS ON: runlen_encode
          CALL RUNLEN_ENCODE(PPFIELD,PPHORIZ_OUT,BUFOUT,PPHORIZ_OUT,    &
     &                     NUM_OUT,RMDI,ICODE,CMESSAGE)
        ! Size of run length encoded data is greater than unpacked
        ! field therefore leave field unpacked.
          if (NUM_OUT  >=  PPHORIZ_OUT) then
            PACKING_TYPE = 0
            DO JJ=1,PPHORIZ_OUT
              BUFOUT(JJ) = PPFIELD(JJ)
            END DO
            NUM_WORDS=PPHORIZ_OUT
            LEN_BUF_WORDS=LENBUF
          else
            NUM_WORDS=NUM_OUT
            LEN_BUF_WORDS=((NUM_WORDS+um_sector_size-1)/um_sector_size)*&
     &      um_sector_size
          endif

        end if
      ELSE  ! No packing required.
#if defined(PRINT84)
        WRITE(6,*)'FROM PPOUT  N_ROWS_OUT  N_COLS_OUT'
        WRITE(6,*)N_ROWS_OUT,N_COLS_OUT
#endif
        DO 1 JJ=1,PPHORIZ_OUT
        BUFOUT(JJ) = PPFIELD(JJ)
    1   CONTINUE
        NUM_WORDS=PPHORIZ_OUT
        LEN_BUF_WORDS=LENBUF
      ENDIF

      IF (VAR_GRID_TYPE >  0) THEN

         ! If we have a variable horizontal grid then we
         ! must add the grid data to the end of the PP record
         ! as "EXTRA DATA".
! DEPENDS ON: extra_variable_grid
         CALL EXTRA_VARIABLE_GRID(BUFOUT(NUM_WORDS+1),LEN_EXTRA         &
     &            ,SROW_OUT,N_ROWS_OUT                                  &
     &            ,WCOL_OUT,N_COLS_OUT)

         ! Adjust output buffer size for WFIO.
         NUM_WORDS = NUM_WORDS + LEN_EXTRA

         IF ((PACKING_TYPE == 1).OR.(PACKING_TYPE == 4)) THEN
           LEN_BUF_WORDS=((NUM_WORDS+um_sector_size-1)/um_sector_size)* &
     &                     um_sector_size ! No of words output
          ELSE
             LEN_BUF_WORDS=LENBUF+LEN_EXTRA
         ENDIF
      ENDIF ! If variable grid data needed adding

      if (mype  ==  current_io_pe) then
      DO JJ=NUM_WORDS+1,LEN_BUF_WORDS
        BUFOUT(JJ)= 0.0
      ENDDO
#if !defined(MPP)
! DEPENDS ON: setpos
      CALL SETPOS(UNITPP,IWA,ICODE)
! DEPENDS ON: buffout
      CALL BUFFOUT(UNITPP,BUFOUT(1),LEN_BUF_WORDS,LEN_IO,IX)
#else
      CALL SETPOS_single(UNITPP,IWA,ICODE)
      CALL BUFFOUT_single(UNITPP,BUFOUT(1),LEN_BUF_WORDS,LEN_IO,IX)

#endif
!     WRITE(6,102) IWA,LEN_BUF_WORDS
  100 FORMAT(//,32X,'   ARRAY BUFOUT AT END OF PPOUT ',//,32(10F8.0/))
  102 FORMAT(' FROM PP_FILE    IWA  LEN_BUF_WORDS ',2I12)
!
      IF (IX /= -1.0.OR.LEN_IO /= LEN_BUF_WORDS) THEN
! DEPENDS ON: ioerror
        CALL IOERROR('Buffer out Data Field',IX,LEN_IO,                 &
     &                LEN_BUF_WORDS)
        CMESSAGE='PPFILE  : I/O error - PP Data Field Output'
        ICODE=7
        RETURN
      ENDIF
      endif ! (mype  ==  current_io_pe)
!     If we have found an error, leave it in icode.  If no error
!     occurred then check if the original input value of icode was
!     non-zero (a previous untrapped error/warning), and copy this
!     back into ICODE before eaving the routine.
      IF (icode  ==  0 .and. ocode  /=  0) then
         icode = ocode
      END IF
  999 CONTINUE
      RETURN
      END SUBROUTINE PP_FILE
#endif
