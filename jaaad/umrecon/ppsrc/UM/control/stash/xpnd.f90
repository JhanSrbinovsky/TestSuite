

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL   SUBROUTINE COEX,COEX2,CMPS,XPND,INSTIN,EXTRIN -----------------
!LL
!LL   PURPOSE:   PACK TO AND UNPACK FROM WGDOS FORMAT
!LL
!LL   (Note that for optimal performance the routines INSTIN and
!LL    EXTRIN are inline expanded (enable option 8 on fpp)
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.2  19/04/93  Code for ne real missing data indicator (TCJ).
!LL   4.2  Nov. 96   T3E migration: WHENNE,WHENEQ replaced by
!LL                   portable fortran code.   S.J.Swarbrick
!LL   4.3  Apr. 97 T3E migration: Calling the functions CRI2IBM/IBM2CRI
!LL                directly replaces the calls to the (now unsupported)
!LL                routines USICTI, USICTC, USSCTC, USSCTI. Ian Edmond
!LL   5.0 05/03/99 Remove DEF,C98_1A. D. Robinson.
!LL   5.1  23/03/00 Add declaration for CRI2IBM. S.D.Mullerworth
!LL   5.2  15/11/00 Allow use in reconfiguration. P.Selwood.
!LL
!LL   5.1 May. 00  The max number of bits per value used in the wgdos
!LL                packed data is 31. This implies that the maximum
!LL                difference from the base value which can be stored
!LL                is 2.147483E+09. A check is now added to prevent
!LL                values greater than this being stored.  Ian Edmond
!LL   5.3  07/06/01 Declare cri2ibm.  A van der Wal
!     5.5  14/02/03 Add Implicit None + appropriate declarations where
!                   necessary. Tidy up some declarations. T.White
!LL   5.5  28/02/03 Insert code for portable data conversion routines
!LL                 to replace Cray-specific CRI2IBM etc.     P.Dando
!     6.0  10/09/03 Conversion of portable data conversion routines
!                   (IEEE2IBM etc) into functions with error return
!                   codes matching those of CRAY routines.   P.Dando
!     6.0  25/09/03 Optimisation for NEC.
!                   Bob Carruthers, Gerrit v.d. Velde, Paul Selwood.
!     6.1  23/08/04 Add FLDCALC def. D Robinson
!     6.1  13/09/04 Fix bug with insertion of zero and missing data
!                   bit maps (filling to end of next 32-bit word)
!                   in CMPS routine when row_length is a multiple of
!                   64.                                    P.Dando
!     6.2  26/10/05 Added function ishift which behaves in a standard
!                   way (as expected by coex) to fix WGDOS errors in
!                   ifc compiled runs.                     T. Edwards
!     6.2  24/11/05 Optimisation for NEC. Two new subroutines added -
!                   CMPS_ALL and XPND_ALL. J-C.Rioual(NEC)/D.Robinson.
!LL
!LL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
!LL   STANDARD B, VERSION 2, DATED 18/01/90
!LL
!LL  Logical component number: S72
!LL
!LL   SYSTEM TASK: P7
!LL   OWNED BY P J SMITH
!LL
!LL   DOCUMENTATION:  ??????
!LLEND-------------------------------------------------------------





!     This is a portable replacement for ishft, which behaves as the UM
!     expects it too when shifting more than 64 bits.






      SUBROUTINE XPND(IX,ICOMP,FIELD,APREC,IBIT,BASE,NOP,RMDI)
!FPP$ NOCONCUR R
!FPP$ EXPAND (EXTRIN)

      Implicit None

!     Subroutine arguments
      INTEGER IX,NOP,ICOMP(NOP+1),IBIT
      REAL FIELD(IX),APREC,BASE,RMDI
      REAL ATEMP(IX)

!     Local variables
      LOGICAL OBTMIN,OBTMIS,OBTZER,OBTMAP
      INTEGER :: I,JWORD,JBIT,JJ,NPOINT,II,IEND
      INTEGER :: IMAP(IX),IMIS(IX),IZERO(IX),IMIN(IX)
      INTEGER :: full_word, iword, idiv, jword_end, jbit_end
      INTEGER IGATH(IX)

!     FUNCTIONS
      INTEGER ishift

!
!     CRAY VERSION C.1.1  16/11/90  P J SMITH
!
!         IX   =  LENGTH OF FIELD
!      ICOMP   =  DATA FOR EXPANSION
!      FIELD   =  FIELD OF EXPANDED DATA
!      APREC   =  PRECISION
!       IBIT   =  NUMBER OF BITS INTO WHICH DATA PACKED
!       BASE   =  REFERENCE (MIN) VALUE FOR PACKED DATA
!        NOP   =  SIZE OF COMP
!       RMDI   =  MISSING DATA INDICATOR
!
!
!     INITIALISE VARIABLES/TEMP ARRAYS
!
      OBTMAP   = .FALSE.
      OBTMIS   = .FALSE.
      OBTMIN   = .FALSE.
      OBTZER   = .FALSE.
      DO I = 1,IX
          IMAP(I)  = 1
          IMIS(I)  = 0
          IMIN(I)  = 0
          IZERO(I) = 0
      END DO
!
!     CHECK IF BITMAP USED FOR ZERO VALUES
!
      IF (IBIT >= 128) THEN
          OBTZER = .TRUE.
          OBTMAP = .TRUE.
          IBIT   = IBIT-128
      ELSE
          OBTMIN = .FALSE.
      ENDIF
!
!     CHECK IF BITMAP USED FOR MINIMUM VALUES (NOT CURRENTLY USED)
!
      IF (IBIT >= 64) THEN
          OBTMIN = .TRUE.
          OBTMAP = .TRUE.
          IBIT   = IBIT - 64
      ELSE
          OBTMIN = .FALSE.
      ENDIF
!
!     CHECK IF BITMAP USED FOR MISSING DATA
!
      IF (IBIT >= 32) THEN
          OBTMIS = .TRUE.
          OBTMAP = .TRUE.
          IBIT   = IBIT - 32
      ELSE
          OBTMIS = .FALSE.
      ENDIF
!
!     SET START POSITION IN ICOMP
!
      JWORD = 1
      JBIT  = 63
!
!     EXTRACT BITMAPS
!
      IF (OBTMIS) THEN
!
!         EXTRACT MISSING DATA BITMAP
!
          DO 10 I=1,IX
              IF (BTEST(ICOMP(JWORD),JBIT)) THEN
                  IMIS(I) = 1
                  IMAP(I) = 0
              ELSE
                  IMIS(I) = 0
              ENDIF
              IF(JBIT >  0) THEN
                  JBIT    = JBIT - 1
              ELSE
                  JBIT    = 63
                  JWORD   = JWORD + 1
              ENDIF
   10     CONTINUE
      ENDIF
      IF (OBTMIN) THEN
!
!         EXTRACT MINIMUM VALUE BITMAP (NOT USED AT PRESENT)
!
          DO 20 I=1,IX
              IF(BTEST(ICOMP(JWORD),JBIT)) THEN
                  IMIN(I) = 1
                  IMAP(I) = 0
              ELSE
                  IMIN(I) = 0
              ENDIF
              IF(JBIT >  0) THEN
                  JBIT    = JBIT - 1
              ELSE
                  JBIT    = 63
                  JWORD   = JWORD + 1
              ENDIF
   20     CONTINUE
      ENDIF
      IF (OBTZER) THEN
!
!         EXTRACT ZERO VALUE BITMAP
!
!         For faster execution we first check if all 64 bits in a word
!         are zero before we test the bits separately.
!
          IEND = MIN(IX,JBIT+1)
          DO 30 I=1,IEND
              IF(BTEST(ICOMP(JWORD),JBIT-I+1)) THEN
                  IZERO(I)= 0
              ELSE
                  IZERO(I)= 1
                  IMAP(I) = 0
              ENDIF
   30     CONTINUE
          JBIT = JBIT - IEND
          IF(JBIT <  0) THEN
              JBIT    = 63
              JWORD   = JWORD + 1
          ENDIF
          DO 35 II=IEND+1,IX,64
              IF (ICOMP(JWORD) == 0) THEN
                  IZERO(II:MIN(IX,II+64-1)) = 1
                  IMAP (II:MIN(IX,II+64-1)) = 0
              ELSE
                  DO 32 I=II,MIN(IX,II+64-1)
                      IF(BTEST(ICOMP(JWORD),JBIT-I+II)) THEN
                          IZERO(I)= 0
                      ELSE
                          IZERO(I)= 1
                          IMAP(I) = 0
                      ENDIF
   32             CONTINUE
              ENDIF
              JBIT = JBIT - ( MIN(IX,II+64-1) - II + 1 )
              IF(JBIT <  0) THEN
                  JBIT    = 63
                  JWORD   = JWORD + 1
              ENDIF
   35     CONTINUE
      ENDIF
!
!     IF BIT MAP USED FIND NUMBER OF POINTS STORED
!     AND RESET POINTERS TO BEGINNING OF 32 BIT BOUNDARY
!
      IF (OBTMAP) THEN
          NPOINT = 0
          DO 40 I=1,IX
              IF (IMAP(I) == 1) NPOINT = NPOINT + 1
   40     CONTINUE
          IF (JBIT /= 63) THEN
              IF (JBIT >= 31) THEN
                  JBIT  = 31
              ELSE
                  JBIT  = 63
                  JWORD = JWORD + 1
              ENDIF
          ENDIF
      ELSE
          NPOINT = IX
      ENDIF
      IF (IBIT >  0) THEN
!
!         UNPACK SCALED VALUES TO TEMP ARRAY
!
          jword_end = jword + (NPOINT*IBIT)/64
          jbit_end  = jbit  - MOD(NPOINT*IBIT,64)
          IF (MOD(IBIT,64) == 0) THEN
              IDIV = 64
          ELSE IF (MOD(IBIT,32) == 0) THEN
              IDIV = 32
          ELSE IF (MOD(IBIT,16) == 0) THEN
              IDIV = 16
          ELSE IF (MOD(IBIT,8) == 0) THEN
              IDIV = 8
          ELSE IF (MOD(IBIT,4) == 0) THEN
              IDIV = 4
          ELSE IF (MOD(IBIT,2) == 0) THEN
              IDIV = 2
          ELSE
              IDIV = 1
          ENDIF
          DO II = 0, 63/IDIV
!CDIR NODEP
          iword = jword
          DO I = II+1,NPOINT,64/IDIV
!
! Get the significant bits from the first word into the
! top of 'full_word'
!
! DEPENDS ON: ishift
              full_word=ishift(icomp(iword), 63-jbit)
!
! Shift the second word to generate spaces for the part from the
! first word, and then combine the two parts in 'full_word'
!
              full_word=ior(full_word,                                  &
! DEPENDS ON: ishift
     &         ishift(icomp(iword+1), -(jbit+1)))
!
! Now shift 'full_word' down so that the required integer is in
! the bottom most bits
!
! DEPENDS ON: ishift
              full_word=ishift(full_word, -(64-ibit))
!
!         ADD DIFFERENCES TO MINIMUM VALUE AND UNSCALE
!
              ATEMP(I)=full_word*APREC+BASE
              iword = iword + IBIT/IDIV
          ENDDO
!
          JBIT = JBIT - IBIT
          IF (JBIT <  0) THEN
              JWORD = JWORD + 1
              JBIT  = JBIT + 64
          ENDIF

          ENDDO
!
! Calculate the position in the source array
!
          jword = jword_end
          jbit  = jbit_end
          IF (JBIT <  0) THEN
              JWORD = JWORD + 1
              JBIT  = JBIT + 64
          END IF
!
!
!         MOVE INTO UNPACKED ARRAY FIELD
!
! FIRST GET GATHER INDEX
!
          JJ  =0
          DO I=1,IX
            IF (IMAP(I) == 1) THEN
              JJ       =JJ+1
              IGATH(JJ)=I
            END IF
          END DO
!
          DO I=1,IX
          FIELD(I)=0.
          END DO

          DO I=1,JJ
          FIELD(IGATH(I)) = ATEMP(I)
          END DO
!
!         IF MINIMUMS BIT MAPPED FILL ZEROS IN FIELD WITH BASE
!
          IF (OBTMIN) THEN
              DO 80 I=1,IX
                  IF(IMIN(I) == 1) FIELD(I) = BASE
   80         CONTINUE
          ENDIF
!
!         IF MISSING DATA BIT MAPPED FILL ZEROS IN FIELD WITH MDI
!
          IF (OBTMIS) THEN
              DO 90 I=1,IX
                  IF(IMIS(I) == 1) FIELD(I) = RMDI
   90         CONTINUE
          ENDIF

      ELSEIF (IBIT == 0) THEN

!
!         ALL POINTS IN ROW HAVE SAME VALUE E.G. POLE ROW
!
          DO 100 I=1,IX
              FIELD(I)=BASE
  100     CONTINUE
!
!         IF MISSING DATA BIT MAPPED FILL ZEROS IN FIELD WITH MDI
!
          IF (OBTMIS) THEN
              DO 110 I=1,IX
                  IF(IMIS(I) == 1) FIELD(I) = RMDI
  110         CONTINUE
          ENDIF
      ENDIF
!
      RETURN

      END SUBROUTINE XPND







