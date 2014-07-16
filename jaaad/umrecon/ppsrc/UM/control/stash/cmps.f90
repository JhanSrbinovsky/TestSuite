

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



      SUBROUTINE CMPS(IX,FIELD,ICOMP,NOP,APREC,IBIT,BASE,RMDI,          &
     &                ICODE,CMESSAGE)
!FPP$ NOCONCUR R
!FPP$ EXPAND (INSTIN)

      Implicit None

!     Subroutine arguments
      INTEGER IX,NOP,ICOMP(IX-1),IBIT
      INTEGER ICODE
      REAL FIELD(IX),APREC,BASE,RMDI
      CHARACTER CMESSAGE*80

!     Local variables
      LOGICAL OBTMIS,OBTZER
      INTEGER IMAP(IX),IMIS(IX),IZERO(IX),ITEMP(IX)
      INTEGER IGATH1(IX),IGATH2(IX)
      INTEGER I,IBASE,IMAX,I2,JBIT,JJ,JWORD,NMIS,NPOINT,NZERO
      INTEGER IDIV, jword_end, jbit_end, II, iword

      Real :: Atemp(Ix)
      Real :: Bprec

!     FUNCTIONS
      INTEGER ishift

!
!     CRAY VERSION C.1.1  16/11/90  P J SMITH
!
!         IX   =  LENGTH OF FIELD
!      FIELD   =  FIELD FOR COMPRESSING
!      ICOMP   =  RETURNED COMPRESSED DATA
!        NOP   =  NUMBER OF WORDS OF COMP FILLED BY THIS CALL
!      APREC   =  PRECISION
!       IBIT   =  NUMBER OF BITS INTO WHICH DATA PACKED
!       BASE   =  REFERENCE (MIN) VALUE FOR PACKED DATA
!       RMDI   =  MISSING DATA INDICATOR
!
!     INITIALISE VARIABLES/TEMP ARRAYS
!
      OBTMIS = .FALSE.
      OBTZER = .FALSE.
      BPREC  = 1./APREC
      DO I = 1,IX
          IMAP(I)  = 1
      END DO
      DO I = 1,IX-1
          ICOMP(I) = 0
      END DO
!
!     SCAN FIELD FOR MISSING DATA AND STORE RESULT IN IMIS,
!     SCAN FIELD FOR ZERO VALUES AND STORE RESULT IN IZERO,
!     SCAN FIELD FOR MINIMUM VALUE (IGNORE MISSING DATA INDICATOR)
!
      BASE  = 999999.
      NMIS  = 0
      NZERO = 0
      JJ    = 0

      DO I = 1,IX
       IF (FIELD(I) /= RMDI) THEN
         IMIS(I)=0
       ELSE
         IMIS(I)=1
       ENDIF
      END DO

      DO I = 1,IX
       IF (FIELD(I) /= 0.0) THEN
         IZERO(I)=0
       ELSE
         IZERO(I)=1
       ENDIF
      END DO

!
! GET NO. OF NON-RMDI POINTS + COMPRESS INDEX TO REMOVE THEM
!
       JJ  =0
       DO I=1,IX
         IF (FIELD(I) /= RMDI) THEN
           JJ        =JJ+1
           IGATH1(JJ)=I
         END IF
       END DO
!
      NMIS=IX-JJ
!
      IF(JJ /= 0)THEN
! REMOVE MISSING DATA
!DIR$ IVDEP
       DO I =1,JJ
       ATEMP(I)=FIELD(IGATH1(I))
       END DO
!
! GET NO. OF NON-ZERO (NON-RMDI) POINTS + COMPRESS INDEX TO REMOVE THEM
!
       NPOINT=0
       DO I  =1,JJ
         IF (ATEMP(I) /= 0.0) THEN
           NPOINT        =NPOINT+1
           IGATH2(NPOINT)=I
         END IF
       END DO
!
       NZERO=JJ-NPOINT
!
! SET BASE VALUE
       DO I =1,JJ
       IF(ATEMP(I) <  BASE) BASE=ATEMP(I)
       END DO
      ENDIF
!
!     CHECK IF A BITMAP FOR MISSING DATA IS NEEDED,
!
      IF (NMIS >  0) THEN
          OBTMIS = .TRUE.
      ELSE
          OBTMIS = .FALSE.
      END IF
!
!     ROUND BASE TO PRECISION REQUIRED
!
      IF (BASE /= 0.0) THEN
          BASE  = BASE*BPREC
          IBASE = NINT(BASE)
          BASE  = IBASE*APREC
      ELSE
          IBASE=0
      END IF
! Fujitsu vectorization directive
!OCL NOVREC
!
!     FIND DIFFERENCE FROM BASE AND SCALE
!     FIND MAXIMUM DIFFERENCE
!
      IMAX = 0
      DO 20 I = 1,JJ
          ATEMP(I) = ATEMP(I)*BPREC
          ITEMP(I) = NINT(ATEMP(I)) - IBASE
          IF(ITEMP(I) <  0) ITEMP(I) = 0
          IF (IMAX <  ITEMP(I)) IMAX = ITEMP(I)
   20 CONTINUE
!
!     FIND NUMBER OF BITS REQUIRED TO STORE MAX DIFFERENCE
!
      IBIT  = 0
      ! Enable a maximum range of 2^5 bits
      IF (IMAX  >   2.147483E+09) THEN
         ICODE = 2
         CMESSAGE='COEX: Unable to WGDOS pack to this accuracy'
         GOTO 999
      ELSE
      IF (IMAX >  0) THEN
          I2    = 1
          DO WHILE(IMAX >= I2)
              IBIT  = IBIT + 1
              I2    = I2*2
          ENDDO
      ENDIF
      ENDIF
!
!     SET START POSITION IN OUTPUT ARRAY
!
      JWORD = 1
      JBIT  = 63
!
!     IF BIT-MAPPING FOR MISSING DATA THEN INSERT BITMAP
!
      IF (OBTMIS) THEN
          DO 30 I = 1,IX
              IF (IMIS(I) == 1) THEN
                  ICOMP(JWORD) = IBSET(ICOMP(JWORD),JBIT)
                  IMAP(I)      = 0
              END IF
              IF (JBIT == 0) THEN
                  JBIT  = 63
                  JWORD = JWORD + 1
              ELSE
                  JBIT  = JBIT - 1
              END IF
   30     CONTINUE
      END IF
!
!     IF WORTHWHILE USE BIT MAP AND COMPRESS OUT ZEROS.
!
      IF (IBIT >  0) THEN
          IF (NZERO >  IX/IBIT) THEN
              OBTZER = .TRUE.
              DO 40 I = 1,IX
                  IF (IZERO(I) == 1) THEN
                      ICOMP(JWORD) = IBCLR(ICOMP(JWORD),JBIT)
                      IMAP(I)      = 0
                  ELSE
                      ICOMP(JWORD) = IBSET(ICOMP(JWORD),JBIT)
                  END IF
                  IF (JBIT == 0) THEN
                      JBIT  = 63
                      JWORD = JWORD + 1
                  ELSE
                      JBIT  = JBIT - 1
                  END IF
   40         CONTINUE
          ELSE
              OBTZER = .FALSE.
          END IF
!
!         IF BIT MAP INCLUDED FILL TO END OF CURRENT 32 BIT WORD
!         AND SET POINTERS TO START OF NEXT 32 BIT BOUNDARY.
!
          IF (OBTZER .OR. OBTMIS) THEN
              DO 50 I = 0,JBIT
                  ICOMP(JWORD) = IBSET(ICOMP(JWORD),I)
   50         CONTINUE
              IF (JBIT /= 63) THEN
                  IF (JBIT >= 31) THEN
                      JBIT = 31
!
! We have set bits in the bottom half of the 64-bit word - clear
! them again
!
! DEPENDS ON: ishift
                      icomp(jword)=ishift(ishift(icomp(jword), -32), 32)
                  ELSE
                      JWORD = JWORD + 1
                      JBIT = 63
                  ENDIF
              ELSE
!
! We have set all the bits in the 64-bit word - clear them again
!
                icomp(jword)=0
                JBIT = 63
              ENDIF
          ELSE
!
! We have set all the bits in the 64-bit word - clear them again
!
              icomp(jword)=0
              JBIT = 63
          END IF
!
!         IF BIT MAPPING ZEROS - COMPRESS OUT UNWANTED ZERO VALUES
!        (OTHERWISE PACK ALL NON RMDI DATA (JJ) )
!
          IF (OBTZER) THEN
!DIR$ IVDEP
           DO I= 1,NPOINT
           ITEMP(I)=ITEMP(IGATH2(I))
           END DO
          ELSE
           NPOINT=JJ
          END IF
!
!         MOVE INTO OUTPUT ARRAY USING MINIMUM NUMBER OF BITS REQUIRED
!
          jword_end = jword + NPOINT*IBIT/64
          jbit_end  = jbit  - MOD(NPOINT*IBIT,64)
          IF (MOD(IBIT,32) == 0) THEN
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
              iword = jword
!CDIR NODEP
              DO I = II+1,NPOINT,64/IDIV
!
! Adjust the position of the bits in the first word, and then
! insert the bits into the word
!
! DEPENDS ON: ishift
                  icomp(iword)=ior(ishift(itemp(i), (jbit+1)-ibit),     &
     &                             icomp(iword))
!
! Position the bits correctly for the second word, and 'or' them in
!
                  icomp(iword+1)=                                       &
! DEPENDS ON: ishift
     &             ior(ishift(itemp(i), min(64, 64-(ibit-(jbit+1)))),   &
     &                 icomp(iword+1))
                  iword = iword + IBIT/IDIV
              ENDDO
!
              JBIT = JBIT - IBIT
              IF (JBIT <  0) THEN
                  JWORD = JWORD + 1
                  JBIT  = JBIT + 64
              END IF
          ENDDO
          jword = jword_end
          jbit  = jbit_end
          IF (JBIT <  0) THEN
              JWORD = JWORD + 1
              JBIT  = JBIT + 64
          END IF

      ELSEIF(IBIT == 0) THEN
!
!         IF BIT MAP INCLUDED FILL TO END OF CURRENT 32 BIT WORD
!
          IF (OBTMIS) THEN
              DO 80 I = 0,JBIT
                  ICOMP(JWORD) = IBSET(ICOMP(JWORD),I)
   80         CONTINUE
          END IF
      END IF
!
!     CALCULATE LENGTH IN 32 BIT WORDS
!
      NOP = JWORD*64 - JBIT - 1
      NOP = (NOP+31)/32
!
!     SET FLAGS TO INDICATE BIT MAPS
!
      IF (OBTZER) IBIT = IBIT + 128
      IF (OBTMIS) IBIT = IBIT + 32


 999  CONTINUE
      RETURN
! Fujitsu vectorization directive
!OCL NOVREC

      END SUBROUTINE CMPS










