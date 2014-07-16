

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




      SUBROUTINE COEX2(FIELD,M,ICOMP,N,IX,IY,NUM,ISC,OCO,RMDI,          &
     &                 ICODE,CMESSAGE)

      Implicit None

!     Subroutine arguments
      INTEGER N,ICOMP(N),M,IX,IY,NUM,ISC
      INTEGER ICODE
      CHARACTER CMESSAGE*80
      REAL FIELD(M),RMDI
      LOGICAL OCO

!     Local variables
      INTEGER JJ,IST,ICX,JCX,NOB,IER2
      INTEGER IC(IY),ICB(IY),NOP(IY),IBIT(IY),ISTART(IY)
      INTEGER JCOMP(IX,IY),IERR1(IY),IERR2(IY),IERR3(IY)
      Integer :: Ip
      Integer :: Ierr
      Real :: Acc
      Real :: Aprec
      Real :: Base(Iy)

!     External functions used




      Integer :: IEEE2IBM
      Integer :: IBM2IEEE


!
!     CRAY VERSION C.1.1  16/11/90  P J SMITH
!     CRAY ROUTINE     64 BIT WORDS - CRAY FULLWORD
!
!                     OCO=.TRUE.                 OCO=.FALSE.
!
!      FIELD   =  FIELD FOR COMPRESSING     RETURNED EXPANDED DATA
!          M   =  SIZE OF FIELD             SIZE OF FIELD
!      ICOMP   =  RETURNED COMPRESSED DATA  FIELD FOR EXPANSION
!          N   =  SIZE OF COMP                 -------
!         IX   =  X DIMENSION OF FIELD      X DIMENSION OF FIELD
!         IY   =  Y DIMENSION OF FIELD      Y DIMENSION OF FIELD
!        NUM   =  TOTAL NO. OF COMPRESSED      -------
!                 (32 BIT) WORDS RETURNED
!        ISC   =  ACCURACY IN POWER OF 2    ACCURACY IN POWER OF 2
!        OCO   =  .TRUE. FOR COMPRESSION    .FALSE. FOR EXPANSION
!        IRC   =  RETURN CODE FOR EACH ROW  RETURN CODE FOR EACH ROW
!
!     INITIALISE TEMP. VARIABLES/ARRAYS
!
      IF (.NOT. OCO) THEN
        DO JJ=1,IY
          NOP(JJ)=0
          IBIT(JJ)=0
        ENDDO
        DO JJ=1,IY
          BASE(JJ)=0.0
          IBIT(JJ)=0.0
        ENDDO
      ENDIF

      DO JJ=1,IY
          DO JCX=1,IX
             JCOMP(JCX,JJ) = 0
          END DO
      END DO
      IF (OCO) THEN
!
!         COMPRESSION OF DATA
!
          IC(1) = 1
          ACC   = ISC
          APREC = 2.**ACC
!
!         PUT PACKED DATA FOR EACH ROW INTO TEMP ARRAY JCOMP
!
!FPP$ CNCALL
          DO 10 JJ = 1,IY
!             IP      = POSITION IN INPUT ARRAY
              IP      = (JJ-1)*IX + 1

! DEPENDS ON: cmps
              CALL CMPS(IX,FIELD(IP),JCOMP(2,JJ),NOP(JJ),APREC,         &
     &         IBIT(JJ),BASE(JJ),RMDI,ICODE,CMESSAGE)

            IF (ICODE  /=  0) GOTO 9990
   10     CONTINUE
!
!         ADD BASE VALUE, NUMBER OF BITS, AND LENGTH
!
!FPP$ CNCALL
          DO 11 JJ = 1,IY





              IERR1(JJ)=IEEE2IBM(3,1, JCOMP(1,JJ),0, BASE(JJ),1,64,32)
              IERR2(JJ)=IEEE2IBM(2,1,JCOMP(1,JJ),32,IBIT(JJ),1,64,16)
              IERR3(JJ)=IEEE2IBM(2,1,JCOMP(1,JJ),48,NOP(JJ),1,64,16)

   11     CONTINUE
!
!         CHECK ROW HEADER AND SET RETURN CODES
!
!FPP$ NOCONCUR
!         DO 12 JJ = 1,IY
!             IF (IERR1(JJ) /= 0) IRC(JJ) = IRC(JJ) + 1
!             IF (IERR2(JJ) /= 0) IRC(JJ) = IRC(JJ) + 2
!             IF (IERR3(JJ) /= 0) IRC(JJ) = IRC(JJ) + 4
!             IF (JCOMP(1,JJ) == 0) THEN
!                 IF (IBIT(JJ) /= 0) IRC(JJ) = IRC(JJ) + 8
!                 IF (NOP(JJ)  /= 0) IRC(JJ) = IRC(JJ) + 16
!             ENDIF
!  12     CONTINUE
!
!         CALCULATE POSITIONS IN OUTPUT ARRAY FOR PACKED ROWS
!         (FIRST PACKED ROW STARTS AT WORD 1; BIT 31)
!
          IC(1)     = 2
          ICB(1)    = -1
          ISTART(1) = 5
!FPP$ NOCONCUR
          DO 20 JJ = 2,IY
              IF (MOD(NOP(JJ-1),2) == 1) THEN
                  IC(JJ ) = IC(JJ-1) + NOP(JJ-1)/2 + 1
                  ICB(JJ) = -ICB(JJ-1)
                  IF (ICB(JJ) >  0) IC(JJ) = IC(JJ) + 1
              ELSE
                  IC(JJ)  = IC(JJ-1) + (NOP(JJ-1)+1)/2 + 1
                  ICB(JJ) = ICB(JJ-1)
              ENDIF
              ISTART(JJ)  = 5
              IF(ICB(JJ) == 1) ISTART(JJ) = 1
   20     CONTINUE
!
!         MOVE TEMP. ARRAY INTO OUTPUT ARRAY
!
!FPP$ NOCONCUR
          DO 30 JJ=1,IY
              NOB  = NOP(JJ)*4 + 8
! CHECK IF PACKED FIELD GREATER THAN UN PACKED FIELD
!             IF(NOB >  IX*8)IRC(JJ)=IRC(JJ)+32
              IST  = ISTART(JJ)
              ICX  = IC(JJ)



              CALL MOVEBYTES(JCOMP(1,JJ),1,NOB,ICOMP(ICX),IST)

   30     CONTINUE
!
!         INSERT TOTAL LENGTH OF THIS FIELD
          NUM = IC(IY)*2 + NOP(IY)
          IF (ICB(IY) <  0) NUM = NUM + 1



          IER2=IEEE2IBM(2,1,ICOMP(1),0,NUM,1,64,32)

!
!         END OF COMPRESSION SECTION
!
      ELSE
!
!         EXPANSION SECTION
!
          ACC   = ISC
          APREC = 2.**ACC
          ICX   = 2
          JCX   = -1
!FPP$ NOCONCUR
          DO 40 JJ = 1,IY
!
!             MOVE PACKED ROWS INTO TEMP ARRAYS
!
              IF (JCX <  0) THEN
!
!                 EXTRACT BASE, NO. BITS, NO 32 BIT WORDS
!





                  IERR=IBM2IEEE(3,1,ICOMP(ICX),32,BASE(JJ),1,64,32)
                  IERR=IBM2IEEE(2,1,ICOMP(ICX+1),0,IBIT(JJ),1,64,16)
                  IERR=IBM2IEEE(2,1,ICOMP(ICX+1),16,NOP(JJ),1,64,16)

!                 SAVE START POSITION OF ROW
                  IC(JJ)     = ICX
                  ISTART(JJ) = 5
              ELSE





                  IERR=IBM2IEEE(3,1,ICOMP(ICX),0,BASE(JJ),1,64,32)
                  IERR=IBM2IEEE(2,1,ICOMP(ICX),32,IBIT(JJ),1,64,16)
                  IERR=IBM2IEEE(2,1,ICOMP(ICX),48,NOP(JJ),1,64,16)

!                 SAVE START POSITION OF ROW
                  IC(JJ)     = ICX
                  ISTART(JJ) = 1
              END IF
!
!             CALCULATE START OF NEXT ROW
!
              IF (MOD(NOP(JJ),2) == 1) THEN
                  ICX   = ICX + NOP(JJ)/2 + 1
                  JCX   = -JCX
                  IF (JCX >  0) ICX = ICX + 1
              ELSE
                  ICX   = ICX + (NOP(JJ)+1)/2 + 1
              END IF
   40     CONTINUE
!
!         MOVE EACH PACKED ROW INTO TEMP ARRAY JCOMP
!
!FPP$ NOCONCUR
          DO 50 JJ = 1,IY
              ICX  = IC(JJ)
              IST  = ISTART(JJ)
              NOB  = NOP(JJ)*4 + 8



              CALL MOVEBYTES(ICOMP(ICX),IST,NOB,JCOMP(1,JJ),1)

   50     CONTINUE
!
!         CALCULATE START OF EACH ROW IN FIELD
!
          ICX = 1
          DO 60 JJ = 1,IY
              IC(JJ) = ICX
              ICX    = ICX + IX
   60     CONTINUE
!
!         UNPACK DATA INTO FIELD
!
!FPP$ CNCALL
          DO 70 JJ=1,IY
              ICX    = IC(JJ)
! DEPENDS ON: xpnd
              CALL XPND(IX,JCOMP(2,JJ),FIELD(ICX),APREC,IBIT(JJ),       &
     &                                            BASE(JJ),NOP(JJ),RMDI)
   70     CONTINUE
      END IF
 9990  CONTINUE
      RETURN

      END SUBROUTINE COEX2

!     This is a portable replacement for ishft, which behaves as the UM
!     expects it too when shifting more than 64 bits.













