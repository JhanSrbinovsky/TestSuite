#if defined(C84_1A) || defined(FLDIO)     || defined(FLDCALC)          \
 || defined(UTILIO) || defined(VAROPSVER) || defined(RECON)
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












      SUBROUTINE EXTRIN(ICOMP,IWORD,ISTART,NBIT,INUM,ISIGN)
!FPP$ NOCONCUR R

      Implicit None

!     Subroutine arguments
      INTEGER IWORD,ICOMP(IWORD),ISTART,NBIT,INUM(*),ISIGN

!     Local variables
      Integer :: Iscomp
      Integer :: Isnum
      Integer :: Num
      Integer :: Ibit

!
!     CRAY VERSION C.1.1  16/11/90  P J SMITH
!
!
!     SUBROUTINE TO EXTRACT AN INTEGER VALUE FROM NBITS OF ICOMP
!     ==========================================================
!
!     ICOMP  = ARRAY FROM WHICH AN NBIT INTEGER IS TO BE RETRIEVED
!     IWORD  = DIMESION OF ICOMP (ENABLES INTEGER TO CROSS WORDS)
!     ISTART = POSITION IN ICOMP WHERE AN NBIT INTEGER STARTS
!     NBIT   = NUMBER OF BITS IN WHICH INTEGER HAS BEEN STORED
!     INUM   = INTEGER EXTRACTED FROM ICOMP
!     ISIGN  = INDICATOR FOR STORAGE IN ICOMP: 0= +VE INTEGER
!                                              1= SIGNED INTEGER
!                                              2= SIGN BIT & +VE INT.
!
!     START POSITIONS IN ICOMP
!     |               |               |               |
!     |       ^       |       ^       |       ^       |       ^       ^
!     6666555555555544444444443333333333222222222211111111110000000000
!     3210987654321098765432109876543210987654321098765432109876543210
!
!     0000000001111111111222222222233333333334444444444555555555566666
!     1234567890123456789012345678901234567890123456789012345678901234
!
!
      INUM(1)= 0
      ISCOMP = 64 - ISTART
      ISNUM  = 64 - NBIT + 1
      NUM    = NBIT
!
!     MOVE INTEGER
!
      IF (ISIGN == 0) THEN
!
!         POSITIVE INTEGER
!
#if defined(CRAY)
          CALL MOVBIT(ICOMP,ISCOMP,NUM,INUM,ISNUM)
#else
          CALL MOVEBITS(ICOMP,ISCOMP,NUM,INUM,ISNUM)
#endif

      ELSEIF(ISIGN == 1) THEN
!
!         SIGNED INTEGER
!
!         SIGN
#if defined(CRAY)
          CALL MOVBIT(ICOMP,ISCOMP,1,INUM,1)
#else
          CALL MOVEBITS(ICOMP,ISCOMP,1,INUM,1)
#endif
          ISCOMP = ISCOMP + 1
          ISNUM  = ISNUM +1
          NUM    = NUM -1
!         INTEGER
#if defined(CRAY)
          CALL MOVBIT(ICOMP,ISCOMP,NUM,INUM,ISNUM)
#else
          CALL MOVEBITS(ICOMP,ISCOMP,NUM,INUM,ISNUM)
#endif
!         SET UNDIFINED IF INUM NEGATIVE
          IF (INUM(1) <  0) THEN
              DO IBIT=NUM,63
                  INUM(1) = IBSET(INUM(1),IBIT)
              END DO
          ENDIF

      ELSEIF(ISIGN == 2) THEN
!
!         SIGN BIT PLUS POSITIVE INTEGER
!
          ISCOMP = ISCOMP + 1
          ISNUM  = ISNUM +1
          NUM    = NUM -1
!         INTEGER
#if defined(CRAY)
          CALL MOVBIT(ICOMP,ISCOMP,NUM,INUM,ISNUM)
#else
          CALL MOVEBITS(ICOMP,ISCOMP,NUM,INUM,ISNUM)
#endif
!         SIGN
          IF (BTEST(ICOMP(1),ISTART)) INUM(1) = -INUM(1)
      ENDIF

      RETURN

      END SUBROUTINE EXTRIN

#if defined(T3E) && defined(MPP)



#endif
#if defined(NEC)
#endif
#endif
