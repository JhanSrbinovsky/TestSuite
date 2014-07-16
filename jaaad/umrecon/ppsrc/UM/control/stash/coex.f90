

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

      SUBROUTINE COEX(FIELD,M,ICOMP,N,IX,IY,NUM,ISC,OCO,RMDI,LWORD,     &
     &                ICODE,CMESSAGE)

      Implicit None

!     Subroutine arguments
      INTEGER N,ICOMP(N),M,IX,IY,NUM,ISC,LWORD
      INTEGER ICODE
      CHARACTER CMESSAGE*80
      REAL FIELD(M),RMDI
      LOGICAL OCO

!     Local variables
      Integer :: Ier
      Integer :: Ocode

!     External functions used




      Integer :: IEEE2IBM
      Integer :: IBM2IEEE


!
!     CRAY VERSION C.1.1  16/11/90  P J SMITH
!
!
!                     OCO=.TRUE.                 OCO=.FALSE.
!
!      FIELD   =  FIELD FOR COMPRESSING     RETURNED EXPANDED DATA
!          M   =  SIZE OF FIELD             SIZE OF FIELD
!      ICOMP   =  RETURNED COMPRESSED DATA  FIELD FOR EXPANSION
!          N   =  SIZE OF COMP                 -------
!         IX   =  X DIMENSION OF FIELD      RETURNED X DIMENSION
!         IY   =  Y DIMENSION OF FIELD      RETURNED Y DIMENSION
!        NUM   =  TOTAL NO. OF COMPRESSED      -------
!                 (32 BIT) WORDS RETURNED
!        ISC   =  ACCURACY IN POWER OF 2
!        OCO   =  .TRUE. FOR COMPRESSION    .FALSE. FOR EXPANSION
!
!        USERS OF THIS ROUTINE SHOULD ENSURE THAT THE ARRAY 'COMP'
!      WILL BE BIG ENOUGH TO HOLD THE COMPRESSED DATA.
!        IF INPUT ARRAY IS ONE DIMENSIONAL PUT IY=1 AND IX=DIMENSION,
!      WHEN USING WITH PRINTFILE FIELDS USE IX=192,IY=121  NOT IX=23232,
!      IY=1 THIS WILL MAKE THE COMPRESSION MORE EFFICIENT.
!      FOR MOST PRINTFILE FIELDS USING AN ACCURACY OF AN EIGHTH (ISC=-3)
!      A COMPRESSION FACTOR OF 4 CAN BE ACHIEVED. FOR DDDFFF FIELDS
!      USERS ARE ADVISED TO SPLIT THE FIELD INTO U AND V
!      COMPONENTS.
!
!     CRAY ROUTINE     64 BIT WORDS - CRAY FULLWORD
!
!    Check if an error has already been encountered, and get out
!    if it has.
      ocode = 0
      IF (icode  >   0) then
         goto 9999
      ELSE IF (icode  <   0)then
         ocode = icode
         icode = 0
      END IF

      IF (OCO) THEN
!     WRITE(6,*) 'CRAY COEX - PACKING'
!
!         COMPRESSION OF DATA
!
!
!         INSERT SCALE FACTOR, COLS AND ROWS INTO HEADER
!





          IER=IEEE2IBM(2,1,ICOMP(1),32,ISC,1,64,32)
          IER=IEEE2IBM(2,1,ICOMP(2),0,IX,1,64,16)
          IER=IEEE2IBM(2,1,ICOMP(2),16,IY,1,64,16)


!
!         CALL AUTOTASKED/DYNAMIC ARRAY PART OF COEX
!
! DEPENDS ON: coex2
          CALL COEX2(FIELD,M,ICOMP,N,IX,IY,NUM,ISC,OCO,RMDI,            &
     &               ICODE,CMESSAGE)
!
!         CHECK RETURN CODES
!
!         IFAIL=0
!         DO JJ = 1,IY
!             IF (IRC(JJ) /= 0) IFAIL=IFAIL+1
!         ENDDO
!         IF(IFAIL >  0)THEN
!         DO JJ = 1,IY
!             IF (IRC(JJ) /= 0) THEN
!              WRITE(6,*)'RETURN CODE',IRC(JJ),'FROM COEX FOR ROW',JJ
!             ENDIF
!         ENDDO
!         ENDIF
!
!         END OF COMPRESSION SECTION
!
      ELSE
!
!         EXPANSION SECTION
!
!
!         EXTRACT SCALE FACTOR, COLS AND ROWS FROM HEADER
!
          IER=IBM2IEEE(2,1,ICOMP(1),32,ISC,1,64,32)
          IER=IBM2IEEE(2,1,ICOMP(2),0,IX,1,64,16)
          IER=IBM2IEEE(2,1,ICOMP(2),16,IY,1,64,16)
!
!         CALL AUTOTASKED/DYNAMIC ARRAY PART OF COEX
!
! DEPENDS ON: coex2
          CALL COEX2(FIELD,M,ICOMP,N,IX,IY,NUM,ISC,OCO,RMDI,            &
     &               ICODE,CMESSAGE)
!
!         CHECK RETURN CODES (NOT YET IMPLEMENTED)
!
!         DO JJ = 1,IY
!             IF (IRC(JJ) /= 0) THEN
!                 WRITE(6,*)' NON-ZERO RETURN CODE FOR ROW ',JJ,IRC(JJ)
!             ENDIF
!         ENDDO
      END IF
!
!
      IF (icode  ==  0 .and. ocode  /=  0) then
         icode = ocode
      END IF

 9999 CONTINUE
      RETURN

      END SUBROUTINE COEX




!     This is a portable replacement for ishft, which behaves as the UM
!     expects it too when shifting more than 64 bits.













