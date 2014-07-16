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













#if defined(T3E) && defined(MPP)



#endif
#if defined(NEC)
      SUBROUTINE XPND_ALL(FIELD,ICOMP64,N,IX,IY,ISC,RMDI,               &
     &                    ICODE,CMESSAGE)

!     Better vectorized version of XPND expanding the whole field at
!     once

      IMPLICIT NONE

!     Subroutine arguments

      INTEGER, INTENT(IN)  :: n, ix, iy, isc
      INTEGER, INTENT(IN)  :: icomp64(n)
      REAL,    INTENT(IN)  :: rmdi
      REAL,    INTENT(OUT) :: field(ix,iy)
      INTEGER              :: icode
      CHARACTER            :: cmessage*80

!     Local variables

      INTEGER :: i, j, nshft, num, iword, ioff, imask, ival, mant, iexp
      INTEGER :: i1, i2, nbits_bmap
      INTEGER, DIMENSION(3*ix) :: itmp
      INTEGER, DIMENSION(ix)   :: idx, imap
      INTEGER, DIMENSION(iy)   :: istart, nop, ibase, nbits
      INTEGER, DIMENSION(iy*(2*ix+2)+4) :: icomp

      REAL                 :: aprec
      REAL, DIMENSION (iy) :: base
      REAL, DIMENSION (ix) :: tmp

      LOGICAL, DIMENSION(iy) :: obtzer, obtmin, obtmis, obtmap

      INTEGER, PARAMETER :: MASK16 = Z'FFFF'
      INTEGER, PARAMETER :: MASK32 = Z'FFFFFFFF'
      INTEGER, PARAMETER :: MASK_MANT_IBM = Z'00FFFFFF'
      INTEGER, PARAMETER :: MASK_EXPT_IBM = Z'7F000000'
      INTEGER, PARAMETER :: MASK_SIGN_IBM = Z'80000000'

      INTEGER, SAVE :: MASK_BITS(0:63), first
      DATA first /1/

      IF(first/=0) THEN
        DO i=0,63
          MASK_BITS(i) = ISHFT(1,63-i)
        ENDDO
        first = 0
      ENDIF

! Scale factor

      aprec = 2.**isc

! All lengths and alignments in WGDOS packing are for 32-bit words,
! so life gets much easier when we treat the packed data as 32-bit
! words.
! We split therefore the 64-bit compressed data into two 32 bit words

      num = ISHFT(icomp64(1),-32) ! Number of 32 bit words

      IF (num > SIZE(icomp)-2) THEN
        ICODE = 2
        CMESSAGE='COEX: Compressed data has too many elements'
        RETURN
      ENDIF

      DO i=1,(num+1)/2
        icomp(2*i-1) = IAND(ISHFT(icomp64(i),-32),MASK32)
        icomp(2*i)   = IAND(icomp64(i),MASK32)
      ENDDO
      ! The following word MUST be 0, it is used during decomposition!
      icomp(num+1) = 0
      icomp(num+2) = 0

! Get start word and length of every row

      istart(1) = 6
      nop(1) = IAND(icomp(5),MASK16)

      DO j=2,iy
        istart(j) = istart(j-1) + nop(j-1) + 2
        nop(j) = IAND(icomp(istart(j)-1),mask16)
        IF (istart(j)+nop(j)-1>num) THEN
          ICODE = 2
          CMESSAGE='COEX: Compressed data inconsistent'
          RETURN
        ENDIF
      ENDDO

! Get base (as a 32-bit IBM floating point number) and number of bits
! for every row and convert IBM floats to native floats
! The routine IBM2IEEE does a rather bad job, so we code it explicitly

      DO j=1,iy
        ibase(j) = icomp(istart(j)-2)
        nbits(j) = IAND(ISHFT(icomp(istart(j)-1),-16),mask16)
      ENDDO

      DO j=1,iy
         mant = IAND(ibase(j),MASK_MANT_IBM)
         iexp = ISHFT(IAND(ibase(j),MASK_EXPT_IBM),-24)-64-6
         base(j) = 16.0**iexp*mant
         if(IAND(ibase(j),MASK_SIGN_IBM) /= 0) base(j) = -base(j)
      ENDDO

! Check if bitmaps are used

      DO j=1,iy
        obtzer(j) = IAND(nbits(j),128) /= 0
        obtmin(j) = IAND(nbits(j),64)  /= 0
        obtmis(j) = IAND(nbits(j),32)  /= 0
        obtmap(j) = obtzer(j) .OR. obtmin(j) .OR. obtmis(j)
        nbits(j) = IAND(nbits(j),31)
      ENDDO


! Decode data row by row

      DO j=1,iy

        ! Care about bitmaps

        imap(:) = 1 ! Data present indicator

        nbits_bmap = 0
        IF(obtmis(j)) nbits_bmap = nbits_bmap + ix
        IF(obtmin(j)) nbits_bmap = nbits_bmap + ix
        IF(obtzer(j)) nbits_bmap = nbits_bmap + ix

        IF(nbits_bmap > 0) THEN
          iword = istart(j)
          DO i1=1,nbits_bmap,64
            ival  = IOR(ISHFT(icomp(iword),32),icomp(iword+1))
            iword = iword+2
            DO i2=0,MIN(nbits_bmap-i1,63)
              itmp(i1+i2) = MERGE(1,0,IAND(ival,MASK_BITS(i2))/=0)
            ENDDO
          ENDDO
          istart(j) = istart(j) + (nbits_bmap+31)/32
        ENDIF

        nbits_bmap = 0

        ! Extract missing data bitmap

        IF(obtmis(j)) THEN
          WHERE(itmp(nbits_bmap+1:nbits_bmap+ix)/=0)
            field(:,j) = rmdi
            imap (:) = 0
          END WHERE
          nbits_bmap = nbits_bmap + ix
        ENDIF

        ! Extract minimum value bitmap

        IF(obtmin(j)) THEN
          WHERE(itmp(nbits_bmap+1:nbits_bmap+ix)/=0)
            field(:,j) = base(j)
            imap (:) = 0
          END WHERE
          nbits_bmap = nbits_bmap + ix
        ENDIF

        ! Extract zero value bitmap

        IF(obtzer(j)) THEN
          WHERE(itmp(nbits_bmap+1:nbits_bmap+ix)==0)
            field(:,j) = 0.
            imap (:) = 0
          END WHERE
          nbits_bmap = nbits_bmap + ix
        ENDIF

        IF(nbits(j)==0) THEN

          ! All points in row have same value

          IF(obtmap(j)) THEN
            WHERE(imap(:)/=0) field(:,j) = base(j)
          ELSE
            field(:,j) = base(j)
          ENDIF

        ELSE

          ! Get number [and index] of values to decode

          IF(obtmap(j)) THEN
            num = 0
            DO i=1,ix
              IF(imap(i) /= 0) THEN
                num = num+1
                idx(num) = i
              ENDIF
            ENDDO
          ELSE
            num = ix
          ENDIF

          imask = ISHFT(1,nbits(j))-1

          ! Decode data

          DO i=1,num

            ! Bit offset to value:
            ioff  = (i-1)*nbits(j)

            ! Number of word in icomp which contains first bit:
            iword = ISHFT(ioff,-5)+istart(j)

            ! We load this word and the following into ival,
            ! this way we don't have to care if a word boundary
            ! is crossed. This requires that ival is a 64 bit word!
            ival  = IOR(ISHFT(icomp(iword),32),icomp(iword+1))

            ! Number of bits we have to shift to the right:
            nshft = 64 - IAND(ioff,31) - nbits(j)

            ! Normally we could now code:
            !   ival = ISHFT(ival,-nshft)
            ! but since vector-shift-by-vector is not
            ! implemented in H/W we have to do:

            if(IAND(nshft,32)/=0) ival = ISHFT(ival,-32)
            if(IAND(nshft,16)/=0) ival = ISHFT(ival,-16)
            if(IAND(nshft, 8)/=0) ival = ISHFT(ival, -8)
            if(IAND(nshft, 4)/=0) ival = ISHFT(ival, -4)
            if(IAND(nshft, 2)/=0) ival = ISHFT(ival, -2)
            if(IAND(nshft, 1)/=0) ival = ISHFT(ival, -1)

            ! Mask ival and calculate decoded value:
            ival = IAND(ival,imask)
            tmp(i) = ival*aprec + base(j)
          ENDDO

          ! Write decoded values to field

          IF(obtmap(j)) THEN
            field(idx(1:num),j) = tmp(1:num)
          ELSE
            field(:,j) = tmp(:)
          ENDIF

        ENDIF

      ENDDO

      END SUBROUTINE XPND_ALL
#endif
#endif
