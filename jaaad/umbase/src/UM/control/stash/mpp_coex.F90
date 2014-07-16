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

      SUBROUTINE MPP_COEX(                                              &
     &  FIELD,M,ICOMP,N,IX,IY,NUM,ISC,OCO,RMDI,                         &
     &  LOCAL_LEVS,GLOBAL_LEVS,PE_FOR_LEVEL,LOCAL_LEVEL,                &
     &  ICODE,CMESSAGE)

      IMPLICIT NONE

      INTEGER                                                           &
     &  M                                                               &
              ! IN: horizontal size of field
     &, N                                                               &
              ! IN: size of compressed data
     &, IX                                                              &
              ! IN: X dimension of field
     &, IY                                                              &
              ! IN: Y dimension of field
     &, NUM                                                             &
              ! OUT: Total number of compressed (32 bit) words output
     &, ISC                                                             &
              ! IN: Accuracy in power of 2

     &, LOCAL_LEVS                                                      &
                    ! IN: number of levels held on this PE
     &, GLOBAL_LEVS ! IN: total number of levels being COEX'd

      LOGICAL                                                           &
     &  OCO   ! IN: .TRUE.=compress  .FALSE.=expand (unsupported)

      INTEGER                                                           &
     &  ICOMP(N,LOCAL_LEVS)                                             &
                            ! OUT: returned compressed data

     &, PE_FOR_LEVEL(GLOBAL_LEVS)                                       &
                                  ! IN: which PE holds each level
     &, LOCAL_LEVEL(GLOBAL_LEVS)  ! IN: local level number on that PE

      REAL                                                              &
     &  FIELD(IX,IY,LOCAL_LEVS)                                         &
                                ! IN: field to compress
     &, RMDI                    ! IN: real missing data indicator
      INTEGER ICODE,OCODE
      CHARACTER*80 CMESSAGE

! COMMON blocks and PARAMETERS
#include "parvars.h"

      INTEGER min_rows_per_pe  ! minimum number of rows per proc.
      PARAMETER(min_rows_per_pe=1)

! Local variables

      INTEGER K,JJ,IERR

      INTEGER                                                           &
     &  total_rows                                                      &
                       ! total number of rows to be COEXd
     &, rows_per_pe                                                     &
                       ! number of rows per PE when split up
     &, n_full_pes                                                      &
                       ! number of PEs doing rows_per_pe rows
!                      ! all the others do rows_per_pe-1 rows
     &, first_gr_index                                                  &
                       ! first global row index I will COEX
     &, n_rows         ! number of rows I will be COEXing

! Functions
#if defined(CRAY)
      INTEGER CRI2IBM
#else
      INTEGER IEEE2IBM
#endif

!    Check if an error has already been encountered, and get out
!    if it has.
      ocode = 0
      IF (icode  >   0) then
         goto 9999
      ELSE IF (icode  <   0)then
         ocode = icode
         icode = 0
      END IF

! Check not been called for expand (not yet supported)

      IF (.NOT. OCO) THEN
        WRITE(6,*) 'MPP_COEX called for expand'
        WRITE(6,*) 'Not yet supported.'
        GOTO 9999
      ENDIF

! Calculate scale factor,columns and rows.

#if defined(CRAY)
      IERR=CRI2IBM(2,1,ICOMP(1,1),32,ISC,1,64,32)
      IERR=CRI2IBM(2,1,ICOMP(2,1),0,IX,1,64,16)
      IERR=CRI2IBM(2,1,ICOMP(2,1),16,IY,1,64,16)
#else
      IERR=IEEE2IBM(2,1,ICOMP(1,1),32,ISC,1,64,32)
      IERR=IEEE2IBM(2,1,ICOMP(2,1),0,IX,1,64,16)
      IERR=IEEE2IBM(2,1,ICOMP(2,1),16,IY,1,64,16)
#endif

! Copy these values to other levels

      DO K=2,LOCAL_LEVS
        ICOMP(1,K)=ICOMP(1,1)
        ICOMP(2,K)=ICOMP(2,1)
      ENDDO

! Calculate mapping
! Distribute each row of every level to a given processor

      total_rows=IY*GLOBAL_LEVS

      rows_per_pe=((total_rows-1)/nproc)+1

      IF (rows_per_pe  >=  min_rows_per_pe) THEN

! However, this can give a few too many rows, so only n_full_pes
! PEs will do rows_per_pe, and all the rest will do rows_per_pe-1

        n_full_pes=nproc-(rows_per_pe*nproc - total_rows)

! Calculate the first global_row_index that I am responsible for
! and the number of rows I will COEX

        IF (mype  <=  n_full_pes) THEN
          first_gr_index=1+(mype*rows_per_pe)
        ELSE
          first_gr_index=1+(n_full_pes*rows_per_pe)+                    &
     &                   (mype-n_full_pes)*(rows_per_pe-1)
        ENDIF

        IF (mype  <   n_full_pes) THEN
          n_rows=rows_per_pe
        ELSE
          n_rows=rows_per_pe-1
        ENDIF

      ELSE ! we force some PE's to do min_rows_per_pe, others none

        IF ((mype*min_rows_per_pe)  >=  total_rows) THEN
          n_rows=0
          first_gr_index=1
        ELSE
          n_rows=MIN(min_rows_per_pe,                                   &
     &               (total_rows-(mype*min_rows_per_pe)))
          first_gr_index=1+(mype*min_rows_per_pe)
        ENDIF

      ENDIF

! DEPENDS ON: mpp_coex2
      CALL MPP_COEX2(                                                   &
     &  FIELD,M,ICOMP,N,IX,IY,NUM,ISC,OCO,RMDI,                         &
     &  LOCAL_LEVS,GLOBAL_LEVS,PE_FOR_LEVEL,LOCAL_LEVEL,                &
     &  first_gr_index,n_rows, icode, cmessage)

      IF (icode  ==  0 .and. ocode  /=  0) then
         icode = ocode
      END IF

 9999 CONTINUE

      RETURN
      END SUBROUTINE MPP_COEX


#endif
#if defined(NEC)
#endif
#endif
