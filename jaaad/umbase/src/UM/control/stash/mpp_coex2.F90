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


      SUBROUTINE MPP_COEX2(                                             &
     &  FIELD,M,ICOMP,N,IX,IY,NUM,ISC,OCO,RMDI,                         &
     &  LOCAL_LEVS,GLOBAL_LEVS,PE_FOR_LEVEL,LOCAL_LEVEL,                &
     &  first_gr_index,n_rows, icode, cmessage)

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
     &, GLOBAL_LEVS                                                     & 
                    ! IN: total number of levels being COEX'd

     &, first_gr_index                                                  &
                       ! IN: first global row index I will COEX
     &, n_rows         ! IN: number of rows I will COEX

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

! COMMON blocks and PARAMETERS
#include "parvars.h"

! Local variables

      INTEGER K,KK,JJ,II,IERR

      INTEGER                                                           &
     &  gri                                                             &
                   ! global row index counter
     &, local_row                                                       &
                   ! local row corresponding to gri
     &, local_lev                                                       &
                   ! local level corresponding to gri
     &, global_lev                                                      &
                   ! global level corresponding to gri
     &, proc                                                            &
                   ! PE corresponding to gri
     &, start_row                                                       &
                   ! first row in chunk to process
     &, last_row                                                        &
                   ! last row this PE is responsible for
     &, n_rows_in_chunk ! number of rows in chunk

      INTEGER                                                           &
     &  IST,ICX,JCX,JBIT,JLEN,NOB,NOC                                   &
     &, IC(IY),ICB(IY),ISTART(IY)                                       &
     &, NOP(IY,LOCAL_LEVS),NOP_PROC(IY)                                 &
     &, IBIT(IY,LOCAL_LEVS),IBIT_PROC                                   &
     &, JCOMP(IX,IY,LOCAL_LEVS),JCOMP_PROC(IX,IY)

      REAL                                                              &
     &  ACC,APREC                                                       &
     &, BASE(IY,LOCAL_LEVS),BASE_PROC                                   &
     &, FIELD_PROC(IX,IY)

      INTEGER ICODE
      CHARACTER*80 CMESSAGE
! Shmem addressing stuff
      INTEGER                                                           &
     &  nvars

      PARAMETER(nvars=3)

      INTEGER                                                           &
     &  address_list(nvars)                                             &
     &, remote_address_list(nvars)

      COMMON /coex_shmem_align/                                         &
     &  address_list

      INTEGER                                                           &
     &  remote_FIELD(IX,IY,GLOBAL_LEVS)                                 &
     &, remote_JCOMP(IX,IY,GLOBAL_LEVS)                                 &
     &, remote_NOP(IY,GLOBAL_LEVS)

      POINTER                                                           &
     & (ptr_remote_FIELD,remote_FIELD),                                 &
     & (ptr_remote_JCOMP,remote_JCOMP),                                 &
     & (ptr_remote_NOP,remote_NOP)

      EQUIVALENCE                                                       &
     & (ptr_remote_FIELD,remote_address_list(1)),                       &
     & (ptr_remote_JCOMP,remote_address_list(2)),                       &
     & (ptr_remote_NOP,remote_address_list(3))

!---End of Shmem addressing stuff

! Functions
#if defined(CRAY)
      INTEGER CRI2IBM
#else
      INTEGER IEEE2IBM
#endif

! Check not been called for expand (not yet supported)

      IF (.NOT. OCO) THEN
        WRITE(6,*) 'MPP_COEX called for expand'
        WRITE(6,*) 'Not yet supported.'
        GOTO 9999
      ENDIF

! Initialise output array (compressed data might not fill it)

      DO K=1,LOCAL_LEVS
        DO JJ=1,IY
          DO II=1,IX
            JCOMP(II,JJ,K)=0.0
          ENDDO
        ENDDO
      ENDDO

! Initialise SHMEM addressing stuff.

      address_list(1)=LOC(FIELD)
      address_list(2)=LOC(JCOMP)
      address_list(3)=LOC(NOP)

      CALL barrier()

! Do the COEX

      ACC=ISC
      APREC=2.0**ACC

      IF (n_rows  /=  0) THEN

        gri=first_gr_index              ! global row index
        global_lev=((gri-1)/IY)+1         ! global level number
        local_row=gri-(global_lev-1)*IY ! local row number
        local_lev=LOCAL_LEVEL(global_lev)  ! local level number
        proc=PE_FOR_LEVEL(global_lev) ! PE this level is stored on

! Get addresses of array on processor PROC

        CALL shmem_get(remote_address_list,address_list,                &
     &                 nvars,proc)

! Implicit relationship by EQUIVALENCE
!        ptr_remote_FIELD=remote_address_list(1)
!        ptr_remote_JCOMP=remote_address_list(2)
!        ptr_remote_NOP=remote_address_list(3)

! Loop over the rows I am responsible for, getting the data into
! FIELD_PROC from the remote FIELD array, COEXing it, and then
! sending the compressed data back.

        start_row=gri
        last_row=gri+n_rows-1

        DO WHILE (start_row  <=  last_row)

! Process a chunk of contiguous rows.
! Either do all the rows to the end of the level - or just up
! to last_row if this is less.

          n_rows_in_chunk=                                              &
     &      MIN(((global_lev*IY)-start_row+1),                          &
     &          (last_row-start_row+1))

          CALL shmem_get(FIELD_PROC,                                    &
     &                   remote_FIELD(1,local_row,local_lev),           &
     &                   IX*n_rows_in_chunk,proc)

!          DO JJ=start_row,start_row+n_rows_in_chunk-1
          DO JJ=1,n_rows_in_chunk

! DEPENDS ON: cmps
            CALL CMPS(IX,FIELD_PROC(1,JJ),JCOMP_PROC(2,JJ),             &
     &                NOP_PROC(JJ),APREC,                               &
     &   IBIT_PROC,BASE_PROC,RMDI,ICODE,CMESSAGE)


! Encode base value, number of bits and length into JCOMP

#if defined(CRAY)
            IERR=CRI2IBM(3,1,JCOMP_PROC(1,JJ),0, BASE_PROC,1,64,32)
            IERR=CRI2IBM(2,1,JCOMP_PROC(1,JJ),32,IBIT_PROC,1,64,16)
            IERR=CRI2IBM(2,1,JCOMP_PROC(1,JJ),48,NOP_PROC(JJ),1,64,16)
#else
            IERR=IEEE2IBM(3,1,JCOMP_PROC(1,JJ),0, BASE_PROC,1,64,32)
            IERR=IEEE2IBM(2,1,JCOMP_PROC(1,JJ),32,IBIT_PROC,1,64,16)
            IERR=IEEE2IBM(2,1,JCOMP_PROC(1,JJ),48,NOP_PROC(JJ),         &
     &                       1,64,16)
#endif

          ENDDO  ! JJ

! Send the compressed data back to the PE from whence it came

          CALL shmem_put(remote_JCOMP(1,local_row,local_lev),           &
     &                   JCOMP_PROC,                                    &
     &                   IX*n_rows_in_chunk,proc)

          CALL shmem_put(remote_NOP(local_row,local_lev),               &
     &                   NOP_PROC,                                      &
     &                   n_rows_in_chunk,proc)

! Increment the row counter, and move to the next level, unless
! we've now finished all our work

          start_row=start_row+n_rows_in_chunk

          IF (start_row  <=  last_row) THEN
            local_row=1
            global_lev=global_lev+1
            local_lev=LOCAL_LEVEL(global_lev)
            proc=PE_FOR_LEVEL(global_lev)

            CALL shmem_get(remote_address_list,address_list,            &
     &                     nvars,proc)

! Implicit relation from EQUIVALENCE statement
!            ptr_remote_FIELD=remote_address_list(1)
!            ptr_remote_JCOMP=remote_address_list(2)
!            ptr_remote_NOP=remote_address_list(3)

          ENDIF

        ENDDO ! DO WHILE over chunks

      ENDIF ! IF n_rows  /=  0

      CALL barrier()  ! wait for all my compressed data to arrive back

      IF (ICODE  /=  0) GOTO 9999
! Calculate positions in output array for packed row
! (First packed row starts at word 1, bit 31)

      DO KK=1,GLOBAL_LEVS

        IF (PE_FOR_LEVEL(KK)  ==  mype) THEN

          K=LOCAL_LEVEL(KK)

          IC(1)=2
          ICB(1)=-1
          ISTART(1)=5

          DO JJ=2,IY
            IF (MOD(NOP(JJ-1,K),2)  ==  1) THEN
              IC(JJ ) = IC(JJ-1) + NOP(JJ-1,K)/2 + 1
              ICB(JJ) = -ICB(JJ-1)
              IF (ICB(JJ) >  0) IC(JJ) = IC(JJ) + 1
            ELSE
              IC(JJ)  = IC(JJ-1) + (NOP(JJ-1,K)+1)/2 + 1
              ICB(JJ) = ICB(JJ-1)
            ENDIF
            ISTART(JJ)  = 5
            IF(ICB(JJ) == 1) ISTART(JJ) = 1
          ENDDO

! Move temporary array into output array

          DO JJ=1,IY
            NOB  = NOP(JJ,K)*4 + 8
            IST  = ISTART(JJ)
            ICX  = IC(JJ)
#if defined(CRAY)
            CALL STRMOV(JCOMP(1,JJ,K),1,NOB,ICOMP(ICX,K),IST)
#else
            CALL MOVEBYTES(JCOMP(1,JJ,K),1,NOB,ICOMP(ICX,K),IST)
#endif

          ENDDO

! Insert total length of this field

          NUM = IC(IY)*2 + NOP(IY,K)
          IF (ICB(IY) <  0) NUM = NUM + 1
#if defined(CRAY)
          IERR=CRI2IBM(2,1,ICOMP(1,K),0,NUM,1,64,32)
#else
          IERR=IEEE2IBM(2,1,ICOMP(1,K),0,NUM,1,64,32)
#endif

        ENDIF ! If this level is on this PE

      ENDDO ! K loop over local levels

 9999 CONTINUE

      RETURN
      END SUBROUTINE MPP_COEX2

#endif
#if defined(NEC)
#endif
#endif
