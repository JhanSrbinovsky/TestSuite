#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(UTILIO) || defined(FLDIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Parallel UM interface to BUFFOUT
!
! Subroutine Interface:
      SUBROUTINE WRITE_MULTI(NFT,D1,ISIZE,LEN_IO,LOCAL_LEN,             &
     &                       LOOKUP,FIXHD12,COMPBUF,                    &
     &                       ADDR_INFO,ICODE,CMESSAGE)
      IMPLICIT NONE
!
! Description:
!  This routine provides an interface to BUFFOUT for the parallel
!  Unified Model. It is used where each process must write out a
!  local section of a global field.
!
! Method:
!  Each processor sends its local part of the global field to PE 0
!  which assembles all the parts, and then writes them to disk.
!  Fields which are compressed to land points are expanded before
!  sending to PE 0, PE 0 then compresses the global field before
!  writing it to disk.
!
! Current Code Owner: Paul Burton
!
! History:
!  Model    Date     Modification history from model version 3.5
!  version
!    3.5    5/1/95   New DECK created for the Parallel Unified
!                    Model. P.Burton + D.Salmond
!    4.1    18/3/96   Simplified communications    P.Burton
!    4.2    18/11/96  Added *CALL AMAXSIZE for IOVARS comdeck
!                     Added atmos_ prefix to landmask fields P.Burton
!    4.2    18/9/96   Modify send/receive maps and change args to
!                     alltoall for GCOM/GCG v1.1     P.Burton
!    4.2    18/10/96  New name for group of processors in scatter_field
!                     P.Burton
!    4.3    18/03/97  Added D1_ADDR argument, and rewrote field
!                     recognition tests. Now handles diagnostic
!                     fields too.                         P.Burton
!    4.4    13/06/97  Use GENERAL_GATHER_FIELD for collecting data
!                     and change decomposition depending on model
!                     type of field being read in            P.Burton
!    4.5    13/01/98  Replace SHMEM COMMON block with dynamic array
!                                                          P.Burton
!LL  5.0    14/9/99  Remove IOSTAT and added ICODE arguments
!LL                   General changes for 5.x dump format
!LL                                                       P.Burton
!LL   5.2  24/07/00  Add new argument to GENERAL_GATHER_FIELD  P.Burton
!     5.3  22/11/01  Enable MPP as the only option for
!                    small executables         E.Leung
!    5.5    02/08/00  Modification for parallelisation of WAM.
!                     Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!     5.5  07/04/03  For SX, a direct copying of PE0 field is done
!                    as single processor is used          E.Leung
!     5.5  05/02/03  Portability changes allowing for big_endian
!                    I/O on little_endian platforms.        P.Dando
!  6.0  17/09/03  Add def for new NEC opt section C96_1C. R Barnes
!
! Subroutine Arguments:

      INTEGER                                                           &
     &  NFT                                                             &
                     !   IN : FORTRAN unit number
     & ,ISIZE                                                           &
                     !   IN : no. of words to write out
     & ,LEN_IO                                                          &
                     !  OUT : no. of words written out
     & ,LOCAL_LEN                                                       &
                     !  OUT : size of the local field written out
     & ,LOOKUP(64)                                                      &
                     !   IN : LOOKUP header from dump
     & ,FIXHD12      !   IN : 12th. element of fixed length header
                     !        required for packing fields

! Required for dimensioning ADDR_INFO
#include "d1_addr.h"

      INTEGER                                                           &
     &  ADDR_INFO(D1_LIST_LEN)   ! IN addressing info about field
      REAL                                                              &
     &  D1(*)                                                           &
                     !   IN : Array to write out
     & ,COMPBUF(*)   !   IN : Workspace for compressing field

      INTEGER                                                           &
     &  ICODE        !  OUT : Return code
      CHARACTER*(80)                                                    &
     &  CMESSAGE     !  OUT : Error message

! Parameters and Common blocks

#include "clookadd.h"
#include "parvars.h"
#include "csmid.h"
#include "decomptp.h"

      REAL buf(ISIZE*2)  ! Buffer for holding data to be written.
!                        ! Factor of two is incase the data is
!                        ! packed (ISIZE is the size on disk)
!DIR$ CACHE_ALIGN buf
      INTEGER                                                           &
     &  I      ! loop counter

      REAL                                                              &
     &  IOSTAT ! return code from buffout_single

      INTEGER                                                           &
     &  orig_decomp                                                     &
                     ! decomposition on entry
     &, new_decomp   ! decomposition to change to

! ------------------------------------------------------------------
      LEN_IO=ISIZE
      LOCAL_LEN=0

      orig_decomp=current_decomp_type
      new_decomp=orig_decomp

#if defined(FLDOP) || defined(MERGE) || defined(PPTOANC)
      new_decomp=decomp_smexe
#else
      IF ((ADDR_INFO(d1_imodl)  ==  ATMOS_IM) .AND.                     &
     &    (orig_decomp  /=  decomp_standard_atmos)) THEN

        new_decomp=decomp_standard_atmos

      ELSEIF ((ADDR_INFO(d1_imodl)  ==  OCEAN_IM) .AND.                 &
     &        (ADDR_INFO(d1_object_type)  ==  prognostic) .AND.         &
     &        (orig_decomp  /=  decomp_standard_ocean)) THEN

        new_decomp=decomp_standard_ocean

      ELSEIF ((ADDR_INFO(d1_imodl)  ==  OCEAN_IM) .AND.                 &
     &        (ADDR_INFO(d1_object_type)  /=  prognostic) .AND.         &
     &        (orig_decomp  /=  decomp_nowrap_ocean)) THEN

        new_decomp=decomp_nowrap_ocean

      ELSEIF ((ADDR_INFO(d1_imodl)  ==  WAVE_IM) .AND.                  &
     &    (orig_decomp  /=  decomp_standard_wave)) THEN

        new_decomp=decomp_standard_wave

      ENDIF
#endif

      IF (new_decomp  /=  orig_decomp) THEN

! DEPENDS ON: change_decomposition
        CALL CHANGE_DECOMPOSITION(new_decomp,icode)

        IF (icode  /=  0) THEN
          WRITE(6,*) 'WTMULT : Error changing to decomposition ',       &
     &      new_decomp
          WRITE(6,*)                                                    &
     &      'Attempting to write field (Model,Section,Item) ',          &
     &      ADDR_INFO(d1_imodl),                                        &
     &      ADDR_INFO(d1_section),                                      &
     &      ADDR_INFO(d1_item)
          ICODE=1
          CMESSAGE='Failure changing decomposition'
          GOTO 9999
        ENDIF

      ENDIF

! Gather the field from the local D1 array to buf


#if defined(UTILIO) || defined(FLDIO)
      do i=1,LOOKUP(LBLREC)
        buf(i)=D1(i)
      enddo
#else
! DEPENDS ON: general_gather_field
      CALL GENERAL_GATHER_FIELD(                                        &
     &  D1,buf,LOCAL_LEN,LOOKUP(LBLREC),1,                              &
     &  ADDR_INFO,0,                                                    &
     &  ICODE,CMESSAGE)

        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'WTMULT : Call to GENERAL_GATHER_FIELD failed'
          WRITE(6,*) 'Return code was ',ICODE
          WRITE(6,*) 'Error message was ',CMESSAGE
          WRITE(6,*) 'Field number ',LOOKUP(ITEM_CODE)
          WRITE(6,*) 'Dimensions ',LOOKUP(LBNPT),' x ',LOOKUP(LBROW)
          WRITE(6,*) 'Grid type ',ADDR_INFO(d1_grid_type)
          WRITE(6,*) 'Field was not written out'

          ICODE=300
          CMESSAGE='Failure to gather field'
          GOTO 9999
        ENDIF
#endif

! ------------------------------------------------------------------
! ------------------------------------------------------------------
! And finally the code to write the global field in array buf
! out to disk.

      IF (mype  ==  0) THEN
!       Does this field need to be compressed?
        IF(MOD((LOOKUP(LBPACK)),10)  ==  2) THEN
          IF(LOOKUP(DATA_TYPE)  ==  1) THEN
! DEPENDS ON: pack21
            CALL PACK21(LOOKUP(LBLREC),buf,                             &
                        COMPBUF)
          ENDIF
        ELSE ! no compression required - just do a copy
          DO i=1,LOOKUP(LBLREC)
            COMPBUF(i)=buf(i)
          ENDDO
        ENDIF

! Now write out the global field

        IF(MOD((LOOKUP(LBPACK)),10)  ==  2) THEN
          IF(LOOKUP(DATA_TYPE)  ==  1) THEN
! Data is packed using CRAY 32 bit method - note that we need to write
! out 2*ISIZE 32 bit words using BUFFO32
            CALL buffo32_single(NFT,COMPBUF,2*ISIZE,LEN_IO,IOSTAT)
! And then halve LEN_IO to satisfy tests against ISIZE
            LEN_IO = LEN_IO/2
          ENDIF
        ELSE
! For non-packed data
          CALL buffout_single(NFT,COMPBUF,ISIZE,LEN_IO,IOSTAT)
        ENDIF
        IF ((IOSTAT  /=  -1.0) .OR. (LEN_IO  /=  ISIZE)) THEN
          WRITE(6,*) 'WTMULT : Error in call to BUFFOUT_SINGLE'
          WRITE(6,*) 'LEN_IO : ',LEN_IO
          WRITE(6,*) 'IOSTAT : ',IOSTAT
          WRITE(6,*)                                                    &
     &      'Attempting to read field (Model,Section,Item) ',           &
     &      ADDR_INFO(d1_imodl),                                        &
     &      ADDR_INFO(d1_section),                                      &
     &      ADDR_INFO(d1_item)
          ICODE=400
          CMESSAGE='Failure writing out field'
          GOTO 9999
        ENDIF

      ENDIF ! am I PE 0 ?


! If the field was compressed for writing on disk, we need to compress
! and expand the field in memory. This ensures the same field exists in
! memory that would exist if this dump was read back in.

      IF(MOD((LOOKUP(LBPACK)),10)  ==  2) THEN
        IF(LOOKUP(DATA_TYPE)  ==  1) THEN
! DEPENDS ON: pack21
          CALL PACK21(LOCAL_LEN,D1,COMPBUF)
! DEPENDS ON: expand21
          CALL EXPAND21(LOCAL_LEN,COMPBUF,D1)
        ENDIF
      ENDIF

      IF (new_decomp  /=  orig_decomp) THEN  ! change back

! DEPENDS ON: change_decomposition
        CALL CHANGE_DECOMPOSITION(orig_decomp,icode)

      ENDIF
 9999 CONTINUE  ! point to jump to if there is a failure

      RETURN
      END SUBROUTINE WRITE_MULTI
#endif
