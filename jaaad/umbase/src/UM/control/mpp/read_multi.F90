#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(UTILIO) || defined(FLUXPROC)

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Parallel UM interface to BUFFIN
!
! Subroutine Interface:
      SUBROUTINE READ_MULTI(NFT,D1,ISIZE,UNPACK_SIZE,LEN_IO,LOCAL_LEN,  &
     &                      LOOKUP,FIXHD12,                             &
     &                      ADDR_INFO,N_LEVELS,                         &
#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)      \
 || defined(VAROPSVER) || defined(MAKEBC)
     &                      expand,                                     &
#endif
     &                      ICODE,CMESSAGE)
      IMPLICIT NONE
!
! Description:
!  This routine provides an interface to BUFFIN for the parallel
!  Unified Model. It is used where each process must read in a
!  local section of a global field.
!
! Method:
!  PE 0 reads in the global field, and then distributes the
!  relevant parts of it to each processor.
!  Fields compressed to land points are expanded by PE 0, and
!  recompressed after being received by the relevant processor.
!
! Current Code Owner: Paul Burton
!
! Subroutine Arguments:

      INTEGER                                                           &
     &  NFT                                                             &
                     !   IN : FORTRAN unit number
     & ,ISIZE                                                           &
                     !   IN : no. of words to be read in (global field)
     & ,UNPACK_SIZE                                                     &
                     !   IN : no. of words after any unpacking (global)
     & ,LEN_IO                                                          &
                     !  OUT : no. of words read in (global field)
     & ,LOCAL_LEN                                                       &
                     !  OUT : no. of words in local field
#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)      \
 || defined(VAROPSVER) || defined(MAKEBC)
     & ,expand                                                          &
                     !   IN : expansion for small exec
#endif
     & ,LOOKUP(64)                                                      &
                     !   IN : LOOKUP header from dump
     & ,FIXHD12                                                         &
                     !   IN : 12th element of fixed length header
     & ,N_LEVELS                                                        &
                     !   IN : Number of levels
     &, ICODE        !   OUT: Return code
! Required for dimensioning ADDR_INFO
#include "d1_addr.h"

      INTEGER                                                           &
     &  ADDR_INFO(D1_LIST_LEN)   ! IN addressing info about field

      CHARACTER*(80)                                                    &
     &  CMESSAGE     !  OUT : Error message

      REAL                                                              &
     &  D1(*)        !  OUT : Array to read data in to

! Parameters and Common blocks

#include "clookadd.h"
#include "parvars.h"


! Local variables

#include "cppxref.h"
#include "decomptp.h"
#include "csmid.h"
      INTEGER                                                           &
     &  info                                                            &
               ! GCOM return code
     &, I                                                               &
               ! loop counter
     &, map(N_LEVELS) ! Which processor holds which level

      INTEGER                                                           &
     &  orig_decomp                                                     &
                     ! decomposition on entry
     &, new_decomp   ! decomposition to change to
      REAL                                                              &
     &  IOSTAT  ! return code from BUFFIN

      REAL buf(UNPACK_SIZE)  ! buffer for reading the field into
!DIR$ CACHE_ALIGN buf

! ------------------------------------------------------------------

      LEN_IO=ISIZE
      LOCAL_LEN=0

      orig_decomp=current_decomp_type
      new_decomp=orig_decomp

#if defined(UTILIO)
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

      ENDIF
#endif

      IF (new_decomp  /=  orig_decomp) THEN

! DEPENDS ON: change_decomposition
        CALL CHANGE_DECOMPOSITION(new_decomp,icode)

        IF (icode  /=  0) THEN
          WRITE(6,*) 'RDMULT : Error changing to decomposition ',       &
     &      new_decomp
          WRITE(6,*)                                                    &
     &      'Attempting to read field (Model,Section,Item) ',           &
     &      ADDR_INFO(d1_imodl),                                        &
     &      ADDR_INFO(d1_section),                                      &
     &      ADDR_INFO(d1_item)
          ICODE=1
          CMESSAGE='Failure changing decomposition'
          GOTO 9999
        ENDIF

      ENDIF

! First thing to do is to read the field in to PE 0

      IF (mype  ==  0) THEN

        IF (MOD((LOOKUP(LBPACK)),10)  ==  2) THEN
! Data is packed using CRAY 32 bit method - note that we need to read
! in 2*ISIZE 32 bit words using BUFFIN32
          CALL buffin32_single(NFT,buf,2*ISIZE,LEN_IO,IOSTAT)
! And then halve LEN_IO to satisfy tests against ISIZE
          LEN_IO = LEN_IO/2
        ELSE
! For non-packed data
          CALL BUFFIN_SINGLE(NFT,buf,ISIZE,LEN_IO,IOSTAT)
        ENDIF

!       Has the data been read in OK?
        IF ((IOSTAT  /=  -1.0) .OR. (LEN_IO  /=  ISIZE)) THEN
          WRITE(6,*) 'RDMULT : Error in call to BUFFIN_SINGLE'
          WRITE(6,*) 'LEN_IO : ',LEN_IO
          WRITE(6,*) 'IOSTAT : ',IOSTAT
          WRITE(6,*)                                                    &
     &      'Attempting to read field (Model,Section,Item) ',           &
     &      ADDR_INFO(d1_imodl),                                        &
     &      ADDR_INFO(d1_section),                                      &
     &      ADDR_INFO(d1_item)
          ICODE=100
          CMESSAGE='Failure reading in field'
          GOTO 9999
        ENDIF

! If it's a compressed REAL field, expand it out
! DEPENDS ON: read_unpack
      CALL READ_UNPACK(buf,isize,unpack_size,lookup,fixhd12,            &
#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)                \
 || defined(VAROPSVER) || defined(MAKEBC)
     & EXPAND,                                                          &
#endif
     & icode,cmessage)

      IF (ICODE /= 0) THEN
        CMESSAGE='Failure unpacking field'
        GOTO 9999
      ENDIF
      ENDIF ! IF (mype  ==  0)

! Now we can distribute it out to the other processes

! For atmosphere zonal ozone fields - set to zonal grid type

      IF (( ADDR_INFO(d1_grid_type)  ==  ppx_atm_ozone) .AND.           &
     &    ( LOOKUP(LBNPT)  ==  1)) THEN

        ADDR_INFO(d1_grid_type)=ppx_atm_tzonal

      ENDIF
! Create map showing which processor contains which level.
! At present, all levels are stored on PE 0
      DO i=1,N_LEVELS
        map(i)=0
      ENDDO


! Now decompose the field in buf to the local D1 arrays

#if defined(UTILIO) || defined(FLUXPROC)
      do i=1,unpack_size
        d1(i)=buf(i)
      enddo
#else
! DEPENDS ON: general_scatter_field
      CALL GENERAL_SCATTER_FIELD(                                       &
     &  D1,buf,LOCAL_LEN,LOOKUP(LBLREC),N_LEVELS,                       &
     &  ADDR_INFO,map,                                                  &
     &  ICODE,CMESSAGE)

      IF (ICODE  ==  1) THEN
        WRITE(6,*) 'RDMULT : Call to GENERAL_SCATTER_FIELD failed'
        WRITE(6,*) 'Return code was ',ICODE
        WRITE(6,*) 'Error message was ',CMESSAGE
        WRITE(6,*) 'Field number ',LOOKUP(ITEM_CODE)
        WRITE(6,*) 'Grid type ',ADDR_INFO(d1_grid_type)
        WRITE(6,*) 'Dimensioned : ',                                    &
     &               LOOKUP(LBNPT),' x ',LOOKUP(LBROW)
        ICODE=300
        CMESSAGE='Failure decomposing field'
        GOTO 9999
      ENDIF

      CALL GC_GSYNC(NPROC, ICODE)
#endif
      IF (new_decomp  /=  orig_decomp) THEN  ! change back

! DEPENDS ON: change_decomposition
        CALL CHANGE_DECOMPOSITION(orig_decomp,icode)

      ENDIF

9999  CONTINUE

      RETURN
      END SUBROUTINE READ_MULTI
#endif
