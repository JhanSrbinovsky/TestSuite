! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT 2000, Met Office, All Rights Reserved.
! Please refer to file $UMDIR/vn$VN/copyright.txt for further details
! *****************************COPYRIGHT*******************************
!
!
!
! Subroutine interface:
      SUBROUTINE READ_UNPACK(buf,isize,unpack_size,lookup,fixhd12,      &
#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)                \
 || defined(VAROPSVER) || defined(MAKEBC) || defined(VOMEXT)           \
 || defined(FRAMES)
     & EXPAND,                                                          &
#endif
     & icode,cmessage)

      IMPLICIT NONE
!
! Description: Unpacking codes from READFL1A & RDMULT1A is
!              combined into this new deck.  This enables
!              readin data to be unpacked for both UM and SX
!
!
! Subroutine Arguments:

      INTEGER                                                           &
     &  PACK_CODE                                                       &
     & ,FIXHD12                                                         &
                     !   IN : 12th element of fixed length header
     & ,LOOKUP(64)                                                      &
                     !   IN : LOOKUP header from dump
     & ,ISIZE                                                           &
                     !   IN :
     & ,UNPACK_SIZE                                                     &
                     !   IN : Size of data when unpacked
     & ,ICODE        !   OUT: return code
#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)                \
 || defined(VAROPSVER) || defined(MAKEBC) || defined(VOMEXT)           \
 || defined(FRAMES)
      INTEGER                                                           &
     &  EXPAND
#endif

      REAL                                                              &
     &  buf(unpack_size) ! INOUT: Holds field to be unpacked

      CHARACTER*80                                                      &
     &  CMESSAGE

!local variables
      INTEGER                                                           &
     & I                                                                &
     & ,num_unpack_values                                               &
                            ! used to detect error
     & ,len_full_word                                                   &
     & ,process_size

#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)                \
 || defined(VAROPSVER) || defined(MAKEBC) || defined(VOMEXT)           \
 || defined(FRAMES)
      INTEGER                                                           &
     &  dimx                                                            &
     &, dimy                                                            &
     &, idum
#endif

#include "cprintst.h"
      REAL                                                              &
     & tempbuf(unpack_size)                                             & 
                            ! Holds field while it's unpacked from buf
     & ,probuf(isize)                                                   &
     & ,extbuf(isize)


#include "clookadd.h"
#include "c_mdi.h"


! ----------------------------------------------------------------

      LEN_FULL_WORD=64
      PACK_CODE=MOD((LOOKUP(LBPACK)),10)

      IF (pack_code == 2) THEN ! 32 Bit CRAY packing
        IF (LOOKUP(DATA_TYPE)  ==  1) THEN ! if it's REAL
! DEPENDS ON: expand32b
          CALL EXPAND32B( LOOKUP(LBLREC) , BUF, FIXHD12 )
        ELSE ! not a real field, so not allowed to be compressed
          WRITE(6,*) 'READ_UNPACK : Error trying to uncompress field'
          WRITE(6,*) 'This field cannot be compressed/uncompressed ',   &
     &               'at it is type ',LOOKUP(DATA_TYPE)
          ICODE=200
          CMESSAGE='Failure uncompressing field'
          GOTO 9999
        ENDIF  ! if it's REAL

#if defined(UTILIO) || defined(FLDIO)

        num_unpack_values=LOOKUP(LBLREC)

#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)                \
 || defined(VAROPSVER) || defined(MAKEBC) || defined(VOMEXT)           \
 || defined(FRAMES)

      ELSE IF (pack_code == 1) THEN ! WGDOS packing
        IF (expand == 1) THEN  ! Field is to be expanded
          idum=0

! DEPENDS ON: coex
          CALL COEX(tempbuf,unpack_size,buf,isize,                      &
     &        dimx,dimy,idum,idum,.FALSE.,rmdi,len_full_word,           &
     &        icode,cmessage)

          DO I=1,UNPACK_SIZE
            buf(I)=tempbuf(I)
          ENDDO
          num_unpack_values=dimx*dimy
        ENDIF

      ELSEIF(pack_code == 4) then

        IF (expand == 1) then
          num_unpack_values=lookup(lbrow) * lookup(lbnpt)
          ! Enabling expansion of pack_code 4 files with extra data
          if(lookup(lbext) >  0) then

            ! process size not lbrow*lbnpt as data is packed
            process_size=lookup(lblrec)-lookup(lbext)

            ! decompose BUF into grid data array & extra data array
            do i=1,process_size
              probuf(i)=buf(i)
            enddo
            do i=1,lookup(lbext)
              extbuf(i)=buf(process_size+i)
            enddo

! DEPENDS ON: runlen_decode
            call runlen_decode(tempbuf,num_unpack_values,probuf,        &
     &                         process_size,rmdi,icode,cmessage)

            ! combine expanded grid data with extra data
            do i=1,num_unpack_values
              buf(i)=tempbuf(i)
            enddo


            do i=1,lookup(lbext)
              buf(num_unpack_values+i)=extbuf(i)
            enddo

          else   ! if no extra data

! DEPENDS ON: runlen_decode
          call runlen_decode(tempbuf,num_unpack_values,buf,             &
     &                        lookup(lblrec),rmdi,icode,cmessage)
          endif
        ENDIF
      ELSE IF (pack_code == 3) THEN ! GRIB packing
#if defined(NECSX6)
! DEPENDS ON: degrib
        CALL DEGRIB(buf,tempbuf,unpack_size,isize,                      &
     &      LOOKUP(1),rmdi,num_unpack_values,len_full_word)
        DO I=1,UNPACK_SIZE
          buf(I)=tempbuf(I)
        ENDDO
#else
        WRITE(6,*) 'Grib unpacking only supported on NEC SX6.'
! DEPENDS ON: ereport
        CALL EREPORT('READ_UNPACK', 1000,                               &
     &   'Grib unpacking only support on NEC SX6')

#endif

#endif
#endif

      ELSE IF (pack_code /= 0) THEN
        icode=6
        cmessage='READ_UNPACK: packing type not supported'
        RETURN
      ENDIF


#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)                \
 || defined(VAROPSVER) || defined(MAKEBC) || defined(VOMEXT)           \
 || defined(FRAMES)
      IF (pack_code /= 0 .AND.                                          &
     &    ((pack_code /= 1 .AND. pack_code /= 4).OR.expand == 1)) THEN
        IF (unpack_size <  num_unpack_values) THEN
          icode=7
          WRITE(6,*)'READ_UNPACK: UNPACK_SIZE = ',unpack_size,          &
     &      ' but NUM_UNPACK_VALUES = ',num_unpack_values
          cmessage='READ_UNPACK: workspace too small for Unpacked Data'
          RETURN
        ENDIF
      ENDIF

! Adjust the value of the data record length if data has been unpacked
        IF (((pack_code == 1.OR.pack_code == 4) .AND. expand == 1) .OR. &
     &      pack_code == 3) THEN
          IF (LOOKUP(lblrec) /= num_unpack_values) THEN
#if !defined(VAROPSVER)
            IF (lookup(lbext) >  0) THEN
              num_unpack_values=num_unpack_values+lookup(lbext)
            ENDIF
            if(printstatus >= prstatus_diag)then
              WRITE(6,*)'Record Length Changed from ',                  &
     &          LOOKUP(lblrec), 'to ',num_unpack_values
            endif
#endif
            LOOKUP(lblrec)=num_unpack_values
          ENDIF
        ENDIF


#endif


 9999 CONTINUE

      RETURN
      END SUBROUTINE READ_UNPACK
