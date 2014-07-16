#if defined(VOMEXT) || defined(FRAMES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Reads in a number of fields from UM format file
!
! Subroutine Interface:

SUBROUTINE READFLDS_SERIAL (NFTIN, NUMBER_OF_FIELDS,              &
                            FIRST_FIELD, LOOKUP, LEN1_LOOKUP_ARG, &
                            D1, FIXHD,                            &
                            EXPAND, ICODE, CMESSAGE)

  IMPLICIT NONE
!
! Description:

!  Reads in NUMBER_OF_FIELDS fields from file on unit NFTIN,
!  starting at field number FIRST_FIELD. The data is returned
!  in the D1 array.
 
!  ==============================================================
!  ========    WARNING    WARNING   WARNING    ==================
!  ==============================================================
!  This routine is set up to read any number of fields but it only
!  works if number of fields is set to 1. Needs investigating why
!  it fails with more than one level.
!  ===============================================================
!  ===============================================================
 
! Current code owner: Glenn Greed
!
! Subroutine Arguments:

  INTEGER :: NFTIN               ! IN: unit number to read data from
  INTEGER :: NUMBER_OF_FIELDS    ! IN: number of fields to read in
  INTEGER :: FIRST_FIELD         ! IN: first field to read in
  INTEGER :: LEN1_LOOKUP_ARG     ! IN: first dimension of LOOKUP table
  INTEGER :: LOOKUP(LEN1_LOOKUP_ARG,*)  ! IN: lookup table starting
                                 !     at field 1
  INTEGER :: FIXHD(*)            ! IN: fixed length header

  INTEGER :: EXPAND        ! IN: (=1 if WGDOS or RLE packed data
                           !      is to be expanded)
  INTEGER :: ICODE         ! OUT: return code

  REAL :: D1(*)            ! OUT: array to return the data in

  CHARACTER*80 :: CMESSAGE ! OUT: Error message if ICODE <> 0

! COMMON blocks and PARAMETERs

#include "c_mdi.h"
#include "clookadd.h"

! Local variables

  INTEGER :: k                ! loop over fields to read in
  INTEGER :: pack_code        ! packing code for field
  INTEGER :: field_start      ! location of field on disk
  INTEGER :: data_size        ! number of words of data on disk
                              ! (including padding for WFIO)
  INTEGER :: data_read_size   ! number of words to read from disk
  INTEGER :: data_full_size   ! number of words after any unpacking
  INTEGER :: len_io           ! number of words read from disk
  INTEGER :: field_item       ! Item number of field
  INTEGER :: field_sect       ! Section number of field
  INTEGER :: field_model      ! Model ID of field
  INTEGER :: i                ! loop index

  REAL    :: a_io             ! Return code from BUFFIN
  
  CHARACTER(LEN=*), PARAMETER :: RoutineName = 'READFLDS_SERIAL'

!--------------------------------------------------------------

  IF (FIXHD(12) < 403) THEN
    WRITE (6,*) 'READFLDS_SERIAL: file created by UM version ', FIXHD(12)
    ICODE = 10
    CMESSAGE = 'READFLDS_SERIAL: Cannot read files from before vn4.3'
! DEPENDS ON: ereport
    CALL Ereport ( RoutineName, Icode, Cmessage )
  END IF


  DO K = FIRST_FIELD, FIRST_FIELD+NUMBER_OF_FIELDS-1

    field_item  = MOD(LOOKUP(42,K),1000)
    field_sect  = (LOOKUP(42,K)-field_item)/1000
    field_model = LOOKUP(45,K)
    pack_code   = MOD((LOOKUP(LBPACK,K)),10)

!   ---------------------------------------
!   Determine location of the field on disk
!   ---------------------------------------

    field_start = LOOKUP(LBEGIN,K) ! position of field in file
    IF (field_start <= 0) THEN
      WRITE (6,*) 'READFLDS_SERIAL: start address =',field_start
      ICODE = 20
      CMESSAGE = 'READFLDS_SERIAL: start address of field not given'
      CALL Ereport ( RoutineName, Icode, Cmessage )
    END IF

!   --------------------
!   Determine data sizes
!   --------------------

! DATA_SIZE : contains the number of words of data used to store the
! field on disk (needs to be halved if 32 bit packing has been used)

    IF (pack_code == 2) THEN
      data_size = (LOOKUP(LBLREC,K)+1)/2
    ELSE
      data_size = LOOKUP(LBLREC,K)
    END IF

! DATA_FULL_SIZE : is the number of words required to store the field
! in memory after any unpacking is done.

! This is to give buf the correct size in RDUNPCK, as
! buf will be the final expanded size of whole field
! including extra data
! WARNING LBEXT - may be -32768 MISSING VALUE !

    IF ((pack_code == 4) .and. (lookup(lbext, k) > 0)) THEN
      data_full_size = max(lookup(lbrow, k)*lookup(lbnpt, k)  &
           +lookup(lbext, k) ,lookup(lblrec,k))
    ELSE
      data_full_size = max(lookup(lbrow, k)*lookup(lbnpt, k)  &
           ,lookup(lblrec,k))
    END IF

    IF ( (lookup(lbrow,k) < 0) .or. (lookup(lbnpt,k) < 0) ) THEN
      data_full_size = lookup(lblrec,k)
    END IF

! DATA_READ_SIZE : contains the number of words to data that need to
! be read in for a field. Each field has extra words of dummy data
! added at the end to ensure each field starts on a disk sector
! boundary. The last field on a dump does not have these extra words
! added

    IF (K /= (FIRST_FIELD+NUMBER_OF_FIELDS-1)) THEN
      data_read_size = LOOKUP(LBNREC,K)
    ELSE
      data_read_size = data_size
    END IF

    IF (data_read_size < 0) THEN
      WRITE (6,*) 'READFLDS_SERIAL: number of words to read =',   &
           data_read_size
      ICODE = 30
      CMESSAGE = 'READFLDS_SERIAL: number of words to read not given'
      CALL Ereport ( RoutineName, Icode, Cmessage )
    END IF

! data_full_size needs to be at least as big as data_read_size since
! it is used to dimension the BUF array (in READ_MULTI ?)

    data_full_size = max(data_full_size, data_read_size)

!   -------------------------------------------
!   Move file pointer to the start of the field
!   -------------------------------------------

! DEPENDS ON: SETPOS
    CALL SETPOS (NFTIN, field_start, ICODE)
    IF (ICODE /= 0) THEN
      WRITE (6,*) 'READFLDS_SERIAL - SETPOS failed to move to ',  &
                   field_start,' on unit ',NFTIN
      WRITE (6,*) 'SETPOS returned error code ',ICODE
      ICODE = 40
      CMESSAGE = 'SETPOS failed in READFLDS_SERIAL. See output.'
      CALL Ereport ( RoutineName, Icode, Cmessage )
    END IF

!   ------------------------
!   Read the field from disk
!   ------------------------

    IF (pack_code == 2) THEN

!     Data is packed using CRAY 32 bit method - note 
!     that we need to read in 2*data_read_size 32 bit
!     words using BUFFIN32_single

      CALL BUFFIN32_single (NFTIN,D1,2*data_read_size,LEN_IO,ICODE)

!     And then halve LEN_IO to satisfy tests against
!     data_read_size

      LEN_IO = LEN_IO/2

    ELSE

!     Data is not 32-bit packed

      CALL BUFFIN_SINGLE (NFTIN,D1,data_read_size,LEN_IO,A_IO)

    END IF

!   --------------------------------
!   Check that data been read in OK?
!   --------------------------------

    IF ((A_IO /= -1.0) .OR. (LEN_IO /= data_read_size)) THEN
      WRITE (6,*) 'READFLDS_SERIAL : Error in call to BUFFIN_SINGLE'
      WRITE (6,*) 'LEN_IO : ',LEN_IO
      WRITE (6,*) 'A_IO : ',A_IO
      WRITE (6,*) 'Attempting to read field (Model,Section,Item) ',  &
                   field_model, field_sect, field_item
      WRITE (6,'(10(E10.5,X))') D1(1:100)
      ICODE = 50
      CMESSAGE = 'READFLDS_SERIAL: Failure reading in field'
      CALL Ereport ( RoutineName, Icode, Cmessage )
    END IF

!   ----------------------------------------------
!   If it's a compressed REAL field, expand it out
!   ----------------------------------------------

! DEPENDS ON: READ_UNPACK
    CALL READ_UNPACK (D1, data_read_size, data_full_size,  &
                      LOOKUP(1,K), FIXHD(12),              &
                      EXPAND,                              &
                      icode, cmessage)

    IF (ICODE /= 0) THEN
      WRITE(6,*)'READFLDS_SERIAL: Failure unpacking field'
      CALL Ereport ( RoutineName, Icode, Cmessage )
    END IF

  END DO ! K : loop over fields

  RETURN
END SUBROUTINE READFLDS_SERIAL
#endif
