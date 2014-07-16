#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines to read requested fields from a UM fieldsfile.


!=======================================================================

SUBROUTINE ReadFld( UMHdr,           &  ! in
                    LDecode,         &  ! in
                    PPHdrMod,        &  ! in
                    Field,           &  ! inout
                    ErrorStatus )       ! inout

! Description:
!   This subroutine reads a given pp-field from a UM fieldsfile.
!   The location of the field in the lookup should already be set within
!   Field, and memory allocated to store the data.
!
! Method:
!   The type of file is ascertained from the pp-header, and the field
!   read in.  If the field is packed and LDecode is true, the field is
!   decompressed.  The data type is checked for compatibility with the
!   program.  If the input field is of a type requiring header
!   modification and PPHdrMod is true, the appropriate modification is
!   done.
!
! Owner: Dave Robinson
!
! History:
! Version Date     Comment
! ------- ----     -------
! 1.0     02/05/03 Original Code.  Sara James
! 6.0     12/09/03 Code implemented into UM. Dave Robinson
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  UM_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
IMPLICIT None

! Subroutine Arguments:
TYPE(UM_Header_type), INTENT(IN) :: UMHdr    ! Fieldsfile header info
LOGICAL, INTENT(IN) :: LDecode               ! Unpack? indicator
LOGICAL, INTENT(IN) :: PPHdrMod              ! Mod of certain headers

TYPE(PP_Field_type), INTENT(INOUT) :: Field  ! Output field
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "ReadFld"
#include "c_mdi.h"

! Local Variables:
INTEGER :: NumCrayWords      ! no of values in an input field
INTEGER :: NumPts            ! The num of points in a data field
INTEGER :: i                 ! local counter
INTEGER :: DataLen           ! Length of a particular field
INTEGER :: InputAddr         ! Word Address in call SETPOS
INTEGER :: Addr              ! Address of a field in the data store
INTEGER :: PackType          ! Packing type N1 of LBPACK
INTEGER :: PackType_I        ! Packing type N1 of LBPACK in loop
INTEGER :: DataType          ! Input data type (integer, real etc.)
INTEGER :: Len_IO            ! Length of data written
INTEGER :: Year, Month, Date           ! Temp variables:
INTEGER :: Hour, Min, Sec, DayNo       !  used in date / time
INTEGER :: EndTimeDays, EndTimeSecs    !  conversion for
INTEGER :: DataTimeDays, DataTimeSecs  !  time-mean &
INTEGER :: FCRangeSecs                 !  accumulated fields
REAL :: Err_IO               ! Error Code
REAL :: AMDI                 ! Missing data indicator
LOGICAL :: Lcal360 = .FALSE. ! 30-day month indicator for time routines
CHARACTER(LEN=80) :: ErrMessage

INTEGER, ALLOCATABLE :: ITempData(:,:)

! End of header --------------------------------------------------------

IF( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

! Decode LBPACK code
PackType  = MOD(Field % Hdr % LBPack,10)

IF ( Field % Hdr % LBNRec == 0 ) THEN
  ! Reading a model type dump - no direct addressing only relative
  IF ( PackType == 2 ) THEN               ! Is the field WGDOS packed
    NumCrayWords = Field % Hdr % LBNRec/2 ! LBNRec is no. 32 bit words
  ELSE
    NumCrayWords = Field % Hdr % LBNRec
  END IF
  NumPts = Field % Hdr % LBLRec ! No of data points
  Addr   = UMHdr % FixHd(160)
  ! Find field posn
  DO i = 1, Field % LookupPos-1
    PackType_i = MOD(UMHdr % Lookup(i) % LBPack,10)
    IF ( PackType_i == 2 ) THEN             ! 32 Bit packed
      DataLen = UMHdr % Lookup(i) % LBLRec/2
    ELSE
      DataLen = UMHdr % Lookup(i) % LBLrec
    END IF
    Addr = Addr + DataLen
  ENDDO
  InputAddr = Addr - 1
ELSE
  ! Reading a PP type file
  NumCrayWords = Field % Hdr % LBLRec
  InputAddr    = Field % Hdr % DataPos
  NumPts       = Field % Hdr % LBUser6
END IF

! DEPENDS ON: setpos
CALL SETPOS( UMHdr % UnitNum, InputAddr, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, ErrorStatus, "Failure in SETPOS" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

IF ( .NOT. LDecode ) THEN
  ! Read data without unpacking or number conversion
! DEPENDS ON: buffin
  CALL BUFFIN( UMHdr % UnitNum, Field % RData, NumCrayWords, &
                              Len_IO, Err_IO )

ELSE
  ! Want real, unpacked data for calculations
  DataType = Field % Hdr % LBUser1
  IF ( DataType == 1 ) THEN

    ! Real field
! DEPENDS ON: buffin
    CALL BUFFIN( UMHdr % UnitNum, Field % RData, NumCrayWords, &
                                Len_IO, Err_IO )
    IF ( PackType > 0 ) THEN        ! Is the field packed?

! #if !defined(NODBG)
!       WRITE(6,*) "Unpacking real data"
! #endif

      AMDI = Field % Hdr % BMDI
      ! Compare with standard RMDI
      IF ( AMDI /= RMDI ) THEN
        WRITE(6,*)" WARNING non-standard MDI in use"
      END IF
! DEPENDS ON: unpackflds
      CALL UnPackFlds( PackType, Field % Hdr % NumCols, &
                       Field % Hdr % NumRows, Field,    &
                       ErrorStatus )
      IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
        CALL EReport( RoutineName, StatusWarning, &
                      "ReadFld - could not unpack input field" )
        ErrorStatus = StatusWarning
        GO TO 9999
      END IF
    END IF

  ELSE IF ( DataType == 2 ) THEN

    ! Integer field
    ALLOCATE( ITempData(Field % Hdr % NumCols,  &
                        Field % Hdr % NumRows) )
    IF ( PackType > 0 ) THEN
! DEPENDS ON: ereport
      CALL EReport( RoutineName, StatusWarning, &
                    "Cannot handle packed integer data fields" )
      ErrorStatus = StatusWarning
      GO TO 9999
    END IF
! DEPENDS ON: buffin
    CALL BUFFIN( UMHdr % UnitNum, ITempData, NumCrayWords, &
                                Len_IO, Err_IO )
    Field % RData = REAL( ITempData )
    DEALLOCATE( ITempData )
    Field % Hdr % LBUser1 = 1   ! Field is now real

  ELSE IF ( DataType == 3 ) THEN  ! Logical

    ! Logical Field
! DEPENDS ON: ereport
    CALL EReport( RoutineName, StatusWarning, &
                  "Cannot handle logical data fields" )
    ErrorStatus = StatusWarning
    GO TO 9999

  ELSE

    WRITE( ErrMessage, '(A24,I2)')   &
                  "Unrecognised data type: ", Field % Hdr % LBUser1
! DEPENDS ON: ereport
    CALL EReport( RoutineName, StatusWarning, ErrMessage )
    ErrorStatus = StatusWarning
    GO TO 9999

  END IF
END IF

! Perform optional header modifications
IF ( PPHdrMod ) THEN
  IF ( Field % Hdr % MO8Type == 58 ) THEN
   ! Modify MetO8 code for 1.5m max/min/instantaneous temperature
    IF      ( Field % Hdr % LBPROC == 4096 ) THEN
      WRITE(6,*) "Changing MetO8 code from 58 to 157"
      Field % Hdr % MO8Type = 157  ! Minimum temperature
    ELSE IF ( Field % Hdr % LBPROC == 8192 ) THEN
      WRITE(6,*) "Changing MetO8 code from 58 to 156"
      Field % Hdr % MO8Type = 156  ! Maximum temperature
    END IF
  END IF
  IF ( Field % Hdr % LBTIM /= 11 ) THEN
   ! Modify date/time for time-mean / accumulated fields
    IF ( Field % Hdr % ValidYear > 0 ) THEN
      WRITE(6,*) "Modifying date / time header"
     ! Move validity time to correct position in header
      Field % Hdr % ValidYear  = Field % Hdr % DataYear
      Field % Hdr % ValidMonth = Field % Hdr % DataMonth
      Field % Hdr % ValidDate  = Field % Hdr % DataDate
      Field % Hdr % ValidHour  = Field % Hdr % DataHour
      Field % Hdr % ValidMin   = Field % Hdr % DataMin
      Field % Hdr % ValidDayNo = Field % Hdr % DataDayNo
     ! Calculate data time : (validity time - FCRange)
      Year  = Field % Hdr % DataYear
      Month = Field % Hdr % DataMonth
      Date  = Field % Hdr % DataDate
      Hour  = Field % Hdr % DataHour
      Min   = Field % Hdr % DataMin
      Sec   = 0
      DayNo = Field % Hdr % DataDayNo
      FCRangeSecs = Field % Hdr % FCRange * 3600
! DEPENDS ON: time2sec
      CALL TIME2SEC( Year, Month, Date, Hour, Min, Sec,  &
                     0, 0, EndTimeDays, EndTimeSecs, Lcal360 )
! DEPENDS ON: time_df
      CALL TIME_DF( EndTimeDays,  EndTimeSecs, 0, FCRangeSecs,  &
                    DataTimeDays, DataTimeSecs )
! DEPENDS ON: sec2time
      CALL SEC2TIME( 0, 0, DataTimeDays, DataTimeSecs,   &
                     Year, Month, Date, Hour, Min, Sec, DayNo, Lcal360 )
     ! Put data time into correct position in header
      Field % Hdr % DataYear  = Year
      Field % Hdr % DataMonth = Month
      Field % Hdr % DataDate  = Date
      Field % Hdr % DataHour  = Hour
      Field % Hdr % DataMin   = Min
      Field % Hdr % DataDayNo = DayNo
    ELSE
! DEPENDS ON: ereport
      CALL EReport( RoutineName, StatusWarning,  &
                    "Missing temporal data - cannot recalculate data time" )
    END IF
  END IF
END IF

9999 CONTINUE

END SUBROUTINE ReadFld

!=======================================================================


#endif
