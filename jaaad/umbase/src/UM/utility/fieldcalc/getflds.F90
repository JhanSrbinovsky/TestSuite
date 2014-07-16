#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines to read requested fields from a UM fieldsfile.

SUBROUTINE GetFlds ( NumFlds,     &  ! in
                     MaxFlds,     &  ! in
                     STCode,      &  ! in
                     MO8Level,    &  ! in
                     FCTime,      &  ! in
                     LBProc,      &  ! in
                     MinsPastHr,  &  ! in
                     Store,       &  ! in
                     PPHdrMod,    &  ! in
                     UMHdr,       &  ! in
                     Fields,      &  ! inout
                     ErrorStatus )   ! inout

! Description:
!   This subroutine reads a given number of fields, matching given
!   criteria from an open UM fieldsfile.  The fields are then store in
!   given locations in the 'Fields' array.
!
! Method:
!   See online documentation.
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
  PP_Field_type,          &
  UM_Header_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
IMPLICIT None

! Subroutine arguments:
INTEGER, INTENT(IN) :: NumFlds
INTEGER, INTENT(IN) :: MaxFlds
INTEGER, INTENT(IN) :: STCode  (NumFlds)
INTEGER, INTENT(IN) :: MO8Level(NumFlds)
INTEGER, INTENT(IN) :: FCTime  (NumFlds)
INTEGER, INTENT(IN) :: LBProc  (NumFlds)
INTEGER, INTENT(IN) :: Store   (NumFlds)
INTEGER, INTENT(IN) :: MinsPastHr (NumFlds)
LOGICAL, INTENT(IN) :: PPHdrMod
TYPE(UM_Header_type), INTENT(IN) :: UMHdr

TYPE(PP_Field_type), INTENT(INOUT) :: Fields(MaxFlds)
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "GetFlds"

! Local Variables:
INTEGER :: i                      ! local counter
INTEGER :: ifld
INTEGER :: ihdr
CHARACTER(LEN=80) :: ErrMessage
LOGICAL :: LDecode = .true.

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

!-------------------------------
! Loop through required fields
!-------------------------------
DO i = 1,NumFlds
  ifld = Store(i)     ! Array index where field should be stored
  ! See if correct field is already there (to reuse fields)
  IF ( (Fields(ifld) % Hdr % STCode   == STCode(i)  )   .AND. &
       (Fields(ifld) % Hdr % MO8Level == MO8Level(i))   .AND. &
       (Fields(ifld) % Hdr % FCRange  == FCTime(i)  )   .AND. &
       (LBProc(i)                     == -1         )   .AND. &
       (Fields(ifld) % Hdr % ValidMin == MinsPastHr(i)) .AND. &
       (ASSOCIATED( Fields(ifld) % RData ))  ) THEN
    CYCLE      ! Already exists - move on to next field
  END IF

  ! Clear space for field
  IF ( ASSOCIATED( Fields(ifld) % RData ) ) THEN
    DEALLOCATE( Fields(ifld) % RData )
    NULLIFY( Fields(ifld) % RData )
  END IF

  !----------------------------------------
  !  Search lookup for the required FIELD
  !----------------------------------------
  ! Search on LBTYP, LBLEV, LBFT
  Fields(ifld) % LookupPos = 0
  DO ihdr = 1,UMHdr % Len2Lookup
    IF ( (    STCode(i) == UMHdr % Lookup(ihdr) % STCode)   .AND. &
         (  MO8Level(i) == UMHdr % Lookup(ihdr) % MO8Level) .AND. &
         (    FCTime(i) == UMHdr % Lookup(ihdr) % FCRange)  .AND. &
         (MinsPastHr(i) == UMHdr % Lookup(ihdr) % ValidMin) ) THEN 
      IF ( (LBProc(i) == -1) .OR. &
           (LBProc(i) == UMHdr % Lookup(ihdr) % LBProc) ) THEN
        Fields(ifld) % LookupPos = ihdr
      END IF
    END IF
  END DO

  IF ( Fields(ifld) % LookupPos == 0 ) THEN
    WRITE( ErrMessage, '(A55,3I5)' )                                &
        "GetFlds: Field not found for STCode,MO8Level,FCRange: ",   &
        STCode(i), MO8Level(i), FCTime(i)
! DEPENDS ON: ereport
    CALL EReport( RoutineName, StatusWarning, ErrMessage )
    ErrorStatus = StatusWarning
    GOTO 9999
  END IF

  Fields(ifld) % Hdr = UMHdr % Lookup( Fields(ifld) % LookupPos )
  ! Allocate space for unpacked field
  ALLOCATE ( Fields(ifld) % RData( Fields(ifld) % Hdr % NumCols,  &
                                   Fields(ifld) % Hdr % NumRows ) )

  !----------------------------------------
  ! End of search.  Start retrieving data
  !----------------------------------------

! DEPENDS ON: readfld
  CALL ReadFld( UMHdr, LDecode, PPHdrMod, Fields(ifld),  &
                ErrorStatus )
  IF ( ErrorStatus /= StatusOK ) THEN
    GO TO 9999
  END IF

END DO

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE GetFlds

!=======================================================================


!=======================================================================


#endif
