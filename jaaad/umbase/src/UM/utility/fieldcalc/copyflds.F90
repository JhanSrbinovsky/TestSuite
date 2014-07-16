#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines for copying fields from one UM fieldsfile to another


!=======================================================================

SUBROUTINE CopyFlds( NumFlds,     &  ! in
                     MaxFlds,     &  ! in
                     STCode,      &  ! in
                     MO8Level,    &  ! in
                     FCTime,      &  ! in
                     LBProc,      &  ! in
                     MinsPastHr,  &  ! in
                     Factor,      &  ! in
                     Store,       &  ! in
                     PPHdrMod,    &  ! in
                     UMHdr_in,    &  ! in
                     UMHdr_out,   &  ! inout
                     Fields,      &  ! inout
                     ErrorStatus )   ! inout

! Description:
!
! Method:
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
  PP_Field_type,          &
  PP_Header_type,         &
  UM_Header_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumFlds
INTEGER, INTENT(IN) :: MaxFlds
INTEGER, INTENT(IN) :: STCode     (NumFlds)
INTEGER, INTENT(IN) :: MO8Level   (NumFlds)
INTEGER, INTENT(IN) :: FCTime     (NumFlds)
INTEGER, INTENT(IN) :: LBProc     (NumFlds)
INTEGER, INTENT(IN) :: Store      (NumFlds)
INTEGER, INTENT(IN) :: MinsPastHr (NumFlds)
REAL,    INTENT(IN) :: Factor     (NumFlds)
LOGICAL, INTENT(IN) :: PPHdrMod
TYPE(UM_Header_type), INTENT(IN) :: UMHdr_in

TYPE(UM_Header_type), INTENT(INOUT) :: UMHdr_out
TYPE(PP_Field_type), INTENT(INOUT) :: Fields(MaxFlds)

INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "CopyFlds"
REAL, PARAMETER :: VSmall = 1.e-9
CHARACTER(LEN=*), PARAMETER :: CPackType = "WGDOS"

! Local Variables:
INTEGER :: NumLookups
INTEGER :: MaxFldSize
INTEGER :: i
INTEGER :: ihdr
INTEGER :: ifld
INTEGER :: istore
INTEGER :: PackType
INTEGER :: NumRows
INTEGER :: NumCols
INTEGER :: DataLen
INTEGER :: NumStore
INTEGER :: NumScl
INTEGER :: CFldCount
INTEGER :: WFldCount
INTEGER :: SFldCount
REAL    :: PackAcc
CHARACTER(LEN=80) :: ErrMessage
LOGICAL :: LDecode
TYPE(PP_Field_type) :: TempField
TYPE(PP_Field_type) :: XpndField
TYPE(PP_Field_type) :: CompField

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

! Allocate temporary field
NumLookups = UMHdr_in % NumFlds
MaxFldSize = MAXVAL( UMHdr_in % Lookup(1:NumLookups) % NumCols ) * &
             MAXVAL( UMHdr_in % Lookup(1:NumLookups) % NumRows )
ALLOCATE( TempField % RData( MaxFldSize, 1) )
NumStore  = 0
NumScl    = 0
CFldCount = 0
WFldCount = 0
SFldCount = 0

!--------------------
! Loop through list
!--------------------
DO i = 1,NumFlds

  IF ( Store(i) > 0 ) THEN
    NumStore = NumStore + 1
  END IF
  IF ( ABS(Factor(i)) > VSmall ) THEN
    NumScl = NumScl + 1
  END IF

  istore = Store(i)
  IF ( (Store(i) > 0) .AND.  &
      ((MO8Level(i) == -1) .OR. (FCTime(i) == -1)) ) THEN
    WRITE(6,*) "Cannot store multiple fields in one position"
    WRITE(6,*) "Store = ", Store, " Level = ", MO8Level, &
               " FCTime = ", FCTime
    istore = 0
  END IF

  !----------------------------------------
  !  Search lookup for the required FIELD
  !----------------------------------------
  ! -1 is a wildcard for MO8Level and FCRange. STCode must be fixed.
  TempField % LookupPos = 0
  DO ihdr = 1, NumLookups
    IF ( STCode(i) == UMHdr_in % Lookup(ihdr) % STCode ) THEN
      IF( ((MO8Level(i) == UMHdr_in % Lookup(ihdr) % MO8Level)    .OR.  &
           (MO8Level(i) == -1)                                )   .AND. &
          ((FCTime(i)   == UMHdr_in % Lookup(ihdr) % FCRange)     .OR.  & 
           (FCTime(i)   == -1)                                )   .AND. & 
          (MinsPastHr(i) == UMHdr_in % Lookup(ihdr) % ValidMin)   .AND. & 
          ((LBProc(i) == UMHdr_in % Lookup(ihdr) % LBProc)        .OR.  &
           (LBProc(i)   == -1)                                ) ) THEN

        !---------------------------------------
        ! Found matching field.  Retrieve data
        !---------------------------------------
        CFldCount = CFldCount + 1
        TempField % LookupPos = ihdr
        TempField % Hdr = UMHdr_in % Lookup( TempField % LookupPos )
        ! Read field into temporary space.
        LDecode = .FALSE.
! DEPENDS ON: readfld
        CALL ReadFld( UMHdr_in, LDecode, PPHdrMod, TempField,   &
                      ErrorStatus )

        ifld = 0
        IF ( (LBProc(i) == -1) .OR. &
             (LBProc(i) == TempField % Hdr % LBProc) ) THEN
          ifld = istore
        END IF

        IF ( ErrorStatus /= StatusOK ) THEN
          WRITE( ErrMessage, '(A36,I4)' ) &
                   "Could not read input field position ", ihdr
! DEPENDS ON: ereport
          CALL EReport( RoutineName, StatusWarning, ErrMessage )
        ELSE

          IF ( (ifld > 0) .OR. (ABS(Factor(i)) > VSmall) ) THEN ! Unpack
            ! Clear space for field
            XpndField % LookupPos = TempField % LookupPos
            XpndField % Hdr       = TempField % Hdr
            PackType = MOD(XpndField % Hdr % LBPack,10)
            PackAcc  = XpndField % Hdr % BACC
            NumCols  = XpndField % Hdr % NumCols
            NumRows  = XpndField % Hdr % NumRows
            DataLen  = XpndField % Hdr % LBLREC
            ALLOCATE( XpndField % RData( NumCols, NumRows ) )

            IF ( XpndField % Hdr % LBUser1 /= 1 ) THEN  ! not real field
              LDecode = .TRUE.
! DEPENDS ON: readfld
              CALL ReadFld( UMHdr_in, LDecode, PPHdrMod, XpndField,   &
                            ErrorStatus )

            ELSE
              XpndField % RData = RESHAPE(SOURCE = TempField % RData,&
                                          SHAPE  = (/NumCols, NumRows/))
              IF ( PackType > 0 ) THEN  ! real packed field
! DEPENDS ON: unpackflds
                CALL UnPackFlds( PackType, NumCols, NumRows, &
                                 XpndField, ErrorStatus )
              END IF
            END IF
            IF ( ErrorStatus /= StatusOK ) THEN
              WRITE(6,*) "Problem with field number ", i
            END IF
          END IF

          IF ( ifld > 0 ) THEN
            ! Keep field for later use
            WFldCount = WFldCount + 1
            IF ( ASSOCIATED( Fields(ifld) % RData ) ) THEN
              DEALLOCATE( Fields(ifld) % RData )
            END IF
            Fields(ifld) % LookupPos = XpndField % LookupPos
            Fields(ifld) % Hdr       = XpndField % Hdr
            IF ( ABS(Factor(i)) > VSmall ) THEN ! Scaling of output reqd
              ALLOCATE( Fields(ifld) % RData( NumCols, NumRows ) )
              Fields(ifld) % RData = XpndField % RData
            ELSE
              Fields(ifld) % RData => XpndField % RData
              NULLIFY( XpndField % RData )
            END IF
          END IF

          IF ( ABS(Factor(i)) > VSmall ) THEN
            ! Scaling of output required
            SFldCount = SFldCount + 1
            XpndField % RData = Factor(i) * XpndField % RData
! DEPENDS ON: pack_single
            CALL Pack_Single( PackAcc, CPackType, XpndField, &
                              CompField, ErrorStatus )
            TempField % Hdr = CompField % Hdr
            TempField % RData( 1:SIZE(CompField % RData), :) = &
                                                    CompField % RData
            DEALLOCATE( XpndField % RData )
            NULLIFY( XpndField % RData )
            DEALLOCATE( CompField % RData )
            NULLIFY( CompField % RData )
          END IF

          ! Write field
! DEPENDS ON: fldout
          CALL FldOut( TempField, UMHdr_out, ErrorStatus )

        END IF

      END IF
    END IF

  END DO

  IF ( TempField % LookupPos == 0 ) THEN
    WRITE( ErrMessage, '(A46,3I5)' )                         &
        "Field not found for STCode, MO8Level, FCRange: ",   &
        STCode(i), MO8Level(i), FCTime(i)
    IF ( Store(i) > 0 ) THEN
! DEPENDS ON: ereport
      CALL EReport( RoutineName, StatusWarning, ErrMessage )
      ErrorStatus = StatusWarning
      GO TO 9999
    ELSE
      WRITE(6,*) ErrMessage
    END IF
  END IF

END DO

WRITE(6,*) CFldCount, " fields copied."
WRITE(6,*) WFldCount, " stored, out of ", NumStore, " requested."
WRITE(6,*) SFldCount, " scaled, out of ", NumScl, " requested."
DEALLOCATE( TempField % RData )

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE CopyFlds

#endif
