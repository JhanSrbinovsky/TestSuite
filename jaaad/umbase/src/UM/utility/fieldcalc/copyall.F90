#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines for copying fields from one UM fieldsfile to another

SUBROUTINE CopyAll( NumFlds,     &  ! in
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
  StatusWarning,          &
  StatusFatal
IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumFlds
INTEGER, INTENT(IN) :: MaxFlds
INTEGER, INTENT(IN) :: STCode     (NumFlds)
INTEGER, INTENT(IN) :: MO8Level   (NumFlds)
INTEGER, INTENT(IN) :: FCTime     (NumFlds)
INTEGER, INTENT(IN) :: LBProc     (NumFlds)
INTEGER, INTENT(IN) :: MinsPastHr (NumFlds)
INTEGER, INTENT(IN) :: Store      (NumFlds)
REAL,    INTENT(IN) :: Factor     (NumFlds)
LOGICAL, INTENT(IN) :: PPHdrMod
TYPE(UM_Header_type), INTENT(IN) :: UMHdr_in

TYPE(UM_Header_type), INTENT(INOUT) :: UMHdr_out
TYPE(PP_Field_type), INTENT(INOUT) :: Fields(MaxFlds)
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "CopyAll"
CHARACTER(LEN=*), PARAMETER :: CPackType = "WGDOS"
REAL, PARAMETER :: VSmall = 1.e-9

! Local Variables:
INTEGER :: i, j
INTEGER :: NumLookups
INTEGER :: MaxFldSize
INTEGER :: PackType
INTEGER :: NumRows
INTEGER :: NumCols
INTEGER :: DataLen
INTEGER :: NumStore
INTEGER :: NumScl
INTEGER :: WFldCount
INTEGER :: SFldCount
INTEGER :: DFldCount
INTEGER :: ifld
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

! Allocate space for field
NumLookups = UMHdr_in % NumFlds
MaxFldSize = MAXVAL( UMHdr_in % Lookup(1:NumLookups) % NumCols ) * &
             MAXVAL( UMHdr_in % Lookup(1:NumLookups) % NumRows )
ALLOCATE( TempField % RData( MaxFldSize, 1 ) )
WRITE(6,*) "Copying ", NumLookups, " Fields"
NULLIFY( XpndField % RData )
NULLIFY( CompField % RData )

! Check there is space in the lookup table
IF ( NumLookups + UMHdr_out % NumFlds > UMHdr_out % Len2Lookup ) THEN
  WRITE(6,*) "Insufficient space allocated for Lookup Table"
  WRITE(6,*) "Increase the value of MaxFldsOut in main prog to correct"
! DEPENDS ON: ereport
  CALL EReport( RoutineName, StatusFatal, &
               "Insufficient space allocated for Lookup Table" )
END IF

NumStore  = 0
NumScl    = 0
WFldCount = 0
SFldCount = 0
DFldCount = 0

! Clear space for saved fields
DO j = 1, NumFlds
  IF ( Store(j) > 0 ) THEN
    NumStore = NumStore + 1
    IF ( ASSOCIATED( Fields(Store(j)) % RData ) ) THEN
      DEALLOCATE( Fields(Store(j)) % RData )
      NULLIFY( Fields(Store(j)) % RData )
    END IF
  END IF
  IF ( ABS(Factor(j)) > VSmall ) THEN
    NumScl = NumScl + 1
  END IF
END DO

WFldCount = 0
SFldCount = 0
DFldCount = 0

DO i = 1, NumLookups
  ! Read input field - do not unpack
  TempField % LookupPos = i
  TempField % Hdr = UMHdr_in % Lookup( TempField % LookupPos )
  LDecode = .FALSE.
! DEPENDS ON: readfld
  CALL ReadFld( UMHdr_in, LDecode, PPHdrMod, TempField, &
                ErrorStatus )
  IF ( ErrorStatus /= StatusOK ) THEN
    GO TO 9999
  END IF

  ! Check if input field is listed
  DO j = 1, NumFlds
    IF ( (STCode  (j) == TempField % Hdr % STCode)      .AND. &
         (MO8Level(j) == TempField % Hdr % MO8Level)    .AND. &
         (FCTime  (j) == TempField % Hdr % FCRange)     .AND. &
         (MinsPastHr(j) == TempField % Hdr % ValidMin)  .AND. &
        ((LBProc  (j) == TempField % Hdr % LBProc)      .OR.  &
         (LBProc  (j) == -1))                         ) THEN

      ifld = 0
      IF ( Store(j) > 0 ) THEN
        WFldCount = WFldCount + 1
        IF ( (LBProc(j) == -1) .OR. &
             (LBProc(j) == TempField % Hdr % LBProc) ) THEN
          ifld = Store(j)
        END IF
      END IF
      IF ( ifld > 0 ) THEN
        IF( ASSOCIATED( Fields(ifld) % RData ) ) THEN
          DFldCount = DFldCount + 1
          WRITE(6,*) "Field ", ifld, " is already found - ",  &
                     "storing first match"
          ifld = 0
        END IF
      END IF

      IF ( (ifld > 0) .OR. (ABS(Factor(j)) > VSmall) ) THEN
        ! Unpack field
        Xpndfield % LookupPos = TempField % LookupPos
        XpndField % Hdr       = TempField % Hdr
        PackType = MOD( XpndField % Hdr % LBPack,10 )
        PackAcc  = XpndField % Hdr % BACC
        NumCols  = XpndField % Hdr % NumCols
        NumRows  = XpndField % Hdr % NumRows
        DataLen  = XpndField % Hdr % LBLREC
        ALLOCATE( XpndField % RData( NumCols, NumRows ) )

        IF ( TempField % Hdr % LBUser1 /= 1 ) THEN  ! not real field
          LDecode = .TRUE.
! DEPENDS ON: readfld
          CALL ReadFld( UMHdr_in, LDecode, PPHdrMod, XpndField,   &
                        ErrorStatus )
        ELSE
          XpndField % RData = RESHAPE( SOURCE = TempField % RData, &
                                       SHAPE  = (/NumCols, NumRows/) )
          IF ( PackType > 0 ) THEN
! DEPENDS ON: unpackflds
            CALL UnPackFlds( PackType, NumCols, NumRows, XpndField, &
                             ErrorStatus )
          END IF
        END IF
        IF ( ErrorStatus /= StatusOK ) THEN
          WRITE(6,*) "Problem with field number ", j
        END IF
      END IF

      IF ( ifld > 0 ) THEN
        ! Keep field for later use
        Fields(ifld) % LookupPos = XpndField % LookupPos  ! = i
        Fields(ifld) % Hdr       = XpndField % Hdr
        IF ( ABS(Factor(j)) > VSmall ) THEN     ! Scaling of output reqd
          ALLOCATE( Fields(ifld) % RData( NumCols, NumRows ) )
          Fields(ifld) % RData = XpndField % RData
        ELSE
          Fields(ifld) % RData => XpndField % RData
          NULLIFY( XpndField % RData )
        END IF
      END IF

      IF ( ABS(Factor(j)) > VSmall ) THEN
        ! Scaling of output required
        SFldCount = SFldCount + 1
        XpndField % RData = Factor(j) * XpndField % RData
! DEPENDS ON: pack_single
        CALL Pack_Single( PackAcc, CPackType, XpndField, CompField, &
                          ErrorStatus )
        TempField % Hdr = CompField % Hdr
        TempField % RData( 1:SIZE( CompField % RData), :) =    &
                                                    CompField % RData
        DEALLOCATE( XpndField % RData )
        NULLIFY( XpndField % RData )
        DEALLOCATE( CompField % RData )
        NULLIFY( CompField % RData )
      END IF

    END IF
  END DO

  ! Write field to output file
! DEPENDS ON: fldout
  CALL FldOut( TempField, UMHdr_out, ErrorStatus )
  IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
    CALL EReport( RoutineName, StatusWarning, &
                  "Could not write field to file" )
  END IF

END DO

WRITE(6,*) WFldCount, " stored, out of ", NumStore, " requested"
IF ( DFldCount /= 0 ) THEN
  WRITE(6,*) "Warning: ", DFldCount, " duplicate fields found"
END IF
WRITE(6,*) SFldCount, " scaled, out of ", NumScl, " requested"
DEALLOCATE( TempField % RData )

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE CopyAll

!=======================================================================


#endif
