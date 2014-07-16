#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program for post-processing of UM output

PROGRAM FieldCalc

! Description:
!   This program performs a range of functions required for the post-
!   processing of UM output.  It can rewrite some or all of the input
!   fields, calculate diagnostics using standard level and model level
!   fields and scale wind fields.
!   The functions it performs are dependent on the contents of the user-
!   defined namelist control file.
!
! Method:
!   After setting up all the input and output files, the program loops
!   through each entry in the namelist file.  Each entry contains an
!   "Action" field, which is checked against the list of possible
!   actions using an IF construct.  Possible actions include reading
!   input fields into memory, using generic subroutines (e.g. Sum, Dif)
!   to manipulate fields already in memory, calculating diagnostics from
!   fields which may or may not be in memory and writing fields out to
!   file.
!   Various checks and constructs have been added to reduce IO and
!   memory overheads.  See online documentation for more details.
!
! Input Files:
!   Input files are defined using the following UNIX environment
!   variables:-
!   FIELDCALC_NL - Namelist file used to control the functions performed
!                  be Fieldcalc.  Essential for the running of the
!                  program.
!   UNIT20   - UM Fieldsfile containing standard level fields at the
!              required forecast ranges.  Essential for most
!              diagnostics.  Failure to open causes a fatal error.
!   UNIT21   - UM Fieldsfile containing model level fields at the
!              required forecast ranges.  Essential for most aviation
!              diagnostics.  Failure to open causes a fatal error.
!   UNIT22   - UM Fieldsfile containing standard level fields at
!              forecast ranges 6-12hrs earlier.  Failure to open will
!              result in a non-fatal error, but all accumulation
!              diagnostics will be unavailable.
!   OROGFILE - UM Fieldsfile containing the model orography field.
!              Failure to open will result in a non-fatal error, but
!              some aviation diagnostics will be unavailable.
!
! Output Files:
!   Output files are defined using the following UNIX environment
!   variables:-
!   UNIT30    - UM Fieldsfile.  Should contain only standard level flds.
!
! Owner: Dave Robinson
!
! History:
! Version Date     Comment
! ------- ----     -------
! 1.0     02/05/03 Original Code.   Sara James
! 6.0     12/09/03 Code implemented into UM. Dave Robinson
! 6.1     02/11/04 Increase CA_NumLevs_Gl from 22 to 24. Dave Robinson
! 6.1     09/06/04 Add new action TOPBASE and call to subroutine
!                  TOPBASE. Dave Robinson
! 6.2     28/07/05 Increase MaxFlds & MaxFldsOut. Increase length for
!                  Namelist filename. D.Robinson
! 6.2     16/01/06 Enable FieldCalc work with 38 or 50 model levels.
!                  D.Robinson
! 6.4     16/01/07 Allow action COPYFLDS to copy fields from the
!                  STDPREV stream as well as MOD and STD
!                  R Sempers (frpz)
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type,          &
  UM_Header_type
USE FldCodes_mod, ONLY:   &
  ST_Urho,   ST_Vrho,     &
  ST_Htheta,              &
  ST_BlkCld, ST_ConCld,   &
  ST_CClBP,  ST_CClTP,    &
  ST_LCClBP, ST_LCClTP,   &
  ST_CPNRT,               &
  ST_Prho,   ST_Ptheta,   &
  ST_GWSU,   ST_GWSV,     &
  ST_Ttheta, ST_Pstar,    &
  ST_Ustd,   ST_Vstd,     &

! Dust
  ST_Dust1,  ST_Dust2,    &
  ST_Dust3,  ST_Dust4,    &
  ST_Dust5,  ST_Dust6,    &
! End Dust

  ST_Orog,   LV_Surface
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning,          &
  StatusFatal,            &
  EndofFile
IMPLICIT None

! Local constants:

CHARACTER(LEN=*), PARAMETER :: ProgName = "FieldCalc"

INTEGER, PARAMETER :: MaxFlds    = 400
INTEGER, PARAMETER :: MaxFldsOut = 10000 ! Max no. fields to be written
INTEGER, PARAMETER :: MaxWrite   = 20
INTEGER, PARAMETER :: NLUnit     = 81

! Internal Field positions: diagnostics requiring model levels only
INTEGER, PARAMETER :: StdPosn  = 20         ! These constants mark
INTEGER, PARAMETER :: UPosn    = 40         ! certain areas of the
INTEGER, PARAMETER :: VPosn    = 80         ! Fields array for certain
INTEGER, PARAMETER :: PPosn    = 120        ! variables.  This makes it
INTEGER, PARAMETER :: TPosn    = 160        ! easier to reuse model
INTEGER, PARAMETER :: ZPosn    = 200        ! fields.

! Dust
INTEGER, PARAMETER :: DST1Posn = 231
INTEGER, PARAMETER :: DST2Posn = 242
INTEGER, PARAMETER :: DST3Posn = 253
INTEGER, PARAMETER :: DST4Posn = 264
INTEGER, PARAMETER :: DST5Posn = 275
INTEGER, PARAMETER :: DST6Posn = 286
! End Dust

INTEGER :: UM_NumLevs  ! Number of model levels
INTEGER, PARAMETER :: Global = 0

! List of model levels that Fieldcalc will work for
INTEGER, PARAMETER :: Model_Levels_38 = 38
INTEGER, PARAMETER :: Model_Levels_50 = 50
INTEGER, PARAMETER :: Model_Levels_70 = 70

! Independent of number of model levels
INTEGER, PARAMETER :: MX_Zsea_Upr   = 20000 !   Upr lim if no.levs /= 38

! Variables for model level ranges required to derive diagnostics.

INTEGER :: TP_NumLevs_Gl   ! TropHeight Fields:
INTEGER :: TP_ZeroLev_Gl   !
INTEGER :: IT_NumLevs_Gl   ! Isotherm Fields:
INTEGER :: IT_ZeroLev_Gl   !
INTEGER :: CT_NumLevs_Gl   ! Contrail Fields:
INTEGER :: CT_ZeroLev_Gl   !
INTEGER :: MX_NumLevs_Gl   ! MaxWind Fields:
INTEGER :: MX_ZeroLev_Gl   !
INTEGER :: CA_NumLevs_Gl   ! CAT Fields:
INTEGER :: CA_ZeroLev_Gl   !
INTEGER :: MW_NumLevs_Gl   ! MtnStress Fields:
INTEGER :: MW_ZeroLev_Gl   !

INTEGER :: WT_NumLevs_Gl   ! WAFC CAT turb:
INTEGER :: WT_ZeroLev_Gl   !
INTEGER :: IC_NumLevs_Gl   ! Icing Fields (original alg):
INTEGER :: IC_ZeroLev_Gl   !
INTEGER :: LI_NumLevs_Gl   ! Icing Fields:
INTEGER :: LI_ZeroLev_Gl   !
INTEGER :: ICT_NumLevs_Gl  ! In (Layer) cloud turb Fields:
INTEGER :: ICT_ZeroLev_Gl  !
INTEGER :: CB_NumLevs_Gl   ! CB Fields:
INTEGER :: CB_ZeroLev_Gl   !

! Number of standard levels required to derive diagnostics.

INTEGER, PARAMETER :: CA_NumStd     = 6   !  CAT Fields
INTEGER, PARAMETER :: MW_NumStd     = 2   !  MtnStress Fields

INTEGER, PARAMETER :: WT_NumStd     = 10  !  WAFC CAT Turb (Shear)

! Note that if the value of WTMW_NumStd is changed, the values set for
! MW_PRef must also be updated, since it's dimension is derived from
! WTMW_NumStd
INTEGER, PARAMETER :: WTMW_NumStd   = 8   !  WAFC CAT Turb (MW)

INTEGER, PARAMETER :: ICT_NumStd    = 5   !  In (Layer) Cloud Turb
INTEGER, PARAMETER :: CB_NumStd     = 4   !  CB Fields


! Additional variables for WAFC_TURB action
INTEGER, PARAMETER :: MW_n_PRef = WTMW_NumStd / 2 ! Number of pressure levels
                                                  !   where MW predictor required

REAL, PARAMETER   :: MW_PRef(MW_n_PRef) =    &
                     ( / 30000.0, 25000.0, 20000.0, 15000.0 / )
                                                  ! Array of pressure levels where
                                                  !   MW predictor required (Pa)

! Dust

! These values are OK for 38 and 50 level models, since the model level
! heights are the same at this level. There is a check at the beginning
! of the action which makes sure either 38 or 50 levels are used.
INTEGER, PARAMETER :: DT_NumLevs_Gl = 10  ! for dust concs
INTEGER, PARAMETER :: DT_ZeroLev_Gl = 0   ! read in model levels 1-10
! End Dust


! Local variables:

INTEGER :: P_Levels                    ! Number of model levels
INTEGER :: GridType                    ! Horizontal Grid Type
INTEGER :: NumFlds                     ! No. fields to perform action on
INTEGER :: NumLevs                     ! Number of Levels
INTEGER :: ZeroLev                     ! Level below lowest required
INTEGER :: i                           ! Loop counter
INTEGER :: ErrorStatus                 ! Program status monitor
INTEGER :: ReadStatus                  ! Namelist read status
REAL :: Zsea, n
CHARACTER(LEN=9)  :: Action     = "DUMMY"
CHARACTER(LEN=7)  :: LevType    = "STD"
CHARACTER(LEN=6)  :: PackType   = "NONE"
CHARACTER(LEN=100):: NLFileName = ""
CHARACTER(LEN=80) :: ErrMessage = ""
CHARACTER(LEN=10) :: RunEnv  = ""
LOGICAL :: PPHdrMod = .FALSE.
TYPE(PP_Field_type)  :: OrogField      ! Orography PP-Field
TYPE(UM_Header_type) :: OrogHdr_in     ! UM Headers:  orography,
TYPE(UM_Header_type) :: StdLevHdr_in   !   standard level input
TYPE(UM_Header_type) :: ModLevHdr_in   !   model level input
TYPE(UM_Header_type) :: StdPrvHdr_in   !   previous std level input
TYPE(UM_Header_type) :: StdLevHdr_out  !   output

INTEGER :: Source  (MaxFlds)           ! Position in Fields array
INTEGER :: Store   (MaxFlds)           ! Position in Fields array
INTEGER :: STCode  (MaxFlds)           ! STASH code     (LBUSER(4))
INTEGER :: MO8Level(MaxFlds)           ! MetO8 Level    (LBLEV)
INTEGER :: FCTime  (MaxFlds)           ! Forecast Range (LBFT)
INTEGER :: LBProc  (MaxFlds)           ! Variable type  (LBPROC)
INTEGER :: USource (MaxFlds)           ! Model level source field
INTEGER :: VSource (MaxFlds)           !   postions in Fields array.
INTEGER :: PSource (MaxFlds)           !   Used in conjuction with UPosn
INTEGER :: TSource (MaxFlds)           !   etc. to minimise IO.
INTEGER :: ZSource (MaxFlds)
INTEGER :: MinsPastHr(MaxFlds)         ! Minutes past whole hour for COPYFLDS
! Dust
INTEGER :: D1Source (MaxFlds)
INTEGER :: D2Source (MaxFlds)
INTEGER :: D3Source (MaxFlds)
INTEGER :: D4Source (MaxFlds)
INTEGER :: D5Source (MaxFlds)
INTEGER :: D6Source (MaxFlds)
! End Dust


REAL    :: Factor  (MaxFlds)           ! Input for diagnostics
REAL    :: PackAcc (MaxFlds)           ! Packing Accuracy

! Additional local variables for icing actions
INTEGER :: NumLyrs                     ! No. layers to output fields for
INTEGER :: LyrMO8L (MaxFlds)           ! LBLEV values assigned to layers
REAL    :: LyrLwrB (MaxFlds)           ! Layer lower boundaris (BRLEV)
REAL    :: LyrUprB (MaxFlds)           ! Layer upper boundaries (BLEV)

TYPE(PP_Header_type) :: NewPPHdr(MaxWrite)
TYPE(PP_Field_type)  :: Fields(MaxFlds)! Input and output fields

! cntl_io required for um_sector_size
#include "cntl_io.h"

Character (Len=80)   :: Cmessage = ' '   ! Error Message
Character (Len=8)    :: c_um_sector_size ! Char string for Env Var

NAMELIST / PPHdrModNL / &      ! Header modification flag namelist
  PPHdrMod
NAMELIST / CommandNL /  &      ! Command namelist - for program control
  Action,               &
  NumFlds,              &
  Source,               &
  Store,                &
  STCode,               &
  MO8Level,             &
  FCTime,               &
  LBProc,               &
  LevType,              &
  PackType,             &
  PackAcc,              &
  MinsPastHr,           &
  Factor,               &
  NumLyrs,              &
  LyrMO8L,              &
  LyrLwrB,              &
  LyrUprB

NAMELIST / NewPPHdrNL / &      ! PP Header namelist - for writing fields
  NewPPHdr

! End of Header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( ProgName, 1 )

! DEPENDS ON: initprintstatus
CALL InitPrintStatus

!-----------------------------------------------------------------------
! 1 PROGRAM SETUP
!-----------------------------------------------------------------------

DO i=1,MaxFlds
  NULLIFY( Fields(i) % RData )
  Fields(i) % ArrayPos  = i
  Fields(i) % LookupPos = 0
END DO
NULLIFY( OrogField % RData )
OrogField % ArrayPos = 9999

!--------------------------------------------
! 1.1 Initialise UM_SECTOR_SIZE from Env Var
!--------------------------------------------

Call FORT_GET_ENV('UM_SECTOR_SIZE',14, c_um_sector_size,8, ErrorStatus)

If ( ErrorStatus /=  0) THEN
  Write(Cmessage,*) 'UM_SECTOR_SIZE has not been set: Default = 2048'
! DEPENDS ON: ereport
  Call Ereport( ProgName, StatusWarning, Cmessage )
  um_sector_size = 2048
Else
  Read( c_um_sector_size, '(I4)') um_sector_size
End If

!-------------------------
! 1.2 Open Namelist File
!-------------------------
CALL FORT_GET_ENV( "FIELDCALC_NL", 12, NLFileName, 100, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( ProgName, StatusFatal, &   ! Cannot run without namelist
                "Cannot read env FIELDCALC_NL" )
END IF
WRITE(6,*) "Namelist File : ", NLFileName
OPEN( UNIT   = NLUnit,        &
      ACCESS = "SEQUENTIAL",  &
      ACTION = "READ",        &
      FILE   = NLFileName,    &
      FORM   = "FORMATTED",   &
      IOSTAT = ErrorStatus,   &
      STATUS = "OLD" )
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( ProgName, StatusFatal, &   ! Cannot run without namelist
               "Cannot open Namelist file" )
END IF

!-----------------------------
! 1.3 Initialise Input Files
!-----------------------------
StdLevHdr_in % UnitNum = 20
StdLevHdr_in % FileNameEnv = "UNIT20"
ModLevHdr_in % UnitNum = 21
ModLevHdr_in % FileNameEnv = "UNIT21"
StdPrvHdr_in % UnitNum = 22
StdPrvHdr_in % FileNameEnv = "UNIT22"

! DEPENDS ON: read_umhdr
CALL Read_UMHdr( StdLevHdr_in, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( ProgName, StatusFatal,     &    ! Fatal error - no input
               "Error Reading Standard Level Header" )
END IF

! DEPENDS ON: read_umhdr
CALL Read_UMHdr( ModLevHdr_in, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( ProgName, StatusFatal,     &    ! Fatal error - no input
               "Error Reading Model Level Header" )
END IF

! DEPENDS ON: read_umhdr
CALL Read_UMHdr( StdPrvHdr_in, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( ProgName, StatusWarning,   &    ! Non-fatal error
               "Error Reading Previous Standard Level Header" )
  ErrorStatus = StatusOK
  NULLIFY( StdPrvHdr_in % Lookup )
END IF

! Get number of model levels from Integer Constants, Word 8
P_Levels = ModLevHdr_in % IntC(8)

! Get grid type from Fixed Header, Word 4
GridType = ModLevHdr_in % FixHd(4)

!---------------------------------------------------------------
! 1.3.1 Initialise model level ranges now that P_Levels is known
!---------------------------------------------------------------

Select Case ( P_Levels )

  Case ( Model_Levels_38 )

    UM_NumLevs    = 38    ! 38 Model Levels

    TP_NumLevs_Gl = 27    ! TropHeight Fields:
    TP_ZeroLev_Gl = 5     !   Model Levels 6-32
    IT_NumLevs_Gl = 32    ! Isotherm Fields:
    IT_ZeroLev_Gl = 0     !   Model Levels 1-32
    CT_NumLevs_Gl = 30    ! Contrail Fields:
    CT_ZeroLev_Gl = 0     !   Model Levels 1-30
    MX_NumLevs_Gl = 32    ! MaxWind Fields:
    MX_ZeroLev_Gl = 0     !   Model Levels 1-32
    CA_NumLevs_Gl = 24    ! CAT Fields:
    CA_ZeroLev_Gl = 5     !   Model Levels 6-29
    MW_NumLevs_Gl = 15    ! MtnStress Fields:
    MW_ZeroLev_Gl = 13    !   Model Levels 14-28

    WT_NumLevs_Gl = 32    ! WAFC CAT turb ...
    WT_ZeroLev_Gl = 0     !   Model Levels 1-32
    IC_NumLevs_Gl = 24    ! Icing Fields (original alg):
    IC_ZeroLev_Gl = 0     !   Model Levels 1-24
    LI_NumLevs_Gl = 33    ! Icing Fields:
    LI_ZeroLev_Gl = 0     !   Model Levels 1-33
    ICT_NumLevs_Gl= 24    ! In (Layer) cloud turb Fields:
    ICT_ZeroLev_Gl= 0     !   Model Levels 1-24
    CB_NumLevs_Gl = 32    ! CB Fields:
    CB_ZeroLev_Gl = 0     !   Model Levels 1-32

  Case ( Model_Levels_50 )

    UM_NumLevs    = 50    ! 50 Model Levels

    TP_NumLevs_Gl = 28    ! TropHeight Fields:
    TP_ZeroLev_Gl = 5     !   Model Levels 6-33
    IT_NumLevs_Gl = 33    ! Isotherm Fields:
    IT_ZeroLev_Gl = 0     !   Model Levels 1-33
    CT_NumLevs_Gl = 30    ! Contrail Fields:
    CT_ZeroLev_Gl = 0     !   Model Levels 1-30
    MX_NumLevs_Gl = 33    ! MaxWind Fields:
    MX_ZeroLev_Gl = 0     !   Model Levels 1-33
    CA_NumLevs_Gl = 24    ! CAT Fields:
    CA_ZeroLev_Gl = 5     !   Model Levels 6-29
    MW_NumLevs_Gl = 15    ! MtnStress Fields:
    MW_ZeroLev_Gl = 13    !   Model Levels 14-28

    WT_NumLevs_Gl = 33    ! WAFC CAT turb ...
    WT_ZeroLev_Gl = 0     !   Model Levels 1-33
    IC_NumLevs_Gl = 24    ! Icing Fields (original alg):
    IC_ZeroLev_Gl = 0     !   Model Levels 1-24
    LI_NumLevs_Gl = 33    ! Icing Fields:
    LI_ZeroLev_Gl = 0     !   Model Levels 1-33
    ICT_NumLevs_Gl= 24    ! In (Layer) cloud turb Fields:
    ICT_ZeroLev_Gl= 0     !   Model Levels 1-24
    CB_NumLevs_Gl = 33    ! CB Fields:
    CB_ZeroLev_Gl = 0     !   Model Levels 1-33

  Case ( Model_Levels_70 )

    UM_NumLevs    = 70    ! 70 Model Levels

!   Note that TP_NumLevs_Gl -> CB_ZeroLev_Gl as for 38 and 50
!   levels are not defined here. Fieldcalc will not derive
!   diagnostics for 70 model levels yet.

!   Defining 70 levels here will enable Fieldcalc
!   be used for filtering 70 level Fieldsfiles
!   for the UK4 model in the operational suite.

  Case Default

! DEPENDS ON: ereport
    CALL EReport( ProgName, StatusFatal,   &  ! Fatal error - no output
                 "Number of model levels not catered for." )

End Select

Write (6,*) ' '
Write (6,*) ' FieldCalc has been set up for a ',UM_NumLevs,    &
            ' level model.'
Write (6,*) ' '

!---------------------------
! 1.4 Read Orography Field
!---------------------------
OrogHdr_in % UnitNum = 23
OrogHdr_in % FileNameEnv = "OROGFILE"
! DEPENDS ON: read_umhdr
CALL Read_UMHdr( OrogHdr_in, ErrorStatus )
IF ( ErrorStatus == StatusOK ) THEN
  NumFlds       =       1
  STCode    (1) = ST_Orog
  MO8Level  (1) = LV_Surface
  FCTime    (1) =       0
  LBProc    (1) =      -1
  MinsPastHr(1) =       0
  Store     (1) =       1

! The argument for MinsPastHr has been set to 0 manually because the
! namelist hasn't been read yet.

! DEPENDS ON: getflds
  CALL GetFlds( NumFlds, MaxFlds, STCode, MO8Level, FCTime,  &
                LBProc, MinsPastHr, Store, PPHdrMod, OrogHdr_in,     &
                OrogField, ErrorStatus )
END IF
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( ProgName, StatusWarning,               &
               "Error Reading Orography Field - &
               &some diagnostics will be unavailable" )
  NULLIFY( OrogField % RData )
END IF
! DEPENDS ON: file_close
CALL File_Close ( OrogHdr_in % UnitNum,                &
                  OrogHdr_in % FileNameEnv,            &
                  LEN_TRIM(OrogHdr_in % FileNameEnv),  &
                  0, 0, ErrorStatus )

!-----------------------------
! 1.5 Initialise output file
!-----------------------------
StdLevHdr_out % UnitNum = 30
StdLevHdr_out % FileNameEnv = "UNIT30"
! DEPENDS ON: new_umhdr
CALL New_UMHdr( StdLevHdr_in, MaxFldsOut, &
                StdLevHdr_out, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( ProgName, StatusFatal,    &    ! Fatal error - no output
               "Error Setting up Output File." )
END IF

!--------------------------
! 1.6 Get other variables
!--------------------------
READ( UNIT = NLUnit,       &
       NML = PPHdrModNL,   &
    IOSTAT = ReadStatus )
IF ( ReadStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( ProgName, ReadStatus, &
               "Error found in PPHdrModNL" )
END IF
WRITE(6,*) "Header modification = ", PPHdrMod
WRITE(6,*) "Global model = ", (GridType == Global)

!-----------------------------------------------------------------------
! 2 MAIN DIAGNOSTIC PROCESSING
!-----------------------------------------------------------------------

DO WHILE ( Action /= "END" )

 !----------------------------
 ! 2.1 Read Command Namelist
 !----------------------------
  NumFlds   = 0
  Source    = Store
  Store(:)  = 0
  LBProc(:) = -1
  MinsPastHr(:) = 00
  Factor(:) = 0.0
  LevType   = "STD"
  PackType  = "NONE"
  READ( UNIT = NLUnit,       &
         NML = CommandNL,    &
      IOSTAT = ReadStatus )
  IF      ( ReadStatus == EndofFile ) THEN
    WRITE(6,*) "Unexpectedly reached end of Namelist file - quitting"
    GO TO 9999
  ELSE IF ( ReadStatus /= StatusOK ) THEN
    ErrorStatus = ReadStatus
! DEPENDS ON: ereport
    CALL EReport( ProgName, ReadStatus, &
                 "Error found in CommandNL" )
  END IF

  WRITE(6,'(A8,A9,A2,I3,A1)') "Action: ", Action, " (", NumFlds, ")"

  IF      ( Action == "CLEAR_ERR" ) THEN
   ! 2.1.1 Clear previous errors
    ErrorStatus = StatusOK
    ErrMessage = ""

  ELSE IF ( Action == "END" )       THEN
   ! 2.1.2 End of namelist
    ErrorStatus = StatusOK
    WRITE(6,*) "All actions completed successfully."

  ELSE IF ( ErrorStatus /= StatusOK ) THEN
   ! 2.1.3 Previous Error - skip actions until error cleared
    WRITE(6,*) "Command ignored due to previous error."
    IF ( Action == "WRITEFLDS" ) THEN
      ! Read the NewPPHdrNL that will come up after WRITEFLDS
      READ( UNIT = NLUnit,      &
            NML  = NewPPHdrNL )     ! Don't worry about errors
    END IF

 !---------------------------------
 ! 2.2 Actions for performing I/O
 !---------------------------------
  ELSE IF ( Action == "COPYALL" )   THEN

   ! 2.2.1 Rewrite entire contents of input file without modification
   !       Save & unpack specified fields
   !       Scale specified fields before writing
    IF      ( LevType == "STD" ) THEN
! DEPENDS ON: copyall
      CALL CopyAll( NumFlds, MaxFlds, STCode, MO8Level, FCTime,      &
                    LBProc, MinsPastHr, Factor, Store, PPHdrMod,     &
                    StdLevHdr_in, StdLevHdr_out, Fields, ErrorStatus )
    ELSE IF ( LevType == "MOD" ) THEN
! DEPENDS ON: copyall
      CALL CopyAll( NumFlds, MaxFlds, STCode, MO8Level, FCTime,      &
                    LBProc, MinsPastHr, Factor, Store, PPHdrMod,     &
                    ModLevHdr_in, StdLevHdr_out, Fields, ErrorStatus )
    ELSE
! DEPENDS ON: ereport
      CALL EReport( ProgName, StatusWarning, &
                    "COPYALL : Unknown level type " // LevType )
      ErrorStatus = StatusWarning
    END IF

  ELSE IF ( Action == "COPYFLDS" )  THEN

   ! 2.2.2 Rewrite specified fields
   !       Save & unpack those with Store /= 0
    IF      ( LevType == "STD" ) THEN
! DEPENDS ON: copyflds
      CALL CopyFlds ( NumFlds, MaxFlds, STCode, MO8Level, FCTime,     &
                      LBProc, MinsPastHr, Factor, Store, PPHdrMod,    &
                      StdLevHdr_in, StdLevHdr_out, Fields,            &
                      ErrorStatus )
    ELSE IF ( LevType == "MOD" ) THEN
! DEPENDS ON: copyflds
      CALL CopyFlds ( NumFlds, MaxFlds, STCode, MO8Level, FCTime,     &
                      LBProc, MinsPastHr, Factor, Store, PPHdrMod,    &
                      ModLevHdr_in, StdLevHdr_out, Fields,            &
                      ErrorStatus )
    ELSE IF ( LevType == "STDPREV" ) THEN
! DEPENDS ON: copyflds
      CALL CopyFlds ( NumFlds, MaxFlds, STCode, MO8Level, FCTime,     &
                      LBProc, MinsPastHr, Factor, Store, PPHdrMod,    &
                      StdPrvHdr_in, StdLevHdr_out, Fields,            &
                      ErrorStatus )
    ELSE
! DEPENDS ON: ereport
      CALL EReport( ProgName, StatusWarning, &
                    "COPYFLDS : Unknown level type " // LevType )
      ErrorStatus = StatusWarning
    END IF

  ELSE IF ( Action == "GETFLDS" )   THEN

   ! 2.2.3 Read fields
    IF      ( LevType == "STD" )     THEN
! DEPENDS ON: getflds
      CALL GetFlds ( NumFlds, MaxFlds, STCode, MO8Level, FCTime,        &
                     LBProc, MinsPastHr, Store, PPHdrMod, StdLevHdr_in, &
                     Fields, ErrorStatus )
    ELSE IF ( LevType == "MOD" )     THEN
! DEPENDS ON: getflds
      CALL GetFlds ( NumFlds, MaxFlds, STCode, MO8Level, FCTime,        &
                     LBProc, MinsPastHr, Store, PPHdrMod, ModLevHdr_in, &
                     Fields, ErrorStatus )
    ELSE IF ( LevType == "STDPREV" ) THEN
      IF ( ASSOCIATED( StdPrvHdr_in % Lookup ) ) THEN
! DEPENDS ON: getflds
        CALL GetFlds ( NumFlds, MaxFlds, STCode, MO8Level, FCTime,        &
                       LBProc, MinsPastHr, Store, PPHdrMod, StdPrvHdr_in, &
                       Fields, ErrorStatus )
      ELSE
        WRITE(6,*) "GETFLDS : Cannot read Previous Standard Level file"
        ErrorStatus = StatusWarning
      END IF
    ELSE
      WRITE(6,*) "GETFLDS : Unknown level type ", LevType
      ErrorStatus = StatusWarning
    END IF

  ELSE IF ( Action == "PACK" )      THEN

   ! 2.2.4 Pack fields
! DEPENDS ON: packflds
    CALL PackFlds( NumFlds, MaxFlds, PackType, Source, PackAcc, &
                   Fields, ErrorStatus )

  ELSE IF ( Action == "WRITEFLDS" ) THEN

   ! 2.2.5 Write fields

! #if !defined(NODBG)
!     WRITE(6,*) "Getting new header info"
! #endif

    NewPPHdr(1:NumFlds) = Fields(Source(1:NumFlds)) % Hdr
    READ( UNIT = NLUnit,     &
           NML = NewPPHdrNL, &
        IOSTAT = ReadStatus )
    IF ( ReadStatus /= StatusOK ) THEN
      ErrorStatus = ReadStatus
! DEPENDS ON: ereport
      CALL EReport( ProgName, ReadStatus, "Error found in NewPPHdrNL" )
    END IF

    Fields(Source(1:NumFlds)) % Hdr = NewPPHdr(1:NumFlds)

! #if !defined(NODBG)
!     WRITE(6,*) Fields(Source(1:NumFlds)) % Hdr
! #endif

! DEPENDS ON: writeflds
    CALL WriteFlds( NumFlds, MaxFlds, MaxFldsOut, PackType, Source, &
                    PackAcc, Fields, StdLevHdr_out, ErrorStatus )

 !---------------------------------
 ! 2.3 General Arithmetic Actions
 !---------------------------------
  ELSE IF ( Action == "SUM" )       THEN

   ! 2.3.1 Sum two or more fields
! DEPENDS ON: sum
    CALL Sum( Fields(Source(1)), Fields(Source(2)),     &
              Fields(Store(1)), ErrorStatus )
    DO i = 3,NumFlds
! DEPENDS ON: sum
      CALL Sum( Fields(Source(i)), Fields(Store(1)),    &
                Fields(Store(1)), ErrorStatus )
    END DO

  ELSE IF ( Action == "DIF" )       THEN

   ! 2.3.2 Take difference of 2 fields
    IF ( NumFlds /= 2 ) THEN
      WRITE(6,*) "Must have 2 fields for differencing"
      ErrorStatus = StatusWarning
    ELSE
! DEPENDS ON: dif
      CALL Dif( Fields(Source(1)), Fields(Source(2)),   &
                Fields(Store(1)), ErrorStatus )
    END IF

  ELSE IF ( Action == "ADD" )       THEN

   ! 2.3.3 Add constant to fields
    DO i = 1,NumFlds
! DEPENDS ON: add
      CALL Add( Factor(1), Fields(Source(i)),   &
                Fields(Store(i)), ErrorStatus )
    END DO

  ELSE IF ( Action == "SCALE" )     THEN

   ! 2.3.4 Multipy fields by a constant
    DO i = 1,NumFlds
! DEPENDS ON: scale
      CALL Scale( Factor(1), Fields(Source(i)),   &
                  Fields(Store (i)), ErrorStatus )
    END DO

  ELSE IF ( Action == "WINDSPEED" ) THEN

   ! 2.3.5 Calculate magnitude of a vector (2 flds) (e.g. - wind speed)
    IF ( NumFlds /= 2 ) THEN
      WRITE(6,*) "Must have a U and a V field for WindSpeed"
      ErrorStatus = StatusWarning
    ELSE
! DEPENDS ON: vecmag
      CALL VecMag( Fields(Source(1)), Fields(Source(2)),   &
                   Fields(Store(1)), ErrorStatus )
    END IF

 !---------------------------------------------------
 ! 2.4 Diagnostic Actions Requiring Standard Levels
 !---------------------------------------------------
  ELSE IF ( Action == "DIVERG" )    THEN

   ! 2.4.1 Calculate the Divergence of 2 fields
    IF ( NumFlds /= 2 ) THEN
      WRITE(6,*) "Must have a U and a V field for Diverg"
      ErrorStatus = StatusWarning
    ELSE
! DEPENDS ON: diverg
      CALL Diverg( Fields(Source(1)), Fields(Source(2)),   &
                   Fields(Store(1)), ErrorStatus )
    END IF

  ELSE IF ( Action == "VORTIC" )    THEN

   ! 2.4.2 Calculate the Vorticity of 2 fields
    IF ( NumFlds /= 2 ) THEN
      WRITE(6,*) "Must have a U and a V field for Vortic"
      ErrorStatus = StatusWarning
    ELSE
! DEPENDS ON: vortic
      CALL Vortic( Fields(Source(1)), Fields(Source(2)),   &
                   Fields(Store(1)), ErrorStatus )
    END IF

  ELSE IF ( Action == "THERMADV" )  THEN

   ! 2.4.3 Thermal Advection
    IF ( NumFlds /= 3 ) THEN
      WRITE(6,*) "Must have U, V and T fields for Thermal Advection"
      ErrorStatus = StatusWarning
    ELSE
! DEPENDS ON: thermadv
      CALL ThermAdv( Fields(Source(1)), Fields(Source(2)),   &
                     Fields(Source(3)),                      &
                     Fields(Store(1)), ErrorStatus )
    END IF

  ELSE IF ( Action == "ICAO_HT" )   THEN

   ! 2.4.4 Convert pressure to ICAO height
    DO i = 1,NumFlds
! DEPENDS ON: icaoheight
      CALL ICAOHeight( Fields(Source(i)),    &
                       Fields(Store(i)), ErrorStatus )
    END DO

  ELSE IF ( Action == "SNOWPROB" )  THEN

   ! 2.4.5 Snow Probability
    IF ( NumFlds /= 2 ) THEN
      WRITE(6,*) "Must have two height fields for Snow Probability"
      ErrorStatus = StatusWarning
    ELSE
! DEPENDS ON: snowprob
      CALL SnowProb( Fields(Source(1)), Fields(Source(2)),  &
                     Fields(Store(1)), ErrorStatus )
    END IF

  ELSE IF ( Action == "PRSYM" )     THEN

   ! 2.4.6 Total Precipitation Code
! DEPENDS ON: prsym
    CALL PrSym( Fields(Source(1)), Fields(Source(2)),     &
                Fields(Source(3)), Fields(Source(4)),     &
                Fields(Source(5)),                        &
                Fields(Store(1)), ErrorStatus )

  ELSE IF ( Action == "WXCODE" )    THEN

   ! 2.4.7 Present Weather Code
! DEPENDS ON: wxcode
    CALL WXCode( Fields(Source( 1)), Fields(Source( 2)),  &
                 Fields(Source( 3)), Fields(Source( 4)),  &
                 Fields(Source( 5)), Fields(Source( 6)),  &
                 Fields(Source( 7)), Fields(Source( 8)),  &
                 Fields(Source( 9)), Fields(Source(10)),  &
                 Fields(Source(11)), Fields(Source(12)),  &
                 Fields(Store(1)), ErrorStatus )
    IF ( NumFlds == 2 ) THEN
! DEPENDS ON: prsym
      CALL PrSym( Fields(Source(1)), Fields(Source(2)),   &
                  Fields(Source(3)), Fields(Source(4)),   &
                  Fields(Source(5)),                      &
                  Fields(Store(2)), ErrorStatus )
    END IF

 !------------------------------------------------
 ! 2.5 Diagnostic Actions Requiring Model Levels
 !------------------------------------------------
  ELSE IF ( Action == "TROPHGHT" )  THEN

   ! 2.5.1 Tropopause Temperature, Pressure and Height
    IF ( P_Levels /= UM_NumLevs ) THEN
      WRITE(6,*) "Can only run TROPHGHT on a 38 or 50 level model"
      ErrorStatus = StatusWarning
    END IF
    NumLevs = TP_NumLevs_Gl
    ZeroLev = TP_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)
    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)

    ! Read theta level pressure
    STCode  (1:NumLevs) = ST_Ptheta
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,           &
                  LBProc, MinsPastHr, PSource, PPHdrMod, ModLevHdr_in,  &
                  Fields, ErrorStatus )
    ! Read theta level temperature
    STCode  (1:NumLevs) = ST_Ttheta
    TSource (1:NumLevs) = TPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,           &
                  LBProc, MinsPastHr, TSource, PPHdrMod, ModLevHdr_in,  &
                  Fields, ErrorStatus )
    ! Calculate theta level heights
    ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)
! DEPENDS ON: calczflds
    CALL CalcZFlds( NumLevs, OrogField,                         &
                    Fields(PSource(1):PSource(NumLevs)),        &
                    Fields(ZSource(1):ZSource(NumLevs)),        &
                    ErrorStatus )

    ! Do calculations
! DEPENDS ON: tropheight
    CALL TropHeight( NumLevs,                                   &
                     Fields(PSource(1):PSource(NumLevs)),       &
                     Fields(TSource(1):TSource(NumLevs)),       &
                     Fields(ZSource(1):ZSource(NumLevs)),       &
                     Fields(Store(1)), Fields(Store(2)),        &
                     Fields(Store(3)), ErrorStatus )

  ELSE IF ( Action == "ISOTHERM" )  THEN

   ! 2.5.2 Isotherm for temperature given in Factor
    IF ( P_Levels /= UM_NumLevs ) THEN
      WRITE(6,*) "Can only run ISOTHERM on a 38 or 50 level model"
      ErrorStatus = StatusWarning
    END IF
    NumLevs = IT_NumLevs_Gl
    ZeroLev = IT_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)
    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)

    ! Read theta level pressure
    STCode  (1:NumLevs) = ST_Ptheta
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, PSource, PPHdrMod, ModLevHdr_in, &
                  Fields, ErrorStatus )
    ! Read theta level temperature
    STCode  (1:NumLevs) = ST_Ttheta
    TSource (1:NumLevs) = TPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, TSource, PPHdrMod, ModLevHdr_in, &
                  Fields, ErrorStatus )
    ! Calculate theta level heights
    ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)
! DEPENDS ON: calczflds
    CALL CalcZFlds( NumLevs, OrogField,                         &
                    Fields(PSource(1):PSource(NumLevs)),        &
                    Fields(ZSource(1):ZSource(NumLevs)),        &
                    ErrorStatus )
    ! Read pstar
    STCode  (1) = ST_Pstar
    MO8Level(1) = LV_Surface
    Source  (1) = StdPosn + 1
! DEPENDS ON: getflds
    CALL GetFlds( 1, MaxFlds, STCode, MO8Level, FCTime,               &
                  LBProc, MinsPastHr, Source, PPHdrMod, ModLevHdr_in, &
                  Fields, ErrorStatus )

    ! Do calculations
! DEPENDS ON: isotherm
    CALL IsoTherm( NumLevs, Factor(1), OrogField, Fields(StdPosn+1), &
                   Fields(PSource(1):PSource(NumLevs)),              &
                   Fields(TSource(1):TSource(NumLevs)),              &
                   Fields(ZSource(1):ZSource(NumLevs)),              &
                   Fields(Store(1)), Fields(Store(2)),               &
                   ErrorStatus )

  ELSE IF ( Action == "CONTRAIL" )  THEN

   ! 2.5.3 Upper and Lower Contrail Limit
    IF ( P_Levels /= UM_NumLevs ) THEN
      WRITE(6,*) "Can only run CONTRAIL on a 38 or 50 level model"
      ErrorStatus = StatusWarning
    END IF
    NumLevs = CT_NumLevs_Gl
    ZeroLev = CT_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)
    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)

    ! Read theta level pressure
    STCode  (1:NumLevs) = ST_Ptheta
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, PSource, PPHdrMod, ModLevHdr_in, &
                  Fields, ErrorStatus )
    ! Read theta level temperature
    STCode  (1:NumLevs) = ST_Ttheta
    TSource (1:NumLevs) = TPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, TSource, PPHdrMod, ModLevHdr_in, &
                  Fields, ErrorStatus )
    ! Calculate theta level heights
    ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)
! DEPENDS ON: calczflds
    CALL CalcZFlds( NumLevs, OrogField,                         &
                    Fields(PSource(1):PSource(NumLevs)),        &
                    Fields(ZSource(1):ZSource(NumLevs)),        &
                    ErrorStatus )
     ! Read pstar
    STCode  (1) = ST_Pstar
    MO8Level(1) = LV_Surface
    Source  (1) = StdPosn + 1
! DEPENDS ON: getflds
    CALL GetFlds( 1, MaxFlds, STCode, MO8Level, FCTime,               &
                  LBProc, MinsPastHr, Source, PPHdrMod, ModLevHdr_in, &
                  Fields, ErrorStatus )

    ! Do calculations
! DEPENDS ON: contrail
    CALL Contrail( NumLevs, Fields(StdPosn + 1),                &
                   Fields(PSource(1):PSource(NumLevs)),         &
                   Fields(TSource(1):TSource(NumLevs)),         &
                   Fields(ZSource(1):ZSource(NumLevs)),         &
                   Fields(Store(1)), Fields(Store(2)),          &
                   ErrorStatus )

  ELSE IF ( Action == "MAXWIND" )   THEN

   ! 2.5.4 Maximum Wind U&V Values and Pressure level on B grid
    IF ( P_Levels == UM_NumLevs ) THEN
      NumLevs = MX_NumLevs_Gl
    ELSE
      ! Find lowest level with Zsea >= MX_Zsea_Upr
      NumLevs = 1
      n = ModLevHdr_in % Len1LevDepC * 6
      DO i = 1, ModLevHdr_in % Len1LevDepC
        Zsea = ModLevHdr_in % LevDepC(i + n)
        IF ( (Zsea < MX_Zsea_Upr) .AND. (Zsea > 0) ) THEN
          NumLevs = NumLevs + 1
        END IF
      END DO
      WRITE(6,*) "MaxWind levels reset to ", NumLevs
    END IF
    ZeroLev = MX_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)

    ! Read U @ 300hPa to use as a template for the B-grid
    STCode   = ST_Ustd
    MO8Level = 300
    Source   = StdPosn + 1
! DEPENDS ON: getflds
    CALL GetFlds( 1, MaxFlds, STCode, MO8Level, FCTime,               &
                  LBProc, MinsPastHr, Source, PPHdrMod, StdLevHdr_in, &
                  Fields, ErrorStatus )

    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)

    ! Read rho level u-component of wind & interpolate to B-grid

    STCode  (1:NumLevs) = ST_Urho
    USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,             &
                    MinsPastHr, USource, PPHdrMod, Fields(StdPosn+1) % Hdr, &
                    ModLevHdr_in, Fields, ErrorStatus )

    ! Read rho level v-component of wind & interpolate to B-grid

    STCode  (1:NumLevs) = ST_Vrho
    VSource (1:NumLevs) = VPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,             &
                    MinsPastHr, VSource, PPHdrMod, Fields(StdPosn+1) % Hdr, &
                    ModLevHdr_in, Fields, ErrorStatus )

    ! Read rho level pressure & interpolate to B-grid

    STCode  (1:NumLevs) = ST_Prho
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,             &
                    MinsPastHr, PSource, PPHdrMod, Fields(StdPosn+1) % Hdr, &
                    ModLevHdr_in, Fields, ErrorStatus )

    ! Do calculations

    ! Store1/2/3 : Max Wind U/V/P

! DEPENDS ON: maxwind
    CALL MaxWind( NumLevs,                                      &
                  Fields(USource(1):USource(NumLevs)),          &
                  Fields(VSource(1):VSource(NumLevs)),          &
                  Fields(PSource(1):PSource(NumLevs)),          &
                  Fields(Store(1)),                             &
                  Fields(Store(2)),                             &
                  Fields(Store(3)),                             &
                  ErrorStatus )

  ELSE IF ( Action == "TOPBASE" )   THEN

    ! 2.5.5 Max wind Top and Base
    ! Action MAXWIND must be called first

    ! Source1/2/3 : Max Wind U/V/P
    ! Store1/2    : Max Wind Base/Top

! DEPENDS ON: topbase
    CALL TopBase( NumLevs,                                      &
                  Fields(USource(1):USource(NumLevs)),          &
                  Fields(VSource(1):VSource(NumLevs)),          &
                  Fields(PSource(1):PSource(NumLevs)),          &
                  Fields(Source(1)),                            &
                  Fields(Source(2)),                            &
                  Fields(Source(3)),                            &
                  Fields(Store(1)),                             &
                  Fields(Store(2)),                             &
                  ErrorStatus )

  ELSE IF ( Action == "CATURB" )    THEN

   ! 2.5.6 Clear Air Turbulence Predictor
    IF ( (P_Levels /= UM_NumLevs ) .OR. &
         (GridType /= Global) ) THEN
      WRITE(6,*) "Can only run CATURB for 38 or 50 level Global field"
      ErrorStatus = StatusWarning
    END IF
    NumLevs = CA_NumLevs_Gl
    ZeroLev = CA_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)

    ! Read standard level fields
    STcode  (1:CA_NumStd) = (/ ST_Ustd, ST_Ustd, ST_Ustd,       &
                               ST_Vstd, ST_Vstd, ST_Vstd /)
    MO8Level(1:CA_NumStd) = (/     300,     250,     200,       &
                                   300,     250,     200 /)
    Source  (1:CA_NumStd) = StdPosn + (/ (i, i=1,CA_NumStd) /)
! DEPENDS ON: getflds
    CALL GetFlds( CA_NumStd, MaxFlds, STCode, MO8Level, FCTime,       &
                  LBProc, MinsPastHr, Source, PPHdrMod, StdLevHdr_in, &
                  Fields, ErrorStatus )

    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)
    ! Read rho level u-component of wind & interpolate to B-grid
    STCode  (1:NumLevs) = ST_Urho
    USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,             &
                    MinsPastHr, USource, PPHdrMod, Fields(StdPosn+1) % Hdr, &
                    ModLevHdr_in, Fields, ErrorStatus )

    ! Read rho level v-component of wind & interpolate to B-grid
    STCode  (1:NumLevs) = ST_Vrho
    VSource (1:NumLevs) = VPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,             &
                    MinsPastHr, VSource, PPHdrMod, Fields(StdPosn+1) % Hdr, &
                    ModLevHdr_in, Fields, ErrorStatus )

    ! Read rho level pressure & interpolate to B-grid
    STCode  (1:NumLevs) = ST_Prho
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,             &
                    MinsPastHr, PSource, PPHdrMod, Fields(StdPosn+1) % Hdr, &
                    ModLevHdr_in, Fields, ErrorStatus )

    ! Calculate rho level heights using B(UV)-grid orography
    ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)
    IF ( ASSOCIATED( OrogField % RData ) ) THEN
      ! Interpolate Orography to B-grid
! DEPENDS ON: ctobgrid
      CALL CtoBgrid( OrogField, Fields(StdPosn+1) % Hdr,        &
                     Fields(StdPosn+10), ErrorStatus )
    ELSE
! DEPENDS ON: ereport
      CALL EReport( ProgName, StatusWarning,                    &
                   "Orography not available - cannot calculate height fields" )
      ErrorStatus = StatusWarning
    END IF
! DEPENDS ON: calczflds
    CALL CalcZFlds( NumLevs, Fields(StdPosn+10),                &
                    Fields(PSource(1):PSource(NumLevs)),        &
                    Fields(ZSource(1):ZSource(NumLevs)),        &
                    ErrorStatus )

    ! Do calculations
    DO i = 1,3
! DEPENDS ON: caturb
      CALL CATurb( NumLevs,                                     &
                   Fields(StdPosn+i), Fields(StdPosn+3+i),      &
                   Fields(USource(1):USource(NumLevs)),         &
                   Fields(VSource(1):VSource(NumLevs)),         &
                   Fields(PSource(1):PSource(NumLevs)),         &
                   Fields(ZSource(1):ZSource(NumLevs)),         &
                   Fields(Store(i)), ErrorStatus )
    END DO
! DEPENDS ON: maxcaturb
    CALL MaxCATurb( Fields(Store(1)), Fields(Store(2)),         &
                    Fields(Store(3)),                           &
                    Fields(Store(4)), Fields(Store(5)),         &
                    ErrorStatus )

  ELSE IF ( Action == "MWTURB" )    THEN

   ! 2.5.7 Mountain Wave Turbulence Predictor
    IF ( (P_Levels /= UM_NumLevs ) .OR. &
         (GridType /= Global) ) THEN
      WRITE(6,*) "Can only run MWTURB for 38 or 50 level Global field"
      ErrorStatus = StatusWarning
    END IF
    NumLevs = MW_NumLevs_Gl
    ZeroLev = MW_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)

    ! Read standard level fields
    Source  (1:MW_NumStd) = StdPosn + (/3,6/)
    STCode  (1:MW_NumStd) = (/ ST_Ustd, ST_Vstd /)
    MO8Level(1:MW_NumStd) = (/     200,     200 /)
! DEPENDS ON: getflds
    CALL GetFlds( MW_NumStd, MaxFlds, STCode, MO8Level, FCTime,       &
                  LBProc, MinsPastHr, Source, PPHdrMod, StdLevHdr_in, &
                  Fields, ErrorStatus )

    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)

    ! Read u-component of gravity wave stress & interpolate to B-grid
    STCode  (1:NumLevs) = ST_GWSU
    USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,             &
                    MinsPastHr, USource, PPHdrMod, Fields(StdPosn+3) % Hdr, &
                    ModLevHdr_in, Fields, ErrorStatus )

    ! Read v-component of gravity wave stress & interpolate to B-grid
    STCode  (1:NumLevs) = ST_GWSV
    VSource (1:NumLevs) = VPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,             &
                    MinsPastHr, VSource, PPHdrMod, Fields(StdPosn+3) % Hdr, &
                    ModLevHdr_in, Fields, ErrorStatus )

    ! Do calculations
! DEPENDS ON: mtnstress
    CALL MtnStress( NumLevs,                                    &
                    Fields(StdPosn+3), Fields(StdPosn+6),       &
                    Fields(USource(1):USource(NumLevs)),        &
                    Fields(Vsource(1):VSource(NumLevs)),        &
                    Fields(Store(1)), ErrorStatus )

  ELSE IF ( Action == "WAFC_TURB" )    THEN

   ! 2.5.8 WAFC CAT turb Predictor (includes shear & mountain wave turb)
    IF ( (P_Levels /= UM_NumLevs) .OR. &
         (GridType /= Global) ) THEN
      WRITE(6,*) "Can only run WAFC_TURB on a 38 or 50 level Global field"
      ErrorStatus = StatusWarning
    END IF

    NumLevs = WT_NumLevs_Gl
    ZeroLev = WT_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)


    !---Calculate shear turbulence predictor (TI1)--------------
    ! This will be calculated for 400, 300, 250, 200, 150 hPa & put into
    ! positions Store(i) to Store(5).

    ! Read standard level fields
    STcode  (1:WT_NumStd) = (/ ST_Ustd, ST_Ustd, ST_Ustd,   &
                               ST_Ustd, ST_Ustd,            &
                               ST_Vstd, ST_Vstd, ST_Vstd,   &
                               ST_Vstd, ST_Vstd /)
    MO8Level(1:WT_NumStd) = (/ 400, 300, 250, 200, 150,    &
                               400, 300, 250, 200, 150 /)
    Source  (1:WT_NumStd) = StdPosn + (/ (i, i=1,WT_NumStd) /)
! DEPENDS ON: getflds
    CALL GetFlds( WT_NumStd, MaxFlds, STCode, MO8Level, FCTime,       &
                  LBProc, MinsPastHr, Source, PPHdrMod, StdLevHdr_in, &
                  Fields, ErrorStatus )
    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)


    ! Read rho level u-component of wind & interpolate to B-grid
    STCode  (1:NumLevs) = ST_Urho
    USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,             &
                    MinsPastHr, USource, PPHdrMod, Fields(StdPosn+1) % Hdr, &
                    ModLevHdr_in, Fields, ErrorStatus )

    ! Read rho level v-component of wind & interpolate to B-grid
    STCode  (1:NumLevs) = ST_Vrho
    VSource (1:NumLevs) = VPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,             &
                    MinsPastHr, VSource, PPHdrMod, Fields(StdPosn+1) % Hdr, &
                    ModLevHdr_in, Fields, ErrorStatus )

    ! Read rho level pressure & interpolate to B-grid
    STCode  (1:NumLevs) = ST_Prho
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
    CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,             &
                    MinsPastHr, PSource, PPHdrMod, Fields(StdPosn+1) % Hdr, &
                    ModLevHdr_in, Fields, ErrorStatus )

    ! Calculate rho level heights using B(UV)-grid orography
    ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)
    IF ( ASSOCIATED( OrogField % RData ) ) THEN
      ! Interpolate Orography to B-grid
! DEPENDS ON: ctobgrid
      CALL CtoBgrid( OrogField, Fields(StdPosn+1) % Hdr,        &
                     Fields(StdPosn+15), ErrorStatus )
    ELSE
! DEPENDS ON: ereport
      CALL EReport( ProgName, StatusWarning,                    &
                   "Orography not available - cannot calculate height fields" )
      ErrorStatus = StatusWarning
    END IF
! DEPENDS ON: calczflds
    CALL CalcZFlds( NumLevs, Fields(StdPosn+15),                &
                    Fields(PSource(1):PSource(NumLevs)),        &
                    Fields(ZSource(1):ZSource(NumLevs)),        &
                    ErrorStatus )

    ! Step through each pressure level & calculate TI1 Index
    DO i = 1, WT_NumStd / 2
! DEPENDS ON:wafc_caturb
      CALL WAFC_CATurb( NumLevs,                                &
                   Fields(StdPosn+i), Fields(StdPosn+5+i),      &
                   Fields(USource(1):USource(NumLevs)),         &
                   Fields(VSource(1):VSource(NumLevs)),         &
                   Fields(PSource(1):PSource(NumLevs)),         &
                   Fields(ZSource(1):ZSource(NumLevs)),         &
                   Fields(Store(i)),ErrorStatus )
    END DO


    !-------------------------------------------------------------
    ! Fields only available to calculate mountain wave part up to T+24
    ! After that just duplicate the shear turbulence fields

    IF (FCTime(1) <= 24) THEN

    !---Calculate mountain wave turbulence predictor --------------
    ! This will be calculated for 300, 250, 200, 150 hPa & put into
    ! positions Store(6) to Store(9).
    ! The standard level fields u and v are re-used from the shear CAT part rather
    ! than re-reading them in. The positions for the model level fields in the
    ! shear CAT part are re-used to read in the required model fields for the MW
    ! predictor calculation.

      NumLevs = MW_NumLevs_Gl
      ZeroLev = MW_ZeroLev_Gl
      FCTime  (1:NumLevs) = FCTime(1)

      !Read model level fields (u and v component of gravity wave stress etc)
      MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)

      ! Read u-component of gravity wave stress & interpolate to B-grid
      STCode  (1:NumLevs) = ST_GWSU
      USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
      CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,             &
                      MinsPastHr, USource, PPHdrMod, Fields(StdPosn+3) % Hdr, &
                      ModLevHdr_in, Fields, ErrorStatus )

      ! Read v-component of gravity wave stress & interpolate to B-grid
      STCode  (1:NumLevs) = ST_GWSV
      VSource (1:NumLevs) = VPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
      CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,             &
                      MinsPastHr, VSource, PPHdrMod, Fields(StdPosn+3) % Hdr, &
                      ModLevHdr_in, Fields, ErrorStatus )

      ! Read pressure & interpolate to B-grid
      STCode  (1:NumLevs) = ST_Ptheta
      PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getb_flds
      CALL GetB_Flds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,             &
                      MinsPastHr, PSource, PPHdrMod, Fields(StdPosn+3) % Hdr, &
                      ModLevHdr_in, Fields, ErrorStatus )

      ! Loop through and do calculations for each reference pressure
      ! level required defined in array MW_PRef(i)
      ! Store results in the  Fields elements after the CAT fields
      ! (i.e. in the locations after Store((WT_NumStd/2)*2)   )

      DO i = 1, MW_n_PRef
! DEPENDS ON: mtnstress_pref
        CALL MtnStress_PRef( NumLevs,                                    &
                             MW_PRef(i),                                 &
                             Fields(StdPosn+1+i), Fields(StdPosn+6+i),   &
                             Fields(USource(1):USource(NumLevs)),        &
                             Fields(Vsource(1):VSource(NumLevs)),        &
                             Fields(Psource(1):PSource(NumLevs)),        &
                             Fields( Store(i+(WT_NumStd/2)) ), ErrorStatus)
      END DO


      ! Combine CAT and MW fields by taking the maximum of the fields at each gridpoint
      ! for each level where both CAT and MW turbulence are produced (i.e.
      ! 300, 250, 200, 150 hPa (400 hPa no produced for MW turbulence)).

      DO i = 1, MW_n_PRef
! DEPENDS ON: maximum
        CALL Maximum( Fields( Store(i+1) ),                         &  !CAT turbulence field
                      Fields( Store( (WT_NumStd/2)+i ) ),           &  !MW turbulence field
                      Fields( Store( (WT_NumStd/2)+MW_n_PRef+i ) ), &  !combined field
                      ErrorStatus )
      END DO


      !Duplicate the combined CAT/MW fields so we have fields for each level
      !labelled with mean and maximum CAT stash codes

      !for 400mb (contains CAT only)

! DEPENDS ON: dup_cat
      CALL Dup_CAT( Fields(Store(1)), &                              !input CAT/MW field
                    Fields(Store((WT_NumStd/2)+(2*MW_n_PRef)+1)), &  !output duplicate field
                    ErrorStatus)

      !Loop through 300mb, 250mb, 200mb, 150mb combined CAT/MW fields
      DO i=1, MW_n_PRef

! DEPENDS ON: dup_cat
        CALL Dup_CAT( Fields(Store((WT_NumStd/2)+MW_n_PRef+i)),       &  !input CAT/MW field
                      Fields(Store((WT_NumStd/2)+(2*MW_n_PRef)+1+i)),  & !output duplicate field
                      ErrorStatus )

      END DO

    ELSE

      !Mountain wave fields not available - just duplicate shear CAT fields
      DO i=1, WT_NumStd/2

! DEPENDS ON: dup_cat
        CALL Dup_CAT( Fields( Store(i) ),                    &  !input CAT field
                      Fields( Store((WT_NumStd/2)+i) ),      &  !output duplicate field
                      ErrorStatus )

      END DO

    ENDIF

  ELSE IF ( Action == "OLD_ICING" )  THEN

   ! 2.5.9 WAFC Icing Predictor (Large scale cloud layer Icing algorithm LIBLKCLD)
   !       BuLK CLouD fraction for cloud in a given temperature range

    IF ( P_Levels /= UM_NumLevs ) THEN
      WRITE(6,*) "Can only run LIBLKCLD on a 38 or 50 level model"
      ErrorStatus = StatusWarning
    END IF

    NumLevs = IC_NumLevs_Gl
    ZeroLev = IC_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)
    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)

    USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
    VSource (1:NumLevs) = VPosn + MO8Level(1:NumLevs)
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
    TSource (1:NumLevs) = TPosn + MO8Level(1:NumLevs)
    ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)

    ! Read theta level bulk cloud fraction
    STCode  (1:NumLevs) = ST_BlkCld
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, USource, PPHdrMod, ModLevHdr_in, &
                  Fields, ErrorStatus )

    ! Read theta level pressure in Pa
    STCode  (1:NumLevs) = ST_Ptheta
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, PSource, PPHdrMod, ModLevHdr_in, &
                  Fields, ErrorStatus )

    ! Read theta level temperature in K
    STCode  (1:NumLevs) = ST_Ttheta
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, TSource, PPHdrMod, ModLevHdr_in, &
                  Fields, ErrorStatus )

    ! Find theta level gridboxes in the given temp range
! DEPENDS ON: temperature_mask
    CALL Temperature_Mask( Factor(1:2), NumLevs,                        &
                   Fields(TSource(1):TSource(NumLevs)),         &
                   Fields(VSource(1):VSource(NumLevs)),         &
                   ErrorStatus )

    ! Calculate theta level cloud fraction for the temp range
! DEPENDS ON: multiply
    CALL Multiply( NumLevs,                                     &
                   Fields(USource(1):USource(NumLevs)),         &
                   Fields(VSource(1):VSource(NumLevs)),         &
                   Fields(ZSource(1):ZSource(NumLevs)),         &
                   ErrorStatus )

    ! Calculate pressure layer cloud fraction for the temp range
! DEPENDS ON: iconplyr
    CALL ICOnPLyr( NumLevs, NumLyrs,                            &
                   LyrMO8L(1:NumLyrs),                          &
                   LyrLwrB(1:NumLyrs),                          &
                   LyrUprB(1:NumLyrs),                          &
                   Fields(ZSource(1):ZSource(NumLevs)),         &
                   Fields(PSource(1):PSource(NumLevs)),         &
                   Fields(Store(1):Store(2*NumLyrs)),           &
                   ErrorStatus )



  ELSE IF ( Action == "ICING" )  THEN

!-----------------------------------------------------------------------------!
!                 Relative Humidity based icing diagnostic                    !
!-----------------------------------------------------------------------------!
! Method:                                                                     !
! 1. Read in pressure, temperature and specific humidity on model levels.     !
!    Calculate saturated mixing ratio from P & T and then combine this with   !
!    specific humidity to get RH (detailed in Rel_humdy subroutine). Then     !
!    determine RH on pressure levels.                                         !
! 2. Read in cloud fraction on model levels. Create mask where cloud is       !
!    present. Determine cloud fraction on pressure levels.                    !
! 3. Create -20 to 0C temperature mask. Convert on to pressure levels.        !
! 4. Multiply the three pressure level arrays to give RH where cloud is       !
!    present and temperature is between 0 and -20C                            !
!-----------------------------------------------------------------------------!

   ! 2.5.9 WAFC Icing Predictor
   !       Relative humidity for cloud in a given temperature range

    IF ( P_Levels /= UM_NumLevs ) THEN
      WRITE(6,*) "Can only run Icing algorithm on a 38 or 50 level model"
      ErrorStatus = StatusWarning
    END IF

    NumLevs = LI_NumLevs_Gl
    ZeroLev = LI_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)
    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)


    USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
    VSource (1:NumLevs) = VPosn + MO8Level(1:NumLevs)
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
    TSource (1:NumLevs) = TPosn + MO8Level(1:NumLevs)
    ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)

!----------------------------------------------------------------
!---- Step 1: Read in pressure,temperature, specific humidity.---
!---- Calculate RH and interpolate on to pressure levels      ---
!----------------------------------------------------------------

    ! Read theta level pressure
    STCode  (1:NumLevs) = ST_Ptheta
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, PSource, PPHdrMod, ModLevHdr_in, &
                  Fields, ErrorStatus )


    ! Read theta level temperature in K
    STCode  (1:NumLevs) = ST_Ttheta
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, TSource, PPHdrMod, ModLevHdr_in, &
                  Fields, ErrorStatus )


    ! Read theta level specific humidity
    STCode  (1:NumLevs) = ST_Htheta
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, USource, PPHdrMod, ModLevHdr_in, &
                  Fields, ErrorStatus )

    ! Calculate saturated mixing ratio on model levels
! DEPENDS ON: sat_mratio
    CALL Sat_MRatio( NumLevs,                                          &
                     Fields(PSource(1):PSource(NumLevs)),              &
                     Fields(TSource(1):TSource(NumLevs)),              &
                     Fields(VSource(1):VSource(NumLevs)),              &
                     ErrorStatus )

    ! Calculate relative humidity on model levels
! DEPENDS ON: rel_humdy
    CALL Rel_humdy( NumLevs,                                           &
                    Fields(USource(1):USource(NumLevs)),               &
                    Fields(VSource(1):VSource(NumLevs)),               &
                    Fields(ZSource(1):ZSource(NumLevs)),               &
                    ErrorStatus )

    ! Nullify bit of the array that's about to be reused
    DO i=1,NumLevs
      DEALLOCATE ( Fields(USource(i)) % RData )
      DEALLOCATE ( Fields(VSource(i)) % RData )
      NULLIFY( Fields(USource(i)) % RData )
      NULLIFY( Fields(VSource(i)) % RData )
    END DO

    ! Calculate relative humidity on pressure levels
! DEPENDS ON: lionplyr
    CALL LIOnPLyr( NumLevs, NumLyrs,                                   &
                   LyrMO8L, LyrLwrB, LyrUprB,                          &
                   Fields(ZSource(1):ZSource(NumLevs)),                &
                   Fields(PSource(1):PSource(NumLevs)),                &
                   Fields(USource(1):USource(NumLyrs)),                &
                   ErrorStatus )

! Nullify bit of the array that's about to be reused
    DO i=1, NumLevs
      DEALLOCATE ( Fields(ZSource(i)) % RData )
      NULLIFY( Fields(ZSource(i)) % RData )
    END DO


!----------------------------------------------------------------
!---- Step 2: Determine temperature on pressure levels        ---
!---- by calulating mean temperature across atmospheric layer ---
!----------------------------------------------------------------

    ! Determine temperature on pressure levels
! DEPENDS ON: lionplyr
    CALL LIOnPLyr( NumLevs, NumLyrs,                                   &
                   LyrMO8L, LyrLwrB, LyrUprB,                          &
                   Fields(TSource(1):TSource(NumLevs)),                &
                   Fields(PSource(1):PSource(NumLevs)),                &
                   Fields(ZSource(1):ZSource(NumLyrs)),                &
                   ErrorStatus )

    ! Nullify bit of the array that's about to be reused
    DO i=1, NumLevs
      DEALLOCATE( Fields(TSource(i)) % RData )
      NULLIFY( Fields(TSource(i)) % RData )
    END DO

    ! Calculate which gridboxes have temperatures between -20 and 0C
! DEPENDS ON: temperature_mask
    CALL Temperature_Mask( Factor(1:2), NumLyrs,                       &
                           Fields(ZSource(1):ZSource(NumLyrs)),        &
                           Fields(TSource(1):TSource(NumLyrs)),        &
                           ErrorStatus )

!----------------------------------------------------------------
!---- Step 3: Read in cloud fraction on model levels         ----
!---- Determine cloud fraction on pressure levels            ----
!----------------------------------------------------------------

    ! Read theta level bulk cloud fraction
    STCode  (1:NumLevs) = ST_BlkCld
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,         &
                  LBProc, MinsPastHr, ZSource, PPHdrMod, ModLevHdr_in,&
                  Fields, ErrorStatus )


    ! Determine cloud fraction on pressure levels
! DEPENDS ON: lionplyr
    CALL LIOnPLyr( NumLevs, NumLyrs,                                  &
                   LyrMO8L, LyrLwrB, LyrUprB,                         &
                   Fields(ZSource(1):ZSource(NumLevs)),               &
                   Fields(PSource(1):PSource(NumLevs)),               &
                   Fields(VSource(1):VSource(NumLyrs)),               &
                   ErrorStatus )

    ! Nullify bit of the array that's about to be reused
    DO i=1, NumLevs
      DEALLOCATE( Fields(ZSource(i)) % RData )
      NULLIFY( Fields(ZSource(i)) % RData )
    END DO


    ! Calculate which gridboxes have cloud present
! DEPENDS ON: cloud_mask
    CALL Cloud_Mask( (/ 0.0, 1.1 /), NumLyrs,                              &
                     Fields(VSource(1):VSource(NumLyrs)),              &
                     Fields(ZSource(1):ZSource(NumLyrs)),              &
                     ErrorStatus )

! Transfer masked cloud mask back into VSource section, has to be done this way
! due to rdata being a pointer.  If we perform field1 = field2 then the field1
! rdata will point to field2 and will not copy the data.
    DO i=1, NumLyrs
      Fields(VSource(i)) % Hdr       = Fields(ZSource(i)) % Hdr
      Fields(VSource(i)) % Rdata     = Fields(ZSource(i)) % Rdata
      Fields(VSource(i)) % LookupPos = Fields(ZSource(i)) % LookupPos
      Fields(VSource(i)) % ArrayPos  = Fields(ZSource(i)) % ArrayPos
      
      DEALLOCATE ( Fields(ZSource(i)) % RData )
      NULLIFY( Fields(ZSource(i)) % RData )
    END DO

!-----------------------------------------------------------------
!---- Step 4: Apply temperature and cloud masks to             ---
!----         relative humidity                                ---
!-----------------------------------------------------------------

    ! Apply temperature mask to relative humidity fields
! DEPENDS ON: multiply
    CALL Multiply( NumLyrs,                                           &
                   Fields(USource(1):USource(NumLyrs)),               &
                   Fields(TSource(1):TSource(NumLyrs)),               &
                   Fields(ZSource(1):ZSource(NumLyrs)),               &
                   ErrorStatus )


    ! Apply cloud mask to relative humidity fields masked by temperature
! DEPENDS ON: multiply
    CALL Multiply( NumLyrs,                                           &
                   Fields(ZSource(1):ZSource(NumLyrs)),               &
                   Fields(VSource(1):VSource(NumLyrs)),               &
                   Fields(Store(1):Store(NumLyrs)),                   &
                   ErrorStatus )


!---------------------------------------------------------------
!---- Step 5: Duplicate icing predictor field and set       ----
!---- fieldscodes in the field headers.                     ----
!---------------------------------------------------------------

    DO i=1, NumLyrs
! DEPENDS ON: dup_ice
      CALL Dup_ice( Fields(Store(i)),                                &
                    Fields(Store(NumLyrs+i)),                        &
                    ErrorStatus )
    END DO


!-----------------------------------------------------------------------------
!--- Output:                                                                 !
!---   Fields(Store(1):Store(NumLyrs) = fields marked as mean icing          !
!---   Fields(Store(NumLyrs+1):Store(2*NumLyrs) = fields marked as max icing !
!-----------------------------------------------------------------------------


  ELSE IF ( Action == "CLD_TURB" )  THEN

   ! 2.5.10 In-cloud turbulence algorithm

    IF ( P_Levels /= UM_NumLevs ) THEN
      WRITE(6,*) "Can only run LIBLKCLD on a 38 or 50 level model"
      ErrorStatus = StatusWarning
    END IF


    !Read in model fields
    NumLevs = ICT_NumLevs_Gl
    ZeroLev = ICT_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)
    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)


    !Read in temprature fields and calculate equivelent potential temperature fields

    ! Read theta level pressure in Pa
    STCode  (1:NumLevs) = ST_Ptheta
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, PSource, PPHdrMod, ModLevHdr_in, &
                  Fields, ErrorStatus )

    ! Read theta level temperature in K
    STCode  (1:NumLevs) = ST_Ttheta
    TSource (1:NumLevs) = TPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, TSource, PPHdrMod, ModLevHdr_in, &
                  Fields, ErrorStatus )


    !Call subroutine to calculate equivalent potential temperature
    !on model levels (put into VSource space in Fields)

    VSource (1:NumLevs) = VPosn + MO8Level(1:NumLevs)
! DEPENDS ON: eq_pot_temp
    CALL Eq_pot_temp( NumLevs,                                 &
                      Fields(TSource(1):TSource(NumLevs)),     &
                      Fields(PSource(1):PSource(NumLevs)),     &
                      Fields(VSource(1):VSource(NumLevs)),     &
                      ErrorStatus )

    ! Read theta level bulk cloud fraction
    STCode  (1:NumLevs) = ST_BlkCld
    USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, USource, PPHdrMod, ModLevHdr_in, &
                  Fields, ErrorStatus )


    ! Calculate theta level heights
    ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)
! DEPENDS ON: calczflds
    CALL CalcZFlds( NumLevs, OrogField,                         &
                    Fields(PSource(1):PSource(NumLevs)),        &
                    Fields(ZSource(1):ZSource(NumLevs)),        &
                    ErrorStatus )


   ! Do calculations

    MO8Level(1:ICT_NumStd) = (/ 700, 600, 500, 400, 300/)

    DO i = 1, ICT_NumStd

! DEPENDS ON: cld_turb
      CALL cld_turb( MO8Level(i),                           &  ! pressure in hPa (in)
                     NumLevs,                               &  ! no. model levels (in)
                     Fields(USource(1):USource(NumLevs)),   &  ! cloud on model level fields (in)
                     Fields(VSource(1):VSource(NumLevs)),   &  ! equivalnt potential temp (in)
                     Fields(PSource(1):PSource(NumLevs)),   &  ! corresponding pressure (in)
                     Fields(ZSource(1):ZSource(NumLevs)),   &  ! corresponding height (in)
                     Fields(Store(i)),                      &  ! output incloud turb predictor
                     ErrorStatus )                             !

      !Duplicate the field and label as the maximum
! DEPENDS ON: dup_ict
      CALL Dup_ICT( Fields(Store(i)),                       &  !input in-cloud turb field
                    Fields(Store(ICT_NumStd+i)),            &  !output duplicate field
                    ErrorStatus )
    END DO

  ELSE IF ( Action == "CB_ACT" )  THEN

   ! 2.5.11 Cumulonimbus (CB) fields

    IF ( P_Levels /= UM_NumLevs ) THEN
      WRITE(6,*) "Can only run LIBLKCLD on a 38 or 50 level model"
      ErrorStatus = StatusWarning
    END IF


! Read in standard level fields

! Fields(10) : Convective Cloud Base Pressure : 5/207
! Fields(11) : Convective Cloud Top  Pressure : 5/208
! Fields(12) : Convective Precipitation Rate  : 5/205
! Fields(13) : Lowest Conv Cloud Top Pressure : 5/223

! Each field is only on one level: level 8888

    STcode  (1:CB_NumStd) = (/  ST_CClBP,  ST_CClTP,            &
                                ST_CPNRT, ST_LCClTP /)
    MO8Level(1:CB_NumStd) = (/      8888,      8888,            &
                                    8888,      8888 /)
    Source  (1:CB_NumStd) = 9 + (/ (i, i=1,CB_NumStd) /)


! DEPENDS ON: getflds
    CALL GetFlds( CB_NumStd, MaxFlds, STCode, MO8Level, FCTime,       &
                  LBProc, MinsPastHr, Source, PPHdrMod, StdLevHdr_in, &
                  Fields, ErrorStatus )


! Read in model level fields

    NumLevs = IC_NumLevs_Gl
    ZeroLev = IC_ZeroLev_Gl
    FCTime  (1:NumLevs) = FCTime(1)
    MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)

!   TSource set later, ZSource not used here
!   TSource (1:NumLevs) = TPosn + MO8Level(1:NumLevs)
!   ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)

    ! Read theta level bulk cloud fraction (0/266)

    STCode  (1:NumLevs) = ST_BlkCld
    USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs,MaxFlds, STCode, MO8Level, FCTime,           &
                  LBProc, MinsPastHr, USource, PPHdrMod, ModLevHdr_in, &
                  Fields, ErrorStatus )

    ! Read Convective Cloud Amount (5/212)

    STCode  (1:NumLevs) = ST_ConCld
    VSource (1:NumLevs) = VPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs,MaxFlds, STCode, MO8Level, FCTime,           &
                  LBProc, MinsPastHr, VSource, PPHdrMod, ModLevHdr_in, &
                  Fields, ErrorStatus )

    ! Read theta level pressure in Pa (0/408)

    STCode  (1:NumLevs) = ST_Ptheta
    PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, PSource, PPHdrMod, ModLevHdr_in, &
                  Fields, ErrorStatus )

    ! Read theta level Temperature in K (16/004)

    STCode  (1:NumLevs) = ST_Ttheta
    TSource (1:NumLevs) = TPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
    CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                  LBProc, MinsPastHr, TSource, PPHdrMod, ModLevHdr_in, &
                  Fields, ErrorStatus )

! Do calculations

! DEPENDS ON: convact
    CALL ConvAct( NumLevs,Fields,Errorstatus)

  ELSE IF ( Action == "DUST_CONC" )  THEN

         ! 3.0.1 Calculate the Surf dust concentration

          IF ( ( P_Levels /= UM_NumLevs )                       .OR. &
               ((UM_NumLevs /= 38) .AND. (UM_NumLevs /= 50))      &
             ) THEN
             WRITE(6,*) "Can only run DUST_CONC on a 38 or 50 level model"
             ErrorStatus = StatusWarning
          END IF

          ! for the dust calc have chosen to read in bottom theta level

          NumLevs = DT_NumLevs_Gl
          ZeroLev = DT_ZeroLev_Gl
          FCTime  (1:NumLevs) = FCTime(1)
          MO8Level(1:NumLevs) = (/ (i+ZeroLev, i=1,NumLevs) /)

          ! Read rho level pressure  (0/407)

          STCode  (1:NumLevs) = ST_Prho
          USource (1:NumLevs) = UPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
          CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                        LBProc, MinsPastHr, USource, PPHdrMod, ModLevHdr_in, &
                        Fields, ErrorStatus )

          ! Calculate rho level heights

          ZSource (1:NumLevs) = ZPosn + MO8Level(1:NumLevs)
! DEPENDS ON: calczflds
          CALL CalcZFlds( NumLevs, OrogField,                   &
                    Fields(USource(1):USource(NumLevs)),        &
                    Fields(ZSource(1):ZSource(NumLevs)),        &
                    ErrorStatus )

          ! Read theta level pressure  (0/408)

          STCode  (1:NumLevs) = ST_Ptheta
          PSource (1:NumLevs) = PPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
          CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                        LBProc, MinsPastHr, PSource, PPHdrMod, ModLevHdr_in, &
                        Fields, ErrorStatus )

          ! Read theta level temperature  (16/004)

          STCode  (1:NumLevs) = ST_Ttheta
          TSource (1:NumLevs) = TPosn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
          CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,          &
                        LBProc, MinsPastHr, TSource, PPHdrMod, ModLevHdr_in, &
                        Fields, ErrorStatus )

          ! Read dust mmr 1 (0/431)

          STCode  (1:NumLevs) = ST_Dust1
          D1Source (1:NumLevs) = DST1Posn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
          CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,           &
                        LBProc, MinsPastHr, D1Source, PPHdrMod, ModLevHdr_in, &
                        Fields, ErrorStatus )

          ! Read dust mmr 2 (0/432)

          STCode  (1:NumLevs) = ST_Dust2
          D2Source (1:NumLevs) = DST2Posn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
          CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,           &
                        LBProc, MinsPastHr, D2Source, PPHdrMod, ModLevHdr_in, &
                        Fields, ErrorStatus )

          ! Read dust mmr 3 (0/433)

          STCode  (1:NumLevs) = ST_Dust3
          D3Source (1:NumLevs) = DST3Posn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
          CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,           &
                        LBProc, MinsPastHr, D3Source, PPHdrMod, ModLevHdr_in, &
                        Fields, ErrorStatus )

          ! Read dust mmr 4 (0/434)

          STCode  (1:NumLevs) = ST_Dust4
          D4Source (1:NumLevs) = DST4Posn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
          CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,           &
                        LBProc, MinsPastHr, D4Source, PPHdrMod, ModLevHdr_in, &
                        Fields, ErrorStatus )

          ! Read dust mmr 5 (0/435)

          STCode  (1:NumLevs) = ST_Dust5
          D5Source (1:NumLevs) = DST5Posn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
          CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,           &
                        LBProc, MinsPastHr, D5Source, PPHdrMod, ModLevHdr_in, &
                        Fields, ErrorStatus )

          ! Read dust mmr 6 (0/436)

          STCode  (1:NumLevs) = ST_Dust6
          D6Source (1:NumLevs) = DST6Posn + MO8Level(1:NumLevs)
! DEPENDS ON: getflds
          CALL GetFlds( NumLevs, MaxFlds, STCode, MO8Level, FCTime,           &
                        LBProc, MinsPastHr, D6Source, PPHdrMod, ModLevHdr_in, &
                        Fields, ErrorStatus )

          ! Calculate surface dust concentration

! DEPENDS ON: dust_con
          CALL DUST_CON (Numlevs, OrogField,                    &
                         Fields(PSource(1):PSource(NumLevs)),   &
                         Fields(TSource(1):TSource(NumLevs)),   &
                         Fields(ZSource(1):ZSource(NumLevs)),   &
                         Fields(D1Source(1):D1Source(NumLevs)), &
                         Fields(D2Source(1):D2Source(NumLevs)), &
                         Fields(D3Source(1):D3Source(NumLevs)), &
                         Fields(D4Source(1):D4Source(NumLevs)), &
                         Fields(D5Source(1):D5Source(NumLevs)), &
                         Fields(D6Source(1):D6Source(NumLevs)), &
                         Fields(Store(1)),Fields(Store(2)),     &
                         ErrorStatus )

  ELSE

   ! Unrecognised Action call
! DEPENDS ON: ereport
    CALL EReport( ProgName, StatusWarning, &
                  "Command not recognised.  Action = " // Action )
    ErrorStatus = StatusWarning

  END IF

  ! Check if current / previous error has caused Action to be aborted
  IF (ErrorStatus /= StatusOK) THEN
!DEPENDS ON: ereport
    CALL EReport( ProgName, StatusWarning, &
                  "Error - skipping to next Action.")
    ErrorStatus = StatusWarning
  END IF

END DO

9999 CONTINUE

!-----------------------------------------------------------------------
! 3 CLOSE INPUT FILES
!-----------------------------------------------------------------------
IF ( ASSOCIATED( OrogField % RData ) ) THEN
  DEALLOCATE( OrogField % RData )
END IF
DO i = 1, MaxFlds
  IF ( ASSOCIATED( Fields(i) % RData ) ) THEN
    DEALLOCATE( Fields(i) % RData )
    NULLIFY( Fields(i) % RData )
  END IF
END DO

CLOSE( NLUnit )
! DEPENDS ON: file_close
CALL File_Close ( StdLevHdr_in % UnitNum,                &
                  StdLevHdr_in % FileNameEnv,            &
                  LEN_TRIM(StdLevHdr_in % FileNameEnv),  &
                  0, 0, ErrorStatus )
! DEPENDS ON: file_close
CALL File_Close ( ModLevHdr_in % UnitNum,                &
                  ModLevHdr_in % FileNameEnv,            &
                  LEN_TRIM(ModLevHdr_in % FileNameEnv),  &
                  0, 0, ErrorStatus )
IF ( ASSOCIATED( StdPrvHdr_in % Lookup ) ) THEN
! DEPENDS ON: file_close
  CALL File_Close ( StdPrvHdr_in % UnitNum,              &
                    StdPrvHdr_in % FileNameEnv,          &
                    LEN_TRIM(StdPrvHdr_in % FileNameEnv),&
                    0, 0, ErrorStatus )
END IF
!-----------------------------------------------------------------------
! 4 CLOSE OUTPUT FILE
!-----------------------------------------------------------------------
! Write out updated lookup table
! DEPENDS ON: writelookup
CALL WriteLookup( StdLevHdr_out, ErrorStatus )
WRITE(6,*) StdLevHdr_out % NumFlds, " fields written."

! DEPENDS ON: file_close
CALL File_Close ( StdLevHdr_out % UnitNum,               &
                  StdLevHdr_out % FileNameEnv,           &
                  LEN_TRIM(StdLevHdr_out % FileNameEnv), &
                  0, 0, ErrorStatus )

! DEPENDS ON: timer
CALL Timer( ProgName, 2 )

END PROGRAM FieldCalc
#endif
