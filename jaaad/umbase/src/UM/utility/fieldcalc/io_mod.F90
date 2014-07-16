#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module containing derived types used for IO in Fieldcalc

MODULE IO_Mod

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

IMPLICIT None

! Global Constants:
INTEGER, PARAMETER :: LenFixHd = 256   ! Length of Fixed_Length_Header
INTEGER, PARAMETER :: LenWord  = 64

! Global Type Definitions:
TYPE PP_Header_type

  INTEGER :: ValidYear    ! 1  LBYR   : Year       \ .
  INTEGER :: ValidMonth   ! 2  LBMON  : Month       \ .
  INTEGER :: ValidDate    ! 3  LBDAT  : Day          \  Validity time
  INTEGER :: ValidHour    ! 4  LBHR   : Hour         / .
  INTEGER :: ValidMin     ! 5  LBMIN  : Minute      / .
  INTEGER :: ValidDayNo   ! 6  LBDAY  : Day number / .
  INTEGER :: DataYear     ! 7  LBYRD  : Year       \ .
  INTEGER :: DataMonth    ! 8  LBMOND : Month       \ .
  INTEGER :: DataDate     ! 9  LBDATD : Day          \ Data time
  INTEGER :: DataHour     ! 10 LBHRD  : Hour         / .
  INTEGER :: DataMin      ! 11 LBMIND : Minute      / .
  INTEGER :: DataDayNo    ! 12 LBDAYD : Day number / .
  INTEGER :: LBTim        ! 13          Time indicator
  INTEGER :: FCRange      ! 14 LBFT   : Forecast range (hrs)
  INTEGER :: LBLRec       ! 15          Length of data record
  INTEGER :: LBCode       ! 16          Grid type code
  INTEGER :: LBHem        ! 17          Hemisphere indicator
  INTEGER :: NumRows      ! 18 LBROW  : Number of rows in grid
  INTEGER :: NumCols      ! 19 LBNPT  : Number of columns in grid
  INTEGER :: LBExt        ! 20          Length of extra data
  INTEGER :: LBPack       ! 21          Packing method indicator
  INTEGER :: LBRel        ! 22          Header release number
  INTEGER :: PPCode       ! 23 LBFC   : Field code
  INTEGER :: LBCFC        ! 24          Second field code
  INTEGER :: LBProc       ! 25          Processing code
  INTEGER :: LBVC         ! 26          Vertical co-ordinate type
  INTEGER :: LBRVC        ! 27          Vert co-ord type for reference
  INTEGER :: LBExp        ! 28          Experiment number
  INTEGER :: DataPos      ! 29 LBEGIN : Word address of data
  INTEGER :: LBNRec       ! 30          Disk length / No. records
  INTEGER :: MO8Proj      ! 31 LBPROJ : MetO8 Projection number
  INTEGER :: MO8Type      ! 32 LBTYP  : MetO8 Field Type
  INTEGER :: MO8Level     ! 33 LBLEV  : MetO8 Level code
  INTEGER :: LBRsvd(4)    ! 34-37       PP-Package use
  INTEGER :: LBSrce       ! 38          Version & Model
  INTEGER :: LBUser1      ! 39          Data type (real, etc.) (UM only)
  INTEGER :: LBUser2      ! 40          Start address in DATA (UM only)
  INTEGER :: LBUser3      ! 41          No. periods for timeseries
  INTEGER :: STCode       ! 42          STASH code (UM only)
  INTEGER :: LBUser5      ! 43          STASH pseudo dimension (UM only)
  INTEGER :: LBUser6      ! 44          Free (used by FIELDCOS)
  INTEGER :: LBUser7      ! 45          Sub-model number (UM only)
  REAL    :: BULev        ! 46          Upper layer boundary
  REAL    :: BHULev       ! 47          Upper layer boundary
  REAL    :: BRsvd(2)     ! 48-49       PP-Package use
  REAL    :: BDatum       ! 50          Datum value
  REAL    :: BAcc         ! 51          Packing accuracy
  REAL    :: RLevel       ! 52 BLEV   : Level (B or Zsea for hybrid)
  REAL    :: RefLevel     ! 53 BRLEV  : Reference level - lower layer bd
  REAL    :: BHLev        ! 54          Level (A or C for hybrid)
  REAL    :: BHRLev       ! 55          Lower layer boundary
  REAL    :: PseudoLat    ! 56 BPLAT  : Real lat of 'pseudo' N Pole
  REAL    :: PseudoLon    ! 57 BPLON  : Real lon of 'pseudo' N Pole
  REAL    :: BGOR         ! 58          Grid orientation
  REAL    :: ZerothLat    ! 59 BZY    : Zeroth latitude
  REAL    :: LatInt       ! 60 BDY    : Latitude interval
  REAL    :: ZerothLon    ! 61 BZX    : Zeroth longitude
  REAL    :: LonInt       ! 62 BDX    : Longitude interval
  REAL    :: BMDI         ! 63          Missing data indicator
  REAL    :: BMKS         ! 64          M.K.S. scaling factor

END TYPE PP_Header_type

TYPE PP_Field_type

  INTEGER              :: LookupPos   ! No of header in lookup
  INTEGER              :: ArrayPos    ! No of header in Fields array
  REAL, POINTER        :: RData(:,:)
  TYPE(PP_Header_type) :: Hdr

END TYPE PP_Field_type

TYPE UM_Header_type

  INTEGER :: LenIntC       ! Length of Integer_Constants array
  INTEGER :: LenRealC      ! Length of Real_Constants array
  INTEGER :: Len1LevDepC   ! 1st dim \ Level_Dependent_Constants
  INTEGER :: Len2LevDepC   ! 2nd dim / array
  INTEGER :: Len1RowDepC   ! 1st dim \ Row_Dependent_Constants
  INTEGER :: Len2RowDepC   ! 2nd dim / array
  INTEGER :: Len1ColDepC   ! 1st dim \ Column_Dependent_Constants
  INTEGER :: Len2ColDepC   ! 2nd dim / array
  INTEGER :: Len1FldsOfC   ! 1st dim \ Fields_Of_Constants
  INTEGER :: Len2FldsOfC   ! 2nd dim / array
  INTEGER :: LenExtraC     ! Length of Extra_Constants array
  INTEGER :: LenHistFile   ! Length of Temp_History_File
  INTEGER :: LenCompFldI1  ! Length of Compressed_Field_Index1
  INTEGER :: LenCompFldI2  ! Length of Compressed_Field_Index2
  INTEGER :: LenCompFldI3  ! Length of Compressed_Field_Index3
  INTEGER :: Len1Lookup    ! 1st dim \ Lookup table
  INTEGER :: Len2Lookup    ! 2nd dim /
  INTEGER :: LenData       ! Length of Data array
  INTEGER :: StartData     ! Position of start of Data array
  INTEGER :: NumFlds       ! Number of Data fields
  INTEGER :: UnitNum       ! Unit number associated with UM dump

  INTEGER, POINTER :: FixHd(:)      ! Fixed_Length_Header
  INTEGER, POINTER :: IntC(:)       ! Integer_Constants array
  INTEGER, POINTER :: CompFldI1(:)  ! Compressed_Field_Index1 array
  INTEGER, POINTER :: CompFldI2(:)  ! Compressed_Field_Index2 array
  INTEGER, POINTER :: CompFldI3(:)  ! Compressed_Field_Index3 array
  TYPE(PP_Header_type), POINTER :: Lookup(:) ! Lookup table

  REAL, POINTER :: RealC(:)
  REAL, POINTER :: LevDepC(:)
  REAL, POINTER :: RowDepC(:)
  REAL, POINTER :: ColDepC(:)
  REAL, POINTER :: FldsOfC(:)
  REAL, POINTER :: ExtraC(:)
  REAL, POINTER :: HistFile(:)

  CHARACTER(LEN=16) :: FileNameEnv

END TYPE UM_Header_type

END MODULE IO_Mod

#endif
