
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Top level routine for reading, handling and writeing GRIB data.

Module Rcf_Grib_Control_Mod

! SUBROUTINE Rcf_Grib_Control
!
! Description: This is the top level routine for reading in GRIB data
!              and generating a UM style intermediary dump containing
!              that data
!
! Method: read in the GRIB records one at a time storing the information
!         in a set of dynamically allocated lists. Using this info, set
!         up the corresponding UM style headers.
!         re-read records from the GRIB file, this time performing
!         simple transformations (if necessary) on the data before
!         writing it out to the intermediary UM dump.
!
! Current Code Owner: Roddy Sharp
!
! History:
! Version   Date      Comment
! -------  --------   -------------------------
!  5.4     12/06/02   Original code. Roddy Sharp (frtz)
!  5.5     21/01/03   Adding 'file close' & 'free unit' R.Sharp
!  6.2     20/06/05   Added initialisation of vertical coordinates
!                     array. Paul Earnshaw (frpe)
!  6.3     22/09/06   Fix problems with deallocating Hdr_Dmy when
!                     compiled under FCM. P.Selwood
!  6.4     22/01/07   Add ozone to vertical level processing
!                     C.Mathison
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
contains
Subroutine Rcf_Grib_Control()

! uses variables and routines from other modules

Use Rcf_Parvars_mod, Only : &
    mype

Use Rcf_GRIB_Block_Params_Mod

Use Rcf_GRIB_FldSort_Mod, Only  :  &
  Rcf_GRIB_FldSort

Use Rcf_StashCodes_Mod

Use Rcf_HeadAddress_Mod

Use Rcf_GRIB_Lookups_Mod               ! Contains cross ref table for
                                       ! Stash => ECMWF parameter Id's
                                       ! and appropriate params

Use Rcf_Grib_Debug_Tools_Mod, Only : &
     Grib_Debug_Print_Basics,        &
     Grib_Debug_Print_Blocks,        &
     Grib_Debug_ListCounts

Use Rcf_Grib_Assign_Mod, Only : &
    Rcf_Grib_Assign

Use Rcf_Grib_Check_Mod, Only : &
    Rcf_Grib_Check

Use Rcf_Grib_SetHdr_Mod, Only : &
    Rcf_Grib_SetHdr

Use Rcf_Grib_Dest_List_Mod, Only : &
    Rcf_Grib_Dest_List

Use Rcf_UMhead_Mod, Only :  &
    LenFixHd,               &
    um_header_type            ! Derived containing UM header info

Use EReport_Mod, Only :     &
    EReport

Use Rcf_WriteUMhdr_Mod, Only :     &
    Rcf_WriteUMhdr

Use SetPos8_Mod, Only :     &
    SetPos8

Use Buffin8_Mod, Only :     &
    BuffIn8

Use Rcf_FreeUMhdr_Mod, Only :     &
    Rcf_FreeUMhdr

Use Rcf_FortranIO_Mod, Only : &
    Rcf_Get_Unit,             &
    RCF_Free_Unit

Use Rcf_Grib_Read_Data_Mod, Only : &
    Rcf_Grib_Read_Data

Use Rcf_Grib_Spcl_Ctl_Mod, Only : &
    Rcf_Grib_Spcl_Ctl

Use Rcf_Grib_Spcl_Hdr_Mod, Only : &
    Rcf_Grib_Spcl_Hdr

Use Rcf_Grid_Type_Mod, Only : &
    Grid_type,              &
    Output_Grid

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Min,           &       ! =1 Minimum output
    PrStatus_Normal,        &       ! =2 Short informative output
    PrStatus_Oper,          &       ! =3 Full informative output
    PrStatus_Diag                   ! =4 Extra Diagnostic output

Implicit None

! Comdecks
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------
!*L------------------ COMDECK LOOKADD ----------------------------------
!LL
!LL Purpose : Contains information about the format
!LL           of the PP header
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   4.0  12/09/95   Change NPERIODS to LBUSER3, BRSVD1 to BULEV,
!LL                   BRSVD2 to BHULEV and definitions for BRLEV and
!LL                   BHRLEV. Corresponding changes made to STWORK1A
!LL                   and PPHEAD1A. (Andrew Brady)
!LL  4.0  12/10/95  Change item 45 from lbuser7 to model_code. RTHBarnes
!LL  5.1  17/04/00    Fixed/Free format. P.Selwood.
!LL  5.2  25/09/00    Add LBCC_xxxx variables for the compressed
!LL                   LBC LOOKUP array                  P.Burton
!LL
!LL Programming standard :
!LL
!LL Logical components covered : F092
!LL
!LL Project task :
!LL
!LL External documentation:
!LL
!LLEND -----------------------------------------------------------------
!

! Validity time
      Integer, Parameter :: LBYR   =1   ! Year
      Integer, Parameter :: LBMON  =2   ! Month
      Integer, Parameter :: LBDAT  =3   ! Day of month
      Integer, Parameter :: LBHR   =4   ! Hour
      Integer, Parameter :: LBMIN  =5   ! Minute
      Integer, Parameter :: LBDAY  =6   ! Day number

! Data time
      Integer, Parameter :: LBYRD  =7   ! Year
      Integer, Parameter :: LBMOND =8   ! Month
      Integer, Parameter :: LBDATD =9   ! Day of month
      Integer, Parameter :: LBHRD  =10  ! Hour
      Integer, Parameter :: LBMIND =11  ! Minute
      Integer, Parameter :: LBDAYD =12  ! Day number

      Integer, Parameter :: LBTIM  =13  ! Time indicator
      Integer, Parameter :: LBFT   =14  ! Forcast period (hours)
      Integer, Parameter :: LBLREC =15  ! Length of data record
      Integer, Parameter :: LBCODE =16  ! Grid type code
      Integer, Parameter :: LBHEM  =17  ! Hemisphere indicator
      Integer, Parameter :: LBROW  =18  ! Number of rows in grid
      Integer, Parameter :: LBNPT  =19  ! Number of points per row
      Integer, Parameter :: LBEXT  =20  ! Length of extra data
      Integer, Parameter :: LBPACK =21  ! Packing method indicator
      Integer, Parameter :: LBREL  =22  ! Header release number
      Integer, Parameter :: LBFC   =23  ! Field code
      Integer, Parameter :: LBCFC  =24  ! Second field code
      Integer, Parameter :: LBPROC =25  ! Processing code
      Integer, Parameter :: LBVC   =26  ! Vertical coordinate type
      Integer, Parameter :: LBRVC  =27  ! Coordinate type for reference
                                        ! level

      Integer, Parameter :: LBEXP  =28  ! Experiment number
      Integer, Parameter :: LBEGIN =29  ! Start record
      Integer, Parameter :: LBNREC =30  ! No of records-Direct access
                                        ! only
      Integer, Parameter :: LBPROJ =31  ! Met-O-8 projection number
      Integer, Parameter :: LBTYP  =32  ! Met-O-8 field type
      Integer, Parameter :: LBLEV  =33  ! Met-O-8 level code
      Integer, Parameter :: LBRSVD1=34  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD2=35  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD3=36  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD4=37  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBSRCE =38  ! =1111 to indicate following
                                        ! apply to UM
      Integer, Parameter :: DATA_TYPE =39  ! Indicator for real/int or
                                           ! timeseries
      Integer, Parameter :: NADDR  =40  ! Start address in DATA_REAL or
                                        ! DATA_INT
      Integer, Parameter :: LBUSER3=41  ! Free for user-defined function
      Integer, Parameter :: ITEM_CODE =42  !Stash code
      Integer, Parameter :: LBPLEV =43  ! Pseudo-level indicator (if
                                        ! defined)
      Integer, Parameter :: LBUSER6=44  ! Free for user-defined function
      Integer, Parameter :: MODEL_CODE =45 ! internal model identifier

      Integer, Parameter :: BULEV  =46  ! Upper level boundary
      Integer, Parameter :: BHULEV =47  ! Upper level boundary
      Integer, Parameter :: BRSVD3 =48  ! Reserved for future PP-package
                                        ! use
      Integer, Parameter :: BRSVD4 =49  ! Reserved for future PP-package
                                        ! use
      Integer, Parameter :: BDATUM =50  ! Datum value
      Integer, Parameter :: BACC   =51  ! (Packed fields) Packing
                                        ! accuracy
      Integer, Parameter :: BLEV   =52  ! Level
      Integer, Parameter :: BRLEV  =53  ! Lower level boundary
      Integer, Parameter :: BHLEV  =54  ! (Hybrid levels) A-level of
                                        ! value
      Integer, Parameter :: BHRLEV =55  ! Lower level boundary
      Integer, Parameter :: BPLAT  =56  ! Real latitude of 'pseudo'
                                        ! N Pole
      Integer, Parameter :: BPLON  =57  ! Real longitude of 'pseudo'
                                        ! N Pole
      Integer, Parameter :: BGOR   =58  ! Grid orientation
      Integer, Parameter :: BZY    =59  ! Zeroth latitude
      Integer, Parameter :: BDY    =60  ! Latitude interval
      Integer, Parameter :: BZX    =61  ! Zeroth longitude
      Integer, Parameter :: BDX    =62  ! Longitude interval
      Integer, Parameter :: BMDI   =63  ! Missing data indicator
      Integer, Parameter :: BMKS   =64  ! M,K,S scaling factor

      Integer, Parameter :: LBCC_LBYR   = 1  ! Year
      Integer, Parameter :: LBCC_LBMON  = 2  ! Month
      Integer, Parameter :: LBCC_LBDAT  = 3  ! Day of the month
      Integer, Parameter :: LBCC_LBHR   = 4  ! Hour
      Integer, Parameter :: LBCC_LBMIN  = 5  ! Minute
      Integer, Parameter :: LBCC_LBDAY  = 6  ! Day number
      Integer, Parameter :: LBCC_LBEGIN = 7  ! Start record
      Integer, Parameter :: LBCC_NADDR  = 8  ! Start address of DATA
! Mapping of MPP_LOOKUP; analogous to mapping in PP header

      Integer, Parameter :: P_NADDR=1    ! Address on local PE
      Integer, Parameter :: P_LBLREC=2   ! Local length of record

!*----------------------------------------------------------------------
! NADDR IS LOCATION IN PP-HEADER (LOOKUP) FOR START POSN OF VARIABLE
! ITEM_CODE is the location in PP header for a code defined as
!           (section number)*1000+item number
! DATA_TYPE is the location in the PP header defining data as REAL or
!           INTEGER.
! LBNPT is the location defining the number of points per row
!
! contains LBLREC (amongst others)

! Local variables

Type (Um_Header_type)            :: Hdr_Dmy
Type (Um_Header_type)            :: Hdr_Itm
! local var to preserve Output_Grid values during 1st stage
Type (Grid_type)                 :: Storage_Grid

Character (Len=*), Parameter     :: RoutineName='Rcf_GRIB_Control'
Character (Len=80)               :: Cmessage(2)   ! used for EReport
Character (Len=20)               :: cFormat
Integer                          :: ErrorStatus   ! used for EReport
Integer                          :: errstat
Integer                          :: I
Integer                          :: Count,Criteria
Integer                          :: dummy1
Integer                          :: dummy2
Integer                          :: disk_address

!-----------------------------------------------------------------------
!  Variables to do specifically with the GRIB record
!-----------------------------------------------------------------------

Type (Grib_Record),Pointer       :: Current       !\ Pointer to current
                                                  !/ grib record

!An array of pointer pairs to form the head and tail of the field lists
Type (List_Marker)               :: Lists(0:grib_max_fields)

Integer, Parameter               :: len_max_c = 256
Character (Len=len_max_c)        :: c_Tmp_Buff(1)
Integer                          :: act_io_len
Integer                          :: pos_in_file

Logical                          :: Order
Logical, Parameter               :: Ascending  = .True.
Logical, Parameter               :: Descending = .False.

!=======================================================================
!  Initialise variables
!=======================================================================

Storage_grid = Output_Grid

Nullify(Current)

! Allocate a dummy header (used to hold values other routines would copy
! from the input dump) and the header for the Intermediate dump
Allocate (Hdr_Dmy % FixHd(LenFixHd))
Allocate (Hdr_Itm % FixHd(LenFixHd))

! Nullify the components of Hdr_Dmy that are not being set anywhere
! to avoid later problems with deallocation.
Nullify (Hdr_Dmy % IntC)
Nullify (Hdr_Dmy % CompFldI1)
Nullify (Hdr_Dmy % CompFldI2)
Nullify (Hdr_Dmy % CompFldI3)
Nullify (Hdr_Dmy % Lookup)
Nullify (Hdr_Dmy % RealC)
Nullify (Hdr_Dmy % LevDepC)
Nullify (Hdr_Dmy % RowDepC)
Nullify (Hdr_Dmy % ColDepC)
Nullify (Hdr_Dmy % FldsOfC)
Nullify (Hdr_Dmy % ExtraC)
Nullify (Hdr_Dmy % HistFile)

! Make sure the list pointers are NULL and the counts are zero
Do I = 0,grib_max_fields
  Nullify( Lists(I) % Begin )
  Nullify( Lists(I) % End )
  Lists(I) % LstCount = 0
End Do

pos_in_file = 3   ! Offset to allow for small rewind in each loop

!=======================================================================
!  Open GRIB File
!=======================================================================

Call Rcf_Get_Unit( Hdr_Dmy % UnitNum )

! DEPENDS ON: file_open
Call File_Open( Hdr_Dmy % UnitNum, 'AINITIAL', 8, 0, 0, errstat)

If ( errstat /= 0 ) Then
  Cmessage(1) = 'Failed to Open GRIB file'
  ErrorStatus = 10
  Call EReport( RoutineName, ErrorStatus, Cmessage(1) )
Else
  If ( PrintStatus >= PrStatus_Diag  ) Then
    If ( mype == 0 ) Then
      Write (6,'(A)') "Opened Grib Data File"
    End If
  End If
End If

!=======================================================================
!  Loop, reading blocks from file, until no more records found
!=======================================================================

Reading: Do                      ! Endless loop until exit supplied

  pos_in_file = pos_in_file - 3  ! Small rewind to allow Index to seek
                                 ! for GRIB properly

!=======================================================================
!  Find the begining of a record by seeking 'GRIB'
!=======================================================================

  ! Set the position within the file to 'pos_in_file'.
  Call SetPos8(Hdr_Dmy % UnitNum, pos_in_file,errstat)
  If ( errstat /= 0 ) Then
    Write(Cmessage(1),'(A,I2)')   'Failed trying to SetPos to ',      &
                                  pos_in_file
    ErrorStatus = 20
    Call EReport( RoutineName, ErrorStatus, Cmessage(1) )
  End If

  ! Clear the buffer and read in section of file to search for 'GRIB'
  c_Tmp_Buff(1) = " "
  Call Buffin8 (Hdr_Dmy % UnitNum, c_Tmp_Buff(1), 1, len_max_c,       &
                act_io_len, errstat)

  If ( errstat > 0 ) Then
    Cmessage(1)    = 'Failed trying to read short block from GRIB file'
    ErrorStatus = 30
    Call EReport( RoutineName, ErrorStatus, Cmessage(1) )
  End If

  I = Index (c_Tmp_Buff(1), "GRIB")   ! Serch the block for text 'GRIB'

!=======================================================================
!        -> _Didn't_ Find 'GRIB' - try again
!=======================================================================

  If (I == 0) Then               ! Index didn't find GRIB in the buffer

    If (PrintStatus >= PrStatus_Diag) Then
      If ( mype == 0 ) Then
        Write (6,*) "Buffer block read contained no 'GRIB' indicator"
      End If
    End If

    ! If data read is smaller than data requested then must have read
    ! until the end of the file. => Exit Do loop
    If (len_max_c > act_io_len) Then
      If (PrintStatus >= PrStatus_Diag) Then
        If ( mype == 0 ) Then
          Write (6,'(2(A,I3),A)') 'len_max_c =', len_max_c,           &
                 ' act_io_len =', act_io_len, 'Exiting from Read loop'
        End If
      End If

      Exit Reading             ! This forces exit from reading the loop
    Else
      ! Force (if some what slow) progress through the file
      pos_in_file = pos_in_file + act_io_len
    End If

!=======================================================================
!        -> _Found_ 'GRIB' - read and decode data
!=======================================================================

  Else                                   ! GRIB WAS found

!=======================================================================
!  Allocate 'GRIB Record' to store header Info
!=======================================================================

    Allocate(Current)

    Nullify(Current % VertCoords)

    Current % Block_1(:) = 0
    Current % Block_2(:) = 0
    Current % Block_3(:) = 0
    Current % Block_4(:) = 0
    Current % Block_R(:) = 0.000

    Current % Start_pos  = 0
    Current % StashCode  = -1
    Current % Data_Type  = Grb_Data_Real
    Current % Num_Fp     = 0
    Current % Num_Vert   = 0
    Current % Num_Bitmap = 0
    Current % Num_Quasi  = 0

!=======================================================================
!  Read binary form of record into CRecord work array
!=======================================================================

    pos_in_file = pos_in_file + I - 1      ! set the pos in the file
    Current % Start_pos = pos_in_file      ! record start pos in header

    Call Rcf_Grib_Read_Data(Hdr_Dmy % UnitNum,Current,FpData,         &
                            LenArrayMax,pos_in_file)

!=======================================================================
!  Test record and assign to correct list
!=======================================================================

    Call Rcf_Grib_Assign(Current, Lists)

   ! Set position in file to end of record read
    pos_in_file = pos_in_file + Current % Block_0(p_Mes_Len)

!=======================================================================
!  End Loop "until no more records found"
!=======================================================================

  End If                               ! Test to find 'GRIB' in buffer

End Do Reading

!=======================================================================
!  Sort the 'Lists' into sensible orders
!=======================================================================

cFormat = "(3A,I2,A)"
Do I = 1, grib_max_fields

  If (Associated(Lists(I) % Begin) .AND.                              &
      Lists(I) % LstCount > 1) Then    ! Checks list has members

    ! For each case, set the order desired, the criteria (currently
    ! any value in Block 1) on which to perform the ordering, and
    ! the output message text.
    Select Case(I)

      ! Pressure level fields ordered away from surface.
      Case (grib_U_field,                                             &
            grib_V_field,                                             &
            grib_W_field,                                             &
            grib_Temp_field,                                          &
            grib_Q_field,                                             &
            grib_ozone_field,                                         &
            grib_NOX_field,                                           &
            grib_CH4_field,                                           &
            grib_CO_field,                                            &
            grib_HCHO_field,                                          &
            grib_GO3_field)

        Order = Descending
        Criteria = p_Lvl_Desc_1     ! 1st Level description parameter
        Write(cMessage(1),cFormat) "Re-ordering ",                    &
                            Trim(Lists(I) % Begin % Desc),            &
                           ". ", Lists(I) % LstCount,  " entries"

      ! Soil levels ordered down from surface 
      Case (grib_Soil_Temp_field,                                     &
            grib_Soil_Moist_field)

        Order = Ascending
        Criteria = p_Lvl_Desc_1     ! 1st Level description parameter
        Write(cMessage(1),cFormat) "Re-ordering ",                    &
                            Trim(Lists(I) % Begin % Desc),            &
                           ". ", Lists(I) % LstCount,  " entries"

      Case Default
        Order = Descending
        Criteria = p_Lvl_Desc_1
        Write(cMessage(1),'(A,I2,A,I2,A)')                            &
                           "List ", I, ", ", Lists(I) % LstCount,     &
        " entries re-ordered using defaults. No Specific rule exists"
    End Select

    Call rcf_Grib_FldSort(Lists(I), Criteria, Order)

  Else
    Write(cMessage(1),*) "List ", I ,                                 &
      " Not re-ordered because it contained ", Lists(I) % LstCount,   &
      " entry"
  End If

  If (PrintStatus >= PrStatus_Diag) Then
    If ( mype == 0 ) Then
      Write (6,'(A)') cMessage(1)
    End If
  End If

End Do ! I over no. of fields

!=======================================================================
! Loop through lists showing count
!=======================================================================

  If ( PrintStatus >= PrStatus_Diag .And. mype == 0 ) Then
    Write (6,*)
    Call Grib_Debug_ListCounts(Lists)
    Write (6,*)
  End If

!=======================================================================
! Loop through lists displaying basic info
!=======================================================================

If ( PrintStatus >= PrStatus_Diag .And. mype == 0 ) Then
  Write (6,*)
  Call Grib_Debug_Print_Basics(Lists)
  Write (6,*)
End If

!=======================================================================
! Perform Basic Data Integrity Checks
!=======================================================================

Call Rcf_Grib_Check(Lists)

!=======================================================================
! Setting up the headers.
!=======================================================================

Call Rcf_Grib_SetHdr(Lists,Output_Grid,Hdr_Dmy,Hdr_Itm)

!=======================================================================
!  Open the Output Dump for the Intermediary Dump
!=======================================================================

Call Rcf_Get_Unit( Hdr_Itm % UnitNum )

If ( PrintStatus >= PrStatus_Diag  ) Then
  If ( mype == 0 ) Then
    Write (6,*) "Opening Intermediary Dump"
  End If
End If
! DEPENDS ON: file_open
Call File_Open( Hdr_Itm % UnitNum, 'RECONTMP', 8, 1, 0, errstat)

If ( errstat /= 0 ) Then
  Cmessage(1)    = 'Failed to Open Intermediary Dump file'
  ErrorStatus = 99
  Call EReport( RoutineName, ErrorStatus, Cmessage(1) )
End If

!=======================================================================
!  Set up addressing info used when writing Data
!=======================================================================

! DEPENDS ON: set_dumpfile_address
Call Set_Dumpfile_Address( Hdr_Itm % FixHd,                           &
                           LenFixHd,                                  &
                           Hdr_Itm % Lookup,                          &
                           Hdr_Itm % Len1Lookup,                      &
                           Hdr_Itm % Len2Lookup,                      &
                           dummy1, dummy2, disk_address)

!=======================================================================
!  Write the Field data
!=======================================================================

count = 0

!Loop across all lists
Do I = 1, grib_max_fields

  If (Associated(Lists(I) % Begin) ) Then

    Current => Lists(I) % Begin
    Do While (Associated(Current))

      count = count + 1

!=======================================================================
!  Read binary form of record into CRecord work array
!=======================================================================

      Call Rcf_Grib_Read_Data(Hdr_Dmy % UnitNum,                 &
                              Current,FpData,                    &
                              Current % Block_0(p_Mes_Len),      &
                              Current % Start_pos)

      ! Call to routine which calls for 'special' handling of given data

      Call Rcf_Grib_Spcl_Ctl(FpData,LgData,Lists,Current,             &
                             Hdr_Dmy,Hdr_Itm,I,count)

      ! Call to Write the field data
      ! .False. Indicates data is held on a single PE
      If (Current % Data_Type == Grb_Data_Real ) Then
! DEPENDS ON: rcf_writflds
        Call Rcf_WritFlds( Hdr_Itm % UnitNum, 1,count,                &
                   Hdr_Itm % Lookup, Hdr_Itm % Len1Lookup, FpData,    &
                   Hdr_Itm % Lookup( LBLREC , count ),                &
                   Hdr_Itm % FixHd,   &
                   ErrorStatus, Cmessage, .False. )

      Else If (Current % Data_Type == Grb_Data_Log ) Then
! DEPENDS ON: rcf_writflds
        Call Rcf_WritFlds( Hdr_Itm % UnitNum, 1,count,                &
                   Hdr_Itm % Lookup, Hdr_Itm % Len1Lookup, LgData,    &
                   Hdr_Itm % Lookup( LBLREC , count ),                &
                   Hdr_Itm % FixHd,   &
                   ErrorStatus, Cmessage, .False. )

      Else
        Cmessage(1) = 'Failed to write to temp dump :Unknown Data Type'
        ErrorStatus = 99
        Call EReport( RoutineName, ErrorStatus, Cmessage(1) )
      End If

      Current => Current % Next

    End Do  ! members of a list

  End If ! Associated(Lists(I) % Begin)

End Do  ! (I) loop over all lists

!=======================================================================
!  Write the Header on the Intermediary Dump
!=======================================================================

Call Rcf_WriteUMhdr( Hdr_Itm )

!=======================================================================
!  Closing the Grib File
!=======================================================================

! DEPENDS ON: file_close
Call File_Close ( Hdr_Dmy % UnitNum, 'AINITIAL', 8, 0, 0, errstat)
Call Rcf_Free_Unit ( Hdr_Dmy % UnitNum )

! DEPENDS ON: file_close
Call File_Close ( Hdr_Itm % UnitNum, 'RECONTMP', 8, 0, 0, errstat)
Call Rcf_Free_Unit ( Hdr_Itm % UnitNum )

!=======================================================================
!  Quick tidy up behind ourselves.
!=======================================================================

! Free up space taken by the headers
Call Rcf_FreeUMhdr ( Hdr_Dmy )
Call Rcf_FreeUMhdr ( Hdr_Itm )

! Deallocate the lists
Do I = 0, grib_max_fields
  Call Rcf_Grib_Dest_List (Lists(I))
End Do

! return Output_Grid values to those read in by namelist
Output_Grid = Storage_Grid

Return

End Subroutine Rcf_Grib_Control
End Module Rcf_Grib_Control_Mod
