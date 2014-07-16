#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Read, from file, and Decode 1 field's worth of GRIB data

Module Rcf_Grib_Read_Data_Mod

! SUBROUTINE Rcf_Grib_Read_Data : Read raw info from file
!
! Description:
!   This routine is used to handle the reading of raw GRIB encoded data
!   in. This data is then passed to DECODE to be decoded.
!
! Method:
!   Recieve position markers and Storage arrays (or Types) from
!   calling routine.
!   Use SETPOS8 to set position in file to start of GRIB record.
!   Read a block of Raw data.
!   Call DECODE to decode raw Data.
!   Store the decoded data
!   Use the decoded length to calculate the next start point and move on
!
! Current Code Owner: Roddy Sharp
!
! History:
! Version   Date      Comment
! -------  --------   -------------------------
!  5.4     09/07/02   Original code. Roddy Sharp (frtz)
!  6.2     20/06/05   Added ability to save vertical coordinate
!                     information from GRIB file. Paul Earnshaw (frpe)
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
Contains
Subroutine Rcf_Grib_Read_Data(Unit_Num ,Current,Data,Len_Data,        &
                              Pos_in_File)

Use Rcf_GRIB_Block_Params_Mod, Only : &
    Grib_Record,     &       ! derived type for storing a GRIB header
    p_B4Undef

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Min,           &       ! =1 Minimum output
    PrStatus_Normal,        &       ! =2 Short informative output
    PrStatus_Oper,          &       ! =3 Full informative output
    PrStatus_Diag                   ! =4 Extra Diagnostic output

Use Rcf_Parvars_mod, Only : &
    mype

Use EReport_Mod, Only :     &
    EReport                         ! Subroutine to generate error msgs

Use SetPos8_Mod, Only :     &
    SetPos8                         ! Subroutine sets position in a file

Use Buffin8_Mod, Only :     &
    BuffIn8                         ! Subroutine to block read from file

Implicit None

! Subroutine arguments

!< Scalar arguments with intent(in):>
Integer, Intent(In)              :: Unit_Num  ! unit number file is on
Integer, Intent(In)              :: Len_Data  ! Length of data to read
Integer, Intent(In)              :: Pos_in_File ! position to start
                                                ! reading data from

!< Array  arguments with intent(InOut):>
Type (Grib_Record),Pointer       :: Current        !\ Pointer to current
                                                   !/ grib record
!< Array  arguments with intent(out):>
Real, Intent(Out)                :: Data(Len_Data) ! Array deocded field
                                                   ! is returned in.

! Local constants
! This was usedfor debugging. It allowed all non debugging messages to
! be turned off.
Integer , Parameter              :: PrStatus_Debug = 0 ! debug
Integer , Parameter              :: GC_MaxRead     = 65536 ! debug

! Local variables

Character (Len=*), Parameter     :: RoutineName='Rcf_Grib_Read_Data'

Character (Len=80)               :: Cmessage(2)   ! used for EReport
Integer                          :: ErrorStatus   ! used for EReport

Integer                          :: Error         ! Error Status
Integer                          :: I, Itemp
Integer                          :: IStart,IFinish
Integer                          :: act_io_len,act_io_len2

!=======================================================================
!  Variables Usedfor calls to 'DECODE'
!=======================================================================

Integer , Parameter    :: LenVertCord = 3000   ! Len of Vert_Coords
Integer , Parameter    :: LenBitmap   =    1   ! Len of bitmap used
Integer , Parameter    :: LenQuasi    =    1   ! Len of array Quasi
Integer , Parameter    :: Word_Size   =   64   ! Word size for machine
                                               ! 64 bit on Cray else 32

! LenWork1 and LenWorkR should be >= to no. of columns in grid
! LenWork2 should be >= twice the no. of rows in grid.
Integer , Parameter    :: LenWork1    =  288   ! Len of Work array 1
Integer , Parameter    :: LenWork2    =  500   ! Len of Work array 2
Integer , Parameter    :: LenWorkR    =  288   ! Len of Work array R
Integer , Parameter    :: LenPosn     =    4   ! Len of array Posn

!Variables used when calling DECODE
Integer                :: Width                ! No. of bits used to
                                               ! encode field els
Integer                :: MsgLvl               ! Message level in
                                               ! DECODE

Integer                :: Words                ! no of words
Integer, Parameter     :: Error_Unit =     6   ! Error unit used by
                                               ! DECODE

Integer                :: Bitmap(LenBitmap)    ! Array for bitmap of
                                               ! non full field data
Integer                :: Quasi(LenQuasi)      ! desc of Quasi-reg grid
Integer                :: I_Record(Len_Data)   ! Int array holding
                                               ! Character array of
                                               !encoded GRIB data
Integer                :: Work_Int1(LenWork1)  ! Work Array
Integer                :: Work_Int2(LenWork2)  ! Work Array

Real                   :: Off                  ! Word Offset - not used

Real                   :: FpWork(Len_Data)      ! Work Array
Real                  :: VertCoords(LenVertCord) ! Vertical coord params
Real                   :: Posn(LenPosn)        ! Not Used- reqd 4 decode
Real                   :: Work_Re1(LenWorkR)   ! Work Array
Character(len=1)       :: C_Record(Len_Data) ! Encoded data read into
                                               ! this using Buffin8

!=======================================================================
!  Read binary form of record into CRecord work array
!=======================================================================
! (parameter GC_MaxRead is used and) buffin called more than once
! because the GC routines seem to be hard wired with a limit
! of 65536 bytes

C_Record(:) = " "                           ! Make sure buffer is clear
act_io_len  = 0
act_io_len2 = 0

!Set MsgLvl used by Decode - almost the reverse of rcf system
Select Case (PrintStatus)
Case (PrStatus_Min,PrStatus_Normal)
  MsgLvl = 3          ! No messages
Case (PrStatus_Oper)
  MsgLvl = 1          ! Errors and Warnings
Case (PrStatus_Diag)
  MsgLvl = 0          ! Errors, Warnings and Notes
End Select
! there is also level 2, 'Errors Only' available.

Itemp       = Int ( ( Len_Data - 1 ) / GC_MaxRead )

Do I = 0, Itemp - 1   ! loop to get the maxiumum no of full GC_MaxRead
                      ! long blocks

  IStart  = ( GC_MaxRead * I ) + 1
  IFinish = ( GC_MaxRead * ( I + 1 ) )

  ! Set Position in file to be read
  Call SetPos8( Unit_Num, Pos_in_File + IStart - 1, Error)

  If ( Error /= 0 ) Then
    Write(Cmessage(1),'(A,I2)') 'Failed trying to SetPos to ',        &
                                 Pos_in_File + IStart - 1
    ErrorStatus = 10
    Call EReport( RoutineName, ErrorStatus, Cmessage(1) )
  End If

 ! Read the block of data into CRecord
  Call Buffin8 (Unit_Num, C_Record(IStart:IFinish),                 &
                GC_MaxRead, 1, act_io_len2, Error)

  If ( Error > 0 ) Then
    Cmessage(1)    = 'Failed trying to read record from GRIB file'
    ErrorStatus = 20
    Call EReport( RoutineName, ErrorStatus, Cmessage(1) )
  End If

 ! Add to running total of 'bytes' read
  act_io_len = act_io_len + act_io_len2

End Do

! Repeat the process for the 'remainder' block of data
IStart  = ( GC_MaxRead * Itemp ) + 1
IFinish =  Len_Data

! Set Position in file to be read
Call SetPos8( Unit_Num, Pos_in_File + IStart - 1, Error)

If ( Error /= 0 ) Then
  Write(Cmessage(1),'(A,I2)') 'Failed trying to SetPos to ',        &
                               Pos_in_File + IStart - 1
  ErrorStatus = 30
  Call EReport( RoutineName, ErrorStatus, Cmessage(1) )
End If

! Read the block of data into CRecord
Call Buffin8 (Unit_Num, C_Record(IStart:IFinish),                 &
              IFinish - IStart +1, 1, act_io_len2, Error)

If ( Error > 0 ) Then
  Cmessage(1)    = 'Failed trying to read record from GRIB file'
  ErrorStatus = 40
  Call EReport( RoutineName, ErrorStatus, Cmessage(1) )
End If

! Add to running total of 'bytes' read
act_io_len = act_io_len + act_io_len2

! Report on progress
If ( PrintStatus >= PrStatus_Diag .And. mype == 0 ) Then
    Write (6,*) "Whilst trying to read full record tried for ",       &
                 Len_Data, "and got", act_io_len
End If

!=======================================================================
!  Prepare for and Call DECODE routine
!=======================================================================
! Initialise some of the arrays used in the call to DECODE
Posn(:)               = 0.00
Words                 = 0
Off                   = 0.00
Error                 = 0
I_Record(:)           = 0
Current % Num_Fp      = 0
Current % Num_vert    = 0
Current % Num_Bitmap  = 0
Current % Num_Quasi   = 0

! Despite not being used DECODE expects this variable to be 0
Current % Block_4(p_B4Undef) = 0

! Transfer the character array that the encoded data is in into
! an Integer array which DECODE expects
! After some initial problems with Transfer it was found to be more
! reliable when 'size' was supplied. careful note must be made of the
! relative sizes of the 2 types. On the Cray Int >= Char(len=1) so
! this should not cause problems.

I_Record(:Len_Data) = Transfer(C_Record(:Len_Data),                   &
                                I_Record(1), Len_Data)

Call Decode(Data, FpWork, Len_Data,         &
            Current % Num_Fp,               &
            VertCoords, LenVertCord,        &  !
            Current % Num_vert,             &
            Bitmap, LenBitmap,              &
            Current % Num_Bitmap,           &
            Quasi, LenQuasi,                &
            Current % Num_Quasi,            &
            Width, Word_Size,               &
            Current % Block_0,              &  !
            Current % Block_1,              &  !
            Current % Block_2,              &  !
            Current % Block_3,              &  !
            Current % Block_4,              &
            Current % Block_R,              &  !
            I_Record, act_io_len, Posn, Words, Off,  &
            Error, Work_Int1, Work_Int2, Work_Re1,  &
            Error_Unit, MsgLvl )

If ( Error /= 0 ) Then
  Write (Cmessage(1),'(A,I2,A)') "Error ",Error," returned by DECODE"

  If (Error == 3) Then       ! DECODE is reporting garbage out
    ErrorStatus = 50
    Cmessage(1) = Cmessage(1) // "-Garbage Out-"
    Call EReport( RoutineName, ErrorStatus, Cmessage(1) )

  Else If (Error == 2) Then  ! Decode is complaining it's an extended
                             ! table 2 (e.g. ECMWF)
    If ( PrintStatus >= PrStatus_Diag ) Then
      ErrorStatus = -55        ! DECODE is reporting warnings
      Cmessage(1) = Cmessage(1) // "-Extended Table 2 (ECMWF ?)-"
      Call EReport( RoutineName, ErrorStatus, Cmessage(1) )
    End If
  Else
    ErrorStatus = -50        ! DECODE is reporting warnings
    Call EReport( RoutineName, ErrorStatus, Cmessage(1) )
  End If
End If

! If vertical coordinates exist, store in grib record
If (Current % Num_vert /= 0) Then
  Allocate(Current % VertCoords(Current % Num_vert))
  Current % VertCoords=VertCoords(1:Current % Num_vert)
End If

Return

End Subroutine Rcf_Grib_Read_Data
End Module Rcf_Grib_Read_Data_Mod
#endif
