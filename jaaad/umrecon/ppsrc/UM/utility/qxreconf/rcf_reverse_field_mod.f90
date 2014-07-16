
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ reverse a N-S grid to give a S-N one

Module Rcf_Reverse_Field_Mod

! Description:
!              Reverse the position of "n_rows" rows for each level of
! the 3D data array passed into the routine such that the first row
! swaps place with the row at position n_rows. The 2nd row swaps place
! with th row at position n_rows-1 and so on. The lookup headers for
! the field are altered to be consistent with the new row ordering.
! Based on the routine - PF_REVERSE at UM4.5
!
! Method:
!
! Current Code Owner: Roddy Sharp
!
! History:
! Version   Date      Comment
! -------  --------   -------------------------
!  5.4     10/07/02   Original code. Roddy Sharp (frtz)
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! -Note- This routine has been overloaded to handle both real and
! logical fields. Please ensure you duplicate 'fixes' performed in one
! routine in the other routine where necessary

Interface Rcf_Reverse_Field
  Module Procedure Rcf_Reverse_Field_Log, Rcf_Reverse_Field_Real
End Interface

Contains

Subroutine Rcf_Reverse_Field_Real(RData,row_length,n_rows,levels,     &
                                  field_no,UM_Hdr)

Use Rcf_UMhead_Mod, Only :  &
    UM_Header_Type            ! Derived containing UM header info

Use Rcf_HeadAddress_Mod, Only  : &
    RC_FirstLat,       &
    RC_LatSpacing

Use EReport_Mod, Only :     &
    EReport

Implicit None
! Subroutine arguments

!< Scalar arguments with intent(in):>
Integer, Intent(In)              :: row_length  ! length of each row
Integer, Intent(In)              :: n_rows      ! no. of rows per level
Integer, Intent(In)              :: levels      ! no. of levels
Integer, Intent(In)              :: field_no    ! pos in lookup of
                                                ! first field

!< Array  arguments with intent(InOut):>
Real,    Intent(InOut)              :: RData(*) ! the data to be reversd
Type (UM_Header_Type),Intent(InOut) :: UM_Hdr   ! UM header info

! Comdecks
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
! defines BDY and BZY

! Local variables

Integer                          :: Int_Val
Integer                          :: top_row
Integer                          :: bottom_row
Integer                          :: point
Integer                          :: lev_shift
Integer                          :: b_row_shift
Integer                          :: t_row_shift
Integer                          :: point_loc
Integer                          :: k
Real                             :: temp
Real                             :: Zeroth_Lat
Real                             :: Lat_Spacing
Real                             :: Real_Val

Integer                          :: ErrorStatus
Character(len=80)                :: cMessage
Character (Len=*), Parameter     :: RoutineName='Rcf_Reverse_Field_Real'

!=======================================================================
! Subroutine Code.
!=======================================================================
!loop over levels
Do k=1,levels

  lev_shift = (k-1) * (row_length * n_rows)

  Do top_row=1,n_rows/2

    bottom_row = n_rows + 1 - top_row

    b_row_shift = lev_shift + ((bottom_row - 1) * row_length)
    t_row_shift = lev_shift + ((top_row    - 1) * row_length)

    Do point=1,row_length

      ! Swap over corresponding point at top and bottom of data.
      temp                      = RData(b_row_shift + point)
      RData(b_row_shift + point) = RData(t_row_shift + point)
      RData(t_row_shift + point) = temp

    End Do ! point

  End Do ! top_row

  ! Now update the Latitude info held in Lookups
  Lat_Spacing = Transfer( UM_Hdr % Lookup( BDY, field_no ) , Real_Val)
  Zeroth_Lat  = Transfer( UM_Hdr % Lookup( BZY, field_no ) , Real_Val)

  ! Transpose Zeroth latitude.
  Zeroth_Lat  = Zeroth_Lat + ( (n_rows + 1) * Lat_Spacing )

  ! Sign of the latitude interval is changed.
  Lat_Spacing = - Lat_Spacing


  ! Put new values back into the Lookups
  UM_Hdr % Lookup( BDY, field_no ) = Transfer(Lat_Spacing , Int_Val)
  UM_Hdr % Lookup( BZY, field_no ) = Transfer(Zeroth_Lat  , Int_Val)

  ! double check the values in the Real Headers
  If (Um_Hdr % RealC (RC_FirstLat) /= ( Zeroth_Lat + Lat_Spacing)) Then
    Um_Hdr % RealC (RC_FirstLat) = ( Zeroth_Lat + Lat_Spacing)
  End If
  If (Um_Hdr % RealC (RC_LatSpacing) /= Abs(Lat_Spacing) ) Then
    !Um_Hdr % RealC (RC_LatSpacing) = Abs(Lat_Spacing)
    Write(cMessage,'(A)') "Latitude Spacing appears to have changed"
    ErrorStatus = 20
    Call EReport( RoutineName, ErrorStatus, Cmessage)

  End If


  ! ***************************************************
  ! At UM4.5 something was done here to re-align u or v
  ! onto the correct spacings.
  ! ***************************************************

End Do ! k over levels


Return

End Subroutine Rcf_Reverse_Field_Real


Subroutine Rcf_Reverse_Field_Log(LData,row_length,n_rows,levels,      &
                                  field_no,UM_Hdr)

Use Rcf_UMhead_Mod, Only :  &
    UM_Header_Type            ! Derived containing UM header info

Use Rcf_HeadAddress_Mod, Only  : &
    RC_FirstLat,       &
    RC_LatSpacing

Use EReport_Mod, Only :     &
    EReport

Implicit None
! Subroutine arguments

!< Scalar arguments with intent(in):>
Integer, Intent(In)              :: row_length  ! length of each row
Integer, Intent(In)              :: n_rows      ! no. of rows per level
Integer, Intent(In)              :: levels      ! no. of levels
Integer, Intent(In)              :: field_no    ! pos in lookup of
                                                ! first field

!< Array  arguments with intent(InOut):>
Logical, Intent(InOut)              :: LData(*) ! the data to be reversd
Type (UM_Header_Type),Intent(InOut) :: UM_Hdr   ! UM header info

! Comdecks
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
! defines BDY and BZY

! Local variables

Integer                          :: Int_Val
Integer                          :: top_row
Integer                          :: bottom_row
Integer                          :: point
Integer                          :: lev_shift
Integer                          :: b_row_shift
Integer                          :: t_row_shift
Integer                          :: point_loc
Integer                          :: k
Real                             :: Zeroth_Lat
Real                             :: Lat_Spacing
Real                             :: Real_Val
Logical                          :: temp

Integer                          :: ErrorStatus
Character(len=80)                :: cMessage
Character (Len=*), Parameter     :: RoutineName='Rcf_Reverse_Field_Log'

!=======================================================================
! Subroutine Code.
!=======================================================================
!loop over levels
Do k=1,levels

  lev_shift = (k-1) * (row_length * n_rows)

  Do top_row=1,n_rows/2

    bottom_row = n_rows + 1 - top_row

    b_row_shift = lev_shift + ((bottom_row - 1) * row_length)
    t_row_shift = lev_shift + ((top_row    - 1) * row_length)

    Do point=1,row_length

      ! Swap over corresponding point at top and bottom of data.
      temp                      = LData(b_row_shift + point)
      LData(b_row_shift + point) = LData(t_row_shift + point)
      LData(t_row_shift + point) = temp

    End Do ! point

  End Do ! top_row

  ! Now update the Latitude info held in Lookups
  Lat_Spacing = Transfer( UM_Hdr % Lookup( BDY, field_no ) , Real_Val)
  Zeroth_Lat  = Transfer( UM_Hdr % Lookup( BZY, field_no ) , Real_Val)

  ! Transpose Zeroth latitude.
  Zeroth_Lat  = Zeroth_Lat + ( (n_rows + 1) * Lat_Spacing )

  ! Sign of the latitude interval is changed.
  Lat_Spacing = - Lat_Spacing


  ! Put new values back into the Lookups
  UM_Hdr % Lookup( BDY, field_no ) = Transfer(Lat_Spacing , Int_Val)
  UM_Hdr % Lookup( BZY, field_no ) = Transfer(Zeroth_Lat  , Int_Val)

  ! double check the values in the Real Headers
  If (Um_Hdr % RealC (RC_FirstLat) /= ( Zeroth_Lat + Lat_Spacing)) Then
    Um_Hdr % RealC (RC_FirstLat) = ( Zeroth_Lat + Lat_Spacing)
  End If
  If (Um_Hdr % RealC (RC_LatSpacing) /= Abs(Lat_Spacing) ) Then
    !Um_Hdr % RealC (RC_LatSpacing) = Abs(Lat_Spacing)
    Write(cMessage,'(A)') "Latitude Spacing appears to have changed"
    ErrorStatus = 20
    Call EReport( RoutineName, ErrorStatus, Cmessage)

  End If


  ! ***************************************************
  ! At UM4.5 something was done here to re-align u or v
  ! onto the correct spacings.
  ! ***************************************************

End Do ! k over levels

Return

End Subroutine Rcf_Reverse_Field_Log

End Module Rcf_Reverse_Field_Mod
