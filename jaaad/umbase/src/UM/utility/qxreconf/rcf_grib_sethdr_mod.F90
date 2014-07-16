#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Set up a UM header describing the GRIB data

Module Rcf_Grib_SetHdr_Mod

! SUBROUTINE Rcf_Grib_SetHdr
!
! Description: Set up the UM style header using the GRIB header
!              information.
!
! Method: Copies information from the GRIB headers into the various
!         derived types used for UM headers and then calls the normal
!         header set up routines.
!
! Current Code Owner: Roddy Sharp
!
! History:
! Version   Date      Comment
! -------  --------   -------------------------
!  5.4     12/06/02   Original code. Roddy Sharp (frtz)
!  6.0     04/08/03   Add fix for reading model level data from GRIB
!                     files. Paul Earnshaw (frpe)
!  6.2     09/03/06   Initialise FixHd in Hdr_Dmy to IMDI.
!                     Prevents downstream errors with dump  R.Sharp
!  6.2     20/06/05   Modified model level code to remove ambiguity.
!                     Paul Earnshaw
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
contains
Subroutine Rcf_Grib_SetHdr(Lists,Output_Grid,Hdr_Dmy,Hdr_Itm)

Use Rcf_Grib_Block_Params_Mod
! Has derived types of List_Marker and Grib_record
! plus also containd most parameter definitions starting with p_

Use Rcf_Grib_Lookups_Mod, Only : &
    grib_max_fields,             &
    grib_U_field,                &
    grib_Soil_Temp_field,        &
    grib_Soil_Moist_field

Use Rcf_Grib_Block_Params_Mod, Only : &
    Tb3_Pressure,                  &
    Tb3_Hybrid

Use Rcf_Grid_Type_Mod, Only : &
    Grid_type

Use Rcf_UMhead_Mod, Only :  &
    LenFixHd,               &
    um_header_type            ! Derived containing UM header info

Use Rcf_HeadAddress_Mod

Use Rcf_headers_Mod, Only : &
    FixHd

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Min,           &       ! =1 Minimum output
    PrStatus_Normal,        &       ! =2 Short informative output
    PrStatus_Oper,          &       ! =3 Full informative output
    PrStatus_Diag                   ! =4 Extra Diagnostic output

Use Rcf_Parvars_mod, Only : &
    mype

Use Rcf_Submodel_Mod, Only :    &
    Atmos_IM

Use EReport_Mod, Only :     &
    EReport

Use Rcf_AllocHdr_Mod, Only :     &
    Rcf_AllocHdr

Use Rcf_Setup_FixHd_Mod, Only :     &
    Rcf_Setup_FixHd

Use Rcf_Generate_Heights_Mod, Only :     &
    height_gen_ecmwf_press,         &
    height_gen_ecmwf_hybrd

Implicit None

! Subroutine arguments

!< Array  arguments with intent(in):>
Type (List_Marker)               :: Lists(0:grib_max_fields)

!< Scalar arguments with intent(InOut):>
Type (Um_Header_type)                    :: Hdr_Dmy
Type (Um_Header_type)                    :: Hdr_Itm
Type (Grid_type)                         :: Output_Grid

! Comdecks
#include "c_mdi.h"
#include "clookadd.h"

! Local variables

Type (Grib_Record), Pointer      :: Current
Type (Grib_Record), Pointer      :: First_Multi  ! first multi level
                                                 ! list for hdr info

Integer                          :: NameLst_Temp(LenFixHd)

Character (Len=*), Parameter     :: RoutineName='Rcf_GRIB_SetHdr'
Character (Len=80)               :: Cmessage(2)   ! used for EReport
Integer                          :: ErrorStatus   ! used for EReport
Integer                          :: I,Count
Integer                          :: Sign

Integer                          :: Fields_tot
Integer                          :: Int_Val  ! Used for mold in transfer
Integer                          :: First_Multi_List
Integer                          :: First_Multi_Type

!=======================================================================
! Preparation for setting up the fixed header.
!=======================================================================

Hdr_Dmy % FixHd (:)     = IMDI

! Count the number of records in all the lists _except_ the misc list
! Also, set pointer to the first entry in the first pressure field
! list encountered.

Fields_tot       = 0
First_Multi_List = 0
Nullify (First_Multi)

Do I = 1, grib_max_fields
  If (Associated(Lists(I) % Begin) ) Then
    Fields_tot = Fields_tot + Lists(I) % LstCount

    If ( (.Not. Associated(First_Multi)) .And.  (                     &
       (Lists(I) % Begin % Block_1(p_Lvl_Type) == Tb3_Pressure)  .OR. &
       (Lists(I) % Begin % Block_1(p_Lvl_Type) == Tb3_Hybrid  )       &
                                                              )) Then
      First_Multi      => Lists(I) % Begin
      First_Multi_List = I
      First_Multi_Type = Lists(I) % Begin % Block_1(p_Lvl_Type)
    End If
  End If
End Do

! Double check that 'First_Multi' was assigned
If (.Not. Associated(First_Multi)) Then
  CMessage(1) = "No multi-level fields found in GRIB data"
  ErrorStatus = 10
  Call EReport( RoutineName, ErrorStatus, Cmessage(1))
End If

! Set Values that would normally have been taken from the start dump
Output_Grid % glob_p_row_length = First_Multi % Block_2(p_Pnts_Prll)
Output_Grid % glob_p_rows       = First_Multi % Block_2(p_Pnts_Merid)
Output_Grid % model_levels      = Lists(First_Multi_List) % LstCount

! Check for Soil Temp Levels before using their count.
If (Associated(Lists(grib_Soil_Temp_field) % Begin)) Then
  Output_Grid % st_levels     = Lists(grib_Soil_Temp_field) % LstCount
End If

! Check for Soil Moisture Levels before using their count.
If (Associated(Lists(grib_Soil_Moist_field) % Begin)) Then
  Output_Grid % sm_levels     = Lists(grib_Soil_Moist_field) % LstCount
End If

Hdr_Dmy % FixHd (FH_RunId)      = 0
Hdr_Dmy % FixHd (FH_ExptNo)     = IMDI
Hdr_Dmy % FixHd (FH_Submodel)   = Atmos_IM

Hdr_Itm % LenIntC               = p_Len_IntC
Hdr_Itm % LenRealC              = p_Len_RealC
Hdr_Itm % Len1LevDepC           = Lists(grib_U_field) % LstCount + 1
Hdr_Itm % Len2LevDepC           = p_Len2_LevDepC
Hdr_Itm % Len1RowDepC           = p_Len1_RowDepC
Hdr_Itm % Len2RowDepC           = p_Len2_RowDepC
Hdr_Itm % Len1ColDepC           = p_Len1_ColDepC
Hdr_Itm % Len2ColDepC           = p_Len2_ColDepC
Hdr_Itm % Len1FldsofC           = p_Len1_FldsofC
Hdr_Itm % Len2FldsofC           = p_Len2_FldsofC
Hdr_Itm % LenExtraC             = p_Len_ExtraC
Hdr_Itm % LenHistFile           = p_Len_HistFile
Hdr_Itm % LenCompFldI1          = p_Len_CompFldI1
Hdr_Itm % LenCompFldI2          = p_Len_CompFldI2
Hdr_Itm % LenCompFldI3          = p_Len_CompFldI3
Hdr_Itm % Len1Lookup            = p_Len1_Lookup
Hdr_Itm % Len2Lookup            = Fields_tot
Hdr_Itm % LenData               = Fields_tot *                        &
                                 First_Multi % Block_2(p_Pnts_Prll) * &
                                 First_Multi % Block_2(p_Pnts_Merid)

! Set up the Initial data Time and Validity Time
Hdr_Dmy % FixHd (FH_DTYear)   = First_Multi % Block_1 (p_Ref_Year)    &
                     + ((First_Multi % Block_1 (p_Ref_Cent) -1) * 100)
Hdr_Dmy % FixHd (FH_DTMonth:FH_DTMinute) =                            &
                        First_Multi % Block_1 (p_Ref_Month:p_Ref_Min)
Hdr_Dmy % FixHd (FH_DTSecond) = 0
Hdr_Dmy % FixHd (FH_DTDayNo)  = 0

Hdr_Dmy % FixHd (FH_VTYear)   = First_Multi % Block_1 (p_Ref_Year)    &
                     + ((First_Multi % Block_1 (p_Ref_Cent) -1) * 100)
Hdr_Dmy % FixHd (FH_VTMonth:FH_VTMinute) =                            &
                        First_Multi % Block_1 (p_Ref_Month:p_Ref_Min)
Hdr_Dmy % FixHd (FH_VTSecond) = 0
Hdr_Dmy % FixHd (FH_VTDayNo)  = 0

If (PrintStatus >= PrStatus_Diag) Then
  If ( mype == 0 ) Then
    Write (6,'(A,I4,2("/",I2),/,A,2(I2,":"),I2,/,A,I3)')            &
     "Date from GRIB data appears to be (yyyy/mm/dd) : ",           &
      Hdr_Dmy % FixHd (FH_DTYear:FH_DTDay),                         &
     "And the Time is : ",Hdr_Dmy % FixHd (FH_DTHour:FH_DTSecond),  &
     "Dayno. is :",Hdr_Dmy % FixHd (FH_DTDayNo)
  End If
End If

! Set the 'global' flag if appropriate
! Check that the grid represents 360 by 180 degrees
If ( ( (First_Multi % Block_2(p_Pnts_Prll) *                          &
        First_Multi % Block_2(p_IncrPrll )  ) >= 359999 ) .AND.       &
     ( (First_Multi % Block_2(p_Pnts_Merid) *                         &
        First_Multi % Block_2(p_IncrMerid)   ) >= 179999 ) )  Then

  Output_Grid % Global = .True.

Else
  ! Areas of the GRIB handling code _assume_ the incoming grid is a
  ! global grid and thus an abort is placed here
  CMessage(1) = "Does not appear to be Global grid"
  ErrorStatus = 20
  Call EReport( RoutineName, ErrorStatus, Cmessage(1))

End If  ! test to see if using global grid

! Call to Set up fixed header.

! First - store Recon Namelist and send dummy blank one to prevent
!         overwriting of data
NameLst_Temp(:) = FixHd(:)
FixHd(:) = IMDI

Call Rcf_Setup_FixHd(Hdr_Dmy, Hdr_Itm)

! Reset values describing the grid which are 'hardwired' in above call
Hdr_Itm % FixHd ( FH_VertCoord )   = FH_VertCoord_Pressure
Hdr_Itm % FixHd ( FH_GridStagger ) = FH_GridStagger_A

! Not forgetting to restore the Recon Namelist
FixHd(:) = NameLst_Temp(:)

! Call to Allocate memory for rest of header info
Call Rcf_AllocHdr (Hdr_Itm)

!=======================================================================
! Initialising Integer Constants
!=======================================================================
! Emulating what was done at 4.5

Hdr_Itm % IntC (:) = IMDI

Hdr_Itm % IntC (IC_XLen)          = Output_Grid % glob_p_row_length
Hdr_Itm % IntC (IC_YLen)          = Output_Grid % glob_p_rows
Hdr_Itm % IntC (IC_PLevels)       = Output_Grid % model_levels
Hdr_Itm % IntC (IC_WetLevels)     = Output_Grid % model_levels
Hdr_Itm % IntC (IC_SoilTLevels)   = Output_Grid % st_levels
Hdr_Itm % IntC (IC_SoilMoistLevs) = Output_Grid % sm_levels

! Why set number of boundary layer levels like this?
! Hardwired to pressure levels!
Hdr_Itm % IntC (IC_BLevels)     = 0
Current => Lists(First_Multi_List) % Begin  ! just to keep following
Do While (Associated(Current))
  If (Current % Block_1 (p_Lvl_Desc_1) > 850 ) Then
    Hdr_Itm % IntC (IC_BLevels) = Hdr_Itm % IntC (IC_BLevels) +1
  End If
  Current => Current % Next
End Do

Hdr_Itm % IntC (IC_MDI)         = IMDI

! Set height_generation_method for Pressure or Hybrid
Select Case (First_Multi_Type)

  Case (Tb3_Pressure)
    Hdr_Itm % IntC (IC_HeightMethod)  = height_gen_ecmwf_press

  Case (Tb3_Hybrid  )
    Hdr_Itm % IntC (IC_HeightMethod)  = height_gen_ecmwf_hybrd

  Case Default
    Cmessage(1)    =                                                  &
         'Multi-level field type not set for height generation method'
    ErrorStatus = 55
    Call EReport( RoutineName, ErrorStatus, Cmessage(1) )

End Select

!=======================================================================
! Initialising Real Constants
!=======================================================================

Hdr_Itm % RealC (:) = RMDI

! First latitude and longitude points in Real Consts Header
Hdr_Itm % RealC (RC_FirstLat)    = First_Multi % Block_2(p_LatGrdPnt1) &
                                   / 1000.000
Hdr_Itm % RealC (RC_FirstLong)   = First_Multi % Block_2(p_LonGrdPnt1) &
                                   / 1000.000

! Latitude and longitude spacing values for Real Const (normally +ve)
Hdr_Itm % RealC (RC_LatSpacing) = First_Multi % Block_2(p_IncrMerid)  &
                                   / 1000.000
Hdr_Itm % RealC (RC_LongSpacing) = First_Multi % Block_2(p_IncrPrll)  &
                                   / 1000.000

! GRIB data holds the Lat and Long of the Southern Pole _iff_ the grid
! is rotated - otherwise both values are zero.
If ( (First_Multi % Block_2(p_S_Pole_Lat) == 0.00 ) .And.             &
     (First_Multi % Block_2(p_S_Pole_Lon) == 0.00 ) ) Then
  Hdr_Itm % RealC (RC_PoleLat)     = 90.0
  Hdr_Itm % RealC (RC_PoleLong)    = 00.0

  If (PrintStatus >= PrStatus_Diag) Then
    If (mype == 0 ) Then
      Write (6,*) "Using default value for North Pole lat and Long"
    End If
  End If
Else
  If (PrintStatus >= PrStatus_Diag) Then
    If (mype == 0 ) Then
      Write (6,*) "Using calculated North Pole lat and Long"
    End If
  End If

  Hdr_Itm % RealC (RC_PoleLat) =                                      &
                  90.0 - (First_Multi % Block_2(p_S_Pole_Lat) / 1000 )
  Hdr_Itm % RealC (RC_PoleLong) =                                     &
                  (First_Multi % Block_2(p_S_Pole_Lon) / 1000 ) + 180.00
  If ( Hdr_Itm % RealC (RC_PoleLong)>= 360.00 ) Then
    Hdr_Itm % RealC (RC_PoleLong)= Hdr_Itm % RealC (RC_PoleLong)- 360.00
  End If

End If

!=======================================================================
!  Initialising Level Dependant Constants
!=======================================================================

Hdr_Itm % LevDepC (:,:) = RMDI
Current => Lists(First_Multi_List) % Begin
Select Case (First_Multi_Type)

  Case (Tb3_Pressure)
    ! *100 to go from HPa to Pa
!cdir novector
    Do I = 1 , Hdr_Itm % Len1LevDepC
      Hdr_Itm % LevDepC (I,LDC_Pressure) =                          &
                          Current % Block_1 ( p_Lvl_Desc_1 ) * 100
      Current => Current % Next
    End Do

  Case (Tb3_Hybrid  )
    ! Model levels keep model level index from ECMWF GRIB file
!cdir novector
    Do I = 1 , Hdr_Itm % Len1LevDepC
      Hdr_Itm % LevDepC (I,LDC_MLIndex) =                           &
                          Current % Block_1 ( p_Lvl_Desc_1 )
      Current => Current % Next
    End Do

  Case Default
    Cmessage(1)    =                                                &
           'Multi-level field type not set to correct value'
    ErrorStatus = 56
    Call EReport( RoutineName, ErrorStatus, Cmessage(1) )

End Select

! Set soil level dependent constants
! Note: rcf_grib_check sets soil index in p_Lvl_Desc_1 and 
!       soil depths in p_Lvl_Desc_2
If (Associated(Lists(grib_Soil_Temp_field) % Begin)) Then
  Current => Lists(grib_Soil_Temp_field) % Begin
!cdir novector
  Do I = 1 , Hdr_Itm % IntC (IC_SoilTLevels)
    ! Level depth stored in cm, convert to m
    Hdr_Itm % LevDepC (I,SoilDepths) =                              &
                        Current % Block_1 ( p_Lvl_Desc_2 ) / 100.
    Current => Current % Next
  End Do
Else If (Associated(Lists(grib_Soil_Moist_field) % Begin)) Then
  Current => Lists(grib_Soil_Moist_field) % Begin
!cdir novector
  Do I = 1 , Hdr_Itm % IntC (IC_SoilMoistLevs)
    ! Level depth stored in cm, convert to m
    Hdr_Itm % LevDepC (I,SoilDepths) =                              &
                        Current % Block_1 ( p_Lvl_Desc_2 ) / 100.
    Current => Current % Next
  End Do
End If

!=======================================================================
!  Initialising Row Dependant Constants
!=======================================================================

Hdr_Itm % RowDepC( :,: ) = RMDI

!=======================================================================
!  Initialising Col Dependant Constants
!=======================================================================

Hdr_Itm % ColDepC( :,: ) = RMDI

!=======================================================================
!  Initialising Addresses and Lengths in Lookup
!=======================================================================

!Clear the decks
Hdr_Itm % Lookup(  1 : 45, : ) = 0
Hdr_Itm % Lookup( 46 : 64, : ) = Transfer( 0.0, int_val)

!Set the data time - twice
Hdr_Itm % Lookup( LBYR   , : ) = Hdr_Dmy % FixHd (FH_DTYear)
Hdr_Itm % Lookup( LBMON  , : ) = Hdr_Dmy % FixHd (FH_DTMonth)
Hdr_Itm % Lookup( LBDAT  , : ) = Hdr_Dmy % FixHd (FH_DTDay)
Hdr_Itm % Lookup( LBHR   , : ) = Hdr_Dmy % FixHd (FH_DTHour)
Hdr_Itm % Lookup( LBMIN  , : ) = Hdr_Dmy % FixHd (FH_DTMinute)
Hdr_Itm % Lookup( LBYRD  , : ) = Hdr_Dmy % FixHd (FH_DTYear)
Hdr_Itm % Lookup( LBMOND , : ) = Hdr_Dmy % FixHd (FH_DTMonth)
Hdr_Itm % Lookup( LBDATD , : ) = Hdr_Dmy % FixHd (FH_DTDay)
Hdr_Itm % Lookup( LBHRD  , : ) = Hdr_Dmy % FixHd (FH_DTHour)
Hdr_Itm % Lookup( LBMIND , : ) = Hdr_Dmy % FixHd (FH_DTminute)

!Set Field dimensions
Hdr_Itm % Lookup( LBNPT  , : ) = Output_Grid % glob_p_row_length
Hdr_Itm % Lookup( LBROW  , : ) = Output_Grid % glob_p_rows
Hdr_Itm % Lookup( LBLREC , : ) = Output_Grid % glob_p_rows *      &
                                   Output_Grid % glob_p_row_length

!Set Data type to 'real' accross the board.
Hdr_Itm % Lookup( DATA_TYPE , : ) = 1

!Set Packing Type to 'no packing' accross the board.
Hdr_Itm % Lookup( LBPACK    , : ) = 0

!Set Header Release Number
Hdr_Itm % Lookup( LBREL     , : ) = 2

!Set Internal Model No.
Hdr_Itm % Lookup( MODEL_CODE , : ) = 1

!Words 46 to 64 are reals stored as integers

! Real Lat and Long of Pseudo North Pole
Hdr_Itm % Lookup( BPLAT , : ) = Transfer (                          &
                             Hdr_Itm % RealC (RC_PoleLat) , int_val )

Hdr_Itm % Lookup( BPLON , : ) = Transfer (                          &
                             Hdr_Itm % RealC (RC_PoleLong), int_val )

! Lat and Long spacing
! calculate sign of latitude spacing (from 'last - first' grid point)
Sign = ( First_Multi % Block_2(p_LatExtrmPt) -                        &
                       First_Multi % Block_2(p_LatGrdPnt1) )
Sign = Sign / Abs(Sign)

Hdr_Itm % Lookup( BDY   , : ) = Transfer (                          &
                    Sign * Hdr_Itm % RealC (RC_LatSpacing), int_val )

Hdr_Itm % Lookup( BDX   , : ) = Transfer (                          &
                          Hdr_Itm % RealC (RC_LongSpacing), int_val )

! Zeroth Lat and Long - (Zeroth not First, Hence subtraction)
Hdr_Itm % Lookup( BZY   , : ) = Transfer(                           &
    Hdr_Itm % RealC (RC_FirstLat) -                                 &
                         ( Sign * Hdr_Itm % RealC (RC_LatSpacing) ) &
                                            , int_val )

Hdr_Itm % Lookup( BZX   , : ) = Transfer (                          &
  Hdr_Itm % RealC (RC_FirstLong) - Hdr_Itm % RealC (RC_LongSpacing) &
                                          , int_val )

! Set the Lookup which holds the stash code

count = 0

!Loop across all lists
Do I = 1, grib_max_fields

  If (Associated(Lists(I) % Begin) ) Then

    Current => Lists(I) % Begin
    Do While (Associated(Current))

      count = count + 1

      Hdr_Itm % Lookup( ITEM_CODE , count ) = Current % StashCode
      Hdr_Itm % Lookup( BLEV      , count ) =                       &
           Transfer (Real(Current % Block_1 (p_Lvl_Desc_1)),int_val)

      Current => Current % Next

    End Do  ! members of a list

  End If ! Associated(Lists(I) % Begin)

End Do ! I over all lists


Return

End Subroutine Rcf_Grib_SetHdr
End Module Rcf_Grib_SetHdr_Mod
#endif
