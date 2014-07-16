#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Sets the output dump fixed header.

Module Rcf_Setup_FixHd_Mod

!  Subroutine Rcf_Setup_FixHd - sets the output dump fixed header.
!
! Description:
!   Sets up the fixed header for the output dump based on namelist
!   information.
!
! Method:
!   UMDP F3 defines the fixed header.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Setup_FixHd( Hdr_In, Hdr_Out )

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

Use Rcf_HeadAddress_Mod     ! Huge amounts of this...

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_UMhead_Mod, Only : &
    Um_header_type,    &
    LenFixHd

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Oper

Use Rcf_Recon_Mod, Only : &
    Reset,                &
    LCal360

Use Rcf_readnl_horizont_mod, Only : &
    Delta_Lambda,                   Delta_Phi,                      &
    Phi_First,                      Phi_Npole,                      &
    IProj

Use Rcf_headers_Mod, ONly : &
    FixHd

Implicit None

! Arguments
Type (Um_header_type), Target      :: Hdr_In
Type (Um_header_type), Target      :: Hdr_Out

! Comdecks
#include "c_mdi.h"

! Local parameters:
Integer, PARAMETER :: no_version = -9999 ! Value returned by function
                                         ! get_umversion if environment
                                         ! variable VN not set
! Local Variables
Character (Len=*), Parameter :: RoutineName='Setup_FixHd'
Character (Len=80)           :: Cmessage
Character (Len=8)            :: c_um_version
Integer, Pointer             :: FixHd_In( : )
Integer, Pointer             :: FixHd_Out( : )
Integer                      :: LLLL
Integer                      :: i
Integer                      :: ipos
Integer                      :: um_version
Integer                      :: um_revision
Integer                      :: Horiz_grid_type
Integer                      :: ErrorStatus

 ! Function & Subroutine calls:
Integer            :: get_um_version

External Date_Time

! Setup pointer to fixhd arrays
FixHd_In  => Hdr_In % FixHd
FixHd_Out => Hdr_Out % FixHd

!-----------------------------------------------------------------
! First Setup all parts of fixed header to IMDI
!-----------------------------------------------------------------
FixHd_Out( : ) = IMDI

!-----------------------------------------------------------------
! Data Set Format Version Number
!-----------------------------------------------------------------
FixHd_Out( FH_Version ) = FH_Version_Value
!-----------------------------------------------------------------
! Times
!-----------------------------------------------------------------

! Current Time
Call Date_Time(FixHd_Out( FH_CTYear  ), FixHd_Out( FH_CTMonth ), &
         FixHd_Out( FH_CTDay   ), FixHd_Out( FH_CTHour  ), &
         FixHd_Out( FH_CTMinute), FixHd_Out( FH_CTSecond) )

! Validity Time - copied from input.
FixHd_Out( FH_VTYear : FH_VTDayNo ) = FixHd_In( FH_VTYear : FH_VTDayNo )
FixHd_Out( FH_VTSecond ) = 0

! Initial Time - either copied from input or from validity (if RESET)
If (RESET) Then
  FixHd_Out( FH_DTYear : FH_DTDayNo ) =                       &
                   FixHd_Out( FH_VTYear : FH_VTDayNo )
  If (PrintStatus >= PrStatus_Oper) Then
    Write (6,'('' ANALYSIS TIME RESET '')')
    Write (6,'('' Analysis:'',7I9)') (FixHd_Out(i), i = FH_DTYear, &
                                                  FH_DTDayNo)
  End If
Else
  FixHd_Out( FH_DTYear : FH_DTDayNo ) =                       &
                   FixHd_In( FH_DTYear : FH_DTDayNo )
End If

!-------------------------------------------------------------------
! Other Fixhd stuff - non-addressing
!-------------------------------------------------------------------

! Atmos or Ocean?
FixHd_Out( FH_SubModel ) = FixHd_In( FH_Submodel )

! Vertical Co-ord type
FixHd_Out( FH_VertCoord ) = FH_VertCoord_CP

! Horizontal grid type
If ( Output_Grid % Global ) Then
  delta_lambda=360./Real( Output_Grid % glob_p_row_length )
  delta_phi=180./Real( Output_Grid % glob_p_rows - 1 )
  horiz_grid_type= FH_HorizGrid_Global
Else If( (phi_first > 89.99 ) .AND.                                    &
         (delta_lambda*Output_Grid % glob_p_row_length > 359.99) .AND. &
         (delta_phi * (Output_Grid % glob_p_rows-1) > 89.99) .AND.     &
         (delta_phi * (Output_Grid % glob_p_rows-1) < 90.01)) Then

! Warn user that this option has not been used for several years
! and output should be checked carefully.

  Cmessage = 'N. hemisphere-only grid - check carefully'
  ErrorStatus = -10
  Call Ereport( RoutineName, ErrorStatus, Cmessage )

  delta_lambda=360./Real(Output_Grid % glob_p_row_length)
  delta_phi=90./Real(Output_Grid % glob_p_rows-1)
  horiz_grid_type=FH_HorizGrid_NH

Else If((phi_first < -89.99) .AND.                                     &
        (delta_lambda*Output_Grid % glob_p_row_length > 359.99) .AND.  &
        (delta_phi * (Output_Grid % glob_p_rows-1) > 89.99) .AND.      &
        (delta_phi * (Output_Grid % glob_p_rows-1) < 90.01)) Then


! Warn user that this option has not been used for several years
! and output should be checked carefully.

  Cmessage = 'S. hemisphere-only grid - check carefully'
  ErrorStatus = -20
  Call Ereport( RoutineName, ErrorStatus, Cmessage )

  delta_lambda=360./Real(Output_Grid % glob_p_row_length)
  delta_phi=90./Real(Output_Grid % glob_p_rows-1)
  horiz_grid_type= FH_HorizGrid_SH

Else If(Output_Grid % glob_p_row_length * delta_lambda > 359.99) Then
  If (phi_npole > 89.99) Then
    horiz_grid_type= FH_HorizGrid_LamWrap
  Else
    horiz_grid_type= FH_HorizGrid_LamWrapEq
  End If
Else
  If (phi_npole > 89.99) Then
    horiz_grid_type= FH_HorizGrid_LamNoWrap
  Else
    horiz_grid_type= FH_HorizGrid_LamNoWrapEq
  End If
End If

FixHd_Out( FH_HorizGrid ) = horiz_grid_type


! Instantaneous Data
FixHd_Out( FH_Dataset      ) = FH_Dataset_InstDump
FixHd_Out( FH_RunId        ) = FixHd_In( FH_RunId        )
FixHd_Out( FH_ExptNo       ) = FixHd_In( FH_ExptNo       )

! Calendar - 360 or 365 day
If (LCal360) Then
  FixHd_Out( FH_CalendarType ) = 2     ! 360
Else
  FixHd_Out( FH_CalendarType ) = 1     ! 365
End If

! Grid Staggering - C for atmos
FixHd_Out( FH_GridStagger ) = FH_GridStagger_C

FixHd_Out( FH_AncilDataId ) = FixHd_In( FH_AncilDataId )

! Projection number - need from horizont namelist
FixHd_Out( FH_ProjNo ) = IProj

! Dump version number - take this form VN env var.
! DEPENDS ON: get_um_version
FixHd_Out( FH_ModelVersion ) = get_um_version()


!------------------------------------------------------------------
! Addressing information
!------------------------------------------------------------------
! Integer constants
ipos= LenFixHd+1
FixHd_Out( FH_IntCStart )=ipos
FixHd_Out( FH_IntCSize  )=hdr_out % LenIntC

! Real constants
ipos=ipos+FixHd_Out( FH_IntCSize )
FixHd_Out( FH_RealCStart )=ipos
FixHd_Out( FH_RealCSize  )=hdr_out % LenRealC

! Level dependent constants
ipos=ipos+FixHd_Out( FH_RealCSize )
LLLL=hdr_out % Len1LevDepC*hdr_out % Len2LevDepC
If (LLLL == 0) Then
  FixHd_Out( FH_LevDepCStart )=IMDI
  FixHd_Out( FH_LevDepCSize1 )=IMDI
  FixHd_Out( FH_LevDepCSize2 )=IMDI
Else
  FixHd_Out( FH_LevDepCStart )=ipos
  FixHd_Out( FH_LevDepCSize1 )=hdr_out % Len1LevDepC
  FixHd_Out( FH_LevDepCSize2 )=hdr_out % Len2LevDepC
EndIf

! Row dependent constants
ipos=ipos+LLLL
LLLL=hdr_out % Len2RowDepC*hdr_out % Len1RowDepC
If (LLLL == 0) Then
  FixHd_Out( FH_RowDepCStart )=IMDI
  FixHd_Out( FH_RowDepCSize1 )=IMDI
  FixHd_Out( FH_RowDepCSize2 )=IMDI
Else
  FixHd_Out( FH_RowDepCStart )=ipos
  FixHd_Out( FH_RowDepCSize1 )=hdr_out % Len1RowDepC
  FixHd_Out( FH_RowDepCSize2 )=hdr_out % Len2RowDepC
EndIf

! Column dependent constants
ipos=ipos+LLLL
LLLL=hdr_out % Len1ColDepC*hdr_out % Len2ColDepC
If (LLLL == 0) Then
  FixHd_Out( FH_ColDepCStart )=IMDI
  FixHd_Out( FH_ColDepCSize1 )=IMDI
  FixHd_Out( FH_ColDepCSize2 )=IMDI
Else
  FixHd_Out( FH_ColDepCStart )=ipos
  FixHd_Out( FH_ColDepCSize1 )=hdr_out % Len1ColDepC
  FixHd_Out( FH_ColDepCSize2 )=hdr_out % Len2ColDepC
EndIf

! Fields of constants
ipos=ipos+LLLL
LLLL=hdr_out % Len1FldsOfC*hdr_out % Len2FldsOfC
If (LLLL == 0) Then
  FixHd_Out( FH_FldsOfCStart )=IMDI
  FixHd_Out( FH_FldsOfCSize1 )=IMDI
  FixHd_Out( FH_FldsOfCSize2 )=IMDI
Else
  FixHd_Out( FH_FldsOfCStart )=ipos
  FixHd_Out( FH_FldsOfCSize1 )=hdr_out % Len1FldsOfC
  FixHd_Out( FH_FldsOfCSize2 )=hdr_out % Len2FldsOfC
EndIf

! Extra constants
ipos=ipos+LLLL
LLLL=hdr_out % LenExtraC
If (LLLL == 0) Then
  FixHd_Out( FH_ExtraCStart )=IMDI
  FixHd_Out( FH_ExtraCSize )=IMDI
Else
  FixHd_Out( FH_ExtraCStart )=ipos
  FixHd_Out( FH_ExtraCSize )=hdr_out % LenExtraC
EndIf

! Temp history record
ipos=ipos+LLLL
LLLL=hdr_out % LenHistFile
If (LLLL == 0) Then
  FixHd_Out( FH_HistStart )=IMDI
  FixHd_Out( FH_HistSize )=IMDI
Else
  FixHd_Out( FH_HistStart )=ipos
  FixHd_Out( FH_HistSize )=hdr_out % LenHistFile
EndIf

! Compress index 1
ipos=ipos+LLLL
LLLL=hdr_out % LenCompFldI1
If (LLLL == 0) Then
  FixHd_Out( FH_CompFldI1Start )=IMDI
  FixHd_Out( FH_CompFldI1Size )=IMDI
Else
  FixHd_Out( FH_CompFldI1Start )=ipos
  FixHd_Out( FH_CompFldI1Size )=hdr_out % LenCompFldI1
EndIf
! Compress index 2
ipos=ipos+LLLL
LLLL=hdr_out % LenCompFldI2
If (LLLL == 0) Then
  FixHd_Out( FH_CompFldI2Start )=IMDI
  FixHd_Out( FH_CompFldI2Size )=IMDI
Else
  FixHd_Out( FH_CompFldI2Start )=ipos
  FixHd_Out( FH_CompFldI2Size )=hdr_out % LenCompFldI2
EndIf
! Compress index 3
ipos=ipos+LLLL
LLLL=hdr_out % LenCompFldI3
If (LLLL == 0) Then
  FixHd_Out( FH_CompFldI3Start )=IMDI
  FixHd_Out( FH_CompFldI3Size )=IMDI
Else
  FixHd_Out( FH_CompFldI3Start )=ipos
  FixHd_Out( FH_CompFldI3Size )=hdr_out % LenCompFldI3
EndIf

! Lookup
ipos=ipos+LLLL
LLLL=hdr_out % Len1Lookup*hdr_out % Len2Lookup
If (LLLL == 0) Then
  FixHd_Out( FH_LookupStart )=IMDI
  FixHd_Out( FH_LookupSize1 )=IMDI
  FixHd_Out( FH_LookupSize2 )=IMDI
Else
  FixHd_Out( FH_LookupStart )=ipos
  FixHd_Out( FH_LookupSize1 )=hdr_out % Len1Lookup
  FixHd_Out( FH_LookupSize2 )=hdr_out % Len2Lookup
!      For dumps, add no of prognostic fields to FixHd_Out
!      In Reconfiguration hdr_out % Len2Lookup = No of prognostic fields
 If (FixHd_Out( FH_DataSet ) == FH_DataSet_InstDump ) Then
   FixHd_Out( FH_NumProgFields ) = hdr_out % Len2Lookup
 EndIf
EndIf

! Model data
ipos=ipos+LLLL
LLLL=hdr_out % LenData
If (LLLL == 0) Then
  FixHd_Out( FH_DataStart )=IMDI
  FixHd_Out( FH_DataSize )=IMDI
Else
          ! make sure the data starts on a sector bndry
  fixhd_Out( FH_DataStart )=    &
       ((ipos+um_sector_size-1)/um_sector_size)*um_sector_size+1
  FixHd_Out( FH_DataSize )=hdr_out % LenData
EndIf

!-----------------------------------------------------------------
! Finally, do we need to overwrite with namelist values?
!-----------------------------------------------------------------
Do i = 1, LenFixHd
  If ( FixHd(i) /= IMDI ) Then
    If ( PrintStatus >= PrStatus_Oper ) Then
      Write (6,*) 'FixHd(',i,') has been reset from ', FixHd_out(i), &
                  ' to ', FixHd(i)
    End If

    FixHd_Out(i) = FixHd(i)
  End If
End Do

Return
End Subroutine Rcf_Setup_FixHd
End Module Rcf_Setup_FixHd_Mod
#endif
