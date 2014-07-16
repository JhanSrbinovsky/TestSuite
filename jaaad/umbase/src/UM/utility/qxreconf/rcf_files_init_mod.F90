#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initialisation of input and output files

Module Rcf_Files_Init_Mod

!  Subroutine Rcf_Files_Init - initialisation of files
!
! Description:
!   Open files to be opened, setup input header and output header
!   default sizes. Set up Input_Grid datatype and copy input
!   dump if rotation required.
!
! Method:
!   Input and output dumps opened with File_Open. Input header
!   setup and read in. Output header default sizes setup from
!   input sizes.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Contains

Subroutine Rcf_Files_Init( hdr_in, hdr_out )

Use Rcf_HeadAddress_Mod, Only : &
    IC_1stConstRho,         IC_BLevels,            IC_TracerLevs,    &
    IC_XLen,                IC_YLen,               IC_NumLandPoints, &
    IC_PLevels,             IC_WetLevels,          IC_NoCloudLevels, &
    IC_SoilTLevels,         IC_SoilMoistLevs,      IC_NumOzoneLevs,  &
    LDC_EtaRho,             LDC_EtaTheta,          RC_ModelTop,      &
    FH_HorizGrid_Global,    FH_HorizGrid,          IC_ConvectLevs,   &
    RC_PoleLat,             RC_PoleLong,           FH_HorizGrid_NH,  &
    FH_HorizGrid_LamWrap,   FH_HorizGrid_LamWrapEQ,FH_HorizGrid_SH,  &
    IC_HeightMethod,        IC_RiverRows,          IC_RiverRowLength


Use Rcf_Recon_Mod, Only : &
    LEN_INTHD,          LEN_REALHD,        LEN2_LEVDEPC, &
    LEN2_ROWDEPC,       LEN2_COLDEPC,      LEN2_FLDDEPC, &
    LEN_EXTCNST,        LEN_DUMPHIST,      GRIB

Use Rcf_Generate_Heights_Mod, Only : &
    height_gen_original,             &
    height_gen_smooth,               &
    height_gen_ECMWF_Press,          &
    height_gen_ECMWF_Hybrd

Use Rcf_Submodel_Mod, Only : &
    Internal_Model_List, &
    N_Internal_Model

Use Rcf_NRecon_Mod, Only : &
    DumpProgLevs,      &
    PrimDataLen

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Grid_Type_Mod, Only : &
    Input_Grid,        &
    Output_Grid

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_FortranIO_Mod, Only : &
    Rcf_Get_Unit

Use Rcf_ReadUMhdr_Mod, Only : &
    Rcf_ReadUMhdr

Use Rcf_Parvars_mod, Only : &
    mype

Use Rcf_cntlatm_mod, Only : &
    model_domain

Use Rcf_FortranIO_Mod, Only : &
    Max_Filename_Len

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

Implicit None

! Arguments
Type (Um_Header_type), Intent(InOut)    :: hdr_in
Type (Um_Header_type), Intent(InOut)    :: hdr_out

! Comdecks
#include "c_mdi.h"
#include "domtyp.h"

! Local variables
Character (Len=*), Parameter            :: RoutineName='Rcf_Files_Init'
Character (Len=80)                      :: Cmessage
Character (Len=Max_Filename_Len)        :: DumpName
Integer                                 :: ErrorStatus
Integer                                 :: err
Integer                                 :: i

!------------------------------------------------------------------
! First open and setup Input Dump
!------------------------------------------------------------------

Call Rcf_Get_Unit( hdr_in % UnitNum )
If ( Grib ) Then
! DEPENDS ON: file_open
  Call File_Open( hdr_in % UnitNum, 'RECONTMP', 8, 0, 0, err)
Else
! DEPENDS ON: file_open
  Call File_Open( hdr_in % UnitNum, 'AINITIAL', 8, 0, 0, err)
End If

If ( err /= 0 ) Then
  Cmessage    = 'Failed to Open Start Dump'
  ErrorStatus = 10
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

If (PrintStatus >= PrStatus_Normal .and. mype == 0) Then

  Call Fort_Get_Env( 'AINITIAL', 8, DumpName, Max_Filename_Len, err)

  write (6,*)
  If ( Grib ) Then
    write (6,*) 'Input GRIB data : ',DumpName(1:len_trim(DumpName))
  Else
    write (6,*) 'Input dump : ',DumpName(1:len_trim(DumpName))
  End If

End If

Call Rcf_ReadUMhdr( hdr_in )

! Question: Do we need to zero IMDIs in the FLH?

!-------------------------------------------------------------------
! Setup Output Dump header sizes (from Input) and open file
!-------------------------------------------------------------------

Call Rcf_Get_Unit( hdr_out % UnitNum )

! DEPENDS ON: file_open
Call File_Open ( hdr_out % UnitNum, 'ASTART', 6, 1, 0, err)

If ( err /= 0 ) Then
  Cmessage    = 'Failed to Open Output Dump'
  ErrorStatus = 20
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

If (PrintStatus >= PrStatus_Normal .and. mype == 0) Then

  Call Fort_Get_Env( 'ASTART', 6, DumpName, Max_Filename_Len, err)
 
  write (6,*)
  write (6,*) 'Output dump : ',DumpName(1:len_trim(DumpName))


End If

! Initially set some values as in input dump
hdr_out % LenIntC      = hdr_in % LenIntC
hdr_out % LenRealC     = hdr_in % LenRealC
hdr_out % Len2LevDepC  = hdr_in % Len2LevDepC
hdr_out % Len2RowDepC  = hdr_in % Len2RowDepC
hdr_out % Len1Lookup   = hdr_in % Len1Lookup
hdr_out % LenCompFldI1 = hdr_in % LenCompFldI1
hdr_out % LenCompFldI2 = hdr_in % LenCompFldI2
hdr_out % LenCompFldI3 = hdr_in % LenCompFldI3
hdr_out % Len2ColDepC  = hdr_in % Len2ColDepC
hdr_out % Len1FldsOfC  = hdr_in % Len1FldsOfC
hdr_out % Len2FldsOfC  = hdr_in % Len2FldsOfC
hdr_out % LenExtraC    = hdr_in % LenExtraC
hdr_out % LenHistFile  = hdr_in % LenHistFile

! Other values may be overwritten from Namelist RECON
If ( LEN_INTHD /= IMDI ) THEN
  hdr_out % LenIntC = LEN_INTHD
End If
If ( LEN_REALHD /= IMDI ) THEN
  hdr_out % LenRealC = LEN_REALHD
End If
If ( LEN2_LEVDEPC /= IMDI ) THEN
  hdr_out % Len2LevDepC = LEN2_LEVDEPC
End If
If ( LEN2_ROWDEPC /= IMDI ) THEN
  hdr_out % Len2RowDepC = LEN2_ROWDEPC
End If
If ( LEN2_COLDEPC /= IMDI ) THEN
  hdr_out % Len2ColDepC = LEN2_COLDEPC
End If
If ( LEN2_FLDDEPC /= IMDI ) THEN
  hdr_out % Len2FldsOfC = LEN2_FLDDEPC
End If
If ( LEN_EXTCNST /= IMDI ) THEN
  hdr_out % LenExtraC = LEN_EXTCNST
End If
If ( LEN_DUMPHIST /= IMDI ) THEN
  hdr_out % LenHistFile = LEN_DUMPHIST
End If

! Check that sizes of integer constants are big enough for vn5.0
!  dumps
If (hdr_out % LenIntC < 46) Then
  ErrorStatus = 30
  Cmessage = 'Length of Integer Constants needs to be at least&
             & 46 for vn5.0 and higher dumps'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! Other values can be deduced from Output_Grid
hdr_out % Len1ColDepC  = Output_Grid % glob_p_row_length
hdr_out % Len1RowDepC  = Output_Grid % glob_p_rows

hdr_out % Len1LevDepC  = Output_Grid % model_levels + 1

! Others have to be deduced from STASH generated addressing.
hdr_out % Len2Lookup   = 0
hdr_out % LenData      = 0
Do i = 1, N_Internal_Model
  hdr_out % Len2Lookup = hdr_out % Len2Lookup +                    &
                         DumpProgLevs( Internal_Model_List( i ) )
  hdr_out % LenData    = hdr_out % LenData +                       &
                         PrimDataLen( Internal_Model_List( i ) )
End Do

!------------------------------------------------------------------
! Set the input grid data type with info from header
! Assume C grid throughout - Ocean however is B grid (??)
!------------------------------------------------------------------
If ( hdr_in % FixHd( FH_HorizGrid ) == FH_HorizGrid_Global ) Then
  Input_Grid % global          =  .TRUE.
  Input_Grid % Rotated         =  .FALSE.
Else
  Input_Grid % global          = .FALSE.
  If ( hdr_in % RealC( RC_PoleLat ) /= 90 .OR.  &
       hdr_in % RealC( RC_PoleLong) /= 0) Then
    Input_Grid % Rotated = .TRUE.
  Else
    Input_Grid % Rotated = .FALSE.
  End If
End If

Input_Grid % glob_p_row_length = hdr_in % IntC( IC_XLen )
Input_Grid % glob_p_rows       = hdr_in % IntC( IC_YLen )

If ( GRIB ) Then    ! GRIB data on A grid ( v rows not 1 short)
  Input_Grid % glob_u_row_length = hdr_in % IntC( IC_XLen )
  Input_Grid % glob_u_rows       = hdr_in % IntC( IC_YLen )
  Input_Grid % glob_v_row_length = hdr_in % IntC( IC_XLen )
  Input_Grid % glob_v_rows       = hdr_in % IntC( IC_YLen )

Else     ! C grid (global or LAM) Atmos
  Input_Grid % glob_u_row_length = hdr_in % IntC( IC_XLen )
  Input_Grid % glob_u_rows       = hdr_in % IntC( IC_YLen )
  Input_Grid % glob_v_row_length = hdr_in % IntC( IC_XLen )
  If (model_domain == mt_bi_cyclic_lam) THEN
    Input_Grid % glob_v_rows       = hdr_in % IntC( IC_YLen )
  Else
    Input_Grid % glob_v_rows       = hdr_in % IntC( IC_YLen ) - 1
  Endif
End If

Input_Grid % glob_land_field   = hdr_in % IntC( IC_NumLandPoints )

Input_Grid % model_levels      = hdr_in % IntC( IC_PLevels )

Input_Grid % wet_levels        = hdr_in % IntC( IC_WetLevels )
Input_Grid % cloud_levels      = hdr_in % IntC( IC_NoCloudLevels )
Input_Grid % st_levels         = hdr_in % IntC( IC_SoilTLevels )
Input_Grid % sm_levels         = hdr_in % IntC( IC_SoilMoistLevs )
Input_Grid % bl_levels         = hdr_in % IntC( IC_BLevels )
Input_Grid % ozone_levels      = hdr_in % IntC( IC_NumOzoneLevs )
Input_Grid % tr_levels         = hdr_in % IntC( IC_TracerLevs )
Input_Grid % conv_levels       = hdr_in % IntC( IC_ConvectLevs )
Input_Grid % z_top_of_model    = hdr_in % RealC( RC_ModelTop )
Input_Grid % first_constant_r_rho_level = &
                                 hdr_in % IntC( IC_1stConstRho)

! Original height generation method is set if Integer Const is
! either 1 or MDI (ie unset)

If ( hdr_in % IntC( IC_HeightMethod ) == IMDI .OR.                 &
     hdr_in % IntC( IC_HeightMethod ) == height_gen_original ) Then

  Input_Grid % height_gen_method = height_gen_original

Else If ( hdr_in % IntC( IC_HeightMethod ) == height_gen_smooth ) Then
  Input_Grid % height_gen_method = height_gen_smooth

Else If (hdr_in%IntC(IC_HeightMethod) == height_gen_ECMWF_Press ) Then
  Input_Grid % height_gen_method = height_gen_ECMWF_Press

Else If (hdr_in%IntC(IC_HeightMethod) == height_gen_ECMWF_Hybrd ) Then
  Input_Grid % height_gen_method = height_gen_ECMWF_Hybrd

Else
  ErrorStatus = 40
  Cmessage = 'Input dump has unknown height generation method'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! If available set the river routing rows/row_length from header.
! Otherwise set to 0 to allow decomposition to work correctly.
If ( hdr_in % IntC( IC_RiverRows ) /= IMDI ) Then
  Input_Grid % glob_r_rows       = hdr_in % IntC( IC_RiverRows )
  Input_Grid % glob_r_row_length = hdr_in % IntC( IC_RiverRowLength )
Else
  Input_Grid % glob_r_rows       = 0
  Input_Grid % glob_r_row_length = 0
End If

Allocate( Input_grid % eta_theta_levels( 0:Input_Grid % model_levels))
Allocate( Input_grid % eta_rho_levels( Input_Grid % model_levels))
Allocate( Input_grid % rh_crit( Input_Grid % model_levels ) )
Allocate( Input_grid % soil_depths( Input_Grid % sm_levels ) )

Input_grid % eta_theta_levels( 0 : Input_Grid % model_levels ) =  &
             hdr_in % LevDepC( 1 : Input_Grid % model_levels + 1, 1)

Input_grid % eta_rho_levels( 1 : Input_Grid % model_levels ) =    &
           hdr_in % LevDepC( 1 : Input_Grid % model_levels, 2)

Input_Grid % rh_crit( 1 : Input_Grid % model_levels ) =           &
    hdr_in % LevDepC( 1 : Input_Grid % model_levels, 3 )

Input_Grid % soil_depths( 1 : Input_Grid % sm_levels ) =          &
        hdr_in % LevDepC( 1 : Input_Grid % sm_levels, 4 )

Input_Grid % glob_p_field = Input_Grid % glob_p_row_length *        &
                            Input_Grid % glob_p_rows
Input_Grid % glob_u_field = Input_Grid % glob_u_row_length *        &
                            Input_Grid % glob_u_rows
Input_Grid % glob_v_field = Input_Grid % glob_v_row_length *        &
                            Input_Grid % glob_v_rows

!--------------------------------------------------------------------
! If input grid is rotated, the winds will be unrotated and
! written back to file, so will copy the file here and
! make sure the input is not overwritten
!--------------------------------------------------------------------
If ( Input_Grid % Rotated .AND. .NOT. GRIB ) Then
  If (mype == 0 ) Then
    Call Shell('cp $AINITIAL $RECONTMP; chmod +rw $RECONTMP',43)
  End If
! DEPENDS ON: file_close
  Call File_Close( hdr_in % UnitNum, 'AINITIAL', 8, 0, 0, err)

! DEPENDS ON: file_open
  Call File_Open( hdr_in % UnitNum, 'RECONTMP', 8, 1, 0, err)

  If ( err /= 0 ) Then
    Cmessage    = 'Failed to Open Copied Start Dump'
    ErrorStatus = 50
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

End If

Return
End Subroutine Rcf_Files_Init
End Module Rcf_Files_Init_Mod
#endif
