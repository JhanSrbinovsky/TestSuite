#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Read the HORIZONT namelist

Module rcf_readnl_horizont_mod

!  Subroutine Rcf_Readnl_Horizont - Read the HORIZONT namelist
!
! Description:
!   Reads the HORIZONT namelist controlling horizontal interpolation
!   and domains.
!
! Method:
!   Data  read and Output_Grid set accordingly.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

 
! Data for namelist reads - mostly for LAM
Real, Save            :: Delta_Lambda
Real, Save            :: Delta_Phi
Real, Save            :: Lambda_First
Real, Save            :: Phi_First
Real, Save            :: Lambda_NPole
Real, Save            :: Phi_NPole
Integer, Save         :: IProj
Integer, Save         :: orog_blend_width      !} For orography blending
Real, Pointer, Save   :: blend_weights(:)      !} zone
Integer, Save         :: Extended_Halo_Size_EW
Integer, Save         :: Extended_Halo_Size_NS
 
! Variable Horizontal Grid
Integer, Parameter  :: max_halo_size = 8
! Comdecks
#include "amaxsize.h"
Real, save          :: lambda_input_p( row_length_max )
Real, save          :: lambda_input_u( row_length_max )
Real, save          :: phi_input_p( rows_max )
Real, save          :: phi_input_v( rows_max )
  
Contains

Subroutine rcf_readnl_horizont( nft, nft_horizgrid )
    
Use Rcf_Interp_Weights_Mod, Only : &
    h_int_method,              &
    bilinear,                  &
    area_weighted,             &
    nearest_neighbour,         &
    smcp_int_method

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Oper

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

Use Rcf_Parvars_Mod, Only : &
    mype

Use Rcf_Recon_Mod, Only : &
    Rimwidtha

Use Rcf_cntlatm_mod, Only : &
    model_domain

Use Rcf_FortranIO_Mod, Only : &
    Rcf_Free_Unit,            &
    Max_Filename_Len

Implicit None

! Arguments
Integer, Intent(In)              :: nft             ! File Unit
Integer, Intent(In)              :: nft_horizgrid   ! File Unit
Integer                          :: status
Logical                          :: l_exist

! Comdecks
#include "c_mdi.h"
#include "domtyp.h"
 
! Local vars/params
Character (Len=*), Parameter     :: RoutineName = 'Rcf_readnl_horizont'
Character (Len=Max_Filename_Len) :: FileName
Character (Len=Max_Filename_Len) :: FileName_horizgrid
Character (Len=80)               :: Cmessage
Integer                          :: ErrorStatus
Integer, Parameter               :: orog_blend_max = 25

! Temp storage for namelist information
Real                             :: orog_blend_weights(orog_blend_max)
Logical                          :: global
Logical                          :: L_VARGRID

Namelist /HORIZONT/ global, L_VARGRID, h_int_method, delta_lambda,     &
                    delta_phi, lambda_first, phi_first, lambda_npole,  &
                    phi_npole, iproj,                                  &
                    orog_blend_width, orog_blend_weights,              &
                    extended_halo_size_EW, extended_halo_size_NS,      &
                    smcp_int_method

! Variable Horizontal Grid
Namelist /HORIZGRID/ lambda_input_p, lambda_input_u,                   &
                     phi_input_p, phi_input_v 
! Set defaults
smcp_int_method  = IMDI
global           = .TRUE.
L_VARGRID        = .FALSE.
h_int_method     = bilinear
delta_lambda     = RMDI
delta_phi        = RMDI
lambda_first     = 0.
phi_first        = -90.
lambda_npole     = 0.
phi_npole        = 90.
iproj            = IMDI
orog_blend_width = 0
orog_blend_weights(:)     = 0.0
extended_halo_size_EW     = 0
extended_halo_size_NS     = 0
lambda_input_p(:)      = 0.0
lambda_input_u(:)      = 0.0
phi_input_p(:)         = 0.0
phi_input_v(:)         = 0.0
 
! Read horizont Namelist
Read( Unit = nft, Nml=horizont )

! Write out namelist for diagnostic
If (PrintStatus >= PrStatus_Oper) Then
  If ( mype == 0 ) Then
    Write ( 6, horizont )
  End If
End If

! ------------------------
! Horizontal Grid Namelist
! ------------------------

! First, find the namelist filename from Env Vars
Call Fort_Get_Env( 'RCF_NAMELIST', 12, FileName,                   &
                    Max_Filename_Len, ErrorStatus )

If (L_VARGRID) then

! Find the filename containing horizontal grid from Env Vars
  Call Fort_Get_Env('VAR_GRID', 8, FileName_horizgrid,              &
                     Max_Filename_Len, ErrorStatus)

  If ( ErrorStatus /= 0 ) Then
    ErrorStatus = 30
    Cmessage =  &
    'Unable to Obtain horizontal grid Filename from Environment'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

  FileName_horizgrid = Trim( FileName_horizgrid )

! Check file exists - Do all this by hand (not using file_open)
! as wish to do this on *all* PEs together
  Inquire( file=FileName_horizgrid, exist=l_exist )

  If ( .Not. l_exist ) Then
    Write (6,*) 'Horizontal grid file: ',FileName_horizgrid 
    ErrorStatus = 40
    Cmessage = ' Horizontal grid Namelist file does not exist!'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If
 
! Read horizgrid Namelist 
! Open the file containing horizontal grid
  Open( Unit=nft_horizgrid, File=FileName_horizgrid, IOstat=status ) 

! Write out namelist for diagnostic
  If ( PrintStatus >= PrStatus_Oper ) Then
    If ( mype == 0 ) Then  
      Write (6,*) 'horizontal grid file: ',FileName_horizgrid
    End If
  End If 
  Read( Unit = nft_horizgrid, Nml=horizgrid )
  
! Set the followings to BMDI for var grids model.
  delta_lambda     = RMDI
  delta_phi        = RMDI
  lambda_first     = RMDI
  phi_first        = RMDI 
  
End If  

! Some checks
If ( h_int_method /= bilinear .AND. h_int_method /= area_weighted ) Then
  Cmessage = 'Unsupported horizontal interpolation method'
  ErrorStatus = 10
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! Can't have a orography blending for a global model.
If (orog_blend_width /= 0 .AND. Global) Then
  Cmessage='orog_blend_width is not zero. This is not right for&
          & a global model'
  ErrorStatus = 20
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! Check that we have don't have too many weights for allocated
! array sizes
If ( orog_blend_width > orog_blend_max) Then
  Write (Cmessage,*) 'Current blending weight limit is ',           &
                      orog_blend_max,'. Please increase orog_blend_max'
  ErrorStatus = 50
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! Allocate and set the internal weights for orography blending
If (orog_blend_width > 0) Then
  Allocate( blend_weights( orog_blend_width ) )

  ! blending zone weights come from namelist
  blend_weights( 1 : orog_blend_width ) = &
                 orog_blend_weights( 1 : orog_blend_width )
End If


! Fill in gaps about Output Grid
Output_Grid % global = global

If (global) Then   ! Atmos C grid global
  Output_Grid % glob_u_row_length = Output_Grid % glob_p_row_length
  Output_Grid % glob_v_row_length = Output_Grid % glob_p_row_length
  Output_Grid % glob_u_rows       = Output_Grid % glob_p_rows
  Output_Grid % glob_v_rows       = Output_Grid % glob_p_rows - 1
  Output_Grid % rotated           = .FALSE.
Else    ! Atmos C grid LAM
! Commented out "correct size" wrapping LAM
!  If ( Output_Grid % glob_p_row_length * Delta_Lambda > 359.99 ) Then
! wrapping LAM
    Output_Grid % glob_u_row_length = Output_Grid % glob_p_row_length
!  Else
! non-wrapping  LAM
!    Output_Grid % glob_u_row_length = Output_Grid % glob_p_row_length-1
!  End If
  Output_Grid % glob_v_row_length = Output_Grid % glob_p_row_length
  Output_Grid % glob_u_rows       = Output_Grid % glob_p_rows
  IF (model_domain == mt_bi_cyclic_lam) THEN
    Output_Grid % glob_v_rows       = Output_Grid % glob_p_rows
  Else
    Output_Grid % glob_v_rows       = Output_Grid % glob_p_rows - 1
  Endif
  If (Phi_NPole /= 90 .OR. Lambda_NPole /= 0 ) Then
    Output_Grid % Rotated         = .TRUE.
  Else
    Output_Grid % Rotated         = .FALSE.
  End If
End If

Output_Grid % glob_p_field = Output_Grid % glob_p_row_length *    &
                             Output_Grid % glob_p_rows
Output_Grid % glob_u_field = Output_Grid % glob_u_row_length *    &
                             Output_Grid % glob_u_rows
Output_Grid % glob_v_field = Output_Grid % glob_v_row_length *    &
                             Output_Grid % glob_v_rows
Output_Grid % glob_r_field = Output_Grid % glob_r_row_length *    &
                             Output_Grid % glob_r_rows

Close( Unit=nft_horizgrid)

Call Rcf_Free_Unit( nft_horizgrid)

Return
End Subroutine rcf_readnl_horizont
End Module rcf_readnl_horizont_mod
#endif
