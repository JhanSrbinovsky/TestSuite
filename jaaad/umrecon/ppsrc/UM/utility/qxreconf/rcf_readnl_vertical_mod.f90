
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ reads the vertical namelist

Module rcf_readnl_vertical_Mod

!  Subroutine Rcf_Readnl_vertical - reads the vertical namelist
!
! Description:
! Module to read in the vertical namelist
! Note that it *must* be called after readnl_recona as
! some sizing is required!
!
! Method:
!  Namelist read in and sets Output_Grid variables appopriately
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


! Some data we need to carry around before it can be assigned
! properly to the output

Integer, Parameter   :: max_levels = 100   ! For temporary array sizing

Contains

Subroutine Rcf_Readnl_Vertical( nft_vertical, nft_vertlevs )

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

Use Rcf_V_Int_Ctl_Mod, Only : &
    v_int_order,          &
    Linear,               &
    Linear_NoEx,          &
    Cubic,                &
    Quintic

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Oper

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Parvars_mod, Only : &
    mype

Use Rcf_Generate_Heights_Mod, Only : &
    height_gen_smooth,               &
    height_gen_original

Implicit None

! Subroutine arguments
Integer, Intent(In)          :: nft_vertical  ! File unit - VERTICAL
Integer, Intent(In)          :: nft_vertlevs  ! File unit - VERTLEVS

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

! local variables/constants
Character (Len=*), Parameter :: RoutineName = 'Rcf_Readnl_Vertical'
Character (Len=80)           :: Cmessage
Integer                      :: ErrorStatus
Integer                      :: levels_theta  ! number of levels of
Integer                      :: levels_rho    ! values read in.
Integer                      :: i             ! looper
Real                         :: last_value    ! used in checking code

! temporary storage for
Real                         :: eta_theta( max_levels  + 1)
Real                         :: eta_rho( max_levels )
Real                         :: rh_crit( max_levels )
Real                         :: soil_depths( max_levels )
Real                         :: z_top_of_model
Integer                      :: first_constant_r_rho_level
Integer                      :: height_method

Namelist /VERTICAL/ rh_crit, soil_depths, v_int_order, height_method

Namelist /VERTLEVS/ z_top_of_model, eta_theta, eta_rho,         &
                    first_constant_r_rho_level

! Set VERTICAL defaults
rh_crit(:)       = RMDI
soil_depths(:)   = 0.0
soil_depths(1:4) = (/0.1, 0.25, 0.65, 2.0/)
v_int_order      = Linear
height_method    = height_gen_smooth

! Set VERTLEVS defaults
eta_theta(:)          = RMDI
eta_rho(:)            = RMDI

! Quick error check to make sure parameter max_levels isn't too small
If ( max_levels < Output_Grid % model_levels ) Then
  ErrorStatus = 10
  Cmessage = 'Internal paramter max_levels is too small!'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! Read VERTICAL Namelist
Read( Unit = nft_vertical, Nml=vertical )

! Read VERTLEVS Namelist
Read( Unit = nft_vertlevs, Nml=vertlevs )

! Write out namelist for diagnostic
If (PrintStatus >= PrStatus_Oper) Then
  If ( mype == 0 ) Then
    Write ( Unit = 6, Nml = vertical)

    Write ( Unit = 6, Fmt = * ) 'Values in VERTLEVS Namelist.'
    Write ( Unit = 6, Fmt = * ) 'z_top_of_model ',z_top_of_model
    Write ( Unit = 6, Fmt = * ) 'first_constant_r_rho_level ',    &
                                 first_constant_r_rho_level
    Write ( Unit = 6, Fmt = * ) 'Eta_Theta'
    Write ( Unit = 6, Fmt = '(5F10.7)' )                          &
          ( eta_theta(i),i=1,Output_Grid % Model_levels+1)
    Write ( Unit = 6, Fmt = * ) 'Eta_Rho'
    Write ( Unit = 6, Fmt = '(5F10.7)' )                          &
          ( eta_rho(i),i=1,Output_Grid % Model_levels)

  End If
End If

! Some sanity checks before we allocate space
If ( Output_Grid % model_levels < 1 ) Then
  Cmessage = 'Model levels for output grid is < 1!'
  ErrorStatus = 20
  Call Ereport( RoutineName, ErrorStatus, Cmessage )

Else If ( Output_Grid % model_levels > 1000 ) Then
  Cmessage = 'Model Levels for output grid is > 1000 - Is this correct?'
  ErrorStatus = -30
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

If ( v_int_order /= Linear .AND. v_int_order /= Linear_NoEx .AND. &
     v_int_order /= Cubic  .AND. v_int_order /= Quintic ) Then
  Cmessage = 'Vertical Interpolation Order not recognised'
  ErrorStatus = 40
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! Check that the number of eta values correspond to the number
! of model levels
levels_theta = 0
levels_rho   = 0
Do i = 1, max_levels
  If (eta_theta(i) /= RMDI) Then
    levels_theta = levels_theta + 1
  End If

  If (eta_rho(i) /= RMDI) Then
    levels_rho = levels_rho + 1
  End If
End Do

If ( levels_theta /= Output_Grid % model_levels + 1) Then
  Cmessage = 'Mismatch in number of model levels and number of &
            &eta_theta values'
  ErrorStatus = 50
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

If ( levels_rho /= Output_Grid % model_levels ) Then
  Cmessage = 'Mismatch in number of model levels and number of &
             &eta_rho values'
  ErrorStatus = 60
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! Check that eta_rho levels and eta_theta levels are monotone
! ascending.
last_value = RMDI
Do i = 1, Output_Grid % model_levels + 1
  If ( eta_theta(i) <= last_value ) Then
    Cmessage = 'Eta_Theta values are not monotone ascending'
    ErrorStatus = 70
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

  last_value = eta_theta(i)
End Do

last_value = RMDI
Do i = 1, Output_Grid % model_levels
  If ( eta_rho(i) <= last_value ) Then
    Cmessage = 'Eta_Rho values are not monotone ascending'
    ErrorStatus = 80
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

  last_value = eta_rho(i)
End Do

  If ( height_method /= height_gen_original .AND.   &
       height_method /= height_gen_smooth ) Then
    ErrorStatus = 90
    Cmessage = 'Unrecognised height generation method'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

! Allocate space and fill up relevant parts of Output_Grid
Allocate (Output_Grid % eta_theta_levels                   &
                       ( 0 : Output_Grid % model_levels) )
Allocate ( Output_Grid % eta_rho_levels( Output_Grid % model_levels ) )
Allocate ( Output_Grid % rh_crit( Output_Grid % model_levels ) )
Allocate ( Output_Grid % soil_depths( Output_Grid % sm_levels ) )

Output_Grid % eta_theta_levels( 0 : Output_Grid % model_levels ) = &
                     eta_theta( 1 : Output_Grid % model_levels + 1)
Output_Grid % eta_rho_levels( 1 : Output_Grid % model_levels ) =   &
                     eta_rho( 1 : Output_Grid % model_levels )
Output_Grid % rh_crit( 1 : Output_Grid % model_levels ) =          &
              rh_crit( 1 : Output_Grid % model_levels )
Output_Grid % soil_depths( 1 : Output_Grid % sm_levels ) =         &
              soil_depths( 1 : Output_Grid % sm_levels )
Output_Grid % z_top_of_model = z_top_of_model
Output_Grid % first_constant_r_rho_level = first_constant_r_rho_level
Output_Grid % height_gen_method = height_method



Return
End Subroutine Rcf_Readnl_Vertical
End Module Rcf_Readnl_Vertical_Mod
