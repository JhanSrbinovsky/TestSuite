#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Generate SI dynamics height field

Module Rcf_generate_heights_mod

!  Subroutine Rcf_Generate_Heights - generate a height field
!
! Description:
!   Generates a height field from eta-theta and eta-rho levels and
!   orography field.
!
! Method:
!   Heights are generated for theta and rho levels on theta points.
!   Either the relevant set is passed back or, futher processing is
!   done eg to average to get heights on u or v points or to get
!   zonal points. Soil depths can also be generated.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.2    6/03/00   Choice of 2 generation schemes. P.Selwood.
!   5.3    5/10/01   Cater for extra rho level at top. D.Robinson
!   5.3   28/01/02   Use heights on theta levels for ozone. D. Robinson
!   5.4    4/09/02   Use top of atmosphere for ozone. P.Selwood
!   5.4   29/08/02   Enable use of bicyclic LAM. C. Roadnight
!   5.4   11/03/02   Remove comment lines from after #include
!                                                  S. Carroll
!   5.4   09/09/02   Add 2 new height gen method parameters. R.Sharp
!   5.5   09/01/03   Add new height generation subroutine. R.Sharp
!   6.0   04/08/03   Add new height generation subroutine for ECMWF
!                    GRIB model level data and modified subroutine
!                    for ECMWF GRIB pressure level data. P.Earnshaw
!   6.2   20/06/05   Removed hardwired ECMWF model level definitions,
!                    now read from GRIB file. Paul Earnshaw
!   6.2   21/10/05   Replace GSYNC with SSYNC. P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


! Global Parameters - for methods of generating height fields

Integer, Parameter    ::  height_gen_original    =  1
Integer, Parameter    ::  height_gen_smooth      =  2
Integer, Parameter    ::  height_gen_ecmwf_press = 10
Integer, Parameter    ::  height_gen_ecmwf_hybrd = 11

Contains

Subroutine Rcf_generate_heights( grid, orography,  gridcode, levelT, &
                                 heights, levelsize )

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Use Ereport_mod, Only : &
    Ereport

Use Rcf_Recon_Mod, ONly : &
    lozone_zonal

Use Rcf_Parvars_Mod, Only : &
    gc_proc_row_group,      &
    gc_all_proc_group,      &
    mype,                   &
    nproc,                  &
    bound,                  &
    BC_CYCLIC

Use Rcf_Scatter_Field_Mod, Only : &
    Rcf_Scatter_Field_Real

Use Rcf_Gather_Field_Mod, Only : &
    Rcf_Gather_Field_Real

Use Rcf_PrintStatus_Mod, Only : &
    LTimer

Implicit None

! Arguments
Type (grid_type), Intent(In)   :: grid
Real, Intent(In)               :: orography( grid % loc_p_field )
Integer, Intent(In)            :: gridcode
Integer, Intent(In)            :: levelT
Integer, Intent(In)            :: levelsize
Real, Intent(Out)              :: heights( levelsize,                 &
                                    0 : grid % model_levels + 1)

! Comdecks
! earth_radius
#include "c_a.h"
#include "cppxref.h"

! Local Data
Real                           :: denom    ! denominator in relax
Real                           :: etk_etbl ! quantaties pulled out of
Real                           :: erk_etbl ! loop for optimisation
Real                           :: r_ref_theta( grid % model_levels )
Real                           :: r_ref_rho  ( grid % model_levels )
Real                           :: r_theta_levels( grid % loc_p_field, &
                                               0 : grid % model_levels )
Real                           :: r_rho_levels( grid % loc_p_field,   &
                                                grid % model_levels + 1)
Real, Allocatable              :: r_level(:,:)
Real, Allocatable              :: u_level(:,:)
Real, Allocatable              :: v_level(:,:)

Integer                        :: istat
Integer                        :: i
Integer                        :: j
Integer                        :: k
Integer                        :: div
Integer                        :: rem
Integer                        :: pe
Integer                        :: ozone_lev  ! local ozone level

Integer                        :: ErrorStatus
Character (Len=*), Parameter   :: RoutineName = 'Rcf_Generate_heights'
Character (Len=80)             :: Cmessage

External gcg_rvecsumr, Timer

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 3)

! First check to see if eta values are set.
If (.NOT. Associated( grid % eta_theta_levels ) .OR. &
    .NOT. Associated( grid % eta_rho_levels ) ) Then
  Cmessage = 'Need to set eta values before heights can be generated'
  ErrorStatus = 10
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! Initialise heights for safety
heights( :, :) = 0.

!------------------------------------------------------------------
! Setup basic theta and rho heights - 3 methods currently
! available for this
!------------------------------------------------------------------

Select Case (grid % height_gen_method)
  Case( height_gen_original)
    Call Rcf_Gen_Height_Original( grid, orography, levelT,      &
                                  r_theta_levels, r_rho_levels )

  Case( height_gen_smooth)
    Call Rcf_Gen_Height_Smooth( grid, orography, levelT,      &
                                r_theta_levels, r_rho_levels )

  Case( height_gen_ecmwf_press)
    Call Rcf_Gen_Height_Pressure( grid, orography, levelT,      &
                                  r_theta_levels, r_rho_levels )

  Case( height_gen_ecmwf_hybrd)
    Call Rcf_Gen_Height_Hybrid( grid, orography, levelT,      &
                                  r_theta_levels, r_rho_levels )

  Case Default
    ErrorStatus = 15
    Cmessage = 'Height generation method is not recognised'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )

End Select

!------------------------------------------------------------------
! Decide which height field we need here
!------------------------------------------------------------------

Select Case ( gridcode )
!------------------------------------------------------------------
! Theta grid
!------------------------------------------------------------------
  Case ( ppx_atm_tall )
    If ( levelT == ppx_theta_level ) Then
      ! We want r_theta_levels
      heights( :, 0 : grid % model_levels ) = r_theta_levels( :, :)
    Else If ( levelT == ppx_rho_level ) Then
      ! Rho levels
      heights( :, 1 : grid % model_levels + 1) = r_rho_levels( :, :)
    Else
      ! Don't know what to do with neither theta nor rho levels
      Write (6,*) 'levelT = ', levelT
      Cmessage = 'LevelT is neither full nor half'
      ErrorStatus = 20
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

!-------------------------------------------------------------------
! U and V grids
!-------------------------------------------------------------------
  Case (ppx_atm_cuall, ppx_atm_cvall )
    Allocate( r_level( grid % glob_p_row_length, grid % glob_p_rows ) )
    Allocate( u_level( grid % glob_u_row_length, grid % glob_u_rows ) )
    Allocate( v_level( grid % glob_v_row_length, grid % glob_v_rows ) )

    ! U and V need to do calculations on P fields (linear average)
    ! So do Gather 1 level per PE and then Scatter again.
    div = grid % model_levels / nproc
    rem = Mod( grid % model_levels, nproc )
    pe = 0

    Do i = 1, div
      Do j = ((i-1) * nproc) + 1, i * nproc
        ! Will gather level j on PE pe
        Call Rcf_Gather_Field_Real( r_rho_levels( :, j), r_level, &
                               grid % loc_p_row_length,      &
                               grid % loc_p_rows,            &
                               grid % glob_p_row_length,     &
                               grid % glob_p_rows, pe,       &
                               gc_all_proc_group )

        pe = pe + 1
        If (pe == nproc) pe = 0
      End Do

      Call Gc_Ssync( nproc, ErrorStatus )

      ! U and V need to be done seperately
      ! Here do U
      If (gridcode == ppx_atm_cuall) Then
        Do k = 1, grid % glob_u_rows
          Do j = 1, grid % glob_u_row_length - 1
            u_level(j,k) = (r_level(j,k) + r_level(j+1,k)) / 2.0
          End Do

          ! If cyclic take the correct average, otherwise just the
          ! nearest point for the `extra' u point
          If (bound(1) == BC_CYCLIC) Then
            u_level( grid % glob_u_row_length,k) = (r_level(1,k) + &
                     r_level( grid % glob_p_row_length,k) ) / 2.0
          Else
            u_level( grid % glob_u_row_length,k) = &
                     r_level( grid % glob_p_row_length,k)
          End If
        End Do

        Call Gc_Ssync( nproc, ErrorStatus )

        Do j = ((i-1) * nproc) + 1, i * nproc

          Call Rcf_Scatter_Field_Real( heights(:, j), u_level,       &
                                  grid % loc_u_row_length,      &
                                  grid % loc_u_rows,            &
                                  grid % glob_u_row_length,     &
                                  grid % glob_u_rows, pe,       &
                                  gc_all_proc_group )


          pe = pe + 1
          If (pe == nproc) pe = 0
        End Do
      End If

      ! Now do V
      If (gridcode == ppx_atm_cvall) Then
        Do k = 1, grid % glob_v_rows - 1
          Do j = 1, grid % glob_v_row_length
            v_level(j,k) = ( r_level(j,k) + r_level(j,k+1) ) /2.0
          End Do
        End Do

        If (bound(2) == BC_CYCLIC) then
          Do j = 1, grid % glob_v_row_length
            v_level(j,grid % glob_v_rows) = &
            (r_level(j,1) + r_level(j,grid % glob_p_rows))/2.0
          End Do
        else
          k = grid % glob_v_rows
          Do j = 1, grid % glob_v_row_length
            v_level(j,grid % glob_v_rows) = &
            (r_level(j,k) + r_level(j,k+1))/2.0
          End Do
        Endif


        Call Gc_Ssync( nproc, ErrorStatus )

        Do j = ((i-1) * nproc) + 1, i * nproc

          Call Rcf_Scatter_Field_Real( heights(:, j), v_level,       &
                                  grid % loc_v_row_length,      &
                                  grid % loc_v_rows,            &
                                  grid % glob_v_row_length,     &
                                  grid % glob_v_rows, pe,       &
                                  gc_all_proc_group )


          pe = pe + 1
          If (pe == nproc) pe = 0
        End Do
      End If
    End Do

!-------------------------------------------------------------------
! There are rem levels left to process. Will do these now.
!-------------------------------------------------------------------
    pe = 0
    Do i = 1, rem
      j = nproc * div + i
      ! Will gather level j on PE pe
      Call Rcf_Gather_Field_Real( r_rho_levels( :, j), r_level, &
                             grid % loc_p_row_length,           &
                             grid % loc_p_rows,                 &
                             grid % glob_p_row_length,          &
                             grid % glob_p_rows, pe,            &
                             gc_all_proc_group )

      pe = pe + 1
    End Do

    Call Gc_Ssync( nproc, ErrorStatus )

    ! U and V need to be done seperately
    ! Here do U
    If (gridcode == ppx_atm_cuall) Then
      If (mype < pe) Then
        Do k = 1, grid % glob_u_rows
          Do j = 1, grid % glob_u_row_length - 1
            u_level(j,k) = (r_level(j,k) + r_level(j+1,k)) / 2.0
          End Do

          ! If cyclic take the correct average, otherwise just the
          ! nearest point for the `extra' u point
          If (bound(1) == BC_CYCLIC) Then
            u_level( grid % glob_u_row_length,k) = (r_level(1,k) + &
                     r_level( grid % glob_p_row_length,k) ) / 2.0
          Else
            u_level( grid % glob_u_row_length,k) = &
                     r_level( grid % glob_p_row_length,k)
          End If
        End Do
      End If

      Call Gc_Ssync( nproc, ErrorStatus )
      pe = 0
      Do i = 1, rem
        j = nproc * div + i

        Call Rcf_Scatter_Field_Real( heights(:, j), u_level,  &
                                grid % loc_u_row_length,      &
                                grid % loc_u_rows,            &
                                grid % glob_u_row_length,     &
                                grid % glob_u_rows, pe,       &
                                gc_all_proc_group )
        pe = pe + 1
        If (pe == nproc) pe = 0
       End Do
    End If

    ! Here do V
    If (gridcode == ppx_atm_cvall) Then
      If (mype < pe) Then
        Do k = 1, grid % glob_v_rows - 1
          Do j = 1, grid % glob_v_row_length
            v_level(j,k) = ( r_level(j,k) + r_level(j,k+1) ) /2.0
          End Do
        End Do

        If (bound(2) == BC_CYCLIC) then
          Do j = 1, grid % glob_v_row_length
            v_level(j,grid % glob_v_rows) = &
            (r_level(j,1) + r_level(j,grid % glob_p_rows))/2.0
          End Do
        else
          k = grid % glob_v_rows
          Do j = 1, grid % glob_v_row_length
            v_level(j,grid % glob_v_rows) = &
            (r_level(j,k) + r_level(j,k+1))/2.0
          End Do
        Endif

      End If

      Call Gc_Ssync( nproc, ErrorStatus )
      pe = 0
      Do i = 1, rem
        j = nproc * div + i

        Call Rcf_Scatter_Field_Real( heights(:, j), v_level,  &
                                grid % loc_v_row_length,      &
                                grid % loc_v_rows,            &
                                grid % glob_v_row_length,     &
                                grid % glob_v_rows, pe,       &
                                gc_all_proc_group )
        pe = pe + 1
        If (pe == nproc) pe = 0
      End Do
    End If

    Deallocate( r_level )
    Deallocate( u_level )
    Deallocate( v_level )

!--------------------------------------------------------------------
! Ozone grid
!-------------------------------------------------------------------
  Case ( ppx_atm_ozone )
    If ( LOZONE_ZONAL ) Then
      Do j = 1, grid % ozone_levels

        ozone_lev = grid % model_levels - grid % ozone_levels + j

        ! Sum rows of r_theta_levels to create heights
        Call gcg_rvecsumr( grid % loc_p_row_length,      &
                           grid % loc_p_row_length,      &
                           1, grid % loc_p_rows,         &
                           r_theta_levels(1,ozone_lev),  &
                           gc_proc_row_group,            &
                           istat, heights(1,j) )

        ! Create the mean
        Do i = 1, grid % loc_p_rows
          heights(i,j) = heights(i,j)/grid % glob_p_row_length
        End Do
      End Do
    Else
      heights( :, 1 : grid % ozone_levels) =                          &
        r_theta_levels( :, grid % model_levels -                      &
                           grid % ozone_levels + 1 :                  &
                           grid % model_levels )
    End If

!--------------------------------------------------------------------
! Deep Soil grid
!--------------------------------------------------------------------
  Case ( ppx_atm_compressed )
    If ( levelT == ppx_soil_level ) Then
    ! Not heights so much as soil level depths...
      Do i = 1, levelsize
        heights( i, 1 : grid % sm_levels) = grid % soil_depths( : )
      End Do

    Else
      Write(6,*) 'Land compressed field, levelT=', levelT
      Cmessage = 'Unsupported grid/level combination for height &
                 &calculation'
      ErrorStatus = 30
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

  Case Default
    Write (6,*) 'Gridcode = ', gridcode
    Cmessage = 'Unsupported grid type for height calculation'
    ErrorStatus = 40
    Call Ereport( RoutineName, ErrorStatus, Cmessage )

End Select

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 4)

Return
End Subroutine Rcf_generate_heights


!**********************************************************************
! Rcf_Gen_Height_Original - subroutine to generate original set of
! heights.
!**********************************************************************
Subroutine Rcf_Gen_Height_Original( grid, orography, levelT,      &
                                    r_theta_levels, r_rho_levels)

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Implicit None

! Arguments
Type (grid_type), Intent(In)   :: grid
Real, Intent(In)               :: orography( grid % loc_p_field )
Integer, Intent(In)            :: levelT
Real, Intent(Out)              :: r_theta_levels( grid % loc_p_field, &
                                                0 : grid % model_levels)
Real, Intent(Out)              :: r_rho_levels( grid % loc_p_field,   &
                                          grid % model_levels + 1)

! Comdecks
! earth radius
#include "c_a.h"
#include "cppxref.h"

! Local Variables
Integer                        :: k       ! looper
Integer                        :: j       ! looper
Real                           :: denom    ! denominator in relax
Real                           :: etk_etbl ! quantaties pulled out of
Real                           :: erk_etbl ! loop for optimisation
Real                           :: r_ref_theta( grid % model_levels )
Real                           :: r_ref_rho  ( grid % model_levels )

!---------------------------------------------------------------------
! Set reference height profile
!---------------------------------------------------------------------
Do k = 1, grid % model_levels
  r_ref_theta(k) = grid % eta_theta_levels(k) * grid % z_top_of_model
  r_ref_rho(k)   = grid % eta_rho_levels(k) * grid % z_top_of_model
End Do

!---------------------------------------------------------------------
! set bottom level, ie: orography
!---------------------------------------------------------------------
Do j = 1, grid % loc_p_field
  r_theta_levels(j,0) =  orography(j) + Earth_radius
End Do

!---------------------------------------------------------------------
! For boundary layer levels set depth to be constant.
!---------------------------------------------------------------------
If (levelT == ppx_theta_level) Then
  Do k = 1, grid % bl_levels
    Do j = 1, grid % loc_p_field
      r_theta_levels(j,k) = r_theta_levels(j,0) + r_ref_theta(k)
    End Do
  End Do
End If

If (levelT == ppx_rho_level) Then
  Do k = 1, grid % bl_levels
    Do j = 1, grid % loc_p_field
      r_rho_levels(j,k)   = r_theta_levels(j,0) + r_ref_rho(k)
    End Do
  End Do

  ! Need theta_levels at bl_levels
  k = grid % bl_levels
  Do j = 1, grid % loc_p_field
    r_theta_levels(j,k)   = r_theta_levels(j,0) + r_ref_theta(k)
  End Do
End If

!---------------------------------------------------------------------
! For constant levels set r to be a constant on the level
!---------------------------------------------------------------------
  Do k = grid % first_constant_r_rho_level, grid % model_levels
    Do j = 1, grid % loc_p_field
      r_theta_levels(j,k) = Earth_radius + r_ref_theta(k)
    End Do
  End Do

  k= grid % first_constant_r_rho_level
  Do j = 1, grid % loc_p_field
    r_rho_levels(j,k)   = Earth_radius + r_ref_rho(k)
  End Do

  Do k = grid % first_constant_r_rho_level, grid % model_levels
    Do j = 1, grid % loc_p_field
      r_rho_levels(j,k)   = Earth_radius + r_ref_rho(k)
    End Do
  End Do

!---------------------------------------------------------------------
! For intermediate levels use linear relaxation to constant value.
! set orographic heights.
!---------------------------------------------------------------------
denom = (grid % eta_rho_levels(grid % first_constant_r_rho_level) -&
         grid % eta_theta_levels(grid % bl_levels) )

If ( levelT == ppx_rho_level) Then
  Do k = grid % bl_levels + 1, grid % first_constant_r_rho_level - 1

    erk_etbl = ( grid % eta_rho_levels(k) -                           &
             grid % eta_theta_levels(grid % bl_levels) )

    Do j = 1, grid % loc_p_field
      r_rho_levels(j,k) =                                             &
           ( r_rho_levels(j,grid % first_constant_r_rho_level) -      &
             r_theta_levels(j,grid % bl_levels) ) *                   &
             erk_etbl /                                               &
             denom                                                    &
               +  r_theta_levels(j,grid % bl_levels)
    End Do
  End Do
End If

If (levelT == ppx_theta_level) Then
  Do k = grid % bl_levels + 1, grid % first_constant_r_rho_level - 1

    etk_etbl = ( grid % eta_theta_levels(k) -                         &
             grid % eta_theta_levels(grid % bl_levels) )

    Do j = 1, grid % loc_p_field
      r_theta_levels(j,k) =                                           &
           ( r_rho_levels(j,grid % first_constant_r_rho_level) -      &
             r_theta_levels(j,grid % bl_levels) ) *                   &
             etk_etbl /                                               &
             denom                                                    &
               +  r_theta_levels(j,grid % bl_levels)
    End Do
  End Do
End If

!---------------------------------------------------------------------
! Derive height for extra rho level at top from top rho & theta levels
!---------------------------------------------------------------------

If ( levelT == ppx_rho_level) Then
  k = grid % model_levels
  Do j = 1, grid % loc_p_field
    r_rho_levels(j,k+1) = r_theta_levels(j,k) +    &
     ( r_theta_levels(j,k) - r_rho_levels(j,k) )
  End Do
End If

Return
End Subroutine Rcf_Gen_Height_Original

!**********************************************************************
! Rcf_Gen_Height_Smooth - subroutine to generate smooth set of
! heights.
!**********************************************************************
Subroutine Rcf_Gen_Height_Smooth( grid, orography, levelT,      &
                                  r_theta_levels, r_rho_levels )

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Implicit None

! Arguments
Type (grid_type), Intent(In)   :: grid
Real, Intent(In)               :: orography( grid % loc_p_field )
Integer, Intent(In)            :: levelT
Real, Intent(Out)              :: r_theta_levels( grid % loc_p_field, &
                                                0 : grid % model_levels)
Real, Intent(Out)              :: r_rho_levels( grid % loc_p_field,   &
                                                grid % model_levels + 1)

! Comdecks
! earth radius
#include "c_a.h"
#include "cppxref.h"

! Local Variables
Integer                        :: k       ! looper
Integer                        :: j       ! looper
Real                           :: r_ref_theta( grid % model_levels )
Real                           :: r_ref_rho  ( grid % model_levels )


!---------------------------------------------------------------------
! Set reference height profile
!---------------------------------------------------------------------
Do k = 1, grid % model_levels
  r_ref_theta(k) = grid % eta_theta_levels(k) * grid % z_top_of_model
  r_ref_rho(k)   = grid % eta_rho_levels(k) * grid % z_top_of_model
End Do

!---------------------------------------------------------------------
! set bottom level, ie: orography
!---------------------------------------------------------------------
Do j = 1, grid % loc_p_field
  r_theta_levels(j,0) =  orography(j) + Earth_radius
End Do

!---------------------------------------------------------------------
! For constant levels set r to be a constant on the level
!---------------------------------------------------------------------
  Do k = grid % first_constant_r_rho_level, grid % model_levels
    Do j = 1, grid % loc_p_field
      r_theta_levels(j,k) = Earth_radius + r_ref_theta(k)
    End Do
  End Do

  Do k = grid % first_constant_r_rho_level, grid % model_levels
    Do j = 1, grid % loc_p_field
      r_rho_levels(j,k)   = Earth_radius + r_ref_rho(k)
    End Do
  End Do

!---------------------------------------------------------------------
! The rest of the levels are set with a quadratic scheme.
!---------------------------------------------------------------------
If ( levelT == ppx_rho_level) Then
  Do k = 1, grid % first_constant_r_rho_level - 1
    Do j = 1, grid % loc_p_field
      r_rho_levels(j,k) = Earth_radius +                               &
               grid % eta_rho_levels(k) * grid % z_top_of_model +      &
               orography(j) * (1.0 - grid % eta_rho_levels(k) /        &
           grid % eta_rho_levels(grid % first_constant_r_rho_level))**2
    End Do
  End Do
End If

If ( levelT == ppx_theta_level) Then
  Do k = 1, grid % first_constant_r_rho_level - 1
    Do j = 1, grid % loc_p_field
      r_theta_levels(j,k) = Earth_radius +                             &
               grid % eta_theta_levels(k) * grid % z_top_of_model +    &
               orography(j) * (1.0 - grid % eta_theta_levels(k) /      &
           grid % eta_rho_levels(grid % first_constant_r_rho_level))**2
    End Do
  End Do
End If

!---------------------------------------------------------------------
! Derive height for extra rho level at top from top rho & theta levels
!---------------------------------------------------------------------

If ( levelT == ppx_rho_level) Then
  k = grid % model_levels
  Do j = 1, grid % loc_p_field
    r_rho_levels(j,k+1) = r_theta_levels(j,k) +    &
     ( r_theta_levels(j,k) - r_rho_levels(j,k) )
  End Do
End If

Return
End Subroutine Rcf_Gen_Height_Smooth

!**********************************************************************
! Rcf_Gen_Height_Pressure - subroutine to generate set of heights
! from ECMWF Pressure levels.
!**********************************************************************
Subroutine Rcf_Gen_Height_Pressure( grid, orography, levelT,      &
                                    r_theta_levels, r_rho_levels  &
                                    )

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Use Rcf_GRIB_T_n_Pstar_H_Interp_Mod, Only :  &
    GRIB_T,             &
    GRIB_Pstar,         &
    GRIB_Levels

Implicit None

! Arguments
Type (grid_type), Intent(In)   :: grid
Real, Intent(In)               :: orography( grid % loc_p_field )
Integer, Intent(In)            :: levelT
Real, Intent(Out)              :: r_theta_levels( grid % loc_p_field, &
                                                0 : grid % model_levels)
Real, Intent(Out)              :: r_rho_levels( grid % loc_p_field,   &
                                          grid % model_levels + 1)

! Comdecks
! earth radius
#include "c_a.h"
! gravitational constant
#include "c_g.h"
! R (gas constant)
#include "c_r_cp.h"
#include "cppxref.h"

! Local Variables
Integer                        :: k       ! looper
Integer                        :: j       ! looper
Real                           :: Roverg,rtemp,logpp

!---------------------------------------------------------------------
! Calculate R/g
!---------------------------------------------------------------------
Roverg = R/g

!---------------------------------------------------------------------
! set bottom level, ie: orography
!---------------------------------------------------------------------
Do j = 1, grid % loc_p_field
  r_theta_levels(j,0) =  orography(j) + Earth_radius
End Do

!---------------------------------------------------------------------
! Calculate heights from Pstar and T
!---------------------------------------------------------------------

!Use hydrostatic approximation and ideal gas equation to derive
!expression for difference in height between two pressure levels.
!Temperature assumed constant between two successive levels.

!Calculate first level heights from pstar
!Temperature is from first level
!Do this because no reference data at surface
Do j = 1, grid % loc_p_field
  logpp = log( GRIB_Pstar % Data(j,1) / GRIB_Levels(1) )
  rtemp = Roverg * GRIB_T % Data(j,1) * logpp
  r_theta_levels(j,1) = r_theta_levels(j,0) + rtemp
  r_rho_levels(j,1)   = r_theta_levels(j,1)
End Do

!Calculate rest of the level heights from previous pressure level
!Temperature between levels is mean of temperatures at each level
Do k = 2, grid % model_levels
  logpp=0.5*Roverg*log( GRIB_Levels(k-1) / GRIB_Levels(k) )
  Do j = 1, grid % loc_p_field
    rtemp = logpp*(GRIB_T%Data(j,k)+GRIB_T%Data(j,k-1))
    r_theta_levels(j,k) = r_theta_levels(j,k-1) + rtemp
    r_rho_levels(j,k)   = r_theta_levels(j,k)
  End Do
End Do

!---------------------------------------------------------------------
! Derive height for extra rho level at top from top rho & theta levels
!---------------------------------------------------------------------

If ( levelT == ppx_rho_level) Then
  k = grid % model_levels
  Do j = 1, grid % loc_p_field
    r_rho_levels(j,k+1) = r_theta_levels(j,k) +    &
     ( r_theta_levels(j,k) - r_rho_levels(j,k-1) )
  End Do
End If

Return
End Subroutine Rcf_Gen_Height_Pressure

!**********************************************************************
! Rcf_Gen_Height_Hybrid - subroutine to generate set of heights
! from ECMWF Model levels.
!**********************************************************************
Subroutine Rcf_Gen_Height_Hybrid( grid, orography, levelT,      &
                                    r_theta_levels, r_rho_levels  &
                                    )

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Use Rcf_GRIB_T_n_Pstar_H_Interp_Mod, Only :  &
    GRIB_T,             &
    GRIB_Pstar,         &
    GRIB_Levels,  &
    ak, bk

Implicit None

! Arguments
Type (grid_type), Intent(In)   :: grid
Real, Intent(In)               :: orography( grid % loc_p_field )
Integer, Intent(In)            :: levelT
Real, Intent(Out)              :: r_theta_levels( grid % loc_p_field, &
                                                0 : grid % model_levels)
Real, Intent(Out)              :: r_rho_levels( grid % loc_p_field,   &
                                          grid % model_levels + 1)

! Comdecks
! earth radius
#include "c_a.h"
! gravitational constant
#include "c_g.h"
! R (gas constant)
#include "c_r_cp.h"
#include "cppxref.h"

! Local Variables
Integer                        :: k       ! looper
Integer                        :: j       ! looper
Real                           :: Roverg,rtemp,t_tmp,p_tmp1,p_tmp2

!---------------------------------------------------------------------
! Calculate R/g
!---------------------------------------------------------------------
Roverg = R/g

!---------------------------------------------------------------------
! set bottom level, ie: orography
!---------------------------------------------------------------------
Do j = 1, grid % loc_p_field
  r_theta_levels(j,0) =  orography(j) + Earth_radius
End Do

!---------------------------------------------------------------------
! Calculate heights from Pstar and T
!---------------------------------------------------------------------

!Use hydrostatic approximation and ideal gas equation to derive
!expression for difference in height between two model levels.

!Use ECMWF definition of model levels to obtain pressure at point
!on model level.

!Temperature assumed constant between two successive model levels.

!In following, GRIB_Levels(k) contains index of level k relative
!to ECMWF model level set but is a real number, hence need to use int
!to make it an integer.

!Calculate first level heights from pstar with temperature from
!first model level. Do this because no temperature data at surface.
Do j = 1, grid % loc_p_field

  !Calculate ratio of pressures of level 1 and surface
        p_tmp1 = ak(int(GRIB_Levels(1)))/GRIB_Pstar%Data(j,1) &
           + bk(int(GRIB_Levels(1)))

  !Calculate thickness of first level layer (note log(p_tmp1) < 0)
  rtemp = - Roverg * GRIB_T % Data(j,1) * log( p_tmp1 )

  !Add thickness of layer to surface height
  r_theta_levels(j,1) = r_theta_levels(j,0) + rtemp

  !Set r_rho equal to r_theta
  r_rho_levels(j,1)   = r_theta_levels(j,1)

End Do

!Calculate rest of the level heights from previous model level
!Temperature between levels is mean of temperatures at each level
Do k = 2, grid % model_levels
  Do j = 1, grid % loc_p_field

    !Calculate pressure at level k-1
    p_tmp1 = ak(int(GRIB_Levels(k-1))) +  &
             bk(int(GRIB_Levels(k-1)))*GRIB_Pstar%Data(j,1)

    !Calculate pressure at level k
    p_tmp2 = ak(int(GRIB_Levels(k))) +  &
             bk(int(GRIB_Levels(k)))*GRIB_Pstar%Data(j,1)

    !Calculate mean temperature across layer
    t_tmp = 0.5*( GRIB_T%Data(j,k) + GRIB_T%Data(j,k-1) )

    !Calculate thickness of layer between levels k and k-1
    rtemp = Roverg * t_tmp * log( p_tmp1 / p_tmp2 )

    !Add thickness of layer to level k-1 height
    r_theta_levels(j,k) = r_theta_levels(j,k-1) + rtemp

    !Set r_rho equal to r_theta
    r_rho_levels(j,k)   = r_theta_levels(j,k)

  End Do
End Do

!---------------------------------------------------------------------
! Derive height for extra rho level at top from top rho & theta levels
!---------------------------------------------------------------------

If ( levelT == ppx_rho_level) Then
  k = grid % model_levels
  Do j = 1, grid % loc_p_field
    r_rho_levels(j,k+1) = r_theta_levels(j,k) +    &
     ( r_theta_levels(j,k) - r_rho_levels(j,k-1) )
  End Do
End If

Return
End Subroutine Rcf_Gen_Height_Hybrid

End Module Rcf_generate_heights_mod
#endif
