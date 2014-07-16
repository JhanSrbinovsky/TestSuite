
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Sets up horizontal interpolation weights

Module Rcf_init_h_interp_mod

!  Subroutine Rcf_init_h_interp - initialises horizontal interp. weights
!
! Description:
!   Sets up interpolation weights based on choice of scheme.
!
! Method:
!   Allocates space for and sets up weights
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Init_H_Interp( grid_in, grid_out, hdr_in, hdr_out )

Use Rcf_Interp_Weights_Mod  ! Almost all of this used.

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Use Rcf_UMhead_Mod, Only :    &
    um_header_type

Use Ereport_mod, Only : &
    Ereport

Use Rcf_Lsm_Mod, Only : &
    cyclic

Use Rcf_Recon_Mod, Only : &
    GRIB

Use rcf_h_int_init_bl_mod, Only : &
    rcf_h_int_init_bl

Use Rcf_PrintStatus_Mod, Only : &
    LTimer

Implicit None

! Arguments
Type (grid_type), Intent(In)      :: grid_in
Type (grid_type), Intent(In)      :: grid_out
Type (um_header_type), Intent(In) :: hdr_in
Type (um_header_type), Intent(In) :: hdr_out

! Local data
Integer                       :: ErrorStatus
Character (Len=*), Parameter  :: RoutineName = 'Init_H_Interp'
Character (Len=80)            :: Cmessage

External Timer

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 3)

! Set the cyclic flag for coastal adjustment
If ( Hdr_Out % Fixhd(4) /= 3 .AND. Hdr_Out % FixHd(4) /= 103 ) Then
  CYCLIC=.TRUE.
Else
  CYCLIC=.FALSE.
End If

!----------------------------------------------------------------
! Test if we know which interpolation scheme we are doing!
!----------------------------------------------------------------
If ( h_int_method == unset ) Then
  Cmessage = 'No interpolation method specified - cannot set weights'
  ErrorStatus = 10
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If


Select Case ( h_int_method )
  Case ( bilinear )
    Allocate( bl_index_b_l( grid_out % glob_p_field, idim ) )
    Allocate( bl_index_b_r( grid_out % glob_p_field, idim ) )
    Allocate( bl_index_nearest( grid_out % glob_p_field ) )
    Allocate( coeff1( grid_out % glob_p_field ) )
    Allocate( coeff2( grid_out % glob_p_field ) )
    Allocate( coeff3( grid_in % glob_p_field ) )
    Allocate( coeff4( grid_in % glob_p_field ) )
    Allocate( weight_t_r( grid_out % glob_p_field, idim ) )
    Allocate( weight_b_r( grid_out % glob_p_field, idim ) )
    Allocate( weight_t_l( grid_out % glob_p_field, idim ) )
    Allocate( weight_b_l( grid_out % glob_p_field, idim ) )

    Call rcf_h_int_init_bl( grid_in, grid_out, hdr_in, hdr_out, GRIB)

  Case ( area_weighted )
    Allocate( aw_index_targ_lhs( grid_out % glob_p_row_length+1, idim ))
    Allocate( aw_index_targ_top( grid_out % glob_p_rows + 1, idim ))
    Allocate( aw_colat_t( grid_out % glob_p_rows + 1, idim ))
    Allocate( aw_long_l( grid_out % glob_p_row_length+1, idim ))
    Allocate( bl_index_b_l( grid_out % glob_p_field, idim ) )
    Allocate( bl_index_b_r( grid_out % glob_p_field, idim ) )
    Allocate( bl_index_nearest( grid_out % glob_p_field ) )
    Allocate( weight_t_r( grid_out % glob_p_field, idim ) )
    Allocate( weight_b_r( grid_out % glob_p_field, idim ) )
    Allocate( weight_t_l( grid_out % glob_p_field, idim ) )
    Allocate( weight_b_l( grid_out % glob_p_field, idim ) )

! DEPENDS ON: h_int_init_aw
    Call h_int_init_aw( icof, idim, grid_out % glob_p_field,         &
                 grid_in % glob_p_rows,                              &
                 grid_out % glob_p_rows, grid_in % glob_p_row_length,&
                 grid_out % glob_p_row_length, grid_in % glob_u_field,&
                 grid_out % glob_u_field,                            &
                 grid_in % glob_u_rows, grid_out % glob_u_rows,      &
                 grid_out % global, GRIB, hdr_in % FixHd,            &
                 hdr_out % FixHd, hdr_in % RealC, hdr_out % RealC,   &
                 aw_area_box, aw_index_targ_lhs, aw_index_targ_top,  &
                 bl_index_b_l, bl_index_b_r, bl_index_nearest,       &
                 aw_colat_t, aw_long_l, weight_t_r, weight_b_r,      &
                 weight_t_l, weight_b_l )



  Case Default
    Cmessage = 'Interpolation method not supported'
    ErrorStatus = 20
    Call Ereport( RoutineName, ErrorStatus, Cmessage )

End Select

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 4)

Return
End Subroutine Rcf_Init_H_Interp

End Module Rcf_init_h_interp_mod
