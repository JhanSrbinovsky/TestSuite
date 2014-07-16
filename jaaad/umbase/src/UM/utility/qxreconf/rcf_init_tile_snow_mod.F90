#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initialisations for snow amount on tiles

Module Rcf_Init_Tile_Snow_Mod

! Subroutine Rcf_Init_Tile_Snow
!
! Description:
!   Initialises snow amount on tiles.
!
! Method:
!   Initialises snow amount on tiles from gridbox mean snow amount.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2   30/11/00   Original code.  R.A.Betts
!   5.3   04/06/01   Bug fix to ensure correct indexing of snow on
!                    land surface tiles
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Init_Tile_Snow( fields_out, field_count_out, hdr_out, &
                               snow_tile )

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_DecompTP_Mod, Only : &
    decomp_rcf_output

Use Rcf_Field_Type_Mod, Only : &
    Field_Type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_mean_snow,       &
    stashcode_prog_sec

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,                &
    PrStatus_Normal

Use Rcf_Parvars_Mod, Only : &
    mype

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field
Use Rcf_Lsm_Mod, Only : &
    local_lsm_out


Use Rcf_Field_Equals_Mod, Only : &
    Rcf_Field_Equals

Implicit None

! Arguments
Type( field_type ), Pointer          :: fields_out(:)
Type( field_type ), Intent( InOut )  :: snow_tile
Type( um_header_type ), Intent( In ) :: hdr_out

Integer, Intent( In )                :: field_count_out

! Local variables
Type( field_type ), Pointer          :: mean_snow
Integer                              :: size
Integer                              :: pos  ! field position
Integer                              :: i    ! looper
Real                                 :: mean_snow_land &
                                       ( snow_tile % level_size, 1)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,*) 'Initialising snow amount on tiles from gb mean snow'
End If

!----------------------------------------------------------------------
! Find snow amount in output fields and read it in
!----------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_mean_snow,            &
                 fields_out, field_count_out, pos)
mean_snow => fields_out(pos)
Call Rcf_Alloc_Field( mean_snow )
Call Rcf_Read_Field( mean_snow, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Copy global mean-snow values into land only variable
!----------------------------------------------------------------------
! DEPENDS ON: to_land_points
  Call To_Land_Points( mean_snow % Data,                  &
                       mean_snow_land,                    &
                       local_lsm_out,                     &
                       mean_snow % level_size,            &
                       size )

! Copy snow amount into tile snow amount (for all psuedo levels) and
! write it out
!----------------------------------------------------------------------
Do i = 1, snow_tile % levels
  snow_tile % Data(:,i) = mean_snow_land(:,1)
End Do

!----------------------------------------------------------------------
! Tidy Up
!----------------------------------------------------------------------
Call Rcf_Dealloc_Field( mean_snow )


End Subroutine Rcf_Init_Tile_Snow

End Module Rcf_Init_Tile_Snow_Mod
#endif
