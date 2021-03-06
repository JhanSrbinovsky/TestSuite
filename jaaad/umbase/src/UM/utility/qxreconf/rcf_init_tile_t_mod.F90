#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initialisations for tile t*, vegetation field

Module Rcf_Init_Tile_T_Mod

! Subroutine Rcf_Init_Tile_T
!
! Description:
!   Initialises surface temperature on tiles.
!
! Method:
!   Specifies TSTAR on each tile to equal the land mean TSTAR.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2   13/10/00   Original code.  P.Selwood
!   5.3      06/01   Add coastal tiling mods:
!                    Set TSTAR on each to that of the land mean
!                    value, rather than the gridbox mean.
!                    Also correct indexing of mean land surface
!                    temperature on land surface tiles
!                                                     Nic Gedney.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Init_Tile_T( fields_out, field_count_out, hdr_out, &
                            tstar_tile )

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_DecompTP_Mod, Only : &
    decomp_rcf_output

Use Rcf_Field_Type_Mod, Only : &
    Field_Type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_tstar,           &
    stashcode_tstar_land,      &
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

Use Rcf_CntlAtm_Mod, Only : &
    L_CTILE

Use Rcf_Lsm_Mod, Only : &
    local_lsm_out

Implicit None

! Arguments
Type( field_type ), Pointer          :: fields_out(:)
Type( field_type ), Intent( InOut )  :: tstar_tile
Type( um_header_type ), Intent( In ) :: hdr_out

Integer, Intent( In )                :: field_count_out

! Local variables
Type( field_type ), Pointer          :: tstar_land
Integer                              :: size
Integer                              :: pos  ! field position
Integer                              :: i    ! looper
Real                                 :: tstar_land_landpts  &
                                       (tstar_tile % level_size,1)


!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,*) 'Initialising T* on tiles from T*'
End If

!----------------------------------------------------------------------
! Find T* temperature in output fields and read it in
!----------------------------------------------------------------------
If (L_CTILE)THEN
  Call Rcf_Locate( stashcode_prog_sec, stashcode_tstar_land,         &
                   fields_out, field_count_out, pos)
Else
  Call Rcf_Locate( stashcode_prog_sec, stashcode_tstar,              &
                   fields_out, field_count_out, pos)
Endif
tstar_land => fields_out(pos)
Call Rcf_Alloc_Field( tstar_land )
Call Rcf_Read_Field( tstar_land, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Copy global land T* values into land only variable
!----------------------------------------------------------------------
! DEPENDS ON: to_land_points
  Call To_Land_Points( tstar_land % Data,                 &
                       tstar_land_landpts,                &
                       local_lsm_out,                     &
                       tstar_land % level_size,           &
                       size )

!----------------------------------------------------------------------
! Copy T* into tile T* (for all psuedo levels) and write it out
!----------------------------------------------------------------------
Do i = 1, tstar_tile % levels
  tstar_tile % Data(:,i) = tstar_land_landpts(:,1)
End Do

!----------------------------------------------------------------------
! Tidy Up
!----------------------------------------------------------------------
Call Rcf_Dealloc_Field( tstar_land )


End Subroutine Rcf_Init_Tile_T

End Module Rcf_Init_Tile_T_Mod
#endif
