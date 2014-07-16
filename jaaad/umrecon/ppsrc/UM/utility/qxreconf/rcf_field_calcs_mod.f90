
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  Performs Source=8 field initialisation calculations

Module Rcf_Field_Calcs_Mod

!  Subroutine Rcf_Field_Calcs - field initialisation calculations.
!
! Description:
!   Some fields may have hard-coded Source=8 initialisation (certain
!   slab fields may be set by UMUI in this way also). These fields
!   are initialised by the calculations in this routine.
!
! Method:
!   Choice of method applied determined by stashcode.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Field_Calcs( fields_in, fields_out, field_count_in, &
                            field_count_out, data_source, hdr_in,  &
                            hdr_out )

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use rcf_data_source_Mod, Only : &
    data_source_type,           &
    Field_Calcs

Use Rcf_Submodel_Mod, Only : &
    Submodel_Ident,      &
    Atmos_IM

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_DecompTP_Mod, Only : &
    decomp_rcf_output

Use Rcf_Alloc_Field_mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

Use Rcf_Write_Field_Mod, Only : &
    Rcf_Write_Field

Use Rcf_Derv_Ice_Temp_Mod, Only : &
    Rcf_Derv_Ice_Temp

Use Rcf_Derv_3D_CCA_Mod, Only : &
    Rcf_Derv_3D_CCA

Use Rcf_Derv_2D_CCA_Mod, Only : &
    Rcf_Derv_2D_CCA

Use Rcf_Derv_Ice_Thick_Mod, Only : &
    Rcf_Derv_Ice_Thick

Use Rcf_Derv_Sea_Ice_Temp_Mod, Only :&
    Rcf_Derv_Sea_Ice_Temp

Use Rcf_PrintStatus_Mod, Only : &
    LTimer

Use Rcf_Init_Canopy_Water_Mod, Only : &
    Rcf_Init_Canopy_Water

Use Rcf_Init_Tile_Snow_Mod, Only : &
    Rcf_Init_Tile_Snow

Use Rcf_Init_Tile_T_Mod, Only : &
    Rcf_Init_Tile_T

Use Rcf_Calc_Gamtot_Mod, Only : &
    Rcf_Calc_Gamtot

Use Rcf_Est_Zw_Mod, Only :      &
    Rcf_Est_Zw

Use Rcf_Est_Sthzw_Mod, Only :      &
    Rcf_Est_Sthzw

Use Rcf_Fit_Fsat_Mod, Only : &
    Rcf_Fit_Fsat

Use Rcf_Calc_Fsat_Mod, Only :    &
    Rcf_Calc_Fsat

Use Rcf_Derv_Ice_Cat_Thick_Mod, Only :    &
    Rcf_Derv_Ice_Cat_Thick

Use Rcf_Derv_Adv_Winds_Mod, Only : &
    Rcf_Derv_Adv_Winds

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_3d_cca,                                     &
    stashcode_ice_conc_cat,    stashcode_ice_thick_cat,   &
    stashcode_ice_temp_cat,    stashcode_ice_snow_depth,  &
    stashcode_icethick,                                   &
    stashcode_cca,             stashcode_can_water_tile,  &
    stashcode_sea_ice_temp,                               &
    stashcode_snow_tile,       stashcode_tstar_tile,      &
    stashcode_gamtot,          stashcode_zw,              &
    stashcode_sthzw,           stashcode_fsat,            &
    stashcode_fwetl,                                      &
    stashcode_a_fsat,          stashcode_c_fsat,          &
    stashcode_a_fwet,          stashcode_c_fwet,          &
    stashcode_prog_sec,                                   &
    stashcode_u_adv,           stashcode_v_adv

Implicit None

! Arguments
Integer,            Intent(In)       :: field_count_in
Integer,            Intent(In)       :: field_count_out
Type( field_type ), Pointer          :: fields_in( : )
Type( field_type ), Pointer          :: fields_out( : )
Type( data_source_type ), Pointer    :: data_source( : )
Type( um_header_type ), Intent(In)   :: hdr_in
Type( um_header_type ), Intent(In)   :: hdr_out
Real , Allocatable                   :: a_fsat(:,:)
Real , Allocatable                   :: c_fsat(:,:)
Real , Allocatable                   :: a_fwet(:,:)
Real , Allocatable                   :: c_fwet(:,:)

! Local variables
Integer                           :: i
Integer                           :: ErrorStatus
Character (Len=*), Parameter      :: RoutineName='Rcf_Field_Calcs'
Character (Len=80)                :: Cmessage
Logical :: L_Calc_Fit_Fsat   !If true then call Rcf_Fit_Fsat


External Timer

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 3)

L_Calc_Fit_Fsat = .true.   ! initialise logical to true

!-----------------------------------------------------------------
! Loop around output fields/sources until a source=field_calcs
! is found
!-----------------------------------------------------------------
Do i = 1, field_count_out
  If ( data_source( i ) % Source == Field_Calcs ) Then

    Call Rcf_Alloc_Field( fields_out( i ) )

    Select Case( fields_out( i ) % stashmaster % model )

    Case ( Atmos_IM )

      ! -----------
      ! Atmos items
      ! -----------

      Select Case( fields_out( i ) % stashmaster % section )

      ! For the moment only section zero fields can be used here.
      Case ( stashcode_prog_sec )

        Select Case( fields_out(i) % stashmaster % item )

        Case( stashcode_3d_cca )
          Call Rcf_Derv_3D_CCA( fields_in, fields_out, field_count_in, &
                                field_count_out, hdr_in, hdr_out ,     &
                                fields_out( i ) )

        Case( stashcode_cca )
          Call Rcf_Derv_2D_CCA( fields_in, fields_out, field_count_in, &
                                field_count_out, hdr_in, hdr_out ,     &
                                fields_out( i ) )

        Case( stashcode_icethick )
          Call Rcf_Derv_Ice_Thick( fields_out, field_count_out,        &
                                   hdr_out, fields_out(i) )

        Case( stashcode_sea_ice_temp )
          Call Rcf_Derv_Sea_Ice_Temp( fields_out, field_count_out,  &
                                      fields_out( i ), hdr_out )

        Case( stashcode_can_water_tile )
          Call Rcf_Init_Canopy_Water( fields_in, field_count_in,       &
                                      hdr_in, fields_out( i ) )

        Case( stashcode_snow_tile )
          Call Rcf_Init_Tile_Snow( fields_out, field_count_out,        &
                                   hdr_out, fields_out( i ) )

        Case( stashcode_ice_temp_cat )
          Call Rcf_Derv_Ice_Temp( fields_out,field_count_out, hdr_out, &
                                 fields_out( i) )

        Case( stashcode_tstar_tile )
          Call Rcf_Init_Tile_T( fields_out, field_count_out, hdr_out,  &
                                fields_out( i ) )

        Case( stashcode_Gamtot )
          Call Rcf_Calc_Gamtot( fields_out, field_count_out,  &
                                   fields_out( i ), hdr_out )
! If calculation required for the LSH fitting variables,
! must first allocate space for them.

        Case( stashcode_a_fsat, stashcode_c_fsat, &
               stashcode_a_fwet, stashcode_c_fwet)

          If(L_Calc_Fit_Fsat)Then

            Allocate (a_fsat &
             (fields_out(i) % level_size, fields_out(i) % levels ) )
            Allocate (c_fsat &
             (fields_out(i) % level_size, fields_out(i) % levels ) )
            Allocate (a_fwet &
             (fields_out(i) % level_size, fields_out(i) % levels ) )
            Allocate (c_fwet &
             (fields_out(i) % level_size, fields_out(i) % levels ) )

            Call Rcf_Fit_Fsat ( fields_out, field_count_out,      &
               fields_out(i),a_fsat,c_fsat,                       &
               a_fwet,c_fwet,hdr_out )

! set to false so Rcf_Fit_Fsat is not called again
            L_Calc_Fit_Fsat = .false.

          End If

          If(.not.L_Calc_Fit_Fsat)Then

            Select Case( fields_out( i ) % stashmaster % item )

              Case( stashcode_a_fsat)
                fields_out(i) % Data (:,:) = a_fsat (:,:)

              Case( stashcode_c_fsat )
                fields_out(i) % Data (:,:) = c_fsat (:,:)

              Case( stashcode_a_fwet)
                fields_out(i) % Data (:,:) = a_fwet (:,:)

              Case( stashcode_c_fwet )
                fields_out(i) % Data (:,:) = c_fwet (:,:)

            End Select

          End If

        Case( stashcode_Sthzw )
           Call Rcf_Est_Sthzw( fields_out, field_count_out,  &
                                   fields_out( i ), hdr_out )

        Case( stashcode_Zw )
          Call Rcf_Est_Zw( fields_out, field_count_out,  &
                                   fields_out( i ), hdr_out )

        Case( stashcode_Fsat )
          Call Rcf_Calc_Fsat( fields_out, field_count_out,  &
                     stashcode_Fsat, fields_out( i ), hdr_out )

        Case( stashcode_Fwetl )
          Call Rcf_Calc_Fsat( fields_out, field_count_out,  &
                     stashcode_Fwetl, fields_out( i ), hdr_out )

        Case( stashcode_ice_conc_cat, stashcode_ice_thick_cat )
          Call Rcf_Derv_Ice_Cat_Thick (                               &
                                fields_out( i ) % stashmaster % item, &
                                fields_out, field_count_out,          &
                                hdr_out, fields_out( i ) )

        Case( stashcode_u_adv, stashcode_v_adv )
          Call Rcf_Derv_Adv_Winds(                                   &
                               fields_out( i ) % stashmaster % item, &
                               fields_out, field_count_out,          &
                               hdr_out, fields_out( i ) )

        Case Default
          ErrorStatus = 30
          Write (Cmessage,*) 'No Field Calculations specified for &
          &section ', fields_out( i ) % stashmaster % section,    &
          &'item', fields_out( i ) % stashmaster % item
          Call Ereport( RoutineName, ErrorStatus, Cmessage )

        End Select

      Case Default
        ErrorStatus = 30
        Write (Cmessage,*) 'No Field Calculations specified for &
        &section ', fields_out( i ) % stashmaster % section
        Call Ereport( RoutineName, ErrorStatus, Cmessage )

      End Select

    Case Default                          ! Selection on Internal Model
      ErrorStatus = 40
      Write (Cmessage,*) 'No Field Calculations specified for &
                    &Submodel ', fields_out( i ) % stashmaster % model
      Call Ereport( RoutineName, ErrorStatus, Cmessage )

    End Select                            ! Selection on Internal Model
!-----------------------------------------------------------------
! Clean up and write out
!-----------------------------------------------------------------
    Call Rcf_Write_Field( fields_out( i ), Hdr_out, decomp_rcf_output )
    Call Rcf_DeAlloc_Field( fields_out( i ) )
  End If

End Do

! Deallocate arrays if they have been allocated
If (Allocated(a_fsat)) Deallocate(a_fsat)
If (Allocated(c_fsat)) Deallocate(c_fsat)
If (Allocated(a_fwet)) Deallocate(a_fwet)
If (Allocated(c_fwet)) Deallocate(c_fwet)


! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 4)

Return
End Subroutine Rcf_Field_Calcs
End Module Rcf_Field_Calcs_Mod
