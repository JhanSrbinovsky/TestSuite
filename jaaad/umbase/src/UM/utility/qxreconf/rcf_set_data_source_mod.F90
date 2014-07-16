#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Sets the data_source array corresponding to the fields array.

Module Rcf_Set_Data_Source_Mod

!  Subroutine Rcf_Set_Data_Source - sets the data_source array.
!
! Description:
!   Allocates the data_source array and sets the elements according
!   to how the data in the the associated fields array should be
!   initialised.
!
! Method:
!   Uses data from the ITEMS namelist.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Set_data_source( data_source, fields_in, fields_out, &
                                field_count_in, field_count_out,    &
                                hdr_in, hdr_out )

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Items_Mod             ! All of it

Use Rcf_Data_Source_Mod   ! Most of It

Use Rcf_Field_Type_Mod, Only : &
    Field_type

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_Parvars_Mod, Only : &
    mype

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Oper

Use Rcf_Submodel_Mod, Only : &
    Atmos_IM,            &
    Submodel_Ident

Use Rcf_Locate_mod, Only : &
    Rcf_Locate

Use Rcf_UMhead_Mod, Only : &
    UM_Header_type

Use Rcf_Recon_Mod, Only : &
    nice

Use Rcf_Stashcodes_Mod   ! So many of these are used the whole lot
                         ! have been included. Anything starting
                         ! with stashcode_ came from here.

Use Rcf_Recon_Mod, Only : &
    GRIB

Implicit None

! Arguments
Type( data_source_type ), Pointer  :: data_source(:)
Type( field_type), Pointer         :: fields_in(:)
Type( field_type ) , Pointer       :: fields_out(:)
Type( um_header_type ), Intent(In) :: hdr_in
Type( um_header_type ), Intent(In) :: hdr_out
Integer                            :: field_count_in
Integer                            :: field_count_out

! Comdecks
#include "cppxref.h"
#include "clookadd.h"

! Local variables
Integer                            :: item_index
Integer                            :: i
Integer                            :: j
Integer                            :: pos
Integer                            :: ErrorStatus
Integer                            :: pos_itc
Character (Len=*), Parameter       :: RoutineName='Rcf_Set_Data_Source'
Character (Len=80)                 :: Cmessage
Character (Len=*), Parameter       :: form="(2i5,4i4,' ',e9.3,' ',50a)"


!---------------------------------------------------------------
! check some of the arrays needed have been set up
!---------------------------------------------------------------

If ( Associated( data_source ) ) Then
  ErrorStatus = -10
  Cmessage = 'data_source is already set up!'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

If (.NOT. Associated( fields_out )  )Then
  ErrorSTatus = 20
  Cmessage = 'Fields have not yet been set up'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

!---------------------------------------------------------------
! Allocate the space
!---------------------------------------------------------------
Allocate( data_source( field_count_out ) )

!---------------------------------------------------------------
! Find any approriate item from namelist and set values
!----------------------------------------------------------------
Do i = 1, field_count_out
  item_index = 0
  Do j = 1, num_items
    If ( fields_out( i ) % stashmaster % item == item_array( j ) .AND.     &
         fields_out( i ) % stashmaster % section == sctn_array( j ) ) Then
      item_index = j
      Exit
    End If
  End Do

  If ( item_index > 0 ) Then       ! Have found the item - fill in gaps

    data_source( i ) % source      = source_array( item_index )
    data_source( i ) % Domain      = area_array  ( item_index )
    data_source( i ) % Ancil_SctnC = upas_array  ( item_index )
    data_source( i ) % Ancil_ItemC = upaa_array  ( item_index )
    data_source( i ) % RConst      = uprc_array  ( item_index )
    data_source( i ) % Ancil_File  = upaf_array  ( item_index )

  Else                             ! No item found - use defaults

    data_source( i ) % source      = Input_Dump
    data_source( i ) % Domain      = Whole_Grid
    data_source( i ) % Ancil_SctnC = 0
    data_source( i ) % Ancil_ItemC = 0
    data_source( i ) % RConst      = 0.0
    data_source( i ) % Ancil_File  = ' '

  End If

 ! Setting Source for items required from input dump but not
  ! available (and known about with relevant code written)
  Call Rcf_Locate( fields_out( i ) % stashmaster % section,         &
                   fields_out( i ) % stashmaster % item,            &
                   fields_in, field_count_in, pos, .TRUE.)

  If ( pos == 0 .AND. data_source( i ) % source == Input_Dump ) Then

  ! ----------------------------------
  ! This item is not in the input dump
  ! ----------------------------------

    Select Case( fields_out( i ) % stashmaster % model )

    Case ( Atmos_IM )                 ! Atmosphere Sub-Model

     ! ----------------
     ! Atmosphere items
     ! ----------------

      Select Case( fields_out( i ) % stashmaster % section )

      Case( stashcode_prog_sec )

        Select Case( fields_out( i ) % stashmaster % item )

        Case( stashcode_3d_cca,          & ! 3D Convective cloud
              stashcode_cca )              ! 2D Convective cloud

          data_source( i ) % source = Field_Calcs

        Case( stashcode_npp_pft_acc,     & ! Carbon accumulation fields
              stashcode_g_lf_pft_acc,    &
              stashcode_g_ph_lf_pft_acc, &
              stashcode_rsp_w_pft_acc,   &
              stashcode_rsp_s_acc,       &
              stashcode_catch_snow,      & ! NLT canopy snow capacity
              stashcode_catch_tile,      & ! Tiled canopy capacity
              stashcode_z0_tile          & ! Tiled roughness length
                                   )
          data_source( i ) % source = set_to_const
          data_source( i ) % rconst = -1.0

        Case( stashcode_can_water_tile,  & ! Vegetation fields
              stashcode_tstar_tile,      &
              stashcode_ice_temp_cat,    &
              stashcode_ice_conc_cat,    & ! ice conc (fraction) cats
              stashcode_ice_thick_cat,   & ! ice thickness categories
              stashcode_Gamtot,          & !\ large-scale hydrology
              stashcode_Zw,              & !/ fields
              stashcode_Sthzw,           & !/
              stashcode_Fsat,            & !\.
              stashcode_Fwetl,           & !/.
              stashcode_a_fsat,         & !\.
              stashcode_c_fsat,         & !/.
              stashcode_a_fwet,         & !\.
              stashcode_c_fwet,         & !/.
              stashcode_snow_tile)

          data_source( i ) % source = Field_Calcs

        Case( stashcode_infil_max_tile,  &
              stashcode_snow_grnd,       & ! Snow beneath NLT canopy
              stashcode_snow_on_ice,     &
              stashcode_ice_snow_depth,  &
              stashcode_sw_down_tile,    &
              stashcode_sw_down,         &
              stashcode_lw_up_diff,      &
              stashcode_mmr_smoke_fr,    &
              stashcode_mmr_smoke_ag,    &
              stashcode_mmr_smoke_cl,    &
              stashcode_biom_surf_em,    &
              stashcode_biom_elev_em,    &
              stashcode_dust1_mmr,       &
              stashcode_dust2_mmr,       &
              stashcode_dust3_mmr,       &
              stashcode_dust4_mmr,       &
              stashcode_dust5_mmr,       &
              stashcode_dust6_mmr,       &
              stashcode_albedo_sice,     &
              stashcode_albedo_land,     &
              stashcode_soilcarb_dpm,    & ! RothC soil C prognostic
              stashcode_soilcarb_rpm,    & ! RothC soil C prognostic
              stashcode_soilcarb_bio,    & ! RothC soil C prognostic
              stashcode_soilcarb_hum,    & ! RothC soil C prognostic
              stashcode_qcf2,            & ! second cloud ice prognostic
              stashcode_qrain,           & ! prognostic rain
              stashcode_qgraup,          & ! prognostic graupel
              stashcode_lcbase,          & ! lowest convective cloud base
              stashcode_3d_ccw           & ! CCW profile sent to radiation
                                        )

          data_source( i ) % source = set_to_zero

        Case( stashcode_rgrain )

          data_source( i ) % source = set_to_const
          data_source( i ) % rconst = 50.0

        Case( stashcode_tstar_land,     &
              stashcode_tstar_sea,      &
              stashcode_tstar_sice )

          data_source( i ) % source = Already_Processed

        End Select                       ! based on STASH item code

        !kdcorbin, 05/10 - added boundary layer section
        Case (stashcode_bl_sec)
              data_source(i) % source = set_to_zero
 
     End Select                         ! based on STASH section code

      !For files not found when reading from GRIB data
      If ( GRIB ) Then

        Select Case( fields_out( i ) % stashmaster % section )

        Case( stashcode_prog_sec )

          Select Case( fields_out( i ) % stashmaster % item )

          Case( stashcode_sea_ice_temp,     &
                stashcode_u_adv,            &
                stashcode_v_adv)

            data_source( i ) % source = Field_Calcs

          Case( stashcode_w_adv,            &
                stashcode_qcf,              &
                stashcode_cc_lwp,           &
                stashcode_unfrozen_soil,    &
                stashcode_frozen_soil,      &
                stashcode_qcl,              &
                stashcode_n_turb_mixlvs,    &
                stashcode_lvl_bse_dp_sc,    &
                stashcode_lvl_top_dp_sc,    &
                stashcode_bl_conv_flag,     &
                stashcode_turb_temp,        &
                stashcode_turb_humid,       &
                stashcode_area_cf,          &
                stashcode_bulk_cf,          &
                stashcode_liquid_cf,        &
                stashcode_frozen_cf,        &
                stashcode_sfc_zonal_cur,    &
                stashcode_sfc_merid_cur,    &
                stashcode_3d_cca,           & ! has to overide source=8
                stashcode_can_water_tile,   &
                stashcode_rho,              & ! rho calc'd after interp
                stashcode_exner,            &
                stashcode_can_conduct,      &
                stashcode_ccb,              &
                stashcode_cct,              &
                stashcode_mean_canopyw,     &
                stashcode_surf_z_curr,      &
                stashcode_surf_m_curr,      &
                stashcode_w)

            data_source( i ) % source = Set_To_Zero

          Case( stashcode_bl_depth )

            data_source( i ) % source = set_to_const
            data_source( i ) % rconst = 500.000

          Case( stashcode_z0 )

            data_source( i ) % source = set_to_const
            data_source( i ) % rconst = 0.500

          End Select                 ! select by item code
        End Select                   ! select by section code

      End If                       ! If GRIB section

    Case Default                 ! Couldn't find a proper Internal Model

      ! The Model type didn't match one of the three defined types.
      Write (Cmessage,*) " Couldn't find a Sub-Model ID type for : ", &
                  " Model ", fields_out( i ) % stashmaster % model,   &
                  " Section ", fields_out( i ) % stashmaster % model, &
                  " Item ",  fields_out( i ) % stashmaster % item
      ErrorStatus = 25
      Call Ereport( RoutineName, ErrorStatus, Cmessage )

    End Select                   ! select by Internal model

    ! Check that Source is now set correctly otherwise, fail
    If ( data_source( i ) % source == Input_Dump ) Then
      Write ( Cmessage, *) 'Section ',                              &
                           fields_out( i ) % stashmaster % section, &
                           'Item ',                                 &
                           fields_out( i ) % stashmaster % item ,   &
                           ' : Required field is not in input dump!'
      ErrorStatus = 30
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

  End If  ! If (pos == 0 ..)

  ! ---------------------------------------------------------------
  ! Some fields, even when found in the input dump, still need
  ! extra work. Mostly for fields that may have differing numbers
  ! of pseudo levels. But also fields needing to be recalculated to
  ! ensure consistancy with others.
  ! ---------------------------------------------------------------

  If ( pos /= 0 .AND. data_source( i ) % source == Input_Dump ) Then

  ! ----------------------------------
  ! This item is in the input dump
  ! ----------------------------------

    Select Case( fields_out( i ) % stashmaster % model )

    ! -----------
    ! Atmos items
    ! -----------

    Case ( Atmos_IM )
      
      Select Case( fields_out( i ) % stashmaster % section )

      Case( stashcode_prog_sec )

          Select Case( fields_out( i ) % stashmaster % item )

          Case( stashcode_icethick)

            If (h_int_active) Then
              data_source( i ) % source = Field_Calcs
            End If

          Case( stashcode_ice_temp_cat,       &
                stashcode_ice_conc_cat,       &
                stashcode_ice_thick_cat )

          !-------------------------------------------------------------
          ! Check number of pseudo levels in input dump is same as
          ! number of sea ice categories (nice) in output dump.
          ! If not -recalculate
          !-------------------------------------------------------------

            Call Rcf_Locate( fields_out( i ) % stashmaster % section, &
                             fields_out( i ) % stashmaster % item,    &
                             fields_in, field_count_in, pos_itc )

            If (fields_in(pos_itc) % levels /= nice) Then
               data_source( i ) % source = Field_Calcs
               If ( mype == 0 .AND.                                   &
                   Printstatus >= PrStatus_Normal ) Then
                Write (6,*) "Setting source of stashcode item ",      &
                  fields_out( i ) % stashmaster % item,               &
                  " to calculations due to difference in NICE"
              EndIf
            EndIf

          Case( stashcode_ice_snow_depth )

          !-------------------------------------------------------------
          ! Check number of pseudo levels in input dump is same as
          ! number of sea ice categories (nice) in output dump.
          ! If not -set to zero
          !-------------------------------------------------------------

            Call Rcf_Locate( fields_out( i ) % stashmaster % section, &
                             fields_out( i ) % stashmaster % item,    &
                             fields_in, field_count_in, pos_itc )

            If (fields_in(pos_itc) % levels /= nice) Then
               data_source( i ) % source = set_to_zero
               If ( mype == 0 .AND.                                   &
                   Printstatus >= PrStatus_Normal ) Then
                Write (6,*) "Setting source of stashcode item ",      &
                  fields_out( i ) % stashmaster % item,               &
                  " to -set-to-zero- due to difference in NICE"
              EndIf
            EndIf

          Case( stashcode_infil_max_tile, &
                stashcode_sw_down_tile )

          !-------------------------------------------------------------
          ! Check number of tiles in input dump is same as
          ! number of tiles in output dump.
          ! If not -set to zero
          !-------------------------------------------------------------

            Call Rcf_Locate( fields_out( i ) % stashmaster % section, &
                             fields_out( i ) % stashmaster % item,    &
                             fields_in, field_count_in, pos_itc )

            If (fields_in(pos_itc) % levels /=                        &
                fields_out( i )    % levels ) Then
               data_source( i ) % source = set_to_zero
               If ( mype == 0 .AND.                                   &
                   Printstatus >= PrStatus_Normal ) Then
                Write (6,*) "Setting source of stashcode item ",      &
                  fields_out( i ) % stashmaster % item,               &
                  " to -set-to-zero- due to difference in no. of tiles"
              EndIf
            EndIf

          Case( stashcode_catch_tile,     &
                stashcode_z0_tile )

          !-------------------------------------------------------------
          ! Check number of tiles in input dump is same as
          ! number of tiles in output dump.
          ! If not -set to user constant -1
          !-------------------------------------------------------------

            Call Rcf_Locate( fields_out( i ) % stashmaster % section, &
                             fields_out( i ) % stashmaster % item,    &
                             fields_in, field_count_in, pos_itc )

            If (fields_in(pos_itc) % levels /=                        &
                fields_out( i )    % levels ) Then
               data_source( i ) % source = set_to_const
               data_source( i ) % rconst = -1.000
               If ( mype == 0 .AND.                                   &
                   Printstatus >= PrStatus_Normal ) Then
                Write (6,*) "Setting source of stashcode item ",      &
                  fields_out( i ) % stashmaster % item,               &
                  " to -1 due to difference in no. of tiles"
              EndIf
            EndIf

          Case( stashcode_can_water_tile, &
                stashcode_tstar_tile,     &
                stashcode_snow_tile)

          !-------------------------------------------------------------
          ! Check number of tiles in input dump is same as
          ! number of tiles in output dump.
          ! If not -set to Field Calcs
          !-------------------------------------------------------------

            Call Rcf_Locate( fields_out( i ) % stashmaster % section, &
                             fields_out( i ) % stashmaster % item,    &
                             fields_in, field_count_in, pos_itc )

            If (fields_in(pos_itc) % levels /=                        &
                fields_out( i )    % levels ) Then
               data_source( i ) % source = Field_Calcs
               If ( mype == 0 .AND.                                   &
                   Printstatus >= PrStatus_Normal ) Then
                Write (6,*) "Setting source of stashcode item ",      &
                  fields_out( i ) % stashmaster % item,               &
                  " to -FielcCalcs- due to difference in no. of tiles"
              EndIf
            EndIf
        End Select                        ! Select on item number
      End Select                          ! Select on section number
    
    Case Default                 ! Couldn't find a proper Internal Model

      ! The Model type didn't match one of the three defined types.
      Write (Cmessage,*) " Couldn't find a Sub-Model ID type for : ",   &
                  " Model ", fields_out( i ) % stashmaster % model,     &
                  " Section ", fields_out( i ) % stashmaster % section, &
                  " Item ",  fields_out( i ) % stashmaster % item
      ErrorStatus = 70
      Call Ereport( RoutineName, ErrorStatus, Cmessage )

    End Select                            ! Select on Internal model

  End If                                  ! Item was in Dump

  ! Check on Rimwidtha if an atmos LBC file for copying
  If ( data_source( i ) % source == Input_Dump .AND.         &
     (fields_out( i ) % stashmaster % grid_type ==ppx_atm_lbc_theta.OR.&
      fields_out( i ) % stashmaster % grid_type == ppx_atm_lbc_u   .OR.&
      fields_out( i ) % stashmaster % grid_type == ppx_atm_lbc_v) ) Then

    If ( hdr_in  % Lookup( lbrow, fields_out( i ) % dump_pos ) /=  &
         hdr_out % Lookup( lbrow, fields_out( i ) % dump_pos) ) Then
      Cmessage = 'Rimwidth needs to be equal for input and output grids&
                 &if LBCs are copied'
      ErrorStatus = 40
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If
  End If

End Do

!------------------------------------------------------------------
! Print out some diagnostics if required....
!------------------------------------------------------------------
If (mype == 0 .AND. PrintStatus >= PrStatus_Oper ) Then
  Write (6,*)
  Write (6,*) 'Data source'

  Do i = 1, field_count_out
    Write (6,form) fields_out( i ) % stashmaster % section,            &
                   fields_out( i ) % stashmaster % item,               &
                   data_source( i ) % source,                          &
                   data_source( i ) % Domain,                          &
                   data_source( i ) % Ancil_SctnC,                     &
                   data_source( i ) % Ancil_ItemC,                     &
                   data_source( i ) % RConst,                          &
   data_source(i) % Ancil_File(1:Len_Trim(data_source(i) % Ancil_File))
  End Do
  Write (6,*)

End If


Return
End Subroutine Rcf_Set_Data_Source
End Module Rcf_Set_Data_Source_Mod
#endif
