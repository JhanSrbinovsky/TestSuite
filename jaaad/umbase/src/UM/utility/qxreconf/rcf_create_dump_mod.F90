#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Dump creation `main loop'

Module Rcf_Create_Dump_Mod

!  Subroutine Rcf_Create_Dump - main dump creation loop.
!
! Description:
!   Creates the output dump based on the input choices.
!
! Method:
!   Data sources are setup and then fields in the output dump are
!   looped over. Each one is set appropriately ( interpolated from
!   input dump, set to 0/MDI/constant/from file/etc). Pre and Post
!   processing is performed as required.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Create_Dump( Hdr_In, Hdr_Out, fields_in, fields_out, &
                            field_count_in, field_count_out )

Use Rcf_Set_Orography_Mod, Only : &
    Rcf_Set_Orography

Use Rcf_Pre_Interp_Transform_Mod, Only : &
    Rcf_Pre_Interp_Transform

Use Rcf_Post_Interp_Transform_Mod, Only : &
    Rcf_Post_Interp_Transform

Use Rcf_Post_Process_Mod, Only :&
    Rcf_Post_Process_Atmos

Use Rcf_Recon_Mod, Only : &
    Trans,            &
    Uars,             &
    GRIB

Use Rcf_FortranIO_Mod, Only : &
    Rcf_Get_Unit,             &
    Rcf_Free_Unit

Use Rcf_FreeUMhdr_mod, Only : &
    Rcf_FreeUMhdr

Use Rcf_UMhead_Mod, Only : &
    UM_header_type

Use Ereport_mod, Only : &
    Ereport

Use Rcf_Data_Source_Mod    ! All of it...

Use rcf_read_field_mod, Only : &
    Rcf_Read_Field

Use rcf_write_field_mod, Only : &
    Rcf_Write_Field

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

Use Rcf_DecompTP_Mod, Only : &
    Decomp_rcf_input,    &
    Decomp_rcf_output

Use Rcf_Grid_Type_Mod, Only : &
    Input_Grid,           &
    Output_Grid

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Alloc_Field_mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

Use Rcf_Locate_mod, Only : &
    Rcf_Locate

Use Rcf_Parvars_Mod, Only : &
    mype

Use Rcf_field_equals_mod, Only : &
    Rcf_field_equals

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_V_Int_Ctl_Mod, Only : &
    v_int_active,         &
    v_int_active_soil

Use rcf_interpolate_mod, Only : &
    rcf_interpolate

Use rcf_set_data_source_mod, Only : &
    Rcf_Set_Data_Source

Use Rcf_Lsm_Mod, Only : &
    local_lsm_out

Use Rcf_aux_file_mod, Only : &
    tracers,                  uars_data,      &
    user_prog,                transplant,     &
    rcf_aux_file

Use Rcf_Field_Calcs_Mod, Only : &
    Rcf_Field_Calcs

Use Rcf_Set_Interp_Flags_Mod, Only :               &
    Rcf_Set_Interp_Flags,                          &
    interp_v_only,            interp_h_only,       &
    interp_all,               interp_no_op,        &
    interp_copy

Use Rcf_Rotate_Mod, Only : &
    Rcf_Rotate,            &
    ToStandard,            &
    FromStandard

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_lsm,             &
    stashcode_orog,            &
    stashcode_icefrac,         &
    stashcode_tstar,           &
    stashcode_tstar_land,      &
    stashcode_tstar_sea,       &
    stashcode_tstar_sice,      &
    stashcode_prog_sec

Use Rcf_GRIB_Interp_TnPstar_Mod, Only : &
    Rcf_GRIB_Interp_TnPstar

Implicit None

! Arguments
Integer, Intent(In)                       :: field_count_in
Integer, Intent(In)                       :: field_count_out
Type (UM_Header_Type), Intent(InOut)      :: Hdr_In
Type (UM_Header_Type), Intent(InOut)      :: Hdr_Out
Type (field_type), Pointer                :: fields_in(:)
Type (field_type), Pointer                :: fields_out(:)

! Comdecks
#include "c_mdi.h"
#include "cppxref.h"

! Local vars
Integer                          :: pos      ! position in fields array
Integer                          :: i        ! looper
Integer                          :: err      ! for file open/close
Integer                          :: IOStatus
Integer                          :: ErrorStatus
Logical                          :: l_exist
Character (Len=80)               :: Cmessage
Character (Len=Max_Filename_Len) :: TracerFileName
Character (Len=*), Parameter     :: RoutineName='Rcf_Create_Dump'

Type (field_type), Target  :: interp_orog   ! interpolated input orog
Type (field_type), Pointer :: orog_out      ! ptr to output orog
Type (field_type), Pointer :: orog_in       ! ptr to intput orog
Type (field_type), Pointer :: input_field   ! ptr to input field
Type (field_type), Pointer :: output_field  ! ptr to output field
Type (UM_Header_Type)      :: Hdr_Aux       ! auxillary file header

Type (data_source_type), Pointer :: data_source(:)

! Formatting
Character (Len=*), Parameter :: form="(a20, i3, ' ( Section', i4, ' )', &
                                      &' ( Stashcode', i4, ' )', ' ', a36)"

External File_Open, File_Close, Fort_Get_Env
!---------------------------------------------------------------
! Initialise some data
!---------------------------------------------------------------
Nullify( interp_orog % data )
Nullify( interp_orog % data_int )
Nullify( interp_orog % data_log )
Nullify( data_source )

! set the source of the data for output dump
Call Rcf_Set_Data_Source( data_source, fields_in, fields_out,  &
                          field_count_in, field_count_out,     &
                          hdr_in, hdr_out )

!-----------------------------------------------------------
! Set up the input, output and interpolated orography fields
!-----------------------------------------------------------
Call Rcf_Set_Orography( fields_in, fields_out, field_count_in,      &
                        field_count_out, hdr_in, hdr_out,           &
                        data_source, orog_in, orog_out, interp_orog )

! If recon from GRIB T & Pstar need to interp'd horizontally
! for height generation
If ( GRIB ) Then
  Call Rcf_GRIB_Interp_TnPstar( fields_in, fields_out, Hdr_In,      &
                                  field_count_in, field_count_out)
End If

!-------------------------------------------------------------
! Land Sea Mask - need to deal with this a little differently
! as data already read in and processed.
!-------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_lsm,                  &
                 fields_out, field_count_out, pos, .TRUE.)
If ( pos /= 0 ) Then
  Call Rcf_Alloc_Field( fields_out ( pos ) )
  fields_out( pos ) % Data_Log(:,1) = local_lsm_out(:)

  Call Rcf_Write_Field( fields_out( pos ), Hdr_Out, decomp_rcf_output )
  Call Rcf_DeAlloc_Field( fields_out( pos ) )
End If

!-------------------------------------------------------------
! Setup interpolation flags
!-------------------------------------------------------------
Call Rcf_Set_Interp_Flags( fields_in, fields_out, field_count_in, &
                           field_count_out, data_source )

!------------------------------------------------------------
! Rotate input winds if so required
!------------------------------------------------------------
If ( Input_Grid % Rotated .AND. h_int_active ) Then
  If (mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
    Write (6,*) 'Rotating input winds'
  End If

  Call Rcf_Rotate( fields_in, field_count_in, Input_Grid, Hdr_In,  &
                   decomp_rcf_input, ToStandard)
End If

!-----------------------------------------------------
! Main loop
!-----------------------------------------------------
Do i = 1, field_count_out ! run through the output fields

  !--------------------------------------------
  ! Which fields have we already dealt with?
  ! Orography
  ! LSM
  ! T*             )  if source == 4
  ! Ice fraction   )  if source == 4
  !--------------------------------------------
  If ( fields_out( i ) % stashmaster % section == stashcode_prog_sec .AND. &
       ( fields_out( i ) % stashmaster % item == stashcode_orog .OR.       &
         fields_out( i ) % stashmaster % item == stashcode_lsm  .OR.       &
         (fields_out( i ) % stashmaster % item == stashcode_tstar .AND.    &
          data_source( i ) % source == Set_to_MDI )              .OR.      &
         (fields_out( i ) % stashmaster % item == stashcode_icefrac .AND.  &
          data_source( i ) % source == Set_to_MDI ) ) .OR.                 &
       data_source( i ) % source == Already_Processed ) Then

    If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
      Write (6,form) 'Already processed ', i,                       &
                      fields_out( i ) % stashmaster % section,      &
                      fields_out( i ) % stashmaster % item,         &
                      fields_out( i ) % stashmaster % name
    End If

    Cycle    ! just skip this iteration of the loop.

  End If

  ! Set the output data space - don't actually use this in the case
  ! of external files used for the data, but will do it anyway to
  ! simplify the algorithm

  output_field => fields_out( i )
  Call Rcf_Alloc_Field( output_field )

  Select Case( Data_Source( i ) % Source )

!---------------------------------------------------------------
! Data interpolated (or copied) from input dump
!---------------------------------------------------------------
    Case( Input_Dump )

      ! set up input fields
      Call Rcf_Locate( fields_out( i ) % stashmaster % section,    &
                       fields_out( i ) % stashmaster % item,       &
                       fields_in, field_count_in, pos )
      input_field => fields_in( pos )
      Call Rcf_Alloc_Field( input_field )

      Call Rcf_Read_Field( input_field, Hdr_In, decomp_rcf_input )

      ! Write out appropriate message
      Select Case( input_field % interp )

        Case( interp_h_only, interp_v_only, interp_all)
          If (mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
            Write (6,form) 'Interpolating Field', i,                &
                           input_field % stashmaster % section,     &
                           input_field % stashmaster % item,        &
                           input_field % stashmaster % name
          End If

        Case( interp_copy )
          If (mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
            Write (6,form) 'Copying Field', i,                      &
                           input_field % stashmaster % section,     &
                           input_field % stashmaster % item,        &
                           input_field % stashmaster % name
          End If

        Case( interp_no_op )
          If (mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
            Write (6,form) 'Skipping Field', i,                     &
                           input_field % stashmaster % section,     &
                           input_field % stashmaster % item,        &
                           input_field % stashmaster % name
          End If

          Call Rcf_DeAlloc_Field( input_field )
          Call Rcf_DeAlloc_Field( output_field )
          Cycle

       End Select

      ! convert fields approriately for interpolation
      Call Rcf_Pre_Interp_Transform( input_field, fields_in,     &
                                     field_count_in, hdr_in, orog_in )

      Call Rcf_Interpolate( input_field, output_field, Input_Grid, &
                              Output_Grid, interp_orog, orog_out )

      ! Convert fields back to original form (if possible) and
      ! perform any simple post-processing required
      Call Rcf_Post_Interp_Transform( output_field, fields_out,  &
                                      field_count_out )

      Call Rcf_DeAlloc_Field( input_field )

!------------------------------------------------------------------
! Data from Ancillary File
!------------------------------------------------------------------
    Case( Ancillary_File )

      If ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
        Write (6,form) 'Ancillary field', i,                       &
                    fields_out( i ) % stashmaster % section,       &
                    fields_out( i ) % stashmaster % item,          &
                    fields_out( i ) % stashmaster % name
      End If

!------------------------------------------------------------------
! Data to be set to zero
!------------------------------------------------------------------
    Case( Set_To_Zero )

      If ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
        Write (6,form) 'Set to Zero, field', i,                    &
                    output_field % stashmaster % section,          &
                    output_field % stashmaster % item,             &
                    output_field % stashmaster % name
      End If

      Select Case ( output_field % stashmaster % data_type )
        Case ( ppx_type_real )
          output_field % Data( :, : ) = 0.0

        Case ( ppx_type_int )
          output_field % Data_Int( : , : ) = 0

        Case ( ppx_type_log )
          output_field % Data_Log( : , : ) = .FALSE.
          If (mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
            Write (6,*) 'Assuming zero means FALSE for logicals'
          End If

        Case Default
          Write (6,*) 'Cannot set fields of this type to zero!'
          Write (6,*) 'ppx_type_? =',&
            output_field % stashmaster % data_type
      End Select

!------------------------------------------------------------------
! Data to be set as missing
!------------------------------------------------------------------
    Case( Set_To_MDI )

      If ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
        Write (6,form) 'Set to MDI ', i,                          &
                    output_field % stashmaster % section,         &
                    output_field % stashmaster % item,            &
                    output_field % stashmaster % name
      End If

      Select Case ( output_field % stashmaster % data_type )
        Case ( ppx_type_real )
          output_field % Data( :, : ) = RMDI

        Case ( ppx_type_int )
          output_field % Data_Int( : , : ) = IMDI

        Case Default
          Write (6,*) 'Cannot set fields of this type to MDI!'
      End Select

!----------------------------------------------------------------
! Tracer data
     !!!kdcorbin, 03/10 - Replaced using ATRACER file with a 
     !!!   user prognostic ancillary file
!----------------------------------------------------------------
    Case( Tracer_File )

      If ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
        Write (6,form) 'Tracer file data ', i,                   &
                    output_field % stashmaster % section,        &
                    output_field % stashmaster % item,           &
                    output_field % stashmaster % name
      End If

      ! Get Tracer File name from env var ATRACER
      !Call Fort_Get_Env( 'ATRACER', 7, TracerFileName, &
      !                   Max_Filename_Len, err )

      !If ( err /= 0 ) Then
      !  Write (6,*) 'Cannot get Tracer File name from Env Var ATRACER '
      !  cmessage = 'Cannot get Tracer File name from Env Var ATRACER'
      !  Call Ereport ( RoutineName, err, cmessage)
      !Endif

      TracerFileName=data_source(i)%Ancil_File

!     Check that Tracer File exists
      Inquire ( file=TracerFileName, exist=l_exist, iostat=IOStatus )

      If ( .not. l_exist ) then
        Write (6,*) 'Tracer File does not exist.'
        Write (6,*) 'File : ',TracerFileName
        ErrorStatus=10
        Cmessage = 'Tracer File does not exist.'
        Call Ereport ( RoutineName, ErrorStatus, Cmessage )
      End If

      Call Rcf_Get_Unit( Hdr_Aux % UnitNum )
! DEPENDS ON: file_open
      !Replaced ATRACER - kdcorbin, 03/10
      Call File_Open( Hdr_Aux % UnitNum, TracerFileName,Max_Filename_Len , &
             0 ,1 ,err )

      !Replaced IMDI with info in last two spots - kdcorbin, 08/10
      Call Rcf_Aux_File( Hdr_Aux, Hdr_Out, fields_out, field_count_out,&
                         tracers, output_field % stashmaster % section,&
                         output_field % stashmaster % item, &
                         output_field%stashmaster%section,  &
                         output_field%stashmaster%item)

! DEPENDS ON: file_close
      !Replace ATRACER with TracerFileName - kdcorbin, 03/10
      Call File_Close( Hdr_Aux % UnitNum,trim(TracerFileName), &
                   Max_Filename_Len ,1 ,0 ,err )
      Call Rcf_Free_Unit( Hdr_Aux % UnitNum)
      Call Rcf_FreeUMhdr( Hdr_Aux )

!----------------------------------------------------------------
! Data to be set to constant from the namelist
!----------------------------------------------------------------
    Case( Set_To_Const )
      If ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
        Write (6,form) 'Set to user const', i,                    &
                    output_field % stashmaster % section,         &
                    output_field % stashmaster % item,            &
                    output_field % stashmaster % name
      End If

      Select Case ( fields_out( i ) % stashmaster % data_type )
        Case ( ppx_type_real )
          output_field % Data( :, : ) = data_source( i ) % RConst

        Case ( ppx_type_int )
          output_field % Data_Int( : , : ) =                         &
                                     Nint( data_source( i ) % RConst )
        Case ( ppx_type_log )
          If ( data_source( i ) % Rconst > 0.5 ) Then
            output_field % Data_Log( : ,: ) = .TRUE.
          Else
            output_field % Data_Log( :, : ) = .FALSE.
          End If

      End Select

!----------------------------------------------------------------
! User prognostics from external dump
!----------------------------------------------------------------
    Case( External_Dump )

      If ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
        Write (6,form) 'User Prognostic', i,                     &
                    output_field % stashmaster % section,        &
                    output_field % stashmaster % item,           &
                    output_field % stashmaster % name
      End If

!     Check that ancillary file exists
      Inquire ( file=data_source( i ) % Ancil_File, &
               exist=l_exist, iostat=IOStatus )

      If ( .not. l_exist ) then
        Write (6,*) 'User Prognostic File does not exist.'
        Write (6,*) 'File : ',data_source( i ) % Ancil_File
        ErrorStatus=10
        Cmessage = 'User Prognostic File does not exist.'
        Call Ereport ( RoutineName, ErrorStatus, Cmessage )
      End If

      Call Rcf_Get_Unit( Hdr_Aux % UnitNum )
! DEPENDS ON: file_open
      Call File_Open( Hdr_Aux % UnitNum, data_source( i ) % Ancil_File,&
                        Max_Filename_Len ,0 ,1 ,err )

      Call Rcf_Aux_File( Hdr_Aux, Hdr_Out, fields_out, field_count_out,&
                     user_prog, data_source( i ) % Ancil_SctnC,        &
                     data_source( i ) % Ancil_ItemC,                   &
                     output_field % stashmaster % section,             &
                     output_field % stashmaster % item)

! DEPENDS ON: file_close
      Call File_Close(Hdr_Aux % UnitNum, data_source( i ) % Ancil_File,&
                        Max_Filename_Len ,1 ,0 ,err )
      Call Rcf_Free_Unit( Hdr_Aux % UnitNum)
      Call Rcf_FreeUMhdr( Hdr_Aux )

!-----------------------------------------------------------------
! Calculations for fields that are missing in the input dump.
! These are skipped now and done in a later loop
!-----------------------------------------------------------------
    Case( Field_Calcs)
      If ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
        Write (6,form) 'Calcs. needed, skip', i,                 &
                    output_field % stashmaster % section,        &
                    output_field % stashmaster % item,           &
                    output_field % stashmaster % name
      End If

!-----------------------------------------------------------------
! Unrecognised source for data - set it to missing
!-----------------------------------------------------------------
    Case Default

      ErrorStatus = -10
      Cmessage = 'Source code not recognised - will set field to MDI'
      Call Ereport( RoutineName, ErrorStatus, Cmessage )

      If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
        Write (6,form) 'Unknown source', i,                      &
                  & output_field % stashmaster % section,        &
                  & output_field % stashmaster % item, 'Set to MDI!'
      End If

      Select Case ( fields_out(i) % stashmaster % data_type )
        Case ( ppx_type_real )
          output_field % Data( :, : ) = RMDI

        Case ( ppx_type_int )
          output_field % Data_Int( : , : ) = IMDI

        Case ( ppx_type_log )
          output_field % Data_Log( :, : ) = .FALSE.
      End Select


  End Select

!-------------------------------------------------------------------
! Write out data if required
!-------------------------------------------------------------------
  If ( Data_Source( i ) % Source /= External_Dump  .AND.  &
       Data_Source( i ) % Source /= Tracer_File    .AND.  &
       Data_Source( i ) % Source /= Field_Calcs    .AND.  &
       Data_Source( i ) % Source /= Ancillary_File ) Then

    Call Rcf_Write_Field( output_field, Hdr_Out, decomp_rcf_output )

  End If

  Call Rcf_DeAlloc_Field( output_field )
End Do

!------------------------------------------------------------------
! Perform post-processing
!------------------------------------------------------------------
Call Rcf_Post_Process_Atmos( fields_in, field_count_in, orog_in,    &
                             fields_out, field_count_out, orog_out, &
                             hdr_in, hdr_out, data_source)

!------------------------------------------------------------------
! Tidy up orography fields
!------------------------------------------------------------------
Call Rcf_DeAlloc_Field( orog_in )
Call Rcf_DeAlloc_Field( orog_out )
Call Rcf_DeAlloc_Field( interp_orog )

!------------------------------------------------------------------
! Need to revisit the source = field_calcs (8) fields
!------------------------------------------------------------------

Call Rcf_Field_Calcs( fields_in, fields_out, field_count_in, &
                      field_count_out, data_source, hdr_in, hdr_out )

!-------------------------------------------------------------------
! Rotate the winds if so required
!-------------------------------------------------------------------
If (h_int_active .AND. Output_Grid % Rotated) Then
  If (mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
    Write (6,*) 'Rotating Output Winds'
  End If

  Call Rcf_Rotate( fields_out, field_count_out, Output_Grid, Hdr_Out, &
                   decomp_rcf_output, FromStandard)
End If

!-------------------------------------------------------------------
! Transplant Data
!-------------------------------------------------------------------
If ( TRANS ) Then

  Call Rcf_Get_Unit( Hdr_Aux % UnitNum )
! DEPENDS ON: file_open
  Call File_Open( Hdr_Aux % UnitNum, 'TRANSP',6 ,0 ,0 ,err )

  Call Rcf_Aux_File( Hdr_Aux, Hdr_Out, fields_out, field_count_out,&
                     transplant, IMDI, IMDI, IMDI, IMDI)

! DEPENDS ON: file_close
  Call File_Close( Hdr_Aux % UnitNum, 'TRANSP',6 ,0 ,0 ,err )
  Call Rcf_Free_Unit( Hdr_Aux % UnitNum)
  Call Rcf_FreeUMhdr( Hdr_Aux )

End If

!-------------------------------------------------------------------
! UARS Data
!-------------------------------------------------------------------
If ( UARS ) Then

  Call Rcf_Get_Unit( Hdr_Aux % UnitNum )
! DEPENDS ON: file_open
  Call File_Open( Hdr_Aux % UnitNum, 'SSU',3 ,0 ,0 ,err )

  Call Rcf_Aux_File( Hdr_Aux, Hdr_Out, fields_out, field_count_out,&
                     uars_data, IMDI, IMDI, IMDI, IMDI)

! DEPENDS ON: file_close
  Call File_Close( Hdr_Aux % UnitNum, 'SSU',3 ,0 ,0 ,err )
  Call Rcf_Free_Unit( Hdr_Aux % UnitNum)
  Call Rcf_FreeUMhdr( Hdr_Aux )

End If

End Subroutine Rcf_Create_Dump

End Module Rcf_Create_Dump_Mod
#endif
