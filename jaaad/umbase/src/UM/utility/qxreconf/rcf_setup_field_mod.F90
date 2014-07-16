#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initialises the field data type for a given dump and grid

Module Rcf_Setup_Field_mod

!  Subroutine Rcf_Setup_Field - sets up a field data type array
!
! Description:
!   Sets up the field array for a given grid and dump and "fills in"
!   all the relevant data parts.
!
! Method:
!   Sizes are calculated based on grid code.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Setup_Field( field, hdr, grid, field_count, title, &
                        local_lsm_size )

Use Rcf_Address_Length_Mod, Only : &
    Rcf_Address_Length

Use Rcf_set_interp_flags_mod, Only : &
    interp_no_op

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Parvars_Mod, Only  :  &
    mype,                 &
    nproc

Use Rcf_UMhead_Mod, Only :    &
    UM_header_type

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Use Ereport_mod, Only :   &
    Ereport

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Oper

Use Rcf_Exppx_Mod, Only :     &
    Rcf_Exppx

Use Rcf_Recon_Mod, Only : &
    Rimwidtha

Use Rcf_Level_Code_Mod, Only : &
    Rcf_Level_Code

Use Rcf_Recon_Mod, Only : &
    Var_Recon


Implicit None

! Arguments
Type (um_header_type), Intent (In)  :: hdr      ! Dump header
Type (grid_type), Intent (In)       :: grid     ! grid for fields
Type (field_type), Pointer          :: field( : )
Integer, Intent (Out)               :: field_count ! no. of fields
Character (Len=*), Intent(In)       :: title
Integer, Intent(In), Optional       :: local_lsm_size
                                       ! size of local land point
                                       ! fields if known

! Comdecks
#include "clookadd.h"
#include "cppxref.h"
#include "c_mdi.h"

! Local data
Integer            :: i, k        ! loopers
Integer            :: ErrorStatus
Integer            :: model
Integer            :: sec_item
Integer            :: item
Integer            :: section
Integer            :: start_level
Integer            :: end_level
Integer            :: levels
Integer            :: land_sea
Integer            :: size
Integer            :: lsm_size
Logical            :: new_field
Character (Len=80) :: Cmessage
Character (Len=80) :: Phrase
Character (Len=*), Parameter :: RoutineName = 'setupfield'


! Note that the grid (will be either input or output) should
! correspond to the one referred to in the header lookups etc.

!-----------------------------------------------------------------
! Initialisation of lsm_size from optional value input
!-----------------------------------------------------------------
If (Present( local_lsm_size ) ) Then
  lsm_size = local_lsm_size
Else
  lsm_size = imdi
EndIf

!--------------------------------------------------------------
! Tests for allowed behaviour/values
!-------------------------------------------------------------
If (Associated(field)) Then
  Cmessage = 'Field sizes already set'
  ErrorStatus = -20
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
  Goto 9999
End If

If (grid % model_levels /= hdr % IntC(8)) Then
  Cmessage = 'Grid and headers disagree on number of levels'
  ErrorStatus = 10
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

!--------------------------------------------------------------
! Need to count the total number of fields for allocation
! of memory.
!---------------------------------------------------------------
field_count = 1
Do i = 1, hdr % Len2Lookup
  If ( i /= 1 ) Then
    ! If not the same field as before
    If ( .NOT. (hdr % Lookup(item_code, i) ==          &
                hdr % Lookup(item_code, i - 1) ) ) Then
      field_count = field_count + 1
    End if
  End If
End Do

! Allocate the space
Allocate( field( field_count ) )

!---------------------------------------------------------------
! Now initialise the fields
!---------------------------------------------------------------

field_count     = 1
new_field = .TRUE.

field( field_count ) % levels   = 1
field( field_count ) % dump_pos = 1

Do i = 1, hdr % Len2Lookup
  If ( i /= 1 ) Then
    If (hdr % Lookup(item_code, i) ==          &
        hdr % Lookup(item_code, i - 1) )  Then

      field( field_count ) % levels = field( field_count ) % levels + 1
      new_field = .FALSE.
    Else
      field_count = field_count + 1
      field( field_count ) % levels = 1
      field( field_count ) % dump_pos =  &
                                 field( field_count - 1) % dump_pos +  &
                                 field(field_count - 1) % levels
      new_field = .TRUE.
    End If
  End If

  If ( new_field ) Then
    Nullify( field( field_count ) % Data )
    Nullify( field( field_count ) % Data_Int )
    Nullify( field( field_count ) % Data_Log )
  End If

  ! Default - nothing doing....
  field( field_count ) % interp =  interp_no_op

! Need STASHmaster record for some information.
  sec_item = hdr % Lookup(item_code, i)
  item     = Mod( hdr % Lookup(item_code, i), 1000)
  section  = (hdr % Lookup(item_code, i) - item) / 1000
  model    = hdr % Lookup(model_code , i)

  field( field_count) % stashmaster => Rcf_Exppx( model, section, item )

! Need number of levels for lbc sizing
  If ( field( field_count ) % stashmaster % grid_type ==           &
                                            ppx_atm_lbc_theta .OR. &
       field( field_count ) % stashmaster % grid_type ==           &
                                            ppx_atm_lbc_u     .OR. &
       field( field_count ) % stashmaster % grid_type ==           &
                                            ppx_atm_lbc_v ) Then

    Call Rcf_Level_Code( field( field_count ) % stashmaster % lb_code, &
                         start_level, grid )
    Call Rcf_Level_Code( field( field_count ) % stashmaster % lt_code, &
                         end_level,   grid )

    levels = (end_level - start_level) + 1
  End If

!------------------------------------------------------------------
! Sizes depend on the grid type - cases will need to be added as
! they are needed. Assume we don't mix our grids in a single
! call to this function (ie some fields B, some C grid etc)
!------------------------------------------------------------------
  Select Case ( field( field_count ) % stashmaster % grid_type )
    Case (ppx_atm_tall, ppx_atm_tsea, ppx_ocn_tall)      ! Theta points
      field( field_count ) % rows            = grid % loc_p_rows
      field( field_count ) % row_len         = grid % loc_p_row_length
      field( field_count ) % level_size      = grid % loc_p_rows *     &
                                               grid % loc_p_row_length

      field( field_count ) % glob_rows       = grid % glob_p_rows
      field( field_count ) % glob_row_len    = grid % glob_p_row_length
      field( field_count ) % glob_level_size = grid % glob_p_rows *    &
                                               grid % glob_p_row_length

    Case (ppx_atm_uall, ppx_atm_usea, ppx_ocn_uall)    ! B grid u points
      field( field_count ) % rows            = grid % loc_u_rows
      field( field_count ) % row_len         = grid % loc_u_row_length
      field( field_count ) % level_size      = grid % loc_u_rows *     &
                                               grid % loc_u_row_length

      field( field_count ) % glob_rows       = grid % glob_u_rows
      field( field_count ) % glob_row_len    = grid % glob_u_row_length
      field( field_count ) % glob_level_size = grid % glob_u_rows *    &
                                               grid % glob_u_row_length

    Case (ppx_atm_cuall)      ! C grid u points
      field( field_count ) % rows            = grid % loc_u_rows
      field( field_count ) % row_len         = grid % loc_u_row_length
      field( field_count ) % level_size      = grid % loc_u_rows *     &
                                               grid % loc_u_row_length

      field( field_count ) % glob_rows       = grid % glob_u_rows
      field( field_count ) % glob_row_len    = grid % glob_u_row_length
      field( field_count ) % glob_level_size = grid % glob_u_rows *    &
                                               grid % glob_u_row_length

    Case (ppx_atm_cvall)      ! C grid v points
      field( field_count ) % rows            = grid % loc_v_rows
      field( field_count ) % row_len         = grid % loc_v_row_length
      field( field_count ) % level_size      = grid % loc_v_rows *     &
                                               grid % loc_v_row_length

      field( field_count ) % glob_rows       = grid % glob_v_rows
      field( field_count ) % glob_row_len    = grid % glob_v_row_length
      field( field_count ) % glob_level_size = grid % glob_v_rows *    &
                                               grid % glob_v_row_length

    Case (ppx_atm_ozone)      ! Ozone grid
      field( field_count ) % rows            = grid % loc_p_rows
      field( field_count ) % glob_rows       = grid % glob_p_rows

      If ( hdr % Lookup(lbnpt,i) == 1) Then
        field( field_count ) % row_len         = 1
        field( field_count ) % glob_row_len    = 1
      Else
        field( field_count ) % row_len         = grid % loc_p_row_length
        field( field_count ) % glob_row_len    = &
                                            grid % glob_p_row_length
      End If

      field( field_count ) % level_size      =  &
              field(field_count) % row_len * field(field_count) % rows
      field( field_count ) % glob_level_size =  &
      field(field_count) % glob_row_len * field(field_count) % glob_rows

    Case (ppx_atm_compressed)     ! Land compressed points
      ! Can only set global size here. Local level_size will need to
      ! be set when sizes are available.
      field( field_count ) % rows            = imdi
      field( field_count ) % row_len         = imdi
      field( field_count ) % level_size      = lsm_size


      field( field_count ) % glob_rows       = imdi
      field( field_count ) % glob_row_len    = imdi
      field( field_count ) % glob_level_size = hdr % Lookup(lblrec,i)

    Case (ppx_atm_river)    ! River routing points
      field( field_count ) % rows            = grid % loc_r_rows
      field( field_count ) % row_len         = grid % loc_r_row_length
      field( field_count ) % level_size      = grid % loc_r_rows *     &
                                               grid % loc_r_row_length

      field( field_count ) % glob_rows       = grid % glob_r_rows
      field( field_count ) % glob_row_len    = grid % glob_r_row_length
      field( field_count ) % glob_level_size = grid % glob_r_rows *    &
                                               grid % glob_r_row_length

    Case (ppx_atm_lbc_theta, ppx_atm_lbc_u, ppx_atm_lbc_v ) ! Atmos LBC

      ! Uses Address_Length to calculate the sizes - this assumes that
      ! both the input and output grid LBCs are the same sizes...

      Call Rcf_Address_Length(                                        &
           field( field_count) % stashmaster % grid_type,             &
           field( field_count ) % stashmaster % halo_type, size )

      field( field_count ) % glob_rows       = imdi
      field( field_count ) % glob_row_len    = imdi
      field( field_count ) % glob_level_size = size * levels

      field( field_count ) % rows            = imdi
      field( field_count ) % row_len         = imdi
      size = field( field_count ) % glob_level_size / nproc
      If (mype == nproc - 1) Then
        size = size + Mod( field( field_count ) % glob_level_size,nproc)
      End If
      field( field_count ) % level_size      = size

    Case Default

      If (section /= 0 .AND. section /= 33 .AND. section /= 34 &
          .AND. .NOT. var_recon) Then
        ! Diagnostic in input dump.
        ! Rcf does not process diagnostics in input dump, so set
        ! relevant field dimensions to imdi in field.  Only 
        ! reconfigure diagnostics for VAR.

        field( field_count ) % glob_rows       = imdi
        field( field_count ) % glob_row_len    = imdi
        field( field_count ) % rows            = imdi
        field( field_count ) % row_len         = imdi
        field( field_count ) % glob_level_size = imdi
        field( field_count ) % level_size      = imdi

      Else

        Write (6,*) 'Unsupported Grid-type ', &
                     field( field_count ) % stashmaster % grid_type
        Cmessage = 'Grid type not yet catered for - size will be &
                   &unset in field data structure'
        ErrorStatus = -10
        Call Ereport( RoutineName, ErrorStatus, Cmessage )

      End If
  End Select
End Do

!------------------------------------------------------------------
! This should have all values calculated correctly
! Only remains to print them out if required.
!------------------------------------------------------------------
If (mype == 0 .AND. PrintStatus >= PrStatus_Oper ) Then
  WRITE(6,'(''  '',/,'' '',A/)')TITLE

  i = 1
  Do k = 1, field_count

    Phrase = field(k) % stashmaster % name
    If ( field(k) % stashmaster % grid_type == ppx_atm_compressed ) Then
      land_sea = 1
    Else
      land_sea = 0
    End If

    i = i + field(k) % levels
    WRITE(6,'('' '',I4,I5,I8,I4,3I6,2x,A36)')                 &
         land_sea, field(k) % levels,                         &
         field(k) % glob_level_size,                          &
         field(k) % stashmaster % data_type,                  &
         field(k) % stashmaster % section,                    &
         field(k) % stashmaster % item,                       &
         field(k) % dump_pos, phrase
  End Do
End if

9999 Continue
End Subroutine Rcf_Setup_Field
End Module Rcf_Setup_Field_Mod
#endif
