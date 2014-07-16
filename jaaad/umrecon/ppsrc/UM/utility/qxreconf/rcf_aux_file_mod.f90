
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reads in auxillary data into dump

Module Rcf_Aux_File_Mod

!  Subroutine Rcf_Aux_File  - reads auxillary data
!
! Description:
!   Reads data from external dumps and incorporates it into the
!   output dump. Four modes are supported
!                Tracer Data
!                User Prognostics
!                UARS Data
!                Area Tranplants
!
! Method:
!   A fields array is set up for the auxillary file and the relevant
!   fields located therein. Data is copied as appropriate for the
!   Mode (ie transplants within an area, all of a user prognostic
!   copied, upper levels only of UARS or Tracers copied).
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


! Parameters describing actions possible
Integer, Parameter       :: tracers   = 1
Integer, Parameter       :: user_prog = 2
Integer, Parameter       :: uars_data = 3
Integer, Parameter       :: transplant= 4

Contains

Subroutine Rcf_Aux_File( Hdr_Aux, Hdr_Out, Fields_Out, Field_Count_Out,&
                         Action, Sctn_code_aux, Item_Code_aux,         & 
                         Sctn_code_out, Item_code_out )

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

Use Rcf_UMhead_Mod, Only : &
    Um_Header_Type

Use Rcf_Field_Type_Mod, Only : &
    Field_Type

Use Rcf_Parvars_mod, Only : &
    mype

Use rcf_trans_mod              ! all of it

Use Rcf_Write_Field_Mod, Only : &
    Rcf_Write_Field

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Setup_Field_mod, Only : &
    Rcf_Setup_Field

Use Rcf_Alloc_Field_mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

use Rcf_Locate_mod, Only : &
    Rcf_Locate

Use Rcf_ReadUMhdr_Mod, Only : &
    Rcf_ReadUMhdr

Use Rcf_HeadAddress_Mod, Only :  &
IC_XLen,                         IC_YLen,    &
FH_DTYear,                       FH_VTYear,  &
FH_DTDayNo,                      FH_VTDayNo,       &
FH_Dataset,                      FH_Dataset_Ancil

Use Rcf_Grid_Type_Mod, Only :&
    Output_Grid

Use Rcf_DecompTP_Mod, Only : &
    decomp_rcf_output

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Lsm_Mod, Only : &
    local_land_out,     &
    local_lsm_out

Use Rcf_Global_To_Local_Mod, Only : &
    Rcf_Global_To_Local_Subdomain

Implicit None

! Arguments
Type( Um_Header_type ), Intent(InOut)  :: Hdr_Aux ! only unit no.,
                                                  ! file already open
Type( Um_Header_type ), Intent(In)     :: Hdr_Out
Type( Field_Type ), Pointer            :: Fields_Out(:)

Integer, Intent(In)                    :: Field_Count_Out
Integer, Intent(In)                    :: sctn_code_aux
Integer, Intent(In)                    :: item_code_Aux
Integer, Intent(In)                    :: sctn_code_out
Integer, Intent(In)                    :: item_code_out ! only for UP
Integer, Intent(In)                    :: Action

! Comdecks
!LL  Comdeck: CPPXREF --------------------------------------------------
!LL
!LL  Purpose: Holds PARAMETER definitions to describe the structure of
!LL           each STASHmaster file record plus some valid entries.
!LL
!LL  Author    Dr T Johns
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
!LL                  1.Removes the limit on primary STASH item numbers.
!LL                  2.Removes the assumption that (section,item)
!LL                    defines the sub-model.
!LL                  3.Thus allows for user-prognostics.
!LL                  Add a PPXREF record for model number.
!LL  4.0   26/07/95  T.Johns.  Add codes for real/int/log data types.
!LL  3.5   10/3/94   Sub-Models project:
!LL                 List of PPXREF addressing codes augmented, in order
!LL                 to include all of the pre_STASH master information
!LL                 in the new PPXREF file.
!LL                 PPXREF_CODELEN increased to 38.
!LL                 PPXREF_IDLEN deleted - no longer relevant.
!LL                   S.J.Swarbrick
!LL  4.1   June 96  Wave model parameters included.
!LL                 ppx_ address parameters adjusted to allow for
!LL                  reading option code as 4x5 digit groups.
!LL                   S.J.Swarbrick
!LL  5.0   29/06/99  Add halo type parameter for new dynamics.
!LL                  New grid codes for LAM boundary conditions
!LL                  D.M. Goddard
!LL  5.1   07/03/00  Fixed/Free format conversion
!LL  5.2   19/09/00  Added ppx_atm_lbc_orog descriptor   P.Burton
!LL  5.3   21/08/01  Added ocean lbc descriptors.   M. J. Bell
!LL  5.3   23/07/01  Add valid pp_lbvc codes referenced in UM. R Rawlins
!LL  5.5   30/01/03  Option code increase from 20 to 30 digits thus
!LL                  requiring option code address range increase by
!LL                  2 so all subsequent addressing codes need to be
!LL                  increased by 2 to make a gap.
!LL                  W Roseblade
!LL
!LL  Logical components covered: C40
!LL
!-----------------------------------------------------------------------
! Primary file record definition
      ! length of ID in a record
      Integer, Parameter :: PPXREF_IDLEN      = 2

      ! total length of characters *WARNING* must be multiple of 4
      ! to avoid overwriting
      Integer, Parameter :: PPXREF_CHARLEN    = 36

      ! number of packing profiles
      Integer, Parameter :: PPXREF_PACK_PROFS = 10

      ! total length of codes = no. of codes (excluding profs)
      ! + pack_profs
      Integer, Parameter :: PPXREF_CODELEN    = 33 + PPXREF_PACK_PROFS

! Derived file record sizes
      ! Assume that an integer is at least 4 bytes long. Wastes some
      ! space on an 8 byte machine.
      ! ppx_charword = 9.
      Integer, Parameter :: PPX_CHARWORD      = ((PPXREF_CHARLEN+3)/4)

      ! read buffer record length
      Integer, Parameter :: PPX_RECORDLEN = PPX_CHARWORD+PPXREF_CODELEN
!
!-----------------------------------------------------------------------
! Addressing codes within PPXREF
      Integer, Parameter ::  ppx_model_number   = 1  ! Model number
                                                     ! address
      Integer, Parameter ::  ppx_section_number = 2  ! Section number
                                                     ! address
      Integer, Parameter ::  ppx_item_number    = 3  ! Item number
                                                     ! address
      Integer, Parameter ::  ppx_version_mask   = 4  ! Version mask
                                                     ! address
      Integer, Parameter ::  ppx_space_code     = 5  ! Space code
                                                     ! address
      Integer, Parameter ::  ppx_timavail_code  = 6  ! Time availability
                                                     !  code  address
      Integer, Parameter ::  ppx_grid_type      = 7  ! Grid type code
                                                     ! address
      Integer, Parameter ::  ppx_lv_code        = 8  ! Level type code
                                                     ! address
      Integer, Parameter ::  ppx_lb_code        = 9  ! First level code
                                                     !  address
      Integer, Parameter ::  ppx_lt_code        =10  ! Last level code
                                                     ! address
      Integer, Parameter ::  ppx_lev_flag       =11  ! Level compression
                                                     !  flag  address
      Integer, Parameter ::  ppx_opt_code       =12  ! Sectional option
                                                     ! code  address
      Integer, Parameter ::  ppx_pt_code        =18  ! Pseudo dimension
                                                     ! type  address
      Integer, Parameter ::  ppx_pf_code        =19  ! First pseudo dim
                                                     ! code  address
      Integer, Parameter ::  ppx_pl_code        =20  ! Last pseudo dim
                                                     ! code  address
      Integer, Parameter ::  ppx_ptr_code       =21  ! Section 0 point-
                                                     ! back code address
      Integer, Parameter ::  ppx_dump_packing   =22  ! Dump packing code
                                                     ! address
      Integer, Parameter ::  ppx_lbvc_code      =23  ! PP LBVC code
                                                     ! address
      Integer, Parameter ::  ppx_rotate_code    =24  ! Rotation code
                                                     ! address
      Integer, Parameter ::  ppx_field_code     =25  ! PP field code
                                                     ! address
      Integer, Parameter ::  ppx_user_code      =26  ! User code address
      Integer, Parameter ::  ppx_meto8_levelcode=27  ! CF level code
                                                     ! address
      Integer, Parameter ::  ppx_meto8_fieldcode=28  ! CF field code
                                                     ! address
      Integer, Parameter ::  ppx_cf_levelcode   =27
      Integer, Parameter ::  ppx_cf_fieldcode   =28
      Integer, Parameter ::  ppx_base_level     =29  ! Base level code
                                                     ! address
      Integer, Parameter ::  ppx_top_level      =30  ! Top level code
                                                     ! address
      Integer, Parameter ::  ppx_ref_lbvc_code  =31  ! Ref level LBVC
                                                     ! code address
      Integer, Parameter ::  ppx_data_type      =32  ! Data type code
                                                     ! address
      Integer, Parameter ::  ppx_halo_type      =33
      Integer, Parameter ::  ppx_packing_acc    =34  ! Packing accuracy
                                                     ! code  address
      Integer, Parameter ::  ppx_pack_acc       =34  ! Must be last:


                                                 ! multiple pack_acc to
                                                 ! fill up remaining
                                                 ! array elements


!-------------------------------------------------------------------
! Valid grid type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_atm_nonstd=0      ! Non-standard atmos
                                                  ! grid
      Integer, Parameter :: ppx_atm_tall=1        ! All T points (atmos)
      Integer, Parameter :: ppx_atm_tland=2       ! Land-only T points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_tsea=3        ! Sea-only T points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_tzonal=4      ! Zonal field at T
                                                  ! points  (atmos)
      Integer, Parameter :: ppx_atm_tmerid=5      ! Merid field at T
                                                  ! points  (atmos)
      Integer, Parameter :: ppx_atm_uall=11       ! All u points (atmos)
      Integer, Parameter :: ppx_atm_uland=12      ! Land-only u points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_usea=13       ! Sea-only u points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_uzonal=14     ! Zonal field at u
                                                  ! points  (atmos)
      Integer, Parameter :: ppx_atm_umerid=15     ! Merid field at u
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_scalar=17     ! Scalar (atmos)
      Integer, Parameter :: ppx_atm_cuall=18      ! All C-grid (u)
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_cvall=19      ! All C-grid (v)
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_compressed=21 ! Compressed land
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_ozone=22      ! Field on ozone
                                                  ! grid (atmos)
      Integer, Parameter :: ppx_atm_river=23      ! River routing
                                                  ! grid (atmos)
      Integer, Parameter :: ppx_atm_rim=25        ! Rim type field
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_theta=26  ! All T points
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_u=27      ! All u points
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_v=28      ! All v points
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_orog=29   ! Orography field
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_ocn_nonstd=30     ! Non-standard ocean
                                                  ! grid
      Integer, Parameter :: ppx_ocn_tcomp=31      ! Compressed T points
                                                  !  (ocean)
      Integer, Parameter :: ppx_ocn_ucomp=32      ! Compressed u points
                                                  !  (ocean)
      Integer, Parameter :: ppx_ocn_tall=36       ! All T points incl.
                                                  ! cyclic  (ocean)
      Integer, Parameter :: ppx_ocn_uall=37       ! All u points incl.
                                                  ! cyclic  (ocean)
      Integer, Parameter :: ppx_ocn_cuall=38      ! All C-grid (u)
                                                  ! points (ocean)
      Integer, Parameter :: ppx_ocn_cvall=39      ! All C-grid (v)
                                                  ! points (ocean)
      Integer, Parameter :: ppx_ocn_tfield=41     ! All non-cyclic T
                                                  ! points  (ocean)
      Integer, Parameter :: ppx_ocn_ufield=42     ! All non-cyclic u
                                                  ! points  (ocean)
      Integer, Parameter :: ppx_ocn_tzonal=43     ! Zonal n-c field at
                                                  ! T points  (ocean)
      Integer, Parameter :: ppx_ocn_uzonal=44     ! Zonal n-c field at
                                                  ! u points (ocean)
      Integer, Parameter :: ppx_ocn_tmerid=45     ! Merid n-c field at
                                                  ! T points  (ocean)
      Integer, Parameter :: ppx_ocn_umerid=46     ! Merid n-c field at
                                                  ! u points  (ocean)
      Integer, Parameter :: ppx_ocn_scalar=47     ! Scalar (ocean)
      Integer, Parameter :: ppx_ocn_rim=51        ! Rim type field
                                                  ! (LAM BCs ocean)
      Integer, Parameter :: ppx_ocn_lbc_theta=52  ! Ocean rim fields
      Integer, Parameter :: ppx_ocn_lbc_u=53      ! on T & U grids
      Integer, Parameter :: ppx_wam_all=60        ! All points (wave
                                                  ! model)
      Integer, Parameter :: ppx_wam_sea=62        ! Sea points only
                                                  ! (wave model)
      Integer, Parameter :: ppx_wam_rim=65        ! Rim type field
                                                  ! (LAM BCs wave)

!--------------------------------------------------------------------
! Valid rotation type codes
!--------------------------------------------------------------------
      Integer, Parameter :: ppx_unrotated=0       ! Unrotated output
                                                  ! field
      Integer, Parameter :: ppx_elf_rotated=1     ! Rotated ELF field

!-------------------------------------------------------------------
! Valid level type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_full_level=1      ! Model full level
      Integer, Parameter :: ppx_half_level=2      ! Model half level
      Integer, Parameter :: ppx_rho_level=1       ! Model rho level
      Integer, Parameter :: ppx_theta_level=2     ! Model theta level
      Integer, Parameter :: ppx_single_level=5    ! Model single level
      Integer, Parameter :: ppx_soil_level=6      ! Deep Soil level

!-------------------------------------------------------------------
! Valid data type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_type_real=1       ! Real data type
      Integer, Parameter :: ppx_type_int=2        ! Integer data type
      Integer, Parameter :: ppx_type_log=3        ! Logical data type

!-------------------------------------------------------------------
! Valid meto8 level type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_meto8_surf=9999   ! MetO8 surface type
                                                  ! code

!-------------------------------------------------------------------
! Valid dump packing codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_pack_off=0        ! Field not packed
                                                  ! (ie. 64 bit)
      Integer, Parameter :: ppx_pack_32=-1        ! Field packed to
                                                  ! 32 bit in  dump
      Integer, Parameter :: ppx_pack_wgdos=1      ! Field packed by
                                                  ! WGDOS method
      Integer, Parameter :: ppx_pack_cfi1=11      ! Field packed using
                                                  ! CFI1  (ocean)

!-------------------------------------------------------------------
! Add valid lbvc codes referenced in model (pp header output labels)
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_lbvc_height  =  1 ! height
      Integer, Parameter :: ppx_lbvc_depth   =  2 ! depth (ocean)
      Integer, Parameter :: ppx_lbvc_pressure=  8 ! pressure
      Integer, Parameter :: ppx_lbvc_theta   = 19 ! potential T
      Integer, Parameter :: ppx_lbvc_hybrid  = 65 ! hybrid height(atmos)
      Integer, Parameter :: ppx_lbvc_PV      = 82 ! potential vorticity
      Integer, Parameter :: ppx_lbvc_surface =129 ! surface

! Local variables
Integer                     :: i, j    ! loopers
Integer                     :: k,l     ! loopers
Integer                     :: pos_out ! field position
Integer                     :: pos_aux ! field position
Integer                     :: lrow1   ! Local positions of
Integer                     :: lrow2   ! trans data validity
Integer                     :: lcol1
Integer                     :: lcol2
Integer                     :: level_base
Integer                     :: level_top
Integer                     :: size
Integer                     :: field_count_aux
Integer                     :: ErrorStatus
Integer                     :: Copy_Count ! counter for no. of x a field
                                          ! is copied

Integer, Parameter          :: st_no_data = -3

Character (Len=20)          :: title
Character (Len=*), Parameter:: RoutineName='Rcf_Aux_File'
Character (Len=80)          :: Cmessage

Type( field_type ), Pointer :: fields_aux(:)

Nullify( fields_aux )
!------------------------------------------------------------------
! Read in the Auxillary file header
!------------------------------------------------------------------
Call Rcf_ReadUMhdr( Hdr_Aux )

!------------------------------------------------------------------
! If uars or tracers check that data time of the aux file matches
! verification time of the output file
!------------------------------------------------------------------
If ( action == tracers .OR. action == uars_data ) Then
  Do i = 0, 6
    If ( Hdr_Aux % FixHd( FH_DTYear + i ) /= &
         Hdr_Out % FixHd( FH_VTYear + i ) ) Then
      Write (6,*) 'Mismatch in date info for auxillary file'
      Write (6,*) 'Aux file Data times = ', (Hdr_Aux % Fixhd( j ), &
                                             j = FH_DTYear, FH_DTDayNo )
      Write (6,*) 'Out file Ver. times = ', (Hdr_Out % Fixhd( j ), &
                                             j = FH_VTYear, FH_VTDayNo )

      Cmessage = 'Date information mismatch between aux and dump files'
      ErrorStatus = 10
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If
  End Do
End If

!-------------------------------------------------------------------
! Check that resolutions of two files match
!-------------------------------------------------------------------
If ( Hdr_Aux % IntC( IC_XLen ) /= Hdr_Out % IntC( IC_XLen ) .OR. &
     Hdr_Aux % IntC( IC_YLen ) /= Hdr_Out % IntC( IC_YLen ) ) Then

  Cmessage = 'Dimensions of AUX file and dump file do not match'
  ErrorStatus = 20
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

!-------------------------------------------------------------------
! Set up field data-types for aux file
! (Grid resolutions should be as for output grid)
!-------------------------------------------------------------------
title = 'Auxillary File'
If ( action == user_prog .AND.                                 &
     Hdr_Aux % FixHd(FH_Dataset)  == FH_Dataset_Ancil ) Then
  ! User prognostic ancillary files have full fields for land-only
  ! fields
  Call Rcf_Setup_Field( fields_aux, Hdr_Aux, Output_Grid,      &
                        field_count_aux, title,                &
                        Output_Grid % loc_p_rows *             &
                        Output_Grid % loc_p_row_length )
Else
  Call Rcf_Setup_Field( fields_aux, Hdr_Aux, Output_Grid,      &
                        field_count_aux, title, local_land_out )
End If

!-------------------------------------------------------------------
! Main data handling/replacement - start with TRANS data
!-------------------------------------------------------------------
If ( action == transplant ) Then
  Do i = 1, num_trans

    If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
      Write (6,*) 'Transplanting data for stashcode ', itemc_array( i )
    End If

    ! Read aux field
    Call Rcf_Locate( sctnc_array( i ), itemc_array( i ),             &
                     fields_aux, field_count_aux, pos_aux)
    Call Rcf_Alloc_Field( fields_aux( pos_aux ) )
    Call Rcf_Read_Field( fields_aux( pos_aux ), Hdr_Aux,             &
                         decomp_rcf_output )

    ! Read dump field
    Call Rcf_Locate( sctnc_array( i ), itemc_array( i ),             & 
                     fields_out, field_count_out, pos_out)

    Call Rcf_Alloc_Field( fields_out( pos_out ) )
    Call Rcf_Read_Field( fields_out( pos_out ), Hdr_Out,             &
                         decomp_rcf_output )

    ! Cannot (yet) do a transplant for a land compressed field
    If (fields_out(pos_out) % stashmaster % grid_type ==             &
                                            ppx_atm_compressed) Then
      Cmessage = 'Cannot transplant a land compressed field'
      ErrorStatus = 30
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

    ! Convert the co-ordinates for replacement with local ones
    Call Rcf_Global_to_Local_Subdomain( .FALSE., .FALSE.,             &
                    fields_out( pos_out ) % stashmaster % grid_type , &
                    mype,                                             &
                    Row2_array( i ), Col2_array( i ),                 &
                    Row1_array( i ), Col1_array( i ),                 &
                    lrow2, lcol2, lrow1, lcol1 )

    ! If the local area is on my pe, do the transplant
    If ( lrow1 /= st_no_data .AND. lrow2 /= st_no_data .AND. &
         lcol1 /= st_no_data .AND. lcol2 /= st_no_data ) Then

      Select Case ( fields_out( pos_out ) % stashmaster % data_type )
        Case (ppx_type_real)
          Do j = lrow1, lrow2
            Do k = lcol1, lcol2
              l = (j-1) * fields_out(pos_out) % row_len + k
              fields_out(pos_out) %                      &
                    data(l,lev1_array(i):lev2_array(i))= &
              fields_aux(pos_aux) % data(l,lev1_array(i):lev2_array(i))
            End Do
          End Do

        Case (ppx_type_int)
          Do j = lrow1, lrow2
            Do k = lcol1, lcol2
              l = (j-1) * fields_out(pos_out) % row_len + k
              fields_out(pos_out) %                          &
                    data_int(l,lev1_array(i):lev2_array(i))= &
              fields_aux(pos_aux) %                          &
                    data_int(l,lev1_array(i):lev2_array(i))
            End Do
          End Do

        Case (ppx_type_log)
          Do j = lrow1, lrow2
            Do k = lcol1, lcol2
              l = (j-1) * fields_out(pos_out) % row_len + k
              fields_out(pos_out) %                          &
                    data_log(l,lev1_array(i):lev2_array(i))= &
              fields_aux(pos_aux) %                          &
                    data_log(l,lev1_array(i):lev2_array(i))
            End Do
          End Do

      End Select
    End If

    ! Write out the field
    Call Rcf_Write_Field( fields_out( pos_out ), Hdr_Out, &
                         decomp_rcf_output)

    Call Rcf_Dealloc_Field( fields_out( pos_out ) )
    Call Rcf_Dealloc_Field( fields_aux( pos_aux ) )

  End Do

Else       ! not transplant....
!----------------------------------------------------------------
! Loop through the auxillary fields for consideration for
! output dump inclusion.
! After problems with an ocean dump containing a time series diagnostic
! feild with identical section and item no. but different size the
! Copy_Count now forces only the first instance in the dump to be used
!----------------------------------------------------------------
  Copy_Count=0
  Do i = 1, field_count_aux
    If ( (action == tracers     .OR.   &
          action == uars_data   .OR.   &
          action == user_prog ) .AND.  &
         ( item_code_aux == fields_aux( i ) % stashmaster % item .AND. &
          sctn_code_aux == fields_aux( i ) % stashmaster % section ) ) Then

      ! Increment the copy count and check its the 1st instance.
      Copy_Count=Copy_Count + 1
      If ( Copy_Count == 1 ) Then

        If ( action == user_prog ) Then
          Call Rcf_Locate( sctn_code_out, item_code_out,               &
                           fields_out, field_count_out, pos_out)
        Else
          Call Rcf_Locate( fields_aux( i ) % stashmaster % section,    &
                           fields_aux( i ) % stashmaster % item,       &
                           fields_out, field_count_out, pos_out )
        End If

        Call Rcf_Alloc_Field( fields_aux( i ) )
        Call Rcf_Read_Field( fields_aux(i), Hdr_Aux, decomp_rcf_output )

        Call Rcf_Alloc_Field( fields_out( pos_out ) )

!---------------------------------------------------------------
! Copy fields if user prog
!---------------------------------------------------------------
        If ( action == user_prog ) Then

          ! If a land only field from ancillary, can compress onto
          ! output field. This assumes a real field only
          If (fields_aux( i ) % stashmaster % grid_type ==             &
                                ppx_atm_compressed .AND.               &
              Hdr_Aux % FixHd( FH_Dataset) == FH_Dataset_Ancil ) Then

            Do j = 1, fields_out( pos_out ) % levels
! DEPENDS ON: to_land_points
              Call To_Land_Points( fields_aux( i ) % Data(:,j),        &
                                   fields_out( pos_out ) % Data(:,j),  &
                                   local_lsm_out,                      &
                                   fields_aux( i ) % level_size,       &
                                   size )
            End Do

          Else    ! not land compressed

          Select Case ( fields_out( pos_out ) % stashmaster % data_type)
            Case (ppx_type_real)
              fields_out( pos_out ) % Data( :, : ) =                   &
                                      fields_aux( i ) % Data( :, : )

            Case (ppx_type_int)
              fields_out( pos_out ) % Data_int( :, : ) =               &
                                     fields_aux( i ) % Data_int( :, : )

            Case (ppx_type_log)
              fields_out( pos_out ) % Data_Log( :, : ) =               &
                                     fields_aux( i ) % Data_Log( :, : )

          End Select
          End If

        Else             ! Must be uars or tracers

          ! If levels don't match for tracers, issue a warning
          If ( action == tracers .AND. fields_out( pos_out ) % levels/=&
                                       fields_aux( i ) % levels) Then
            ErrorStatus = -40
            Cmessage = 'Not all tracer levels have been initialised!'
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If

          !kdcorbin, 05/10 - commenting out reading of fields
          !  as they were read in from auxillary file.
          !Call Rcf_Read_Field( fields_out( pos_out ), Hdr_Out,         &
          !                                           decomp_rcf_output)

        ! Only copy the aux levels over the top of the output levels
          level_top  = fields_out( pos_out ) % levels
          level_base = fields_out( pos_out ) % levels -                &
                       fields_aux( i ) % levels + 1

          Select Case ( fields_out( pos_out ) % stashmaster % data_type)
            Case (ppx_type_real)
            fields_out( pos_out ) % Data( :, level_base : level_top) = &
                                      fields_aux( i ) % Data( :, : )

            Case (ppx_type_int)
              fields_out(pos_out) % Data_int( :,level_base:level_top)= &
                                    fields_aux( i ) % Data_int( :, : )

            Case (ppx_type_log)
              fields_out(pos_out) % Data_Log( :,level_base:level_top)= &
                                    fields_aux( i ) % Data_Log( :, : )

          End Select

        End If

        Call Rcf_Write_Field( fields_out( pos_out ), Hdr_Out,          &
                              decomp_rcf_output )

        Call Rcf_Dealloc_Field( fields_out( pos_out ) )
        Call Rcf_Dealloc_Field( fields_aux( i ) )

      Else  ! Copy_Count /= 1, must have already copied field over

        ErrorStatus = -50
        Cmessage="Was about to overwrite user_prog data with 2nd field."
        Call EReport ( RoutineName, ErrorStatus, Cmessage )

      End If ! Copy_Count = 1
    End If
  End Do

End If

! Clean up the fields
Deallocate( fields_aux )
Nullify( fields_aux )

Return
End Subroutine Rcf_Aux_File
End Module Rcf_Aux_File_Mod
