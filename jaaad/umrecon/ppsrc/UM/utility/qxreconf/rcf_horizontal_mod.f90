
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Performs horizontal interpolation and related tasks

Module Rcf_horizontal_Mod

!  Subroutine Rcf_horizontal - horizontal interpolation
!
! Description:
!   This module contains a wrapper subroutine for horizontal
!   interpolation. The interpolation is done level-by-level
!   but on a vertically decomposed grid. Thus we need to gather
!   data by levels onto the compute PEs and then rescatter it onto
!   the output grids. This won't be the fastest routine in the world,
!   but should be reliable and be easily verifiable for
!   bit-reproducibility etc.
!
! Method:
!   Note that land compressed fields are uncompressed and recompressed,
!   copying is done if that is all that is required, coastal
!   adjustment is performed and polar row averaging is done.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_horizontal( field_in, field_out, grid_in, grid_out )

Use Ereport_mod, Only : &
    Ereport

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Parvars_Mod, Only : &
    nproc,              &
    mype,               &
    gc_all_proc_group

Use Rcf_DecompTP_Mod, Only : &
    decomp_rcf_input,    &
    decomp_rcf_output

Use Rcf_select_weights_mod, Only : &
    Rcf_select_weights

Use Rcf_average_polar_mod, Only : &
    Rcf_average_polar

Use Rcf_H_Int_Ctl_Mod, Only : &
    Rcf_H_Int_Ctl

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal,        &
    LTimer

Use Rcf_Lsm_Mod, Only : &
    n_coastal_points,                n_land_points_unres, &
    n_sea_points_unres,              coast_index_in,      &
    coast_index_out,                 index_targ_land,     &
    index_targ_sea,                  land_unres_index,    &
    sea_unres_index,                 local_lsm_in,        &
    local_lsm_out,                   glob_lsm_out,        &
    cyclic

Use Rcf_Recon_Mod, Only : &
    Lspiral_S

Use Rcf_Set_Interp_Flags_Mod, Only : &
    interp_h_only,                   interp_v_only,       &
    interp_all,                      interp_copy,         &
    interp_no_op

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_prog_sec,    &
    stashcode_tstar,       &
    stashcode_vol_smc_wilt,&
    stashcode_vol_smc_cri, &
    stashcode_vol_smc_sat

Use Rcf_Change_Decomposition_Mod, Only : &
    Rcf_Change_Decomposition

Use Rcf_Scatter_Zonal_Field_Mod, Only : &
    Rcf_Scatter_Zonal_Field

Use Rcf_Scatter_Field_Mod, Only : &
    Rcf_Scatter_Field

Use Rcf_Gather_Field_Mod, Only : &
    Rcf_Gather_Field

Use Rcf_Gather_Zonal_Field_Mod, Only : &
    Rcf_Gather_Zonal_Field

Use Rcf_Interp_Weights_Mod  ! All of it

Implicit None

! Derived Type Arguments
Type (field_type), Target, Intent(InOut) :: field_in  ! Input data field
Type (field_type), Target, Intent(InOut) :: field_out !Output data field
Type (grid_type), Intent(In)     :: grid_in   ! Input grid sizes
Type (grid_type), Intent(In)     :: grid_out  ! Output grid sizes

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

! Local data
Integer             :: div         ! for calculation of partition
Integer             :: rem         ! for calculation of partition
Integer             :: pe          ! PE from which to gather/scatter
Integer             :: i,j         ! Looping
Integer             :: size        ! check value for lsm
Integer             :: stat        ! GCOM status
Integer             :: field_averaged ! 1/0 status of field
Integer             :: orig_h_int_method ! Stores original interp method
                                         ! since soil can use different one.
Integer             :: sea_points_unres_tmp  ! tmp stuff for coast aj
Integer             :: land_points_unres_tmp ! tmp stuff for coast aj
Integer             :: lsm_tmp( grid_out % glob_p_rows * &
                                grid_out % glob_p_row_length )
Integer             :: Index_targ_land_tmp( grid_out % glob_p_rows * &
                                          grid_out % glob_p_row_length )
Integer             :: Index_targ_sea_tmp( grid_out % glob_p_rows * &
                                          grid_out % glob_p_row_length )
Integer             :: maxdim      ! a size parameter for coast aj
Integer             :: ErrorStatus
Logical             :: averaged    ! return from polar_average
Character (Len=*), Parameter :: RoutineName = 'Rcf_horzontal'
Character (Len=80)           :: Cmessage

                                       ! Single level after scatter
Real, Allocatable   :: level_field_in( : )

                                        ! Single level after interp.
Real, Allocatable   :: level_field_out( : )

Type (field_type), Target     :: field_in_tmp
Type (field_type), Target     :: field_out_tmp
Type (field_type), Pointer    :: ptr_field_in
Type (field_type), Pointer    :: ptr_field_out

! Pointers for choice of data etc
Integer, Pointer              :: ptr_bl_index_b_l(:)
Integer, Pointer              :: ptr_bl_index_b_r(:)
Integer, Pointer              :: ptr_aw_index_targ_lhs(:)
Integer, Pointer              :: ptr_aw_index_targ_top(:)
Real, Pointer                 :: ptr_aw_colat_t(:)
Real, Pointer                 :: ptr_aw_long_l(:)
Real, Pointer                 :: ptr_weight_b_l(:)
Real, Pointer                 :: ptr_weight_b_r(:)
Real, Pointer                 :: ptr_weight_t_l(:)
Real, Pointer                 :: ptr_weight_t_r(:)

External Intf_Coast_AJ, Timer

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 3)

! If input and output datatypes are different, we have a problem.
If (  field_in % stashmaster % data_type /=                     &
     field_out % stashmaster % data_type ) Then
  Cmessage = 'Input and Output field datatypes differ!'
  ErrorStatus = 10
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! Make sure that field_averaged is initialised
field_averaged = 0

!-------------------------------------------------------------------
! Is interpolation activated? If not, copy data across and that's
! all we will do.
!-------------------------------------------------------------------
Select Case( field_in % interp )
  Case( interp_copy, interp_v_only )     ! copy data

  ! Sizes should be the same, but will check...
  If ( field_in % level_size /= field_out % level_size ) Then
    write (6,*) "Aborting due to mismatch in local datasizes "
    write (6,*) "Input dump data size =",field_in % level_size
    write (6,*) "Output dump data size =",field_out % level_size
    Cmessage = 'No interpolation required, but input and output &
               &data fields are different sizes'
    ErrorStatus = 20
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

  Select Case( field_in % stashmaster % data_type )
    Case ( ppx_type_real )
      field_out % Data(:,:) = field_in % Data(:,:)

    Case ( ppx_type_int )
      If ( Associated( field_in %  Data ) ) Then
        field_out % Data(:,:) = field_in % Data(:,:)
      Else
        field_out % Data_Int(:,:) = field_in % Data_Int(:,:)
      End If

    Case ( ppx_type_log )
      If ( Associated( field_in % Data ) ) Then
        field_out % Data(:,:) = field_in % Data(:,:)
      Else
        field_out % Data_Log(:,:) = field_in % Data_Log(:,:)
      End If

    Case Default
      Cmessage = 'Unsupported Data-Type'
      ErrorStatus = -30
      Call Ereport( RoutineName, ErrorStatus, Cmessage )

  End Select

  Case( interp_all, interp_h_only )

!---------------------------------------------------------------
! If data is compressed onto land points, need to expand it
! To do this neatly requires some fiddling with pointers.
! ptr_field_in and ptr_field_out are the ones to work with. Set
! up here.
!--------------------------------------------------------------
  If ( field_in % stashmaster % grid_type == ppx_atm_compressed) Then
    Allocate( field_in_tmp % Data( grid_in % loc_p_rows *           &
                             grid_in % loc_p_row_length ,           &
                             field_in % levels ) )

    Allocate( field_out_tmp % Data( grid_out % loc_p_rows *         &
                              grid_out % loc_p_row_length,          &
                              field_out % levels ) )

    field_in_tmp % levels       = field_in % levels
    field_in_tmp % rows         = grid_in % loc_p_rows
    field_in_tmp % row_len      = grid_in % loc_p_row_length
    field_in_tmp % level_size   = field_in_tmp % rows *              &
                                  field_in_tmp % row_len
    field_in_tmp % glob_rows    = grid_in % glob_p_rows
    field_in_tmp % glob_row_len = grid_in % glob_p_row_length
    field_in_tmp % glob_level_size = field_in_tmp % glob_rows *      &
                                     field_in_tmp % glob_row_len
    field_in_tmp % stashmaster => field_in % stashmaster

    field_out_tmp % levels       = field_out % levels
    field_out_tmp % rows         = grid_out % loc_p_rows
    field_out_tmp % row_len      = grid_out % loc_p_row_length
    field_out_tmp % level_size   = field_out_tmp % rows *              &
                                   field_out_tmp % row_len
    field_out_tmp % glob_rows    = grid_out % glob_p_rows
    field_out_tmp % glob_row_len = grid_out % glob_p_row_length
    field_out_tmp % glob_level_size = field_out_tmp % glob_rows *      &
                                      field_out_tmp % glob_row_len
    field_out_tmp % stashmaster => field_out % stashmaster

    ! Expand level by level
    Do i = 1, field_in % levels
! DEPENDS ON: from_land_points
      Call From_Land_Points( field_in_tmp % Data(1,i),                &
                             field_in % Data(1,i),                    &
                             local_lsm_in,  field_in_tmp % level_size,&
                             field_in % level_size )
    End Do

    ptr_field_in  => field_in_tmp
    ptr_field_out => field_out_tmp
  Else
    ptr_field_in  => field_in
    ptr_field_out => field_out
  End If

!------------------------------------------------------------------
! Can now allocate space for local levels for interpolation
!------------------------------------------------------------------
  Allocate( level_field_in( ptr_field_in % glob_level_size ) )
  Allocate( level_field_out( ptr_field_out % glob_level_size ) )

!-------------------------------------------------------------------
! Need to work out which weights we wish to use. This depends on
! the type of field that is being interpolated.
!-------------------------------------------------------------------

! For soil moisture fields we can have nearest neighbour interpolation.
! This is assuming that bilinear is being used as rcf_select_weights
! enforces this condition.
  orig_h_int_method = h_int_method
  If ( field_in % stashmaster % section  == stashcode_prog_sec      .AND.  &
       ( field_in % stashmaster % item   == stashcode_vol_smc_wilt  .OR.   &
         field_in % stashmaster % item   == stashcode_vol_smc_cri   .OR.   &
         field_in % stashmaster % item   == stashcode_vol_smc_sat ) .AND.  &
       smcp_int_method == nearest_neighbour ) Then
    h_int_method = smcp_int_method
  End If

  Call Rcf_select_weights( ptr_bl_index_b_l, ptr_bl_index_b_r, &
                           ptr_weight_b_l, ptr_weight_b_r,     &
                           ptr_weight_t_l, ptr_weight_t_r,     &
                           ptr_aw_index_targ_lhs,              &
                           ptr_aw_index_targ_top,              &
                           ptr_aw_colat_t, ptr_aw_long_l,      &
                           field_in % stashmaster % grid_type, &
                           field_in % stashmaster % section,   &
                           field_in % stashmaster % item )

!--------------------------------------------------------------------
! We need to gather <nproc> levels onto pes.
!--------------------------------------------------------------------

  div = field_in % levels / nproc
  rem = Mod( field_in % levels, nproc )
  pe = 0


  Do i = 1, div

    Call Rcf_Change_Decomposition( decomp_rcf_input )

    Do j = ((i-1) * nproc) + 1, i * nproc
      ! Will gather level j on pe

      If (ptr_field_in % glob_row_len == 1) Then ! Zonal Data
        Call Rcf_Gather_Zonal_Field( ptr_field_in % Data(:,j),         &
                                     level_field_in,                   &
                                     ptr_field_in % level_size,        &
                                     ptr_field_in % glob_level_size, 1,&
                                     ppx_atm_tzonal, pe )

      Else

        Call Rcf_Gather_Field( ptr_field_in % Data(:,j),              &
                               level_field_in,                        &
                               ptr_field_in % row_len,                &
                               ptr_field_in % rows,                   &
                               ptr_field_in % glob_row_len,           &
                               ptr_field_in % glob_rows, pe,          &
                               gc_all_proc_group )
      End If

      pe = pe + 1
      if (pe == nproc) pe = 0
    End Do

!-------------------------------------------------------------------
! All PEs are currently full of level data - so do the interpolation
!-------------------------------------------------------------------
    Call Gc_Ssync( nproc, ErrorStatus )

    Call Rcf_H_Int_Ctl( ptr_field_out % glob_level_size,                 &
                        ptr_field_in % glob_row_len,                     &
                        ptr_field_out % glob_row_len,                    &
                        ptr_field_in % glob_rows,                        &
                        ptr_field_out % glob_rows,                       &
                        grid_out % global, ptr_aw_index_targ_lhs,        &
                        ptr_aw_index_targ_top, ptr_bl_index_b_l,         &
                        ptr_bl_index_b_r, ptr_aw_colat_t, ptr_aw_long_l, &
                        level_field_in, ptr_weight_t_r, ptr_weight_b_r,  &
                        ptr_weight_t_l, ptr_weight_b_l, level_field_out )
  
!-------------------------------------------------------------------
! Coastal Adjustment for land-only, sea-only and T*  fields
!-------------------------------------------------------------------
! DEPENDS ON: timer
    If (LTimer) Call Timer( 'Coastal adjustment', 103)
  
    If (ptr_field_out % stashmaster % grid_type == ppx_atm_compressed.OR.   &
        ptr_field_out % stashmaster % grid_type == ppx_atm_tsea .OR.        &
        ( ptr_field_out % stashmaster % section == stashcode_prog_sec .AND. &
          ptr_field_out % stashmaster % item    == stashcode_tstar )) Then
      Do j = 1, n_coastal_points
        level_field_out( coast_index_out( j ) ) = &
                         level_field_in( coast_index_in( j ) )
      End Do
  
      If (LSpiral_S) Then     ! Spiral adjustment
        maxdim = min(grid_out % glob_p_rows, grid_out % glob_p_row_length)
        Do j = 1, grid_out % glob_p_rows * grid_out % glob_p_row_length
          index_targ_sea_tmp( j )  = index_targ_sea( j )
          index_targ_land_tmp( j ) = index_targ_land( j )
          If ( glob_lsm_out( j ) ) Then
            lsm_tmp( j ) = 1
          Else
            lsm_tmp( j ) = 0
          End If
        End Do
  
        sea_points_unres_tmp  = n_sea_points_unres
        land_points_unres_tmp = n_land_points_unres
  
        If (ptr_field_out % stashmaster % grid_type /=                   &
                                          ppx_atm_compressed) Then
! Only do the sea coastal points if the field isn't land only.
! DEPENDS ON: intf_coast_aj
          Call Intf_Coast_AJ( lsm_tmp, index_targ_sea_tmp,               &
               sea_points_unres_tmp, grid_out % glob_p_rows,             &
               grid_out % glob_p_row_length, level_field_out, 0, cyclic, &
               maxdim )
        End If

! DEPENDS ON: intf_coast_aj
        Call Intf_Coast_AJ( lsm_tmp, index_targ_land_tmp,              &
             land_points_unres_tmp, grid_out % glob_p_rows,            &
             grid_out % glob_p_row_length, level_field_out, 1, cyclic, &
             maxdim )

      Else    ! Non-spiral adjustment

        Do j = 1, n_land_points_unres
          level_field_out( index_targ_land( j ) ) = &
                           level_field_out( land_unres_index( j ) )
        End Do

        Do j = 1, n_sea_points_unres
          level_field_out( index_targ_sea( j ) ) = &
                         level_field_out( sea_unres_index( j ) )
        End Do

      End If ! Spiral
    End If

! DEPENDS ON: timer
    If (LTimer) Call Timer( 'Coastal adjustment', 104)

!-----------------------------------------------------------------
! Average the polar rows
!-----------------------------------------------------------------
    If ( ptr_field_out % stashmaster % grid_type <= 3 ) Then
      Call Rcf_Average_polar(level_field_out, ptr_field_out % glob_rows, &
                             ptr_field_out % glob_row_len,               &
                             grid_out % global, averaged )

      If ( averaged .AND. PrintStatus >= PrStatus_Normal) Then
        field_averaged = 1
      End If
    End If

!-------------------------------------------------------------------
! And re-scatter the data back to original PEs
!------------------------------------------------------------------
    Call Rcf_Change_Decomposition( decomp_rcf_output )
    Call Gc_Ssync( nproc, ErrorStatus)
    Do j = ((i-1) * nproc) + 1, i * nproc

      If (ptr_field_out % glob_row_len == 1) Then ! Zonal Data
        Call Rcf_Scatter_Zonal_Field( ptr_field_out % Data(:,j),        &
                                  level_field_out,                      &
                                  ptr_field_out % level_size,           &
                                  ptr_field_out % glob_level_size, 1,   &
                                  ppx_atm_tzonal, pe )
      Else
        Call Rcf_Scatter_Field( ptr_field_out % Data(:,j),              &
                                level_field_out,                        &
                                ptr_field_out % row_len,                &
                                ptr_field_out % rows,                   &
                                ptr_field_out % glob_row_len,           &
                                ptr_field_out % glob_rows, pe,          &
                                gc_all_proc_group )
      End If
      pe = pe + 1
      If (pe == nproc) pe = 0
    End Do
  End Do

!-------------------------------------------------------------------
! There are rem levels left to process. Will do these now.
!-------------------------------------------------------------------
  Call Rcf_Change_Decomposition( decomp_rcf_input )
  pe = 0
  Do i = 1, rem
    j = nproc * div + i

    If (ptr_field_in % glob_row_len == 1) Then ! Zonal Data
      Call Rcf_Gather_Zonal_Field( ptr_field_in % Data(:,j),            &
                                   level_field_in,                      &
                                   ptr_field_in % level_size,           &
                                   ptr_field_in % glob_level_size, 1,   &
                                   ppx_atm_tzonal, pe )
    Else
      Call Rcf_Gather_Field( ptr_field_in % Data(:,j),                  &
                             level_field_in,                            &
                             ptr_field_in % row_len,                    &
                             ptr_field_in % rows,                       &
                             ptr_field_in % glob_row_len,               &
                             ptr_field_in % glob_rows, pe,              &
                             gc_all_proc_group )
    End If

    pe = pe + 1
  End Do

  Call Gc_Ssync( nproc, ErrorStatus )

  If (mype < pe) Then
    Call Rcf_H_Int_Ctl( ptr_field_out % glob_level_size,                 &
                        ptr_field_in % glob_row_len,                     &
                        ptr_field_out % glob_row_len,                    &
                        ptr_field_in % glob_rows,                        &
                        ptr_field_out % glob_rows,                       &
                        grid_out % global, ptr_aw_index_targ_lhs,        &
                        ptr_aw_index_targ_top, ptr_bl_index_b_l,         &
                        ptr_bl_index_b_r, ptr_aw_colat_t, ptr_aw_long_l, &
                        level_field_in, ptr_weight_t_r, ptr_weight_b_r,  &
                        ptr_weight_t_l, ptr_weight_b_l, level_field_out )

!-------------------------------------------------------------------
! Coastal Adjustment for land-only, sea-only and T*  fields
!-------------------------------------------------------------------
! DEPENDS ON: timer
    If (LTimer) Call Timer( 'Coastal adjustment', 103)
    If (ptr_field_out % stashmaster % grid_type == ppx_atm_compressed.OR.   &
        ptr_field_out % stashmaster % grid_type == ppx_atm_tsea .OR.        &
        ( ptr_field_out % stashmaster % section == stashcode_prog_sec .AND. &
          ptr_field_out % stashmaster % item    == stashcode_tstar )) Then
      Do j = 1, n_coastal_points
        level_field_out( coast_index_out( j ) ) = &
                         level_field_in( coast_index_in( j ) )
      End Do

      If (LSpiral_S) Then     ! Spiral adjustment
        maxdim = min( grid_out % glob_p_rows,grid_out % glob_p_row_length)
        Do j = 1, grid_out % glob_p_rows * grid_out % glob_p_row_length
          index_targ_sea_tmp( j )  = index_targ_sea( j )
          index_targ_land_tmp( j ) = index_targ_land( j )
          If ( glob_lsm_out( j ) ) Then
            lsm_tmp( j ) = 1
          Else
            lsm_tmp( j ) = 0
          End If
        End Do

        sea_points_unres_tmp  = n_sea_points_unres
        land_points_unres_tmp = n_land_points_unres

        If (ptr_field_out % stashmaster % grid_type /=                   &
                                          ppx_atm_compressed) Then
! Only do the sea coastal points if the field isn't land only.
! DEPENDS ON: intf_coast_aj
          Call Intf_Coast_AJ( lsm_tmp, index_targ_sea_tmp,               &
               sea_points_unres_tmp, grid_out % glob_p_rows,             &
               grid_out % glob_p_row_length, level_field_out, 0, cyclic, &
               maxdim )
        End If

! DEPENDS ON: intf_coast_aj
        Call Intf_Coast_AJ( lsm_tmp, index_targ_land_tmp,              &
             land_points_unres_tmp, grid_out % glob_p_rows,            &
             grid_out % glob_p_row_length, level_field_out, 1, cyclic, &
             maxdim )

      Else    ! Non-spiral adjustment

        Do j = 1, n_land_points_unres
          level_field_out( index_targ_land( j ) ) = &
                           level_field_out( land_unres_index( j ) )
        End Do

        Do j = 1, n_sea_points_unres
          level_field_out( index_targ_sea( j ) ) = &
                           level_field_out( sea_unres_index( j ) )
        End Do

      End If ! Spiral
    End If

! DEPENDS ON: timer
    If (LTimer) Call Timer( 'Coastal adjustment', 104)

!-----------------------------------------------------------------
! Average the polar rows
!-----------------------------------------------------------------
    If ( ptr_field_out % stashmaster % grid_type <= 3 ) Then
      Call Rcf_Average_polar( level_field_out, ptr_field_out % glob_rows,&
                       ptr_field_out % glob_row_len,                     &
                       grid_out % global, averaged )

      If ( averaged .AND. PrintStatus >= PrStatus_Normal) Then
        field_averaged = 1
      End If
    End If

  End If

!-------------------------------------------------------------------
! And re-scatter data
!-------------------------------------------------------------------

  Call Rcf_Change_Decomposition( decomp_rcf_output )
  Call Gc_Ssync( nproc, ErrorStatus )
  pe = 0
  Do i = 1, rem
    j = nproc * div + i
    If (ptr_field_out % glob_row_len == 1) Then ! Zonal Data
      Call Rcf_Scatter_Zonal_Field( ptr_field_out % Data(:,j),        &
                                level_field_out,                      &
                                ptr_field_out % level_size,           &
                                ptr_field_out % glob_level_size, 1,   &
                                ppx_atm_tzonal , pe )
    Else
      Call Rcf_Scatter_Field( ptr_field_out % Data(:,j),              &
                              level_field_out,                        &
                              ptr_field_out % row_len,                &
                              ptr_field_out % rows,                   &
                              ptr_field_out % glob_row_len,           &
                              ptr_field_out % glob_rows, pe,          &
                              gc_all_proc_group )
    End If
    pe = pe + 1
  End Do

!----------------------------------------------------------------
! Print out a message if required for polar averaged fields
!----------------------------------------------------------------

  Call GC_IMAX( 1, nproc, stat, field_averaged )

  If ( field_averaged == 1 .AND. mype == 0 .AND.                       &
       PrintStatus >= PrStatus_Normal ) Then
    Write (6,*) 'Interpolation has required averaging of polar rows ', &
                'for section: ', ptr_field_out % stashmaster % section,&
                ', item: ', ptr_field_out % stashmaster % item
  End If

!----------------------------------------------------------------
! Deallocate levels
!----------------------------------------------------------------
  Deallocate( level_field_in )
  Deallocate( level_field_out )

!-----------------------------------------------------------------
! Reverse process above, take fields back to land points
!-----------------------------------------------------------------

  If (field_out % stashmaster % grid_type == ppx_atm_compressed) Then
    Do i = 1, field_out % levels
! DEPENDS ON: to_land_points
      Call To_Land_Points( field_out_tmp % Data(1,i),                 &
                           field_out % Data(1,i),                     &
                           local_lsm_out, field_out_tmp % level_size, &
                           size )

      If (size /= field_out % level_size ) Then
        Cmessage = 'Recompression onto land points - sizes mismatch'
        ErrorStatus = 60
        Call Ereport( RoutineName, ErrorStatus, Cmessage )
      End If

    End Do

  ! Release temp. memory
    Deallocate( field_out_tmp % Data )
    Deallocate( field_in_tmp % Data )
  End If

  Nullify( ptr_field_out )
  Nullify( ptr_field_in  )

 ! Earlier check if using soil could have changed this.  Set it back for
 ! future fields. 
  h_int_method = orig_h_int_method

  Case( interp_no_op)
    ! do nothing

End Select

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 4)

Return
End Subroutine Rcf_horizontal

End Module Rcf_horizontal_Mod
