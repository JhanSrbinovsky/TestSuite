
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! + Scatters any type of field from one processor to many processors

Module Rcf_General_Scatter_Field_Mod

!  Subroutine Rcf_General_Scatter_Field
!
! Description:
!   Chooses scattering method for a field based on grid_code.
!
! Method:
!   Land compressed fields are uncompressed, scattered and recompressed
!   LBCs and zonal fields have specialist subroutines
!   Other fields are "normal" and can be dealt with by Rcf_Scatter_Field
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.5   06/01/03   River routing support. P.Selwood.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.



Contains

SUBROUTINE Rcf_General_Scatter_Field ( LOCAL_FIELD,  GLOBAL_FIELD , &
                                       LOCAL_SIZE,   GLOBAL_SIZE ,  &
                                       STASH_RECORD, SCATTER_PE )

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Lsm_Mod, Only :   &
    Local_Land_Field,     &
    Local_Atmos_Landmask, &
    Glob_Atmos_Landmask

Use Rcf_Parvars_Mod, Only :                    &
    blsizep,                blsizeu,           &
    blsizev,                glsizeu,           &
    glsizev,                glsize,            &
    blsizer,                glsizer,           &
    lasize,                 gc_all_proc_group, &
    fld_type_p,             fld_type_u,        &
    fld_type_v,             fld_type_r,        &
    mype

Use Rcf_Scatter_Atmos_LBCs_Mod, Only : &
    Rcf_Scatter_Atmos_LBCs

Use Rcf_Scatter_Ocean_LBCs_Mod, Only : &
    Rcf_Scatter_Ocean_LBCs

Use Rcf_Scatter_Zonal_Field_Mod, Only : &
    Rcf_Scatter_Zonal_Field

Use Rcf_Scatter_Field_Mod, Only : &
    Rcf_Scatter_Field

Use Rcf_Global_To_Local_Mod, Only : &
    Rcf_Get_Fld_Type

Use Rcf_Ppx_Info_Mod, Only : &
    STM_Record_Type

IMPLICIT NONE

! Subroutine arguments:
Integer, Intent(In)   :: Global_Size  ! size of global_field
Integer, Intent(In)   :: Scatter_PE   ! PE on which global field lives
Integer, Intent(Out)  :: Local_Size   ! size of local_field

Real, Intent(In)      :: Global_Field( Global_Size ) ! field to scatter
Real, Intent(Out)     :: Local_Field( * )

Type( STM_Record_Type ), Intent(In) :: Stash_Record

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
Integer               :: dummy       ! ignored argument
Integer               :: fld_type    ! P or U or V
Integer               :: local_x     ! local x size
Integer               :: local_y     ! local y size
Integer               :: global_x    ! global x size
Integer               :: global_y    ! global y size
Integer               :: ErrorStatus ! error code
Character (Len=*), Parameter   :: RoutineName='Rcf_General_Scatter_&
                                              &Field'
Character (Len=80)             :: Cmessage

Real                  :: buf_expand( glsize(1) * glsize(2) )
Real                  :: buf_expand_local( lasize(1) * lasize(2) )

!===================================================================

! Choose action depending on grid code
Select Case( Stash_Record % grid_type )

!-------------------------------------------------------------------
! Land compressed fields
!-------------------------------------------------------------------
  Case( ppx_atm_compressed )

    IF (mype .EQ. SCATTER_PE) THEN

! DEPENDS ON: from_land_points
      CALL FROM_LAND_POINTS(buf_expand, GLOBAL_FIELD, &
                            glob_atmos_landmask,      &
                            glsize(1)*glsize(2),dummy)

! SCATTER_PE now contains the expanded version of the full field

    ENDIF

! Now scatter this to all the other processors, putting the local
! part of the field into the array buf_expand_local

    CALL Rcf_Scatter_Field(buf_expand_local , buf_expand,       &
                           lasize(1) ,        lasize(2),        &
                           glsize(1) ,        glsize(2),        &
                           SCATTER_PE,        GC_ALL_PROC_GROUP )

! Pack the local field down to local land points and put
! the packed field into LOCAL_FIELD

! DEPENDS ON: to_land_points
    CALL TO_LAND_POINTS(buf_expand_local,LOCAL_FIELD, &
                        local_atmos_landmask,         &
                        lasize(1)*lasize(2), dummy)

    local_size = local_land_field

!-------------------------------------------------------------------
! Atmosphere Lateral boundary fields
!-------------------------------------------------------------------
  Case ( ppx_atm_lbc_theta, ppx_atm_lbc_u, ppx_atm_lbc_v )

    CALL Rcf_Scatter_atmos_LBCS( Global_Field, Global_Size, &
                                 Local_Field,  Local_Size,  &
                                 Stash_Record, Scatter_PE)


!-------------------------------------------------------------------
! Ocean Lateral boundary fields
!-------------------------------------------------------------------
  Case (ppx_ocn_rim)

    CALL Rcf_Scatter_Ocean_LBCS( Global_Field, Global_Size, &
                                 Local_Field,  Local_Size,  &
                                 Stash_Record, Scatter_PE)

!-------------------------------------------------------------------
! Zonal fields
!-------------------------------------------------------------------
  Case ( ppx_atm_tzonal, ppx_atm_uzonal, &
         ppx_ocn_tzonal, ppx_ocn_uzonal )

    If ( Stash_Record % grid_type == ppx_atm_tzonal .OR. &
         Stash_Record % grid_type == ppx_ocn_tzonal ) Then    ! P grid
      global_y = glsize(2)
      local_y  = blsizep(2)

    Else                ! U grid
      global_y = glsizeu(2)
      local_y  = blsizeu(2)
    End If

    local_size = local_y

    Call Rcf_Scatter_Zonal_Field( Local_Field, Global_Field,           &
                                  local_y,     global_y,               &
                                  1,          Stash_Record % grid_type,&
                                  Scatter_PE )

!-------------------------------------------------------------------
! Normal fields
!-------------------------------------------------------------------
  Case &
!     atmosphere grids
       ( ppx_atm_tall,   &! Atmos T points
         ppx_atm_tland,  &! Atmos T land points
         ppx_atm_tsea,   &! Atmos T sea points
         ppx_atm_uall,   &! Atmos U points
         ppx_atm_uland,  &! Atmos U land points
         ppx_atm_usea,   &! Atmos U sea points
         ppx_atm_cuall,  &! Atmos C grid U pts
         ppx_atm_cvall,  &! Atmos C grid V pts
         ppx_atm_ozone,  &! Atmos ozone field
         ppx_atm_river,  &! Atmos river routing field
!     ocean grids
         ppx_ocn_tcomp,  &! Ocean "Compressed" T
         ppx_ocn_ucomp,  &! Ocean "Compressed" u
         ppx_ocn_tall,   &! Ocean T points (cyc)
         ppx_ocn_uall,   &! Ocean U points (cyc)
         ppx_ocn_cuall,  &! Ocean C grid U pts
         ppx_ocn_cvall,  &! Ocean C grid V pts
         ppx_ocn_tfield, &! Ocean T points
         ppx_ocn_ufield)! Ocean U points


    fld_type = Rcf_Get_Fld_Type(Stash_Record % grid_type)

    IF (fld_type .EQ. fld_type_p) THEN
      global_x = glsize(1)
      global_y = glsize(2)
      local_x  = blsizep(1)
      local_y  = blsizep(2)
    ELSE IF (fld_type .EQ. fld_type_u) THEN
      global_x = glsizeu(1)
      global_y = glsizeu(2)
      local_x  = blsizeu(1)
      local_y  = blsizeu(2)
    ELSE IF (fld_type .EQ. fld_type_v) THEN
      global_x = glsizev(1)
      global_y = glsizev(2)
      local_x  = blsizev(1)
      local_y  = blsizev(2)
    ELSE IF (fld_type .EQ. fld_type_r) THEN
      global_x = glsizer(1)
      global_y = glsizer(2)
      local_x  = blsizer(1)
      local_y  = blsizer(2)
    ENDIF

    local_size = local_x * local_y

    Call Rcf_Scatter_Field( Local_Field, Global_Field,     &
                            local_x,     local_y,          &
                            global_x,    global_y,         &
                            Scatter_PE,  GC_ALL_PROC_GROUP )

!-------------------------------------------------------------------
! Any other type of field
!-------------------------------------------------------------------
  Case Default

    ErrorStatus = 10
    Cmessage='Field type not recognized for Scatter'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )

End Select

RETURN
END Subroutine Rcf_General_Scatter_Field
End Module Rcf_General_Scatter_Field_Mod