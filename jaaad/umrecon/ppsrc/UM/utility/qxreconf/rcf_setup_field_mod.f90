
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
!*L------------------ COMDECK LOOKADD ----------------------------------
!LL
!LL Purpose : Contains information about the format
!LL           of the PP header
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   4.0  12/09/95   Change NPERIODS to LBUSER3, BRSVD1 to BULEV,
!LL                   BRSVD2 to BHULEV and definitions for BRLEV and
!LL                   BHRLEV. Corresponding changes made to STWORK1A
!LL                   and PPHEAD1A. (Andrew Brady)
!LL  4.0  12/10/95  Change item 45 from lbuser7 to model_code. RTHBarnes
!LL  5.1  17/04/00    Fixed/Free format. P.Selwood.
!LL  5.2  25/09/00    Add LBCC_xxxx variables for the compressed
!LL                   LBC LOOKUP array                  P.Burton
!LL
!LL Programming standard :
!LL
!LL Logical components covered : F092
!LL
!LL Project task :
!LL
!LL External documentation:
!LL
!LLEND -----------------------------------------------------------------
!

! Validity time
      Integer, Parameter :: LBYR   =1   ! Year
      Integer, Parameter :: LBMON  =2   ! Month
      Integer, Parameter :: LBDAT  =3   ! Day of month
      Integer, Parameter :: LBHR   =4   ! Hour
      Integer, Parameter :: LBMIN  =5   ! Minute
      Integer, Parameter :: LBDAY  =6   ! Day number

! Data time
      Integer, Parameter :: LBYRD  =7   ! Year
      Integer, Parameter :: LBMOND =8   ! Month
      Integer, Parameter :: LBDATD =9   ! Day of month
      Integer, Parameter :: LBHRD  =10  ! Hour
      Integer, Parameter :: LBMIND =11  ! Minute
      Integer, Parameter :: LBDAYD =12  ! Day number

      Integer, Parameter :: LBTIM  =13  ! Time indicator
      Integer, Parameter :: LBFT   =14  ! Forcast period (hours)
      Integer, Parameter :: LBLREC =15  ! Length of data record
      Integer, Parameter :: LBCODE =16  ! Grid type code
      Integer, Parameter :: LBHEM  =17  ! Hemisphere indicator
      Integer, Parameter :: LBROW  =18  ! Number of rows in grid
      Integer, Parameter :: LBNPT  =19  ! Number of points per row
      Integer, Parameter :: LBEXT  =20  ! Length of extra data
      Integer, Parameter :: LBPACK =21  ! Packing method indicator
      Integer, Parameter :: LBREL  =22  ! Header release number
      Integer, Parameter :: LBFC   =23  ! Field code
      Integer, Parameter :: LBCFC  =24  ! Second field code
      Integer, Parameter :: LBPROC =25  ! Processing code
      Integer, Parameter :: LBVC   =26  ! Vertical coordinate type
      Integer, Parameter :: LBRVC  =27  ! Coordinate type for reference
                                        ! level

      Integer, Parameter :: LBEXP  =28  ! Experiment number
      Integer, Parameter :: LBEGIN =29  ! Start record
      Integer, Parameter :: LBNREC =30  ! No of records-Direct access
                                        ! only
      Integer, Parameter :: LBPROJ =31  ! Met-O-8 projection number
      Integer, Parameter :: LBTYP  =32  ! Met-O-8 field type
      Integer, Parameter :: LBLEV  =33  ! Met-O-8 level code
      Integer, Parameter :: LBRSVD1=34  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD2=35  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD3=36  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD4=37  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBSRCE =38  ! =1111 to indicate following
                                        ! apply to UM
      Integer, Parameter :: DATA_TYPE =39  ! Indicator for real/int or
                                           ! timeseries
      Integer, Parameter :: NADDR  =40  ! Start address in DATA_REAL or
                                        ! DATA_INT
      Integer, Parameter :: LBUSER3=41  ! Free for user-defined function
      Integer, Parameter :: ITEM_CODE =42  !Stash code
      Integer, Parameter :: LBPLEV =43  ! Pseudo-level indicator (if
                                        ! defined)
      Integer, Parameter :: LBUSER6=44  ! Free for user-defined function
      Integer, Parameter :: MODEL_CODE =45 ! internal model identifier

      Integer, Parameter :: BULEV  =46  ! Upper level boundary
      Integer, Parameter :: BHULEV =47  ! Upper level boundary
      Integer, Parameter :: BRSVD3 =48  ! Reserved for future PP-package
                                        ! use
      Integer, Parameter :: BRSVD4 =49  ! Reserved for future PP-package
                                        ! use
      Integer, Parameter :: BDATUM =50  ! Datum value
      Integer, Parameter :: BACC   =51  ! (Packed fields) Packing
                                        ! accuracy
      Integer, Parameter :: BLEV   =52  ! Level
      Integer, Parameter :: BRLEV  =53  ! Lower level boundary
      Integer, Parameter :: BHLEV  =54  ! (Hybrid levels) A-level of
                                        ! value
      Integer, Parameter :: BHRLEV =55  ! Lower level boundary
      Integer, Parameter :: BPLAT  =56  ! Real latitude of 'pseudo'
                                        ! N Pole
      Integer, Parameter :: BPLON  =57  ! Real longitude of 'pseudo'
                                        ! N Pole
      Integer, Parameter :: BGOR   =58  ! Grid orientation
      Integer, Parameter :: BZY    =59  ! Zeroth latitude
      Integer, Parameter :: BDY    =60  ! Latitude interval
      Integer, Parameter :: BZX    =61  ! Zeroth longitude
      Integer, Parameter :: BDX    =62  ! Longitude interval
      Integer, Parameter :: BMDI   =63  ! Missing data indicator
      Integer, Parameter :: BMKS   =64  ! M,K,S scaling factor

      Integer, Parameter :: LBCC_LBYR   = 1  ! Year
      Integer, Parameter :: LBCC_LBMON  = 2  ! Month
      Integer, Parameter :: LBCC_LBDAT  = 3  ! Day of the month
      Integer, Parameter :: LBCC_LBHR   = 4  ! Hour
      Integer, Parameter :: LBCC_LBMIN  = 5  ! Minute
      Integer, Parameter :: LBCC_LBDAY  = 6  ! Day number
      Integer, Parameter :: LBCC_LBEGIN = 7  ! Start record
      Integer, Parameter :: LBCC_NADDR  = 8  ! Start address of DATA
! Mapping of MPP_LOOKUP; analogous to mapping in PP header

      Integer, Parameter :: P_NADDR=1    ! Address on local PE
      Integer, Parameter :: P_LBLREC=2   ! Local length of record

!*----------------------------------------------------------------------
! NADDR IS LOCATION IN PP-HEADER (LOOKUP) FOR START POSN OF VARIABLE
! ITEM_CODE is the location in PP header for a code defined as
!           (section number)*1000+item number
! DATA_TYPE is the location in the PP header defining data as REAL or
!           INTEGER.
! LBNPT is the location defining the number of points per row
!
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
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

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
