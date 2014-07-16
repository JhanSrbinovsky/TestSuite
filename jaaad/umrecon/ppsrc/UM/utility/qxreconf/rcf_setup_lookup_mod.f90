
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Sets up the output dump lookup tables

Module Rcf_Setup_Lookup_Mod

!  Subroutine Rcf_Setup_Lookup - sets up lookups for output dump
!
! Description:
!   The lookup headers are filled in - with calculated addressing
!   and other namelist derived variables
!
! Method:
!    UMDP F3 defines the lookups.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Setup_Lookup( Hdr_In, Hdr_Out )

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Exppx_mod, Only : &
    Rcf_Exppx

Use Rcf_Ppx_Info_mod, Only : &
    STM_Record_Type,         &
    NSectP

Use Rcf_Submodel_Mod, Only : &
    Internal_model_list, &
    N_Internal_model

Use Rcf_NRecon_Mod, Only : &
    Recondat_Node,         &
    RecondatList,          &
    DumpProgLevs

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Model_Mod, Only : &
    Necf

Use Rcf_HeadAddress_Mod, Only :&
    FH_HorizGrid,         RC_LatSpacing,        RC_LongSpacing, &
    RC_PoleLat,           RC_PoleLong,          FH_VertCoord,   &
    FH_VertCoord_CP,      RC_FirstLat,          LDC_ZseaTheta,  &
    LDC_CkTheta,          LDC_ZseaRho,          LDC_CkRho

Use Rcf_generate_heights_mod, Only : &
    height_gen_original,             &
    height_gen_smooth

Use Rcf_Readnl_horizont_Mod, Only : &
    Iproj

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

Use Rcf_Items_Mod, Only : &
    num_items,          area_array,             item_array

Use Rcf_Recon_Mod, Only : &
    Lcal360,          &
    Lozone_zonal,     &
    Dump_Pack,        &
    Rimwidtha,        &
    Var_Recon

Use Rcf_Level_Code_Mod, Only : &
    Rcf_Level_Code

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_prog_sec,        &
    stashcode_ozone,           &
    stashcode_ozone_tracer,    &
    stashcode_o3_prod_loss,    &
    stashcode_o3_p_l_vmr,      &
    stashcode_o3_vmr,          &
    stashcode_o3_p_l_temp,     &
    stashcode_o3_temp,         &
    stashcode_o3_p_l_colo3,    &
    stashcode_o3_colo3,        &
    stashcode_u

Implicit None

! Arguments
Type (Um_Header_type), Intent(In) :: Hdr_In
Type (Um_Header_type), Target     :: Hdr_Out

! Comdecks
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
!
! Description:
!   Contains model id used in LBSRCE in UM lookup tables
!   1111 is used to define unified model
!
! Current Code Owner: D.M. Goddard
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   5.3  30/01/02   Original code. D.M. Goddard
!
! Declarations:

! Global parameters:
      INTEGER, PARAMETER :: model_id = 1111

!- End of COMDECK declaration

! Local vars.
Integer                :: int_val = 1      ! used for transfers
Integer, Pointer       :: Lookup(:,:)
Type (STM_Record_Type) :: STM_Record
Integer                :: i, j, jj, k ! Loopers
Integer                :: icount     ! Counter for calc. Lookup(40)
Integer                :: sec_item
Integer                :: item
Integer                :: section
Integer                :: model
Integer                :: k_out
Integer                :: whole      ! Whole packing indicator -
                                     ! Lookup(21) - see UMDP F3
Integer                :: n1, n2, n3 ! Packing and compression parts
                                     ! of above
Integer                :: length
Integer                :: start_address
Integer                :: n_levels
Integer                :: n_plevels
Integer                :: lev          ! level
Integer                :: bot_lev      ! bottom level for field
Integer                :: lblev_val    ! value for lblev
Integer                :: area
Integer                :: ppxref_grid_type
Integer                :: ErrorStatus
Integer                :: Area_Expand( Hdr_Out % Len2Lookup )
Character (Len=*), Parameter :: RoutineName='Rcf_Setup_Lookup'
Character (Len=80)     :: Cmessage
Real                   :: depth
Real                   :: level( Hdr_Out % Len1LevDepC )
Real                   :: zsea  !} height defining constants
Real                   :: Ck    !}
Real                   :: riv_bzx
Real                   :: riv_bzy

! Function & Subroutine calls:
Integer            :: get_um_version_id

!------------------------------------------------------------------
! First set all elements of lookup to 0 - except dates which
! come from the input lookup (there must be 1).  If reconfiguring
! for VAR then need to make sure date and time is from an 
! instantaneous field (problem when time accumulated fields are at
! the beginning) therefore find where u velocity is.
!------------------------------------------------------------------
Lookup => Hdr_Out % Lookup
k = 1
If (Var_Recon) Then
  Do i = 1, Hdr_In % Len2Lookup
    If (Hdr_In % Lookup(item_code, i) == stashcode_u) Then
      k = i
      Exit
    End If
  End Do
End If

Do i = lbyr, lbdayd
  Lookup( i, : ) = Hdr_In % Lookup( i, k )
End Do
Lookup( 13 : 45, : )  = 0
Lookup( 46 : 64, : ) = Transfer( 0.0, int_val)

K_OUT=1

! 5: Loop through prognostic items and initialise Lookup
Do J=1,N_INTERNAL_MODEL
  Do JJ=0,NSectP
    recondat_node => RecondatList(J,JJ)
    Do While (Associated(recondat_node % recondat_info)) 
      ! 5.1: Extract addressing and number of levels from linked list
      SEC_ITEM      = recondat_node % recondat_info % sec_item
      N_LEVELS      = recondat_node % recondat_info % rlevs
      LENGTH        = recondat_node % recondat_info % len
      START_ADDRESS = recondat_node % recondat_info % raddress
      N_PLEVELS     = recondat_node % recondat_info % rplevs

      area = 1
      Do k = 1, num_items
        If ( sec_item == item_array(k) ) Then
          area = area_array(k)
        End If
      End Do

      ICOUNT=0
      Do  K=K_OUT,K_OUT+N_LEVELS*N_PLEVELS-1

        area_expand(k) = area

        ! 5.3.1: Initialise Lookup
        Lookup( item_code, K ) = SEC_ITEM
        Lookup( model_code,K )=INTERNAL_MODEL_LIST(J)

        item    = Mod(sec_item,1000)
        section = (Lookup(item_code,K) - item)/1000
        model   = internal_model_list(j)
        STM_record = Rcf_Exppx( model, section, item )

        ! 5.3.2: Set addressing information
        Lookup( lblrec,K)=LENGTH/(N_LEVELS*N_PLEVELS)

        Lookup( naddr, K)=START_ADDRESS+ICOUNT
        ICOUNT=ICOUNT+Lookup( lblrec,K)

        !-----------------------------------------------------------
        !        Calculate levels from level dependent constants
        !        Set levels for multi-level fields only
        !        Levels not set for single level fields
        !-----------------------------------------------------------
        If (N_LEVELS > 1) Then

          ! Set the level
          Call Rcf_Level_Code( STM_record % lb_code, bot_lev,          &
          Output_Grid)

          ! Set level number. For ozone we only use top of the
          ! atmosphere levels. Otherwise we count from the surface.
          lev = Mod(K - K_OUT, N_levels) + bot_lev

          If (lev == 0) Then
            lblev_val = 9999                      ! surface
          Else
            lblev_val = lev
          End If

          If ( section == stashcode_prog_sec .AND.                     & 
               item    == stashcode_ozone) Then
            lev = lev + Output_Grid % model_levels -                   &
            Output_Grid % ozone_levels
          End If
          
          LOOKUP( lblev,K ) = lblev_val

          If ( Output_Grid % height_gen_method == height_gen_smooth) Then

            If (STM_record % lbvc_code == 9 .OR.                  &
            STM_record % lbvc_code ==65 ) Then  ! Hybrid/Eta levels

              If (STM_record % lv_code == ppx_theta_level) Then

                If (lev == Output_Grid % model_levels ) Then ! Top level

                  ! Level is current theta level
                  zsea = Hdr_Out % LevDepC( lev+1, LDC_ZseaTheta )
                  Ck   = Hdr_Out % LevDepC( lev+1, LDC_CkTheta )
                  Lookup( blev, K )  = Transfer( zsea, int_val )
                  Lookup( bhlev, K ) = Transfer( Ck, int_val )

                  ! Lower boundary is rho level below
                  zsea = Hdr_Out % LevDepC( lev, LDC_ZseaRho )
                  Ck   = Hdr_Out % LevDepC( lev, LDC_CkRho )
                  Lookup( brlev, K )  = Transfer( zsea, int_val )
                  Lookup( bhrlev, K ) = Transfer( Ck, int_val )

                  ! Upper boundary is rho level above - calculated
                  zsea  = 2.0 * Hdr_Out % LevDepC( lev+1, LDC_ZseaTheta) &
                  - Hdr_Out % LevDepC( lev, LDC_ZseaRho )
                  Ck    = 0.0    ! above 1st constant rho level
                  Lookup( bulev, K )  = Transfer( zsea, int_val )
                  Lookup( bhulev, K ) = Transfer( Ck, int_val )

                Else If ( lev > Output_Grid % model_levels ) Then

                  Lookup( bulev, lev )  = Transfer( 0., int_val )
                  Lookup( blev, lev )   = Transfer( 0., int_val )
                  Lookup( brlev, lev )  = Transfer( 0., int_val )
                  Lookup( bhulev, lev ) = Transfer( 0., int_val )
                  Lookup( bhlev, lev )  = Transfer( 0., int_val )
                  Lookup( bhrlev, lev ) = Transfer( 0., int_val )

                Else

                  ! Level is current theta level
                  zsea = Hdr_Out % LevDepC( lev+1, LDC_ZseaTheta )
                  Ck   = Hdr_Out % LevDepC( lev+1, LDC_CkTheta )
                  Lookup( blev, K )  = Transfer( zsea, int_val )
                  Lookup( bhlev, K ) = Transfer( Ck, int_val )

                  ! Upper boundary is rho level above
                  zsea = Hdr_Out % LevDepC( lev+1, LDC_ZseaRho )
                  Ck   = Hdr_Out % LevDepC( lev+1, LDC_CkRho )
                  Lookup( bulev, K )  = Transfer( zsea, int_val )
                  Lookup( bhulev, K ) = Transfer( Ck, int_val )

                  ! Set the lower boundary
                  If (lev <= 1) Then
                    ! Lowest theta level has a lower boundary of the
                    ! physical surface, as defined in the physics.
                    zsea = 0.0
                    Ck   = 1.0
                  Else
                    ! Lower boundary is rho level below
                    zsea = Hdr_Out % LevDepC( lev, LDC_ZseaRho )
                    Ck   = Hdr_Out % LevDepC( lev, LDC_CkRho )
                  End If
                  Lookup( brlev, K )  = Transfer( zsea, int_val )
                  Lookup( bhrlev, K ) = Transfer( Ck, int_val )
                End If

              Else If (STM_record % lv_code == ppx_rho_level ) Then
                If ( lev == Output_Grid % model_levels + 1) Then

                  ! Level is top rho level - calculate
                  zsea = 2.0 * Hdr_Out % LevDepC( lev, LDC_ZseaTheta ) - &
                  Hdr_Out % LevDepC( lev - 1, LDC_ZseaRho )
                  Ck   = 0.0      ! Above 1st const rho by definition
                  Lookup( blev, K )  = Transfer( zsea, int_val )
                  Lookup( bhlev, K ) = Transfer( Ck, int_val )

                  ! Upper level boundary same as level
                  Lookup( bulev, K )  = Transfer( zsea, int_val )
                  Lookup( bhulev, K ) = Transfer( Ck, int_val )

                  ! Lower boundary is top theta level
                  zsea = Hdr_Out % LevDepC( lev, LDC_ZseaTheta )
                  Ck   = Hdr_Out % LevDepC( lev, LDC_CkTheta )
                  Lookup( brlev, K )  = Transfer( zsea, int_val )
                  Lookup( bhrlev, K ) = Transfer( Ck, int_val )

                Else

                  ! Level above is Theta level
                  zsea = Hdr_Out % LevDepC( lev+1, LDC_ZSeaTheta )
                  Ck   = Hdr_Out % LevDepC( lev+1, LDC_CkTheta )
                  Lookup( bulev, K )  = Transfer( zsea, int_val )
                  Lookup( bhulev, K ) = Transfer( Ck, int_val )

                  ! Level is Rho level
                  zsea = Hdr_Out % LevDepC( lev, LDC_ZSeaRho )
                  Ck   = Hdr_Out % LevDepC( lev, LDC_CkRho )
                  Lookup( blev, K )  = Transfer( zsea, int_val )
                  Lookup( bhlev, K ) = Transfer( Ck, int_val )

                  ! Level below is Theta level
                  If (lev <=0) Then  ! orography
                    zsea = 0.0
                    Ck   = 1.0
                  Else
                    zsea = Hdr_Out % LevDepC( lev, LDC_ZSeaTheta )
                    Ck   = Hdr_Out % LevDepC( lev, LDC_CkTheta )
                  End If
                  Lookup( brlev, K )  = Transfer( zsea, int_val )
                  Lookup( bhrlev, K ) = Transfer( Ck, int_val )

                End If    ! number of levels
              End If      ! rho levels
            End If        ! Hybrid Levels

          Else If (Output_Grid % height_gen_method ==                  &
          height_gen_original) Then

            ! Alternative lbvc code=65 enabled in case adopted in future
            If (STM_record % lbvc_code == 9 .or. &
            STM_record % lbvc_code ==65 )  Then   ! Hybrid/Eta levels
              Lookup( bhlev, K )  = Transfer( RMDI, int_val )
              Lookup( bhulev, K ) = Transfer( RMDI, int_val )
              Lookup( bhrlev, K ) = Transfer( RMDI, int_val )

              If (STM_record % lv_code == ppx_theta_level) Then

                If ( lev == Output_Grid % model_levels ) Then ! Top level
                  Lookup( bulev, K ) = Transfer(                           &
                  2.0 * Output_Grid % eta_theta_levels( lev ) -    &
                  Output_Grid % eta_rho_levels( lev ), int_val )
                  Lookup( blev, K ) = Transfer(                            &
                  Output_Grid % eta_theta_levels( Lev), int_val )
                  Lookup( brlev, K ) = Transfer(                           &
                  Output_Grid % eta_rho_levels  ( Lev), int_val )
                Else If ( lev > Output_Grid % model_levels ) Then
                  Lookup( bulev, K ) = Transfer( 0., int_val )
                  Lookup( blev, K  ) = Transfer( 0., int_val )
                  Lookup( brlev,K )  = Transfer( 0., int_val )
                Else
                  Lookup( bulev, K ) = Transfer(                           &
                  Output_Grid % eta_rho_levels( Lev + 1), int_val )
                  Lookup( blev, K ) = Transfer(                            &
                  Output_Grid % eta_theta_levels( Lev), int_val )

                  If (lev <= 0) Then
                    Lookup( brlev, K ) = Transfer( 0., int_val )
                  Else
                    Lookup( brlev, K ) = Transfer(                         &
                    Output_Grid % eta_rho_levels  ( Lev), int_val )
                  End If
                End If

              Else If ( STM_record % lv_code == ppx_rho_level ) Then

                If ( lev == Output_Grid % model_levels + 1) Then
                  Lookup( bulev, K ) = Transfer(                           &
                  2.0 * Output_Grid % eta_theta_levels(lev - 1) -  &
                  Output_Grid % eta_rho_levels( lev - 1), int_val )
                  Lookup( blev, K  ) = Lookup( bulev, K )
                  Lookup( brlev,K )  = Transfer(                           &
                  Output_Grid % eta_theta_levels( lev - 1), int_val)
                Else
                  Lookup( bulev, K ) = Transfer(                           &
                  Output_Grid % eta_theta_levels( Lev), int_val )
                  Lookup( blev, K ) = Transfer(                            &
                  Output_Grid % eta_rho_levels  ( Lev), int_val )

                  If ( lev <= 0 ) Then   ! bottom level
                    Lookup( brlev, K ) = Transfer( 1.0, int_val )
                  Else
                    Lookup( brlev, K ) = Transfer(                         &
                    Output_Grid % eta_theta_levels( Lev - 1), int_val)
                  End If
                End If
              Else
                Write (6,*) 'ERROR:- lv_code=',STM_record % lv_code
                Write (6,*) 'model, section, item = ', model, section, item
                ErrorStatus = 10
                Cmessage = 'lv_code not right from STASHmaster'
                Call Ereport(  RoutineName, ErrorStatus, Cmessage )

              End If

            End If   ! Hybrid levels
          End If   ! Height Gen methods
        End If    ! N_LEVELS > 1

        !--------------------------------------------------------
        ! Set pseudo level number if variable is on pseudo levels
        !--------------------------------------------------------
        If (STM_Record % pt_code > 0) Then
          Lookup(lbplev,K) = k - k_out + 1
        Endif

      End Do ! K=K_OUT...
      K_OUT=K_OUT+(N_LEVELS*N_PLEVELS)
      recondat_node => recondat_node % next
    End Do ! While associated
  End Do ! All Sections
End Do ! All internal/sub models
!-------------------------------------------------------------------
! Initialise LOOKUP fields from PPXREF
!-------------------------------------------------------------------

Do K=1, Hdr_Out % Len2Lookup
  ITEM=MOD(Lookup(item_code,K),1000)
  SECTION=(Lookup(item_code,K)-ITEM)/1000
  MODEL=Lookup(model_code,K)
  STM_Record = Rcf_Exppx( model, section, item)

  If ( Hdr_Out % FIXHD( FH_HorizGrid ) < 100) Then
    Lookup( lbcode,K)= 1
    Lookup( lbhem, K)= Hdr_Out % FixHd( FH_HorizGrid )
  Else
    Lookup( lbcode,K)= 101 !100 added for non-standard polar axis
    Lookup( lbhem, K)= Hdr_Out % FixHd( FH_HorizGrid ) - 100
  End If

  ! LBCs have an lbhem value of 99
  If (STM_Record % grid_type == ppx_atm_lbc_theta .OR.          &
  STM_Record % grid_type == ppx_atm_lbc_u      .OR.         &
  STM_Record % grid_type == ppx_atm_lbc_v) Then
    Lookup( lbhem, K)= 99
  End If

  Lookup( lbext,K )  = 0 ! No extra data
  Lookup( lbrel,K )  = 2 ! Header release number currently 2
  Lookup( lbfc,K )   = STM_Record % field_code
  Lookup( lbvc,K )   = STM_Record % lbvc_code
  Lookup( lbegin,K ) = 0
  Lookup( lbnrec,K ) = 0
  Lookup( lbproj,K ) = IPROJ
  Lookup( lbtyp,K )  = STM_Record % cf_fieldcode

  If (Lookup( lblev,K )  ==  0) Then
    Lookup( lblev,K ) = STM_Record % cf_levelcode
  End If

  ! DEPENDS ON: get_um_version_id
  Lookup( lbsrce,K )=get_um_version_id(model_id)

  If (Lookup( data_type,K )  ==  0) Then
    Lookup( data_type,K ) = STM_Record % data_type
  End If

  If (Lookup( lbpack,K )  ==  0) Then
    Lookup( lbpack,K ) = STM_Record % dump_packing
    If (DUMP_PACK.eq.2 .or. DUMP_PACK.eq.3 ) Then
      ! Do not pack data ; Override packing indicator from PPXREF
      N1 = 0   !   No packing
      Lookup( lbpack,K ) = (Lookup( lbpack,K )/10)*10 + N1
    End If
  End If

  Lookup( bmdi,K ) = Transfer( rmdi, int_val )
  Lookup( bmks,K ) = Transfer( 1.0,  int_val )
End Do

!-------------------------------------------------------------------
! Change LOOKUP to allow for change in horizontal dimensions
!-------------------------------------------------------------------

Do K=1, Hdr_Out % Len2Lookup

  ITEM=MOD(Lookup( item_code,K ),1000)
  SECTION=(Lookup( item_code,K )-ITEM)/1000
  MODEL=Lookup( model_code,K )
  STM_Record = Rcf_Exppx( model, section, item)

  If (AREA_Expand(K) == 1) Then

    ! Get N2 and N3 from whole value of LBPACK
    WHOLE=Lookup( lbpack,K )
    N2=MOD(INT(WHOLE/10),10)
    N3=MOD(INT(WHOLE/100),10)
    
    ! Find grid type to calculate grid info
    PPXREF_GRID_TYPE = STM_Record % grid_type

    If (N2 == 2 .AND. N3 == 1) Then
      Lookup( lbrow,K ) = 0
      Lookup( lbnpt,K ) = 0
    ElseIf (N2 == 1 .AND. N3 == 1) Then
      Lookup( lbrow,K ) = 0
      Lookup( lbnpt,K ) = 0
    Else
      ! Only deal with C-grid
      Select Case( ppxref_grid_type )
      Case (ppx_atm_cuall)
        Lookup( lbrow,K ) = Output_Grid % Glob_u_rows
        Lookup( lbnpt,K ) = Output_Grid % Glob_u_row_length

      Case (ppx_atm_cvall)
        Lookup( lbrow,K ) = Output_Grid % Glob_v_rows
        Lookup( lbnpt,K ) = Output_Grid % Glob_v_row_length

      Case (ppx_atm_lbc_theta) ! LBC theta points
        Lookup( lbrow,K ) = rimwidtha
        Lookup( lbnpt,K ) = Output_Grid % Glob_p_row_length

      Case (ppx_atm_lbc_u) ! LBC U points
        Lookup( lbrow,K ) = rimwidtha
        Lookup( lbnpt,K ) = Output_Grid % Glob_u_row_length - 1
                            ! Note this is actual rather than
                            ! oversized and only exists for LAM

      Case (ppx_atm_lbc_v) ! LBC V points
        Lookup( lbrow,K ) = rimwidtha
        Lookup( lbnpt,K ) = Output_Grid % Glob_v_row_length

      Case (ppx_atm_river)
        Lookup( lbrow,K ) = Output_Grid % Glob_r_rows
        Lookup( lbnpt,K ) = Output_Grid % Glob_r_row_length

      Case Default
        Lookup( lbrow,K ) = Output_Grid % Glob_p_rows
        If (LOZONE_ZONAL .AND. SECTION == 0    &     
                         .AND. (  ITEM == stashcode_ozone        &   !STASH = 60 
                         .OR.     ITEM == stashcode_o3_prod_loss &   !STASH = 481
                         .OR.     ITEM == stashcode_o3_p_l_vmr   &   !STASH = 482
                         .OR.     ITEM == stashcode_o3_vmr       &   !STASH = 483
                         .OR.     ITEM == stashcode_o3_p_l_temp  &   !STASH = 484
                         .OR.     ITEM == stashcode_o3_temp      &   !STASH = 485
                         .OR.     ITEM == stashcode_o3_p_l_colo3 &   !STASH = 486
                         .OR.     ITEM == stashcode_o3_colo3)    &   !STASH = 487
                                             ) Then

          Lookup( lbnpt,K ) = 1
        Else
          Lookup( lbnpt,K ) = Output_Grid % Glob_p_row_length
        End If
      End Select
    End If

    Lookup( bplat,K ) = Transfer( Hdr_Out % RealC( RC_PoleLat  ),    &
                                  int_val)
    Lookup( bplon,K ) = Transfer( Hdr_Out % RealC( RC_PoleLong ),    &
                                  int_val)
    Lookup( bgor, K ) = Transfer( 0.                            ,    &
                                  int_val)

    Lookup( bdy,K ) = Transfer( Hdr_Out % RealC(2), int_val )
    Lookup( bzy,K ) = Transfer( Hdr_Out % RealC(3) -               &
                                 Hdr_Out % RealC(2), int_val )

    Lookup( bzx,K ) = Transfer( Hdr_Out % RealC(4) -                  &
                             Hdr_Out % RealC(1), int_val)
    If (Lookup( lbnpt,K ) ==  1) Then
      Lookup( bdx,K ) = Transfer( 360., int_val )
    Else
      Lookup( bdx,K ) = Transfer( Hdr_Out % RealC( RC_LongSpacing),   &
                                  int_val )
    End If

    If (PPXREF_GRID_TYPE == ppx_atm_cuall) Then

      Lookup( bzy,K ) = Transfer( Hdr_Out % RealC(3) -           &
                                 Hdr_Out % RealC(2), int_val )
      Lookup( bzx,K ) = Transfer( Hdr_Out % RealC(4) -           &
                                 Hdr_Out % RealC(1) * .5, int_val )
    Else If (PPXREF_GRID_TYPE == ppx_atm_cvall) Then

      Lookup( bzy,K ) = Transfer( Hdr_Out % RealC(3) -           &
                                 Hdr_Out % RealC(2) *.5, int_val )
      Lookup( bzx,K ) = Transfer( Hdr_Out % RealC(4) -           &
                                 Hdr_Out % RealC(1), int_val )
    Else If (PPXREF_GRID_TYPE == ppx_atm_river) Then
      Lookup( bdy,K ) = Transfer( - 180.0/Output_Grid % Glob_r_rows, &
                            int_val)
      riv_bzy = - Hdr_Out % RealC(3) &
                  - 0.5*(- 180.0/Output_Grid % Glob_r_rows)
      Lookup( bzy,K ) = Transfer( riv_bzy,int_val)
      Lookup( bdx,K ) =                                          &
                  Transfer( 360.0/Output_Grid % Glob_r_row_length, &
                       int_val)
      riv_bzx = Hdr_Out % RealC(4) &
                    - 0.5*(360.0/Output_Grid % Glob_r_row_length)
      Lookup( bzx,K ) = Transfer( riv_bzx, int_val )
    End If
  End If

  ! Set lookup 13 from LCAL360. prefromed in UI prior to vn 3.5
  If (LCAL360) Then
    Lookup( lbtim,K ) = 2
  Else
    Lookup( lbtim,K ) = 1
  End If

End Do

! Set BZY and BZX to RMDI in variable resolution header 
  If (Hdr_Out %RealC(1)< 0) Then
    Do K=1, Hdr_Out % Len2Lookup 
      Lookup( bzy,K ) = Transfer( rmdi, int_val )
      Lookup( bzx,K ) = Transfer( rmdi, int_val )
    End Do
  End If

Return
End Subroutine Rcf_Setup_Lookup
End Module Rcf_Setup_Lookup_Mod

