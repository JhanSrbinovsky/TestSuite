
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Parallel RCF: Transform from global to local co-ordinates:

Module Rcf_Global_To_Local_Mod

!  Subroutine Global_To_Local_Subdomain : global subdomain boundaries
!                                         to local ones
!  Subroutine global_to_local_RC: converts global row,column co-ords to
!                                 processor co-ordinates plus local
!                                 co-ordinates within the processor.
!   Function Rcf_Get_Fld_Type : Determines P, U or V type of field
!
! Description:
!   Takes a global definition of a subdomain region (in terms of
!   model gridpoints) and translates it into local numbers.
!   This effectively means local co-ordinates of the region of the
!   subdomain which intersects with this processor's area.
!
! Method:
!   Use the datastart variable in PARVARS to see if the requested
!   subdomain intersects with this processor's area, if it does
!   then use datastart to convert to local co-ordinate and do a bit
!   of logic using MAX and MIN to ensure the local co-ordinates
!   actually lie within the local area  Then make any corrections
!   necessary to account for a subdomain which crosses over the
!   0 longitude line. Finally, if L_include_halos is set to
!   .TRUE. - include any relevant halo regions.
!
! Derived from UM 4.5 code
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

SUBROUTINE RCF_GLOBAL_TO_LOCAL_SUBDOMAIN(                           &
                                   L_include_halosEW,               &
                                   L_include_halosNS,               &
                                   grid_code,procid,                &
                                   global_north_in,global_east_in,  &
                                   global_south_in,global_west_in,  &
                                   local_north,local_east,          &
                                   local_south,local_west)
Use Rcf_Parvars_Mod

IMPLICIT NONE

! Subroutine arguments:
Logical, Intent(In)   :: L_include_halosEW  ! If true, include East-West
                                            ! halos in local region
Logical, Intent(In)   :: L_include_halosNS  ! if true incl. North-South
                                            ! halos in local region
Integer, Intent(In)   :: grid_code          ! STASH grid type
Integer, Intent(In)   :: procid             ! process result wanted for
Integer, Intent(In)   :: global_north_in    ! global northern boundary
Integer, Intent(In)   :: global_east_in     ! global eastern boundary
Integer, Intent(In)   :: global_south_in    ! global southern boundary
Integer, Intent(In)   :: global_west_in     ! global western boundary

Integer, Intent(Out)  :: local_north        ! local northern boundary
Integer, Intent(Out)  :: local_south        ! local sothern boundary
Integer, Intent(Out)  :: local_east         ! local eastern boundary
Integer, Intent(Out)  :: local_west         ! local westernboundary

! Parameters and Common blocks
Integer, Parameter    :: st_no_data = -3    ! Magic number

! Local variables
! Copies of the input arguments, that can be modified for
! wrap-around calculations
Integer   :: global_north, global_east, global_south, global_west
Integer   :: fld_type           ! is field on P or U or V grid?
Integer   :: row_len_nh         ! row length when halos are removed
Integer   :: nrows_nh           ! number of rows when halos are removed
Integer   :: first_global_pt_EW ! global point number of first and last
Integer   :: last_global_pt_EW  ! local points in local area
Integer   :: first_global_pt_NS ! in the East-West and
Integer   :: last_global_pt_NS  ! North-South directions

! Logicals indicating if this processor contains part of a
! subdomain
Logical   :: NS_intersect, EW_intersect
Logical   :: wrap      ! set to .TRUE. if the subdomain passes over the
                       ! the 0 degree longitude line
Logical   :: fullfield ! if the field is NOT a subdomain

! ------------------------------------------------------------------

! Copy the global_in variables into local variables

global_north=global_north_in
global_east=global_east_in
global_south=global_south_in
global_west=global_west_in

! Find out if the data is on a mass or velocity grid

fld_type = Rcf_Get_Fld_Type(grid_code)

IF (fld_type .EQ. fld_type_unknown) THEN
  WRITE(6,*) 'GLOBAL_TO_LOCAL_SUBDOMAIN encountered ',     &
             'field with gridtype code ',grid_code
  WRITE(6,*) 'Unable to process this field.'
  local_north=st_no_data
  local_south=st_no_data
  local_east=st_no_data
  local_west=st_no_data
  GOTO 9999
ENDIF

! Set up logical indicating if this is a full field, or just
! a subdomain

fullfield= ((( global_west .EQ. 1 ) .AND.            &
             ( global_south .EQ. 1 )) .AND.          &
           (((fld_type .EQ. fld_type_p ) .AND.       &
             ( global_north .EQ. glsize(2) ) .AND.   &
             ( global_west  .EQ. glsize(1) )) .OR.   &
            ((fld_type .EQ. fld_type_u ) .AND.       &
             ( global_north .EQ. glsizeu(2) ) .AND.  &
             ( global_west  .EQ. glsizeu(1))) .OR.   &
            ((fld_type .EQ. fld_type_v ) .AND.       &
             ( global_north .EQ. glsizev(2) ) .AND.  &
             ( global_west  .EQ. glsizev(1) ))))

! If this is a fullfield (ie. not a subdomain) the local addressing
! is easy:

IF (fullfield) THEN

  IF (L_include_halosNS) THEN
    local_north = g_lasize(2,procid)
    local_south = 1
  ELSE
    local_north = g_lasize(2,procid)-Offy
    local_south = 1+Offy
  ENDIF
  IF (L_include_halosEW) THEN
    local_west=1
    local_east=g_lasize(1,procid)
  ELSE
    local_west=1+Offx
    local_east=g_lasize(1,procid)-Offx
  ENDIF

ELSE ! a subdomain requires some careful analysis:

  IF (fld_type .EQ. fld_type_p) THEN
    row_len_nh=g_blsizep(1,procid)
    nrows_nh=g_blsizep(2,procid)
  ELSE IF (fld_type .EQ. fld_type_u) THEN
    row_len_nh=g_blsizeu(1,procid)
    nrows_nh=g_blsizeu(2,procid)
  ELSE
    row_len_nh=g_blsizev(1,procid)
    nrows_nh=g_blsizev(2,procid)
  ENDIF

! Set up variables giving the global point numbers of the
! start and end of this processor's subdomain

  first_global_pt_EW=g_datastart(1,procid)
  last_global_pt_EW=first_global_pt_EW+row_len_nh-1

  first_global_pt_NS=g_datastart(2,procid)
  last_global_pt_NS=first_global_pt_NS+nrows_nh-1

! If global_east is greater than the global row length, this
! indicates a wrap around - but not in the format this code
! wants - where it expects a wrap around to be indicated by
! the east column being less than the west column.

  IF (global_east .LT. global_west) THEN
    wrap=.TRUE.
  ELSEIF (global_east .GT. glsize(1)) THEN
    wrap=.TRUE.
    global_east=global_east-glsize(1)
  ELSE
    wrap=.FALSE.
  ENDIF

  EW_intersect =                                         &
         (( .NOT. wrap) .AND.                            &
          ((global_east .GE. first_global_pt_EW) .AND.   &
           (global_west .LE. last_global_pt_EW)))        &
         .OR.                                            &
         ((wrap) .AND.                                   &
          ((global_west .LE. last_global_pt_EW) .OR.     &
           (global_east .GE. first_global_pt_EW)))

  NS_intersect =                                         &
         ((global_south .LE. last_global_pt_NS) .AND.    &
          (global_north .GE. first_global_pt_NS))

  IF (NS_intersect) THEN

    IF ((global_south .GE. first_global_pt_NS) .AND.     &
        (global_south .LE. last_global_pt_NS)) THEN
! This processor contains the NS start of the subarea
      local_south=global_south-first_global_pt_NS+Offy+1
    ELSE
! This processor is to the North of the start of the subarea
      local_south=1+Offy
    ENDIF

    IF ((global_north .GE. first_global_pt_NS) .AND.     &
        (global_north .LE. last_global_pt_NS)) THEN
! This processor contains the NS end of the subarea
      local_north=global_north-first_global_pt_NS+Offy+1
    ELSE
! This processor is to the South of the subarea
      local_north=nrows_nh+Offy
    ENDIF

  ELSE

    local_north=st_no_data
    local_south=st_no_data

  ENDIF

  IF (EW_intersect) THEN

    IF ((global_west .GE. first_global_pt_EW) .AND.     &
        (global_west .LE. last_global_pt_EW)) THEN
! This processor contains the EW start of the subarea
      local_west=global_west-first_global_pt_EW+Offx+1
    ELSE
! This processor is to the right of the start of the subarea
      local_west=1+Offx
    ENDIF

    IF ((global_east .GE. first_global_pt_EW) .AND.     &
        (global_east .LE. last_global_pt_EW)) THEN
! This processor contains the EW end of the subarea
      local_east=global_east-first_global_pt_EW+Offx+1
    ELSE
! This processor is to the left of the end of the subarea
      local_east=Offx+row_len_nh
    ENDIF

  ELSE

    local_east=st_no_data
    local_west=st_no_data

  ENDIF

ENDIF ! is this a fullfield?

 9999 CONTINUE

RETURN
End Subroutine Rcf_Global_To_Local_Subdomain

!--------------------------------------------------------------------
! This routine is included for completeness
!--------------------------------------------------------------------

! Subroutine Interface:
SUBROUTINE Rcf_Global_To_Local_RC( grid_code,                      &
                                   global_column_in , global_row,  &
                                   processor_x , processor_y,      &
                                   local_column, local_row)

Use Rcf_Parvars_Mod

IMPLICIT NONE
! Description:
! Takes a global co-ordinate, in model gridpoints, and returns
! the processor co-ordinate of the processor containing that
! point, and the local co-ordinates of the point on that processor.
!
!
! Current code owner : Paul Burton
!
! History:
!  Model    Date     Modification history from model version 4.2
!  version
!  4.2      17 /09/96 New deck created for mpp code.  P.Burton
!  4.3      13/03/97  Various bug fixes               P.Burton
!  4.4      18/06/97  Check that row number is valid      P.Burton
!           06/10/97  Set correct row length and n_rows
!                     in dowhile loop.                    P.Burton
!  5.1      14/03/00  New reconfiguration                 P.Selwood

! Subroutine arguments:

Integer, Intent(In)  :: grid_code          ! STASH grid type code
Integer, Intent(In)  :: global_column_in   ! global column number
Integer, Intent(In)  :: global_row         ! global row number
Integer, Intent(Out) :: processor_x        ! processor X (EW) co-ord.
                                           !       (0->nproc_x)
Integer, Intent(Out) :: processor_y        ! processor Y (NS) co-ord.
                                           ! (0->nproc_y)
Integer, Intent(Out) :: local_column       ! local column no. on proc.
Integer, Intent(Out) :: local_row          ! local row number on proc.

! Parameters and COMMON blocks
Integer, Parameter   :: st_no_data = -3


! Local variables

Integer   :: global_column ! modified version of global_column_in which
!                          ! takes account of domains wrapping over
!                          ! 0 degree longitude
Integer   :: fld_type      ! field stored on P grid or U grid?
Integer   :: row_len_nh,nrows_nh ! row_len and n_rows when halos removed
Integer   :: proc  ! loop counter for loop over processors

! global column and row numbers delimiting a processors area
Integer   :: start_col, end_col, start_row, end_row

! ------------------------------------------------------------------

! Find out if the data is on a mass or velocity grid

fld_type=Rcf_Get_Fld_Type(grid_code)

IF (fld_type .EQ. fld_type_unknown) THEN
  WRITE(6,*) 'GLOBAL_TO_LOCAL_RC encountered ',       &
             'field with gridtype code ',grid_code
  WRITE(6,*) 'Unable to process this field.'
  processor_x=st_no_data
  processor_y=st_no_data
  local_column=st_no_data
  local_row=st_no_data
  GOTO 9999
ENDIF

! If global_column_in is more than the global row length, perform
! a wrap around to ensure it falls within the global bounds

IF (global_column_in .GT. glsize(1)) THEN
  global_column=MOD(global_column_in+1,glsize(1))-1
ELSE
  global_column=global_column_in
ENDIF

IF ((global_column .LT. 1) .OR.                        &
    (global_row .LT. 1) .OR.                           &
    (global_row .GT. glsize(2))) THEN

  WRITE(6,*) 'GLOBAL_TO_LOCAL_RC encountered ',            &
             'impossible global row/column co-ordinates ', &
             'row: ',global_row,' column: ',global_column

  processor_x=st_no_data
  processor_y=st_no_data
  local_column=st_no_data
  local_row=st_no_data

ENDIF

! Make a first guess at the processor co-ordinates

processor_x=MIN(global_column/(glsize(1)/gridsize(1)), nproc_x-1)
processor_y=MIN(global_row/(glsize(2)/gridsize(2)), nproc_y-1)

proc=processor_x+processor_y*gridsize(1)

IF (fld_type .EQ. fld_type_p) THEN
  row_len_nh=g_blsizep(1,proc)
  nrows_nh=g_blsizep(2,proc)
ELSE IF (fld_type .EQ. fld_type_u) THEN
  row_len_nh=g_blsizeu(1,proc)
  nrows_nh=g_blsizeu(2,proc)
ELSE
  row_len_nh=g_blsizev(1,proc)
  nrows_nh=g_blsizev(2,proc)
ENDIF

start_col=g_datastart(1,proc)
end_col=start_col+row_len_nh-1
start_row=g_datastart(2,proc)
end_row=start_row+nrows_nh-1

! Now iterate around these processors until we hit the right one

DO WHILE  (((global_column .LT. start_col) .OR.     &
            (global_column .GT. end_col  ))         &
        .OR.                                        &
           ((global_row .LT. start_row) .OR.        &
            (global_row .GT. end_row)))


  IF (global_column .LT. start_col) THEN
    processor_x=processor_x-1
  ELSEIF (global_column .GT. end_col) THEN
    processor_x=processor_x+1
  ENDIF

  IF (global_row .LT. start_row) THEN
    processor_y=processor_y-1
  ELSEIF (global_row .GT. end_row) THEN
    processor_y=processor_y+1
  ENDIF

  proc=processor_x+processor_y*gridsize(1)

IF (fld_type .EQ. fld_type_p) THEN
  row_len_nh=g_blsizep(1,proc)
  nrows_nh=g_blsizep(2,proc)
ELSE IF (fld_type .EQ. fld_type_u) THEN
  row_len_nh=g_blsizeu(1,proc)
  nrows_nh=g_blsizeu(2,proc)
ELSE
  row_len_nh=g_blsizev(1,proc)
  nrows_nh=g_blsizev(2,proc)
ENDIF
  start_col=g_datastart(1,proc)
  end_col=start_col+row_len_nh-1
  start_row=g_datastart(2,proc)
  end_row=start_row+nrows_nh-1

ENDDO

! Now we have the processor co-ordinates, we can calculate the
! local co-ordinates.

local_column=Offx+global_column-start_col+1
local_row=Offy+global_row-start_row+1

 9999 CONTINUE

RETURN
END Subroutine Rcf_Global_To_Local_RC

!-------------------------------------------------------------------
!
! Rcf only version
!
!-------------------------------------------------------------------

! Function Interface
INTEGER FUNCTION RCF_GET_FLD_TYPE (grid_type_code)

Use Rcf_Parvars_Mod
IMPLICIT NONE

!
! Description:
! Takes a STASH grid type code, and returns which type of
! grid this is - mass or wind grid.
!
! Current code owner : Paul Burton
!
! History:
!  Model    Date     Modification history from model version 4.2
!  version
!  4.2      21/11/96 New deck created for mpp code.  P.Burton
!  5.5      06/01/03 River routing support. P.Selwood.
!
! Subroutine arguments:

INTEGER, Intent(In)   :: grid_type_code     ! IN : STASH grid type code

! Parameters

! Comdeck
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

IF (((grid_type_code .GE. ppx_atm_tall)    .AND.    &
     (grid_type_code .LE. ppx_atm_tmerid)) .OR.     &
     (grid_type_code .EQ. ppx_atm_ozone)   .OR.     &
     (grid_type_code .EQ. ppx_atm_compressed) .OR.  &
     (grid_type_code .EQ. ppx_ocn_tall)    .OR.     &
     (grid_type_code .EQ. ppx_ocn_cuall)  .OR.      &
     (grid_type_code .EQ. ppx_ocn_tfield)  .OR.     &
     (grid_type_code .EQ. ppx_ocn_tzonal)  .OR.     &
     (grid_type_code .EQ. ppx_ocn_tmerid)) THEN
  Rcf_Get_Fld_Type=fld_type_p
ELSEIF                                              &
   (((grid_type_code .GE. ppx_atm_uall)    .AND.    &
     (grid_type_code .LE. ppx_atm_umerid)) .OR.     &
     (grid_type_code .EQ. ppx_atm_cuall)   .OR.     &
     (grid_type_code .EQ. ppx_ocn_uall)    .OR.     &
     (grid_type_code .EQ. ppx_ocn_ufield)  .OR.     &
     (grid_type_code .EQ. ppx_ocn_uzonal)  .OR.     &
     (grid_type_code .EQ. ppx_ocn_umerid)) THEN
  Rcf_Get_Fld_Type=fld_type_u
ELSEIF                                              &
   (((grid_type_code .EQ. ppx_atm_cvall)   .OR.     &
     (grid_type_code .EQ. ppx_ocn_cvall))) THEN
  Rcf_Get_Fld_Type=fld_type_v
ELSEIF                                              &
   (grid_type_code .EQ. ppx_atm_river) THEN
  Rcf_Get_Fld_Type=fld_type_r
ELSE
  Rcf_Get_Fld_Type=fld_type_unknown
ENDIF

RETURN

END Function Rcf_Get_Fld_Type

End Module Rcf_Global_To_Local_Mod



