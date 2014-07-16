#if defined(FRAMES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program FRAMES : Top-level program to create boundary dataset
!                   from model analyses/dumps.
!
PROGRAM FRAMES

  IMPLICIT NONE
!
! Description : Generates FRAMES data for MakeBC.
!
! Method : 
!    The top level routine FRAMES set up the program to
!    generate the Frames data in the LBC_Frames routine.
!
! Code Owner : Thomas Green
!
! Code Description :
! Language : FORTRAN 90
! Code standards : as in UMDP3 
 
  INTEGER :: internal_model   !  Internal Model Identifier

  CHARACTER (LEN=*), PARAMETER :: RoutineName='FRAMES'

! Yet to test if any of the following can be removed.

#include "parparm.h"
#include "typsize.h"
#include "csubmodl.h"
#include "cntl_io.h"
#include "typstsz.h"
#include "cprintst.h"
#include "chsunits.h"
#include "cntlall.h"

! DEPENDS ON: Timer
  CALL Timer ( RoutineName, 1 )

! DEPENDS ON: InitPrintStatus
  CALL InitPrintStatus

  WRITE (6,*) ' ######################'
  WRITE (6,*) ' Running Frames Utility'
  WRITE (6,*) ' ######################'
  WRITE (6,*) ' '

!     Only Atmosphere Model catered for
  n_internal_model = 1
  internal_model = 1
  internal_model_index(internal_model) = 1

! --------------------------------------------------
! Get list of LBC variables Frames data required for
! --------------------------------------------------

  CALL Get_LBC_List

! ---------------------
! Get required LBC area
! ---------------------

  CALL Get_LBC_Areas
 
! -----------------------------------------------
! Determine LBC grid dimensions for selected area
! -----------------------------------------------

  CALL Get_LBC_Grid

! ------------------------
! Generate the Frames data
! ------------------------

  CALL LBC_Frames

! DEPENDS ON: Timer
  CALL Timer ( RoutineName, 2 )

  WRITE (6,*) ' '
  WRITE (6,*) ' ##################################'
  WRITE (6,*) ' FRAMES program completed normally.'
  WRITE (6,*) ' ##################################'

  STOP   !  required ?

CONTAINS
!     =============================================

  SUBROUTINE LBC_Frames

!
! Description :
!       Main subroutine of Frames program to generate the data for MakeBC
!
! Method :
!       Order of processing in routine :
!       1. Determine total number of lookup headers in input fieldsfiles.
!       2. Read in all lookup headers.
!       3. Determine grid sizes for LBC area selected.
!       4. Call LBC_Interp_Coeffs to work out data points required on
!          input grid.
!       5. Set up a mask (1/0) corresponding to data points required on
!          input grid (1=required) for MakeBC. Resulting mask will be in
!          the shape of a frame.
!       6. Determine sub-area containing the frame.
!       7. Get forecast time (ie. T+3) of input data required.
!       8. Open output file for farmes data and set up headers.
!       9. Determine which input fieldsfiles to read data from.
!      10. Read in data.
!      11. Extract sub-area data from full-area data.
!      12. Apply mask to set all data in sub_area that is not required
!          to RMDI.
!      13. Set up ookup header for output file.
!      14. Write out sub-area data.
!      15. Write out headers to output file.
!
! Code Owner : Thomas Green
!
!   ----------------------------------------------------------------------
!   ----------------------------------------------------------------------
!   This subroutine was never meant to be as long as it stands. There is 
!   plenty of scope to break it down into smaller routines. It was decided
!   to leave this until after the long term requirements of this program
!   have been reviewed.
!   ----------------------------------------------------------------------
!   ----------------------------------------------------------------------

    USE LBC_Frames_mod

    IMPLICIT NONE

!   -------------------------
!   Input grid specifications
!   -------------------------

    INTEGER :: src_row_len
    INTEGER :: src_rows
    INTEGER :: src_model_levels
    INTEGER :: src_wet_levels

    REAL    :: src_delta_lat
    REAL    :: src_delta_long
    REAL    :: src_first_lat
    REAL    :: src_first_long
    REAL    :: src_pole_lat
    REAL    :: src_pole_long

    LOGICAL :: src_cyclic
    LOGICAL :: src_rotated

!   -------------------------
!   LBC grid specifications
!   -------------------------

    INTEGER :: lbc_row_len
    INTEGER :: lbc_rows
    REAL    :: lbc_delta_lat
    REAL    :: lbc_delta_long
    REAL    :: lbc_first_lat
    REAL    :: lbc_first_long
    REAL    :: lbc_pole_lat
    REAL    :: lbc_pole_long
    INTEGER :: rimwidth
    INTEGER :: lbc_halo_x
    INTEGER :: lbc_halo_y

    INTEGER :: lbc_size
    INTEGER :: max_lbc_size = -1

    INTEGER, ALLOCATABLE :: lbc_index_bl   (:)
    INTEGER, ALLOCATABLE :: lbc_index_br   (:)
    REAL,    ALLOCATABLE :: lbc_weights_bl (:)
    REAL,    ALLOCATABLE :: lbc_weights_br (:)
    REAL,    ALLOCATABLE :: lbc_weights_tl (:)
    REAL,    ALLOCATABLE :: lbc_weights_tr (:)
    REAL,    ALLOCATABLE :: lbc_coeff1     (:)
    REAL,    ALLOCATABLE :: lbc_coeff2     (:)
    REAL,    ALLOCATABLE :: lambda_in      (:)
    REAL,    ALLOCATABLE :: phi_in         (:)

    INTEGER, ALLOCATABLE :: mask (:,:)
    INTEGER, ALLOCATABLE :: col_mask (:)

    INTEGER, ALLOCATABLE :: indata (:,:)
    INTEGER, ALLOCATABLE :: sub_area_data (:,:)


    INTEGER :: i,j,ij,ijmin

    INTEGER :: jn, js, iw, ie   !  Boundaries for sub-area
    INTEGER :: jnv
    INTEGER :: i_uv
    INTEGER :: fld_type
    INTEGER :: lbc_fld_type
    INTEGER :: lbc_halo_type
    REAL    :: percentage
    INTEGER :: npts_mask
    INTEGER :: npts_full_area
    INTEGER :: npts_sub_area
    INTEGER :: npts_sub_area_v
    INTEGER :: sub_area_row_len
    INTEGER :: sub_area_rows
    INTEGER :: sub_area_rows_v

    INTEGER :: ifld
    LOGICAL :: L_GM
    LOGICAL :: print_sub_area_data
    LOGICAL :: l_apply_mask
    LOGICAL :: l_u_field
    LOGICAL :: l_v_field
    REAL :: mask_value

    INTEGER :: n_lookups (0:20)
    INTEGER :: LookUp_Start_Address(20)
    INTEGER :: i_FF, n_FF
    INTEGER :: ihead

    INTEGER :: len_env    = 6  !  Length of env. variable
    INTEGER :: env_var_0  = 0  !  Indicator that filename is in env var
    INTEGER :: env_var_1  = 1  !  Indicator that filename is in variable
    INTEGER :: read_only  = 0  !  Input Fieldsfiles  - read only
    INTEGER :: read_write = 1  !  Frames Fieldsfiles - read & write
    INTEGER :: keep_file  = 0  !  Keep Frames FFs at end of job
    CHARACTER*6 env            !  Env Variable for input fieldsfiles

    INTEGER :: FF_unit_no         ! Unit No for FFs (31-)
    INTEGER :: frames_unit_no = 20

    INTEGER icode            !  Error code
    CHARACTER*80 cmessage    !  Error Message
    CHARACTER (len=100) :: FileName
    CHARACTER (Len=*), PARAMETER :: RoutineName = 'LBC_Frames'

    INTEGER :: ReturnCode
    REAL    :: a_io
    INTEGER :: len_io

    INTEGER :: row_length
    INTEGER :: rows
    INTEGER :: nlevs

!   ----------------------------------
!   Header sizes for Input Fieldsfiles
!   ----------------------------------

    INTEGER :: Len_FixHd = 256
!     Integer :: FixHd(Len_FixHd)  ! Compiler doesn't like this line ?
    INTEGER :: FixHd (256)

    INTEGER :: Len_IHead = 46
    INTEGER :: Len_RHead = 38
    INTEGER :: IntHead (46)
    REAL    :: RealHead(38)

    INTEGER :: Len1_LevDepC
    INTEGER :: Len2_LevDepC
    INTEGER :: Len_LevDepC
    REAL, ALLOCATABLE :: LevDepC (:,:)

    INTEGER :: Len1_Lookup
    INTEGER :: Len2_Lookup
    INTEGER :: Len_Lookup
    INTEGER, ALLOCATABLE :: Lookup (:,:)
    INTEGER, ALLOCATABLE :: Lookup_FF_No (:)

!   -----------------------------------
!   Header sizes for Output Frames File
!   -----------------------------------

    INTEGER :: Frames_FixHd (256)
    INTEGER :: Frames_IntHd (46)
    REAL    :: Frames_RealHd(38)
    REAL, ALLOCATABLE :: Frames_LevDepC (:,:)
    INTEGER :: Frames_LookUp(64,4096)

!   ---------------------------------
!   Data arrays for Input/Output data
!   ---------------------------------

    REAL, ALLOCATABLE :: Data_FF       (:,:,:)
    REAL, ALLOCATABLE :: Data_Frames   (:,:,:)


    INTEGER :: ipt(Num_LBC_Vars)
    INTEGER :: Stash_Code
    INTEGER :: jvar
    INTEGER :: k
    INTEGER :: expand = 1

    INTEGER :: int_val = 1      ! used for transfers
    REAL :: frames_first_lat
    REAL :: frames_first_long
    REAL :: zeroth_lat
    REAL :: zeroth_long

    INTEGER :: k0 = 0
    INTEGER :: kk
    INTEGER :: Start_Address
    INTEGER :: FC_Time
    CHARACTER (Len=5) :: C_FC_Time

!   --------------------------------------------------------
!   FRAMES has not been tested with variable resolution yet.
!   Hardwire so VR switched off for now.
!   --------------------------------------------------------
    LOGICAL :: l_var_lbc = .FALSE.


#include "cmaxsize.h"
#include "cintfa.h"
#include "c_mdi.h"

#include "cprintst.h"

#include "parvars.h"
#include "parlbcs.h"

!   ----------------------------------
!   Work out how many FFs to be opened
!   ----------------------------------

    n_FF = 0
    DO i_FF = 1, 2

      ff_unit_no = 30 + i_FF

      env = 'FILE  '
      WRITE (env(5:6),'(I2)') ff_unit_no

      FileName = 'NotSet'

      CALL Fort_Get_Env (env, 6, FileName, 100, ICode)

!      write (6,*) ' ICode ',ICode
!      write (6,*) ' Filename ',Filename

      IF (Icode == -1) THEN

! Return Code is -1 if Env Var is not found
        WRITE (6,*) ' Env Var ',env,' not found.'
        CYCLE

      END IF

      IF (Icode == 0) THEN

! Env Var found
        WRITE (6,*) ' Env Var ',env,' found.'
        WRITE (6,*) ' Filename ',Filename
        n_FF = n_FF + 1

      END IF

    END DO
    WRITE (6,*) ' n_FF = ',n_FF

    IF (n_FF == 0) THEN
      ICode = 13
      CMessage = 'No input FieldsFiles available ??'
! DEPENDS ON: Ereport
      CALL Ereport ( RoutineName, icode, CMessage )
    END IF

!   --------------------
!   Open the fieldsfiles
!   --------------------

    n_lookups (:) = 0

    DO i_FF = 1, n_FF

      ff_unit_no = 30 + i_FF

      env = 'FILE  '
      WRITE (env(5:6),'(I2)') ff_unit_no
      WRITE (6,*) ' '

!     Test env var exists

! DEPENDS ON: File_Open
      CALL File_Open (ff_unit_no,                                &
         env,len_env,read_only,env_var_0,icode)
      IF (icode /= 0) THEN
        CMessage = 'Error in opening fieldfiles'
! DEPENDS ON: Ereport
        CALL Ereport ( RoutineName, icode, CMessage )
      END IF

!     ------------------------
!     Read in the fixed header
!     ------------------------

! DEPENDS ON: Read_FLH
      CALL Read_FLH ( ff_unit_no,fixhd,len_fixhd,icode,cmessage )
      IF (icode > 0) THEN
        CMessage = 'Error in Read_FLH'
! DEPENDS ON: Ereport
        CALL Ereport ( RoutineName, icode, CMessage )
      END IF

      IF (fixhd(5) == 3) THEN
        WRITE (6,*) 'This is a FieldsFile'
      ELSE
        CMessage= 'Invalid input file type'
        Returncode = 10
! DEPENDS ON: Ereport
        CALL Ereport( RoutineName, ReturnCode, CMessage )
      END IF

!     -------------------------------
!     Move to start of Integer Header
!     -------------------------------

! DEPENDS ON: SetPos
      CALL SetPos ( FF_unit_no, FixHd(100)-1, ICode )

      IF (ICode /= 0) THEN
        CMessage = 'Error in SetPos for Integer Header'
! DEPENDS ON: Ereport
        CALL Ereport ( RoutineName, icode, CMessage )
      END IF

!     ----------------------
!     Read in Integer Header
!     ----------------------

! DEPENDS ON: BuffIn
      CALL BuffIn ( FF_unit_no,                                  &
                    IntHead, Len_IHead, Len_io, a_io )

      IF (A_IO /= -1.0 .OR. LEN_IO /= len_ihead) THEN
! DEPENDS ON: IOERROR
        CALL IOERROR ('buffer in of Integer header',             &
                       A_IO, LEN_IO, len_ihead)
        CMESSAGE = 'I/O ERROR with buffin of Integer Header'
        ICode    = 11
! DEPENDS ON: Ereport
        CALL EReport (RoutineName, ICode, CMessage)
      END IF

!     ---------------------------------------
!     Get grid dimensions from Integer Header
!     ---------------------------------------

      src_row_len = IntHead (6)
      src_rows    = IntHead (7)
      src_model_levels  = IntHead (8)
      src_wet_levels    = IntHead (9)

      WRITE (6,*) ' Model grid : ',src_row_len,' x ',src_rows
      WRITE (6,*) ' No of model levels in FFs ',src_model_levels
      WRITE (6,*) ' No of wet   levels in FFs ',src_wet_levels

!     ----------------------------
!     Move to start of Real Header
!     ----------------------------

! DEPENDS ON: SetPos
      CALL SetPos ( FF_unit_no, FixHd(105)-1, ICode )

      IF (ICode /= 0) THEN
        CMessage = 'Error in SetPos for Real Header'
! DEPENDS ON: Ereport
        CALL Ereport ( RoutineName, icode, CMessage )
      END IF

!     -------------------
!     Read in Real Header
!     -------------------

! DEPENDS ON: BuffIn
      CALL BuffIn ( FF_unit_no,                                &
                    RealHead, Len_RHead, Len_io, a_io )

      IF (A_IO /= -1.0 .OR. LEN_IO /= len_rhead) THEN
! DEPENDS ON: IOERROR
        CALL IOERROR ('buffer in of Real header',              &
                       A_IO, LEN_IO, len_rhead)
        CMESSAGE = 'I/O ERROR with buffin of Real Header'
        ICode    = 11
! DEPENDS ON: Ereport
        CALL EReport (RoutineName, ICode, CMessage)
      END IF

!     ----------------------------------------
!     Get grid information from Integer Header
!     ----------------------------------------

      src_delta_long = RealHead (1)
      src_delta_lat  = RealHead (2)
      src_first_lat  = RealHead (3)
      src_first_long = RealHead (4)
      src_pole_lat   = RealHead (5)
      src_pole_long  = RealHead (6)

      WRITE (6,*) ' src_delta_lat  ',src_delta_lat
      WRITE (6,*) ' src_delta_long ',src_delta_long
      WRITE (6,*) ' src_first_lat  ',src_first_lat
      WRITE (6,*) ' src_first_long ',src_first_long
      WRITE (6,*) ' src_pole_lat   ',src_pole_lat
      WRITE (6,*) ' src_pole_long  ',src_pole_long

!     ------------------------------------------
!     Move to start of Level Dependent Constants
!     ------------------------------------------

! DEPENDS ON: SetPos
      CALL SetPos ( FF_unit_no, FixHd(110)-1, ICode )

      IF (ICode /= 0) THEN
        CMessage = 'Error in SetPos for LevDepC'
! DEPENDS ON: Ereport
        CALL Ereport ( RoutineName, icode, CMessage )
      END IF

!     ---------------------------------
!     Read in Level Dependent Constants
!     ---------------------------------

      Len1_LevDepC = FixHd(111)
      Len2_LevDepC = FixHd(112)
      Len_LevDepC  = Len1_LevDepC * Len2_LevDepC
!     write (6,*) ' Len1_LevDepC ',Len1_LevDepC
!     write (6,*) ' Len2_LevDepC ',Len2_LevDepC

      ALLOCATE ( LevDepC (Len1_LevDepC, Len2_LevDepC) )

! DEPENDS ON: BuffIn
      CALL BuffIn ( FF_unit_no,                               &
                    LevDepC, Len_LevDepC, Len_io, a_io )

      IF (A_IO /= -1.0 .OR. LEN_IO /= Len_LevDepC) THEN
! DEPENDS ON: IOERROR
        CALL IOERROR ('buffer in of LevDepC header',          &
                       A_IO, LEN_IO, len_LevDepC)
        CMESSAGE = 'I/O ERROR with buffin of LevDepC '
        ICode    = 12
! DEPENDS ON: Ereport
        CALL EReport (RoutineName, ICode, CMessage)
      END IF

!     ---------------------------------------
!     Get lookup dimensions from fixed header
!     ---------------------------------------

      LookUp_Start_Address(i_FF) = FixHd(150)
      Len1_Lookup = FixHd(151)
      Len2_Lookup = FixHd(152)
      Len_Lookup  = Len1_Lookup * Len2_Lookup
      WRITE (6,*) ' len1_lookup ',len1_lookup,  &
         ' len2_lookup ',len2_lookup,  &
         ' len_lookup  ',len_lookup

!     -------------------------------
!     Allocate space for LookUp table
!     -------------------------------

      ALLOCATE ( LookUp (Len1_Lookup, Len2_Lookup) )

!     -----------------------------
!     Move to start of LookUp table
!     -----------------------------

! DEPENDS ON: SetPos
      CALL SetPos ( FF_unit_no, FixHd(150)-1, ICode )

      IF (ICode /= 0) THEN
        CMessage = 'Error in SetPos for LookUp table'
! DEPENDS ON: Ereport
        CALL Ereport ( RoutineName, icode, CMessage )
      END IF

!     --------------------
!     Read in LookUP table
!     --------------------

! DEPENDS ON: BuffIn
      CALL BuffIn ( FF_unit_no,                                 &
                    LookUp, Len_Lookup, Len_io, a_io )

      IF (A_IO /= -1.0 .OR. LEN_IO /= len_lookup) THEN
! DEPENDS ON: IOERROR
        CALL IOERROR ('buffer in of lookup header',             &
                       A_IO, LEN_IO, len_lookup)
        CMESSAGE = 'I/O ERROR with buffin of Lookup Header'
        ICode    = 11
! DEPENDS ON: Ereport
        CALL EReport (RoutineName, ICode, CMessage)
      END IF

!     -----------------------------
!     Find number of LookUp headers
!     -----------------------------

      DO i = 1, FixHd(152)
        IF ( Lookup (42,i) == -99 ) THEN
          EXIT
        ELSE
          n_lookups(i_FF) = n_lookups(i_FF) +1
        END IF
      END DO

!     Update total no of LookUp Entries
      n_lookups(0) = n_lookups(0) + n_lookups(i_FF)

      WRITE (6,*) ' n_lookups ',n_lookups(i_FF)

!     Deallocate LookUp Table for File i_FF
      DEALLOCATE (LookUp)

    END DO   !  Loop over i_FF

    WRITE (6,*) ' total n_lookups ',n_lookups(0)

!   ----------------------------------------
!   Total number of LookUp entries now known
!   Allocate space for all LookUp entries
!   ----------------------------------------

    ALLOCATE ( LookUp (64, n_lookups(0) ) )
    ALLOCATE ( LookUp_FF_No (n_lookups(0)) )

    ihead = 0

    DO i_FF = 1, n_FF

      ff_unit_no = 30 + i_FF

      env = 'FILE  '
      WRITE (env(5:6),'(I2)') ff_unit_no
      WRITE (6,*) ' '

!     -----------------------------
!     Move to start of LookUp table
!     -----------------------------

! DEPENDS ON: SetPos
      CALL SetPos ( FF_unit_no, LookUp_Start_Address(i_FF)-1, ICode )

      IF (ICode /= 0) THEN
        CMessage = 'Error in SetPos for LookUp table'
! DEPENDS ON: Ereport
        CALL Ereport ( RoutineName, icode, CMessage )
      END IF

!     --------------------
!     Read in LookUP table
!     --------------------

      Len_LookUp = 64 * n_lookups(i_FF)

! DEPENDS ON: BuffIn
      CALL BuffIn ( FF_unit_no, LookUp (1,ihead+1),       &
                    Len_Lookup, Len_io, a_io )

      IF (A_IO /= -1.0 .OR. LEN_IO /= len_lookup) THEN
! DEPENDS ON: IOERROR
        CALL IOERROR ('buffer in of lookup header',       &
                       A_IO, LEN_IO, len_lookup)
        CMESSAGE = 'I/O ERROR with buffin of Lookup Header'
        ICode    = 11
! DEPENDS ON: Ereport
        CALL EReport (RoutineName, ICode, CMessage)
      END IF

      DO i = 1, n_lookups (i_FF)
        LookUp_FF_No (ihead+i) = i_FF
      END DO

      ihead = ihead + n_lookups (i_FF)

    END DO

!   ------------------------------
!   Initialise Src_Grid in parlbcs
!   ------------------------------

!   Do Ifld = 1, NFld_Max
    DO Ifld = 1, 3

      IF     (Ifld == fld_type_p) THEN   ! p-grid

!         Src_Grid (Ifld,1) = realhd(rh_baselat)
!         Src_Grid (Ifld,2) = realhd(rh_baselong)

        Src_Grid (Ifld,1) = RealHead(3)
        Src_Grid (Ifld,2) = RealHead(4)

      ELSEIF (Ifld == fld_type_u) THEN   ! u-grid

!         Src_Grid (Ifld,1) = realhd(rh_baselat)
!         Src_Grid (Ifld,2) = realhd(rh_baselong) +
!    &                        0.5 * realhd(rh_deltaEW)

        Src_Grid (Ifld,1) = RealHead(3)
        Src_Grid (Ifld,2) = RealHead(4) + 0.5 * RealHead(1)

      ELSEIF (Ifld == fld_type_v) THEN   ! v-grid

!         Src_Grid (Ifld,1) = realhd(rh_baselat) +
!    &                        0.5 * realhd(rh_deltaNS)
!         Src_Grid (Ifld,2) = realhd(rh_baselong)

        Src_Grid (Ifld,1) = RealHead(3) + 0.5 * RealHead(2)
        Src_Grid (Ifld,2) = RealHead(4)

      END IF

      WRITE (6,*) ' ifld ',ifld,                     &
         ' src_grid_1 ',Src_Grid (Ifld,1),  &
         ' src_grid_2 ',Src_Grid (Ifld,2)

    END DO

    Src_Cyclic = (FixHd(4) < 3)  !  If Global or NH or SH
    Src_Rotated= (RealHead(5)/= 90.0 .OR. RealHead(6) /= 0.0)

    WRITE (6,*) ' Src_Cyclic  : ',Src_Cyclic
    WRITE (6,*) ' Src_Rotated : ',Src_Rotated

!   ----------------------------------------------------------
!   Get the STASHmaster information required for LBC variables
!   ----------------------------------------------------------

    CALL Get_StashMaster_Information ( src_model_levels,  &
                                       src_wet_levels )

!   ----------------------------
!   Allocate/Initialise the mask
!   ----------------------------

    ALLOCATE ( mask (src_row_len, src_rows) )
    mask(:,:) = 0

    DO ifld = 1, Num_LBC_Vars

      lbc_fld_type  = lbc_variables (ifld) % fld_type
      lbc_halo_type = lbc_variables (ifld) % halo_type

      WRITE (6,*) ' ################################ '
      WRITE (6,*) ' ifld ',ifld,' LBC Stash Code ',   &
         LBC_Variables(ifld) % Stashcode

      WRITE (6,*) ' fld_type ',lbc_fld_type,          &
         ' halo_type ',lbc_halo_type

      WRITE (6,*) ' lbc_fld_type  ',lbc_fld_type
      WRITE (6,*) ' lbc_halo_type ',lbc_halo_type

      lbc_size = lbc_Interp_lenrima ( lbc_fld_type, lbc_halo_type)
      WRITE (6,*) ' lbc_size from LBC_Interp_Lenrima ',lbc_size
!     write (6,*) ' test 0,0 ',lbc_Interp_lenrima (0,0)

      IF ( (lbc_fld_type == fld_type_u) .OR.          &
           (lbc_fld_type == fld_type_v) ) THEN
        i_uv = 1
      ELSE
        i_uv = 0
      END IF

      WRITE (6,*) ' lbc_size  ',lbc_size
      WRITE (6,*) ' i_uv      ',i_uv

      ALLOCATE ( lbc_index_bl   (lbc_size) )
      ALLOCATE ( lbc_index_br   (lbc_size) )
      ALLOCATE ( lbc_weights_bl (lbc_size) )
      ALLOCATE ( lbc_weights_br (lbc_size) )
      ALLOCATE ( lbc_weights_tl (lbc_size) )
      ALLOCATE ( lbc_weights_tr (lbc_size) )
      ALLOCATE ( lbc_coeff1     (lbc_size) )
      ALLOCATE ( lbc_coeff2     (lbc_size) )

      lbc_row_len   = intf_row_length(1)
      lbc_rows      = intf_p_rows(1)
      lbc_delta_lat = intf_nsspace(1)
      lbc_delta_long= intf_ewspace(1)
      lbc_first_lat = intf_firstlat(1)
      lbc_first_long= intf_firstlong(1)
      lbc_pole_lat  = intf_polelat(1)
      lbc_pole_long = intf_polelong(1)
      rimwidth      = intfwidtha(1)
      lbc_halo_x    = intf_exthalo_ew(1)
      lbc_halo_y    = intf_exthalo_ns(1)

!     write (6,*) ' lbc_row_len    ',lbc_row_len
!     write (6,*) ' lbc_rows       ',lbc_rows
!     write (6,*) ' lbc_delta_lat  ',lbc_delta_lat
!     write (6,*) ' lbc_delta_long ',lbc_delta_long
!     write (6,*) ' lbc_first_lat  ',lbc_first_lat
!     write (6,*) ' lbc_first_long ',lbc_first_long
!     write (6,*) ' lbc_pole_lat   ',lbc_pole_lat
!     write (6,*) ' lbc_pole_long  ',lbc_pole_long
!     write (6,*) ' rimwidth       ',rimwidth
!     write (6,*) ' lbc_halo_x     ',lbc_halo_x
!     write (6,*) ' lbc_halo_y     ',lbc_halo_y

      ALLOCATE ( lambda_in (1 - lbc_halo_x : lbc_row_len + lbc_halo_x) )
      ALLOCATE ( phi_in    (1 - lbc_halo_y : lbc_rows    + lbc_halo_y) )

!     -------------------------------------------------------------------------
!     Call UM routine LBC_Interp_Coeffs. The required output arrays are 
!     LBC_INDEX_BL and LBC_INDEX_BR. These are the data points on the input grid
!     required for interpolation to the LBC points. NB. 4 points are required,
!     the other two points are TL=(BL+row_length) and TR=(BR+row_length)
!     -------------------------------------------------------------------------

! DEPENDS ON: lbc_interp_coeffs
      CALL lbc_interp_coeffs (            &
             lbc_size,                    &
             src_row_len,                 &
             src_rows,                    &
             src_delta_lat,               &
             src_delta_long,              &
             src_grid (lbc_fld_type, 1),  & ! first lat
             src_grid (lbc_fld_type, 2),  & ! first long
             src_pole_lat,                &              
             src_pole_long,               &
             src_cyclic,                  & 
             src_rotated,                 &
             lbc_row_len,                 &
             lbc_rows,                    &
             l_var_lbc,                   & 
             lambda_in,                   &
             phi_in,                      &
             lbc_delta_lat,               &
             lbc_delta_long,              &
             lbc_first_lat,               &
             lbc_first_long,              &
             lbc_pole_lat,                &
             lbc_pole_long,               &
             rimwidth,                    &
             lbc_halo_x,                  & 
             lbc_halo_y,                  &
             lbc_index_bl,                &
             lbc_index_br,                &
             lbc_weights_tr,              &
             lbc_weights_br,              &
             lbc_weights_bl,              &
             lbc_weights_tl,              &
             lbc_coeff1,                  &
             lbc_coeff2,                  &
             i_uv                         &
         )


!      if (ifld<4) then 

!        do ij = 1,lbc_size

!        write (6,*) ' ifld ',ifld
!        do ij = 24310,24450
!        write (6,*) ' ij ',ij,' ind_bl ',lbc_index_bl(ij),
!    &                         ' ind_br ',lbc_index_br(ij)
!        write (6,*) ' ij ',ij,' wts_bl ',lbc_weights_bl(ij),
!    &                         ' wts_br ',lbc_weights_br(ij)
!        write (6,*) ' ij ',ij,' wts_tl ',lbc_weights_tl(ij),
!    &                         ' wts_tr ',lbc_weights_tr(ij)
!        enddo

!      endif 

      max_lbc_size = MAX (lbc_size, max_lbc_size)


 
      WRITE (6,*) ' ################'
      WRITE (6,*) ' Setting the mask'
      WRITE (6,*) ' ################'

!     ---------------
!     Set up the mask
!     ---------------

      DO ij = 1, lbc_size
!      Do ij = 1, 10

!        write (6,*) ' ij ',ij,' lbc_index_bl ',lbc_index_bl(ij)

!       ----------------------------------------
!       For BL point, determine row/col for mask
!       ----------------------------------------

        j = INT ( (lbc_index_bl(ij)-1)/src_row_len ) + 1
        i = MOD ( (lbc_index_bl(ij)-1), src_row_len ) + 1

!        write (6,*) ' ij ',ij,' i ',i,' j ',j

!        if (ifld==2) then
!        write (6,*) ' ij ',ij,' lbc_index_bl ',lbc_index_bl(ij)
!        write (6,*) ' ij ',ij,' i ',i,' j ',j
!        endif

!       ---------------------------------------
!       (I,J)   is the point corrsponding to BL
!       (I,J+1) is the point corrsponding to TL
!       ---------------------------------------

        mask (i,j)   = 1
        mask (i,j+1) = 1

!       --------------------------------------------------
!       Extend mask to include extra point below (j-1), on 
!       LH side (i-1) and above (j+2). See notes below.
!       --------------------------------------------------

        mask (i  ,j-1) = 1
        mask (i-1,j-1) = 1
        mask (i-1,j  ) = 1
        mask (i-1,j+1) = 1
        mask (i-1,j+2) = 1
        mask (i  ,j+2) = 1

!        write (6,*) ' ij ',ij,' lbc_index_br ',lbc_index_br(ij)

!       ----------------------------------------
!       For BR point, determine row/col for mask
!       ----------------------------------------

        j = INT ( (lbc_index_br(ij)-1)/src_row_len ) + 1
        i = MOD ( (lbc_index_br(ij)-1), src_row_len ) + 1

!        write (6,*) ' ij ',ij,' i ',i,' j ',j

!        if (ifld==2) then
!        write (6,*) ' ij ',ij,' lbc_index_br ',lbc_index_br(ij)
!        write (6,*) ' ij ',ij,' i ',i,' j ',j
!        endif

!       ---------------------------------------
!       (I,J)   is the point corrsponding to BR
!       (I,J+1) is the point corrsponding to TR
!       ---------------------------------------

        mask (i,j)   = 1
        mask (i,j+1) = 1

!       --------------------------------------------------
!       Extend mask to include extra point below (j-1), on
!       RH side (i+1) and above (j+2). See notes below.
!       --------------------------------------------------

        mask (i  ,j-1) = 1
        mask (i+1,j-1) = 1
        mask (i+1,j  ) = 1
        mask (i+1,j+1) = 1
        mask (i+1,j+2) = 1
        mask (i  ,j+2) = 1

!       -------------------------------------------------------------------
!       NOTES on Mask Extension

!       The mask extension may not be necessary. When testing MakeBC using 
!       input Frames data, it was (occassionally) picking up RMDI data one
!       point outside the FRAME so the rim has been extended by one point
!       both on the inside and outside the rim as a safety measure until it
!       has been fully investigated.
!       -------------------------------------------------------------------

      END DO

      DEALLOCATE ( lbc_index_bl   )
      DEALLOCATE ( lbc_index_br   )
      DEALLOCATE ( lbc_weights_bl )
      DEALLOCATE ( lbc_weights_br )
      DEALLOCATE ( lbc_weights_tl )
      DEALLOCATE ( lbc_weights_tr )
      DEALLOCATE ( lbc_coeff1     )
      DEALLOCATE ( lbc_coeff2     )
      DEALLOCATE ( lambda_in      )
      DEALLOCATE ( phi_in         )


    END DO   !  ifld

!   ----------------------------------------------------------------
!   Set up column mask ; this records whether any of the mask points
!   have been set in each column.
!   ----------------------------------------------------------------

    ALLOCATE ( col_mask (src_row_len) )
    col_mask(:) = 0

    DO i=1,src_row_len
      DO j=1,src_rows
        IF (mask(i,j)==1) THEN
          col_mask(i) = 1
          EXIT
        END IF
      END DO
    END DO

!   -------------------------------
!   Count no of data points in MASK
!   -------------------------------

    npts_mask = 0
    DO j=1,src_rows
      DO i=1,src_row_len
        IF (mask(i,j) == 1 ) THEN
          npts_mask = npts_mask + 1
        END IF
      END DO
    END DO

!   --------------
!   Print out MASK
!   --------------

    IF (PrintStatus == PrStatus_Diag) THEN

      DO i=1,src_row_len,100
        ijmin=MIN(i+99,src_row_len)
        WRITE (6,*) ' ijmin ',ijmin
        DO j=1,src_rows
          WRITE (6,10) (mask(ij,j),ij=i,ijmin)
10        FORMAT (' ',100I1)
        END DO
        WRITE (6,*) ' '
      END DO

    END IF

!   ------------------
!   Print out col mask
!   ------------------

    IF (PrintStatus == PrStatus_Diag) THEN

      DO i=1,src_row_len,50
        ijmin=MIN(i+49,src_row_len)
        WRITE (6,*) ' ijmin ',ijmin
        WRITE (6,10) (col_mask(ij),ij=i,ijmin)
        WRITE (6,*) ' '
      END DO

    END IF

!   ------------------------------------------------
!   Check if LBC grid crosses the Greenwich Meridian
!   ------------------------------------------------

    L_GM = ( col_mask(1) == 1 .AND. col_mask(src_row_len) == 1 ) 
    IF (L_GM) THEN
      WRITE (6,*) ' ======================================= '
      WRITE (6,*) ' LBC grid crosses the Greenwich Meridian '
      WRITE (6,*) ' ======================================= '
    END IF

!   --------------------------------------------
!   Determine JS : Southern boundary of Sub-Area
!   --------------------------------------------

    js = -1
    South : DO j = 1, src_rows
      DO i = 1, src_row_len
        IF (mask(i,j) == 1) THEN
          JS = J
          EXIT South
        END IF
      END DO
    END DO  South

    IF (JS == -1) THEN
      CMESSAGE = ' JS = -1 ; Southern Boundary for Mask not set.'
      ICODE = 61
! DEPENDS ON: Ereport
      CALL EReport (RoutineName, ICode, CMessage)
    END IF

    WRITE (6,*) ' Southern Row for sub-area   : ',js

!   --------------------------------------------
!   Determine JN : Northern boundary of Sub-Area
!   --------------------------------------------

    jn = -1
    North : DO j = src_rows, 1, -1
      DO i = 1, src_row_len
        IF (mask(i,j) == 1) THEN
          JN = J
          EXIT North
        END IF
      END DO
    END DO  North

    IF (JN == -1) THEN
      CMESSAGE = ' JN = -1 ; Northern Boundary for Mask not set.'
      ICODE = 62
! DEPENDS ON: Ereport
      CALL EReport (RoutineName, ICode, CMessage)
    END IF

    WRITE (6,*) ' Northern Row for sub-area - v   : ',jn

    jnv = jn
    jn  = jn + 1

    WRITE (6,*) ' Northern Row for sub-area       : ',jn

!   -------------------------------------------
!   Determine IW : Western boundary of Sub-Area
!   -------------------------------------------

    IF (L_GM) THEN

      iw = -1
      DO i = src_row_len, 1, -1
        IF (col_mask(i) == 1) THEN
          iw = i
        ELSE
          EXIT
        END IF
      END DO

    ELSE

      iw = -1
      DO i = 1, src_row_len
        IF (col_mask(i) == 1) THEN
          iw = i
          EXIT 
        END IF
      END DO

    END IF

    IF (IW == -1) THEN
      CMESSAGE = ' IW = -1 ; Western Boundary for Mask not set.'
      ICODE = 63
! DEPENDS ON: Ereport
      CALL EReport (RoutineName, ICode, CMessage)
    END IF

    WRITE (6,*) ' Western Column for sub-area : ',iw

!   -------------------------------------------
!   Determine IE : Eastern boundary of Sub-Area
!   -------------------------------------------

    IF (L_GM) THEN

      ie = -1
      DO i = 1, src_row_len
        IF (col_mask(i) == 1) THEN
          ie = i
        ELSE
          EXIT
        END IF
      END DO

    ELSE

      ie = -1
      DO i = src_row_len, 1, -1
        IF (col_mask(i) == 1) THEN
          ie = i
          EXIT
        END IF
      END DO

    END IF

    IF (IE == -1) THEN
      CMESSAGE = ' IE = -1 ; Eastern Boundary for Mask not set.'
      ICODE = 64
! DEPENDS ON: Ereport
      CALL EReport (RoutineName, ICode, CMessage)
    END IF

    WRITE (6,*) ' Eastern Column for sub-area : ',ie

!   --------------------------------
!   Sub-Area has now been determined
!   Boundaries are JS, JN, IW and IE
!   --------------------------------

!   -------------------------------------------------
!   Next section needs checking that it was set up to
!   validate this program.
!   Could protect by PrintStatus for now.
!   -------------------------------------------------

    ALLOCATE (indata (src_row_len, src_rows) )

    DO j = 1, src_rows
      DO i = 1, src_row_len
        indata (i,j) = i*1000 + j
      END DO
    END DO

    IF (L_GM) THEN
      ALLOCATE (sub_area_data ( iw:ie+src_row_len, js:jn ) )
    ELSE
      ALLOCATE (sub_area_data ( iw:ie, js:jn ) )
    END IF

    IF (L_GM) THEN

      DO j = js, jn
        DO i = iw, src_row_len
          IF (mask(i,j) == 1) THEN
            sub_area_data (i,j) = indata (i,j)
          ELSE
            sub_area_data (i,j) = 0
          END IF
        END DO
      END DO

      DO j = js, jn
        DO i = 1, ie
          IF (mask(i,j) == 1) THEN
            sub_area_data (i+src_row_len,j) = indata (i,j)
          ELSE
            sub_area_data (i+src_row_len,j) = 0
          END IF
        END DO
      END DO

    ELSE

      DO j = js, jn
        DO i = iw, ie
          IF (mask(i,j) == 1) THEN
            sub_area_data (i,j) = indata (i,j)
          ELSE
            sub_area_data (i,j) = 0
          END IF
        END DO
      END DO

    END IF

    print_sub_area_data = .FALSE.
    IF (print_sub_area_data) THEN
      IF (L_GM) THEN
        WRITE (6,*) ' '
        DO i=iw,ie+src_row_len,10
          ijmin=MIN(i+9,ie+src_row_len)
          WRITE (6,*) ' ijmin ',ijmin
          WRITE (6,"('     ',10I7)") (ij,ij=i,ijmin)
          DO j=js,jn
            WRITE (6,"(' ',I4,10I7)") j,(sub_area_data(ij,j),ij=i,ijmin)
          END DO
          WRITE (6,*) ' '
        END DO

      ELSE

        WRITE (6,*) ' '
        DO i=iw,ie,10
          ijmin=MIN(i+9,ie)
          WRITE (6,*) ' ijmin ',ijmin
          WRITE (6,"('     ',10I7)") (ij,ij=i,ijmin)
          DO j=js,jn
            WRITE (6,"(' ',I4,10I7)") j,(sub_area_data(ij,j),ij=i,ijmin)
          END DO
          WRITE (6,*) ' '
        END DO

      END IF
    END IF

    DEALLOCATE (sub_area_data)
    DEALLOCATE (indata)

!   -------------------------------------------
!   End of section that may no longer be needed
!   -------------------------------------------


!   --------------------------
!   Statistics on data volumes
!   --------------------------

    npts_full_area = src_row_len * src_rows

    sub_area_rows = (jn-js+1)
    sub_area_rows_v = sub_area_rows - 1

    IF (L_GM) THEN
      sub_area_row_len = (ie+src_row_len-iw+1)
    ELSE
      sub_area_row_len = (ie-iw+1)
    END IF

    npts_sub_area    = sub_area_row_len * sub_area_rows
    npts_sub_area_v  = sub_area_row_len * sub_area_rows_v

    WRITE (6,*) ' '
    WRITE (6,*) ' DRIVING MODEL GRID '
    WRITE (6,*) ' Row Length ',src_row_len
    WRITE (6,*) ' Rows       ',src_rows
    WRITE (6,*) ' Levels     ',src_model_levels
    WRITE (6,*) ' Points (one level) ',npts_full_area
    WRITE (6,*) ' '
    WRITE (6,*) ' LBC GRID'
    WRITE (6,*) ' Max No of points in LBC grid ',max_lbc_size
    WRITE (6,*) ' '
    WRITE (6,*) ' SUB_AREA GRID'
    WRITE (6,*) ' Southern Row on driving grid ',js
    WRITE (6,*) ' Northern Row on driving grid ',jn
    WRITE (6,*) ' Western Col  on driving grid ',iw
    WRITE (6,*) ' Eastern Col  on driving grid ',ie
    WRITE (6,*) ' Row Length   ',sub_area_row_len
    WRITE (6,*) ' Rows         ',sub_area_rows
    WRITE (6,*) ' No of points ',npts_sub_area
    percentage =  &
       ( float(npts_sub_area) / float(npts_full_area) ) * 100.0
    WRITE (6,*) ' Percentage of driving grid ',percentage,' %'
    WRITE (6,*) ' '
    WRITE (6,*) ' LBC DATA POINTS ON DRIVING GRID '
    WRITE (6,*) ' No of points ',npts_mask
    percentage =  &
       ( float(npts_mask) / float(npts_full_area) ) * 100.0
    WRITE (6,*) ' Percentage of full     grid ',percentage,' %'
    percentage =  &
       ( float(npts_mask) / float(npts_sub_area) ) * 100.0
    WRITE (6,*) ' Percentage of sub_area grid ',percentage,' %'

!   -----------------
!   Get forecast time
!   -----------------

    CALL Fort_Get_Env ('FCTIME', 6, C_FC_Time, 5, ICode)
    READ(C_FC_Time,'(I4)') FC_Time

!   -------------------
!   Open the Frames FFs
!   -------------------

    CALL Fort_Get_Env ('UNIT20', 6, FileName, 100, ICode)
    len_env = LEN_TRIM(FileName)

    WRITE (6,*) ' Frames FFs : ',                   &
       FileName(1:LEN_TRIM(FileName))

! DEPENDS ON: File_Open
    CALL File_Open (frames_unit_no,                 &
       FileName, len_env, read_write, env_var_1, icode)

    WRITE (6,*) ' Frames FFs opened : icode ',icode

    IF (icode /= 0) THEN
      CMessage = 'Error in opening Frames FFs'
! DEPENDS ON: Ereport
      CALL Ereport ( RoutineName, icode, CMessage )
    END IF

!   What do the next 2 lines do ???

    FileName = 'NotSet'
    CALL Fort_Get_Env (env, 6, FileName, 100, ICode)

!   -----------------------------
!   Set up headers for Frames FFs
!   -----------------------------

    Frames_FixHD(:)  = IMDI
    Frames_IntHd(:)  = IMDI
    Frames_RealHd(:) = RMDI
    Frames_LookUp(1:45,:)  = IMDI
    Frames_Lookup(46:64,:) = TRANSFER ( IMDI, int_val)

    Frames_FixHD(1)  = 20
    Frames_FixHd(2)  = 1
    Frames_FixHd(3)  = 5
    Frames_FixHd(4)  = 0
    Frames_FixHd(5)  = 3
    Frames_FixHd(8)  = 1
    Frames_FixHd(9)  = 3
    Frames_Fixhd(12) = 601

    DO I = 21, 41
      Frames_Fixhd(I) = FixHd(I)
    END DO

    Frames_FixHd(100) = 257
    Frames_Fixhd(101) = 46
    Frames_FixHd(105) = Frames_FixHd(100) + Frames_FixHd(101)
    Frames_FixHd(106) = 38
    Frames_FixHd(110) = Frames_FixHd(105) + Frames_FixHd(106)
    Frames_FixHd(111) = Len1_LevDepC
    Frames_FixHd(112) = Len2_LevDepC
    Frames_FixHd(150) = Frames_FixHd(110) +                   &
       Frames_FixHd(111) * Frames_FixHd(112)
    Frames_FixHd(151) = 64
    Frames_FixHd(152) = 4096
    Frames_FixHd(160) = Frames_FixHd(150) +                   &
       Frames_FixHd(151) * Frames_FixHd(152)

    Frames_FixHd(160) = ((Frames_FixHd(160) + 2048 - 1)/      &
       2048)*2048

    Frames_FixHd(160) = Frames_FixHd(160) + 1

    Start_Address = Frames_FixHd(160) - 1

    Frames_FixHd(161) = 0 ! ?????????


    Frames_IntHd(6)  = IntHead(6)    !  # cols
    Frames_IntHd(7)  = IntHead(7)    !  # rows
    Frames_IntHd(8)  = IntHead(8)    !  Model levels
    Frames_IntHd(9)  = IntHead(9)    !  Wet levels
!     Frames_IntHd(12) = IntHead(12)   !  Tracer levels
    Frames_IntHd(12) = 0             !  Tracer levels
    Frames_IntHd(17) = IntHead(17)
    Frames_IntHd(24) = IntHead(24)

    Frames_RealHd(1)  = RealHead(1)
    Frames_RealHd(2)  = RealHead(2)
    Frames_RealHd(3)  = RealHead(3)
    Frames_RealHd(4)  = RealHead(4)
    Frames_RealHd(5)  = RealHead(5)
    Frames_RealHd(6)  = RealHead(6)
    Frames_RealHd(16) = RealHead(16)

!   ------------------------------------------------
!   Set up Frames_LevDepC to be same size as LevDepC
!   ------------------------------------------------

    ALLOCATE ( Frames_LevDepC (Len1_LevDepC, Len2_LevDepC) )

    Frames_LevDepC (:,:) = LevDepC (:,:)

!     do j=1,len2_levdepc
!     Do i=1,len1_levdepc  
!     write (6,*) ' i ',i,' j ',j,
!    &            ' LevDepC ',LevDepC (i,j),
!    &            ' Frames_LevDepC ',Frames_LevDepC(i,j)
!     enddo
!     enddo

    Frames_Lookup(  1:45, : ) = -99
    Frames_Lookup( 46:64, : ) = TRANSFER ( -99.0, int_val)

!   ------------------------------
!   Set up pointers to data fields
!   ------------------------------

    ipt(:)=0
    DO jvar = 1, Num_LBC_Vars

      Stash_Code = LBC_Variables (jvar) % Sect0_StashCode

      DO i = 1, n_lookups(0)
        IF (LookUp(42,i) == Stash_Code   .AND.  &
           LookUp(14,i) == FC_Time    ) THEN 
          ipt(jvar) = i
          EXIT
        END IF
      END DO

    END DO

    DO jvar = 1, Num_LBC_Vars
      WRITE (6,*) 'var ',jvar,' ipt ',ipt(jvar)
    END DO

    DO jvar = 1, Num_LBC_Vars
      IF (ipt(jvar) == 0) THEN
        ICode = 101
        CMessage = ' Pointer(s) not set. '
! DEPENDS ON: Ereport
        CALL Ereport ( RoutineName, icode, CMessage )
      END IF
    END DO

!   -----------------------------
!   Allocate space for input data
!   -----------------------------

    ALLOCATE ( Data_FF ( src_row_len, src_rows, src_model_levels+1 ) )

!   ------------------------------
!   Allocate space for FRAMES data
!   ------------------------------

    IF (L_GM) THEN
      ALLOCATE ( Data_Frames   ( iw:ie+src_row_len, js:jn,         &
                                 src_model_levels+1 )       )
    ELSE
      ALLOCATE ( Data_Frames   ( iw:ie, js:jn, src_model_levels+1 ) )
    END IF

!   ---------------------------------
!   Set up unit number for input data
!   ---------------------------------

    FF_Unit_No = 30 + LookUp_FF_No (ipt(1))
    WRITE (6,*) ' FF_Unit_No ',FF_Unit_No

!   -----------------------
!   Loop over LBC variables
!   -----------------------

    DO jvar = 1, Num_LBC_Vars

      fld_type  = lbc_variables (jvar) % fld_type
      l_u_field = (fld_type == fld_type_u)
      l_v_field = (fld_type == fld_type_v)
      nlevs     = lbc_variables (jvar) % n_levels

!     ----------------------
!     Read in input data
!     Data stored in Data_FF
!     ----------------------

      DO k = 1, nlevs

!       ----------------------------------------------------
!       WARNING : This routine must be called for each level
!       It will fail if trying to read in NLEVS levels.
!       Requires investigating.
!       ----------------------------------------------------

! DEPENDS ON: ReadFlds_Serial
        CALL ReadFlds_Serial (FF_unit_no,       &
                              1,                &
                              ipt(jvar)+k-1,    &
                              lookup,           &
                              64,               & 
                              Data_FF(1,1,k),   &
                              fixhd,            &
                              expand,           &
                              icode,            &
                              CMessage)

        IF (icode /= 0) THEN
          WRITE (6,*) ' Error in ReadFlds_Serial ?'
! DEPENDS ON: Ereport
          CALL Ereport ( RoutineName, icode, CMessage )
        END IF

      END DO  ! k

!     ------------
!     Extract data
!     ------------

      nlevs = lbc_variables (jvar) % n_levels

      DO k = 1, nlevs

        IF (L_GM) THEN

          DO j = js, jn

            DO i = iw, src_row_len

!               write (6,*) ' i ',i,' j ',j,
!    &                      ' lat ',
!    &                      src_first_lat +(j-1)*src_delta_lat,
!    &                      ' lon ',
!    &                      src_first_long+(i-1)*src_delta_long,
!    &                      ' data_ff ',data_ff(i,j,1)

              Data_Frames (i,j,k) = Data_FF(i,j,k)

            END DO

            DO i = 1, ie

!               write (6,*) ' i ',i,' j ',j,
!    &                      ' lat ',
!    &                      src_first_lat +(j-1)*src_delta_lat,
!    &                      ' lon ',
!    &                      src_first_long+(i-1)*src_delta_long,
!    &                      ' data_ff ',data_ff(i,j,1)

              Data_Frames (i+src_row_len,j,k) = Data_FF(i,j,k)

            END DO

          END DO

        ELSE

          DO j = js, jn
            DO i = iw, ie

!               write (6,*) ' i ',i,' j ',j,
!    &                      ' lat ',
!    &                      src_first_lat +(j-1)*src_delta_lat,
!    &                      ' lon ',
!    &                      src_first_long+(i-1)*src_delta_long,
!    &                      ' data_ff ',data_ff(i,j,1)

              Data_Frames (i,j,k) = Data_FF(i,j,k)

            END DO
          END DO

        END IF

      END DO

!     ------------------------------
!     Frames Data now in Data_Frames
!     ------------------------------

!     ------------------------------------------------------------
!     Apply mask - all data in sub-area not required for MakeBC is
!     to MASK_VALUE.
!     ------------------------------------------------------------

      l_apply_mask = .TRUE.
!     l_apply_mask = .false.

!     ---------------------------------------------------------
!     Yet to be validated : is RMDI the best value to use for 
!     packing the data prior to transmiision to external sites.
!     ---------------------------------------------------------

      mask_value = RMDI
!     mask_value = 0

      IF (l_apply_mask) THEN

        nlevs = lbc_variables (jvar) % n_levels

        DO k = 1, nlevs

          IF (L_GM) THEN

            DO j = js, jn
              DO i = iw, src_row_len

                IF ( Mask(i,j) == 0 ) THEN
                  Data_Frames (i,j,k) = Mask_Value
                END IF

              END DO
            END DO

            DO j = js, jn
              DO i = 1, ie

                IF ( Mask(i,j) == 0 ) THEN
                  Data_Frames (i+src_row_len,j,k) = Mask_Value
                END IF

              END DO
            END DO


          ELSE

            DO j = js, jn
              DO i = iw, ie

                IF ( Mask(i,j) == 0 ) THEN
                  Data_Frames (i,j,k) = Mask_Value
                END IF

              END DO
            END DO

          END IF

        END DO

      END IF

!     --------------------------
!     Set up FRAMES lookup table
!     --------------------------

!       --------------------------------------------------------
!       Yet to do :
!       - Use variables in CLOOKUP for position in LookUp header
!       - use UM_SECTOR_SIZE instead of hardwired 2048
!       - remove hard-wired UM version 6.1 in Word 38
!       --------------------------------------------------------

!     write (6,*) ' k0 ',k0

      DO k = 1, nlevs

        kk = k0 + k

        DO i = 1, 12
          Frames_Lookup (i,kk) = Lookup (i,ipt(jvar)+k-1)
        END DO

        Frames_Lookup (13,kk) = Lookup (13,ipt(jvar)+k-1)
        Frames_Lookup (14,kk) = Lookup (14,ipt(jvar)+k-1)
        IF (l_v_field) THEN
          Frames_Lookup (15,kk) = npts_sub_area_v
        ELSE
          Frames_Lookup (15,kk) = npts_sub_area
        END IF
        Frames_Lookup (16,kk) = Lookup (16,ipt(jvar)+k-1)

        IF (l_v_field) THEN
          Frames_Lookup (18,kk) = sub_area_rows_v
        ELSE
          Frames_Lookup (18,kk) = sub_area_rows
        END IF

        Frames_Lookup (19,kk) = sub_area_row_len
        Frames_Lookup (20,kk) = 0
        Frames_Lookup (21,kk) = 0
        Frames_Lookup (22,kk) = 2
        Frames_Lookup (23,kk) = Lookup (23,ipt(jvar)+k-1)
        Frames_Lookup (29,kk) = Start_Address

!       write (6,*) ' LookUp : Field No ',kk,
!    &  ' start_address ',Start_Address

        Frames_Lookup (30,kk) = ((Frames_Lookup (15,kk) + 2048 - 1)/ &
           2048)*2048

        Frames_Lookup (32,kk) = Lookup (32,ipt(jvar)+k-1)
        Frames_Lookup (33,kk) = Lookup (33,ipt(jvar)+k-1)
        Frames_Lookup (38,kk) = 6011111
        Frames_Lookup (39,kk) = 1
        Frames_Lookup (40,kk) = Start_Address
        Frames_Lookup (41,kk) = 0
        Frames_Lookup (42,kk) = Lookup (42,ipt(jvar)+k-1)
        Frames_Lookup (43,kk) = Lookup (43,ipt(jvar)+k-1)
        Frames_Lookup (44,kk) = Lookup (44,ipt(jvar)+k-1)
        Frames_Lookup (45,kk) = 1

        frames_first_lat  = src_first_lat  + (js-1)*src_delta_lat
        frames_first_long = src_first_long + (iw-1)*src_delta_long

        zeroth_lat  = frames_first_lat  - src_delta_lat
        zeroth_long = frames_first_long - src_delta_long

        IF ( l_u_field ) THEN
          zeroth_long = zeroth_long + 0.5 * src_delta_long
        END IF

        IF ( l_v_field ) THEN
! Bug in UM FFs ??, match for now.
!         zeroth_lat = zeroth_lat + 0.5 * src_delta_lat
          zeroth_lat = zeroth_lat - 0.5 * src_delta_lat
        END IF

!       write (6,*) ' frames_first_lat  ',frames_first_lat
!       write (6,*) ' frames_first_long ',frames_first_long
!       write (6,*) ' zeroth_lat        ',zeroth_lat
!       write (6,*) ' zeroth_long       ',zeroth_long

        Frames_Lookup (56,kk) = TRANSFER ( Frames_RealHd(5), int_val)
        Frames_Lookup (57,kk) = TRANSFER ( Frames_RealHd(6), int_val)

        Frames_Lookup (58,kk) = TRANSFER ( 0.0, int_val)
        Frames_Lookup (59,kk) = TRANSFER ( zeroth_lat, int_val)
        Frames_Lookup (60,kk) = TRANSFER ( src_delta_lat, int_val)
        Frames_Lookup (61,kk) = TRANSFER ( zeroth_long, int_val)
        Frames_Lookup (62,kk) = TRANSFER ( src_delta_long, int_val)
        Frames_Lookup (63,kk) = TRANSFER ( rmdi, int_val)
        Frames_Lookup (64,kk) = TRANSFER ( 0.0, int_val)

        Start_Address = Start_Address + Frames_Lookup (30,kk)

      END DO   !  k

!     -----------------------------------
!     Write extracted data to Frames FFs
!     ----------------------------------

      DO k = 1, nlevs

!       write (6,*) ' k0 + k ',k0+k

! DEPENDS ON: setpos
        CALL SetPos ( frames_unit_no, Frames_Lookup(29,k0+k), icode )

        IF (ICode /= 0) THEN
          CMessage = 'Error in SetPos for data in Frame FieldsFiles'
! DEPENDS ON: Ereport
          CALL Ereport ( RoutineName, icode, CMessage )
        END IF

! DEPENDS ON: buffout
        CALL buffout (frames_unit_no,                   &
                      Data_Frames(iw,js,k),             &
                      Frames_Lookup(30,k0+k),           &
                      len_io, a_io)

        IF (A_IO /= -1.0 .OR. LEN_IO /= Frames_Lookup(30,k0+k) ) THEN
! DEPENDS ON: IOERROR
          CALL IOERROR ('buffout of Frames data',       & 
                         A_IO, LEN_IO, Frames_Lookup(30,k0+k) )
          CMESSAGE = 'I/O ERROR with buffout of Frames data'
          ICode    = 11
! DEPENDS ON: Ereport
          CALL EReport (RoutineName, ICode, CMessage)
        END IF

      END DO

      k0 = k0 + nlevs

    END DO  !  jvar

!   -----------------
!   Write out headers
!   -----------------

!   -----
!   FixHd
!   -----

! DEPENDS ON: setpos
    CALL SetPos ( frames_unit_no, 0, icode )

    IF (ICode /= 0) THEN
      CMessage = 'Error in SetPos for Frames Fixed Header'
! DEPENDS ON: Ereport
      CALL Ereport ( RoutineName, icode, CMessage )
    END IF

! DEPENDS ON: buffout
    CALL buffout (frames_unit_no,Frames_FixHd,256,len_io,a_io)

    IF (A_IO /= -1.0 .OR. LEN_IO /= len_fixhd) THEN
! DEPENDS ON: IOERROR
      CALL IOERROR ('buffout of Frames fixed header',         &
                     A_IO, LEN_IO, len_fixhd)
      CMESSAGE = 'I/O ERROR with buffin of Integer Header'
      ICode    = 11
! DEPENDS ON: Ereport
      CALL EReport (RoutineName, ICode, CMessage)
    END IF

!   -----
!   IntHd
!   -----

! DEPENDS ON: setpos
    CALL SetPos ( frames_unit_no, Frames_FixHd(100)-1, icode )

    IF (ICode /= 0) THEN
      CMessage = 'Error in SetPos for Frames Integer Header'
! DEPENDS ON: Ereport
      CALL Ereport ( RoutineName, icode, CMessage )
    END IF

! DEPENDS ON: buffout
    CALL buffout (frames_unit_no,                       &
                  Frames_IntHd,                         &
                  Frames_FixHd(101),                    &
                  len_io, a_io)

    IF (A_IO /= -1.0 .OR. LEN_IO /= len_ihead) THEN
! DEPENDS ON: IOERROR
      CALL IOERROR ('buffout of Frames Integer header',      &
                     A_IO, LEN_IO, len_ihead)
      CMESSAGE = 'I/O ERROR with buffout of Frames Integer Header'
      ICode    = 11
! DEPENDS ON: Ereport
      CALL EReport (RoutineName, ICode, CMessage)
    END IF

!   ------
!   RealHd
!   ------

! DEPENDS ON: setpos
    CALL SetPos ( frames_unit_no, Frames_FixHd(105)-1, icode )

    IF (ICode /= 0) THEN
      CMessage = 'Error in SetPos for Frames Real Header'
! DEPENDS ON: Ereport
      CALL Ereport ( RoutineName, icode, CMessage )
    END IF

! DEPENDS ON: buffout
    CALL buffout (frames_unit_no,                        &
                  Frames_RealHd,                         &
                  Frames_FixHd(106),                     &
                  len_io, a_io)

    IF (A_IO /= -1.0 .OR. LEN_IO /= len_rhead) THEN
! DEPENDS ON: IOERROR
      CALL IOERROR ('buffout of Frames Real Header',     &
                     A_IO, LEN_IO, len_rhead)
      CMESSAGE = 'I/O ERROR with buffout of Frames Real Header'
      ICode    = 11
! DEPENDS ON: Ereport
      CALL EReport (RoutineName, ICode, CMessage)
    END IF

!   -------
!   LevDepC
!   -------

! DEPENDS ON: setpos
    CALL SetPos ( frames_unit_no, Frames_FixHd(110)-1, icode )

    IF (ICode /= 0) THEN
      CMessage = 'Error in SetPos for Frames LevDepC'
! DEPENDS ON: Ereport
      CALL Ereport ( RoutineName, icode, CMessage )
    END IF

! DEPENDS ON: buffout
    CALL buffout (frames_unit_no,                         &
                  Frames_LevDepC,                         &
                  Frames_FixHd(111)*Frames_Fixhd(112),    &
                  len_io, a_io)

    IF (A_IO /= -1.0 .OR. LEN_IO /= len_LevDepC) THEN
! DEPENDS ON: IOERROR
      CALL IOERROR ('buffout of Frames LevDepC header',   &
                     A_IO, LEN_IO, len_lookup)
      CMESSAGE = 'I/O ERROR with buffout of Frames LevDepC'
      ICode    = 11
! DEPENDS ON: Ereport
      CALL EReport (RoutineName, ICode, CMessage)
    END IF

!   ------
!   Lookup
!   ------

! DEPENDS ON: setpos
    CALL SetPos ( frames_unit_no, Frames_FixHd(150)-1, icode )

    IF (ICode /= 0) THEN
      CMessage = 'Error in SetPos for Frames LookUp'
! DEPENDS ON: Ereport
      CALL Ereport ( RoutineName, icode, CMessage )
    END IF

! DEPENDS ON: buffout
    CALL buffout (frames_unit_no,                          &
                  Frames_Lookup,                           &
                  Frames_FixHd(151)*Frames_FixHd(152),     &
                  len_io, a_io)

    IF (A_IO /= -1.0 .OR.                                  &
        LEN_IO /= Frames_FixHd(151)*Frames_FixHd(152) ) THEN
! DEPENDS ON: IOERROR
      CALL IOERROR ('buffout of Frames LookUp',            &
                     A_IO, LEN_IO,                         &
                     Frames_FixHd(151)*Frames_FixHd(152) )
      CMESSAGE = 'I/O ERROR with buffout of Frames LookUp'
      ICode    = 11
! DEPENDS ON: Ereport
      CALL EReport (RoutineName, ICode, CMessage)
    END IF


!   --------------------
!   Close the Frames FFs
!   --------------------

! DEPENDS ON: File_Close
    CALL File_Close (frames_unit_no,                    &
       FileName,len_env,env_var_1,keep_file,icode)
    IF (icode /= 0) THEN
      CMessage = 'Error in closing fieldsfiles'
! DEPENDS ON: Ereport
      CALL Ereport ( RoutineName, Icode, Cmessage )
    END IF

    DEALLOCATE (data_ff)
    DEALLOCATE (data_frames)
    DEALLOCATE (mask)
    DEALLOCATE (col_mask)
    DEALLOCATE (levdepc)
    DEALLOCATE (frames_levdepc)

    RETURN
  END SUBROUTINE LBC_Frames

! =============================================
! =============================================

  SUBROUTINE Get_StashMaster_Information ( src_model_levels,  &
                                           src_wet_levels )

!
! Description :
!       Get STASHmaster information Set up list of LBC variables for which Frames
!       data is required.
!
! Method :
!       For each LBC variable required, set up -
!       - section 32 STASH code
!       - corresponding section 0 STASH code
!       in LBC_Frames module.
!
!       Currently hard-wired to first 13 variables
!       in Section 32.
!
! Code Owner : Thomas Green

    USE LBC_Frames_mod 

    IMPLICIT NONE

    INTEGER :: src_model_levels
    INTEGER :: src_wet_levels

    INTEGER :: ErrorStatus
    CHARACTER (Len=80) :: CMessage
    CHARACTER (Len=*), PARAMETER ::  RoutineName =      &
       'Get_StashM_Info'

#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"

    INTEGER :: grid_type
    INTEGER :: level_type
    INTEGER :: first_level
    INTEGER :: last_level
    INTEGER :: i
    INTEGER :: row_number
    INTEGER :: item_code

    INTEGER :: top_level, bottom_level

    INTEGER, PARAMETER :: Unit_No_StashM  = 22
    INTEGER, PARAMETER :: Unit_No_UStashM = 2

    INTEGER, PARAMETER :: Sect_0  = 0
    INTEGER, PARAMETER :: Sect_32 = 32

    INTEGER get_fld_type
    INTEGER exppxi

    ErrorStatus = 0
    CMessage = ' '

    WRITE (6,*) ' ##############################'
    WRITE (6,*) ' Reading Read_Stashmaster_Files'
    WRITE (6,*) ' ##############################'

!   -------------------------------------------------
!   Determine number of records in Atmosphere SM File
!   -------------------------------------------------

    ppXrecs = 1
! DEPENDS ON: Hdppxrf
    CALL Hdppxrf (Unit_No_StashM, 'STASHmaster_A', ppxRecs, &
                  ErrorStatus, CMessage)

    WRITE (6,*) ' ppxRecs ',ppxrecs

!   Check errorstatus

!   ---------------------------
!   Read in STASHmaster records
!   ---------------------------

    row_number = 0
! DEPENDS ON: Getppx
    CALL Getppx ( Unit_No_StashM, Unit_No_UStashM,  &
                 'STASHmaster_A', row_number,       &
#include "argppx.h"
                  ErrorStatus, CMessage )

    WRITE (6,*) ' row_number ',row_number

!   Check errorstatus

    WRITE (6,*) ' ###############################'
    WRITE (6,*) ' Getting_Stashmaster_Information'
    WRITE (6,*) ' ###############################'

    DO i = 1, Num_LBC_Vars

!     -------------     
!     Get grid type
!     -------------

! DEPENDS ON: exppxi 
      grid_type = exppxi (                       &
                  1, Sect_32, i, ppx_grid_type,  &
#include "argppx.h"
                  ErrorStatus, CMessage )

!     -------------
!     Get halo type
!     -------------

! DEPENDS ON: exppxi
      lbc_variables (i) % halo_type = exppxi (   &
                  1, Sect_32, i, ppx_halo_type,  &
#include "argppx.h"
                  ErrorStatus, CMessage )

!     --------------------------------
!     Derive field type from grid type
!     --------------------------------

! DEPENDS ON: get_fld_type 
      lbc_variables (i) % fld_type = get_fld_type (grid_type)

!       write (6,*) i,
!    &  ' fld_type ',lbc_variables (i) % fld_type,
!    &  ' halo_type ',lbc_variables (i) % halo_type

      item_code = LBC_Variables (i)  % Sect0_Stashcode

!     --------------
!     Get level type
!     --------------

! DEPENDS ON: exppxi
      level_type = exppxi (                            &
                   1, Sect_0, item_code, ppx_lv_code,  &
#include "argppx.h"
                   ErrorStatus, CMessage )

!     ---------------
!     Get first level
!     ---------------


! DEPENDS ON: exppxi
      first_level = exppxi (                           &
                    1, Sect_0, item_code, ppx_lb_code, &
#include "argppx.h"
                    ErrorStatus, CMessage )

!     --------------
!     Get last level
!     --------------

! DEPENDS ON: exppxi
      last_level = exppxi (                            &
                   1, Sect_0, item_code, ppx_lt_code,  &
#include "argppx.h"
                   ErrorStatus, CMessage )

!     ---------------------------
!     Derive bottom and top level
!     ---------------------------

!       write (6,*) ' i ',i,' level_type  ',level_type,
!    &                      ' first_level ',first_level,
!    &                      ' last_level  ',last_level

      IF (level_type /= 5) THEN

!         Call LevCod (first_level, bottom_level, ErrorStatus, CMessage)
!         Call LevCod (last_level,  top_level,    ErrorStatus, CMessage)

!         LevCod too complicated for this program
!         Just set up what is required for now

        SELECT CASE (first_level )

        CASE (1)   ! First model level
          bottom_level = 1
        CASE (38)  ! Surface level on Charney-Phillips theta grid
          bottom_level = 0
        CASE Default

          WRITE (CMessage,*) ' First level ',first_level,   &
                             ' not catered for.'
          ErrorStatus = 10
! DEPENDS ON: Ereport
          CALL Ereport (RoutineName, ErrorStatus, CMessage)

        END SELECT

        SELECT CASE (last_level )

        CASE (2) 
!             top_level = Model_Levels
          top_level = src_model_levels
        CASE (3)  
!             top_level = Wet_Levels
          top_level = src_wet_levels
        CASE (19)
!             top_level = Model_Levels + 1
          top_level = src_model_levels + 1

        CASE Default

          WRITE (CMessage,*) ' Last level ',last_level,    &
                             ' not catered for.'
          ErrorStatus = 20
! DEPENDS ON: Ereport
          CALL Ereport (RoutineName, ErrorStatus, CMessage)

        END SELECT

      ELSE
        bottom_level = 1
        top_level    = 1
      END IF

      lbc_variables (i) % n_levels = top_level - bottom_level + 1
      lbc_variables (i) % bot_lev  = bottom_level
      lbc_variables (i) % top_lev  = top_level 

      WRITE (6,*) ' i ',i,' nlevels ',lbc_variables (i) % n_levels,  &
         ' bot_lev ',lbc_variables (i) % bot_lev,   & 
         ' top_lev ',lbc_variables (i) % top_lev

    END DO

  END SUBROUTINE Get_StashMaster_Information

! =============================================
! =============================================

  SUBROUTINE Get_LBC_Grid
!
! Description :
!       Determine grid information for LBC area
!
! Method :
!       For the LBC area, determine :
!       - the LBC grid size for each combination of
!         field and halo type.
!       - the grid size for the interpolated LBC grid
!         for each combination of field and halo type.
!
! Code Owner : Thomas Green

    IMPLICIT NONE

#include "cmaxsize.h"
#include "cintfa.h"
#include "parvars.h"
#include "parlbcs.h"

    INTEGER :: ihalo
    INTEGER :: ifld

    WRITE (6,*) ' ####################'
    WRITE (6,*) ' calling Get_LBC_Grid'
    WRITE (6,*) ' ####################'

    WRITE (6,*)
    WRITE (6,*) ' intf_row_length(1) ',intf_row_length(1)
    WRITE (6,*) ' intf_p_rows(1)     ',intf_p_rows(1)
    WRITE (6,*) ' intfwidtha(1)      ',intfwidtha(1)
    WRITE (6,*) ' intf_exthalo_ew(1) ',intf_exthalo_ew(1)
    WRITE (6,*) ' intf_exthalo_ns(1) ',intf_exthalo_ns(1)

! DEPENDS ON: LBC_Grid_Sizes
    CALL LBC_Grid_Sizes (1)

    WRITE (6,*) ' lbc_global_lenrima '
    DO ifld=1,nfld_max
      DO ihalo=1,nhalo_max
        WRITE (6,*) ' ihalo ',ihalo,' ifld ',ifld,     &
           lbc_global_lenrima(ifld,ihalo)
      END DO
    END DO

    WRITE (6,*) ' lbc_interp_lenrima '
    DO ifld=1,nfld_max
      DO ihalo=1,nhalo_max
        WRITE (6,*) ' ihalo ',ihalo,' ifld ',ifld,     &
           lbc_interp_lenrima(ifld,ihalo)
      END DO
    END DO

  END SUBROUTINE Get_LBC_Grid

! =============================================
! =============================================

  SUBROUTINE Get_LBC_List

!
! Description :
!       Set up list of LBC variables for which Frames
!       data is required.
!
! Method :
!       For each LBC variable required, set up -
!       - section 32 STASH code
!       - corresponding section 0 STASH code
!       in LBC_Frames module.
!
!       Currently set up for first 13 variables
!       in Section 32.
!
! Code Owner : Thomas Green


    USE LBC_Frames_mod

    IMPLICIT NONE

#include "csubmodl.h"
#include "cmaxsize.h"
#include "parvars.h"
#include "parlbcs.h"

    INTEGER, PARAMETER :: Sect_32 = 32

!   ------------------------------------------------
!   STASH Codes for corresponding Section 0 variable
!   ------------------------------------------------

!   (This may not be the final place for these
!    declarations - decide after MakeBC/Frames review)

    INTEGER, PARAMETER :: StashCode_orog    = 33
    INTEGER, PARAMETER :: StashCode_u       = 2
    INTEGER, PARAMETER :: StashCode_v       = 3
    INTEGER, PARAMETER :: StashCode_w       = 150
    INTEGER, PARAMETER :: StashCode_density = 253
    INTEGER, PARAMETER :: StashCode_theta   = 4
    INTEGER, PARAMETER :: StashCode_q       = 10
    INTEGER, PARAMETER :: StashCode_qcl     = 254
    INTEGER, PARAMETER :: StashCode_qcf     = 12
    INTEGER, PARAMETER :: StashCode_exner   = 255
    INTEGER, PARAMETER :: StashCode_u_adv   = 256
    INTEGER, PARAMETER :: StashCode_v_adv   = 257
    INTEGER, PARAMETER :: StashCode_w_adv   = 258

    INTEGER :: I

    WRITE (6,*) ' ####################'
    WRITE (6,*) ' calling Get_LBC_List'
    WRITE (6,*) ' ####################'

    LBC_Variables (:) % Model   = Atmos_IM
    LBC_Variables (:) % Section = Sect_32

!   -----------------------
!   Stashcode in Section 32
!   -----------------------

    LBC_Variables (1)  % Stashcode = lbc_stashcode_orog
    LBC_Variables (2)  % Stashcode = lbc_stashcode_u
    LBC_Variables (3)  % Stashcode = lbc_stashcode_v
    LBC_Variables (4)  % Stashcode = lbc_stashcode_w
    LBC_Variables (5)  % Stashcode = lbc_stashcode_density
    LBC_Variables (6)  % Stashcode = lbc_stashcode_theta
    LBC_Variables (7)  % Stashcode = lbc_stashcode_q
    LBC_Variables (8)  % Stashcode = lbc_stashcode_qcl
    LBC_Variables (9)  % Stashcode = lbc_stashcode_qcf
    LBC_Variables (10) % Stashcode = lbc_stashcode_exner
    LBC_Variables (11) % Stashcode = lbc_stashcode_u_adv
    LBC_Variables (12) % Stashcode = lbc_stashcode_v_adv
    LBC_Variables (13) % Stashcode = lbc_stashcode_w_adv

    WRITE (6,*) ' LBC_Var_StashCode '
    DO i = 1, Num_LBC_Vars
      WRITE (6,*) i,LBC_Variables (i) % Stashcode
    END DO

!   ------------------------------------
!   Corresponding Stashcode in Section 0
!   ------------------------------------

    LBC_Variables (1)  % Sect0_Stashcode = StashCode_orog
    LBC_Variables (2)  % Sect0_Stashcode = StashCode_u
    LBC_Variables (3)  % Sect0_Stashcode = StashCode_v
    LBC_Variables (4)  % Sect0_Stashcode = StashCode_w
    LBC_Variables (5)  % Sect0_Stashcode = StashCode_density
    LBC_Variables (6)  % Sect0_Stashcode = StashCode_theta
    LBC_Variables (7)  % Sect0_Stashcode = StashCode_q
    LBC_Variables (8)  % Sect0_Stashcode = StashCode_qcl
    LBC_Variables (9)  % Sect0_Stashcode = StashCode_qcf
    LBC_Variables (10) % Sect0_Stashcode = StashCode_exner
    LBC_Variables (11) % Sect0_Stashcode = StashCode_u_adv
    LBC_Variables (12) % Sect0_Stashcode = StashCode_v_adv
    LBC_Variables (13) % Sect0_Stashcode = StashCode_w_adv

    WRITE (6,*) ' LBC_Var_Sect0_StashCode'
    DO i = 1, Num_LBC_Vars
      WRITE (6,*) i,LBC_Variables (i) % Sect0_Stashcode
    END DO


  END SUBROUTINE Get_LBC_List

! =============================================
! =============================================

  SUBROUTINE Get_LBC_Areas

!
! Description :
!       Get LBC Area for which Frames data required
!
! Method :
!       Read in INTFCNSTA namelist and initialise
!       variables stored in TYPE LBC_GRID.
!       NB. Frames only works for one LBC area.
!
! Code Owner : Thomas Green


    IMPLICIT NONE

#include "cmaxsize.h"
#include "cintfa.h"
#include "cnaminfa.h"
#include "cprintst.h"

    TYPE Grid

       INTEGER :: Rows
       INTEGER :: Row_Length
       INTEGER :: RimWidth
       INTEGER :: HaloX
       INTEGER :: HaloY

       REAL :: First_Lat
       REAL :: First_Long
       REAL :: DeltaX
       REAL :: DeltaY
       REAL :: PoleLat
       REAL :: PoleLong

    END TYPE Grid

    TYPE (Grid) :: LBC_Grid

    INTEGER :: n_intf_a = 1
    INTEGER :: jintf
    INTEGER :: j

    WRITE (6,*) ' ####################'
    WRITE (6,*) ' calling Get_LBC_Area'
    WRITE (6,*) ' ####################'


    DO JINTF = 1,MAX_N_INTF_A
      A_INTF_START_HR(JINTF) = 0
      A_INTF_END_HR  (JINTF) = 0
      A_INTF_FREQ_HR (JINTF) = 0
      A_INTF_FREQ_MN (JINTF) = 0
      A_INTF_FREQ_SC (JINTF) = 0
      INTF_METH_LEV_CALC(JINTF) = 5
      INTF_MAX_SIG_HLEV(JINTF)  = 0
      INTF_MIN_PRS_HLEV(JINTF)  = 0
      INTF_PACK(JINTF) = 1
      LBC_ND(JINTF)    = 1
      INTF_PACK(JINTF) = 1
      INTF_VERTLEVS(JINTF)=' '
    END DO
    LBC_Q_MIN = 0.0

!   ---------------------------------------
!   Read INTFCNSTA namelist to get LBC grid
!   Same namelist used in UM and MakeBC.
!   ---------------------------------------

    REWIND (5)
    READ(5,INTFCNSTA)

!   Print out INTFCNTL namelist variables

    IF (PrintStatus >= PrStatus_Oper ) THEN

      WRITE (6,*) ' '
      WRITE (6,*) ' Namelist INTFCNSTA read in'
      DO jintf = 1, n_intf_a
        WRITE (6,*)' '
        WRITE (6,*) ' For area ',jintf
        WRITE (6,*) ' intf_p_rows      ',INTF_P_ROWS(JINTF)
        WRITE (6,*) ' intf_row_length  ',INTF_ROW_LENGTH(JINTF)
!         WRITE (6,*) ' intf_p_levels    ',INTF_P_LEVELS(JINTF)
!         WRITE (6,*) ' intf_q_levels    ',INTF_Q_LEVELS(JINTF)
!         WRITE (6,*) ' intf_tr_levels   ',INTF_TR_LEVELS(JINTF)
        WRITE (6,*) ' intf_firstlat    ',INTF_FIRSTLAT(JINTF)
        WRITE (6,*) ' intf_firstlong   ',INTF_FIRSTLONG(JINTF)
        WRITE (6,*) ' intf_nsspace     ',INTF_NSSPACE(JINTF)
        WRITE (6,*) ' intf_ewspace     ',INTF_EWSPACE(JINTF)
        WRITE (6,*) ' intf_polelat     ',INTF_POLELAT(JINTF)
        WRITE (6,*) ' intf_polelong    ',INTF_POLELONG(JINTF)
        WRITE (6,*) ' intfwidtha       ',INTFWIDTHA(JINTF)
        WRITE (6,*) ' intf_ExtHalo_NS  ',INTF_EXTHALO_NS(JINTF)
        WRITE (6,*) ' intf_ExtHalo_EW  ',INTF_EXTHALO_EW(JINTF)
      END DO
      WRITE (6,*) ' '
    END IF

    jintf = 1
    LBC_Grid % Rows       = Intf_P_Rows(jintf)
    LBC_Grid % Row_Length = Intf_Row_Length(jintf)
    LBC_Grid % First_Lat  = Intf_FirstLat(jintf)
    LBC_Grid % First_Long = Intf_FirstLong(jintf)
    LBC_Grid % DeltaX     = Intf_EWSpace(jintf)
    LBC_Grid % DeltaY     = Intf_NSSpace(jintf)
    LBC_Grid % PoleLat    = Intf_PoleLat(jintf)
    LBC_Grid % PoleLong   = Intf_PoleLong(jintf)
    LBC_Grid % RimWidth   = Intfwidtha(jintf)
    LBC_Grid % HaloX      = Intf_Exthalo_NS(jintf)
    LBC_Grid % HaloY      = Intf_Exthalo_EW(jintf)

    WRITE (6,*) ' LBC_Grid % Rows       ',LBC_Grid % Rows
    WRITE (6,*) ' LBC_Grid % Row_Length ',LBC_Grid % Row_Length
    WRITE (6,*) ' LBC_Grid % First_Lat  ',LBC_Grid % First_Lat
    WRITE (6,*) ' LBC_Grid % First_Long ',LBC_Grid % First_Long
    WRITE (6,*) ' LBC_Grid % DeltaX     ',LBC_Grid % DeltaX
    WRITE (6,*) ' LBC_Grid % DeltaY     ',LBC_Grid % DeltaY
    WRITE (6,*) ' LBC_Grid % PoleLat    ',LBC_Grid % PoleLat
    WRITE (6,*) ' LBC_Grid % PoleLong   ',LBC_Grid % PoleLong
    WRITE (6,*) ' LBC_Grid % RimWidth   ',LBC_Grid % RimWidth
    WRITE (6,*) ' LBC_Grid % HaloX      ',LBC_Grid % HaloX
    WRITE (6,*) ' LBC_Grid % HaloY      ',LBC_Grid % HaloY

  END SUBROUTINE Get_LBC_Areas

END PROGRAM FRAMES
#endif
