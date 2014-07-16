
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine: UM_INDEX-----------------------------------------------
!LL
!LL  Purpose: Calculate addresses of component arrays within a
!LL           series of super arrays, made up of combinations of arrays
!LL           that require dynamic allocation. Lengths of the super
!LL           arrays are calculated and passed into U_MODEL for
!LL           dynamic allocation, reducing the no. of arguments needed
!LL           to be passed between top-level routines.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 6.1.5A
!LL
!LL  Model            Modification history:
!LL version  date
!LL  3.2   30/03/93  Introduced as new DECK to allow dynamic allocation
!LL                  of main data arrays in U_MODEL.
!LL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
!LL                  1.Removes the limit on primary STASH item numbers.
!LL                  2.Removes the assumption that (section,item)
!LL                    defines the sub-model.
!LL                  3.Thus allows for user-prognostics.
!LL                  re-dimension PP_XREF and add INDEX_PPXREF.
!LL  3.5   29/03/95  mpp code : Land point fields are allocated P_FIELD
!LL                  amount of space in D1       P.Burton
!LL  3.5   Apr. 95   Sub-Models project.
!LL                  STASH super array modified in accordance with
!LL                  internal model separation scheme for diagnostics.
!LL                  *CALL CSUBMODL introduced.
!LL                  S.J.Swarbrick
!LL  4.0   06/09/95  Added atmos/ocean stash superarrays. K Rogers
!LL  4.1  15/03/96  Introduce Wave sub-model.  RTHBarnes.
!LL  4.1   04/12/95  Increased  A_IXPTR to accomodate 2 extra prognostic
!LL                  arrays J.Smith
!LL  4.1   26/04/96  Increased A_IXPTR to allow for 12 Sulphur Cycle
!LL                  prognostics and ancillaries        MJWoodage
!LL  4.2  11/10/96   Enable atmos-ocean coupling for mpp.
!LL                  (1): Coupled fields. Change to 'global' sizes
!LL                  instead of local. R.Rawlins
!LL  4.2   11/10/96  Enable atmos-ocean coupling for mpp.
!LL                  (2): Swap D1 memory. Add copies of D1 for atmos and
!LL                  ocean. R.Rawlins
!LL  4.3   26/03/97  Added HadCM2 sulphate loading pattern.  Will Ingram
!LL  4.4   01/07/97  Added padding to place the D1 array on a
!LL                  cache line boundary.
!LL                    Author: Bob Carruthers, Cray Research.
!LL  4.4   05/08/97  Allowed prognostic variable CCA to be 3D. JMG
!LL  4.4   11/10/97  Rename AO_D1_MEMORY to L_AO_D1_MEMORY. D Robinson
!LL  4.4   13/10/97  Initialise LEN_A/O/W_SPSTS. D. Robinson.
!LL  4.5   29/07/98  Use U_FIELD_INTFA to set up A_IXINF(21).
!LL                  Compute new A_IXINF(22-25). D. Robinson.
!LL  4.5   04/03/98  Increase IXPTR to allow for 1 new NH3 var in
!LL                  S Cycle and 3 new soot vars.     M Woodage
!LL  4.5   08/05/98  Increase A_IXPTR by 16 to increase maximum number
!LL                  of multi-level user ancillaries to 20
!LL                  Author D.M. Goddard
!LL  4.5   13/05/98  Added RHcrit variable to D1 pointers. S. Cusack
!LL  4.5   15/07/98  Added 3D CO2 to D1 pointers. C.D.Jones
!LL  4.5   17/08/98  Remove pointers for JSOIL_FLDS and JVEG_FLDS.
!LL                  D. Robinson.
!    4.5   30/03/98  Reserve space for land mask in ocean
!                    stash array. R. Hill
!LL  5.0   28/04/99  Remove references to FLOORFLDSA    P.Burton
!LL  5.0   20/05/99 New pointers for new dynamics data arrays added
!                   Redundant pointers removed.
!                   D.M. Goddard
!LL  5.0   21/05/99 Change variable names P.Burton
!    5.0   21/05/99  (1) STASH sizes now held separately.
!                    (2) Atmosphere STASH array changes. R.Rawlins
!LL  5.1  22/02/00  Add PARVARS for TYPSIZE                 P.Burton
!  5.1  28/04/00   Extra pressure level above top theta level
!                  p_top_theta_level removed                A.Malcolm
!    5.1   07/02/00 Extend super array ao_ixcpl() in order to include
!                   the extra arrays xva() & yva()
!    5.2   10/11/00 Extend a_ixinf to include lbc_eta_theta (26) and
!                   lbc_eta_rho (27). Rename MAX_INTF_P_LEVELS to
!                   MAX_INTF_MODEL_LEVELS. D.Robinson
!
!    5.3   19/06/01 Add tropopause-based ozone array. Dave Tan
!    5.3   01/10/01 Extend A_IXCON by 2 arrays for checkerboard
!                   radiation. S Cusack
!    5.3   15/10/01 Include more ocean boundary data. M J Bell
!    5.3   18/07/01 1) Remove unneeded allocation of 3 local array
!                   images of D1.
!                   2) Tidy: remove superfluous #if defined (mpp) and
!                   add test for writing out info messages.
!                   3) Remove innocuous array out of bounds message
!                   for a_spptr. R Rawlins
!    5.4  21/03/02 Remove comment on same line as #include
!                                               S. Carroll
!    5.5   17/02/03 Set pointers for Wave interface code. D.Holmes-Bell
!    5.5   17/02/03 Upgrade Wave model from 4.1 to 5.5, & incorporate
!                   parallel mods. D.Holmes-Bell
!    5.5  05/08/00  Modification for parallelisation of WAM.
!                   Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!    5.5  05/02/03  Add biomass aerosol arrays           P Davison
!    5.5  28/02/03  Set up Super array of atmos lat. and long.
!                   values for river routing if atmos only model
!                                           . C.Bunton
!    5.5  12/02/03  Include code for mineral dust scheme. S Woodward
!    6.0  30/07/03  Include cloud fraction lbcs. Damian Wilson
!    6.1  01/09/04  End of TIMER(UM_SETUP) moved.  R.Barnes
!    6.1  20/08/03  Add STOCHEM pointers.  C.Johnson
!    6.1  18/08/04  Set up A_IXINF(28-33). D Robinson.
!
!  6.1   04/08/04   Add dynamic arrays for diffusion coefficients
!                                                      Terry Davies
!    6.2  23/11/05  Removed all references to the wavemodel.
!                   T.Edwards
!    6.2  02/03/06  Added pointer for direct surface PAR flux.
!                                        M.G. Sanderson
!  6.2   14/02/06   Modify dynamic arrays for diffusion coefficients
!  6.2   14/02/06   Add dynamic arrays for variable resolution
!                                                      Terry Davies
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: C0
!LL
!LL  Project task: C0
!LL
!LL  External documentation: On-line UM document C1 - The top-level
!LL                          dynamic allocation
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!


!LL  Subroutine: UM_INDEX_A---------------------------------------------
!LL
!LL  Purpose: Calculate addresses and lengths of atmosphere super arrays
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 6.1.5A
!LL
!LL  Model            Modification history:
!LL version  date
!LL  3.2   30/03/93  Introduced as new DECK to allow dynamic allocation
!LL                  of main data arrays in U_MODEL.
!LL  3.2  22/11/93  Add 3 more A_IXPTR values for aerosol pointers.
!LL                                                   R.T.H.Barnes.
!LL  3.4  29/09/94  Add 6 more A_IXPTR values for multi-level murk and
!LL                  user ancillaries.                R.T.H.Barnes.
!LL  3.4    18/5/94 Add extra u field to A_SPCON (argcona). J Thomson
!LL  4.0  06/09/95  Added atmos stash superarray. K Rogers
!LL  5.0  09/06/99  New addressing for C-P C-grid derived constants for
!LL                 atmosphere, A_IXCON().  M L Gallani.
!LL  5.1  22/02/00  Add PARVARS for TYPSIZE                 P.Burton
!LL  5.1  26/04/00  Added LBC variables to A_IXPTR variables  P.Burton
!    5.2   01/09/00 - Remove levels from LBC pointers
!                   - Add orography LBC
!                   - Add extra space for LBC LOOKUPs
!                                                  P.Burton
!LL  5.2  13/09/00  Removed 3 aerosol pointers. P.Selwood
!LL  5.2  09/03/01  Replace eta_levels by zseak and Ck arrays for
!LL                 passing into STASH super-array. R Rawlins
!    5.4  02/09/02  Added sea mask as 12th element of a_ixsts K.Williams
!    5.5  03/02/03  Add 3 moisture variables and lbcs. R.M.Forbes
!    6.2  11/05/05  Correct jch4_stoch pointer. M.G. Sanderson.
!    6.2  01/10/04  Add murk aerosol lbc.   R.M.Forbes
!    6.2  05/01/06   Add true_latitude to A_IXCON. Yongming Tang
!    6.2  14/11/04  Add JTR_UKCA UKCA tracer pointer.   R Barnes
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: C0
!LL
!LL  Project task: C0
!LL
!LL  External documentation: On-line UM document C1 - The top-level
!LL                          dynamic allocation
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE UM_INDEX_A(                                            &
! ARGSZSPA super arrays lengths (atmosphere)
     &  A_SPDUM_LEN,A_SPPTR_LEN,A_SPCON_LEN,A_SPINF_LEN,A_SPANC_LEN,    &
     &  A_SPBND_LEN,A_SPSTS_LEN,                                        &
! ARGSZSPA end
     &              ICODE,CMESSAGE)
!
!*----------------------------------------------------------------------
      USE CSENARIO_MOD
      IMPLICIT NONE
!
!  Subroutines called
!
!  Local variables
!
      INTEGER ICODE             ! Work - Internal return code
      CHARACTER*80  CMESSAGE    ! Work - Internal error message
!
!  Configuration-dependent sizes for dynamic arrays
!
! ------------------------ Comdeck PARVARS -------------------------
! Parameters and common blocks required by the mpp-UM
! Added new comdeck AMAXSIZE required for new arrays in PARCOMM
! Add non-mpp option
!                                                      P.Burton
!========================== COMDECK PARPARM ====================
!   Description:
!
!   This COMDECK contains PARAMETERs for the mpp-UM
!
!   Two sets of parameters are set up -
!     i)  for the mpp-UM itself.
!     ii) for the interface to the Message Passing Software.
!
      !=================================================================
      ! Parameters needed for the mpp-UM
      !=================================================================
      ! maximum number of spatial dimensions
      INTEGER,PARAMETER:: Ndim_max = 3 ! 3d data

      ! number of different halo types
      INTEGER,PARAMETER:: NHalo_max = 3 ! for N.D. atmos. model

      INTEGER,PARAMETER:: halo_type_single   = 1
      INTEGER,PARAMETER:: halo_type_extended = 2
      INTEGER,PARAMETER:: halo_type_no_halo  = 3

! FLDTYPE definitions for the different field types recognised on the
! decomposition
      INTEGER,PARAMETER:: Nfld_max=7 ! maximum number of field types
      INTEGER,PARAMETER:: fld_type_p=1       ! grid on P points
      INTEGER,PARAMETER:: fld_type_u=2       ! grid on U points
      INTEGER,PARAMETER:: fld_type_v=3       ! grid on V points
      INTEGER,PARAMETER:: fld_type_comp_wave  = 4
                              ! Compressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_full_wave  = 5
                              ! Uncompressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_rim_wave   = 6
                              ! Boundary data for WAM Wave Field
      INTEGER,PARAMETER:: fld_type_r=7       ! grid on river points
      INTEGER,PARAMETER:: fld_type_unknown=-1! non-standard grid
! FLDTYPE end

      ! Used in addressing to indicate if calculation is for a local or
      ! global (ie. disk dump) size

      INTEGER,PARAMETER:: local_data=1
      INTEGER,PARAMETER:: global_dump_data=2

      ! maximum permitted size of a halo
      INTEGER,PARAMETER:: Max_Halo_Size=10

      !=================================================================
      ! Parameters needed for the Message Passing Software
      !=================================================================
      INTEGER,PARAMETER:: Maxproc = 512 ! Max number of processors

      ! Processor addresses in the neighbour array
      INTEGER,PARAMETER:: PNorth   = 1
      INTEGER,PARAMETER:: PEast    = 2
      INTEGER,PARAMETER:: PSouth   = 3
      INTEGER,PARAMETER:: PWest    = 4

      ! Value in neighbour array if the domain has  no neighbour in this
      ! direction. Otherwise the value will be the tid of the neighbor
      INTEGER,PARAMETER:: NoDomain = -1

      INTEGER,PARAMETER:: BC_STATIC   = 1 ! Static boundary conditions
      INTEGER,PARAMETER:: BC_CYCLIC   = 2 ! Cyclic boundary conditions
      INTEGER,PARAMETER:: BC_OVERPOLE = 3 ! Transfer over pole
! PARPARM end
!====================== COMDECK AMAXSIZE ========================
! Description
!   This comdeck provides parameters giving the maximum likely sizes
!   of key UM resolution variables, useful for sizing static arrays.
!
!   History:
!   Model    Date     Modification history
!  version
!   4.2      18/11/96 New comdeck created.  P.Burton
!   4.3      24/01/97 Define MaxFieldSize to be a quarter of the
!                     SHMEM common block size.         P.Burton
!   4.4      3/7/97   Add MaxFieldSizeMes. Deborah Salmond
!   4.5     12/01/98  Added new variables, and changed sizes to
!                     correspond to global hi-res forecast - current
!                     largest configuration.                P.Burton
!                     Changed MAX_SHMEM_COMMON_SIZE to 3000000
!                     required for operational data assimilation.
!                                                           P.Burton
!   5.0     29/04/99  Changed variable names:
!                       P_ROWS_MAX -> ROWS_MAX
!                       P_LEVELS_MAX -> MODEL_LEVELS_MAX
!                       Q_LEVELS_MAX -> WET_LEVELS_MAX
!                       MaxHaloSize -> MaxHaloArea
!                     Removed variable:
!                       HALO_MAX (use PARPARM Max_Halo_Size instead)
!    5.0   29/04/99  Remove mpp #define
!    5.3   05/12/01  Remove MaxFieldSize, MaxFieldSizeMes and
!                    Max3DFieldSize.  S.Cusack
!    5.5   22/01/03  Increase ROW_LENGTH_MAX and HORIZ_DIM_MAX
!                    from 432 to 548. D Robinson.
!    6.1   31/08/04  Allow up to 100 levels.  R.Barnes
!    6.2   13/02/06  Increase max values of row_length and
!                    rows to cope with FOAM high res, as well
!                    as Global N320 and NAE.  M Martin.
!    6.2   24/11/05  Use max function for horiz_dim_max. R Barnes
!    6.2     11/01/06 Remove max_shmem_common_size here and
!                     in rdobs2.   Camilla Mathison/R Barnes
!

! Maximum sector size for I/O
      INTEGER,PARAMETER:: IO_SECTOR_SIZE_MAX=4096
      INTEGER,PARAMETER:: ROW_LENGTH_MAX   = 840 ! Maximum row length
      INTEGER,PARAMETER:: ROWS_MAX         = 600 ! Max no of rows

      ! MAX(ROW_LENGTH_MAX,ROWS_MAX)
      INTEGER,PARAMETER:: HORIZ_DIM_MAX=MAX(ROW_LENGTH_MAX,ROWS_MAX)

      INTEGER,PARAMETER:: MODEL_LEVELS_MAX = 100 ! Max no of total levels
      INTEGER,PARAMETER:: WET_LEVELS_MAX   = 100 ! Max no of wet levels
      INTEGER, PARAMETER :: Max2DFieldSize = ROW_LENGTH_MAX*ROWS_MAX +  &
     &  IO_SECTOR_SIZE_MAX
      INTEGER, PARAMETER :: MaxHaloArea    = HORIZ_DIM_MAX*Max_Halo_Size
!========================== COMDECK PARCOMM ====================
!
! *** NOTE : This comdeck requires comdeck PARPARM to be *CALLed
!            first.
!
!   Description:
!
!   This COMDECK contains COMMON blocks for the mpp-UM
!
!
!   Two COMMON blocks are defined:
!     i)  UM_PARVAR holds information required by the
!         Parallel Unified Model itself
!     ii) MP_PARVAR holds information required by the interface to
!         the Message Passing Software used by the PUM
!
!   Key concepts used in the inline documentation are:
!     o global data - the entire data domain processed by the UM
!     o LOCAL data - the fragment of the global data which is
!       stored by this particular process
!     o PERSONAL data - the fragment of the LOCAL data which is
!       updated by this particular process
!     o HALO data - a halo around the PERSONAL data which forms
!       the LOCAL data
!
!     Acronyms used:
!     LPG - Logical Process Grid, this is the grid of logical
!           processors; each logical processor handles one of the
!           decomposed parts of the global data. It does not
!           necessarily represent a physical grid of processors.
!
!   History:
!
!   4.1      27/1/96  New comdeck based on second section of
!                     old PARVARS.   P.Burton
!   4.2     19/08/96  Removed some unused variables, and added
!                     current_decomp_type variable to allow use
!                     of flexible decompositions.
!                     Added nproc_max to indicate the max. number
!                     of processors used for mpp-UM
!                                                      P.Burton
!   5.0     12/04/99  - Added halosize array to common block
!                     - Added halo_i and halo_j to common block
!                     - Added fld_type dimension to glsize
!                     - Added halo type dimension to lasize
!                     - Added fld_type dimension to lasize
!                     - Replaced blsizep/blsizeu by blsize with
!                       extra fld_type dimension
!                     - Replace attop etc. with at_extremity
!                     - Added g_pe_index to common block
!                                                      P.Burton
!   5.1     22/05/00  Removed DATA statement and put in BLKDATA
!                                                      P.Burton
!   5.1     26/01/00  - Renamed g_pe_index -> g_pe_index_EW
!                     - Added g_pe_index_NS
!                                                     P.Burton
!   5.2     02/08/00  Added g_at_extremity        P.Burton
!   5.3     14/09/01  Added sb_model_domain   P.Burton
!   5.5     06/08/00  Modification for parallelisation of WAM.
!                     Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!   5.5     30/01/03  Generalised datastart   P.Selwood.
!   5.5     07/02/03  SX now uses PARCOMM instead of SXCOMM    E.Leung
!   6.0     18/09/03  F90-fy continuation lines.               P.Dando
!   6.2     23/11/05  Removed all references to the wavemodel.
!                     T.Edwards
!
! -------------------- COMMON BLOCKS --------------------
!
! =======================================================
! Common block for the Parallel Unified Model
! =======================================================

      INTEGER :: first_comp_pe       ! top left pe in LPG
      INTEGER :: last_comp_pe        ! bottom right pe in LPG
      INTEGER :: current_decomp_type ! current decomposition type
      INTEGER :: Offx                ! standard halo size in East-West
      INTEGER :: Offy                ! standard halo size in North-South
      INTEGER :: halo_i              ! extended halo size in East-West
      INTEGER :: halo_j              ! extended halo size in North-South
      INTEGER :: halosize(Ndim_max,NHalo_max) ! available halo sizes
      INTEGER :: glsize(Ndim_max,Nfld_max) ! global data size
      INTEGER :: lasize(Ndim_max,Nfld_max,NHalo_max) ! local data size
      INTEGER :: blsize(Ndim_max,Nfld_max) ! personal data size

      ! position of personal data in global data (in terms of standard
      ! Fortran array notation
      INTEGER :: datastart(Ndim_max)

      ! Generalised version of datastart for *all* fieldtypes
      INTEGER :: datastart_f(Ndim_max,Nfld_max)

      INTEGER :: gridsize(Ndim_max)  ! size of the LPG in each dimension

      ! position of this process in the LPG 0,1,2,...,nproc_x-1 etc.
      INTEGER :: gridpos(Ndim_max)

      INTEGER :: sb_model_domain


      ! logicals indicating if a processor is at the edge of the LPG
      LOGICAL :: at_extremity(4)

      COMMON /UM_PARVAR/                                                &
     &  first_comp_pe, last_comp_pe, current_decomp_type, Offx, Offy,   &
     &  halo_i, halo_j, halosize, glsize, lasize, blsize, datastart,    &
     &  datastart_f, gridsize, gridpos                                  &
     &,                 at_extremity,sb_model_domain

      ! Common block for the Message Passing Software

      ! type of boundary (cyclic or static) in each direction
      INTEGER :: bound(Ndim_max)

      ! global copy of local data size
      INTEGER :: g_lasize(Ndim_max,Nfld_max,NHalo_max,0:maxproc)

      ! global copy of personal data area
      INTEGER :: g_blsize(Ndim_max,Nfld_max,0:maxproc)

      ! global copy of datastart
      INTEGER :: g_datastart(Ndim_max,0:maxproc)

      ! global copy of datastart_f
      INTEGER :: g_datastart_f(Ndim_max,Nfld_max,0:maxproc)

      INTEGER :: g_gridpos(Ndim_max,0:maxproc) ! global copy of gridpos

      ! Which processor column a given point is in: 0 -> nproc_x-1
      INTEGER :: g_pe_index_EW(1-Max_Halo_Size:                         &
     &  ROW_LENGTH_MAX+Max_Halo_Size)

      ! Which processor row a given point is in: 0 -> nproc_y-1
      INTEGER :: g_pe_index_NS(1-Max_Halo_Size:ROWS_MAX+Max_Halo_Size)

      INTEGER :: nproc      ! number of processors in current decomp
      INTEGER :: mype      ! number of this processor (starting from 0)
      INTEGER :: nproc_max  ! maximum number of processors
      INTEGER :: nproc_x    ! number of processors in x-direction
      INTEGER :: nproc_y    ! number of processors in y-direction

      ! array with the tids of the four neighbours in the horizontal
      ! plane
      INTEGER :: neighbour(4)

      INTEGER :: gc_proc_row_group  ! GID for procs along a proc row
      INTEGER :: gc_proc_col_group  ! GID for procs along a proc col
      INTEGER :: gc_all_proc_group  ! GID for all procs

      ! at_extremity for each processor
      LOGICAL :: g_at_extremity(4,0:maxproc)

      COMMON /MP_PARVAR/                                                &
     &  bound,                                                          &
     &  g_lasize,g_blsize,g_datastart, g_datastart_f, g_gridpos,        &
     &  g_pe_index_EW,g_pe_index_NS,                                    &
     &  nproc_max,nproc_x,nproc_y,                                      &
     &  neighbour,gc_proc_row_group,                                    &
     &  gc_proc_col_group, gc_all_proc_group                            &
     &  ,nproc,mype                                                     &
     &, g_at_extremity

! PARCOMM end
! TYPSIZE start
!   Description:
!     This file contains sizes needed for dynamic allocation of
!   main data arrays within the model. Sizes read in from the user
!   interface via NAMELISTs are passed by /COMMON/. Other control
!   sizes that are fundamental in the definition of data structures
!   are assigned by PARAMETER statements.
!
      ! All sizes
      ! Not dependent on sub-model
      ! DATA IN NAMLST#x MEMBER OF THE JOB LIBRARY
      ! atmos START
      ! Main sizes of fields for each submodel
      ! Grid-related sizes for ATMOSPHERE submodel.
      INTEGER:: ROW_LENGTH           ! IN: No of points per local row
      INTEGER:: global_ROW_LENGTH    ! IN: Points per global row
      INTEGER:: ROWS                 ! IN: No of local (theta) rows
      INTEGER:: global_ROWS          ! IN: No of global (theta) rows
      INTEGER:: MODEL_LEVELS         ! IN: No of model levels
      INTEGER:: LAND_FIELD           ! IN: No of land points in field
      INTEGER:: NTILES               ! IN: No of land surface tiles

      ! Physics-related sizes for ATMOSPHERE submodel
      INTEGER:: WET_LEVELS          ! IN: No of moist-levels
      INTEGER:: CLOUD_LEVELS        ! IN: No of cloud-levels
      INTEGER:: ST_LEVELS           ! IN: No of soil temperature levels
      INTEGER:: SM_LEVELS           ! IN: No of soil moisture levels
      INTEGER:: BL_LEVELS           ! IN: No of boundary-layer-levels
      INTEGER :: OZONE_LEVELS       ! IN: No of ozone-levels
      INTEGER :: TPPS_OZONE_LEVELS  ! IN: No of tropopause-ozone-levels
      INTEGER :: RIVER_ROWS         ! IN: No of rows for river routing
      INTEGER :: RIVER_ROW_LENGTH   ! IN: Row length for river routing
      ! Dynamics-related sizes for ATMOSPHERE submodel

      INTEGER:: TR_LEVELS            ! IN: No of tracer-levels
      INTEGER:: TR_VARS              ! IN: No of passive tracers
      INTEGER:: TR_UKCA              ! IN: No of UKCA tracers

      ! Dynamics output diagnostic-related sizes for ATMOSPHERE submodel
      INTEGER:: THETA_PV_P_LEVS   ! IN: No of levels requested for pvort

      ! For Small executables
      INTEGER:: TOT_LEVELS
      ! Assimilation-related sizes for ATMOSPHERE submodel
      INTEGER :: N_AOBS           ! IN: No. of atmos observation types

      ! Grid related sizes for data structure
      ! Data structure sizes for ATMOSPHERE submodel
      INTEGER:: A_PROG_LOOKUP     ! IN: No of prognostic fields
      INTEGER:: A_PROG_LEN        ! IN: Total length of prog fields
      INTEGER:: A_LEN_INTHD       ! IN: Length of INTEGER header
      INTEGER:: A_LEN_REALHD      ! IN: Length of REAL header
      INTEGER:: A_LEN2_LEVDEPC    ! IN: No of LEVEL-dependent arrays
      INTEGER:: A_LEN2_ROWDEPC    ! IN: No of ROW-dependent arrays
      INTEGER:: A_LEN2_COLDEPC    ! IN: No of COLUMN-dependent arrays
      INTEGER:: A_LEN2_FLDDEPC    ! IN: No of FIELD arrays
      INTEGER:: A_LEN_EXTCNST     ! IN: No of EXTRA scalar constants
      INTEGER:: A_LEN_CFI1        ! IN: Length of compressed fld index 1
      INTEGER:: A_LEN_CFI2        ! IN: Length of compressed fld index 2
      INTEGER:: A_LEN_CFI3        ! IN: Length of compressed fld index 3
      ! atmos end

      ! OCEAN start
! TYPOCPAR start
!  History:
!  Version   Date     Comment
!  -------   ----     -------
!    4.4   15.06.97   Add free surface scalar R.Lenton
!     5.1   07.01.00   Invert_ocean to false with New Dynamics. JC Thil.
!     5.4   29.08.02   Add N_STRAIT and N_STRAIT_CLM. D. Storkey
!     5.5   03.01.03   Remove typocbas. R. Hill
!
      ! Grid related sizes for OCEAN model
      INTEGER::LSEG!IN:Max no of sets of start/end indices for vorticity
      INTEGER:: NISLE            ! IN: No of islands
      INTEGER:: ISEGM            ! IN: Max no of island segments per box
      INTEGER :: N_STRAIT       ! IN: No of pairs of Strait exchange pts
      INTEGER :: N_STRAIT_CLM   ! IN: No of Strait pts set by climate
                                !     values
      INTEGER:: O_LEN_COMPRESSED ! IN: No of ocean points in 3D field
      INTEGER:: LSEGC            ! IN: No of island basins for mead calc
      ! No of start/end indicies for the free surface solution
      INTEGER:: LSEGFS ! IN

      ! Fourier filtering for OCEAN submodel
      INTEGER :: LSEGF    ! IN: max. no of sets of indices for filtering
      INTEGER :: JFRST    ! IN: first J row of T to be filtered

      ! filtering is done on T with a low pass cut off to make the zonal
      ! dimension of the box filtered effectively the same as that of
      ! the boxes on row JFT0
      INTEGER :: JFT0     ! IN:

      INTEGER :: JFT1     ! IN: last J row of T in SH to be filtered
      INTEGER :: JFT2     ! IN: first J row of T in NH to be filtered
      INTEGER :: JFU0     ! IN: same function as JFT0 but for U,V
      INTEGER :: JFU1     ! IN: last J row of U,V in SH to be filtered
      INTEGER :: JFU2     ! IN: first J row of U,V in NH to be filtered

      ! Variables derived from those above
      INTEGER :: IMU     ! IN: total number of U,V grid boxes zonally
      INTEGER :: IMTP1   ! IN: IMT+1
      INTEGER :: IMTM1   ! IN: IMT-1
      INTEGER :: IMTM2   ! IN: IMT-2
      INTEGER :: IMUM1   ! IN: IMU-1
      INTEGER :: IMUM2   ! IN: IMU-2
      INTEGER :: JMTP1   ! IN: JMT+1
      INTEGER :: JMTM1   ! IN: JMT-1
      INTEGER :: JMTM2   ! IN: JMT-2
      INTEGER :: JSCAN   ! IN: JMTM2+1
      INTEGER :: KMP1    ! IN: KM+1
      INTEGER :: KMP2    ! IN: KM+2
      INTEGER :: KMM1    ! IN: KM-1
      INTEGER :: NSLAB   ! IN: no of words in one slab
      INTEGER :: JSKPT   ! IN: no of rows of T and U,V not filtered in
      INTEGER :: JSKPU   ! IN: low and mid latitudes + 1
      INTEGER :: NJTBFT  ! IN: no of J rows to be filtered on T
      INTEGER :: NJTBFU  ! IN: no of J rows to be filtered on U,V
      INTEGER :: IMTKM   ! IN: IMT*KM
      INTEGER :: NTMIN2  ! IN: maximum of NT or 2
      INTEGER :: IMTD2   ! IN: IMT/2
      INTEGER :: LQMSUM  ! IN: IMTD2*(IMT-IMTD2)
      INTEGER :: LHSUM   ! IN: IMT*IMTP1/2
      INTEGER :: IMTX8   ! IN: IMT*8
      INTEGER :: IMTIMT  ! IN: IMT*IMT
      INTEGER :: IMROT   ! X dimension for Coriolis array
      INTEGER :: JMROT   ! Y dimension for Coriolis array
      INTEGER :: IMBC    ! No of columns in boundary field array
      INTEGER :: JMBC    ! No of rows in boundary field array
      INTEGER :: KMBC    ! No of levels in boundary field array
      INTEGER :: NTBC    ! No of tracers in boundary field array
      INTEGER :: JMMD    ! No of rows for mead diagnostic basin indices
      INTEGER :: LDIV    ! No of divisions mead basin indices

      ! Grid-related switches for OCEAN submodel
      LOGICAL :: CYCLIC_OCEAN        ! IN: TRUE if CYCLIC E-W boundary
      LOGICAL :: GLOBAL_OCEAN        ! IN: TRUE if global domain

      ! TRUE if ocean grid NS-inverted cf atmos
      LOGICAL, PARAMETER :: INVERT_OCEAN=.FALSE.
      ! User interface limit for tracers
      ! Max no. tracers in STASHMASTER
      INTEGER, PARAMETER :: O_MAX_TRACERS=20
! COMOCPAR start
      COMMON /COMOCPAR/ GLOBAL_OCEAN, CYCLIC_OCEAN,                     &
     &  LSEG,NISLE,ISEGM,N_STRAIT,N_STRAIT_CLM,O_LEN_COMPRESSED,LSEGC,  &
     &  LSEGFS,LSEGF,JFRST,JFT0,JFT1,JFT2,JFU0,JFU1,JFU2,IMU,IMTP1,     &
     &  IMTM1,IMTM2,IMUM1,IMUM2,JMTP1,JMTM1,JMTM2,JSCAN,KMP1,KMP2,KMM1, &
     &  NSLAB,JSKPT,JSKPU,NJTBFT,NJTBFU,IMTKM,NTMIN2,                   &
     &  IMTD2,LQMSUM,LHSUM,IMTX8,IMTIMT,                                &
     &  IMROT,JMROT,IMBC,JMBC,KMBC,NTBC,JMMD,LDIV
! COMOCPAR end
! TYPOCPAR end
! TYPOCBAS Physics-related sizes for OCEAN submodel
      INTEGER ::  NT ! IN: No of ocean tracers (inc T,S)
      ! Grid related sizes for OCEAN model
      INTEGER :: IMT  ! IN: No of points per row (incl wrap)
      INTEGER :: JMT  ! IN: No of tracer rows
      INTEGER :: KM   ! IN: No of tracer levels
      INTEGER :: NT_UI     ! Copy of NT
      INTEGER :: IMT_UI    ! Copy of IMT
      INTEGER :: JMT_UI    ! Copy of JMT
      INTEGER :: KM_UI     ! Copy of KM
      INTEGER :: NICE      ! IN: No. of sea ice thickness categories
! COMOCBAS start
      COMMON /COMOCBAS/ NT_UI, IMT_UI, JMT_UI, KM_UI                    &
     &        ,NT, IMT, JMT, KM                                         &
     &                 ,NICE
! COMOCBAS end
! TYPOCBAS end
! TYPOASZ sizes for dynamic allocation of ocean assim.
! 5.2 11/08/00  JO_NMAX_OBS_ICE introduced. JO_MAX_OBS_VAL and
!               JO_NMAX_OBS_ICE put into COMOCASZ. M J Bell
      INTEGER :: JO_MAX_OBS_VAL !max number of values in OBS array

      ! max no of flds reqd at once in sea ice assim
      INTEGER :: JO_NMAX_OBS_ICE

      ! length of climate/covariances array
      INTEGER, PARAMETER:: JO_LEN_COV = 1

      ! max number of columns in climate grid
      INTEGER,PARAMETER:: JO_MAX_COLS_C = 1

      ! max number of rows in climate grid
      INTEGER,PARAMETER:: JO_MAX_ROWS_C = 1

      ! max number of levels in climate grid
      INTEGER,PARAMETER:: JO_MAX_LEVS_C = 1

! COMOCASZ start
      COMMON /COMOCASZ/ JO_MAX_OBS_VAL, JO_NMAX_OBS_ICE
! COMOCASZ end
! TYPOASZ end
      ! OCEAN end

      !  WAVE sub-model start
      ! WAVE end

      ! Grid related sizes for COUPLING between atmos and OCEAN
      ! submodels [For mpp, sizes are 'global' values over all
      ! PEs.Also needed for river routing]
      INTEGER:: AOCPL_IMT                ! Ocean rowlength
      INTEGER:: AOCPL_JMT                ! Ocean no. of rows
      INTEGER:: AOCPL_ROW_LENGTH         ! Atmos rowlength
      INTEGER:: AOCPL_P_ROWS             ! Atmos no. of p rows

      COMMON/SIZE_AOCPL/                                                &
     &  AOCPL_IMT, AOCPL_JMT, AOCPL_ROW_LENGTH, AOCPL_P_ROWS

      ! Data structure sizes for ATMOSPHERE ANCILLARY file control
      ! routines
      INTEGER :: NANCIL_LOOKUPSA  ! IN: Max no of fields to be read

      ! Data structure sizes for ATMOSPHERE INTERFACE file control
      ! routines
      INTEGER::N_INTF_A          ! No of atmosphere interface areas
      INTEGER::MAX_INTF_MODEL_LEVELS ! Max no of model levs in all areas
      INTEGER::MAX_LBCROW_LENGTH ! Max no of lbc row length in all areas
      INTEGER::MAX_LBCROWS ! Max no of lbc rows in all areas
      INTEGER::TOT_LEN_INTFA_P   ! Total length of interface p grids.
      INTEGER::TOT_LEN_INTFA_U    ! Total length of interface u grids.
      INTEGER::U_FIELD_INTFA      ! Length of Model U field (= U_FIELD)

      !  Data structure sizes for ATMOSPHERE BOUNDARY file control
      ! routines
! PARAMETERs defining the RIM_TYPE characteristics

!  History:
!  Date      Vn     Modification
!  31/10/01  5.3    Reset Nrima_max to 1. Change rima_type_orog
!                   from 2 to 1. D. Robinson

      INTEGER, PARAMETER :: Nrima_max = 1

      INTEGER, PARAMETER :: rima_type_norm=1  ! Normal field
      INTEGER, PARAMETER :: rima_type_orog=1  ! Orography field

! At 5.3 rima_type_orog=2 => rima_type_orog=1. This means that
! Orography LBCs has the same rim type as all prognostic LBCs. The
! value was changed rather than remove all occurences from the code to
! retain the functionality if it ever needs to be restored and to
! simplify the number of changes required.
      INTEGER :: RIMWIDTHA(Nrima_max)
      INTEGER :: NRIM_TIMESA      ! IN: Max no of timelevels in rim flds

      ! Data structure sizes for atmos & OCEAN BOUNDARY file control
      !routines

      ! Sizes applicable to all configurations (DUMPS/FIELDSFILE)

      INTEGER :: PP_LEN_INTHD   ! IN: Length of PP file integer header
      INTEGER :: PP_LEN_REALHD  ! IN: Length of PP file real    header

      ! Other sizes passed from namelist into common blocks
      COMMON/NLSIZES/                                                   &
     &  ROW_LENGTH,global_ROW_LENGTH,ROWS,global_ROWS,                  &
     &  LAND_FIELD,MODEL_LEVELS,WET_LEVELS,                             &
     &  NTILES,                                                         &
     &  CLOUD_LEVELS,TR_LEVELS,ST_LEVELS,SM_LEVELS,BL_LEVELS,           &
     &  OZONE_LEVELS,TPPS_OZONE_LEVELS,TR_VARS,TR_UKCA,                 &
     &  RIVER_ROWS, RIVER_ROW_LENGTH,                                   &

     &  THETA_PV_P_LEVS, N_AOBS,                                        &

     &  A_PROG_LOOKUP,A_PROG_LEN,                                       &
     &  A_LEN_INTHD,A_LEN_REALHD,                                       &
     &  A_LEN2_LEVDEPC,A_LEN2_ROWDEPC,A_LEN2_COLDEPC,                   &
     &  A_LEN2_FLDDEPC,A_LEN_EXTCNST,                                   &
     &  A_LEN_CFI1,A_LEN_CFI2,A_LEN_CFI3,                               &

     &  NANCIL_LOOKUPSA,                                                &

     &  N_INTF_A, MAX_INTF_MODEL_LEVELS, MAX_LBCROW_LENGTH,             &
     &  MAX_LBCROWS, TOT_LEN_INTFA_P,                                   &
     &  TOT_LEN_INTFA_U, U_FIELD_INTFA,                                 &

     &  RIMWIDTHA, NRIM_TIMESA,                                         &

     &  PP_LEN_INTHD,PP_LEN_REALHD

      !-----------------------------------------------------------------
      ! data in STASHC#x member of the job library

      ! Data structure sizes for ATMOSPHERE submodel (config dependent)
      INTEGER:: A_LEN2_LOOKUP   ! IN: Total no of fields (incl diags)
      INTEGER:: A_LEN_DATA      ! IN: Total no of words of data
      INTEGER:: A_LEN_D1        ! IN: Total no of words in atmos D1
      ! Data structure sizes for SLAB submodel (config dependent)
      INTEGER:: S_LEN2_LOOKUP   !IN: Tot no of fields (incl diags)
      INTEGER:: S_LEN_DATA      !IN: Tot no of words of data
      ! Data structure sizes for OCEAN submodel (config dependent)
      INTEGER:: O_LEN2_LOOKUP     ! IN: Total no of fields (incl diags)
      INTEGER:: O_LEN_DATA        ! IN: Total no of words of data
      INTEGER:: O_LEN_DUALDATA    ! IN: Words of data at 2 time levels
      INTEGER:: O_LEN_D1          ! IN: Total no of words in ocean D1
      ! Data structure sizes for WAVE submodel (config dependent)
      INTEGER:: W_LEN2_LOOKUP     ! IN: Total no of fields (incl diags)
      INTEGER:: W_LEN_DATA        ! IN: Total no of words of data
      INTEGER:: W_LEN_D1          ! IN: Total no of words in atmos D1

      ! Size of main data array for this configuration

      INTEGER:: LEN_TOT             ! IN: Length of D1 array
      INTEGER:: N_OBJ_D1_MAX         ! IN: No of objects in D1 array

      COMMON/STSIZES/                                                   &
     &  S_LEN2_LOOKUP,S_LEN_DATA,                                       &
     &  A_LEN2_LOOKUP,A_LEN_DATA,A_LEN_D1,                              &
     &  O_LEN2_LOOKUP,O_LEN_DATA,O_LEN_DUALDATA,O_LEN_D1,               &
     &  W_LEN2_LOOKUP,W_LEN_DATA,W_LEN_D1,                              &
     &  LEN_TOT,N_OBJ_D1_MAX
      ! global (ie. dump version) of *_LEN_DATA
      INTEGER:: global_A_LEN_DATA
      INTEGER:: global_O_LEN_DATA
      INTEGER :: global_W_LEN_DATA ! global (ie. dump version) of
                                   !                        W_LEN_DATA

      COMMON /MPP_STSIZES_extra/ global_A_LEN_DATA,global_O_LEN_DATA    &
     &      ,global_W_LEN_DATA
      ! Sizes of Stash Auxillary Arrays and associated index arrays
      ! Initialised in UMINDEX and UMINDEX_A/O/W
      INTEGER:: LEN_A_IXSTS
      INTEGER:: LEN_A_SPSTS
      INTEGER:: LEN_O_IXSTS
      INTEGER:: LEN_O_SPSTS
      INTEGER:: LEN_W_IXSTS
      INTEGER:: LEN_W_SPSTS

      COMMON /DSIZE_STS/                                                &
     &  LEN_A_IXSTS, LEN_A_SPSTS,LEN_O_IXSTS, LEN_O_SPSTS,              &
     &  LEN_W_IXSTS, LEN_W_SPSTS
!     From 4.5, the number of land points is computed for each
!     PE before the addressing section. All prognostics on land
!     points in the D1 space are now dimensioned by the local
!     no of land points rather than the global no of land points.

      INTEGER:: global_land_field    !  Global no of land points
      INTEGER:: local_land_field     !  Local no of land points
      COMMON /mpp_landpts/ global_land_field,local_land_field
      ! ----------------------------------------------------------------
      ! extra variables not passed through user interface

      ! fundamental data sizes :
      ! Fundamental parameter  sizes of data structure
      ! Sizes applicable to all configurations (HISTORY FILE)

      ! Length of history file in dump
      INTEGER, PARAMETER :: LEN_DUMPHIST = 0

      ! Sizes applicable to all configurations (DUMPS/FIELDSFILE)
      ! Length of dump fixed header
      INTEGER, PARAMETER :: LEN_FIXHD = 256

      ! Size of a single LOOKUP header
      INTEGER, PARAMETER :: LEN1_LOOKUP  = 64
      INTEGER, PARAMETER :: MPP_LEN1_LOOKUP= 2

      ! Size of compressed LBC LOOKUP (only used internally and
      ! contains just the items which change between each set of LBCs
      INTEGER, PARAMETER :: LEN1_LBC_COMP_LOOKUP = 8

      ! Sizes applicable to all configurations (STASH)
      ! Moved to <typstsz/typstsz.h>

      INTEGER:: INTF_LEN2_LEVDEPC !1st dim of interface out lev dep cons
      INTEGER:: INTF_LEN2_ROWDEPC !2nd dim of interface out Row dep cons
      INTEGER:: INTF_LEN2_COLDEPC !2nd dim of interface out Col dep cons
      
      COMMON /DSIZE/                                                    &
     &  INTF_LEN2_LEVDEPC,INTF_LEN2_ROWDEPC,INTF_LEN2_COLDEPC
      ! sub-model atmosphere   :
      ! Data structure sizes derived from grid size
      INTEGER:: A_LEN1_LEVDEPC ! IN: 1st dim of level  dep const
      INTEGER:: A_LEN1_ROWDEPC ! IN: 1st dim of row    dep const
      INTEGER:: A_LEN1_COLDEPC ! IN: 1st dim of column dep const
      INTEGER:: A_LEN1_FLDDEPC ! IN: 1st dim of field  dep const

      ! Data structure sizes for ATMOSPHERE INTERFACE file control
      ! routines
      INTEGER :: INTF_LOOKUPSA        ! No of interface lookups.
      INTEGER :: OLD_INTF_LOOKUPSA    ! No of interface lookups
                                      ! for old (4.5) LBCs.
      COMMON /DSIZE_A/                                                  &
     &  A_LEN1_LEVDEPC,A_LEN1_FLDDEPC,A_LEN1_ROWDEPC,A_LEN1_COLDEPC,    &
     &  INTF_LOOKUPSA,OLD_INTF_LOOKUPSA

      ! sub-model atmosphere   : derived sizes
      ! derived from model grid/levels. Arakawa B-grid

                                  ! Size of fields on THETA grid:
      INTEGER :: THETA_FIELD_SIZE     ! IN: with no halos
      INTEGER :: THETA_OFF_SIZE       ! IN: with simple halos
      INTEGER :: THETA_HALO_SIZE      ! IN: with extended halos

                                  ! Size of fields on U grid:
      INTEGER :: U_FIELD_SIZE         ! IN: with no halos
      INTEGER :: U_OFF_SIZE           ! IN: with simple halos
      INTEGER :: U_HALO_SIZE          ! IN: with extended halos

                                  ! Size of fields on V grid
      INTEGER :: V_FIELD_SIZE         ! IN: with no halos
      INTEGER :: V_OFF_SIZE           ! IN: with simple halos
      INTEGER :: V_HALO_SIZE          ! IN: with extended halos

      INTEGER :: N_ROWS               ! IN: No of V-rows
      INTEGER :: N_CCA_LEV            ! IN: No of CCA Levels
      COMMON/DRSIZE_A/                                                  &
     &  N_ROWS,N_CCA_LEV,THETA_FIELD_SIZE,U_FIELD_SIZE,V_FIELD_SIZE,    &
     &  THETA_OFF_SIZE,U_OFF_SIZE,V_OFF_SIZE,                           &
     &  THETA_HALO_SIZE,U_HALO_SIZE,V_HALO_SIZE
      ! boundary updating      : derived values
      ! Variables describing the Atmosphere Lateral Boundary Conditions
      ! Local (per processor) information

      ! Size of atmos LBC for given field type, halo type and rimwidth
      ! type
      INTEGER:: LENRIMA(Nfld_max,NHalo_max,Nrima_max)

      ! Size of given side (PNorth,PEast,PSouth and PWest), field type,
                        ! halo type and rimwidth type
      INTEGER:: LBC_SIZEA(4,Nfld_max,NHalo_max,Nrima_max)

                        ! Start of a given side within the LBC
      INTEGER:: LBC_STARTA(4,Nfld_max,NHalo_max,Nrima_max)

      ! Start of a given side within the LBC on a given processor
      INTEGER:: g_LBC_STARTA(4,Nfld_max,NHalo_max,Nrima_max,0:Maxproc-1)

      ! and global (within the file) information

      ! Size of atmos LBC on disk for given field type, halo type and
      ! rimwidth type
      INTEGER:: global_LENRIMA(Nfld_max,NHalo_max,Nrima_max)

                        ! Size of given side, field type and halo type
      INTEGER:: global_LBC_SIZEA(4,Nfld_max,NHalo_max,Nrima_max)

                        ! Start of a given side within the LBC
      INTEGER:: global_LBC_STARTA(4,Nfld_max,NHalo_max,Nrima_max)

      ! Variables describing the Ocean Lateral Boundary Conditions
      INTEGER:: LENRIMO                ! Size of ocean LBC (theta)
      INTEGER:: LENRIMO_U              ! Size of ocean LBC (velocity)

      ! Variables that may be needed for vn5.2 but have not yet been
      ! dealt with at vn5.1
      INTEGER:: RIMFLDSA
      INTEGER:: RIMFLDSO
      INTEGER:: global_LENRIMDATA_A
      INTEGER:: LENRIMDATA_A
      INTEGER:: LENRIMDATA_O
      INTEGER:: BOUNDFLDS
      INTEGER:: RIM_LOOKUPSA
      INTEGER:: RIM_LOOKUPSO
      INTEGER:: BOUND_LOOKUPSA
      INTEGER:: BOUND_LOOKUPSO
      INTEGER:: BOUND_LOOKUPSW
      INTEGER :: RIM_LOOKUPSW
      INTEGER :: LENRIMDATA_W
      INTEGER :: global_LENRIMDATA_W
      COMMON/DRSIZ_BO/                                                  &
      ! Atmosphere variables
     &  LENRIMA, LBC_SIZEA, LBC_STARTA, g_LBC_STARTA,                   &
     &  global_LENRIMA,global_LBC_SIZEA,global_LBC_STARTA,              &
      ! Wave model variables
     & RIM_LOOKUPSW, LENRIMDATA_W, global_LENRIMDATA_W,                 &
      ! Ocean variables
     &  LENRIMO, LENRIMO_U,                                             &
      ! Variables still to be dealt with
     &  RIMFLDSA,RIMFLDSO,BOUNDFLDS,RIM_LOOKUPSA,RIM_LOOKUPSO,          &
     &  BOUND_LOOKUPSA,BOUND_LOOKUPSO,BOUND_LOOKUPSW,                   &
     &  global_LENRIMDATA_A,                                            &
     &  LENRIMDATA_A,LENRIMDATA_O
! TYPSIZE end
!
! Parameters for constants arrays
! CMAXSIZE maximum sizes for dimensioning arrays
! of model constants whose sizes are configuration dependent. This
! allows constants to be read in from a NAMELIST file and maintain
! the flexibility of dynamic allocation for primary variables. The
! maximum sizes should agree with the maximum sizes implicit in the
! front-end User Interface.

!
!  Model            Modification history:
! version  Date
! 3.2  26/03/93  New COMDECK. Author R.Rawlins
! 3.4  06/08/94: Parameter MAX_NO_OF_SEGS used to dimension addresses
!                in macro-tasked calls to SWRAD, LWRAD & CONVECT.
!                Authors: A.Dickinson, D.Salmond, Reviewer: R.Barnes
! 3.5  22/05/95  Add MAX_N_INTF. D. Robinson
! 4.5  29/07/98  Increase MAX_N_INTF/MAX_N_INTF_A to 8. D. Robinson.
! 5.0  20/04/99  Changes for conversion to C-P C dynamics grid.
!                R. Rawlins
!  6.1   04/08/04  Add diffusion variable max_power     Terry Davies
! 6.2  25/12/05  Add max_updiff_levels/max_sponge_width   Terry Davies

      INTEGER,PARAMETER::max_model_levels = 100 ! Maximum no. of levels

      ! Max levels in boundary layer
      INTEGER,PARAMETER:: max_bl_levels = max_model_levels

      ! Max size of alpha_Cd
      INTEGER,PARAMETER :: max_number_alpha_cds = max_bl_levels

      ! Max no. of levels for pvort output
      INTEGER,PARAMETER :: MAX_REQ_THPV_LEVS = max_model_levels

      ! Max no. 1-2-1 rows in polar filter
      INTEGER,PARAMETER ::  max_121_rows =  8
      ! 0 is used for horizontal diffusion pointer

      ! Max no. of levels (from top) to apply upper level diffusion
      INTEGER,PARAMETER ::  max_updiff_levels = 10

      ! Max size of any sponge zones
      INTEGER,PARAMETER ::  max_sponge_width = 10

      ! Max size of look-up tables for searches
      INTEGER,PARAMETER ::  max_look = 2048

      ! Max no. of atmos interface areas
      INTEGER,PARAMETER :: MAX_N_INTF_A =  8

      ! Max no. of points in LBC      
      INTEGER,PARAMETER :: max_intf_lbcrow_length = 1000
      INTEGER,PARAMETER :: max_intf_lbcrows = 1000
        
      ! Max no. of atmos interface levels
      INTEGER,PARAMETER :: MAX_INTF_LEVELS = max_model_levels

      ! Maximum number of physics segments
      INTEGER,PARAMETER :: MAX_NO_OF_SEGS = 200
      ! MAX_N_INTF/MAX_N_INTF_A to be sorted out in later version
      ! Max no. of interface areas
      INTEGER, PARAMETER :: MAX_N_INTF =  8
! CMAXSIZE end
! CCONSTS start
! Description:
!    This file contains declarations for derived constants within
!   the atmospheric model. Where necessary PARAMETERS are defined to
!   dimension these constants. All constants are placed in the common
!   block CDERIVED, except hardwired constants, e.g. ETA_SPLIT and LENs.
!   file CMAXSIZE must be included first
!
!   The derived constants are calculated in the routine SETCONA1.
!
! PA, WI      <- programmer of some or all of previous code or changes
!
!  Model            Modification history from model version 3.0:
! version  Date
! 3.2   26/03/93  Remove resolution dependent variables for dynamic
!                 allocation. R.Rawlins
!   4.0   20/07/95  Sizes of tables expanded for version 3A
!                   of the radiation in sections 1 and 2. J.M. Edwards
!   5.0   21/06/99  Remove obsolete constants, for C-P C-grid dynamics.
!                   M L Gallani
!   5.1   25/02/00  Replace Data ETA_SPLIT by hard-wire h_split in
!                   Readlsta. Also remove obsolete references to
!                   radiation tables. R Rawlins
! Logical component: F011

      ! No of cloud types ie low/med/high
      INTEGER, PARAMETER :: NUM_CLOUD_TYPES = 3

      ! derived constants:
      INTEGER :: LOW_BOT_LEVEL      ! Bottom level of lowest cloud type
      INTEGER :: LOW_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: MED_BOT_LEVEL      ! Bottom   "    "  med      "    "
      INTEGER :: MED_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: HIGH_BOT_LEVEL     ! Bottom   "    "  top      "    "
      INTEGER :: HIGH_TOP_LEVEL     ! Top      "    "   "       "    "

      ! height values to split model levels into l/m/h cloud
      REAL ::    h_split(NUM_CLOUD_TYPES+1)

      LOGICAL :: ELF                ! T if atmosphere model on LAM grid

      ! Constants for dynamics output independent of resolution but
      ! dependent on choice of levels for output.
      REAL :: REQ_THETA_PV_LEVS(MAX_REQ_THPV_LEVS)

      COMMON /CDERIVED/                                                 &
     &  h_split,LOW_BOT_LEVEL,LOW_TOP_LEVEL,MED_BOT_LEVEL,MED_TOP_LEVEL,&
     &  HIGH_BOT_LEVEL, HIGH_TOP_LEVEL,ELF,REQ_THETA_PV_LEVS
! CCONSTS end
!  Ancillary file parameters for ancillary length calculations
! CONANC start
! History:
! Version  Date     Comment
! -------  ----     -------
!  3.4  04/10/94  Increase no.of datasets for murk & user ancil. RTHB
!  4.1  02/04/96  Introduce wave sub-model.  RTHBarnes.
!  4.1  26/02/96  Increase NANCIL_DATASETSA. D. Robinson.
!  4.3   18/3/97  Increase NANCIL_DATASETSA to 19.  William Ingram.
!  4.4   12/9/97  Increase NANCIL_DATASETSA to 22.  R.A.Betts
!  5.3   19/06/01 Increase NANCIL_DATASETSA to 25 (from 24) for
!                 tropopause-based ozone               Dave Tan
!  5.5   30/12/02 Increase NANCIL_DATASETSA to 28. Dave Robinson.
!  6.1   07/04/04 Increase NANCIL_DATASETSA to 29. Andy Jones
!  6.1   08/11/04 Increase NANCIL_DATASETSA to 32 for RR files. R.Sharp
!  7.3   31/05/10 Increase NANCIL_DATASETSA to 67 from 47 for
!                 tracer fluxes - kdcorbin
!
! To be used in conjunction with file TYPANC, defining dimensions
! of headers from ancillary files. A separate file is needed to
! allow calculation of super array sizes in UM_INDEX.
      INTEGER, PARAMETER :: NANCIL_DATASETSA=67
      INTEGER, PARAMETER :: NANCIL_DATASETSO=10
      INTEGER, PARAMETER :: NANCIL_DATASETSW=1
      INTEGER, PARAMETER :: NANCIL_DATASETS=NANCIL_DATASETSA+           &
     &  NANCIL_DATASETSO+NANCIL_DATASETSW
! CONANC end
!
!  Super array sizes for dynamic allocation in U_MODEL
!
! TYPSZSPA super arrays lengths (atmosphere)
      INTEGER :: A_SPDUM_LEN
      INTEGER :: A_SPPTR_LEN
      INTEGER :: A_SPCON_LEN
      INTEGER :: A_SPINF_LEN
      INTEGER :: A_SPANC_LEN
      INTEGER :: A_SPBND_LEN
      INTEGER :: A_SPSTS_LEN
! TYPSZSPA end
!
!
!  Addresses of arrays in super arrays.
!
!---------------------Start of SPINDEX--------------------------------
!L
!L --------------- D1    array  -----------------------------------
!L (now includes extra copies of atmos/ocean D1 for mpp)
!L
      INTEGER     IXD1_LEN       ! No. of arrays
      PARAMETER(  IXD1_LEN =     4)
      INTEGER     IXD1           ! Addresses of arrays
      COMMON/    CIXD1/      IXD1(IXD1_LEN)
!L
!L --------------- Dump headers--------------------------
!L
!L atmos
      INTEGER   A_IXDUM_LEN       ! No. of arrays
      PARAMETER(A_IXDUM_LEN = 14)
      INTEGER   A_IXDUM           ! Addresses of arrays
      COMMON/  CA_IXDUM/ A_IXDUM(A_IXDUM_LEN)
      INTEGER   A_IXSTS_LEN
      PARAMETER(A_IXSTS_LEN = 12)
      INTEGER   A_IXSTS
      COMMON/  CA_IXSTS/ A_IXSTS(A_IXSTS_LEN)
!L
!L --------------- STASH arrays -----------------------------------
!L
      INTEGER     IXSTS_LEN       ! No. of arrays
      PARAMETER(  IXSTS_LEN = 14)
      INTEGER     IXSTS           ! Addresses of arrays
      COMMON/    CIXSTS/      IXSTS(IXSTS_LEN)
!L
!L --------------- Pointers in D1 array and row/level dependent ---
!L --------------- constants
!L atmos
      INTEGER   A_IXPTR_LEN       ! No. of arrays
      PARAMETER(A_IXPTR_LEN = 184)  !CABLE_LAI
      INTEGER   A_IXPTR           ! Addresses of arrays
      COMMON/  CA_IXPTR/ A_IXPTR(A_IXPTR_LEN)
!L
!L --------------- Pre-calculated arrays of constants -------------
!L
!L atmos
      INTEGER   A_IXCON_LEN       ! No. of arrays
! A_IXCON_LEN is so low b/c the bulk of constants that previously were
! stored in um_index.a are now allocated directly in SETCONA 
      PARAMETER(A_IXCON_LEN = 3)
      INTEGER   A_IXCON           ! Addresses of arrays
      COMMON/  CA_IXCON/ A_IXCON(A_IXCON_LEN)
!L
!L --------------- Headers for output interface datasets (boundary
!L                 conditions out)
!L atmos
      INTEGER   A_IXINF_LEN       ! No. of arrays
      PARAMETER(A_IXINF_LEN = 35)
      INTEGER   A_IXINF           ! Addresses of arrays
      COMMON/  CA_IXINF/ A_IXINF(A_IXINF_LEN)
!L
!L --------------- Headers for ancillary files -------------------
!L
!L atmos
      INTEGER   A_IXANC_LEN       ! No. of arrays
      PARAMETER(A_IXANC_LEN = 4)
      INTEGER   A_IXANC           ! Addresses of arrays
      COMMON/  CA_IXANC/ A_IXANC(A_IXANC_LEN)
!L
!L --------------- Headers from input boundary files -------------
!L
!L NOT SUB-MODEL DEPENDENT
      INTEGER     IXBND_LEN       ! No. of arrays
      PARAMETER(  IXBND_LEN = 1)
      INTEGER     IXBND           ! Addresses of arrays
      COMMON/    CIXBND/ IXBND(IXBND_LEN)
!L
!L atmos
      INTEGER   A_IXBND_LEN       ! No. of arrays
      PARAMETER(A_IXBND_LEN = 5)
      INTEGER   A_IXBND           ! Addresses of arrays
      COMMON/  CA_IXBND/ A_IXBND(A_IXBND_LEN)
!L
!L --------------- Constant arrays needed for atmosphere-ocean----
!L --------------- coupling
!L
      INTEGER   AO_IXCPL_LEN      ! No. of arrays
      PARAMETER(AO_IXCPL_LEN = 10)
      INTEGER   AO_IXCPL          ! Addresses of arrays
      COMMON/ CAO_IXCPL/ AO_IXCPL(AO_IXCPL_LEN)
!L
!-------------------End of SPINDEX------------------------------------
!
!L----------------------------------------------------------------------
!L
      ICODE=0

!L
!L 1.0 Calculate size of atmosphere stash array
!L

! [Note that size of a_ixsts must be the same as o_ixsts and any other
!  submodels, since lower level STASH routines assume the existence of
!  _ixsts(i) even if not used.]

      a_ixsts(1) = 1                              ! zseak_rho
      a_ixsts(2) = a_ixsts(1) + model_levels      ! Ck_rho
      a_ixsts(3) = a_ixsts(2) + model_levels      ! zseak_theta
      a_ixsts(4) = a_ixsts(3) + model_levels + 1  ! Ck_theta
      a_ixsts(5) = a_ixsts(4) + model_levels + 1  ! dummy
      a_ixsts(6) = a_ixsts(5) + 1                 ! dummy
      a_ixsts(7) = a_ixsts(6) + 1                 ! p     pointer in D1
      a_ixsts(8) = a_ixsts(7) + 1                 ! pstar pointer in D1
      a_ixsts(9) = a_ixsts(8) + 1                 ! cos_v_latitude
                                                  !         (no halos)
      a_ixsts(10)= a_ixsts(9) + row_length*n_rows ! cos_theta_latitude
                                                  !         (no halos)
      a_ixsts(11)= a_ixsts(10)+ row_length*rows   ! landsea mask
      a_ixsts(12)= a_ixsts(11)+ row_length*rows   ! land fraction
      a_spsts_len= a_ixsts(12)+ row_length*rows
!     Store A_SPSTS_LEN in TYPSIZE
      LEN_A_SPSTS = A_SPSTS_LEN

!L
!L 1.1   DUMP     super array
!L
!L          super array addresses
      A_IXDUM(1) =1
      A_IXDUM(2) =A_IXDUM(1) + LEN_FIXHD
      A_IXDUM(3) =A_IXDUM(2) + A_LEN_INTHD
      A_IXDUM(4) =A_IXDUM(3) + A_LEN_CFI1+1
      A_IXDUM(5) =A_IXDUM(4) + A_LEN_CFI2+1
      A_IXDUM(6) =A_IXDUM(5) + A_LEN_CFI3+1
      A_IXDUM(7) =A_IXDUM(6) + A_LEN_REALHD
      A_IXDUM(8) =A_IXDUM(7) + A_LEN1_LEVDEPC*A_LEN2_LEVDEPC+1
      A_IXDUM(9) =A_IXDUM(8) + A_LEN1_ROWDEPC*A_LEN2_ROWDEPC+1
      A_IXDUM(10)=A_IXDUM(9) + A_LEN1_COLDEPC*A_LEN2_COLDEPC+1
      A_IXDUM(11)=A_IXDUM(10)+ A_LEN1_FLDDEPC*A_LEN2_FLDDEPC+1
      A_IXDUM(12)=A_IXDUM(11)+ A_LEN_EXTCNST+1
      A_IXDUM(13)=A_IXDUM(12)+ LEN_DUMPHIST+1
      A_IXDUM(14)=A_IXDUM(13)+ LEN1_LOOKUP*A_LEN2_LOOKUP

!L
!L          super array length
      A_SPDUM_LEN  =A_IXDUM(14)+ MPP_LEN1_LOOKUP*A_LEN2_LOOKUP
      A_SPDUM_LEN  =A_SPDUM_LEN -1
!L
!L
!L 1.2   Pointers super array
!L
!L          super array addresses

! Comment at end of each line corresponds to the matching
! pointer in ARGPTRA. eg A_IXPTR(3) in ARTPTRA = jtheta in ARGPTRA

! For each line : A_IXPTR(n+1) = A_IXPTR(n) + n_levs
! where n_levs in the no of levels for pointer n.

      A_IXPTR(1) =1                                  ! ju
      A_IXPTR(2) =A_IXPTR(1) + MODEL_LEVELS          ! jv
      A_IXPTR(3) =A_IXPTR(2) + MODEL_LEVELS          ! jw
      A_IXPTR(4) =A_IXPTR(3) + MODEL_LEVELS+1        ! jrho
      A_IXPTR(5) =A_IXPTR(4) + MODEL_LEVELS          ! jtheta
      A_IXPTR(6) =A_IXPTR(5) + MODEL_LEVELS          ! jq
      A_IXPTR(7) =A_IXPTR(6) + WET_LEVELS            ! jqcl
      A_IXPTR(8) =A_IXPTR(7) + WET_LEVELS            ! jqcf
      A_IXPTR(9) =A_IXPTR(8) + WET_LEVELS            ! jexner_rho_lvls
      A_IXPTR(10)=A_IXPTR(9) + MODEL_LEVELS+1         ! ju_adv
      A_IXPTR(11)=A_IXPTR(10)+ MODEL_LEVELS          ! jv_adv
      A_IXPTR(12)=A_IXPTR(11)+ MODEL_LEVELS          ! jw_adv
      A_IXPTR(13)=A_IXPTR(12)+ MODEL_LEVELS+1        ! jp
      A_IXPTR(14)=A_IXPTR(13)+ MODEL_LEVELS+1         ! jp_theta_levels
      A_IXPTR(15)=A_IXPTR(14)+ MODEL_LEVELS          ! jexner_theta_lvls
      A_IXPTR(16)=A_IXPTR(15)+ MODEL_LEVELS          ! jcca
      A_IXPTR(17)=A_IXPTR(16)+ N_CCA_LEV             ! jcf_area
      A_IXPTR(18)=A_IXPTR(17)+ WET_LEVELS            ! jcf_bulk
      A_IXPTR(19)=A_IXPTR(18)+ WET_LEVELS            ! jcf_liquid
      A_IXPTR(20)=A_IXPTR(19)+ WET_LEVELS            ! jcf_frozen
      A_IXPTR(21)=A_IXPTR(20)+ WET_LEVELS            ! j_deep_soil_temp
      A_IXPTR(22)=A_IXPTR(21)+ ST_LEVELS             ! jsmcl
      A_IXPTR(23)=A_IXPTR(22)+ SM_LEVELS             ! jsthu
      A_IXPTR(24)=A_IXPTR(23)+ SM_LEVELS             ! jsthf
      A_IXPTR(25)=A_IXPTR(24)+ SM_LEVELS             ! jsw_incs
      A_IXPTR(26)=A_IXPTR(25)+ MODEL_LEVELS+2        ! jlw_incs
      A_IXPTR(27)=A_IXPTR(26)+ MODEL_LEVELS+1        ! jozone
      A_IXPTR(28)=A_IXPTR(27)+ OZONE_LEVELS          ! jtracer
      A_IXPTR(29)=A_IXPTR(28)+ TR_LEVELS*(TR_VARS+1) ! jmurk
      A_IXPTR(30)=A_IXPTR(29)+ MODEL_LEVELS          ! jmurk_source
      A_IXPTR(31)=A_IXPTR(30)+ MODEL_LEVELS          ! jso2
      A_IXPTR(32)=A_IXPTR(31)+ MODEL_LEVELS          ! jdms
      A_IXPTR(33)=A_IXPTR(32)+ MODEL_LEVELS          ! jso4_aitken
      A_IXPTR(34)=A_IXPTR(33)+ MODEL_LEVELS          ! jso4_accu
      A_IXPTR(35)=A_IXPTR(34)+ MODEL_LEVELS          ! jso4_diss
      A_IXPTR(36)=A_IXPTR(35)+ MODEL_LEVELS          ! jh2o2
      A_IXPTR(37)=A_IXPTR(36)+ MODEL_LEVELS          ! jnh3
      A_IXPTR(38)=A_IXPTR(37)+ MODEL_LEVELS          ! jsoot_new
      A_IXPTR(39)=A_IXPTR(38)+ MODEL_LEVELS          ! jsoot_agd
      A_IXPTR(40)=A_IXPTR(39)+ MODEL_LEVELS          ! jsoot_cld
      A_IXPTR(41)=A_IXPTR(40)+ MODEL_LEVELS          ! jso2_natem
      A_IXPTR(42)=A_IXPTR(41)+ MODEL_LEVELS          ! joh
      A_IXPTR(43)=A_IXPTR(42)+ MODEL_LEVELS          ! jho2
      A_IXPTR(44)=A_IXPTR(43)+ MODEL_LEVELS          ! jh2o2_limit
      A_IXPTR(45)=A_IXPTR(44)+ MODEL_LEVELS          ! jo3_chem
      A_IXPTR(46)=A_IXPTR(45)+ MODEL_LEVELS          ! jhadcm2_so4
      A_IXPTR(47)=A_IXPTR(46)+ NSULPAT               ! jco2
      A_IXPTR(48)=A_IXPTR(47)+ MODEL_LEVELS          ! juser_mult1
      A_IXPTR(49)=A_IXPTR(48)+ MODEL_LEVELS          ! juser_mult2
      A_IXPTR(50)=A_IXPTR(49)+ MODEL_LEVELS          ! juser_mult3
      A_IXPTR(51)=A_IXPTR(50)+ MODEL_LEVELS          ! juser_mult4
      A_IXPTR(52)=A_IXPTR(51)+ MODEL_LEVELS          ! juser_mult5
      A_IXPTR(53)=A_IXPTR(52)+ MODEL_LEVELS          ! juser_mult6
      A_IXPTR(54)=A_IXPTR(53)+ MODEL_LEVELS          ! juser_mult7
      A_IXPTR(55)=A_IXPTR(54)+ MODEL_LEVELS          ! juser_mult8
      A_IXPTR(56)=A_IXPTR(55)+ MODEL_LEVELS          ! juser_mult9
      A_IXPTR(57)=A_IXPTR(56)+ MODEL_LEVELS          ! juser_mult10
      A_IXPTR(58)=A_IXPTR(57)+ MODEL_LEVELS          ! juser_mult11
      A_IXPTR(59)=A_IXPTR(58)+ MODEL_LEVELS          ! juser_mult12
      A_IXPTR(60)=A_IXPTR(59)+ MODEL_LEVELS          ! juser_mult13
      A_IXPTR(61)=A_IXPTR(60)+ MODEL_LEVELS          ! juser_mult14
      A_IXPTR(62)=A_IXPTR(61)+ MODEL_LEVELS          ! juser_mult15
      A_IXPTR(63)=A_IXPTR(62)+ MODEL_LEVELS          ! juser_mult16
      A_IXPTR(64)=A_IXPTR(63)+ MODEL_LEVELS          ! juser_mult17
      A_IXPTR(65)=A_IXPTR(64)+ MODEL_LEVELS          ! juser_mult18
      A_IXPTR(66)=A_IXPTR(65)+ MODEL_LEVELS          ! juser_mult19
      A_IXPTR(67)=A_IXPTR(66)+ MODEL_LEVELS          ! juser_mult20
      A_IXPTR(68)=A_IXPTR(67)+ MODEL_LEVELS          ! j_orog_lbc
      A_IXPTR(69)=A_IXPTR(68)+ 1                     ! ju_lbc
      A_IXPTR(70)=A_IXPTR(69)+ 1                     ! ju_lbc_tend
      A_IXPTR(71)=A_IXPTR(70)+ 1                     ! jv_lbc
      A_IXPTR(72)=A_IXPTR(71)+ 1                     ! jv_lbc_tend
      A_IXPTR(73)=A_IXPTR(72)+ 1                     ! jw_lbc
      A_IXPTR(74)=A_IXPTR(73)+ 1                     ! jw_lbc_tend
      A_IXPTR(75)=A_IXPTR(74)+ 1                     ! jrho_lbc
      A_IXPTR(76)=A_IXPTR(75)+ 1                     ! jrho_lbc_tend
      A_IXPTR(77)=A_IXPTR(76)+ 1                     ! jtheta_lbc
      A_IXPTR(78)=A_IXPTR(77)+ 1                     ! jtheta_lbc_tend
      A_IXPTR(79)=A_IXPTR(78)+ 1                     ! jq_lbc
      A_IXPTR(80)=A_IXPTR(79)+ 1                     ! jq_lbc_tend
      A_IXPTR(81)=A_IXPTR(80)+ 1                     ! jqcl_lbc
      A_IXPTR(82)=A_IXPTR(81)+ 1                     ! jqcl_lbc_tend
      A_IXPTR(83)=A_IXPTR(82)+ 1                     ! jqcf_lbc
      A_IXPTR(84)=A_IXPTR(83)+ 1                     ! jqcf_lbc_tend
      A_IXPTR(85)=A_IXPTR(84)+ 1                     ! jexner_lbc
      A_IXPTR(86)=A_IXPTR(85)+ 1                     ! jexner_lbc_tend
      A_IXPTR(87)=A_IXPTR(86)+ 1                     ! ju_adv_lbc
      A_IXPTR(88)=A_IXPTR(87)+ 1                     ! ju_adv_lbc_tend
      A_IXPTR(89)=A_IXPTR(88)+ 1                     ! jv_adv_lbc
      A_IXPTR(90)=A_IXPTR(89)+ 1                     ! jv_adv_lbc_tend
      A_IXPTR(91)=A_IXPTR(90)+ 1                     ! jw_adv_lbc
      A_IXPTR(92)=A_IXPTR(91)+ 1                     ! jw_adv_lbc_tend
      A_IXPTR(93)=A_IXPTR(92)+ 1                     ! jtracer_lbc
      A_IXPTR(94)=A_IXPTR(93)+ MAX(TR_VARS,1)        ! jtracer_lbc_tend
!     tracer_lbc (93) and tracer_lbc_tend (94) have lengths of TR_VARS
!     which may be zero so use MAX(TR_VARS,1) for lengths.

      A_IXPTR(95)=A_IXPTR(94)+ MAX(TR_VARS,1)        ! jtppsozone

!     tropopause ozone (95) has TPPS_OZONE_LEVELS which may be zero
!     so use MAX(TPPS_OZONE_LEVELS,1) for length.

      A_IXPTR(96)=A_IXPTR(95)+ MAX(TPPS_OZONE_LEVELS,1) ! jbmass_new
      A_IXPTR(97)=A_IXPTR(96)+ MODEL_LEVELS             ! jbmass_agd
      A_IXPTR(98)=A_IXPTR(97)+ MODEL_LEVELS             ! jbmass_cld
      A_IXPTR(99)=A_IXPTR(98)+ MODEL_LEVELS          ! jqcf2
      A_IXPTR(100)=A_IXPTR(99)  + WET_LEVELS         ! jqrain
      A_IXPTR(101)=A_IXPTR(100) + WET_LEVELS         ! jqgraup
      A_IXPTR(102)=A_IXPTR(101) + WET_LEVELS         ! jqcf2_lbc
      A_IXPTR(103)=A_IXPTR(102) + 1                  ! jqcf2_lbc_tend
      A_IXPTR(104)=A_IXPTR(103) + 1                  ! jqrain_lbc
      A_IXPTR(105)=A_IXPTR(104) + 1                  ! jqrain_lbc_tend
      A_IXPTR(106)=A_IXPTR(105) + 1                  ! jqgraup_lbc
      A_IXPTR(107)=A_IXPTR(106) + 1                  ! jqgraup_lbc_tend
      A_IXPTR(108)=A_IXPTR(107) + 1                  ! jdust_div1
      A_IXPTR(109)=A_IXPTR(108) + MODEL_LEVELS       ! jdust_div2
      A_IXPTR(110)=A_IXPTR(109) + MODEL_LEVELS       ! jdust_div3
      A_IXPTR(111)=A_IXPTR(110) + MODEL_LEVELS       ! jdust_div4
      A_IXPTR(112)=A_IXPTR(111) + MODEL_LEVELS       ! jdust_div5
      A_IXPTR(113)=A_IXPTR(112) + MODEL_LEVELS       ! jdust_div6
      A_IXPTR(114)=A_IXPTR(113) + MODEL_LEVELS     ! jcf_bulk_lbc
      A_IXPTR(115)=A_IXPTR(114) + 1                ! jcf_bulk_lbc_tend
      A_IXPTR(116)=A_IXPTR(115) + 1                ! jcf_liquid_lbc
      A_IXPTR(117)=A_IXPTR(116) + 1                ! jcf_liquid_lbc_tend
      A_IXPTR(118)=A_IXPTR(117) + 1                ! jcf_frozen_lbc
      A_IXPTR(119)=A_IXPTR(118) + 1                ! jcf_frozen_lbc_tend
      A_IXPTR(120)=A_IXPTR(119) + 1                ! jch4_stoch
      A_IXPTR(121)=A_IXPTR(120) + MODEL_LEVELS     ! jo3_stoch
      A_IXPTR(122)=A_IXPTR(121) + MODEL_LEVELS     ! jdirpar
      A_IXPTR(123)=A_IXPTR(122) + 1                     ! jtr_ukca
      A_IXPTR(124)=A_IXPTR(123) + TR_LEVELS*(TR_UKCA+1) ! jmurk_lbc
      A_IXPTR(125)=A_IXPTR(124) + 1                  ! jmurk_lbc_tend
      A_IXPTR(126)=A_IXPTR(125) + 1                  ! jOH_UKCA
      A_IXPTR(127)=A_IXPTR(126) + MODEL_LEVELS       ! jHO2_UKCA
      A_IXPTR(128)=A_IXPTR(127) + MODEL_LEVELS       ! jH2O2_UKCA
      A_IXPTR(129)=A_IXPTR(128) + MODEL_LEVELS       ! jO3_UKCA
      A_IXPTR(130)=A_IXPTR(129) + MODEL_LEVELS       ! jlcbase
      A_IXPTR(131)=A_IXPTR(130) + 1                  ! jccw_rad 
      A_IXPTR(132)=A_IXPTR(131) + WET_LEVELS         ! jozone_tracer
      A_IXPTR(133)=A_IXPTR(132) + MODEL_LEVELS       ! jO3_prod_loss
      A_IXPTR(134)=A_IXPTR(133) + MODEL_LEVELS       ! jO3_P_L_VMR
      A_IXPTR(135)=A_IXPTR(134) + MODEL_LEVELS       ! jO3_VMR
      A_IXPTR(136)=A_IXPTR(135) + MODEL_LEVELS       ! jO3_P_L_temp
      A_IXPTR(137)=A_IXPTR(136) + MODEL_LEVELS       ! jO3_temp
      A_IXPTR(138)=A_IXPTR(137) + MODEL_LEVELS       ! jO3_P_L_colO3
      A_IXPTR(139)=A_IXPTR(138) + MODEL_LEVELS       ! jO3_colO3
      A_IXPTR(140)=A_IXPTR(139) + MODEL_LEVELS       ! jarclbiog_bg
      A_IXPTR(141)=A_IXPTR(140) + MODEL_LEVELS       ! jarclbiom_fr
      A_IXPTR(142)=A_IXPTR(141) + MODEL_LEVELS       ! jarclbiom_ag
      A_IXPTR(143)=A_IXPTR(142) + MODEL_LEVELS       ! jarclbiom_ic
      A_IXPTR(144)=A_IXPTR(143) + MODEL_LEVELS       ! jarclblck_fr
      A_IXPTR(145)=A_IXPTR(144) + MODEL_LEVELS       ! jarclblck_ag
      A_IXPTR(146)=A_IXPTR(145) + MODEL_LEVELS       ! jarclsslt_fi
      A_IXPTR(147)=A_IXPTR(146) + MODEL_LEVELS       ! jarclsslt_je
      A_IXPTR(148)=A_IXPTR(147) + MODEL_LEVELS       ! jarclsulp_ac
      A_IXPTR(149)=A_IXPTR(148) + MODEL_LEVELS       ! jarclsulp_ak
      A_IXPTR(150)=A_IXPTR(149) + MODEL_LEVELS       ! jarclsulp_di
      A_IXPTR(151)=A_IXPTR(150) + MODEL_LEVELS       ! jarcldust_b1
      A_IXPTR(152)=A_IXPTR(151) + MODEL_LEVELS       ! jarcldust_b2
      A_IXPTR(153)=A_IXPTR(152) + MODEL_LEVELS       ! jarcldust_b3
      A_IXPTR(154)=A_IXPTR(153) + MODEL_LEVELS       ! jarcldust_b4
      A_IXPTR(155)=A_IXPTR(154) + MODEL_LEVELS       ! jarcldust_b5
      A_IXPTR(156)=A_IXPTR(155) + MODEL_LEVELS       ! jarcldust_b6
      A_IXPTR(157)=A_IXPTR(156) + MODEL_LEVELS       ! jarclocff_fr
      A_IXPTR(158)=A_IXPTR(157) + MODEL_LEVELS       ! jarclocff_ag
      A_IXPTR(159)=A_IXPTR(158) + MODEL_LEVELS       ! jarclocff_ic
      A_IXPTR(160)=A_IXPTR(159) + MODEL_LEVELS       ! jarcldlta_dl
      A_IXPTR(161)=A_IXPTR(160) + MODEL_LEVELS       ! jocff_new
      A_IXPTR(162)=A_IXPTR(161) + MODEL_LEVELS       ! jocff_agd
      A_IXPTR(163)=A_IXPTR(162) + MODEL_LEVELS       ! jocff_cld
      A_IXPTR(164)=A_IXPTR(163) + SM_LEVELS          ! jtsoil_tile 
      A_IXPTR(165)=A_IXPTR(164) + SM_LEVELS          ! jsmcl_tile 
!!! Not used MRD
!!!      A_IXPTR(166)=A_IXPTR(165) + SM_LEVELS          ! jsthu_tile 
      A_IXPTR(166)=A_IXPTR(165) + SM_LEVELS          ! jsthf_tile 
      A_IXPTR(167)=A_IXPTR(166) + 3                  ! jsnow_depth3l 
      A_IXPTR(168)=A_IXPTR(167) + 3                  ! jsnow_mass3l 
      A_IXPTR(169)=A_IXPTR(168) + 3                  ! jsnow_tmp3l 
      A_IXPTR(170)=A_IXPTR(169) + 3                  ! jsnow_rho3l 
      A_IXPTR(171)=A_IXPTR(170) + 1                  ! jsnow_rho1l  
      A_IXPTR(172)=A_IXPTR(171) + 1                  ! jsnow_age 
      A_IXPTR(173)=A_IXPTR(172) + 1                  ! jsnow_flg3l 
! Lestevens Sept 2012 - adding new prognostics
      A_IXPTR(174)=A_IXPTR(173) + 10                 ! jcpool_tile
      A_IXPTR(175)=A_IXPTR(174) + 10                 ! jnpool_tile
      A_IXPTR(176)=A_IXPTR(175) + 12                 ! jppool_tile
      A_IXPTR(177)=A_IXPTR(176) + 1                  ! jsoil_order
      A_IXPTR(178)=A_IXPTR(177) + 1                  ! jnidep
      A_IXPTR(179)=A_IXPTR(178) + 1                  ! jnifix
      A_IXPTR(180)=A_IXPTR(179) + 1                  ! jpwea
      A_IXPTR(181)=A_IXPTR(180) + 1                  ! jpdust
      A_IXPTR(182)=A_IXPTR(181) + 1                  ! jglai
      A_IXPTR(183)=A_IXPTR(182) + 1                  ! jphenphase
!CABLE_LAI
      A_IXPTR(184)=A_IXPTR(183) + 1                  ! jCABLE_LAI 
!     Super array length. If length of last variable could be zero, use
!     max(nlevs,1) to avoid out of bounds warning messages

!CABLE_LAI: updated starting index + model_levels
      A_SPPTR_LEN  =A_IXPTR(184) + MODEL_LEVELS
      A_SPPTR_LEN  =A_SPPTR_LEN -1
!L
!L
!L 1.3   Derived constants super array
!L
!L          super array addresses
!   The size of the increment on each line is actually the size of
!   the array referenced by the previous line,
!   e.g. (model_levels+1) is the size of eta_theta_levels.
      A_IXCON(1) = 1                         ! land_index
      A_IXCON(2) =A_IXCON(1) + land_field  ! land_ice_index
      A_IXCON(3) =A_IXCON(2) + land_field  ! soil_index
      ! other length info is now dealt w/ by ALLOCATE staments in SETCONA

      A_SPCON_LEN = A_IXCON(3) + land_field
!     A_SPCON_LEN  =A_SPCON_LEN -1
! sza: Above statement cause array out side of bounds if land_field = 0
! SWarbrick Sean suggested following, also fix typlndm.h in subroutine initdump  
      IF (A_SPCON_LEN > 1 ) A_SPCON_LEN = A_SPCON_LEN -1
!L
!L
!L 1.4   Interface output (boundary conditions) super array
!L
!L          super array addresses
      A_IXINF(1) =1
      A_IXINF(2) =A_IXINF(1) + LEN_FIXHD*N_INTF_A
      A_IXINF(3) =A_IXINF(2) + PP_LEN_INTHD*N_INTF_A
      A_IXINF(4) =A_IXINF(3) + LEN1_LOOKUP*INTF_LOOKUPSA*N_INTF_A
      A_IXINF(5) =A_IXINF(4) + PP_LEN_REALHD*N_INTF_A
      A_IXINF(6) =A_IXINF(5) +                                          &
     &            MAX_INTF_MODEL_LEVELS * INTF_LEN2_LEVDEPC * N_INTF_A
      A_IXINF(7) =A_IXINF(6) + TOT_LEN_INTFA_P
      A_IXINF(8) =A_IXINF(7) + TOT_LEN_INTFA_P
      A_IXINF(9) =A_IXINF(8) + TOT_LEN_INTFA_U
      A_IXINF(10)=A_IXINF(9) + TOT_LEN_INTFA_U
      A_IXINF(11)=A_IXINF(10)+ TOT_LEN_INTFA_P
      A_IXINF(12)=A_IXINF(11)+ TOT_LEN_INTFA_P
      A_IXINF(13)=A_IXINF(12)+ TOT_LEN_INTFA_P
      A_IXINF(14)=A_IXINF(13)+ TOT_LEN_INTFA_P
      A_IXINF(15)=A_IXINF(14)+ TOT_LEN_INTFA_U
      A_IXINF(16)=A_IXINF(15)+ TOT_LEN_INTFA_U
      A_IXINF(17)=A_IXINF(16)+ TOT_LEN_INTFA_U
      A_IXINF(18)=A_IXINF(17)+ TOT_LEN_INTFA_U
      A_IXINF(19)=A_IXINF(18)+ TOT_LEN_INTFA_U
      A_IXINF(20)=A_IXINF(19)+ TOT_LEN_INTFA_U
      A_IXINF(21)=A_IXINF(20)+ U_FIELD_INTFA
      A_IXINF(22)=A_IXINF(21)+ U_FIELD_INTFA
      A_IXINF(23)=A_IXINF(22)+ (MAX_INTF_MODEL_LEVELS+1)*N_INTF_A
      A_IXINF(24)=A_IXINF(23)+ (MAX_INTF_MODEL_LEVELS+1)*N_INTF_A
      A_IXINF(25)=A_IXINF(24)+  MAX_INTF_MODEL_LEVELS   *N_INTF_A
      A_IXINF(26)=A_IXINF(25)+  MAX_INTF_MODEL_LEVELS   *N_INTF_A
      A_IXINF(27)=A_IXINF(26)+ (MAX_INTF_MODEL_LEVELS+1)*N_INTF_A
!L
      A_IXINF(28)=A_IXINF(27)+ (MAX_INTF_MODEL_LEVELS*N_INTF_A)
      A_IXINF(29)=A_IXINF(28)+ TOT_LEN_INTFA_U
      A_IXINF(30)=A_IXINF(29)+ TOT_LEN_INTFA_U
      A_IXINF(31)=A_IXINF(30)+ TOT_LEN_INTFA_U
      A_IXINF(32)=A_IXINF(31)+ TOT_LEN_INTFA_U
      A_IXINF(33)=A_IXINF(32)+ TOT_LEN_INTFA_U
      A_IXINF(34)=A_IXINF(33)+ TOT_LEN_INTFA_U
      A_IXINF(35)=A_IXINF(34)+                                          &
     &            MAX_LBCROWS * INTF_LEN2_ROWDEPC * N_INTF_A
!L          super array length
!            [MAX( ,1) to prevent harmless out-of-bounds warnings
!                 when no boundary condition interfaces requested.]
      A_SPINF_LEN  =A_IXINF(35)+                                        &
     &         MAX(MAX_LBCROW_LENGTH * INTF_LEN2_COLDEPC * N_INTF_A ,1)

      A_SPINF_LEN  =A_SPINF_LEN -1

!L
!L 1.5   Ancillary file super array
!L
!L          super array addresses
      A_IXANC(1) =1
      A_IXANC(2) =A_IXANC(1) + LEN_FIXHD*NANCIL_DATASETSA
      A_IXANC(3) =A_IXANC(2) + A_LEN_INTHD*NANCIL_DATASETSA
      A_IXANC(4) =A_IXANC(3) + LEN1_LOOKUP*NANCIL_LOOKUPSA
!L
!L          super array length
      A_SPANC_LEN  =A_IXANC(4)+ A_LEN_REALHD*NANCIL_DATASETSA
      A_SPANC_LEN  =A_SPANC_LEN -1
!L
!L
!L 1.6   Input boundary constants super array
!L
!L          super array addresses
      A_IXBND(1) =1
      A_IXBND(2) =A_IXBND(1) + LEN_FIXHD
      A_IXBND(3) =A_IXBND(2) + A_LEN_INTHD
      A_IXBND(4) =A_IXBND(3) + LEN1_LOOKUP*RIM_LOOKUPSA
      A_IXBND(5) =A_IXBND(4) + LEN1_LBC_COMP_LOOKUP*BOUND_LOOKUPSA
!L
!L          super array length
       A_SPBND_LEN  =A_IXBND(5)+ A_LEN_REALHD
      A_SPBND_LEN  =A_SPBND_LEN -1
!L
!L----------------------------------------------------------------------
!L
      RETURN
      END SUBROUTINE UM_INDEX_A
