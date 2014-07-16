
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculates the Horizontal Interpolation Coefficients for LBCs.
!
! Subroutine Interface:

      Subroutine LBC_Interp_Coeffs (                                    &
     &     lbc_size                                                     &
     &,    src_row_len                                                  &
     &,    src_rows                                                     &
     &,    src_delta_lat                                                &
     &,    src_delta_long                                               &
     &,    src_first_lat                                                &
     &,    src_first_long                                               &
     &,    src_pole_lat                                                 &
     &,    src_pole_long                                                &
     &,    src_cyclic                                                   &
     &,    src_rotated                                                  &
     &,    lbc_row_len                                                  &
     &,    lbc_rows                                                     &
     &,    l_var_lbc                                                    &
     &,    lambda_in                                                    &
     &,    phi_in                                                       &
     &,    lbc_delta_lat                                                &
     &,    lbc_delta_long                                               &
     &,    lbc_first_lat                                                &
     &,    lbc_first_long                                               &
     &,    lbc_pole_lat                                                 &
     &,    lbc_pole_long                                                &
     &,    rimwidth                                                     &
     &,    lbc_halo_x                                                   &
     &,    lbc_halo_y                                                   &
     &,    lbc_index_bl                                                 &
     &,    lbc_index_br                                                 &
     &,    lbc_weights_tr                                               &
     &,    lbc_weights_br                                               &
     &,    lbc_weights_bl                                               &
     &,    lbc_weights_tl                                               &
     &,    lbc_coeff1                                                   &
     &,    lbc_coeff2                                                   &
     &,    i_uv                                                         &
     & )

      Implicit NONE
!
! Description:
!   Calculates the Horizontal Interpolation Coefficients for LBCs.
!
! Method:
!   1. Calculate lat/longs of lbc points on rotated grid.
!   2. Call EQTOLL to get corresponding true lat/longs.
!   3. Call W_COEFF to calculate coefficients to rotate the winds.
!   4. Set up lat/longs for the source grid.
!   5. Call H_INT_CO to calculate the horizontal interpolation
!      coefficients.
!
! Current Code Owner: Dave Robinson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2    13/11/00   Original code. Dave Robinson
!   5.3    09/10/01   New argument, src_cyclic. Correct calculation of
!                     lat/long for Northern side. Correct w_coeff
!                     argument list. D Robinson
!   5.5    06/08/02   New arguments - src_rotated, src_pole_lat and
!                     src_pole_long. Add call to LLTOEQ. D Robinson
!   6.1    28/07/04   Remove surplus checks on longitude ranges.
!                     D Robinson
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables :

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
!*L------------------ COMDECK CINTFA ----------------------------------
!L CMAXSIZE should be called first.
!
!   Contains Variables, Headers and Index blocks for control of
!   generation of boundary information for the limited area model.
!
!   Interfaces to all other models are handled by STASH, and there is
!   no explicit coding written for them in the model.
!
!Interface variables initialised through INTFCNSTA
!namelist read in the interface control routine INTF_CTL.
!L
!L 29/07/98  CINTF comdeck renamed to CINTFA. New arrays LBC_STREAM_A
!L           and LBC_UNIT_NO_A added. INTF_AK/BK/AKH/BKH removed - now
!L           in ARGINFA/TYPINFA. D. Robinson.
!L 10/11/00  5.2 Add Intf_ExtHalo_NS, Intf_ExtHalo_EW, Intf_RimW_Orog,
!L               LBC_ND, LBC_Z_TOP_MODEL, LBC_BL_LEVELS,
!L               LBC_FIRST_R_RHO and LBC_Q_MIN. D.Robinson
!  22/10/01  5.3 Remove Intf_RimW_Orog. D. Robinson
!L 18/09/01  5.3 Add A_INTF_FREQ_MN,A_INTF_FREQ_SC. Peter Clark
!L
      INTEGER                                                           &
     &  INTF_ROW_LENGTH                                                 &
                         ! Interface field row length
     & ,INTF_P_ROWS                                                     &
                         ! Interface field no of rows
     & ,INTF_P_LEVELS                                                   &
                         ! Interface field no of levels
     & ,INTF_Q_LEVELS                                                   &
                         ! Interface field no of wet levels
     & ,INTF_TR_LEVELS                                                  &
                         ! Interface field no of tracer levels
     & ,INTFWIDTHA                                                      &
                         ! Width of interface zone (atmosphere)
     & ,Intf_ExtHalo_NS                                                 &
                         ! Extended Halo in NS direction
     & ,Intf_ExtHalo_EW                                                 &
                         ! Extended Halo in EW direction
     & ,LBC_ND                                                          &
                         ! LBCs for old UM (0) or ND (1)
     & ,A_INTF_START_HR                                                 &
                         ! ) Start and End time in
     & ,A_INTF_FREQ_HR                                                  &
                         ! ) hours, Frequency in h,m,s for which
     & ,A_INTF_FREQ_MN                                                  &
                         ! ) atmosphere interface data
     & ,A_INTF_FREQ_SC                                                  &
                         ! ) is to be generated.
     & ,A_INTF_END_HR                                                   &
                         ! )
     & ,LEN_INTFA_P                                                     &
                         ! Length of interface p field
     & ,LEN_INTFA_U                                                     &
                         ! Length of interface u field
     & ,LEN_INTFA_DATA                                                  &
                         ! Length of interface data
     & ,INTF_PACK                                                       &
                         ! Packing Indicator for boundary data
     & ,LBC_STREAM_A                                                    &
                         ! Output streams in UMUI
     & ,LBC_UNIT_NO_A                                                   &
                         ! Unit Nos for Atmos Boundary Dataset
     & ,LBC_FIRST_R_RHO                                                 &
                         ! First rho level at which height is constant
     & ,LBC_BL_LEVELS                                                   &
                         ! No of Boundary Layer levels
!
! Following 3 variables not in common ; in namelist
     & ,INTF_METH_LEV_CALC(MAX_N_INTF_A)                                &
!                              !Method of calculating Eta level (ETAK)
!                              !from layers (ETAH)
     & ,INTF_MAX_SIG_HLEV(MAX_N_INTF_A)                                 &
!                              !level below which sigma coordinates used
     & ,INTF_MIN_PRS_HLEV(MAX_N_INTF_A)
!                              !level above which pressure coordinates

      REAL                                                              &
     &  INTF_EWSPACE                                                    &
                         ! E-W grid spacing (degrees)
     & ,INTF_NSSPACE                                                    &
                         ! N-S grid spacing (degrees)
     & ,INTF_FIRSTLAT                                                   &
                         ! Latitude of first row (degrees)
     & ,INTF_FIRSTLONG                                                  &
                         ! Longitude of first row (degrees)
     & ,INTF_POLELAT                                                    &
                         ! Real latitude of coordinate pole (degrees)
     & ,INTF_POLELONG                                                   &
                         ! Real longitude of coordinate pole (degrees)
     & ,LBC_Z_TOP_MODEL                                                 &
                         ! Height of top of model
     & ,LBC_Q_MIN                                                       &
                         ! Minimum value for q
!
! VarRes grid spacing
      , LAMBDA_INTF_P(MAX_INTF_LBCROW_LENGTH, MAX_N_INTF_A)             &
      , LAMBDA_INTF_U(MAX_INTF_LBCROW_LENGTH, MAX_N_INTF_A)             &    
      , PHI_INTF_P(MAX_INTF_LBCROWS, MAX_N_INTF_A)                      &
      , PHI_INTF_V(MAX_INTF_LBCROWS, MAX_N_INTF_A)                      &
      , LAMBDA_LBC_P(MAX_INTF_LBCROW_LENGTH)                            &
      , LAMBDA_LBC_U(MAX_INTF_LBCROW_LENGTH)                            &    
      , PHI_LBC_P(MAX_INTF_LBCROWS)                                     &
      , PHI_LBC_V(MAX_INTF_LBCROWS)

! Following variable not in common ; in namelist
      REAL INTF_ETAH(MAX_INTF_LEVELS+1,MAX_N_INTF_A)
!                          !Eta values at model layer boundaries ETAKH

      LOGICAL                                                           &
     &  INTF_VERT_INTERP                                                &
                         ! Switch to request vertical interpolation
     & ,LNEWBND          ! True for initialising new boundary data file

! Switch for variable resolution LBC output
      LOGICAL  INTF_L_VAR_LBC(MAX_N_INTF_A)

! Files for VERTLEVS namelist     
      Character(Len=80) :: INTF_VERTLEVS

! Files for HorzGrid namelist  
      Character(Len=80) :: INTF_HorzGrid(MAX_N_INTF_A)
!*----------------------------------------------------------------------
      COMMON /INTFCTL_ATMOS/                                            &
     &  INTF_EWSPACE(MAX_N_INTF_A)    ,INTF_NSSPACE(MAX_N_INTF_A)       &
     & ,INTF_FIRSTLAT(MAX_N_INTF_A)   ,INTF_FIRSTLONG(MAX_N_INTF_A)     &
     & ,INTF_POLELAT(MAX_N_INTF_A)    ,INTF_POLELONG(MAX_N_INTF_A)      &
     & ,INTF_ROW_LENGTH(MAX_N_INTF_A) ,INTF_P_ROWS(MAX_N_INTF_A)        &
     & ,INTF_P_LEVELS(MAX_N_INTF_A)   ,INTF_Q_LEVELS(MAX_N_INTF_A)      &
     & ,INTF_TR_LEVELS(MAX_N_INTF_A)  ,INTFWIDTHA(MAX_N_INTF_A)         &
     & ,Intf_ExtHalo_NS(Max_N_Intf_A) ,Intf_ExtHalo_EW(Max_N_Intf_A)    &
     & ,LBC_ND(Max_N_Intf_A)                                            &
     & ,A_INTF_START_HR(MAX_N_INTF_A) ,A_INTF_FREQ_HR(MAX_N_INTF_A)     &
     & ,A_INTF_FREQ_MN(MAX_N_INTF_A)  ,A_INTF_FREQ_SC(MAX_N_INTF_A)     &
     & ,A_INTF_END_HR(MAX_N_INTF_A)   ,LEN_INTFA_P(MAX_N_INTF_A)        &
     & ,LEN_INTFA_U(MAX_N_INTF_A)     ,LEN_INTFA_DATA(MAX_N_INTF_A)     &
     & ,LNEWBND(MAX_N_INTF_A)         ,INTF_VERT_INTERP(MAX_N_INTF_A)   &
     & ,INTF_PACK(MAX_N_INTF_A)       ,LBC_STREAM_A(MAX_N_INTF_A)       &
     & ,LBC_UNIT_NO_A(MAX_N_INTF_A)   ,LBC_FIRST_R_RHO(MAX_N_INTF_A)    &
     & ,LBC_BL_LEVELS(MAX_N_INTF_A)   ,LBC_Z_TOP_MODEL(MAX_N_INTF_A)    &
     & ,INTF_VERTLEVS(MAX_N_INTF_A)   ,LBC_Q_MIN                        &
     & ,INTF_L_VAR_LBC                ,INTF_HORZGRID                    &
     & ,LAMBDA_INTF_P                 ,LAMBDA_INTF_U                    &
     & ,PHI_INTF_P                    ,PHI_INTF_V
!---------------------------------------------------------------------
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

! Subroutine arguments

      Integer :: lbc_size

!     Source Grid
      Integer :: src_row_len
      Integer :: src_rows
      Real    :: src_delta_lat
      Real    :: src_delta_long
      Real    :: src_first_lat
      Real    :: src_first_long
      Real    :: src_pole_lat
      Real    :: src_pole_long
      Logical :: src_cyclic
      Logical :: src_rotated

!     LBC Grid
      Integer :: lbc_row_len
      Integer :: lbc_rows
      Integer :: lbc_halo_x
      Integer :: lbc_halo_y
      Real    :: lbc_delta_lat
      Real    :: lbc_delta_long
      Real    :: lbc_first_lat
      Real    :: lbc_first_long
      Real    :: lbc_pole_lat
      Real    :: lbc_pole_long
      
      Logical l_var_lbc
! Input VarRes grid info in degrees 
      Real    :: Lambda_in ( 1-lbc_halo_x: lbc_row_len+lbc_halo_x )  
      Real    :: Phi_in ( 1-lbc_halo_y: lbc_rows+lbc_halo_y )
      
      Integer :: RimWidth

      Integer, dimension (lbc_size) :: lbc_index_bl
      Integer, dimension (lbc_size) :: lbc_index_br
      Real,    dimension (lbc_size) :: lbc_weights_tr
      Real,    dimension (lbc_size) :: lbc_weights_br
      Real,    dimension (lbc_size) :: lbc_weights_tl
      Real,    dimension (lbc_size) :: lbc_weights_bl
      Real,    dimension (lbc_size) :: lbc_coeff1
      Real,    dimension (lbc_size) :: lbc_coeff2

      Integer :: i_uv   !  If 1 => Wind field (u or v). Otherwise 0.

! Local parameters:

      Character (Len=*), Parameter :: RoutineName = 'LBC_Interp_Coeffs'

! Local scalars:

      Integer :: ipt         ! LBC point number
      Integer :: row,pt      ! Loop indices for row and point
      Integer :: iside       ! Loop index for LBC sides
      Integer :: lbc_len     ! Computed no of lbc points
      Integer :: ErrorStatus ! Error Code

      Character (Len=80) :: CMessage

! Local dynamic arrays:

      Real, dimension (:), allocatable :: lambda_lbc
      Real, dimension (:), allocatable :: phi_lbc
      Real, dimension (:), allocatable :: lambda_targ
      Real, dimension (:), allocatable :: phi_targ
      Real, dimension (:), allocatable :: lambda_source
      Real, dimension (:), allocatable :: phi_source
!
! Function & Subroutine calls:
      External EqToLL, H_Int_Co, LLToEq, W_Coeff
!
!- End of header

      ErrorStatus = 0
      CMessage = ' '

! --------------------
! Allocate work arrays
! --------------------

      allocate ( lambda_targ (lbc_size) )
      allocate ( phi_targ    (lbc_size) )
      allocate ( lambda_lbc  (lbc_size) )
      allocate ( phi_lbc     (lbc_size) )

! ---------------------------
! Get lat/Longs of LBC points
! ---------------------------
      If(l_var_lbc) Then 
        
        ipt=0  
   
        Do ISide = 1, 4 
 
          If (ISide == PSouth) Then    !  Southern Boundary (SP)  
 
            Do Row = 1 - lbc_halo_y, RimWidth + i_uv  
              Do Pt = 1 - lbc_halo_x, lbc_row_len + lbc_halo_x   
   
                ipt = ipt + 1  
                Lambda_lbc(ipt) = Lambda_in(pt) 
                Phi_lbc(ipt)    = Phi_in(row)  

              End Do  
            End Do 
  
          End If  !  South  
  
          If (ISide == PEast) Then    !  Eastern Boundary   
  
            Do Row = RimWidth + 1, lbc_rows - RimWidth  
              Do Pt = lbc_row_len - i_uv - RimWidth + 1,                &
     &              lbc_row_len + lbc_halo_x   
  
                ipt = ipt + 1  
                Lambda_lbc(ipt) = Lambda_in(pt) 
                Phi_lbc(ipt)    = Phi_in(row)      

              End Do           
            End Do             
   
          End If  ! East   
 
          If (ISide == PNorth) Then   !  Northern Boundary (NP)  
            Do Row = lbc_rows - RimWidth + 1 - i_uv,                    &
     &             lbc_rows + lbc_halo_y  
              Do Pt = 1-lbc_halo_x, lbc_row_len+lbc_halo_x  
  
              ipt = ipt + 1 
                Lambda_lbc(ipt) = Lambda_in(pt) 
                Phi_lbc(ipt)    = phi_in(row)             

              End Do              
            End Do          
  
          End If  ! North 
  
          If (ISide == PWest) Then   !  Western Boundary  
  
            Do Row = RimWidth + 1, lbc_rows - RimWidth 
              Do Pt = 1 - lbc_halo_x, RimWidth + i_uv 
    
                ipt = ipt + 1  
                Lambda_lbc(ipt) = Lambda_in(pt) 
                Phi_lbc(ipt)    = Phi_in(row) 

              End Do 
            End Do 
  
          End If  ! West 
  
        End Do  !  ISide 

      Else     ! regular grid
      
      ipt=0

      Do ISide = 1, 4

        If (ISide == PSouth) Then    !  Southern Boundary (SP)

          Do Row = 1 - lbc_halo_y, RimWidth + i_uv
            Do Pt = 1 - lbc_halo_x, lbc_row_len + lbc_halo_x

            ipt = ipt + 1
            Lambda_lbc(ipt) = lbc_first_long + (pt-1)  * lbc_delta_long
            Phi_lbc(ipt)    = lbc_first_lat  + (row-1) * lbc_delta_lat

            End Do
          End Do

        End If  !  South

        If (ISide == PEast) Then    !  Eastern Boundary

          Do Row = RimWidth + 1, lbc_rows - RimWidth
            Do Pt = lbc_row_len - i_uv - RimWidth + 1,                  &
     &              lbc_row_len + lbc_halo_x

            ipt = ipt + 1
            Lambda_lbc(ipt) = lbc_first_long + (pt-1)  * lbc_delta_long
            Phi_lbc(ipt)    = lbc_first_lat  + (row-1) * lbc_delta_lat

            End Do
          End Do

        End If  ! East

        If (ISide == PNorth) Then   !  Northern Boundary (NP)
          Do Row = lbc_rows - RimWidth + 1 - i_uv,                      &
     &             lbc_rows + lbc_halo_y
            Do Pt = 1-lbc_halo_x, lbc_row_len+lbc_halo_x

            ipt = ipt + 1
            Lambda_lbc(ipt) = lbc_first_long + (pt-1)  * lbc_delta_long
            Phi_lbc (ipt)   = lbc_first_lat  + (row-1) * lbc_delta_lat

            End Do
          End Do

        End If  ! North

        If (ISide == PWest) Then   !  Western Boundary

          Do Row = RimWidth + 1, lbc_rows - RimWidth
            Do Pt = 1 - lbc_halo_x, RimWidth + i_uv

            ipt = ipt + 1
            Lambda_lbc(ipt) = lbc_first_long + (pt-1)  * lbc_delta_long
            Phi_lbc(ipt)    = lbc_first_lat  + (row-1) * lbc_delta_lat

            End Do
          End Do

        End If  ! West

      End Do  !  ISide
      
      End If  ! variable resolution
      
      lbc_len = ipt

! -----------------------------------------------------
! Check no of lbc points computed ; must match lbc_size
! -----------------------------------------------------

      If (lbc_len /= lbc_size) Then
        write (6,*) ' Mismatch in number of LBC points.'
        write (6,*) ' Expected no of lbc points ',lbc_size
        write (6,*) ' Computed no of lbc points ',lbc_len
        write (CMessage,*) ' Mismatch in number of LBC points.'
        ErrorStatus = 10
! DEPENDS ON: ereport
        Call Ereport ( RoutineName, ErrorStatus, CMessage)
      End If

! -------------------------------------
! Get true lat and longs for lbc points
! -------------------------------------

! DEPENDS ON: eqtoll
      Call EqToLL (                                                     &
     &   phi_lbc                                                        &
     &,  lambda_lbc                                                     &
     &,  phi_targ                                                       &
     &,  lambda_targ                                                    &
     &,  lbc_pole_lat                                                   &
     &,  lbc_pole_long                                                  &
     &,  lbc_len                                                        &
     & )

! ----------------------------------------------
! Calculate the coefficients to rotate the winds
! ----------------------------------------------

      If (i_uv == 1) Then

! DEPENDS ON: w_coeff
        Call W_Coeff (                                                  &
     &       lbc_coeff1                                                 &
     &,      lbc_coeff2                                                 &
     &,      lambda_targ                                                &
     &,      lambda_lbc                                                 &
     &,      lbc_pole_lat                                               &
     &,      lbc_pole_long                                              &
     &,      lbc_len                                                    &
     & )

      End If

! -----------------------------------------------------------
! For rotated model grids, get lat/longs w.r.t the model grid
! -----------------------------------------------------------

      If (src_rotated) Then

! DEPENDS ON: lltoeq
        Call LLToEq (                                                   &
     &       phi_targ                                                   &
     &,      lambda_targ                                                &
     &,      phi_targ                                                   &
     &,      lambda_targ                                                &
     &,      src_pole_lat                                               &
     &,      src_pole_long                                              &
     &,      lbc_len                                                    &
     & )

      End If

! -----------------------------------
! Calculate lat/longs for source grid
! -----------------------------------

      allocate ( lambda_source(src_row_len) )
      allocate ( phi_source   (src_rows)    )

      Do pt=1,src_row_len
        Lambda_Source(pt) = src_first_long + src_delta_long * (pt-1)
      Enddo
      Do row=1,src_rows
        Phi_Source(row)   = src_first_lat  + src_delta_lat  * (row-1)
      Enddo

! --------------------------------------------------
! Calculate the horizontal interpolation cofficients
! --------------------------------------------------

! DEPENDS ON: h_int_co
      Call H_Int_Co (                                                   &
     &   lbc_index_bl                                                   &
     &,  lbc_index_br                                                   &
     &,  lbc_weights_tr                                                 &
     &,  lbc_weights_br                                                 &
     &,  lbc_weights_tl                                                 &
     &,  lbc_weights_bl                                                 &
     &,  Lambda_Source                                                  &
     &,  Phi_Source                                                     &
     &,  Lambda_Targ                                                    &
     &,  Phi_Targ                                                       &
     &,  src_row_len                                                    &
     &,  src_rows                                                       &
     &,  lbc_len                                                        &
     &,  src_cyclic                                                     &
     & )

      deallocate (phi_lbc)
      deallocate (lambda_lbc)
      deallocate (Lambda_Source)
      deallocate (Phi_Source)
      deallocate (Lambda_Targ)
      deallocate (Phi_Targ)

      Return
      END SUBROUTINE LBC_Interp_Coeffs
