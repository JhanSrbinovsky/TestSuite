
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Diagnostic RHcrit Calculation Scheme.
! Subroutine Interface:
      SUBROUTINE LS_CALC_RHCRIT(                                        &
!      Pressure related fields
     &  p_layer_centres                                                 &
!      Array dimensions
     &, LEVELS, row_length, rows, BL_LEVELS, global_row_length          &
!      Prognostic Fields
     &, T, Q, QCF, LAND, LAND_FRAC, ICE_FRAC                            &
!      Logical control
     &, l_mixing_ratio                                                  &
!      Output
     &, RHCPT)
!
      IMPLICIT NONE
!
!     Purpose: To calculate the critical relative humidity in every
!              grid-box.
!
!     Method : The critical relative humidity of a certain grid-box is
!            determined from the variance in a 3*3 set of boxes centred
!            on it. A fit, dependent on pressure, relates the variance
!            of the 3*3 region to the variance within the one grid-box.
!            This variance is converted to a critical relative humidity
!            in a straightforward fashion.
!            Some points in the 3*3 region may be excluded from the
!            variance calculation in the BL layers, if their
!            surfaces do not 'match'. The criterion for matching is that
!            land and sea-ice match, but that open ocean does not match
!            with either of these.
!            In all layers, points in the 3*3 region which lie outside 2
!            std devs of the mean are excluded in a second iteration of
!            the main calculation.
!            The best estimate of the standard deviation in a sample
!            size of n elements in which the mean is calculated from the
!            sample uses (n-1) in the denominator. This subroutine uses
!            n in the denominator because errors in the estimate of the
!            std dev from using such small sample sizes gave rise to
!            concerns over the std dev being too large occasionally,
!            and using n in the denominator rather than (n-1) addresses
!            this issue indirectly. Addressing the issue of large
!            estimates of std dev was considered more important than
!            using the 'best' estimate of the std dev.
!
! Current Owner of Code: S. Cusack / OWNER OF LS CLOUD SECTION
!
! History:
! Version   Date     Comment
!  5.1    09/03/00   Rewritten for New Dynamics.  A.C. Bushell
!                    Based on original code (part of HadAM4 package) at
!                    VN4.5 : subroutine RHCRIT_CALC.   S. Cusack
!
!  5.2    09/09/00   Changes needed to implement critical relative
!                    humidity scheme (fix polar row values).
!                    A.C. Bushell.
!  5.3    19/10/01   Use appropriate gcg routines.   S. Cusack
!
!  5.4    04/09/02   Add gcg_rvecsumf to External list. A.C. Bushell
!
!  6.1    16/08/04   New argument to subroutine added. Changes to the
!                    outlier rejection calculation.         S. Cusack
!  6.4    14/08/06   Use mixing ratio formulation.  Damian Wilson
! Code Description:
!   Language: FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
! Declarations:
!
!  Global Variables:----------------------------------------------------
!
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------
!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

!*----------------------------------------------------------------------
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
!
!  Subroutine arguments
!-----------------------------------------------------------------------
! IN variables
!-----------------------------------------------------------------------
      INTEGER                                                           &
                        !, INTENT(IN)
     &  LEVELS                                                          &
!       No. of levels being processed.
     &, BL_LEVELS                                                       &
!       No. of BL levels being processed.
     &, row_length, rows                                                &
!       Horizontal dimensions of arrays being processed by scheme.
     &, global_row_length
!       Length of full global row (ie. counting all relevant PEs).
!
      LOGICAL                                                           &
                        !, INTENT(IN)
     &  LAND(0:row_length+1,0:rows+1)                                   &
                                          ! The model land mask
     &, l_mixing_ratio                    ! Use mixing ratios
!
      REAL                                                              &
                        !, INTENT(IN)
     &  ICE_FRAC(0:row_length+1,0:rows+1)                               &
!       Ice fraction.
     &, LAND_FRAC(0:row_length+1,0:rows+1)                              &
     &, p_layer_centres(0:row_length+1,0:rows+1,0:LEVELS)               &
!       pressure at all points, on theta levels (Pa).
!       NB: Zero level of array is surface pressure.
     &, Q(0:row_length+1,0:rows+1,LEVELS)                               &
!       Total water content, Q+QCL (kg per kg)
     &, QCF(0:row_length+1,0:rows+1,LEVELS)                             &
!       Cloud ice content at processed levels (kg water per kg air).
     &, T(0:row_length+1,0:rows+1,LEVELS)
!       Liquid/frozen water temperature (K)
!-----------------------------------------------------------------------
! OUT variables
!-----------------------------------------------------------------------
      REAL                                                              &
                        !, INTENT(OUT)
     &  RHCPT(row_length,rows,LEVELS)
!       Critical relative humidity at every gridpoint.
!
!  Local Parameters and other physical constants------------------------
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
! RHCCON2B start
!     Comdeck for use with the RHcrit parametrization of the large-scale
!     cloud scheme, A09_2B.
!     Four constants are used to specify the variable (a function of
!     pressure) which relates the variability of the saturation variable
!     in one box to the variability over 9 climate grid boxes. The
!     variability of the saturation variable in one box is required to
!     specify RHcrit.
!     Note that the constants
!       RHC_CON1=0.522, RHC_CON2=0.122, RHC_CON3=2.5E3, RHC_CON4=1.75E4
!     are only suitable for use in the 2.5*3.75 degrees climate model:
!     these constants depend upon the size of a grid-box.
!
!     The fit is of the form:
!          A=RHC_CON1+RHC_CON2*(p-RHC_CON4)/(RHC_CON3+abs(p-RHC_CON4))
!       where p is the pressure at the layer midpoint.
!     Then,
!           sigma(s) = A * sigma(s,9)
!        where sigma(s) is the std dev of the saturation variable s in
!     one grid-box, and sigma(s,9) is the std dev over 9 boxes.
!
!     RHC_MIN and RHC_MAX are user defined limits on the values of RHc.
!
!     S. Cusack   02-09-98
!
      REAL, PARAMETER :: RHC_CON1 = 0.522
      REAL, PARAMETER :: RHC_CON2 = 0.122
      REAL, PARAMETER :: RHC_CON3 = 2.5E3
      REAL, PARAMETER :: RHC_CON4 = 1.75E4
      REAL, PARAMETER :: RHC_MIN = 0.3
      REAL, PARAMETER :: RHC_MAX = 0.98

! RHCCON2B end
!
      REAL INV9, LS, LSRCP, ERCPR
      PARAMETER ( INV9 = 1./9.                                          &
     &          , LS = LC+LF                                            &
     &          , LSRCP = (LC+LF)/CP                                    &
     &          , ERCPR = EPSILON/(CP*R))
!
!  Local scalars--------------------------------------------------------
!
      INTEGER                                                           &
     &    I, J, K, J8                                                   &
                            ! Simple loop variables
     & ,  COUNT                                                         &
                            ! Counter of points
     & ,  IM1,IP1, JM1,JP1                                              &
                            ! (I-1),(I+1),(J-1),(J+1)
     & ,  ij_field                                                      &
                            ! Number of non-halo points in arrays
     & ,  i_length                                                      &
                            ! Row length for polar row adjustments
     & ,  i_start                                                       &
                            ! Row start point for polar row adjustments
     & ,  ISTAT             ! Status (error code) indicator
!
      REAL                                                              &
     &    MEAN_SUPSAT                                                   &
                            ! MEAN RH OF 3*3 REGION
     & ,  TOT_VAR                                                       &
                            ! TOTAL VARIANCE OF 3*3 REGION
     & ,  SUPSAT_SD_1                                                   & 
                            ! STANDARD DEVIATION OF 'R.H.' IN GRID-BOX
     & ,  SUPSAT_SD_3                                                   & 
                            ! RESOLVED STD DEV OF 'R.H.' IN 3*3 REGION
     & ,  ROOT_6                                                        &
                            ! =sqrt(6.)
     & ,  LATHT                                                         &
                            ! =Lc if T>Tm, ELSE = Lc+Lf
     & ,  TWO_SIGMA                                                     &
                          ! Two times sigma
     & ,  SUPSAT1                                                       &
                            ! Temporary variable
     & ,  SUPSAT2                                                       &
                            ! Temporary variable
     & ,  SUPSAT3                                                       &
                            ! Temporary variable
     & ,  SUPSAT4                                                       &
                            ! Temporary variable
     & ,  SUPSAT5                                                       &
                            ! Temporary variable
     & ,  SUPSAT6                                                       &
                            ! Temporary variable
     & ,  SUPSAT7                                                       &
                            ! Temporary variable
     & ,  SUPSAT8                                                       &
                            ! Temporary variable
     & ,  r_row_length      ! Reciprocal of number of points on FI row)
!
!  Local dynamic arrays-------------------------------------------------
      INTEGER                                                           &
     &   ICOUNT(row_length,rows)       ! Counter of points
!
      LOGICAL                                                           &
     &   OCEAN(0:row_length+1,0:rows+1)! Those points which are not
!                     land, and have a sea-ice fraction less than 0.25.
!
      REAL                                                              &
     &   TL(0:row_length+1,0:rows+1,LEVELS)                             &
                                            ! Conserved temperature
!                                                      (P292.1, UMDP29)
     & , QT(0:row_length+1,0:rows+1,LEVELS)                             &
                                            ! Conserved WATER
!                                                      (P292.2, UMDP29)
     & , QST(0:row_length+1,0:rows+1)                                   &
                                        ! SATURATION VAPOUR PRESSURE
     & , P_GRAD(row_length,rows,LEVELS)                                 &
                                        ! TERM WHICH RELATES
!                                         RH_SD_3 TO RH_SD_1
     & , SUPSAT(0:row_length+1,0:rows+1)                                & 
                                        ! 'RELATIVE HUMIDITY' OF GRIDBOX
     & , AL(0:row_length+1,0:rows+1)                                    &
                                        ! Defined in P292.6 in UMDP 29
     & , SURF_MULT(row_length,rows,8)                                   &
                                        ! Multiplier to take into
!                                         account surface matching.
     & , RHCPT_MEAN(LEVELS)             ! Mean of first interior row.
!
!  External subroutine calls: ------------------------------------------
      EXTERNAL QSAT
      EXTERNAL gcg_rvecsumr, gcg_rvecsumf
!- End of Header
!
! ==Main Block==--------------------------------------------------------
      ROOT_6=SQRT(6.)
!
! Levels_do1:
      DO K=1,LEVELS
! Rows_do1:
        DO J=1,rows
! Rowlen_do1:
          DO I=1,row_length
            P_GRAD(I,J,K) = p_layer_centres(I,J,K) - RHC_CON4
            P_GRAD(I,J,K) = RHC_CON1 +                                  &
     &     (RHC_CON2 * P_GRAD(I,J,K) / (RHC_CON3 + ABS(P_GRAD(I,J,K))) )
!
          END DO ! Rowlen_do1
        END DO ! Rows_do1
      END DO ! Levels_do1
!
! Ocean points defined now as not land and where ice fraction LT 0.25
! Rows_do2:
      DO J=0,rows+1
! Rowlen_do2:
        DO I=0,row_length+1
          OCEAN(I,J) = (LAND_FRAC(I,J) <  0.5) .AND.                    &
     &                                    (ICE_FRAC(I,J) <  2.5E-1)
        END DO ! Rowlen_do2
      END DO ! Rows_do2
!
! A real no. is now assigned to every neighbouring point of every point
! on the grid, if their surfaces match it has the value one, else it is
! zero.
! Eight_do1:
      DO J8=1,8
! Rows_do3:
        DO J=1,rows
! Rowlen_do3:
          DO I=1,row_length
            SURF_MULT(I,J,J8)=0.
          END DO ! Rowlen_do3
        END DO ! Rows_do3
      END DO ! Eight_do1
!
! Rows_do4:
      DO J=1,rows
        JM1 = J - 1
        JP1 = J + 1
! Rowlen_do4:
        DO I=1,row_length
          ICOUNT(I,J)=1
          IM1 = I - 1
          IP1 = I + 1
!
          IF ( (OCEAN(I,J) .AND. OCEAN(IM1,JM1))  .OR.                  &
     &         (.NOT.OCEAN(I,J) .AND. .NOT.OCEAN(IM1,JM1)) ) THEN
            SURF_MULT(I,J,1) = 1.
            ICOUNT(I,J) = ICOUNT(I,J) + 1
          ENDIF
!
          IF ( (OCEAN(I,J) .AND. OCEAN(I,JM1))  .OR.                    &
     &         (.NOT.OCEAN(I,J) .AND. .NOT.OCEAN(I,JM1)) ) THEN
            SURF_MULT(I,J,2) = 1.
            ICOUNT(I,J) = ICOUNT(I,J) + 1
          ENDIF
!
          IF ( (OCEAN(I,J) .AND. OCEAN(IP1,JM1))  .OR.                  &
     &         (.NOT.OCEAN(I,J) .AND. .NOT.OCEAN(IP1,JM1)) ) THEN
            SURF_MULT(I,J,3) = 1.
            ICOUNT(I,J) = ICOUNT(I,J) + 1
          ENDIF
!
          IF ( (OCEAN(I,J) .AND. OCEAN(IM1,J))  .OR.                    &
     &         (.NOT.OCEAN(I,J) .AND. .NOT.OCEAN(IM1,J)) ) THEN
            SURF_MULT(I,J,4) = 1.
            ICOUNT(I,J) = ICOUNT(I,J) + 1
          ENDIF
!
          IF ( (OCEAN(I,J) .AND. OCEAN(IP1,J))  .OR.                    &
     &         (.NOT.OCEAN(I,J) .AND. .NOT.OCEAN(IP1,J)) ) THEN
            SURF_MULT(I,J,5) = 1.
            ICOUNT(I,J) = ICOUNT(I,J) + 1
          ENDIF
!
          IF ( (OCEAN(I,J) .AND. OCEAN(IM1,JP1))  .OR.                  &
     &         (.NOT.OCEAN(I,J) .AND. .NOT.OCEAN(IM1,JP1)) ) THEN
            SURF_MULT(I,J,6) = 1.
            ICOUNT(I,J) = ICOUNT(I,J) + 1
          ENDIF
!
          IF ( (OCEAN(I,J) .AND. OCEAN(I,JP1))  .OR.                    &
     &         (.NOT.OCEAN(I,J) .AND. .NOT.OCEAN(I,JP1)) ) THEN
            SURF_MULT(I,J,7) = 1.
            ICOUNT(I,J) = ICOUNT(I,J) + 1
          ENDIF
!
          IF ( (OCEAN(I,J) .AND. OCEAN(IP1,JP1))  .OR.                  &
     &         (.NOT.OCEAN(I,J) .AND. .NOT.OCEAN(IP1,JP1)) ) THEN
            SURF_MULT(I,J,8) = 1.
            ICOUNT(I,J) = ICOUNT(I,J) + 1
          ENDIF
!
        END DO ! Rowlen_do4
      END DO ! Rows_do4
!
! An initial sweep is done for all grid-cells, obtaining an initial
! estimate of the variance of the 3*3 grid.
!
! Levels_do5:
      DO K=1,BL_LEVELS
! Rows_do5a:
        DO J=0,rows+1
! Rowlen_do5a:
          DO I=0,row_length+1
!   Calculate Tl and QT as in P292.1, P292.2 in UMDP 29.
!  (Assumes version 3A onwards of Section 4)
            TL(I,J,K) = T(I,J,K) - LSRCP*QCF(I,J,K)
            QT(I,J,K) = Q(I,J,K) + QCF(I,J,K)
!
          END DO ! Rowlen_do5a
        END DO ! Rows_do5a
!
! DEPENDS ON: qsat_mix
        call qsat_mix(qst(0,0),tl(0,0,K),p_layer_centres(0,0,K),        &
     &            (row_length+2)*(rows+2),l_mixing_ratio)
! Rows_do5b:
        DO J=0,rows+1
! Rowlen_do5b:
          DO I=0,row_length+1
            IF (TL(I,J,K)  >   TM) THEN
              LATHT = LC/TL(I,J,K)
            ELSE
              LATHT = LS/TL(I,J,K)
            ENDIF
            AL(I,J) = 1./ (1.+LATHT*LATHT*ERCPR*QST(I,J))
! SUPSAT given by P292.3 of UMDP 29.
            SUPSAT(I,J) = AL(I,J)*(QT(I,J,K) - QST(I,J))
          END DO ! Rowlen_do5b
        END DO ! Rows_do5b
!
! Rows_do5c:
        DO J=1,rows
          JM1 = J - 1
          JP1 = J + 1
! Rowlen_do5c:
          DO I=1,row_length
            IM1 = I - 1
            IP1 = I + 1
            SUPSAT1 = SURF_MULT(I,J,1) * SUPSAT(IM1,JM1)
            SUPSAT2 = SURF_MULT(I,J,2) * SUPSAT(I,JM1)
            SUPSAT3 = SURF_MULT(I,J,3) * SUPSAT(IP1,JM1)
            SUPSAT4 = SURF_MULT(I,J,4) * SUPSAT(IM1,J)
            SUPSAT5 = SURF_MULT(I,J,5) * SUPSAT(IP1,J)
            SUPSAT6 = SURF_MULT(I,J,6) * SUPSAT(IM1,JP1)
            SUPSAT7 = SURF_MULT(I,J,7) * SUPSAT(I,JP1)
            SUPSAT8 = SURF_MULT(I,J,8) * SUPSAT(IP1,JP1)
            MEAN_SUPSAT=(SUPSAT1 + SUPSAT2 + SUPSAT3 + SUPSAT4          &
     &                 + SUPSAT5 + SUPSAT6 + SUPSAT7 + SUPSAT8          &
     &                 + SUPSAT(I,J)) / ICOUNT(I,J)
            TOT_VAR = SUPSAT1*SUPSAT1 + SUPSAT2*SUPSAT2                 &
     &              + SUPSAT3*SUPSAT3 + SUPSAT4*SUPSAT4                 &
     &              + SUPSAT5*SUPSAT5 + SUPSAT6*SUPSAT6                 &
     &              + SUPSAT7*SUPSAT7 + SUPSAT8*SUPSAT8                 &
     &              + SUPSAT(I,J)*SUPSAT(I,J)                           &
     &              - ICOUNT(I,J)*MEAN_SUPSAT*MEAN_SUPSAT
            TOT_VAR = ABS(TOT_VAR)
!
!  Now remove the statistical outliers from the 3*3 region, so that
!  sigma, and hence RHcrit, is not biased by extreme values.
!  Points outside 2*sigma of the mean are considered outliers and are
!  rejected.
            IF (ICOUNT(I,J)  >   1) THEN
              TWO_SIGMA = 2. * SQRT(TOT_VAR/ICOUNT(I,J))
            ELSE
              TWO_SIGMA = QST(I,J) * 0.01
            ENDIF
            COUNT=1
!
            IF (ABS(SUPSAT(IM1,JM1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
              SUPSAT1 = 0.
            ELSE IF (SURF_MULT(I,J,1)  >   0.5) THEN
              COUNT = COUNT + 1
            ENDIF
!
            IF (ABS(SUPSAT(I,JM1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
              SUPSAT2 = 0.
            ELSE IF (SURF_MULT(I,J,2)  >   0.5) THEN
              COUNT = COUNT + 1
            ENDIF
!
            IF (ABS(SUPSAT(IP1,JM1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
              SUPSAT3 = 0.
            ELSE IF (SURF_MULT(I,J,3)  >   0.5) THEN
              COUNT = COUNT + 1
            ENDIF
!
            IF (ABS(SUPSAT(IM1,J)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
              SUPSAT4 = 0.
            ELSE IF (SURF_MULT(I,J,4)  >   0.5) THEN
              COUNT = COUNT + 1
            ENDIF
!
            IF (ABS(SUPSAT(IP1,J)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
              SUPSAT5 = 0.
            ELSE IF (SURF_MULT(I,J,5)  >   0.5) THEN
              COUNT = COUNT + 1
            ENDIF
!
            IF (ABS(SUPSAT(IM1,JP1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
              SUPSAT6 = 0.
            ELSE IF (SURF_MULT(I,J,6)  >   0.5) THEN
              COUNT = COUNT + 1
            ENDIF
!
            IF (ABS(SUPSAT(I,JP1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
              SUPSAT7 = 0.
            ELSE IF (SURF_MULT(I,J,7)  >   0.5) THEN
              COUNT = COUNT + 1
            ENDIF
!
            IF (ABS(SUPSAT(IP1,JP1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
              SUPSAT8 = 0.
            ELSE IF (SURF_MULT(I,J,8)  >   0.5) THEN
              COUNT = COUNT + 1
            ENDIF
!
            IF (COUNT >  1) THEN
              MEAN_SUPSAT=(SUPSAT1 + SUPSAT2 + SUPSAT3 + SUPSAT4        &
     &                   + SUPSAT5 + SUPSAT6 + SUPSAT7 + SUPSAT8        &
     &                   + SUPSAT(I,J)) / COUNT
              TOT_VAR = SUPSAT1*SUPSAT1 + SUPSAT2*SUPSAT2               &
     &                + SUPSAT3*SUPSAT3 + SUPSAT4*SUPSAT4               &
     &                + SUPSAT5*SUPSAT5 + SUPSAT6*SUPSAT6               &
     &                + SUPSAT7*SUPSAT7 + SUPSAT8*SUPSAT8               &
     &                + SUPSAT(I,J)*SUPSAT(I,J)                         &
     &                - COUNT*MEAN_SUPSAT*MEAN_SUPSAT
              TOT_VAR = ABS(TOT_VAR)
              SUPSAT_SD_3 = SQRT(TOT_VAR/COUNT)
            ELSE
!           Limit the 3*3 grid variance when scatter is large.
              SUPSAT_SD_3 = QST(I,J)*0.01
            ENDIF
!
! Try to detect if the central point (i,j) is an outlier, as can happen
! in a GPS situation. If so, set the variance to a small value.
            IF (COUNT >  2) THEN
              MEAN_SUPSAT=(SUPSAT1 + SUPSAT2 + SUPSAT3 + SUPSAT4        &
     &         + SUPSAT5 + SUPSAT6 + SUPSAT7 + SUPSAT8) / (COUNT-1.0)
              TOT_VAR = SUPSAT1*SUPSAT1 + SUPSAT2*SUPSAT2               &
     &                + SUPSAT3*SUPSAT3 + SUPSAT4*SUPSAT4               &
     &                + SUPSAT5*SUPSAT5 + SUPSAT6*SUPSAT6               &
     &                + SUPSAT7*SUPSAT7 + SUPSAT8*SUPSAT8               &
     &                - COUNT*MEAN_SUPSAT*MEAN_SUPSAT
              IF (ABS(SUPSAT(I,J)-MEAN_SUPSAT) >                        &
     &                    (2.0*SQRT(ABS(TOT_VAR)/COUNT))) THEN
                SUPSAT_SD_3 = QST(I,J)*0.01
              ENDIF
            ENDIF
!
! P_GRAD determines the relation between 3*3 and sub-grid variance.
            SUPSAT_SD_1 = P_GRAD(I,J,K) * SUPSAT_SD_3
! RHCPT defined from P292.14 in UMDP 29
            RHCPT(I,J,K) = 1. - (ROOT_6*SUPSAT_SD_1 /(AL(I,J)*QST(I,J)))
! RHcrit is now limited to lie between a range defined in RHCCON2B
            RHCPT(I,J,K) = MAX(RHCPT(I,J,K),RHC_MIN)
            RHCPT(I,J,K) = MIN(RHCPT(I,J,K),RHC_MAX)
          END DO ! Rowlen_do5c
        END DO ! Rows_do5c
!
      END DO ! Levels_do5
!
! The same calculations as above are performed, but the 'surface match'
! criterion is now dropped (atmosphere less influenced by surface at
! greater heights).
! Levels_if6:
      IF (LEVELS  >   BL_LEVELS) THEN
!
! An initial sweep is done for all grid-cells, obtaining an initial
! estimate of the variance of the 3*3 grid.
!
! Levels_do6:
        DO K=(BL_LEVELS+1),LEVELS
! Rows_do6a:
          DO J=0,rows+1
! Rowlen_do6a:
            DO I=0,row_length+1
!   Calculate Tl and QT as in P292.1, P292.2 in UMDP 29.
!  (Assumes version 3A onwards of Section 4)
              TL(I,J,K) = T(I,J,K) - LSRCP*QCF(I,J,K)
              QT(I,J,K) = Q(I,J,K) + QCF(I,J,K)
!
            END DO ! Rowlen_do6a
          END DO ! Rows_do6a
!
! DEPENDS ON: qsat_mix
          call qsat_mix(qst(0,0),tl(0,0,K),p_layer_centres(0,0,K),      &
     &              (row_length+2)*(rows+2),l_mixing_ratio)
!
! Rows_do6b:
          DO J=0,rows+1
! Rowlen_do6b:
            DO I=0,row_length+1
              IF (TL(I,J,K)  >   TM) THEN
                LATHT = LC/TL(I,J,K)
              ELSE
                LATHT = LS/TL(I,J,K)
              ENDIF
              AL(I,J) = 1./ (1.+LATHT*LATHT*ERCPR*QST(I,J))
! SUPSAT given by P292.3 of UMDP 29.
              SUPSAT(I,J) = AL(I,J)*(QT(I,J,K) - QST(I,J))
            END DO ! Rowlen_do6b
          END DO ! Rows_do6b
!
! Rows_do6c:
          DO J=1,rows
            JM1 = J - 1
            JP1 = J + 1
! Rowlen_do6c:
            DO I=1,row_length
              IM1 = I - 1
              IP1 = I + 1
              SUPSAT1 = SUPSAT(IM1,JM1)
              SUPSAT2 = SUPSAT(I,JM1)
              SUPSAT3 = SUPSAT(IP1,JM1)
              SUPSAT4 = SUPSAT(IM1,J)
              SUPSAT5 = SUPSAT(IP1,J)
              SUPSAT6 = SUPSAT(IM1,JP1)
              SUPSAT7 = SUPSAT(I,JP1)
              SUPSAT8 = SUPSAT(IP1,JP1)
              MEAN_SUPSAT=(SUPSAT1 + SUPSAT2 + SUPSAT3 + SUPSAT4        &
     &                   + SUPSAT5 + SUPSAT6 + SUPSAT7 + SUPSAT8        &
     &                   + SUPSAT(I,J)) * INV9
              TOT_VAR = SUPSAT1*SUPSAT1 + SUPSAT2*SUPSAT2               &
     &                + SUPSAT3*SUPSAT3 + SUPSAT4*SUPSAT4               &
     &                + SUPSAT5*SUPSAT5 + SUPSAT6*SUPSAT6               &
     &                + SUPSAT7*SUPSAT7 + SUPSAT8*SUPSAT8               &
     &                + SUPSAT(I,J)*SUPSAT(I,J)                         &
     &                - 9. * MEAN_SUPSAT*MEAN_SUPSAT
              TOT_VAR = ABS(TOT_VAR)
!
!  Now remove the statistical outliers from the 3*3 region, so that
!  sigma, and hence RHcrit, is not biased by extreme values.
!  Points outside 2*sigma of the mean are considered outliers and are
!  rejected.
              TWO_SIGMA = 0.67*SQRT(TOT_VAR)  ! =2*SQRT(TOT_VAR/9)
              COUNT=1
!
            IF (ABS(SUPSAT(IM1,JM1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
                SUPSAT1=0.
              ELSE
                COUNT=COUNT+1
              ENDIF
!
            IF (ABS(SUPSAT(I,JM1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
                SUPSAT2=0.
              ELSE
                COUNT=COUNT+1
              ENDIF
!
            IF (ABS(SUPSAT(IP1,JM1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
                SUPSAT3=0.
              ELSE
                COUNT=COUNT+1
              ENDIF
!
            IF (ABS(SUPSAT(IM1,J)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
                SUPSAT4=0.
              ELSE
                COUNT=COUNT+1
              ENDIF
!
            IF (ABS(SUPSAT(IP1,J)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
                SUPSAT5=0.
              ELSE
                COUNT=COUNT+1
              ENDIF
!
            IF (ABS(SUPSAT(IM1,JP1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
                SUPSAT6=0.
              ELSE
                COUNT=COUNT+1
              ENDIF
!
            IF (ABS(SUPSAT(I,JP1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
                SUPSAT7=0.
              ELSE
                COUNT=COUNT+1
              ENDIF
!
            IF (ABS(SUPSAT(IP1,JP1)-MEAN_SUPSAT)  >   TWO_SIGMA) THEN
                SUPSAT8=0.
              ELSE
                COUNT=COUNT+1
              ENDIF
!
              IF (COUNT  >   1) THEN
                MEAN_SUPSAT=(SUPSAT1 + SUPSAT2 + SUPSAT3 + SUPSAT4      &
     &                     + SUPSAT5 + SUPSAT6 + SUPSAT7 + SUPSAT8      &
     &                     + SUPSAT(I,J)) / COUNT
                TOT_VAR = SUPSAT1*SUPSAT1 + SUPSAT2*SUPSAT2             &
     &                  + SUPSAT3*SUPSAT3 + SUPSAT4*SUPSAT4             &
     &                  + SUPSAT5*SUPSAT5 + SUPSAT6*SUPSAT6             &
     &                  + SUPSAT7*SUPSAT7 + SUPSAT8*SUPSAT8             &
     &                  + SUPSAT(I,J)*SUPSAT(I,J)                       &
     &                  - COUNT*MEAN_SUPSAT*MEAN_SUPSAT
                TOT_VAR = ABS(TOT_VAR)
                SUPSAT_SD_3 = SQRT(TOT_VAR/COUNT)
              ELSE
!             Limit the 3*3 grid variance when scatter is large.
                SUPSAT_SD_3 = QST(I,J) * 0.01
              ENDIF
!
! Try to detect if the central point (i,j) is an outlier, as can happen
! in a GPS situation. If so, set the variance to a small value.
            IF (COUNT >  2) THEN
              MEAN_SUPSAT=(SUPSAT1 + SUPSAT2 + SUPSAT3 + SUPSAT4        &
     &         + SUPSAT5 + SUPSAT6 + SUPSAT7 + SUPSAT8) / (COUNT-1.0)
              TOT_VAR = SUPSAT1*SUPSAT1 + SUPSAT2*SUPSAT2               &
     &                + SUPSAT3*SUPSAT3 + SUPSAT4*SUPSAT4               &
     &                + SUPSAT5*SUPSAT5 + SUPSAT6*SUPSAT6               &
     &                + SUPSAT7*SUPSAT7 + SUPSAT8*SUPSAT8               &
     &                - COUNT*MEAN_SUPSAT*MEAN_SUPSAT
              IF (ABS(SUPSAT(I,J)-MEAN_SUPSAT) >                        &
     &                    (2.0*SQRT(ABS(TOT_VAR)/COUNT))) THEN
                SUPSAT_SD_3 = QST(I,J)*0.01
              ENDIF
            ENDIF
!
! P_GRAD determines the relation between 3*3 and sub-grid variance.
              SUPSAT_SD_1 = P_GRAD(I,J,K) * SUPSAT_SD_3
! RHCPT defined from P292.14 in UMDP 29
              RHCPT(I,J,K) = 1.-(ROOT_6*SUPSAT_SD_1 /(AL(I,J)*QST(I,J)))
! RHcrit is now limited to lie between a range defined in RHCCON2B
              RHCPT(I,J,K) = MAX(RHCPT(I,J,K),RHC_MIN)
              RHCPT(I,J,K) = MIN(RHCPT(I,J,K),RHC_MAX)
            END DO ! Rowlen_do6c
          END DO ! Rows_do6c
!
        END DO  ! Levels_do6
      ENDIF  ! Levels_if6
!
! Tidy up at South Pole : Pole is mean of first interior row.
! SouthPole_if1:
      IF (at_extremity(PSouth)) THEN
!
        ij_field = row_length * rows
        i_length = row_length
! Start point of first interior (ie. non-polar) row
        i_start  = row_length + 1
!
        r_row_length = 1. / global_row_length
!       Number of points in sum should be global_row_length: might
!       need to modify r_row_length if i_length lt row_length.
!       r_row_length = r_row_length * row_length / i_length
!
! Sum over points in PEs in order along first interior row
! (gc_proc_row_group is group ID for rows of PEs, here only PSouth).
        CALL gcg_rvecsumr(ij_field, i_length, i_start, LEVELS, RHCPT,   &
     &                    gc_proc_row_group, ISTAT, RHCPT_MEAN)
!
! Levels_do7:
        DO K=1,LEVELS
          RHCPT_MEAN(K) = RHCPT_MEAN(K) * r_row_length
! Rowlen_do7:
          DO I=1,row_length
            RHCPT(I,1,K) = RHCPT_MEAN(K)
          END DO ! Rowlen_do7
        END DO  ! Levels_do7
!
      ENDIF  ! SouthPole_if1
!
! Tidy up at North Pole : Pole is mean of first interior row.
! NorthPole_if1:
      IF (at_extremity(PNorth)) THEN
!
        ij_field = row_length * rows
        i_length = row_length
! Start point of first interior (ie. non-polar) row
        i_start  = (rows - 2) * row_length + 1
!
        r_row_length = 1. / global_row_length
!       Number of points in sum should be global_row_length: might
!       need to modify r_row_length if i_length lt row_length.
!       r_row_length = r_row_length * row_length / i_length
!
! Sum over points in PEs in order along first interior row
! (gc_proc_row_group is group ID for rows of PEs, here only PNorth).
        CALL gcg_rvecsumr(ij_field, i_length, i_start, LEVELS, RHCPT,   &
     &                    gc_proc_row_group, ISTAT, RHCPT_MEAN)
!
! Levels_do8:
        DO K=1,LEVELS
          RHCPT_MEAN(K) = RHCPT_MEAN(K) * r_row_length
! Rowlen_do8:
          DO I=1,row_length
            RHCPT(I,rows,K) = RHCPT_MEAN(K)
          END DO ! Rowlen_do8
        END DO  ! Levels_do8
!
      ENDIF  ! NorthPole_if1
!
      RETURN
      END SUBROUTINE LS_CALC_RHCRIT
