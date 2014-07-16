
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculate simple test diagnostics based on a simple analytic formula
!
! Subroutine Interface:
      SUBROUTINE Testdiag(                                              &
     &  p_field,v_field,rows,n_rows,row_length                          &
     & ,EW_SPACE,NS_SPACE,FIRST_LAT,FIRST_LONG,PHI_POLE,LAMBDA_POLE     &
     & ,ELF                                                             &
     & ,PRESS_LEVELS_LIST,NO_PRESS_LEVELS                               &
     & ,MODEL_LEVELS_LIST,NO_MODEL_LEVELS,FORECAST_HRS                  &
     & ,DIAG1,DIAG2,DIAG3,DIAG4                                         &
     & ,qdia1,qdia2,qdia3,qdia4)
      IMPLICIT NONE
!
! Description:

! Calculate simple test diagnostics based on a simple analytic formula:
!
!    VALUE=A*(LATITUDE+90.)+B*LONGITUDE+C*LEVEL+D*FORECAST_HRS
!    where A=1.0, B=1.0E2, C=1.0E3, D=1.0E4
!    and (LAT,LONG) are in degrees, actual position (rotated for LAM),
!    LEVEL is either the model level (real number) or
!                        pressure level (mb),
!    FORECAST_HRS in T+hours after assimilation time.
!    Theses diagnostics are to be used for checking output procedures
!    for various post-processing routes.
!    Four diagnostics are supported:
!    1. single-level FIELD (LEVEL=0.) at V points of Arakawa C grid
!    2. single-level FIELD (LEVEL=0.) at P points of Arakawa C grid
!    3. multi -level FIELD (LEVEL=press level) at P points of C grid
!    4. multi -level FIELD (LEVEL=model level) at P points of C grid
!
! Method:
!      Calculate lat,long offsets for mpp
!   1. Calculate FIRST DIAGNOSTIC  (V GRID SINGLE LEVEL)
!   2. Calculate ACTUAL LATITUDES, LONGITUDES FOR P FIELDS (DIAG 2-4)
!   3. Calculate SECOND DIAGNOSTIC (P GRID SINGLE LEVEL)
!   4. Calculate THIRD  DIAGNOSTIC (P GRID PRESSURE LEVELS)
!   5. Calculate FOURTH DIAGNOSTIC (P GRID MODEL    LEVELS)
!
!    DOCUMENTATION:  UM Doc Paper D7
!
! Current Code Owner: R Rawlins
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 5.0  29/06/99   New deck based on TESTDI1A, developed to provide
!                 functionality for mpp and upgrade for
!                 C-P/C dynamics grid. R. Rawlins
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Global variables (*CALLed COMDECKs etc...):

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER                                                           &
     &  p_field                                                         &
                         !IN   Horizontal field size p points
     &, v_field                                                         &
                         !IN   Horizontal field size v points
     &, rows                                                            &
                         !IN   No. of rows for p field
     &, n_rows                                                          &
                         !IN   No. of rows for v field
     &, row_length                                                      &
                         !IN   No. of points per row
     &, NO_MODEL_LEVELS                                                 &
                         !IN   model levels for output
     &, NO_PRESS_LEVELS                                                 
                         !IN   press levels for output

! UM6.5 - MODEL_ANALYSIS_HRS chnged to real -
!                   requires FORECAST_HRS changed to REAL also
      REAL FORECAST_HRS     !IN   FORECAST HOURS T+0, etc

      LOGICAL                                                           &
     &  ELF                                                             &
                         !IN  TRUE IF MODEL IS LAM WITH ROTATED GRID
     & ,qdia1                                                           &
                         !IN  STASHflag for DIAG1
     & ,qdia2                                                           &
                         !IN  STASHflag for DIAG2
     & ,qdia3                                                           &
                         !IN  STASHflag for DIAG3
     & ,qdia4            !IN  STASHflag for DIAG4

      REAL                                                              &
     &  EW_SPACE                                                        &
                         !IN  DELTA LONGITUDE (DEGREES)
     &, NS_SPACE                                                        &
                         !IN  DELTA  LATITUDE (DEGREES)
     &, FIRST_LAT                                                       &
                         !IN  latitude  of first p row (degrees)
     &, FIRST_LONG                                                      &
                         !IN  longitude of first p col (degrees)
     &, PHI_POLE                                                        &
                         !IN  latitude  of the pseudo pole
     &, LAMBDA_POLE      !IN  longitude of the pseudo pole

!   Array  arguments with intent(in):
      REAL                                                              &
     &  MODEL_LEVELS_LIST(NO_MODEL_LEVELS)                              &
                                           !IN LEVELS list (for DIAG3)
     &, PRESS_LEVELS_LIST(NO_PRESS_LEVELS) !IN LEVELS list (for DIAG4)

!   Scalar arguments with intent(InOut):

!   Array  arguments with intent(InOut):

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
      REAL                                                              &
     &  DIAG1(v_field)                                                  &
                                        !OUT DIAGNOSTIC 1
     &, DIAG2(p_field)                                                  &
                                        !OUT DIAGNOSTIC 2
     &, DIAG3(p_field,NO_PRESS_LEVELS)                                  &
                                        !OUT DIAGNOSTIC 3
     &, DIAG4(p_field,NO_MODEL_LEVELS)  !OUT DIAGNOSTIC 4


! Local parameters:
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
!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

!*----------------------------------------------------------------------
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
      REAL                                                              &
     &  A,B,C,D  ! COEFFICIENTS FOR CALCULATING VALUES OF FIELD
!
      PARAMETER(A=1.0,B=1.0E2,C=1.0E3,D=1.0E4)

! Local scalars:
      INTEGER                                                           &
     &  I,J,K                                                           &
                     !  LOOP COUNTERS
     & ,L                                                               &
                     !  LOOP INDEX
     & ,OFFSETX                                                         &
                     !  INDEX OFFSETs (for calculating mpp lat,longs)
     & ,OFFSETY

! Local dynamic arrays:
      REAL                                                              &
     &  LATITUDE(p_field)                                               &
                          ! latitude  in degrees
     &, LONGITUDE(p_field)                                              &
                          ! longitude in degrees
     &, LAT(p_field)                                                    &
                          ! latitude  in degrees on equatorial grid
     &, LONG(p_field)     ! longitude in degrees on equatorial grid

! Function & Subroutine calls:
      External                                                          &
     & LLTOEQ

!- End of header


!
!  Calculate lat,long offsets for mpp
!

      OFFSETX=datastart(1)- Offx - 1
      OFFSETY=datastart(2)- Offy - 1

!-----------------------------------------------------------------------
!   1. CALCULATE FIRST DIAGNOSTIC  (V GRID SINGLE LEVEL)
!-----------------------------------------------------------------------
      IF(qdia1) THEN
!
!   1a. FIND EQUATORIAL LATITUDES,LONGITUDES
!
         DO J=1,n_rows
           DO I=1,row_length
             L= I + (J-1)*row_length
             LAT (L)=  FIRST_LAT  - NS_SPACE*(J+OFFSETY-0.5)
             LONG(L)=  FIRST_LONG + EW_SPACE*(I+OFFSETX-0.5)
           ENDDO
         ENDDO
!
!   1b. CONVERT TO ACTUAL LATITUDE,LONGITUDE IF ELF GRID
!
         IF(ELF) THEN
! DEPENDS ON: eqtoll
           CALL EQTOLL(LAT,LONG,LATITUDE,LONGITUDE,PHI_POLE,LAMBDA_POLE,&
     &                 v_field)
         ELSE
           DO I=1,v_field
             LATITUDE(I) =LAT(I)
             LONGITUDE(I)=LONG(I)
           ENDDO
         ENDIF
!
!   1c. CALCULATE VALUE FROM ANALYTIC FUNCTION
!

         DO I=1,v_field
           DIAG1(I)=A*(LATITUDE(I)+90.0) + B*LONGITUDE(I) +             &
     &              D*FORECAST_HRS
         ENDDO

      ENDIF               ! END OF qdia1 TEST

!-----------------------------------------------------------------------
!   2. CALCULATE ACTUAL LATITUDES, LONGITUDES FOR P FIELDS (DIAG 2-4)
!-----------------------------------------------------------------------
      IF(qdia2.OR.qdia3.OR.qdia4) THEN
!
!   2a. FIND EQUATORIAL LATITUDES,LONGITUDES
!
         DO J=1,rows
           DO I=1,row_length
             L= I + (J-1)*row_length
             LAT (L)=  FIRST_LAT  - NS_SPACE*(J+OFFSETY-1)
             LONG(L)=  FIRST_LONG + EW_SPACE*(I+OFFSETX-1)
           ENDDO
         ENDDO
!
!   2b. CONVERT TO ACTUAL LATITUDE,LONGITUDE IF ELF GRID
!
         IF(ELF) THEN
! DEPENDS ON: eqtoll
           CALL EQTOLL(LAT,LONG,LATITUDE,LONGITUDE,PHI_POLE,LAMBDA_POLE,&
     &                 p_field)
         ELSE
           DO I=1,p_field
             LATITUDE(I) =LAT(I)
             LONGITUDE(I)=LONG(I)
           ENDDO
         ENDIF
      ENDIF                      ! END OF qdia2-4 TEST
!-----------------------------------------------------------------------
!   3. CALCULATE SECOND DIAGNOSTIC (P GRID SINGLE LEVEL)
!-----------------------------------------------------------------------
      IF(qdia2) THEN

         DO I=1,p_field
           DIAG2(I)=A*(LATITUDE(I)+90.0) + B*LONGITUDE(I) +             &
     &              D*FORECAST_HRS
         ENDDO

      ENDIF               ! END OF qdia2 TEST
!-----------------------------------------------------------------------
!   4. CALCULATE THIRD  DIAGNOSTIC (P GRID PRESSURE LEVELS)
!-----------------------------------------------------------------------
      IF(qdia3) THEN

         DO K=1,NO_PRESS_LEVELS
           DO I=1,p_field
             DIAG3(I,K)=A*(LATITUDE(I)+90.0)   + B*LONGITUDE(I) +       &
     &                  C*PRESS_LEVELS_LIST(K) + D*FORECAST_HRS
           ENDDO
         ENDDO

      ENDIF               ! END OF qdia3 TEST
!-----------------------------------------------------------------------
!   5. CALCULATE FOURTH DIAGNOSTIC (P GRID MODEL    LEVELS)
!-----------------------------------------------------------------------
      IF(qdia4) THEN

         DO K=1,NO_MODEL_LEVELS
           DO I=1,p_field
             DIAG4(I,K)=A*(LATITUDE(I)+90.0)   + B*LONGITUDE(I) +       &
     &                  C*MODEL_LEVELS_LIST(K) + D*FORECAST_HRS
           ENDDO
         ENDDO

      ENDIF               ! END OF qdia4 TEST

      RETURN
      END SUBROUTINE Testdiag
!=======================================================================
