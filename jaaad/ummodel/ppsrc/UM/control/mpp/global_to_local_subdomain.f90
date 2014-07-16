

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Parallel UM : Transform from global to local co-ordinates:
! GLOBAL_TO_LOCAL_SUBDOMAIN: converts global subdomain boundaries
!                            to local subdomain boundaries
! GLOBAL_TO_LOCAL_RC: converts global row,column co-ordinates to
!                     processor co-ordinates plus local
!                     co-ordinates within the processor.
!
! Subroutine Interface:
      SUBROUTINE GLOBAL_TO_LOCAL_SUBDOMAIN(                             &
     &                              L_include_halosEW,                  &
     &                              L_include_halosNS,                  &
     &                              grid_code,halo_type,procid,         &
     &                              global_start_row_in,                &
     &                              global_end_col_in,                  &
     &                              global_end_row_in,                  &
     &                              global_start_col_in,                &
     &                              local_start_row,local_end_col,      &
     &                              local_end_row,local_start_col)
      IMPLICIT NONE
!
! Description:
! Takes a global definition of a subdomain region (in terms of
! model gridpoints) and translates it into local numbers.
! This effectively means local co-ordinates of the region of the
! subdomain which intersects with this processor's area.
!
! Method:
! Use the datastart variable in PARVARS to see if the requested
! subdomain intersects with this processor's area, if it does
! then use datastart to convert to local co-ordinate and do a bit
! of logic using MAX and MIN to ensure the local co-ordinates
! actually lie within the local area  Then make any corrections
! necessary to account for a subdomain which crosses over the
! 0 longitude line. Finally, if L_include_halos is set to
! .TRUE. - include any relevant halo regions.
!
! Current code owner : Paul Burton
!
! History:
!  Model    Date     Modification history from model version 4.2
!  version
!  4.2      03/09/96 New deck created for mpp code.  P.Burton
!  4.3      13/03/97 Various bug fixes               P.Burton
!  4.4      12/06/97 Another bug fix                 P.Burton
!  5.0      9/6/99   Added halo_type argument
!                    Change references to PARVARS variables to include
!                      the new dimensions
!                    Changed logic to allow for indexing from South
!                                                          P.Burton
!  5.1      28/01/00 Changed variable names to deal with start+end
!                    rather than North,South,East+West which get
!                    confusing with grids which have different
!                    orderings.                           P.Burton
!   5.1     25/01/00 Corrections for uv points on B grid. R.Rawlins
!   5.3     22/11/01 Enable mpp as the only option for
!                    small executables         E.Leung
!  5.3      21/08/01 Added ocean boundary data grid types M. J. Bell
!  5.5      10/08/00 Modification for parallelisation of WAM.
!                    Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!  6.0  17/09/03  Add def for new NEC opt section c96_1c. R Barnes
!
! Subroutine arguments:

      LOGICAL                                                           &
     &  L_include_halosEW                                               &
                           ! IN : include East-West halos in local
!                          !      region if set to .TRUE.
     &, L_include_halosNS  ! IN : include North-South halos in local
!                          !      region if set to .TRUE.
      INTEGER                                                           &

     &  grid_code                                                       &
                         ! IN : STASH grid type of field
     &, halo_type                                                       &
                         ! IN : which type of halo has the field got
     &, procid                                                          &
                         ! IN : processor to produce result for
     &, global_start_row_in                                             &
                            ! IN : first row of global subdomain
     &, global_end_col_in                                               &
                            ! IN : last column of global subdomain
     &, global_end_row_in                                               &
                            ! IN : last row of global subdomain
     &, global_start_col_in                                             &
                            ! IN : first column of global subdomain

     &, local_start_row                                                 &
                            ! OUT : first row of local subdomain
     &, local_end_col                                                   &
                            ! OUT : last column of local subdomain
     &, local_end_row                                                   &
                            ! OUT : last row of local subdomain
     &, local_start_col     ! OUT : first column of local subdomain

! Parameters and Common blocks

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
!LL  Comdeck: STERR ----------------------------------------------------
!LL
!LL  Purpose: PARAMETER names for STASH processing error codes;
!LL           fatal errors have positive codes, warnings negative.
!LL
!LL  Author:   S.Tett
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.3  16/09/93  Add st_illegal_weight error code.
!LL                   Added st_no_data for mpp code
!LL                   (means a processor does not contain any data
!LL                    for a given subdomain)                 P.Burton
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: D70
!LL
!LL  Project task: D7
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage handling and diagnostic
!LL                                 system (STASH)
!
! Warning codes
!
         integer st_upper_less_lower ! warning code for bad domain
         parameter(st_upper_less_lower=-1)

         integer st_not_supported ! warning code for unsupported routine
         parameter(st_not_supported=-2)
         integer st_no_data,st_nd ! indicates no data on a processor
         parameter(st_no_data=-3,st_nd=-3)
!
! Error codes
!
         integer st_bad_array_param ! error code for dodgy array params
         parameter(st_bad_array_param=1)

         integer st_bad_address     ! error code for address violation
         parameter(st_bad_address=2)

         integer st_unknown ! error code for unknown option
         parameter(st_unknown=3)

         integer st_bad_wraparound ! error code for illegal wraparound
         parameter(st_bad_wraparound=4)

         integer st_illegal_weight ! error code for illegal weighting
         parameter(st_illegal_weight=9)

         integer unknown_weight ! error code for an unknown weight
         parameter(unknown_weight=10)

         integer unknown_mask ! error code for an unknown mask
         parameter(unknown_mask=11)

         integer unknown_processing ! error code for unknown processing
         parameter(unknown_processing=12)

         integer nonsense ! error code for general nonsense request
         parameter(nonsense=13)

! Local variables
      INTEGER                                                           &
! Copies of the input arguments, that can be modified for
! wrap-around calculations
     &  global_start_row,global_end_col                                 &
     &, global_end_row,global_start_col                                 &
     &, fld_type                                                        &
                  ! is field on P or U grid?
     &, row_len_nh                                                      &
                      ! row length when halos are removed
     &, nrows_nh                                                        &
                      ! number of rows when halos are removed
     &, first_global_pt_EW                                              &
                           ! global point number of first and last
     &, last_global_pt_EW                                               &
                           ! local points in local area
     &, first_global_pt_NS                                              &
                           ! in the East-West and
     &, last_global_pt_NS  ! North-South directions

      LOGICAL                                                           &
! Logicals indicating if this processor contains part of a
! subdomain
     &  NS_intersect,EW_intersect                                       &
     &, wrap                                                            &
             ! set to .TRUE. if the subdomain passes over the
!            ! the 0 degree longitude line
     &, fullfield ! if the field is NOT a subdomain

      INTEGER GET_FLD_TYPE  ! function

! ------------------------------------------------------------------

! Copy the global_in variables into local variables

      global_start_row=global_start_row_in
      global_end_col=global_end_col_in
      global_end_row=global_end_row_in
      global_start_col=global_start_col_in

! Find out if the data is on a mass or velocity grid

! DEPENDS ON: get_fld_type
      fld_type=GET_FLD_TYPE(grid_code)

      IF (fld_type  ==  fld_type_unknown) THEN
        WRITE(6,*) 'GLOBAL_TO_LOCAL_SUBDOMAIN encountered ',            &
     &    'field with gridtype code ',grid_code
        WRITE(6,*) 'Unable to process this field.'
        local_start_row=st_no_data
        local_end_col=st_no_data
        local_end_row=st_no_data
        local_start_col=st_no_data
        GOTO 9999
      ENDIF

! Set up logical indicating if this is a full field, or just
! a subdomain

      fullfield= ((global_start_col  ==  1) .AND.                       &
     &            (global_end_col  ==  glsize(1,fld_type)) .AND.        &
     &            (global_start_row  ==  1) .AND.                       &
     &            (global_end_row  ==  glsize(2,fld_type)))

! If this is a fullfield (ie. not a subdomain) the local addressing
! is easy:

      IF (fullfield) THEN

        IF (L_include_halosNS) THEN
          local_start_row=1
          local_end_row=g_lasize(2,fld_type,halo_type,procid)
        ELSE
          local_start_row=1+halosize(2,halo_type)
          local_end_row=g_lasize(2,fld_type,halo_type,procid)-          &
     &                  halosize(2,halo_type)
        ENDIF
        IF (L_include_halosEW) THEN
          local_start_col=1
          local_end_col=g_lasize(1,fld_type,halo_type,procid)
        ELSE
          local_start_col=1+halosize(1,halo_type)
          local_end_col=g_lasize(1,fld_type,halo_type,procid)-          &
     &                 halosize(1,halo_type)
        ENDIF

      ELSE ! a subdomain requires some careful analysis:

        row_len_nh=g_blsize(1,fld_type,procid)
        nrows_nh=g_blsize(2,fld_type,procid)

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

        IF (global_end_col  <   global_start_col) THEN
          wrap=.TRUE.
        ELSEIF (global_end_col  >   glsize(1,fld_type)) THEN
          wrap=.TRUE.
          global_end_col=global_end_col-glsize(1,fld_type)
        ELSE
          wrap=.FALSE.
        ENDIF

        EW_intersect =                                                  &
     &    (( .NOT. wrap) .AND.                                          &
     &     ((global_end_col  >=  first_global_pt_EW) .AND.              &
     &      (global_start_col  <=  last_global_pt_EW)))                 &
     &    .OR.                                                          &
     &    ((wrap) .AND.                                                 &
     &     ((global_start_col  <=  last_global_pt_EW) .OR.              &
     &      (global_end_col  >=  first_global_pt_EW)))

        NS_intersect =                                                  &
     &    ((global_start_row  <=  last_global_pt_NS) .AND.              &
     &     (global_end_row  >=  first_global_pt_NS))

        IF (NS_intersect) THEN

          IF ((global_start_row  >=  first_global_pt_NS) .AND.          &
     &        (global_start_row  <=  last_global_pt_NS)) THEN
! This processor contains the NS start of the subarea
            local_start_row=global_start_row-first_global_pt_NS+        &
     &                    halosize(2,halo_type)+1
          ELSE
! This processor is to the North of the start of the subarea
            local_start_row=1+halosize(2,halo_type)
          ENDIF

          IF ((global_end_row  >=  first_global_pt_NS) .AND.            &
     &        (global_end_row  <=  last_global_pt_NS)) THEN
! This processor contains the NS end of the subarea
            local_end_row=global_end_row-first_global_pt_NS+            &
     &                    halosize(2,halo_type)+1
          ELSE
! This processor is to the South of the subarea
            local_end_row=halosize(2,halo_type)+nrows_nh
          ENDIF

        ELSE

          local_start_row=st_no_data
          local_end_row=st_no_data

        ENDIF

        IF (EW_intersect) THEN

          IF ((global_start_col  >=  first_global_pt_EW) .AND.          &
     &        (global_start_col  <=  last_global_pt_EW)) THEN
! This processor contains the EW start of the subarea
            local_start_col=global_start_col-first_global_pt_EW+        &
     &                   halosize(1,halo_type)+1
          ELSE
! This processor is to the right of the start of the subarea
            local_start_col=1+halosize(1,halo_type)
          ENDIF

          IF ((global_end_col  >=  first_global_pt_EW) .AND.            &
     &        (global_end_col  <=  last_global_pt_EW)) THEN
! This processor contains the EW end of the subarea
            local_end_col=global_end_col-first_global_pt_EW+            &
     &                   halosize(1,halo_type)+1
          ELSE
! This processor is to the left of the end of the subarea
            local_end_col=halosize(1,halo_type)+row_len_nh
          ENDIF

        ELSE

          local_start_col=st_no_data
          local_end_col=st_no_data

        ENDIF

      ENDIF ! is this a fullfield?

 9999 CONTINUE

      RETURN
      END SUBROUTINE GLOBAL_TO_LOCAL_SUBDOMAIN

! Subroutine Interface:

! Function Interface

