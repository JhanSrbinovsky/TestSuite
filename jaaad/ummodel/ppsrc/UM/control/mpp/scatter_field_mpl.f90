
! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT Met Office. ALL RIGHTS RESERVED.
! For further details please refer to the file copyright.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ scatters a field from one processor to many processors
!
! SUBROUTINE interface:
SUBROUTINE scatter_field_mpl(  local_field,     global_field,           &
                               local_row_len,   local_rows,             &
                               global_row_len,  global_rows,            &
                               grid_type,       halo_type,              &
                               scatter_pe,      proc_group,             &
                               icode,           cmessage)

USE mpl, ONLY: &
    MPL_REAL

IMPLICIT NONE

!
! Description:
!  Takes a model field which is stored entirely on one processor
!  and distributes it over a group of processors using the
!  standard UM decomposition.
!
! Method:
!   Copies each processors data into contiguous space in a 1D
!   array and then uses mpl_scatterv to do the collective operation.
!
!
! Current code owner: Paul Selwood
!
! SUBROUTINE arguments:

INTEGER, INTENT(IN) :: local_row_len  ! length of rows in local part of field
INTEGER, INTENT(IN) :: local_rows     ! number of rows in local part of field
INTEGER, INTENT(IN) :: global_row_len ! length of rows in global field
INTEGER, INTENT(IN) :: global_rows    ! number of rows in global field
INTEGER, INTENT(IN) :: grid_type      ! type (p,u or v) of grid
INTEGER, INTENT(IN) :: halo_type      ! halo type (hence width) of grid
INTEGER, INTENT(IN) :: scatter_pe     ! processor to scatter global
                                        ! field from
INTEGER, INTENT(IN) :: proc_group     ! group id of processors involved
                                        ! here
INTEGER, INTENT(OUT) :: icode         ! return code

REAL, INTENT(OUT) :: local_field(local_row_len*local_rows)
                                          ! local part of field
REAL, INTENT(INOUT) ::  global_field(global_row_len*global_rows)
                                        ! (on pe scatter_pe) global field
  
CHARACTER*80 :: cmessage         ! out error message

! parameters and common blocks

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

! local variables

REAL :: global_field_copy( global_row_len * global_rows )
                         ! Space for re-ordered copy of global field

INTEGER :: info          ! Status return from MPL/GCOM
INTEGER :: i,j, iproc    ! Loopers
INTEGER :: igpos         ! Position pointer
INTEGER :: igpos1        ! Position pointer
INTEGER :: ipos          ! Position pointer
INTEGER :: my_comm       ! Current communicator

INTEGER :: counts(0:nproc-1)               ! Counts for Scatter
INTEGER :: displacements(0:nproc-1)        ! Displacements for Scatter

INTEGER, SAVE :: old_global_row_len = -1234
INTEGER, SAVE :: old_global_rows    = -1234
INTEGER, SAVE :: old_proc_group     = -1234
INTEGER, SAVE :: old_scatter_pe     = -1234
INTEGER, SAVE :: old_decomp         = -1234
INTEGER, SAVE :: old_grid_type      = -1234
INTEGER, SAVE :: old_halo_type      = -1234
INTEGER, SAVE :: recv_type          = -1234        ! receive data type

CHARACTER (LEN=*), PARAMETER :: routinename='scatter_field_mpl'
!-------------------------------------------------------

! 0.0 Error checking. 
! Do we recognise the field type?
IF (grid_type  ==  fld_type_unknown) THEN
  WRITE (6,*) 'SCATTER_FIELD_MPL : bad field type'
  WRITE (6,*) 'field will not be scattered.'
  WRITE (cmessage,*) 'SCATTER_FIELD_MPL : bad field type'
  icode = 1
  ! DEPENDS ON: ereport
  CALL ereport(routinename, icode, cmessage)
END IF

! We can only do scatters with global processor groups with this
! version currently
IF (proc_group /= gc_all_proc_group) THEN
  WRITE (cmessage, *) 'Can only scatter in an all processor group'
  icode = 2
  ! DEPENDS ON: ereport
  CALL ereport(routinename, icode, cmessage)
END IF
  
! 1.0 Setup MPI datatypes. Can we use the same as last time round?
IF ((global_row_len  /=  old_global_row_len) .OR.                     &
    (global_rows     /=  old_global_rows   ) .OR.                     &
    (proc_group      /=  old_proc_group    ) .OR.                     &
    (scatter_pe      /=  old_scatter_pe    ) .OR.                     &
    (grid_type       /=  old_grid_type     ) .OR.                     &
    (halo_type       /=  old_halo_type     ) .OR.                     &
    (current_decomp_type  /=  old_decomp  )) THEN

  ! Need to create a new receive type
  ! First get rid of the old one if this isn't the first call
  IF (recv_type /= -1234) THEN
    CALL MPL_type_free(recv_type, info)
  END IF
 
  ! New type is a 2D block with a local row length + 2 * halo stride.
  ! This is used to describe the data size and layout of the local 
  ! data without including halos via an MPI datatype to simplify 
  ! communications and reduce copies. Will need to specify the real
  ! start of the data in the communications call.
  CALL MPL_type_vector(blsize(2, grid_type),  blsize(1, grid_type),   &
                       local_row_len,  MPL_REAL, recv_type, info)
  CALL MPL_type_commit(recv_type, info)

END IF  ! New MPI types

! 2.0 set up the scattering array, counts and displacements
!     Note that global data does not include halos.
!     This copy takes the conventional global field ordering and copies
!     into a per-PE contiguous data layout (all of PE0 data,
!     followed by all of PE1 data, followed by ...).
  IF  (mype == scatter_pe) THEN
    ipos = 1 ! start of array
    DO iproc = 0, nproc-1
      DO j = 1, g_blsize(2, grid_type, iproc)
        igpos1 =   (g_datastart_f(2,grid_type,iproc) + j - 2) *        &
                    global_row_len - 1 +                               &
                    g_datastart_f(1,grid_type, iproc)
        DO i = 1, g_blsize(1, grid_type, iproc)

          igpos =  i + igpos1
          global_field_copy(ipos) = global_field(igpos)
          ipos = ipos + 1

        END DO
      END DO
    END DO

    ! Make sure we've got the full array
    IF (ipos /= global_rows * global_row_len + 1) THEN
      WRITE(cmessage, *) 'Addressing gone wrong in field copy'
      icode = 3
      ! DEPENDS ON: ereport
      CALL ereport(routinename, icode, cmessage)
    END IF
  
  END IF

  DO iproc = 0, nproc-1
    counts(iproc) = g_blsize(1, grid_type, iproc) *                   &
                    g_blsize(2, grid_type, iproc)

    IF (iproc == 0) THEN
      displacements(iproc) = 0
    ELSE
      displacements(iproc)=SUM(counts(0:iproc-1))
    END IF

  END DO ! nprocs

! 3.0 Do the exchange of data
CALL gc_get_communicator(my_comm, info)

! ipos is the location in the recv array of the first datapoint
ipos = halosize(2,halo_type) * local_row_len                         &
     + halosize(1,halo_type) + 1
CALL MPL_scatterv(global_field_copy, counts, displacements,          &
                  MPL_REAL, local_field(ipos),                       &
                  1, recv_type, scatter_pe, my_comm, info)


RETURN
END SUBROUTINE scatter_field_mpl