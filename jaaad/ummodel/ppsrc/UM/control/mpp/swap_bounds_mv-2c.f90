
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Direct MPI version of swap_bounds to deal with multiple variables
! at once

SUBROUTINE SWAP_BOUNDS_MV(  &
   input_fields,            &        ! Fields to be swapped
   n_multi,                 &        ! The number of Fields to swap
   row_length,              &        ! field size
   halo_x, halo_y)                   ! halos

USE MPL, ONLY :             &
         MPL_REAL,          &
         MPL_STATUS_SIZE

USE SWAPABLE_FIELD_MOD, ONLY : &
    SWAPABLE_FIELD_POINTER_TYPE

IMPLICIT NONE

!  This code performs the same process as swap_bounds but
!  can swap multiple variables at once and hence utilise bandwidth
!  better. This is the MPL version.
!
!  Note that it can only deal with multiple variables with the same
!  row_length and same halo size.
!
! Purpose:
!   This subroutine takes care of all boundary swapping and
!   extending of arrays at the global boundaries. Data is swapped
!   across the poles for any non-zero halo size in the y direction
!   if it is a global model.
!
! Implementation
!   The logic flow is non-trivial!
!   The across pole differencing in particular must be handled carefully!
!   The basic idea is to copy the data to be transferred into a
!   send_buffer array, which is sent to the receiving processor
!   where it arrives in receive_buffer. The receiving processor
!   then copies this data into the appropriate halo region.
!   The East/West halos are done first, then the North/South
!   halos.
!
!   Note that due to the fact that pointers to the data are used,
!   addressing is (1:row_length+2*halo_x, 1:rows+2*halo_y) rather
!   than (1-halo_x:row_length+halo_x, 1-halo_y:rows+halo_y) as used
!   elsewhere. This is unavoidable as pointers change the addressing
!   mode.
!
! Author: Paul Selwood
! Current code owner: Paul Selwood
!
! Comdecks
! DOMTYP contains different model domain types
!
! Author : P.Burton
! History:
! Version  Date      Comment.
! 5.0      15/04/99  New comdeck
! 5.2      15/11/00  add bi_cyclic_lam domain   A. Malcolm

      INTEGER,PARAMETER:: mt_global        = 1
      INTEGER,PARAMETER:: mt_lam           = 2
      INTEGER,PARAMETER:: mt_cyclic_lam    = 3
      INTEGER,PARAMETER:: mt_bi_cyclic_lam = 4
      INTEGER,PARAMETER:: mt_single_column = 5
! DOMTYP end
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

! Arguments:

INTEGER, INTENT(IN)  :: n_multi      ! number of fields to be swapped
INTEGER, INTENT(IN)  :: row_length   ! number of points on row (no halos)
INTEGER, INTENT(IN)  :: halo_x       ! sixze of "i" halo
INTEGER, INTENT(IN)  :: halo_y       ! sixze of "j" halo

! Fields to swap
TYPE(swapable_field_pointer_type), TARGET :: input_fields(n_multi)

! Local scalar variables
INTEGER :: levels             ! number of levels in field
INTEGER :: rows               ! number of rows in field (no halos)
INTEGER :: field_type    ! The grid type of the field (u,v,p)
INTEGER :: i,j,k         ! Spatial loop counters
INTEGER :: info          ! GCOM return code
INTEGER :: length
INTEGER :: i_field
INTEGER :: ierror        ! MPI return code
INTEGER :: nreq_r_ew     ! number of MPI recv requests ew
INTEGER :: nreq_s_ew     ! number of MPI send requests ew
INTEGER :: nreq_r_ns     ! number of MPI recv requests ns
INTEGER :: nreq_s_ns     ! number of MPI send requests ns
INTEGER :: full_row_length       ! length of row including halos
INTEGER :: max_full_rows         ! max no of rows including halo rows
INTEGER :: EW_Halo_Size          ! size of EW halo region
INTEGER :: NS_Halo_Size          ! size of NS halo region
INTEGER :: west_halo_source      ! first column of data for west halo
INTEGER :: east_halo_source      ! first column of data for east halo
INTEGER :: south_halo_source
INTEGER :: north_halo_source
INTEGER :: half_full_row_length  ! half of the total EW dimension
INTEGER :: half_row_length       ! half of the data (no halo) EW dim
INTEGER :: index_2_start         ! start address of index2 in buffers
INTEGER :: max_levels            ! maximum levels used over all fields
INTEGER :: max_rows              ! maximum rows used over all fields
INTEGER :: buffer_size           ! size of send/recv buffers
INTEGER :: MY_COMM               ! Communicator

LOGICAL :: l_vector      ! TRUE if a vector field


! Local arrays
INTEGER :: istat(MPL_STATUS_SIZE,4)  ! MPI status
INTEGER :: ireq_r_ew(4)              ! MPI requests
INTEGER :: ireq_s_ew(4)              ! MPI requests
INTEGER :: ireq_r_ns(4)              ! MPI requests
INTEGER :: ireq_s_ns(4)              ! MPI requests
INTEGER :: north_off(n_multi)        ! Offsets to use when copying data
INTEGER :: south_off(n_multi)        ! to send around poles

REAL, POINTER :: field(:, :, :)

! Send and receive buffers to be allocated dynamically
REAL, ALLOCATABLE :: send_buffer(:)
REAL, ALLOCATABLE :: receive_buffer(:)

LOGICAL :: change_sign(n_multi)     ! .TRUE. if sign change across pole


! Statement functions for addressing the buffer arrays
INTEGER :: EW_address
INTEGER :: NS_address
 
! Variables used in statement functions

INTEGER :: row
INTEGER :: point
INTEGER :: halo
INTEGER :: level
INTEGER :: index
INTEGER :: j_field

EW_address(row,halo,level,index,j_field)=                                 &
   n_multi*(index-1)*(max_levels*EW_Halo_Size) +                          &
   (j_field-1)*(max_levels*EW_Halo_Size) +                                &
   (level-1)*EW_Halo_Size +                                               &
   (halo-1)*max_full_rows +                                               &
   row

NS_address(point,halo,level,index,j_field)=                               &
   n_multi*(index-1)*(max_levels*NS_Halo_Size) +                          &
   (j_field-1)*(max_levels*NS_Halo_Size) +                                &
   (level-1)*NS_Halo_Size +                                               &
   (halo-1)*full_row_length +                                             &
   point


!------------------------------------------------------------------
! 0.0 Check if there is anything to do
IF (((halo_x == 0) .AND. (halo_y == 0)) .OR. n_multi == 0) GOTO 9999

!------------------------------------------------------------------
! 1.0 Initialise variables

! Maximum rows and levels
max_rows   = 0
max_levels = 0
DO i_field=1, n_multi
  max_levels= max(input_fields(i_field) % levels, max_levels)
  max_rows  = max(input_fields(i_field) % rows,   max_rows)
END DO

full_row_length = row_length + 2*halo_x
max_full_rows   = max_rows + 2*halo_y

half_full_row_length = full_row_length/2
half_row_length      = row_length/2

EW_Halo_Size         = max_full_rows   * halo_x
NS_Halo_Size         = full_row_length * halo_y
 
 
! Allocate buffers for communication
buffer_size =  n_multi * 2*max_levels *                                   &
               (max(max_rows*halo_x, row_length*halo_y) +                 &
               (2 * halo_x * halo_y))

ALLOCATE( send_buffer(buffer_size) )
ALLOCATE( receive_buffer(buffer_size) )

! Get communicator we will be using from GCOM
CALL GC_GET_COMMUNICATOR(MY_COMM,IERROR)

!------------------------------------------------------------------
! 2.0 East-West communications
!---------------------------------------
! 2.1 Simple case of only one processor
!     in the East-West direction

IF (halo_x > 0) THEN

nreq_s_ew = 0
nreq_s_ns = 0
IF (nproc_x == 1) THEN ! only 1 processor East-West

  IF (bound(1) == BC_CYCLIC) THEN        ! cyclic boundary conditions
    west_halo_source=row_length+1        ! copy from opposite end
    east_halo_source=halo_x+1            ! of each row

    DO i_field=1, n_multi
      field       => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      rows        = input_fields(i_field) % rows
!CDIR NOVECTOR
      DO I=1,halo_x
        DO K=1,levels
          DO J=1,rows + (2 * halo_y)

            ! Fill Western halo
            field(I,J,K) = field(west_halo_source+I-1,J,K)

            ! Fill Eastern halo
            field(row_length+halo_x+I,J,K) =                              &
                            field(east_halo_source+I-1,J,K)

          END DO ! J
        END DO ! K
      END DO ! I

    END DO ! loop over fields

  END IF   !  bound(1) == BC_CYCLIC

!---------------------------------------
! 2.1 Now the more common case of having
!     a number of processors in the
!     East-West direction

ELSE ! If there is more than 1 processor East-West

!---------------------------------------
! 2.1.1 Copy the data into send_buffer

  DO i_field=1, n_multi
    field       => input_fields(i_field) % field
    levels      = input_fields(i_field) % levels
    rows        = input_fields(i_field) % rows
  
    DO K=1,levels
      DO J=1,rows + (2*halo_y)
        DO I=1,halo_x

          ! Copy stuff from the Western side of the grid
          send_buffer(EW_address(J,I,K,1,i_field)) =                      &
                      field(halo_x+I,J,K)

          ! Copy stuff from the Eastern side of the grid
          send_buffer(EW_address(J,I,K,2,i_field)) =                      &
                      field(row_length+I,J,K)

        END DO ! I
      END DO ! J
    END DO ! K

  END DO ! loop over fields

!---------------------------------------
! 2.1.2 Send and receive the data

!---------------------------------------
! 2.1.2.1 Special case of 2 processors
!         East-West - both sides are
!         sent to the same processor
!         (if cyclic BC)

  IF ((nproc_x == 2) .AND. (bound(1) == BC_CYCLIC)) THEN

    length = 2 * n_multi * EW_halo_size * max_levels

    CALL GC_RSEND(1,length,neighbour(PEast),info,                         &
                   receive_buffer, send_buffer)

    CALL GC_RRECV(1,length,neighbour(PWest),info,                         &
                   receive_buffer, send_buffer)

!---------------------------------------
! 2.1.2.2 More general case when there
!         are more than 2 processors
!         in the EW direction. Each
!         halo can be sent seperately
!         as there is no danger of them
!         both being sent to the same
!         processor.

  ELSE ! more than 2 processors East-West

    ! index_2_start points to the start of the second index
    ! within the buffer arrays. The first index contains
    ! data sent from the Western side, the second index
    ! contains data sent from the Eastern side

    index_2_start=EW_address(1,1,1,2,1)

    length = n_multi * EW_halo_size * max_levels

    nreq_r_ew=0
    IF (neighbour(PEast) /= NoDomain) THEN

      ! Receive from East
      nreq_r_ew=nreq_r_ew+1
      CALL MPL_IRECV(receive_buffer,                                      &
                     length, MPL_REAL, neighbour(PEast), 2,               &
                     MY_COMM, ireq_r_ew(nreq_r_ew), ierror)

    END IF

    IF (neighbour(PWest) /= NoDomain) THEN

      ! Receive from West
      nreq_r_ew=nreq_r_ew+1
      CALL MPL_IRECV(receive_buffer(index_2_start),                       &
                     length, MPL_REAL, neighbour(PWest), 3,               &
                     MY_COMM, ireq_r_ew(nreq_r_ew), ierror)

    END IF

    nreq_s_ew=0
    IF (neighbour(PEast) /= NoDomain) THEN
      
      ! Send East
      nreq_s_ew=nreq_s_ew+1
      CALL MPL_ISEND(send_buffer(index_2_start),                          &
                     length, MPL_REAL, neighbour(PEast), 3,               &
                     MY_COMM, ireq_s_ew(nreq_s_ew), ierror)
    END IF

    IF (neighbour(PWest) /= NoDomain) THEN
    
      ! Send West
      nreq_s_ew=nreq_s_ew+1
      CALL MPL_ISEND(send_buffer,                                         &
                     length, MPL_REAL, neighbour(PWest), 2,               &
                     MY_COMM, ireq_s_ew(nreq_s_ew), ierror)
    
    END IF
      
    CALL MPL_WAITALL( nreq_r_ew, ireq_r_ew, istat, ierror )

  END IF ! test on numbers of processors East-West

  CALL MPL_WAITALL( nreq_s_ew, ireq_s_ew, istat, ierror )

!---------------------------------------
! 2.1.2 Fill the halos with data

  IF (neighbour(PEast) /= NoDomain) THEN

    ! unpack data from receive_buffer into field

    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      rows        = input_fields(i_field) % rows 

      DO K=1,levels
        DO J=1,rows + (2*halo_y)
          DO I=1,halo_x
            field(row_length+halo_x+I,J,K)=                               &
               receive_buffer(EW_address(J,I,K,1,i_field))
          END DO
        END DO
      END DO
    END DO ! loop over fields

  ELSEIF (sb_Model_domain == mt_global) THEN
       ! No neighbour to my East (ie. at edge of the domain
       ! and it's not a cyclic boundary condition)
       ! Just copy data from last column
! NOTE: This block of code should never be executed. It's being left
!       in so that it can be used in the future by changing the logic
!       in the ELSE statement above.

    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      rows        = input_fields(i_field) % rows 

      DO K=1,levels
        DO J=1,rows+ (2*halo_y)
          DO I=1,halo_x
            field(row_length+halo_x+I,J,K)=field(row_length+halo_x,J,K)
          END DO
        END DO
      END DO
    END DO ! loop over fields

  END IF ! IF (neighbour(PEast) /= NoDomain)

  IF (neighbour(PWest) /= NoDomain) THEN

    ! unpack data from receive_buffer into field

    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      rows        = input_fields(i_field) % rows

      DO K=1,levels
        DO J=1,rows+ (2*halo_y)
          DO I=1,halo_x
            field(I,J,K)=                                                 &
               receive_buffer(EW_address(J,I,K,2,i_field))
          END DO
        END DO
      END DO
    END DO ! loop over fields

  ELSEIF (sb_Model_domain == mt_global) THEN
       ! No neighbour to my West (ie. at edge of the domain
       ! and it's not a cyclic boundary condition)
       ! Just copy data from first column
! NOTE: This block of code should never be executed. It's being left
!       in so that it can be used in the future by changing the logic
!       in the ELSE statement above.

    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      rows        = input_fields(i_field) % rows

      DO K=1,levels
        DO J=1,rows+ (2*halo_y)
          DO I=1,halo_x
            field(I,J,K)=field(1+halo_x,J,K)
          END DO
        END DO
      END DO
    END DO ! loop over fields
 
  END IF ! IF (neighbour(PWest) /= NoDomain)

END IF ! IF (nproc_x == 1)

END IF ! halo_x > 0

!------------------------------------------------------------------
! 3.0 North-South communications
!  section of code for cyclic N-S boundary conditions
IF (halo_y > 0) THEN

IF(bound(2) == BC_CYCLIC) THEN
  
  IF (nproc_y == 1) THEN ! only 1 processor north-south
  ! 2.1 Simple case of only one processor
  
    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      rows        = input_fields(i_field) % rows

      south_halo_source=rows+1              ! copy from opposite end
      north_halo_source=1                   ! of each column

      DO K=1,levels
        DO J=1,halo_y
          DO I=1 ,row_length+ (2*halo_x)
  
            ! Fill southern halo
            field(I,J,K)=field(I,south_halo_source+J,K)
  
            ! Fill northern halo
            field(I,rows+halo_y+J,K)=field(I,north_halo_source+J,K)
  
          END DO ! I
        END DO ! J
      END DO ! K
    END DO ! loop over fields
  
  !---------------------------------------
  ! 2.1 Now the more common case of having
  !     a number of processors in the
  !     North-South direction
  
  ELSE ! If there is more than 1 processor north_south
  
  !---------------------------------------
  ! 2.1.1 Copy the data into buf_send
  
    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      rows        = input_fields(i_field) % rows

      DO K=1,levels
        DO J=1,halo_y
          DO I=1, row_length + (2*halo_x)
  
            ! Copy stuff from the southern side of the grid
            send_buffer(NS_address(I,J,K,1,i_field)) =                    &
               field(I,halo_y+J,K)
  
            ! Copy stuff from the northern side of the grid
            send_buffer(NS_address(I,J,K,2,i_field)) =                    &
               field(I,rows+J,K)
  
          END DO ! I
        END DO ! J
      END DO ! K
    END DO ! loop over fields
  
  !---------------------------------------
  ! 2.1.2 Send and receive the data
  
  !---------------------------------------
  ! 2.1.2.1 Special case of 2 processors
  !         north-south - both sides are
  !         sent to the same processor
  !         as cyclic BC
  
   IF ( nproc_y == 2 ) THEN
  
     length=2*n_multi*NS_halo_size*max_levels
  
      CALL GC_RSEND(1,length,neighbour(PNorth),info,                      &
                     receive_buffer, send_buffer)
  
      CALL GC_RRECV(1,length,neighbour(PSouth),info,                      &
                     receive_buffer, send_buffer)
  
  !---------------------------------------
  ! 2.1.2.2 More general case when there
  !         are more than 2 processors
  !         in the NS direction. Each
  !         halo can be sent seperately
  !         as there is no danger of them
  !         both being sent to the same
  !         processor.
  
    ELSE ! more than 2 processors North-South
  
      ! index_2_start points to the start of the second index
      ! within the buffer arrays. The first index contains
      ! data sent from the southern side, the second index
      ! contains data sent from the northern side
  
      index_2_start=NS_address(1,1,1,2,1)
  
      length = n_multi * ns_halo_size * max_levels
  
      IF (neighbour(PSouth) /= NoDomain) THEN
  
        ! Send south
        CALL GC_RSEND(2,length,neighbour(Psouth),info,                    &
                       receive_buffer, send_buffer)
      END IF
  
      IF (neighbour(Pnorth) /= NoDomain) THEN
  
        ! Send north
        CALL GC_RSEND(3,length,neighbour(Pnorth),info,                    &
                       receive_buffer(index_2_start),                     &
                       send_buffer(index_2_start))
      END IF
  
      IF (neighbour(Pnorth) /= NoDomain) THEN
  
        ! Receive from north
        CALL GC_RRECV(2,length,neighbour(Pnorth),info,                    &
                       receive_buffer, send_buffer)
      END IF
  
      IF (neighbour(Psouth) /= NoDomain) THEN
  
        ! Receive from south
        CALL GC_RRECV(3,length,neighbour(Psouth),info,                    &
                       receive_buffer(index_2_start),                     &
                       send_buffer(index_2_start))
  
      END IF
  
    END IF ! test on numbers of processors north-south
  
  !---------------------------------------
  ! 2.1.2 Fill the halos with data
  
    IF (neighbour(Pnorth) /= NoDomain) THEN
  
      ! unpack data from receive_buffer into field
  
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        rows        = input_fields(i_field) % rows

        DO K=1,levels
          DO J=1,halo_y
            DO I=1,row_length+ (2*halo_x)
              field(I,J+halo_y+rows,K)=                                   &
              receive_buffer(NS_address(I,J,K,1,i_field))
                             
            END DO
          END DO
        END DO
      END DO ! loop over fields
  
    END IF ! IF (neighbour(Pnorth) /= NoDomain)
  
    IF (neighbour(Psouth) /= NoDomain) THEN
  
      ! unpack data from receive_buffer into field
  
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        DO K=1,levels
          DO J=1,halo_y
            DO I=1,row_length+ (2*halo_x)
              field(I,J,K)=                                               &
              receive_buffer(NS_address(I,J,K,2,i_field))
            END DO
          END DO
        END DO
      END DO ! loop over fields
  
  
    END IF ! IF (neighbour(Psouth) /= NoDomain)
  
  END IF ! IF (nproc_y == 1)
  
ELSE                 !!! bc_cyclic in NS

  ! Set up some variables
  
  ! Set up the offsets. When copying data that is to be passed over
  ! the pole, on wind (u or v) grid, then copy data one row away
  ! from the pole
  
  DO i_field=1, n_multi
    north_off(i_field)=0
    south_off(i_field)=0
    field_type = input_fields(i_field) % field_type
    IF ((sb_Model_domain == mt_global) .AND.                              &
            ((field_type == fld_type_p) .OR.                              &
             (field_type == fld_type_u))) THEN
  
      IF (at_extremity(PNorth)) north_off(i_field)=1
      IF (at_extremity(PSouth)) south_off(i_field)=1
   
    END IF
  
  ! Set up the sign factor. If l_vector is true and data has been passed
  ! over the poles in a global model, then the variables must change
  ! sign
  
    l_vector = input_fields(i_field) % vector
    IF (.NOT. ((l_vector) .AND. (sb_Model_domain == mt_global))) THEN
      change_sign(i_field)=.FALSE.
    ELSE
      change_sign(i_field)=.TRUE.
    END IF
  END DO
  
  !---------------------------------------
  ! 3.1 Copy data into the send_buffer
  !     But not if:
  !       - Not a global model and the data is at the North/South edge
  !       - A global model at the North/South edge but only 1 processor
  !         in the East-West direction.
  
  IF (.NOT. (at_extremity(PSouth) .AND.                                   &
              ((nproc_x == 1)  .OR.                                       &
              (sb_Model_domain /= mt_global))))                           &
     THEN
  
    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      DO K=1,levels
        DO J=1,halo_y
          DO I=1,row_length+ (2*halo_x)
  
            send_buffer(NS_address(I,J,K,1,i_field)) =                    &
               field(I,J+halo_y+south_off(i_field),K)
  
          END DO
        END DO
      END DO
    END DO ! loop over fields
  END IF
  
  IF (.NOT. (at_extremity(PNorth) .AND.                                   &
              ((nproc_x == 1) .OR.                                        &
              (sb_Model_domain /= mt_global))))                           &
     THEN
  
    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      rows        = input_fields(i_field) % rows
      DO K=1,levels
        DO J=1,halo_y
          DO I=1,row_length+ (2*halo_x)
  
            send_buffer(NS_address(I,J,K,2,i_field)) =                    &
               field(I,rows-north_off(i_field)+J,K)
  
          END DO
        END DO
      END DO
    END DO ! loop over fields
  END IF
  
  !---------------------------------------
  ! 3.2 Send and receive the data
  
  !---------------------------------------
  ! 3.2.1 The special case where nproc_y=1
  !       Both buffers are sent to the
  !       same processor as one message
  
  IF ((nproc_y == 1) .AND. (sb_Model_domain == mt_global) .AND.           &
       (nproc_x > 1)) THEN
  
    length = 2 * n_multi * ns_halo_size * max_levels
  
    CALL GC_RSEND(10,length,neighbour(PSouth),info,                       &
                     receive_buffer, send_buffer)
  
    CALL GC_RRECV(10,length,neighbour(PSouth),info,                       &
                   receive_buffer, send_buffer)
  END IF
  
  !---------------------------------------
  ! 3.2.2 The more general case, where
  !       each buffer is sent to a
  !       different processor
  
  index_2_start=NS_address(1,1,1,2,1)
  
  nreq_r_ns=0
  
  length = n_multi * ns_halo_size * max_levels
  
  IF (at_extremity(PSouth)) THEN
  
    IF ((neighbour(PSouth) /= NoDomain) .AND.                             &
         (neighbour(PSouth) /= mype)) THEN
  
       nreq_r_ns=nreq_r_ns+1
       CALL MPL_IRECV(receive_buffer(index_2_start),                      &
         length, MPL_REAL, neighbour(PSouth), 11, MY_COMM,                &
         ireq_r_ns(nreq_r_ns), ierror)
  
    END IF
  
  ELSE ! not at the South
  
     nreq_r_ns=nreq_r_ns+1
     CALL MPL_IRECV(receive_buffer(index_2_start),                        &
       length, MPL_REAL, neighbour(PSouth), 14, MY_COMM,                  &
       ireq_r_ns(nreq_r_ns), ierror)
  
  END IF
  
  IF (at_extremity(PNorth)) THEN
  
    IF ((neighbour(PNorth) /= NoDomain) .AND.                             &
         (neighbour(PNorth) /= mype)) THEN
  
       nreq_r_ns=nreq_r_ns+1
       CALL MPL_IRECV(receive_buffer,                                     &
         length, MPL_REAL, neighbour(PNorth), 13, MY_COMM,                &
         ireq_r_ns(nreq_r_ns), ierror)
  
    END IF
  
  ELSE
  
       nreq_r_ns=nreq_r_ns+1
       CALL MPL_IRECV(receive_buffer,                                     &
         length, MPL_REAL, neighbour(PNorth), 12, MY_COMM,                &
         ireq_r_ns(nreq_r_ns), ierror)
  
  END IF
  
  nreq_s_ns=0
  IF (at_extremity(PSouth)) THEN
  
    IF ((neighbour(PSouth) /= NoDomain) .AND.                             &
         (neighbour(PSouth) /= mype)) THEN
  
       nreq_s_ns=nreq_s_ns+1
       CALL MPL_ISEND(send_buffer,                                        &
         length, MPL_REAL, neighbour(PSouth), 11, MY_COMM,                &
         ireq_s_ns(nreq_s_ns),ierror)
  
    END IF
  
  ELSE ! not at the South
  
       nreq_s_ns=nreq_s_ns+1
       CALL MPL_ISEND(send_buffer,                                        &
         length, MPL_REAL, neighbour(PSouth), 12, MY_COMM,                &
         ireq_s_ns(nreq_s_ns),ierror)
  
  END IF
  
  IF (at_extremity(PNorth)) THEN
  
    IF ((neighbour(PNorth) /= NoDomain) .AND.                             &
         (neighbour(PNorth) /= mype)) THEN
  
     nreq_s_ns=nreq_s_ns+1
     CALL MPL_ISEND(send_buffer(index_2_start),                           &
       length, MPL_REAL, neighbour(PNorth), 13, MY_COMM,                  &
       ireq_s_ns(nreq_s_ns),ierror)
  
    END IF
  
  ELSE ! not at the North
  
     nreq_s_ns=nreq_s_ns+1
     CALL MPL_ISEND(send_buffer(index_2_start),                           &
       length, MPL_REAL, neighbour(PNorth), 14, MY_COMM,                  &
       ireq_s_ns(nreq_s_ns),ierror)
  
  END IF
  
  CALL MPL_WAITALL( nreq_r_ns, ireq_r_ns, istat, ierror )
  !---------------------------------------
  ! 3.3 Fill the halos with data
  
  !---------------------------------------
  ! 3.3.1 Southern halo
  
  IF (at_extremity(PSouth)) THEN
  
    IF (neighbour(PSouth) == NoDomain) THEN
      IF (sb_Model_domain /= mt_lam) THEN
  
      ! Just copy adjacent rows into halo area
  ! NOTE: This block of code should never be executed. It's being left
  !       in so that it can be used in the future by changing the logic
  !       in the IF statement above.
  
        DO i_field=1, n_multi
          field      => input_fields(i_field) % field
          levels      = input_fields(i_field) % levels
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,row_length+(2*halo_x)
                field(I,J,K)=field(I,halo_y+1,K)
              END DO
            END DO
          END DO
        END DO ! loop over fields
      END IF ! IF (sb_Model_domain /= mt_lam)
  
    ELSEIF (neighbour(PSouth) == mype) THEN
    ! Local across pole difference
  
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        IF (change_sign(i_field)) THEN
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,half_full_row_length
                field(half_row_length+halo_x+I,1-J+halo_y,K)=             &
                   -field(I+halo_x,J+halo_y+south_off(i_field),K)
                field(I,1-J+halo_y,K)=                                    &
                   -field(half_row_length+I,J+halo_y+                     &
                   south_off(i_field),K)
              END DO
            END DO
          END DO
        ELSE ! don't change sign
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,half_full_row_length
                field(half_row_length+halo_x+I,1-J+halo_y,K)=             &
                   field(I+halo_x,J+halo_y+south_off(i_field),K)
                field(I,1-J+halo_y,K)=                                    &
                   field(half_row_length+I,J+halo_y+                      &
                   south_off(i_field),K)
              END DO
            END DO
          END DO
        END IF ! IF (change_sign(i_field))
      END DO ! loop over fields
  
    ELSE ! data is receive_buffer(index 2)
  
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        IF (change_sign(i_field)) THEN
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,row_length+(2*halo_x)
                field(I,1-J+halo_y,K)=                                    &
                -receive_buffer(NS_address(I,J,K,2,i_field))
              END DO
            END DO
          END DO
        ELSE ! don't change sign
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,row_length+(2*halo_x)
                field(I,1-J+halo_y,K)=                                    &
                 receive_buffer(NS_address(I,J,K,2,i_field))
              END DO
            END DO
          END DO
        END IF !  IF (change_sign(i_field))
      END DO ! loop over fields
  
    END IF ! What type of South extremity
  
  ELSE ! IF (at_extremity(PSouth)
  
    ! not at a South extremity
  
    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      DO K=1,levels
        DO J=1,halo_y
          DO I=1,row_length+(2*halo_x)
            field(I,J,K)=                                                 &
               receive_buffer(NS_address(I,J,K,2,i_field))
          END DO
        END DO
      END DO
    END DO ! loop over fields
  
  END IF ! IF (at_extremity(PSouth)
  
  
  !---------------------------------------
  ! 3.3.2 Northern halo
  
  IF (at_extremity(PNorth)) THEN
  
    IF (neighbour(PNorth) == NoDomain) THEN
      IF (sb_Model_domain /= mt_lam) THEN
      ! Just copy adjacent rows into halo area
  ! NOTE: This block of code should never be executed. It's being left
  !       in so that it can be used in the future by changing the logic
  !       in the IF statement above.
  
        DO i_field=1, n_multi
          field      => input_fields(i_field) % field
          levels      = input_fields(i_field) % levels
          rows        = input_fields(i_field) % rows
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,row_length+(2*halo_x)
                field(I,rows+halo_y+J,K)=field(I,rows+halo_y,K)
              END DO
            END DO
          END DO
        END DO ! loop over fields
      END IF
  
    ELSEIF (neighbour(PNorth) == mype) THEN
    ! Local across pole difference
  
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        rows        = input_fields(i_field) % rows
        IF (change_sign(i_field)) THEN
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,half_full_row_length
  
                field(half_row_length+halo_x+I,rows+halo_y+J,K)=          &
                   -field(I+halo_x,rows+halo_y-J+1-north_off(i_field),K)
                field(I,rows+J+halo_y,K)=                                 &
                   -field(half_row_length+I,                              &
                        rows-J+halo_y+1-north_off(i_field),K)
  
              END DO
            END DO
          END DO
        ELSE ! don't change sign
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,half_full_row_length
   
                field(half_row_length+halo_x+I,rows+halo_y+J,K)=          &
                   field(I+halo_x,rows+halo_y-J+1-north_off(i_field),K)
                field(I,rows+J+halo_y,K)=                                 &
                   field(half_row_length+I,                               &
                       rows-J+1+halo_y-north_off(i_field),K)
  
              END DO
            END DO
          END DO
        END IF ! IF (change_sign(i_field))
      END DO ! loop over fields
  
    ELSE ! data is receive_buffer(index 1)
  
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        rows        = input_fields(i_field) % rows
        IF (change_sign(i_field)) THEN
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,row_length+(2*halo_x)
                field(I,rows+halo_y+J,K)=                                 &
                   -receive_buffer(                                       &
                   NS_address(I,halo_y-J+1,K,1,i_field))
              END DO
            END DO
          END DO
        ELSE ! don't change sign
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,row_length+(2*halo_x)
                field(I,rows+halo_y+J,K)=                                 &
                   receive_buffer(                                        &
                   NS_address(I,halo_y-J+1,K,1,i_field))
              END DO
            END DO
          END DO
        END IF ! IF (change_sign(i_field))
      END DO ! loop over fields
  
    END IF ! What type of North extremity
  
  ELSE ! IF (at_extremity(PNorth)
  
    ! not at a North extremity
  
    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      rows        = input_fields(i_field) % rows
      DO K=1,levels
        DO J=1,halo_y
          DO I=1,row_length+ (2*halo_x)
            field(I,rows+halo_y+J,K)=                                     &
               receive_buffer(NS_address(I,J,K,1,i_field))
          END DO
        END DO
      END DO
    END DO ! loop over fields
  END IF ! IF (at_extremity(PNorth)
END IF ! BC_CYCLIC  for NS
  
CALL MPL_WAITALL( nreq_s_ns, ireq_s_ns, istat, ierror )

END IF ! halo_y > 0
  
9999 CONTINUE

IF (ALLOCATED(send_buffer))    DEALLOCATE(send_buffer)
IF (ALLOCATED(receive_buffer)) DEALLOCATE(receive_buffer)

RETURN
END SUBROUTINE SWAP_BOUNDS_MV
