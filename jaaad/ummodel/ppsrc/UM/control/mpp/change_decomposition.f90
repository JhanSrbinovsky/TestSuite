

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Parallel UM: Select a new decomposition
!
! Subroutine Interface:
      SUBROUTINE CHANGE_DECOMPOSITION(new_decomp,icode)

      IMPLICIT NONE

!
! Description:
! Sets up the PARVARS common blocks with the correct information for
! decomposition new_decomp
!
! Method:
! If new_decomp is already the current decomposition, exit and do
! nothing.
! If decomposition new_decomp has not been initialised, print a
! message and exit, with icode=-1.
! Otherwise, copy the information from the decomp_db arrays in the
! DECOMPDB comdeck into the PARVARS comdecks arrays.
!
! Current Code Owner : Paul Burton
!
! History:
!  Model    Date     Modification history from model version 4.2
!  version
!  4.2      22/8/96  New deck created for mpp code.  P.Burton
!  4.3      17/02/97 Changed ICODE to a positive error no. P.Burton
!  5.0      12/04/99 - added new halosize PARVARS array and halo_i/j
!                    - glsize now has extra Nfld_max dimension
!                    - lasize now has extra NHalo_max and Nfld_max
!                      dimensions
!                    - blsizep/u replace by blsize with new
!                      Nfld_max dimension
!                    - attop etc replaced with at_extremity
!                                                         P.Burton
!  5.1      02/02/00 Code for new g_pe_index variables   P.Burton
!  5.2      02/08/00 Code for new g_at_extremity variable  P.Burton
!  5.3      14/09/01 Added sb_model_domain variable    P.Burton
!  5.3      22/11/01 Enable mpp as the only option for
!                    small executables         E.Leung
!  5.5      06/08/00 Modification for parallelisation of WAM
!                    Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!  5.5      30/01/03 River routing support. P.Selwood
!  6.0  17/09/03  Add def for new NEC opt section c96_1c. R Barnes
!  6.2      23/11/05  Removed all references to the wavemodel.
!                     T.Edwards
!
! Subroutine arguments:

      INTEGER                                                           &
     &  new_decomp                                                      &
                     ! IN : new decomposition to use
     &, icode        ! OUT: return code (-1 is failure)

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
! DECOMPTP comdeck
!
! Description
!
! Magic numbers indicating decomposition types.
! These numbers are used to index the arrays defined in the
! DECOMPDB comdeck, and are required as an argument to
! the CHANGE_DECOMPOSITION subroutine.
!
! Current code owner : P.Burton
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 4.2       19/08/96  Original code.   P.Burton
! 4.3       17/02/97  Added new ocean decomposition decomp_nowrap_ocean
!                     which does not contain extra wrap points at
!                     start and end of row.                  P.Burton
! 5.5       04/08/00  Modification for parallelisation of WAM
!                   Author:  Bob Carruthers, Cray UK Inc(D.Holmes-Bell)

! Magic Numbers indicating decomposition types

      INTEGER                                                           &
     &  max_decomps                                                     &
                               ! maximum number of decompositions
     &, decomp_unset                                                    &
                               ! no decomposition selected
     &, decomp_standard_atmos                                           &
                               ! standard 2D atmosphere
!                              ! decomposition
     &, decomp_standard_ocean                                           &
                               ! standard 1D ocean decomposition
     &, decomp_nowrap_ocean                                             &
                               ! 1D ocean without extra wrap-around
!                              ! points at ends of each row
     &, decomp_smexe                                                    &
     &, decomp_standard_wave   ! standard 1D WAM Wave Model
                               ! decomposition

      PARAMETER (                                                       &
     &  max_decomps=5                                                   &
     &, decomp_unset=-1                                                 &
     &, decomp_standard_atmos=1                                         &
     &, decomp_standard_ocean=2                                         &
     &, decomp_nowrap_ocean=3                                           &
     &, decomp_smexe=4                                                  &
     &, decomp_standard_wave=5)

! End of DECOMPTP comdeck
! DECOMPDB comdeck
!
! Description:
!
! DECOMPDB comdeck (Decomposition Database) contains information
! describing the various decompositions used by the mpp-UM
! The CHANGE_DECOMPOSITION subroutine can be used to select
! a particular decomposition (which copies the appropriate
! decomposition information into the PARVARS common block).
!
! Requires comdeck PARVARS to be *CALLed before it.
!
! Current code owner : P.Burton
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 4.2       19/08/96  Original code.   P.Burton
! 5.0       12/04/99  - added dimension Nfld_max to decomp_db_glsize
!                     - added dimension NHalo_max to decomp_db_halosize
!                     - added dimension NHalo_max to decomp_db_g_lasize
!                     - added dimension Nfld_max to
!                       decomp_db_g_lasize
!                     - replace blsizep/u variables by blsize with new
!                       Nfld_max dimension
!                                                            P.Burton
! 5.1       27/01/00  - changed g_pe_index to g_pe_index_EW
!                     - added g_pe_index_NS
!                                                  P.Burton
! 5.3       27/01/00  - Moved data statement into blkdata. A Van der Wal
! 5.3       14/09/01  Added model_domain variable.  P.Burton
! 5.5       30/01/03  Generalised datastart. P.Selwood.

! Common blocks containing information about each decomposition
! (For description of variables see the PARVARS comdeck)

      INTEGER :: decomp_db_bound(Ndim_max,max_decomps)                  &
     &, decomp_db_sb_model_domain(max_decomps)
      INTEGER :: decomp_db_glsize(Ndim_max,Nfld_max,max_decomps)
      INTEGER :: decomp_db_gridsize(Ndim_max,max_decomps)
      INTEGER :: decomp_db_g_lasize(Ndim_max,Nfld_max,                  &
     &                     NHalo_max,0:maxproc,max_decomps)
      INTEGER :: decomp_db_g_blsize(Ndim_max,Nfld_max,0:maxproc,        &
     &                          max_decomps)
      INTEGER :: decomp_db_g_datastart(Ndim_max,0:maxproc,max_decomps)
      INTEGER :: decomp_db_g_datastart_f(Ndim_max,Nfld_max,0:maxproc,   &
     &                                   max_decomps)
      INTEGER :: decomp_db_g_gridpos(Ndim_max,0:maxproc,max_decomps)
      INTEGER :: decomp_db_g_pe_index_EW(1-Max_Halo_Size:               &
     &  ROW_LENGTH_MAX+Max_Halo_Size, max_decomps)
      INTEGER :: decomp_db_g_pe_index_NS(1-Max_Halo_Size:               &
     &  ROWS_MAX+Max_Halo_Size, max_decomps)
      INTEGER :: decomp_db_halosize(Ndim_max,NHalo_max,max_decomps)
      INTEGER :: decomp_db_neighbour(4,max_decomps)
      INTEGER :: decomp_db_first_comp_pe(max_decomps)
      INTEGER :: decomp_db_last_comp_pe(max_decomps)
      INTEGER :: decomp_db_nproc(max_decomps)
      INTEGER :: decomp_db_gc_proc_row_group(max_decomps)
      INTEGER :: decomp_db_gc_proc_col_group(max_decomps)
      INTEGER :: decomp_db_gc_all_proc_group(max_decomps)

      ! indicates if a decomposition has been initialised

      LOGICAL :: decomp_db_set(max_decomps)

      COMMON /DECOMP_DATABASE/                                          &
     &  decomp_db_bound , decomp_db_sb_model_domain , decomp_db_glsize  &
     & ,decomp_db_g_lasize , decomp_db_gridsize,                        &
     &  decomp_db_g_blsize,                                             &
     &  decomp_db_g_datastart , decomp_db_g_datastart_f,                &
     &  decomp_db_g_gridpos,                                            &
     &  decomp_db_g_pe_index_EW,decomp_db_g_pe_index_NS,                &
     &  decomp_db_halosize , decomp_db_neighbour,                       &
     &  decomp_db_first_comp_pe , decomp_db_last_comp_pe,               &
     &  decomp_db_nproc,                                                &
     &  decomp_db_gc_proc_row_group , decomp_db_gc_proc_col_group,      &
     &  decomp_db_gc_all_proc_group,                                    &
     &  decomp_db_set

! End of DECOMPDB comdeck

! Local variables
      INTEGER ineb,idim,ifld,iproc,ihalo,ipt,iside


! ------------------------------------------------------------------


! Check that the new_decomp argument is sensible
      IF ((new_decomp  >   max_decomps) .OR.                            &
     &   ((new_decomp  <   1) .AND. (new_decomp  /=  decomp_unset)))    &
     & THEN
        IF (mype  ==  0) THEN
          WRITE(6,*) 'Error: Cannot change to decomposition ',          &
     &               new_decomp
          WRITE(6,*) 'This decomposition does not exist'
          WRITE(6,*) 'Exiting.'
        ENDIF
        icode=1
        GOTO 999
      ENDIF

! Check if this is already the current decomposition

      IF (new_decomp  ==  current_decomp_type) GOTO 999

! Check to see if setting decomposition to unset

      IF (new_decomp  ==  decomp_unset) THEN
        current_decomp_type = decomp_unset
        GOTO 999
      ENDIF

! Check if this decomposition has been initialised

      IF ( .NOT. decomp_db_set(new_decomp) ) THEN
        IF (mype  ==  0) THEN
          WRITE(6,*) 'Error : Attempt to select uninitialised ',        &
     &               'decomposition ',new_decomp
          WRITE(6,*) 'Exiting.'
        ENDIF
        icode=1
        GOTO 999
      ENDIF

! Now we can copy the information into PARVARS

      first_comp_pe=decomp_db_first_comp_pe(new_decomp)
      last_comp_pe=decomp_db_last_comp_pe(new_decomp)

      nproc=decomp_db_nproc(new_decomp)
      nproc_x=decomp_db_gridsize(1,new_decomp)
      nproc_y=decomp_db_gridsize(2,new_decomp)

      sb_model_domain=decomp_db_sb_model_domain(new_decomp)

      DO ihalo=1,NHalo_max
        DO idim=1,Ndim_max
          halosize(idim,ihalo)=                                         &
     &      decomp_db_halosize(idim,ihalo,new_decomp)
        ENDDO
      ENDDO

      Offx=decomp_db_halosize(1,halo_type_single,new_decomp)
      Offy=decomp_db_halosize(2,halo_type_single,new_decomp)

      halo_i=decomp_db_halosize(1,halo_type_extended,new_decomp)
      halo_j=decomp_db_halosize(2,halo_type_extended,new_decomp)

      gc_proc_row_group=decomp_db_gc_proc_row_group(new_decomp)
      gc_proc_col_group=decomp_db_gc_proc_col_group(new_decomp)
      gc_all_proc_group=decomp_db_gc_all_proc_group(new_decomp)

      DO ineb=1,4
        neighbour(ineb)=decomp_db_neighbour(ineb,new_decomp)
      ENDDO

      DO idim=1,Ndim_max
        bound(idim)=decomp_db_bound(idim,new_decomp)
        gridsize(idim)=decomp_db_gridsize(idim,new_decomp)

        DO ifld=1,Nfld_max
          glsize(idim,ifld)=decomp_db_glsize(idim,ifld,new_decomp)
        ENDDO
      ENDDO

      DO iproc=first_comp_pe,last_comp_pe
        DO idim=1,Ndim_max
          g_datastart(idim,iproc)=                                      &
     &      decomp_db_g_datastart(idim,iproc,new_decomp)
          g_gridpos(idim,iproc)=                                        &
     &      decomp_db_g_gridpos(idim,iproc,new_decomp)
        ENDDO
        g_at_extremity(PNorth,iproc)=                                   &
     &   (g_gridpos(2,iproc)  ==  (gridsize(2)-1))
        g_at_extremity(PSouth,iproc)=(g_gridpos(2,iproc)  ==  0)
        g_at_extremity(PEast,iproc)=                                    &
     &   (g_gridpos(1,iproc)  ==  (gridsize(1)-1))
        g_at_extremity(PWest,iproc)=(g_gridpos(1,iproc)  ==  0)

        DO ihalo=1,NHalo_max
          DO ifld=1,Nfld_max
            DO idim=1,Ndim_max
              g_lasize(idim,ifld,ihalo,iproc)=                          &
     &          decomp_db_g_lasize(idim,ifld,ihalo,iproc,new_decomp)
            ENDDO
          ENDDO
        ENDDO

        DO ifld=1,Nfld_max
          DO idim=1,Ndim_max
            g_blsize(idim,ifld,iproc)=                                  &
     &          decomp_db_g_blsize(idim,ifld,iproc,new_decomp)

            g_datastart_f(idim,ifld,iproc)=                             &
     &          decomp_db_g_datastart_f(idim,ifld,iproc,new_decomp)
          ENDDO
        ENDDO
      ENDDO

      DO idim=1,Ndim_max
        DO ifld=1,Nfld_max

          DO ihalo=1,NHalo_max
            lasize(idim,ifld,ihalo)=                                    &
     &        g_lasize(idim,ifld,ihalo,mype)
          ENDDO ! ihalo

          blsize(idim,ifld)=g_blsize(idim,ifld,mype)
          datastart_f(idim,ifld)=g_datastart_f(idim,ifld,mype)

        ENDDO ! ifld

        datastart(idim)=g_datastart(idim,mype)
        gridpos(idim)=g_gridpos(idim,mype)

      ENDDO ! idim


      DO iside=1,4
        at_extremity(iside)=g_at_extremity(iside,mype)
      ENDDO

      DO ipt=1-decomp_db_halosize(1,halo_type_extended,new_decomp),     &
     &       decomp_db_glsize(1,fld_type_p,new_decomp)+                 &
     &         decomp_db_halosize(1,halo_type_extended,new_decomp)
        g_pe_index_EW(ipt)=decomp_db_g_pe_index_EW(ipt,new_decomp)
      ENDDO

      DO ipt=1-decomp_db_halosize(2,halo_type_extended,new_decomp),     &
     &       decomp_db_glsize(2,fld_type_p,new_decomp)+                 &
     &         decomp_db_halosize(2,halo_type_extended,new_decomp)
        g_pe_index_NS(ipt)=decomp_db_g_pe_index_NS(ipt,new_decomp)
      ENDDO

      current_decomp_type=new_decomp

 999  CONTINUE

      RETURN
      END SUBROUTINE CHANGE_DECOMPOSITION

