
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine TURB_Smagorinsky
!
          Subroutine TURB_Smagorinsky(                                  &
     &                     rows, row_length, n_rows                     &
     &,                    model_levels                                 &
     &,                    r_theta_levels, r_rho_levels                 &
     &,                    u, v, w, visc, shear                         &
     &,                    z0, RNEUTML, timestep                        &
     &,                    diff_factor, mix_factor, max_diff            &
     &,                    cos_theta_latitude                           &
     &,                    cos_v_latitude                               &
     &,                    delta_lambda, delta_phi)

! Description: Calculates coefficients for use in subgrid turbulence
!              scheme
!
! Method:
!
! Current code owner: Carol Halliwell
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 6.2      05/01/06  New code
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
      Implicit None

!*L------------------COMDECK C_A----------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Replace variable A by more meaningful name for
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Convert to Fixed/Free format. P. Selwood

      ! Mean radius of Earth in m.
      Real, Parameter  :: Earth_Radius = 6371000.

!*----------------------------------------------------------------------
! for earth radius
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
! halo information
! C_VKMAN start
      REAL,PARAMETER:: VKMAN=0.4 ! Von Karman's constant
! C_VKMAN end
!Von Karman constant

! Variables with Intent (In)

      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of v rows.
     &, model_levels     ! number of model levels.

      Real                                                              &
     &  u(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)       &
     &, v(1-offx:row_length+offx, 1-offy:n_rows+offy, model_levels)     &
     &, w(1-offx:row_length+offx, 1-offy:rows+offy, 0:model_levels)     &
     &, z0(row_length,rows)                                             &
                             ! roughness length
     &, timestep                                                        &
     &, diff_factor                                                     &
                    ! factor between 0 and 1 multiplied by max diff
!                   ! coeff for numerical stability to obtain max
!                   ! diffusion coeff for run
     &, max_diff                                                        &
                    ! max diffusion coeff for run
     &, mix_factor                                                      &
                    ! lambda0 = mix_factor * gridlength
     &, r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                    1-halo_j:rows+halo_j, 0:model_levels)         &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                    1-halo_j:rows+halo_j, model_levels)           &
     &, delta_lambda                                                    &
     &, delta_phi                                                       &
     &, cos_theta_latitude (1-Offx:row_length+Offx,                     &
     &                                           1-Offy:rows+Offy)      &
     &, cos_v_latitude (1-Offx:row_length+Offx,1-Offy:n_rows+Offy)

!Local parameters

      Integer                                                           &
     & i,j,k

      Real                                                              &
     &  delta_x, delta_y                                                &
     &, delta_z(row_length, rows, model_levels)                         &
     &, delta_zn(row_length+offx, rows+offy, model_levels)              &
     &, delta_zn_u(row_length, rows, model_levels)                      &
     &, delta_zn_v(row_length, rows, model_levels)                      &
     &, rdz(row_length, rows, model_levels)                             &
     &, rdzn_u(1-offx:row_length+offx, rows, model_levels)              &
     &, rdzn_v(row_length, 1-offy:rows+offy, model_levels)              &
     &, smallp                                                          &
                   ! A small number
     &, RNEUTML(row_length, rows, model_levels)                         &
!               ! mixing length scale (m) (lambda)
     &, RNEUTML_SQ(row_length, rows, model_levels)                      &
!               ! square of RNEUTML
     &, RMLMAX                                                          &
!               ! basic mixing length (lambda0)
     &, z(row_length, rows, model_levels)                               &
     &, shear(row_length, rows, model_levels)

      Real                                                              &
     &  dx_rho(row_length, rows, model_levels)                          &
!               ! dx on rho points
     &, dx_theta_u(row_length, rows, model_levels)                      &
!               ! dx above u points on theta levels
     &, dx_rho_centre(row_length, rows, model_levels)                   &
!               ! dx in centre of grid box on rho levels
     &, dy_rho(row_length, rows, model_levels)                          &
!               ! dy on rho points
     &, dy_theta_v(row_length, rows, model_levels)                      &
!               ! dy above v points on theta levels
     &, dy_rho_centre(row_length, rows, model_levels)                   &
!               ! dy in centre of grid box on rho levels
     &, cx_rho(row_length, rows, model_levels)                          &
!               ! reciprocal of dx_rho
     &, cx_theta_u(1-offx:row_length+offx, rows, model_levels)          &
!               ! reciprocal of dx_theta_u
     &, cx_rho_centre(1-offx:row_length+offx, 1-offy:rows+offy          &
     &                                          , model_levels)         &
!               ! reciprocal of dx_rho_centre
     &, cy_rho(row_length, rows, model_levels)                          &
!               ! reciprocal of dy_rho
     &, cy_theta_v(row_length, 1-offy:rows+offy, model_levels)          &
!               ! reciprocal of dy_theta_v
     &, cy_rho_centre(1-offx:row_length+offx, 1-offy:rows+offy          &
     &                                          , model_levels)         &
!               ! reciprocal of dy_rho_centre
     &, cx2_rho(1-offx:row_length+offx, rows, model_levels)             &
!               ! square of cx_rho
     &, cy2_rho(row_length, 1-offy:rows+offy, model_levels)
!               ! square of cy_rho

      Parameter (smallp=1.e-14)

      Real                                                              &
     &  ssq11                                                           &
                      ! ii component of shear stress
     &, ssq22                                                           &
                      ! jj component of shear stress
     &, ssq33                                                           &
                      ! kk component of shear stress
     &, ssq13                                                           &
                      ! ik component of shear stress
     &, ssq23                                                           &
                      ! jk component of shear stress
     &, ssq12                                                           &
                      ! ij component of shear stress
     &, ssq(row_length, rows, model_levels)                             &
                                     ! Half squared strain rate
     &, sum(row_length, rows, model_levels)   ! sum of Sij (ssqij)

! Variables with Intent (Out)

      Real                                                              &
     &  visc(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j           &
     &                  , model_levels)   ! OUT: lambda^2*S

!----------------------------------------------------------------------
! Initialise variables
!----------------------------------------------------------------------
      Do k= 1, model_levels

        Do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
! viscosity is set to zero at the top of the domain.
            visc(i,j,k) = 0.0
          End Do
        End Do

        Do j = 1, rows
          Do i = 1, row_length
            sum(i,j,k) = 0.0
            ssq(i,j,k) = 0.0
          End Do
        End Do

      End Do

      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            z(i,j,k) = r_theta_levels(i,j,k)                            &
     &                        - r_theta_levels(i,j,0)
           End Do
         End Do
      End Do
!----------------------------------------------------------------------
! Calculate grid spacings
!----------------------------------------------------------------------
! Horizontal grid spacings used in calculation of rate of strain term
! delta_lambda and delta_phi are in radians
!
      delta_x = Earth_radius * delta_lambda
      delta_y = Earth_radius * delta_phi

      Do k = 1, model_levels

        Do j = 1, rows
          Do i = 1, row_length

            dx_rho(i,j,k) =                                             &
     &               r_rho_levels(i,j,k)                                &
     &             * delta_lambda * cos_theta_latitude(i,j)

            dx_theta_u(i,j,k) =                                         &
     &         0.5*(r_theta_levels(i,j,k) + r_theta_levels(i+1,j,k))    &
     &       * delta_lambda * cos_theta_latitude(i,j)

            dy_rho(i,j,k) =                                             &
     &               r_rho_levels(i,j,k) * delta_phi

            dy_theta_v(i,j,k) =                                         &
     &         0.5*(r_theta_levels(i,j,k) + r_theta_levels(i,j+1,k))    &
     &       * delta_phi

            dy_rho_centre(i,j,k) = 0.25*(                               &
     &         r_rho_levels(i,j+1,k) + r_rho_levels(i+1,j+1,k)          &
     &        +r_rho_levels(i,j,k) + r_rho_levels(i+1,j,k) )            &
     &        * delta_phi

            cx_rho(i,j,k) = 1./dx_rho(i,j,k)
            cx_theta_u(i,j,k) =  1./dx_theta_u(i,j,k)
            cy_rho(i,j,k) =  1./dy_rho(i,j,k)
            cy_theta_v(i,j,k) =  1./dy_theta_v(i,j,k)
            cy_rho_centre(i,j,k) = 1./dy_rho_centre(i,j,k)
            cx2_rho(i,j,k) = cx_rho(i,j,k)*cx_rho(i,j,k)
            cy2_rho(i,j,k) =  cy_rho(i,j,k)*cy_rho(i,j,k)
          End Do
        End Do

        Do j = 1, n_rows
          Do i = 1, row_length
            dx_rho_centre(i,j,k) = 0.25*(                               &
     &         r_rho_levels(i,j+1,k) + r_rho_levels(i+1,j+1,k)          &
     &        +r_rho_levels(i,j,k) + r_rho_levels(i+1,j,k) )            &
     &                * delta_lambda * cos_v_latitude(i,j)
            cx_rho_centre(i,j,k) =  1./dx_rho_centre(i,j,k)
          End Do
        End Do

      End Do   !k

      Do k = 1, model_levels
        Do i = 1, row_length
          dx_rho_centre(i,rows,k) = dx_rho_centre(i,n_rows,k)
          cx_rho_centre(i,rows,k) =  1./dx_rho_centre(i,n_rows,k)
        End Do
      End Do

! Vertical grid spacings used in calculation of rate of strain term

      Do k = 1, model_levels

        Do j = 1, rows
          Do i = 1, row_length
            delta_z(i,j,k) = r_theta_levels(i,j,k)                      &
     &                -r_theta_levels(i,j,k-1)
            rdz(i,j,k) = 1./delta_z(i,j,k)
          End Do
        End Do

        Do j = 1, rows+offy
          Do i = 1, row_length+offx
            If (k  ==  1) then
              delta_zn(i,j,k) = r_rho_levels(i,j,k)                     &
     &                         -r_theta_levels(i,j,0)
            Else
              delta_zn(i,j,k) = r_rho_levels(i,j,k)                     &
     &                        - r_rho_levels(i,j,k-1)
            End If
          End Do
        End Do

      End Do

      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            delta_zn_u(i,j,k) =                                         &
     &             0.5*(delta_zn(i,j,k)+delta_zn(i+1,j,k))
            delta_zn_v(i,j,k) =                                         &
     &             0.5*(delta_zn(i,j,k)+delta_zn(i,j+1,k))
            rdzn_u(i,j,k) = 1./delta_zn_u(i,j,k)
            rdzn_v(i,j,k) = 1./delta_zn_v(i,j,k)
          End Do
        End Do
      End Do

! DEPENDS ON: swap_bounds
      call Swap_Bounds(cx2_rho, row_length, rows, model_levels,         &
     &                 offx, 0, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
      call Swap_Bounds(cx_theta_u, row_length, rows, model_levels,      &
     &                 offx, 0, fld_type_u, .false.)

! DEPENDS ON: swap_bounds
      call Swap_Bounds(cx_rho_centre, row_length, rows, model_levels,   &
     &                 offx, offy, fld_type_v, .false.)

! DEPENDS ON: swap_bounds
      call Swap_Bounds(cy2_rho, row_length, rows, model_levels,         &
     &                 0, offy, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
      call Swap_Bounds(cy_theta_v, row_length, rows, model_levels,      &
     &                 0, offy, fld_type_v, .false.)

! DEPENDS ON: swap_bounds
      call Swap_Bounds(cy_rho_centre, row_length, rows, model_levels,   &
     &                 offx, offy, fld_type_v, .false.)

! DEPENDS ON: swap_bounds
      call Swap_Bounds(rdzn_u, row_length, rows, model_levels,          &
     &                 offx, 0, fld_type_u, .false.)

! DEPENDS ON: swap_bounds
      call Swap_Bounds(rdzn_v, row_length, rows, model_levels,          &
     &                 0, offy, fld_type_v, .false.)

!--------------------------------------------------------
! As in the LEM code
! _Now calculate half-squared strain rate SSQ on w-points
! _CX=1./DX, CY=1./DY, RDZ(K)=1./DZ(K),
!       RDZN(K) =1./DZN(K)
! _SSQ= 0.5*^^DU_I/DX_J+DU_J/DX_I^^**2
! _SSQIJ= (DU_I/DX_J+DU_J/DX_I)**2
! _Hence SSQ= SUM(I,J) {0.5*(SSQIJ)}
! _Note that for a simple shear S, SSQ=S**2
!   (as in the notation of Mason and Callen 1986)
!--------------------------------------------------------

      Do k = 1, model_levels - 1

        Do j = 1, rows

          Do i = 1, row_length

          ssq11 =                                                       &
     &     cx2_rho(i,j,k+1)*(u(i,j,k+1)-u(i-1,j,k+1))**2 +              &
     &     cx2_rho(i,j,k)*(u(i,j,k)-u(i-1,j,k))**2
          ssq22 =                                                       &
     &     cy2_rho(i,j,k+1)*(v(i,j,k+1)-v(i,j-1,k+1))**2+               &
     &     cy2_rho(i,j,k)*(v(i,j,k)-v(i,j-1,k))**2
          ssq33 =                                                       &
     &       ((w(i,j,k)-w(i,j,k-1))*rdz(i,j,k))**2 +                    &
     &       ((w(i,j,k+1)-w(i,j,k))*rdz(i,j,k+1))**2
          ssq13= (                                                      &
     &       ((u(i,j,k+1)-u(i,j,k))*rdzn_u(i,j,k+1)+                    &
     &       (w(i+1,j,k)-w(i,j,k))*cx_theta_u(i,j,k))**2+               &
     &       ((u(i-1,j,k+1)-u(i-1,j,k))*rdzn_u(i-1,j,k+1)               &
     &       +(w(i,j,k)-w(i-1,j,k))*cx_theta_u(i-1,j,k))**2             &
     &             )*0.5      ! _averaging ssq13 over 2 points
          ssq23 = (                                                     &
     &       ((w(i,j,k)-w(i,j-1,k))*cy_theta_v(i,j-1,k)                 &
     &      +(v(i,j-1,k+1)-v(i,j-1,k))*rdzn_v(i,j-1,k+1))**2+           &
     &       ((w(i,j+1,k)-w(i,j,k))*cy_theta_v(i,j,k)+                  &
     &       (v(i,j,k+1)-v(i,j,k))*rdzn_v(i,j,k+1))**2                  &
     &            )*0.5        ! _averaging ssq23 over 2 points
          ssq12 = (  (  (                                               &
     &       ((u(i-1,j,k)-u(i-1,j-1,k))*cy_rho_centre(i-1,j-1,k)        &
     &      +(v(i,j-1,k)-v(i-1,j-1,k))*cx_rho_centre(i-1,j-1,k))**2 +   &
     &       ((u(i-1,j+1,k)-u(i-1,j,k))*cy_rho_centre(i-1,j,k)          &
     &      +(v(i,j,k)-v(i-1,j,k))*cx_rho_centre(i-1,j,k))**2           &
     &              )  +  (                                             &
     &       ((u(i,j,k)-u(i,j-1,k))*cy_rho_centre(i,j-1,k)              &
     &      +(v(i+1,j-1,k)-v(i,j-1,k))*cx_rho_centre(i,j-1,k))**2 +     &
     &       ((u(i,j+1,k)-u(i,j,k))*cy_rho_centre(i,j,k)                &
     &      +(v(i+1,j,k)-v(i,j,k))*cx_rho_centre(i,j,k))**2             &
     &        ) )  +  ( (                                               &
     &       ((u(i-1,j,k+1)-u(i-1,j-1,k+1))*cy_rho_centre(i-1,j-1,k+1)  &
     &      +(v(i,j-1,k+1)-v(i-1,j-1,k+1))                              &
     &        *cx_rho_centre(i-1,j-1,k+1))**2+                          &
     &       ((u(i-1,j+1,k+1)-u(i-1,j,k+1))*cy_rho_centre(i-1,j,k+1)    &
     &      +(v(i,j,k+1)-v(i-1,j,k+1))*cx_rho_centre(i-1,j,k+1))**2     &
     &              )  +  (                                             &
     &       ((u(i,j,k+1)-u(i,j-1,k+1))*cy_rho_centre(i,j-1,k+1)        &
     &      +(v(i+1,j-1,k+1)-v(i,j-1,k+1))                              &
     &       *cx_rho_centre(i,j-1,k+1))**2+                             &
     &       ((u(i,j+1,k+1)-u(i,j,k+1))*cy_rho_centre(i,j,k+1)          &
     &      +(v(i+1,j,k+1)-v(i,j,k+1))*cx_rho_centre(i,j,k+1))**2       &
     &         ) ) )*0.125      ! _averaging ssq12 over 8 points

          ssq(i,j,k)=ssq11+ssq22+ssq33+ssq13+ssq23+ssq12+smallp
          sum(i,j,k) = sum(i,j,k) + ssq(i,j,k)

          End Do   ! on i

        End Do    ! on j

      End Do ! on k

      Do k = 1,model_levels-1
        Do j = 1, rows
          Do i = 1, row_length
           shear(i,j,k) = sqrt(sum(i,j,k)) ! already been *0.5 in ssq
          End Do
        End Do
      End Do

      RMLMAX = mix_factor * MAX(delta_x, delta_y)

      Do k = 1, model_levels-1
        Do j = 1, rows
          Do i = 1, row_length

            RNEUTML(i,j,K)=SQRT(1./                                     &
     &     (1./(VKMAN*(z(i,j,k)+z0(i,j)))**2+1./RMLMAX**2) )
            RNEUTML_SQ(i,j,K)=RNEUTML(i,j,K)*RNEUTML(i,j,K)
            visc(i,j,k) = RNEUTML_SQ(i,j,k)*shear(i,j,k)
!
! visc is now lambda^2*S and still needs to be multiplied by
! stability functions (done in ATMSTEP)
!
          End Do
        End Do
      End Do
!
! maximum diffusion coefficient allowed in this run
!
      max_diff = diff_factor*delta_x*delta_x/(8.0*timestep)

      return
      END SUBROUTINE TURB_Smagorinsky

