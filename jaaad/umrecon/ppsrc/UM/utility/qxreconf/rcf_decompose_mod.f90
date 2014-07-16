
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Decomposition of grid for reconfiguration

Module Rcf_Decompose_Mod

!  Subroutine Rcf_Decompose - performs grid decomposition.
!
! Description:
! This routine performs a 2D decomposition - taking the global X and Y
! data sizes and decomposing across nproc_EW processors in the X
! direction and nproc_NS processors  in the Y direction.
!
! Method:
!   Based on UM4.5 code but now uses F90 data-structures.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Decompose( grid, nproc_EW, nproc_NS, decomp )

Use Rcf_Parvars_Mod      ! Most of this is used
Use Rcf_DecompTP_Mod     ! Most of this is used
Use Rcf_DecompDB_Mod     ! Most of this is used

Use Rcf_Lsm_Mod, Only : &
    Glob_LSM_Out,       &
    Glob_LSM_In,        &
    Local_LSM_Out,      &
    Local_LSM_In

Use Rcf_Set_Neighbour_Mod, Only : &
    Rcf_Set_Neighbour

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Use Rcf_Recon_Mod, Only : &
    GRIB

Implicit None

! Arguments
Type( grid_type ), Intent( InOut )  :: grid
Integer, Intent( In )               :: nproc_EW
Integer, Intent( In )               :: nproc_NS
Integer, Intent( In )               :: decomp

! Comdecks
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! GC - General Communication primitives package. For use on
! multiprocessor shared memory and message passing systems.
!
!
! LICENSING TERMS
!
!  GC is provided free of charge. Unless otherwise agreed with SINTEF,
!  use and redistribution in source and binary forms are permitted
!  provided that
!
!      (1) source distributions retain all comments appearing within
!          this file header, and
!
!      (2) distributions including binaries display the following
!          acknowledgement:
!
!              "This product includes software developed by SINTEF.",
!
!          in the documentation or other materials provided with the
!          distribution and in all advertising materials mentioning
!          features or use of this software.
!
!  The name of SINTEF may not be used to endorse or promote products
!  derived from this software without specific prior written
!  permission.  SINTEF disclaims any warranty that this software will
!  be fit for any specific purposes. In no event shall SINTEF be liable
!  for any loss of performance or for indirect or consequential damage
!  or direct or indirect injury of any kind. In no case shall SINTEF
!  be liable for any representation or warranty make to any third party
!  by the users of this software.
!
!
! Fortran header file. PLEASE use the parameter variables in user
! routines calling GC and NOT the numeric values. The latter are
! subject to change without further notice.
!
!---------------------------------------------- ------------------------
! $Id: gps0h501,v 1.6 2000/04/17 10:05:47 t11ps Exp $
! (C) Jorn Amundsen, Roar Skaalin, SINTEF Industrial Mathematics.

!    4.4   30/09/97  Added code to permit the SHMEM/NAM timeout
!                    value to be set from a shell variable.
!                      Author: Bob Carruthers  Cray Research.
!    5.1   17/04/00  Fixed/Free format. P.Selwood.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!     GC general options
      INTEGER, PARAMETER :: GC_OK         =     0
      INTEGER, PARAMETER :: GC_FAIL       =    -1
      INTEGER, PARAMETER :: GC_NONE       =     0
      INTEGER, PARAMETER :: GC_ANY        =    -1
      INTEGER, PARAMETER :: GC_DONTCARE   =    -1
      INTEGER, PARAMETER :: GC_SHM_DIR    =     1
      INTEGER, PARAMETER :: GC_SHM_SAFE   =     2
      INTEGER, PARAMETER :: GC_NAM_TIMEOUT=     4
      INTEGER, PARAMETER :: GC_SHM_GET    = -9999
      INTEGER, PARAMETER :: GC_SHM_PUT    = -9998
      INTEGER, PARAMETER :: GC_USE_GET    = -9999
      INTEGER, PARAMETER :: GC_USE_PUT    = -9998

!     GC functions
      INTEGER GC_COMLEN, GC_ISIZE, GC_RSIZE, GC_ME, GC_NPROC

!     GC groups (GCG) support
      INTEGER, PARAMETER :: GC_ALLGROUP = 0
      INTEGER, PARAMETER :: GCG_ALL = GC_ALLGROUP

!     GC groups (GCG) functions
      INTEGER GCG_ME

!     GC reserved message tags
      INTEGER, PARAMETER :: GC_MTAG_LOW   = 999999901
      INTEGER, PARAMETER :: GC_MTAG_HIGH  = 999999999

!     GCG_RALLETOALLE index parameters
      INTEGER, PARAMETER :: S_DESTINATION_PE = 1
      INTEGER, PARAMETER :: S_BASE_ADDRESS_IN_SEND_ARRAY = 2
      INTEGER, PARAMETER :: S_NUMBER_OF_ELEMENTS_IN_ITEM = 3
      INTEGER, PARAMETER :: S_STRIDE_IN_SEND_ARRAY = 4
      INTEGER, PARAMETER :: S_ELEMENT_LENGTH = 5
      INTEGER, PARAMETER :: S_BASE_ADDRESS_IN_RECV_ARRAY = 6
      INTEGER, PARAMETER :: S_STRIDE_IN_RECV_ARRAY = 7

      INTEGER, PARAMETER :: R_SOURCE_PE = 1
      INTEGER, PARAMETER :: R_BASE_ADDRESS_IN_RECV_ARRAY = 2
      INTEGER, PARAMETER :: R_NUMBER_OF_ELEMENTS_IN_ITEM = 3
      INTEGER, PARAMETER :: R_STRIDE_IN_RECV_ARRAY = 4
      INTEGER, PARAMETER :: R_ELEMENT_LENGTH = 5
      INTEGER, PARAMETER :: R_BASE_ADDRESS_IN_SEND_ARRAY = 6
      INTEGER, PARAMETER :: R_STRIDE_IN_SEND_ARRAY = 7

! Local variables
Integer                     :: iproc
Integer                     :: irest
Integer                     :: jrest
Integer                     :: info
Integer                     :: in_atm_decomp


! ------------------------------------------------------------------
decomp_DB (decomp) % halosize (1) = 0
decomp_DB (decomp) % halosize (2) = 0
decomp_DB (decomp) % halosize (3) = 0

! ------------------------------------------------------------------
! 1.0 Set up global data size
! ------------------------------------------------------------------

decomp_DB (decomp) % glsize (1) = grid % glob_p_row_length
decomp_DB (decomp) % glsize (2) = grid % glob_p_rows
decomp_DB (decomp) % glsize (3) = grid % model_levels

decomp_DB (decomp) % glsizeu (1) = grid % glob_u_row_length
decomp_DB (decomp) % glsizeu (2) = grid % glob_u_rows
decomp_DB (decomp) % glsizeu (3) = grid % model_levels

decomp_DB (decomp) % glsizev (1) = grid % glob_v_row_length
decomp_DB (decomp) % glsizev (2) = grid % glob_v_rows
decomp_DB (decomp) % glsizev (3) = grid % model_levels

decomp_DB (decomp) % glsizer (1) = grid % glob_r_row_length
decomp_DB (decomp) % glsizer (2) = grid % glob_r_rows
decomp_DB (decomp) % glsizer (3) = 1

! ------------------------------------------------------------------
! 2.0 Calculate decomposition
! ------------------------------------------------------------------


! select processors to use for the data decomposition
decomp_DB (decomp) % nproc          = nproc_EW * nproc_NS
decomp_DB (decomp) % first_comp_pe  = 0
decomp_DB (decomp) % last_comp_pe   =  decomp_DB (decomp) % nproc -1

!     Set the grid size

decomp_DB (decomp) % gridsize (1) = nproc_EW
decomp_DB (decomp) % gridsize (2) = nproc_NS
decomp_DB (decomp) % gridsize (3) = 1

! Calculate the local data shape of each processor.
DO iproc=decomp_DB (decomp) % first_comp_pe , &
         decomp_DB (decomp) % last_comp_pe
!       ! Loop over all processors in this decomposition

  decomp_DB (decomp) % g_gridpos (3,iproc) = 0

  decomp_DB (decomp) % g_gridpos (2,iproc) = &
    iproc / decomp_DB (decomp) % gridsize (1)

  decomp_DB (decomp) % g_gridpos (1,iproc) = &
    iproc - decomp_DB (decomp) % g_gridpos (2,iproc)* &
            decomp_DB (decomp) % gridsize (1)

!  Find the number of grid points (blsizep) in each domain and starting
!  points in the global domain (datastart) We first try to divide
!  the total number equally among the processors. The rest is
!  distributed one by one to first processor in each direction.

! The X (East-West) direction:

  decomp_DB (decomp) % g_blsizep (1,iproc) = &
    decomp_DB (decomp) % glsize (1) / &
    decomp_DB (decomp) % gridsize (1)
  irest = decomp_DB (decomp) % glsize (1)- &
          decomp_DB (decomp) % g_blsizep (1,iproc)* &
          decomp_DB (decomp) % gridsize (1)
  decomp_DB (decomp) % g_datastart (1,iproc) = &
    decomp_DB (decomp) % g_gridpos (1,iproc)* &
    decomp_DB (decomp) % g_blsizep (1,iproc) + 1

  IF (decomp_DB (decomp) % g_gridpos (1,iproc) .LT. &
      irest) THEN
    decomp_DB (decomp) % g_blsizep (1,iproc) = &
      decomp_DB (decomp) % g_blsizep (1,iproc)+1
    decomp_DB (decomp) % g_datastart (1,iproc) = &
      decomp_DB (decomp) % g_datastart (1,iproc) + &
      decomp_DB (decomp) % g_gridpos (1,iproc)
  ELSE
    decomp_DB (decomp) % g_datastart (1,iproc) = &
      decomp_DB (decomp) % g_datastart (1,iproc) + &
      irest
  ENDIF

  decomp_DB (decomp) % g_lasize (1,iproc)= &
    decomp_DB (decomp) % g_blsizep (1,iproc) + &
    2*decomp_DB (decomp) % halosize (1)

! X direction for river routing

  decomp_DB (decomp) % g_blsizer (1,iproc) =        &
    decomp_DB (decomp) % glsizer (1) /              &
    decomp_DB (decomp) % gridsize (1)
  irest = decomp_DB (decomp) % glsizer (1) -        &
          decomp_DB (decomp) % g_blsizer (1,iproc)* &
          decomp_DB (decomp) % gridsize (1)
  decomp_DB (decomp) % g_datastartr (1,iproc) =     &
    decomp_DB (decomp) % g_gridpos (1,iproc) *      &
    decomp_DB (decomp) % g_blsizer (1,iproc) + 1

  IF (decomp_DB (decomp) % g_gridpos (1,iproc) .LT. &
      irest) THEN
    decomp_DB (decomp) % g_blsizer (1,iproc) =      &
      decomp_DB (decomp) % g_blsizer (1,iproc) + 1

    decomp_DB (decomp) % g_datastartr (1,iproc) =   &
      decomp_DB (decomp) % g_datastartr (1,iproc) + &
      decomp_DB (decomp) % g_gridpos (1,iproc)
  ELSE
    decomp_DB (decomp) % g_datastartr (1,iproc) =   &
      decomp_DB (decomp) % g_datastartr (1,iproc) + &
      irest
  ENDIF

! The Y (North-South) direction

  decomp_DB (decomp) % g_blsizep (2,iproc) = &
    decomp_DB (decomp) % glsize (2) / &
    decomp_DB (decomp) % gridsize (2)
  jrest = decomp_DB (decomp) % glsize (2)- &
          decomp_DB (decomp) % g_blsizep (2,iproc)* &
          decomp_DB (decomp) % gridsize (2)
  decomp_DB (decomp) % g_datastart (2,iproc) = &
    decomp_DB (decomp) % g_gridpos (2,iproc)* &
    decomp_DB (decomp) % g_blsizep (2,iproc) + 1

  IF (decomp_DB (decomp) % g_gridpos (2,iproc) .LT. &
      jrest) THEN
    decomp_DB (decomp) % g_blsizep (2,iproc) = &
      decomp_DB (decomp) % g_blsizep (2,iproc)+1
    decomp_DB (decomp) % g_datastart (2,iproc) = &
      decomp_DB (decomp) % g_datastart (2,iproc) + &
      decomp_DB (decomp) % g_gridpos (2,iproc)
  ELSE
    decomp_DB (decomp) % g_datastart (2,iproc) = &
      decomp_DB (decomp) % g_datastart (2,iproc) +  jrest
  ENDIF

  decomp_DB (decomp) % g_lasize (2,iproc)= &
    decomp_DB (decomp) % g_blsizep (2,iproc) + &
    2*decomp_DB (decomp) % halosize (2)

! Y decomp for rivers

  decomp_DB (decomp) % g_blsizer (2,iproc) =        &
    decomp_DB (decomp) % glsizer (2) /              &
    decomp_DB (decomp) % gridsize (2)
  jrest = decomp_DB (decomp) % glsizer (2)-         &
          decomp_DB (decomp) % g_blsizer (2,iproc)* &
          decomp_DB (decomp) % gridsize (2)
  decomp_DB (decomp) % g_datastartr (2,iproc) =     &
    decomp_DB (decomp) % g_gridpos (2,iproc) *      &
    decomp_DB (decomp) % g_blsizer (2,iproc) + 1

  IF (decomp_DB (decomp) % g_gridpos (2,iproc) .LT. &
      jrest) THEN
    decomp_DB (decomp) % g_blsizer (2,iproc) =      &
      decomp_DB (decomp) % g_blsizer (2,iproc) + 1

    decomp_DB (decomp) % g_datastartr (2,iproc) =   &
      decomp_DB (decomp) % g_datastartr (2,iproc) + &
      decomp_DB (decomp) % g_gridpos (2,iproc)
  ELSE
    decomp_DB (decomp) % g_datastartr (2,iproc) =   &
      decomp_DB (decomp) % g_datastartr (2,iproc) +  jrest
  ENDIF

!  The Z (vertical) direction (no decomposition):

  decomp_DB (decomp) % g_datastart (3,iproc) = 1
  decomp_DB (decomp) % g_blsizep (3,iproc)   =  grid % model_levels
  decomp_DB (decomp) % g_lasize (3,iproc)    =  grid % model_levels
  decomp_DB (decomp) % g_blsizer(3,iproc)    = 1

!  Commented out stuff only required for "correct" sizing of non-wrap
!  One less u value at the east of a LAM if non-wrapping
!  If (( grid % glob_u_row_length ==                   &
!        grid % glob_p_row_length - 1 ).AND.           &
!      decomp_DB (decomp) % g_gridpos (1,iproc) == &
!      decomp_DB (decomp) % gridsize (1) - 1 ) Then
!    decomp_DB (decomp) % g_blsizeu (1,iproc) = &
!        decomp_DB (decomp) % g_blsizep (1,iproc) - 1
!  Else
    decomp_DB (decomp) % g_blsizeu (1,iproc) = &
        decomp_DB (decomp) % g_blsizep (1,iproc)
!  End If

  decomp_DB (decomp) % g_blsizeu (2,iproc) = &
    decomp_DB (decomp) % g_blsizep (2,iproc)

  decomp_DB (decomp) % g_blsizeu (3,iproc) = &
    decomp_DB (decomp) % g_blsizep (3,iproc)

! One less v row at the north
   decomp_DB (decomp) % g_blsizev (1,iproc) = &
     decomp_DB (decomp) % g_blsizep (1,iproc)

  If ( (decomp_DB (decomp) % g_gridpos (2,iproc) == &
       decomp_DB (decomp) % gridsize (2) - 1) .AND. &
     grid % glob_p_rows /=  grid % glob_v_rows ) Then
    decomp_DB (decomp) % g_blsizev (2,iproc) = &
        decomp_DB (decomp) % g_blsizep (2,iproc) -1
  Else
    decomp_DB (decomp) % g_blsizev (2,iproc) = &
        decomp_DB (decomp) % g_blsizep (2,iproc)
  End If

  decomp_DB (decomp) % g_blsizev (3,iproc) = &
    decomp_DB (decomp) % g_blsizep (3,iproc)

ENDDO ! loop over processors


! ------------------------------------------------------------------
! 3.0 Set boundary conditions
! ------------------------------------------------------------------

! if global or wrapping LAM
If ( grid % glob_p_row_length ==  grid % glob_u_row_length ) Then
  decomp_DB (decomp) % bound (1) = BC_CYCLIC ! Cyclic East-West bdy
Else
  decomp_DB (decomp) % bound (1) = BC_STATIC ! No East-West wrap around
Endif

If ( grid % glob_p_rows ==  grid % glob_v_rows .AND. .NOT.GRIB ) Then
  decomp_DB (decomp) % bound (2) = BC_CYCLIC ! cyclic N-S
Else
  decomp_DB (decomp) % bound (2) = BC_STATIC !No North-South wrap around
Endif
decomp_DB (decomp) % bound (3) = BC_STATIC ! No vertical wrap around

CALL Rcf_Set_Neighbour (decomp)

! ------------------------------------------------------------------
! 4.0 Return the new data sizes and exit subroutine
! ------------------------------------------------------------------

! Set up the GCOM groups:

! 1) Group of all processors on my row

IF ( decomp_DB (decomp) % gridsize (2) .EQ. 1) THEN
 decomp_DB (decomp) % gc_proc_row_group =GCG_ALL
ELSE
  CALL GCG_SPLIT(mype,nproc_max, &
    decomp_DB (decomp) % g_gridpos (2,mype), info, &
    decomp_DB (decomp) % gc_proc_row_group )
ENDIF

! 2) Group of all processors on my column

IF ( decomp_DB (decomp) % gridsize (1) .EQ. 1)  THEN
  decomp_DB (decomp) % gc_proc_col_group =GCG_ALL
ELSE
  CALL GCG_SPLIT(mype,nproc_max, &
    decomp_DB (decomp) % g_gridpos (1,mype), info, &
    decomp_DB (decomp) % gc_proc_col_group )
ENDIF

! 3) Group of all processors in the atmosphere model
IF (decomp_DB (decomp) % nproc  .EQ. nproc_max) &
 THEN
  decomp_DB (decomp) % gc_all_proc_group =GCG_ALL
ELSE
  IF ((mype .GE. decomp_DB (decomp) % first_comp_pe ) .AND. &
     (mype .LE. decomp_DB (decomp) % last_comp_pe ) ) &
   THEN
    in_atm_decomp=1
  ELSE
    in_atm_decomp=0
  ENDIF

  CALL GCG_SPLIT(mype,nproc_max,in_atm_decomp,info, &
    decomp_DB (decomp) % gc_all_proc_group )
ENDIF

! Set logical indicating this decomposition has been initialised
! and is now ready for use

decomp_DB (decomp) % set =.TRUE.

! And return the new horizontal dimensions

grid % loc_p_row_length = decomp_DB (decomp) % g_blsizep (1,mype)
grid % loc_p_rows       = decomp_DB (decomp) % g_blsizep (2,mype)
grid % loc_p_field      = grid % loc_p_row_length * grid % loc_p_rows
grid % loc_u_row_length = decomp_DB (decomp) % g_blsizeu (1,mype)
grid % loc_u_rows       = decomp_DB (decomp) % g_blsizeu (2,mype)
grid % loc_u_field      = grid % loc_u_row_length * grid % loc_u_rows
grid % loc_v_row_length = decomp_DB (decomp) % g_blsizev (1,mype)
grid % loc_v_rows       = decomp_DB (decomp) % g_blsizev (2,mype)
grid % loc_v_field      = grid % loc_v_row_length * grid % loc_v_rows
grid % loc_r_row_length = decomp_DB (decomp) % g_blsizer (1,mype)
grid % loc_r_rows       = decomp_DB (decomp) % g_blsizer (2,mype)
grid % loc_r_field      = grid % loc_r_row_length * grid % loc_r_rows

! Finally allocate space for land-sea-mask
If (decomp == decomp_rcf_output) Then
  Allocate( glob_lsm_out(grid % glob_p_row_length * grid % glob_p_rows))
  Allocate( local_lsm_out( grid % loc_p_row_length * grid % loc_p_rows))

Else If (decomp == decomp_rcf_input) Then
  Allocate( glob_lsm_in(grid % glob_p_row_length * grid % glob_p_rows))
  Allocate( local_lsm_in( grid % loc_p_row_length * grid % loc_p_rows))
End If


Return
End Subroutine Rcf_Decompose
End Module Rcf_Decompose_Mod
