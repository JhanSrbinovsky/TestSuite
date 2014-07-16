#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C)
#if defined(MPP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Parallel UM: Perform 2D data decomposition
!
! Subroutine interface:
      SUBROUTINE DECOMPOSE_ATMOS(global_row_len, global_n_rows,         &
     &                           tot_levels,                            &
     &                           river_rows, river_row_length,          &
     &                           model_type,                            &
     &                           nproc_EW, nproc_NS,                    &
     &                           extended_halo_EW,                      &
     &                           extended_halo_NS,                      &
     &                           rimwidth, nrima_max,                   &
     &                           local_row_len, local_n_rows)
      IMPLICIT NONE
!
! Description:
! This routine performs a 2D decomposition - taking the global X
! (global_row_len) and Y (global_n_rows) data sizes and decomposing
! across nproc_EW processors in the X direction and nproc_NS processors
! in the Y direction.
! The local data size is returned via local_row_len and local_n_rows.
! These values will include a data halo for boundary updates.
!
! Method:
! The local data sizes are calculated and stored in the COMMON block
! DECOMPDB. The boundary conditions are set (cyclic in East/West
! direction if *DEF,GLOBAL).
!
! Current Code Owner: Paul Burton
!
! History:
!  Model    Date     Modification history from model version 3.5
!  version
!    3.5    1/3/95   New DECK created for the Parallel Unified
!                    Model. P.Burton + R.Skaalin
!    4.1    18/3/96  Added first/last_comp_pe variable.  P.Burton
!    4.2   19/08/96  Changed name to DECOMPOSE_ATMOS.
!                    Changed argument list to allow a standard
!                    interface to all decomposition routines.
!                    Changed decomposition description variables to
!                    the decomp_db* form, from the DECOMPDB comdeck
!                    to allow flexible decompositions.
!                    Added code to initialise GCOM groups
!                    Changed LAM model EW BCs to cyclic
!    5.0   12/04/99  - Check for valid nproc, decomposition
!                    - Set up extended halo sizes
!                    - Change decomposition algorithm for ND:
!                      - Even (or 1) processor EW
!                      - Symmetric distribution of excess points
!                    - Added model_type argument
!                                                       P.Burton
!    5.1   27/01/00  - Changed g_pe_index to g_pe_index_EW
!                      Added g_pe_index_NS
!                                                  P.Burton
!    5.3   10/09/01  Corrected loop bounds in calculation of
!                    decomp_db_g_pe_index_NS        P.Burton
!    5.3   14/09/01  Added model_domain variable    P.Burton
!    5.3   14/12/01  Added check on halo size       P.Burton
!    5.3   14/12/01  Output errors to PE0           P.Burton
!    5.3   27/07/01  Rimwidth added to perform an additional check
!                    that a sensible configuration has been chosen.
!                                                       Z. Gardner
! 5.3      15/09/01  add bi_cyclic_lam domain             A. Malcolm
!    5.5   15/01/03  River routing support. P.Selwood.
! 5.5      06/02/03  correct size check for mes           A. Malcolm
!    5.5   26/03/03  Correct error trap in N-S and problem when
!                    size=halo_ns                        A. Malcolm
!  6.0  17/09/03  Add def for new NEC opt section C96_1C. R Barnes
!  6.2    01/02/06  Fixed out of bounds reference.  T.Edwards
!  6.2   02/03/05  revise size check for mes            A. Malcolm
!
! Subroutine Arguments:

      INTEGER                                                           &

     &  global_row_len,                                                 &
                          ! IN  :number of E-W points of entire model
     &  global_n_rows,                                                  &
                          ! IN  :number of P rows of entire model
     &  tot_levels,                                                     &
                          ! IN  :total number of levels
     &  river_rows,                                                     &
                          ! IN  :number of rows in river routing model
     &  river_row_length,                                               &
                          ! IN  :number of E-W points for river model
     &  model_type,                                                     &
                          ! IN  : type (Global,LAM etc) of model
     &  nproc_EW,                                                       &
                          ! IN  : number of processors East-West
     &  nproc_NS,                                                       &
                          ! IN  : number of processors North-South
     &  extended_halo_EW,                                               &
                          ! IN  : size of extended EW halo
     &  extended_halo_NS,                                               &
                          ! IN  : size of extended NS halo
     &  nrima_max,                                                      &
                              ! IN : size of rimwidth
     &  rimwidth(nrima_max),                                            &
                             ! IN : size of blending region in the bcs
     &  local_row_len,                                                  &
                          ! OUT :number of E-W points of this process
     &  local_n_rows      ! OUT :number of rows of this process

! Parameters and Common blocks

#include "parvars.h"
#include "decomptp.h"
#include "decompdb.h"
#include "gccom.h"
#include "domtyp.h"

! Local variables
      INTEGER                                                           &
     &  iproc                                                           &
     &, iproc_x                                                         &
     &, iproc_y                                                         &
     &, ifld                                                            &
     &, ihalo                                                           &
     &, idim                                                            &
     &, ipt                                                             &
     &, irest                                                           &
     &, start                                                           &
     &, size1                                                           &
     &, size                                                            &
     &, prow_N                                                          &
     &, prow_S                                                          &
     &, info                                                            &
     &, in_atm_decomp                                                   &
     &, max_NS                                                          &
     &, max_EW


      LOGICAL                                                           &
     &  at_north

      INTEGER                                                           &
     &  size_x(0:nproc_EW-1)                                            &
     &, size_y(0:nproc_NS-1)                                            &
     &, start_x(0:nproc_EW-1)                                           &
     &, start_y(0:nproc_NS-1)

! For river routing
      INTEGER                                                           &
     &  rsize_x(0:nproc_EW-1)                                           &
     &, rsize_y(0:nproc_NS-1)                                           &
     &, rstart_x(0:nproc_EW-1)                                          &
     &, rstart_y(0:nproc_NS-1)

! Error reporting
      INTEGER       ErrorStatus ! =0 normal exit; >0 error exit
      CHARACTER*256 Cmessage    ! Error message
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='DECOMPOSE_ATMOS')


! ------------------------------------------------------------------
! 0.0 Check for valid decomposition
! ------------------------------------------------------------------

      IF (mype  ==  0) THEN
      IF (nproc_EW*nproc_NS  >   Maxproc) THEN
        ErrorStatus=1
        WRITE(Cmessage,                                                 &
     &    '("Cannot run with decomposition ",I3," x ",I3,               &
     &      " (",I3,") processors. ",                                   &
     &      "Maxproc is ",I3," processors.")')                          &
     &       nproc_EW,nproc_NS,nproc_EW*nproc_NS,Maxproc
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

      IF ((nproc_EW  /=  1) .AND. (MOD(nproc_EW,2)  /=  0)) THEN
        ErrorStatus=2
        WRITE(Cmessage,                                                 &
     &    '("Cannot run with an odd (",I3,") number of processors ",    &
     &      "in the East-West direction.")') nproc_EW
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

      IF (MOD(global_row_len,2)  /=  0) THEN
        ErrorStatus=3
        WRITE(Cmessage,                                                 &
     &    '("Cannot run with an odd (",I3,") number of points ",        &
     &    "in the East-West direction.")') global_row_len
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

      IF (extended_halo_EW  >   Max_Halo_Size) THEN
        ErrorStatus=4
        WRITE(Cmessage,                                                 &
     &    '("East-West extended halo size (",I2,") is too large.",      &
     &      "The maximum permitted size is Max_Halo_Size=",I2)')        &
     &    extended_halo_EW,Max_Halo_Size
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

      IF (extended_halo_NS  >   Max_Halo_Size) THEN
        ErrorStatus=4
        WRITE(Cmessage,                                                 &
     &    '("North-South extended halo size (",I2,") is too large.",    &
     &      "The maximum permitted size is Max_Halo_Size=",I2)')        &
     &    extended_halo_NS,Max_Halo_Size
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF
      ENDIF

! ------------------------------------------------------------------
      decomp_db_sb_model_domain(decomp_standard_atmos)=model_type
      decomp_db_halosize(1,halo_type_single,decomp_standard_atmos) = 1
      decomp_db_halosize(2,halo_type_single,decomp_standard_atmos) = 1
      decomp_db_halosize(3,halo_type_single,decomp_standard_atmos) = 0

      decomp_db_halosize(1,halo_type_extended,decomp_standard_atmos) =  &
     &  extended_halo_EW
      decomp_db_halosize(2,halo_type_extended,decomp_standard_atmos) =  &
     &  extended_halo_NS
      decomp_db_halosize(3,halo_type_extended,decomp_standard_atmos) =  &
     &  0

      decomp_db_halosize(1,halo_type_no_halo,decomp_standard_atmos) = 0
      decomp_db_halosize(2,halo_type_no_halo,decomp_standard_atmos) = 0
      decomp_db_halosize(3,halo_type_no_halo,decomp_standard_atmos) = 0


! ------------------------------------------------------------------
! 1.0 Set up global data size
! ------------------------------------------------------------------

      decomp_db_glsize(1,fld_type_p,decomp_standard_atmos) =            &
     &  global_row_len
      decomp_db_glsize(2,fld_type_p,decomp_standard_atmos) =            &
     &  global_n_rows
      decomp_db_glsize(3,fld_type_p,decomp_standard_atmos) =            &
     &  tot_levels

      decomp_db_glsize(1,fld_type_u,decomp_standard_atmos) =            &
     &  global_row_len
      decomp_db_glsize(2,fld_type_u,decomp_standard_atmos) =            &
     &  global_n_rows
      decomp_db_glsize(3,fld_type_u,decomp_standard_atmos) =            &
     &  tot_levels

      decomp_db_glsize(1,fld_type_v,decomp_standard_atmos) =            &
     &  global_row_len
      If (model_type  /=  mt_bi_cyclic_lam ) then
      decomp_db_glsize(2,fld_type_v,decomp_standard_atmos) =            &
     &  global_n_rows-1
      Else
        decomp_db_glsize(2,fld_type_v,decomp_standard_atmos) =          &
     &  global_n_rows
      Endif
      decomp_db_glsize(3,fld_type_v,decomp_standard_atmos) =            &
     &  tot_levels

      decomp_db_glsize(1,fld_type_r,decomp_standard_atmos) =            &
     &  river_row_length
      decomp_db_glsize(2,fld_type_r,decomp_standard_atmos) =            &
     &  river_rows
      decomp_db_glsize(3,fld_type_r,decomp_standard_atmos) =            &
     &  1

! ------------------------------------------------------------------
! 2.0 Calculate decomposition
! ------------------------------------------------------------------


! select processors to use for the data decomposition
      decomp_db_nproc(decomp_standard_atmos)=nproc_EW*nproc_NS
      decomp_db_first_comp_pe(decomp_standard_atmos) = 0
      decomp_db_last_comp_pe(decomp_standard_atmos) =                   &
     &  decomp_db_nproc(decomp_standard_atmos)-1

!     Set the grid size

      decomp_db_gridsize(1,decomp_standard_atmos) = nproc_EW
      decomp_db_gridsize(2,decomp_standard_atmos) = nproc_NS
      decomp_db_gridsize(3,decomp_standard_atmos) = 1


! Work out the decomposition in the East-West direction. As far as
! possible each processor has the same local row length. However, if
! this is not possible, the extra points are distributed symetrically
! such that each processor has the same number of points as the
! processor on the opposite side of the globe.

      start=1
      size=global_row_len/nproc_EW  ! local data size on each processor
                                    ! assuming nproc_EW divides exactly
                                    ! into global_row_len.
      irest=global_row_len-(size*nproc_EW)
                                    ! If it doesn't divide exactly then
                                    ! irest contains the number of left
                                    ! over points that need to be
                                    ! allocated to processors

! Check the domains are big enough for the extended halos
      IF ((mype  ==  0) .AND. (size  <   extended_halo_EW)) THEN
        ErrorStatus=5
        WRITE(Cmessage,                                                 &
     &    '("Too many processors in the East-West direction ",          &
     &    "(",I3,") to support the extended halo size ",                &
     &    "(",I3,"). Try running with ",I3," processors.")')            &
     &    nproc_EW,extended_halo_EW,(global_row_len/extended_halo_EW)
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF


      DO iproc=1,nproc_EW
        start_x(iproc-1)=start

        IF (iproc  <=  nproc_EW/2) THEN
          ! processor in the first half of the row
          IF (iproc  <=  (irest/2)) THEN
            size_x(iproc-1)=size+1 ! gets one of the extra points
          ELSE
            size_x(iproc-1)=size   ! gets the standard row length
          ENDIF
        ELSE
          ! processor in the second half of the row
          IF (iproc-(nproc_EW/2)  <=  (irest/2)) THEN
            size_x(iproc-1)=size+1 ! gets one of the extra points
          ELSE
            size_x(iproc-1)=size   ! gets the standard row length
          ENDIF
        ENDIF

        start=start+size_x(iproc-1)

      ENDDO

! Work out the decomposition in the East-West direction for the
! river routing grid. As far as possible all processors are given
! the same number of points. If this is not possible the extra "n"
! points are given to the first "n" processors in the LPG.
      start=1
      size=river_row_length/nproc_EW  !local data size on each processor
                                      !assuming nproc_EW divides exactly
                                      !into global_row_len.
      irest=river_row_length-(size*nproc_EW)
                                      !If it doesn't divide exactly then
                                      !irest contains the number of left
                                      !over points that need to be
                                      !allocated to processors

      DO iproc=0, nproc_EW-1
        rstart_x( iproc ) = start

        IF (iproc < irest) THEN
          rsize_x( iproc ) = size + 1
        ELSE
          rsize_x( iproc ) = size
        END IF

        start = start + rsize_x( iproc )
      END DO

! Work out the decomposition in the North-South direction. As far as
! possible each processor has the same number of rows. However, if this
! is not possible, the extra rows are distributed thus:
! - an extra row is given to the Northern most processor
! - the remaining extra rows are distributed symetrically around the
!   equator, starting at the processor(s) closest to the equator.

      start=1
      size=global_n_rows/nproc_NS  ! local data size on each processor
                                   ! assuming nproc_NS divides exactly
                                   ! into global_n_rows
      irest=global_n_rows-(size*nproc_NS)
                                   ! If it doesn't divide exactly then
                                   ! irest contains the number of left
                                   ! over points that need to be
                                   ! allocated to processors

! Check the domains are big enough for the extended halos
      IF ((mype  ==  0) .AND. (size  <=  extended_halo_NS)) THEN
        ErrorStatus=5
        WRITE(Cmessage,                                                 &
     &    '("Too many processors in the North-South direction ",        &
     &    "(",I3,") to support the extended halo size ",                &
     &    "(",I3,"). Try running with ",I3," processors.")')            &
     &    nproc_NS,extended_halo_NS,                                    &
     &    (global_n_rows / (extended_halo_NS+1) )
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF


      DO iproc=1,nproc_NS
        size_y(iproc-1)=size
      ENDDO

      IF (irest  >=  1) THEN
        ! give Northern most processors an extra row
        size_y(nproc_NS-1)=size+1
        irest=irest-1
      ENDIF

! Set up pointers to processor rows to which we will add extra rows
! to. These start around the equator, and will work out towards
! the poles.

      IF (MOD(nproc_NS,2)  ==  0) THEN  ! Even number of NS processors
        prow_S=nproc_NS/2
        prow_N=prow_S+1
      ELSE  ! Odd number of NS processors
        prow_S=(nproc_NS/2)+1
        prow_N=prow_S
      ENDIF

      DO WHILE (irest  >=  1)

        IF (prow_N  ==  prow_S) THEN
          size_y(prow_N-1)=size+1
          irest=irest-1
        ELSE
          size_y(prow_S-1)=size+1
          irest=irest-1
          IF (irest  >=  1) THEN
            size_y(prow_N-1)=size+1
            irest=irest-1
          ENDIF
        ENDIF

        prow_S=MAX(1,prow_S-1)
        prow_N=MIN(nproc_NS,prow_N+1)

      ENDDO

      DO iproc=1,nproc_NS
        start_y(iproc-1)=start
        start=start+size_y(iproc-1)
      ENDDO

! Work out the decomposition in the North-South direction for the
! river routing grid. As far as possible all processors are given
! the same number of rows. If this is not possible the extra "n" rows
! are given to the first "n" processors in the LPG.
      start=1
      size=river_rows/nproc_NS  ! local data size on each processor
                                ! assuming nproc_NS divides exactly
                                ! into global_row_len.
      irest=river_rows - (size*nproc_NS)
                                ! If it doesn't divide exactly then
                                ! irest contains the number of left
                                ! over points that need to be
                                ! allocated to processors

      DO iproc=0, nproc_NS-1
        rstart_y( iproc ) = start

        IF (iproc < irest) THEN
          rsize_y( iproc ) = size + 1
        ELSE
          rsize_y( iproc ) = size
        END IF

        start = start + rsize_y( iproc )
      END DO

! Set the local data shape and offsets of each processor

      DO iproc_y=0,nproc_NS-1

        IF (iproc_y  ==  (nproc_NS-1)) THEN
          at_north=.TRUE.  ! Doing the nothernmost processors
        ELSE
         at_north=.FALSE. ! Not doing the northernmost processors
        ENDIF

        DO iproc_x=0,nproc_EW-1

          iproc=decomp_db_first_comp_pe(decomp_standard_atmos)+         &
     &          iproc_x+(iproc_y*nproc_EW)

! Set the position in the Logical Processor Grid

          decomp_db_g_gridpos(1,iproc,decomp_standard_atmos)=iproc_x
          decomp_db_g_gridpos(2,iproc,decomp_standard_atmos)=iproc_y
          decomp_db_g_gridpos(3,iproc,decomp_standard_atmos)=0

! Set the number of local datapoints (blsize) on the processor

!  Fields on P grid:
          decomp_db_g_blsize(1,fld_type_p,iproc,                        &
     &                       decomp_standard_atmos)=                    &
     &      size_x(iproc_x)
          decomp_db_g_blsize(2,fld_type_p,iproc,                        &
     &                       decomp_standard_atmos)=                    &
     &      size_y(iproc_y)
          decomp_db_g_blsize(3,fld_type_p,iproc,                        &
     &                       decomp_standard_atmos)=                    &
     &      tot_levels

! Fields on U grid:
          decomp_db_g_blsize(1,fld_type_u,iproc,                        &
     &                       decomp_standard_atmos)=                    &
     &      size_x(iproc_x)
          decomp_db_g_blsize(2,fld_type_u,iproc,                        &
     &                       decomp_standard_atmos)=                    &
     &      size_y(iproc_y)
          decomp_db_g_blsize(3,fld_type_u,iproc,                        &
     &                       decomp_standard_atmos)=                    &
     &      tot_levels

! Fields on V grid:
          decomp_db_g_blsize(1,fld_type_v,iproc,                        &
     &                       decomp_standard_atmos)=                    &
     &      size_x(iproc_x)
          IF (at_north) THEN
            decomp_db_g_blsize(2,fld_type_v,iproc,                      &
     &                         decomp_standard_atmos)=                  &
     &      size_y(iproc_y)-1
            If (model_type  ==  mt_bi_cyclic_lam ) then
              decomp_db_g_blsize(2,fld_type_v,iproc,                    &
     &                           decomp_standard_atmos)=                &
     &        size_y(iproc_y)
            endif
          ELSE
            decomp_db_g_blsize(2,fld_type_v,iproc,                      &
     &                         decomp_standard_atmos)=                  &
     &      size_y(iproc_y)
          ENDIF
          decomp_db_g_blsize(3,fld_type_v,iproc,                        &
     &                       decomp_standard_atmos)=                    &
     &      tot_levels

!  Fields on R (river) grid:
          decomp_db_g_blsize(1,fld_type_r,iproc,                        &
     &                       decomp_standard_atmos)=                    &
     &      rsize_x(iproc_x)
          decomp_db_g_blsize(2,fld_type_r,iproc,                        &
     &                       decomp_standard_atmos)=                    &
     &      rsize_y(iproc_y)
          decomp_db_g_blsize(3,fld_type_r,iproc,                        &
     &                       decomp_standard_atmos)=                    &
     &      1

! Set the number of points including the halos on the processor

          DO ihalo=1,NHalo_max
            DO ifld=1,Nfld_max
              DO idim=1,Ndim_max
                decomp_db_g_lasize(idim,ifld,ihalo,iproc,               &
     &                             decomp_standard_atmos)=              &
     &            decomp_db_g_blsize(idim,ifld,iproc,                   &
     &                               decomp_standard_atmos)+            &
     &            2*decomp_db_halosize(idim,ihalo,                      &
     &                               decomp_standard_atmos)
              ENDDO  ! idim
            ENDDO  ! ifld
          ENDDO  ! ihalo

! Set the starting point in the global domain

          decomp_db_g_datastart(1,iproc,decomp_standard_atmos)=         &
     &      start_x(iproc_x)
          decomp_db_g_datastart(2,iproc,decomp_standard_atmos)=         &
     &      start_y(iproc_y)
          decomp_db_g_datastart(3,iproc,decomp_standard_atmos)=1

          decomp_db_g_datastart_f(1,fld_type_p,iproc,                   &
     &                decomp_standard_atmos)=start_x(iproc_x)
          decomp_db_g_datastart_f(2,fld_type_p,iproc,                   &
     &                decomp_standard_atmos)=start_y(iproc_y)
          decomp_db_g_datastart_f(3,fld_type_p,iproc,                   &
     &                decomp_standard_atmos)=1

          decomp_db_g_datastart_f(1,fld_type_u,iproc,                   &
     &                decomp_standard_atmos)=start_x(iproc_x)
          decomp_db_g_datastart_f(2,fld_type_u,iproc,                   &
     &                decomp_standard_atmos)=start_y(iproc_y)
          decomp_db_g_datastart_f(3,fld_type_u,iproc,                   &
     &                decomp_standard_atmos)=1

          decomp_db_g_datastart_f(1,fld_type_v,iproc,                   &
     &                decomp_standard_atmos)=start_x(iproc_x)
          decomp_db_g_datastart_f(2,fld_type_v,iproc,                   &
     &                decomp_standard_atmos)=start_y(iproc_y)
          decomp_db_g_datastart_f(3,fld_type_v,iproc,                   &
     &                decomp_standard_atmos)=1

          decomp_db_g_datastart_f(1,fld_type_r,iproc,                   &
     &                decomp_standard_atmos)=rstart_x(iproc_x)
          decomp_db_g_datastart_f(2,fld_type_r,iproc,                   &
     &                decomp_standard_atmos)=rstart_y(iproc_y)
          decomp_db_g_datastart_f(3,fld_type_r,iproc,                   &
     &                decomp_standard_atmos)=1

        ENDDO ! iproc_x
      ENDDO ! iproc_y

! Set up the pe_index_EW array - for each point along a global row
! it indicates the PE index (along the processor row) which
! contains that point

      DO iproc_x=0,nproc_EW-1
        DO ipt=decomp_db_g_datastart(1,iproc_x,decomp_standard_atmos),  &
     &         decomp_db_g_datastart(1,iproc_x,decomp_standard_atmos)+  &
     &         size_x(iproc_x)
          decomp_db_g_pe_index_EW(ipt,decomp_standard_atmos)=iproc_x
        ENDDO
      ENDDO

! And fill in the halos at either end

      DO ipt=1-extended_halo_EW,0
        decomp_db_g_pe_index_EW(ipt,decomp_standard_atmos)=0
      ENDDO

      DO ipt=global_row_len+1,global_row_len+extended_halo_EW
        decomp_db_g_pe_index_EW(ipt,decomp_standard_atmos)=nproc_x-1
      ENDDO

! Now set up the pe_index_NS_array - for each point along a global
! North-South column it indicates the PE index (along the processor
! column) which contains that point

      DO iproc_y=0,nproc_NS-1
        DO ipt=decomp_db_g_datastart(2,iproc_y*nproc_EW,                &
     &                               decomp_standard_atmos),            &
     &         decomp_db_g_datastart(2,iproc_y*nproc_EW,                &
     &                               decomp_standard_atmos)+            &
     &         size_y(iproc_y)
          decomp_db_g_pe_index_NS(ipt,decomp_standard_atmos)=iproc_y
        ENDDO
      ENDDO

! And fill in the halos at either end

      DO ipt=1-extended_halo_NS,0
        decomp_db_g_pe_index_NS(ipt,decomp_standard_atmos)=0
      ENDDO

      DO ipt=global_n_rows+1,global_n_rows+extended_halo_NS
        decomp_db_g_pe_index_NS(ipt,decomp_standard_atmos)=nproc_y-1
      ENDDO



! ------------------------------------------------------------------
! 3.0 Set boundary conditions
! ------------------------------------------------------------------


      IF (model_type  ==  mt_global) THEN
        decomp_db_bound(1,decomp_standard_atmos) = BC_CYCLIC
        decomp_db_bound(2,decomp_standard_atmos) = BC_OVERPOLE
      ELSEIF (model_type  ==  mt_lam) THEN
        decomp_db_bound(1,decomp_standard_atmos) = BC_STATIC
        decomp_db_bound(2,decomp_standard_atmos) = BC_STATIC
      ELSEIF (model_type  ==  mt_cyclic_lam) THEN
        decomp_db_bound(1,decomp_standard_atmos) = BC_CYCLIC
        decomp_db_bound(2,decomp_standard_atmos) = BC_STATIC
      ELSEIF (model_type  ==  mt_bi_cyclic_lam) THEN
        decomp_db_bound(1,decomp_standard_atmos) = BC_CYCLIC
        decomp_db_bound(2,decomp_standard_atmos) = BC_CYCLIC
         ELSE
           ErrorStatus=10
           WRITE(Cmessage,                                              &
     &       '("Unrecognised model_type: ",I3)') model_type
! DEPENDS ON: ereport
           CALL Ereport(RoutineName,ErrorStatus,Cmessage)
       ENDIF
       decomp_db_bound(3,decomp_standard_atmos) = BC_STATIC

! DEPENDS ON: set_neighbour
      CALL SET_NEIGHBOUR(                                               &
     &  decomp_standard_atmos)

! ------------------------------------------------------------------
! 4.0 Return the new data sizes and exit subroutine
! ------------------------------------------------------------------

! Set up the GCOM groups:

! 1) Group of all processors on my row

      IF ( decomp_db_gridsize(2,decomp_standard_atmos)  ==  1)          &
     & THEN
       decomp_db_gc_proc_row_group(decomp_standard_atmos)=GCG_ALL
      ELSE
        CALL GCG_SPLIT(mype,nproc_max,                                  &
     &    decomp_db_g_gridpos(2,mype,decomp_standard_atmos),            &
     &    info,                                                         &
     &    decomp_db_gc_proc_row_group(decomp_standard_atmos))
      ENDIF

! 2) Group of all processors on my column

      IF ( decomp_db_gridsize(1,decomp_standard_atmos)  ==  1)          &
     & THEN
        decomp_db_gc_proc_col_group(decomp_standard_atmos)=GCG_ALL
      ELSE
        CALL GCG_SPLIT(mype,nproc_max,                                  &
     &    decomp_db_g_gridpos(1,mype,decomp_standard_atmos),            &
     &    info,                                                         &
     &    decomp_db_gc_proc_col_group(decomp_standard_atmos))
      ENDIF

! 3) Group of all processors in the atmosphere model
      IF (decomp_db_nproc(decomp_standard_atmos)  ==  nproc_max)        &
     & THEN
        decomp_db_gc_all_proc_group(decomp_standard_atmos)=GCG_ALL
      ELSE
        IF ((mype  >=  decomp_db_first_comp_pe(decomp_standard_atmos))  &
     &    .AND.                                                         &
     &     (mype  <=  decomp_db_last_comp_pe(decomp_standard_atmos)) )  &
     &   THEN
          in_atm_decomp=1
        ELSE
          in_atm_decomp=0
        ENDIF

        CALL GCG_SPLIT(mype,nproc_max,in_atm_decomp,info,               &
     &    decomp_db_gc_all_proc_group(decomp_standard_atmos))
      ENDIF

! Set logical indicating this decomposition has been initialised
! and is now ready for use

      decomp_db_set(decomp_standard_atmos)=.TRUE.

! And return the new horizontal dimensions

      local_row_len=decomp_db_g_blsize(1,fld_type_p,mype,               &
     &                                 decomp_standard_atmos)
      local_n_rows=decomp_db_g_blsize(2,fld_type_p,mype,                &
     &                                decomp_standard_atmos)

!Check that the decomposition is valid for the mes.  The halos and the
!rimwidth must be contained within a single processor - they cannot
!cross over into another one as the indexing for the lbcs will go
!awry.

      If (model_type  ==  mt_lam) then
        if ((mype  ==  nproc_EW*(nproc_NS-1)) .or.                      &
     &      (mype  ==  nproc_EW-1) ) then   !ie nw corner or se corner
          Do idim = 1, nrima_max
      if(nproc_ns  == 1)then
            size=max(2*rimwidth(idim) +1,                               &
     &               extended_halo_NS + rimwidth(idim) +2 )
      else
        size=extended_halo_NS + rimwidth(idim) +2
      endif
      size1=extended_halo_NS + rimwidth(idim) +2
            If(size  >    local_n_rows ) then
              max_NS = global_n_rows /size1
              Errorstatus = 4
              Write(Cmessage,                                           &
     &         '("Too many processors in the North-South direction.",   &
     &           "The maximum permitted is ",I2)') max_NS
! DEPENDS ON: ereport
              CALL Ereport(RoutineName,ErrorStatus,Cmessage)
            End If
          End Do
        Else If (mype  ==  nproc_EW*nproc_NS-1 .or.                     &
     &           mype  ==  0 ) then        ! ie sw or ne corner
          Do idim = 1, nrima_max
      if(nproc_EW  == 1)then
            size=max(2*rimwidth(idim) +1,                               &
     &               extended_halo_EW + rimwidth(idim) +2 )
      else
        size=extended_halo_EW + rimwidth(idim) +2
      endif
      size1=extended_halo_EW + rimwidth(idim) +2
            If(size  >   local_row_len ) then
              max_EW = global_row_len / size1
              Errorstatus = 4
              Write(Cmessage,                                           &
     &         '("Too many processors in the East-West direction.",     &
     &           "The maximum permitted is ",I2)') 2* int(max_EW/2)
! DEPENDS ON: ereport
              CALL Ereport(RoutineName,ErrorStatus,Cmessage)
            End If
          End Do
        End If
      End If
      RETURN

      END SUBROUTINE DECOMPOSE_ATMOS

#endif
#endif
