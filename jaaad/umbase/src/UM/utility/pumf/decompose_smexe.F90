! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT 2000, Met Office, All Rights Reserved.
! Please refer to file $UMDIR/vn$VN/copyright.txt for further details
! *****************************COPYRIGHT*******************************

!+ Parallel UM: Perform 2D data decomposition

! Subroutine interface:
      SUBROUTINE DECOMPOSE_SMEXE(global_row_len, global_n_rows,         &
     &                           extended_halo_EW,extended_halo_NS,     &
     &                           TOT_LEVELS)
      IMPLICIT NONE

! Description:
! This routine performs a 2D decomposition - taking the global X
! (global_row_len) and Y (global_n_rows) data sizes and decomposing
! across 1 processor in both the X & Y direction for small executables.
! Local data size should be the same as Global field, as no halo
! is included.

! Method:
! Various data is stored in the COMMON block DECOMPDB for later use.


! Subroutine Arguments:

      INTEGER                                                           &

     &  global_row_len,                                                 &
                          ! IN  :number of E-W points of entire model
     &  global_n_rows,                                                  &
                          ! IN  :number of P rows of entire model
     &  extended_halo_EW,                                               &
                          ! IN  :EW halo size
     &  extended_halo_NS,                                               &
                          ! IN  :NS halo size
     &  tot_levels        ! IN  :total number of levels

! Parameters and Common blocks

#include "parvars.h"
#include "decomptp.h"
#include "decompdb.h"
#include "gccom.h"
#include "domtyp.h"




! Local variables
      INTEGER                                                           &
     &  nproc_EW                                                        &
     &, nproc_NS                                                        &
     &, iproc                                                           &
     &, iproc_x                                                         &
     &, iproc_y                                                         &
     &, ifld                                                            &
     &, ihalo                                                           &
     &, idim                                                            &
     &, ipt                                                             &
     &, irest                                                           &
     &, start                                                           &
     &, size                                                            &
     &, prow_N                                                          &
     &, prow_S                                                          &
     &, info                                                            &
     &, in_atm_decomp

      PARAMETER (  nproc_EW=1                                           &
     &,            nproc_NS=1)

      LOGICAL                                                           &
     &  at_north

      INTEGER                                                           &
     &  size_x(0:nproc_EW-1)                                            &
     &, size_y(0:nproc_NS-1)                                            &
     &, start_x(0:nproc_EW-1)                                           &
     &, start_y(0:nproc_NS-1)

! Error reporting
      INTEGER       ErrorStatus ! =0 normal exit; >0 error exit
      CHARACTER*256 Cmessage    ! Error message
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='DECOMPOSE_SMEXE')

! ------------------------------------------------------------------
! 0.0 Check for valid decomposition
! ------------------------------------------------------------------

! No need to check for valid decomposition as nproc_EW, nproc_NS
! extended_halo_EW & extended_halo_NS are set

! ------------------------------------------------------------------
      decomp_db_halosize(1,halo_type_single,decomp_smexe) = 1
      decomp_db_halosize(2,halo_type_single,decomp_smexe) = 1
      decomp_db_halosize(3,halo_type_single,decomp_smexe) = 0

      decomp_db_halosize(1,halo_type_extended,decomp_smexe) =           &
     &  extended_halo_EW
      decomp_db_halosize(2,halo_type_extended,decomp_smexe) =           &
     &  extended_halo_NS
      decomp_db_halosize(3,halo_type_extended,decomp_smexe) =           &
     &  0

      decomp_db_halosize(1,halo_type_no_halo,decomp_smexe) = 0
      decomp_db_halosize(2,halo_type_no_halo,decomp_smexe) = 0
      decomp_db_halosize(3,halo_type_no_halo,decomp_smexe) = 0


! ------------------------------------------------------------------
! 1.0 Set up global data size
! ------------------------------------------------------------------

      decomp_db_glsize(1,fld_type_p,decomp_smexe) =                     &
     &  global_row_len
      decomp_db_glsize(2,fld_type_p,decomp_smexe) =                     &
     &  global_n_rows
      decomp_db_glsize(3,fld_type_p,decomp_smexe) =                     &
     &  tot_levels

      decomp_db_glsize(1,fld_type_u,decomp_smexe) =                     &
     &  global_row_len
      decomp_db_glsize(2,fld_type_u,decomp_smexe) =                     &
     &  global_n_rows
      decomp_db_glsize(3,fld_type_u,decomp_smexe) =                     &
     &  tot_levels

      decomp_db_glsize(1,fld_type_v,decomp_smexe) =                     &
     &  global_row_len
      decomp_db_glsize(2,fld_type_v,decomp_smexe) =                     &
     &  GLOBAL_N_ROWS
      decomp_db_glsize(3,fld_type_v,decomp_smexe) =                     &
     &  tot_levels

! ------------------------------------------------------------------
! 2.0 Calculate decomposition
! ------------------------------------------------------------------


! select processors to use for the data decomposition
      decomp_db_nproc(decomp_smexe)=nproc_EW*nproc_NS
      decomp_db_first_comp_pe(decomp_smexe) = 0
      decomp_db_last_comp_pe(decomp_smexe) =                            &
     &  decomp_db_nproc(decomp_smexe)-1

!     Set the grid size

      decomp_db_gridsize(1,decomp_smexe) = nproc_EW
      decomp_db_gridsize(2,decomp_smexe) = nproc_NS
      decomp_db_gridsize(3,decomp_smexe) = 1


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

! Set the local data shape and offsets of each processor

      DO iproc_y=0,nproc_NS-1

        IF (iproc_y  ==  (nproc_NS-1)) THEN
          at_north=.TRUE.  ! Doing the nothernmost processors
        ELSE
         at_north=.FALSE. ! Not doing the northernmost processors
        ENDIF

        DO iproc_x=0,nproc_EW-1

          iproc=decomp_db_first_comp_pe(decomp_smexe)+                  &
     &          iproc_x+(iproc_y*nproc_EW)

! Set the position in the Logical Processor Grid

          decomp_db_g_gridpos(1,iproc,decomp_smexe)=iproc_x
          decomp_db_g_gridpos(2,iproc,decomp_smexe)=iproc_y
          decomp_db_g_gridpos(3,iproc,decomp_smexe)=0

! Set the number of local datapoints (blsize) on the processor

!  Fields on P grid:
          decomp_db_g_blsize(1,fld_type_p,iproc,                        &
     &                       decomp_smexe)=                             &
     &      size_x(iproc_x)
          decomp_db_g_blsize(2,fld_type_p,iproc,                        &
     &                       decomp_smexe)=                             &
     &      size_y(iproc_y)
          decomp_db_g_blsize(3,fld_type_p,iproc,                        &
     &                       decomp_smexe)=                             &
     &      tot_levels

! Fields on U grid:
          decomp_db_g_blsize(1,fld_type_u,iproc,                        &
     &                       decomp_smexe)=                             &
     &      size_x(iproc_x)
          decomp_db_g_blsize(2,fld_type_u,iproc,                        &
     &                       decomp_smexe)=                             &
     &      size_y(iproc_y)
          decomp_db_g_blsize(3,fld_type_u,iproc,                        &
     &                       decomp_smexe)=                             &
     &      tot_levels

#if defined(MAKEBC)
! Fields on V grid:
          decomp_db_g_blsize(1,fld_type_v,iproc,                        &
     &                       decomp_smexe)=                             &
     &      size_x(iproc_x)
          IF (at_north) THEN
            decomp_db_g_blsize(2,fld_type_v,iproc,                      &
     &                         decomp_smexe)=                           &
     &      size_y(iproc_y)-1
          ELSE
! debug - subtract 1 from size_y(iproc_y) to account for one
! less row in v fields
            decomp_db_g_blsize(2,fld_type_v,iproc,                      &
     &                         decomp_smexe)=                           &
     &      size_y(iproc_y)-1
! debug - subtract 1 from size_y(iproc_y) to account for one
! less row in v fields
          ENDIF
          decomp_db_g_blsize(3,fld_type_v,iproc,                        &
     &                       decomp_smexe)=                             &
     &      tot_levels
#else
! Fields on V grid:
          decomp_db_g_blsize(1,fld_type_v,iproc,                        &
     &                       decomp_smexe)=                             &
     &      size_x(iproc_x)
          IF (at_north) THEN
            decomp_db_g_blsize(2,fld_type_v,iproc,                      &
     &                         decomp_smexe)=                           &
     &      size_y(iproc_y)
          ELSE
            decomp_db_g_blsize(2,fld_type_v,iproc,                      &
     &                         decomp_smexe)=                           &
     &      size_y(iproc_y)
          ENDIF
          decomp_db_g_blsize(3,fld_type_v,iproc,                        &
     &                       decomp_smexe)=                             &
     &      tot_levels
#endif

! Set the number of points including the halos on the processor

          DO ihalo=1,NHalo_max
            DO ifld=1,Nfld_max
              DO idim=1,Ndim_max
                decomp_db_g_lasize(idim,ifld,ihalo,iproc,               &
     &                             decomp_smexe)=                       &
     &            decomp_db_g_blsize(idim,ifld,iproc,                   &
     &                               decomp_smexe)+                     &
     &            2*decomp_db_halosize(idim,ihalo,                      &
     &                               decomp_smexe)
              ENDDO  ! idim
            ENDDO  ! ifld
          ENDDO  ! ihalo

! Set the starting point in the global domain

          decomp_db_g_datastart(1,iproc,decomp_smexe)=                  &
     &      start_x(iproc_x)
          decomp_db_g_datastart(2,iproc,decomp_smexe)=                  &
     &      start_y(iproc_y)
          decomp_db_g_datastart(3,iproc,decomp_smexe)=1

          decomp_db_g_datastart_f(1,fld_type_p,iproc,                   &
     &                decomp_smexe)=start_x(iproc_x)
          decomp_db_g_datastart_f(2,fld_type_p,iproc,                   &
     &                decomp_smexe)=start_y(iproc_y)
          decomp_db_g_datastart_f(3,fld_type_p,iproc,                   &
     &                decomp_smexe)=1

          decomp_db_g_datastart_f(1,fld_type_u,iproc,                   &
     &                decomp_smexe)=start_x(iproc_x)
          decomp_db_g_datastart_f(2,fld_type_u,iproc,                   &
     &                decomp_smexe)=start_y(iproc_y)
          decomp_db_g_datastart_f(3,fld_type_u,iproc,                   &
     &                decomp_smexe)=1

          decomp_db_g_datastart_f(1,fld_type_v,iproc,                   &
     &                decomp_smexe)=start_x(iproc_x)
          decomp_db_g_datastart_f(2,fld_type_v,iproc,                   &
     &                decomp_smexe)=start_y(iproc_y)
          decomp_db_g_datastart_f(3,fld_type_v,iproc,                   &
     &                decomp_smexe)=1

          decomp_db_g_datastart_f(1,fld_type_r,iproc,                   &
     &                decomp_smexe)=0
          decomp_db_g_datastart_f(2,fld_type_r,iproc,                   &
     &                decomp_smexe)=0
          decomp_db_g_datastart_f(3,fld_type_r,iproc,                   &
     &                decomp_smexe)=1
        ENDDO ! iproc_x
      ENDDO ! iproc_y

! Set up the pe_index_EW array - for each point along a global row
! it indicates the PE index (along the processor row) which
! contains that point

      DO iproc_x=0,nproc_EW-1
        DO ipt=decomp_db_g_datastart(1,iproc_x,decomp_smexe),           &
     &         decomp_db_g_datastart(1,iproc_x,decomp_smexe)+           &
     &         size_x(iproc_x)
          decomp_db_g_pe_index_EW(ipt,decomp_smexe)=iproc_x
        ENDDO
      ENDDO

! And fill in the halos at either end

      DO ipt=1-extended_halo_EW,1
        decomp_db_g_pe_index_EW(ipt,decomp_smexe)=0
      ENDDO

      DO ipt=global_row_len+1,global_row_len+1+extended_halo_EW
        decomp_db_g_pe_index_EW(ipt,decomp_smexe)=nproc_x-1
      ENDDO

! Now set up the pe_index_NS_array - for each point along a global
! North-South column it indicates the PE index (along the processor
! column) which contains that point

      DO iproc_y=0,nproc_NS-1
        DO ipt=decomp_db_g_datastart(2,iproc_y*nproc_EW,                &
     &                               decomp_smexe),                     &
     &         decomp_db_g_datastart(2,iproc_y*nproc_EW,                &
     &                               decomp_smexe)+                     &
     &         size_y(iproc_y)
          decomp_db_g_pe_index_NS(ipt,decomp_smexe)=iproc_y
        ENDDO
      ENDDO

! And fill in the halos at either end

      DO ipt=1-extended_halo_NS,1
        decomp_db_g_pe_index_NS(ipt,decomp_smexe)=0
      ENDDO

      DO ipt=global_n_rows+1,global_n_rows+1+extended_halo_NS
        decomp_db_g_pe_index_NS(ipt,decomp_smexe)=nproc_y-1
      ENDDO

! ------------------------------------------------------------------
! 3.0 Return the new data sizes and exit subroutine
! ------------------------------------------------------------------

! Set up the GCOM groups:

! 1) Group of all processors on my row

      IF ( decomp_db_gridsize(2,decomp_smexe)  ==  1)                   &
     & THEN
       decomp_db_gc_proc_row_group(decomp_smexe)=GCG_ALL
      ELSE
        CALL GCG_SPLIT(mype,nproc_max,                                  &
     &    decomp_db_g_gridpos(2,mype,decomp_smexe),                     &
     &    info,                                                         &
     &    decomp_db_gc_proc_row_group(decomp_smexe))
      ENDIF

! 2) Group of all processors on my column

      IF ( decomp_db_gridsize(1,decomp_smexe)  ==  1)                   &
     & THEN
        decomp_db_gc_proc_col_group(decomp_smexe)=GCG_ALL
      ELSE
        CALL GCG_SPLIT(mype,nproc_max,                                  &
     &    decomp_db_g_gridpos(1,mype,decomp_smexe),                     &
     &    info,                                                         &
     &    decomp_db_gc_proc_col_group(decomp_smexe))
      ENDIF

! 3) Group of all processors in the atmosphere model
      IF (decomp_db_nproc(decomp_smexe)  ==  nproc_max)                 &
     & THEN
        decomp_db_gc_all_proc_group(decomp_smexe)=GCG_ALL
      ELSE
        IF ((mype  >=  decomp_db_first_comp_pe(decomp_smexe))           &
     &    .AND.                                                         &
     &     (mype  <=  decomp_db_last_comp_pe(decomp_smexe)) )           &
     &   THEN
          in_atm_decomp=1
        ELSE
          in_atm_decomp=0
        ENDIF

        CALL GCG_SPLIT(mype,nproc_max,in_atm_decomp,info,               &
     &    decomp_db_gc_all_proc_group(decomp_smexe))
      ENDIF

! Set logical indicating this decomposition has been initialised
! and is now ready for use

      decomp_db_set(decomp_smexe)=.TRUE.

      RETURN

      END SUBROUTINE DECOMPOSE_SMEXE
