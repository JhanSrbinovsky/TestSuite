
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine  Q_Pos_Ctl

      subroutine Q_Pos_Ctl                                              &
     &          (q, row_length, rows, wet_model_levels,                 &
     &           global_row_length, global_rows,                        &
     &           me, n_proc, off_x, off_y,                              &
     &           gc_all_proc_group,                                     &
     &           model_domain,                                          &
     &           halo_type, l_q_pos_local, qlimit)

!
!
! Method:
!         Each level of q is gathered onto a processor. (Gather_Field)
!         The q_pos routine is called on each level.
!         The results are then returned back to the original processors.
!                                                       (Scatter_field)
!
! Original Progammer: Andrew J Malcolm
! Current code owner: Andrew J Malcolm
!
! History:
! Version    Date      Comment
! ----     -------     -------
! 6.1      24/04/04    Rewrite from scratch to correct myriad problems
!                                     Andy Malcolm
! 6.2      06/12/05    Add multi-level gather/scatter. P.Selwood
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

      implicit none

! Arguments with Intent IN. ie: Input variables.

      integer                                                           &
     &  row_length                                                      &
                        ! number of points on a row
     &, rows                                                            &
                        ! number of rows of data
     &, wet_model_levels                                                &
     &, global_row_length                                               &
     &, global_rows                                                     &
     &, me                                                              &
                   ! processor number
     &, n_proc                                                          &
     &, off_x                                                           &
                  ! Size of halo in i direction.
     &, off_y                                                           &
                  ! Size of halo in j direction.
     &, gc_all_proc_group                                               &
     &, halo_type                                                       &
     &, model_domain

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

      logical                                                           &
     &  l_q_pos_local   ! logical : true to do Method 2,
                        !         : false to do Method 1

      real                                                              &
     &  qlimit

! Arguments with Intent IN/OUT.

      real                                                              &
     &  q(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      wet_model_levels)

! Local Variables.

      real                                                              &
     &  q_global (global_row_length, global_rows,                       &
     &            (wet_model_levels/n_proc)+1)

      integer                                                           &
     & k

      INTEGER                                                           &
     &  MAP(wet_model_levels)                                           &
                                    ! processor number for level
     &, N_LEVS_ON_PROC(0:N_PROC-1)                                      &
                                   ! number of levels on each processor
     &, MAX_LEVS_PER_CPU                                                &
                                    ! max number of levs per cpu
     &, LEV_ON_GAT_PE(WET_MODEL_LEVELS)                                 &
                                        ! for shmem only
     &, icode                       ! return code from comms.

      CHARACTER*80                                                      &
     &  CMESSAGE         ! OUT error message

! Calculate the mapping of which processor each level will go to
      DO K=0,N_PROC-1
        N_LEVS_ON_PROC(K)=0
      ENDDO

      DO K=1, wet_model_levels
        MAP(K)                 = MOD((K-1),N_PROC)
        N_LEVS_ON_PROC(MAP(K)) = N_LEVS_ON_PROC(MAP(K)) + 1
        LEV_ON_GAT_PE(K)       = N_LEVS_ON_PROC(MAP(K))
      ENDDO

! Distribute Q over the processors
      MAX_LEVS_PER_CPU = WET_MODEL_LEVELS / N_PROC + 1

! DEPENDS ON: gather_field_ml
      CALL GATHER_FIELD_ML(                                             &
     &                      Q,                    Q_GLOBAL,             &
     &                      ROW_LENGTH + 2*OFF_X, ROWS + 2*OFF_Y,       &
     &                      WET_MODEL_LEVELS,                           &
     &                      GLOBAL_ROW_LENGTH,    GLOBAL_ROWS,          &
     &                      MAX_LEVS_PER_CPU,                           &
     &                      MAP,                  LEV_ON_GAT_PE,        &
     &                      FLD_TYPE_P,           HALO_TYPE )

! and call Q_POS with the whole levels on this processor

! DEPENDS ON: q_pos
      call q_pos(                                                       &
     &           q_global, global_row_length, global_rows,              &
     &           N_LEVS_ON_PROC(me),                                    &
     &           l_q_pos_local, qlimit ,                                &
     &           model_domain)

! and scatter back levels to original processors

! DEPENDS ON: scatter_field_ml
      CALL SCATTER_FIELD_ML(                                            &
     &                      Q,                    Q_GLOBAL,             &
     &                      ROW_LENGTH + 2*OFF_X, ROWS + 2*OFF_Y,       &
     &                      WET_MODEL_LEVELS,                           &
     &                      GLOBAL_ROW_LENGTH,    GLOBAL_ROWS,          &
     &                      MAX_LEVS_PER_CPU,     MAP,                  &
     &                      FLD_TYPE_P,           HALO_TYPE )



      RETURN
      END SUBROUTINE Q_Pos_Ctl


!----------------------------------------------------------------------
! Subroutine Q_POS - Does the work.
!----------------------------------------------------------------------

