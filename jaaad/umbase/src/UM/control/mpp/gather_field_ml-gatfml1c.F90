#if defined(C96_1C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Takes 1 or more levels of a model field that have been decomposed
!  over a group of processors, and gathers the data together so that
!  one complete global level is contained on one processor.  Processors
!  can hold one or more such levels, as determined by the 'local_level'
!  array, which gives the index for each level on the processor to
!  which it is sent.  For one level per PE, the setting of the values
!  local_level(1...local_levs) should all be one.  Successive ones
!  obviously range from 1 upwards.
!
! Method:
!  This routine copies all local data (all levels) into a single
!  array which is then sent to the gathering CPU where it is unpacked.
!
! Current Code Owner: Paul Selwood
!
!
! Subroutine Interface:
      Subroutine Gather_Field_ML(                                       &
     &    local_field,    global_field,                                 &
     &    local_row_len,  local_rows,   local_levs,                     &
     &    global_row_len, global_rows,  global_levs,                    &
     &    pe_for_level,   local_level,                                  &
     &    fld_type,       halo_type)

      Use MPL, Only :           &
               MPL_REAL,        &
               MPL_STATUS_SIZE

      Implicit None
!
! Subroutine Arguments:

      Integer, Intent(In) :: local_row_len   ! local field row length
      Integer, Intent(In) :: local_rows      ! local field rows
      Integer, Intent(In) :: local_levs      ! local field levels
      Integer, Intent(In) :: global_row_len  ! global field row length
      Integer, Intent(In) :: global_rows     ! global field rows
      Integer, Intent(In) :: global_levs     ! global field levels
      Integer, Intent(In) :: fld_type        ! field type of grid
      Integer, Intent(In) :: halo_type       ! halo type of field
      Integer, Intent(In) :: pe_for_level(local_levs)  ! PE to gather
                                                       ! each level to
      Integer, Intent(In) :: local_level(local_levs)   ! Not used. Give
                                                       ! same arguments
                                                       ! as C96_1B

      ! Decomposed data
      Real, Intent(In)    :: local_field( local_row_len,                &
     &                                    local_rows, local_levs )
      ! Gathered data
      Real, Intent(Out)   :: global_field( global_row_len,              &
     &                                     global_rows, global_levs )

! Parameters and Common blocks

#include "parvars.h"

! Local variables

      Integer :: i       ! loop index  - cols
      Integer :: j       ! loop index  - rows
      Integer :: k       ! loop index  - levels
      Integer :: iproc   ! loop index  - processors
      Integer :: halo_x  ! halo size - x
      Integer :: halo_y  ! halo size - y
      Integer :: local_row_len_nh  ! row length without halo
      Integer :: local_rows_nh     ! row count without halo
      Integer :: pos               ! position in array

      Integer :: levs_to_send(0 : nproc-1) ! num of levels to each PE
      Integer :: kpos(0 : nproc-1)         ! array position
      Integer :: send_size(0 : nproc-1)    ! size of sent data
      Integer :: recv_size(0 : nproc-1)    ! size of received data

      Integer :: ierr                      ! error code
      Integer :: status(MPL_STATUS_SIZE)   ! MPL status
      Integer :: MY_COMM                   ! MPL communicator


! Array to hold all local data - note contains space for halos
! that won't be used
      Real :: local_buffer(local_row_len * local_rows * global_levs,    &
     &                     0 : nproc-1)

! Array to hold received data - much too big but best I can
! think of
      Real :: recv_buff(global_row_len * global_rows * global_levs,     &
     &                  0 : nproc -1)


!-------------------------------------------------------
! 0) Setup  - get communicator from GCOM
!-------------------------------------------------------
      Call gc_get_communicator(MY_COMM, ierr)

!-------------------------------------------------------
! 1) Copy data from multiple levels into 1 buffer
!    for sending - this reduces the number of messages
!    sent.
!-------------------------------------------------------
      halo_x = halosize(1, halo_type)
      halo_y = halosize(2, halo_type)
! Non halo sizes
      local_row_len_nh = local_row_len - (2 * halo_x)
      local_rows_nh = local_rows - (2 * halo_y)

      Do iproc = 0, nproc - 1
        kpos(iproc) = 0
      End Do


! Copy local_field into local_buffer with halos stripped off
! and with contiguous levels (ie 1 message per PE)
      Do k = 1, local_levs
        Do j = 1+halo_y, local_rows - halo_y
          Do i = 1+halo_x, local_row_len - halo_x
            pos = i - halo_x +                                          &
     &            (j - halo_y - 1) * local_row_len_nh +                 &
     &            kpos(pe_for_level(k)) * local_rows_nh *               &
     &                                    local_row_len_nh
            local_buffer(pos,pe_for_level(k)) =                         &
     &      local_field(i,j,k)
          End Do
        End Do
        kpos(pe_for_level(k)) = kpos(pe_for_level(k)) + 1
      End Do


!-------------------------------------------------------
! 2) Find sizes for send/recv and do the communications
!    Use MPL_Sendrecv to pair up comms.
!-------------------------------------------------------
! Find sizes to send and receive
      levs_to_send(:) = 0
      Do k = 1, local_levs
        levs_to_send(pe_for_level(k)) = levs_to_send(pe_for_level(k))+1
      End Do

      Do iproc = 0, nproc - 1
        send_size(iproc) = local_row_len_nh * local_rows_nh *           &
     &                     levs_to_send(iproc)
        recv_size(iproc) = g_blsize(1,fld_type,iproc) *                 &
     &                     g_blsize(2,fld_type,iproc) *                 &
     &                     levs_to_send(mype)
      End Do

! Do communications using MPL
      Do iproc = 0, nproc - 1
        Call MPL_Sendrecv( local_buffer(1,iproc), send_size(iproc),     &
     &                     MPL_REAL, iproc, 999, recv_buff(1,iproc),    &
     &                     recv_size(iproc), MPL_REAL, iproc, 999,      &
     &                     MY_COMM, status, ierr)

      End Do

!-------------------------------------------------------
! 3) Copy data from recv buffers into final data field
!-------------------------------------------------------
      Do iproc = 0, nproc - 1
        Do k = 1, levs_to_send(mype)
          Do j = 1, g_blsize(2,fld_type,iproc)
            Do i = 1, g_blsize(1,fld_type,iproc)
              global_field( g_datastart_f(1,fld_type,iproc) + i - 1,    &
     &                      g_datastart_f(2,fld_type,iproc) + j - 1 ,k) &
     &        = recv_buff(i + (j - 1) * g_blsize(1,fld_type,iproc) +    &
     &                    (k - 1) * g_blsize(1,fld_type,iproc) *        &
     &                    g_blsize(2,fld_type,iproc), iproc)
            End Do
          End Do
        End Do
      End Do

      Return
      End Subroutine Gather_Field_ML
#endif
