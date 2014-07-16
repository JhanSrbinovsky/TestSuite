#if defined(C96_1C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
! This routine is a direct inverse of GATHER_field_ML. It takes
! full fields (selected levels - chosen by way of a map) and
! scatters them to decomposed data.
!
! Method:
!  This routine copies all global data (selected levels) into a single
!  array which is then sent to all CPUs where it is unpacked into
!  local data. MPL used for speed.
!
! Current Code Owner: Paul Selwood
!
! Subroutine Interface:
      Subroutine Scatter_Field_ML(                                      &
     &    local_field,    global_field,                                 &
     &    local_row_len,  local_rows,    local_levs,                    &
     &    global_row_len, global_rows,   global_levs,                   &
     &    pe_for_level,                                                 &
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
      Integer, Intent(In) :: pe_for_level(local_levs)  ! PE to scatter
                                                       ! level from

      ! Scattered data
      Real, Intent(Out)   :: local_field( local_row_len,                &
     &                                    local_rows, local_levs )
      ! Original data
      Real, Intent(In)    :: global_field( global_row_len,              &
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
      Integer :: local_row_len_nh  ! local row length without halos
      Integer :: local_rows_nh     ! local rows without halos
      Integer :: pos               ! buffer position

      Integer :: levs_to_send(0 : nproc-1) ! num of levs to send
      Integer :: kpos(0 : nproc-1)         ! buffer position
      Integer :: send_size(0 : nproc-1)    ! size to send
      Integer :: recv_size(0 : nproc-1)    ! size to receive

      Integer :: ierr                      ! error flag
      Integer :: status(MPL_STATUS_SIZE)   ! MPL status
      Integer :: MY_COMM                   ! MPL Communicator


! Array to hold all local data - note contains space for halos
! that won't be used
      Real :: local_buffer(local_row_len * local_rows * global_levs,    &
     &                     0 : nproc-1)

! Array to hold received data - much too big but can't think of
! how to make it better sized at the moment.
      Real :: send_buff(global_row_len * global_rows * global_levs,     &
     &                  0 : nproc -1)


!-------------------------------------------------------
! 0) Calculate a few sizes I'll need later
!-------------------------------------------------------
      halo_x = halosize(1, halo_type)
      halo_y = halosize(2, halo_type)

! Non halo sizes
      local_row_len_nh = local_row_len - (2 * halo_x)
      local_rows_nh = local_rows - (2 * halo_y)

! Find sizes to send and receive
      levs_to_send(:) = 0
      Do k = 1, local_levs
        levs_to_send(pe_for_level(k)) = levs_to_send(pe_for_level(k))+1
      End Do

!-------------------------------------------------------
! 0) Setup  - get communicator from GCOM
!-------------------------------------------------------
      Call gc_get_communicator(MY_COMM, ierr)

!-------------------------------------------------------
! 1) Copy data from global fields into send buffer
!-------------------------------------------------------
      Do iproc = 0, nproc - 1
        Do k = 1, levs_to_send(mype)
          Do j = 1, g_blsize(2,fld_type,iproc)
            Do i = 1, g_blsize(1,fld_type,iproc)
              send_buff(i + (j - 1) * g_blsize(1,fld_type,iproc) +      &
     &                    (k - 1) * g_blsize(1,fld_type,iproc) *        &
     &                    g_blsize(2,fld_type,iproc), iproc)            &
     &        = global_field( g_datastart_f(1,fld_type,iproc) + i - 1,  &
     &                      g_datastart_f(2,fld_type,iproc) + j - 1 ,k)
            End Do
          End Do
        End Do
      End Do

!-------------------------------------------------------
! 2) Find sizes for send/recv and do the communications
!    Use MPL_Sendrecv to pair up comms.
!-------------------------------------------------------
      Do iproc = 0, nproc - 1
        recv_size(iproc) = local_row_len_nh * local_rows_nh *           &
     &                     levs_to_send(iproc)
        send_size(iproc) = g_blsize(1,fld_type,iproc) *                 &
     &                     g_blsize(2,fld_type,iproc) *                 &
     &                     levs_to_send(mype)
      End Do

! Do communications using MPL directly
      Do iproc = 0, nproc - 1
        Call MPL_Sendrecv( send_buff(1,iproc), send_size(iproc),        &
     &                     MPL_Real, iproc, 999, local_buffer(1,iproc), &
     &                     recv_size(iproc), MPL_Real, iproc, 999,      &
     &                     MY_COMM, status, ierr)

      End Do

!-------------------------------------------------------
! 3) Copy data from received buffer into proper
!    decomposed data locations
!-------------------------------------------------------
      Do iproc = 0, nproc - 1
        kpos(iproc) = 0
      End Do


! Copy local_buffer (no halos) into local_field (with halos)
! Need to get levels right too.
      Do k = 1, local_levs
        Do j = 1+halo_y, local_rows - halo_y
          Do i = 1+halo_x, local_row_len - halo_x
            pos = i - halo_x +                                          &
     &            (j - halo_y - 1) * local_row_len_nh +                 &
     &            kpos(pe_for_level(k)) * local_rows_nh *               &
     &                                    local_row_len_nh

            local_field(i,j,k)                                          &
     &      = local_buffer(pos,pe_for_level(k))
          End Do
        End Do
        kpos(pe_for_level(k)) = kpos(pe_for_level(k)) + 1
      End Do


      Return
      END SUBROUTINE Scatter_Field_ML
#endif
