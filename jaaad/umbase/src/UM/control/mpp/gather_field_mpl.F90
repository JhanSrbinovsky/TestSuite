#if defined(C96_1C) 
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Gathers a field from many processors to one processor
!
! Subroutine Interface:
SUBROUTINE gather_field_mpl(local_field,    global_field,             &
                            local_row_len,  local_rows,               &
                            global_row_len, global_rows,              &
                            grid_type,      halo_type,                &
                            gather_pe,      proc_group,               &
                            icode,          cmessage)

USE mpl, ONLY: &
    MPL_REAL

IMPLICIT NONE

!
! Description:
!  Takes a model field that has been decomposed over a group of
!  processors, and gathers the data together so that one processor
!  contains the entire global field.
!
! Method:
!  Uses mpl_gatherv to gather into a contiguous 1D array and then
!  copy data to final array location
!
! Current Code Owner: Paul Selwood
!
! Subroutine Arguments:

INTEGER, INTENT(IN) :: local_row_len   ! length of rows in local part of field
INTEGER, INTENT(IN) :: local_rows      ! number of rows in local part of field
INTEGER, INTENT(IN) :: global_row_len  ! length of rows in global field
INTEGER, INTENT(IN) :: global_rows     ! number of rows in global field
INTEGER, INTENT(IN) :: grid_type       ! type (p,u or v) of grid
INTEGER, INTENT(IN) :: halo_type       ! halo type (hence width) of grid
INTEGER, INTENT(IN) :: gather_pe       ! processor to gather global field to
INTEGER, INTENT(IN) :: proc_group      ! group id of processors involved here
INTEGER, INTENT(OUT) :: icode          ! out return code

REAL, INTENT(IN)     :: local_field(local_row_len*local_rows)
                                       ! local part of field
REAL, INTENT(INOUT)  :: global_field(global_row_len*global_rows)
                                       ! (on pe gather_pe) global field

CHARACTER (LEN=80), INTENT(OUT)  :: cmessage  ! error message

! Parameters and Common blocks

#include "parvars.h"

! Local variables

REAL :: global_field_copy( global_row_len * global_rows )
                        ! Space for re-ordered copy of global_field

INTEGER :: info          ! Status return from MPL/GCOM
INTEGER :: i,j, iproc    ! Loopers
INTEGER :: igpos         ! Position pointer
INTEGER :: igpos1        ! Position pointer
INTEGER :: ipos          ! Position pointer
INTEGER :: my_comm       ! Current communicator

INTEGER :: counts(0:nproc-1)               ! Counts for Gather
INTEGER :: displacements(0:nproc-1)        ! Displacements for Gather

INTEGER, SAVE :: old_global_row_len = -1234
INTEGER, SAVE :: old_global_rows    = -1234
INTEGER, SAVE :: old_proc_group     = -1234
INTEGER, SAVE :: old_gather_pe      = -1234
INTEGER, SAVE :: old_decomp         = -1234
INTEGER, SAVE :: old_grid_type      = -1234
INTEGER, SAVE :: old_halo_type      = -1234
INTEGER, SAVE :: send_type          = -1234

CHARACTER (LEN=*), PARAMETER :: routinename='gather_field_mpl'

!-------------------------------------------------------
IF (grid_type  ==  fld_type_unknown) THEN
  WRITE (6,*) 'GATHER_FIELD_MPL : bad field type'
  WRITE (6,*) 'field will not be gathered.'
  WRITE (cmessage,*) 'GATHER_FIELD_MPL : bad field type'
  icode = 1
  ! DEPENDS ON: ereport
  CALL ereport(routinename, icode, cmessage)
END IF

! We can only do gathers with global processor groups with this
! version currently
IF (proc_group /= gc_all_proc_group) THEN
  WRITE (cmessage, *) 'Can only gather in an all processor group'
  icode = 2
  ! DEPENDS ON: ereport
  CALL ereport(routinename, icode, cmessage)
END IF

! 1.0 Setup MPI datatypes. Can we use the same as last time?
!     last time round?

IF ((global_row_len  /=  old_global_row_len) .OR.                 &
    (global_rows     /=  old_global_rows   ) .OR.                 &
    (proc_group      /=  old_proc_group    ) .OR.                 &
    (gather_pe       /=  old_gather_pe     ) .OR.                 &
    (grid_type       /=  old_grid_type     ) .OR.                 &
    (halo_type       /=  old_halo_type     ) .OR.                 &
    (current_decomp_type  /=  old_decomp  )) THEN


  ! Need to create new send type
  ! Must get rid of old one if this isn't the first call
  IF (send_type /= -1234) THEN
    CALL MPL_type_free(send_type, info)
  END IF

  ! New type is a 2D block with a local row length + 2 * halo stride.
  ! This is used to describe the data size and layout of the local 
  ! data without including halos via an MPI datatype to simplify 
  ! communications and reduce copies. Will need to specify the real
  ! start of the data in the communications call.
  CALL MPL_type_vector(blsize(2, grid_type),  blsize(1, grid_type),   &
                       local_row_len,  MPL_REAL, send_type, info)
  CALL MPL_type_commit(send_type, info)

END IF  ! New MPI types

! 2.0 Setup the gathering counts and displacements
DO iproc = 0, nproc-1
  counts(iproc) = g_blsize(1, grid_type, iproc) *                    &
                  g_blsize(2, grid_type, iproc)

  IF (iproc == 0) THEN
    displacements(iproc) = 0
  ELSE
    displacements(iproc)=SUM(counts(0:iproc-1))
  END IF

END DO ! nprocs

! 3.0 do the exchange of data
CALL gc_get_communicator(my_comm, info)

! ipos is the position in the send array of the first data point
ipos = halosize(2,halo_type) * local_row_len +                      &
       halosize(1,halo_type) + 1
CALL MPL_gatherv(local_field(ipos), 1, send_type,                   &
                 global_field_copy, counts, displacements, MPL_REAL,&
                 gather_pe, my_comm, info)


! 4.0 Copy data out of global_field_copy into global_field
!     Note that global data does not include halos.
!     This copy takes the existing data layout (all of PE0 data,
!     followed by all of PE1 data, followed by ...) and copies into 
!     the conventional global field ordering.
IF  (mype == gather_pe) THEN
  ipos = 1 ! start of array
  DO iproc = 0, nproc-1
    DO j = 1, g_blsize(2, grid_type, iproc)
      igpos1 =   (g_datastart_f(2,grid_type,iproc) + j - 2) *        &
                  global_row_len - 1 +                               &
                  g_datastart_f(1,grid_type, iproc)
      DO i = 1, g_blsize(1, grid_type, iproc)

        igpos =  i + igpos1
        global_field(igpos) = global_field_copy(ipos)
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


RETURN
END SUBROUTINE gather_field_mpl
#endif
