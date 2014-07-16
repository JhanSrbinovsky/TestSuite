#if defined(C96_1A) || defined(MAKEBC)
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
!  This routine uses multiple calls to GATHER_FIELD to achieve
!  the required result.
!
! Current Code Owner: Paul Selwood
!
! History:
!  Model    Date      Modification history:
!  version
!    6.2    06/12/05  New code - alternative to optimised versions.
!                       Author: P.Selwood
!
! Subroutine Interface:
      Subroutine Gather_Field_ML(                                       &
     &    local_field,    global_field,                                 &
     &    local_row_len,  local_rows,   local_levs,                     &
     &    global_row_len, global_rows,  global_levs,                    &
     &    pe_for_level,   local_level,                                  &
     &    fld_type,       halo_type)

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
      Integer            :: k                ! loop index  - levels
      Integer            :: icode            ! error code
      Integer            :: n_levs_on_proc(0 : nproc-1)
      Character (Len=80) :: cmessage         ! error message
      Character (Len=*), Parameter  :: routinename = 'Gather_Field_ML'


!-------------------------------------------------------
! 1) Loop over levels gathering each separately.
!-------------------------------------------------------
      icode             = 0
      n_levs_on_proc(:) = 0

      Do k = 1, local_levs
        n_levs_on_proc(pe_for_level(k)) =                               &
     &                 n_levs_on_proc(pe_for_level(k)) + 1

! DEPENDS ON: gather_field
        Call Gather_Field(                                              &
     &              local_field(:,:,k),                                 &
     &              global_field(:,:,n_levs_on_proc(pe_for_level(k))),  &
     &              local_row_len,                                      &
     &              local_rows,                                         &
     &              global_row_len,                                     &
     &              global_rows,                                        &
     &              fld_type,                                           &
     &              halo_type,                                          &
     &              pe_for_level(k),                                    &
     &              gc_all_proc_group,                                  &
     &              icode,                                              &
     &              cmessage )

        If (icode /= 0) Then
! DEPENDS ON: ereport
          Call Ereport(routinename, icode, cmessage)
        End If
      End Do

      Return
      End Subroutine Gather_Field_ML
#endif
