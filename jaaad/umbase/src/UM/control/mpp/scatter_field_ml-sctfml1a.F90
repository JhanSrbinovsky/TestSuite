#if defined(C96_1A) || defined(C96_1B)
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
! This version is the portable one. It justs scatters by looping
! over levels using SCATTER_FIELD
!
! Current Code Owner: Paul Selwood
!
! History:
!  Model    Date      Modification history:
!  version
!    6.2    25/06/04  New code - portable version
!                       Author: P.Selwood
!
! Subroutine Interface:
      Subroutine Scatter_Field_ML(                                      &
     &    local_field,    global_field,                                 &
     &    local_row_len,  local_rows,    local_levs,                    &
     &    global_row_len, global_rows,   global_levs,                   &
     &    pe_for_level,                                                 &
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

      Integer :: k       ! loop index  - levels
      Integer :: n_levs_on_proc(0 : nproc-1)   ! levs on PE

      Integer :: icode                         ! error code
      Character (Len=80) :: cmessage           ! error message
      Character (Len=*), Parameter  :: routinename='Scatter_Field_ML'

!-------------------------------------------------------
! 1) Loop over levels and scatter
!-------------------------------------------------------
      icode             = 0
      n_levs_on_proc(:) = 0

      Do k = 1, local_levs

        n_levs_on_proc(pe_for_level(k)) =                               &
     &                 n_levs_on_proc(pe_for_level(k)) + 1

! DEPENDS ON: scatter_field
        Call Scatter_Field(                                             &
     &               local_field(:,:,k),                                &
     &               global_field(:,:,n_levs_on_proc(pe_for_level(k))), &
     &               local_row_len,                                     &
     &               local_rows,                                        &
     &               global_row_len,                                    &
     &               global_rows,                                       &
     &               fld_type,                                          &
     &               halo_type,                                         &
     &               pe_for_level(k),                                   &
     &               gc_all_proc_group,                                 &
     &               icode,                                             &
     &               cmessage )

        If (icode /= 0) Then
! DEPENDS ON: ereport
          Call Ereport(routinename, icode, cmessage)
        End If
      End Do

      Return
      End Subroutine Scatter_Field_ML
#endif
