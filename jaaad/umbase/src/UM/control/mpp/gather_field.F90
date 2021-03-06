#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
 || defined(UTILIO) || defined(FLDIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Gathers a field from many processors to one processor
!
! Subroutine Interface:
SUBROUTINE gather_field(local_field,    global_field,     &
                        local_row_len,  local_rows,       &
                        global_row_len, global_rows,      &
                        grid_type,      halo_type,        &
                        gather_pe,      proc_group,       &
                        icode,          cmessage)

IMPLICIT NONE

! 
! Description: 
! Interface to potentially 2 methods of taking a field that is decomposed
! over a group of processors and gathering the data so a single processor
! contains the entire global field.
!
! Method:
! For C96_1A or 1B, only the GCOM method of gathering is available.
! For C96_1C GCOM and MPL versions are available. The choice is made
! on a user definable threshold that is given in the UMUI.
!
! Current Code Owner: Paul Selwood
!
! Subroutine Arguments:
INTEGER :: local_row_len   ! length of rows in local part of field
INTEGER :: local_rows      ! number of rows in local part of field
INTEGER :: global_row_len  ! length of rows in global field
INTEGER :: global_rows     ! number of rows in global field
INTEGER :: grid_type       ! type (p,u or v) of grid
INTEGER :: halo_type       ! halo type (hence width) of grid
INTEGER :: gather_pe       ! processor to gather global field to
INTEGER :: proc_group      ! group id of processors involved here
INTEGER :: icode           ! out return code

REAL    :: local_field(local_row_len*local_rows)
                                         ! local part of field
REAL    :: global_field(global_row_len*global_rows)
                                         ! (on pe gather_pe) global field

CHARACTER (LEN=80)  :: cmessage  ! error message



! Parameters and Common Blocks
#include "mppconf.h"
#include "parvars.h"

#if defined(C96_1A) || defined(C96_1B) \
 || defined(UTILIO) || defined(FLDIO) || defined(MAKEBC)

! Only GCOM version available
! DEPENDS ON: gather_field_gcom
CALL gather_field_gcom(local_field,    global_field,       &
                       local_row_len,  local_rows,         &
                       global_row_len, global_rows,        &
                       grid_type,      halo_type,          &
                       gather_pe,      proc_group,         &
                       icode,          cmessage)


#else

! Use GCOM version if number of processors less than or equal to 
! threshold
IF (nproc <= gcom_coll_limit) THEN 

! DEPENDS ON: gather_field_gcom
  CALL gather_field_gcom(local_field,    global_field,     &
                         local_row_len,  local_rows,       &
                         global_row_len, global_rows,      &
                         grid_type,      halo_type,        &
                         gather_pe,      proc_group,       &
                         icode,          cmessage)


ELSE

! DEPENDS ON: gather_field_mpl
  CALL gather_field_mpl(local_field,    global_field,      &
                        local_row_len,  local_rows,        &
                        global_row_len, global_rows,       &
                        grid_type,      halo_type,         &
                        gather_pe,      proc_group,        &
                        icode,          cmessage)


END IF
#endif

RETURN
END SUBROUTINE gather_field
#endif
