#if defined(OASIS3)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE OASIS3_PUT64(data_64,l_cols,l_rows,ierr64,ft           &
     ,s_cols,s_rows,put_step,halo_type_no_halo,gc_all_proc_group,mype)
  !
  !  Description: Buffer routine between OASIS 32 bit put
  !               (send) operations and UM 64 bit data.
  !               Also performs gather if the option to 
  !               communicate with OASIS through a master
  !               PE is selected. 
  !
  !  Author:      R. Hill
  !  Current Code Owner : R. Hill
  !
  !===========================================================
  
  USE oasis3_atm_data_mod
  USE oasis3_atmos_init

  IMPLICIT NONE

#include "cntlatm.h"

  INTEGER :: l_cols, l_rows     ! Local domain sizes
  INTEGER :: s_cols, s_rows     ! Send buffer sizes
  REAL    :: data_64(l_cols,l_rows,1)

  INTEGER :: I,J, N
  INTEGER :: info64
  INTEGER :: halo_type_no_halo
  INTEGER :: gc_all_proc_group
  INTEGER :: mype

  INTEGER :: ierr64, ft
  INTEGER :: ierror
  INTEGER :: info
  CHARACTER*(80) :: cmessage

  LOGICAL :: put_step ! Indicates genuine coupling step
  ! for performance purposes. 

  REAL ,ALLOCATABLE ,DIMENSION(:,:)  :: data_send

  ALLOCATE (DATA_SEND(s_cols,s_rows))

  IF (put_step) THEN

    IF (L_COUPLE_MASTER) THEN
      ! We now have to gather the data on our master PE
      ! But we'll only bother if this is a genuine coupling
      ! timestep - i.e. iff the data will actually be used. 

      ! 1st: Gather all onto 1 PE in the send buffer.
      ! DEPENDS ON: gather_field
      CALL GATHER_FIELD(data_64(1,1,1),data_send,                 &
                  l_cols,                                         &
                  l_rows,                                         &
                  s_cols,                                         &
                  s_rows,                                         &
                  ft,halo_type_no_halo,                           &
                  OASIS_CNTLPE,GC_ALL_PROC_GROUP,info,cmessage)

    ELSE
      ! Move the data to send direct to our send buffer.
      data_send(1:l_cols,1:l_rows)=data_64(1:s_cols,1:s_rows,1)

    END IF

  ENDIF  ! put_step = true.

  ! If we're putting through a master only, call the put
  ! routine on the master PE. Otherwise everybody calls it.
  IF ((L_COUPLE_MASTER.AND.MYPE == OASIS_CNTLPE).OR.                &
          (.NOT.L_COUPLE_MASTER)) THEN

    CALL OASIS3_PUT32(data_send,s_cols,s_rows,ierr64)

  END IF

  DEALLOCATE (DATA_SEND)

  RETURN
END SUBROUTINE OASIS3_PUT64
#endif
