#if defined(OASIS3)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE OASIS3_GET64(data_64,l_cols,l_rows,oasis_info,ft &
   ,r_cols,r_rows,get_step,halo_type_no_halo   &
   ,gc_all_proc_group,mype)  

  !
  !  Description: Buffer routine between OASIS 32 bit get
  !               (receive) operations and UM 64 bit data.
  !               Also performs scatter if the option to 
  !               communicate with OASIS through a master
  !               PE is selected. 
  !
  !  Author:      R. Hill
  !  Current Code Owner : R. Hill
  !===========================================================
  USE oasis3_atm_data_mod
  USE oasis3_atmos_init

  IMPLICIT NONE

#include "cntlatm.h"

  INTEGER :: l_cols, l_rows  ! Local domain size
  INTEGER :: r_cols, r_rows  ! Receive buffer size
  INTEGER :: oasis_info

  REAL    :: data_64(l_cols,l_rows,1)

  INTEGER :: ft,halo_type_no_halo
  INTEGER :: gc_all_proc_group
  INTEGER :: mype

  INTEGER :: info
  REAL ,ALLOCATABLE ,DIMENSION(:,:) :: data_recv

  LOGICAL :: get_step
  CHARACTER*(80) :: cmessage

  ! If we're coupling through a master PE we need to allocate
  ! the array size for incoming data accordingly
  ALLOCATE (DATA_RECV(r_cols,r_rows))
  DATA_RECV(:,:) = 0.0

  IF ((L_COUPLE_MASTER.AND.MYPE == oasis_CNTLPE).OR.            &
     (.NOT.L_COUPLE_MASTER)) THEN
    CALL OASIS3_GET32(DATA_RECV,r_cols,r_rows, oasis_info)
  END IF

  ! Having just got the incoming data, we may have to scatter
  ! it to our component decomposition if the put/get
  ! was done through a master PE.

  IF (GET_STEP) THEN

    IF (L_COUPLE_MASTER) THEN
      ! In this case the master PE must broadcast any
      ! oasis return code to all PEs since they  may need this
      ! info for subsequent operations (e.g to decide
      ! whether to do polar row averaging etc etc).

      ! Now we can scatter this data according to
      ! the true MPP decomposition.
      ! DEPENDS ON: scatter_field
      CALL SCATTER_FIELD(data_64,DATA_RECV,                &
         l_cols,                                               &
         l_rows,                                               &
         r_cols,                                               &
         r_rows,                                               &
         ft,halo_type_no_halo,                                 &
         OASIS_CNTLPE,GC_ALL_PROC_GROUP,info,cmessage)
    ELSE

      ! Nothing to do here beyond copying the receive buffer 
      ! straight to our local domain. 
      data_64(1:l_cols,1:l_rows,1)=DATA_RECV(1:r_cols,1:r_rows)

    END IF
  END IF  ! GET_STEP = true

  DEALLOCATE (DATA_RECV)
  ! Now we're ready to proceed with our incoming coupling data
  ! distributed to local domains. 

  RETURN
END SUBROUTINE OASIS3_GET64
#endif
