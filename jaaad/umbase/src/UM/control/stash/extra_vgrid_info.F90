#if defined(C84_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      SUBROUTINE EXTRA_VGRID_INFO(BUF_EXTRA,LEN_EXTRA                   &
     &               ,EXTRA_ARR,VECTOR_SIZE,vector_code)

      IMPLICIT NONE

#include "parvars.h"

      INTEGER :: LEN_EXTRA     ! IN len extra data in this vector
      INTEGER :: VECTOR_SIZE   ! OUT the size of  vector being created
      INTEGER :: VECTOR_CODE   ! IN indicates type of vector

      REAL :: BUF_EXTRA(LEN_EXTRA+1)  ! OUT the EXTRA_DATA vector:
                                      ! initial integer code +
                                      ! the extra data.
      REAL :: EXTRA_ARR (LEN_EXTRA)   ! IN incoming data to be added
                                      ! to the EXTRA_DATA vector
      INTEGER :: I  ! Local index


      ! Set up the initial Integer value in the extra data vector
! DEPENDS ON: stuff_int
      CALL STUFF_INT(BUF_EXTRA(1),(1000*LEN_EXTRA)+VECTOR_CODE)

      ! One could temporarily ise this version for test reasons
      !BUF_EXTRA(1) = 1000.*LEN_EXTRA+VECTOR_CODE


      ! Now at last we can stuff in the extra data itself.
      DO I = 1, LEN_EXTRA
         BUF_EXTRA(I+1) = EXTRA_ARR(I)
      ENDDO

      VECTOR_SIZE = LEN_EXTRA+1

      RETURN
      END SUBROUTINE EXTRA_VGRID_INFO
#endif
