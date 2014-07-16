#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Constructs STASH index array
!
! Subroutine Interface:

      SUBROUTINE SINDX(NRECS)
      IMPLICIT NONE

! Description:
!   The STASH list (LIST_S) is input to this routine. The output from
!   the routine is the STASH index (INDX_S) consistent with LIST_S.
!   Called by STPROC, DUPLIC.
!
! Method:
!   Before this routine is executed, LIST_S has been ordered by
!   model, section, item, input code. In this routine, INDX_S is
!   set up as follows:
!   INDX_S(1,m,s,i)= position of 1st. occurrence of m,s,i in LIST_S
!   INDX_S(2,m,s,i)= no. of occurrences of m,s,i in LIST_S
!
! Current code owner:  S.J.Swarbrick
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Mar. 95    Original code.  S.J.Swarbrick
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
!  Global variables:

#include "csubmodl.h"
#include "version.h"
#include "cstash.h"
#include "stextend.h"
#include "stparam.h"

!  Subroutine Arguments:

!    Scalar Argument with intent(in):

      INTEGER NRECS                 ! No. of records in LIST_S

!  Local variables:

      INTEGER II
      INTEGER IITM
      INTEGER IM
      INTEGER IREC
      INTEGER IS
      INTEGER ISEC
      INTEGER LSTART
      INTEGER MODL

!  Local arrays:

!- End of Header -------------------------------------------------------


! Initialise Index array

      DO IM    = 1,N_INTERNAL_MODEL_MAX
        DO IS  = 0,NSECTP
          DO II= 1,NITEMP
           INDX_S(1,IM,IS,II)=0
           INDX_S(2,IM,IS,II)=0
          END DO
        END DO
      END DO

! Set up index

      IF(NRECS >= 1) THEN

      MODL=LIST_S(st_model_code  ,1)
      ISEC=LIST_S(st_sect_no_code,1)
      IITM=LIST_S(st_item_code   ,1)

      LSTART=1

      INDX_S(2,MODL,ISEC,IITM) = 1
      INDX_S(1,MODL,ISEC,IITM) = 1

      IF(NRECS >= 2) THEN            ! More than one record in LIST_S

        DO IREC=2,NRECS
          IF((LIST_S(st_model_code  ,IREC) == MODL).AND.                &
                                                         ! Same model,
     &       (LIST_S(st_sect_no_code,IREC) == ISEC).AND.                &
                                                         ! sec,item,
     &       (LIST_S(st_item_code   ,IREC) == IITM))THEN ! as before

            LSTART=LSTART+1

            INDX_S(2,MODL,ISEC,IITM) = INDX_S(2,MODL,ISEC,IITM)+1

          ELSE   ! New model, section, item

            MODL=LIST_S(st_model_code  ,IREC)
            ISEC=LIST_S(st_sect_no_code,IREC)
            IITM=LIST_S(st_item_code   ,IREC)

            LSTART=LSTART+1

            INDX_S(1,MODL,ISEC,IITM) = LSTART
            INDX_S(2,MODL,ISEC,IITM) = 1

          END IF
        END DO

      ELSE      ! Only one record

        INDX_S(1,LIST_S(st_model_code  ,1),                             &
     &           LIST_S(st_sect_no_code,1),                             &
     &           LIST_S(st_item_code   ,1)) = 1

        INDX_S(2,LIST_S(st_model_code  ,1),                             &
     &           LIST_S(st_sect_no_code,1),                             &
     &           LIST_S(st_item_code   ,1)) = 1

      END IF
      END IF

      RETURN
      END SUBROUTINE SINDX
#endif
