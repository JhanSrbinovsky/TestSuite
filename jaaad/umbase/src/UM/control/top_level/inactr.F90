#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Add inactive records to STASH list, when space is required
!
! Subroutine Interface:

      SUBROUTINE INACTR(                                                &
#include "argppx.h"
     &                   NRECS,ErrorStatus,CMESSAGE)
      IMPLICIT NONE

! Description:
!
! Method:
!
! Current code owner:  S.J.Swarbrick
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Mar. 95    Original code.  S.J.Swarbrick
!   4.0     Oct. 95                    S.J.Swarbrick
!   4.1     Apr. 96      Add ErrorStatus arguments & general
!                         improvements             S.J.Swarbrick
!   4.5     Jul. 98    Clarify error message. S.D.Mullerworth
!LL  5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!    6.0  15/12/03  Update IOPN inline with
!                   stashmaster 30 digit option codes.    M.Hughes
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Global variables:
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "parparm.h"
#include "typsize.h"
#include "cstash.h"
#include "stextend.h"
#include "model.h"
#include "stparam.h"

! Subroutine arguments:
!   Scalar arguments with intent(in):
      INTEGER NRECS
!   Scalar arguments with intent(out):
      CHARACTER*80 CMESSAGE

! ErrorStatus:
      INTEGER ErrorStatus

! Local scalars
      LOGICAL LMASK
      LOGICAL LIMPLIED
      LOGICAL LDUM
      INTEGER I
      INTEGER ITEM
      INTEGER ISEC
      INTEGER Im_ident
      INTEGER Im_index
      INTEGER LBVC

! Function and subroutine calls:
      INTEGER  EXPPXI
      EXTERNAL EXPPXI,IMPLIED,TSTMSK,ADDIN

!- End of Header ------------------------------------------------------


      DO Im_ident=1,N_INTERNAL_MODEL_MAX
      Im_index   = INTERNAL_MODEL_INDEX(Im_ident)
      IF (Im_index >  0) THEN
       DO ISEC   =0,PPXREF_SECTIONS
       DO ITEM   =1,PPXREF_ITEMS

          IF ((INDX_S(2,Im_ident,ISEC,ITEM) == 0).AND.                  &
     &        (PPXPTR(  Im_index,ISEC,ITEM) /= 0))    THEN

! No requests for this diag; check whether it is an implied diag.,
!                        or one for which space is always required

! DEPENDS ON: exppxi
        VMSK    = EXPPXI(Im_ident   ,ISEC   ,ITEM  ,ppx_version_mask,   &
#include "argppx.h"
     &                                        ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        ISPACE  = EXPPXI(Im_ident   ,ISEC   ,ITEM  ,ppx_space_code  ,   &
#include "argppx.h"
     &                                        ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        ILEV    = EXPPXI(Im_ident   ,ISEC   ,ITEM  ,ppx_lv_code     ,   &
#include "argppx.h"
     &                                        ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IBOT    = EXPPXI(Im_ident   ,ISEC   ,ITEM  ,ppx_lb_code     ,   &
#include "argppx.h"
     &                                        ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        ITOP    = EXPPXI(Im_ident   ,ISEC   ,ITEM  ,ppx_lt_code     ,   &
#include "argppx.h"
     &                                        ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IFLAG   = EXPPXI(Im_ident   ,ISEC   ,ITEM  ,ppx_lev_flag    ,   &
#include "argppx.h"
     &                                        ErrorStatus,CMESSAGE)
        DO I=1,6

! DEPENDS ON: exppxi
        IOPN(I) = EXPPXI(Im_ident   ,ISEC   ,ITEM  ,ppx_opt_code+I-1,   &
#include "argppx.h"
     &                                        ErrorStatus,CMESSAGE)
        END DO
! DEPENDS ON: exppxi
        LBVC    = EXPPXI(Im_ident   ,ISEC   ,ITEM  ,ppx_lbvc_code   ,   &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)

! Check whether this diag is implied by any of
!       the other diags in LIST_S

            IF(ISPACE == 6) THEN
! DEPENDS ON: implied
              CALL IMPLIED                                              &
     &       (Im_ident,ISEC,ITEM,LIMPLIED,ErrorStatus,CMESSAGE)
            ELSE
              LIMPLIED=.FALSE.
            END IF

            IF((ISPACE == 1).OR.LIMPLIED)THEN
! Check availability of diag
! DEPENDS ON: tstmsk
              CALL TSTMSK                                               &
     &       (Im_ident,ISEC,LMASK,LDUM,ErrorStatus,CMESSAGE)
              IF(LMASK.AND.(NRECS <  NRECDP)) THEN
! Diag to be included
                NRECS=NRECS+1
! Add diag to LIST_S
! DEPENDS ON: addin
                CALL ADDIN                                              &
     &         (NRECS,ITEM,ISEC,Im_ident,LBVC,ErrorStatus,CMESSAGE)
              ELSE IF (NRECS >= NRECDP) THEN
                WRITE(6,*)'ERROR, INACTR: TOO MANY S_LIST ENTRIES '     &
     &        ,'CANNOT ADD ENTRIES FOR ARRAYS REQUIRED BY THE MODEL'
              END IF
            END IF

          END IF   ! INDX_S

       END DO    ! Items
       END DO    ! Sections
      END IF    ! Im_index>0
      END DO    ! Models

      RETURN
      END SUBROUTINE INACTR


!+Find whether ST_list entry (Im_ident,ISEC,ITEM) is an implied diag
! Subroutine Interface:



!+Add diagnostic to the STASH list (LIST_S)
! Subroutine Interface:

#endif
