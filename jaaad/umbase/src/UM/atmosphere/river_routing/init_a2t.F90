#if defined(ATMOS) && defined(A26_1A) && !defined(OCEAN)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      SUBROUTINE INIT_A2T (                                             &
#include "argd1.h"
#include "argsts.h"
#include "argduma.h"
#include "argptra.h"
#include "argatcpl.h"
     &                      ICODE,CMESSAGE )

!  Routine: INIT_A2T ------------------------------------------------
!
!  Purpose: Initialises the Lat and Long values of ATMOS for
!           regridding to TRIP river routing grid (based on INITA2O)
!           NB this will need to be extended to pick up the values
!           for the river routing grid from the ancil header of one
!           file later so will leave in redundant code for guidance.
!
!  Author:   C.Bunton         Date: 28 Feb 2003
!
!  Model            Modification history from model version 5.5:
! version  date
!   6.0   12/09/03  Change DEF from A20 to A26. D. Robinson
!   6.0   31/10/03  Remove CONTROL def requirement. S.D.Mullerworth
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!  -------------------------------------------------------------------

      IMPLICIT NONE
!
#include "csubmodl.h"
#include "cntlatm.h"
#include "parparm.h"
#include "typsize.h"
#include "typd1.h"
#include "typsts.h"
#include "typduma.h"
#include "typptra.h"
#include "typatcpl.h"
!
      INTEGER ICODE                    ! OUT - Error return code
      CHARACTER*(*) CMESSAGE           ! OUT - Error return message
!*------------------------------------------------------------------
!  Common blocks
!
!  Local variables
!
      INTEGER                                                           &
     &       I,J

      ICODE=0

! Pick up the ATMOS x coords from the dump header
      XUA(0)=A_REALHD(4)-0.5*A_REALHD(1)
      DO I=1,AOCPL_ROW_LENGTH
        XPA(I)=A_REALHD(4)+(I-1)*A_REALHD(1)
        XVA(I)=A_REALHD(4)+(I-1)*A_REALHD(1)
        XUA(I)=A_REALHD(4)+(I-0.5)*A_REALHD(1)
      ENDDO
      XPA(AOCPL_ROW_LENGTH+1)=A_REALHD(4)+AOCPL_ROW_LENGTH*A_REALHD(1)
      XVA(AOCPL_ROW_LENGTH+1)=A_REALHD(4)+AOCPL_ROW_LENGTH*A_REALHD(1)

! Pick up the ATMOS y coords from the dump header
      DO J=1,AOCPL_P_ROWS
! JCT  S->N now  YTA(J)=A_REALHD(3)-(J-1)*A_REALHD(2)
        YPA(J)=A_REALHD(3)+(J-1)*A_REALHD(2)
        YUA(J)=A_REALHD(3)+(J-1)*A_REALHD(2)
      ENDDO
      DO J=0,AOCPL_P_ROWS
! JCT  S->N now  YUA(J)=A_REALHD(3)-(J-0.5)*A_REALHD(2)
        YVA(J)=A_REALHD(3)+(J-0.5)*A_REALHD(2)
      ENDDO
!
      RETURN
!----------------------------------------------------------------------
      END SUBROUTINE INIT_A2T
#endif
