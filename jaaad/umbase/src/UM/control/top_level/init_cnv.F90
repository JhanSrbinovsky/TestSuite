#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Ensure conv.cloud cover & liquid water path zero if no conv.cloud.
!
! Subroutine Interface:
      SUBROUTINE INIT_CNV(                                              &
#include "argd1.h"
#include "argptra.h"
     &                    ICODE,CMESSAGE)

      IMPLICIT NONE
!
! Description:
!   Sets conv.cloud cover and liquid water path to zero at timestep 0
!   if conv.cloud base and top are zero (no conv.cloud is present).
!
! Method:
!   For full model field tests conv.cloud base and top, and if either
!   is zero sets conv.cloud cover and liquid water path to zero.
!   This consistency check at timestep 0 is needed as interpolation in
!   the reconfiguration can give rise to inconsistent fields.
!
! Current Code Owner: R.T.H.Barnes (FR)
!
! History:
! Version  Date         Comment
! -------  ----         -------
!  3.3  25/02/94  New routine. R.T.H.Barnes.
!  4.4  26/09/97  Changes to initialisation of CCA to allow for
!                 radiative representation of cld anvils. Julie Gregory
!  5.0  30/06/99  Change loop over horizontal points to use C-P C-grid
!                 dynamics variables instead of P_FIELD. M.L.Gallani
!LL  5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: <appropriate code>
! System Task:              <appropriate code>
!
! Declarations: these are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Global variables (*CALLed common blocks etc.):
#include "cmaxsize.h"
#include "parparm.h"
#include "typsize.h"
#include "typd1.h"
#include "typptra.h"

! Subroutine arguments
!   Scalar arguments with Intent(In):
!   Array  arguments with Intent(In):
!   Scalar arguments with Intent(InOut):
!   Array  arguments with Intent(InOut):
!   Scalar arguments with Intent(Out):
!   Array  arguments with Intent(Out):

!   ErrorStatus <Delete this & the next 2 lines if ErrorStatus not used>
      INTEGER      ICODE                ! Error flag (0 = OK)
      CHARACTER*80 CMESSAGE             ! Error message if ICODE>0

! Local parameters:
! Local scalars:
      INTEGER   I,K,IJ,J ! Loop counters over ROWS,ROW_LENGTH,N_CCA_LEV

! Local dynamic arrays:

! Function & Subroutine calls:
!     External - NONE

!- End of header

! 1.0 Ensure that conv.cloud cover and liquid water path are zero
!      when there is zero conv.cloud base and top.
      WRITE(6,*)'INIT_CNV:resets conv.cld cover zero for base/top zero'
      DO  J = 1,ROWS
        DO i=1, row_length
          ij = i + (j-1)*row_length
          IF ( ID1(JCCB+IJ-1) == 0 .OR. ID1(JCCT+IJ-1) == 0 ) THEN
            DO  K = 1,N_CCA_LEV
              D1(JCCA(K)+IJ-1) = 0.0
            END DO
            D1(JCCLWP+IJ-1) = 0.0
          END IF
        END DO ! i
      END DO ! j

      RETURN
      END SUBROUTINE INIT_CNV
#endif
