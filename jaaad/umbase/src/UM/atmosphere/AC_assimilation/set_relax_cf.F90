#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE SET_RELAX_CF -------------------------------------------
!LL
!LL  Purpose :
!LL
!LL  For use on Cray Y-MP
!LL  For Cray - Global  ; Enable defs GLOBAL
!LL
!LL  Created from INITAC. D Robinson 2/3/92
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.2  8/7/93      Eliminate QA FORTRAN complaints    S Bell
!    4.2 25/11/96: t3e mods Stuart Bell
!LL
!LL  4.1  4/09/96:  Port to CRAY T3E  Deborah Salmond
!LL  4.4 28/11/97:  Correct unsafe code and remove bounds check
!LL                 messages. Rick Rawlins
!LL  4.5 26/01/98:  Correct unsafe code to avoid oob access to
!LL                 RELAX_CF                          P.Burton
!   5.2 12/12/00:  change attop,base to elements of at_extremity
!                  -amend glsize for extra dimension  B Macpherson
!   5.3 12/07/01:  correct use of at_extremity for S->N ND grid order
!                  disable latitudinal variation in global model
!                  B Macpherson
!   6.0 19/06/03:  Remove non-MPP parts of code. T. White
!  for ease, set nudging coeff to S.Hem. value at all rows
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
!*L  Arguments:---------------------------------------------------
      SUBROUTINE SET_RELAX_CF (JGROUP, N_ROWS_LOCAL,                    &
     &                       RELAX_CF_LOCAL, LWIND,                     &
     &           TIMESTEP, TIMESTEP_NO, ITER_NO,                        &
     &           ICODE, CMESSAGE)

!     Do not multi-task this routine.
!FPP$ NOCONCUR R

      IMPLICIT NONE
#include "acparm.h"

      INTEGER JGROUP               !IN identifier for group
      INTEGER N_ROWS_LOCAL
      REAL    RELAX_CF_LOCAL(N_ROWS_LOCAL)
      LOGICAL LWIND                !IN switch to identify wind grid
      REAL    TIMESTEP             !IN timestep
      INTEGER TIMESTEP_NO          !IN timestep number
      INTEGER ITER_NO              !IN iteration count
      INTEGER ICODE                !OUT error code and message
      CHARACTER*256 CMESSAGE
!*---------------------------------------------------------------------
!     AC Comdecks
#include "comacp.h"
#include "comacdg.h"
#include "comag.h"
#include "commg.h"
#include "parvars.h"
#include "mppac.h"
      INTEGER N_ROWS
      REAL RELAX_CF(P_ROWS_MAX)
!L---------------------------------------------------------------------
!     UM Constant Comdecks
#include "c_pi.h"
#include "c_r_cp.h"
!-----------------------------------------------------------------------
      INTEGER JROW
#if defined(GLOBAL)
      INTEGER IR1,IR2,IR3,IR4,IR5,INR1,INR2,INR3,INR4,INR5
      REAL TROPLATN, TROPLATS
      REAL GRADN,GRADS,NUDGEC
      REAL NUDGE_NH
      REAL NUDGE_TR
      REAL NUDGE_SH
#else
      REAL NUDGE_LAM
#endif
!-----------------------------------------------------------------------
!     RELAXATION COEFFICIENTS
!     -----------------------
!
!     these are now input as nudging (N in 3.23 of TN27) and
!     converted to relaxation coefficients in RELAX_CF (using 3.23)

      N_ROWS=glsize(2,fld_type_p)
#if defined(GLOBAL)


      NUDGE_NH  = DEF_NUDGE_NH (GROUP_INDEX(JGROUP))
      NUDGE_TR  = DEF_NUDGE_TR (GROUP_INDEX(JGROUP))
      NUDGE_SH  = DEF_NUDGE_SH (GROUP_INDEX(JGROUP))

      DO JROW=1,N_ROWS
        RELAX_CF(JROW) = NUDGE_SH
      ENDDO
#else

      NUDGE_LAM = DEF_NUDGE_LAM(GROUP_INDEX(JGROUP))
      DO JROW=1,N_ROWS
        RELAX_CF(JROW) = NUDGE_LAM
      ENDDO
#endif

!     convert the nudging coefficients to relaxation coeffs
!     using eqn 3.23 of TN27 (lamda=N*dT/(1+N*dT) )

      DO JROW=1,N_ROWS

        RELAX_CF(JROW) =                                                &
     &  RELAX_CF(JROW)*TIMESTEP/(1+RELAX_CF(JROW)*TIMESTEP)

      ENDDO
      IF (N_ROWS+1  <=  P_ROWS_MAX) THEN
        RELAX_CF(N_ROWS+1)=RELAX_CF(N_ROWS)
      ENDIF

      IF (LDIAGAC .AND. TIMESTEP_NO == 1 .AND. ITER_NO == 1) THEN
      if(mype == 0)then

        PRINT '(/,A,I3)', ' RELAX_CF for Group No',JGROUP
        PRINT '(10F10.3)', (RELAX_CF(JROW),JROW=1,N_ROWS)

      ENDIF
      endif
      if(at_extremity(PSouth)) then
         relax_cf_local(1)=relax_cf(datastart(2))
      else
         relax_cf_local(1)=relax_cf(datastart(2)-1)
      endif
      do jrow=2,n_rows_local
      relax_cf_local(jrow)=relax_cf(jrow+datastart(2)-2)
      enddo
      if(at_extremity(PNorth)) then
        relax_cf_local(n_rows_local)=                                   &
     &    relax_cf(n_rows_local+datastart(2)-3)
      else
        relax_cf_local(n_rows_local)=                                   &
     &    relax_cf(n_rows_local+datastart(2)-2)
      endif

      RETURN
      END SUBROUTINE SET_RELAX_CF
#endif
