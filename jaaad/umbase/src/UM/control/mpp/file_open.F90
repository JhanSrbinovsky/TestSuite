#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(RECON) \
 || defined(UTILIO) || defined(FLDIO) || defined(FLUXPROC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Parallel UM version of FILE_OPEN
!
! Subroutine Interface:
      SUBROUTINE FILE_OPEN(NFTIN,ENV,NENV,READ_WRITE,ENV_VAR,ERR)
#if defined(RECON)
      Use Rcf_Parvars_Mod
#endif

      IMPLICIT NONE
!
! Description:
!  This routine provides an interface to FILE_OPEN for the Parallel
!  Unified Model.
!
! Method:
!  The C FILE_OPEN is renamed OPEN_SINGLE under *DEF,MPP. This
!  routine causes OPEN_SINGLE to be called by PE 0 only.
!
! Current Code Owner: Paul Burton
!
! History:
!  Model    Date     Modification history from model version 3.5
!  version
!    3.5    5/1/95   New DECK created for the Parallel Unified
!                    Model. A.Dickinson + D.Salmond
!    4.3    7/3/96   Broadcast return code to all PEs  P.Burton
!LL  5.1    10/04/00 New reconfiguration support.
!    5.3   22/11/01   Enable MPP as the only option for
!                     small executables         E.Leung
!    5.5   25/04/03   Add defs FLUXPROC           E.Leung
!  6.0  17/09/03  Add def for new NEC opt section C96_1C. R Barnes
!
! Subroutine Arguments:

      INTEGER                                                           &
     & NFTIN                                                            &
                     !IN Unit number for I/O
     &,NENV                                                             &
                     !IN Length of ENV
     &,READ_WRITE                                                       &
                     !IN =0 read only, <> 0 read and write
     &,ENV_VAR                                                          &
                     !IN =0 file name stored in environment var
                     !   <>0 file name specified explicitly
     &,ERR           !OUT =0 file OPENED
                     !   <>0 file NOT OPENED because of error

      CHARACTER*(*)                                                     &
     & ENV           !IN Environment name or explicit file name

      INTEGER info
! Parameters and Common blocks

#if !defined(RECON)
#include "parvars.h"
#endif

! ------------------------------------------------------------------

      ERR=0
      IF (mype  ==  0) THEN    ! only PE 0 does any I/O
        CALL OPEN_SINGLE(NFTIN,ENV,NENV,READ_WRITE,ENV_VAR,ERR)
      ENDIF
      CALL GC_IBCAST(1,1,0,nproc,info,ERR)

      RETURN
      END SUBROUTINE FILE_OPEN

#endif
