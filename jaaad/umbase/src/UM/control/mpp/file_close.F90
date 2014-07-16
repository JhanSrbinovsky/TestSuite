#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
 || defined(RECON) \
 || defined(UTILIO) || defined(FLDIO)

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Parallel UM version of FILE_CLOSE
!
! Subroutine Interface:
      SUBROUTINE FILE_CLOSE(NFTIN,ENV,NENV,ENV_VAR,DELETE,ERR)

#if defined(RECON)
      Use Rcf_Parvars_Mod
#endif
      IMPLICIT NONE
!
! Description:
!  This routine provides an interface to FILE_CLOSE for the Parallel
!  Unified Model.
!
! Method:
!  The C FILE_CLOSE is renamed CLOSE_SINGLE under *DEF,MPP. This
!  routine causes CLOSE_SINGLE to be called by PE 0 only.
!
! Current Code Owner: Paul Burton
!
! History:
!  Model    Date     Modification history from model version 3.5
!  version
!    3.5    5/1/95   New DECK created for the Parallel Unified
!                    Model. A.Dickinson + D.Salmond
!    4.2    15/10/96 Added ERR argument to match non-MPP code
!                    P.Burton
!LL  5.1    10/04/00 New reconfiguration support.
!    5.3    07/12/01 Initialise ERR. Adam Clayton
!
!    5.3   22/11/01   Enable MPP as the only option for
!                     small executables         E.Leung
!    6.0   17/09/03   Add def for new NEC opt section C96_1C. R Barnes
!    6.0    Close and flush buffers when closing file.
!           JC Rioual, P. Selwood.
!    6.0   02/07/03   Add FLDIO defs      E.Leung
! Subroutine Arguments:

      INTEGER                                                           &
     & NFTIN                                                            &
                     !IN Unit number for I/O
     &,NENV                                                             &
                     !IN Length of ENV
     &,DELETE                                                           &
                     !IN =0, do not delete file
                     !   <>0, delete file
     &,ENV_VAR                                                          &
                     !IN =0 file name stored in environment var
                     !   <>0 file name specified explicitly
     &,ERR      ! OUT return code (0=no error)

      CHARACTER*(*)                                                     &
     & ENV          !IN Environment name of file

! Parameters and Common blocks

#if !defined(RECON)
#include "parvars.h"
#endif

! ------------------------------------------------------------------

      ERR = 0
      IF (mype  ==  0) THEN    ! only PE 0 does any I/O
! Flush buffers if not small execs
#if !defined(UTILIO) && !defined(RECON)
! DEPENDS ON: test_and_flush_pp
        Call test_and_flush_pp(nftin)
#endif

        CALL CLOSE_SINGLE(NFTIN,ENV,NENV,ENV_VAR,DELETE,ERR)
      ENDIF

      RETURN
      END SUBROUTINE FILE_CLOSE

#endif
