#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(RECON) \
 || defined(UTILIO) || defined(FLDIO) || defined(FLUXPROC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Parallel UM version of SETPOS
!
! Subroutine Interface:
      SUBROUTINE SETPOS(NFT,IPOS,ICODE)

#if defined(RECON)
      Use Rcf_Parvars_Mod
#endif
      IMPLICIT NONE
!
! Description:
!  This routine provides an interface to SETPOS for the Parallel
!  Unified Model.
!
! Method:
!  The C SETPOS is renamed SETPOS_SINGLE under *DEF,MPP. This
!  routine causes SETPOS_SINGLE to be called by PE 0 only.
!
! Current Code Owner: Paul Burton
!
! History:
!  Model    Date     Modification history from model version 3.5
!  version
!    3.5    5/1/95   New DECK created for the Parallel Unified
!                    Model. A.Dickinson + D.Salmond
!    4.1    21/05/96   Added ICODE argument   P.Burton
!    5.1    10/04/00 New reconfiguration support.
!    5.3    17/10/01 ICODE initialised. Adam Clayton
!    5.3    22/11/01 Enable MPP as the only option for
!                    small executables         E.Leung
!    5.5    25/04/03 Add defs FLUXPROC           E.Leung
!  6.0  17/09/03  Add def for new NEC opt section C96_1C. R Barnes
!
! Subroutine Arguments:

      INTEGER NFT,                                                      &
                     ! IN : Fortran unit number
     &        IPOS                                                      &
                     ! IN : Position in file
     &,       ICODE  ! OUT : Return code

! Parameters and Common blocks

#if !defined(RECON)
#include "parvars.h"
#endif

! ------------------------------------------------------------------

      ICODE = 0

      IF (mype  ==  0) THEN    ! only PE 0 does any I/O
        CALL SETPOS_SINGLE(NFT,IPOS,ICODE)
      ENDIF

      RETURN
      END SUBROUTINE SETPOS

#endif
