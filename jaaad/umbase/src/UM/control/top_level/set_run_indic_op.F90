#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine: SET_RUN_INDIC_OP-----------------------------------
!LL
!LL  Purpose: Interface routine required to set RUN_INDIC_OP in dump
!LL           header. Called by INITIAL.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 6.1.5A
!LL
!LL  Model            Modification history:
!LL version  date
!LL  3.2   27/03/93  New routine required by dynamic allocation
!LL                  to interface with dump header array.
!LL  5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: C0
!LL
!LL  Project task: C0
!LL
!LL  External documentation: On-line UM document C1 - The top-level
!LL                          dynamic allocation
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE SET_RUN_INDIC_OP(                                      &
#include "argduma.h"
     &              ICODE,CMESSAGE)
!
!*----------------------------------------------------------------------
      IMPLICIT NONE
!
!  Subroutines called: none
!
!  Arguments
!
!  Configuration-dependent sizes and arrays
!
#include "parparm.h"
#include "typsize.h"
#include "typduma.h"
!
#include "chsunits.h"
#include "ihisto.h"
!
      INTEGER ICODE             ! Work - Internal return code
      CHARACTER*256 CMESSAGE    ! Work - Internal error message
!
!  Local variables: none
!
      A_FIXHD(6) = RUN_INDIC_OP
!
      RETURN
      END SUBROUTINE SET_RUN_INDIC_OP
#endif
