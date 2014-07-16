#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine: JOBCTL ------------------------------------------------
!LL
!LL  Purpose: Outputs job release request to output processing server
!LL           task at selected timesteps.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Author:   T.C.Johns
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL  3.1  3/02/93 : added comdeck CHSUNITS to define NUNITS for i/o
!LL
!LL   3.1  05/02/93    Portable Fortran unit no assigns
!LL                    Author: A. Dickinson    Reviewer: R. Stratton
!LL 3.4  16/6/94 : Change CHARACTER*(*) to CHARACTER*(80) N.Farnon
!LL  3.5  25/04/95  Sub_models stage 1: make generic. RTHBarnes.
!LL  4.3  31-01-97  Parallelise writes to the pipe.   LCWiles
!LL 4.4  09/10/97   Change the closes on unit 8 to flushes
!LL                   Author: Bob Carruthers, Cray Research
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!LL  6.2  05/05/06  Print out message when job released. D. Robinson
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: C0
!LL
!LL  Project task: C0
!LL
!LL  External documentation:
!LL    On-line UM document C0 - The top-level control system
!LL    On-line UM document - Automated output processing system
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE JOBCTL(I_AO,ICODE,CMESSAGE)
!
!*----------------------------------------------------------------------
      IMPLICIT NONE
!
      INTEGER I_AO               ! Internal model indicator
!                                ! 1 - Atmosphere; 2 - Ocean
      INTEGER ICODE              ! Error return code
      CHARACTER*(80) CMESSAGE    ! Error return message
!
!  Common blocks
!
#include "cmaxsize.h"
#include "csubmodl.h"
#include "chsunits.h"
#include "chistory.h"
#include "ccontrol.h"
#include "ctime.h"
#include "cprintst.h"
!
#if defined(MPP)
#include "parvars.h"
#endif
!
!
!  Local variables
!
! External subroutines
!
      EXTERNAL GET_FILE
      INTEGER      I                                                    &
                                 ! Loop counter
     & ,MODJOBREL                ! modulus of jobrel_step
      CHARACTER*2  CDIGIT(10)    ! Character array of digits
      CHARACTER*80 CREQUEST      ! Template processing request
      CHARACTER*80 FILENAME
      DATA CDIGIT                                                       &
     & / "1 ","2 ","3 ","4 ","5 ","6 ","7 ","8 ","9 ","10" /
      DATA CREQUEST                                                     &
     & / "%%%  _JOB_MEMBER_   REL NET" /
!
!L----------------------------------------------------------------------
! Get name of pipe
      CALL GET_FILE(8,FILENAME,80,ICODE)
!L 1. Set character for internal model type
!L
      IF (I_AO == 1) THEN
        CREQUEST(5:5) = "A"
      ELSE IF (I_AO == 2) THEN
        CREQUEST(5:5) = "O"
      ELSE
        ICODE = I_AO
        CMESSAGE = 'JOBCTL1: invalid internal model indicator'
        GO TO 9999
      END IF

!L 2. Job Release - loop over release timesteps
!L
        DO I=1,10
          IF (JOBREL_STEPim(I,I_AO) <  0) THEN
            MODJOBREL=ABS(JOBREL_STEPim(I,I_AO))
            IF (MOD(STEPim(I_AO),MODJOBREL) == 0) THEN
              CREQUEST(18:19) = CDIGIT(I)
              IF (mype == 0) THEN
                IF (PrintStatus >= PrStatus_Oper) THEN
                  WRITE(6,*) 'JobCtl : Timestep ',STEPim(I_AO),         &
     &            ' Job No ',I,' released. CREQUEST : ',                &
     &              CREQUEST (1:Len_Trim(CREQUEST))
                END IF
              END IF
#if defined(MPP)
              IF(mype == 0) THEN
                WRITE(8,*) CREQUEST
#if defined(T3E)
! DEPENDS ON: um_fort_flush
                call um_fort_flush(8, icode)
#else
                CLOSE(8)
                OPEN(8,FILE=FILENAME)
#endif
              ENDIF
#else
              WRITE(8,*) CREQUEST
              CLOSE(8)
              OPEN(8,FILE=FILENAME)
#endif
            ENDIF
          ELSE
            IF (JOBREL_STEPim(I,I_AO) == STEPim(I_AO)) THEN
              CREQUEST(18:19) = CDIGIT(I)
              IF (mype == 0) THEN
                IF (PrintStatus >= PrStatus_Oper) THEN
                  WRITE(6,*) 'JobCtl : Timestep ',STEPim(I_AO),         &
     &            ' Job No ',I,' released. CREQUEST : ',                &
     &              CREQUEST (1:Len_Trim(CREQUEST))
                END IF
              END IF
#if defined(MPP)
              IF(mype == 0) THEN
                WRITE(8,*) CREQUEST
#if defined(T3E)
! DEPENDS ON: um_fort_flush
                call um_fort_flush(8, icode)
#else
                CLOSE(8)
                OPEN(8,FILE=FILENAME)
#endif
              ENDIF
#else
              WRITE(8,*) CREQUEST
              CLOSE(8)
              OPEN(8,FILE=FILENAME)
#endif
            ENDIF
          ENDIF
        ENDDO
!
 9999 CONTINUE
      RETURN
      END SUBROUTINE JOBCTL
#endif
