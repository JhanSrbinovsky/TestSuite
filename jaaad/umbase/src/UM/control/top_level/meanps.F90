#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL
!LL    Subroutine:
!LL    MEANPS
!LL
!LL    Purpose:
!LL    To mean partial sums. Sums obtained from the D1 array, put there
!LL    by ACUMPS.
!LL
!LL    Tested under compiler cft77
!LL    Tested under OS version:
!LL    UNICOS 5.1
!LL
!LL T.J., D.R.  <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.1  19/02/93  Use FIXHD(12) not FIXHD(1) as Version no in P21BITS
!LL   3.1   25/01/93 : Corrected LBPACK after change to definition.
!LL 3.4  16/6/94 : Change CHARACTER*(*) to CHARACTER*(80) N.Farnon
!     4.1  18/06/96  Changes to cope with changes in STASH addressing
!                    Author D.M. Goddard.
!LL   4.2  27/11/96  MPP changes for T3E.  Using READDUMP to
!LL                  read in partial mean files. K Rogers
!LL   4.3  23/01/97  Use MPP_LOOKUP adresses when on MPP
!LL        17/04/97  And pass D1_ADDR to UM_READDUMP
!LL                  S.D.Mullerworth
!LL   4.3  10/04/97  Add READHDR argument to READDUMP. K Rogers
!LL   5.0  3/6/99    Remove ICODE and CMESSAGE from UM_READDUMP
!LL                                                        P.Burton
!LL   5.1  24/01/00  Optimised to mean only tagged diagnostics, and
!LL                  to get data from D1 rather than from dump.
!LL                  S.D.Mullerworth
!LL   5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!LL
!LL    Programming standard:
!LL    UM Doc Paper 3
!LL
!LL    Logical system components covered:
!LL    C5
!LL
!LL    Project tasks:
!LL    C5,C51,C52
!LL
!LL    External documentation:
!LL    On-line UM document C5 - Control of means calculations
!LL
!*L    Interface and arguments:
      SUBROUTINE MEANPS(                                                &
     &  N_OBJS_D1,D1_ADDR                                               &
     &  ,LEN_DATA,D1,LD1,ID1,                                           &
#include "argsts.h"
     &  MEANING_PERIOD                                                  &
     &  )
!
      IMPLICIT NONE
#include "d1_addr.h"
!
      INTEGER                                                           &
     &  N_OBJS_D1

      INTEGER                                                           &
     &  D1_ADDR(D1_LIST_LEN,N_OBJS_D1)

      INTEGER                                                           &
     &  LEN_DATA,                                                       &
                                ! IN Length of model data
     &  MEANING_PERIOD          ! IN Meaning period (in multiples
                                !             of restart frequency)
!
      INTEGER                                                           &
     &  ID1(LEN_DATA)           ! IN Integer equivalence of data block
!
      REAL                                                              &
     &  D1(LEN_DATA)            ! IN/OUT Real equivalence of data block
                                !    containing meaned fields
!
      LOGICAL                                                           &
     &  LD1(LEN_DATA)           ! IN Logical equivalence of data block
!
!      Common blocks
!
#include "c_mdi.h"
#include "csubmodl.h"
#include "parparm.h"
#include "typsize.h"
#include "typsts.h"
#include "stparam.h"
!
!*L
!
!      Cray specific functions  UNIT,LENGTH

!
!      Local variables
!
      INTEGER                                                           &
     &  J,K                                                             &
                                ! Loop indices
     &  ,TAG                                                            &
                                ! Climate mean tag
     &  ,ADDRESS                                                        &
     &  ,PTD1
!
      REAL                                                              &
     &  FACTOR                                                          &
                                ! Meaning period (real)
     &  ,RFACTOR                 ! Reciprocal of FACTOR


!
!      Constants
!
      REAL                                                              &
     &       ONE                  ! 1.0


! Calculate divisor
      ONE=1.0
      FACTOR=MEANING_PERIOD
      RFACTOR=ONE/FACTOR

!----------------------------------------------------------------------
!     Loop through STASH list and process climate mean fields
!     NOTE: D1 contains partial sums put there by preceding ACUMPS call
!----------------------------------------------------------------------

      DO K=1,TOTITEMS
        TAG=STLIST(st_macrotag,K)/1000
        PTD1=STLIST(st_d1pos,K)
        IF(TAG /= 0.AND.STLIST(s_modl,k) == D1_ADDR(d1_imodl,PTD1))THEN
! Object is tagged for climate meaning and in relevant internal model
          ADDRESS=D1_ADDR(d1_address,PTD1) ! local address
! Divide whole field by FACTOR - except for RMDI
          DO J=ADDRESS,ADDRESS+D1_ADDR(d1_length,PTD1)-1
            IF (D1(J) /= RMDI)THEN
              D1(J)=D1(J)*RFACTOR
            ENDIF
          ENDDO
        ENDIF
      ENDDO

!**********************************************************************
!     End of loop over STASH list
!**********************************************************************

 999  CONTINUE
      RETURN
      END SUBROUTINE MEANPS
#endif
