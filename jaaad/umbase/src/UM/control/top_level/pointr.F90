#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Change pointer system for "child" diags to final STASH list version
!
! Subroutine Interface:

      SUBROUTINE POINTR(NRECS)
      IMPLICIT NONE

! Description:
!   The stash list with preliminary pointer system, and the stash
!   index, are input to this routine. The stash list with final
!   pointer system is output.
!   Called by STPROC.
!
!   Fuller explanation:
!   Any diag in the stash list which has a processing code in the range
!   2-7 (i.e., accumulate, time mean, append time series, max, min,
!   trajectory) has one or more "child records". A child is another
!   diag, with the same m,s,i. The output from the parent diag is used
!   as input to the child diag, which is then processed to produce
!   further output. Each child record has an entry which points to its
!   parent record - the st_input_code entry. In routine PRELIM, a
!   preliminary pointer system is set up, involving the use
!   of the "extra entry", NELEMP+1. In each record, entry NELEMP+1 is
!   set to the current value of NRECS, i.e., the position of that
!   record in the prelim stash list. The value of the st_input_code
!   entry for a child record is set to the negative of the NRECS value
!   for its parent. Note that, in the prelim stash list, the children
!   of a particular parent appear immediately after the parent.
!   So, after PRELIM, each record in the stash list identifies
!   itself by its NELEMP+1 entry, and each child record identifies
!   its parent by its st_input_code entry. The final position of
!   each record in the stash list is given by the INDX_S array.
!   This subroutine therefore changes the  st_input_code entry of
!   each child record so that it agrees with INDX_S.
!   The NELEMP+1 entry is then no longer relevant.
!
! Method:
!   Uses INDX_S array to identify parent records (i.e., diagnostics
!   which have more than one entry in the stash list).
!
! Current code owner:  S.J.Swarbrick
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Apr. 95    Original code.  S.J.Swarbrick
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
#include "version.h"
#include "cstash.h"
#include "stextend.h"
#include "stparam.h"

! Subroutine arguments:

!   Scalar arguments with intent(in):

      INTEGER NRECS    ! No. of records in stash list

! Local scalars

      INTEGER MODL  ! Loop counter for internal models
      INTEGER ISEC  ! Do. sections
      INTEGER IITM  ! Do. items
      INTEGER ISTR  ! Position of parent record in stash list
      INTEGER IEND  ! Position of final child record in stash list
      INTEGER I1
      INTEGER I2
      INTEGER I3

!- End of Header ----------------------------------------------------

! Loop over models, section, items

      DO MODL=1,N_INTERNAL_MODEL_MAX
      DO ISEC=0,NSECTP
      DO IITM=1,NITEMP

! Examine INDX_S entry to find out whether there are child record(s)

        IF(INDX_S(2,MODL,ISEC,IITM) >= 2) THEN

          ISTR=     INDX_S(1,MODL,ISEC,IITM)
          IEND=ISTR+INDX_S(2,MODL,ISEC,IITM)-1

          DO I1=ISTR,IEND-1
            DO I2=I1+1,IEND
              IF(LIST_S(st_input_code,I2) ==                            &
     &          -LIST_S(NELEMP+1     ,I1))    THEN
                 LIST_S(st_input_code,I2)=-I1-NRECS
              END IF
            END DO
          END DO

          DO I3=ISTR,IEND
            IF(LIST_S(st_input_code,I3) <  0) THEN
               LIST_S(st_input_code,I3)=                                &
     &         LIST_S(st_input_code,I3)+NRECS
            END IF
          END DO

        END IF

      END DO  ! Items
      END DO  ! Sections
      END DO  ! Internal models

      RETURN
      END SUBROUTINE POINTR

!- End of subroutine code ---------------------------------------------
#endif
