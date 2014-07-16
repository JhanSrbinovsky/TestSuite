#if defined(CONTROL) || defined(SETUP) || defined(COMB)                \
 || defined(PICK) || defined(HPRT)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: INITCHST
!LL
!LL  Purpose: To set integer areas of history file namelist to
!LL           zero and set character areas to blank
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.0
!LL
!LL  Author:   A.Sangster
!LL
!LL  Code version no: 1           Date: 20 January 1990
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL  3.5  06/04/95  Sub-Models stage 1: revise History and Control file
!LL                 contents.  RTHBarnes.
!LL  4.3  17/02/97  Further initialisations needed for namelist reads
!LL                 to have initialised fields             L C Wiles
!LL
!LL  Programming standard: UM Doc Paper 3, draft version 3 (15/1/90)
!LL
!LL  Logical components covered:
!LL
!LL  Project task: H
!LL
!LL  Documentation:  Unified Model Documentation Paper
!LL                  H- History Bricks
!LL
!
!*L  Interface and arguments:
!
      SUBROUTINE INITCHST
!
      IMPLICIT NONE
!*
!
!L Common blocks
!
#include "csubmodl.h"
#include "chsunits.h"
#include "chistory.h"
!
!*L EXTERNAL subroutines called
!     None
!*
!
! Local variables
!
      INTEGER I  ! Loop counter
      INTEGER J  ! Loop counter
!L
!L 1. Set common block area to zero or blank
!L
! for NLIHISTO - integer overall model variables
      do  i = 1,6
      MODEL_DATA_TIME(i)=-32768
      end do
      RUN_MEANCTL_RESTART=0
      RUN_INDIC_OP=0
      do  i = 1,6
      RUN_RESUBMIT_TARGET(i)=-32768
      end do
      do  i = 20,nunits
      FT_LASTFIELD(i)=-32768
      end do
! for NLCHISTO - character overall model variables
      RUN_HIST_TYPE='          '
      RUN_COMPCODE='              '
      RUN_LAST_MEAN='              '
      RUN_MEANS_TO_DO=' '
      RUN_OCEAN_FIRST=' '
      RUN_TYPE='Reconfig'
      RUN_JOB_NAME='        '
      RUN_ID='     '
      RUN_RESUBMIT=' '
      RUN_RESUBMIT_Q='            '
      RUN_RESUBMIT_TIME='                    '
      RUN_RESUBMIT_CPU='      '
      RUN_RESUBMIT_MEMORY='      '
      RUN_RESUBMIT_PRTY='  '
      RUN_RESUBMIT_JOBNAME='        '
      do  i = 20,nunits
      FT_ACTIVE(i)=' '
      end do
! for NLIHISTG - integer generic model variables
      do  i = 1,n_internal_model_max
      H_STEPim(i)=0
      H_GROUPim(i)=0
      MEAN_OFFSETim(i)=-32768
      OFFSET_DUMPSim(i)=-32768
      end do
!!    MEAN_NUMBERim(1)=4
!!    do  i = 2,n_internal_model_max
!!    MEAN_NUMBERim(i)=0
!!    end do
      do  i = 1,n_internal_model_max
      do  j = 1,4
      RUN_MEANCTL_INDICim(j,i)=1
      end do
      end do
! for NLCHISTG - character generic model variables
      do  i = 1,n_internal_model_max
      END_DUMPim(i)='              '
      RESTARTim(i)='                                                    &
     &                          '
      SAFEDMPim(i)='              '
      NEWSAFEim(i)='              '
      LASTATMim(i)='              '
      CURRATMim(i)='              '
      LASTDMPim(i)='              '
      end do
!
      do  i = 1,nunits
         MODEL_FT_UNIT(I)=' '
      end do
!
!L
!L 2. Return
!L
      RETURN
      END SUBROUTINE INITCHST
#endif
