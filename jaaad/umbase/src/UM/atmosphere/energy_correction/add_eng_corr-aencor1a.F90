#if defined(A14_1B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE ADD_ENG_CORR--------------------------------------
!LL
!LL  PURPOSE : PART OF ENERGY CORRECTION SUITE OF ROUTINES
!LL            - TO ADD IN TEMPERATURE CORRECTION TO
!LL              GLOBAL TEMPERATURE FIELD SO TO
!LL              CONSERVE TOTAL ENERGY GLOBALLY
!LL
!LL  NOT SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  MODEL            MODIFICATION HISTORY:
!LL VERSION  DATE
!LL   5.1   16/03/00  Altered for new dynamics
!LL   5.2   04/09/00  Add stash diagnostics. R.A.Stratton
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL
!LL  DOCUMENTATION :
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE ADD_ENG_CORR (ENERGY_CORR,T,row_length,rows,           &
     &                         model_levels,TSTEP,at_extremity,         &
#include "argsts.h"
     &                         STASHwork14)
!
      IMPLICIT NONE
!
#include "c_r_cp.h"
#include "csubmodl.h"
#include "typsts.h"

!----------------------------------------------------------------------
! VECTOR LENGTHS
!----------------------------------------------------------------------
!
      INTEGER                                                           &
     &     row_length                                                   &
                            ! IN row length
     &,    rows                                                         &
                            ! IN rows
     &,    model_levels     ! IN NUMBER OF LEVELS IN VERTICAL

      Logical                                                           &
     & at_extremity(4)      ! indicates whether PE at poles
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!----------------------------------------------------------------------
!
      REAL ENERGY_CORR         ! IN ENERGY CORRECTION
!
      REAL TSTEP               ! IN TIMESTEP
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE IN AND OUT
!----------------------------------------------------------------------
!
! INOUT sum of temperature increments
      REAL T(row_length,rows,model_levels)
!  diagnostics  out
      REAL                                                              &
     & STASHwork14(*)   ! STASH workspace
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!----------------------------------------------------------------------
      REAL                                                              &
     & work(row_length,rows,model_levels)

      INTEGER I,J,K                                                     &
                                   ! LOOP COUNTERs
     &,  icode                                                          &
                        ! return code
     &,im_index         ! model index

      character*80                                                      &
     &  cmessage        ! return message
      character(*)                                                      &
     &  RoutineName
      parameter (RoutineName='Add_Eng_corr')
!----------------------------------------------------------------------
! EXTERNAL SUBROUTINE CALLS
!----------------------------------------------------------------------
      External copydiag_3d,ereport
!
!*---------------------------------------------------------------------
!
      icode=0      ! initialise return code to zero

!----------------------------------------------------------------------
! CORRECT TEMPERATURE FOR ERROR IN ENERGY BUDGET OF THE
! PREVIOUS DAY
!----------------------------------------------------------------------
!
      DO K=1,model_levels
        DO J=1,rows
          DO I=1,row_length
            T(I,J,K) = T(I,J,K) + ENERGY_CORR*TSTEP
       END DO
      END DO
      END DO

#if !defined(SCMA)
!-----------------------------------------------------------------------
! Stash diagnostics
!-----------------------------------------------------------------------
! 14 181  T increment on model levels

      im_index = internal_model_index(atmos_im)

      If (sf(181,14)) then
        DO K=1,model_levels
          DO J=1,rows
            DO I=1,row_length
              work(I,J,K) = ENERGY_CORR*TSTEP
            END DO
          END DO
        END DO
! DEPENDS ON: copydiag_3d
        call copydiag_3d (stashwork14(si(181,14,im_index)),             &
     &       work,row_length,rows,model_levels,0,0,0,0, at_extremity,   &
     &       stlist(1,stindex(1,181,14,im_index)),len_stlist,           &
     &       stash_levels,num_stash_levels+1,atmos_im,14,181,           &
     &       icode,cmessage)
        if (icode >  0) then
          cmessage=":error in copydiag_3d(item 181)"//cmessage
        endif

      ENDIF

! Note old 14201 would be expensive to calculate here as it is
! now defined as cv*dt*(column integral of dry mass) and is therefore
! calculated and output from section 30.

      if (icode /= 0) then
! DEPENDS ON: ereport
         call ereport(RoutineName,icode,cmessage)
      endif
#endif
!----------------------------------------------------------------------
! ADD ENERGY CORRECTION INTO SUM OF DIABATIC FLUXES
!----------------------------------------------------------------------
! From UM 5.1 Moved elsewhere in code

      RETURN
      END SUBROUTINE ADD_ENG_CORR
#endif
