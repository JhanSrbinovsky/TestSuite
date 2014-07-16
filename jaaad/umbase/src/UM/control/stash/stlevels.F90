#if defined(C84_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine STLEVELS -----------------------------------------------
!LL
!LL  Purpose: Generate a level index from STASHrecord and level_lists
!LL           and number of levels tailored to a particular diagnostic.
!LL           Also set levels and pseudo-levels information for encoding
!LL           PPheader details.  (This subroutine based on a merger
!LL           between GEN_INDEX and PP_COMPUTE_LEVEL).
!LL                  New subroutine STLEVELS is based on GEN_INDEX and
!LL                  PP_COMPUTE_LEVEL with merged functionality.
!LL           A general note as levels list is an integer
!LL           real values are multiplied by a 1000.0.
!LL           When computing the real value of the level for the
!LL           pp header it is necessary to divide by a 1000.0.
!LL           Levels that are affected by this are theta, pressure and
!LL           height. S. Anderson.
!LL
!LL  Author:   T.Johns
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.1  14/01/93  Set pseudo_level to 0 not IMDI if no pseudo-level.
!LL        29/01/93  Correct FORMAT statements.
!LL  3.1     14/01/93 Include PV levels as levels divided by 1000.
!LL   3.2  19/04/93  Correct roundoff error for LEVEL type conversion.
!LL                  1.0E-10 is added after REAL divide by 1000 (TCJ).
!LL  4.0  14/12/95  Correct long-standing error in input levels range
!LL                 to output levels list conversion.  RTHBarnes.
!    4.4  02/12/96 Time mean timeseries added R A Stratton.
!    5.0  10/06/99 Remove references to ak,bk,ak_lev,bk_lev, which are
!                  not used in downstream processing. Rick Rawlins
!    5.3  23/07/01 Replace lbvc magic no.s by parameters. R Rawlins
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered : C4?
!LL
!LL  Project task: C4
!LL
!LL  External documentation : UMDP no C4
!LL
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE STLEVELS(stash_control,stash_control_size,             &
     &     stash_levels,num_stash_levels,num_level_lists,               &
     &     stash_pseudo_levels,num_stash_pseudo,num_pseudo_lists,       &
     &     max_stash_levs,num_levs_in,num_levs_out,index_size,          &
     &     index_lev,level_list,                                        &
     &     lbvcl,level,pseudo_level,                                    &
     &     icode,cmessage)
!
      IMPLICIT NONE
!
      INTEGER                                                           &
     &       stash_control_size,                                        &
                                 ! IN size of stash control record
     &       stash_control(stash_control_size),                         &
                                               ! IN  stash control
     &       num_stash_levels,                                          &
                                 ! IN max. no of hts for a levels list
     &       num_level_lists,                                           &
                                 ! IN max. no of level lists
     &       stash_levels(num_stash_levels+1,num_level_lists),          &
                                                               ! IN
!                                !    lookup table for level lists
     &       num_stash_pseudo,num_pseudo_lists,                         &
                                               ! IN dims of pseudo_levs
     &       stash_pseudo_levels(num_stash_pseudo+1,num_pseudo_lists),  &
!                                ! IN lookup table for pseudo-lev lists
     &       max_stash_levs,                                            &
                                 ! IN max. no of output levels
     &       num_levs_in,                                               &
                                 ! OUT no of levels in input data
     &       num_levs_out,                                              &
                                 ! OUT no of levels in output data
     &       index_size,                                                &
                                 ! OUT no of levels in levels index
     &       index_lev(max_stash_levs),                                 &
                                        ! OUT index of output level
!                                               relative to input level
     &       level_list(max_stash_levs),                                &
                                         ! OUT value of model level
     &       pseudo_level(max_stash_levs),                              &
                                           ! OUT Value of pseudo levels
     &       lbvcl,                                                     &
                                 ! IN  vertical coordinate PP code
     &       icode               ! OUT error code
      REAL                                                              &
     &       level(max_stash_levs)  ! OUT Value of output levels (real)
      CHARACTER*(*)                                                     &
     &       cmessage            ! OUT error message
!*----------------------------------------------------------------------
! Parameters
!
#include "sterr.h"
#include "stparam.h"
#include "cppxref.h"
!
! Local variables
!
      INTEGER                                                           &
     &       index_pseudo_lev(max_stash_levs),                          &
                                               ! Pseudo-level 1D index
     &       num_pseudo_in,num_pseudo_out,                              &
                                               ! Number of pseudo levs
     &       k2,ml,kl,                                                  &
                                       ! loop counts
     &       NI,NO,                                                     &
                                       ! Number In/Out
     &       indx1,                                                     &
                                       ! index count
     &       ilev,                                                      &
                                       ! Integer level/pseudo-level
     &       what_mean,what_proc       ! Meaning and processing code
!
! First compute the index for physical levels
!
      IF(STASH_CONTROL(st_input_bottom) <  0) THEN ! Input LEVELS list
        NI=-STASH_CONTROL(st_input_bottom)
        NUM_LEVS_IN=STASH_LEVELS(1,NI)
        IF(STASH_CONTROL(st_output_bottom) <  0) THEN ! LEVELS LIST out
          NO=-STASH_CONTROL(st_output_bottom)
          NUM_LEVS_OUT=STASH_LEVELS(1,NO)
          INDX1=0
          DO ML=1,NUM_LEVS_OUT
            ilev=STASH_LEVELS(ML+1,NO)    !  Level required
            DO KL=1,NUM_LEVS_IN
              IF(STASH_LEVELS(KL+1,NI) == ilev) THEN
                INDX1=INDX1+1
                INDEX_LEV(INDX1)=KL   ! Relative position of Input to Ou
                level_list(indx1)=ilev
                GOTO 400
              ENDIF
            ENDDO
            ICODE=nonsense
            WRITE(CMESSAGE,101) 'Output level ',ilev,                   &
     &                          ' not found in input levels list'
            GOTO 999
 400        CONTINUE
            ENDDO
        ELSE           !  Output as a Level range
          NUM_LEVS_OUT=STASH_CONTROL(st_output_top)-                    &
     &                 STASH_CONTROL(st_output_bottom)+1
          ilev=STASH_CONTROL(st_output_bottom) !1st output model level
          DO KL=1,NUM_LEVS_IN
            IF(STASH_LEVELS(KL+1,NI) == ilev) THEN
              INDEX_LEV(1)=KL ! Relative posn of Input to the 1st level
              level_list(1)=ilev
              GOTO 401
            ENDIF
          ENDDO
          ICODE=nonsense
          WRITE(CMESSAGE,101) 'Output bottom model level ',ilev,        &
     &                        ' not found in input levels list'
          GOTO 999
 401      CONTINUE
          DO KL=2,NUM_LEVS_OUT
            INDEX_LEV(KL)=INDEX_LEV(KL-1)+1
            level_list(kl)=level_list(kl-1)+1
          ENDDO
        ENDIF
      ELSEIF(STASH_CONTROL(st_input_bottom) == 100) THEN !Special level
          NUM_LEVS_IN=1
          NUM_LEVS_OUT=1
          INDEX_LEV(1)=1
          level_list(1)=1 ! could be worth setting to some nonsense no.
      ELSE     !  Input is Model level range
        NUM_LEVS_IN=STASH_CONTROL(st_input_top)-                        &
     &              STASH_CONTROL(st_input_bottom)+1
        IF(STASH_CONTROL(st_output_bottom) <  0) THEN ! LEVELS LIST out
          NO=-STASH_CONTROL(st_output_bottom)
          NUM_LEVS_OUT=STASH_LEVELS(1,NO)
          INDX1=0
          DO ML=1,NUM_LEVS_OUT
            ilev=STASH_LEVELS(ML+1,NO)    ! Output level reqd
            DO KL=1,NUM_LEVS_IN
              IF((STASH_CONTROL(st_input_bottom)+KL-1) == ilev) THEN
                INDX1=INDX1+1
                INDEX_LEV(INDX1)=KL   ! Relative posn of output to inpt
                level_list(INDX1)=ilev
                GOTO 402
              ENDIF
            ENDDO
            ICODE=nonsense
            WRITE(CMESSAGE,101) 'Output model level ',ilev,             &
     &                          ' not in input model level range'
            GOTO 999
 402        CONTINUE
          ENDDO
        ELSE     !   Output as model level range
! Do some consistency checks here to ensure valid processing request
! output bottom should be greater or equal to input bottom
          IF (stash_control(st_output_bottom) <                         &
     &       stash_control(st_input_bottom)) THEN
            icode=nonsense
            write(cmessage,103)'bad level spec, bot input>output',      &
     &       stash_control(st_input_bottom),                            &
     &       stash_control(st_output_bottom)
            goto 999 ! jump to error
          ELSEIF (stash_control(st_output_top) >                        &
     &         stash_control(st_input_top)) THEN
            icode=nonsense
            write(cmessage,103)'bad level spec, top input<output',      &
     &        stash_control(st_input_top),                              &
     &        stash_control(st_output_top)
              goto 999 ! jump to error
          ENDIF
          NUM_LEVS_OUT=STASH_CONTROL(st_output_top)-                    &
     &                 STASH_CONTROL(st_output_bottom)+1
          INDEX_LEV(1)=STASH_CONTROL(st_output_bottom)-                 &
     &                 STASH_CONTROL(st_input_bottom)+1
          level_list(1)=stash_control(st_output_bottom)
          DO kl=2,NUM_LEVS_OUT
            INDEX_LEV(kl)=INDEX_LEV(kl-1)+1
            level_list(kl)=level_list(kl-1)+1
          ENDDO
        ENDIF
      ENDIF
      index_size=num_levs_out
      IF (num_levs_out >  num_levs_in) THEN   ! things very badly wrong
        icode=nonsense
        write(cmessage,103)'asking for num_levs_out>num_levs_in',       &
     &   num_levs_out,num_levs_in
        goto 999 ! jump to return
      ENDIF
!
! Next, compute actual (physical) levels for encoding PPheaders
!
      IF (STASH_CONTROL(st_output_bottom) <  0) THEN ! Levels List ?
        NO=-STASH_CONTROL(st_output_bottom)     ! Index of Levels list

          ! Remove scaling (by factor 1000) of vertical level coord
          ! for certain types of STASH output [originally needed to
          ! store in an intermediary integer array]
          IF( LBVCL  ==  ppx_lbvc_height   .OR.                         &
                                                !  height levels
     &        LBVCL  ==  ppx_lbvc_pressure .OR.                         &
                                                ! pressure levels
     &        LBVCL  ==  ppx_lbvc_theta    .OR.                         &
                                                ! theta levels
     &        LBVCL  ==  ppx_lbvc_PV ) THEN     ! potential vorticity


          DO ML=1,NUM_LEVS_OUT
            LEVEL(ML)=REAL(STASH_LEVELS(ML+1,NO))*0.001+1.0E-10
          ENDDO
        ELSE
          DO ML=1,NUM_LEVS_OUT
            LEVEL(ML)=REAL(STASH_LEVELS(ML+1,NO))
          ENDDO
        ENDIF
      ELSEIF (STASH_CONTROL(st_output_bottom) == st_special_code) THEN
       ! Special level.
       ! The LEVEL array is not used by the model except to construct pp
       ! header items at output. The value of -1.0 is set as a flag for
       ! special levels so that routine PP_HEAD will insert the lbvc
       ! item in STASHmaster record.
        DO ML=1,NUM_LEVS_OUT
          LEVEL(ML)=-1.0
        ENDDO
      ELSE
        DO ML=1,NUM_LEVS_OUT
          LEVEL(ML)=REAL(STASH_CONTROL(st_output_bottom)+ML-1)
        ENDDO
      ENDIF
!
!
! Now reset the number of output levels to 1 if vertical compression is
! to be done in SPATIAL.  NB: index_lev and level_list need to be filled
! with values corresponding to the full range of levels processed.
!
      what_proc=STASH_CONTROL(st_proc_no_code)
      what_mean=(STASH_CONTROL(st_gridpoint_code)/block_size)*block_size
      IF(what_mean == vert_mean_base .OR. what_mean == global_mean_base &
     &   .OR. what_proc == st_time_series_code                          &
     &   .OR. what_proc == st_time_series_mean                          &
     &   .OR. what_proc == st_append_traj_code) num_levs_out=1
!
! Next compute the index for pseudo levels, if there are any
!
      IF(STASH_CONTROL(st_pseudo_in) >  0) THEN ! Input PSEUDO_LEVELS
        NI=STASH_CONTROL(st_pseudo_in)
        num_pseudo_in=STASH_PSEUDO_LEVELS(1,NI)
        IF(STASH_CONTROL(st_pseudo_out) >  0) THEN ! Output PSEUDO_LEVS
          NO=STASH_CONTROL(st_pseudo_out)
          num_pseudo_out=STASH_PSEUDO_LEVELS(1,NO)
          INDX1=0
          DO ML=1,NUM_PSEUDO_OUT
            ilev=STASH_PSEUDO_LEVELS(ML+1,NO)   !  Level required
            DO KL=1,NUM_PSEUDO_IN
              IF(STASH_PSEUDO_LEVELS(KL+1,NI) == ilev) THEN
                INDX1=INDX1+1
                INDEX_PSEUDO_LEV(INDX1)=KL
                pseudo_level(indx1)=ilev
                GOTO 500
              ENDIF
            ENDDO
            ICODE=nonsense
            WRITE(CMESSAGE,101) 'Output pseudo level ',ilev,            &
     &                          ' not found in input levels list'
            GOTO 999
 500        CONTINUE
          ENDDO
        ELSE  ! Illegal combination
          ICODE=nonsense
          WRITE(CMESSAGE,101) 'Input pseudo level list ',NI,            &
     &         ' has illegal output pseudo levels list'
          GOTO 999
        ENDIF
      ELSE  ! Only levels lists are supported for pseudo levels
        num_pseudo_out=0
      ENDIF
!
! Next expand the separate indexes and physical levels arrays into
! combined arrays if necessary, taking care not to overwrite earlier
! parts of the arrays.  If no pseudo-levels, set pseudo-level to 0.
!
      IF (num_pseudo_out >  0) THEN
        DO K2=num_pseudo_out,1,-1
          DO ML=1,num_levs_out
            INDEX_LEV(ML+(K2-1)*num_levs_out)=                          &
     &        (INDEX_PSEUDO_LEV(K2)-1)*num_levs_in+INDEX_LEV(ML)
            level(ML+(K2-1)*num_levs_out)=level(ML)
          ENDDO
          DO ML=num_levs_out,1,-1
            pseudo_level(ML+(K2-1)*num_levs_out)=pseudo_level(K2)
          ENDDO
        ENDDO
        num_levs_out=num_levs_out*num_pseudo_out
      ELSE
        DO ML=1,num_levs_out
          pseudo_level(ML)=0
        ENDDO
      ENDIF
!
999   CONTINUE ! jump here for error return
 101  FORMAT('STLEVELS : ',a,i6,a)
 103  FORMAT('STLEVELS : >> FATAL ERROR <<',a,2i5)
      RETURN
      END SUBROUTINE STLEVELS

#endif
