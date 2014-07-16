#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+
! Subroutine Interface:

      SUBROUTINE INPUTL(NRECS,                                          &
#include "argppx.h"
     &                   NLEVELS,ErrorStatus,CMESSAGE)
      IMPLICIT NONE
! Description:
!
! Method:
!
! Current code owner:  S.J.Swarbrick
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Apr. 95    Original code.  S.J.Swarbrick
!   4.1     Apr. 96      Allow for prognostics with pseudo-levels:
!                        LevFlag=0 now implies input on all available
!                        levels and pseudo-levels. Also other additions
!                        for wave-model grid.
!                                                          S.J.Swarbrick
!   4.2     06/09/96   MPP code : Use local size for calculating
!                      horizontal dimension.    P.Burton
!
!   4.3     30/01/97   Ensure that domain decomposition is consistent
!                      with submodel. R.Rawlins
!LL 5.0     23/06/99   Added halo_type argument to GLOBAL_TO_LOCAL
!LL                    ADDRLN routines.                   P.Burton
!LL  5.1  22/02/00  Add PARVARS for TYPSIZE                 P.Burton
!   5.1     17/04/00   Removed ErrorStatus argument from ADDRLN
!                                                      P.Burton
!   5.1     15/05/00   S-N ordering consistency correction. R Rawlins
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!LL  5.5  08/08/00  Modification for parallelisation of WAM.
!                   Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!LL  6.2  23/11/05  Removed all references to the wavemodel.
!LL                 T.Edwards
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
#include "ppxlook.h"
#include "parvars.h"
#include "typsize.h"
#include "cstash.h"
#include "stextend.h"
#include "model.h"
#include "stparam.h"
#if defined(MPP)
#include "decomptp.h"
#endif

! Subroutine arguments:
!   Scalar arguments with intent(in):
      INTEGER NRECS

!   Scalar arguments with intent(out):
      INTEGER NLEVELS   ! total no. of sets of stash levels

!   Scalar arguments with intent(out):

! ErrorStatus:
      INTEGER ErrorStatus

! Local scalars:
      LOGICAL MODEL_LEV
      CHARACTER*80 CMESSAGE
      LOGICAL LADD
      LOGICAL LDUPLL
      INTEGER I,IL,ILIN
      INTEGER ISTART,IEND
      INTEGER MODL
      INTEGER ISEC
      INTEGER IITM
      INTEGER IP_IN
      INTEGER IX1,IX2
      INTEGER IY1,IY2
      INTEGER IZ_IN
      INTEGER LEN_IN
      INTEGER LEN_PRIMIN
      INTEGER NDUPLL
      INTEGER NLEVIN
      INTEGER LENO
      INTEGER IPF,IPL
#if defined(MPP)
! local versions of the global subdomain limits
      INTEGER local_IX1,local_IX2,local_IY1,local_IY2
      INTEGER                                                           &
     &        orig_decomp                                               &
                               ! MPP decomposition before start
     &       ,decomp_type                                               &
                               ! decomposition type
     &       ,sm_ident         ! submodel identifier
#endif

! Function and subroutine calls:
      LOGICAL  DISCT_LEV
      INTEGER  EXPPXI
      EXTERNAL EXPPXI,LEVSRT,LLTORC,ADDRLN,OCNVOL
#if defined(MPP)
      External CHANGE_DECOMPOSITION,GLOBAL_TO_LOCAL_SUBDOMAIN
#endif

!- End of Header ----------------------------------------------------

#if defined(MPP)
      orig_decomp = current_decomp_type
#endif

      DO MODL=1,N_INTERNAL_MODEL_MAX
#if defined(MPP)
!
!    Ensure that domain decomposition is consistent with submodel
!
      sm_ident = SUBMODEL_PARTITION_INDEX(MODL)
      IF(sm_ident == atmos_sm) THEN
         decomp_type = decomp_standard_atmos
      ELSEIF(sm_ident == ocean_sm) THEN
         decomp_type = decomp_standard_ocean
      ELSE                            ! No decomposition defined:
         decomp_type = orig_decomp    !  return to original
      ENDIF

! DEPENDS ON: change_decomposition
      CALL CHANGE_DECOMPOSITION(decomp_type,ErrorStatus)

      IF(ErrorStatus >  0) THEN
         CMESSAGE='INPUTL: ERROR in changing MPP decomposition'
         write(6,*) CMESSAGE
         GOTO 999
      ENDIF
#endif
      DO ISEC=0,PPXREF_SECTIONS
      DO IITM=1,PPXREF_ITEMS
        IF(INDX_S(2,MODL,ISEC,IITM) >= 1) THEN
! At least one stash rec
! DEPENDS ON: exppxi
        IGP     = EXPPXI(MODL,ISEC,IITM,ppx_grid_type     ,             &
#include "argppx.h"
     &                                         ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        ILEV    = EXPPXI(MODL,ISEC,IITM,ppx_lv_code       ,             &
#include "argppx.h"
     &                                         ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IFLAG   = EXPPXI(MODL,ISEC,IITM,ppx_lev_flag      ,             &
#include "argppx.h"
     &                                         ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IPSEUDO = EXPPXI(MODL,ISEC,IITM,ppx_pt_code       ,             &
#include "argppx.h"
     &                                         ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IPFIRST = EXPPXI(MODL,ISEC,IITM,ppx_pf_code       ,             &
#include "argppx.h"
     &                                         ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IPLAST  = EXPPXI(MODL,ISEC,IITM,ppx_pl_code       ,             &
#include "argppx.h"
     &                                         ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        HALO_TYPE  = EXPPXI(MODL,ISEC,IITM,ppx_halo_type,               &
#include "argppx.h"
     &                                         ErrorStatus,CMESSAGE)

          ISTART=       INDX_S(1,MODL,ISEC,IITM)   ! Pos of 1st rec
          IEND  =ISTART+INDX_S(2,MODL,ISEC,IITM)-1 ! Pos of last rec
! Diagnostics with input on levels list (IFLAG=1),
!  rather than on all possible levels
          IF((IFLAG  == 1   ).AND.                                      &
                                                    !Input on lev list
     &       (ISTART == IEND).AND.                                      &
                                                    !Only 1 stash rec
     &       (LIST_S(st_output_bottom,ISTART) <  0))                    &
                                                    !Output on lev list
     &        THEN
! Only one stash record for this m,s,i - output levels list is
!  the same as the input levels list
            LIST_S(st_input_bottom ,ISTART)=                            &
     &      LIST_S(st_output_bottom,ISTART)
          ELSE IF (IFLAG == 1.AND.ILEV /= 5) THEN
! Input on levels list & more than one stash request -
!  construct input levels list
            NLEVELS=NLEVELS+1
            IF (NLEVELS >  NLEVLSTSP) THEN
              WRITE(6,*) 'ERROR IN ROUTINE INPUTL:'
              WRITE(6,*) 'TOO MANY STASH LEVELS LISTS REQUESTED ',      &
     &                   'ARRAYS WILL BE OVERWRITTEN'
              WRITE(6,*) 'REDUCE NUMBER OF LEVELS LISTS'
              ErrorStatus=1
              GO TO 999
            END IF
! Construct input levels list: this is the combined list of all
!  the output levels for all the stash requests for this m,s,i
            NLEVIN=1
! Set levels list type - real or integer
! DEPENDS ON: disct_lev
            MODEL_LEV=DISCT_LEV(ILEV,ErrorStatus,CMESSAGE)
            IF (.NOT.MODEL_LEV) THEN
! Non-model levels - real
              LLISTTY(NLEVELS)='R'
! Model levels - integer
            ELSE
              LLISTTY(NLEVELS)='I'
            END IF
! Loop over stash recs for this m,s,i
            DO I=ISTART,IEND
! Pointer for input level list
              LIST_S(st_input_bottom,I)=-NLEVELS
              IF (LIST_S(st_output_bottom,I) <  0) THEN
! There is an output levels list:
!  For each of the levels in the output levels lists for the stash
!   record I, find out whether this level is already present in the
!   input levels list NLEVELS constructed so far.
!   If it is, set LADD=F. Otherwise, LADD=T.
! Loop over output levels and check each one
                DO IL=2,LEVLST_S(1,-LIST_S(st_output_bottom,I))+1
                  LADD=.TRUE.
                  IF(NLEVIN >  1) THEN
                    DO ILIN=2,NLEVIN
                      IF(LIST_S(st_output_top,I) /= 1) THEN
! Non-model levels: real
                        IF(RLEVLST_S(IL,-LIST_S(st_output_bottom,I))    &
     &                   ==                                             &
     &                     RLEVLST_S(ILIN,NLEVELS)) LADD=.FALSE.
                      ELSE
! Model levels: integer
                        IF( LEVLST_S(IL,-LIST_S(st_output_bottom,I))    &
     &                   ==                                             &
     &                      LEVLST_S(ILIN,NLEVELS)) LADD=.FALSE.
                      END IF
                    END DO
                  END IF

! If LADD=T, add level 'IL' from stash record 'I' output levels list
!  to input levels list NLEVELS
                  IF (LADD) THEN
                    NLEVIN=NLEVIN+1
                    IF(LIST_S(st_output_top,I) /= 1) THEN
                      RLEVLST_S(NLEVIN,NLEVELS)=                        &
     &                RLEVLST_S(IL,-LIST_S(st_output_bottom,I))
                    ELSE
                      LEVLST_S(NLEVIN,NLEVELS)=                         &
     &                LEVLST_S(IL,-LIST_S(st_output_bottom,I))
                    END IF
                  END IF
                END DO     ! Loop over levels

              ELSE
! Contiguous range of model levels for output, rather than list
!  Compare output levels range for stash record I with input levs
!   range NLEVELS. Any of the output levels not already present
!   in the input range is added to the input list.
                DO IL=LIST_S(st_output_bottom,I),                       &
     &                LIST_S(st_output_top   ,I)
                  LADD=.TRUE.
                  DO ILIN=2,NLEVIN
                    IF(IL == LEVLST_S(ILIN,NLEVELS)) LADD=.FALSE.
                  END DO
                  IF(LADD) THEN
                    NLEVIN=NLEVIN+1
                    LEVLST_S(NLEVIN,NLEVELS)=IL
                  END IF
                END DO
              END IF   !  Levels list/range
            END DO     !  Loop over stash recs

! Record no. of levels in input list just constructed
            LEVLST_S(1,NLEVELS)=NLEVIN-1

            IF(NLEVIN-1 == 0) THEN
              WRITE(6,*) 'ORDINARY LEVEL'
              WRITE(6,*) 'ISEC=',ISEC
              WRITE(6,*) 'IITM=',IITM
              WRITE(6,*) 'NLEVELS=',NLEVELS
              DO I=ISTART,IEND
                WRITE(6,*) 'I=',I
                WRITE(6,*) 'LIST_S(st_output_bottom)=',                 &
     &                      LIST_S(st_output_bottom,I)
                IF (LIST_S(st_output_bottom,I) <  0) THEN
                  DO IL=2,LEVLST_S(1,-LIST_S(st_output_bottom,I))+1
                    WRITE(6,*) 'IL=',IL
                    WRITE(6,*)                                          &
     &              'LEVLST',LEVLST_S(IL,-LIST_S(st_output_bottom,I))
                  END DO
                ELSE
                  WRITE(6,*)                                            &
     &            'LIST_S(st_output_top=',LIST_S(st_output_top,I)
                END IF
              END DO
            END IF
! Sort levels list
! DEPENDS ON: levsrt
            CALL LEVSRT(LLISTTY(  NLEVELS), LEVLST_S(1,NLEVELS),        &
     &                 LEVLST_S(2,NLEVELS),RLEVLST_S(2,NLEVELS))

! Determine whether this levels list is a duplicate of another list
! DEPENDS ON: duplevl
            CALL DUPLEVL(NLEVELS,LDUPLL,NDUPLL)
            IF (LDUPLL) THEN
! Duplicate list at NDUPLL - reset pointer and reduce NLEVELS by 1
              NLEVELS=NLEVELS-1
              DO I=ISTART,IEND
                LIST_S(st_input_bottom,I)=-NDUPLL
              END DO
            END IF
          END IF   !Levels lists

! Pseudo levels lists
          IF((IFLAG  == 1   ).AND.                                      &
     &      ((ISTART == IEND).OR.(IPSEUDO == 0)) ) THEN
! Either no pseudo levels or only one request:
! Input pseudo levels list equals output list
            LIST_S(st_pseudo_in,ISTART)=LIST_S(st_pseudo_out,ISTART)
          ELSE IF (IFLAG == 1) THEN
! Input pseudo levels list with more than one request
            NPSLISTS=NPSLISTS+1
            IF(NPSLISTS >  NPSLISTP) THEN
              WRITE(6,*) 'ERROR IN ROUTINE INPUTL:'
              WRITE(6,*)                                                &
     &       'TOO MANY STASH PSEUDO LEVELS LISTS REQUESTED ',           &
     &       'ARRAYS WILL BE OVERWRITTEN'
              WRITE(6,*) 'REDUCE NUMBER OF PSEUDO LISTS'
              ErrorStatus=1
              GO TO 999
            END IF
! Construct input pseudo list: combined list of all output
!  pseudo levels for all stash requests for this m,s,i
            NLEVIN=0
            DO I=ISTART,IEND
              LIST_S(st_pseudo_in,I)=NPSLISTS
              DO IL=1,LENPLST(LIST_S(st_pseudo_out,I))
                LADD=.TRUE.
                IF(NLEVIN >  0) THEN
                  DO ILIN=1,NLEVIN
                    IF( PSLIST_D(IL,LIST_S(st_pseudo_out,I)) ==         &
     &                  PSLIST_D(ILIN,NPSLISTS)) LADD=.FALSE.
                  END DO
                END IF
                IF(LADD) THEN
                  NLEVIN=NLEVIN+1
                  PSLIST_D(NLEVIN,NPSLISTS)=                            &
     &            PSLIST_D(IL,LIST_S(st_pseudo_out,I))
                END IF
              END DO
            END DO
            LENPLST(NPSLISTS)=NLEVIN

            IF(NLEVIN == 0) THEN
              WRITE(6,*) 'PSEUDO LEVEL'
              WRITE(6,*) 'ISEC=',ISEC
              WRITE(6,*) 'IITM=',IITM
              WRITE(6,*) 'NPSLISTS=',NPSLISTS
              DO I=ISTART,IEND
                WRITE(6,*) 'I=',I
                WRITE(6,*) 'LENPLST=',LENPLST(LIST_S(st_pseudo_out,I))
                DO IL=1,LENPLST(LIST_S(st_pseudo_out,I))
                  WRITE(6,*) 'IL=',IL
                  WRITE(6,*)                                            &
     &            'PSLIST_D',PSLIST_D(IL,LIST_S(st_pseudo_out,I))
                END DO
              END DO
            END IF
! Sort input pseudo levels list
! DEPENDS ON: levsrt
            CALL LEVSRT('I',LENPLST(NPSLISTS),PSLIST_D(1,NPSLISTS),     &
     &                  REAL(PSLIST_D(1,NPSLISTS)))
! Find out if duplicate
! DEPENDS ON: duppsll
            CALL DUPPSLL(LDUPLL,NDUPLL)
            IF(LDUPLL) THEN
! Duplicate pseudo list at NDUPLL
              NPSLISTS=NPSLISTS-1
              DO I=ISTART,IEND
                LIST_S(st_pseudo_in,I)=NDUPLL
              END DO
            END IF
          ELSE IF (IFLAG == 0.AND.IPSEUDO /= 0) THEN
! Input pseudo levels list contains all possible pseudo levels for
!  this diagnostic
            NPSLISTS=NPSLISTS+1
            DO I=ISTART,IEND
              LIST_S(st_pseudo_in,I)=NPSLISTS
! Decode first & last pseudo level codes from stash master
! DEPENDS ON: pslevcod
              CALL PSLEVCOD(IPFIRST,IPF,'F',ErrorStatus,CMESSAGE)
! DEPENDS ON: pslevcod
              CALL PSLEVCOD(IPLAST ,IPL,'L',ErrorStatus,CMESSAGE)
! Construct list
              DO NLEVIN = IPF,IPL
                PSLIST_D(NLEVIN,NPSLISTS)=NLEVIN
              END DO
            END DO
            LENPLST(NPSLISTS)=IPL-IPF+1
          END IF   ! Pseudo levels

! Calculate horizontal factor for input length
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,90,-90,0,360,IY1,IY2,IX1,IX2)

#if defined(MPP)
! Convert from global to local subdomain limits
! DEPENDS ON: global_to_local_subdomain
        CALL GLOBAL_TO_LOCAL_SUBDOMAIN( .TRUE., .TRUE.,                 &
     &                                  IGP,halo_type,mype,             &
     &                                  IY1,IX2,IY2,IX1,                &
     &                                  local_IY1,local_IX2,            &
     &                                  local_IY2,local_IX1)
        IX1=local_IX1
        IX2=local_IX2
        IY1=local_IY1
        IY2=local_IY2
#endif
! All sub-model grids: atmos/ocean/wave now ordered S->N
           LEN_IN=(IX2-IX1+1)*(IY2-IY1+1)

! Calculate vertical levels factor for input length
          IF(ILEV /= 5) THEN
! More than one level
            IF(LIST_S(st_input_bottom,ISTART) <  0) THEN
! Level list
              IZ_IN=LEVLST_S(1,-LIST_S(st_input_bottom,ISTART))
            ELSE
! Range of model levs
              IZ_IN=LIST_S(st_input_top   ,ISTART)-                     &
     &              LIST_S(st_input_bottom,ISTART)+1
            END IF
          ELSE
! Single level input
            IZ_IN=1
          END IF

! Calculate pseudo levels factor for input length
          IF(IPSEUDO /= 0) THEN
            IP_IN=LENPLST(LIST_S(st_pseudo_in,ISTART))
          ELSE
            IP_IN=1
          END IF

! Calculate input length for this diag. and store in LIST_S
! Input_code <  0 means that a diag already processed into D1 is being
!   reprocessed, so input length of child diag equals output length of
!   parent.
! Otherwise, the input len is given by the product of the appropriate
!   x-,y-,z-, and p-dimensions.
          DO I=ISTART,IEND
            IF(LIST_S(st_input_code  ,I) >= 0) THEN
               LIST_S(st_input_length,I)=LEN_IN*IZ_IN*IP_IN
            ELSE
               LIST_S(st_input_length ,I)=                              &
     &         LIST_S(st_output_length,-LIST_S(st_input_code,I))
            END IF
 ! Store model no. in last element of LIST_S - for ADDRES
            LIST_S(NELEMP+1,I)=MODL
          END DO

! Recalculate input length for non-primary (length unchanged for
! most cases) and store in IN_S array.
          IF (ISEC /= 0) THEN
            IF ((IGP /= 31).AND.(IGP /= 32))THEN
! DEPENDS ON: addrln
              CALL ADDRLN(IGP,halo_type,LEN_PRIMIN,local_data)
              IN_S(2,MODL,ISEC,IITM)=LEN_PRIMIN*IZ_IN*IP_IN
            ELSE
! DEPENDS ON: ocnvol
              CALL OCNVOL(LENO,LIST_S(st_input_bottom,ISTART),          &
     &                         LIST_S(st_input_top   ,ISTART))
              IN_S(2,MODL,ISEC,IITM)=LENO
            END IF
          END IF

        END IF ! At least one stash record for m,s,i

      END DO   ! Items
      END DO   ! Sections
      END DO   ! Models
#if defined(MPP)
!
! DEPENDS ON: change_decomposition
      CALL CHANGE_DECOMPOSITION(orig_decomp,ErrorStatus)

      IF(ErrorStatus >  0) THEN
         CMESSAGE='INPUTL: ERROR in original MPP decomposition'
         write(6,*) CMESSAGE
         GOTO 999
      ENDIF
#endif

 999  CONTINUE

      RETURN
      END SUBROUTINE INPUTL

!- End of subroutine code -------------------------------------------


!+Determine whether a levels list is a duplicate of another levels list
! Subroutine Interface:


!- End of subroutine code ----------------------------------------


!+Determine whether a pseudo lev list is a duplicate of another one
! Subroutine Interface:


!- End of subroutine code ---------------------------------------
#endif
