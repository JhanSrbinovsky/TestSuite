#if defined(A31_1A)
#if defined(ATMOS) && !defined(GLOBAL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine READ_ATMOS_LBCS
!
! Purpose : Reads in all the Atmos LBCs
!
! Author :  P.Burton
!
! Model                Modification history from model version 5.2
! version    Date
! 5.2        28/09/00  New deck at 5.2                    P.Burton
! 5.3   12/10/01  extended haloes for exner in lbcs    Terry Davies
! 5.5   03/02/03  Include qcf2,qrain,qgraup lbcs. R.M.Forbes
! 6.0   29/07/03  Include cloud fractions in lbcs. Damian Wilson
! 6.2   05/01/05  Change mcr/pc2 lbc reads to be dependent on _lbc
!                 logicals only.   R.Forbes
! 6.2   01/10/04  Include murk aerosol lbcs.  R.M.Forbes
!
! ---------------------------------------------------------

      SUBROUTINE READ_ATMOS_LBCS(                                       &
     &  LENRIM, global_LENRIM,                                          &
     &  MODEL_LEVELS,WET_LEVELS,                                        &
     &  TR_VARS,TR_LEVELS,                                              &
     &  L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, l_pc2,                   &
     &  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, l_pc2_lbc,   &
     &  L_murk, L_murk_lbc,                                             &
     &  LEN1_LOOKUP,LEN_FIXHD,NFTIN, N_LOOKUPS,LOOKUP, FIXHD,           &
     &  U_LBC,V_LBC,W_LBC,RHO_LBC,THETA_LBC,Q_LBC,QCL_LBC,QCF_LBC,      &
     &  QCF2_LBC, QRAIN_LBC, QGRAUP_LBC,                                &
     &  CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC,                      &
     &  EXNER_LBC, U_ADV_LBC, V_ADV_LBC, W_ADV_LBC,                     &
     &  MURK_LBC, TRACER_LBCS,                                          &
#include "argppx.h"
     &  ICODE,CMESSAGE)

      IMPLICIT NONE

! Parameters required for argument declarations
#include "parparm.h"

! Arguments:

      INTEGER                                                           &
     &  LENRIM(Nfld_max,NHalo_Max)                                      &
                                    ! IN : Size of a level of LBC
     &, global_LENRIM(Nfld_max,NHalo_Max)                               &
                                    ! IN : Full (disk) size of a
                                    !      level of LBC
     &, MODEL_LEVELS                                                    &
                                    ! IN : Number of model levels
     &, WET_LEVELS                                                      &
                                    ! IN : Number of wet model levels
     &, TR_VARS                                                         &
                                    ! IN : Number of tracers
     &, TR_LEVELS                                                       &
                                    ! IN : Number of tracer levels
     &, LEN1_LOOKUP                                                     &
                                    ! IN : Size of LOOKUP header
     &, LEN_FIXHD                                                       &
                                    ! IN : Size of fixed length header
     &, NFTIN                                                           &
                                    ! IN : Unit for LBC file
     &, N_LOOKUPS                                                       &
                                    ! IN  : Number of lookup headers
     &, LOOKUP(LEN1_LOOKUP,N_LOOKUPS)                                   &
                                    ! IN : LOOKUP headers
     &, FIXHD(LEN_FIXHD)            ! IN : Fixed header

      Logical, Intent (In) ::                                           &
     &  L_mcr_qcf2                                                      &
                          ! true if second cloud ice active
     &, L_mcr_qrain                                                     &
                          ! true if rain active
     &, L_mcr_qgraup                                                    &
                          ! true if graupel active
     &, L_mcr_qcf2_lbc                                                  &
                          ! true if second cloud ice lbcs active
     &, L_mcr_qrain_lbc                                                 &
                          ! true if rain lbcs active
     &, L_mcr_qgraup_lbc                                                &
                          ! true if graupel lbcs active
     &, L_pc2                                                           &
                          ! true if prognostic cloud fractions used
     &, L_pc2_lbc                                                       &
                          ! true if cloud fractions in lbcs
     &, L_murk                                                          &
                          ! true if prognostic murk aerosol active
     &, L_murk_lbc        ! true if murk aerosol in lbcs

      REAL                                                              &
     &  U_LBC(LENRIM(fld_type_u,halo_type_extended),                    &
     &        MODEL_LEVELS)                                             &
                                    ! OUT : U LBC
     &, V_LBC(LENRIM(fld_type_v,halo_type_extended),                    &
     &        MODEL_LEVELS)                                             &
                                    ! OUT : V LBC
     &, W_LBC(LENRIM(fld_type_p,halo_type_extended),                    &
     &        0:MODEL_LEVELS)                                           &
                                    ! OUT : V LBC
     &, RHO_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
     &          MODEL_LEVELS)                                           &
                                    ! OUT : Rho LBC
     &, THETA_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &            MODEL_LEVELS)                                         &
                                    ! OUT : Theta LBC
     &, Q_LBC(LENRIM(fld_type_p,halo_type_extended),                    &
     &        WET_LEVELS)                                               &
                                    ! OUT : Q LBC
     &, QCL_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
     &          WET_LEVELS)                                             &
                                    ! OUT : QCL LBC
     &, QCF_LBC(LENRIM(fld_type_p,halo_type_extended),                  &
     &          WET_LEVELS)                                             &
                                    ! OUT : QCL LBC
     &, QCF2_LBC(LENRIM(fld_type_p,halo_type_extended),                 &
     &          WET_LEVELS)                                             &
                                    ! OUT : QCF2 LBC
     &, QRAIN_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &          WET_LEVELS)                                             &
                                    ! OUT : QRAIN LBC
     &, QGRAUP_LBC(LENRIM(fld_type_p,halo_type_extended),               &
     &          WET_LEVELS)                                             &
                                    ! OUT : QGRAUP LBC
     &, CF_BULK_LBC(LENRIM(fld_type_p,halo_type_extended),              &
     &          WET_LEVELS)                                             &
                                    ! OUT : CF_BULK LBC
     &, CF_LIQUID_LBC(LENRIM(fld_type_p,halo_type_extended),            &
     &          WET_LEVELS)                                             &
                                    ! OUT : CF_LIQUID LBC
     &, CF_FROZEN_LBC(LENRIM(fld_type_p,halo_type_extended),            &
     &          WET_LEVELS)                                             &
                                    ! OUT : CF_FROZEN LBC
     &, EXNER_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &            MODEL_LEVELS+1)                                       &
                                    ! OUT : Exner LBC
     &, U_ADV_LBC(LENRIM(fld_type_u,halo_type_extended),                &
     &            MODEL_LEVELS)                                         &
                                    ! OUT : U_ADV LBC
     &, V_ADV_LBC(LENRIM(fld_type_v,halo_type_extended),                &
     &            MODEL_LEVELS)                                         &
                                    ! OUT : V_ADV LBC
     &, W_ADV_LBC(LENRIM(fld_type_p,halo_type_extended),                &
     &            0:MODEL_LEVELS)                                       &
                                    ! OUT : W LBC
     &, MURK_LBC(LENRIM(fld_type_p,halo_type_single),                   &
     &           MODEL_LEVELS)                                          &
                                       ! OUT : MURK LBC
     &, TRACER_LBCS(LENRIM(fld_type_p,halo_type_extended),              &
     &              TR_LEVELS,TR_VARS)
                                    ! OUT : Tracer LBCs

#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "clookadd.h"

      INTEGER                                                           &
     &  ICODE                           ! Error code

      CHARACTER*(80)                                                    &
     &  CMESSAGE                        ! Error message

! Local variables

      INTEGER                                                           &
     &  tracer                                                          &
                 ! loop counter over tracer variables
     &, lbc_num  ! number of active lbcs (used for microphysics lbcs)

! ---------------------------------------------------------

! Check the field about to be read is U_LBC - if not, something
! has gone wrong.
      IF ( LOOKUP(ITEM_CODE,1)  /=  31002) THEN
        WRITE(6,*) 'Expected to find U_LBC (31002) in LBC file '
        WRITE(6,*) 'But found ',LOOKUP(ITEM_CODE,1),' instead.'
        ICODE=1
        CMESSAGE='READ_ATMOS_LBCS : Found wrong field in LBC file'
        GOTO 9999
      ENDIF

! Now read in all the fields

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,1),LEN1_LOOKUP,U_LBC,            &
     &              global_LENRIM(fld_type_u,halo_type_extended)*       &
     &                MODEL_LEVELS,                                     &
     &              FIXHD,                                              &
#include "argppx.h"
     &              ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading U_LBC'
        GOTO 9999
      ENDIF


! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,2),LEN1_LOOKUP,V_LBC,            &
     &              global_LENRIM(fld_type_v,halo_type_extended)*       &
     &                MODEL_LEVELS,                                     &
     &              FIXHD,                                              &
#include "argppx.h"
     &              ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading V_LBC'
        GOTO 9999
      ENDIF

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,3),LEN1_LOOKUP,W_LBC,            &
     &              global_LENRIM(fld_type_p,halo_type_extended)*       &
     &                (MODEL_LEVELS+1),                                 &
     &              FIXHD,                                              &
#include "argppx.h"
     &              ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading W_LBC'
        GOTO 9999
      ENDIF

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,4),LEN1_LOOKUP,RHO_LBC,          &
     &              global_LENRIM(fld_type_p,halo_type_extended)*       &
     &                MODEL_LEVELS,                                     &
     &              FIXHD,                                              &
#include "argppx.h"
     &              ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading RHO_LBC'
        GOTO 9999
      ENDIF

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,5),LEN1_LOOKUP,THETA_LBC,        &
     &              global_LENRIM(fld_type_p,halo_type_extended)*       &
     &                MODEL_LEVELS,                                     &
     &              FIXHD,                                              &
#include "argppx.h"
     &              ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading THETA_LBC'
        GOTO 9999
      ENDIF

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,6),LEN1_LOOKUP,Q_LBC,            &
     &              global_LENRIM(fld_type_p,halo_type_extended)*       &
     &                WET_LEVELS,                                       &
     &              FIXHD,                                              &
#include "argppx.h"
     &              ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading Q_LBC'
        GOTO 9999
      ENDIF

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,7),LEN1_LOOKUP,QCL_LBC,          &
     &              global_LENRIM(fld_type_p,halo_type_extended)*       &
     &                WET_LEVELS,                                       &
     &              FIXHD,                                              &
#include "argppx.h"
     &              ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading QCL_LBC'
        GOTO 9999
      ENDIF

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,8),LEN1_LOOKUP,QCF_LBC,          &
     &              global_LENRIM(fld_type_p,halo_type_extended)*       &
     &                WET_LEVELS,                                       &
     &              FIXHD,                                              &
#include "argppx.h"
     &              ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading QCF_LBC'
        GOTO 9999
      ENDIF

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,9),LEN1_LOOKUP,EXNER_LBC,        &
     &              global_LENRIM(fld_type_p,halo_type_extended)*       &
     &                (MODEL_LEVELS+1),                                 &
     &              FIXHD,                                              &
#include "argppx.h"
     &              ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading EXNER_LBC'
        GOTO 9999
      ENDIF

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,10),LEN1_LOOKUP,U_ADV_LBC,       &
     &              global_LENRIM(fld_type_u,halo_type_extended)*       &
     &                MODEL_LEVELS,                                     &
     &              FIXHD,                                              &
#include "argppx.h"
     &              ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading U_ADV_LBC'
        GOTO 9999
      ENDIF


! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,11),LEN1_LOOKUP,V_ADV_LBC,       &
     &              global_LENRIM(fld_type_v,halo_type_extended)*       &
     &                MODEL_LEVELS,                                     &
     &              FIXHD,                                              &
#include "argppx.h"
     &              ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading V_ADV_LBC'
        GOTO 9999
      ENDIF

! DEPENDS ON: readflds
      CALL READFLDS(NFTIN,1,1,LOOKUP(1,12),LEN1_LOOKUP,W_ADV_LBC,       &
     &              global_LENRIM(fld_type_p,halo_type_extended)*       &
     &                (MODEL_LEVELS+1),                                 &
     &              FIXHD,                                              &
#include "argppx.h"
     &              ICODE,CMESSAGE)

      IF (ICODE  >   0) THEN
        WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading W_ADV_LBC'
        GOTO 9999
      ENDIF

      ! set number of active lbcs for lookup table
      lbc_num = 12

      ! ----------------------------------------------------------------
      ! Read ice crystal (qcf2) lbcs if present in lbc file
      ! ----------------------------------------------------------------

      IF (L_mcr_qcf2_lbc) THEN

        ! qcf2 is in input lbc file
        ! Increment number of lbcs
        lbc_num = lbc_num + 1
! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
     &                LEN1_LOOKUP,QCF2_LBC,                             &
     &                global_LENRIM(fld_type_p,halo_type_extended)*     &
     &                WET_LEVELS,                                       &
     &                FIXHD,                                            &
#include "argppx.h"
     &                ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading QCF2_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_mcr_qcf2) THEN

        ! qcf2 prognostic is active but qcf2 lbcs are not in input file
        ! set lbc array to zero
        QCF2_LBC(:,:) = 0.0

      ENDIF  ! on L_mcr_qcf2, L_mcr_qcf2_lbc

      ! ----------------------------------------------------------------
      ! Read rain lbcs if present in lbc file
      ! ----------------------------------------------------------------

      IF (L_mcr_qrain_lbc) THEN

        ! qrain is in input lbc file
        ! Increment number of lbcs
        lbc_num = lbc_num + 1
! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
     &                LEN1_LOOKUP,QRAIN_LBC,                            &
     &                global_LENRIM(fld_type_p,halo_type_extended)*     &
     &                WET_LEVELS,                                       &
     &                FIXHD,                                            &
#include "argppx.h"
     &                ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading QRAIN_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_mcr_qrain) THEN

        ! qrain is active but qrain lbcs are not in input file
        ! set lbc array to zero
        QRAIN_LBC(:,:) = 0.0

      ENDIF  ! on L_mcr_qrain, L_mcr_qrain_lbc

      ! ----------------------------------------------------------------
      ! Read graupel lbcs if present in lbc file
      ! ----------------------------------------------------------------

      IF (L_mcr_qgraup_lbc) THEN

        ! qgraup is in input lbc file
        ! Increment number of lbcs
        lbc_num = lbc_num + 1
! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
     &                LEN1_LOOKUP,QGRAUP_LBC,                           &
     &                global_LENRIM(fld_type_p,halo_type_extended)*     &
     &                WET_LEVELS,                                       &
     &                FIXHD,                                            &
#include "argppx.h"
     &                ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading QGRAUP_LBC'
          GOTO 9999
        ENDIF

      ELSE  IF (L_mcr_qgraup) THEN

        ! qgraup is active but qgraup lbcs are not in input file
        ! set lbc array to zero
        QGRAUP_LBC(:,:) = 0.0

      ENDIF  ! on L_mcr_qgraup, L_mcr_qgraup_lbc

      ! ----------------------------------------------------------------
      ! Read cloud fraction lbcs present in lbc file
      ! ----------------------------------------------------------------

      ! Set cloud fraction lbcs
      IF (L_pc2_lbc) THEN

        ! Cloud fractions variables are in input lbc file
        ! Increment number of lbcs
        lbc_num = lbc_num + 1
! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
     &                LEN1_LOOKUP,CF_BULK_LBC,                          &
     &                global_LENRIM(fld_type_p,halo_type_extended)*     &
     &                WET_LEVELS,                                       &
     &                FIXHD,                                            &
#include "argppx.h"
     &                ICODE,CMESSAGE)
        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading CF_BULK_LBC'
          GOTO 9999
        ENDIF

       lbc_num = lbc_num + 1
! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
     &                LEN1_LOOKUP,CF_LIQUID_LBC,                        &
     &                global_LENRIM(fld_type_p,halo_type_extended)*     &
     &                WET_LEVELS,                                       &
     &                FIXHD,                                            &
#include "argppx.h"
     &                ICODE,CMESSAGE)
        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading CF_LIQUID_LBC'
          GOTO 9999
        ENDIF

        lbc_num = lbc_num + 1
! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
     &                LEN1_LOOKUP,CF_FROZEN_LBC,                        &
     &                global_LENRIM(fld_type_p,halo_type_extended)*     &
     &                WET_LEVELS,                                       &
     &                FIXHD,                                            &
#include "argppx.h"
     &                ICODE,CMESSAGE)
        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading CF_FROZEN_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_pc2) THEN

        ! Cloud fraction prognostics are active but cloud fraction lbcs
        ! are not in input file so set the lbc arrays to zero.
        CF_BULK_LBC(:,:)   = 0.0
        CF_LIQUID_LBC(:,:) = 0.0
        CF_FROZEN_LBC(:,:) = 0.0

      ENDIF  ! on L_pc2, L_pc2_lbc

      ! ----------------------------------------------------------------
      ! Read murk aerosol lbcs if expected in lbc file
      ! ----------------------------------------------------------------

      IF (L_murk_lbc) THEN   ! murk lbcs are in lbc file

        ! Increment number of lbcs
        lbc_num = lbc_num + 1

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num),                      &
     &                LEN1_LOOKUP,MURK_LBC,                             &
     &                global_LENRIM(fld_type_p,halo_type_single)*       &
     &                MODEL_LEVELS,                                     &
     &                FIXHD,                                            &
#include "argppx.h"
     &                ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading MURK_LBC'
          GOTO 9999
        ENDIF

      ELSE IF (L_murk) THEN

        ! murk prognostic is active but murk lbcs are not in input file
        ! set lbc array to zero
        MURK_LBC(:,:) = 0.0

      ENDIF ! on L_murk_lbc, L_murk

      ! ----------------------------------------------------------------
      ! Read tracer lbcs if expected in lbc file
      ! ----------------------------------------------------------------

      DO tracer=1,TR_VARS

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP(1,lbc_num+tracer),               &
     &                LEN1_LOOKUP, TRACER_LBCS(1,1,tracer),             &
     &                global_LENRIM(fld_type_p,halo_type_extended)*     &
     &                  TR_LEVELS,                                      &
     &                FIXHD,                                            &
#include "argppx.h"
     &                ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'READ_ATMOS_LBCS : Problem reading Tracer LBC ',   &
     &               tracer
          GOTO 9999
        ENDIF

      ENDDO  ! tracer


 9999 CONTINUE
      RETURN
      END SUBROUTINE READ_ATMOS_LBCS
#endif
#endif
