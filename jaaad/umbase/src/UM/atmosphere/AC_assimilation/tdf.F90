#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Main routine for TDF scheme.

      SUBROUTINE TDF (                                                  &
#include "argd1.h"
#include "argptra.h"
#include "argduma.h"
#include "argppx.h"
     &                 Weight,                                          &
                                     ! in
     &                 L_MiddleCall,                                    &
                                     ! in
     &                 L_LastCall,                                      &
                                     ! in
     &                 D1_TDF )      ! inout

! Description:
!
!   Main routine for Temporal Digital Filtering (TDF) scheme.
!
!   Fields filtered:
!
!     Always:     u, v, w, rho, theta, exner
!     Optionally: q, qCL, qCF, u_adv, v_adv, w_adv
!
!   Filtering of the humidity fields q, qCL and qCF is requested via
!   the logicals L_TDF_FilterQ, L_TDF_FilterQCL and L_TDF_FilterQCF in
!   the temporal filtering namelist RUN_TFilt. By default, all three
!   are set to .FALSE..
!
!   For the advected wind components u_adv, v_adv and w_adv there are
!   three options:
!
!     1. Leave alone.
!     2. Filter directly.
!     3. Set to filtered versions of non-advected winds.
!
!   The choice is specified via the namelist variable TDF_AdvWindOpt.
!   By default, the first option is chosen.
!
!   If either qCL or qCF is filtered, or the model is being run
!   without physics, the resulting TDF dump will generally be
!   unsuitable for starting runs involving physics unless problems
!   and inconsistences are removed from the cloud variables. The
!   required modifications are requested by setting the logical
!   L_TDF_ModifyCloud to .TRUE..
!
! Method:
!
!   On each call, for each field to be filtered, add a fraction of the
!   main D1 field to the corresponding field in the TDF field array.
!
!   On middle call, write out ordinary model dump.
!
!   On final call, overwrite fields that have been modified.
!
!
! Current Code Owner: Adam Clayton.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   13/04/00   Original code. Adam Clayton
!   5.2   13/10/00   Revise cloud modification option. Adam Clayton
!   5.3   08/06/01   Declare N_INTERNAL_MODEL before uisng.
!                    A van der Wal
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! Declarations:

      IMPLICIT NONE

! Common blocks:

#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
#include "typptra.h"
#include "typduma.h"
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "clookadd.h"
#include "cprintst.h"
#include "ctfilt.h"

! Subroutine arguments:

      REAL,    INTENT(IN) :: Weight       ! TDF weight.
      LOGICAL, INTENT(IN) :: L_MiddleCall ! .TRUE. for middle call.
      LOGICAL, INTENT(IN) :: L_LastCall   ! .TRUE. for last   call.

      REAL, INTENT(INOUT) :: D1_TDF(D1_TFilt_len) ! TDF field array.

! Local variables:

      LOGICAL :: NeedFieldNum

      INTEGER :: ICode,                                                 &
     &           i, j, k,                                               &
     &           BufLen,                                                &
     &           SubModel_num,                                          &
     &           number_of_data_words_in_memory,                        &
     &           number_of_data_words_on_disk,                          &
     &           disk_address,                                          &
     &           FirstLev,                                              &
     &           Addr,                                                  &
     &           Addr1, Addr2,                                          &
     &           Addr3, Addr4,                                          &
     &           IM_index,                                              &
     &           row_u_adv,                                             &
     &           row_v_adv,                                             &
     &           row_w_adv

      REAL    :: lower_limit

      INTEGER :: NumLevs    (18),                                       &
     &           FieldLen   (18),                                       &
     &           STASH_code (18),                                       &
     &           Addr_D1_TDF(18),                                       &
     &           Addr_D1    (12),                                       &
     &           FieldNum   (18)

      LOGICAL :: FilterField(12)

      CHARACTER(17) ::                                                  &
     &           FieldDesc  (18)

      LOGICAL, SAVE ::                                                  &
     &           FirstWrite ! .TRUE. if D1_TDF yet to be written to.
      DATA       FirstWrite / .TRUE. /

      CHARACTER(80) :: CMessage

      CHARACTER(*)  :: RoutineName
      PARAMETER      ( RoutineName='TDF' )

! External subroutines called:

      EXTERNAL FILE_OPEN,  EReport,              UM_WRITDUMP,           &
     &         FILE_CLOSE, SET_DUMPFILE_ADDRESS, WRITFLDS,              &
     &         READFLDS

!- End of header ------------------------------------------------------

!----------------------------------------------------------------------
! [1]: Info for fields than may be modified.
!----------------------------------------------------------------------

      FieldDesc = (/ 'u                ', 'v                ',          &
     &               'w                ', 'rho              ',          &
     &               'theta            ', 'exner            ',          &
     &               'q                ', 'qCL              ',          &
     &               'qCF              ', 'u_adv            ',          &
     &               'v_adv            ', 'w_adv            ',          &
     &               'qCL              ', 'qCF              ',          &
     &               'area_cloud_frac  ', 'bulk_cloud_frac  ',          &
     &               'cloud_frac_liquid', 'cloud_frac_frozen' /)

      NumLevs(:)      = model_levels
      NumLevs(3)      = model_levels+1 ! w
      NumLevs(6)      = model_levels+1 ! exner
      NumLevs(7:9)    = wet_levels     ! q, qCL, qCF
      NumLevs(12)     = model_levels+1 ! w_adv
      NumLevs(13:14)  = wet_levels     ! qCL, qCF in TDF dump
      NumLevs(15:18)  = wet_levels     ! cloud fractions in TDF dump

      FieldLen(1)     = u_off_size       *  model_levels
      FieldLen(2)     = v_off_size       *  model_levels
      FieldLen(3)     = theta_off_size   * (model_levels+1)
      FieldLen(4:5)   = theta_off_size   *  model_levels
      FieldLen(6)     = theta_off_size   * (model_levels+1)
      FieldLen(7:9)   = theta_halo_size  *  wet_levels
      FieldLen(10)    = u_halo_size      *  model_levels
      FieldLen(11)    = v_halo_size      *  model_levels
      FieldLen(12)    = theta_halo_size  * (model_levels+1)
      FieldLen(13:14) = theta_halo_size  *  wet_levels
      FieldLen(15:18) = theta_field_size *  wet_levels

      Addr_D1 = (/ JU(1),     JV(1),                                    &
     &             JW(0),     JRHO(1),                                  &
     &             JTHETA(1), JEXNER_RHO_LEVELS(1),                     &
     &             JQ(1),     JQCL(1),                                  &
     &             JQCL(1),   JU_ADV(1),                                &
     &             JV_ADV(1), JW_ADV(0) /)

      STASH_code = (/ 2,   3,   150, 253, 4,  255, 10,  254, 12,        &
     &                256, 257, 258, 254, 12, 265, 266, 267, 268 /)

! Filter field?

      FilterField(1:6) = .TRUE.
      FilterField(7)   = L_TDF_FilterQ
      FilterField(8)   = L_TDF_FilterQCL
      FilterField(9)   = L_TDF_FilterQCF

      IF (TDF_AdvWindOpt == DirectFilt) THEN
        FilterField(10:12) = .TRUE.
      ELSE
        FilterField(10:12) = .FALSE.
      END IF

! Find field addresses in TDF field array.

      Addr = 1

      DO i = 1, 12
        IF ( FilterField(i) ) THEN
          Addr_D1_TDF(i) = Addr
          Addr           = Addr + FieldLen(i)
        END IF
      END DO

      IF ( TDF_AdvWindOpt == CopyFiltUVW ) THEN
        Addr_D1_TDF(10) = Addr_D1_TDF(1) ! u_adv -> u
        Addr_D1_TDF(11) = Addr_D1_TDF(2) ! v_adv -> v
        Addr_D1_TDF(12) = Addr_D1_TDF(3) ! w_adv -> w
      END IF

      Addr = 1 ! Reuse workspace for modification of cloud variables.

      DO i = 13, 18
          Addr_D1_TDF(i) = Addr
          Addr           = Addr + FieldLen(i)
      END DO

!----------------------------------------------------------------------
! [2]: Accumulate fields in TDF field array.
!----------------------------------------------------------------------

      IF (ABS(Weight) < Weight_LL) THEN

        IF (PrintStatus >= PrStatus_Normal) THEN
          WRITE (6,*) ''
          WRITE (6,'(A,F13.10,A)') 'TDF: Weight (', Weight,             &
     &                             ') too small to bother with.'
          WRITE (6,*) ''
        END IF

      ELSE

        DO i = 1, 12
          IF ( FilterField(i) ) THEN

            Addr1 = Addr_D1_TDF(i)
            Addr2 = Addr_D1_TDF(i) + FieldLen(i) - 1
            Addr3 = Addr_D1(i)
            Addr4 = Addr_D1(i)     + FieldLen(i) - 1

            IF (FirstWrite) THEN
              D1_TDF(Addr1:Addr2) = D1(Addr3:Addr4) * Weight
            ELSE
              D1_TDF(Addr1:Addr2) = D1_TDF(Addr1:Addr2) +               &
     &                              D1(Addr3:Addr4) * Weight
            END IF

          END IF
        END DO

        FirstWrite = .FALSE.

      END IF ! (ABS(Weight) < Weight_LL)

!----------------------------------------------------------------------
! [3]: On middle call, write model dump. Fields to be modified will be
!      overwritten on the last call.
!----------------------------------------------------------------------

      IF (L_MiddleCall) THEN

! DEPENDS ON: file_open
        CALL FILE_OPEN(TDF_unit, 'TDF_dump', 8, 1, 0, ICode)
        IF (ICode > 0) THEN
          CMessage = 'Error opening TDF dump.'
! DEPENDS ON: ereport
          CALL EReport (RoutineName, ICode, CMessage)
        END IF

        SubModel_num = SUBMODEL_FOR_SM(ATMOS_SM)

        BufLen = 0
        DO i = 1, A_LEN2_LOOKUP
          BufLen = MAX(BufLen, A_LOOKUP(LBLREC, i))
        END DO

        IF (PrintStatus >= PrStatus_Normal) THEN
          WRITE (6,*) ''
          WRITE (6,*) 'TDF: Writing TDF dump on unit ', TDF_unit
          WRITE (6,*) '     Fields to be filtered will be'              &
     &                  //' overwritten later.'
          WRITE (6,*) ''
        END IF

! DEPENDS ON: um_writdump
        CALL UM_WRITDUMP (                                              &
     &    TDF_unit,                                                     &
                                                        ! in
     &    A_FIXHD,      LEN_FIXHD,                                      &
                                                        ! in
     &    A_INTHD,      A_LEN_INTHD,                                    &
                                                        ! in
     &    A_REALHD,     A_LEN_REALHD,                                   &
                                                        ! in
     &    A_LEVDEPC,    A_LEN1_LEVDEPC, A_LEN2_LEVDEPC,                 &
                                                        ! in
     &    A_ROWDEPC,    A_LEN1_ROWDEPC, A_LEN2_ROWDEPC,                 &
                                                        ! in
     &    A_COLDEPC,    A_LEN1_COLDEPC, A_LEN2_COLDEPC,                 &
                                                        ! in
     &    A_FLDDEPC,    A_LEN1_FLDDEPC, A_LEN2_FLDDEPC,                 &
                                                        ! in
     &    A_EXTCNST,    A_LEN_EXTCNST,                                  &
                                                        ! in
     &    A_DUMPHIST,   LEN_DUMPHIST,                                   &
                                                        ! in
     &    A_CFI1,       A_LEN_CFI1,                                     &
                                                        ! in
     &    A_CFI2,       A_LEN_CFI2,                                     &
                                                        ! in
     &    A_CFI3,       A_LEN_CFI3,                                     &
                                                        ! in
     &    A_LOOKUP,     LEN1_LOOKUP,    A_LEN2_LOOKUP,                  &
                                                        ! in
#if defined(MPP)
     &    A_MPP_LOOKUP, MPP_LEN1_LOOKUP,                                &
                                                        ! in
#endif
     &    BufLen,                                                       &
                                                        ! in
#include "argppx.h"
     &    ATMOS_SM,                                                     &
                                                        ! in
     &    NO_OBJ_D1  (SubModel_num),                                    &
                                                        ! in
     &    D1_ADDR(1,1,SubModel_num),                                    &
                                                        ! in
     &    A_LEN_DATA,                                                   &
                                                        ! in
     &    D1 )                                          ! in

! DEPENDS ON: file_close
        CALL FILE_CLOSE(TDF_unit, 'TDF_dump', 8, 0, 0, ICode)
        IF (ICode > 0) THEN
          CMessage = 'Error closing TDF dump.'
! DEPENDS ON: ereport
          CALL EReport (RoutineName, ICode, CMessage)
        END IF

      END IF  ! L_MiddleCall

!----------------------------------------------------------------------
! [4]: On last call, overwrite fields we chose to modify.
!----------------------------------------------------------------------

      IF (L_LastCall) THEN

! Enforce lower limits for humidities.

        DO i = 7, 9

          IF (i == 7) THEN
            lower_limit = q_min ! q
          ELSE
            lower_limit = 0.0   ! qCL, qCF
          END IF

          IF ( FilterField(i) ) THEN
            DO j = Addr_D1_TDF(i), ADDR_D1_TDF(i) + FieldLen(i) - 1
              IF (D1_TDF(j) < lower_limit) D1_TDF(j) = lower_limit
            END DO
          END IF

        END DO

! If required, reset stratospheric humidities.

        IF (L_TDF_CallStratQ) THEN

          IF (PrintStatus >= PrStatus_Normal) THEN
            WRITE (6,*) ''
            WRITE (6,*) 'TDF: Currently no StratQ for new dynamics.'
            WRITE (6,*) ''
          END IF

        END IF

! Open TDF dump.

! DEPENDS ON: file_open
        CALL FILE_OPEN(TDF_unit, 'TDF_dump', 8, 1, 0, ICode)
        IF (ICode > 0) THEN
          CMessage = 'Error opening TDF dump.'
! DEPENDS ON: ereport
          CALL EReport (RoutineName, ICode, CMessage)
        END IF

! DEPENDS ON: set_dumpfile_address
        CALL SET_DUMPFILE_ADDRESS (                                     &
     &    A_FIXHD,     LEN_FIXHD,                                       &
                                          ! in
     &    A_LOOKUP,                                                     &
                                          ! inout
     &    LEN1_LOOKUP, A_LEN2_LOOKUP,                                   &
                                          ! in
     &    number_of_data_words_in_memory,                               &
                                          ! out
     &    number_of_data_words_on_disk,                                 &
                                          ! out
     &    disk_address )                  ! out

        BufLen = 0
        DO i = 1, A_LEN2_LOOKUP
          BufLen = MAX(BUFLEN, A_LOOKUP(LBLREC, i))
        END DO

! Find field numbers in dump.

        FieldNum(:) = 0

        DO j = 1, A_LEN2_LOOKUP
          DO i = 1, 18

            IF (i == 3 .OR. i == 12) THEN ! w, w_adv
              FirstLev = 9999
            ELSE
              FirstLev = 1
            END IF

            IF ( A_LOOKUP(ITEM_CODE, j) == STASH_code(i) .AND.          &
     &           A_LOOKUP(LBLEV,     j) == FirstLev )  THEN
              FieldNum(i) = j
            END IF

          END DO
        END DO

        NeedFieldNum = ( (i <= 12 .AND. FilterField(i)) .OR.            &
     &                   (i >= 13 .AND. L_TDF_ModifyCloud) )

        DO i = 1, 18
          IF (NeedFieldNum .AND. FieldNum(i) == 0) THEN
            ICode = 1
            WRITE (CMessage, FMT='(A,A,A)') 'Could not find ',          &
     &        FieldDesc(i), ' in TDF dump.'
! DEPENDS ON: ereport
            CALL EReport (RoutineName, ICode, CMessage)
          END IF
        END DO

! Overwrite fields we chose to filter.

        IF (TDF_AdvWindOpt == CopyFiltUVW) THEN

         ! Advected wind fields have larger halos than their
         ! non-advected counterparts, so temporarily change their
         ! PPXI halo entries to those of the normal wind fields so
         ! we can correctly overwrite the advected fields in the
         ! dump with the filtered versions of the normal wind fields.

          IM_index  = INTERNAL_MODEL_INDEX(A_IM)
          row_u_adv = PPXPTR(IM_index, 0, 256)
          row_v_adv = PPXPTR(IM_index, 0, 257)
          row_w_adv = PPXPTR(IM_index, 0, 258)

          PPXI(row_u_adv, ppx_halo_type) = halo_type_single
          PPXI(row_v_adv, ppx_halo_type) = halo_type_single
          PPXI(row_w_adv, ppx_halo_type) = halo_type_single

        END IF

        WRITE (6,*) ''

        DO i = 1, 12

          IF ( FilterField(i) ) THEN

            IF (PrintStatus >= PrStatus_Normal)                         &
     &        WRITE (6,*) 'TDF: Overwriting ', FieldDesc(i)

! DEPENDS ON: writflds
            CALL WRITFLDS (TDF_unit,                                    &
                                                   ! in
     &                     NumLevs(i),                                  &
                                                   ! in
     &                     FieldNum(i),                                 &
                                                   ! in
     &                     A_LOOKUP,                                    &
                                                   ! in
     &                     LEN1_LOOKUP,                                 &
                                                   ! in
     &                     D1_TDF(Addr_D1_TDF(i)),                      &
                                                   ! in
     &                     BufLen,                                      &
                                                   ! in
     &                     A_FIXHD,                                     &
                                                   ! in
#include "argppx.h"
     &                     ICode,                                       &
                                                   ! out
     &                     CMessage)               ! out

            IF (ICode > 0) THEN
              WRITE (CMessage, FMT='(A,A)') 'Error overwriting ',       &
     &                                      FieldDesc(i)
! DEPENDS ON: ereport
              CALL EReport (RoutineName, ICode, CMessage)
            END IF

          END IF ! ( FilterField(i) )

        END DO ! i

        IF (TDF_AdvWindOpt == CopyFiltUVW) THEN

         ! Restore PPXI halo entries for advected wind fields.
          PPXI(row_u_adv, ppx_halo_type) = halo_type_extended
          PPXI(row_v_adv, ppx_halo_type) = halo_type_extended
          PPXI(row_w_adv, ppx_halo_type) = halo_type_extended

        END IF

! If required, remove problems and inconsistencies from cloud
! variables in TDF dump so that it becomes suitable for starting
! forecasts including physics.

        IF (L_TDF_ModifyCloud) THEN

         ! Read cloud variables from TDF dump.
          DO i = 13, 18

! DEPENDS ON: readflds
              CALL READFLDS (TDF_unit,                                  &
                                                     ! in
     &                       NumLevs(i),                                &
                                                     ! in
     &                       FieldNum(i),                               &
                                                     ! in
     &                       A_LOOKUP,                                  &
                                                     ! in
     &                       LEN1_LOOKUP,                               &
                                                     ! in
     &                       D1_TDF(Addr_D1_TDF(i)),                    &
                                                     ! out
     &                       BufLen,                                    &
                                                     ! in
     &                       A_FIXHD,                                   &
                                                     ! in
#include "argppx.h"
     &                       ICode,                                     &
                                                     ! out
     &                       CMessage)               ! out

              IF (ICode > 0) THEN
                WRITE (CMessage, FMT='(A,A)') 'Error reading ',         &
     &                                        FieldDesc(i)
! DEPENDS ON: ereport
                CALL EReport (RoutineName, ICode, CMessage)
              END IF

          END DO

         ! Modify cloud variables.
! DEPENDS ON: modifycloudvars
          CALL ModifyCloudVars (                                        &
     &      D1_TDF(Addr_D1_TDF(13)),                                    &
                                      ! inout, qCL
     &      D1_TDF(Addr_D1_TDF(14)),                                    &
                                      ! inout, qCF
     &      D1_TDF(Addr_D1_TDF(15)),                                    &
                                      ! inout, area_cloud_fraction
     &      D1_TDF(Addr_D1_TDF(16)),                                    &
                                      ! inout, bulk_cloud_fraction
     &      D1_TDF(Addr_D1_TDF(17)),                                    &
                                      ! inout, cloud_fraction_liquid
     &      D1_TDF(Addr_D1_TDF(18)) ) ! inout, cloud_fraction_frozen

         ! Write modified cloud variables to TDF dump.
          DO i = 13, 18

              IF (PrintStatus >= PrStatus_Normal)                       &
     &          WRITE (6,*) 'TDF: Overwriting ', FieldDesc(i)

! DEPENDS ON: writflds
              CALL WRITFLDS (TDF_unit,                                  &
                                                     ! in
     &                       NumLevs(i),                                &
                                                     ! in
     &                       FieldNum(i),                               &
                                                     ! in
     &                       A_LOOKUP,                                  &
                                                     ! in
     &                       LEN1_LOOKUP,                               &
                                                     ! in
     &                       D1_TDF(Addr_D1_TDF(i)),                    &
                                                     ! in
     &                       BufLen,                                    &
                                                     ! in
     &                       A_FIXHD,                                   &
                                                     ! in
#include "argppx.h"
     &                       ICode,                                     &
                                                     ! out
     &                       CMessage)               ! out

              IF (ICode > 0) THEN
                WRITE (CMessage, FMT='(A,A)') 'Error overwriting ',     &
     &                                        FieldDesc(i)
! DEPENDS ON: ereport
                CALL EReport (RoutineName, ICode, CMessage)
              END IF

          END DO

        END IF ! (L_TDF_ModifyCloud)

! Close TDF dump.

! DEPENDS ON: file_close
        CALL FILE_CLOSE(TDF_unit, 'TDF_dump', 8, 0, 0, ICode)
        IF (ICode > 0) THEN
          CMessage = 'Error closing TDF dump.'
! DEPENDS ON: ereport
          CALL EReport (RoutineName, ICode, CMessage)
        END IF

        IF (PrintStatus >= PrStatus_Normal) THEN
          WRITE (6,*) ''
          WRITE (6,*) 'TDF: TDF dump completed successfully.'
          WRITE (6,*) ''
        END IF

      END IF  ! L_LastCall


      RETURN
      END SUBROUTINE TDF


!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
!+ Remove problems and inconsistencies from cloud variables.

#endif
