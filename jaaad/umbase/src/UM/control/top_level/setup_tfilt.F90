#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Temporal filtering set-up routine.

      SUBROUTINE Setup_TFilt

! Description:
!
!   Temporal filtering set-up routine.
!
!   1. Obtain values for namelist variables.
!
!   2. Calculate work array dimensions.
!
!
! Current Code Owner: Adam Clayton.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! Declarations:

      IMPLICIT NONE

! Common blocks:

#include "parvars.h"
#include "typsize.h"
#include "nstypes.h"
#include "ctfilt.h"
#include "clookadd.h"
#include "chsunits.h"
#include "cntlall.h"
#include "cmaxsize.h"
#include "cruntimc.h"

! Local variables:

      INTEGER :: ICode,                                                 &
     &           D1_IAU_len,                                            &
                                       ! Space required for IAU fields.
     &           Code,                                                  &
     &           LenIO,                                                 &
     &           FieldNum,                                              &
     &           Lookup_len,                                            &
     &           i,                                                     &
     &           LocFldLen,                                             &
     &           pack_code

      REAL :: A_IO

      INTEGER, ALLOCATABLE :: IAU_lookup(:,:)

      CHARACTER(80) :: CMessage

      CHARACTER(*)  :: RoutineName
      PARAMETER      ( RoutineName='Setup_TFilt' )

! External subroutines called:

      EXTERNAL EReport, FILE_OPEN, READ_FLH, SETPOS, BUFFIN, IOERROR

!- End of header ------------------------------------------------------

!----------------------------------------------------------------------
! [1]: Set up q_min.
!----------------------------------------------------------------------

      IF (L_QPOS) THEN
        q_min = qlimit ! qlimit specified via UMUI
      ELSE
        q_min = 1.0E-8
      END IF

! Add parameter for ozone
      oz_min = 1.0E-8 

!----------------------------------------------------------------------
! [2]: Obtain values for namelist variables.
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! [2.1]: Set up defaults.
!----------------------------------------------------------------------

      L_IAU = .FALSE.
      L_TDF = .FALSE.

      IAU_StartMin = 0
      TDF_StartMin = 0

      IAU_EndMin = 0
      TDF_EndMin = 0

      IAU_FilterType = 'Uniform   '
      TDF_FilterType = 'Uniform   '

      IAU_ApexMin = 0
      TDF_ApexMin = 0

      IAU_Cutoff_period = 0.0
      TDF_Cutoff_period = 0.0

      IAU_SBE_period = 0.0
      TDF_SBE_period = 0.0

      L_IAU_RemoveSS = .TRUE.

      L_IAU_CallStratQ = .TRUE.
      L_TDF_CallStratQ = .FALSE.

      IAU_LL_strat_pv = 2.5E-06
      IAU_UL_strat_p  = 4.0E+04
      IAU_LL_trop_p   = 1.0E+04

      L_IAU_Diags = .FALSE.

      L_IAU_LowMem = .FALSE.

      L_IAU_CalcCloudIncs = .FALSE.
      L_IAU_IncrementIce  = .FALSE.
      L_IAU_ScaleCloud    = .FALSE.
      L_IAU_CalcExnerIncs = .FALSE.
      L_IAU_CalcThetaIncs = .FALSE.
      L_IAU_CalcRhoIncs   = .FALSE.

      L_IAU_IncTStar = .FALSE.

      L_IAU_ResetPoles = .FALSE.

      L_IAU_DumpTS0State = .FALSE.
      L_IAU_UPPER_THETA  = .FALSE.
      L_TDF_FilterQ   = .FALSE.
      L_TDF_FilterQCL = .FALSE.
      L_TDF_FilterQCF = .FALSE.

      L_TDF_ModifyCloud = .FALSE.

      TDF_AdvWindOpt = LeaveAlone
      L_IAU_SetOzoneMin = .FALSE.
!----------------------------------------------------------------------
! [2.2]: Read in namelist (which should already be open on unit 5).
!----------------------------------------------------------------------

      REWIND (UNIT=5)

      READ (UNIT=5, NML=RUN_TFilt, IOSTAT=ICode)
      IF (ICode > 0) THEN
        CMessage='Error reading RUN_TFilt namelist.'
! DEPENDS ON: ereport
        CALL EReport (RoutineName, ICode, CMessage)
      END IF

      REWIND (UNIT=5)

!----------------------------------------------------------------------
! [2.3]: Check for errors/inconsistencies etc.
!----------------------------------------------------------------------

      ! Turn off IAU if in FASTRUN mode:
      IF (RUN_ASSIM_MODE  ==  "FastRun   ") THEN
        IF (L_IAU) THEN
          RUN_ASSIM_MODE = "NoIAU     "
          L_IAU = .FALSE.
        ELSE
          RUN_ASSIM_MODE = "None      "
        END IF
      END IF

      ! No need to modify cloud variables if running with physics and
      ! not filtering qCL or qCF:
      IF ( L_Physics .AND. .NOT.L_TDF_FilterQCL                         &
     &               .AND. .NOT.L_TDF_FilterQCF )                       &
     &  L_TDF_ModifyCloud = .FALSE.

      ! If there is just one IAU update time, might as well use the low
      ! memory option:
      IF (IAU_StartMin == IAU_EndMin) L_IAU_LowMem = .TRUE.

      ! Check for overlapping filtering periods:
      IF (L_IAU .AND. L_TDF) THEN

        IF (TDF_StartMin  <   IAU_EndMin) THEN
          IF (TDF_EndMin  >=  IAU_StartMin) THEN

            ICode    = 1
            CMessage = 'IAU/TDF filtering periods must not overlap.'
! DEPENDS ON: ereport
            CALL EReport (RoutineName, ICode, CMessage)

          END IF
        END IF

      END IF

!----------------------------------------------------------------------
! [3]: Calculate work array dimensions.
!----------------------------------------------------------------------

      ! Defaults:
      D1_TFilt_len   = 1
      D1_IAU_k4_len  = 1
      IAU_Len1Lookup = 1
      IAU_Len2Lookup = 1
      L_Pack_D1_IAU  = .FALSE.

      IF (L_TDF) THEN

        D1_TFilt_len =   u_off_size     *  model_levels                 &
                                                           ! u
     &                 + v_off_size     *  model_levels                 &
                                                           ! v
     &                 + theta_off_size * (model_levels+1)              &
                                                           ! w
     &                 + theta_off_size *  model_levels                 &
                                                           ! rho
     &                 + theta_off_size *  model_levels                 &
                                                           ! theta
     &                 + theta_off_size * (model_levels+1) ! exner

        IF (L_TDF_FilterQ)                                              &
     &    D1_TFilt_len = D1_TFilt_len + theta_halo_size * wet_levels

        IF (L_TDF_FilterQCL)                                            &
     &    D1_TFilt_len = D1_TFilt_len + theta_halo_size * wet_levels

        IF (L_TDF_FilterQCF)                                            &
     &    D1_TFilt_len = D1_TFilt_len + theta_halo_size * wet_levels

        IF (TDF_AdvWindOpt == DirectFilt)                               &
     &    D1_TFilt_len = D1_TFilt_len                                   &
     &                   + u_halo_size     *  model_levels              &
     &                   + v_halo_size     *  model_levels              &
     &                   + theta_halo_size * (model_levels+1)

      END IF

      IF (L_IAU) THEN

        ! Get local field lengths for IAU fields:
        DO i = 1, IAU_NumFldCodes

          Code = IAU_FldCodes(i)

          ! Default:
          IAU_LocFldLens(i) = theta_field_size

          ! Exceptions:
          IF (Code == 2) IAU_LocFldLens(i) = u_field_size ! u
          IF (Code == 3) IAU_LocFldLens(i) = v_field_size ! v

          IF (Code == 4 .AND. L_IAU_CalcThetaIncs)                      &
     &      IAU_LocFldLens(i) = 0 ! Don't read theta

          IF (Code == 253 .AND. L_IAU_CalcRhoIncs)                      &
     &      IAU_LocFldLens(i) = 0 ! Don't read rho

          IF (Code == 255 .AND. L_IAU_CalcExnerIncs)                    &
     &      IAU_LocFldLens(i) = 0 ! Don't read exner

          IF (Code == 407 .AND. .NOT.L_IAU_CalcExnerIncs)               &
     &      IAU_LocFldLens(i) = 0 ! Don't read p

        END DO

        ! Open IAU increment file.
        ! (It is left open until the end of the IAU insertion period.)
! DEPENDS ON: file_open
        CALL FILE_OPEN (IAU_unit, 'IAU_inc', 7, 0, 0, ICode)
        IF (ICode > 0) THEN
          CMessage = 'Error opening IAU increment file.'
! DEPENDS ON: ereport
          CALL EReport (RoutineName, ICode, CMessage)
        END IF

        ! Read in fixed-length header:
! DEPENDS ON: read_flh
        CALL READ_FLH (IAU_unit, IAU_FixHd, LEN_FIXHD,                  &
     &                 ICode,    CMessage)
        IF (ICode > 0) THEN
! DEPENDS ON: ereport
          CALL EReport (RoutineName, ICode, CMessage)
        END IF

        IAU_Len1Lookup = IAU_FixHd(151)
        IAU_Len2Lookup = IAU_FixHd(152)

        IF (IAU_Len2Lookup == 0) THEN

          ICode    = 1
          CMessage = 'No fields in IAU increment file.'
! DEPENDS ON: ereport
          CALL EReport (RoutineName, ICode, CMessage)

        END IF

        ! Calculate required length of work array:
        IF (.NOT.L_IAU_LowMem) THEN

          ! Read in lookup tables from IAU increment file:
          ALLOCATE (IAU_lookup(IAU_Len1Lookup, IAU_Len2Lookup))

! DEPENDS ON: setpos
          CALL SETPOS (IAU_unit, IAU_FixHd(150)-1, ICode)
          IF (ICode > 0) THEN
            CMessage = 'SETPOS error moving to start of lookup'         &
     &               //' tables in IAU increment file.'
! DEPENDS ON: ereport
            CALL EReport (RoutineName, ICode, CMessage)
          END IF

          Lookup_len = IAU_Len1Lookup * IAU_Len2Lookup

! DEPENDS ON: buffin
          CALL BUFFIN (IAU_unit, IAU_lookup(1,1), Lookup_len,           &
     &                 LenIO,    A_IO)
          IF ( (A_IO  /= -1.0      ) .OR.                               &
     &         (LenIO /= Lookup_len) ) THEN
! DEPENDS ON: ioerror
            CALL IOERROR ('Buffer in of IAU inc lookups',               &
     &                    A_IO, LenIO, Lookup_len)
            ICode    = NINT(A_IO) + 1
            CMessage = 'Error reading IAU increment lookup tables.'
! DEPENDS ON: ereport
            CALL EReport (RoutineName, ICode, CMessage)
          END IF

          D1_IAU_len = 0
          L_Pack_D1_IAU = .TRUE.

          ! Loop over fields in IAU increment file:
          DO FieldNum = 1, IAU_Len2Lookup

            Code = IAU_lookup(ITEM_CODE, FieldNum)

            ! Get local field length:
            LocFldLen = 0
            DO i = 1, IAU_NumFldCodes
              IF (IAU_FldCodes(i) == Code)                              &
     &          LocFldLen = IAU_LocFldLens(i)
            END DO

            D1_IAU_len = D1_IAU_len + LocFldLen

            ! If any of the increment fields are not 32-bit packed, do
            ! not pack IAU D1 array:
            IF (LocFldLen > 0) THEN
              pack_code = MOD(IAU_lookup(LBPACK, FieldNum), 10)
              IF (pack_code /= 2) L_Pack_D1_IAU = .FALSE.
            END IF

          END DO

          DEALLOCATE (IAU_lookup)

          ! Do not use packed D1 array if also running TDF scheme:
          IF (L_TDF) L_Pack_D1_IAU = .FALSE.

          IF (L_Pack_D1_IAU) THEN
            D1_IAU_k4_len = D1_IAU_len
          ELSE
            D1_TFilt_len = MAX(D1_TFilt_len, D1_IAU_len)
          END IF

        END IF ! (.NOT.L_IAU_LowMem)

      END IF ! (L_IAU)


      RETURN
      END SUBROUTINE Setup_TFilt
#endif
