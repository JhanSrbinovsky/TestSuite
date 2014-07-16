#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Combined control routine for IAU and TDF schemes.

      SUBROUTINE TFilt_cntl (                                           &
      U, V, W, U_adv, V_adv, W_adv,                                     &
      Theta, Exner_rho_levels, Rho,                                     &
      Q, Qcl, Qcf, Murk, OZONE_TRACER,                                  &
      Deep_soil_temp,                                                   &
      P, P_theta_levels, Exner_theta_levels,                            &
      snodep,                                                           &
      cf_area, cf_bulk, cf_liquid, cf_frozen,                           &
      Pstar, Tstar, Tstar_tile,                                         &
#include "argduma.h"
#include "argcona.h"
#include "arglndm.h"
#include "argppx.h"
     &                        l_mixing_ratio,                           &
                                              ! In. Use mixing ratios
     &                        IAU_lookup,                               &
                                          ! inout
     &                        D1_IAU_k4,                                &
                                          ! inout
     &                        D1_TFilt )  ! inout

! Description:
!
!   Combined control routine for Incremental Analysis Update (IAU) and
!   Temporal Digital Filtering (TDF) schemes.
!
!   1. Check to see if this is the first call, and if so carry out some
!      basic checks and call a routine to calculate the filter weights.
!
!   2. If required, call the relevant filtering routine.
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
#include "cmaxsize.h"
#include "typduma.h"
#include "typcona.h"
#include "typlndm.h"
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "clookadd.h"
#include "chsunits.h"
#include "cprintst.h"
#include "cntlall.h"
#include "ctfilt.h"
#include "ctime.h"
#include "cntlatm.h"
#include "c_kinds.h"

! Subroutine arguments:

      Real, Intent (InOut) ::                                           &
     &  u(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)       &
     &, v(1-offx:row_length+offx, 1-offy:n_rows+offy, model_levels)     &
     &, w(1-offx:row_length+offx, 1-offy:rows+offy, 0:model_levels)     &
     &, u_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &           model_levels)                                          &
     &, v_adv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,       &
     &           model_levels)                                          &
     &, w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &           0:model_levels)                                        &
     &, theta(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)   &
     &, exner_rho_levels(1-offx:row_length+offx, 1-offy:rows+offy,      &
     &            model_levels+1)                                       &
     &, rho(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)     &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &           wet_levels)                                            &
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &           wet_levels)                                            &
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &           wet_levels)                                            &
     &, murk(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)    &
! Add cariolle specific parameters for ozone tracer  
     &, OZONE_TRACER(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &              model_levels)                            



      Real, Intent (InOut) ::                                           &
     &  p(1-offx:row_length+offx, 1-offy:rows+offy,                     &
     &      model_levels+1)                                             &
     &, p_theta_levels(1-offx:row_length+offx,                          &
     &                   1-offy:rows+offy, model_levels)                &
     &, exner_theta_levels(1-offx:row_length+offx,                      &
     &                     1-offy:rows+offy, model_levels)              &
     &, snodep(theta_field_size)                                        &
     &, cf_area(row_length, rows, wet_levels)                           &
     &, cf_bulk(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &           wet_levels)                                            &
     &, cf_liquid(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,     &
     &             wet_levels)                                          &
     &, cf_frozen(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,     &
     &             wet_levels)                                          &
     &, pstar(row_length,rows)                                          &
     &, Tstar(row_length*rows)                                          &
     &, Tstar_tile(land_points,ntiles)                                  &
     &, deep_soil_temp(land_points,sm_levels)

      INTEGER,      INTENT(INOUT) ::                                    &
                                     ! Lookup tables of IAU inc file

     &  IAU_lookup(IAU_Len1Lookup, IAU_Len2Lookup)

      REAL(KIND=real32), INTENT(INOUT) ::                               &
                                          ! Array for packed IAU fields

     &  D1_IAU_k4(D1_IAU_k4_len)

      REAL,         INTENT(INOUT) ::                                    &
                                     ! Array for unpacked IAU/TDF fields

     &  D1_TFilt(D1_TFilt_len)

      Logical, intent(in) ::                                            &
     &  l_mixing_ratio               ! Use mixing ratios (if available)

! Local variables:

      INTEGER :: i,                                                     &
     &           ICode,                                                 &
     &           TS_len_secs,                                           &
     &           NumWeights,                                            &
     &           WeightNum,                                             &
     &           ElapsedDays,                                           &
     &           ElapsedSecs,                                           &
     &           IAU_VTSecs,                                            &
     &           IAU_VTMins,                                            &
     &           Lookup_len,                                            &
     &           LenIO,                                                 &
     &           Sec,                                                   &
     &           Min,                                                   &
     &           ss,                                                    &
     &           MiddleSec,                                             &
     &           CallNumber

      REAL ::    A_IO,                                                  &
     &           Weight,                                                &
     &           Timestep

      LOGICAL :: L_CallIAU,                                             &
     &           L_CallTDF,                                             &
     &           L_MiddleCall,                                          &
     &           L_LastCall

      CHARACTER(6)  :: Min_str
      CHARACTER(80) :: CMessage

      CHARACTER(*)  :: RoutineName
      PARAMETER      ( RoutineName='TFilt_cntl' )

!- End of header ------------------------------------------------------

      TS_len_secs = INT(SECS_PER_STEPIM(A_IM))
      Timestep = secs_per_stepim(atmos_im)

!----------------------------------------------------------------------
! [1]: On first call, check inputs, calculate weights, etc.
!----------------------------------------------------------------------

      IF (STEPIM(A_IM) == 0) THEN

!----------------------------------------------------------------------
! [1.1]: IAU scheme.
!----------------------------------------------------------------------

        IF (L_IAU) THEN

! Calculate filter Weights:

! DEPENDS ON: calc_tfiltwts
          CALL Calc_TFiltWts ( IAU_StartMin*60,                         &
                                                  ! in
     &                         IAU_EndMin*60,                           &
                                                  ! in
     &                         TS_len_secs,                             &
                                                  ! in
     &                         IAU_ApexMin*60,                          &
                                                  ! in
     &                         IAU_Cutoff_period,                       &
                                                  ! in
     &                         IAU_SBE_period,                          &
                                                  ! in
     &                         IAU_FilterType,                          &
                                                  ! in
     &                         MaxNumWeights,                           &
                                                  ! in
     &                         IAU_Weights,                             &
                                                  ! out
     &                         NumWeights )       ! out

! Read in lookup tables from IAU increment file:

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
     &                  A_IO, LenIO, Lookup_len)
            ICode    = NINT(A_IO) + 1
            CMessage = 'Error reading IAU increment lookup tables.'
! DEPENDS ON: ereport
            CALL EReport (RoutineName, ICode, CMessage)
          END IF

! Basic checks:

          IF (IAU_FixHd(4) /= A_FIXHD(4)) THEN
            ICode = 1
            CMessage = 'Increment has wrong grid type.'
! DEPENDS ON: ereport
            CALL EReport (RoutineName, ICode, CMessage)
          END IF

! DEPENDS ON: time2sec
          CALL TIME2SEC (                                               &
     &      IAU_lookup(LBYR,1),  IAU_lookup(LBMON,1),                   &
                                                      ! in
     &      IAU_lookup(LBDAT,1), IAU_lookup(LBHR,1),                    &
                                                      ! in
     &      IAU_lookup(LBMIN,1), 0,                                     &
                                                      ! in
     &      BASIS_TIME_DAYS,     BASIS_TIME_SECS,                       &
                                                      ! in
     &      ElapsedDays,         ElapsedSecs,                           &
                                                      ! out
     &      LCAL360)                                  ! in

          IAU_VTSecs = ElapsedSecs + ElapsedDays * 86400
          IAU_VTMins = IAU_VTSecs / 60

          IF ( (IAU_VTMins < IAU_StartMin) .OR.                         &
     &         (IAU_VTMins > IAU_EndMin) ) THEN
            ICode = 1
            CMessage = 'Increment validity time outside IAU period.'
! DEPENDS ON: ereport
            CALL EReport (RoutineName, ICode, CMessage)
          END IF

! Print details:

          IF (PrintStatus >= PrStatus_Normal) THEN

            WRITE (6,*) ''
            WRITE (6,*) 'TFilt_cntl: Incremental Analysis Update (IAU)' &
     &                           //' details.'
            WRITE (6,*) '            <><><><><><><><><><><><><><><><><' &
     &                           //'><><><><>'
            WRITE (6,*) ''
            WRITE (6,*) '            Start/end minutes:       ',        &
     &                               IAU_StartMin, '/', IAU_EndMin

            WRITE (6,*) '            Calculate cloud incs?:   ',        &
     &        L_IAU_CalcCloudIncs

            WRITE (6,*) '            Calculate qcf  incs?:   ',        &
     &        L_IAU_IncrementIce

            WRITE (6,*) '            Scale cloud incs?:       ',        &
     &        L_IAU_ScaleCloud

            WRITE (6,*) '            Calculate exner incs?:   ',        &
     &        L_IAU_CalcExnerIncs

            WRITE (6,*) '            Calculate theta incs?:   ',        &
     &        L_IAU_CalcThetaIncs

            WRITE (6,*) '            Limit upper level theta incs ?:  ',&
     &        L_IAU_UPPER_THETA

            WRITE (6,*) '            Calculate rho   incs?:   ',        &
     &        L_IAU_CalcRhoIncs

            WRITE (6,*) '            Call StratQ?:            ',        &
     &        L_IAU_CallStratQ

            WRITE (6,*) '            Remove supersaturation?: ',        &
     &        L_IAU_RemoveSS

            WRITE (6,*) '            Increment TStar/TSoil?:  ',        &
     &        L_IAU_IncTStar

            IF (Model_domain == mt_global) THEN
              WRITE (6,*) '            Reset polar rows?:       ',      &
     &          L_IAU_ResetPoles
            END IF

            WRITE (6,*) '            Low memory version?:     ',        &
     &        L_IAU_LowMem

            IF (.NOT.L_IAU_LowMem) THEN
              WRITE (6,*) '            Pack IAU D1 array?:      ',      &
     &          L_Pack_D1_IAU
            END IF

            WRITE (6,*) '            Filter type:             ',        &
     &                               IAU_FilterType

            IF (IAU_FilterType == 'Triangular')                         &
     &        WRITE (6,*) '            Apex minute:             ',      &
     &                                 IAU_ApexMin

            IF (IAU_FilterType == 'Unwindowed' .OR.                     &
     &          IAU_FilterType == 'LancWin   ')                         &
     &        WRITE (6,*) '            Cut-off period:          ',      &
     &                                 IAU_Cutoff_period, ' hrs'

            IF (IAU_FilterType == 'Dolph     ')                         &
     &        WRITE (6,*) '            Stop band edge period:   ',      &
     &                                 IAU_SBE_period, ' hrs'

            WRITE (6,*) ''
            WRITE (6,*) '            Filter Weights:'
            WRITE (6,*) ''
            WRITE (6,*) '              Min:Sec  Weight'
            WRITE (6,*) '              -------  ------'

            DO WeightNum = 1, NumWeights
              Sec = (WeightNum - 1) * TS_len_secs + IAU_StartMin * 60
              Min = Sec / 60
              ss  = ABS(Sec - Min * 60)
              WRITE (6,'(A,I6,A,I2.2,A,F13.10)') '            ',        &
     &          Min, ':', ss, '  ', IAU_Weights(WeightNum)
            END DO

            WRITE (6,*) ''

          END IF ! (PrintStatus >= PrStatus_Normal)

      END IF ! (L_IAU)

!----------------------------------------------------------------------
! [1.2]: TDF scheme.
!----------------------------------------------------------------------

        IF (L_TDF) THEN

! Calculate filter Weights:

! DEPENDS ON: calc_tfiltwts
          CALL Calc_TFiltWts ( TDF_StartMin*60,                         &
                                                  ! in
     &                         TDF_EndMin*60,                           &
                                                  ! in
     &                         TS_len_secs,                             &
                                                  ! in
     &                         TDF_ApexMin*60,                          &
                                                  ! in
     &                         TDF_Cutoff_period,                       &
                                                  ! in
     &                         TDF_SBE_period,                          &
                                                  ! in
     &                         TDF_FilterType,                          &
                                                  ! in
     &                         MaxNumWeights,                           &
                                                  ! in
     &                         TDF_Weights,                             &
                                                  ! out
     &                         NumWeights )       ! out

! Basic checks:

          i = (NumWeights/2) * 2

          IF (i == NumWeights) THEN
            ICode = 1
            CMessage = 'TDF filtering period must be centred on a'      &
     &               //' timestep.'
! DEPENDS ON: ereport
            CALL EReport (RoutineName, ICode, CMessage)
          END IF

! Print details:

          IF (PrintStatus >= PrStatus_Normal) THEN

            WRITE (6,*) ''
            WRITE (6,*) 'TFilt_cntl: Temporal Digital Filtering (TDF)'  &
     &                           //' details.'
            WRITE (6,*) '            <><><><><><><><><><><><><><><><>'  &
     &                           //'<><><><><'
            WRITE (6,*) ''
            WRITE (6,*) '            Start/end minutes:     ',          &
     &                               TDF_StartMin, '/', TDF_EndMin

            IF (L_TDF_FilterQ) THEN
              WRITE (6,*) '            Filtering q?:          Yes'
            ELSE
              WRITE (6,*) '            Filtering q?:          No'
            END IF

            IF (L_TDF_FilterQCL) THEN
              WRITE (6,*) '            Filtering qCL?:        Yes'
            ELSE
              WRITE (6,*) '            Filtering qCL?:        No'
            END IF

            IF (L_TDF_FilterQCF) THEN
              WRITE (6,*) '            Filtering qCF?:        Yes'
            ELSE
              WRITE (6,*) '            Filtering qCF?:        No'
            END IF

            IF (L_TDF_ModifyCloud) THEN
              WRITE (6,*) '            Modifying cloud?:      Yes'
            ELSE
              WRITE (6,*) '            Modifying cloud?:      No'
            END IF

            IF (TDF_AdvWindOpt == LeaveAlone)                           &
     &        WRITE (6,*) '            Advected winds option: '         &
     &                          //'Leave alone'

            IF (TDF_AdvWindOpt == DirectFilt)                           &
     &        WRITE (6,*) '            Advected winds option: '         &
     &                          //'Filter directly'

            IF (TDF_AdvWindOpt == CopyFiltUVW)                          &
     &        WRITE (6,*) '            Advected winds option: '         &
     &                          //'Set to filtered non-advected winds'

            WRITE (6,*) '            Filter type:           ',          &
     &                               TDF_FilterType

            IF (TDF_FilterType == 'Triangular')                         &
     &        WRITE (6,*) '            Apex minute:           ',        &
     &                                 TDF_ApexMin

            IF (TDF_FilterType == 'Unwindowed' .OR.                     &
     &          TDF_FilterType == 'LancWin   ')                         &
     &        WRITE (6,*) '            Cut-off period:        ',        &
     &                                 TDF_Cutoff_period, ' hrs'

            IF (TDF_FilterType == 'Dolph     ')                         &
     &        WRITE (6,*) '            Stop band edge period: ',        &
     &                                 TDF_SBE_period, ' hrs'

            WRITE (6,*) ''
            WRITE (6,*) '            Filter Weights:'
            WRITE (6,*) ''
            WRITE (6,*) '              Min:Sec  Weight'
            WRITE (6,*) '              -------  ------'

            DO WeightNum = 1, NumWeights
              Sec = (WeightNum - 1) * TS_len_secs + TDF_StartMin * 60
              Min = Sec / 60
              ss  = ABS(Sec - Min * 60)
              WRITE (6,'(A,I6,A,I2.2,A,F13.10)') '            ',        &
     &          Min, ':', ss, '  ', TDF_Weights(WeightNum)
            END DO

            WRITE (6,*) ''

          END IF ! (PrintStatus >= PrStatus_Normal)

        END IF ! (L_TDF)

      END IF ! (STEPIM(A_IM) == 0)

!----------------------------------------------------------------------
! [2]: If required, call main IAU or TDF routine.
!----------------------------------------------------------------------

      Sec =  STEPIM(A_IM) * TS_len_secs

      L_CallIAU = L_IAU .AND. (Sec >= IAU_StartMin * 60)                &
     &                  .AND. (Sec <= IAU_EndMin   * 60)

      ! Cater for runs restarted during the IAU insertion period:
      IF (IAU_StartMin < 0 .AND. Sec == 0) L_CallIAU = .FALSE.

      L_CallTDF = L_TDF .AND. (Sec >= TDF_StartMin * 60)                &
     &                  .AND. (Sec <= TDF_EndMin   * 60)

      IF (L_CallIAU) THEN

        L_LastCall  = (Sec == IAU_EndMin * 60)

        CallNumber  = 1 + (Sec - IAU_StartMin * 60) / TS_len_secs
        Weight      = IAU_Weights(CallNumber)

        Min = Sec / 60
        WRITE (Min_str,*) Min
        ss  = ABS(Sec - Min * 60)

        IF (PrintStatus >= PrStatus_Normal) THEN
          WRITE (6,*) ''
          WRITE (6,*) 'TFilt_cntl: Calling main IAU routine.'
          WRITE (6,*) '           '
          WRITE (6,'(A,F13.10)')                                        &
     &               '               Weight:  ', Weight
          WRITE (6,'(A,A,A,I2.2)')                                      &
     &               '               Min:Sec: ', TRIM(Min_str), ':', ss
          WRITE (6,*) ''
        END IF

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER ('IAU     ', 3)

! DEPENDS ON: iau
        CALL IAU (                                                      &
#include "argppx.h"
#include "argcona.h"
#include "arglndm.h"
     &             Weight, Timestep,                                    &
                                                            ! in
     &             L_LastCall,                                          &
                                                            ! in
     &             l_mixing_ratio,                                      &
                                                            ! in
! IAU work arrays:
     &             IAU_lookup,                                          &
                                                            ! inout
     &             D1_IAU_k4,                                           &
                                                            ! inout
     &             D1_TFilt,                                            &
                                                            ! inout
! Model fields:   
     &             U,   V,   W,   U_ADV,   V_ADV,   W_ADV,              &
                                                            ! inout
     &             THETA, EXNER_RHO_LEVELS, RHO,                        &
                                                            ! inout
     &             Q,   QCL,   QCF,                                     &
                                                            ! inout
     &             MURK,   DEEP_SOIL_TEMP,   TSTAR,   TSTAR_TILE,       &
                                                            ! inout
     &             P,   PSTAR,   P_THETA_LEVELS,   EXNER_THETA_LEVELS,  &
                                                            ! inout
     &             SNODEP,                                              &
                                                             !inout
     &             CF_AREA,   CF_BULK,   CF_LIQUID, CF_FROZEN,          &
                                                            !inout
     &             OZONE_TRACER)

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER ('IAU     ', 4)

      END IF ! (L_CallIAU)

      IF (L_CallTDF) THEN

        L_LastCall   = (Sec == TDF_EndMin * 60)

        MiddleSec    = (TDF_StartMin + TDF_EndMin) * 60 / 2
        L_MiddleCall = (Sec == MiddleSec)

        CallNumber   = 1 + (Sec - TDF_StartMin * 60) / TS_len_secs
        Weight       = TDF_Weights(CallNumber)

        Min = Sec / 60
        WRITE (Min_str,*) Min
        ss  = ABS(Sec - Min * 60)

        IF (PrintStatus >= PrStatus_Normal) THEN
          WRITE (6,*) ''
          WRITE (6,*) 'TFilt_cntl: Calling main TDF routine.'
          WRITE (6,*) '           '
          WRITE (6,'(A,F13.10)')                                        &
     &               '               Weight:  ', Weight
          WRITE (6,'(A,A,A,I2.2)')                                      &
     &               '               Min:Sec: ', TRIM(Min_str), ':', ss
          WRITE (6,*) ''
        END IF

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER ('TDF     ', 3)

! DEPENDS ON: tdf
        CALL TDF (                                                      &
     &             U,   V,   W,   Rho,   Theta,   Exner_rho_levels,     &
                                  ! inout
     &             Q,   Qcl,   Qcf,   U_adv,   V_adv,   W_adv,          &
                                  ! inout
#include "argduma.h"
#include "argppx.h"
     &             Weight,                                              &
                                 ! in
     &             L_MiddleCall,                                        &
                                 ! in
     &             L_LastCall,                                          &
                                 ! in
     &             D1_TFilt )    ! inout

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER ('TDF     ', 4)

      END IF ! (L_CallTDF)


      RETURN
      END SUBROUTINE TFilt_cntl
#endif
