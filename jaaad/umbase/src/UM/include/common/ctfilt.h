! ----------------------- header file: CTFilt  -----------------------
!
! Description:
!
!   Parameters and variables for Incremental Analysis Update (IAU) and
!   Temporal Digital Filtering (TDF) schemes.
!
!
! Current Code Owner: Adam Clayton.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! Declarations:

      INTEGER :: TDF_unit,                                              &
                             ! Unit number for TDF output dump.
     &           IAU_unit    ! Unit number for IAU increment.

      PARAMETER (TDF_unit = 107)
      PARAMETER (IAU_unit = 108)

      INTEGER :: LeaveAlone,                                            &
                             ! \  Options for filtering
     &           DirectFilt,                                            &
                             !  > of advected wind
     &           CopyFiltUVW ! /  components.

      PARAMETER (LeaveAlone  = 1)
      PARAMETER (DirectFilt  = 2)
      PARAMETER (CopyFiltUVW = 3)

      REAL    :: q_CC_tol,                                              &
                             ! qCL/qCF tolerance for cloud clearing.
     &           Weight_LL   ! Lower limit for absolute value of
                             ! filter weight.

      PARAMETER (q_CC_tol  = 1.0E-12  )
      PARAMETER (Weight_LL = 0.0000001)

      REAL    :: q_min       ! Minimum value allowed for q after
                             ! addition of q increments.

      REAL    :: oz_min      ! Minimum value allowed for ozone after
                             ! addition of ozone increments

      INTEGER :: MaxNumWeights ! Maximum number of filter weights.

      PARAMETER (MaxNumWeights = 1000)

      REAL    :: IAU_Weights(MaxNumWeights),                            &
                                             ! IAU filter weights array.
     &           TDF_Weights(MaxNumWeights)  ! TDF filter weights array.

      INTEGER :: IAU_FixHd(LEN_FIXHD) ! Fixed-length header of
                                      ! IAU increment dump.

      INTEGER :: D1_TFilt_len,                                          &
                                 ! D1 length for unpacked fields
     &           D1_IAU_k4_len,                                         &
                                 ! D1 length for packed IAU fields
     &           IAU_Len1Lookup,                                        &
                                 ! \ IAU increment
     &           IAU_Len2Lookup  ! / lookup dimensions

      LOGICAL :: L_Pack_D1_IAU   ! Hold IAU fields in packed form?

      INTEGER :: IAU_NumFldCodes
      PARAMETER (IAU_NumFldCodes = 13)

      ! Codes of fields that may be read in from the IAU file:
      INTEGER :: IAU_FldCodes(IAU_NumFldCodes)
      PARAMETER (IAU_FldCodes = (/ 2,   3,   4,   10,  12, 90,          &
     &                             150, 253, 254, 255, 407, 18001,480 /))

      ! Local lengths of fields (0 if field not required):
      INTEGER :: IAU_LocFldLens(IAU_NumFldCodes)

      ! Field descriptions:
      CHARACTER(7) :: IAU_FldDescs(IAU_NumFldCodes)
      PARAMETER  (IAU_FldDescs = (/ 'u      ', 'v      ', 'theta  ',    &
     &                              'q      ', 'qCF    ', 'aerosol',    &
     &                              'w      ', 'rho    ', 'qCL    ',    &
     &                              'exner  ', 'p      ', 'qT     ',     &
     &                              'ozone  ' /))

  ! IAU namelist variables:
  ! -----------------------

      LOGICAL L_IAU               ! Activate IAU scheme?

      INTEGER IAU_StartMin,                                             &
                                  ! Start minute of filtering period.
     &        IAU_EndMin,                                               &
                                  ! End   minute of filtering period.
     &        IAU_ApexMin         ! Apex minute for triangular filter.

      REAL    IAU_Cutoff_period,                                        &
                                  ! Filter cut-off period in hours.
     &        IAU_SBE_period      ! Stop band edge period in hours.

      LOGICAL L_IAU_CalcExnerIncs ! Use p increments to calculate
                                  ! exner increments?

      LOGICAL L_IAU_CalcThetaIncs ! Calculate theta increments using
                                  ! exner and q increments?

      LOGICAL L_IAU_CalcRhoIncs   ! Calculate rho increments using
                                  ! exner, theta and (possibly) q
                                  ! increments?

      LOGICAL L_IAU_IncTStar      ! If set, add level-one temperature
                                  ! increments to surface temperature
                                  ! and top-level soil temperature.

      LOGICAL L_IAU_ResetPoles    ! If set, reset polar rows of
                                  ! relevant increment fields to
                                  ! their mean values.

      LOGICAL L_IAU_RemoveSS      ! Remove supersaturation wrt water?

      LOGICAL L_IAU_CallStratQ    ! Reset stratospheric humidities at
                                  ! end of IAU insertion period?

      LOGICAL L_IAU_Diags         ! If set, write out IAU diagnostics.

      LOGICAL L_IAU_LowMem        ! If set, activate low memory (but
                                  ! high IO) version of the IAU code.

      LOGICAL L_IAU_DumpTS0State  ! If set, write out model state
                                  ! immediately after timestep-zero
                                  ! call to TFilt_cntl.

      LOGICAL L_IAU_CalcCloudIncs ! If set, calculate q, qcl & Cl
                                  ! from q or qT increments.

      LOGICAL L_IAU_IncrementIce  ! If this and L_IAU_CalcCloudIncs,
                                  ! calculate qcf & Cf

      LOGICAL L_IAU_ScaleCloud    ! If set, scale qcl, qcf, Cl & Cf
                                  ! increments to be in physical bounds
                                  
      LOGICAL L_IAU_UPPER_THETA   ! If set, then constrain the upper theta 
                                  ! increments.

      LOGICAL L_IAU_SetOzoneMin  ! If set, reset ozone to oz_min        
                                 ! in IAU if ozone was negative

      CHARACTER*10                                                      &
     &        IAU_FilterType      ! Filter type.

      REAL    IAU_LL_strat_pv,                                          &
                                  ! Lower-limit for strat ABS(PV).
     &        IAU_UL_strat_p,                                           &
                                  ! Upper-limit for strat pressure.
     &        IAU_LL_trop_p       ! Lower-limit for trop  pressure.

  ! TDF namelist variables:
  ! -----------------------

      LOGICAL L_TDF               ! Activate TDF scheme?

      INTEGER TDF_StartMin,                                             &
                                  ! Start minute of filtering period.
     &        TDF_EndMin,                                               &
                                  ! End   minute of filtering period.
     &        TDF_ApexMin         ! Apex minute for triangular filter.

      REAL    TDF_Cutoff_period,                                        &
                                  ! Filter cut-off period in hours.
     &        TDF_SBE_period      ! Stop band edge period in hours.

      LOGICAL L_TDF_FilterQ,                                            &
                                  ! Filter q?
     &        L_TDF_FilterQCL,                                          &
                                  ! Filter qCL?
     &        L_TDF_FilterQCF     ! Filter qCF?

      LOGICAL L_TDF_ModifyCloud   ! If set, modify cloud variables in
                                  ! TDF dump so that it becomes suitable
                                  ! for starting forecasts including
                                  ! physics.

      INTEGER TDF_AdvWindOpt      ! Filtering option for advected winds.

      LOGICAL L_TDF_CallStratQ    ! Reset stratospheric humidities in
                                  ! TDF dump?

      CHARACTER*10                                                      &
     &        TDF_FilterType      ! Filter type.


      COMMON / CTFilt /                                                 &
     & q_min, oz_min, IAU_Weights, TDF_Weights, IAU_FixHd,              &
     & D1_TFilt_len, D1_IAU_k4_len, IAU_Len1Lookup, IAU_Len2Lookup,     &
     & L_Pack_D1_IAU, IAU_LocFldLens,                                   &
     & L_IAU, IAU_StartMin, IAU_EndMin,                                 &
     & IAU_ApexMin, IAU_Cutoff_period, IAU_SBE_period,                  &
     & L_IAU_CalcCloudIncs, L_IAU_IncrementIce, L_IAU_ScaleCloud,       &
     & L_IAU_CalcExnerIncs, L_IAU_CalcThetaIncs, L_IAU_CalcRhoIncs,     &
     & L_IAU_IncTStar, L_IAU_ResetPoles, L_IAU_RemoveSS,                &
     & L_IAU_CallStratQ, L_IAU_Diags, L_IAU_LowMem, L_IAU_DumpTS0State, &
     & L_IAU_UPPER_THETA,                                               &
     & IAU_LL_strat_pv, IAU_UL_strat_p, IAU_LL_trop_p,                  &
     & L_TDF, TDF_StartMin, TDF_EndMin,                                 &
     & TDF_ApexMin, TDF_Cutoff_period, TDF_SBE_period,                  &
     & L_TDF_FilterQ, L_TDF_FilterQCL, L_TDF_FilterQCF,                 &
     & L_TDF_ModifyCloud, TDF_AdvWindOpt,                               &
     & L_TDF_CallStratQ, L_IAU_SetOzoneMin,                             &
     ! Character variables at the end of common block
     & IAU_FilterType, TDF_FilterType

      NAMELIST / RUN_TFilt /                                            &
     & L_IAU, IAU_StartMin, IAU_EndMin,                                 &
     & IAU_ApexMin, IAU_Cutoff_period, IAU_SBE_period,                  &
     & L_IAU_CalcCloudIncs, L_IAU_IncrementIce, L_IAU_ScaleCloud,       &
     & L_IAU_CalcExnerIncs, L_IAU_CalcThetaIncs, L_IAU_CalcRhoIncs,     &
     & L_IAU_IncTStar, L_IAU_ResetPoles, L_IAU_RemoveSS,                &
     & L_IAU_CallStratQ, L_IAU_Diags, L_IAU_LowMem, L_IAU_DumpTS0State, &
     & L_IAU_UPPER_THETA,                                               &
     & IAU_LL_strat_pv, IAU_UL_strat_p, IAU_LL_trop_p,                  &
     & L_TDF, TDF_StartMin, TDF_EndMin,                                 &
     & TDF_ApexMin, TDF_Cutoff_period, TDF_SBE_period,                  &
     & L_TDF_FilterQ, L_TDF_FilterQCL, L_TDF_FilterQCF,                 &
     & L_TDF_ModifyCloud, TDF_AdvWindOpt,                               &
     & L_TDF_CallStratQ, L_IAU_SetOzoneMin,                             &
     & IAU_FilterType, TDF_FilterType
