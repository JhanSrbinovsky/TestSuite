#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
      ! Critical Richardson number, where Z0M_EFF=Z0M.
      ! Linear interpolation between RIB=0 and RI_CRIT
      REAL,PARAMETER:: RI_CRIT=0.5

      ! Tunable parameters in calculation of explicit orographic stresss
      REAL,PARAMETER:: ALPHA    = 12.0                                  &
                                             ! Tunable parameter for form
     &,                BETA     = 1.0                                   &
                                             ! drag calculation.
     &,                FD_DECAY = 0.3333                                &
                                             ! Height scale factors for
     &,                MAX_HT_SCALE  = 300.0                            &
                                             ! stress profiles
     &,                MIN_HT_SCALE  = 100.0
      LOGICAL,PARAMETER:: L_LOWHILL=.FALSE.  ! Set to .TRUE. for Wood &
!                                            ! Mason (1993) drag formula
!                                            ! Set to .FALSE. for steep
!                                            ! hill expression

!*----------------------------------------------------------------------
#endif
