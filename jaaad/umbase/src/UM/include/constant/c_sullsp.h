#if defined(A04_3B)
!--------------------COMDECK C_SULLSP-----------------------------
!
! LARGE SCAL PPN SCAVENGING COEFFS FOR SULPHUR CYCLE
!
      ! parameters for S Cycle scavenging
      ! Provisional values. May be set to 0 for no scavenging
      REAL,PARAMETER:: KRAIN_SO2    = 0.5E-4
      REAL,PARAMETER:: KSNOW_SO2    = 0.0
      REAL,PARAMETER:: KRAIN_NH3    = 0.5E-4
      REAL,PARAMETER:: KSNOW_NH3    = 0.0
      REAL,PARAMETER:: KRAIN_SO4AIT = 0.0 ! SO4_AITKEN
      REAL,PARAMETER:: KSNOW_SO4AIT = 0.0 ! SO4_AITKEN
      REAL,PARAMETER:: KRAIN_SO4ACC = 0.0 ! SO4_ACCU
      REAL,PARAMETER:: KSNOW_SO4ACC = 0.0 ! SO4_ACCU

      REAL,PARAMETER:: KRAIN_SO4DIS = 0.0 ! SO4_DISS
      ! N.B. Not same as SCONSCV value

      REAL,PARAMETER:: KSNOW_SO4DIS = 0.0 ! SO4_DISS

! C_SULCON end
#endif
