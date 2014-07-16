! NVEGPARM start

! Surface and vegetation parameters
      REAL :: ALBSNC_NVG(NNVG)  ! Snow-covered albedo.
      REAL :: ALBSNF_NVG(NNVG)  ! Snow-free albedo.
      REAL :: CATCH_NVG(NNVG)   ! Canopy capacity (kg/m2).
      REAL :: GS_NVG(NNVG)          ! Surface conductance (m/s).
      REAL :: INFIL_NVG(NNVG)       ! Infiltration enhancement factor.
      REAL :: ROOTD_NVG(NNVG)   ! Rootdepth (m).
      REAL :: Z0_NVG(NNVG)          ! Roughness length (m).
      REAL :: CH_NVG(NNVG)         ! "Canopy" heat capacity (J/K/m2)
      REAL :: VF_NVG(NNVG)         ! Fractional "canopy" coverage

      COMMON  /RUN_BLVEG/ALBSNC_NVG,ALBSNF_NVG,CATCH_NVG,GS_NVG,        &
     &  INFIL_NVG,ROOTD_NVG,Z0_NVG,CH_NVG,VF_NVG

! NVEGPARM end
