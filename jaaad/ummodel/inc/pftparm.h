! PFTPARM defines Surface parameters for each Plant Functional Type

      REAL:: ALBSNC_MAX(NPFT)  ! Snow-covered albedo for large LAI.
      REAL:: ALBSNC_MIN(NPFT)  ! Snow-covered albedo for zero LAI.
      REAL:: ALBSNF_MAX(NPFT)  ! Snow-free albedo for large LAI.
      REAL:: DZ0V_DH(NPFT)     ! Rate of change of vegetation
                               ! roughness length with height.
      REAL:: CATCH0(NPFT)      ! Minimum canopy capacity (kg/m2).
      REAL:: DCATCH_DLAI(NPFT) ! Rate of change of canopy capacity
                               ! with LAI.
      REAL:: INFIL_F(NPFT)     ! Infiltration enhancement factor.
      REAL:: KEXT(NPFT)        ! Light extinction coefficient.
      REAL:: ROOTD_FT(NPFT)    ! Rootdepth (m).
      !----------------------------------------------------------------
      !                     BT    NT   C3G   C4G    S
      !----------------------------------------------------------------
      COMMON  /RUN_PFT/ALBSNC_MAX,ALBSNC_MIN,ALBSNF_MAX,DZ0V_DH,        &
     &  CATCH0,DCATCH_DLAI,INFIL_F,KEXT,ROOTD_FT

! PFTPARM end
