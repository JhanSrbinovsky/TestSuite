! ROCGM; for ocean assimilation; documented in CODG
      INTEGER :: I_TYP_GRD_M
      INTEGER :: NO_COLS_M
      INTEGER :: NO_ROWS_M
      INTEGER :: NO_LEVS_M
      INTEGER :: J_1_M
      INTEGER :: J_NROWS_M
      INTEGER :: JST_M
      INTEGER :: JFIN_M
      INTEGER :: NO_ROWS_GL_M
      INTEGER :: NO_ROWS_OFF_M
      INTEGER :: NO_COLS_HALO_M
      INTEGER :: NO_ROWS_HALO_M
      INTEGER :: MY_PE_M
      INTEGER :: NO_PE_M

      REAL :: POL_LAT_M
      REAL :: POL_LON_M
      REAL :: RLAT_N_M
      REAL :: RLON_W_M
      REAL :: RLAT_S_M
      REAL :: RLON_E_M
      REAL :: D_FIX_LAT_M
      REAL :: D_FIX_LON_M
      REAL :: RLAT_S_GL_M
      REAL :: RLAT_M(MAX_ROWS_M)
      REAL :: RLON_M(MAX_COLS_M)
      REAL :: DEP_LEVS_M(MAX_LEVS_M)

      LOGICAL :: LL_CYC_M
      LOGICAL :: LL_GLBL_M
! ROCGM end
