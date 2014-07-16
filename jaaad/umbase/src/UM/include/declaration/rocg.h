! ROCG; for ocean assimilation; documented in CODG
      INTEGER :: I_TYP_GRD
      INTEGER :: NO_COLS
      INTEGER :: NO_ROWS
      INTEGER :: NO_GRD_LEVS
      INTEGER :: J_1
      INTEGER :: J_NROWS
      INTEGER :: JST
      INTEGER :: JFIN
      INTEGER :: NO_ROWS_GL
      INTEGER :: NO_ROWS_OFF
      INTEGER :: NO_COLS_HALO
      INTEGER :: NO_ROWS_HALO
      INTEGER :: MY_PE
      INTEGER :: NO_PE

      REAL :: POL_LAT
      REAL :: POL_LON
      REAL :: RLAT_N
      REAL :: RLON_W
      REAL :: RLAT_S
      REAL :: RLON_E
      REAL :: RLAT_S_GL
      REAL :: D_FIX_LAT
      REAL :: D_FIX_LON
      REAL :: RLAT(NO_ROWS)
      REAL :: RLON(NO_COLS)
      REAL :: DEP_GRD_LEVS(NO_GRD_LEVS)

      LOGICAL :: LL_CYC
      LOGICAL :: LL_GLBL
! ROCG end
