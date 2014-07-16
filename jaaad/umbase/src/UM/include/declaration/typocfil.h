! TYPOCFIL start
      LOGICAL:: INITDN        ! Switch for initialising FILTER constants
      INTEGER:: ICBASE(IMTP1) ! Coeffs required for filtering routine
      INTEGER:: IDBASE(IMTP1) !                      "
      INTEGER:: IND(IMTX8)    !                      "
      REAL   :: COSSAV(LQMSUM)!                      "
      REAL   :: DENMSV(LHSUM) !                      "
      REAL   :: COSNPI(IMT)   !                      "

      ! Number of master pes I deal with as a slave process (U & V).
      INTEGER :: MAST_CNT_U(JMT)

      ! Number of master pes I deal with as a slave process (T & S).
      INTEGER :: MAST_CNT_T(JMT)

      ! No of master pes I deal with as a slave process (Free surf).
      INTEGER :: MAST_CNT_F

      ! Define the various indices in relation to master PEs
      INTEGER:: MAST_PE_U(KM*LSEGF,JMT) ! The master PEs (U/V)
      INTEGER:: MAST_PE_T(KM*LSEGF,JMT) ! The master PEs (T/S)
      INTEGER:: MAST_PE_F(LSEGF*JMT)    ! The master PEs (Fsurf)
      INTEGER:: MAST_K_U(KM*LSEGF,JMT)  ! The vert levels (U/V)
      INTEGER:: MAST_K_T(KM*LSEGF,JMT)  ! The vert levels (T/S)
      INTEGER:: MAST_SEG_U(KM*LSEGF,JMT)! The segments (U/V)
      INTEGER:: MAST_SEG_T(KM*LSEGF,JMT)! The segments (T/S)
      INTEGER:: MAST_SEG_F(LSEGF*JMT)   ! The segments (Fsurf)
      INTEGER:: MAST_ROW_F(LSEGF*JMT)   ! The rows (Free surf)

      ! Number of slave pes I deal with as a master process (U & V)
      INTEGER :: SLAV_CNT_U(JMT)

      ! Number of slave pes I deal with as a master process (T & S)
      INTEGER :: SLAV_CNT_T(JMT)

      ! Number of slave pes I deal with as a master process (Free surf)
      INTEGER :: SLAV_CNT_F

      INTEGER :: MAX_ROW_INDEX        ! Maximum number of rows on any pe
      ! The following use of JMT in dimensioning assumes an approximate
      ! regular load balance. If irregular load balancing is used, the
      ! only truely safe approach is to replace JMT with JMT_GLOBAL.
      INTEGER :: MIN_K_U(JMT,0:O_NPROC-1)
      INTEGER :: MIN_K_T(JMT,0:O_NPROC-1)
      INTEGER :: MAX_K_U(JMT,0:O_NPROC-1)
      INTEGER :: MAX_K_T(JMT,0:O_NPROC-1)
      LOGICAL :: MYSLAVE(JMT,0:O_NPROC-1)

      INTEGER :: MAST_CNT_ZTD            ! How many masters I work for
      INTEGER :: MAST_PE_ZTD(LSEGF*JMT_GLOBAL)
      INTEGER :: MAST_J_ZTD(LSEGF*JMT_GLOBAL)
      INTEGER :: MAST_L_ZTD(LSEGF*JMT_GLOBAL)
      INTEGER :: SLAV_PE_ZTD(LSEGF*JMT_GLOBAL)
      INTEGER :: SLAV_J_ZTD(LSEGF*JMT_GLOBAL)
      INTEGER :: SLAV_L_ZTD(LSEGF*JMT_GLOBAL)
      INTEGER :: ROW_CNT_FSF
      INTEGER :: SLAV_CNT_ZTD            ! How many slaves work for me
#include "comocfil.h"
! TYPOCFIL end
