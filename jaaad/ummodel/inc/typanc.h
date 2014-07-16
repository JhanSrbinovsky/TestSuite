! TYPANC start
!L
!L This COMDECK needs COMDECK TYPSIZE *CALLed first
!L
#include "conanc.h"
#if defined(ATMOS)
      ! Headers from  ancillary datasets

      INTEGER :: FTNANCILA
      INTEGER :: F
      INTEGER :: LOOKUP_START_ANCILA

      INTEGER :: FIXHD_ANCILA(LEN_FIXHD,NANCIL_DATASETSA)
      INTEGER :: INTHD_ANCILA(A_LEN_INTHD,NANCIL_DATASETSA)
      INTEGER :: LOOKUP_ANCILA(LEN1_LOOKUP,NANCIL_LOOKUPSA)

      REAL :: REALHD_ANCILA(A_LEN_REALHD,NANCIL_DATASETSA)

      COMMON/ANCILHDA/                                                  &
     &  FTNANCILA(NANCIL_DATASETSA),                                    &
     &  LOOKUP_START_ANCILA(NANCIL_DATASETSA)
#endif
