! TROPFILT start

      REAL :: COS_TEMP

      INTEGER :: J_TO_FILTER_TEMP
      INTEGER :: SEG_TO_FILTER_TEMP
      INTEGER :: WORKING_FOR
      INTEGER :: SEG_START_TEMP
      INTEGER :: SEG_START_TEMP_EXTRA
      INTEGER :: SEG_LENGTH_TEMP
      INTEGER :: SEG_LENGTH_TEMP_EXTRA

#if defined(SHMEM)
       ! Arrays used in t3e mpp mode only for shmem calls

      REAL :: ETAA_TEMP(400)
      REAL :: UBTA_TEMP(400)
      REAL :: VBTA_TEMP(400)

      COMMON /SHMEM_TROPFL/                                             &
     &  ETAA_TEMP, UBTA_TEMP, VBTA_TEMP, COS_TEMP, J_TO_FILTER_TEMP,    &
     &  SEG_TO_FILTER_TEMP, WORKING_FOR, SEG_START_TEMP,                &
     &  SEG_START_TEMP_EXTRA, SEG_LENGTH_TEMP, SEG_LENGTH_TEMP_EXTRA

! TROPFILT end
#endif
