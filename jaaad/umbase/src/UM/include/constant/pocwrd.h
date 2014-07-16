! POCWRD contains number of bits per word on this machine
#if defined(CRAY)
      INTEGER, PARAMETER :: NO_BITS_WRD = 64
#else
      INTEGER, PARAMETER :: NO_BITS_WRD = 32
#endif
! POCWRD end
