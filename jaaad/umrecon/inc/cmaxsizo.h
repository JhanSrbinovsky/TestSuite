#if defined(OCEAN)
! CMAXSIZO contains maximum sizes for dimensioning arrays
! of model constants whose sizes are configuration dependent. This
! allows constants to be read in from a NAMELIST file and maintain
! the flexibility of dynamic allocation for primary variables. The
! maximum sizes should agree with the maximum sizes implicit in the
! front-end User Interface.
      ! Max no. of ocean interface areas
      INTEGER,PARAMETER :: MAX_N_INTF_O = 4

      ! Max no. of ocean interface levels
      INTEGER,PARAMETER :: MAX_INTF_LEVELS_O = 40
! CMAXSIZO end
#endif
