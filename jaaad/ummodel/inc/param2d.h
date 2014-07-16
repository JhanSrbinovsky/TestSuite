!
!     parameter param2d
!     Parameters for 2-D model/photolysis
!
!     Parameter    Description
!     ---------    -----------
!     nolat        Number of 2-D model latitudes
!     nolev        Number of 2-D model levels
!     nlphot       Number of 2-D model sub-levels used in the calculation
!                  of photolysis rates (approx. 1km spacing).
!     ntphot       Number of times per day that data from the 2-D model
!                  photolysis scheme is stored (time 3 is noon).
!     jpphio       Fortran i/o unit to write out and read in photolysis
!                  rates
!
!     ----------------------------------------------------------------
!
      integer nolat,nolev,nlphot,ntphot,jpphio,jpphin
!
      parameter(nolat=19,nolev=17,nlphot=51,ntphot=3)
      parameter(jpphio=84,jpphin=58)
