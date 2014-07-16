!
!  Parameters for the photolysis scheme
!
      INTEGER NIN,JPLEV,JPLEVP1,JPWAV,JPLO,JPHI,JPTEM
      INTEGER JPO3P,JPS90,JPCHI,JPCHIN,NTAB,NSTD
      REAL SZAMAX,TMIN,TMAX,O3MIN,O3MAX
!
!     Number of altitudes in STDTO3 file
!      PARAMETER (NIN=37)
      PARAMETER (NIN=44)
!     
!     Number of levels in jtable.
      PARAMETER (JPLEVP1=64,JPLEV=JPLEVP1-1)
!     Number of zenith angles in jtable.
      PARAMETER (JPCHI=20)
!     Number of za greater than 90 deg in jtable.
      PARAMETER (JPS90=5)
!     Index of 90 degrees in jtable.
      PARAMETER (JPCHIN=JPCHI-JPS90)
!     Maximum zenith angle in jtable (degrees)
      PARAMETER (SZAMAX=98.0)
!     Number of wavelength intervals.
      PARAMETER (JPWAV=203)
!     Range of wavelength intervals to use.
      PARAMETER (JPLO=46,JPHI=203)
!     Number of temperatures in jtable
      PARAMETER (JPTEM=3)
!     Min and max temperatures
      PARAMETER (TMIN=200.0, TMAX=250.0)
!     Number of O3 profiles in jtable
      PARAMETER (JPO3P=5)
!     Min and max O3 profile factors
      PARAMETER (O3MIN=0.3, O3MAX=2.0)
!
!     Fortran channel for writing/reading photolysis table
      PARAMETER (NTAB=80)
!     Fortran channel for reading STDTO3
      PARAMETER (NSTD=72)
!
