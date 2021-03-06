!C_DUSTGEN_CAM.........................................................
! Description: Contains parameters for mineral dust generation
!
! This version includes values of the constants that are currently used 
! in the NWP dust setup in the Crisis Area Models, which is based on
! an old configuration of the HadCM3 dust code. In the future, we hope
! to unify the code used in  HadGEM and NWP.
!
! Current Code Owner: Stephanie Woodward
! 
!----------------------------------------------------------------------
!
!
!     ! A,B and C in calculation of threshold friction velocity
!     ! Use value of A suitable for Iverson and White T*t(dry) 
!     ! and value of C tuned for use in the Southern Asia CAM 
      REAL, PARAMETER :: UST_A = 0.05
      REAL, PARAMETER :: UST_B = 0
      REAL, PARAMETER :: UST_C = -0.10

      REAL, PARAMETER :: HORIZ_C = 2.61 ! C in horizontal flux calc.
      REAL, PARAMETER :: VERT_A = 13.4 ! A in vertical flux calc
      REAL, PARAMETER :: VERT_B = -6. ! B in vertical flux calc

      REAL,PARAMETER :: RHOP = 2.65E+3  ! density of a dust particle
      REAL,PARAMETER :: Z0B = 0.0003    !roughness length for bare soil
      REAL Z0S(NDIV)  ! smooth roughness len (calc.d from part size)
      REAL DREP(NDIV) ! representative particle diameter
      REAL DMAX(NDIV) ! max diameter of particles in each div.
      REAL DMIN(NDIV) ! min diameter of particles in each div.
                         ! note that by using two arrays here we can set
                         ! up overlapping divisions, however this means
                         ! we have to be careful to make them consistent
!
      DATA Z0S/ .374894E-08, .118552E-07, .374894E-07, .118552E-06,     &
     &           .374894E-06, .118552E-05/
      DATA DREP/ .112468E-06, .355656E-06, .112468E-05, .355656E-05,    &
     &           .112468E-04, .355656E-04/
      DATA DMAX/2.0E-7,6.32456E-7,2.0E-6,                               &
     &          6.32456E-6,2.0E-5,6.32456E-5/
      DATA DMIN/6.32456E-8,2.0E-7,6.32456E-7,                           &
     &          2.0E-6,6.32456E-6,2.0E-5/
!.......................................................................
