!C_DUSTSCAV.............................................................
! Description: Contains  mineral dust scavenging coefficients
! Current Code Owner: Stephanie Woodward
!
! History:
! Version  Date     Comment
! -------  ----     -------
!  1      12/02/03  Original Code.   Stephanie Woodward
!
! Mineral dust scavenging coefficients
!
      REAL :: KRAIN_DUST(NDIV)  ! scav. coeff. for rain
      REAL :: KSNOW_DUST(NDIV)  ! scav. coeff. for snow
!
      DATA KRAIN_DUST/2.e-5,2.e-5,3.e-5,6.e-5,4.e-4,4.e-4/
      DATA KSNOW_DUST/2.e-5,2.e-5,3.e-5,6.e-5,4.e-4,4.e-4/
!.......................................................................
