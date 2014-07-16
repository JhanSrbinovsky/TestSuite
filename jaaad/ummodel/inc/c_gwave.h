!
!  Description: This comdeck defines the constants for the 3B and 4A
!               versions of the Gravity Wave Drag Code. These are
!               tuneable parameters but are unlikely to be changed.
!
!  History:
!  Version    Date     Comment
!  -------    ----     -------
!    3.4     18/10/94  Original Version    J.R. Mitchell
!    4.3      7/03/97  Remove KAY_LEE (now set in RUNCNST) S.Webster
!    4.5     03/08/98  Add GAMMA_SATN (Used in 06_3B). D. Robinson
!    5.0     14/07/99  Remove redundant switches/variables: only keep
!                      s-version 06_3B parameters. R Rawlins
!    5.2     15/11/00  Set parameters for the 4A scheme.
!                                                  Stuart Webster.
!    5.3     16/10/01  Partition 3B and 4A scheme parameters.
!                                                  Stuart Webster.
!     5.3     21/09/01  Set parameters for the spectral (middle
!                       atmosphere) gravity wave forcing scheme.
!                       Warner and McIntyre, J. Atm. Sci., 2001,
!                       Scaife et al., Geophys. Res. Lett., 2000
!                       Scaife et al, J. Atm. Sci., 2002 give
!                       further details.
!                                                    Adam Scaife
!    5.4     Alters c_gwave.h include file to change from launch level
!            to launch eta to make less model level dependent.
!                                                    Adam Scaife
!    5.5     25/02/03  Remove 3B GWD parameter settings. Stuart Webster
!
!    6.2     16/06/05  Move CCL parameter to gw_ussp. Andrew Bushell
!
!
      ! Number of standard deviations above the mean orography of top
      !  of sub-grid mountains
      REAL,PARAMETER :: NSIGMA = 2.5

      ! Switch to determine which wave saturation test is used
      INTEGER,PARAMETER :: Amplitude_saturation = 1
      INTEGER,PARAMETER :: Stress_saturation    = 0

! SPECTRAL GRAVITY WAVE SCHEME PARAMETERS:

! LMINL = 1/max wavelength at launch
! LSTAR = 1/characteristic wavelength
! ETALAUNCH = eta for model launch

      REAL,PARAMETER:: LMINL = 1./20000
      REAL,PARAMETER:: LSTAR = 1./4300
      REAL,PARAMETER:: ETALAUNCH = 0.045
