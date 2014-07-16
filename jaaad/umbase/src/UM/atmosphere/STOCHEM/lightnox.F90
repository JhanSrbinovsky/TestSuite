#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE LIGHTNOX(cc_base_level,cc_top_level,t,acr,orog,        &
     &  land_fraction,totlnoxem,lnox,nflashes,z_top_of_model,           &
     &  first_constant_r_rho_level)
! ----------------------------------------------------------------------
! Purpose:
! Calculates production of NO by lightning.
!
! Method:
! Uses cloud height to calculate number of flashes per minute.
! (Price and Rind, 1992).
!
! Original Programmer: Michael Sanderson
!
! Current code owner: Michael Sanderson
!
! History:
! Date        Version     Comment
! -------     -------     -------
! 01/08/02      4.5         Created. Michael Sanderson
! 05/08/02      5.0         ND version  Colin Johnson
! 05/10/04      5.5         Vectorised code. M.G. Sanderson
! 01/03/06      6.2         Now uses Price and Rind (1992)
!                           formulation.  M.G. Sanderson

! Variables:
! Input
!  cc_base_level- convective cloud top level number
!  cc_top_level - convective cloud top level number
!  T        - 3D temperature array (K)
!  ACR      - Convective precipitation (kg m-2 s-1 == mm s-1)
!  OROG     - Orography (m)
!  LTDAT    - Number of STOCHEM grid rows per processor
!  Land_Fraction   - Land fraction
! From modules
!  Stochem_Grid_Area - Area of a grid cell on STOCHEM grid
!  Stochem_Advection_Step  - STOCHEM time step (s)
! Output
!  TOTLNOXEM- Total NOx produced by lightning
!  LNOX     - NOx produced by lightning (vmr s-1) on STOCHEM grid
! Local
!  Prof_Alt - 0-16 km in metres,
!  TROPICN, - Indices of tropical region on STOCHEM grid
!   TROPICS
!  NO_PER_CGFLASH - Molecules NO produced per cloud-to-ground (CG) flash
!                   intracloud flashes produce 10 times less NO.
!  MIN_CLOUD_DEPTH - Minimum cloud depth for lightning to occur;
!                    set to 5 km (Stockwell et al. JGR 1999)
!  AA,BB,CC,DD,EE - Coefficients of quartic polynomial for calculating
!             fraction of flashes that are CG from cold cloud depth.
!  CtoG_Frac - Fraction of flashes that are cloud-to-ground
!  ETA_0C   - Eta coordinate of zero C isotherm
!  NO_PROD  - molecules NO produced per timestep
!  Cloud_Depth - Cold cloud depth (depth of cloud from 0 C level to
!                cloud top).
!  N_PROFILE- Fraction of total NO produced in given eta level
!  NOX      - NOx produced by lightning (molecules NO per timestep) on
!             UM grid
! Functions
!  PROF_TO_ETA-Calculates N profiles for eta levels from fixed altitudes
!
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!-----------------------------------------------------------------------

      INTEGER, INTENT(IN)  :: first_constant_r_rho_level
      INTEGER, DIMENSION(nlonpe,nlatpe),INTENT(IN) :: cc_base_level
      INTEGER, DIMENSION(nlonpe,nlatpe),INTENT(IN) :: cc_top_level

      REAL, INTENT(IN)                          :: z_top_of_model
      REAL, DIMENSION(nlonpe,nlatpe),INTENT(IN) :: acr
      REAL, DIMENSION(nlonpe,nlatpe),INTENT(IN) :: orog
      REAL, DIMENSION(nlonpe,nlatpe),INTENT(IN) :: land_fraction
      REAL, DIMENSION(nlonpe,nlatpe,nmetlev),INTENT(IN) :: t

      REAL, INTENT(OUT) :: totlnoxem
      REAL, DIMENSION(nlnpe,nlpe,nlev), INTENT(OUT) :: lnox
      REAL, DIMENSION(nlnpe,nlpe,2),    INTENT(OUT) :: nflashes

      INTEGER :: i
      INTEGER :: j
      INTEGER :: k
      INTEGER :: l
      INTEGER :: m_storm

! Tropics assumed to be +/- 30 deg latitude
      REAL, PARAMETER :: tropics = 60.0
      REAL, PARAMETER :: tropicn = 120.0

! Per CG flash; IC is 10 times smaller
      REAL, PARAMETER :: no_per_cgflash = 6.7e26

! Min and max cloud depths (km) for lightning
      REAL, PARAMETER :: min_cloud_depth = 5.5,                         &
     &                   max_cloud_depth = 14.0

! Area of 1 deg square of ISCCP database (m2).
      REAL, PARAMETER :: cell_area = 6.25e10

! Fraction of lightning flashes that are cloud-to-ground calculated from
! quartic polynomial function of cold cloud depth with coefficients
! aa, bb, cc, dd, ee.

      REAL, PARAMETER :: aa = 0.021
      REAL, PARAMETER :: bb = -0.648
      REAL, PARAMETER :: cc = 7.49
      REAL, PARAMETER :: dd = -36.54
      REAL, PARAMETER :: ee = 64.09
      REAL, DIMENSION(0:nlkms), PARAMETER ::                            &
     &  prof_alt = (/ (1.0e03*i,i=0,nlkms) /)

      REAL :: nfm
      REAL :: nfc
      REAL :: tot_num_flashes
      REAL :: rlat
      REAL :: eta_0c
      REAL :: no_prod
      REAL :: ctog_frac
      REAL :: cloud_depth
      REAL :: cloud_height
      REAL :: eta_cloudbase
      REAL :: eta_cloudtop
      REAL, DIMENSION(nlev) :: n_profile
      REAL, DIMENSION(0:nlkms), SAVE :: eta_prof_alt
      REAL, DIMENSION(0:nlkms) :: ht_in_eta
      REAL, DIMENSION(nlonpe,nlatpe,nlev) :: nox
      REAL, DIMENSION(nlonpe,nlatpe,2) :: um_flashes

      CHARACTER(LEN=72) :: cmessage
      LOGICAL :: first=.TRUE.

      nox = 0.0
      um_flashes = 0.0

! ETA_PROF_ALT contains the eta coordinates of the NOx
!   production profiles at z = 0, 1, 2, ... 16 km.
! Approximate as zero orography is assumed.
      IF (first) THEN
        DO i=0,nlkms
! DEPENDS ON: rtoeta
          eta_prof_alt(i) = RTOETA(prof_alt(i)+earth_radius,0.0,        &
     &     z_top_of_model,first_constant_r_rho_level)
        END DO
        first = .FALSE.
      END IF

      DO j = 1, rowspe
        DO i = 1, nlonpe
          IF (acr(i,j) > 1.0e-8) THEN ! If no convective rain no lightni

! Calculate cold cloud depth.
! First, find height of 0 C isotherm in eta coordinates (ETA_0C)
            IF (t(i,j,1) < zerodegc) THEN
              eta_0c = eta_theta(2)
            ELSE
              DO k = 2, nmetlev
                IF (t(i,j,k-1)>=zerodegc .AND. t(i,j,k)<zerodegc) EXIT
              END DO
              IF (k > nmetlev-1) k = nmetlev-1
! Linear interpolation
              eta_0c = eta_theta(k-1)+((eta_theta(k)-eta_theta(k-1))/   &
     &          (t(i,j,k-1)-t(i,j,k)))*(t(i,j,k-1)-zerodegc)
            END IF

            IF (eta_0c > 1.0 .OR. eta_0c < 0.0) THEN
              cmessage = 'ETA_0C out of range'
              write(6,*) cmessage,eta_0c
! DEPENDS ON: ereport
              CALL EREPORT('LIGHTNOX',1,cmessage)
            END IF

! Find cloud base eta
            eta_cloudbase = eta_rho(cc_base_level(i,j))

! Find cloud top eta
            eta_cloudtop = eta_rho(cc_top_level(i,j))

! Calculate cloud depth from 0 C level to cloud top in km
! DEPENDS ON: etator
            cloud_depth = 1.0e-3*(ETATOR(eta_cloudtop,orog(i,j),        &
     &        z_top_of_model,first_constant_r_rho_level)-               &
! DEPENDS ON: etator
     &        ETATOR(eta_0c,orog(i,j),                                  &
     &        z_top_of_model,first_constant_r_rho_level))
            IF (cloud_depth < min_cloud_depth) CYCLE
            IF (cloud_depth > max_cloud_depth)                          &
     &        cloud_depth = max_cloud_depth

! Calculate fraction of lightning flashes that are cloud-to-ground
            ctog_frac = 1.0 / (ee + cloud_depth * (dd + cloud_depth *   &
     &        (cc + cloud_depth * (bb + cloud_depth * aa))))

! Now calculate absolute cloud height in km.
! DEPENDS ON: etator
            cloud_height = 1.0e-3*(ETATOR(eta_cloudtop,orog(i,j),       &
     &        z_top_of_model,first_constant_r_rho_level)-               &
! DEPENDS ON: etator
     &        ETATOR(eta_cloudbase,orog(i,j),                           &
     &        z_top_of_model,first_constant_r_rho_level))

! Identify type of storm (maritime or continental)
! Assumes midlat maritime storms are same as tropical ones
            rlat = latm(j+lobound-1)
            IF (land_fraction(i,j) > 0.5) THEN
              tot_num_flashes = 3.44e-5 * (cloud_height**4.9)
              IF (rlat > tropicn .AND. rlat <= tropics) THEN  ! In the t
                m_storm = 3       ! Tropical continental storm
              ELSE
                m_storm = 1       ! Midlatitude continental storm
              END IF
            ELSE
              tot_num_flashes = 6.4e-4 * (cloud_height**1.73)
              m_storm = 2         ! Tropical maritime storm
            END IF

! Convert to flashes m-2 s-1
            tot_num_flashes = tot_num_flashes / (cell_area * 60.0)

! Calculate NO Production (molecules NO m-2 s-1)
            no_prod = tot_num_flashes * no_per_cgflash * (ctog_frac  +  &
     &        (1.0-ctog_frac)*0.1)

! Now distribute the NO produced in the column using the appropriate
! profile.  The profile is strictly a mass distribution.
! Scale the heights of the profiles to match the new cloud height

            ht_in_eta(0) = eta_prof_alt(0)
            DO k = 1, nlkms
              ht_in_eta(k) = (eta_prof_alt(k) *                         &
     &          eta_cloudtop / eta_prof_alt(nlkms))
            END DO

! Redistribute the appropriate profile onto the vertical STOCHEM grid
! DEPENDS ON: prof_to_nlev
            CALL PROF_TO_NLEV(m_storm,ht_in_eta,n_profile)

! Multiply the NO production fractions by the total NO produced
! and the lightning NOx multiplier
            nox(i,j,:) = n_profile(:) * no_prod  ! * lightmult
!
! Store the number of cloud-to-ground and intercloud flashes
            um_flashes(i,j,1) = tot_num_flashes * ctog_frac
            um_flashes(i,j,2) = tot_num_flashes * (1.0 - ctog_frac)
          END IF
        END DO
      END DO

! Convert emission and number of flashes from UM grid to STOCHEM grid
! Array nox contains NO produced m-2 s-1 on each level

! DEPENDS ON: met2data
      CALL MET2DATA(lnox,nox,nlev,nlev)
! DEPENDS ON: met2data
      CALL MET2DATA(nflashes,um_flashes,2,2)

! Convert NO production from molecules m-2 s-1 to vmr NO s-1

      DO j = 1, nlpe
        DO k = 1, nlev
          lnox(:,j,k) = lnox(:,j,k) * stochem_grid_area(j+ltdat-1) /    &
     &      lmolec
        END DO
      END DO

      totlnoxem = SUM(lnox)

! Convert number of flashes from flashes m-2 s-1 to
! flashes m-2 timestep-1

      DO j = 1, nlpe
        nflashes(:,j,1) = nflashes(:,j,1) * stochem_advection_step
        nflashes(:,j,2) = nflashes(:,j,2) * stochem_advection_step
      END DO

      END SUBROUTINE LIGHTNOX
#endif
