#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Current code owner: M.G. Sanderson
!
!  Model
! Version   Date          Modification history
!   4.5   08/04/02   Created for H2 only. M.G. Sanderson
!   4.5   09/05/02   Extended for other species. M.G. Sanderson
!   6.1   26/11/04   Vectorised version of dry deposition code.
!                    M.G. Sanderson
!   6.2   03/02/05   Now treats HNO3/H2O2/ROOH separately, as the
!                    surface resistances for these species are
!                    invariant.  M.G. Sanderson.
!   6.2   01/03/06   Ozone cuticular resistance varies with LAI.
!                    M.G. Sanderson
!
!    Programming standard:
!
      SUBROUTINE SURF_DDEP_RES(t0,p0,rh,soilmc,so4_vd,gsf,stcon,t0tile, &
     &  lai_ft,canwc,rc,o3_stom_frac)
!
!-----------------------------------------------------------------------
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!-----------------------------------------------------------------------

      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) ::                     &
     &  t0                                                              &
                  ! Surface temperature (K)
     & ,p0                                                              &
                  ! Surface pressure (Pa)
     & ,rh                                                              &
                  ! Relative humidity (fraction)
     & ,soilmc                                                          &
                  ! Soil moisture content (Fraction by volume)
     &, so4_vd    ! Aerosol deposition velocity (m s-1).

      REAL, DIMENSION(nlonpe,nlatpe,ntype), INTENT(IN) ::               &
     &  gsf                                                             &
                  ! Global surface fractions
     & ,t0tile    ! Surface temperature on tiles (K)

      REAL, DIMENSION(nlonpe,nlatpe,npft), INTENT(IN) ::                &
     &  stcon                                                           &
                  ! Stomatal conductance (m s-1)
     & ,lai_ft                                                          &
                  ! Leaf area index (m2 leaf m-2)
     & ,canwc     ! Canopy water content (mm)

! Surface resistance on tiles (s m-1).
      REAL, DIMENSION(nlonpe,nlatpe,ntype,nc), INTENT(OUT) :: rc
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(OUT) :: o3_stom_frac
!
! Local variables
!
      LOGICAL :: first = .true.

      INTEGER ::                                                        &
     &  i                                                               &
                  ! Loop count over longitudes
     & ,k                                                               &
                  ! Loop count over latitudes
     & ,j                                                               &
                  ! Loop count over species that deposit
     & ,n                                                               &
                  ! Loop count over tiles
     & ,nj        ! Species depositing
!
      REAL ::                                                           &
     &  sm                                                              &
                            ! Soil moisture content of gridbox.
     & ,rr                                                              &
                            ! General temporary store for resistances.
     & ,ts                                                              &
                            ! Temperature of a particular tile.
     & ,y                                                               &
                            ! Latitude in degrees.
     & ,f                                                               &
                            ! Factor to modify CH4 uptake fluxes
     & ,r_cuticle                                                       &
                            ! Cuticular resistance
     & ,mml = 1.008e5                                                   &
                            ! Factor to convert methane flux to dry dep
     & ,r_wet_o3 = 500.0                                                &
                            ! Wet soil surface resistance for ozone (s m
     & ,cuticle_o3 = 5000.0 ! Constant for calculation of O3 cuticular r
!
! MML - Used to convert methane flux in ug m-2 h-1 to dry dep vel in
!  m s-1; MML=3600*0.016*1.0E9*1.75E-6, where 0.016=RMM methane (kg),
!  1.0E9 converts ug -> kg, 1.75E-6 = assumed CH4 vmr
!
      REAL, DIMENSION(nlonpe,npft) :: r_cut_o3 ! Cuticular Resistance fo
      REAL, DIMENSION(nlonpe,npft,5) :: r_stom ! Stomatal resistance
      REAL, DIMENSION(ntype,nc), SAVE :: rsurf ! Standard surface resist

      REAL :: tundra_s_lim = 150.0 ! Southern limit of tundra (reverse g
!
! Following arrays used to set up surface resistance array rsurf.
!
      REAL,PARAMETER,DIMENSION(ntype) :: zero =                         &
     &  (/ (r_null,i=1,ntype) /)
      REAL,PARAMETER,DIMENSION(ntype) :: rooh =                         &
     &  (/ 30.0, (10.0,i=1,ntype-1) /)
      REAL,PARAMETER,DIMENSION(ntype) :: hno3 =                         &
     &  (/ (10.0,i=1,ntype) /)
      REAL,PARAMETER,DIMENSION(ntype) :: aerosol =                      &
     &  (/ (r_null, i=1,npft+1), 1000.0, r_null, 20000.0 /)
!
! Hydrogen - linear dependence on soil moisture (except savannah)
      REAL, DIMENSION(npft) :: h2dd_c = (/ 0.00197, 0.00197,            &
     &  0.00177, 1.2346, 0.0001 /)
      REAL, DIMENSION(npft) :: h2dd_m = (/ -0.00419, -0.00419,          &
     &  -0.00414, -0.472, 0.0 /)
      REAL :: h2dd_q = 0.27 ! Quadratic term for H2 loss to savannah
!
! Resistances for Tundra if different to standard value in rsurf().
      REAL, DIMENSION(5) :: r_tundra =                                  &
     &  (/ 1200.0, 25000.0, 800.0, 3850.0, 1100.0/)
!           NO2      CO      O3      H2     PAN
!
! CH4 uptake fluxes, in ug m-2 hr-1
      REAL, DIMENSION(ntype) :: ch4_up_flux =                           &
     &  (/ 39.5, 50.0, 30.0, 37.0, 27.5, 0.0, 0.0, 27.5, 0.0 /)
!
! CH4 loss to tundra - Cubic polynomial fit to data. N.B. Loss flux is
! in units of ug(CH4) m-2 s-1
      REAL, DIMENSION(4) :: ch4dd_tun = (/ -4.757e-6, 4.0288e-3,        &
     &                                     -1.13592, 106.636 /)
!
! HNO3 dry dep to ice; quadratic dependence
      REAL, DIMENSION(3) :: hno3dd_ice = (/ -13.57, 6841.9, -857410.6 /)
!
! SO2 dry dep to snow/ice; quadratic dependence
      REAL, DIMENSION(3) :: so2dd_ice = (/ 0.0001, 0.003308, 0.1637 /)
!
! Diffusion correction for stomatal conductance
      REAL, DIMENSION(5) ::                                             &
     &  dif = (/ 1.6, 1.6, 2.6, 1.9, 0.97 /)
!                NO2   O3  PAN  SO2  NH3
!
      LOGICAL, DIMENSION(nlonpe,ntype) ::                               &
     &  todo     ! True if tile fraction > 0.0
      LOGICAL, DIMENSION(nlonpe,nlatpe) ::                              &
     &  microb   ! True if T > 5 C and RH > 40%
!                  (i.e. microbes in soil are active).
#include "trif.h"
!
! Set up standard resistance array rsurf on first call only
      IF (first) THEN
        DO n = 1, ntype
          DO i = 1, nc
            rsurf(n,i) = r_null
          END DO
        END DO
!
! Standard surface resistances (s m-1). Values are for 9 tiles in
! order: Broadleaved trees, Needleleaf trees, C3 Grass, C4 Grass,
! Shrub, Urban, Water, Bare Soil, Ice.
!
        rsurf(:,i_no2)=(/225.,225.,400.,400.,600.,1200.,2600.,1200.,    &
     &    3500. /)
        rsurf(:,i_co)=(/3700.,7300.,4550.,1960.,r_null,r_null,r_null,   &
     &    r_null,r_null /)
        rsurf(:,i_ch4)=zero
        rsurf(:,i_o3)=(/200.,200.,200.,200.,400.,800.,2200.,800.,       &
     &    2500. /)
        rsurf(:,i_h2)=zero
        rsurf(:,i_pan)=(/500.,500.,500.,500.,500.,r_null,12500.,        &
     &    500.0, 12500. /)
        rsurf(:,i_ch3ooh)=rooh
        rsurf(:,i_so2)=(/100.,100.,150.,350.,400.,400.,10.0,700.0,      &
     &    r_null /)
        rsurf(:,i_sa)=aerosol
        rsurf(:,i_ammsul)=aerosol
        rsurf(:,i_naer)=aerosol
        rsurf(:,i_msa)=aerosol
        rsurf(:,i_orgnit)=aerosol
        rsurf(:,i_hno3)=hno3
        rsurf(:,i_h2o2)=hno3
        rsurf(:,i_c3h7ooh)=rooh
        rsurf(:,i_c2h5ooh)=rooh
        rsurf(:,i_c4h9ooh)=rooh
        rsurf(:,i_isopooh)=rooh
        rsurf(:,i_mvkooh)=rooh
        rsurf(:,i_nh3)=hno3

        first = .false.

      END IF
!
      o3_stom_frac = 0.0
!
! Set flag for microbial action
!
      DO k = 1, rowspe
        DO i = 1, nlonpe
          microb(i,k) = (rh(i,k) > 0.40 .AND. t0(i,k) > 278.0)
        END DO
      END DO
!
! Set surface resistances to standard values. rsurf is the resistance of
! the soil, rock, water etc. Set all tiles to standard values. These
! values will be modified below as necessary. Extra terms for vegetated
! tiles (stomatal, cuticular) will be added if required. Loop over all
! parts of array rc to ensure all of it is assigned a value.
!
      DO n = 1, ntype
        DO k = 1, rowspe
          DO i = 1, nlonpe
            IF (gsf(i,k,n) > 0.0) THEN
              DO j = 1, nc
                rc(i,k,n,j) = rsurf(n,j)
              END DO
            ELSE
              DO j = 1, nc
                rc(i,k,n,j) = r_null
              END DO
            END IF
          END DO
        END DO
      END DO
!
! Now begin assigning specific surface resistances.
!
! Aerosols. Assume vds is valid for all land types and aerosol types
!
      DO n = 1, npft+1
        DO k = 1, rowspe
          DO i = 1, nlonpe
            IF (so4_vd(i,k) > 0.0 .AND. gsf(i,k,n) > 0.0) THEN
              rr = 1.0 / so4_vd(i,k)
              rc(i,k,n,i_sa)     = rr
              rc(i,k,n,i_ammsul) = rr
              rc(i,k,n,i_naer)   = rr
              rc(i,k,n,i_orgnit) = rr
              rc(i,k,n,i_msa)    = rr
            END IF
          END DO
        END DO
      END DO
      n = soil
      DO k = 1, rowspe
        DO i = 1, nlonpe
          IF (so4_vd(i,k) > 0.0 .AND. gsf(i,k,n) > 0.0) THEN
            rr = 1.0 / so4_vd(i,k)
            rc(i,k,n,i_sa)     = rr
            rc(i,k,n,i_ammsul) = rr
            rc(i,k,n,i_naer)   = rr
            rc(i,k,n,i_orgnit) = rr
            rc(i,k,n,i_msa)    = rr
          END IF
        END DO
      END DO
!
! Tundra; some species have different standard values. Assume bare soil
! is also tundra in this region. Tundra is defined as shrub/bare soil
! north of 60 degrees north. Note that the grid runs S -> N, so zero
! is the S pole and 180 is the N pole.
!
      DO k = 1, rowspe
        n = npft
        y = latm(k+lobound-1)
        IF (y > tundra_s_lim) THEN
          DO i = 1, nlonpe
            IF (gsf(i,k,n) > 0.0) THEN
              rc(i,k,n,i_no2) = r_tundra(1)         ! NO2
              rc(i,k,n,i_o3)  = r_tundra(3)         ! O3
              rc(i,k,n,i_pan) = r_tundra(5)         ! PAN
            END IF
            IF (gsf(i,k,soil) > 0.0) THEN
              rc(i,k,soil,i_no2) = r_tundra(1)      ! NO2
              rc(i,k,soil,i_o3)  = r_tundra(3)      ! O3
              rc(i,k,soil,i_pan) = r_tundra(5)      ! PAN
            END IF
          END DO
!
! Only assign values for H2, CO if microbes are active.
          DO i = 1, nlonpe
            IF (microb(i,k)) THEN
              IF (gsf(i,k,n) > 0.0) THEN
                rc(i,k,n,i_co)  = r_tundra(2)       ! CO
                rc(i,k,n,i_h2)  = r_tundra(4)       ! H2
              END IF
              IF (gsf(i,k,soil) > 0.0) THEN
                rc(i,k,soil,i_co)  = r_tundra(2)    ! CO
                rc(i,k,soil,i_h2)  = r_tundra(4)    ! H2
              END IF
            END IF
          END DO
        END IF
      END DO
!
! O3: Change land deposition values if surface is wet; soil moisture
! value of 0.3 fairly arbitrary.
!
      DO k = 1, rowspe
        DO i = 1, nlonpe
          IF (soilmc(i,k) > 0.3) THEN
            DO n = 1, npft
              IF (gsf(i,k,n) > 0.0) THEN
                rc(i,k,n,i_o3) = r_wet_o3
              END IF
            END DO
            IF (gsf(i,k,soil) > 0.0) THEN
              rc(i,k,soil,i_o3) = r_wet_o3
            END IF
          END IF
        END DO
      END DO
!
! All standard resistances have been applied. Now calculate surface
! resistances that depend on temperature, soil moisture and humidity.
! Stomatal and cuticular resistances will also be calculated here.
!
      DO k = 1, rowspe
        y = latm(k+lobound-1) - 90.0
!
! Set logical for surface type fractions > 0.0.
        DO n = 1, ntype
          DO i = 1, nlonpe
            todo(i,n) = (gsf(i,k,n) > 0.0)
          END DO
        END DO
!
! Calculate stomatal resistances for this latitude.
! Assumes that there is a mesophyll resistance for NO2 equal to 50% of
! the stomatal resistance.
!
        DO n = 1, npft
          DO i = 1, nlonpe
            IF (todo(i,n) .AND. stcon(i,k,n) > glmin(n)) THEN
              r_stom(i,n,1) = 1.5 * dif(1) / stcon(i,k,n)
              r_stom(i,n,2) = dif(2) / stcon(i,k,n)
              r_stom(i,n,3) = dif(3) / stcon(i,k,n)
              r_stom(i,n,4) = dif(4) / stcon(i,k,n)
              r_stom(i,n,5) = dif(5) / stcon(i,k,n)
            ELSE
              r_stom(i,n,1) = r_null
              r_stom(i,n,2) = r_null
              r_stom(i,n,3) = r_null
              r_stom(i,n,4) = r_null
              r_stom(i,n,5) = r_null
            END IF
          END DO
        END DO
!
! Set cuticular resistance for ozone.
!
        DO n = 1, npft
          DO i = 1, nlonpe
            IF (todo(i,n) .AND. lai_ft(i,k,n) > 0.0) THEN
              r_cut_o3(i,n) = cuticle_o3 / lai_ft(i,k,n)
            ELSE
              r_cut_o3(i,n) = r_null
            END IF
          END DO
        END DO
!
! --------------------------------------------------------------------
! Broadleaved Forests
        n = 1
!
! H2 dry dep vel has linear dependence on soil moisture
! Limit sm to avoid excessively high deposition velocities
        DO i = 1, nlonpe
          IF (todo(i,n) .AND. microb(i,k)) THEN
            sm = MAX(soilmc(i,k),0.1)
            rc(i,k,n,i_h2) = 1.0 / (h2dd_m(n) * sm + h2dd_c(n))
          END IF
        END DO
!
! CH4: Calculate an uptake flux initially. The uptake flux depends on
! soil moisture, based on results of Reay et al. (2001).
        DO i = 1, nlonpe
          IF (todo(i,n) .AND. microb(i,k)) THEN
            sm = soilmc(i,k)
            IF (sm < 0.16) THEN
              f = sm / 0.16
            ELSE IF (sm > 0.30) THEN
              f = (0.50 - sm) / 0.20
            ELSE
              f = 1.0
            END IF
            f = MAX(f,0.0)
            rc(i,k,n,i_ch4) = ch4_up_flux(n) * f
          END IF
        END DO
!
! NO2: Calculate plant deposition terms.
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = (1.0/r_stom(i,n,1)) + (1.0/rc(i,k,n,i_no2))
            rc(i,k,n,i_no2) = 1.0 / rr
          END IF
        END DO
!
! O3: Calculate plant deposition terms.
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = (1.0/r_stom(i,n,2)) + (1.0/r_cut_o3(i,n)) +            &
     &        (1.0/rc(i,k,n,i_o3))
            rc(i,k,n,i_o3) = 1.0 / rr
            o3_stom_frac(i,k) = o3_stom_frac(i,k) +                     &
     &        gsf(i,k,n) / (rr * r_stom(i,n,2))
          END IF
        END DO
!
! PAN: Calculate plant deposition terms.
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = (1.0/r_stom(i,n,3)) + (1.0/rc(i,k,n,i_pan))
            rc(i,k,n,i_pan) = 1.0 / rr
          END IF
        END DO
!
! SO2: Calculate plant deposition terms.
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = rh(i,k)
            ts = t0tile(i,k,n)
            IF (ts < 268.0) THEN
              r_cuticle = 500.0
            ELSE IF (ts >= 268.0 .AND. ts < 272.0) THEN
              r_cuticle = 200.0
            ELSE
! If CANWC > 0.25 mm assume cuticle is wet (Wyers/Erisman Atm Env 1998)
              IF (canwc(i,k,n) > 0.25) THEN
                r_cuticle = 1.0
              ELSE
                IF (rr >= 0.813) THEN
                  r_cuticle = 5.8e11 * EXP(-27.8 * rr)
                ELSE
                  r_cuticle = 2.5e4 * EXP(-6.93 * rr)
                END IF
              END IF
            END IF
            rr = (1.0/r_stom(i,n,4)) + (1.0/r_cuticle) +                &
     &        (1.0/rc(i,k,n,i_so2))
            rc(i,k,n,i_so2) = 1.0 / rr
          END IF
        END DO
!
! NH3: Calculate plant deposition terms.
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = rh(i,k)
            ts = t0tile(i,k,n)
            IF (ts < 268.0) THEN
              r_cuticle = 1000.0
            ELSE IF (ts >= 268.0 .AND. ts < 272.0) THEN
              r_cuticle = 200.0
            ELSE
! If CANWC > 0.25 mm assume cuticle is wet (Wyers/Erisman Atm Env 1998)
              IF (canwc(i,k,n) > 0.25) THEN
                r_cuticle = 10.0
              ELSE
                r_cuticle = 5.0 *                                       &
     &            LOG10(ts-zerodegc+2.0) * EXP((1.0-rr)/0.12)
              END IF
            END IF
            rr = (1.0/r_stom(i,n,5)) + (1.0/r_cuticle) +                &
     &        (1.0/rc(i,k,n,i_nh3))
            rc(i,k,n,i_nh3) = 1.0 / rr
          END IF
        END DO
!
! --------------------------------------------------------------------
! Needleleaved Forests
        n = 2
!
! H2 dry dep vel has linear dependence on soil moisture
! Limit sm to avoid excessively high deposition velocities
        DO i = 1, nlonpe
          IF (todo(i,n) .AND. microb(i,k)) THEN
            sm = MAX(soilmc(i,k),0.1)
            rc(i,k,n,i_h2) = 1.0 / (h2dd_m(n) * sm + h2dd_c(n))
          END IF
        END DO
!
! CH4: Calculate an uptake flux initially. The uptake flux depends on
! soil moisture, based on results of Reay et al. (2001).
        DO i = 1, nlonpe
          IF (todo(i,n) .AND. microb(i,k)) THEN
            sm = soilmc(i,k)
            IF (sm < 0.16) THEN
              f = sm / 0.16
            ELSE IF (sm > 0.30) THEN
              f = (0.50 - sm) / 0.20
            ELSE
              f = 1.0
            END IF
            f = MAX(f,0.0)
            rc(i,k,n,i_ch4) = ch4_up_flux(n) * f
          END IF
        END DO
!
! NO2: Calculate plant deposition terms.
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = (1.0/r_stom(i,n,1)) + (1.0/rc(i,k,n,i_no2))
            rc(i,k,n,i_no2) = 1.0 / rr
          END IF
        END DO
!
! O3:
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = (1.0/r_stom(i,n,2)) + (1.0/r_cut_o3(i,n)) +            &
     &        (1.0/rc(i,k,n,i_o3))
            rc(i,k,n,i_o3) = 1.0 / rr
            o3_stom_frac(i,k) = o3_stom_frac(i,k) +                     &
     &        gsf(i,k,n) / (rr * r_stom(i,n,2))
          END IF
        END DO
!
! PAN:
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = (1.0/r_stom(i,n,3)) + (1.0/rc(i,k,n,i_pan))
            rc(i,k,n,i_pan) = 1.0 / rr
          END IF
        END DO
!
! SO2:
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = rh(i,k)
            ts = t0tile(i,k,n)
            IF (ts < 268.0) THEN
              r_cuticle = 500.0
            ELSE IF (ts >= 268.0 .AND. ts < 272.0) THEN
              r_cuticle = 200.0
            ELSE
! If CANWC > 0.25 mm assume cuticle is wet (Wyers/Erisman Atm Env 1998)
              IF (canwc(i,k,n) > 0.25) THEN
                r_cuticle = 1.0
              ELSE
                IF (rr >= 0.813) THEN
                  r_cuticle = 5.8e11 * EXP(-27.8 * rr)
                ELSE
                  r_cuticle = 2.5e4 * EXP(-6.93 * rr)
                END IF
              END IF
            END IF
            rr = (1.0/r_stom(i,n,4)) + (1.0/r_cuticle) +                &
     &        (1.0/rc(i,k,n,i_so2))
            rc(i,k,n,i_so2) = 1.0 / rr
          END IF
        END DO
!
! NH3:
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = rh(i,k)
            ts = t0tile(i,k,n)
            IF (ts < 268.0) THEN
              r_cuticle = 1000.0
            ELSE IF (ts >= 268.0 .AND. ts < 272.0) THEN
              r_cuticle = 200.0
            ELSE
! If CANWC > 0.25 mm assume cuticle is wet (Wyers/Erisman Atm Env 1998)
              IF (canwc(i,k,n) > 0.25) THEN
                r_cuticle = 10.0
              ELSE
                r_cuticle = 5.0 *                                       &
     &            LOG10(ts-zerodegc+2.0) * EXP((1.0-rr)/0.12)
              END IF
            END IF
            rr = (1.0/r_stom(i,n,5)) + (1.0/r_cuticle) +                &
     &        (1.0/rc(i,k,n,i_nh3))
            rc(i,k,n,i_nh3) = 1.0 / rr
          END IF
        END DO
!
! --------------------------------------------------------------------
! C3 grasslands
        n = 3
!
! H2 dry dep vel has linear dependence on soil moisture
! Limit sm to avoid excessively high deposition velocities
        DO i = 1, nlonpe
          IF (todo(i,n) .AND. microb(i,k)) THEN
            sm = MAX(soilmc(i,k),0.1)
            rc(i,k,n,i_h2) = 1.0 / (h2dd_m(n) * sm + h2dd_c(n))
          END IF
        END DO
!
! CH4: Calculate an uptake flux initially. The uptake flux depends on
! soil moisture, based on results of Reay et al. (2001).
        DO i = 1, nlonpe
          IF (todo(i,n) .AND. microb(i,k)) THEN
            sm = soilmc(i,k)
            IF (sm < 0.16) THEN
              f = sm / 0.16
            ELSE IF (sm > 0.30) THEN
              f = (0.50 - sm) / 0.20
            ELSE
              f = 1.0
            END IF
            f = MAX(f,0.0)
            rc(i,k,n,i_ch4) = ch4_up_flux(n) * f
          END IF
        END DO
!
! NO2: Calculate plant deposition terms.
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = (1.0/r_stom(i,n,1)) + (1.0/rc(i,k,n,i_no2))
            rc(i,k,n,i_no2) = 1.0 / rr
          END IF
        END DO
!
! O3:
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = (1.0/r_stom(i,n,2)) + (1.0/r_cut_o3(i,n)) +            &
     &        (1.0/rc(i,k,n,i_o3))
            rc(i,k,n,i_o3) = 1.0 / rr
            o3_stom_frac(i,k) = o3_stom_frac(i,k) +                     &
     &        gsf(i,k,n) / (rr * r_stom(i,n,2))
          END IF
        END DO
!
! PAN:
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = (1.0/r_stom(i,n,3)) + (1.0/rc(i,k,n,i_pan))
            rc(i,k,n,i_pan) = 1.0 / rr
          END IF
        END DO
!
! SO2:
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = rh(i,k)
            ts = t0tile(i,k,n)
            IF (ts < 268.0) THEN
              r_cuticle = 500.0
            ELSE IF (ts >= 268.0 .AND. ts < 272.0) THEN
              r_cuticle = 200.0
            ELSE
! If CANWC > 0.25 mm assume cuticle is wet (Wyers/Erisman Atm Env 1998)
              IF (canwc(i,k,n) > 0.25) THEN
                r_cuticle = 1.0
              ELSE
                IF (rr >= 0.813) THEN
                  r_cuticle = 5.8e11 * EXP(-27.8 * rr)
                ELSE
                  r_cuticle = 2.5e4 * EXP(-6.93 * rr)
                END IF
              END IF
            END IF
            rr = (1.0/r_stom(i,n,4)) + (1.0/r_cuticle) +                &
     &        (1.0/rc(i,k,n,i_so2))
            rc(i,k,n,i_so2) = 1.0 / rr
          END IF
        END DO
!
! NH3:
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = rh(i,k)
            ts = t0tile(i,k,n)
            IF (ts < 268.0) THEN
              r_cuticle = 1000.0
            ELSE IF (ts >= 268.0 .AND. ts < 272.0) THEN
              r_cuticle = 200.0
            ELSE
! If CANWC > 0.25 mm assume cuticle is wet (Wyers/Erisman Atm Env 1998)
              IF (canwc(i,k,n) > 0.25) THEN
                r_cuticle = 10.0
              ELSE
                r_cuticle = 5.0 *                                       &
     &            LOG10(ts-zerodegc+2.0) * EXP((1.0-rr)/0.12)
              END IF
            END IF
            rr = (1.0/r_stom(i,n,5)) + (1.0/r_cuticle) +                &
     &        (1.0/rc(i,k,n,i_nh3))
            rc(i,k,n,i_nh3) = 1.0 / rr
          END IF
        END DO
!
! --------------------------------------------------------------------
! C4 grasslands
        n = 4
!
! H2 dry dep has quadratic-log dependence on soil moisture
        DO i = 1, nlonpe
          IF (todo(i,n) .AND. microb(i,k)) THEN
            sm = LOG(MAX(soilmc(i,k),0.1))
            rr = (h2dd_c(4) + sm*(h2dd_m(4) + sm*h2dd_q)) * 1.0e-4
            IF (rr > 0.00131) rr = 0.00131 ! Conrad/Seiler Max value
            rc(i,k,n,i_h2) = 1.0 / rr
          END IF
        END DO
!
! CH4: Calculate an uptake flux initially. The uptake flux depends on
! soil moisture, based on results of Reay et al. (2001).
        DO i = 1, nlonpe
          IF (todo(i,n) .AND. microb(i,k)) THEN
            sm = soilmc(i,k)
            IF (sm < 0.16) THEN
              f = sm / 0.16
            ELSE IF (sm > 0.30) THEN
              f = (0.50 - sm) / 0.20
            ELSE
              f = 1.0
            END IF
            f = MAX(f,0.0)
            rc(i,k,n,i_ch4) = ch4_up_flux(n) * f
          END IF
        END DO
!
! NO2: Calculate plant deposition terms.
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = (1.0/r_stom(i,n,1)) + (1.0/rc(i,k,n,i_no2))
            rc(i,k,n,i_no2) = 1.0 / rr
          END IF
        END DO
!
! O3:
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = (1.0/r_stom(i,n,2)) + (1.0/r_cut_o3(i,n)) +            &
     &        (1.0/rc(i,k,n,i_o3))
            rc(i,k,n,i_o3) = 1.0 / rr
            o3_stom_frac(i,k) = o3_stom_frac(i,k) +                     &
     &        gsf(i,k,n) / (rr * r_stom(i,n,2))
          END IF
        END DO
!
! PAN:
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = (1.0/r_stom(i,n,3)) + (1.0/rc(i,k,n,i_pan))
            rc(i,k,n,i_pan) = 1.0 / rr
          END IF
        END DO
!
! SO2:
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = rh(i,k)
            ts = t0tile(i,k,n)
            IF (ts < 268.0) THEN
              r_cuticle = 500.0
            ELSE IF (ts >= 268.0 .AND. ts < 272.0) THEN
              r_cuticle = 200.0
            ELSE
! If CANWC > 0.25 mm assume cuticle is wet (Wyers/Erisman Atm Env 1998)
              IF (canwc(i,k,n) > 0.25) THEN
                r_cuticle = 1.0
              ELSE
                IF (rr >= 0.813) THEN
                  r_cuticle = 5.8e11 * EXP(-27.8 * rr)
                ELSE
                  r_cuticle = 2.5e4 * EXP(-6.93 * rr)
                END IF
              END IF
            END IF
            rr = (1.0/r_stom(i,n,4)) + (1.0/r_cuticle) +                &
     &        (1.0/rc(i,k,n,i_so2))
            rc(i,k,n,i_so2) = 1.0 / rr
          END IF
        END DO
!
! NH3:
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = rh(i,k)
            ts = t0tile(i,k,n)
            IF (ts < 268.0) THEN
              r_cuticle = 1000.0
            ELSE IF (ts >= 268.0 .AND. ts < 272.0) THEN
              r_cuticle = 200.0
            ELSE
! If CANWC > 0.25 mm assume cuticle is wet (Wyers/Erisman Atm Env 1998)
              IF (canwc(i,k,n) > 0.25) THEN
                r_cuticle = 10.0
              ELSE
                r_cuticle = 5.0 *                                       &
     &            LOG10(ts-zerodegc+2.0) * EXP((1.0-rr)/0.12)
              END IF
            END IF
            rr = (1.0/r_stom(i,n,5)) + (1.0/r_cuticle) +                &
     &        (1.0/rc(i,k,n,i_nh3))
            rc(i,k,n,i_nh3) = 1.0 / rr
          END IF
        END DO
!
! --------------------------------------------------------------------
! Shrub
        n = 5
!
! Shrub: H2 dry dep velocity has no dependence on soil moisture
        rr = 1.0 / h2dd_c(npft)
        DO i = 1, nlonpe
          IF (todo(i,n) .AND. microb(i,k)) THEN
            rc(i,k,n,i_h2) = rr
          END IF
        END DO
!
! Tundra. Assume tile 5, 'Shrub', is tundra if latitude is north of 60.0
! CH4 uptake flux has cubic dependence on temperature.
!
        IF (y > tundra_s_lim) THEN
          DO i = 1, nlonpe
            IF (todo(i,n) .AND. microb(i,k)) THEN
              ts = t0tile(i,k,n)
              rr = ch4dd_tun(4) + ts * (ch4dd_tun(3) +                  &
     &          ts * (ch4dd_tun(2) + ts * ch4dd_tun(1)))
              rr = rr * 3600.0 ! Convert from s-1 to h-1
              rc(i,k,n,i_ch4) = MAX(rr,0.0)
            END IF
          END DO
        ELSE
!
! CH4: Calculate an uptake flux initially. The uptake flux depends on
! soil moisture, based on results of Reay et al. (2001).
          DO i = 1, nlonpe
            IF (todo(i,n) .AND. microb(i,k)) THEN
              sm = soilmc(i,k)
              IF (sm < 0.16) THEN
                f = sm / 0.16
              ELSE IF (sm > 0.30) THEN
                f = (0.50 - sm) / 0.20
              ELSE
                f = 1.0
              END IF
              f = MAX(f,0.0)
              rc(i,k,n,i_ch4) = ch4_up_flux(n) * f
            END IF
          END DO
        END IF
!
! NO2: Calculate plant deposition terms.
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = (1.0/r_stom(i,n,1)) + (1.0/rc(i,k,n,i_no2))
            rc(i,k,n,i_no2) = 1.0 / rr
          END IF
        END DO
!
! O3:
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = (1.0/r_stom(i,n,2)) + (1.0/r_cut_o3(i,n)) +            &
     &        (1.0/rc(i,k,n,i_o3))
            rc(i,k,n,i_o3) = 1.0 / rr
            o3_stom_frac(i,k) = o3_stom_frac(i,k) +                     &
     &        gsf(i,k,n) / (rr * r_stom(i,n,2))
          END IF
        END DO
!
! PAN:
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = (1.0/r_stom(i,n,3)) + (1.0/rc(i,k,n,i_pan))
            rc(i,k,n,i_pan) = 1.0 / rr
          END IF
        END DO
!
! SO2:
!
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = rh(i,k)
            ts = t0tile(i,k,n)
            IF (ts < 268.0) THEN
              r_cuticle = 500.0
            ELSE IF (ts >= 268.0 .AND. ts < 272.0) THEN
              r_cuticle = 200.0
            ELSE
! If CANWC > 0.25 mm assume cuticle is wet (Wyers/Erisman Atm Env 1998)
              IF (canwc(i,k,n) > 0.25) THEN
                r_cuticle = 1.0
              ELSE
                IF (rr >= 0.813) THEN
                  r_cuticle = 5.8e11 * EXP(-27.8 * rr)
                ELSE
                  r_cuticle = 2.5e4 * EXP(-6.93 * rr)
                END IF
              END IF
            END IF
            rr = (1.0/r_stom(i,n,4)) + (1.0/r_cuticle) +                &
     &        (1.0/rc(i,k,n,i_so2))
            rc(i,k,n,i_so2) = 1.0 / rr
          END IF
        END DO
!
! NH3:
        DO i = 1, nlonpe
          IF (todo(i,n)) THEN
            rr = rh(i,k)
            ts = t0tile(i,k,n)
            IF (ts < 268.0) THEN
              r_cuticle = 1000.0
            ELSE IF (ts >= 268.0 .AND. ts < 272.0) THEN
              r_cuticle = 200.0
            ELSE
! If CANWC > 0.25 mm assume cuticle is wet (Wyers/Erisman Atm Env 1998)
              IF (canwc(i,k,n) > 0.25) THEN
                r_cuticle = 10.0
              ELSE
                r_cuticle = 5.0 *                                       &
     &            LOG10(ts-zerodegc+2.0) * EXP((1.0-rr)/0.12)
              END IF
            END IF
            rr = (1.0/r_stom(i,n,5)) + (1.0/r_cuticle) +                &
     &        (1.0/rc(i,k,n,i_nh3))
            rc(i,k,n,i_nh3) = 1.0 / rr
          END IF
        END DO
!
! --------------------------------------------------------------------
! Bare soil
        n = soil
! H2 and CH4:
! Tundra. Assume Bare Soil is tundra if latitude is north of 60.0 N.
! CH4 uptake flux has cubic dependence on temperature.
! H2 resistance already assigned.
        IF (y > tundra_s_lim) THEN
          DO i = 1, nlonpe
            IF (todo(i,n) .AND. microb(i,k)) THEN
              ts = t0tile(i,k,n)
              rr = ch4dd_tun(4) + ts * (ch4dd_tun(3) +                  &
     &          ts * (ch4dd_tun(2) + ts * ch4dd_tun(1)))
              rr = rr * 3600.0 ! Convert from s-1 to h-1
              rc(i,k,n,i_ch4) = MAX(rr,0.0)
            END IF
          END DO
        ELSE
!
          DO i = 1, nlonpe
            IF (todo(i,n) .AND. microb(i,k)) THEN
              sm = MAX(soilmc(i,k),0.1)
              rc(i,k,n,i_h2) = 1.0 / (h2dd_m(3) * sm + h2dd_c(3))
!
              sm = soilmc(i,k)
              IF (sm < 0.16) THEN
                f = sm / 0.16
              ELSE IF (sm > 0.30) THEN
                f = (0.50 - sm) / 0.20
              ELSE
                f = 1.0
              END IF
              f = MAX(f,0.0)
              rc(i,k,n,i_ch4) = ch4_up_flux(n) * f
            END IF
          END DO
        END IF
!
! --------------------------------------------------------------------
! Convert CH4 uptake fluxes (ug m-2 h-1) to resistance (s m-1).
        DO n = 1, npft
          DO i = 1, nlonpe
            IF (todo(i,n) .AND. microb(i,k)) THEN
              rr = rc(i,k,n,i_ch4)
              IF (rr > 0.0) THEN
                rc(i,k,n,i_ch4) = p0(i,k) * mml /                       &
     &            (rgc * t0tile(i,k,n) * rr)
              ELSE
                rc(i,k,n,i_ch4) = r_null
              END IF
            END IF
          END DO
        END DO
        n = soil
        DO i = 1, nlonpe
          IF (todo(i,n) .AND. microb(i,k)) THEN
            rr = rc(i,k,n,i_ch4)
            IF (rr > 0.0) THEN
              rc(i,k,n,i_ch4) = p0(i,k) * mml /                         &
     &          (rgc * t0tile(i,k,n) * rr)
            ELSE
              rc(i,k,n,i_ch4) = r_null
            END IF
          END IF
        END DO
!
      END DO  ! DO k = 1, rowspe
!
! Calculate resistances for SO2 and HNO3 deposition to ice, which
! depend on temperature. Ensure resistance for HNO3 does not fall
! below 10 s m-1.
!
      n = ntype
      DO k = 1, rowspe
        DO i = 1, nlonpe
          IF (gsf(i,k,n) > 0.0) THEN
            ts = t0tile(i,k,n)
            rr = hno3dd_ice(3) + ts * (hno3dd_ice(2) +                  &
     &        ts * hno3dd_ice(1))
            rc(i,k,n,i_hno3) = MAX(rr,10.0)
            rc(i,k,n,i_so2) = 1.0 / (so2dd_ice(1) + so2dd_ice(2) *      &
     &        EXP(so2dd_ice(3) * (ts-zerodegc)))
          END IF
        END DO
      END DO
!
      END SUBROUTINE SURF_DDEP_RES
#endif
