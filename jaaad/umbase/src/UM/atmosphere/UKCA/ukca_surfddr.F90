#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Current code owner: M.G. Sanderson
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
      SUBROUTINE UKCA_SURF_DD_RES(row_length, rows, sinlat,            &
        t0, p0, q, soilmc, gsf, stcon, t0tile,                         &
        lai_ft, canwc, rc, o3_stom_frac)

      USE ASAD_MOD,              ONLY: ndepd, nldepd, speci
      USE UKCA_CONSTANTS,        ONLY: rsurf, rmol
      IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "parvars.h"
#include "c_0_dg_c.h"
#include "trif.h"

!     Null resistance for deposition (1/r_null ~ 0)

      REAL, PARAMETER :: r_null = 1.0e50

      INTEGER, INTENT(IN) :: row_length, rows
      REAL, DIMENSION(row_length,rows), INTENT(IN) :: sinlat
                      ! Sine(latitude)                          
      REAL, DIMENSION(row_length,rows), INTENT(IN) :: t0  
                      ! Surface temperature (K)                 
      REAL, DIMENSION(row_length,rows), INTENT(IN) :: p0 
                      ! Surface pressure (Pa)                   
      REAL, DIMENSION(row_length,rows), INTENT(IN) :: q     
                      ! Specific humidity (kg H2O / kg air)     
      REAL, DIMENSION(row_length,rows), INTENT(IN) :: soilmc    
                      ! Soil moisture content (Fraction by volume)

      REAL, DIMENSION(row_length,rows,ntype), INTENT(IN) :: gsf
                      ! Global surface fractions                
      REAL, DIMENSION(row_length,rows,ntype), INTENT(IN) :: t0tile
                      ! Surface temperature on tiles (K)

      REAL, DIMENSION(row_length,rows,npft), INTENT(IN) ::  stcon
                      ! Stomatal conductance (m s-1)            
      REAL, DIMENSION(row_length,rows,npft), INTENT(IN) :: lai_ft
                      ! Leaf area index (m2 leaf m-2)           
      REAL, DIMENSION(row_length,rows,npft), INTENT(IN) :: canwc   
                      ! Canopy water content (mm)

!     Surface resistance on tiles (s m-1).

      REAL, DIMENSION(row_length,rows,ntype,jpdd), INTENT(OUT) :: rc
      REAL, DIMENSION(row_length,rows), INTENT(OUT) :: o3_stom_frac

!     Local variables

      LOGICAL :: first = .true.

      INTEGER :: i         ! Loop count over longitudes
      INTEGER :: k         ! Loop count over latitudes               
      INTEGER :: j         ! Loop count over species that deposit    
      INTEGER :: n         ! Loop count over tiles

      REAL :: sm            ! Soil moisture content of gridbox
      REAL :: rr            ! General temporary store for resistances.
      REAL :: ts            ! Temperature of a particular tile.
      REAL :: f             ! Factor to modify CH4 uptake fluxes
!     REAL :: r_cuticle_so2 ! Cuticular resistance for SO2
!     REAL :: r_cuticle_nh3 ! Cuticular resistance for NH3
      REAL :: mml = 1.008e5 ! Factor to convert methane flux to dry dep vel.
      REAL :: r_wet_o3 = 500.0  ! Wet soil surface resistance for ozone (s m-1)
      REAL :: cuticle_o3 = 5000.0 ! Constant for caln of O3 cuticular resistance
      REAL :: tundra_s_limit = 0.866 ! Southern limit for tundra (SIN(60))

!     MML - Used to convert methane flux in ug m-2 h-1 to dry dep vel in
!     m s-1; MML=3600*0.016*1.0E9*1.75E-6, where 0.016=RMM methane (kg),
!     1.0E9 converts ug -> kg, 1.75E-6 = assumed CH4 vmr

      REAL, DIMENSION(row_length) :: qs            ! Saturation mass mixing ratio
      REAL, DIMENSION(row_length,rows) :: rh       ! Relative humidity
      REAL, DIMENSION(row_length,npft) :: r_cut_o3 ! Cuticular Resistance for O3
      REAL, DIMENSION(row_length,rows,npft,5) :: r_stom ! Stomatal resistance
! transferred to UKCA_CONSTANTS as allocatable
!      REAL, DIMENSION(ntype,jpdd), SAVE :: rsurf   ! Standard surface resistances

!     Following arrays used to set up surface resistance array rsurf.

      REAL,PARAMETER,DIMENSION(ntype) :: zero =                        &
        (/ (r_null,i=1,ntype) /)
      REAL,PARAMETER,DIMENSION(ntype) :: rooh =                        &
        (/ 30.0, (10.0,i=1,ntype-1) /)
      REAL,PARAMETER,DIMENSION(ntype) :: hno3 =                        &
        (/ (10.0,i=1,ntype) /)
      REAL,PARAMETER,DIMENSION(ntype) :: aerosol =                     &
        (/ (r_null, i=1,npft+1), 1000.0, r_null, 20000.0 /)

!     Hydrogen - linear dependence on soil moisture (except savannah)

      REAL, DIMENSION(npft) :: h2dd_c = (/ 0.00197, 0.00197,  &
        0.00177, 1.2346, 0.0001 /)
      REAL, DIMENSION(npft) :: h2dd_m = (/ -0.00419, -0.00419,&
        -0.00414, -0.472, 0.0 /)
      REAL :: h2dd_q = 0.27 ! Quadratic term for H2 loss to savannah

!     Resistances for Tundra if different to standard value in rsurf().

      REAL, DIMENSION(5) :: r_tundra =                                 &
        (/ 1200.0, 25000.0, 800.0, 3850.0, 1100.0/)           
!           NO2      CO      O3      H2     PAN

!     CH4 uptake fluxes, in ug m-2 hr-1

      REAL, DIMENSION(ntype) :: ch4_up_flux =                          &
        (/ 39.5, 50.0, 30.0, 37.0, 27.5, 0.0, 0.0, 27.5, 0.0 /)

!     CH4 loss to tundra - Cubic polynomial fit to data. 
!     N.B. Loss flux is in units of ug(CH4) m-2 s-1

      REAL, DIMENSION(4) :: ch4dd_tun = (/ -4.757e-6, 4.0288e-3,       &
                                           -1.13592, 106.636 /)

!     HNO3 dry dep to ice; quadratic dependence

      REAL, DIMENSION(3) :: hno3dd_ice = (/ -13.57, 6841.9, -857410.6 /)

!     SO2 dry dep to snow/ice; quadratic dependence

      REAL, DIMENSION(3) :: so2dd_ice = (/ 0.0001, 0.003308, 0.1637 /)

!     Diffusion correction for stomatal conductance

      REAL, DIMENSION(5) ::                                            &
        dif = (/ 1.6, 1.6, 2.6, 1.9, 0.97 /)
!                NO2   O3  PAN  SO2  NH3

      LOGICAL, DIMENSION(row_length,rows,ntype) ::                     &
        todo     ! True if tile fraction > 0.0
      LOGICAL, DIMENSION(row_length,rows) ::                           &
        microb   ! True if T > 5 C and RH > 40%
!                  (i.e. microbes in soil are active).


! Set up standard resistance array rsurf on first call only
      IF (first) THEN
        ALLOCATE(rsurf(ntype,jpdd))
        rsurf(:,:) = r_null

! Standard surface resistances (s m-1). Values are for 9 tiles in
! order: Broadleaved trees, Needleleaf trees, C3 Grass, C4 Grass,
! Shrub, Urban, Water, Bare Soil, Ice.
!
        DO n = 1, ndepd
          SELECT CASE (speci(nldepd(n)))
            CASE ('O3        ','O3S       ')
        rsurf(:,n)=(/200.,200.,200.,200.,400.,800.,2200.,800.,         &
          2500. /)
            CASE ('NO2       ')
        rsurf(:,n)=(/225.,225.,400.,400.,600.,1200.,2600.,1200.,       &
          3500. /)
            CASE ('HNO3      ','HONO2     ')
        rsurf(:,n)=hno3
            CASE ('H2O2      ','HOOH      ','HO2H      ')
        rsurf(:,n)=hno3
            CASE ('CH3OOH    ','MeOOH     ','C2H5OOH   ','EtOOH     ', &
                  'n_C3H7OOH ','i_C3H7OOH ','n-PrOOH   ','i-PrOOH   ', &
                  'MeCOCH2OOH','ISOOH     ','MACROOH   ','MeCO3H    ')
        rsurf(:,n)=rooh
            CASE ('PAN       ','PPAN      ','MPAN      ')
        rsurf(:,n)=(/500.,500.,500.,500.,500.,r_null,12500.,           &
          500.0, 12500. /)
            CASE ('NH3       ')
        rsurf(:,n)=hno3
            CASE ('CO        ')
        rsurf(:,n)=(/3700.,7300.,4550.,1960.,4550.0,r_null,r_null,     &
          4550.0,r_null /)  ! Shrub+bare soil set to C3 grass (guess)
            CASE ('CH4       ')
        rsurf(:,n)=zero
            CASE ('H2        ')
        rsurf(:,n)=zero
            CASE ('SO2       ')
        rsurf(:,n)=(/100.,100.,150.,350.,400.,400.,10.0,700.0,         &
          r_null /)
          END SELECT
       END DO

       first = .false.
!      icounter=0
        
      END IF

      o3_stom_frac = 0.0

!     rco3global = 0.0
!     gsfglobal = 0.0


!     Set logical for surface types

      DO n = 1, ntype
        DO k = 1, rows
          DO i = 1, row_length
            todo(i,k,n) = (gsf(i,k,n) > 0.0)
          END DO
        END DO
      END DO

!     Calculate relative humidity. First, obtain saturated mixing ratio
!     (kg H2O / kg air)

      DO k = 1, rows
! DEPENDS ON: qsat
        CALL QSAT(qs,t0(:,k),p0(:,k),row_length)
!       Calculate relative humidity and set logical for microbial action.
        DO i = 1, row_length
          rh(i,k) = q(i,k) / qs(i)
          rh(i,k) = MIN(rh(i,k),1.0)
          microb(i,k) = (rh(i,k) > 0.4 .AND. t0(i,k) > 278.0)
        END DO
      END DO

!     Set surface resistances to standard values. rsurf is the 
!     resistance of the soil, rock, water etc. Set all tiles to 
!     standard values. These values will be modified below as 
!     necessary. Extra terms for vegetated tiles (stomatal, cuticular) 
!     will be added if required. Loop over all parts of array rc to 
!     ensure all of it is assigned a value.

      DO n = 1, ntype
        DO k = 1, rows
          DO i = 1, row_length
            IF (todo(i,k,n)) THEN
              DO j = 1, ndepd
                rc(i,k,n,j) = rsurf(n,j)
              END DO
            ELSE
              DO j = 1, ndepd
                rc(i,k,n,j) = r_null
              END DO
            END IF
          END DO
        END DO
      END DO

!     Calculate stomatal resistances

      DO n = 1, npft
        DO k = 1, rows
          DO i = 1, row_length
            IF (todo(i,k,n) .AND. stcon(i,k,n) > glmin(n)) THEN
              r_stom(i,k,n,1) = 1.5 * dif(1) / stcon(i,k,n) ! NO2
              r_stom(i,k,n,2) = dif(2) / stcon(i,k,n)       ! O3
              r_stom(i,k,n,3) = dif(3) / stcon(i,k,n)       ! PAN
              r_stom(i,k,n,4) = dif(4) / stcon(i,k,n)       ! SO2
              r_stom(i,k,n,5) = dif(5) / stcon(i,k,n)       ! NH3
            ELSE
              r_stom(i,k,n,1) = r_null
              r_stom(i,k,n,2) = r_null
              r_stom(i,k,n,3) = r_null
              r_stom(i,k,n,4) = r_null
              r_stom(i,k,n,5) = r_null
            END IF
          END DO
        END DO
      END DO

!     Now begin assigning specific surface resistances.

      DO j = 1, ndepd

!       O3: Change land deposition values if surface is wet; 
!       soil moisture value of 0.3 fairly arbitrary.

        IF (speci(nldepd(j)) == 'O3        ') THEN
          DO k = 1, rows
            DO i = 1, row_length
              IF (soilmc(i,k) > 0.3) THEN
                DO n = 1, npft
                  IF (todo(i,k,n)) THEN
                    rc(i,k,n,j) = r_wet_o3
                  END IF
                END DO
                IF (todo(i,k,soil)) THEN
                  rc(i,k,soil,j) = r_wet_o3
                END IF
              END IF
            END DO

!           Change values for tundra regions

            n = npft
            DO i = 1, row_length
              IF (sinlat(i,k) > tundra_s_limit) THEN
                IF (todo(i,k,n)) rc(i,k,n,j)  = r_tundra(3)
                IF (todo(i,k,soil)) rc(i,k,soil,j)  = r_tundra(3)
              END IF
            END DO

!           Cuticular resistance for ozone.

            DO n = 1, npft
              DO i = 1, row_length
                IF (todo(i,k,n) .AND. lai_ft(i,k,n) > 0.0) THEN
                  r_cut_o3(i,n) = cuticle_o3 / lai_ft(i,k,n)
                ELSE
                  r_cut_o3(i,n) = r_null
                END IF
              END DO
            END DO

!           Calculate plant deposition terms.

            DO n = 1, npft
              DO i = 1, row_length
                IF (todo(i,k,n)) THEN
                  rr = (1.0/r_stom(i,k,n,2)) +                         &
                       (1.0/r_cut_o3(i,n)) +                           &
                       (1.0/rc(i,k,n,j))
                  rc(i,k,n,j) = 1.0 / rr
                  o3_stom_frac(i,k) = o3_stom_frac(i,k) +              &
                    gsf(i,k,n) / (rr * r_stom(i,k,n,2))
                END IF
              END DO
            END DO
          END DO

!       NO2

        ELSE IF (speci(nldepd(j)) == 'NO2       ') THEN

          DO k = 1, rows

!           Change values for tundra regions

            n = npft
            DO i = 1, row_length
              IF (sinlat(i,k) > tundra_s_limit) THEN
                IF (todo(i,k,n)) rc(i,k,n,j)  = r_tundra(1)
                IF (todo(i,k,soil)) rc(i,k,soil,j)  = r_tundra(1)
              END IF
            END DO

!           Calculate plant deposition terms.

            DO n = 1, npft
              DO i = 1, row_length
                IF (todo(i,k,n)) THEN
                  rr = (1.0/r_stom(i,k,n,2)) + (1.0/rc(i,k,n,j))
                  rc(i,k,n,j) = 1.0 / rr
                END IF
              END DO
            END DO
          END DO

!       PAN

        ELSE IF (speci(nldepd(j)) == 'PAN       ') THEN

          DO k = 1, rows

!           Change values for tundra regions

            n = npft
            DO i = 1, row_length
              IF (sinlat(i,k) > tundra_s_limit) THEN
                IF (todo(i,k,n)) rc(i,k,n,j)  = r_tundra(5)
                IF (todo(i,k,soil)) rc(i,k,soil,j)  = r_tundra(5)
              END IF
            END DO

!           Calculate plant deposition terms.

            DO n = 1, npft
              DO i = 1, row_length
                IF (todo(i,k,n)) THEN
                  rr = (1.0/r_stom(i,k,n,2)) + (1.0/rc(i,k,n,j))
                  rc(i,k,n,j) = 1.0 / rr
                END IF
              END DO
            END DO
          END DO

!       H2 dry dep vel has linear dependence on soil moisture
!       Limit sm to avoid excessively high deposition velocities

        ELSE IF (speci(nldepd(j)) == 'H2        ') THEN
          DO k = 1, rows
            DO n = 1, npft-2
              DO i = 1, row_length
                IF (todo(i,k,n) .AND. microb(i,k)) THEN
                  sm = MAX(soilmc(i,k),0.1)
                  rc(i,k,n,j) = 1.0 / (h2dd_m(n) * sm + h2dd_c(n))
                END IF
              END DO
            END DO

!           C4 grass: H2 dry dep has quadratic-log dependence on 
!           soil moisture

            n = npft - 1
            DO i = 1, row_length
              IF (todo(i,k,n) .AND. microb(i,k)) THEN
                sm = LOG(MAX(soilmc(i,k),0.1))
                rr = (h2dd_c(4) + sm*(h2dd_m(4) + sm*h2dd_q)) * 1.0e-4
                IF (rr > 0.00131) rr = 0.00131 ! Conrad/Seiler Max value
                rc(i,k,n,j) = 1.0 / rr
              END IF
            END DO

!           Shrub: H2 dry dep velocity has no dependence on 
!           soil moisture

            n = npft
            rr = 1.0 / h2dd_c(n)
            DO i = 1, row_length
              IF (sinlat(i,k) > tundra_s_limit) THEN
                IF (todo(i,k,n) .AND. microb(i,k)) THEN
                  rc(i,k,n,j) = rr
                  rc(i,k,soil,j) = rr
                END IF
              END IF
            END DO

!           Bare soil

            DO i = 1, row_length
              IF (sinlat(i,k) <= tundra_s_limit) THEN
                IF (todo(i,k,n) .AND. microb(i,k)) THEN
                  sm = MAX(soilmc(i,k),0.1)
                  rc(i,k,n,j) = 1.0 / (h2dd_m(3) * sm + h2dd_c(3))
                  rc(i,k,soil,j) = 1.0 / (h2dd_m(3) * sm + h2dd_c(3))
                END IF
              END IF
            END DO
          END DO

!       CH4

        ELSE IF (speci(nldepd(j)) == 'CH4       ') THEN

!        CH4: Calculate an uptake flux initially. The uptake 
!        flux depends on soil moisture, based on results of 
!        Reay et al. (2001). Do the PFTs Broadleaf, Needleleaf, 
!        C3 and C4 grasses first.

          DO k = 1, rows
            DO i = 1, row_length
              IF (microb(i,k)) THEN
                sm = soilmc(i,k)
                IF (sm < 0.16) THEN
                  f = sm / 0.16
                ELSE IF (sm > 0.30) THEN
                  f = (0.50 - sm) / 0.20
                ELSE
                  f = 1.0
                END IF
                f = MAX(f,0.0)
                IF (todo(i,k,1)) rc(i,k,1,j) = ch4_up_flux(1) * f
                IF (todo(i,k,2)) rc(i,k,2,j) = ch4_up_flux(2) * f
                IF (todo(i,k,3)) rc(i,k,3,j) = ch4_up_flux(3) * f
                IF (todo(i,k,4)) rc(i,k,4,j) = ch4_up_flux(4) * f
              END IF
            END DO
          END DO

!         Now do shrub and bare soil, assumed to be tundra if 
!         latitude > 60N; Ctherwise, calculate an uptake flux 
!         initially. The uptake flux depends on soil moisture, 
!         based on results of Reay et al. (2001).

          n = npft
          DO k = 1, rows
            DO i = 1, row_length
              IF (microb(i,k)) THEN
                IF (sinlat(i,k) > tundra_s_limit) THEN
                  ts = t0tile(i,k,n)
                  rr = ch4dd_tun(4) + ts * (ch4dd_tun(3) +             &
                    ts * (ch4dd_tun(2) + ts * ch4dd_tun(1)))
                  rr = rr * 3600.0 ! Convert from s-1 to h-1
                  IF (todo(i,k,n)) rc(i,k,n,j) = MAX(rr,0.0)
                  IF (todo(i,k,soil)) rc(i,k,soil,j) = MAX(rr,0.0)
                ELSE
                  sm = soilmc(i,k)
                  IF (sm < 0.16) THEN
                    f = sm / 0.16
                  ELSE IF (sm > 0.30) THEN
                    f = (0.50 - sm) / 0.20
                  ELSE
                    f = 1.0
                  END IF
                  f = MAX(f,0.0)
                  IF (todo(i,k,n)) rc(i,k,n,j) = ch4_up_flux(n) * f
                  IF (todo(i,k,soil)) rc(i,k,soil,j) =                 &
                    ch4_up_flux(soil) * f
                END IF
              END IF
            END DO
          END DO

!         Convert CH4 uptake fluxes (ug m-2 h-1) to 
!         resistance (s m-1).

          DO k = 1, rows
            DO n = 1, npft
              DO i = 1, row_length
                IF (todo(i,k,n) .AND. microb(i,k)) THEN
                  rr = rc(i,k,n,j)
                  IF (rr > 0.0) THEN
                    rc(i,k,n,j) = p0(i,k) * mml /                      &
                      (rmol * t0tile(i,k,n) * rr)
                  ELSE
                    rc(i,k,n,j) = r_null
                  END IF
                END IF
              END DO
            END DO
            n = soil
            DO i = 1, row_length
              IF (todo(i,k,n) .AND. microb(i,k)) THEN
                rr = rc(i,k,n,j)
                IF (rr > 0.0) THEN
                  rc(i,k,n,j) = p0(i,k) * mml /                        &
                    (rmol * t0tile(i,k,n) * rr)
                ELSE
                  rc(i,k,n,j) = r_null
                END IF
              END IF
            END DO
          END DO
 
        ELSE IF (speci(nldepd(j)) == 'CO        ') THEN

!         Only assign values for CO if microbes are active.

          n = npft
          DO k = 1, rows
            DO i = 1, row_length
              IF (sinlat(i,k) > tundra_s_limit .AND. microb(i,k)) THEN
                IF (todo(i,k,n)) THEN
                  rc(i,k,n,j)  = r_tundra(2)       ! CO
                END IF
                IF (todo(i,k,soil)) THEN
                  rc(i,k,soil,j)  = r_tundra(2)    ! CO
                END IF
              END IF
            END DO
          END DO

!       HONO2

!       Calculate resistances for HONO2 deposition to ice, which
!       depend on temperature. Ensure resistance for HONO2 does not fall
!       below 10 s m-1.
 
        ELSE IF (speci(nldepd(j)) == 'HNO3      ' .OR.                 &
          speci(nldepd(j)) == 'HONO2     ') THEN

          n = ntype
          DO k = 1, rows
            DO i = 1, row_length
              IF (todo(i,k,n)) THEN

!               Limit temperature to a minimum of 252K. Curve used 
!               only used data between 255K and 273K.

                ts = MAX(t0tile(i,k,n), 252.0)
                rr = hno3dd_ice(3) + ts * (hno3dd_ice(2) +             &
                     ts * hno3dd_ice(1))
                rc(i,k,n,j) = MAX(rr,10.0)
              END IF
            END DO
          END DO

        END IF ! End of IF (speci == species name)

      END DO  ! End of DO j = 1, ndepd

      END SUBROUTINE UKCA_SURF_DD_RES
#endif
