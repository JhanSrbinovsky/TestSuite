#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT*******************************
!
! Description:
!  Main driver routine for chemistry
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!   Called from UKCA_MAIN1.
!
! Current code owner: Colin Johnson/Olaf Morgenstern
!                     Fiona O'Connor
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!------------------------------------------------------------------
!
      SUBROUTINE UKCA_CHEMISTRY_CTL(i_month, i_day_number, i_hour,     &
                      i_minute, secs_per_step,                         &
                      p_levelsda, ntracers,                            &
                      ndiags,                                          &
                      nfluxdiags,                                      &
                      sinlat,                                          &
                      coslat,                                          &
                      true_longitude,                                  &
                      pres, temp, q,                                   &
                      qcf, qcl, rh,                                    &
                      p_layer_boundaries,                              &
                      r_theta_levels,                                  &
                      cos_zenith_angle,                                &
                      tracer,                                          &
                      user_diags,                                      &
                      mflux_diags,                                     &
                      t_surf, dzl, z0m, u_s,                           &
                      drain, crain,                                    &
                      cloud_frac,                                      &
                      fastj_dj,                                        &
                      volume, mass,                                    &
                      land_points, land_index,                         &
                      tile_pts, tile_index, tile_frac,                 &
                      zbl, surf_hf, seaice_frac, stcon,                &
                      soilmc_lp, fland, laift_lp, canhtft_lp,          &
                      z0tile_lp, t0tile_lp, canwctile_lp,              &
                      nbimol, nphotf, n_opl, nabund,                   &
                      abundance, phot, bimol, termol,                  &
                      hetero, drydep, wetdep,                          &
                      ozonebud,                                        &
                      pv_at_theta,                                     &
                      theta,                                           &
                      um_ozone3d,                                      &
                      delso2_wet_h2o2,                                 &
                      delso2_wet_o3,                                   &
                      delso2_dry_oh,                                   &
                      delso2_drydep,                                   &
                      delso2_wetdep,                                   &
                      so4_sa                                           &
                      )

      USE ASAD_MOD
      USE UKCA_D1_DEFS
      USE UKCA_CSPECIES
      USE UKCA_tropopause
      USE UKCA_strat_update
      USE      UKCA_phot2d,          ONLY: UKCA_PHOTIN, UKCA_CURVE,    &
                                           UKCA_INPR2D,                &
                                           nolev, nlphot, ntphot, pjin
      IMPLICIT NONE

#include "parvars.h"
#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "typsize.h"
#include "c_v_m.h"
#include "c_pi.h"
#include "c_g.h"
#include "c_sulchm.h"
#include "csubmodl.h"
#include "typsts.h"

        INTEGER, INTENT(IN) :: p_levelsda        ! no. of model levels
        INTEGER, INTENT(IN) :: ntracers          ! no. of tracers
        INTEGER, INTENT(IN) :: ndiags            ! no. of diagnostics
        INTEGER, INTENT(IN) :: nfluxdiags        ! no. of flux diagnostics
        INTEGER, INTENT(IN) :: i_month           ! month
        INTEGER, INTENT(IN) :: i_day_number      ! day
        INTEGER, INTENT(IN) :: i_hour            ! hour
        INTEGER, INTENT(IN) :: i_minute          ! minute
        INTEGER, INTENT(IN) :: nbimol            ! No bimolecular fluxes
        INTEGER, INTENT(IN) :: nphotf            ! No photolytic fluxes
        INTEGER, INTENT(IN) :: n_opl             ! No O3 budget terms
        INTEGER, INTENT(IN) :: nabund            ! No abundance terms

!       Variables for interactive dry deposition scheme
        INTEGER, INTENT(IN) :: land_points
        INTEGER, INTENT(IN) :: land_index(land_points)
        INTEGER, INTENT(IN) :: tile_pts(ntype)
        INTEGER, INTENT(IN) :: tile_index(land_points,ntype)

        REAL, INTENT(IN) :: secs_per_step                      ! time step
        REAL, INTENT(IN) :: sinlat(row_length, rows)           ! sin(latitude)
        REAL, INTENT(IN) :: coslat(row_length, rows)           ! cos(latitude)
        REAL, INTENT(IN) :: true_longitude(row_length,rows)    ! longitude
        REAL, INTENT(IN) :: pres(row_length,rows,p_levelsda)   ! pressure
        REAL, INTENT(IN) :: p_layer_boundaries(row_length,rows,    &
                                               0:p_levelsda)   ! pressure
        REAL, INTENT(IN) :: r_theta_levels(row_length,rows,0:p_levelsda)
        REAL, INTENT(IN) :: temp(row_length,rows,p_levelsda)   ! actual temp
        REAL, INTENT(IN) :: dzl(row_length, rows, bl_levels)   ! thickness
        REAL, INTENT(IN) :: u_s(row_length, rows)              ! ustar
        REAL, INTENT(IN) :: z0m(row_length, rows)              ! roughness
        REAL, INTENT(IN) :: t_surf(row_length, rows)           ! surface temp
        REAL, INTENT(IN) :: drain(row_length,rows, p_levelsda) ! 3-D LS rain
        REAL, INTENT(IN) :: crain(row_length,rows, p_levelsda) ! 3-D convec
        REAL, INTENT(IN) :: cos_zenith_angle(row_length, rows) ! cosine of ZA
        REAL, INTENT(IN) :: volume(row_length,rows,p_levelsda) ! cell vol.
        REAL, INTENT(IN) :: mass(row_length, rows, p_levelsda) ! cell mass
        REAL, INTENT(IN) :: pv_at_theta(row_length, rows, p_levelsda) ! PV
        REAL, INTENT(IN) :: theta(row_length, rows, p_levelsda)    ! theta
        REAL, INTENT(IN) :: um_ozone3d(row_length, rows, p_levelsda) ! O3
        REAL, INTENT(IN) :: qcf(row_length, rows, wet_levels)  ! qcf
        REAL, INTENT(IN) :: qcl(row_length, rows, wet_levels)  ! qcl
        REAL, INTENT(IN) :: rh(row_length, rows, wet_levels)    ! RH frac
        REAL, INTENT(IN) :: cloud_frac(row_length, rows, wet_levels)
        REAL, INTENT(IN) :: so4_sa(row_length,rows,p_levelsda)  ! aerosol
!                                                        surface area

!       Variables for interactive dry deposition scheme

        REAL, INTENT(IN) :: tile_frac(land_points,ntype)
        REAL, INTENT(IN) :: zbl(row_length,rows)
        REAL, INTENT(IN) :: surf_hf(row_length,rows)
        REAL, INTENT(IN) :: seaice_frac(row_length,rows)
        REAL, INTENT(IN) :: stcon(row_length,rows,npft)
        REAL, INTENT(IN) :: soilmc_lp(land_points)
        REAL, INTENT(IN) :: fland(land_points)
        REAL, INTENT(IN) :: laift_lp(land_points,npft)
        REAL, INTENT(IN) :: canhtft_lp(land_points,npft)
        REAL, INTENT(IN) :: z0tile_lp(land_points,ntype)
        REAL, INTENT(IN) :: t0tile_lp(land_points,ntype)
        REAL, INTENT(IN) :: canwctile_lp(land_points,ntype)

        REAL, INTENT(INOUT) :: fastj_dj(row_length,rows,p_levelsda,jppj)
        REAL, INTENT(INOUT) :: q(row_length,rows,p_levelsda)   ! water vapo
        REAL, INTENT(INOUT) :: tracer(row_length,rows,              &
                                      p_levelsda,ntracers)     ! tracer MMR

        REAL, INTENT(OUT) :: user_diags(row_length,rows,            &
                                        p_levelsda,ndiags)     ! chem diags
        REAL, INTENT(OUT) :: mflux_diags(row_length,rows,           &
                                        p_levelsda,nfluxdiags) ! flux diags
! budget variables:
      REAL, INTENT(OUT) :: abundance(row_length,rows,p_levelsda,nabund)
      REAL, INTENT(OUT) :: phot(row_length,rows,p_levelsda,nphotf)
      REAL, INTENT(OUT) :: bimol(row_length,rows,p_levelsda,nbimol)
      REAL, INTENT(OUT) :: termol(row_length,rows,p_levelsda)
      REAL, INTENT(OUT) :: hetero(row_length,rows,p_levelsda)
      REAL, INTENT(OUT) :: drydep(row_length,rows,p_levelsda)
      REAL, INTENT(OUT) :: wetdep(row_length,rows,p_levelsda)
      REAL, INTENT(OUT) :: ozonebud(row_length,rows,p_levelsda,n_opl)

! SO2 increments
      REAL, INTENT(INOUT) :: delSO2_wet_H2O2(row_length,rows,p_levelsda)
      REAL, INTENT(INOUT) :: delSO2_wet_O3(row_length, rows, p_levelsda)
      REAL, INTENT(INOUT) :: delSO2_dry_OH(row_length, rows, p_levelsda)
      REAL, INTENT(INOUT) :: delSO2_drydep(row_length, rows, p_levelsda)
      REAL, INTENT(INOUT) :: delSO2_wetdep(row_length, rows, p_levelsda)


!       Local variables

        INTEGER, SAVE :: nr            ! no of rxns for BE
        INTEGER, SAVE :: n_be_calls    ! no of call to BE solver
        INTEGER, SAVE :: first_row
        INTEGER, SAVE :: first_column

        INTEGER :: i             ! Loop variable
        INTEGER :: j             ! loop variable
        INTEGER :: js            ! loop variable
        INTEGER :: k             ! loop variable
        INTEGER :: l             ! loop variable
        INTEGER :: ll            ! loop variable
        INTEGER :: m             ! loop variable
        INTEGER :: sp            ! loop variable
        INTEGER :: n_pnts        ! no. of pts in 2D passed to CDRIVE
        INTEGER :: icnt

        CHARACTER(LEN=72)                 :: cmessage  ! Error message

        REAL, PARAMETER :: fxb = 23.45 * Pi_Over_180 ! tropic of capricorn
        REAL, PARAMETER :: fxc = 24.0/ pi
        REAL, PARAMETER :: DTS = 300.0               ! B. Euler timestep

!       Pressure level above which top boundary conditions are applied
        REAL, PARAMETER :: p_above = 7000.           ! Pa

        REAL :: tgmt               ! GMT time (decimal represent
        REAL :: declin             ! declination
        REAL :: total_water        ! Total water
        REAL :: const              ! constant
        REAL, PARAMETER :: limit = 1200.   ! seconds. For timesteps
                                           ! greater than this we halve
                                           ! the chemical timestep
                                           ! (for IMPACT).

! SO2 increments in molecules/cm^3
      REAL :: delta_SO2_wetox_H2O2_total(theta_field_size)
      REAL :: delta_SO2_wetox_O3_total(theta_field_size)
      REAL :: delta_SO2_dryox_OH_total(theta_field_size)

        REAL, ALLOCATABLE :: BE_rc(:,:)
        REAL :: zftr(theta_field_size,jpctr) ! 1-D array of tracers
        REAL :: zp  (theta_field_size) ! 1-D pressure
        REAL :: zt  (theta_field_size) ! 1-D temperature
        REAL :: zclw(theta_field_size)       ! 1-D cloud liquid water
        REAL :: zfcloud(theta_field_size)    ! 1-D cloud fraction
        REAL :: cdot(theta_field_size,jpctr) ! 1-D chem. tendency
        REAl :: pq(row_length,rows)         ! 2-D water vapour mmr
        REAL :: pjinda(rows, ntphot, jppj)  ! PJIN at one level
        REAL :: zprt(row_length, rows, jppj)   ! 2-D photolysis rates
        REAL :: BE_zprt(theta_field_size, jppj) ! 2-D photolysis rates
        REAL :: wks(theta_field_size, jppj)
        REAL :: tloc  (row_length, rows)     ! local time
        REAL :: daylen(row_length, rows)     ! local daylength
        REAL :: cs_hour_ang(row_length, rows)! cosine hour angle
        REAL :: tanlat(row_length, rows)     ! tangens of latitude
        REAL :: zdryrt(row_length, rows, jpdd)  ! dry dep rate
        REAL :: zdryrt2(theta_field_size, jpdd) ! dry dep rate
        REAL :: zwetrt(row_length, rows, p_levelsda, jpdw) ! wet dep rate
        REAL :: zwetrt2(theta_field_size, jpdw)            ! wet dep rat
        REAL :: zfrdiss2(theta_field_size,jpdw,2)    ! dissolved fraction
        REAL :: zfrdiss(row_length, rows, p_levelsda, jpdw, 2)
        REAL :: kp_nh(row_length, rows, p_levelsda) ! Dissociation const
        REAL :: kp_nh2(theta_field_size)            ! Dissociation const
        REAL :: nlev_in_bl(row_length, rows)        ! No levs in bl
        REAL :: nlev_in_bl2(theta_field_size)       ! No levs in bl
        REAL :: ozonecol(row_length, rows, p_levelsda) ! for strat chem
        REAL :: BE_tnd(theta_field_size)            ! total no density
        REAL :: BE_h2o(theta_field_size)            ! water vapour concn
        REAL :: BE_o2(theta_field_size)             ! oxygen concn
        REAL :: BE_vol(theta_field_size)            ! gridbox volume
        REAL :: BE_wetrt(theta_field_size,jpspec)   ! wet dep rates (s-1)
        REAL :: BE_dryrt(theta_field_size,jpspec)   ! dry dep rates (s-1)
        REAL :: BE_deprt(theta_field_size,jpspec)   ! dep rates (s-1)
        REAL :: BE_frdiss(theta_field_size,jpspec,2) ! dissolved fraction
        REAL :: BE_y (theta_field_size,jpspec)
        REAL :: k_dms(theta_field_size,5)           ! dms rate coeffs
        REAL :: sflux(theta_field_size,nfluxdiags)  ! rxn fluxes summed
        REAL :: pr2d(nolev)             ! 2D model level pressures
        REAL :: pr2dj(nlphot)           ! 2D photolysis level pres
        REAL :: sc_3d(row_length,rows,p_levelsda,jpspec) ! 3-D sc
        REAL, SAVE :: first_lat
        REAL, SAVE :: dellat
        REAL, SAVE :: first_lon
        REAL, SAVE :: dellon

        LOGICAL, SAVE :: firstcall = .true.

! Variables for heterogeneous chemistry
        REAL, ALLOCATABLE :: shno3_3d(:,:,:)
        LOGICAL :: stratflag(theta_field_size)

        COMMON /STREQ/ PR2D

        n_pnts = rows * row_length
        nr = nr_therm + nr_phot

        IF (firstcall) then

!         Check that theta_field_size >= n_pnts

          IF (theta_field_size /= n_pnts) THEN
            cmessage='theta_field_size not equal to n_pnts'
          WRITE(6,*) cmessage
! DEPENDS ON: ereport
            CALL EREPORT('UKCA_CHEMISTRY_CTL',n_pnts,cmessage)
          ENDIF

!         Check that nfluxdiags >= nr+jpdd+jpdw

          IF (L_ukca_BEflux .AND. nfluxdiags < nr+jpdd+jpdw) THEN
            cmessage='nfluxdiags needs to be increased'
            WRITE(6,*) cmessage
! DEPENDS ON: ereport
            CALL EREPORT('UKCA_CHEMISTRY_CTL',nr+jpdd+jpdw,cmessage)
        ENDIF

! Determine where the domain is in latitude
          first_lat = ASIN(sinlat(1,1))
          dellat = ASIN(sinlat(1,2)) - first_lat
          first_row = INT((first_lat + 0.5*pi) / dellat + 0.00001) + 1

! Determine where the domain is in longitude
          first_lon = true_longitude(1,2)
          dellon = true_longitude(2,2) - first_lon
          first_column = INT(first_lon / dellon + 0.00001) + 1

!         Update ASAD timestep variables

          dtime = secs_per_step
          kcdt = INT(secs_per_step)

!         Backward Euler timestep variables

          n_be_calls = INT(dtime/dts)

!         Check whether water vapour is advective tracer. Then,
!         check whether UM and ASAD advective tracers correspond
!         to each other.

          IF ((L_ukca_advh2o .AND. (ntracers+1 /= jpctr))               &
          .OR.                                                          &
         (.NOT.(L_ukca_advh2o) .AND. (ntracers /= jpctr)))              &
          THEN
            cmessage=' ntracers inconsistent with ASAD.'
            WRITE (6,*) cmessage,ntracers,jpctr
! DEPENDS ON: ereport
            CALL EREPORT('UKCA_CHEMISTRY_CTL',1,cmessage)
          ENDIF

          IF ((L_ukca_advh2o) .AND. (n_h2o == 0)) THEN
            cmessage='No tracer for advected water vapour'
! DEPENDS ON: ereport
            CALL ereport('UKCA_CHEMISTRY_CTL',4,cmessage)
          ENDIF

          firstcall = .false.

        END IF  ! of initialization of chemistry subroutine

!       Calculate local time as function of longitude

        tgmt = real(i_hour) + real(i_minute)/60.0                      &
                            + secs_per_step * 0.5 / 3600.0
!        IF (tgmt < 0.) tgmt = tgmt + 24.

        tloc = tgmt + 24.0 * true_longitude/pi/2.0
        WHERE (tloc > 24.0) tloc = tloc - 24.0

!       Calculate Declination Angle and Daylength for each row for
!       current day of the year.
!       Ensure COS of HOUR ANGLE does not exceed + or - 1, and set DAY
!       LENGTH of 1st & last rows to be the same as that of adjacent rows
!       to avoid possible problems at the poles (tan(90)=infinity).

        daylen = 0.0
        DO i=1,rows
          DO j=1,row_length
            IF (ABS(coslat(j,i)) < 1e-10) THEN
              IF (sinlat(j,i) >= 0.0) THEN
                tanlat(j,i) = 1.0e20
              ELSE
                tanlat(j,i) = -1.0e20
              ENDIF
            ELSE
              tanlat(j,i)=sinlat(j,i)/coslat(j,i)
            END IF
          END DO
        END DO

        declin = fxb * SIN(pi_over_180*(266.0+i_day_number))
        cs_hour_ang = -tanlat * TAN(declin)
        WHERE (cs_hour_ang < -1.0) cs_hour_ang = -1.0
        WHERE (cs_hour_ang >  1.0) cs_hour_ang =  1.0
        daylen = fxc * acos(cs_hour_ang)

!       Call routine to calculate dry deposition rates.

        zdryrt = 0.0
        IF (ndepd /= 0) THEN

          IF (L_ukca_intdd) THEN           ! Call interactive dry dep

! DEPENDS ON: ukca_drydep_ctl
            CALL UKCA_DRYDEP_CTL(row_length, rows, bl_levels,          &
              land_points, land_index, tile_pts, tile_index,           &
              secs_per_step, sinlat, tile_frac, t_surf,                &
              p_layer_boundaries(:,:,0), dzl, zbl, surf_hf, u_s,       &
              q, stcon, soilmc_lp, fland, seaice_frac, laift_lp,       &
              canhtft_lp, z0tile_lp, t0tile_lp, canwctile_lp,          &
              nlev_in_bl, zdryrt)

          ELSE                             ! Call prescribed dry dep

! DEPENDS ON: ukca_ddeprt
            CALL UKCA_DDEPRT(daylen, tloc, n_pnts, dzl, bl_levels,     &
                                     z0m, u_s, t_surf,                 &
                                     sinlat, i_month,                  &
                                     1, n_pnts,                        &
                                     zdryrt)

          ENDIF
        ENDIF

!       Call routine to calculate wet deposition rates.

        zwetrt = 0.0
        IF (ndepw /= 0) THEN

! DEPENDS ON: ukca_wdeprt
          CALL UKCA_WDEPRT(drain,crain,n_pnts,p_LEVELSDA,              &
                               temp,SINLAT,                            &
                               SECS_PER_STEP,                          &
                               1,n_pnts, ZWETRT)
        ENDIF

! Calculate dissolved fraction
      IF (L_ukca_aerchem) THEN
! DEPENDS ON: ukca_fracdiss
        CALL UKCA_FRACDISS(n_pnts,p_levelsda,wet_levels, temp,rh,qcl,  &
                           zfrdiss, kp_nh)
      ENDIF

!       Call routine to read in 2D photolysis rates once per day and
!       interplate to model latitude and levels.

        IF (L_ukca_phot2d .AND. i_hour == 0 .AND. i_minute == 0) THEN
! DEPENDS ON: ukca_phot2d
          CALL UKCA_PHOTIN(i_day_number, row_length, rows,             &
                           p_levelsda, n_pnts,                         &
                           first_row, global_row_length, global_rows,  &
                           reshape(sinlat,(/theta_field_size/)),       &
                           pres, jppj)
        ELSE IF (L_ukca_fastj) THEN
! DEPENDS ON: ukca_phot2d
          CALL UKCA_INPR2D(pr2d,pr2dj)
        ENDIF

!       Calculate ozone column for stratospheric photolysis

        IF (L_ukca_strat .OR. L_ukca_stratcfc .OR. L_ukca_strattrop)    &
          THEN
          IF (n_ro3 > 0) THEN
            CALL UKCA_CALC_OZONECOL(p_levelsda, rows, row_length,       &
                           p_layer_boundaries, pres,                    &
                           tracer(:,:,:,n_ox)/c_o3,                     &
                           ozonecol)
          ENDIF

! Calculate total chlorine and total bromine before chemistry

          IF (nn_cl > 0)  THEN
            CALL conserve(row_length, rows, p_levelsda, ntracers,       &
                    advt, tracer, pres, um_ozone3d, .TRUE.)
          ENDIF
        ENDIF

! if heterogeneous chemistry is selected, allocate solid HNO3 array
      IF (L_ukca_het_psc) THEN
        ALLOCATE(shno3_3d(row_length, rows, p_levelsda))
        shno3_3d = 0.
      ENDIF

!       Initialize budget variables

      IF (l_ukca_budget2) THEN
        abundance = 0.
        bimol = 0.
        termol = 0.
        hetero = 0.
        drydep = 0.
        wetdep = 0.
        phot = 0.
      END IF

      IF (L_ukca_BEflux) THEN
        DO l = 1, nfluxdiags
          DO k = 1, p_levelsda
            DO i = 1, rows
              DO j = 1, row_length
                mflux_diags(j,i,k,l) = 0.0
              END DO
            END DO
          END DO
        END DO
      ENDIF

      DO k=1,p_levelsda

! Copy water vapour and ice field into 1-D arrays
        DO i=1,rows
          DO j=1,row_length
            l = (i-1) * row_length + j
            pq(j,i) = q(j,i,k)/c_h2o
            IF (L_ukca_het_psc) sph2o(l) = qcf(j,i,k)/c_h2o
          END DO
        END DO

!       Put pressure, temperature and tracer mmr into 1-D arrays
!       for use in ASAD chemical solver

        zclw(:) = 0.0
        zfcloud(:) = 0.0

        DO i=1,rows
          DO j=1,row_length
            l = (i-1) * row_length + j
            zp  (l  ) = pres(j,i,k)
            zt  (l  ) = temp(j,i,k)
            IF(L_UKCA_aerchem) THEN
              IF(k <= wet_levels) THEN
                zclw(l) = qcl(j,i,k)
                zfcloud(l) = cloud_frac(j,i,k)
              ENDIF
            ENDIF
        ENDDO
      ENDDO

!       Convert mmr into vmr for tracers
!       If water is advective tracer then the water vapour
!       variable needs to be nudged into the array of tracers.
!       Otherwise this is not necessary. MAKE SURE H2O IS NEITHER
!       THE FIRST NOR THE LAST ADVECTIVE TRACER.

        IF (L_ukca_advh2o) THEN
          DO js = 1, jpctr
            DO i=1,rows
              DO j=1,row_length
                l = (i-1) * row_length + j
                IF (js < n_h2o) THEN
                  zftr(l,js) = tracer(j,i,k,js) / c_species(js)
                ELSE IF (js == n_h2o) THEN
                  zftr(l,n_h2o) = pq(j,i)
                ELSE IF (js > n_h2o) THEN
                  zftr(l,js) = tracer(j,i,k,js-1) / c_species(js)
                END IF
              END DO
            END DO
          END DO
        ELSE
          DO js=1,jpctr
            DO i=1,rows
              DO j=1,row_length
                l = (i-1) * row_length + j
                zftr(l,js) = tracer(j,i,k,js) / c_species(js)
              END DO
            END DO
          END DO
        END IF

!       Interplate 2D photolysis rates as function of longitude
!       per model step.

        IF (L_ukca_phot2d) THEN

!         Interplate 2D photolysis rates as function of longitude
!         per model step.

          pjinda = pjin(1:rows,k,:,:)
! DEPENDS ON: ukca_phot2d
          CALL UKCA_CURVE(pjinda, reshape(tloc,(/n_pnts/)),            &
                          reshape(daylen,(/n_pnts/)), n_pnts,          &
                          p_levelsda, rows, row_length, wks)
          zprt = reshape(wks,(/row_length, rows, jppj/))

!         Call STRAT_PHOTOL to update photolysis with stratospheric
!         rates from SLIMCAT (Lary scheme)

          IF ((L_ukca_strat .OR. L_ukca_strattrop .OR.                  &
               L_ukca_stratcfc) .AND. (k > 54))  THEN
            CALL ukca_strat_photol(pres(:,:,k), temp(:,:,k),            &
                            ozonecol(:,:,k), cos_zenith_angle,          &
                            jppj, zprt)
          ENDIF
        ELSE IF (L_ukca_fastj) THEN
          zprt(:,:,:) = fastj_dj(:,:,k,:)
        ENDIF    ! L_ukca_phot2d

!       Call ASAD routines to do chemistry integration
!       In lowest levels choose half the dynamical timestep for
!       chemistry. If dynamical timestep > 20 min, use half and
!       quarter of dynamical timestep for chemistry.

        IF (.NOT.(L_ukca_trop .OR. L_ukca_aerchem .OR.                  &
                  L_ukca_tropisop)) THEN

! fill stratospheric flag indicator and SO4 surface area
          DO i=1,rows
            DO j=1,row_length
              l = (i-1) * row_length + j
              stratflag(l) = (.NOT. L_troposphere(i,j,k))
              za(l) = so4_sa(j,i,k)
            ENDDO
          ENDDO

          IF (method == 3) THEN
! timestep for Newton-Raphson solver: 1 hour. Note that in this case
! UKCA_CHEMISTRY_CTL is only called every 2nd/3rd dynamical timestep
! (for a 30 / 20 minutes stepsize), and dtime should equal 1 h.
! In case of non-convergence Newton-Raphson contains automatic stepsize
! refinement.
            ncsteps = 1
            cdt = dtime
          ELSEIF (method == 1) THEN
! use about 15 or 10 minutes, depending on dynamical timestep
            IF (secs_per_step < limit) THEN
              cdt = dtime
              ncsteps = 1
            ELSE
              cdt = 0.5*dtime
              ncsteps = 2
            END IF
          END IF

! DEPENDS ON: asad_cdrive
          CALL ASAD_CDRIVE(cdot,zftr,zp,zt,pq,k,zdryrt,zwetrt,          &
                      zprt,n_pnts)

          IF (L_ukca_het_psc) THEN
! Save MMR of NAT PSC particles into 3-D array for PSC sedimentation.
! Note that sphno3 is NAT in number density of HNO3.
            IF (ANY(sphno3(1:theta_field_size) > 0.)) THEN
              DO i=1,rows
                DO j=1,row_length
                  l =  (i-1) * row_length + j
                  shno3_3d(j,i,k) = sphno3(l) / tnd(l) * c_hono2
                END DO
              END DO
            ELSE
              shno3_3d(:,:,k) = 0.
            END IF
          END IF

          IF (L_ukca_budget2) THEN  ! do budget calculation

! DEPENDS ON: ukca_budget2
            CALL ukca_budget2( rows, row_length, n_pnts,                &
                        k, p_levelsda,                                  &
                        secs_per_step,                                  &
                        zdryrt, zwetrt,                                 &
                        zprt,                                           &
                        nabund, nbimol, nphotf,                         &
                        abundance,                                      &
                        bimol,                                          &
                        termol,                                         &
                        hetero, drydep,                                 &
                        wetdep,                                         &
                        phot, volume, pres,                             &
                        n_opl, n_o3, ozonebud,                          &
                        global_rows, global_row_length, first_row,      &
                        first_column, .TRUE.                            &
                     )
          END IF          ! L_ukca_budget2

! Rescale bromine and chlorine tracers to guarantee conservation of total
! chlorine, bromine, and hydrogen over timestep. Only makes sense if at least
! chlorine chemistry is present.
          IF (nn_cl > 0) THEN
            CALL conserve(row_length, rows, p_levelsda, ntracers,       &
                          advt, tracer, pres, um_ozone3d,.false.)
          END IF

!         Write species concentration array in vmr.

          DO i=1,n_pnts
            sc(i,k,:) = y(i,:)/tnd(i)
          END DO

!         Calculate mmrs from vmrs

          DO i=1,rows
            DO j=1,row_length
              l = (i-1) * row_length + j

!             Bring results back from vmr to mmr. Treat water
!             vapour separately if necessary.

              IF (L_ukca_h2o_feedback) THEN
                tracer(J,I,K,:) = zftr(l,:    )*c_species
                q(j,i,k  ) = zftr(l,n_h2o)*c_h2o
              ELSE
                IF (n_h2o < jpctr)                                      &
                  tracer(j,i,k,n_h2o+1:jpctr) = zftr(l,n_h2o+1:jpctr) * &
                                               c_species(n_h2o+1:jpctr)
                IF (n_h2o > 1)                                          &
                  tracer(j,i,k,1:n_h2o-1) = zftr(l,1:n_h2o-1) *         &
                                               c_species(1:n_h2o-1)
              END IF

! Copy SS and family member VMRs into diagnostic fields for family
! chemistry
              IF (L_ukca_family) THEN
                user_diags(j,i,k,1) = y(l,nn_oh)/tnd(l)    ! OH vmr
                user_diags(j,i,k,2) = y(l,nn_ho2)/tnd(l)   ! HO2 vmr
                user_diags(j,i,k,3) = y(l,nn_no)/tnd(l)    ! NO vmr
                user_diags(j,i,k,4) = y(l,nn_no2)/tnd(l)   ! NO2 vmr
                user_diags(j,i,k,5) = y(l,nn_o1d)/tnd(l)   ! O1D vmr
                user_diags(j,i,k,6) = y(l,nn_meoo)/tnd(l)  ! MeOO vmr
                IF (nn_meco3 > 0)                                       &
                  user_diags(j,i,k,7) = user_diags(j,i,k,7)             &
                    + y(l,nn_meco3)/tnd(l)
                IF (nn_etoo > 0)                                        &
                  user_diags(j,i,k,7) = user_diags(j,i,k,7)             &
                    + y(l,nn_etoo)/tnd(l)
                IF (nn_etco3 > 0)                                       &
                  user_diags(j,i,k,7) = user_diags(j,i,k,7)             &
                    + y(l,nn_etco3)/tnd(l)
                IF (nn_nproo > 0)                                       &
                  user_diags(j,i,k,7) = user_diags(j,i,k,7)             &
                    + y(l,nn_nproo)/tnd(l)
               IF (nn_nprooh > 0)                                      &
                  user_diags(j,i,k,7) = user_diags(j,i,k,7)             &
                    + y(l,nn_nprooh)/tnd(l)
                user_diags(j,i,k,8) = y(l,nn_o3)/tnd(l)       ! O3 vmr
                user_diags(j,i,k,9) = y(l,nn_o3p)/tnd(l)      ! O(3P) vmr
                IF (nn_n > 0)                                           &
                  user_diags(j,i,k,10) = y(l,nn_n)/tnd(l)     ! N vmr
                IF (nn_cl > 0)                                          &
                  user_diags(j,i,k,11) = y(l,nn_cl)/tnd(l)    ! Cl vmr
                IF (nn_clo > 0)                                         &
                  user_diags(j,i,k,12) = y(l,nn_clo)/tnd(l)   ! ClO vmr
                IF (nn_cl2o2 > 0)                                       &
                  user_diags(j,i,k,13) = y(l,nn_cl2o2)/tnd(l) ! Cl2O2 vmr
                IF (nn_meo > 0)                                         &
                  user_diags(j,i,k,14) = y(l,nn_meo)/tnd(l)   ! MeO vmr
                IF (nn_br > 0)                                          &
                  user_diags(j,i,k,15) = y(l,nn_br)/tnd(l)    ! Br vmr
                IF (nn_bro > 0)                                         &
                  user_diags(j,i,k,16) = y(l,nn_bro)/tnd(l)   ! BrO vmr

! copy ozone into RO3 tracer
                IF (n_ro3 > 0)                                          &
                  tracer(j,i,k,n_ro3) = sc(l,k,nn_o3) * c_o3
              END IF
            END DO
          END DO

        ELSE  ! Backward Euler with non-families

!         Calculate total number  density, o2, h2o, and tracer
!         concentrations for Backward Euler solver

          const = Boltzmann*1.0e6
          DO i=1,rows
            DO j=1,row_length
              l = (i-1) * row_length + j
              nlev_in_bl2(l) = nlev_in_bl(j,i)
!              BE_tnd(l) = zp(l) / (const * zt(l) )
              BE_tnd(l) = zp(l) / (1.3806e-23 * 1.0e6 * zt(l) )
              BE_o2(l)  = 0.2095 * BE_tnd(l)
              BE_h2o(l) = pq(j,i) * BE_tnd(l)
              BE_vol(l) = volume(j,i,k) * 1.0e6 ! m3 -> cm3
              zwetrt2(l,:) = zwetrt(j,i,k,:)
              zdryrt2(l,:) = zdryrt(j,i,:)
              BE_zprt(l,:) = zprt(j,i,:)
              zftr(l,:) = zftr(l,:) * BE_tnd(l)
              IF (L_ukca_aerchem) THEN
                zfrdiss2(l,:,:) = zfrdiss(j,i,k,:,:)
                kp_nh2(l) = kp_nh(j,i,k)
              ENDIF
            END DO
          END DO

!         Assign wet and dry deposition rates to species

! DEPENDS ON: ukca_be_wetdep
          CALL UKCA_BE_WETDEP(n_pnts, zwetrt2, be_wetrt)

! DEPENDS ON: ukca_be_drydep
          CALL UKCA_BE_DRYDEP(k, n_pnts, nlev_in_bl2, zdryrt2, be_dryrt)

! Assign fractional dissociation
          IF(L_ukca_aerchem) THEN
! DEPENDS ON: ukca_be_wetdep
            CALL UKCA_BE_WETDEP(n_pnts, zfrdiss2(1:n_pnts,:,1),         &
                                  be_frdiss(1:n_pnts,:,1))
! DEPENDS ON: ukca_be_wetdep
            CALL UKCA_BE_WETDEP(n_pnts, zfrdiss2(1:n_pnts,:,2),         &
                                  be_frdiss(1:n_pnts,:,2))
          ENDIF

!         Calculate reaction rate coefficients

          ALLOCATE(BE_rc(theta_field_size,nr_therm))
! DEPENDS ON: ukca_chemco
          CALL UKCA_CHEMCO(nr_therm, n_pnts, zt(1:n_pnts),              &
                           BE_tnd(1:n_pnts), BE_h2o(1:n_pnts),          &
                           BE_o2(1:n_pnts), zclw(1:n_pnts),             &
                           zfcloud(1:n_pnts), BE_frdiss(1:n_pnts,:,:),  &
                           k_dms(1:n_pnts,:), BE_rc(1:n_pnts,:))

!         Assign tracer concentrations to species concentrations

          BE_y(1:n_pnts,1:jpspec) = 0.0
        DO i = 1,jpspec
          DO j = 1,jpctr
            IF (speci(i) == advt(j)) THEN
                BE_y(1:n_pnts,i) = zftr(1:n_pnts,j)
            EXIT
            ENDIF
          ENDDO
          ENDDO

          IF (L_ukca_aerchem) THEN
            delta_SO2_wetox_H2O2_total(:) = 0.0
            delta_SO2_dryox_OH_total(:)   = 0.0
            delta_SO2_wetox_O3_total(:)   = 0.0
          ELSE
            delSO2_wet_H2O2(:,:,:) = 0.0
            delSO2_wet_O3(:,:,:)   = 0.0
            delSO2_dry_OH(:,:,:)   = 0.0
            delSO2_drydep(:,:,:)   = 0.0
            delSO2_wetdep(:,:,:)   = 0.0
          ENDIF

!         Call Backward Euler solver
!         N.B. Emissions already added, via call to TR_MIX from
!         UKCA_EMISSION_CTL

          IF (L_ukca_aerchem) THEN
! DEPENDS ON: ukca_deriv_aero
            CALL UKCA_DERIV_AERO(nr_therm, nr_phot, n_be_calls,         &
                          n_pnts, nfluxdiags, BE_rc(1:n_pnts,:),        &
                          BE_wetrt(1:n_pnts,:), BE_dryrt(1:n_pnts,:),   &
                          BE_zprt(1:n_pnts,:), k_dms(1:n_pnts,:),       &
                          BE_h2o(1:n_pnts), BE_tnd(1:n_pnts),           &
                          BE_o2(1:n_pnts), BE_vol(1:n_pnts),            &
                          dts, ldepd, ldepw, L_ukca_BEflux,             &
                          BE_y(1:n_pnts,:), sflux(1:n_pnts,:),          &
                          delta_SO2_wetox_H2O2_total(1:n_pnts),         &
                          delta_SO2_wetox_O3_total(1:n_pnts),           &
                          delta_SO2_dryox_OH_total(1:n_pnts) )
          ELSE

! DEPENDS ON: ukca_deriv
          CALL UKCA_DERIV(nr, n_be_calls, n_pnts, nfluxdiags,          &
                       BE_rc(1:n_pnts,:), BE_wetrt(1:n_pnts,:),        &
                       BE_dryrt(1:n_pnts,:), BE_zprt(1:n_pnts,:),      &
                       BE_h2o(1:n_pnts), BE_tnd(1:n_pnts),             &
                       BE_o2(1:n_pnts), BE_vol(1:n_pnts),              &
                       dts, ldepd, ldepw, L_ukca_BEflux,               &
                       BE_y(1:n_pnts,:), sflux(1:n_pnts,:))
          ENDIF
          DEALLOCATE(BE_rc)

!         Assign species concentrations to tracer concentrations

        DO j = 1,jpctr
          DO i = 1,jpspec
            IF (advt(j) == speci(i)) THEN
              zftr(1:n_pnts,j) = BE_y(1:n_pnts,i)
              EXIT
            ENDIF
          ENDDO
        ENDDO

!         Write to species array in vmr and convert tracer
!         concentrations back to mmr


          DO js = 1, jpspec
            DO l = 1, n_pnts
              sc(l,k,js) = BE_y(l,js) / BE_tnd(l)
            END DO
          END DO

!         Convert tracers back to mmr

          DO js = 1, jpctr
            DO i = 1, rows
              DO j = 1, row_length
                l = (i-1) * row_length + j
                tracer(j,i,k,js) = zftr(l,js)                          &
                                 * c_species(js) / BE_tnd(l)
              END DO
            END DO
          END DO

         DO i = 1, rows
            DO j = 1, row_length
              l = (i-1) * row_length + j
              user_diags(j,i,k,1)  = BE_y(l,nn_oh)*c_oh/BE_tnd(l)
              user_diags(j,i,k,2)  = BE_y(l,nn_ho2)*c_ho2/BE_tnd(l)
              user_diags(j,i,k,3)  = BE_y(l,nn_o3p)*c_o3p/BE_tnd(l)
              user_diags(j,i,k,4)  = BE_y(l,nn_o1d)*c_o1d/BE_tnd(l)
              user_diags(j,i,k,5)  = BE_y(l,nn_meoo)*c_meoo/BE_tnd(l)
              user_diags(j,i,k,6)  = BE_y(l,nn_etoo)*c_etoo/BE_tnd(l)
              user_diags(j,i,k,7)  = BE_y(l,nn_meco3)*c_meco3/BE_tnd(l)
              user_diags(j,i,k,8)  = BE_y(l,nn_nproo)*c_proo/BE_tnd(l)
              user_diags(j,i,k,9)  = BE_y(l,nn_iproo)*c_proo/BE_tnd(l)
              user_diags(j,i,k,10) = BE_y(l,nn_etco3)*c_etco3/BE_tnd(l)
              user_diags(j,i,k,11) = BE_y(l,nn_mecoch2oo)              &
                                   *c_mecoch2oo/BE_tnd(l)

              IF (L_troposphere(j,i,k)) THEN
              user_diags(j,i,k,16) = zftr(l,n_ch4)*volume(j,i,k)     &
                           *1.0e6/avogadro   ! trop ch4 burden in moles
                user_diags(j,i,k,17) = zftr(l,n_o3)*volume(j,i,k)      &
                           *1.0e6/avogadro   ! trop o3 burden in moles
              ELSE
                user_diags(j,i,k,16) = 0.0
                user_diags(j,i,k,17) = 0.0
              ENDIF
              user_diags(j,i,k,18) = theta_trop(j,i)
              user_diags(j,i,k,19) = pv_trop(j,i)
              user_diags(j,i,k,20) = p_tropopause(j,i)

              IF (L_ukca_aerchem) THEN
                delSO2_wet_H2O2(j,i,k)=delta_SO2_wetox_H2O2_total(l)
                delSO2_wet_O3(j,i,k) = delta_SO2_wetox_O3_total(l)
                delSO2_dry_OH(j,i,k) = delta_SO2_dryox_OH_total(l)
                delSO2_drydep(j,i,k) = BE_dryrt(l,nn_so2)               &
                        * BE_y(l,nn_so2) * DTS
                delSO2_wetdep(j,i,k) = BE_wetrt(l,nn_so2)               &
                        * BE_y(l,nn_so2) * DTS
              ENDIF
            END DO
          END DO

!         Set mean fluxes above tropopause to zero

          IF (L_ukca_BEflux) THEN
          DO js = 1, nfluxdiags
            DO i=1,rows
              DO j=1,row_length
                IF (L_troposphere(j,i,k)) THEN  ! in troposphere
                  l = (i-1) * row_length + j
                  mflux_diags(j,i,k,js) = sflux(l,js)               &
                                          /dtime  ! flux in moles/s
                  ELSE
                  mflux_diags(j,i,k,js) = 0.0   ! in stratosphere
                END IF
              END DO
            END DO
          END DO
        ENDIF     ! L_ukca_BEflux

        END IF
      END DO            ! level loop (k)

      IF (L_ukca_budget2)  THEN

! Do final budget calculations for ASAD, summing over processors
! DEPENDS ON: ukca_budget2
        CALL ukca_budget2( rows, row_length, n_pnts, k, p_levelsda,     &
                      secs_per_step,                                    &
                      zdryrt, zwetrt, zprt,                             &
                      nabund, nbimol, nphotf,                           &
                      abundance, bimol,                                 &
                      termol, hetero,                                   &
                      drydep, wetdep,                                   &
                      phot, volume, pres,                               &
                      n_opl, n_o3, ozonebud,                            &
                      global_rows, global_row_length, first_row,        &
                      first_column, .FALSE. )
      ENDIF

      IF (L_ukca_strat .OR. L_ukca_stratcfc .OR. L_ukca_strattrop)      &
         THEN
        IF (L_ukca_het_psc) THEN
! Do NAT PSC sedimentation

! take NAT out of gasphase again
          tracer(:,:,:,n_hono2) = tracer(:,:,:,n_hono2) - shno3_3d

! DEPENDS ON: ukca_sediment
          CALL ukca_sediment(rows, row_length, p_levelsda, shno3_3d,    &
                   qcf, r_theta_levels, mass, secs_per_step)

! add solid-phase HNO3 back to gasphase HNO3
          tracer(:,:,:,n_hono2) = tracer(:,:,:,n_hono2) + shno3_3d
        END IF

! Again calculate ozone column, this time from actual post-chemistry
! ozone, for diagnostic purposes
        IF (ndiags > 19) THEN
          IF (n_ro3 > 0) THEN
            CALL UKCA_CALC_OZONECOL(p_levelsda, rows, row_length,       &
                               p_layer_boundaries, pres,                &
                               tracer(:,:,:,n_ro3)/c_o3,                &
                               user_diags(:,:,:,20))
          ELSE
            CALL UKCA_CALC_OZONECOL(p_levelsda, rows, row_length,       &
                               p_layer_boundaries, pres,                &
                               tracer(:,:,:,n_ox)/c_o3,                 &
                               user_diags(:,:,:,20))
          END IF   ! ro3

! Copy NAT MMR into user_diagostics
          IF ((L_ukca_het_psc) .AND. (ndiags > 18))                     &
            user_diags(:,:,:,19)=shno3_3d
        END IF

      ELSE     ! tropospheric chemistry
! Call routine to overwrite O3, CH4 and NOy species once per day
! above level defined by p_above. Only for tropospheric chemistry
!
! DEPENDS ON: ukca_stratf
        CALL UKCA_STRATF(i_day_number, row_length,rows, p_levelsda,    &
                      n_pnts, first_row, global_row_length,            &
                      global_rows, jpctr, sinlat, pres,                &
                      um_ozone3d, p_above,                             &
                      tracer(1:row_length,1:rows,1:p_levelsda,1:jpctr))

        DO js = 1, jpspec
          DO k = 1, p_levelsda
            DO i = 1, rows
              DO j = 1, row_length
                l = (i-1) * row_length + j
                sc_3d(j,i,k,js) = sc(l,k,js)
              END DO
            END DO
          END DO
        END DO

      ENDIF     ! L_ukca_strat etc

      IF (L_ukca_het_psc) DEALLOCATE(shno3_3d)

      END SUBROUTINE UKCA_CHEMISTRY_CTL
!
#endif
