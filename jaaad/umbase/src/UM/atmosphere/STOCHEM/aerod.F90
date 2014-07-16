#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE AEROD(t0,p0,hf,bl,tauu,tauv,canht,gsf,                 &
     &  z0,ra,rq,so4_vd)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods :
!-     Calculates aerodynamic and quasi-laminar resistances.
!-     Returns RA, RQ, SO4_VD
!
! Current Code Owner: M.G. Sanderson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   4.0   08/10/97   Created (as a function). W.J. Collins.
!   6.1   03/12/04   Extensively rewritten for new dry deposition
!                    scheme. Now a subroutine. M.G. Sanderson
!   6.2   28/03/06   Minor changes for vn6.2  M.G. Sanderson
!
!----------------------------------------------------------------------
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!----------------------------------------------------------------------
!
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) ::                     &
     &  t0                                                              &
                             ! Surface temperature (K)
     & ,p0                                                              &
                             ! Surface pressure (Pa)
     & ,hf                                                              &
                             ! Surface heat flux (W m-2)
     & ,bl                                                              &
                             ! Boundary layer height (m)
     & ,tauu                                                            &
                             ! Wind stress in u-direction (N m-2)
     & ,tauv                 ! Wind stress in v-direction (N m-2)
      REAL, DIMENSION(nlonpe,nlatpe,npft), INTENT(IN) ::                &
     &  canht                ! Canopy height of vegetation (m)
      REAL, DIMENSION(nlonpe,nlatpe,ntype), INTENT(IN) ::               &
     &  gsf                  ! Surface tile fractions
      REAL, DIMENSION(nlonpe,nlatpe,ntype), INTENT(INOUT) ::            &
     &  z0                   ! Roughness length on tiles (m)
      REAL, DIMENSION(nlonpe,nlatpe,ntype), INTENT(OUT) ::              &
     &  ra                   ! Aerodynamic resistance (s m-1)
      REAL, DIMENSION(nlonpe,nlatpe,nc), INTENT(OUT) ::                 &
     &  rq                   ! Quasi-laminar resistance (s m-1)
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(OUT) ::                    &
     &  so4_vd               ! Sulphate aerosol dep velocity (m s-1)
!
#include "c_vkman.h"
!
      LOGICAL :: first = .true.
!
      INTEGER :: i, j, k, n     ! Loop counts
!
      REAL ::                                                           &
     &  cp = 1004.67                                                    &
                             ! Heat capacity of air under const p (J K-1
     & ,pr = 0.72                                                       &
                             ! Prandtl number of air
     & ,vair296 = 1830.0e-8                                             &
                             ! Dynamic viscosity of air at 296 K (kg m-1
     & ,twoth = 2.0 / 3.0                                               &
                             ! 2/3
     & ,zref = 50.0          ! Reference height for dry deposition (m).
!
      REAL ::                                                           &
     &  b0                                                              &
                             ! Temporary store
     & ,b1                                                              &
                             ! Temporary store
     & ,bl_l                                                            &
                             ! Boundary layer height / Monin-Obukhov len
     & ,sc_pr                ! Ratio of Schmidt and Prandtl numbers
!
      REAL, SAVE :: reus
      REAL, DIMENSION(nc), SAVE :: d0
!
      REAL, DIMENSION(nlonpe,nlatpe) ::                                 &
     &  ustar                                                           &
                  ! Friction velocity (m s-1)
     & ,l                                                               &
                  ! Monin-Obukhov length (m)
     & ,rho_air                                                         &
                  ! Density of air at surface (kg m-3)
     & ,kva       ! Kinematic velocity of air (m2 s-1)

      REAL, DIMENSION(nlonpe,nlatpe,npft) ::                            &
     &  d         ! Zero-plane displacement (m)

      REAL, DIMENSION(nlonpe,nlatpe,ntype) ::                           &
     &  z         ! Reference height (m)

      REAL, DIMENSION(nlonpe,nlatpe,ntype) ::                           &
     &  psi       ! Businger function (dimensionless)
!
! Assign diffusion coefficients, units m2 s-1. Set to -1 unless species
! dry deposits. If no value found in literature, D0 calculated using:
!  D(X) = D(H2O) * SQRT[RMM(H2O)/RMM(X)], where X is the species in
! question and D(H2O) = 2.08 x 10^-5 m2 s-1 (Marrero & Mason, J Phys
! Chem Ref Dat, 1972).
!
      IF (first) THEN
        d0 = -1.0
        d0(i_no2) = 1.4e-5
        d0(i_o3) = 1.4e-5
        d0(i_pan) = 0.31e-5
        d0(i_so2) = 1.2e-5
        d0(i_nh3) = 2.08e-5
        d0(i_co) = 1.86e-5
        d0(i_ch4) = 5.74e-5
        d0(i_h2) = 6.7e-5
        d0(i_hno3) = 1.2e-5
        d0(i_h2o2) = 1.46e-5
        d0(i_ch3ooh) = 1.27e-5
        d0(i_c3h7ooh) = 1.01e-5
        d0(i_c2h5ooh) = 1.12e-5
        d0(i_c4h9ooh) = 0.93e-5
        d0(i_isopooh) = 0.81e-5
        d0(i_mvkooh) = 0.81e-5
!
! Values for aerosols are set to avoid problems in calculation
! of RQ below. RQ for aerosols explicitly set below.
        d0(i_sa) = 1.0e-5
        d0(i_msa) = 1.0e-5
        d0(i_ammsul) = 1.0e-5
        d0(i_naer) = 1.0e-5
        d0(i_orgnit) = 1.0e-5
!
! Set up other constants
!
        reus = -cp / (vkman * g)
!
        first = .false.
!
      END IF
!
! Calculate resistance parameters
!   1. RA
!
      ra = r_null
      DO k = 1, rowspe
        DO i = 1, nlonpe
!
! Calculate density of air in kg m-3.
          rho_air(i,k) = p0(i,k) * mair / (rgc * t0(i,k))
!
! Calculate friction velocity USTAR; ensure is non-zero
          ustar(i,k) = sqrt(sqrt(tauu(i,k)**2 + tauv(i,k)**2) /         &
     &      rho_air(i,k))
          IF (ustar(i,k) == 0.0) ustar(i,k) = 1.0e-2
        END DO
      END DO
!
! Set zero-plane displacement to 0.7 of canopy height (Smith et al., 200
!
      DO n = 1, npft
        DO k = 1, rowspe
          DO i = 1, nlonpe
            IF (canht(i,k,n) > 0.0) THEN
              d(i,k,n) = canht(i,k,n) * 0.70
            ELSE
              d(i,k,n) = 0.0
            END IF
          END DO
        END DO
      END DO
!
! Ensure undefined parts of Z0 are set to 0.001 to avoid any
! floating-point exceptions. Initialise z to reference height.
!
      DO n = 1, ntype
        DO k = 1, rowspe
          DO i = 1, nlonpe
            IF (z0(i,k,n) <= 0.0) z0(i,k,n) = 0.001
            z(i,k,n) = zref
          END DO
        END DO
      END DO
!
! Set Z to height above surface minus zero plane displacement
! for vegetated tiles.
!
      DO n = 1, npft
        DO k = 1, rowspe
          DO i = 1, nlonpe
            z(i,k,n) = zref - d(i,k,n)
          END DO
        END DO
      END DO
!
! Calculate kinematic viscosity of air, KVA; = dynamic viscosity / densi
! Formula from Kaye and Laby, 1952, p.43.
!
      DO k = 1, rowspe
        DO i = 1, nlonpe
          kva(i,k) = (vair296 - 4.83e-8 * (296.0-t0(i,k))) /            &
     &      rho_air(i,k)
!
! Calculate roughness length over oceans using Charnock formula in
! form given by Hummelshoj et al., 1992.
!
          IF (gsf(i,k,nwater) > 0.0) THEN
            z0(i,k,nwater) = (kva(i,k) / (9.1 * ustar(i,k))) +          &
     &        0.016 * ustar(i,k) * ustar(i,k) / g
          END IF
        END DO
      END DO
!
! Calculate Monin-Obukhov length L
!
      DO k = 1, rowspe
        DO i = 1, nlonpe
          IF (hf(i,k) /= 0.0) THEN
! Stable or unstable b.l.
            l(i,k) = reus * t0(i,k) * rho_air(i,k) * (ustar(i,k)**3) /  &
     &        hf(i,k)
          ELSE
! Neutral b.l.
            l(i,k) = 10000.0
          END IF
        END DO
      END DO
      IF (ANY(l == 0.0)) WRITE(6,*) 'L=0:'
!
! Calculate Businger functions (PSI(ZETA))
!
      DO k = 1, rowspe
        DO i = 1, nlonpe
          IF (l(i,k) > 0.0) THEN
! Stable b.l.
            DO n = 1, ntype
              psi(i,k,n) = -5.0 * (z(i,k,n) - z0(i,k,n)) / l(i,k)
            END DO
          ELSE
! Unstable b.l.
            DO n = 1, ntype
              b1 = (1.0 - 16.0 * (z(i,k,n) / l(i,k)))**0.5
              b0 = (1.0 - 16.0 * (z0(i,k,n) / l(i,k)))**0.5
              psi(i,k,n) = 2.0 * LOG((1.0+b1) / (1.0+b0))
            END DO
          END IF
        END DO
      END DO
!
! Calculate aerodynamic resistance (RA) for all species
!
      DO n = 1, ntype
        DO k = 1, rowspe
          DO i = 1, nlonpe
            ra(i,k,n) = (LOG(z(i,k,n) / z0(i,k,n)) - psi(i,k,n)) /      &
     &        (vkman * ustar(i,k))
          END DO
        END DO
      END DO
!
!   2. RQ
!
! Calculate Schmidt number divided by Prandtl number for each species
! Schmidt number = Kinematic viscosity of air divided by molecular diffu
! Calculate quasi-laminar boundary layer resistance (Hicks et al., 1987)
!
      rq = r_null
!
      DO j = 1, n_spec_deposit
        n = dep_spec(j)
        DO k = 1, rowspe
          DO i = 1, nlonpe
            sc_pr = kva(i,k) / (pr * d0(n))
            rq(i,k,n) = (sc_pr**twoth) / (vkman * ustar(i,k))
          END DO
        END DO
      END DO
!
! Set RQ for aerosols to 1.0
      DO k = 1, rowspe
        DO i = 1, nlonpe
          rq(i,k,i_sa) = 1.0
          rq(i,k,i_msa) = 1.0
          rq(i,k,i_orgnit) = 1.0
          rq(i,k,i_naer) = 1.0
          rq(i,k,i_ammsul) = 1.0
        END DO
      END DO
!
!   3. SO4_VD
! Calculate surface deposition velocity term for sulphate particles
! using simple parameterisation of Wesely in the form given by
! Zhang et al. (2001) Atmos Environ, 35, 549-560.
! These values will be used for all aerosols (SA, AMMSUL, ORGNIT, NAER,
! and MSA).
!
      so4_vd = 1.0 / r_null
!
      DO k = 1, rowspe
        DO i = 1, nlonpe
          so4_vd(i,k) = 0.0
          bl_l = bl(i,k) / l(i,k)
          IF (l(i,k) >= 0.0) so4_vd(i,k) = 0.002 * ustar(i,k)
          IF (bl_l < -30.0) so4_vd(i,k) = 0.0009 * ustar(i,k) *         &
     &      ((-bl_l)**twoth)
          IF (bl_l >= -30.0 .AND. l(i,k) < 0.0) so4_vd(i,k) =           &
     &      0.002 * ustar(i,k) *  (1.0 + (-300.0/l(i,k))**twoth)
        END DO
      END DO
!
!     WRITE(6,*) 'Aerodynamic Resistance'
!     WRITE(6,*) ra
!     WRITE(6,*) '-----------------------------------------------'
!     WRITE(6,*) 'Quasi-laminar Resistance'
!     WRITE(6,*) rq
!     WRITE(6,*) '-----------------------------------------------'
!     WRITE(6,*) 'Sulphate deposition velocity
!     WRITE(6,*) so4_dv
!     WRITE(6,*) '-----------------------------------------------'
!
      END SUBROUTINE AEROD
#endif
