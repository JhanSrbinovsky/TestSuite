#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE EMCALC_NATURAL_VOC(month,day,totpar,dirpar,            &
     &  gsf,t0tile,lai_ft,isopem,terpem)
!
! Purpose:
! Calculates emission of isoprene and monoterpenes, given PAR flux densi
!   and temperature. Based on model developed by
!   Guenther et al. JGR 1995, V.100(D5) 8873-8892.
!
! Method:
!
! History:
! Version   Date     Comment
! ---      -----     -------
! 5.5      02/08/04  Created. M.G. Sanderson
! 5.5      24/09/04  Added tile fraction test. M.G. Sanderson
! 6.2      01/03/06  Improved vectorisation. M.G. Sanderson
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM, ONLY : mcarb
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!-----------------------------------------------------------------------
!
! Arguments with Intent IN. ie: Input variables.

      INTEGER, INTENT(IN) ::                                            &
     &  month                 ! Month
      REAL, INTENT(IN) ::                                               &
     &  day                   ! Day of month + hour as fraction
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) ::                     &
     &  totpar,                                                         &
                              ! Total surface PAR flux, W m-2
     &  dirpar                ! Direct component of PAR flux, W m-2
      REAL, DIMENSION(nlonpe,nlatpe,ntype), INTENT(IN) ::               &
     &  t0tile,                                                         &
                              ! Surface temperature of each tile
     &  gsf                   ! Global surface type fractions
      REAL, DIMENSION(nlonpe,nlatpe,npft), INTENT(IN) ::                &
     &  lai_ft                ! Leaf area index m2 m-2
!
! Arguments with Intent OUT

      REAL, DIMENSION(nlnpe,nlpe), INTENT(OUT) ::                       &
     &  isopem,                                                         &
                              ! Isoprene emission, vmr grid square-1 s-1
     &  terpem                ! Terpene emission, vmr grid square-1 s-1

! Local variables
      INTEGER :: i, k, n, jj
      REAL, DIMENSION(nlonpe) :: cza
      REAL, DIMENSION(nlonpe,nlatpe,2) :: voce
      REAL, DIMENSION(nlnpe,nlpe,2) :: vocse
      REAL, DIMENSION(nlonpe) :: pbeam, pscat, fi, ft
!
      REAL :: lai_sun, lai_shade, q1, qsun, qshade
      REAL :: al, icf, tcf, ct, denom, v, tt
      REAL :: ifac, tfac, clsun, clshade
      REAL :: x, y
      REAL :: lw, gamma, ss
!
! Parameters for G95 model (Guenther et al., 1995).
      REAL :: alpha = 0.0027, cli = 1.066
      REAL :: ct1 = 9.5e4, ct2 = 2.3e5, ts = 303.0, tmax = 314.0
      REAL :: beta = 0.09
      REAL :: cos_mlsa = 0.5 ! Cos(Mean Leaf-Sun Angle); assumed 60 deg.
      REAL :: w2photon = 0.219
      REAL :: red_fac = 0.6875  ! Reduce emissions to 550 Tg/yr
!
! Isoprene emission factors, in ug C g-1 h-1, from Guenther et al. (1995
!
      REAL, DIMENSION(npft) :: iefac = (/                               &
     & 24.0,                                                            &
               ! Broadleaf trees
     &  8.0,                                                            &
               ! Needleleaf trees
     & 16.0,                                                            &
               ! C3 Grass
     & 16.0,                                                            &
               ! C4 grass
     & 16.0 /) ! Shrub
!
! Terpene emission factors, in ug C g-1 h-1, from Guenther et al. (1995)
!
      REAL, DIMENSION(npft) :: tefac = (/                               &
     &  0.4,                                                            &
               ! Broadleaf trees
     &  2.4,                                                            &
               ! Needleleaf trees
     &  0.8,                                                            &
               ! C3 Grass
     &  0.8,                                                            &
               ! C4 grass
     &  0.8 /) ! Shrub
!
! Specific leaf weights in g m-2, from Guenther et al. (1995)
!
      REAL, DIMENSION(npft) :: slw = (/                                 &
     & 125.0,                                                           &
               ! Broadleaf trees
     & 150.0,                                                           &
               ! Needleleaf trees
     & 125.0,                                                           &
               ! C3 grass
     & 125.0,                                                           &
               ! C4 grass
     & 125.0/) ! Shrub
!
! Factors to convert isoprene and terpene emissions from ug C m-2 h-1
! to molecules m-2 s-1. Factor of 5 & 10 needed as isoprene has 5 carbon
! atoms (C5H8), and terpenes 10 (C10H16).
!
      icf = na * 1.0e-9 / (5.0 * mcarb * 3600.0)
      tcf = na * 1.0e-9 / (10.0 * mcarb * 3600.0)
!
      voce = 0.0
      vocse = 0.0
      isopem = 0.0
      terpem = 0.0
      pbeam = 0.0
      pscat = 0.0
!
! DEPENDS ON: secs
      ss = SECS(day,month,1)
      DO k = 1, rowspe
        jj = k + lobound - 1
        y = latm_half(jj) - 90.0
!
! Calculate cosines of solar zenith angles along longitude row
        DO i = 1, nlonpe
          x = longm_half(i+lnbound-1)
! DEPENDS ON: zen
          cza(i) = COS(ZEN(ss,y,x))
!
! Calculate fluxes of direct and diffuse PAR in umol photons m-2 s-1
          pbeam(i) = dirpar(i,k) / w2photon
          pscat(i) = (totpar(i,k) - dirpar(i,k)) / w2photon
        END DO
!
        fi = 0.0
        ft = 0.0
!
        DO n = 1, npft
          lw = slw(n)
          ifac = iefac(n)
          tfac = tefac(n)
          DO i = 1, nlonpe
            tt = t0tile(i,k,n)
            al = lai_ft(i,k,n)
!
! Calculate area of sunlit and shaded leaves
            IF (cza(i) > 0.0 .AND. al > 0.0 .AND.                       &
     &        gsf(i,k,n) > 0.0) THEN
              lai_sun = (1.0 - EXP(-0.5*al / cza(i))) * cza(i)          &
     &          / cos_mlsa
              lai_shade = al - lai_sun
!
! Calculate flux density of PAR on a sunlit leaf, QSUN, and on a
! shaded leaf, QSHADE
              q1 = 0.07*pbeam(i) * (1.1 - 0.1*al) * EXP(-cza(i))
              qshade = pscat(i) * exp(-0.5 * (al**0.7)) + q1
              qsun = pbeam(i) * (cos_mlsa / cza(i)) + qshade
!
              v = alpha * qsun
              clsun = cli * v / sqrt(1.0 + v*v)
!
              v = alpha * qshade
              clshade = cli * v / sqrt(1.0 + v*v)
!
              denom = rgc * ts * tt
              ct = exp((ct1*(tt-ts) / denom)) /                         &
     &          (1.0 + exp((ct2 * (tt - tmax) / denom)))
!
! Sum isoprene emissions from sunlit and shaded leaves
!
              fi(i) = fi(i) + lw * ct * ifac * gsf(i,k,n) *             &
     &          ((lai_sun * clsun) + (lai_shade * clshade))
!
            END IF
!
! Calculate temperature factor for monoterpene emissions
!
            gamma = exp(beta * (tt - ts))
!
! Sum terpene emissions from leaves
!
            ft(i) = ft(i) + al * lw * gamma * tfac * gsf(i,k,n)
          END DO
        END DO
!
! Store emissions in units of ugC m-2 h-1
!
        DO i = 1, nlonpe
          voce(i,k,1) = fi(i)
          voce(i,k,2) = ft(i)
        END DO
!
      END DO
!
! Convert emission data from UM grid to STOCHEM grid
!
! DEPENDS ON: met2data
      CALL MET2DATA(vocse,voce,2,2,.FALSE.)
!
! Convert emissions from ugC m-2 h-1 to vmr grid cell-1 s-1
!
      DO k = 1, nlpe
        DO i = 1, nlnpe
          isopem(i,k) = vocse(i,k,1) * icf * red_fac *                  &
     &      stochem_grid_area(k+ltdat-1) / lmolec
          terpem(i,k) = vocse(i,k,2) * tcf * red_fac *                  &
     &      stochem_grid_area(k+ltdat-1) / lmolec
        END DO
      END DO
!
      END SUBROUTINE EMCALC_NATURAL_VOC
#endif
