#if defined(A06_4A)
! *********************************************************************
! (c) CROWN COPYRIGHT 2001, Met Office, All Rights Reserved.
! Please refer to file $UMDIR/vn$VN/copyright.txt for further details
! *********************************************************************
!
!+ Gravity Wave (Ultra-Simple Spectral Parametrization) Scheme.
! Subroutine Interface:
      SUBROUTINE GW_USSP(LEVELS, MODEL_DOMAIN, ROWS, NROWS,             &
     &  OFF_X, OFF_Y, HALO_I, HALO_J, ROW_LENGTH,                       &
     &  R_RHO_LEVELS, R_THETA_LEVELS, P_LAYER_BOUNDARIES,               &
     &  R_U, R_V, L_USSP_OPAQUE,                                        &
     &  SIN_THETA_LONGITUDE, SIN_THETA_LATITUDE,                        &
     &  THETA, RHO, TIMESTEP, U, V, AT_EXTREMITY,                       &
     &  GWSPEC_EFLUX,GWSPEC_SFLUX,GWSPEC_WFLUX,GWSPEC_NFLUX,            &
     &  GWSPEC_EWACC,GWSPEC_NSACC,ETA_THETA_LEVELS,                     &
     &  GWSPEC_EFLUX_ON, GWSPEC_SFLUX_ON, GWSPEC_WFLUX_ON,              &
     &  GWSPEC_NFLUX_ON, GWSPEC_EWACC_ON, GWSPEC_NSACC_ON)
!
! purpose:    This subroutine calculates the vertical momentum flux
!             divergence due to gravity waves as parametrised by the
!             Warner and McIntyre Ultra Simple Spectral gravity wave
!             Parametrization adapted for use in the UM.
!
! Current code owner: A.C. Bushell (from July 2003).
!
! History:
! Version   Date    Comment
! ----   -------   --------
! 5.3   08/10/2001 This deck introduced.  A.A. Scaife
! 5.4   05/09/02 Add diagnostics for spectral (non-orographic) gravity
!                wave forcing.       Adam Scaife
!
! 5.5   16/07/2003 Allow parametrized gravity waves out of top of model.
!                                         A.C. Bushell
!
! 6.0   29/12/2003 Optimization changes for NEC machine (re-order main
!                  loop to improve vectorization).  A.C. Bushell
!
! 6.2   08/06/2005 Further optimization changes for NEC machine.
!                  Vectorization by expanding scalar variables and
!                  breaking loops. A_A, B_B, DDU, MCRIT, MGUESS, MXL,
!                  MXI are expanded in four dimensions. Local variable
!                  MTEST is expanded in the row direction and the Newton
!                  search direction (jwhile) only.
!                                                   JC. Rioual
!
! 6.2   16/06/2005 Changes to facilitate resolution tuning of scheme.
!                  Matching code with new UMDP for USSP sheme.
!                                                   A.C. Bushell
!
! 6.4   19/05/2006 Algorithm changes to speed up scheme (pre-selection).
!                  Correction of diagnostic flux output.
!                                                   A.C. Bushell
!
! code description:
!   language: fortran 77 + cray extensions
!   this code is written to umdp3 programming standards.
!   documentation: Unified Model Documentation Paper 34 (Non-Orog GW)

      IMPLICIT NONE

!--------------------
! Global Parameters
!--------------------
#include "c_a.h"
#include "c_g.h"
#include "c_pi.h"
#include "c_r_cp.h"
#include "c_omega.h"
#include "fldtype.h"
#include "domtyp.h"
!----------------
! Local constants
!----------------
#include "c_gwave.h"
! ----------------------------------------------------------------------+-------
!     Subroutine arguments of GW_USSP
! ----------------------------------------------------------------------+-------
      INTEGER                                                           &
     &  LEVELS                                                          &
                             !IN Number of model levels
     &, ROWS                                                            &
                             !IN Number of rows for u field
     &, NROWS                                                           &
                             !IN Number of rows for v field
     &, OFF_X                                                           &
                             !IN offset longitude
     &, OFF_Y                                                           &
                             !IN offset latitude
     &, HALO_I                                                          &
                             !IN Halo in longitude
     &, HALO_J                                                          &
                             !IN Halo in latitude
     &, ROW_LENGTH                                                      &
                             !IN Number of grid points in row
     &, MODEL_DOMAIN         !IN Model type (global, LAM etc)
      REAL                                                              &
     &  SIN_THETA_LONGITUDE(ROW_LENGTH,ROWS)                            &
                                              !Grid point longitudes
     &, SIN_THETA_LATITUDE(ROW_LENGTH,ROWS)                             &
                                              !P-GRID Latitudes
     &, GWSPEC_EFLUX(ROW_LENGTH,ROWS,LEVELS)                            &
                                              ! Fp in each of 4
     &, GWSPEC_SFLUX(ROW_LENGTH,NROWS,LEVELS)                           &
                                              ! azimuths for diags.
     &, GWSPEC_WFLUX(ROW_LENGTH,ROWS,LEVELS)                            &
                                              !
     &, GWSPEC_NFLUX(ROW_LENGTH,NROWS,LEVELS)                           &
                                              !
     &, GWSPEC_EWACC(ROW_LENGTH,ROWS,LEVELS)                            &
                                              !
     &, GWSPEC_NSACC(ROW_LENGTH,NROWS,LEVELS)                           &
                                              !
     &, ETA_THETA_LEVELS(0:LEVELS)                                      &
     &, THETA(ROW_LENGTH,ROWS,LEVELS)                                   &
                             !IN Primary model array for theta
     &, RHO(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y,LEVELS)         &
                             !IN Primary model array for Density
                             !x (radius earth)^2.
     &, TIMESTEP                                                        &
                             !IN Timestep
     &, U(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y,LEVELS)           &
                             !INOUT Primary model array for U field
     &, V(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:NROWS+OFF_Y,LEVELS)          &
                             !INOUT Primary model array for V field
     &, R_RHO_LEVELS(1-HALO_I:ROW_LENGTH+HALO_I,1-HALO_J:ROWS           &
     &               +HALO_J,LEVELS)                                    &
                   !Distance of rho levels from Earth centre.
     &, R_THETA_LEVELS(1-HALO_I:ROW_LENGTH+HALO_I,1-HALO_J:             &
     &                 ROWS+HALO_J,0:LEVELS)                            &
                   !Distance of theta levels from Earth centre.
     &, P_LAYER_BOUNDARIES(ROW_LENGTH,ROWS,0:LEVELS)                    &
     &, R_U(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y,LEVELS)         &
     &, R_V(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:NROWS+OFF_Y,LEVELS)
      LOGICAL                                                           &
     &  L_USSP_OPAQUE   !IN  Switch for Opaque Upper boundary condition
!
! ----------------------------------------------------------------------+-------
!     Local parameters
! ----------------------------------------------------------------------+-------
!
!     Max number of directions, typically four.
      Integer, Parameter :: IDIR           = 4
!
!     N processor address
      Integer, Parameter :: PNORTH         = 1
!
!     E processor address
      Integer, Parameter :: PEAST          = 2
!
!     S processor address
      Integer, Parameter :: PSOUTH         = 3
!
!     W processor address
      Integer, Parameter :: PWEST          = 4
!
!     Maximum number of iterations of Newton Raphson DO (While) loop
      Integer, Parameter :: MAXWHILE       = 9
!
      Real, Parameter :: R_EARTH_RADIUS_SQ =                            &
     &                                  1. / (Earth_Radius*Earth_Radius)
!
!     Reciprocal of middle atmosphere mean scale height for pressure,
!     normally assumed to be around 7km
      Real, Parameter :: RSCALE_H          = G / (R * 239.145)
!
!     Parameter beta in the launch spectrum total energy equation
      Real, Parameter :: BETA_E0           = 1.0227987125E-1
!
!     Azimuthal sector for launch spectrum integral Delta Phi / 2
      Real, Parameter :: DDPHIR2           = PI / IDIR
!
!     Parameter p in B_0(p) for launch spectrum intrinsic frequency
!     NOTE: This parameter determines the intrinsic frequency spectrum
!           shape and hence the integral form in 4.1, which is strictly
!           valid only for p > 1. !!IF contemplating changes BE WARNED!!
      Real, Parameter :: PSAT              = 5.0 / 3.0
!
!     Psat - 1
      Real, Parameter :: PSATM1            = PSAT - 1.0
!
!     2 - Psat
      Real, Parameter :: TWOMPSAT          = 2.0 - PSAT
!
!     Wavenumber at peak in spectrum
      Real, Parameter ::  MSTAR            = 2.0 * PI * LSTAR
!
!     Reciprocal of mstar (m) and mstar^2 (m^2)
      Real, Parameter ::  RMSTAR           = 1. / MSTAR
      Real, Parameter ::  RMSTARSQ         = RMSTAR * RMSTAR
!
!     Minimum vertical wavenumber at launch (/m )
      Real, Parameter ::  MMINL            = 2.0 * PI * LMINL
!
!     Normalised minimum vertical wavenumber at launch
      Real, Parameter ::  MNLMIN           = LMINL / LSTAR
!
!     Equatorial planetary vorticity gradient parameter B_eq (/m /s )
      Real, Parameter :: BETA_EQ_RMSTAR    = 2.3E-11 * RMSTAR
!
!     Power s of vertical wavenumber spectrum A_0(s,t) at low m
      Real, Parameter :: SS                = 1.0
!
!     s + 1, s - 1
      Real, Parameter :: SSP1              = SS + 1.0
!
!     Power t=t_sat of vertical wavenumber spectrum at large m due to
!     saturation by gravity wave breaking (and shape of chopping fn)
      Real, Parameter :: TT                = 3.0
!
!     t - 1, t - 2, 1 / (t-2), (t-3) / (t-2), 2 - t
      Real, Parameter :: TTM1              = TT - 1.0
      Real, Parameter :: TTM2              = TT - 2.0
      Real, Parameter :: RTTM2             = 1.0 / TTM2
      Real, Parameter :: TTRAT             = (TT - 3.0) * RTTM2
      Real, Parameter :: TWOMTT            = 2.0 - TT
!
!     s + t, 1 / (s+t)
      Real, Parameter :: SSPTT             = SS + TT
      Real, Parameter :: RSSPTT            = 1.0 / SSPTT
!
!     Weight for (n+1)th guess at mNlX in iteration solution
      Real, Parameter :: MWEIGHT           = 0.8
!
!     Strength coefficient constant for Launch spectrum (CCL / A0)
      Real, Parameter :: CCL0 = 3.41910625e-9
!
! ----------------------------------------------------------------------+-------
!     Security parameters
! ----------------------------------------------------------------------+-------
!
!     Minimum allowed value of buoyancy frequency squared
      Real, Parameter ::  SQNMIN           = 1.0E-4
!
!     Minimum allowed non-zero value of Curvature Coefficient A
      Real, Parameter ::  ASECP            =  1.0E-20
      Real, Parameter ::  ASECN            = -(ASECP)
!
! ----------------------------------------------------------------------+-------
!     Local Constants (a) Physical
! ----------------------------------------------------------------------+-------
      REAL                                                              &
     &  A0_R_SP1TM1                                                     &
                           !A0/(s+1)(t-1) normalisation factor for the
!                          ! launch spectrum vertical wavenumber
     &, A0_R_1MT                                                        &
                           !-A0/(t-1) norm factor for launch spectrum m
     &, TAIL_CHOP2B                                                     &
                           ! Integral segment [(s+1) * mNLmin**(1-t)]
     &, HEAD_CHOP2A                                                     &
                           ! Integral segment [(t-1) * mNLmin**(s+1)]
     &, COSPHI(4)                                                       &
                           !cos(phi_j)
     &, SINPHI(4)          !sin(phi_j)
!
! ----------------------------------------------------------------------+-------
!     Local variables (scalars) used in GW_USSP
!     Some effectively expanded to workspace (using vector registers)
! ----------------------------------------------------------------------+-------
      INTEGER                                                           &
     &  I                                                               &
                            !longitude
     &, J                                                               &
                            !latitude
     &, K                                                               &
                            !level
     &, ILAUNCH                                                         &
                            !Minimum level for launch level (all points)
     &, LAUNCHLEV                                                       &
                            !Launch level at specific point
     &, NNI, NNJ, NNJD                                                  &
                            !Index values
     &, NCHOP2                                                          &
                            !Number of spectra with low m intersect
!    &, NSPIRA              !Number spiral A solution points
     &, JDIR                                                            &
                            !direction
     &, JLEV                                                            &
                            !level
     &, JWHILE              !Counter for while loop
      REAL                                                              &
     &  F_F                                                             &
                            !Inertial frequency at current latitude
                            !(rad s^-1)
     &, MNLY                                                            &
                            ! High wavenumber intersection point
     &, MKILL                                                           &
                            ! Wavenumber where Doppler transformed
!                             spectrum is reduced to zero
     &, OMINRNBV                                                        &
                            ! omega_min(launch) / N (k)
     &, CCS0                                                            &
                            ! Constant component of saturation spectrum
     &, FMINUS                                                          &
                            ! Minimum range value of function f
     &, FTERM                                                           &
                            ! Intermediate  value of function f
     &, FPLUS                                                           &
                            ! Maximum range value of function f
     &, GMINUS                                                          &
                            ! Minimum range value of function g
     &, GTERM                                                           &
                            ! Intermediate  value of function g
     &, GPLUS               ! Maximum range value of function g
      LOGICAL                                                           &
     &  GWSPEC_EFLUX_ON                                                 &
                            !Switches for diagnostics of
     &, GWSPEC_SFLUX_ON                                                 &
                            !Fp in each of 4 azimuths
     &, GWSPEC_WFLUX_ON                                                 &
                            !
     &, GWSPEC_NFLUX_ON                                                 &
                            !
     &, GWSPEC_EWACC_ON                                                 &
                            !and accelerations
     &, GWSPEC_NSACC_ON                                                 &
                            !
     &, AT_EXTREMITY(4)                                                 &
                            !Edge of domain indicator
     &, L_CHOP2             ! Indicates spectra with low m intersect
! ----------------------------------------------------------------------+-------
!     Local variables (dynamic arrays) used in GW_USSP
! ----------------------------------------------------------------------+-------
!     INTEGER
!    &  NVIEW(LEVELS)       ! check output array
!    &, NVIEW2(LEVELS)      ! check output array
      REAL                                                              &
     &  DDU_a(ROW_LENGTH,ROWS,LEVELS,IDIR)                              &
!                            Delta U=udotk(launch)-udotk(jlev)
     &, FPFAC(ROW_LENGTH,ROWS,LEVELS,IDIR)                              &
                                           ! Record of total flux
     &, FPTOT(ROW_LENGTH,ROWS,LEVELS,IDIR)                              &
                                           ! Pseudomomentum flux
!                            integrated over
                            !azimuthal sector (kg m^-1 s^-2)
     &, G_G(ROW_LENGTH,ROWS,LEVELS,IDIR)                                &
                                           !Wave-induced force per unit
                            !mass due to azimuthal sectors (m s^-2)
     &, G_X(ROW_LENGTH,ROWS,LEVELS)                                     &
                                       !Zonal component of wave-induced
                            !force at longitude, level (m s^-2)
                            !Note that in our notation G_X is
                            !equivalent to DU_DT, the zonal
                            !wind tendency (m s^-2)
     &, G_Y(ROW_LENGTH,ROWS,LEVELS)                                     &
                                   !Meridional component of wave-induced
                            !force at longitude, level (m s^-2)
                            !Note that in our notation G_Y is
                            !equivalent to DV_DT, the meridional
                            !wind tendency (m s^-2)
     &, ACOEFF(ROW_LENGTH,ROWS,LEVELS,IDIR)                             &
!                            Coefficient A in intersect point equation
     &, CURVATURE(ROW_LENGTH,ROWS,LEVELS,IDIR)                          &
!                            Term (A / B) in intersect point equation
     &, ATTENUATION(ROW_LENGTH,ROWS,LEVELS,IDIR)                        &
!                            Coefficient B in intersect point equation
     &, INTERCEPT1(ROW_LENGTH,ROWS,LEVELS,IDIR)                         &
!                            Chop function B*[1 + (A/B)]^(t-2)
     &, MINTERCEPT(ROW_LENGTH,ROWS,LEVELS,IDIR)                         &
!                            Chop function B*[1 + (A/B)*mNlmin]^(t-2)
     &, MGUESS_a(ROW_LENGTH,ROWS,LEVELS+1,IDIR)                         &
!                            Starting value of vertical wavenumber for
!                            crossing point search
     &, MNLX(IDIR*ROWS*ROW_LENGTH,0:MAXWHILE)                           &
                                              ! Intersect mNlX estimates
     &, INDEXI(IDIR*ROWS*ROW_LENGTH)                                    &
                                      ! I location of chop type points
     &, INDEXJ(IDIR*ROWS*ROW_LENGTH)                                    &
                                      ! J location of chop type points
     &, INDEXJD(IDIR*ROWS*ROW_LENGTH)                                   &
                                      ! JDIR location of chop type pnts
     &, ATTE_C(IDIR*ROWS*ROW_LENGTH)                                    &
                                      ! Compressed attenuation array
     &, CURV_C(IDIR*ROWS*ROW_LENGTH)                                    &
                                      ! Compressed curvature array
     &, WGTN(IDIR*ROWS*ROW_LENGTH)                                      &
                                      ! Weighting of n term in iter
     &, OMIN(ROW_LENGTH,ROWS)                                           &
                               ! Either f_f or the equatorial minimum
!                            frequency, whichever is less  (rad s^-1)
     &, NBV(ROW_LENGTH,ROWS,LEVELS)                                     &
                                           ! Buoyancy frequency
!                             [Brunt Vaisala frequency] on half-levels
                            !(rad s^-1)
     &, RHOCL(ROW_LENGTH,ROWS,IDIR)                                     &
                                           ! [Rho . Cl]_klaunch
     &, RHOCSK(ROW_LENGTH,ROWS,LEVELS)                                  &
                                           ! [Rho(z) . Csat(z)]_k
     &, RHO_TH(ROW_LENGTH,ROWS,LEVELS)                                  &
                                           ! Rho on theta levels
     &, UDOTK(ROW_LENGTH,ROWS,LEVELS+1,IDIR)                            &
                                             ! Component of wind
!                            in phi_jdir direction (m s^-1)
     &, UONP(ROW_LENGTH,ROWS,LEVELS)                                    &
                            !U and V on Rho grid for local use
     &, VONP(ROW_LENGTH,ROWS,LEVELS)                                    &
     &, UONP_HALO(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y,          &
     &            LEVELS)                                               &
                            !U and V with halo for local use
     &, VONP_HALO(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y,          &
     &            LEVELS)                                               &
     &, RHONT(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y,              &
     &        LEVELS)                                                   &
                            !Rho on Theta grid.
     &, UINC(ROW_LENGTH,ROWS,LEVELS)                                    &
                                         ! Increments on rho grid.
     &, VINC(ROW_LENGTH,NROWS,LEVELS)                                   &
                                          !
     &, WORK_HALO(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y,          &
     &              0:LEVELS)
      LOGICAL                                                           &
     &  L_FTHENG(IDIR*ROWS*ROW_LENGTH) ! Indicate dir of spiral solution
!
!  External subroutine calls: ------------------------------------------+-------
      EXTERNAL P_TO_T, U_TO_P, V_TO_P, P_TO_U, P_TO_V
!- End of Header
!
! ==Main Block==--------------------------------------------------------+-------
! --------------------------------
!     Local Constants (a) Physical
! --------------------------------
      DATA COSPHI/0,-1,0,1/
      DATA SINPHI/1,0,-1,0/
!
      FMINUS = MNLMIN**SSPTT
      TAIL_CHOP2B = SSP1 / (MNLMIN**TTM1)
      HEAD_CHOP2A = TTM1 * (MNLMIN**SSP1)
      A0_R_SP1TM1 = 1.0 / ( SS + TT - HEAD_CHOP2A )
      A0_R_1MT    = (-(SSP1) ) * A0_R_SP1TM1
      CCS0 = BETA_E0 * SIN(DDPHIR2) * PSATM1 / (PI * TWOMPSAT)
!
! ----------------------------------------------------------------------+-------
! Find model level for launch
! ----------------------------------------------------------------------+-------
!     Minimum value if variable height sources are used
      ILAUNCH=1
      DO K=2,LEVELS
        IF (ETA_THETA_LEVELS(K) >  ETALAUNCH.AND.                       &
     &      ETA_THETA_LEVELS(K-1) <  ETALAUNCH) THEN
          IF ((ETA_THETA_LEVELS(K)-ETALAUNCH) <                         &
     &        (ETALAUNCH-ETA_THETA_LEVELS(K-1))) THEN
            ILAUNCH=K
          ELSE
            ILAUNCH=K-1
          ENDIF
        ENDIF
      ENDDO
      LAUNCHLEV = ILAUNCH
!
! ----------------------------------------------------------------------+-------
!     Interpolate RHO onto T grid and U,V onto P grid
! ----------------------------------------------------------------------+-------
! DEPENDS ON: p_to_t
      CALL P_TO_T (ROW_LENGTH, ROWS, HALO_I, HALO_J,                    &
     &             OFF_X,OFF_Y,LEVELS-1,R_THETA_LEVELS,                 &
     &             R_RHO_LEVELS,RHO,RHONT)
!     Extrapolate topmost level to be a scale height from level below
      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        RHONT(I,J,LEVELS)=RHONT(I,J,LEVELS-1)                           &
     &   *EXP(-(R_THETA_LEVELS(I,J,LEVELS)-                             &
     &          R_THETA_LEVELS(I,J,LEVELS-1)) * RSCALE_H)
       ENDDO
      ENDDO
! DEPENDS ON: u_to_p
      CALL U_TO_P (U,ROW_LENGTH,ROWS,LEVELS,                            &
     &            OFF_X,OFF_Y,MODEL_DOMAIN,                             &
     &            AT_EXTREMITY,UONP)
! DEPENDS ON: v_to_p
      CALL V_TO_P (V,ROW_LENGTH,ROWS,NROWS,LEVELS,                      &
     &            OFF_X,OFF_Y,MODEL_DOMAIN,                             &
     &            AT_EXTREMITY,VONP)
! ---------------------------
!     Set polar winds to zero
! ---------------------------
      If(model_domain  ==  mt_global) Then
        If (at_extremity(psouth) ) Then
          VONP(:,1,:) = 0.0
          UONP(:,1,:) = 0.0
        Endif
        If (at_extremity(pnorth) ) Then
          VONP(:,rows,:) = 0.0
          UONP(:,rows,:) = 0.0
        Endif
      Endif
!
! ----------------------------------------------------------------------+-------
!     Initialize local arrays to zero
! ----------------------------------------------------------------------+-------
! ---------------------------------
!     Set winds with haloes to zero
! ---------------------------------
      Levels_do1: DO K=1,LEVELS
        Rows_do1: DO J=1-OFF_Y,ROWS+OFF_Y
          Row_length_do1: DO I=1-OFF_X,ROW_LENGTH+OFF_X
            UONP_HALO(I,J,K) = 0.
            VONP_HALO(I,J,K) = 0.
          END DO  Row_length_do1
        END DO  Rows_do1
      END DO  Levels_do1
!
      Levels_do1a: DO K=1,LEVELS
        Rows_do1a: DO J=1,ROWS
          Row_length_do1a: DO I=1,ROW_LENGTH
! ----------------------------------------------------------------------+-------
!           Zero vertical divergence of pseudomomentum flux
! ----------------------------------------------------------------------+-------
            G_X(I,J,K) = 0.0
            G_Y(I,J,K) = 0.0
          END DO  Row_length_do1a
        END DO  Rows_do1a
!       NVIEW(K) = 0
!       NVIEW2(K) = 0
      END DO  Levels_do1a
!
! ----------------------------------------------------------------------+-------
! 1.0   Set variables that are to be defined on all model levels
! ----------------------------------------------------------------------+-------
!
      Levels_do2: DO JLEV=2,(LEVELS - 1)
! ----------------------------------------------------------------------+-------
! 1.1   Density, buoyancy frequency and altitude for middle levels
! ----------------------------------------------------------------------+-------
        Rows_do2: DO J=1,ROWS
          Row_length_do2: DO I=1,ROW_LENGTH
            RHO_TH(I,J,JLEV) = RHONT(I,J,JLEV) * r_earth_radius_sq
!           Buoyancy (Brunt-Vaisala) frequency calculation
            NBV(I,J,JLEV) = ( g*(THETA(I,J,JLEV+1)-THETA(I,J,JLEV-1))/  &
     &                          (THETA(I,J,JLEV) *                      &
     &       (R_THETA_LEVELS(I,J,JLEV+1)-R_THETA_LEVELS(I,J,JLEV-1))) )
            NBV(I,J,JLEV) = MAX(NBV(I,J,JLEV), SQNMIN)
            NBV(I,J,JLEV) = SQRT(NBV(I,J,JLEV))
          END DO  Row_length_do2
        END DO  Rows_do2
      END DO  Levels_do2
!
! ----------------------------------------------------------------------+-------
!     Density, Buoyancy (Brunt-Vaisala) frequency at top and bottom
! ----------------------------------------------------------------------+-------
      Rows_do3: DO J=1,ROWS
        Row_length_do3: DO I=1,ROW_LENGTH
!       Set rho and theta level value of Z_TH.
          RHO_TH(I,J,1)      = RHONT(I,J,1) * r_earth_radius_sq
          NBV(I,J,1)         = NBV(I,J,2)
          RHO_TH(I,J,LEVELS) = RHONT(I,J,LEVELS) * r_earth_radius_sq
          NBV(I,J,LEVELS)    = NBV(I,J,LEVELS-1)
        END DO  Row_length_do3
      END DO  Rows_do3
!
      Levels_do4: DO JLEV=LEVELS,2,-1
! ----------------------------------------------------------------------+-------
! 1.2   Set buoyancy frequency constant up to 1km altitude
! ----------------------------------------------------------------------+-------
        Rows_do4: DO J=1,ROWS
          Row_length_do4: DO I=1,ROW_LENGTH
            IF ( (R_THETA_LEVELS(I,J,JLEV) - EARTH_RADIUS) <  1.0E3)    &
     &      NBV(I,J,JLEV-1) = NBV(I,J,JLEV)
          END DO  Row_length_do4
        END DO  Rows_do4
      END DO  Levels_do4
!
! ----------------------------------------------------------------------+-------
! 1.3  Compute component of wind U in each wave-propagation direction.
!      U is the dot product of (u,v) with k_0 but n.b. UDOTK is half-way
!      between rho levels.
! ----------------------------------------------------------------------+-------
      IDir_do1: DO JDIR=1,IDIR
!
        Levels_do5: DO JLEV=1,LEVELS-1
          Rows_do5: DO J=1,ROWS
            Row_length_do5: DO I=1,ROW_LENGTH
!           Assume theta levels are half way between rho levels.
              UDOTK(I,J,JLEV,JDIR) =                                    &
     &        0.5*(UONP(I,J,JLEV) + UONP(I,J,JLEV+1))*COSPHI(JDIR) +    &
     &        0.5*(VONP(I,J,JLEV) + VONP(I,J,JLEV+1))*SINPHI(JDIR)
            END DO  Row_length_do5
          END DO  Rows_do5
        END DO  Levels_do5
!
! ----------------------------------------------------------------------+-------
!      Set wind component for top level, to be equal to that on the top
!      Rho level, and total flux of horizontal pseudomomentum at bottom
! ----------------------------------------------------------------------+-------
        Rows_do5a: DO J=1,ROWS
          Row_length_do5a: DO I=1,ROW_LENGTH
            UDOTK(I,J,LEVELS,JDIR) =                                    &
     &       UONP(I,J,LEVELS)*COSPHI(JDIR) + VONP(I,J,LEVELS)*SINPHI(JDIR)
            FPTOT(I,J,1,JDIR) = 0.0
          END DO  Row_length_do5a
        END DO  Rows_do5a
!
! ----------------------------------------------------------------------+-------
! 2.0  Initialize variables that need to be defined up to launch level
! ----------------------------------------------------------------------+-------
!
        Levels_do6: DO JLEV=1,ILAUNCH
          Rows_do6: DO J=1,ROWS
            Row_length_do6: DO I=1,ROW_LENGTH
! ----------------------------------------------------------------------+-------
!           Vertical divergence of total flux of horizontal pseudomom
! ----------------------------------------------------------------------+-------
              G_G(I,J,JLEV,JDIR) = 0.0
            END DO  Row_length_do6
          END DO  Rows_do6
        END DO  Levels_do6
!
! ----------------------------------------------------------------------+-------
! 3.0  Initialize gravity wave spectrum variables for the launch level
! ----------------------------------------------------------------------+-------
!
        Rows_do7: DO J=1,ROWS
          Row_length_do7: DO I=1,ROW_LENGTH
! ----------------------------------------------------------------------+-------
!         Value at launch of (m*^2 * Total vertical flux of horizontal
!         pseudomomentum ... UMDP 34: 1.14)
! ----------------------------------------------------------------------+-------
            RHOCL(I,J,JDIR) = RHO_TH(I,J,ILAUNCH) * CCL0
          END DO  Row_length_do7
        END DO  Rows_do7
      END DO  IDir_do1
!
! ----------------------------------------------------------------------+-------
! 3.1 Compute minimum intrinsic frequency OMIN from inertial frequency
!     squared F_F. See UMDP 34, eqn. (A.7).
!-----------------------------------------------------------------------+-------
      Rows_do7a: DO J=1,ROWS
        Row_length_do7a: DO I=1,ROW_LENGTH
          F_F = (TWO_OMEGA * SIN_THETA_LATITUDE(I,J))**2
          OMIN(I,J) = NBV(I,J,ILAUNCH) * BETA_EQ_RMSTAR
          OMIN(I,J) = MAX( OMIN(I,J), F_F )
          OMIN(I,J) = SQRT(OMIN(I,J))
        END DO  Row_length_do7a
      END DO  Rows_do7a
!
! ----------------------------------------------------------------------+-------
! 4.0 Calculations carried out at levels from Launch to Lid
! ----------------------------------------------------------------------+-------
!
      Levels_do8: DO JLEV=2,LEVELS
        Rows_do8: DO J=1,ROWS
          Row_length_do8: DO I=1,ROW_LENGTH
            IDir_do1a: DO JDIR=1,IDIR
! ----------------------------------------------------------------------+-------
!           Total vertical flux of horizontal pseudomomentum (analytic
!           integral under curve = rho_l * C_l / m*^2 ... UMDP 34: 1.14)
! ----------------------------------------------------------------------+-------
              IF (JLEV == LAUNCHLEV)  THEN
                FPTOT(I,J,JLEV,JDIR) = RMSTARSQ * RHOCL(I,J,JDIR)
              ELSE
                FPTOT(I,J,JLEV,JDIR) = 0.0
              END IF
              FPTOT(I,J,JLEV,JDIR) = FPTOT(I,J,JLEV-1,JDIR) +           &
     &                               FPTOT(I,J,JLEV,JDIR)
            END DO  IDir_do1a
          END DO  Row_length_do8
        END DO  Rows_do8
      END DO  Levels_do8
!
! ----------------------------------------------------------------------+-------
! 4.1 Compute [rho(z) . C(z)]_k for combined spectrum ...(UMDP 34: 1.16)
! ----------------------------------------------------------------------+-------
!
!     IF (ABS(PSATM1) >= 0.1) THEN
!     For current setting of parameter psat this test is always true
        Levels_do8a: DO JLEV=ILAUNCH,LEVELS
          Rows_do8a: DO J=1,ROWS
            Row_length_do8a: DO I=1,ROW_LENGTH
              OMINRNBV = OMIN(I,J) / NBV(I,J,JLEV)
              RHOCSK(I,J,JLEV) = RHO_TH(I,J,JLEV) * CCS0 *              &
     &         (NBV(I,J,JLEV))**2 * (OMINRNBV**PSATM1) *                &
     &         (1.0 - (OMINRNBV**TWOMPSAT)) / (1.0 - (OMINRNBV**PSATM1))
            END DO  Row_length_do8a
          END DO  Rows_do8a
        END DO  Levels_do8a
!     ELSE
!     Require a different functional form for normalisation factor B0
!             BBS = 1.0 / ALOG(NBV(I,J,JLEV) / OMIN(I,J))
!     END IF
! ----------------------------------------------------------------------+-------
!     Loop over directions and levels and calculate horizontal component
!     of the vertical flux of pseudomomentum for each azimuthal
!     direction and for each altitude
! ----------------------------------------------------------------------+-------
      Levels_do92: DO JLEV=ILAUNCH+1,LEVELS
        L_CHOP2 = .FALSE.
        IDir_do2a: DO JDIR=1,IDIR
          Rows_do92: DO J=1,ROWS
            Row_length_do92: DO I=1,ROW_LENGTH
! ----------------------------------------------------------------------+-------
!     Initialise MGUESS (start point for iterative searches if needed)
! ----------------------------------------------------------------------+-------
              MGUESS_a(I,J,JLEV,JDIR) = 0.0
              Fptot_if1: IF (FPTOT(I,J,JLEV-1,JDIR) >  0.0) THEN
! ----------------------------------------------------------------------+-------
! 4.2       Calculate variables that define the Chop Type Cases.
! ----------------------------------------------------------------------+-------
              DDU_a(I,J,JLEV,JDIR) = UDOTK(I,J,ILAUNCH,JDIR)            &
     &                                - UDOTK(I,J,JLEV,JDIR)
! ----------------------------------------------------------------------+-------
!             UMDP 34: 1.23 coefficient B
! ----------------------------------------------------------------------+-------
              ATTENUATION(I,J,JLEV,JDIR) =                              &
     &           ( RHOCSK(I,J,JLEV) / RHOCL(I,J,JDIR) ) *               &
     &           ( NBV(I,J,ILAUNCH) / NBV(I,J,JLEV) )**TTM1
! ----------------------------------------------------------------------+-------
!             UMDP 34: 1.22 coefficient A = (A/B) * B
! ----------------------------------------------------------------------+-------
              CURVATURE(I,J,JLEV,JDIR)  =  DDU_a(I,J,JLEV,JDIR) * MSTAR &
     &                                     / NBV(I,J,ILAUNCH)
!
              ACOEFF(I,J,JLEV,JDIR)     =  CURVATURE(I,J,JLEV,JDIR) *   &
     &                                   ATTENUATION(I,J,JLEV,JDIR)
!
              MINTERCEPT(I,J,JLEV,JDIR) = ATTENUATION(I,J,JLEV,JDIR)    &
     &         * ( (MNLMIN * CURVATURE(I,J,JLEV,JDIR)) + 1.0 )**TTM2
!
              Curv_if1: IF (CURVATURE(I,J,JLEV,JDIR) <  ASECN)  THEN
! ----------------------------------------------------------------------+-------
!             Negative Doppler Shift : factor will hit zero (kill point)
! ----------------------------------------------------------------------+-------
                MKILL = 1.0 / ABS(CURVATURE(I,J,JLEV,JDIR))
!
                Mkill_if1: IF (MKILL <= MNLMIN )  THEN
! ----------------------------------------------------------------------+-------
!               Chop Type IV : No flux propagates
! ----------------------------------------------------------------------+-------
                  FPTOT(I,J,JLEV,JDIR) = 0.0
                ELSE
                  IF (MKILL >  1.0)  THEN
                  INTERCEPT1(I,J,JLEV,JDIR) = ATTENUATION(I,J,JLEV,JDIR)&
     &              * ( 1.0 + CURVATURE(I,J,JLEV,JDIR) )**TTM2
                  ELSE
! ----------------------------------------------------------------------+-------
!                 Doppler factor minimum (kill point) situated below
!                 mstar in the low-m part of the launch spectrum
! ----------------------------------------------------------------------+-------
                    INTERCEPT1(I,J,JLEV,JDIR) = 0.0
                  END IF
!
                  Lowend_if1: IF (INTERCEPT1(I,J,JLEV,JDIR) >= 1.0) THEN
! ----------------------------------------------------------------------+-------
!                 Chop Type I: Intersection in high wavenumber part only
! ----------------------------------------------------------------------+-------
                    MNLY = ( ATTENUATION(I,J,JLEV,JDIR)**TTRAT -        &
     &              ATTENUATION(I,J,JLEV,JDIR)) / ACOEFF(I,J,JLEV,JDIR)
!
                    FPTOT(I,J,JLEV,JDIR) = FPTOT(I,J,JLEV,JDIR) *       &
     &     (1.0 - (A0_R_1MT * CURVATURE(I,J,JLEV,JDIR) * MNLY**TWOMTT))
                  ELSE
                    IF (MINTERCEPT(I,J,JLEV,JDIR) <= FMINUS)  THEN
! ----------------------------------------------------------------------+-------
!                 Chop Type IIb: Low wavenumber intersect only below min
! ----------------------------------------------------------------------+-------
                      FPTOT(I,J,JLEV,JDIR) = FPTOT(I,J,JLEV,JDIR) *     &
     &            A0_R_SP1TM1 * TAIL_CHOP2B * MINTERCEPT(I,J,JLEV,JDIR) &
     &                 * ( (MNLMIN * CURVATURE(I,J,JLEV,JDIR)) + 1.0 )
                    ELSE
! ----------------------------------------------------------------------+-------
!                 Chop Type IIa: Low wavenumber intersect only
! ----------------------------------------------------------------------+-------
                      L_CHOP2 = .TRUE.
                      MGUESS_a(I,J,JLEV,JDIR) = MIN(MKILL, 1.0)
                      FPFAC(I,J,JLEV,JDIR) = FPTOT(I,J,JLEV,JDIR) * A0_R_SP1TM1
                      FPTOT(I,J,JLEV,JDIR) = 0.0
                    END IF
                  END IF  Lowend_if1
!
                END IF  Mkill_if1
!
              ELSE IF (CURVATURE(I,J,JLEV,JDIR) >  ASECP)  THEN
! ----------------------------------------------------------------------+-------
!             Positive Doppler Shift : non-zero factor (no kill point)
! ----------------------------------------------------------------------+-------
                INTERCEPT1(I,J,JLEV,JDIR) = ATTENUATION(I,J,JLEV,JDIR)  &
     &            * ( 1.0 + CURVATURE(I,J,JLEV,JDIR) )**TTM2
!
                Chop3_if1: IF (INTERCEPT1(I,J,JLEV,JDIR) <  1.0)  THEN
! ----------------------------------------------------------------------+-------
!               Chop Type III: Intersection in both wavenumber parts
! ----------------------------------------------------------------------+-------
                  FPFAC(I,J,JLEV,JDIR) = FPTOT(I,J,JLEV,JDIR) * A0_R_SP1TM1
!
! ----------------------------------------------------------------------+-------
!                 First find intersect in high wavenumber part
!                 UMDP 34: 1.25
! ----------------------------------------------------------------------+-------
                  MNLY = ( ATTENUATION(I,J,JLEV,JDIR)**TTRAT -          &
     &              ATTENUATION(I,J,JLEV,JDIR)) / ACOEFF(I,J,JLEV,JDIR)
!
                  FPTOT(I,J,JLEV,JDIR) = FPTOT(I,J,JLEV,JDIR) *         &
     &              A0_R_1MT * CURVATURE(I,J,JLEV,JDIR) * MNLY**TWOMTT
!
! ----------------------------------------------------------------------+-------
!                 Then find intersect in low wavenumber part to reckon
!                 its flux contribution for addition when available
! ----------------------------------------------------------------------+-------
                  IF (MINTERCEPT(I,J,JLEV,JDIR) <= FMINUS)  THEN
! ----------------------------------------------------------------------+-------
!                 Chop Type IIIb: Low wavenumber intersect below min
! ----------------------------------------------------------------------+-------
                    FPTOT(I,J,JLEV,JDIR) = FPTOT(I,J,JLEV,JDIR) +       &
     &                                   ( FPFAC(I,J,JLEV,JDIR) *       &
     &                 TAIL_CHOP2B *  MINTERCEPT(I,J,JLEV,JDIR) *       &
     &                ( (MNLMIN * CURVATURE(I,J,JLEV,JDIR)) + 1.0 ) )
                  ELSE
! ----------------------------------------------------------------------+-------
!                 Chop Type IIIa: Low wavenumber intersect
! ----------------------------------------------------------------------+-------
                    L_CHOP2 = .TRUE.
                    MGUESS_a(I,J,JLEV,JDIR) = 1.0
                  END IF
!
! ----------------------------------------------------------------------+-------
!               ELSE Chop Type 0: No intersection (spectrum unaltered)
! ----------------------------------------------------------------------+-------
                END IF  Chop3_if1
              ELSE
! ----------------------------------------------------------------------+-------
!             Negligible Doppler shift
! ----------------------------------------------------------------------+-------
!               Strictly this is analytic solution mNLX.  UMDP 34: 1.27
                MNLY = ATTENUATION(I,J,JLEV,JDIR)**RSSPTT
                IF (MNLY <= MNLMIN)  THEN
! ----------------------------------------------------------------------+-------
!               Chop Type IIb: Low wavenumber intersect only below min
! ----------------------------------------------------------------------+-------
                  FPTOT(I,J,JLEV,JDIR) = FPTOT(I,J,JLEV,JDIR) *         &
     &           A0_R_SP1TM1 * TAIL_CHOP2B * ATTENUATION(I,J,JLEV,JDIR)
                ELSE
                  IF (MNLY <  1.0)  FPTOT(I,J,JLEV,JDIR) =              &
! ----------------------------------------------------------------------+-------
!               Chop Type IIc: Low wavenumber intersect only (analytic)
! ----------------------------------------------------------------------+-------
     &              FPTOT(I,J,JLEV,JDIR) * A0_R_SP1TM1 *                &
     &               ( (SSPTT * (MNLY**SSP1)) - HEAD_CHOP2A )
! ----------------------------------------------------------------------+-------
!                 ELSE Chop Type 0: No intersection (spectrum unaltered)
! ----------------------------------------------------------------------+-------
                END IF
              END IF  Curv_if1
!
              END IF  Fptot_if1
            END DO  Row_length_do92
          END DO  Rows_do92
        END DO  IDir_do2a
!
        Lchop2_if1: IF (L_CHOP2)  THEN
! ----------------------------------------------------------------------+-------
!       Process low wavenumber contribution: evaluate intersect mNX
! ----------------------------------------------------------------------+-------
          NCHOP2 = 0
!
          IDir_do2b: DO JDIR=1,IDIR
            Rows_do93: DO J=1,ROWS
              Row_length_do93: DO I=1,ROW_LENGTH
                IF (MGUESS_a(I,J,JLEV,JDIR) >  0.0)  THEN
                  NCHOP2 = NCHOP2 + 1
!
                  INDEXJD(NCHOP2) = JDIR
                  INDEXJ(NCHOP2)  = J
                  INDEXI(NCHOP2)  = I
                END IF
              END DO  Row_length_do93
            END DO  Rows_do93
          END DO  IDir_do2b
!         NVIEW(JLEV) = NCHOP2
!         NSPIRA = 0
!
          Nchop2_do1: DO I=1,NCHOP2
! ----------------------------------------------------------------------+-------
!               Chop Type IIa : / Full solution required for mNlX
!          or   Chop Type IIIa: \
! ----------------------------------------------------------------------+-------
            NNJD = INDEXJD(I)
            NNJ  = INDEXJ(I)
            NNI  = INDEXI(I)
!
            FPLUS  = MGUESS_a(NNI,NNJ,JLEV,NNJD)**SSPTT
            GPLUS  = INTERCEPT1(NNI,NNJ,JLEV,NNJD)
!           FMINUS = MNLMIN**SSPTT    Defined as a constant
            GMINUS = MINTERCEPT(NNI,NNJ,JLEV,NNJD)
            ATTE_C(I) = ATTENUATION(NNI,NNJ,JLEV,NNJD)
            CURV_C(I) = CURVATURE(NNI,NNJ,JLEV,NNJD)
!
            FTERM = ( ((FMINUS / ATTE_C(I))**RTTM2) - 1.0 ) / CURV_C(I)
            GTERM = GMINUS**RSSPTT
            L_FTHENG(I) = .FALSE.
!
            Curv_if2: IF (CURVATURE(NNI,NNJ,JLEV,NNJD) >  ASECP)  THEN
! ----------------------------------------------------------------------+-------
!           Positive Doppler Shift
! ----------------------------------------------------------------------+-------
              WGTN(I) = 0.0
            ELSE
! ----------------------------------------------------------------------+-------
!           Negative Doppler Shift
! ----------------------------------------------------------------------+-------
              WGTN(I) = 1.0 - MWEIGHT
!
              IF (FPLUS <= GMINUS  .AND.  GPLUS >  FMINUS)  THEN
                FTERM = (((FPLUS / ATTE_C(I))**RTTM2) - 1.0)/ CURV_C(I)
                GTERM = GPLUS**RSSPTT
                L_FTHENG(I) = (GTERM  <   FTERM)
!
              ELSE IF (FPLUS >  GMINUS  .AND.  GPLUS <= FMINUS)  THEN
                L_FTHENG(I) = (GTERM >= FTERM)
!
              ELSE IF (FPLUS <= GMINUS  .AND.  GPLUS <= FMINUS)  THEN
                L_FTHENG(I) = .TRUE.
!
!             ELSE Use default settings
              END IF
            END IF  Curv_if2
!
            IF (L_FTHENG(I))  THEN
!             NSPIRA = NSPIRA + 1
              MNLX(I,0) = FTERM
            ELSE
              MNLX(I,0) = GTERM
            END IF
          END DO  Nchop2_do1
!         NVIEW2(JLEV) = NSPIRA
!
          Jwhile_do2: DO JWHILE=0,MAXWHILE-1
            Nchop2_do2: DO I=1,NCHOP2
!
              IF (L_FTHENG(I))  THEN
! ----------------------------------------------------------------------+-------
!           Obtain m_n+1 from g_n+1  = f_n (m_n)
! ----------------------------------------------------------------------+-------
                MNLX(I,JWHILE+1) = (                                    &
     &           (((MNLX(I,JWHILE)**SSPTT) / ATTE_C(I))**RTTM2) - 1.0 ) &
     &           / CURV_C(I)
              ELSE
! ----------------------------------------------------------------------+-------
!           Obtain m_n+1 from f_n+1  = g_n (m_n)
! ----------------------------------------------------------------------+-------
                MNLX(I,JWHILE+1) = ( (ATTE_C(I) *                       &
     &          ((1.0 + (CURV_C(I) * MNLX(I,JWHILE)))**TTM2))**RSSPTT )
              END IF
!
              MNLX(I,JWHILE+1) = ((1.0 - WGTN(I)) * MNLX(I,JWHILE+1)) + &
     &                                  (WGTN(I)  * MNLX(I,JWHILE))
!
            END DO  Nchop2_do2
          END DO  Jwhile_do2
!
!CDIR NODEP
          Nchop2_do3: DO I=1,NCHOP2
            NNJD = INDEXJD(I)
            NNJ  = INDEXJ(I)
            NNI  = INDEXI(I)
!
            FPTOT(NNI,NNJ,JLEV,NNJD) = FPTOT(NNI,NNJ,JLEV,NNJD) +       &
     &     (FPFAC(NNI,NNJ,JLEV,NNJD) * ( ((MNLX(I,MAXWHILE)**SSP1) *    &
     &      ( SSPTT + (SSP1 * MNLX(I,MAXWHILE) * CURV_C(I)) )) - HEAD_CHOP2A ))
!
          END DO  Nchop2_do3
        END IF  Lchop2_if1
!
        IDir_do2c: DO JDIR=1,IDIR
          Rows_do10: DO J=1,ROWS
            Row_length_do10: DO I=1,ROW_LENGTH
!-----------------------------------------------------------------------+-------
!         Now correct pseudomomentum flux in the evolved spectrum if the
!         new value is non-physical (pseudomomentum flux cannot increase
!         with altitude)
!-----------------------------------------------------------------------+-------
              FPTOT(I,J,JLEV,JDIR) =                                    &
     &          MIN(FPTOT(I,J,JLEV,JDIR), FPTOT(I,J,JLEV-1,JDIR))
!
            END DO  Row_length_do10
          END DO  Rows_do10
        END DO  IDir_do2c
!
      END DO  Levels_do92
!
! ----------------------------------------------------------------------+-------
! 4.5   If choosing Opaque Upper Boundary set fluxes to zero at top
! ----------------------------------------------------------------------+-------
      IF (L_USSP_OPAQUE) THEN
        IDir_do3: DO JDIR=1,IDIR
          Rows_do12: DO J=1,ROWS
            Row_length_do12: DO I=1,ROW_LENGTH
              FPTOT(I,J,LEVELS,JDIR) =  0.0
            END DO  Row_length_do12
          END DO  Rows_do12
        END DO  IDir_do3
      ENDIF
!
! ----------------------------------------------------------------------+-------
! 5.0   Compute vertical divergence of pseudomomentum flux.
! ----------------------------------------------------------------------+-------
      IDir_do4: DO JDIR=1,IDIR
        Levels_do14: DO JLEV=ILAUNCH+1,LEVELS
          Rows_do14: DO J=1,ROWS
            Row_length_do14: DO I=1,ROW_LENGTH
!           Pseudomomentum flux
              G_G(I,J,JLEV,JDIR) =                                      &
     &         g * (FPTOT(I,J,JLEV,JDIR) - FPTOT(I,J,JLEV-1,JDIR)) /    &
     &          (P_LAYER_BOUNDARIES(I,J,JLEV) - P_LAYER_BOUNDARIES(I,J,JLEV-1))
              G_X(I,J,JLEV) = G_X(I,J,JLEV) + G_G(I,J,JLEV,JDIR) * COSPHI(JDIR)
              G_Y(I,J,JLEV) = G_Y(I,J,JLEV) + G_G(I,J,JLEV,JDIR) * SINPHI(JDIR)
            END DO  Row_length_do14
          END DO  Rows_do14
        END DO  Levels_do14
      END DO  IDir_do4
!
! ----------------------------------------------------------------------+-------
!     WRITE(6,*) ' Level of Launching : ',ILAUNCH
!     WRITE(6,*) ' Number of Specials : ',NVIEW
!     WRITE(6,*) ' Number of SpiralAs : ',NVIEW2
! ----------------------------------------------------------------------+-------
!
      Levels_do15: DO JLEV=ILAUNCH,LEVELS
! ----------------------------------------------------------------------+-------
! 5.1   Wind and temperature increments from wave dissipation
! ----------------------------------------------------------------------+-------
        Rows_do15: DO J=1,ROWS
          Row_length_do15: DO I=1,ROW_LENGTH
            UONP_HALO(I,J,JLEV) = TIMESTEP * G_X(I,J,JLEV)
            VONP_HALO(I,J,JLEV) = TIMESTEP * G_Y(I,J,JLEV)
          END DO  Row_length_do15
        END DO  Rows_do15
      END DO  Levels_do15
!
! ----------------------------------------------------------------------+-------
!     Put U,V increments onto U,V grid after initialising increments.
! ----------------------------------------------------------------------+-------
      UINC(:,:,:)=0.
      VINC(:,:,:)=0.
! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS(UONP_HALO(1-off_x,1-off_y,1)                     &
     &                 ,ROW_LENGTH,ROWS,LEVELS,                         &
     &                OFF_X,OFF_Y,FLD_TYPE_P,.FALSE.)
! DEPENDS ON: fill_external_halos
      Call FILL_EXTERNAL_HALOS(UONP_HALO(1-off_x,1-off_y,1),            &
     &                         row_length,rows,levels,off_x,off_y)
! DEPENDS ON: p_to_u
      CALL P_TO_U(UONP_HALO(1-off_x,1-off_y,1),ROW_LENGTH,ROWS,LEVELS,  &
     &            OFF_X,OFF_Y,UINC)
!
! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS(VONP_HALO(1-off_x,1-off_y,1)                     &
     &                ,ROW_LENGTH,ROWS,LEVELS,                          &
     &                OFF_X,OFF_Y,FLD_TYPE_P,.FALSE.)
! DEPENDS ON: fill_external_halos
      Call FILL_EXTERNAL_HALOS(VONP_HALO(1-off_x,1-off_y,1),            &
     &                         row_length,rows,levels,off_x,off_y)
! DEPENDS ON: p_to_v
      CALL P_TO_V(VONP_HALO(1-off_x,1-off_y,1),                         &
     &            ROW_LENGTH,ROWS,NROWS,LEVELS,                         &
     &            OFF_X,OFF_Y,VINC)
!
! ----------------------------------------------------------------------+-------
! Add increments to wind and temperature
! ----------------------------------------------------------------------+-------
      DO K=ILAUNCH,LEVELS
       DO J=1,ROWS
        DO I=1,ROW_LENGTH
          R_U(I,J,K)=R_U(I,J,K)+UINC(I,J,K)
        ENDDO
       ENDDO
       DO J=1,NROWS
        DO I=1,ROW_LENGTH
          R_V(I,J,K)=R_V(I,J,K)+VINC(I,J,K)
        ENDDO
       ENDDO
      ENDDO
!
! ----------------------------------------------------------------------+-------
!     Diagnostic output : Northward Flux of Horiz Pseudomom (on theta levels)
! ----------------------------------------------------------------------+-------
      GWspec_Flux_if1n: IF(GWSPEC_NFLUX_ON) THEN
        Levels_do20n: DO K=1,LEVELS
          Rows_do20n: DO J=1,ROWS
            Row_length_do20n: DO I=1,ROW_LENGTH
              WORK_HALO(I,J,K) = FPTOT(I,J,K,1)
            END DO  Row_length_do20n
          END DO  Rows_do20n
        END DO  Levels_do20n
!
! DEPENDS ON: swap_bounds
        CALL SWAP_BOUNDS(                                               &
     &    WORK_HALO(1-off_x,1-off_y,1), ROW_LENGTH, ROWS,               &
     &    LEVELS,OFF_X,OFF_Y,FLD_TYPE_P,.FALSE.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,1),          &
     &                           row_length,rows,levels,off_x,off_y)
! DEPENDS ON: p_to_v
        CALL P_TO_V(                                                    &
     &    WORK_HALO(1-off_x,1-off_y,1), ROW_LENGTH, ROWS, NROWS,        &
     &    LEVELS, OFF_X, OFF_Y, GWSPEC_NFLUX(1,1,1))
      END IF  GWspec_Flux_if1n
!
! ----------------------------------------------------------------------+-------
!     Diagnostic output : Westward Flux of Horiz Pseudomom (on theta levels)
! ----------------------------------------------------------------------+-------
      GWspec_Flux_if1w: IF(GWSPEC_WFLUX_ON) THEN
        Levels_do20w: DO K=1,LEVELS
          Rows_do20w: DO J=1,ROWS
            Row_length_do20w: DO I=1,ROW_LENGTH
              WORK_HALO(I,J,K) = FPTOT(I,J,K,2)
            END DO  Row_length_do20w
          END DO  Rows_do20w
        END DO  Levels_do20w
!
! DEPENDS ON: swap_bounds
        CALL SWAP_BOUNDS(                                               &
     &    WORK_HALO(1-off_x,1-off_y,1), ROW_LENGTH, ROWS,               &
     &    LEVELS, OFF_X, OFF_Y, FLD_TYPE_P, .FALSE.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,1),          &
     &                           row_length,rows,levels,off_x,off_y)
! DEPENDS ON: p_to_u
        CALL P_TO_U(                                                    &
     &    WORK_HALO(1-off_x,1-off_y,1), ROW_LENGTH, ROWS,               &
     &    LEVELS, OFF_X, OFF_Y, GWSPEC_WFLUX(1,1,1))
      END IF  GWspec_Flux_if1w
!
! ----------------------------------------------------------------------+-------
!     Diagnostic output : Southward Flux of Horiz Pseudomom (on theta levels)
! ----------------------------------------------------------------------+-------
      GWspec_Flux_if1s: IF(GWSPEC_SFLUX_ON) THEN
        Levels_do20s: DO K=1,LEVELS
          Rows_do20s: DO J=1,ROWS
            Row_length_do20s: DO I=1,ROW_LENGTH
              WORK_HALO(I,J,K) = FPTOT(I,J,K,3)
            END DO  Row_length_do20s
          END DO  Rows_do20s
        END DO  Levels_do20s
!
! DEPENDS ON: swap_bounds
        CALL SWAP_BOUNDS(                                               &
     &    WORK_HALO(1-off_x,1-off_y,1), ROW_LENGTH, ROWS,               &
     &    LEVELS, OFF_X, OFF_Y, FLD_TYPE_P, .FALSE.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,1),          &
     &                           row_length,rows,levels,off_x,off_y)
! DEPENDS ON: p_to_v
        CALL P_TO_V(                                                    &
     &    WORK_HALO(1-off_x,1-off_y,1), ROW_LENGTH, ROWS, NROWS,        &
     &    LEVELS, OFF_X, OFF_Y, GWSPEC_SFLUX(1,1,1))
      END IF  GWspec_Flux_if1s
!
! ----------------------------------------------------------------------+-------
!     Diagnostic output : Eastward Flux of Horiz Pseudomom (on theta levels)
! ----------------------------------------------------------------------+-------
      GWspec_Flux_if1e: IF(GWSPEC_EFLUX_ON) THEN
        Levels_do20e: DO K=1,LEVELS
          Rows_do20e: DO J=1,ROWS
            Row_length_do20e: DO I=1,ROW_LENGTH
              WORK_HALO(I,J,K) = FPTOT(I,J,K,4)
            END DO  Row_length_do20e
          END DO  Rows_do20e
        END DO  Levels_do20e
!
! DEPENDS ON: swap_bounds
        CALL SWAP_BOUNDS(                                               &
     &    WORK_HALO(1-off_x,1-off_y,1), ROW_LENGTH, ROWS,               &
     &    LEVELS, OFF_X, OFF_Y, FLD_TYPE_P, .FALSE.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,1),          &
     &                           row_length,rows,levels,off_x,off_y)
! DEPENDS ON: p_to_u
        CALL P_TO_U(                                                    &
     &    WORK_HALO(1-off_x,1-off_y,1), ROW_LENGTH, ROWS,               &
     &    LEVELS, OFF_X, OFF_Y, GWSPEC_EFLUX(1,1,1))
      END IF  GWspec_Flux_if1e
!
! ----------------------------------------------------------------------+-------
!     Diagnostic output : Acceleration of Zonal Wind (on rho levels)
! ----------------------------------------------------------------------+-------
      IF(GWSPEC_EWACC_ON) THEN
       DO K=1,LEVELS
        DO J=1,ROWS
         DO I=1,ROW_LENGTH
          WORK_HALO(I,J,K) = G_X(I,J,K)
         ENDDO
        ENDDO
       ENDDO
! DEPENDS ON: swap_bounds
       CALL SWAP_BOUNDS(                                                &
     &      WORK_HALO(1-off_x,1-off_y,1), ROW_LENGTH, ROWS,             &
     &      LEVELS, OFF_X, OFF_Y, FLD_TYPE_P, .FALSE.)
! DEPENDS ON: fill_external_halos
       CALL FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,1),           &
     &                   row_length,rows,levels,off_x,off_y)
! DEPENDS ON: p_to_u
       CALL P_TO_U(                                                     &
     &      WORK_HALO(1-off_x,1-off_y,1), ROW_LENGTH, ROWS,             &
     &      LEVELS, OFF_X, OFF_Y, GWSPEC_EWACC)
      ENDIF
! ----------------------------------------------------------------------+-------
!     Diagnostic output : Acceleration of Meridional Wind (on rho levels)
! ----------------------------------------------------------------------+-------
      IF(GWSPEC_NSACC_ON) THEN
       DO K=1,LEVELS
        DO J=1,NROWS
         DO I=1,ROW_LENGTH
          WORK_HALO(I,J,K) = G_Y(I,J,K)
         ENDDO
        ENDDO
       ENDDO
! DEPENDS ON: swap_bounds
       CALL SWAP_BOUNDS(                                                &
     &      WORK_HALO(1-off_x,1-off_y,1), ROW_LENGTH, ROWS,             &
     &      LEVELS, OFF_X, OFF_Y, FLD_TYPE_P, .FALSE.)
! DEPENDS ON: fill_external_halos
       CALL FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,1),           &
     &                   row_length,rows,levels,off_x,off_y)
! DEPENDS ON: p_to_v
       CALL P_TO_V(                                                     &
     &      WORK_HALO(1-off_x,1-off_y,1), ROW_LENGTH, ROWS,             &
     &      NROWS,                                                      &
     &      LEVELS,OFF_X,OFF_Y,GWSPEC_NSACC)
      ENDIF
      RETURN
!
      END SUBROUTINE GW_USSP
#endif
