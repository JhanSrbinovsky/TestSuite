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
! Purpose: Main chemistry driver routine
!
!     Chemistry driver routine. If necessary, model concentrations are
!     converted from vmr to number density. If the chemistry is to be
!     "process-split" from the transport, the chemistry tendencies are
!     integrated by the chosen method and an average chemistry tendency
!     over the model timestep is returned. If the chemistry is not to
!     be integrated, instantaneous chemistry tendencies are returned.
!     The tracer tendencies returned to the calling routine are the
!     final values held in the chemistry.
!
!     Note: It is important for conservation that the concentrations
!     passed to this routine are non-negative.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from UKCA_CHEMISTRY_CTL
!
! Current code owner: Glenn Carver/Colin Johnson
!                     Oliver Wild
!
!     Arguments:
!        cdot        - Tracer tendencies due to chemistry.
!        ftr         - Tracer concentrations.
!        pp          - Pressure (Nm-2).
!        pt          - Temperature (K).
!        pq          - Water vapor field (vmr).
!        nlev        - Model level
!        dryrt       - Dry deposition rates (s-1)
!        wetrt       - Wet deposition rates (s-1)
!        n_points    - No. of points calculations be done.
!
!     Method
!     ------
!     Since the heterogeneous rates may depend on species
!     concentrations, the call to hetero is made inside the
!     loop over chemistry sub-steps.
!
!     Photolysis rates may be computed at frequency set by the
!     variable nfphot and so the call to photol is also inside
!     the chemical sub-step loop.
!
!     Externals
!     ---------
!     asad_bimol  - Calculates bimolecular rate coefficients.
!     asad_trimol - Calculates trimolecular rate coefficients.
!     ukca_photol - Calculates photolysis rate coefficients.
!     ukca_drydep - Calculates dry deposition rates.
!     ukca_wetdep - Calculates wet deposition rates.
!     asad_emissn - Calculates emission rates.
!     asad_totnud - Calculates total number densities.
!     asad_impact - IMPACT time integration scheme.
!     asad_hetero - Calculates heterogeneous rate coefficients.
!     asad_posthet- Performs housekeeping after heterogeneous chemistry.
!     asad_ftoy   - Partitions families.
!     asad_diffun - Calculates chemistry tendencies.
!
!     Local variables
!     ---------------
!     ifam       Family index of in/out species.
!     itr        Tracer index of in/out species.
!     iodd       Number of odd atoms in in/out species.
!     gfirst     .true. when the species need to be
!                initialised on the first chemical step.
!     gphot      .true. if photol needs to be called.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_CDRIVE(cdot,ftr,pp,pt,pq,nlev,dryrt,           &
                               wetrt,prt,n_points,stratflag)
        USE ASAD_MOD
        USE UKCA_HETERO_MOD
        IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "parvars.h"
#include "typsize.h"

! Subroutine interface
        INTEGER, INTENT(IN) :: n_points               ! No of points
        INTEGER, INTENT(IN) :: nlev                   ! Model level

        REAL, INTENT(IN) :: prt(n_points,jppj)        ! Photolysis rates
        REAL, INTENT(IN) :: pq(n_points)              ! Water vapour
        REAL, INTENT(IN) :: dryrt(n_points,jpdd)      ! Dry dep rates
        REAL, INTENT(IN) :: wetrt(n_points,model_levels,jpdw) ! Wet dep rates
        REAL, INTENT(IN) :: pp(theta_field_size)       ! Pressure
        REAL, INTENT(IN) :: pt(theta_field_size)       ! Temperature
        LOGICAL, INTENT(IN) :: stratflag(theta_field_size) ! Strat indicator

        REAL, INTENT(INOUT) :: ftr(theta_field_size,jpctr) ! Tracer concs

        REAL, INTENT(OUT) :: cdot(theta_field_size,jpctr)  ! Tracer tendencies

!       Local variables

        INTEGER :: jtr                                ! Loop variable
        INTEGER :: jl                                 ! Loop variable
        INTEGER :: js                                 ! Loop variable
        INTEGER :: nl
        INTEGER :: ifam
        INTEGER :: itr
        INTEGER :: iodd

        LOGICAL :: gfirst
        LOGICAL :: gphot

        CHARACTER(len=72) :: cmessage          ! Error message

!       1.  Initialise variables and arrays

!       1.1   Clear tendencies to avoid contributions from levels
!             on which no chemistry is performed

        DO jtr = 1, jpctr
          DO jl = 1, n_points
            cdot(jl,jtr) = 0.0
          ENDDO
        ENDDO

!       1.2  Copy pressure and temperature to common

        DO jl = 1, n_points
          p(jl) = pp(jl)
          t(jl) = pt(jl)
        ENDDO

!       1.2.1 Copy water vapor to common

        DO jl = 1, n_points
          wp(jl) = pq(jl)
        ENDDO

!       2.  Calculate total number densities

! DEPENDS ON: asad_totnud
        CALL asad_totnud(n_points)

!       3.  Read model tracer concentrations into working array,
!           and if necessary, convert vmr to number densities

        IF ( lvmr ) THEN
          do jtr = 1, jpctr
            do jl  = 1, n_points
              ftr(jl,jtr) = ftr(jl,jtr) * tnd(jl)
              f(jl,jtr)   = ftr(jl,jtr)
            ENDDO
          ENDDO
        ELSE
          do jtr = 1, jpctr
            do jl  = 1, n_points
              f(jl,jtr)   = ftr(jl,jtr)
            ENDDO
          ENDDO
        ENDIF

!       4.  Calculate reaction rate coefficients
!           --------- -------- ---- ------------

! DEPENDS ON: asad_bimol
        CALL asad_bimol (n_points)
! DEPENDS ON: asad_trimol
        CALL asad_trimol(n_points)

!       5.  Calculate deposition and emission rates
!           --------- ---------- --- -------- -----

! DEPENDS ON: ukca_wetdep
        IF ( ndepw /= 0 ) CALL UKCA_WETDEP(nlev, wetrt, n_points)
! DEPENDS ON: ukca_drydep
        IF ( ndepd /= 0 ) CALL UKCA_DRYDEP(nlev, dryrt, n_points)
! DEPENDS ON: asad_emissn
        IF ( nemit /= 0 ) CALL asad_emissn

!       6.  Integrate chemistry by chosen method. Otherwise,
!           simply calculate tendencies due to chemistry
!           ------ --------- ---------- --- -- ---------

        gphot = .true.
        IF ( method /= 0 ) THEN

          DO jsubs = 1, ncsteps

            gfirst = jsubs  ==  1
            gphot  = gfirst
            IF ( nfphot /= 0 .AND. .NOT.gfirst )                       &
                          gphot = mod(jsubs-1,nfphot) == 0

!           ---------------------------------------------------
!           NON-STIFF integrators take values in the range 1-9.
!           ---------------------------------------------------

!           6.1  IMPACT integration: first compute heterogeneous
!                and photolysis rates, species and tendencies.
!                ===============================================

            nl = n_points
            IF ( method == 1 ) THEN
! DEPENDS ON: ukca_photol
              IF ( gphot ) CALL UKCA_PHOTOL(prt,n_points)
              IF (L_ukca_het_psc)  CALL ukca_hetero(n_points,stratflag)
! DEPENDS ON: asad_ftoy
              CALL asad_ftoy( gfirst, nitfg, n_points )
! DEPENDS ON: asad_diffun
              CALL asad_diffun( nl )
! DEPENDS ON: asad_jac
              CALL asad_jac( n_points )
! DEPENDS ON: asad_impact
              CALL asad_impact( n_points )

!             6.2.  Quasi-steady state scheme.
!             ================================

            ELSEIF ( method == 2 ) THEN
              cmessage='QSSA not in UM6.5 build'
! DEPENDS ON: ereport
              CALL EREPORT('ASAD_CDRIVE',1,cmessage)

!              6.3   Sparse Newton-Raphson solver
!              ==================================

            ELSEIF ( method == 3 ) THEN
! DEPENDS ON: ukca_photol
              IF ( gphot ) CALL UKCA_photol(prt,n_points)
              IF (L_ukca_het_psc)  CALL UKCA_HETERO(n_points,stratflag)
! DEPENDS ON: asad_ftoy
              CALL asad_ftoy( gfirst, nitfg, n_points )
! DEPENDS ON: asad_spmjpdriv
              CALL ASAD_SPMJPDRIV(nlev,mype,n_points)

!              6.5   Backward Euler solver
!              ===========================

            ELSEIF ( method == 5 ) then
! DEPENDS ON: ukca_photol
              IF ( gphot ) CALL UKCA_PHOTOL(prt,n_points)
              IF (L_ukca_het_psc)  CALL UKCA_HETERO(n_points,stratflag)
! DEPENDS ON: asad_ftoy
              CALL ASAD_FTOY( gfirst, nitfg, n_points )
! DEPENDS ON: asad_bedriv
              CALL ASAD_BEDRIV(nlev,mype,n_points)

!           -------------------------------------------------
!           STIFF integrators take values in the range 10-19.
!           -------------------------------------------------

!           6.10  NAG BDF stiff integrator.

            ELSEIF ( method == 10 ) THEN
              cmessage='NAG BDF not in UM6.5 build'
! DEPENDS ON: ereport
              CALL EREPORT('ASAD_CDRIVE',1,cmessage)

!             6.11  SVODE ODE stiff integrator from NETLIB.

            ELSEIF ( method == 11 ) THEN
              cmessage='SVODE not in UM6.5 build'
! DEPENDS ON: ereport
              CALL EREPORT('ASAD_CDRIVE',1,cmessage)

            ENDIF

!           6.12  Do any final work for the heterogeneous
!                 chemistry before the end of the time loop.
!                 =========================================

!            IF ( L_ukca_het_psc ) CALL ukca_solidphase(n_points)
! DEPENDS ON: asad_posthet
            IF ( L_ukca_het_psc ) CALL asad_posthet

          ENDDO        ! End of looping over chemical timesteps

        ELSE           ! Method is equal to zero

!       6.99  Not integrating: just compute tendencies.

          method = 0
! DEPENDS ON: ukca_photol
          IF ( gphot ) CALL ukca_photol(prt,n_points)
          IF ( L_ukca_het_psc )  CALL ukca_hetero(n_points,stratflag)
! DEPENDS ON: asad_ftoy
          CALL asad_ftoy( .true., nit0, n_points )
! DEPENDS ON: asad_diffun
          CALL asad_diffun( nl )
          IF ( L_ukca_het_psc ) CALL ukca_solidphase( n_points )
! DEPENDS ON: asad_posthet
          IF ( L_ukca_het_psc ) CALL asad_posthet

        ENDIF        ! End of IF statement for method


!       7.  Determine concentrations and tendencies to be returned to
!           the model depending on whether or not the chemistry has
!           been integrated  -- -- ------- -- --- --- --------- ---
!           ---- ----------

!       7.1  Obtain model family concentrations by subtracting
!            concentrations of in/out species from ASAD families

        DO js = 1, jpspec
          IF ( ctype(js) == jpif ) THEN
            ifam = moffam(js)
            itr  = madvtr(js)
            iodd = nodd(js)
            DO jl = 1, n_points
              IF ( linfam(jl,itr) ) f(jl,ifam) =                       &
                                    f(jl,ifam) - iodd*f(jl,itr)
            ENDDO
          ENDIF
        ENDDO

!       7.2  Returned values of concentration and chemical tendency

        DO jtr = 1, jpctr
          IF ( method /= 0 ) THEN
            DO jl = 1, n_points
              cdot(jl,jtr) = ( f(jl,jtr)-ftr(jl,jtr)) / (cdt*ncsteps)
              ftr(jl,jtr)  = f(jl,jtr)
            ENDDO
          ELSE
            DO jl = 1, n_points
              cdot(jl,jtr) = fdot(jl,jtr)
              ftr(jl,jtr)  = f(jl,jtr)
            ENDDO
          ENDIF
        ENDDO


!       8.  If necessary, convert from number densities back to vmr
!           -- ---------- ------- ---- ------ --------- ---- -- ---

        IF ( lvmr ) THEN
          DO jtr = 1, jpctr
            DO jl = 1, n_points
              ftr(jl,jtr)  = ftr(jl,jtr)  / tnd(jl)
              cdot(jl,jtr) = cdot(jl,jtr) / tnd(jl)
            ENDDO
          ENDDO
        ENDIF

        RETURN
        END SUBROUTINE ASAD_CDRIVE
#endif
