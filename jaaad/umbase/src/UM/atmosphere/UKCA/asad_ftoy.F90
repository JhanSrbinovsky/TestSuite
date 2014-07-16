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
! Purpose: Sets the species concentrations prior to calculating the
!          chemistry tendencies. This includes partitioning families.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CDRIVE and ASAD_IMPACT
!
! Current code owner: Glenn Carver/Colin Johnson
!                     Oliver Wild
!
!     Interface
!     ---------
!     Must always be called before diffun.
!
!     Arguments:
!        ofirst     Should be .true. if this is first call to ftoy
!                   during the current model/dynamical timestep.
!        iter       Maximum number of iterations.
!
!     Method
!     ------
!    After initialisation, the family member concentrations and
!    steady state species are iterated until they converge to a
!    solution. Firstly, the ratios for family members and steady
!    state concentrations for the non-family steady state species
!    are computed. Then, the family member concentrations are
!    evaluated using the ratios. A convergence test is performed
!    for all species at all spatial points. If convergence is not
!    achieved for any species at any point, the iterations continue
!    until either convergence is achieved or the maximum number of
!    iterations is reached.
!
!    In/out species are determined to be in or out of the family
!    on the first call to ftoy during a model timestep according
!    to their lifetime compared with some threshold value. If the
!    species is found to be 'in' the family, its concentration
!    is added to the family f concentration before the sum
!    is partitioned amongst the members. The in/out test is
!    performed at all spatial points. This species then stays in
!    the family for the rest of the call to cdrive. The species is
!    subtracted from the family at the end of cdrive.
!
!    If the species self reacts, the steady state concentration
!    is found by solving a quadratic. Otherwise, it is simply
!    found by dividing the production rate by the loss rate.
!    The self reacting term, qa, has previously been calculated
!    in fyself. Note qa is +ve.
!
!    The ratio of a minor family member to the major member is
!    determined by dividing the steady state concentration of the
!    minor member from fystst by the family member concentration:
!                 R1m = Y1/Ym, R2m = Y2/Ym ....
!    The ratio of the major member to the family as a whole is
!    determined from the minor member ratios weighted by the number
!    of odd atoms in each species:
!                 Rmf = 1 / ( Nm + N1*R1m +N2*R2m + .....)
!    The concentrations are then found from:
!                 Ym = Rmf*Z, Y1 = R1m*Ym, Y2 = R2m*Ym, .....
!    where Z is the remaining family concentration after the
!    appropriate in/out species concentrations have been subtracted.
!
!    Although ASAD supports the user of deposition and emissions,
!    we do not include these terms in the calculation of the ratios
!    even though they represent a loss or production from or in to
!    the family. This is because, family members are assumed to be
!    in steady state, generally taken to be a 'chemical' steady
!    state, one that is controlled by the concentrations of other
!    species. If for instance, a strong emission source was present,
!    for NO, this would give a high bogus value for the NO ratio,
!    which would be wrong because by assuming steady state we assume
!    that the extra NO has already reacted to reach a new equilibrium
!    with its family members.
!
!    Externals
!    ---------
!    fyinit         Sets initial species concentrations.
!    fyself         Calculates self-reacting terms.
!    fyfixr         Calculates family member concentrations
!                   using ratios computed on a previous call.
!    fystst         Calculates concentration of a species
!                   in steady state.
!    prls           Calculates production & loss terms.
!
!    Local variables
!    ---------------
!    ifam           Index of family to which species belongs.
!    imaj           Index of major member of family to which
!                   species belongs.
!    itr            Index of model tracer to which species
!                   corresponds.
!    gconv          .true. if convergence has been achieved.
!    gonce          .true. for only the first call to this routine.
!    zthresh        Threshold value of loss rate to determine
!                   whether in/out species are in or out of family.
!    zy             Value of y on previous iteration.
!    ilstmin        List of species which are either steady state,
!                   'SS', family members ('FM' and 'FT') but EXCLUDING
!                   major family species.
!    istmin         No. of entries in ilstmin.
!    ilft           List of species of type 'FT'.
!    ift            No. of entries in ilft.
!    zb, zc, zd     Variables used in the quadratic. Note that
!                   zb is really '-b' in the quadratic.
!    zb             The loss rate less self reacting terms.
!    zc             The production rate.
!    zd             "b^2 - 4*a*c"
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_FTOY(ofirst,iter,n_points)

        USE ASAD_MOD
        IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "parvars.h"
#include "typsize.h"

        INTEGER, INTENT(IN)    :: n_points  ! No of spatial points
        INTEGER, INTENT(INOUT) :: iter      ! Max no of iterations

        LOGICAL, INTENT(IN) :: ofirst

!       Local variables

        INTEGER       :: j           ! Loop variable
        INTEGER       :: jit         ! Loop variable
        INTEGER       :: jl          ! Loop variable
        INTEGER       :: js          ! Loop variable
        INTEGER       :: ifam
        INTEGER       :: imaj
        INTEGER       :: nl
        INTEGER       :: istart
        INTEGER       :: iend
        INTEGER       :: iodd
        INTEGER       :: itr

        INTEGER, SAVE :: istmin
        INTEGER, SAVE :: ift
!        INTEGER       :: ilstmin(jpspec)      ! Now in ASAD_MOD
!        INTEGER       :: ilft(jpspec)         ! Now in ASAD_MOD

        CHARACTER (LEN=72) :: cmessage     ! Error message

        REAL          :: zthresh
        REAL          :: sl
        REAL          :: zy(theta_field_size,jpspec)
        REAL          :: zb(theta_field_size)
        REAL          :: zc(theta_field_size)
        REAL          :: zd(theta_field_size)

        LOGICAL       :: gconv
        LOGICAL, SAVE :: gonce  = .true.
        LOGICAL, SAVE :: gdepem = .false.

!       1.  Initialise.
!           -----------

        IF ( gonce ) THEN
          gonce  = .false.
          istmin = 0
          ift    = 0
          DO j = 1, jpspec
            ilstmin(j) = 0
            ilft(j) = 0
          ENDDO

          DO j = 1, nstst
            js = nlstst(j)
            IF ( ctype(js) == jpfm ) THEN
              ifam = moffam(js)
              imaj = majors(ifam)
              IF ( imaj /= js ) THEN
                istmin          = istmin + 1
                ilstmin(istmin) = js
              ENDIF
            ELSEIF ( ctype(js) == jpif ) THEN
              istmin          = istmin + 1
              ilstmin(istmin) = js
              ift             = ift + 1
              ilft(ift)       = js
            ELSEIF ( ctype(js) == jpna ) THEN
              istmin          = istmin + 1
              ilstmin(istmin) = js
            ELSE
              WRITE (6,*) '**** ASAD ERROR in ASAD_FTOY!! '
              WRITE (6,*) 'ASAD_FTOY found an unexpected species type'
              WRITE (6,*) 'in the species list nlstst ', nlstst
              cmessage='Found unexpected species type'
! DEPENDS ON: ereport
              CALL EREPORT('ASAD_FTOY',js,cmessage)
            ENDIF
          ENDDO
        ENDIF   ! End of IF (gonce) statement

!       1.1 Set concentrations/initialise species

! DEPENDS ON: asad_fyinit
        CALL asad_FYINIT(ofirst,n_points)

!       1.2 If there are no family or steady species then exit

        IF ( nstst == 0 ) RETURN

!       1.3 Initialise local variables and do sanity checks.

        nl = n_points

        DO j = 1, nstst
          js = nlstst(j)
          DO jl = 1, n_points
            zy(jl,js) = 0.0
          ENDDO
        ENDDO

        IF (ofirst .AND. iter < 5) THEN
          WRITE (6,*) 'WARNING:ITER TOO LOW ON FIRST CALL TO ASAD_FTOY'
          WRITE (6,*) 'RESETTING TO 5'
          iter = 5
        END IF

        zthresh = 2.0 / cdt

!       2.  Calculate self-reacting terms
!           --------- ------------- -----

! DEPENDS ON: asad_fyself
        IF ( ofirst ) CALL ASAD_FYSELF(n_points)

!       3.  Calculate family members using previous ratios
!           --------- ------ ------- ----- -------- ------

        IF (iter == 0) THEN
! DEPENDS ON: asad_fyfixr
          CALL ASAD_FYFIXR(n_points)
          RETURN
        ENDIF

!       4.  Iterate family members (including in/out species)
!           ------- ------ ------- --------------------------
!           and non-model steady-state species
!           --- --------- ------------ -------

        DO jit = 1, iter

!         check to see if a species has gone zero with a nonzero
!         production.

          IF ( jit > 1 ) THEN
            DO j = 1, nstst
              js = nlstst(j)
              DO jl = 1, n_points
                IF (abs(y(jl,js)) < peps .AND.                         &
                    prod(jl,js) > peps ) THEN
                  y(jl,js) = peps
                ENDIF
              ENDDO
            ENDDO
          ENDIF

!         4.1 Calculate production and loss terms.
!             Note that the effect of deposition (wet & dry) and
!             emissions are not included in the ratios calculation.
!             Regardless of whether the user has them turned on or not.
!             See method above.

! DEPENDS ON: asad_prls
          CALL ASAD_PRLS( nl, istmin, ilstmin, gdepem )

!         4.2 Initialise ratios

          DO j = 1, nstst
            js = nlstst(j)
            DO jl = 1, n_points
             ratio(jl,js) = 0.0
            ENDDO
          ENDDO

!         4.3  Compute species values for species in steady state
!              (minor members of family, 'SS' and 'FT' species)
!              N.B. because of the way asad computes slos, the y
!              value can go slightly negative during the quadratic
!              due to loss of precision. We forcibly fix it here.
!              Also note that we cannot permit y=0.0 during the ftoy
!              iteration because of the need to get 'sl'.

          DO j = 1, istmin
            js = ilstmin(j)
            DO jl = 1, n_points
              IF ( y(jl,js) < peps ) THEN
                sl = 0.0
              ELSE
                sl = slos(jl,js) / y(jl,js)
              ENDIF
              IF ( qa(jl,js) > peps ) THEN
                zb(jl) = sl - qa(jl,js) * y(jl,js)
                zc(jl) = prod(jl,js)
                zd(jl) = zb(jl)*zb(jl) + 4.0*qa(jl,js)*zc(jl)
                IF ( zd(jl) > 0.0 ) THEN
                  y(jl,js) = (zb(jl) - sqrt(zd(jl)))/(-2.0*qa(jl,js))
                  if ( y(jl,js)  <   0.0 ) y(jl,js) = 10.0*peps
                ELSE
                  y(jl,js) = zb(jl) / ( -2.0 * qa(jl,js) )
                ENDIF
              ELSE IF (qa(jl,js) <= peps .AND. sl > peps ) THEN
                y(jl,js) = prod(jl,js) / sl
              ELSE
                y(jl,js) = 0.0
              END IF
            ENDDO
          ENDDO

!         4.4  Now compute ratios for minor species members.
!              ** could possibly take out 1/y(imaj) and convert to '*'

          istart = nlmajmin(3)
          iend   = nlmajmin(4)
          DO j = nlmajmin(3), nlmajmin(4)
            js = nlmajmin(j)
            iodd = nodd(js)
            ifam = moffam(js)
            imaj = 0
            IF ( ifam /= 0 ) imaj = majors(ifam)
            DO jl = 1, n_points
              IF ( y(jl,imaj) > peps ) THEN
                ratio(jl,js)   = y(jl,js) / y(jl,imaj)
                ratio(jl,imaj) = ratio(jl,imaj) +                      &
                                 iodd * ratio(jl,js)
              ELSE
                ratio(jl,js) = 0.0
              END IF
            ENDDO
          ENDDO

!         4.5  Now compute ratios from species of type 'FT' if
!              applicable. If this is the first call, we also set
!              whether the species is in or out of the family.

          DO j = 1, ift
            js = ilft(j)
            iodd = nodd(js)
            ifam = moffam(js)
            itr  = madvtr(js)
            imaj = majors(ifam)
            IF ( ofirst .and. jit == 1 ) THEN
              DO jl = 1, n_points
                IF ( y(jl,js) > peps ) then
                  sl = slos(jl,js) / y(jl,js)
                ELSE
                  sl = 0.0
                ENDIF
                linfam(jl,itr) = sl  >   zthresh
                IF ( linfam(jl,itr) )                                  &
                     f(jl,ifam) = f(jl,ifam) + iodd*f(jl,itr)
              ENDDO
            END IF

            DO jl = 1, n_points
              if ( linfam(jl,itr) ) then
                if ( y(jl,imaj)  >   peps ) then
                  ratio(jl,js)   = y(jl,js) / y(jl,imaj)
                  ratio(jl,imaj) = ratio(jl,imaj) +                    &
                                   iodd * ratio(jl,js)
                else
                  ratio(jl,js) = 0.0
                endif
              else
                y(jl,js) = f(jl,itr)
              endif
            ENDDO
          ENDDO

!         4.6  Finally compute ratio of major species to family and
!              hence set the minor family members.

          istart = nlmajmin(1)
          iend   = nlmajmin(2)
          DO j = istart, iend
            js = nlmajmin(j)
            iodd = nodd(js)
            ifam = moffam(js)
            DO jl = 1, n_points
              ratio(jl,js) = 1.0 / (iodd + ratio(jl,js))
              y(jl,js)     = f(jl,ifam) * ratio(jl,js)
            ENDDO
          ENDDO

          istart = nlmajmin(3)
          iend   = nlmajmin(4)
          DO j = istart, iend
            js = nlmajmin(j)
            ifam = moffam(js)
            imaj = majors(ifam)
            IF ( ctype(js) /= jpif ) THEN
              DO jl = 1, n_points
                y(jl,js) = y(jl,imaj) * ratio(jl,js)
              ENDDO
            else
              itr = madvtr(js)
              DO jl = 1, n_points
                if ( linfam(jl,itr) ) y(jl,js) =                       &
                     y(jl,imaj) * ratio(jl,js)
              ENDDO
            ENDIF
          ENDDO

!         4.7 Test for convergence.

          gconv = .true.
          DO j = 1,nstst
            js = nlstst(j)
              DO jl = 1, n_points
                IF ( ABS(y(jl,js)-zy(jl,js)) >  ftol*y(jl,js)          &
                .AND. y(jl,js) >  pmintnd(jl) ) gconv=.false.
              ENDDO
            IF ( .NOT. gconv ) EXIT
          ENDDO
          IF (gconv) RETURN

          DO j = 1, nstst
            js = nlstst(j)
            DO jl = 1, n_points
              zy(jl,js) = y(jl,js)
            ENDDO
          ENDDO

        ENDDO    ! End of iterations


!       9. Convergence achieved or max. iterations reached.

        RETURN
        END SUBROUTINE ASAD_FTOY
#endif
