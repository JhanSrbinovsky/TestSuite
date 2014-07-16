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
! Purpose: sets or initialises species concentrations
!
!     Sets the individual species concentrations y from those of
!     the model tracers/families f. If the species concentrations
!     are later to be iterated in ftoy, then they are merely
!     initialised in this routine.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_FTOY
!
! Current code owner: Glenn Carver/Colin Johnson
!                     Oliver Wild
!
!     Arguments:
!        ofirst     Should be .true. if this is first call to FTOY
!                   for the chemistry time loop ie. ofirst should be
!                   true only when starting a new dynamical timestep.
!
!     Method
!     ------
!     If this is the first call to ftoy for the chemical time loop
!     ie. we're starting a new model/dynamical timestep, then we need
!     to initialise the species array, y. For subsequent calls, all we
!     need to do is copy the species of type 'TR' and 'FT' from the f
!     (family/tracer) array to their places in the species array, y.
!
!     When ofirst is .true., family are initialised by setting the
!     major family member to the family concentration and the minor
!     family members to a small number. Steady state species are also
!     set to a small number. Constant species (species type 'CT')
!     are also only set when .ofirst. is true, since these species are
!     never integrated (see diffun.f). Constant field species (type
!     'CF') are allowed to change every chemical time step however.
!
!     Local variables
!     ---------------
!     zcnst          VMR of fixed concentration species.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_FYINIT(ofirst,n_points)

        USE ASAD_MOD
        IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"

        INTEGER, INTENT(IN) :: n_points   ! No of spatial points

        LOGICAL, INTENT(IN) :: ofirst     ! True on first call

!       Local variables

        INTEGER, SAVE :: iss
        INTEGER, SAVE :: ict
        INTEGER, SAVE :: iftr
        INTEGER, SAVE :: icf
        INTEGER, SAVE :: istep = 0

        INTEGER       :: j           ! Loop variable
        INTEGER       :: jl          ! Loop variable
        INTEGER       :: js          ! Index
        INTEGER       :: istart
        INTEGER       :: iend
        INTEGER       :: ifam
        INTEGER       :: itr
        INTEGER       :: kstep
        INTEGER       :: inl

        REAL          :: peps10
        REAL          :: zcnst

        LOGICAL :: gonce = .true.
        CHARACTER(len=72) :: cmessage

!       1.  Species initialised only at start of dynamical step.
!           ------- ----------- ---- -- ----- -- --------- -----

        IF ( gonce ) THEN

!         Build lists for this routine. This lot should really be a
!         common variable - next version!

          gonce = .false.
          iss = 0
          ict = 0
          icf = 0
          iftr = 0
          DO js = 1, jpspec
            IF ( ctype(js) == jpna ) THEN
              iss       = iss + 1
              ilss(iss) = js
            end IF
            IF ( ctype(js) == jpco ) THEN
              ict       = ict + 1
              ilct(ict) = js
            end IF
            IF ( ctype(js) == jpcf ) THEN
              icf       = icf + 1
              ilcf(icf) = js
            end IF
            IF (ctype(js) == jpif .or.                                 &
                ctype(js) == jpsp ) THEN
              iftr        = iftr + 1
              ilftr(iftr) = js
            end IF
          ENDDO
        ENDIF           ! End of IF (gonce) statement

        IF ( ofirst ) THEN
          peps10 = 10.0 * peps

!         1.1 Ordinary family members (major then minor )

          istart = nlmajmin(1)
          iend   = nlmajmin(2)
          DO j = istart, iend
            js = nlmajmin(j)
            ifam = moffam(js)
            DO jl = 1, n_points
              y(jl,js) = f(jl,ifam)
            ENDDO
          ENDDO

          istart = nlmajmin(3)
          iend   = nlmajmin(4)
          DO j = istart, iend
            js = nlmajmin(j)
            DO jl = 1, n_points
              y(jl,js) = peps10
            ENDDO
          ENDDO

!         1.2 Non-model species in steady state

          DO j = 1, iss
            js = ilss(j)
            DO jl = 1, n_points
              y(jl,js) = peps10
            ENDDO
          ENDDO

!         1.3 Species set to constants (unaffected by remainder of ftoy

          DO j = 1, ict
            js = ilct(j)
            IF ( speci(js) == 'N2        ' ) THEN
              zcnst = fn2
            ELSE IF ( speci(js) == 'O2        ' ) THEN
              zcnst = fo2
            ELSE IF ( speci(js) == 'CO2       ' ) THEN
              zcnst = fco2
            ELSE IF ( speci(js) == 'H2        ' ) THEN
              zcnst = fh2
            ELSE IF ( speci(js) == 'CH4       ' ) THEN
              zcnst = fch4
            ELSE
              cmessage=' Value not supplied for '//speci(js)//          &
                       ',setting value to zero'
! DEPENDS ON: ereport
              CALL EREPORT('ASAD_FYINIT',-1,cmessage)
              zcnst = 0.0
            ENDIF

            DO jl = 1, n_points
              y(jl,js) = zcnst * tnd(jl)
            ENDDO
          ENDDO
        ENDIF

!       2.   Initialise rest of species (potentially each call).
!            ---------- ---- -- ------- ------------ ---- ------

!       2.1  Species of type 'TR' and 'FT' are copied from tracers
!            on every call to this routine.

        DO j = 1, iftr
          js = ilftr(j)
          itr = madvtr(js)
          DO jl = 1, n_points
            y(jl,js) = f(jl,itr)
          ENDDO
        ENDDO

!       2.2  Constant but set by user, type 'CF' are only
!            initialised on first call per chemical step (this
!            routine may be called more than once per chemical
!            step).

!       IF ( istep /= jsubs ) THEN
!         istep = jsubs
!         kstep = jsubs
!         inl = n_points
!         DO j = 1, icf
!           js = ilcf(j)
!           CALL INICNT( speci(js), y(1,js), inl )
!         ENDDO
!       ENDIF
!
!      *****       corrected to call inicnt not only on first call
!                  per chemical step.

        kstep = jsubs
        inl = n_points
        DO j = 1, icf
          js = ilcf(j)
! DEPENDS ON: asad_inicnt
          CALL ASAD_INICNT( speci(js), y(1,js), inl )
        ENDDO

        RETURN
        END SUBROUTINE ASAD_FYINIT
#endif
