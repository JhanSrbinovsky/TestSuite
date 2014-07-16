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
!  Description:
!    Symbolic Backward-Euler solver for ASAD chemical integration
!    package.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Current Code Owner:       Olaf Morgenstern/Glenn Carver
!                            Colin Johnson
!
!
!     ASAD: bedriv                Version: bedriv.f 1.1 30/08/07
!
!     Purpose.
!     --------
!     To organise the integration of the chemical rate equations using
!     the Backward Euler implicit integrator.
!
!     Interface
!     ---------
!     Called from chemistry driver routine cdrive.
!
!     Method.
!     -------
!     Solves the non-linear system of equations via a Backward Euler
!     technique. NO3 and N2O5 are solved in a coupled way. See
!
!     Hertel et al., Test of two numerical schemes for use in atmospheric
!     transport-chemistry models, Atm. Env., 27A, 16, 2591-2611, 1993.
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
      SUBROUTINE asad_bedriv(nslon,nslat, n_points)

      USE ASAD_MOD
      IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "parvars.h"
#include "typsize.h"

! Subroutine interface
      INTEGER, INTENT(IN) :: n_points
      INTEGER, INTENT(IN) :: nslon
      INTEGER, INTENT(IN) :: nslat

! Local variables
      INTEGER, PARAMETER :: prodmax = 50 ! maximum number of reactions that
                                         ! produce or destroy any species
! Number of iterations in BE solver
      INTEGER, PARAMETER :: besteps = 8

! Number of SS/TR species
      INTEGER, SAVE :: knspec

      INTEGER :: jit
      INTEGER :: j
      INTEGER :: jp
      INTEGER :: jr
      INTEGER :: js
      INTEGER :: k
      INTEGER :: k1
      INTEGER :: i
      INTEGER :: if
      INTEGER :: jr1
      INTEGER :: jr2
      INTEGER :: jr3
      INTEGER :: jr4
      INTEGER :: jr5

      INTEGER, SAVE :: in2o5=0   ! position of N2O5 in species array
      INTEGER, SAVE :: ino3=0    ! position of NO3 in species array
      INTEGER, SAVE :: ino2=0    ! position of NO2 in species array
      INTEGER, SAVE :: pn2o5=0   ! position of N2O5 + h nu in reaction array
      INTEGER, SAVE :: tn2o5=0   ! position of N2O5 + M in reaction array
      INTEGER, SAVE :: tno2no3=0 ! position of NO2 + NO3 + M in reaction array

      LOGICAL, SAVE :: first = .TRUE.

! Position of SS/TR species
      INTEGER :: kspec(jpspec)

! Total production, loss and intermediate terms in budget calculation
      REAL, DIMENSION(theta_field_size) :: pdn
      REAL, DIMENSION(theta_field_size) :: l
      REAL, DIMENSION(theta_field_size) :: r1
      REAL, DIMENSION(theta_field_size) :: r2
      REAL, DIMENSION(theta_field_size) :: l1
      REAL, DIMENSION(theta_field_size) :: l2
      REAL, DIMENSION(theta_field_size) :: l3

! List of reactions, for all species, that create or destroy it
      INTEGER :: uprodreac(jpspec, prodmax)  ! unimol. prod reactions
      INTEGER :: ulossreac(jpspec, prodmax)  ! unimol. loss reactions
      INTEGER :: bprodreac(jpspec, prodmax)  ! bi,-termol. prod reactions
      INTEGER :: blossreac(jpspec, prodmax)  ! bi. loss reactions
      INTEGER :: blosspartner(jpspec, prodmax) ! partner in loss reaction

! Number of reactions that produce or cost a species
      INTEGER :: nuintprod(jpspec)
      INTEGER :: nbintprod(jpspec)
      INTEGER :: nfracprod(jpspec)
      INTEGER :: nuloss(jpspec)
      INTEGER :: nbloss(jpspec)

! Fractional production coefficient
      REAL :: fraction(jpspec, prodmax)

! Number and pointer of tracer variables in species array
      INTEGER :: iftr, ilftra(jpctr)

! k mod 5, precomputed
      INTEGER, SAVE :: mod5(0:prodmax)

      IF (first) THEN
        first = .FALSE.
! find N2O5, NO3, and NO2 species in species array
        DO js = 1,jpspec
          SELECT CASE (speci(js))
            CASE ('N2O5      ')
              in2o5 = js
            CASE ('NO3       ')
              ino3  = js
            CASE ('NO2       ')
              ino2  = js
          END SELECT
        END DO

! find N2O5 + M and NO2 + NO3 + M in termolecular reactions array
        DO js = 1,jptk
          SELECT CASE (spt(js,1))
            CASE ('N2O5      ')
              tn2o5   = ntrkx(js)
            CASE ('NO2       ')
              IF (spt(js,2) == 'NO3       ') tno2no3 = ntrkx(js)
            CASE ('NO3       ')
              IF (spt(js,2) == 'NO2       ') tno2no3 = ntrkx(js)
          END SELECT
        END DO

! find N2O5 photolysis in photolysis array
        DO js = 1,jppj
          IF (spj(js,1) == 'N2O5      ') pn2o5 = nprkx(js)
        END DO

! Find SS species
        knspec = 0
        DO js=1,jpspec
          IF (ctype(js) == jpna) THEN
            knspec = knspec + 1
            kspec(knspec) = js
          END IF
        END DO
! Find TR species
        DO js=1,jpspec
          IF (ctype(js) == jpsp) THEN
            knspec = knspec + 1
            kspec(knspec) = js
          END IF
        END DO

        IF (mype == 0) THEN
          WRITE(6,*) 'N2O5=',in2o5,' NO3 ',ino3,' NO2 ',ino2
          WRITE(6,*) 'tn2o5=',tn2o5,'tno2no3=',tno2no3,               &
            ' pn2o5 ',pn2o5
        END IF

! Fill production and destruction arrays
        uprodreac = 0
        ulossreac = 0
        bprodreac = 0
        blossreac = 0
        fraction = 1.
        nuintprod = 0
        nbintprod = 0
        nfracprod= 0
        nuloss = 0
        nbloss = 0
        DO jr=1,jpnr
! Do production first, for reactions with non-fractional products
          IF (nfrpx(jr) == 0) THEN
            DO jp=3,jpmsp
              js = nspi(jr,jp)
              IF (js > 0) THEN
                IF (jr <= nuni) THEN
                  nuintprod(js) = nuintprod(js) + 1
                  IF (nuintprod(js) > prodmax)                          &
! DEPENDS ON: ereport
                    CALL ereport('ASAD_BEDRIV',1,                       &
                         'Increase prodmax.')
                  uprodreac(js,nuintprod(js)) = jr
                ELSE
                  nbintprod(js) = nbintprod(js) + 1
                  IF (nbintprod(js) > prodmax)                          &
! DEPENDS ON: ereport
                    CALL ereport('ASAD_BEDRIV',1,                       &
                         'Increase prodmax.')
                  bprodreac(js,nbintprod(js)) = jr
                END IF
              END IF
            END DO
          END IF
! Do loss
          DO jp=1,2
            js = nspi(jr,jp)
            IF (js > 0) THEN
              IF (jr <= nuni) THEN
                nuloss(js) = nuloss(js) + 1
                IF (nuloss(js) > prodmax)                               &
! DEPENDS ON: ereport
                  CALL ereport('ASAD_BEDRIV',1,                         &
                     'Increase prodmax.')
                ulossreac(js,nuloss(js)) = jr
              ELSE
                nbloss(js) = nbloss(js) + 1
                IF (nbloss(js) > prodmax)                               &
! DEPENDS ON: ereport
                  CALL ereport('ASAD_BEDRIV',1,                         &
                     'Increase prodmax.')
                blossreac(js,nbloss(js)) = jr
! Remember which partner the loss reaction is done with
                IF (jp == 1) THEN
                  blosspartner(js,nbloss(js)) = nspi(jr,2)
                ELSE
                  blosspartner(js,nbloss(js)) = nspi(jr,1)
                END IF
              END IF
            END IF
          END DO
        END DO

! Do production by reactions with fractional products
        nfracprod = nbintprod
        DO i=1,nnfrp
          js = ntabfp(i,1) ! species index
          jr = ntabfp(i,2) ! reaction index
          if = ntabfp(i,3) ! pointer to fraction index
          nfracprod(js) = nfracprod(js) + 1
          IF (nfracprod(js) > prodmax)                                  &
! DEPENDS ON: ereport
            CALL ereport('ASAD_BEDRIV',1,                               &
                 'Increase prodmax.')
          bprodreac(js,nfracprod(js)) = jr
          fraction(js,nfracprod(js)) = frpb(if)
        END DO
        IF (mype == 0) THEN
          WRITE(6,*) 'NUINTPROD = ',nuintprod
          WRITE(6,*) 'NBINTPROD = ',nbintprod
          WRITE(6,*) 'NFRACPROD = ',nfracprod
          WRITE(6,*) 'NULOSS = ',nuloss
          WRITE(6,*) 'NBLOSS = ',nbloss
          WRITE(6,*) 'UPRODREAC = ',uprodreac(1,:)
          WRITE(6,*) 'BPRODREAC = ',bprodreac(1,:)
          WRITE(6,*) 'ULOSSREAC = ',ulossreac(1,:)
          WRITE(6,*) 'BLOSSREAC = ',blossreac(1,:)
        END IF

! fill tracer pointer
        iftr = 0
        DO js=1,jpspec
          IF ((ctype(js) == jpif) .OR. (ctype(js) == jpsp)) THEN
            iftr = iftr + 1
            ilftra(iftr) = js
          END IF
        END DO

! Precompute k mod 5
        DO k=0,prodmax
          mod5(k) = mod(k,5)
        END DO
      END IF  ! End of initialization

! find parameters needed in chemistry, special reactions etc
! Assign sensible values to species array y
! DEPENDS ON: asad_ftoy
      CALL asad_ftoy(.false.,1, n_points)
      IF (nstst  /=  0) THEN 
! DEPENDS ON: asad_steady
        CALL asad_steady( n_points )
      ENDIF

! Save previous state of species array
      ydot = y
!
!  Start Loop - perform backward Euler iteration
      DO jit=1,besteps

! Loop over species excluding unchanged ones
        DO j=1,knspec
          js = kspec(j)  ! js is now species number

          IF (js  /=  in2o5) THEN

! calculate production. Unimolecular reactions
! group terms in groups of 5.
! Number of unimolecular production terms for species js
            k1 = nuintprod(js)
            SELECT CASE (mod5(k1))
              CASE (4)
                jr1 = uprodreac(js,k1  )
                jr2 = uprodreac(js,k1-1)
                jr3 = uprodreac(js,k1-2)
                jr4 = uprodreac(js,k1-3)
                pdn = rk(:,jr1) * y(:,nspi(jr1,1))                      &
                  + rk(:,jr2) * y(:,nspi(jr2,1))                        &
                  + rk(:,jr3) * y(:,nspi(jr3,1))                        &
                  + rk(:,jr4) * y(:,nspi(jr4,1))
              CASE (3)
                jr1 = uprodreac(js,k1  )
                jr2 = uprodreac(js,k1-1)
                jr3 = uprodreac(js,k1-3)
                pdn = rk(:,jr1) * y(:,nspi(jr1,1))                      &
                  + rk(:,jr2) * y(:,nspi(jr2,1))                        &
                  + rk(:,jr3) * y(:,nspi(jr3,1))
              CASE (2)
                jr1 = uprodreac(js,k1  )
                jr2 = uprodreac(js,k1-1)
                pdn = rk(:,jr1) * y(:,nspi(jr1,1))                      &
                  + rk(:,jr2) * y(:,nspi(jr2,1))
              CASE (1)
                jr1 = uprodreac(js,k1  )
                pdn = rk(:,jr1) * y(:,nspi(jr1,1))
              CASE (0)
                pdn = 0.
            END SELECT
            DO k = 1,k1-4,5
              jr1 = uprodreac(js,k  )
              jr2 = uprodreac(js,k+1)
              jr3 = uprodreac(js,k+2)
              jr4 = uprodreac(js,k+3)
              jr5 = uprodreac(js,k+4)
              pdn = pdn + rk(:,jr1) * y(:,nspi(jr1,1))                  &
                    + rk(:,jr2) * y(:,nspi(jr2,1))                      &
                    + rk(:,jr3) * y(:,nspi(jr3,1))                      &
                    + rk(:,jr4) * y(:,nspi(jr4,1))                      &
                    + rk(:,jr5) * y(:,nspi(jr5,1))
            END DO

! Production by bimolecular reactions, integer products
            k = nbintprod(js)
            DO WHILE (k >= 5)
              jr1 = bprodreac(js,k  )
              jr2 = bprodreac(js,k-1)
              jr3 = bprodreac(js,k-2)
              jr4 = bprodreac(js,k-3)
              jr5 = bprodreac(js,k-4)
              pdn = pdn + rk(:,jr1)*y(:,nspi(jr1,1))*y(:,nspi(jr1,2))   &
                    + rk(:,jr2)*y(:,nspi(jr2,1))*y(:,nspi(jr2,2))       &
                    + rk(:,jr3)*y(:,nspi(jr3,1))*y(:,nspi(jr3,2))       &
                    + rk(:,jr4)*y(:,nspi(jr4,1))*y(:,nspi(jr4,2))       &
                    + rk(:,jr5)*y(:,nspi(jr5,1))*y(:,nspi(jr5,2))
              k = k - 5
            END DO
            SELECT CASE (k)
              CASE (4)
                jr1 = bprodreac(js,1)
                jr2 = bprodreac(js,2)
                jr3 = bprodreac(js,3)
                jr4 = bprodreac(js,4)
                pdn = pdn + rk(:,jr1)*y(:,nspi(jr1,1))*y(:,nspi(jr1,2)) &
                      + rk(:,jr2)*y(:,nspi(jr2,1))*y(:,nspi(jr2,2))     &
                      + rk(:,jr3)*y(:,nspi(jr3,1))*y(:,nspi(jr3,2))     &
                      + rk(:,jr4)*y(:,nspi(jr4,1))*y(:,nspi(jr4,2))
              CASE (3)
                jr1 = bprodreac(js,1)
                jr2 = bprodreac(js,2)
                jr3 = bprodreac(js,3)
                pdn = pdn + rk(:,jr1)*y(:,nspi(jr1,1))*y(:,nspi(jr1,2)) &
                      + rk(:,jr2)*y(:,nspi(jr2,1))*y(:,nspi(jr2,2))     &
                      + rk(:,jr3)*y(:,nspi(jr3,1))*y(:,nspi(jr3,2))
              CASE (2)
                jr1 = bprodreac(js,1)
                jr2 = bprodreac(js,2)
                pdn = pdn + rk(:,jr1)*y(:,nspi(jr1,1))*y(:,nspi(jr1,2)) &
                      + rk(:,jr2)*y(:,nspi(jr2,1))*y(:,nspi(jr2,2))
              CASE (1)
                jr1 = bprodreac(js,1)
                pdn = pdn + rk(:,jr1)*y(:,nspi(jr1,1))*y(:,nspi(jr1,2))
            END SELECT

! Fractional products here. Note that only bimolecular reactions
! can have fractional products at the moment.
            DO k=nbintprod(js) + 1, nfracprod(js)
              jr = bprodreac(js,k)
              pdn = pdn + rk(:,jr) * y(:,nspi(jr,1)) *                  &
                      y(:,nspi(jr,2)) * fraction(js,k)
            END DO

! Calculate loss. Unroll loss loops in groups of 5, to speed up
! calculation.
! Unimolecular loss first
            k1 = nuloss(js)
            SELECT CASE (mod5(k1))
              CASE (4)
                l = rk(:,ulossreac(js,k1  ))+rk(:,ulossreac(js,k1-1))   &
                  + rk(:,ulossreac(js,k1-2))+rk(:,ulossreac(js,k1-3))
              CASE (3)
                l = rk(:,ulossreac(js,k1  ))+rk(:,ulossreac(js,k1-1))   &
                  + rk(:,ulossreac(js,k1-2))
              CASE (2)
                l = rk(:,ulossreac(js,k1  ))+rk(:,ulossreac(js,k1-1))
              CASE (1)
                l = rk(:,ulossreac(js,k1))
              CASE (0)
                l = 0.
            END SELECT
            DO k=1,k1-4,5
              l = l + rk(:,ulossreac(js,k  ))+rk(:,ulossreac(js,k+1))   &
                    + rk(:,ulossreac(js,k+2))+rk(:,ulossreac(js,k+3))   &
                    + rk(:,ulossreac(js,k+4))
            END DO

! bimolecular loss
            k = nbloss(js)
            DO WHILE (k >= 5)
              l = l+rk(:,blossreac(js,k  ))*y(:,blosspartner(js,k  ))   &
                   +rk(:,blossreac(js,k-1))*y(:,blosspartner(js,k-1))   &
                   +rk(:,blossreac(js,k-2))*y(:,blosspartner(js,k-2))   &
                   +rk(:,blossreac(js,k-3))*y(:,blosspartner(js,k-3))   &
                   +rk(:,blossreac(js,k-4))*y(:,blosspartner(js,k-4))
              k = k - 5
            END DO
            SELECT CASE (k)
              CASE (4)
               l = l+rk(:,blossreac(js,4))*y(:,blosspartner(js,4))      &
                    +rk(:,blossreac(js,3))*y(:,blosspartner(js,3))      &
                    +rk(:,blossreac(js,2))*y(:,blosspartner(js,2))      &
                    +rk(:,blossreac(js,1))*y(:,blosspartner(js,1))
              CASE (3)
               l = l+rk(:,blossreac(js,3))*y(:,blosspartner(js,3))      &
                    +rk(:,blossreac(js,2))*y(:,blosspartner(js,2))      &
                    +rk(:,blossreac(js,1))*y(:,blosspartner(js,1))
              CASE (2)
               l = l+rk(:,blossreac(js,2))*y(:,blosspartner(js,2))      &
                    +rk(:,blossreac(js,1))*y(:,blosspartner(js,1))
              CASE (1)
               l = l+rk(:,blossreac(js,1))*y(:,blosspartner(js,1))
            END SELECT

! add dry and wet deposition terms
            IF (ldepd(js)) l = l + dpd(:,js)
            IF (ldepw(js)) l = l + dpw(:,js)

! add emission terms
            IF (lemit(js)) pdn = pdn + emr(:,js)

! calculate new concentration. Note that family chemistry is not supported
! by the BE solver

! do steady-state species first
            IF (ctype(js) == jpna) THEN

              y(:,js) = pdn/l

! Do tracer here. Regular case first (i.e., tracer # NO3 or N2O5)

            ELSEIF (js  /=  ino3) THEN

              y(:,js) = (ydot(:,js) + cdt*pdn)/(1.0 + cdt * l)

            ELSE

! do NO3 and N2O5 tracers combined

! R1 =  N2O5 + h nu and N2O5 + M
              r1 = rk(:,pn2o5) + rk(:,tn2o5)

! Pdn = all production terms for NO3 except N2O5 + h nu and N2O5 + M
              pdn = pdn - r1 * y(:,in2o5)
! L = all loss terms for NO3 (unchanged)


! R2 = NO2 + NO3 + M
              r2 = rk(:,tno2no3) * y(:,ino2)

! L1 = all loss terms of N2O5
              l1 = 0.
! Calculate loss
              DO k=1,nuloss(in2o5)
                l1 = l1 + rk(:,ulossreac(in2o5,k))
              END DO
              DO k=1,nbloss(in2o5)
                l1 = l1 + rk(:,blossreac(in2o5,k)) *                    &
                          y(:,blosspartner(in2o5,k))
              END DO

! add dry and wet deposition terms
              IF (ldepd(in2o5)) l1 = l1 + dpd(:,in2o5)
              IF (ldepw(in2o5)) l1 = l1 + dpw(:,in2o5)

! L2 = 1 + cdt * L
              l2 = 1.0 + cdt * l

! L3 = 1 + cdt * L1
              l3 = 1.0 + cdt * l1

! New value for NO3
              y(:,ino3) = (l3 * (ydot(:,js) + p*cdt) +                  &
                           r1 * cdt * ydot(:,in2o5)) /                  &
                          (l3 * l2 - r1 * r2 * cdt * cdt)

! New value for N2O5
              y(:,in2o5) = (ydot(:,in2o5) + r2*cdt*y(:,js))/l3

            END IF    ! distinction between SS and TR species
          END IF      ! exclude N2O5 tracer
        END DO        ! species loop
      END DO          ! BE iteration loop

! update tracer information
      DO j = 1, jpctr
        f(:,j) = y(:,ilftra(j))
      END DO

! finish BE step
      RETURN
      END SUBROUTINE asad_bedriv
#endif
