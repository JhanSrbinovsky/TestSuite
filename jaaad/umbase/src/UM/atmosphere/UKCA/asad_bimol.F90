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
! Purpose: To calculate bimolecular rate coefficients using 
!          data from the bimolecular ratefile ratb.d or
!          ukca_chem1 module
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CDRIVE 
!
!     Method
!     ------
!     Bimolecular rate coefficients are calculated from the Ahrenius
!     expression: k(T) = k * (T/300)^a * exp(-b/T) where T is the
!     temperature (in K). The parameters k, a and b are taken from
!     the bimolecular ratefile ratb.d.
!
!     The reactions CO + OH -> H + CO2 and OH + HONO2 are pressure
!     or density dependent and therefore need to be calculated
!     separately. However, this code should never need changing even
!     if neither of these reaction are included in the chemistry.
!
!     The reactions HO2+MeCO3, MeOO+MeOO and OH+C3H8 have temperature
!     dependent branching ratios. Therefore, they are calculated
!     separately.
!
!     Local Variables
!     ---------------
!     ih2o           Tracer index for H2O if it's an advective tracer
!     iohco          Reaction index for OH + CO.
!     iohhno3        Reaction index for OH + HONO2.
!     iohc3h8a       Reaction index for OH+C3H8 -> n-PrOO + H2O
!     iohc3h8b       Reaction index for OH+C3H8 -> i-PrOO + H2O
!     imeoomeooa     Reaction index for MeOO+MeOO -> HO2+HCHO+HO2+HCHO
!     imeoomeoob     Reaction index for MeOO+MeOO -> MeOH+HCHO
!     imeooho2a      Reaction index for MeOO+HO2 -> MeOOH
!     imeooho2b      Reaction index for MeOO+HO2 -> HCHO
!     z1,z2,z3,z4    Variables used to calculate pressure
!                    dependent rate coefficients.
!     ratioa2b       Branching ratio of branch A to branch B
!
!
! Current code owner: Glenn Carver/Colin Johnson
!                     Oliver Wild
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_BIMOL( n_points )

        USE ASAD_MOD
        IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"

        INTEGER, INTENT(IN) :: n_points

!       Local variables

        INTEGER :: jl
        INTEGER :: iohco
        INTEGER :: iohhno3
        INTEGER :: iho2
        INTEGER :: jtr
        INTEGER :: j
        INTEGER :: jr
        INTEGER :: iohc3h8a
        INTEGER :: iohc3h8b
        INTEGER :: imeoomeooa
        INTEGER :: imeoomeoob
        INTEGER :: imeooho2a
        INTEGER :: imeooho2b
        INTEGER :: asad_findreaction
        INTEGER, SAVE :: ih2o

        REAL :: ab1
        REAL :: ab2
        REAL :: ab3
        REAL :: z1
        REAL :: z2
        REAL :: z3
        REAL :: z4
        REAL :: ratioa2b
        REAL :: ratiob2total

        CHARACTER*10 :: r1
        CHARACTER*10 :: r2
        CHARACTER*10 :: prods(jpspb)

        LOGICAL, SAVE :: first = .true.

!       1. Calculate bimolecular rate coefficients
!          --------- ----------- ---- ------------

!       Compute intermediate results

        DO jl = 1, n_points
          t300(jl) = t(jl) / 300.0
        END DO

!       Check if H2O is an advected tracer

        IF (first) then
          ih2o = 0
          DO jtr = 1, jpctr
            if ( advt(jtr)  ==  'H2O       ' ) ih2o = jtr
          ENDDO
          first = .false.
        ENDIF

        iohco      = 0
        iohhno3    = 0
        iho2       = 0
        iohc3h8a   = 0
        iohc3h8b   = 0
        imeoomeooa = 0
        imeoomeoob = 0
        imeooho2a  = 0
        imeooho2b  = 0

!       Look for the reactions which need special treatment.

        r1 = '          '
        r2 = r1
        prods(:) = r1
        r1 = 'OH        '
        r2 = 'C3H8      '
        prods(1) = 'n-PrOO    '
        prods(2) = 'H2O       '
! DEPENDS ON: asad_findreaction
        iohc3h8a = asad_findreaction( r1, r2, prods, 2, spb, nbrkx,    &
                               jpbk+1, jpspb )
        prods(1) = 'i-PrOO    '
        prods(2) = 'H2O       '
! DEPENDS ON: asad_findreaction
        iohc3h8b = asad_findreaction( r1, r2, prods, 2, spb, nbrkx,    &
                               jpbk+1, jpspb )
        r1 = 'MeOO      '
        r2 = 'MeOO      '
        prods(1) = 'MeOH      '
        prods(2) = 'HCHO      '
! DEPENDS ON: asad_findreaction
        imeoomeooa = asad_findreaction( r1, r2, prods, 2, spb, nbrkx,  &
                                 jpbk+1, jpspb )

!! WARNING!!  must use 4 product version of ratb.d for this reaction!

        prods(1) = 'HO2       '
        prods(2) = 'HCHO      '
        prods(3) = 'HO2       '
        prods(4) = 'HCHO      '
! DEPENDS ON: asad_findreaction
        imeoomeoob = asad_findreaction( r1, r2, prods, 4, spb, nbrkx,  &
                                 jpbk+1, jpspb )
        r1 = 'HO2       '
        r2 = 'MeOO      '
        prods(1) = 'MeOOH     '
! DEPENDS ON: asad_findreaction
        imeooho2a = asad_findreaction( r1, r2, prods, 1, spb, nbrkx,   &
                                jpbk+1, jpspb )
        prods(1) = 'HCHO      '
! DEPENDS ON: asad_findreaction
        imeooho2b = asad_findreaction( r1, r2, prods, 1, spb, nbrkx,   &
                                jpbk+1, jpspb )

!       1.2  Compute rates

        DO j = 1, jpbk
          jr = nbrkx(j)

          IF ( ( spb(j,1) == 'OH     '.and.spb(j,2) == 'CO     ' ).or. &
             ( spb(j,1) == 'CO     '.and.spb(j,2) == 'OH     ' ) )     &
            iohco = jr
          IF ( ( spb(j,1) == 'OH     '.and.spb(j,2) == 'HONO2  ' ).or. &
             (  spb(j,1) == 'HONO2  '.and.spb(j,2) == 'OH     ' ) )    &
            iohhno3 = jr
          IF ( spb(j,1) == 'HO2    '.and.spb(j,2) == 'HO2    ' )       &
            iho2 = jr


          ab1 = ab(j,1)
          ab2 = ab(j,2)
          ab3 = ab(j,3)

          DO jl = 1, n_points
            IF ( ab2 == 0.0 .and. ab3 == 0.0 ) THEN
              rk(jl,jr) = ab1
            ELSE IF ( ab2 == 0.0 ) THEN
              rk(jl,jr) = ab1 * exp( -ab3 / t(jl) )
            ELSE IF ( ab3 == 0.0 ) THEN
              rk(jl,jr) = ab1 * t300(jl)**ab2
            ELSE
              rk(jl,jr) = ab1 * t300(jl)**ab2 * exp( -ab3 / t(jl) )
            ENDIF
          END DO
        END DO  ! end of loop over jpbk

!       2. Dependent reactions.
!          --------- ----------

        DO jl = 1, n_points

! OH + CO; updated with IUPAC March 2005 (Paul Young)
! GC: note the original March IUPAC summary had a mistake for this
! reaction. The values given in the datasheet are correct.
! k = k' . (1 + [N2]/4.2E19).. we use TND below instead of [N2] but 
! Paul suggests the reaction would probably go with [O2] anyway.

          IF ( iohco /= 0 )                                            &
            rk(jl,iohco)=rk(jl,iohco)*(1.0 + tnd(jl)/4.2e19)

! OH + HONO2; no change with IUPAC Nov 2003 (Paul Young)
          IF ( iohhno3 /= 0 ) THEN
            z1 = 2.4e-14 * exp(460.0/t(jl))
            z3 = 6.5e-34 * exp(1335.0/t(jl))
            z4 = 2.7e-17 * exp(2199.0/t(jl))
            z2 = z3*tnd(jl) / ( 1.0+z3*tnd(jl)/z4 )
            rk(jl,iohhno3) = z1 + z2
          ENDIF

! HO2 + HO2; no change with IUPAC Nov 2003 (Paul Young)
          IF ( ih2o /= 0 .and. iho2 /= 0) THEN
!           water is an advected tracer
            rk(jl,iho2) = rk(jl,iho2) *                                &
                     ( 1.0 + 1.4E-21*f(jl,ih2o)*exp(2200.0/t(jl)) )
          ELSE IF (ih2o == 0 .and. iho2 /= 0) then
!           use model water concentration
            rk(jl,iho2) = rk(jl,iho2) *                                &
            ( 1.0 + 1.4E-21*wp(jl)*tnd(jl)*exp(2200./t(jl)) )
          ENDIF
        END DO       ! end of jl loop

!       3. Temperature-Dependent branching ratios
!          ----------- --------- --------- ------
!       rk above was calculated using the total rate coefficients.
!       Here, rk is reduced according to the branching ratio.

        DO jl=1,n_points

!         OH + C3H8 -> n-PrOO ... Branch A
!         OH + C3H8 -> i-PrOO ... Branch B
          IF ( iohc3h8a /= 0 ) THEN
            ratioa2b        = 226.0 * t(jl)**(-0.64)*exp(-816.0/t(jl))
            rk(jl,iohc3h8a) = rk(jl,iohc3h8a)*(ratioa2b/(ratioa2b+1.0))
          ENDIF

          IF (iohc3h8b /= 0) THEN
           ratioa2b        = 226.0*t(jl)**(-0.64)*exp(-816.0/t(jl))
           rk(jl,iohc3h8b) = rk(jl,iohc3h8b)/(ratioa2b+1.0)
          ENDIF

!         MeOO + MeOO -> MeOH + HCHO       ... Branch A
!         MeOO + MeOO -> 2HO2 + 2HCHO      ... Branch B
          IF ( imeoomeooa /= 0 ) THEN
             ratiob2total      = 1.0/(1.0+(EXP(1300.0/t(jl)))/33.0)
             rk(jl,imeoomeooa) = rk(jl,imeoomeooa)*(1.0-ratiob2total)
          ENDIF

          IF ( imeoomeoob /= 0 ) THEN
            ratiob2total      = 1.0/(1.0+(EXP(1300.0/t(jl)))/33.0)
            rk(jl,imeoomeoob) = rk(jl,imeoomeoob)*(ratiob2total)
          ENDIF

!         Added from IUPAC March 2005 - Paul Young March, 2005
!         MeOO + HO2 -> MeOOH     ... Branch A
!         MeOO + HO2 -> HCHO      ... Branch B

          IF ( imeooho2a /= 0 ) THEN
            ratiob2total     = 1.0 / (1.0 + 498.0*EXP(-1160.0/t(jl)))
            rk(jl,imeooho2a) = rk(jl,imeooho2a)*(1.0-ratiob2total)
          ENDIF

          IF ( imeooho2b /= 0 ) THEN
           ratiob2total     = 1.0 / (1.0 + 498.0*EXP(-1160.0/t(jl)))
           rk(jl,imeooho2b) = rk(jl,imeooho2b)*(ratiob2total)
          ENDIF

        ENDDO

        RETURN
        END SUBROUTINE ASAD_BIMOL
#endif
