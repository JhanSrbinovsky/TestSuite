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
! Purpose: Calculates trimolecular rate coefficients
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CDRIVE
!
! Current code owner: Glenn Carver/Colin Johnson
!                     Oliver Wild
!
!     Method
!     ------
!     See the IUPAC reference material on their website for details on
!     calculation of termolecular rates.
!     http://www.iupac-kinetic.ch.cam.ac.uk/
!
!     Local variables
!     ---------------
!     zo         Low pressure limit to rate*density
!     zi         High pressure limit to rate
!     zr         Ratio of zo/zi
!     iho2       Reaction index for HO2+HO2+M
!     ih2o       Array index for advected tracer H2O
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_TRIMOL(n_points)

        USE ASAD_MOD
        IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"

        INTEGER, INTENT(IN) :: n_points

!       Local variables

        INTEGER, SAVE :: ih2o
        INTEGER       :: iho2
        INTEGER       :: j           ! Loop variable
        INTEGER       :: jl          ! Loop variable
        INTEGER       :: jtr         ! Loop variable
        INTEGER       :: jr          ! Index

        REAL :: zo
        REAL :: zi
        REAL :: zfc
        REAL :: zr

        LOGICAL, SAVE :: first = .true.

!       1.  Calculate trimolecular rate coefficients
!           --------- ------------ ---- ------------

        iho2  = 0

!       Check if H2O is an advected tracer
        IF (first) THEN
          ih2o = 0
          DO jtr = 1, jpctr
            if ( advt(jtr)  ==  'H2O       ' ) ih2o = jtr
          ENDDO
          first = .false.
        ENDIF

        DO j = 1, jptk
          jr = ntrkx(j)

          IF ( spt(j,1) == 'HO2    '.and.spt(j,2) == 'HO2    ' )       &
               iho2 = jr

          DO jl = 1, n_points
            zo = at(j,2) * t300(jl)**at(j,3) *                         &
                           exp( -at(j,4)/t(jl) ) * tnd(jl)
            zi = at(j,5) * t300(jl)**at(j,6) * exp( -at(j,7)/t(jl) )
            IF ( zo < peps ) THEN
              rk(jl,jr) = zi
            ELSE IF ( zi < peps ) THEN
              rk(jl,jr) = zo
            ELSE
              IF( at(j,1) <= 1.0 ) THEN
                zfc = at(j,1)
              ELSE
                zfc = exp( -t(jl)/at(j,1) )
              END IF
              zr = zo / zi
              rk(jl,jr) = (zo/(1.0+zr)) *                              &
                           zfc**( 1.0 / (1.0+((alog10(zr))**2)) )
            END IF
          END DO
        ENDDO              ! end of loop over jptk

!       2. Dependent reactions.
!          --------- ----------

!       HO2 + HO2 [+ M]
        IF (ih2o /= 0 .and. iho2 /= 0 ) THEN
!         h2o is an advected tracer
          DO jl = 1, n_points
            rk(jl,iho2) = rk(jl,iho2) *                                &
            ( 1.0 + 1.4E-21*f(jl,ih2o)*exp(2200./t(jl)) )
          ENDDO
        ELSE IF (ih2o == 0 .AND. iho2 /= 0) THEN
!         use modelled water concentration
          DO jl = 1, n_points
            rk(jl,iho2) = rk(jl,iho2) *                                &
            ( 1.0 + 1.4E-21*wp(jl)*tnd(jl)*exp(2200./t(jl)) )
          ENDDO
        ENDIF

        RETURN
        END SUBROUTINE ASAD_TRIMOL
#endif
