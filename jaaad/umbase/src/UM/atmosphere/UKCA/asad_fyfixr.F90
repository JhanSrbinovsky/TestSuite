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
! Purpose: Calculates family concentrations from fixed ratios, calculated
!          during a previous call to asad_ftoy.
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
!     Interface
!     ---------
!     Called from ftoy if the number of iterations requested equals
!     zero.
!
!     Method
!     ------
!     The family is partitioned amongst the members
!     using ratios calculated on a previous call to ftoy.
!     The concentrations are found from:
!                  Ym = Rmf*Z, Y1 = R1m*Ym, Y2 = R2m*Ym, .....
!
!     Local variables
!     ---------------
!     ifam           Index of family to which species belongs.
!     imaj           Index of major member of family to which
!                    species belongs.
!     itr            Index of model tracer to which species
!                    corresponds.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_FYFIXR(n_points)

        USE ASAD_MOD
        IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"

        INTEGER, INTENT(IN) :: n_points

!       Local variables

        INTEGER :: istart
        INTEGER :: iend
        INTEGER :: j       ! Loop variable
        INTEGER :: jl      ! Loop variable
        INTEGER :: js      ! Index
        INTEGER :: ifam    ! Index
        INTEGER :: imaj    ! Index
        INTEGER :: itr     ! Index

!       1.  Calculate major species of family
!           --------- ----- ------- -- ------

        istart = nlmajmin(1)
        iend   = nlmajmin(2)
        DO j = istart, iend
          js   = nlmajmin(j)
          ifam = moffam(js)
          DO jl = 1, n_points
            y(jl,js) = f(jl,ifam) * ratio(jl,js)
          ENDDO
        ENDDO

!       3.  Calculate minor species of family
!           --------- ----- ------- -- ------

        istart = nlmajmin(3)
        iend   = nlmajmin(4)
        DO j = istart, iend
          js   = nlmajmin(j)
          ifam = moffam(js)
          imaj = majors(ifam)
          IF ( ctype(js) /= jpif ) THEN
            DO jl = 1, n_points
              y(jl,js) = y(jl,imaj) * ratio(jl,js)
            ENDDO
          ELSE
            itr = madvtr(js)
            DO jl = 1, n_points
              IF ( linfam(jl,itr) ) y(jl,js) =                         &
                                    y(jl,imaj) * ratio(jl,js)
            ENDDO
          ENDIF
        ENDDO

        RETURN
        END SUBROUTINE ASAD_FYFIXR
#endif
