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
! Purpose: Calculates the total number density in gridboxes
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
!     The total number density at a point is given by: p=nkt
!     p-pressure , n-number density, k-boltzmann's constant,
!     t-temperature.
!
!     local variables
!     ---------------
!     zboltz     boltzmann's constant
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_TOTNUD(n_points)

        USE ASAD_MOD,             ONLY: tnd, p, t, pmintnd, pmin
        USE UKCA_CONSTANTS,       ONLY: zboltz
        IMPLICIT NONE


        INTEGER, INTENT(IN) :: n_points

!       Local variables

        INTEGER :: jl

        REAL :: zb

        zb = zboltz*1.0e6

!       1. Total number density (1e6 converts numbers to /cm**3).
!          ----- ------ ------- ---- -------- ------- -- --------

        DO jl = 1, n_points
          tnd(jl)     = p(jl) / ( zb * t(jl) )
          pmintnd(jl) = pmin * tnd(jl)
        ENDDO

        RETURN
        END SUBROUTINE ASAD_TOTNUD
#endif
