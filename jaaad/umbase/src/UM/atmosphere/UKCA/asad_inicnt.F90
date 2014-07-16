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
! Purpose: In this case, ASAD will treat the species as a constant
!     but will call this routine so that the user may set the
!     values differently at each gridpoint for example.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_FYINIT
!
! Current code owner: Glenn Carver/Colin Johnson
!                     Oliver Wild
!
!     Interface
!     On entry, the following will be set:
!              species - character name of species to set.
!                        Will be the same as listed in chch.d file
!              klen    - length of array, y.
!
!     On exit, the following must be set:
!              y       - Array of points to set for the species.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_INICNT( species, y, klen )

        USE ASAD_MOD,    ONLY: wp, tnd
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: klen      ! No of spatial points

        CHARACTER (LEN=10), INTENT(IN)  :: species  ! Species char strng

        REAL, INTENT(OUT)   :: y(klen)   ! Species concentration

!       Local variables

        INTEGER :: jl

        CHARACTER (LEN=72) :: cmessage

!       1.  Copy water into ASAD array.

        IF ( species(1:3) /= 'H2O' ) THEN
           cmessage= 'Expected species H2O but got '//species
! DEPENDS ON: ereport
           CALL EREPORT('ASAD_INICNT',124,cmessage)
        ENDIF

        DO jl = 1, klen
          y(jl) = wp(jl)*tnd(jl)
        ENDDO

        RETURN
        END SUBROUTINE ASAD_INICNT
#endif
