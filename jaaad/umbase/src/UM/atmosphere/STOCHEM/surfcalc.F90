#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE SURFCALC(si,sa,p0,t0,va,seaicefr,snowfr,psurf,         &
     & tsurf,vaero)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Calculate snow fraction on Eulerian grid from
!-                         snow amount data.  Interpolate other fields
!-                         onto the Eulerian grid from the Met. data.
!-
!-   Inputs  : SI,SA,P0,T0,VA
!-   Outputs : SEAICEFR,SNOWFR,PSURF,TSURF,VAERO
!-   Controls:
!-
!
! Current Code Owner: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  4.5    27/02/98  Created.  C.E. Johnson
!  4.5    24/09/98  Does P0, T0 and VA too.  W.J. Collins
!  6.1    22/10/04  No change
!
!-
!VVV  V2.2  SURFCALC 20/X/99
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!-----------------------------------------------------------------------
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: si
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: sa
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: t0
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: p0
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: va

      REAL, DIMENSION(nlnpe,nlpe),     INTENT(OUT) :: seaicefr
      REAL, DIMENSION(nlnpe,nlpe),     INTENT(OUT) :: snowfr
      REAL, DIMENSION(nlnpe,nlpe),     INTENT(OUT) :: psurf
      REAL, DIMENSION(nlnpe,nlpe),     INTENT(OUT) :: tsurf
      REAL, DIMENSION(nlnpe,nlpe),     INTENT(OUT) :: vaero

! The snow and ice fractions are only available once per day
! DEPENDS ON: met2data
      CALL MET2DATA(seaicefr,si,1,1)
! DEPENDS ON: met2data
      CALL MET2DATA(snowfr,sa,1,1)
! DEPENDS ON: met2data
      CALL MET2DATA(psurf,p0,1,1)
! DEPENDS ON: met2data
      CALL MET2DATA(tsurf,t0,1,1)
! DEPENDS ON: met2data
      CALL MET2DATA(vaero,va,1,1)

! Convert snow depth (kg/m2 or mm water equiv.) to snow cover
! fraction (Essery).
! Linear function of snow amount up to 10 kg/m2.
      WHERE(snowfr /= 0.0) snowfr=min(snowfr*0.1,1.0)

      END SUBROUTINE SURFCALC
#endif
