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
!    Computes steady-state species for Newton-Raphson integrator.
!    Part of the ASAD chemical solver.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Current Code Owner:       Olaf Morgenstern/Oliver Wild
!                            Colin Johnson
!
!     ASAD: ycn                      Version: steady.f 4.1 20/7/07
!
!     Purpose
!     -------
!
!  Routine to explicitly define steady state expressions - prevents
!  generality, but removes need for sluggish iteration round family
!  stuff.
!
!  Important Notes:
!     1) Ordering of calculations is important - need to avoid feedbacks!
!     2) Needs to be rewritten whenever reaction numbering is changed.
!
!  Additions:
!     To improve generality, reactions involving steady state species
!  are selected in 'setsteady' and loaded into nss* integer arrays,
!  which are then used in this routine. This causes a very slight
!  increase in CPU time, but removes the need to rewrite the routine
!  whenever a new species is added or a reaction changed.
!
!                                            Oliver   (3 Feb 1998)
! We add more general terms for the steady-state species. It is assumed
! that O(1D), O(3P), H and N can be put in steady state.
!
!
!     Method
!     ------
!
!   We sum up production and loss terms for SS species, and divide.
!   Moreover, corresponding terms in the Jacobian are calculated that
!   account for the dependence of steady-state variables on tracer
!   variables.
!
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
      SUBROUTINE asad_steady( kl )

      USE ASAD_MOD
      IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "parvars.h"
#include "typsize.h"

! Subroutine interface
      INTEGER, INTENT(IN) :: kl

! Local variables
      INTEGER, PARAMETER :: n_o3=1
      INTEGER, PARAMETER :: n_oh=2
      INTEGER, PARAMETER :: n_ho2=3
      INTEGER, PARAMETER :: n_no=4

      INTEGER :: jl
      INTEGER :: jl1
      INTEGER :: jl2
      INTEGER :: jr
      INTEGER :: ix
      INTEGER :: i
      INTEGER :: j

      REAL :: ssnum(theta_field_size)
      REAL :: ssden(theta_field_size)

! Add here derivatives w.r.t. O3, OH, and HO2, and NO of numerator and
! denominator of steady state species
      REAL :: dssnum(theta_field_size,4)
      REAL :: dssden(theta_field_size,4)
!
! Set up loops correctly
      jl1 = 1
      jl2 = kl
      IF (kl == 1) THEN
        jl1 = jl1+jlst
        jl2 = jl1
      END IF
!
! Loop through steady state species

      DO ix = 1,nsst
        ssnum = 0.
        ssden = 0.
        dssden = 0.
        dssnum = 0.

! Production terms
        DO jr = 1,nsspt(ix)
          i = nsspi(ix,jr)
          IF (i <= nuni) THEN
            ssnum(jl1:jl2) = ssnum(jl1:jl2) +                           &
                rk(jl1:jl2,i)*y(jl1:jl2,nspi(i,1))

            IF ((ix < 5) .AND. (nspi(i,1) == nspo3 ))                   &
! add terms to derivative for d(j[O3])/d[O3] = j_o3
              dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +             &
                                       rk(jl1:jl2,i)
            IF ((ix < 5) .AND. (nspi(i,1) == nspno ))                   &
! add terms to derivative for d(j[NO])/d[NO] = j_no
              dssnum(jl1:jl2,n_no) = dssnum(jl1:jl2,n_no) +             &
                                       rk(jl1:jl2,i)
          ELSE
            ssnum(jl1:jl2) = ssnum(jl1:jl2) +                           &
                rk(jl1:jl2,i)*y(jl1:jl2,nspi(i,1))*y(jl1:jl2,nspi(i,2))
            IF (ix < 5) THEN

! add terms for derivative w.r.t. ozone.
              IF (nspi(i,1) == nspo1d)                                  &
                dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +           &
                     rk(jl1:jl2,i)                                      &
                    *y(jl1:jl2,nspi(i,2))*deriv(jl1:jl2,1,n_o3)

              IF (nspi(i,2) == nspo1d)                                  &
                  dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +         &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))*deriv(jl1:jl2,1,n_o3)

              IF (nspi(i,1) == nspo3p)                                  &
                  dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +         &
                     rk(jl1:jl2,i)                                      &
                    *y(jl1:jl2,nspi(i,2))*deriv(jl1:jl2,2,n_o3)

              IF (nspi(i,2) == nspo3p)                                  &
                  dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +         &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))*deriv(jl1:jl2,2,n_o3)

              IF (nspi(i,1) == nspo3)                                   &
                dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,2))

              IF (nspi(i,2) == nspo3)                                   &
                dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))

! add terms for derivative w.r.t OH
              IF (nspi(i,1) == nspo3p)                                  &
! add terms to derivates for d(a[A][B])
                dssnum(jl1:jl2,n_oh) = dssnum(jl1:jl2,n_oh) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,2))*deriv(jl1:jl2,2,n_oh)

              IF (nspi(i,2) == nspo3p)                                  &
! add terms to derivates for d(a[O2][O1D])/d[O3] and b[N2][O1D]
                dssnum(jl1:jl2,n_oh) = dssnum(jl1:jl2,n_oh) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))*deriv(jl1:jl2,2,n_oh)

              IF (nspi(i,1) == nspoh)                                   &
! add terms to derivates for d(a[O2][O1D])/d[O3] and b[N2][O1D]
                dssnum(jl1:jl2,n_oh) = dssnum(jl1:jl2,n_oh) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,2))

              IF (nspi(i,2) == nspoh)                                   &
! add terms to derivates for d(a[O2][O1D])/d[O3] and b[N2][O1D]
                dssnum(jl1:jl2,n_oh) = dssnum(jl1:jl2,n_oh) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))

! add terms for derivative w.r.t HO2
              IF (nspi(i,1) == nspo3p)                                  &
                dssnum(jl1:jl2,n_ho2) = dssnum(jl1:jl2,n_ho2) +         &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,2))*deriv(jl1:jl2,2,n_ho2)

              IF (nspi(i,2) == nspo3p)                                  &
                dssnum(jl1:jl2,n_ho2) = dssnum(jl1:jl2,n_ho2) +         &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))*deriv(jl1:jl2,2,n_ho2)

              IF (nspi(i,1) == nspho2)                                  &
                dssnum(jl1:jl2,n_ho2) = dssnum(jl1:jl2,n_ho2) +         &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,2))

              IF (nspi(i,2) == nspho2)                                  &
                dssnum(jl1:jl2,n_ho2) = dssnum(jl1:jl2,n_ho2) +         &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))

! add terms for derivative w.r.t NO
              IF (nspi(i,1) == nspno)                                   &
                dssnum(jl1:jl2,n_no) = dssnum(jl1:jl2,n_no) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,2))

              IF (nspi(i,2) == nspno)                                   &
                dssnum(jl1:jl2,n_no) = dssnum(jl1:jl2,n_no) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))

            END IF
          END IF
        END DO
!
! Destruction terms
        DO jr = 1,nssrt(ix)
          i = nssri(ix,jr)
          j = nssrx(ix,jr)
          IF (i <= nuni) THEN
            ssden(jl1:jl2) = ssden(jl1:jl2) + rk(jl1:jl2,i)
          ELSE
            ssden(jl1:jl2) = ssden(jl1:jl2) +                           &
              rk(jl1:jl2,i) * y(jl1:jl2,nspi(i,j))
            IF (ix < 5) THEN
              IF (nspi(i,j) == nspo3 )                                  &
                dssden(jl1:jl2,n_o3 ) = dssden(jl1:jl2,n_o3 ) +         &
                  rk(jl1:jl2,i)
              IF (nspi(i,j) == nspoh )                                  &
                dssden(jl1:jl2,n_oh ) = dssden(jl1:jl2,n_oh ) +         &
                  rk(jl1:jl2,i)
              IF (nspi(i,j) == nspho2)                                  &
                dssden(jl1:jl2,n_ho2) = dssden(jl1:jl2,n_ho2) +         &
                  rk(jl1:jl2,i)
              IF (nspi(i,j) == nspno )                                  &
                dssden(jl1:jl2,n_no ) = dssden(jl1:jl2,n_no ) +         &
                  rk(jl1:jl2,i)
            END IF
          END IF
        END DO
!
! Steady state and derivatives of steady state
        y(jl1:jl2,nssi(ix)) = ssnum(jl1:jl2)/ssden(jl1:jl2)
        IF (ix < 5) THEN
          DO jr =1,4
            deriv(jl1:jl2,ix,jr) =                                      &
                (ssden(jl1:jl2)*dssnum(jl1:jl2,jr) -                    &
                 ssnum(jl1:jl2)*dssden(jl1:jl2,jr))/                    &
                 (ssden(jl1:jl2) * ssden(jl1:jl2))
          END DO
        END IF
      END DO

! rescale deriv to mean [O3]/[O] * d[O]/d[O3], where [O] = [O(1D)] or [O(3P)]
! rescale derivative of O(1D)
      DO jl=jl1,jl2
        IF (y(jl,nspo1d) > peps) THEN
          deriv(jl,1,n_o3) = deriv(jl,1,n_o3 )*y(jl,nspo3 )/y(jl,nspo1d)
          deriv(jl,1,n_oh) = deriv(jl,1,n_oh )*y(jl,nspoh )/y(jl,nspo1d)
          deriv(jl,1,n_ho2)= deriv(jl,1,n_ho2)*y(jl,nspho2)/y(jl,nspo1d)
          deriv(jl,1,n_no )= deriv(jl,1,n_no )*y(jl,nspno )/y(jl,nspo1d)
        ELSE
          deriv(jl,1,:) = 1.
        ENDIF
        IF (nspo3p > 0) THEN
          IF (y(jl,nspo3p) > peps) THEN
            deriv(jl,2,n_o3 )= deriv(jl,2,n_o3 )*y(jl,nspo3 )/          &
                               y(jl,nspo3p)
            deriv(jl,2,n_oh )= deriv(jl,2,n_oh )*y(jl,nspoh )/          &
                               y(jl,nspo3p)
            deriv(jl,2,n_ho2)= deriv(jl,2,n_ho2)*y(jl,nspho2)/          &
                               y(jl,nspo3p)
            deriv(jl,2,n_no )= deriv(jl,2,n_no )*y(jl,nspno )/          &
                               y(jl,nspo3p)
          ELSE
            deriv(jl,2,:) = 1.
          ENDIF
        END IF
        IF (nsph > 0) THEN
        IF (y(jl,nsph  ) > peps) THEN
          deriv(jl,4,n_o3) = deriv(jl,4,n_o3 )*y(jl,nspo3 )/y(jl,nsph  )
          deriv(jl,4,n_oh) = deriv(jl,4,n_oh )*y(jl,nspoh )/y(jl,nsph  )
          deriv(jl,4,n_ho2)= deriv(jl,4,n_ho2)*y(jl,nspho2)/y(jl,nsph  )
          deriv(jl,4,n_no )= deriv(jl,4,n_no )*y(jl,nspno )/y(jl,nsph  )
        ELSE
          deriv(jl,4,:) = 1.
        ENDIF
        END IF
       END DO
!
      RETURN
      END SUBROUTINE asad_steady
#endif
