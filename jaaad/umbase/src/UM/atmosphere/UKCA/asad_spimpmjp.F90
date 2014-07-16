#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
!
! Copyright (c) 2008, Regents of the University of California
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer. 
!     * Redistributions in binary form must reproduce the above
!       copyright notice, this list of conditions and the following
!       disclaimer in the documentation and/or other materials provided
!       with the distribution. 
!     * Neither the name of the University of California, Irvine nor the
!       names of its contributors may be used to endorse or promote
!       products derived from this software without specific prior
!       written permission.
!
!       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
!       IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!       TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!       PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!       OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!       EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!       PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!       PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!       NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!       SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Driver for fully implicit ODE integrator.
!    Part of the ASAD chemical solver
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Current Code Owner:       Olaf Morgenstern/Oliver Wild
!                            Colin Johnson
!
!   Called from asad_spmjpdriv
!
!     Purpose.
!     --------
!     To organise the integration of the chemical rate equations using
!     the MJP implicit integrator.
!
!     Interface
!     ---------
!     Called from chemistry driver routine via the mjpdriv driver.
!
!     Method.
!     -------
!     Solves the non-linear system of equations via a Newton-Raphson
!     technique. The full Jacobian 'fj' is constructed in 'fuljac', and
!     the net change in family concentration is calculated in 'linslv'.
!     The first few iterations (currently 7) are controlled to prevent
!     very rapid initial changes leading to divergence. Convergence is
!     determined by checking that the total concentration error 'errxo'
!     is less than tolerance 'raferr' or that the total rate error 'errpl'
!     is less than tolerance 'rafpml'. If convergence is not achieved
!     after 'nrsteps', or if divergence is encountered, the routine resets
!     the family concentrations to their initial values, and exits with a
!     non-zero value for ndxraf.
!
!     Global variables
!     ----------------
!     rafpml - tolerance - set to  1.0E-10 in input file
!     rafmin - limit for first few iterations, 0.1 in input file
!     rafmax - limit for first few iterations, 1.0E+04 in input file
!     raferr - tolerance (again?) - set to 1.0E-06 in input file
!     rafbig - maximum concentration above which divergence is assumed
!     rafeps - small non-zero concentration
!
!     Local variables
!     ---------------
!     ifi           Number of ftoy iterations.
!     zf            Values of f at start of chemistry step.
!     zprf          Value of f on previous Newton-Rhapson iteration.
!     ndxraf        Convergence exit code
!     damp1         Damping factor to apply to the first iteration
!     deltt         Reciprocal time step length  (1./cdt)
!
!  Changes for whole-atmosphere chemistry:
!  a. Increase rafeps to sqrt(peps)
!  b. Disactivate crash if too many negatives occur. They are now
!     fine during the iteration.
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
      SUBROUTINE asad_spimpmjp(ndxraf,nslon,nslat, n_points)

      USE ASAD_MOD
      USE ASAD_SPARSE_VARS, ONLY : spfj, pointer2, spfuljac,            &
                              splinslv2, spresolv2

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
      INTEGER, INTENT(OUT):: ndxraf

! Local variables
      INTEGER :: jtr
      INTEGER :: jit
      INTEGER :: ifi
      INTEGER :: i
      INTEGER :: itr
      INTEGER :: nl
      INTEGER :: jl
      INTEGER :: kr
      INTEGER :: j
      INTEGER :: ip

      REAL :: ztmp
      REAL :: rafpml
      REAL :: rafmin
      REAL :: rafmax
      REAL :: raferr
      REAL :: rafbig
      REAL :: rafeps
      REAL :: deltt
      REAL :: damp1
      REAL :: errxo
      REAL :: errpl

      LOGICAL :: errfl80
      LOGICAL :: ltrig

      REAL :: zf(theta_field_size,jpctr)
      REAL :: xoo(theta_field_size,jpctr)
      REAL :: fxo(theta_field_size,jpctr)
      REAL :: zsum(jpctr)
      REAL :: tmprc(theta_field_size,jpctr)
      REAL :: bx(jpctr)
!
      COMMON /trig/ltrig
!
      rafpml=1.0e-10
      rafmin=1.0e-01
      rafmax=1.0e+04
      raferr=ptol*10.   !   care...
!
      rafbig=1.0/SQRT(peps)
      rafeps=SQRT(peps)
!
      ndxraf = 0
      errxo = 1.
      lsvjac = .FALSE.
      deltt = 1./cdt
      damp1 = 0.5
!
      nl = n_points
      jlp = 0
!
!  Save values of f at start of step and make linearised first guess
      zf = f
      WHERE (f<rafeps) f = rafeps
      CALL spfuljac(n_points)
      DO jtr=1,jpctr
        ip = pointer2(jtr,jtr)
        DO jl=1,n_points
          IF  (spfj(jl,ip) > 0.) THEN
            f(jl,jtr) = zf(jl,jtr) + cdt*fdot(jl,jtr)
          ELSE
            f(jl,jtr) = zf(jl,jtr) + (cdt*fdot(jl,jtr))                 &
                                     /(1.0-cdt*spfj(jl,ip))
          END IF
          IF (f(jl,jtr) < rafeps) f(jl,jtr) = rafeps
        END DO
      END DO

!  Start Loop - ensure mixing ratios positive and non-zero on entry
      DO jit=1,nrsteps
        ifi = 0
        IF(jit == 1) ifi=nitnr
!
! DEPENDS ON: asad_ftoy
        CALL asad_ftoy(.FALSE.,ifi, n_points)
        IF (nstst /= 0) THEN
! DEPENDS ON: asad_steady
          if(ifi == 0) CALL asad_steady( nl )
        END IF
!
        IF (ltrig) THEN
          DO jl=1,n_points
            WRITE(6,*) 'Level ',jl
            IF (jit == 1) THEN
             WRITE(6,"(1x,i2,20(1x,1pG12.4))")                          &
                         jit-1,(zf(jl,jtr),jtr=1,jpctr)
             WRITE(6,"(1x,i2,20(1x,1pG12.4))") jit-1,                   &
                         ((zf(jl,jtr)+cdt*fdot(jl,jtr)),jtr=1,jpctr)
            END IF
            WRITE(6,"(1x,i2,20(1x,1pG12.4))")                           &
                            jit-1,(f(jl,jtr),jtr=1,jpctr),              &
                            (y(jl,i),i=1,2)
          END DO
        END IF
!
! DEPENDS ON: asad_diffun
        CALL asad_diffun( nl )
!
!  Temporary prod+loss array
        tmprc(1:n_points,:) = prod(1:n_points,1:jpctr) +                &
                              slos(1:n_points,1:jpctr)
!
        IF (errxo < raferr) THEN
          IF(jit >= 51) THEN               ! Set to 50 if want LTRIG  !!
            ndxraf = 4
            WRITE (6,                                                   &
      "('Convergence problems (',i3,1x,'iter) at lon=',i3,'  lat=',i3)")&
              jit, nslon, nslat
          END IF
          GOTO 9999
        END IF

        IF (jit == nrsteps) THEN
          WRITE(6,                                                      &
          "('Convergence not achieved in impmjp (iter',i3,')  lon=',"// &
          "i3,'  lat=',i3,'; halving step')") jit, nslon, nslat
          ndxraf = 2
          f = zf
          IF(jit >= 51) THEN               ! Set to 50 if want LTRIG  !!
            ndxraf = 4
            WRITE (6,                                                   &
      "('Convergence problems (',i3,1x,'iter) at lon=',i3,'  lat=',i3)")&
              jit, nslon, nslat
          END IF
          GOTO 9999
        END IF
!
!  Calculate fxo
        fxo = (f - zf)*deltt - fdot
!
!  Test for convergence
        errfl80 = .FALSE.
        errpl = 0.
        DO jl=1,n_points
          DO jtr=1,jpctr
            IF(ABS(tmprc(jl,jtr)) > rafeps)                             &
                   errpl=MAX(errpl,ABS(fxo(jl,jtr)/tmprc(jl,jtr)))
            IF (f(jl,jtr) > rafbig)  THEN
              errfl80 = .TRUE.
            END IF
          END DO
        END DO
        IF (errfl80) THEN
          ndxraf = 3
          f = zf
          IF(jit >= 51) THEN               ! Set to 50 if want LTRIG  !!
            ndxraf = 4
            WRITE (6,                                                   &
      "('Convergence problems (',i3,1x,'iter) at lon=',i3,'  lat=',i3)")&
              jit, nslon, nslat
          END IF
          GOTO 9999
        END IF
        IF (errpl < rafpml) THEN
          IF(jit >= 51) THEN               ! Set to 50 if want LTRIG  !!
            ndxraf = 4
            WRITE (6,                                                   &
      "('Convergence problems (',i3,1x,'iter) at lon=',i3,'  lat=',i3)")&
              jit, nslon, nslat
          END IF
          GOTO 9999
        END IF
!
!  Fill in the Jacobian, or just solve if lsvjac = .true.
        IF (lsvjac) THEN
! Sparse resolve routine
          CALL spresolv2(fxo,xoo,n_points,rafeps)
        ELSE
          CALL spfuljac(n_points)
!
          IF (ltrig) THEN
            WRITE(60,*) 'Iteration ',jit
            DO jl=1,n_points
              WRITE(60,*) 'Level ',jl
              fj(jl,:,:) = 0.
              DO jtr = 1,jpctr
                DO itr = 1,jpctr
                  IF (pointer2(jtr,itr) > 0)                            &
                    fj(jl,jtr,itr) = spfj(jl,pointer2(jtr,itr))
                END DO
              END DO
              DO jtr=1,jpctr
                WRITE(60,"(1x,i2,20(1x,1pG12.4))") jtr,                 &
                   (fj(jl,jtr,itr),itr=1,jpctr)
              END DO
            END DO
          END IF

          CALL splinslv2(fxo,xoo,n_points,rafeps,rafbig)

          IF(ltrig) THEN
            DO jl=1,n_points
              WRITE(60,*) jl
              WRITE(60,"(1x,a3,20(1x,1pG12.4))")                        &
                  'fxo',(fxo(jl,jtr),jtr=1,jpctr)
              WRITE(60,"(1x,a3,20(1x,1pG12.4))")                        &
                  'fdt',(fdot(jl,jtr),jtr=1,jpctr)
              WRITE(60,"(1x,a3,20(1x,1pG12.4))")                        &
                  'del',((f(jl,jtr)-zf(jl,jtr))*deltt,jtr=1,jpctr)
              WRITE(60,"(1x,a3,20(1x,1pG12.4))")                        &
                  'f  ',(f(jl,jtr),jtr=1,jpctr)
              WRITE(60,"(1x,a3,20(1x,1pG12.4))")                        &
                  'xoo',(xoo(jl,jtr),jtr=1,jpctr)
            END DO
! DEPENDS ON: asad_fuljac
            CALL asad_fuljac(n_points)
            WRITE(60,*) jit
            DO jl=1,n_points
              DO jtr=1,jpctr
                zsum(jtr) = 0.
                DO itr=1,jpctr
                  zsum(jtr)=zsum(jtr)+fj(jl,jtr,itr)*xoo(jl,itr)
                END DO
                WRITE(60,"(1x,i2,20(1x,1pG12.4))") jtr,                 &
                        (fj(jl,jtr,itr)*xoo(jl,itr),itr=1,jpctr)
              END DO
              WRITE(60,"(1x,a3,20(1x,1pG12.4))") 'sum',                 &
                (zsum(jtr),jtr=1,jpctr)
            END DO
          END IF
!
        END IF
!
        errxo=0.
        ndxraf=0
        DO jtr=1,jpctr
          DO jl=1,n_points
!  Filter for negatives
!  Special kick for troublesome convergence
            IF (jit == 1) xoo(jl,jtr) = damp1*xoo(jl,jtr)
!  Calculate error
            xoo(jl,jtr) = MIN(MAX(xoo(jl,jtr),-rafbig),rafbig)
            IF (ABS(xoo(jl,jtr)) > 1.e-16)                              &
                errxo = MAX(errxo,ABS(xoo(jl,jtr)/                      &
                   MAX(f(jl,jtr),rafeps)))
!  New mixing ratios
            ztmp = f(jl,jtr) + xoo(jl,jtr)
!  Put limit on MAXimum correction for first few (6) iterations
            IF (jit < 7)                                                &
              ztmp = MAX(rafmin*f(jl,jtr),MIN(rafmax*f(jl,jtr),ztmp))
!  Filter negatives and zeros
            IF (ztmp == 0.) ztmp = rafeps
            IF (ztmp < 0.) THEN
              ztmp = rafeps
              ndxraf = ndxraf+1
            END IF
!  Final mixing ratios
            f(jl,jtr) = ztmp
          END DO
        END DO
! if 5 or more negatives, drop out and halve step
        IF(ndxraf > 1000) THEN
          WRITE(6,                                                      &
            "(1x,'Too many negatives (>5)  in impmjp (iter',i3,"//      &
            "')  lon=',i3,'  lat=',i3,'; halving step')")               &
            jit, nslon, nslat
          ndxraf=2
          f = zf
          IF(jit >= 51) THEN               ! Set to 50 if want LTRIG  !!
            ndxraf = 4
            WRITE (6,"(1x,'Convergence problems (',i3,1x,"//            &
             "'iter) at lon=',i3,'  lat=',i3)")                         &
              jit, nslon, nslat
          END IF
          GOTO 9999
        END IF

        ndxraf=0
      END DO
!
!  Failure to Converge - reset f's and exit
      WRITE(6,                                                        &
          "('Convergence not achieved in impmjp (iter',i3,')  lon=',"// &
          "i3,'  lat=',i3,'; halving step')") jit, nslon, nslat
      ndxraf = 2
      f = zf
      IF(jit >= 51) THEN               ! Set to 50 if want LTRIG  !!
         ndxraf = 4
         WRITE (6,                                                    &
      "('Convergence problems (',i3,1x,'iter) at lon=',i3,'  lat=',i3)")&
              jit, nslon, nslat
      END IF

 9999 CONTINUE
      RETURN
      END SUBROUTINE asad_spimpmjp
#endif
