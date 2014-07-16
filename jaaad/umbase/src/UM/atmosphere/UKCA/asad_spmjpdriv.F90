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
!   Called from ASAD_CDRIVE
!
!
!     MJPDRIV  - Driver for MJP fully implicit ODE integrator.
!
!     Michael Prather            Earth System Science
!     Oliver Wild                University of California, Irvine
!
!     ASAD: mjpdriv              Version: mjpdriv.f 1.0 04/17/97
!
!     Purpose.
!     --------
!     To organise the integration of the chemical rate equations using
!     the MJP implicit integrator.
!
!     Interface
!     ---------
!     Called from chemistry driver routine *cdrive*.
!
!     This routine assumes that all the bi-,tri-, phot- and het-
!     reaction rates have been computed prior to this routine.
!     It is also assumed that the species array, y, has been set
!     from the array, f, passed to ASAD, and that constant species
!     have been set. This can be done by calling the routine fyinit.
!
!     Method.
!     -------
!     This routine calls the MJP integrator once for each gridpoint
!     of the one-dimensional arrays passed to it, setting 'jlst' to
!     the current level in the same way as the SVODE integrator.
!     If convergence isn't achieved, the time step is halved, and the
!     integrator called again - this is continued until either
!     convergence is achieved or the minimum time step length is
!     encountered (currently 1.E-05 seconds).
!
!     Local variables
!     ---------------
!     ncst    -  Stores number of basic chemical steps 'ncsteps'
!     ctrd    -  Stores basic chemical time step length 'ctd'
!     zf      -  Stores family concentrations at beginning of call
!     ndxraf  -  Error code from the integrator:
!                   0 = successful return
!                   1 = negatives encountered
!                   2 = convergence failure after 'nrsteps' iterations
!                   3 = convergence failure due to divergence - 'NaN's
!                   4 = convergence failure (as '2') but set debugging
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
      SUBROUTINE asad_spmjpdriv(nslon,nslat,n_points)

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
      LOGICAL :: ltrig
      INTEGER :: ndxraf
      INTEGER :: ncst
      INTEGER :: iredo
      INTEGER :: jl
      INTEGER :: jtr
      INTEGER :: i
      INTEGER :: nl

      REAL :: ctrd
      REAL :: zf(theta_field_size,jpctr)
      COMMON /trig/ltrig
!
      ncst = ncsteps
      ctrd = cdt
      ltrig=.FALSE.
!
      nl = n_points
! DEPENDS ON: asad_diffun
      CALL asad_diffun( nl )
!
      iredo = 1
      zf(1:n_points,:)=f(1:n_points,:)
!
! Start iterations here.
      i = 1
      DO WHILE (i <= iredo)
! DEPENDS ON: asad_spimpmjp
        CALL asad_spimpmjp(ndxraf,nslon,nslat, n_points)

!  Debug slow convergence systems - switch this on in 'impmjp'
        IF (ndxraf == 4) THEN
          IF (ltrig) THEN
            CLOSE(60)
! DEPENDS ON: ereport
            CALL ereport('ASAD_MJPDRIV',1,                              &
              'You need to debug slow-converging system.')
          END IF
          OPEN (60,FILE='jac')
          ltrig=.TRUE.
          f(1:n_points,:)=zf(1:n_points,:)
! DEPENDS ON: asad_ftoy
          CALL asad_ftoy( .FALSE., nitfg, n_points )
! DEPENDS ON: asad_diffun
          CALL asad_diffun( nl )
          i = 1
        ELSE
!
!  Reset for failed convergence
          IF (ndxraf > 1) THEN
            ncsteps = ncsteps*2
            cdt = cdt/2.
            iredo = iredo*2
            IF(cdt < 1.0e-05) THEN
! DEPENDS ON: ereport
              CALL ereport('ASAD_MJPDRIV',2,' Time step now too short')
            END IF
            f(1:n_points,:)=zf(1:n_points,:)
! Drop out at some point - if 3 successive halvings fail
            IF (iredo >= 16) THEN
              WRITE(6,"(' Reseting column after',i3,' iterations')")    &
                    iredo
              WRITE(6,*) 'NO CONVERGENCE ',nslon,nslat
              ncsteps = ncst
              cdt = ctrd
              GOTO 9999
            END IF
! DEPENDS ON: asad_ftoy
            CALL asad_ftoy( .FALSE., nitfg, n_points )
! DEPENDS ON: asad_diffun
            CALL asad_diffun( nl )
            i = 1
          ELSE
            i = i + 1
          END IF
        END IF
      END DO

      IF(iredo > 1) THEN
        IF (iredo > 2) WRITE(6,"('   No. iterations =',i2)") iredo
        ncsteps = ncst
        cdt = ctrd
      END IF

 9999 CONTINUE
      RETURN
      END SUBROUTINE asad_spmjpdriv
#endif


