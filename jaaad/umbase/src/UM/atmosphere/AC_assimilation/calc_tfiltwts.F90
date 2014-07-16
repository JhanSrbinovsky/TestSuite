#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculate weights for a temporal filtering scheme.

      SUBROUTINE Calc_TFiltWts ( StartSec,                              &
                                                ! in
     &                           EndSec,                                &
                                                ! in
     &                           TS_len_secs,                           &
                                                ! in
     &                           ApexSec,                               &
                                                ! in
     &                           Cutoff_period,                         &
                                                ! in
     &                           SBE_period,                            &
                                                ! in
     &                           FilterType,                            &
                                                ! in
     &                           MaxNumWeights,                         &
                                                ! in
     &                           Weights,                               &
                                                ! out
     &                           NumWeights )   ! out

! Description:
!
!   Calculate weights for an Incremental Analysis Update (IAU) or
!   Temporal Digital Filtering (TDF) scheme.
!
!   FilterType must be one of:
!
!     'Uniform   ':  Uniform weights.
!
!     'Triangular':  Triangular weights. ApexSec is then the second
!                    at which the triangle's apex is to occur.
!                    ApexSec does not have to fall on a timestep.
!
!                    To avoid zero weights at the beginning and end of
!                    the filter period, the `base' of the triangle
!                    extends one timestep either side.
!
!     'Unwindowed':  Weights as determined by a straightforward Fourier
!                    analysis with Cutoff_period the nominal cut-off
!                    period in hours. For details, see Lynch and Huang
!                    (1992), MWR, 120, 1019-1034.
!
!     'LancWin   ':  As above, but with a Lanczos window designed to
!                    reduce Gibbs' oscillations in the pass band.
!
!     'Dolph     ':  Dolph filter, as descibed by Lynch (1997),
!                    MWR, 124, 655-660. SBE_period is the stop band
!                    edge period in hours.
!
!   For each filter type, some of the subroutine arguments will be
!   irrelevant:
!
!     Filter type   Irrelevant arguments
!     -----------   --------------------
!
!     'Uniform   '  ApexSec, Cutoff_period, SBE_period.
!     'Triangular'           Cutoff_period, SBE_period.
!     'Unwindowed'  ApexSec,                SBE_period.
!     'LancWin   '  ApexSec,                SBE_period.
!     'Dolph     '  ApexSec, Cutoff_period.
!
!
! Current Code Owner: Adam Clayton.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   13/04/00   Original code. Adam Clayton
!   5.5   29/01/03   Arguments in seconds rather than minutes.
!   5.5              Adam Clayton
!   6.1   20/05/04   Section 2 rearranged to remove SX6 warning.
!   6.1              Adam Clayton
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! Declarations:

      IMPLICIT NONE

! Common blocks:

#include "c_pi.h"

! Subroutine arguments:

      INTEGER,       INTENT(IN) ::                                      &
      
     &  StartSec,                                                       &
                       ! Start second of filter.
     &  EndSec,                                                         &
                       ! End   second of filter.
     &  TS_len_secs,                                                    &
                       ! Timestep length in seconds.
     &  ApexSec,                                                        &
                       ! Max weight second for triangular filter.
     &  MaxNumWeights  ! Maximum number of weights.

      REAL,          INTENT(IN) ::                                      &
      
     &  Cutoff_period,                                                  &
                       ! Cut-off period in hours.
     &  SBE_period     ! Stop band edge period in hours.

      CHARACTER(10), INTENT(IN) ::                                      &
      
     &  FilterType     ! Filter type.

      REAL,          INTENT(OUT) ::                                     &
      
     &  Weights(MaxNumWeights) ! Filter weights.

      INTEGER,       INTENT(OUT) ::                                     &
      
     &  NumWeights     ! Number of weights required for filter.

! Local variables:

      INTEGER :: m, n, j, index,                                        &
     &           ICode,                                                 &
     &           WeightNum,                                             &
                               ! Weight number.
     &           Sec,                                                   &
                               ! Second.
     &           BaseStartSec,                                          &
                               ! Start second of triangle base.
     &           BaseEndSec    ! End   second of triangle base.

      REAL ::    X, Sum,                                                &
     &           ThetaC,                                                &
                               ! Digital cut-off frequency.
     &           ThetaS,                                                &
                               ! Digital SBE     frequency.
     &           ThetaI,                                                &
     &           ThetaJ,                                                &
     &           X_zero,                                                &
                               ! Main lobe width.
     &           r,                                                     &
                               ! Ripple ratio.
     &           ChebPoly,                                              &
                               ! Chebyshev polynomial.
     &           SumWeights    ! Sum of weights.

      CHARACTER(80) :: CMessage

      CHARACTER(*)  :: RoutineName
      PARAMETER      ( RoutineName='Calc_TFiltWts' )

! External subroutines called:

      EXTERNAL EReport

!- End of header ------------------------------------------------------

!----------------------------------------------------------------------
! [1]: Check inputs and set up some variables.
!----------------------------------------------------------------------

      IF (FilterType /= 'Uniform   ' .AND.                              &
     &    FilterType /= 'Triangular' .AND.                              &
     &    FilterType /= 'Unwindowed' .AND.                              &
     &    FilterType /= 'LancWin   ' .AND.                              &
     &    FilterType /= 'Dolph     ') THEN

        ICode = 1
        WRITE (6,*) ''
        WRITE (CMessage, FMT='(A,A)')                                   &
     &    'Invalid filter type: ', FilterType
        WRITE (6,*) 'Calc_TFiltWts: ', CMessage
        WRITE (6,*) '               Valid types: Uniform, Triangular,'
        WRITE (6,*) '                            Unwindowed, LancWin,'
        WRITE (6,*) '                            Dolph.'
        WRITE (6,*) ''
! DEPENDS ON: ereport
        CALL EReport (RoutineName, ICode, CMessage)

      END IF

      IF (EndSec < StartSec) THEN
        ICode    = 1
        CMessage = 'End second before start second.'
! DEPENDS ON: ereport
        CALL EReport (RoutineName, ICode, CMessage)
      END IF

      IF (MOD(StartSec, TS_len_secs) /= 0 .OR.                          &
     &    MOD(EndSec,   TS_len_secs) /= 0) THEN

        ICode    = 1
        CMessage = 'Filter must start and end on a timestep.'
        WRITE (6,*) ''
        WRITE (6,*) 'Calc_TFiltWts: ', CMessage
        WRITE (6,*) ''
        WRITE (6,*) '  StartSec:     ', StartSec
        WRITE (6,*) '  EndSec:       ', EndSec
        WRITE (6,*) '  TS_len_secs:  ', TS_len_secs
        WRITE (6,*) ''
! DEPENDS ON: ereport
        CALL EReport (RoutineName, ICode, CMessage)

      END IF

      NumWeights = (EndSec - StartSec)/TS_len_secs + 1

      IF (NumWeights > MaxNumWeights) THEN
        ICode = 1
        WRITE (CMessage,*)                                              &
     &    'Number of weights required for filter (', NumWeights,        &
     &    ') exceeds maximum (', MaxNumWeights, ')'
! DEPENDS ON: ereport
        CALL EReport (RoutineName, ICode, CMessage)
      END IF

      IF (NumWeights == 1) THEN
        Weights(1) = 1.0
        RETURN
      END IF

      IF (MOD(NumWeights, 2) /= 1) THEN

        ICode    = 1
        CMessage = 'Filter must be centred on a timestep.'
        WRITE (6,*) ''
        WRITE (6,*) 'Calc_TFiltWts: ', CMessage
        WRITE (6,*) ''
        WRITE (6,*) '  StartSec:     ', StartSec
        WRITE (6,*) '  EndSec:       ', EndSec
        WRITE (6,*) '  TS_len_secs:  ', TS_len_secs
        WRITE (6,*) ''
! DEPENDS ON: ereport
        CALL EReport (RoutineName, ICode, CMessage)

      END IF

      IF (FilterType  ==  'Triangular') THEN

        BaseStartSec = StartSec - TS_len_secs
        BaseEndSec   = EndSec   + TS_len_secs

        IF (ApexSec  <   StartSec .OR.                                  &
     &      ApexSec  >   EndSec) THEN
          ICode    = 1
          CMessage = 'Invalid apex second.'
          WRITE (6,*) ''
          WRITE (6,*) 'Calc_TFiltWts: ', CMessage
          WRITE (6,*) ''
          WRITE (6,*) '  ApexSec:  ', ApexSec
          WRITE (6,*) '  StartSec: ', StartSec
          WRITE (6,*) '  EndSec:   ', EndSec
          WRITE (6,*) ''
! DEPENDS ON: ereport
          CALL EReport (RoutineName, ICode, CMessage)
        END IF

      END IF

      IF (FilterType == 'Unwindowed' .OR.                               &
     &    FilterType == 'LancWin   ') THEN

       ! Digital cut-off frequency:
        ThetaC = 2.0 * PI * REAL(TS_len_secs)                           &
     &         / (Cutoff_period * 3600.0)

        IF (ThetaC > PI) THEN
          ICode    = 1
          CMessage = 'Cut-off period too small.'
! DEPENDS ON: ereport
          CALL EReport (RoutineName, ICode, CMessage)
        END IF

        IF (FilterType == 'LancWin   ') m = (NumWeights - 1)/2

      END IF

      IF (FilterType == 'Dolph     ') THEN

       ! Digital stop band edge frequency:
        ThetaS = 2.0 * PI * REAL(TS_len_secs)                           &
     &         / (SBE_period * 3600.0)

        IF (ThetaS > PI) THEN
          ICode    = 1
          CMessage = 'Stop band edge period too small.'
! DEPENDS ON: ereport
          CALL EReport (RoutineName, ICode, CMessage)
        END IF

       ! Main lobe width:
        X_zero = 1.0 / ( COS(ThetaS/2.0) )

       ! Ripple ratio:
! DEPENDS ON: chebpoly
        r = 1.0 / ChebPoly(X_zero, NumWeights-1)

       m = (NumWeights - 1)/2
       n = NumWeights

      END IF

!----------------------------------------------------------------------
! [2]: Calculate un-normalised filter weights.
!----------------------------------------------------------------------

      Weights(1:MaxNumWeights) = 0.0

!----------------------------------------------------------------------
! [2.1]: Uniform filter.
!----------------------------------------------------------------------

      IF (FilterType == 'Uniform   ') THEN
        DO WeightNum = 1, NumWeights
          Weights(WeightNum) = 1.0
        END DO
      END IF

!----------------------------------------------------------------------
! [2.2]: Triangular filter.
!----------------------------------------------------------------------

      IF (FilterType == 'Triangular') THEN
        DO WeightNum = 1, NumWeights

          Sec = StartSec + (WeightNum - 1) * TS_len_secs

          IF (Sec <= ApexSec) THEN
            Weights(WeightNum) = REAL(Sec     - BaseStartSec) /         &
     &                           REAL(ApexSec - BaseStartSec)
          ELSE
            Weights(WeightNum) = REAL(BaseEndSec - Sec)       /         &
     &                           REAL(BaseEndSec - ApexSec)
          END IF

        END DO
      END IF

!----------------------------------------------------------------------
! [2.3]: Unwindowed and Lanczos-windowed filters.
!----------------------------------------------------------------------

      IF (FilterType == 'Unwindowed' .OR.                               &
     &    FilterType == 'LancWin   ') THEN
        DO WeightNum = 1, NumWeights

          Index = WeightNum - (NumWeights + 1) / 2

          IF (Index /= 0) THEN
            Weights(WeightNum) = SIN(REAL(Index) * ThetaC) /            &
     &                              (REAL(Index) * PI)
          ELSE
            Weights(WeightNum) = ThetaC / PI
          END IF

        END DO
      END IF

      IF (FilterType == 'LancWin   ') THEN
        DO WeightNum = 1, NumWeights

          Index = WeightNum - (NumWeights + 1) / 2

          IF (Index /= 0)                                               &
     &      Weights(WeightNum) = Weights(WeightNum) *                   &
     &                           SIN(Index * PI / REAL(m+1)) /          &
     &                              (Index * PI / REAL(m+1))

        END DO
      END IF

!----------------------------------------------------------------------
! [2.4]: Dolph filter.
!----------------------------------------------------------------------

      IF (FilterType == 'Dolph     ') THEN
        DO WeightNum = 1, NumWeights

          Index = WeightNum - (NumWeights + 1)/2

          ThetaI = 2.0 * PI * Index / REAL(n)

          Sum = 0.0

          DO j = 1, m
            ThetaJ = 2.0 * PI * j / REAL(n)
            X      = X_zero * COS(ThetaJ / 2.0)
! DEPENDS ON: chebpoly
            Sum    = Sum + ChebPoly(X, 2*m) * COS(j * ThetaI)
          END DO

          Weights(WeightNum) = (2.0 * r * Sum + 1.0) / REAL(n)

        END DO
      END IF

!----------------------------------------------------------------------
! [3]: Normalise weights.
!----------------------------------------------------------------------

      SumWeights = 0.0

      DO WeightNum = 1, NumWeights
        SumWeights = SumWeights + Weights(WeightNum)
      END DO

      DO WeightNum = 1, NumWeights
        Weights(WeightNum) = Weights(WeightNum) / SumWeights
      END DO


      RETURN
      END SUBROUTINE Calc_TFiltWts


!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
!+ Chebyshev polynomial function.

#endif
