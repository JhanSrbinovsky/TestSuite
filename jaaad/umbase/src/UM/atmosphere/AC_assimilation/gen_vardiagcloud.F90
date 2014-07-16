#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Diagnostic liquid and frozen cloud scheme for deriving IAU increments

SUBROUTINE Gen_VarDiagCloud(field_size,     & ! in
                            p_theta_levels, & ! in
                            RHc,            & ! in
                            IncrementIce,   & ! in
                            qT,             & ! inout
                            CMessage,       & ! inout
                            ICode,          & ! inout
                            CL,             & ! out   (optional)
                            qCL,            & ! inout (optional)
                            dCLdp,          & ! out   (optional)
                            dCLdqT,         & ! out   (optional)
                            dCLdT,          & ! out   (optional)
                            dqCLdp,         & ! out   (optional)
                            dqCLdqT,        & ! out   (optional)
                            dqCLdT,         & ! out   (optional)
                            CF,             & ! out   (optional)
                            qCF,            & ! inout (optional)
                            qCFmax,         & ! inout (optional)
                            dCFdp,          & ! out   (optional)
                            dCFdqcf,        & ! out   (optional)
                            dCFdT,          & ! out   (optional)
                            dqCFdqCL,       & ! out   (optional)
                            dqCFdT,         & ! out   (optional)
                            T,              & ! in    (optional)
                            BGqcl,          & ! in    (optional)
                            BGqcf,          & ! in    (optional)
                            BGT,            & ! in    (optional)
                            TL)               ! inout (optional)


! Description:
!
!   Documentation located in VAR Scientific Documentation Paper 61
!
!   NOTE 1: Only request gradients when calling from VAR or OPS and
!           the temperature input is T; not TL.
!
!   NOTE 2: This code is intended to mirror the subroutine of the same
!           name held in GenMod_Utilities. Differences between UM and
!           GEN versions are commented in-line. These are primarily
!           differences in CALLs to subroutines for calculating
!           saturated vapour pressure and a commenting-out of linear
!           code in the UM-only version. At some future time, it may be
!           possible for the UM, VAR & OPS all to access this code from
!           a single source file.
!
! Current Code Owner: Martin Sharpe.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.1   17/08/04   Original code. Martin Sharpe
!   6.4   13/12/06   Include code for ice incrementing. Martin Sharpe
!
! Code Description:
!   Language: FORTRAN 90
!
! Declarations:

! When copying this subroutine into GenMod_Utilities, add the
! following USE statements:
!
! USE GenMod_UMConstants, ONLY : &
!     Cp,                        &
!     Epsilon,                   &
!     Lc

IMPLICIT NONE

! Common blocks:
! When copying this subroutine into GenMod_Utilities, remove these
! "#include" statements, and insert the following interface:
!
! INCLUDE 'Var_LTinterp.interface'

#include "parvars.h"
#include "c_lheat.h"
#include "c_r_cp.h"

! Subroutine arguments:

INTEGER,       INTENT(IN)              :: &
    field_size

REAL,          INTENT(IN)              :: &
    p_theta_levels (field_size),          &
    RHc            (field_size)

LOGICAL,       INTENT(IN)              :: &
    IncrementIce

REAL,          INTENT(IN),    OPTIONAL :: &
    T         (:),                        &
    BGqcl     (:),                        &
    BGqcf     (:),                        &
    BGT       (:)

INTEGER,       INTENT(INOUT)           :: &
    ICode

REAL,          INTENT(INOUT)           :: &
    qT         (field_size)

REAL,          INTENT(INOUT), OPTIONAL :: &
    qCl            (:),                   &
    qCF            (:),                   &
    qCFmax         (:),                   &
    TL             (:)

CHARACTER(len=*), INTENT(INOUT)        :: &
    CMessage

REAL, INTENT(OUT),   OPTIONAL          :: &
    CL             (:),                   &
    dCLdp          (:),                   &
    dCLdqT         (:),                   &
    dCLdT          (:),                   &
    dqCLdp         (:),                   &
    dqCLdqT        (:),                   &
    dqCLdT         (:),                   &
    CF             (:),                   &
    dCFdp          (:),                   &
    dCFdqcf        (:),                   &
    dCFdT          (:),                   &
    dqCFdqCL       (:),                   &
    dqCFdT         (:)

! Local variables:

INTEGER           ::          &
    calc,                     &
    loop,                     &
    maxLoops,                 &
    statL,                    &
    statF

LOGICAL           ::          &
    oscillate(field_size),    &
    converged(field_size)

! When copying this subroutine into GenMod_Utilities,
! replace declaration of p_temp & qSat_temp as (field_size,1,1)
REAL              ::          &
    A,                        &
    accuracy0,                &
    accuracy    (field_size), &
    aL          (field_size), &
    dqcf        (field_size), &
    delta       (field_size), &
    p_temp      (field_size), &
    qCL_calc    (field_size), &
    qCL_1       (field_size), &
    qCL_2       (field_size), &
    qCL_3       (field_size), &
    qSatW       (field_size), &
    qSatW_minus (field_size), &
    qSatW_plus  (field_size), &
    qSat_temp   (field_size), &
    Qc          (field_size), &
    QN          (field_size), &
    Tmod        (field_size), &
    Tmod_temp   (field_size)

REAL, ALLOCATABLE ::   &
    alphaL        (:), &
    dPsatdT_minus (:), &
    dPsatdT_plus  (:), &
    f_of_T        (:), &
    f_of_T_plus   (:), &
    fM            (:), &
    fPlus         (:), &
    fPlusM        (:), &
    qCL_4         (:), &
    qCFmin        (:), &
    qCFmax0       (:), &
    qCFin         (:), &
    qCF_1         (:), &
    qCF_2         (:), &
    qCF_3         (:), &
    qCL2_1        (:), &
    qCL2_2        (:), &
    qCl2_3        (:), &
    B             (:), &
    C             (:), &
    D             (:), &
    E             (:), &
    F             (:), &
    G             (:), &
    GCx           (:), &
    H             (:), &
    I             (:), &
    J             (:), &
    K             (:), &
    L             (:), &
    M             (:), &
    N             (:), &
    O             (:), &
    V             (:), &
    W             (:), &
    Y             (:), &
    Z             (:), &
    multiplier    (:)

REAL :: test(field_size)

LOGICAL           :: assignPseudoBG
LOGICAL           :: assign_qCFmax

CHARACTER(len=*), PARAMETER :: RoutineName='Gen_VarDiagCloud'

!- End of header ------------------------------------------------------

ICode = 0
statL = 0
statF = 0

! Check that one temperature vector is passed in
IF (.NOT. PRESENT(T) .AND. .NOT. PRESENT(TL)) THEN

  ICode    = 1
  CMessage = 'Error from Gen_VarDiagCloud: No temperature data'

ELSE IF (PRESENT(T) .AND. PRESENT(TL)) THEN

  ICode    = 1
  CMessage = 'Error from Gen_VarDiagCloud: T & TL supplied'

ELSE

!----------------------------------------------------------------------
! [1]: Set fixed and initial values.
!----------------------------------------------------------------------

  ! Fractional accuracy to which cloud water calculation
  accuracy  = 0.001
  accuracy0 = accuracy(1)

  ! Initial qcf increment
  dqcf      = 0.0

  ! Determine whether to calculate initial or final values
  IF (PRESENT(BGqcl) .AND. PRESENT(BGqcf) .AND.                   &
      PRESENT(qCL)   .AND. PRESENT(qCF)) THEN

    assignPseudoBG = (ALL(BGqcl == qCL)                   .AND.   &
                      ALL(BGqcf == qCF)                   .AND.   &
                      IncrementIce)

  ELSE

    assignPseudoBG = .FALSE.

  END IF

  ! Allocate extra array if qCFmax is to be calculated
  IF (PRESENT(qCFmax)) THEN

    assign_qCFmax = (ALL(qT == qCFmax))
    IF (assign_qCFmax) ALLOCATE ( qCL_4 (field_size) )

  ELSE

    assign_qCFmax = .FALSE.

  END IF

  ! Setup maximum number of iterations for cloud water calculation
  ! If calculating incremented qcf, allocate arrays, setup initial
  !   values and calculate fM & fMPlus
  IF (.NOT. AssignPseudoBG .AND. IncrementIce   .AND.             &
            PRESENT(BGqcl) .AND. PRESENT(BGqcf) .AND.             &
            PRESENT(qCL)   .AND. PRESENT(qCF)   .AND.             &
            .NOT. assign_qCFmax) THEN

    maxLoops = 150

    ALLOCATE ( f           (field_size) )
    ALLOCATE ( f_of_T      (field_size) )
    ALLOCATE ( f_of_T_plus (field_size) )
    ALLOCATE ( fM          (field_size) )
    ALLOCATE ( fPlus       (field_size) )
    ALLOCATE ( fPlusM      (field_size) )
    ALLOCATE ( M           (field_size) )
    ALLOCATE ( qCFin       (field_size) )
    ALLOCATE ( qCFmin      (field_size) )
    ALLOCATE ( qCFmax0     (field_size) )
    ALLOCATE ( qCF_1       (field_size) )
    ALLOCATE ( qCF_2       (field_size) )
    ALLOCATE ( qCF_3       (field_size) )
    ALLOCATE ( qCL2_1      (field_size) )
    ALLOCATE ( qCL2_2      (field_size) )
    ALLOCATE ( qCL2_3      (field_size) )

    qCFmin    = 0.0
    qCFin     = qCF
    qCFmax0   = qCFmax
    WHERE (qcfMax /= 0.0)
      ! Ensure that qCF is in the allowed range
      qCF       = MAX(MIN(qcfMax,qCFin),qcfMin)
      dqcf      = qCF - qCFin
    ELSEWHERE
      ! qCF incrementing not possible
      dqcf = 0.0
    ENDWHERE
    test      = 1.0
    qCF_1     = qCF
    qCF_2     = qCF
    qCF_3     = qCF
    qCL2_1    = 0.0
    qCL2_2    = 0.0
    qCL2_3    = 0.0

    WHERE (BGqcl /= 0.0 .OR. BGqcf /= 0.0)

      f   = BGqcf / (BGqcf + BGqcl)

    ELSEWHERE

      f   = 1.0

    ENDWHERE

    f_of_T      = 0.5 * (1.0 - TANH((BGT - 260.65) * 0.18))
    f_of_T_plus = 0.5 * (1.0 - TANH((T   - 260.65) * 0.18))

    WHERE (BGqcf == 0.0)

      f = 0.0
      fPlus = 0.0
      M = 0.0

      WHERE(f_of_T_plus - f_of_T > 0.0)

        fPlus = f_of_T_plus - f_of_T
        M     = 1.0

      ENDWHERE

    ELSEWHERE

      M     = (BGqcf / (BGqcf + qCL)) / f
      fPlus = f + MIN((MIN(1.0, 1.0 / M) - f), f) *               &
                    (f_of_T_plus - f_of_T)

    ENDWHERE


    WHERE (BGqcl == 0.0 .AND. BGqcf == 0.0)

      f     = 0.0
      fPlus = 0.0
      M     = 0.0

      WHERE (f_of_T_plus - f_of_T > 0.0)

        fPlus = f_of_T_plus - f_of_T
        M     = 1.0

      ENDWHERE

    ELSEWHERE (BGqcl == 0.0)

      f     = 1.0
      fPlus = 1.0
      ! Where qCL is too small to make M /= 1.0 or machine
      !   precision results in unphysical values
      WHERE (f / (BGqcf / (BGqcf + qCL)) <= 1.0 .OR.              &
             (BGqcf / (BGqcf + qCL)) / f >= 1.0 .OR.              &
             (BGqcf / (BGqcf + qCL)) - f >= 0.0 .OR.              &
             f - (BGqcf / (BGqcf + qCL)) <= 0.0)

        M     = 1.0

      ENDWHERE

    ENDWHERE

    fM     = MIN(MAX(f * M, 0.0), 1.0)
    fPlusM = MIN(MAX(fPlus * M, 0.0), 1.0)

    WHERE ((1.0 - fM) * (1.0 - fM) == 0.0 .OR. qCL == 0.0   .OR.  &
           (1.0 - fPlusM) * (1.0 - fPlusM) == 0.0           .OR.  &
           BGqCF == BGqCF + accuracy * qCL * (1.0-accuracy))

      ! Ensure that qcf calculation does not cause divide by zero
      !   or imply increments where initial qcl is very small
      test = HUGE(A)
      qCF=qCFin
      dqcf = 0.0

    ENDWHERE

  ELSE

    maxLoops = 50

  END IF

  IF (PRESENT(TL)) THEN

! DEPENDS ON: qsat_wat
    CALL QSAT_WAT (qSatW, TL, p_theta_levels, field_size)

! When copying this subroutine into GenMod_Utilities replace this
! call with the following:
!
!    CALL Var_LTinterp(                                            &
!        RESHAPE(TL,(/field_size,1,1/)), ! in                      &
!        qSat_temp,                      ! inout                   &
!        WaterOnly=.TRUE.)               ! in
!    qSatW = (Epsilon * RESHAPE(qSat_temp,(/field_size/))) /       &
!            p_theta_levels

    Tmod = TL

  ELSE IF (PRESENT(T)) THEN

! DEPENDS ON: qsat_wat
    CALL QSAT_WAT (qSatW_plus,     T+0.05,                        &
                    p_theta_levels, field_size)
! DEPENDS ON: qsat_wat
    CALL QSAT_WAT (qSatW_minus,    T-0.05,                        &
                    p_theta_levels, field_size)

    aL = 1.0 /                                                    &
        (1.0 + Lc * ((qSatW_plus - qSatW_minus) / 0.1) / Cp)

! When copying this subroutine into GenMod_Utilities replace these
! calls & the definition of aL with the following:
!
!    ALLOCATE ( alphaL        (field_size) )
!    ALLOCATE ( dPsatdT_plus  (field_size) )
!    ALLOCATE ( dPsatdT_minus (field_size) )
!
!    CALL Var_LTinterp(                                            &
!        RESHAPE(T+0.05,(/field_size,1,1/)), ! in                  &
!        qSat_temp,                          ! inout               &
!        p_temp,                             ! inout               &
!        .TRUE.,                             ! in                  &
!        WaterOnly=.TRUE.)                   ! in
!    qSatW_plus = Epsilon * RESHAPE(qSat_temp,(/field_size/)) /    &
!                 p_theta_levels
!    dPsatdT_plus = RESHAPE(p_temp,(/field_size/))
!
!    CALL Var_LTinterp(                                            &
!        RESHAPE(T-0.05,(/field_size,1,1/)), ! in                  &
!        qSat_temp,                          ! inout               &
!        p_temp,                             ! inout               &
!        .TRUE.,                             ! in                  &
!        WaterOnly=.TRUE.)                   ! in
!    qSatW_minus = Epsilon * RESHAPE(qSat_temp,(/field_size/)) /   &
!                  p_theta_levels
!    dPsatdT_minus = RESHAPE(p_temp,(/field_size/))
!
!    alphaL = (qSatW_plus - qSatW_minus) / 0.1
!    aL     = 1.0 / (1.0 + Lc * alphaL / Cp)

    Tmod = T

  END IF

!----------------------------------------------------------------------
! [2]: Calculate qcl & derived quantities
!----------------------------------------------------------------------

  converged = .FALSE.
  loop      = 0
  oscillate = .FALSE.
  qCL_calc  = 0.0
  qCL_1     = 0.0
  qCL_2     = 0.0
  qCL_3     = 0.0
  IF (assign_qCFmax) qCL_4 = 0.0

  DO WHILE (COUNT(converged) + COUNT(oscillate) < field_size .AND.    &
             loop < maxLoops)

    loop = loop + 1
    calc = field_size - COUNT(converged)

    Tmod_temp(1:calc) = PACK(Tmod, .NOT. converged)

    IF (PRESENT(T)) THEN

      p_temp   (1:calc) = PACK(p_theta_levels, .NOT. converged)
! DEPENDS ON: qsat_wat
      CALL QSAT_WAT (qSat_temp(1:calc),                           &
                      Tmod_temp(1:calc),                          &
                      p_temp   (1:calc),                          &
                      calc)
      qSatW = UNPACK(qSat_temp(1:calc), .NOT. converged, qSatW)

! When copying this subroutine into GenMod_Utilities,
! replace this conditional code with the following:
!
!      CALL Var_LTinterp(                                          &
!          RESHAPE(Tmod_temp(1:calc),(/calc,1,1/)), ! in           &
!          qSat_temp(1:calc,:,:),                   ! inout        &
!          WaterOnly=.TRUE.)                        ! in
!      qSatW = UNPACK(RESHAPE(qSat_temp(1:calc,:,:),(/calc/)),     &
!                           .NOT. converged, qSatW)
!      WHERE (.NOT. converged)
!
!        qSatW = Epsilon * qSatW / p_theta_levels
!
!      ENDWHERE

    ELSE IF (PRESENT(TL)) THEN

      p_temp   (1:calc) = PACK(p_theta_levels, .NOT. converged)
! DEPENDS ON: qsat_wat
      CALL QSAT_WAT (qSat_temp(1:calc),                           &
                      Tmod_temp(1:calc)+0.05,                     &
                      p_temp   (1:calc),                          &
                      calc)
      qSatW_plus = UNPACK(qSat_temp(1:calc),                      &
            .NOT. converged, qSatW_plus)
! DEPENDS ON: qsat_wat
      CALL QSAT_WAT (qSat_temp(1:calc),                           &
                      Tmod_temp(1:calc)-0.05,                     &
                      p_temp   (1:calc),                          &
                      calc)
      qSatW_minus = UNPACK(qSat_temp(1:calc),                     &
          .NOT. converged, qSatW_minus)

      WHERE (.NOT. converged)

        aL = 1.0 /                                                &
            (1.0 + Lc * ((qSatW_plus - qSatW_minus) / 0.1) / Cp)

      ENDWHERE

! When copying this subroutine into GenMod_Utilities,
! replace this conditional code with the following:
!
!      CALL Var_LTinterp(                                          &
!          RESHAPE(Tmod_temp(1:calc)+0.05,(/calc,1,1/)), ! in      &
!          qSat_temp(1:calc,:,:),                        ! inout   &
!          p_temp   (1:calc,:,:),                        ! inout   &
!          .TRUE.,                                       ! in      &
!          WaterOnly=.TRUE.)                             ! in
!      qSatW_plus  = UNPACK(Epsilon *                              &
!          RESHAPE(qSat_temp(1:calc,:,:),(/calc/)),                &
!           .NOT. converged, qSatW_plus)
!
!      CALL Var_LTinterp(                                          &
!          RESHAPE(Tmod_temp(1:calc)-0.05,(/calc,1,1/)), ! in      &
!          qSat_temp(1:calc,:,:),                        ! inout   &
!          p_temp   (1:calc,:,:),                        ! inout   &
!          .TRUE.,                                       ! in      &
!          WaterOnly=.TRUE.)                             ! in
!      qSatW_minus = UNPACK(Epsilon *                              &
!          RESHAPE(qSat_temp(1:calc,:,:),(/calc/)),                &
!          .NOT. converged, qSatW_minus)
!
!      WHERE (.NOT. converged)
!
!        qSatW_plus  = qSatW_plus  / p_theta_levels
!        qSatW_minus = qSatW_minus / p_theta_levels
!        aL = 1.0 /
!            (1.0 + Lc * ((qSatW_plus - qSatW_minus) / 0.1) / Cp)
!
!      ENDWHERE

    END IF

    oscillate = .FALSE.

    ! When assigning qCFmax, set qT such that we find the point
    !   where qT = qcl
    IF (assign_qCFmax) THEN
      IF (loop == 0) THEN
        qT  = 0.0
      ELSE
        WHERE (qCL_calc > 0.0)
          qT = qCL_calc
        ELSEWHERE
          qT = 0.0
        ENDWHERE
      END IF
    END IF

    WHERE(.NOT. converged)

      Qc = aL * ((qT-dqcf) - qSatW)
      delta    = qSatW * aL * (1.0 - RHc) / 4.0
      QN       = Qc / delta
      QN       = MAX(QN, -ALOG(HUGE(QN(1)))*0.99)
      qCL_calc = Qc + delta * ALOG( EXP(-QN) + 1.0 )

    ENDWHERE

    IF (IncrementIce) THEN

      WHERE(.NOT. converged)

        qCL_calc = MAX(qCL_calc, 0.0)

      ENDWHERE

    END IF

    IF (PRESENT(TL)) THEN

      WHERE(.NOT. converged)

        Tmod = TL + Lc * qCL_calc / Cp

      ENDWHERE

    ELSE IF (PRESENT(T)) THEN

      WHERE(.NOT. converged)

       Tmod = T  - Lc * qCL_calc / Cp

      ENDWHERE

    ENDIF

    ! Apply convergence criteria
    IF (IncrementIce .AND. .NOT. assign_qCFmax) THEN

      WHERE     (.NOT. converged    .AND.                         &
                 (ABS(qCL_1 - qCL_calc) <                         &
                     accuracy * MIN(qCL_calc, qCL_1) .OR.         &
                  qCL_calc == 0.0))

        converged = .TRUE.

      ELSEWHERE (.NOT. converged   .AND.                          &
                 qCL_2 == qCL_calc .AND.                          &
                 qCL_3 == qCL_1)

        oscillate = .TRUE.

      ENDWHERE

      WHERE(oscillate)

        qCL_calc  = (qCL_calc + qCL_1) * 0.5
        oscillate = .FALSE.
        converged = .TRUE.

      ENDWHERE

    ELSE IF (assign_qCFmax) THEN

      WHERE (.NOT. converged .AND.                                &
             ((qCL_calc == qCL_1 .AND.                            &
               qCL_1    == qCL_2 .AND.                            &
               qCL_2    == qCL_3         ) .OR.                   &
              (qCL_2    == qCL_calc .AND.                         &
               qCL_3    == qCL_1         ) .OR.                   &
              (qCL_3    == qCL_calc      ) .OR.                   &
              (qCL_4    == qCL_calc      )))

        converged = .TRUE.
        WHERE (qcfMax < qCL_calc)

          qcfMax = 0.0

        ELSEWHERE

          qcfMax = MAX(qcfMax - qT, BGqcf) * (1.0 + accuracy)

        ENDWHERE

      ENDWHERE

      IF (loop == maxLoops) THEN

        WHERE (.NOT. converged)

          converged = .TRUE.
          WHERE (qcfMax < qCL_calc)

            qcfMax = 0.0

          ELSEWHERE

            qcfMax = MAX(qcfMax - qT, BGqcf) * (1.0 + accuracy)

          ENDWHERE

        ENDWHERE

      ENDIF


    ELSE

      WHERE     (.NOT. converged   .AND.                          &
                 ABS(qCL_1 - qCL_calc) <                          &
                 accuracy * MIN(qCL_calc, qCL_1))

        converged = .TRUE.

      ELSEWHERE (.NOT. converged   .AND.                          &
                 qCL_2 == qCL_calc .AND.                          &
                 qCL_3 == qCL_1)

        oscillate = .TRUE.

      ENDWHERE

    ENDIF

    IF (assign_qCFmax) THEN

      WHERE(.NOT. converged)

        qCL_4 = qCL_3

      ENDWHERE

    END IF

    WHERE(.NOT. converged)

      qCL_3    = qCL_2
      qCL_2    = qCL_1
      qCL_1    = qCL_calc

    ENDWHERE

    ! Control iteration for calculating incremented cloud water
    IF (.NOT. AssignPseudoBG .AND. IncrementIce   .AND.           &
              PRESENT(BGqcl) .AND. PRESENT(BGqcf) .AND.           &
              PRESENT(qCL)   .AND. PRESENT(qCF)   .AND.           &
        .NOT. (assign_qCFmax .OR. loop == maxloops)) THEN

      ! test == HUGE(A) serves as an indicator that the qcf
      !   calculation is converged
      IF (ANY(converged)) THEN
        WHERE (converged .AND. loop /= maxloops)

          WHERE(test /= HUGE(A))

            WHERE (qcfMax == 0.0)

              qCF  = qCFin
              test = HUGE(A)
              converged = .FALSE.

            ELSEWHERE

              ! Set test = the next value of qcf
              test = qCL_calc * (fPlusM / (1.0 - fPlusM))

            ENDWHERE

          ELSEWHERE(test == HUGE(A) .AND.                           &
                    qCL_calc == 0.0 .AND. qCF == 0.0)

            ! Ensure that qcf increments are zero if qCL_calc is zero
            qCF = qCFin
            WHERE (qCFin /= 0.0)

              converged = .FALSE.

            ENDWHERE

          ENDWHERE

          WHERE (test == HUGE(A))

            ! Do nothing

          ELSEWHERE ((test == 0.0 .AND. qCf == 0.0) .OR.            &
                 (test == test +                                    &
                          accuracy0 * qCL_calc * (1.0 - accuracy0)  &
                  .AND. qCF  == qCF  +                              &
                          accuracy0 * qCL_calc * (1.0 - accuracy0)  &
                  .AND. test /= 0.0 .AND. qCF /= 0.0))

            ! Ensure that qcf increments are zero where qcl is tiny
            WHERE (qCFin /= 0.0)

              converged = .FALSE.

            ENDWHERE
            qCF  = qCFin
            test = HUGE(A)

          ELSEWHERE (                                               &
                 (qCL_calc == qCL_calc +                            &
                          accuracy0 * test * (1.0 - accuracy0)      &
                  .AND. qCL2_1  == qCL2_1  +                        &
                          accuracy0 * test * (1.0 - accuracy0)      &
                  .AND. qCL_calc /= 0.0 .AND. qCL2_1 /= 0.0))

            ! qcf calculation is converged if qcf in cannot be seen
            !   by qcl
            test = HUGE(A)

          ELSEWHERE

            qCF = qCL_calc * (fPlusM / (1.0 - fPlusM))
            WHERE((qCF == qCF_1 .AND. qCF /= 0.0 .AND.              &
                   qCF /= qcfMax0)               .OR.               &
                  ((ABS(qCL_calc - qCL2_1) < accuracy0 *            &
                    MIN(qCL2_1, qCL_calc) .OR.                      &
                    qCL_calc == qCL2_1)          .AND.              &
                   (ABS(qCF - qCF_1) < accuracy0 *                  &
                    MIN(qCF_1, qCF) .OR. qCF == qCF_1)))

              ! qcf calculation is converged if qcf and qcl changes
              !   are within predefined fractional accuracy
              test=HUGE(A)

            ELSEWHERE (qCL_calc /= 0.0    .AND. qCL2_1 /= 0.0 .AND. &
                       (qCL_calc == qCL2_2                          &
                  .OR. qCL_calc == qCL2_1))

              ! qcl calculation is oscillating
              WHERE(qcl_calc == qcl_1)

                ! qcf is converged if oscillation is null
                test=HUGE(A)
                qCL_calc = 0.5*(qCL_calc + qCL2_1)
                qCF = qCL_calc * (fPlusM / (1.0 - fPlusM))

              ELSEWHERE

                ! tighten accuracy requirement
                accuracy  = accuracy * 0.1
                converged = .FALSE.
                qCF_3  = qCF_2
                qCF_2  = qCF_1
                qCF_1  = qCF
                qCL2_3 = qCL2_2
                qCL2_2 = qCL2_1
                qCL2_1 = qCL_calc

              ENDWHERE

            ELSEWHERE((ABS(qCFmax - qCFmin) <                       &
                       accuracy * MIN(qCFmax, qCFmin)   .OR.        &
                       qCFmax == qCFmin) .OR.                       &
                       (qCFmax < qCFmax0 .AND. qcfMin /= 0.0 .AND.  &
                        ABS(qCFMax / (1.0 + Accuracy) -             &
                            qCFmin * (1.0 + Accuracy)) <            &
                        accuracy * (qCFmax / (1.0 + Accuracy) +     &
                                    qCFmin * (1.0 + Accuracy))))

              ! Bounds on qcf calculation are within accuracy
              WHERE (qCF > qCFmax .OR. qCF < qCFmin)

                qCF = 0.5*(qCFmax+qCFmin)

              ENDWHERE
              qCF_3     = qCF
              qCF_2     = qCF
              qCF_1     = qCF
              qCL2_3    = qCL2_2
              qCL2_2    = qCL2_1
              qCL2_1    = qCL_calc
              converged = .FALSE.
              accuracy  = 0.1 * accuracy

            ELSEWHERE(((qCF < qCF_1 .AND. qCF_1 > qCF_2) .OR.       &
                       (qCF > qCF_1 .AND. qCF_1 < qCF_2)     ).AND. &
                      ABS(qCF - qCF_1) > ABS(qCF_1 - qCF_2)       )

              ! qcf solution diverging
              qCF = MAX(MIN(qCF,qCFmax),qCFmin)
              qCF = qCF_2 + (qCF_1 - qCF_2) * (qCF_1 - qCF_2) /     &
                         ((qCF_1 - qCF_2) + (qCF_1 - qCF ))
              WHERE (qCF < qCF_1)

               qCFmax=MIN(qCFmax,qCF_1*(1.0+Accuracy))
               qcfMin=MAX(qcfMin,qCF_2/(1.0+Accuracy))

              ENDWHERE
              WHERE (qCF > qCF_1)

                qCFmax=MIN(qCFmax,qCF_2*(1.0+Accuracy))
                qCFmin=MAX(qCFmin,qCF_1/(1.0+Accuracy))

              ENDWHERE
              qCF_3     = qCF
              qCF_2     = qCF
              qCF_1     = qCF
              qCL2_3    = qCL2_2
              qCL2_2    = qCL2_1
              qCL2_1    = qCL_calc
              converged = .FALSE.

            ELSEWHERE (((qCF < qCF_1 .AND. qCF_1 > qCF_2) .OR.      &
                        (qCF > qCF_1 .AND. qCF_1 < qCF_2)))

              ! qcf solution converging
              WHERE (qCF < qCF_1)

                qcfMax=MIN(qcfMax,qCF_1*(1.0+Accuracy))
                qcfMin=MAX(qcfMin,qCF_2/(1.0+Accuracy))

              ENDWHERE
              WHERE (qCF > qCF_1)

                qcfMAX=MIN(qcfMax,qCF_2*(1.0+Accuracy))
                qcfMin=MAX(qcfMin,qCF_1/(1.0+Accuracy))

              ENDWHERE
              qCF = qCF_2 + (qCF_1 - qCF_2) * (qCF_1 - qCF_2) /     &
                         ((qCF_1 - qCF_2) + (qCF_1 - qCF ))
              qCF_3     = qCF
              qCF_2     = qCF
              qCF_1     = qCF
              qCL2_3    = qCL2_2
              qCL2_2    = qCL2_1
              qCL2_1    = qCL_calc
              converged = .FALSE.

            ELSEWHERE(((qCF < qCF_1 .AND. qCF_1 < qCF_2).OR.        &
                      (qCF > qCF_1 .AND. qCF_1 > qCF_2)))

              ! qcf solution increasing or decreasing
              WHERE (qCF < qCF_1)

                qCFmax=MIN(qCFmax,qCF_1*(1.0+Accuracy))

              ELSEWHERE (qCF > qCF_1)

                qCFmin=MAX(qCFmin,qCF_1/(1.0+Accuracy))

              ENDWHERE
              qCF_3     = qCF_2
              qCF_2     = qCF_1
              qCF_1     = qCF
              qCL2_3    = qCL2_2
              qCL2_2    = qCL2_1
              qCL2_1    = qCL_calc
              converged = .FALSE.

            ELSEWHERE(qCL2_1 == qCL2_3 .AND.                        &
                      qCL_calc == qCL2_2 .AND.                      &
                      (qCFmax == qCFMin .AND. qCFmax /= qCFmax0     &
                                        .AND. qCFmax /= 0.0         &
                                        .AND. qCF == qCFmax))

              ! qcf solution converged to machine accuracy
              test   = HUGE(A)
              qCF_3  = qCF
              qCF_2  = qCF
              qCF_1  = qCF
              qCL2_3 = qCL2_2
              qCL2_2 = qCL2_1
              qCL2_1 = qCL_calc

            ELSEWHERE

              WHERE(.NOT. (fPlusM == 0.0) .AND.                     &
                    .NOT. (qCF /= qCF_1 .AND. qCF_1 == qCF_2))

                ! Preceding logic should make this section unused
                test      = HUGE(A)
                qCF       = qCFin
                converged = .FALSE.

              ELSEWHERE

                converged = .FALSE.
                WHERE (qCF_1 < qcfMin .OR. qCF_1 > qcfMax)

                  ! qcf solution ouside known bounds
                  qCF    = 0.5*(qcfMin+qcfMax)
                  qCF_3  = qCF
                  qCF_2  = qCF
                  qCF_1  = qCF
                  qcfMax = qcfmax0
                  qcfMin = 0.0
                  qCF_1  = qCF
                  qCF_2  = qCF
                  qCF_3  = qCF
                  qcl2_1 = qcl
                  qcl2_2 = qcl
                  qcl2_3 = qcl

                ELSEWHERE (qcfmax <= qcfMin)

                  ! Computational precision issue or convergence
                  qcfMax = qcfmax0
                  qcfMin = 0.0
                  qCF_3  = qCF
                  qCF_2  = qCF
                  qCF_1  = qCF

                ELSEWHERE (qCF < qcfMin)

                  qCFmax=MIN(qCFmax,qCF_1*(1.0+Accuracy))
                  WHERE(qCF_3 < qCF_2 .AND. qCF_2 == qCF_1 .AND.    &
                      ABS(qCF_1 - qCF_3) < accuracy * (qCF_1+qCF_3))

                    ! Need greater accuracy in order to converge
                    accuracy = accuracy*0.1
                    qCF      = 0.5*(qCF_3+qCF_1)
                    qCF_3    = qCF
                    qCF_2    = qCF
                    qCF_1    = qCF

                  ELSEWHERE(qCF_3 < qCF_2 .AND. qCF_2 == qCF_1)

                    ! Reset qCF to know mid-point
                    qCF   = 0.5*(qCF_3+qCF_1)
                    qCF_3 = qCF
                    qCF_2 = qCF
                    qCF_1 = qCF

                  ELSEWHERE

                    ! Reset qCF to know mid-point
                    qCF   = 0.5*(qcfMin+qCF_1)
                    qCF_3 = qCF_1
                    qCF_2 = qCF
                    qCF_1 = qCF

                  ENDWHERE

                ELSEWHERE (qCF > qcfMax)

                  qCFmin=MAX(qcfMin,qCF_1/(1.0+Accuracy))
                  WHERE(qCF_3 > qCF_2 .AND. qCF_2 == qCF_1 .AND.    &
                      ABS(qCF_1 - qCF_3) < accuracy * (qCF_1+qCF_3))

                    ! Need greater accuracy in order to converge
                    accuracy = accuracy*0.1
                    qCF      = 0.5*(qCF_3+qCF_1)
                    qCF_3    = qCF
                    qCF_2    = qCF
                    qCF_1    = qCF

                  ELSEWHERE(qCF_3 > qCF_2 .AND. qCF_2 == qCF_1)

                    ! Reset qCF to know mid-point
                    qCF=0.5*(qCF_3+qCF_1)
                    qCF_3 = qCF
                    qCF_2 = qCF
                    qCF_1 = qCF

                  ELSEWHERE

                    ! Reset qCF to know mid-point
                    qCF=0.5*(qcfMax+qCF_1)
                    qCF_3 = qCF_1
                    qCF_2 = qCF
                    qCF_1 = qCF

                  ENDWHERE

                ELSEWHERE

                  ! Set up succession of previous qCF values
                  qCF_3 = qCF_2
                  qCF_2 = qCF_1
                  qCF_1 = qCF

                ENDWHERE
                ! Set up succession of previous qCL_calc values
                qCL2_3 = qCL2_2
                qCL2_2 = qCL2_1
                qCL2_1 = qCL_calc

              ENDWHERE

            ENDWHERE

          ENDWHERE

          ! Assign qCF increment to dqcf
          dqcf = qCF - qCFin

        ENDWHERE

        IF (loop > maxLoops - 50) THEN

          WHERE(.NOT. converged)

            ! Bail out and just do liquid incrementing if necessary
            qcf   = qcfin
            test  = HUGE(A)
            dqcf  = 0.0

          ENDWHERE

          IF (statF == 0) THEN

            statF = COUNT(.NOT. converged)

          END IF

        END IF
      END IF
    END IF

  END DO

  IF (loop == maxLoops .OR. statF /= 0) THEN

    statL    = COUNT(.NOT. converged)
    ICode    = -1

    IF (PRESENT(TL))THEN
      IF (COUNT(TL < 183.10 .AND. Tmod_temp > 183.10 .AND. &
          .NOT. converged) == COUNT (.NOT. converged)) THEN
        WRITE(CMessage,'(A,I8,A)') 'Gen_VarDiagCloud: all ', &
            COUNT(.NOT.(converged)),                         &
            ' unconverged points at first call due to T~183.15K'
      ELSE
        WRITE(CMessage, '(A,I3,A,I8,A,I8,A,I8)') 'Gen_VarDiagCloud:', loop,  &
        ' loops unconverged, liquid=', statL, ' ice=', statF, ' / ', field_size
      END IF
    ELSE
      WRITE(CMessage, '(A,I3,A,I8,A,I8,A,I8)') 'Gen_VarDiagCloud:', loop,  &
      ' loops unconverged, liquid=', statL, ' ice=', statF, ' / ', field_size
    END IF

  ENDIF

  !----------------------------------------------------------------------
  ! [3]: Assign output values
  !----------------------------------------------------------------------


  IF (PRESENT(TL)) THEN

    TL = Tmod

  END IF

  IF (PRESENT(qCL)) THEN

    qCL = qCL_calc

  END IF

  ! Calculate qcf, given calculated qcl and inputs of qcl & qcf
  IF (assignPseudoBG .AND. .NOT. assign_qCFmax) THEN

    ALLOCATE ( f (field_size)  )
    ALLOCATE ( M (field_size)  )
    ALLOCATE ( fM (field_size) )

    WHERE (BGqcl /= 0.0 .OR. BGqcf /= 0.0)

      f   = BGqcf / (BGqcf + BGqcl)

    ELSEWHERE

      f   = 1.0

    ENDWHERE

    WHERE (BGqcf == 0.0)

      f = 0.0
      M = 0.0

    ELSEWHERE

      M   = (BGqcf / (BGqcf + qCL)) / f

    ENDWHERE

    WHERE (BGqcl == 0.0 .AND. BGqcf == 0.0)

      f = 0.0
      M = 0.0

    ELSEWHERE (BGqcl == 0.0)

      f = 1.0
      ! Where qCL is too small to make M /= 1.0 or machine
      !   precision results in unphysical values
      WHERE (f / (BGqcf / (BGqcf + qCL)) <= 1.0 .OR.              &
             (BGqcf / (BGqcf + qCL)) / f >= 1.0 .OR.              &
             (BGqcf / (BGqcf + qCL)) - f >= 0.0 .OR.              &
             f - (BGqcf / (BGqcf + qCL)) <= 0.0)

        M = 1.0

      ENDWHERE

    ENDWHERE

    fM  = MIN(MAX(f * M, 0.0), 1.0)

    WHERE ((1.0 - fM) * (1.0 - fM) == 0.0 .OR. qCL == 0.0)

      qCF = BGqcf

    ELSEWHERE

      qCF = qCL * (fM / (1.0 - fM))

    ENDWHERE

  ELSE IF (.NOT. AssignPseudoBG .AND. IncrementIce   .AND.        &
                 PRESENT(BGqcl) .AND. PRESENT(BGqcf) .AND.        &
                 PRESENT(qCL)   .AND. PRESENT(qCF)) THEN
    ! NULL

  END IF

  ! Calculate Cl according to scheme in use
  IF (PRESENT(Cl) .AND. .NOT. IncrementIce) THEN

    QN       = Qc / (2.0 * delta)
    QN       = MIN(MAX(QN,                                        &
                       -ALOG(HUGE(QN(1)))*0.99),                  &
                   ALOG(HUGE(QN(1)))*0.99)
    Cl = ( 1.0 + (exp(QN) - exp(-QN)) /                           &
                 (exp(QN) + exp(-QN))  ) / 2.0

    IF (PRESENT(qCL)) THEN

      WHERE(qCL == 0.0)

        Cl = 0.0

      ENDWHERE
      WHERE(Cl == 0.0)

        qCL = 0.0

      ENDWHERE

    END IF

  ELSE IF (PRESENT(Cl) .AND. IncrementIce) THEN

    Cl = 1.0 - EXP(-qCL / delta)

    IF (PRESENT(qCL)) THEN

      WHERE(qCL == 0.0)

        Cl = 0.0

      ENDWHERE
      WHERE(Cl == 0.0)

        qCL = 0.0

      ENDWHERE

    END IF

  END IF

  ! Calculate Cf
  IF (PRESENT(Cf) .AND. IncrementIce) THEN

! DEPENDS ON: qsat_wat
    CALL QSAT_WAT (qSatW_plus,     T,                             &
                    p_theta_levels, field_size)

! When copying this subroutine into GenMod_Utilities,
! replace this call with the following:
!
!    CALL Var_LTinterp(                                            &
!        RESHAPE(T,(/field_size,1,1/)), ! in                       &
!        qSat_temp,                     ! inout                    &
!        WaterOnly=.TRUE.) ! in
!    ! Note that qSatW_plus is being reused to refer to qSatW
!    qSatW_plus = RESHAPE(qSat_temp,(/field_size/))
!    qSatW_plus = Epsilon * qSatW_plus / p_theta_levels


    Cf = 1.0 - EXP(-4.0*qCF/((1.0-RHc)*qSatW_plus))

    IF (PRESENT(qCF)) THEN

      WHERE(qCF == 0.0)

        Cf = 0.0

      ENDWHERE
      WHERE(Cf == 0.0)

        qCF = 0.0

      ENDWHERE

    END IF

  END IF
! When copying this subroutine into GenMod_Utilities,
! uncomment this section; to make linearisation code available to Var:

!   IF ((PRESENT(dCLdp)    .OR.                                     &
!        PRESENT(dCLdqT)   .OR.                                     &
!        PRESENT(dCLdT)    .OR.                                     &
!        PRESENT(dqCLdp)   .OR.                                     &
!        PRESENT(dqCLdqT)  .OR.                                     &
!        PRESENT(dqCLdT) ) .AND.                                    &
!        PRESENT(T)        .AND. .NOT. IncrementIce    ) THEN
!
!     ! Liquid-only scheme linearisation
!
!     ALLOCATE ( B (field_size) )
!     ALLOCATE ( C (field_size) )
!     ALLOCATE ( D (field_size) )
!     ALLOCATE ( E (field_size) )
!     ALLOCATE ( F (field_size) )
!     ALLOCATE ( G (field_size) )
!     ALLOCATE ( H (field_size) )
!     ALLOCATE ( I (field_size) )
!     ALLOCATE ( J (field_size) )
!     ALLOCATE ( K (field_size) )
!     ALLOCATE ( L (field_size) )
!     ALLOCATE ( M (field_size) )
!     ALLOCATE ( N (field_size) )
!     ALLOCATE ( O (field_size) )
!     ALLOCATE ( multiplier (field_size) )
!
!     A  = -Lc / Cp
!     QN = Qc / delta
!     QN = MAX(QN, -ALOG(HUGE(QN(1)))*0.99)
!     B  = 1.0 - EXP(-QN) / (EXP(-QN) + 1.0)
!     C  = ALOG(EXP(-QN) + 1.0) +                                   &
!          Qc * EXP(-QN) / (delta * (EXP(-QN) + 1.0))
!     D  = -aL
!
!     CALL Var_LTinterp(                                            &
!         RESHAPE(Tmod+0.05,(/field_size,1,1/)),                    & ! in
!         qSat_temp,                                                & ! inout
!         p_temp,                                                   & ! inout
!         .TRUE.,                                                   & ! in
!         WaterOnly=.TRUE.) ! in
!     qSatW_plus = RESHAPE(qSat_temp,(/field_size/))
!
!     CALL Var_LTinterp(                                            &
!         RESHAPE(Tmod-0.05,(/field_size,1,1/)),                    & ! in
!         qSat_temp,                                                & ! inout
!         p_temp,                                                   & ! inout
!         .TRUE.,                                                   & ! in
!         WaterOnly=.TRUE.) ! in
!     qSatW_minus = RESHAPE(qSat_temp,(/field_size/))
!
!     E  = Epsilon * ((qSatW_plus - qSatW_minus) / 0.1) /           &
!          p_theta_levels
!     F  = 0.25 * aL * (1.0 - RHc)
!     QN       = Qc / (2.0 * delta)
!     QN       = MIN(MAX(QN,                                        &
!                        -ALOG(HUGE(QN(1)))*0.99),                  &
!                    ALOG(HUGE(QN(1)))*0.99)
!     G = (1.0 - tanh(QN) * tanh(QN)) / (4.0 * delta)
!     H  = -Qc * G / delta
!     I  = -D
!     J  = qT - qSatW
!     K  = -qsatW / p_theta_levels
!     L  = qSatW * (1.0 - RHc) / 4.0
!     M  = - Lc / (Cp *                                             &
!          (1.0 + Lc * alphaL / Cp) * (1.0 + Lc * alphaL / Cp))
!     N  = Epsilon * ((dPsatdT_plus - dPsatdT_minus) / 0.1)         &
!          / p_theta_levels
!     O  = -alphaL / p_theta_levels
!
!     multiplier = 1.0 / (1.0 - A * E * (B * D + C * F))
!
!     IF (PRESENT(dCLdp)) THEN
!
!       dCLdp = multiplier * G * (K * (D - F * QN) + M * O *        &
!           (A * E * (C + B * QN) * (D * L - F * J) - L * QN + J))
!
!     END IF
!
!     IF (PRESENT(dCLdqT)) THEN
!
!       dCLdqT = multiplier * G * I *                               &
!           (1.0 - A * E * F * (C + B * QN))
!
!     END IF
!
!     IF (PRESENT(dCLdT)) THEN
!
!       dCldT = multiplier * G * (E * (D - F * QN) + M * N *        &
!           (A * E * (C + B * QN) * (D * L - F * J) - L * QN + J))
!
!     END IF
!
!     IF (PRESENT(dqCLdp)) THEN
!
!       dqCLdp = multiplier * (K * (B * D + C * F) + M * O *      &
!           (B * J + C * L))
!     END IF
!
!     IF (PRESENT(dqCLdqT)) THEN
!
!       dqCLdqT = multiplier * B * I
!
!     END IF
!
!     IF (PRESENT(dqCLdT)) THEN
!
!       dqCLdT = multiplier * (E * (B * D + C * F) + M * N *      &
!           (B * J + C * L))
!
!     END IF
!
!   ELSE IF ((PRESENT(dCLdp)    .OR.                                &
!             PRESENT(dCLdqT)   .OR.                                &
!             PRESENT(dCLdT)    .OR.                                &
!             PRESENT(dqCLdp)   .OR.                                &
!             PRESENT(dqCLdqT)  .OR.                                &
!             PRESENT(dqCLdT)   .OR.                                &
!             PRESENT(dCFdp)    .OR.                                &
!             PRESENT(dCFdqcf)  .OR.                                &
!             PRESENT(dCFdT)    .OR.                                &
!             PRESENT(dqCFdqCL) .OR.                                &
!             PRESENT(dqCFdT) ) .AND.                               &
!             PRESENT(T)        .AND. IncrementIce) THEN
!
!     ! Liquid & ice scheme linearisation
!
!     ALLOCATE ( B (field_size) )
!     ALLOCATE ( C (field_size) )
!     ALLOCATE ( D (field_size) )
!     ALLOCATE ( E (field_size) )
!     ALLOCATE ( G (field_size) )
!     ALLOCATE ( GCx (field_size) )
!     ALLOCATE ( H (field_size) )
!     ALLOCATE ( I (field_size) )
!     ALLOCATE ( J (field_size) )
!     ALLOCATE ( K (field_size) )
!     ALLOCATE ( L (field_size) )
!     ALLOCATE ( N (field_size) )
!     ALLOCATE ( O (field_size) )
!     ALLOCATE ( V (field_size) )
!     ALLOCATE ( W (field_size) )
!     ALLOCATE ( Y (field_size) )
!     ALLOCATE ( Z (field_size) )
!     ALLOCATE ( multiplier (field_size) )
!
!     A  = -Lc / Cp
!     QN = Qc / delta
!     QN = MAX(QN, -ALOG(HUGE(QN(1)))*0.99)
!     B  = 1.0 - EXP(-QN) / (EXP(-QN) + 1.0)
!     C  = ALOG(EXP(-QN) + 1.0) +                                   &
!          Qc * EXP(-QN) / (delta * (EXP(-QN) + 1.0))
!     D  = -aL
!
!     CALL Var_LTinterp(                                            &
!         RESHAPE(Tmod+0.05,(/field_size,1,1/)),                    & ! in
!         qSat_temp,                                                & ! inout
!         p_temp,                                                   & ! inout
!         .TRUE.,                                                   & ! in
!         WaterOnly=.TRUE.) ! in
!     qSatW_plus = RESHAPE(qSat_temp,(/field_size/))
!
!     CALL Var_LTinterp(                                            &
!         RESHAPE(Tmod-0.05,(/field_size,1,1/)),                    & ! in
!         qSat_temp,                                                & ! inout
!         p_temp,                                                   & ! inout
!         .TRUE.,                                                   & ! in
!         WaterOnly=.TRUE.) ! in
!     qSatW_minus = RESHAPE(qSat_temp,(/field_size/))
!
!     E = aL
!     G = aL * (1.0 - RHc) * 0.25
!     H = Epsilon * ((qSatW_plus - qSatW_minus) / 0.1) /            &
!          p_theta_levels
!     I = -1.0
!     K = qT - qSatW
!     L = qSatW * (1.0 - RHc) * 0.25
!     N = -qsatW / p_theta_levels
!     V = - Lc / (Cp *                                              &
!          (1.0 + Lc * alphaL / Cp) * (1.0 + Lc * alphaL / Cp))
!     W = -alphaL / p_theta_levels
!     Y = Epsilon * ((dPsatdT_plus - dPsatdT_minus) / 0.1)          &
!          / p_theta_levels
!
!     WHERE (BGqcf == 0.0 .OR.                                      &
!            ((1.0 - fM) * (1.0 - fM) == 0.0 .OR. qCL == 0.0))
!
!       ! Equivalent to BGqcf == 0.0 rule for f & M and
!       ! linearisation of qcf calculation not available,
!       ! since there is no gradient in the qcl scheme, or f ~ 1
!       GCx = 0.0
!       O   = 0.0
!       J   = 0.0
!
!     ELSEWHERE
!
!       J = fM / (1.0 - fM)
!       O = qCL * M  /((1.0 - fM) * (1.0 - fM))
!
!       WHERE (M == 0.0)
!         GCx = MIN((1.0 - (BGqcf / (BGqcf + BGqcl))),              &
!                 (BGqcf / (BGqcf + BGqcl)))
!       ELSEWHERE
!         GCx = MIN((MIN(1.0 / M, 1.0) - (BGqcf / (BGqcf + BGqcl))),&
!                 (BGqcf / (BGqcf + BGqcl)))
!       ENDWHERE
!
!     ENDWHERE
!
!     Multiplier = 1.0 / (1.0 - A * H * (B * D + C * G) -           &
!                         B * E * I * J)
!     Z = GCx * (-1.0 / (2.0 * cosh(0.18 * T - 260.65) *            &
!                              cosh(0.18 * T - 260.65)))
!
!     IF (PRESENT(dqCLdp)) THEN
!
!       dqCLdp = Multiplier * ((B * K + C * L) * W * V +            &
!                              (B * D + C * G) * N)
!
!     END IF
!
!     IF (PRESENT(dqCLdqT)) THEN
!
!       dqCLdqT = Multiplier * B * E
!
!     END IF
!
!     IF (PRESENT(dqCLdT)) THEN
!
!       dqCLdT = Multiplier * ((B * D + C * G) * H +                &
!                              (B * K + C * L) * Y * V +            &
!                              E * I * B * O * Z)
!
!     END IF
!
!     IF (PRESENT(dqCFdqCL)) THEN
!
!       dqCFdqCL = J
!
!     END IF
!
!     IF (PRESENT(dqCFdT)) THEN
!
!       dqCFdT = O*Z
!
!     END IF
!
!     ! Reuse qCL_1 & qCL_2 as a temporary variable
!     qCL_1 = EXP(-qCL / delta) / delta
!     qCL_2 = -qCL * EXP(-qCL / delta) / (delta * delta)
!
!     IF (PRESENT(dCLdqT)) THEN
!
!       IF (PRESENT(dqCLdqT)) THEN
!
!         dCLdqT = dqCLdqT * qCL_1 +                                &
!             (Multiplier * (G * E * H * A * B)) * qCL_2
!
!       ELSE
!
!         ICode    = 1
!         CMessage =                                                &
!         'Error from Gen_VarDiagCloud: add dqCLdqT to arg list in CALL'
!
!       END IF
!
!     END IF
!
!     IF (PRESENT(dCLdT)) THEN
!
!       IF (PRESENT(dqCLdT)) THEN
!
!         dCLdT = dqCLdT * qCL_1 +                                  &
!             (Multiplier *                                         &
!              (Y * V * (L + B * (A * H * (G * K - D * L) -         &
!                                 J * E * L * I))           +       &
!              H * G * (1.0 + B * E * I * (O * Z * A - J)))) * qCL_2
!
!       ELSE
!
!         ICode    = 1
!         CMessage =                                                &
!         'Error from Gen_VarDiagCloud: add dqCLdT to arg list in CALL'
!
!       END IF
!
!     END IF
!
!     IF (PRESENT(dCLdp)) THEN
!
!       IF (PRESENT(dqCLdp)) THEN
!
!         dCLdp = dqCLdp * qCL_1 +                                  &
!             (Multiplier *                                         &
!              (W * V * (L + B * (A * H * (G * K - D * L) -         &
!               J * E * L * I)))) * qCL_2
!
!       ELSE
!
!         ICode    = 1
!         CMessage =                                                &
!         'Error from Gen_VarDiagCloud: add dqCLdp to arg list in CALL'
!
!       END IF
!
!     END IF
!
!     CALL Var_LTinterp(                                            &
!         RESHAPE(T,(/field_size,1,1/)),                            & ! in
!         qSat_temp,                                                & ! inout
!         p_temp,                                                   & ! inout
!         .TRUE.,                                                   & ! in
!         WaterOnly=.TRUE.) ! in
!     qSatW_plus = Epsilon * RESHAPE(qSat_temp,(/field_size/)) /    &
!                  p_theta_levels
!     dPsatdT_plus = RESHAPE(p_temp,(/field_size/))
!
!
!     IF (PRESENT(dCFdqCF)) THEN
!       dCFdqCF = -4.0 * EXP((-4.0 * qCF) /                         &
!                            ((1.0 - RHc) * qSatW_plus))            &
!                 / ((1.0 - RHc) * qSatW_plus)
!
!     END IF
!
!     ! Reuse dqcf as a temporary variable
!     dqcf = -4.0 * qCF * EXP((-4.0 * qCF) /                        &
!                             ((1.0 - RHc) * qSatW_plus))           &
!            / ((1.0 - RHc) * qSatW_plus * qSatW_plus)
!
!     IF (PRESENT(dCFdT)) THEN
!
!       dCFdT = dqcf * (Epsilon * dPsatdT_plus / p_theta_levels)
!
!     END IF
!
!     IF (PRESENT(dCFdp)) THEN
!
!       dCFdp = - dqcf * (qSatW_plus / p_theta_levels)
!
!     END IF
!
!   END IF

!----------------------------------------------------------------------
! [4]: Deallocate arrays
!----------------------------------------------------------------------

  IF  (ALLOCATED(alphaL)       ) DEALLOCATE ( alphaL        )
  IF  (ALLOCATED(dPsatdT_minus)) DEALLOCATE ( dPsatdT_minus )
  IF  (ALLOCATED(dPsatdT_plus) ) DEALLOCATE ( dPsatdT_plus  )
  IF  (ALLOCATED(fM)           ) DEALLOCATE ( fM            )
  IF  (ALLOCATED(fPlus)        ) DEALLOCATE ( fPlus         )
  IF  (ALLOCATED(fplusM)       ) DEALLOCATE ( fPlusM        )
  IF  (ALLOCATED(f_of_T)       ) DEALLOCATE ( f_of_T        )
  IF  (ALLOCATED(f_of_T_Plus)  ) DEALLOCATE ( f_of_T_plus   )
  IF  (ALLOCATED(multiplier)   ) DEALLOCATE ( multiplier    )
  IF  (ALLOCATED(qCL_4)        ) DEALLOCATE ( qCL_4         )
  IF  (ALLOCATED(qCFmax0)      ) DEALLOCATE ( qCFmax0       )
  IF  (ALLOCATED(qCFmin)       ) DEALLOCATE ( qCFmin        )
  IF  (ALLOCATED(qCFin)        ) DEALLOCATE ( qCFin         )
  IF  (ALLOCATED(qCF_1)        ) DEALLOCATE ( qCF_1         )
  IF  (ALLOCATED(qCF_2)        ) DEALLOCATE ( qCF_2         )
  IF  (ALLOCATED(qCF_3)        ) DEALLOCATE ( qCF_3         )
  IF  (ALLOCATED(qCL2_1)       ) DEALLOCATE ( qCL2_1        )
  IF  (ALLOCATED(qCL2_2)       ) DEALLOCATE ( qCL2_2        )
  IF  (ALLOCATED(qCL2_3)       ) DEALLOCATE ( qCL2_3        )
  IF  (ALLOCATED(B)            ) DEALLOCATE ( B             )
  IF  (ALLOCATED(C)            ) DEALLOCATE ( C             )
  IF  (ALLOCATED(D)            ) DEALLOCATE ( D             )
  IF  (ALLOCATED(E)            ) DEALLOCATE ( E             )
  IF  (ALLOCATED(F)            ) DEALLOCATE ( F             )
  IF  (ALLOCATED(G)            ) DEALLOCATE ( G             )
  IF  (ALLOCATED(GCx)          ) DEALLOCATE ( GCx           )
  IF  (ALLOCATED(H)            ) DEALLOCATE ( H             )
  IF  (ALLOCATED(I)            ) DEALLOCATE ( I             )
  IF  (ALLOCATED(J)            ) DEALLOCATE ( J             )
  IF  (ALLOCATED(K)            ) DEALLOCATE ( K             )
  IF  (ALLOCATED(L)            ) DEALLOCATE ( L             )
  IF  (ALLOCATED(M)            ) DEALLOCATE ( M             )
  IF  (ALLOCATED(N)            ) DEALLOCATE ( N             )
  IF  (ALLOCATED(O)            ) DEALLOCATE ( O             )
  IF  (ALLOCATED(V)            ) DEALLOCATE ( V             )
  IF  (ALLOCATED(W)            ) DEALLOCATE ( W             )
  IF  (ALLOCATED(Y)            ) DEALLOCATE ( Y             )
  IF  (ALLOCATED(Z)            ) DEALLOCATE ( Z             )

ENDIF

RETURN
END SUBROUTINE Gen_VarDiagCloud

#endif
