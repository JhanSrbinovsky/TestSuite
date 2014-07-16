#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculate basic statistics for a horizontal field.

      SUBROUTINE FieldStats (                                           &
#include "argcona.h"
#include "arglndm.h"
     &                        Field_len,                                &
                                                  ! in
     &                        Field,                                    &
                                                  ! in
     &                        grid_type,                                &
                                                  ! in
     &                        halo_type,                                &
                                                  ! in
     &                        Global_max,                               &
                                                  ! out
     &                        Global_min,                               &
                                                  ! out
     &                        Global_mean,                              &
                                                  ! out
     &                        Global_RMS )        ! out

      Use trignometric_mod, Only : cos_theta_latitude, cos_v_latitude

! Description:
!
!   Calculate basic statistics for a horizontal field on u, v, theta
!   or land points.
!
!
! Current Code Owner: Adam Clayton.
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 5.2       14/12/00 Original code. Adam Clayton
! 6.2       21/10/05 Replace GSYNC with SSYNC. P.Selwood
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! Declarations:

      IMPLICIT NONE

! Common blocks:

#include "cppxref.h"
#include "parvars.h"
#include "typsize.h"
#include "cmaxsize.h"
#include "typcona.h"
#include "typlndm.h"

! Subroutine arguments:

      INTEGER, INTENT(IN)  :: Field_len
      REAL,    INTENT(IN)  :: Field(Field_len) ! Horizontal field
      INTEGER, INTENT(IN)  :: grid_type
      INTEGER, INTENT(IN)  :: halo_type
      REAL,    INTENT(OUT) :: Global_max
      REAL,    INTENT(OUT) :: Global_min
      REAL,    INTENT(OUT) :: Global_mean      ! Area-weighted mean
      REAL,    INTENT(OUT) :: Global_RMS       ! Area-weighted RMS

! Local parameters:

      CHARACTER(*)  :: RoutineName
      PARAMETER      ( RoutineName='FieldStats' )

! Local variables:

      INTEGER :: i, j, k, p
      INTEGER :: ICode
      INTEGER :: xm, ym
      INTEGER :: xh, yh
      INTEGER :: Expected_len
      INTEGER :: Addr

      LOGICAL :: ThetaRows
      LOGICAL :: VPoints
      LOGICAL :: LandPoints

      REAL :: Value
      REAL :: Weight

      REAL :: Local_max
      REAL :: Local_min
      REAL :: Local_sum
      REAL :: Local_sumsq
      REAL :: Local_sumwts

      REAL :: Global_sum
      REAL :: Global_sumsq
      REAL :: Global_sumwts

      REAL :: ReshapedField(row_length,rows)
      REAL :: WeightedField(row_length,rows)
      REAL :: FieldWeights (row_length,rows)

      REAL :: LocalStats(nproc,5)

      CHARACTER(80) :: CMessage

! External subroutines called:

      EXTERNAL EReport

!- End of header ------------------------------------------------------

!----------------------------------------------------------------------
! [1]: Do some checks based on grid and halo type.
!----------------------------------------------------------------------

      ThetaRows  = .FALSE.
      VPoints    = .FALSE.
      LandPoints = .FALSE.

      IF (grid_type == ppx_atm_tall .OR.                                &
     &    grid_type == ppx_atm_cuall) THEN
        ThetaRows = .TRUE.
        xm = row_length
        ym = rows
      ELSE IF (grid_type == ppx_atm_cvall) THEN
        VPoints = .TRUE.
        xm = row_length
        ym = n_rows
      ELSE IF (grid_type == ppx_atm_compressed) THEN
        LandPoints = .TRUE.
        Expected_len = land_field
      ELSE
        ICode = 1
        CMessage = 'Grid type not catered for.'
! DEPENDS ON: ereport
        CALL EReport ( RoutineName, ICode, CMessage )
      END IF

      IF (.NOT.LandPoints) THEN

        IF (halo_type == halo_type_single) THEN
          xh = offx
          yh = offy
        ELSE IF (halo_type == halo_type_extended) THEN
          xh = halo_i
          yh = halo_j
        ELSE IF (halo_type == halo_type_no_halo) THEN
          xh = 0
          yh = 0
        ELSE
          ICode = 1
          CMessage = 'Invalid halo type.'
! DEPENDS ON: ereport
          CALL EReport ( RoutineName, ICode, CMessage )
        END IF

        Expected_len = (xm + 2*xh) * (ym + 2*yh)

      END IF

      IF (Field_len /= Expected_len) THEN
        ICode = 1
        CMessage = 'Supplied field length incorrect.'
! DEPENDS ON: ereport
        CALL EReport ( RoutineName, ICode, CMessage )
      END IF

!----------------------------------------------------------------------
! [2]: Get local statistics for fields on u, v, or theta points.
!----------------------------------------------------------------------

      IF (.NOT.LandPoints) THEN

        ReshapedField(:,:) = 0.0
        WeightedField(:,:) = 0.0
        FieldWeights (:,:) = 0.0

        Addr = 1

        DO j = 1 - yh, ym + yh
          DO i = 1 - xh, xm + xh
            IF (i >=1 .AND. i <= xm .AND.                               &
     &          j >=1 .AND. j <= ym) THEN
              ReshapedField(i,j) = Field(Addr)
            END IF
            Addr = Addr + 1
          END DO
        END DO

        IF (ThetaRows) THEN
          FieldWeights(1:xm,1:ym) = cos_theta_latitude(1:xm,1:ym)
        ELSE
          FieldWeights(1:xm,1:ym) = cos_v_latitude    (1:xm,1:ym)
        END IF

        WeightedField(1:xm,1:ym) = ReshapedField(1:xm,1:ym)             &
     &                           * FieldWeights (1:xm,1:ym)

        Local_max    = MAXVAL( ReshapedField(:,:) )
        Local_min    = MINVAL( ReshapedField(:,:) )
        Local_sum    = SUM   ( WeightedField(:,:) )
        Local_sumsq  = SUM   ( ReshapedField(:,:)                       &
     &                       * WeightedField(:,:) )
        Local_sumwts = SUM   ( FieldWeights (:,:) )

      END IF

!----------------------------------------------------------------------
! [3]: Get local statistics for fields on land points.
!----------------------------------------------------------------------

      IF (LandPoints) THEN

        Local_max    = -HUGE(1.0)
        Local_min    =  HUGE(1.0)
        Local_sum    = 0.0
        Local_sumsq  = 0.0
        Local_sumwts = 0.0

        DO p = 1, land_field
          i = 1 + MOD(land_index(p)-1, row_length)
          j = 1 + (land_index(p)-1) / row_length
          Value        = Field(p)
          Weight       = cos_theta_latitude(i,j)
          Local_max    = MAX(Local_max, Value)
          Local_min    = MIN(Local_min, Value)
          Local_sum    = Local_sum    + Weight * Value
          Local_sumsq  = Local_sumsq  + Weight * Value**2
          Local_sumwts = Local_sumwts + Weight
        END DO

      END IF

!----------------------------------------------------------------------
! [4]: Get global statistics.
!----------------------------------------------------------------------

      LocalStats(:,:) = 0.0

      LocalStats(mype+1,:) = (/ Local_max,                              &
     &                          Local_min,                              &
     &                          Local_sum,                              &
     &                          Local_sumsq,                            &
     &                          Local_sumwts /)

      CALL GC_SSYNC ( nproc, ICode )

      CALL GCG_RSUMR ( nproc*5,                                         &
                                          ! in
     &                 gc_all_proc_group,                               &
                                          ! in
     &                 ICode,                                           &
                                          ! out
     &                 LocalStats )       ! inout

      Global_max    = MAXVAL(LocalStats(:,1))
      Global_min    = MINVAL(LocalStats(:,2))
      Global_sum    = SUM   (LocalStats(:,3))
      Global_sumsq  = SUM   (LocalStats(:,4))
      Global_sumwts = SUM   (LocalStats(:,5))

      Global_mean =      Global_sum   / Global_sumwts
      Global_RMS  = SQRT(Global_sumsq / Global_sumwts)


      RETURN
      END SUBROUTINE FieldStats
#endif
