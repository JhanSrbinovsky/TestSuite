#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
SUBROUTINE Fill_Pressure( Inf1,         &  !in
                          Infm,         &  !in
                          Outf,         &  !in
                          Pfields,      &  !inout
                          ErrorStatus)     !inout

! Description:
!   Interpolate field values within Cb mask area.
!
! Method:
! Performs linear interpolation horizontally E-W and vertically N-S 
! to fill in the field values (pressure or heights of Cb bases and tops
! or horizontal extent of Cb) at points within the Cb mask that have 
! missing values. 
!
! Owner: Dave Jerrett
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6
!
!
!
USE IO_mod, ONLY:        &
    PP_Header_type,      &
    PP_Field_type
USE Err_Mod, ONLY:       &
    StatusOK

IMPLICIT None


INTEGER, INTENT(IN) :: Inf1                        ! Input/output field index
INTEGER, INTENT(IN) :: Infm                        ! Mask field index
INTEGER, INTENT(IN) :: Outf                        ! Dummy field index (can
                                                   !  also be used for output)
TYPE(PP_Field_type), INTENT(INOUT) :: Pfields(270) ! Fields array
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Fill_Pressure"
#include "c_mdi.h"

INTEGER :: i, j, k                                  ! loop counters
REAL    :: a                                        ! holds adjusted field value  

!  Variables p & q contain the input field values in a 3x3 box around the point
!  of interest (p5/q5)
!
!  p1,p2,p3,
!  p4,p5,p6,
!  p7,p8,p9,
!
!  q1,q2,q3,
!  q4,q5,q6,
!  q7,q8,q9
!
! Only p5, q3, q5, q6, q9 are used in this routine
!
REAL    :: p5, q3, q5, q6, q9


! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF



! Copy input field to dummy field ready for checking the field and making 
! adjustments.
Pfields(Outf)%RData(:,:)=Pfields(Inf1)%RData(:,:)


!---- Check through the field and adjust according to neighbouring points ----
!
! If input field has a missing value at point of interest but mask field 
! indicates Cb, check if there is an input field value at the southerly, 
! northery or easterly neighbouring gridpoints. If so, set value at point 
! of interest to be this value. Repeat check through field a further 9 times.   
!
! Optimisation comments:
! kk Exchanging I and J loop avoids bank conflicts AND allows vectorization.

DO k = 1, 10
  DO j = 1, PFields(1) % Hdr % NumRows
    DO i = 1, PFields(1) % Hdr % NumCols

      p5 = Pfields(Infm) % RData(i,j)
      q3 = Pfields(Outf) % RData(i+1,j+1)
      q5 = Pfields(Outf) % RData(i,j)
      q6 = Pfields(Outf) % RData(i+1,j)
      q9 = Pfields(Outf) % RData(i+1,j-1)

      IF ( (q5 == RMDI ) .AND. ( p5 /= RMDI)) THEN
        a=RMDI
        IF (q9 /= RMDI )  THEN
          a=q9
        END IF
        IF (q3 /= RMDI ) THEN
          a=q3
        END IF
        IF (q6 /= RMDI ) THEN
          a=q6
        END IF
        Pfields(Outf)%RData(i,j)=a
      END IF
    END DO
  END DO
END DO


! Copy adjusted field back to Pfields(Inf1) ready for output 
Pfields(Inf1) % Rdata(:,:) = Pfields(Outf) % RData(:,:)


9999 CONTINUE

END SUBROUTINE Fill_Pressure
#endif
