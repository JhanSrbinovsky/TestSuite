#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine used in Fieldcalc CB actions

!=======================================================================
SUBROUTINE ConvAct( NumLevs,       &  ! in
                    Pfields,       &  ! inout
                    ErrorStatus )     ! inout
!
! Description:
!   Calculates CB Action fields
!
! Method:
!  Subroutine ConvAct takes fields of convective cloud combined with the 
!  field of convective precipitation rate to produce a mask of Cumulonimbus 
!  activity. This Cb mask is then compared to convective cloud top and base
!  pressures to give Cb cloud top and base pressures. These fields can later 
!  be converted to heights (in Kft) using the ICAO_HT action.
! 
!
!  Contents of Pfields:
!    Pfields(1)= Will contain Cb base pressure
!    Pfields(2)= Will contain Cb top pressure
!    Pfields(3)= In this subroutine, used to contain CB mask 
!                (after this subroutine has finished this will
!                later be used to contain the Cb base height)
!    Pfields(4)= Used as a dummy/spare field space for adjusting fields 
!                (after this subroutine has finished this will
!                later be used to contain Cb top height)
!    Pfields(5)= Will contain embedded Cb base pressure
!    Pfields(6)= Will contain embedded Cb top pressure
!    Pfields(7)= In this subroutine, used to contain embedded CB mask
!                (after this subroutine has finished this 
!                later be used to contain embedded Cb base height)
!    Pfields(8)= Used as a spare field space for adjusting fields
!                (after this subroutine has finished this will 
!                later be used to contain embedded Cb top height)       
!    Pfields(9)= Will contain Cb horizontal extent (as index)
!
!    PFields on Input :
!    Pfields(10) : Convective Cloud Base Pressure (5/207)
!    Pfields(11) : Convective Cloud Top  Pressure (5/208)
!    Pfields(12) : Convective Precipitation Rate  (5/205)
!    Pfields(13) : Lowest Conv Cloud Top Pressure (5/223)
!
!    Pfields(40-)  : Bulk Cloud Fraction     (0/266)
!    Pfields(80-)  : Convective Cloud Amount (5/212)
!    Pfields(120-) : Theta level Pressure    (0/408)
!    Pfields(160-) : Theta level Temperature (16/004)
!
!Debi Turp: I don't think Pfields(10), Pfields(11), Pfields(13) are used 
!in this subroutine so shouldn't be referenced? 
!
!
! Owner: Dave Jerrett
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6
!
 
USE IO_mod, ONLY:        &
    PP_Header_type,      &
    PP_Field_type
USE Err_Mod, ONLY:       &
    StatusOK
USE FldCodes_mod, ONLY:  &
ST_P_CBB,   MO8_P_CBB,  PP_P_CBB,  VC_P_CBB, &
ST_P_CBT,   MO8_P_CBT,  PP_P_CBT,  VC_P_CBT, &
ST_I_CBB,   MO8_I_CBB,  PP_I_CBB,  VC_I_CBB, &
ST_I_CBT,   MO8_I_CBT,  PP_I_CBT,  VC_I_CBT, &
ST_P_ECBB,  MO8_P_ECBB, PP_P_ECBB, VC_P_ECBB,&
ST_P_ECBT,  MO8_P_ECBT, PP_P_ECBT, VC_P_ECBT,&
ST_I_ECBB,  MO8_I_ECBB, PP_I_ECBB, VC_I_ECBB,&
ST_I_ECBT,  MO8_I_ECBT, PP_I_ECBT, VC_I_ECBT,&
ST_CBHorE,  MO8_CBHorE, PP_CBHorE, VC_CBHorE,&
ST_CPNRT,   MO8_CPNRT,  PP_CPNRT,  VC_CPNRT

IMPLICIT None

INTEGER, INTENT(IN)                :: Numlevs         !No. of levels
TYPE(PP_Field_type), INTENT(INOUT) :: Pfields(270)    !Fields array
INTEGER, INTENT(INOUT)             :: Errorstatus


! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "ConvAct"
#include "c_mdi.h"


! Local Variables:
INTEGER :: i, j, k              ! Loop counters
INTEGER :: n                    ! Debi Turp: I don't think this is used?
INTEGER :: pcnt                 ! Holds index when determining Cb horiz extent
INTEGER :: ii, jj               ! Loop counters in optimisation loops
INTEGER :: z                    ! Loop counter for field index 

! Parameters used in connection with the subroutine Fill_Pressure 
INTEGER :: In_field             ! index for input field to Fill_Pressure
INTEGER :: In_mask              ! index for input mask field to Fill_Pressure
INTEGER :: Out_field            ! index for input field to Fill_Pressure
INTEGER :: Nofill               ! flag, set to 1 to call Fill_Pressure to
                                !   fill in values in area of mask

! Variables minnm, minnmb, diffn, diffnb below used in to check cloud tops
! are above cloud bases 
!
REAL :: diffn                   ! difference of CB base pressure
                                !  and CB top pressure
REAL :: diffnb                  ! difference of embedded CB base pressure  
                                !  and embedded CB top pressure
REAL :: minnm                   ! indicates CB base is above CB top if negative
REAL :: minnmb                  ! indicates embedded CB base is above embedded 
                                !  CB top if negative

!  Variables p & q below are used to contain field values in a 3x3 box  
!  around the point of interest (which is p5/q5)
!
!  p1,p2,p3,
!  p4,p5,p6,
!  p7,p8,p9,
!
!  q1,q2,q3,
!  q4,q5,q6,
!  q7,q8,q9
!
! Debi Turp: I don't think p10 and q10 are used - should be deleted?
!
REAL :: p1, p2, p3, p4, p5, p6, p7, p8, p9, p10
REAL :: q1, q2, q3, q4, q5, q6, q7, q8, q9, q10


! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF



!Debi Turp: why copy Pfields(12) into Pfields(10) below - why not use 
!Pfields(12) directly later in the subroutine?

Pfields(10) % LookupPos  = Pfields(12) % LookupPos
Pfields(10) % ArrayPos   = Pfields(12) % ArrayPos
Pfields(10) % Hdr        = Pfields(12) % Hdr
Pfields(10) % RData(:,:) = Pfields(12) % RData(:,:)

Pfields % Hdr % BMDI    = RMDI
Pfields(21) % Hdr % MO8Level           = 8888



!---- Allocate output fields and set fieldcodes etc. in headers ----

DO z = 1,9
  IF(associated(Pfields(z) % Rdata))then
    deallocate(Pfields(z) % Rdata)
  END IF
  Pfields(z) % LookupPos  = Pfields(21) % LookupPos
  Pfields(z) % ArrayPos   = Pfields(21) % ArrayPos
  Pfields(z) % Hdr        = Pfields(21) % Hdr
END DO

Pfields(1) % Hdr % STCode  = ST_P_CBB
Pfields(1) % Hdr % MO8Type = MO8_P_CBB
Pfields(1) % Hdr % PPcode  = PP_P_CBB
Pfields(1) % Hdr % LBVC    = VC_P_CBB
ALLOCATE( Pfields(1) % RData(Pfields(21) % Hdr % NumCols,     &
                             Pfields(21) % Hdr % NumRows) )
Pfields(1) % RData(:,:) = RMDI

Pfields(2) % Hdr % STCode  = ST_P_CBT
Pfields(2) % Hdr % MO8Type = MO8_P_CBT
Pfields(2) % Hdr % PPCode  = PP_P_CBT
Pfields(2) % Hdr % LBVC    = VC_P_CBT
ALLOCATE( Pfields(2) % RData(Pfields(21) % Hdr % NumCols,     &
                             Pfields(21) % Hdr % NumRows) )
Pfields(2) % RData(:,:) = RMDI

Pfields(3) % Hdr % STCode  = ST_I_CBB
Pfields(3) % Hdr % MO8Type = MO8_I_CBB
Pfields(3) % Hdr % PPCode  = PP_I_CBB
Pfields(3) % Hdr % LBVC    = VC_I_CBB
ALLOCATE( Pfields(3) % RData(Pfields(21) % Hdr % NumCols,     &
                             Pfields(21) % Hdr % NumRows) )
Pfields(3) % RData(:,:) = RMDI

Pfields(4) % Hdr % STCode  = ST_I_CBT
Pfields(4) % Hdr % MO8Type = MO8_I_CBT
Pfields(4) % Hdr % PPCode  = PP_I_CBT
Pfields(4) % Hdr % LBVC    = VC_I_CBT
ALLOCATE( Pfields(4) % RData(Pfields(21) % Hdr % NumCols,     &
                             Pfields(21) % Hdr % NumRows) )
Pfields(4) % RData(:,:) = RMDI

Pfields(5) % Hdr % STCode  = ST_P_ECBB
Pfields(5) % Hdr % MO8Type = MO8_P_ECBB
Pfields(5) % Hdr % PPCode  = PP_P_ECBB
Pfields(5) % Hdr % LBVC    = VC_P_ECBB
ALLOCATE( Pfields(5) % RData(Pfields(21) % Hdr % NumCols,     &
                             Pfields(21) % Hdr % NumRows) )
Pfields(5) % RData(:,:) = RMDI

Pfields(6) % Hdr % STCode  = ST_P_ECBT
Pfields(6) % Hdr % MO8Type = MO8_P_ECBT
Pfields(6) % Hdr % PPCode  = PP_P_ECBT
Pfields(6) % Hdr % LBVC    = VC_P_ECBT
ALLOCATE( Pfields(6) % RData(Pfields(21) % Hdr % NumCols,     &
                             Pfields(21) % Hdr % NumRows) )
Pfields(6) % RData(:,:) = RMDI

Pfields(7) % Hdr % STCode  = ST_I_ECBB
Pfields(7) % Hdr % MO8Type = MO8_I_ECBB
Pfields(7) % Hdr % PPCode  = PP_I_ECBB
Pfields(7) % Hdr % LBVC    = VC_I_ECBB
ALLOCATE( Pfields(7) % RData(Pfields(21) % Hdr % NumCols,     &
                             Pfields(21) % Hdr % NumRows) )
Pfields(7) % RData(:,:) = RMDI

Pfields(8) % Hdr % STCode  = ST_I_ECBT
Pfields(8) % Hdr % MO8Type = MO8_I_ECBT
Pfields(8) % Hdr % PPCode  = PP_I_ECBT
Pfields(8) % Hdr % LBVC    = VC_I_ECBT
ALLOCATE( Pfields(8) % RData(Pfields(21) % Hdr % NumCols,     &
                             Pfields(21) % Hdr % NumRows) )
Pfields(8) % RData(:,:) = RMDI

Pfields(9) % Hdr % STCode  = ST_CBHorE
Pfields(9) % Hdr % MO8Type = MO8_CBHorE
Pfields(9) % Hdr % PPCode  = PP_CBHorE
Pfields(9) % Hdr % LBVC    = VC_CBHorE
ALLOCATE( Pfields(9) % RData(Pfields(21) % Hdr % NumCols,     &
                             Pfields(21) % Hdr % NumRows) )
Pfields(9) % RData(:,:) = RMDI



!---- Generate initial Cb and layer cloud masks as follows:
!
! If convective PPTN rate greater than .00025 then Cb. Put into Pfields(3).
! If layer cloud > .95 on any model level then layer cloud. Put into Pfields(7).
!   This mask will be used to generate embedded Cb mask later.
!  
!Debi Turp: why not use Pfields(12) instead of Pfields(10) below?
 
DO k = 1, NumLevs
  DO j = 1, PFields(1) % Hdr % NumRows
    DO i = 1, PFields(1) % Hdr % NumCols
      IF ( Pfields(10)%RData(i,j)   > 0.00025    )  THEN
         Pfields(3)%RData(i,j) = 1.
      ELSE
         Pfields(3)%RData(i,j) = RMDI
      END IF
      IF  (  Pfields(40+k) % RData(i,j)   > 0.95  )  THEN
         Pfields(7) % RData(i,j) = 1.
      END IF
    END DO
  END DO
END DO


!---  Adjust mask fields using subroutine FILL_N_DSPEC to fill in small gaps 
!     and smooth edges. Subroutine iterates 30 times.
!
! Note: For Cb mask, need to copy mask field into Pfields(4) first
!       as this is used in the iteration process in FILL_N_DSPEC
!       but don't need to do this for the layer cloud mask.

! Cb mask
Pfields(4) % RData(:,:) = Pfields(3) % RData(:,:)

! DEPENDS ON: fill_n_dspec
Call FILL_N_DSPEC(30,Pfields,3,4,ErrorStatus)

Pfields(4) % RData(:,:) = Pfields(3) % RData(:,:)


! Embedded Cb mask
! DEPENDS ON: fill_n_dspec
Call FILL_N_DSPEC(30,Pfields,7,8,ErrorStatus)

Pfields(8) % RData(:,:) = Pfields(7) % RData(:,:)
Pfields(7) % RData(:,:) = RMDI



!---- Identify embedded Cbs and generate embedded Cb mask
!
! There are no output fields of layer cloud only with which to identify
! embedded Cbs. In order to classify Cbs as embedded the following
! criteria are used. A cell which has a Cb will be considered as
! embedded if at least one of the 8 adjacent cells is not a Cb AND is
! not clear of any cloud.
!
! field 3 contains Cb mask
! field 8 contains layer cloud mask
! Output field 7 will contain the embedded Cb mask

DO j = 2, PFields(1) % Hdr % NumRows  -1
  DO i = 2, PFields(1) % Hdr % NumCols -1
    p1 = Pfields(3) % RData(i-1,j+1)
    p2 = Pfields(3) % RData(i,j+1)
    p3 = Pfields(3) % RData(i+1,j+1)
    p4 = Pfields(3) % RData(i-1,j)
    p5 = Pfields(3) % RData(i,j)
    p6 = Pfields(3) % RData(i+1,j)
    p7 = Pfields(3) % RData(i-1,j-1)
    p8 = Pfields(3) % RData(i,j-1)
    p9 = Pfields(3) % RData(i+1,j-1)
    q1 = Pfields(8) % RData(i-1,j+1)
    q2 = Pfields(8) % RData(i,j+1)
    q3 = Pfields(8) % RData(i+1,j+1)
    q4 = Pfields(8) % RData(i-1,j)
    q5 = Pfields(8) % RData(i,j)
    q6 = Pfields(8) % RData(i+1,j)
    q7 = Pfields(8) % RData(i-1,j-1)
    q8 = Pfields(8) % RData(i,j-1)
    q9 = Pfields(8) % RData(i+1,j-1)
    IF ((((q1 > 0.) .AND. (p1 == RMDI))  .OR. &
         ((q2 > 0.) .AND. (p2 == RMDI))  .OR. &
         ((q3 > 0.) .AND. (p3 == RMDI))  .OR. &
         ((q4 > 0.) .AND. (p4 == RMDI))  .OR. &
         ((q6 > 0.) .AND. (p6 == RMDI))  .OR. &
         ((q7 > 0.) .AND. (p7 == RMDI))  .OR. &
         ((q8 > 0.) .AND. (p8 == RMDI))  .OR. &
         ((q9 > 0.) .AND. (p9 == RMDI)) ) .AND. (p5 == 1.) ) THEN
      Pfields(7) % RData(i,j) = 1.
    END IF
  END DO
END DO



!---- Recheck embedded Cb mask ----
!
! If cell is marked as having an embedded Cb, check the surrounding 
! 8 cells in turn - if they contain Cbs, mark them as embedded Cbs too.
!
!Comment on optimisation below:
!kk Optimisation exchanged I and J loop to avoid bank conflicts
!kk In the I loop, there are dependencies in Pfields(7)%RData which
!kk prevent vectorization. Running this loop with stride 3 allows
!kk vectorization, but slightly changes the algorithm .

DO k=1,50
  DO jj=0,2
    DO j = jj+2, PFields(1) % Hdr % NumRows  -1, 3
      DO ii=0,2
        DO i = ii+2, PFields(1) % Hdr % NumCols -1, 3
          p1=Pfields(3) % RData(i-1,j+1)
          p2=Pfields(3) % RData(i  ,j+1)
          p3=Pfields(3) % RData(i+1,j+1)
          p4=Pfields(3) % RData(i-1,j  )
          p6=Pfields(3) % RData(i+1,j  )
          p7=Pfields(3) % RData(i-1,j-1)
          p8=Pfields(3) % RData(i  ,j-1)
          p9=Pfields(3) % RData(i+1,j-1)

          q5=Pfields(7) % RData(i  ,j  )
          IF ( q5 ==1. ) THEN
            IF ( p1 == 1. ) THEN
              Pfields(7) % RData(i-1,j+1) = 1.
            END IF
            IF ( p2 == 1. ) THEN
              Pfields(7) % RData(i,j+1) = 1.
            END IF
            IF ( p3 == 1. ) THEN
              Pfields(7) % RData(i+1,j+1) = 1.
            END IF
            IF ( p4 == 1. ) THEN
              Pfields(7) % RData(i-1,j) = 1.
            END IF
            IF ( p6 == 1. ) THEN
              Pfields(7) % RData(i+1,j) = 1.
            END IF
            IF ( p7 == 1. ) THEN
              Pfields(7) % RData(i-1,j-1) = 1.
            END IF
            IF ( p8 == 1. ) THEN
              Pfields(7) % RData(i,j-1) = 1.
            END IF
            IF ( p9 == 1. ) THEN
              Pfields(7) % RData(i+1,j-1) = 1.
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO
END DO




!---- Identify pressure at base and top of isolated and embedded Cbs ----
! Use:
! Pfields(3)  containing Cb mask
! Pfields(7)  containing Cb mask
! Pfields(80+k)  containing convective cloud fraction at level k
! Pfields(120+k) containing pressure at level k
!
! Output:
! Pfields(1) will contain Cb base pressure
! Pfields(2) will contain Cb top pressure
! Pfields(5) will contain embedded Cb base pressure
! Pfields(6) will contain embedded Cb top pressure
!
! Do this by looping over model levels bottom to top. 
!
! For cloud bases:
!  If Cb mask indicates Cb in this cell, and convective cloud amount at this 
!  level > 0, and no base has yet been identified, set the Cb base to this 
!  level. Repeat for embedded Cb using embedded Cb mask. 
!
! For cloud tops:
!  If Cb mask indicates Cb in this cell, and convective cloud amount at this 
!  level > 0, and cloud base has been set and is lower than or at this level
!  then set the cloud top to be this level. This is overwritten each time the 
!  model level increases and these conditions are met, until the conditions  
!  aren't met anymore. Repeat for embedded Cb using embedded Cb mask. 
!

DO k= 1, numlevs
  DO j = 1, PFields(1) % Hdr % NumRows
    DO i = 1, PFields(1) % Hdr % NumCols
!  find cloud Bases
      IF ( ( Pfields(80+k) % RData(i,j)  /= 0. )  .AND.  &
           ( Pfields(3) % RData(i,j)     /= RMDI )  .AND.  &
           ( Pfields(1) % RData(i,j)     == RMDI )  )  THEN
        Pfields(1) % RData(i,j) = Pfields(120+k) % RData(i,j)
      END IF

      IF ( ( Pfields(80+k) % RData(i,j)  /= 0. )  .AND.  &
           ( Pfields(7) % RData(i,j)     /= RMDI )  .AND.  &
           ( Pfields(5) % RData(i,j)     == RMDI )  )  THEN
        Pfields(5) % RData(i,j) = Pfields(120+k) % RData(i,j)
      END IF

! Find cloud tops
      IF ( ( Pfields(80+k) % RData(i,j)  /= 0. ) .AND.  &
           ( Pfields(3) % RData(i,j)     /= RMDI ) .AND.  &
           ( Pfields(1) % RData(i,j) >= Pfields(120+k) % RData(i,j))) THEN
        Pfields(2) % RData(i,j) = Pfields(120+k) % RData(i,j)
      END IF

      IF ( ( Pfields(80+k) % RData(i,j)  /= 0. ) .AND.  &
           ( Pfields(7) % RData(i,j)     /= RMDI ) .AND.  &
           ( Pfields(5) % RData(i,j) >= Pfields(120+k) % RData(i,j))) THEN
        Pfields(6) % RData(i,j) = Pfields(120+k) % RData(i,j)
      END IF
    END DO
  END DO
END DO



!---- Test that Cb / embedded Cb tops are above bases ----
!
! Do this by calculating the difference in pressure between the cloud base
!  and top. If this is negative, the cloud base is above cloud top. In this
!  case set the flag minnm (for isolated Cb) or minnmb (for embedded Cb) 
!  to be equal to the difference. This is then used to reset all values in the
!  Cb base pressure or embedded Cb base pressure to the min value. (I.e. if 
!  any one point is wrong, all values in the field are set to the min value.)

minnm=1.
minnmb=1.

DO j = 1, PFields(1) % Hdr % NumRows
   DO i = 1, PFields(1) % Hdr % NumCols
     diffn = Pfields(1)%RData(i,j) - Pfields(2)%RData(i,j)
     IF (diffn < 0.) THEN 
       minnm = diffn
     END IF
     diffnb = Pfields(5)%RData(i,j) - Pfields(6)%RData(i,j)
     IF (diffnb < 0.) THEN
       minnmb = diffnb
     END IF
  END DO
END DO

DO j = 1, PFields(1) % Hdr % NumRows
   DO i = 1, PFields(1) % Hdr % NumCols
     IF (minnm < 0.) THEN
       Pfields(1)%RData(i,j) = minnm
     END IF
     IF (minnmb < 0.) THEN
       Pfields(5)%RData(i,j) = minnmb
     END IF
  END DO
END DO



!---- Calculate Cb horizontal extent ----
! Input:
!   Pfields(3)    contains Cb mask
!   Pfields(80+k) contains convective cloud fraction at level k
!
! Output:
!   Pfields(9)    containing Cb horizontal extent, given as index 
!                 with value 0, 1, 2, or 3                      
!
!  For each gridpoint:
!  1) Set Pfields(9)=largest convective cloud fraction in the column
!  2) Reset Pfields(9) as an index as follows:
!      0  < max conv cloud fraction <  .01 then  Pfields(9)=0
!    .01 <= max conv cloud fraction <= .5  then  Pfields(9)=1
!     .5  < max conv cloud fraction <= .75 then  Pfields(9)=2
!     .75 < max conv cloud fraction  then  Pfields(9)=3
! 

DO k = 1, NumLevs
  DO j = 1, PFields(1) % Hdr % NumRows
    DO i = 1, PFields(1) % Hdr % NumCols
      IF ( (Pfields(80+k) % RData(i,j) > Pfields(9) % RData(i,j)) .AND. &
           (Pfields(3) % RData(i,j) > 0. )  ) THEN
        Pfields(9) % RData(i,j) = Pfields(80+k) % RData(i,j)
      END IF
    END DO
  END DO
END DO

DO j = 1, PFields(1) % Hdr % NumRows
  DO i = 1, PFields(1) % Hdr % NumCols
    pcnt = RMDI
    IF ( ( Pfields(9) % RData(i,j) > 0.0 ) .AND. &
         ( Pfields(9) % RData(i,j) < .01) ) THEN
      pcnt = 0.
    END IF
    IF ( ( Pfields(9) % RData(i,j) >= .01) .AND. &
         ( Pfields(9) % RData(i,j) <= .5) ) THEN
      pcnt = 1.
    END IF
    IF ( ( Pfields(9) % RData(i,j) > .5) .AND. &
         ( Pfields(9) % RData(i,j) <= .75) ) THEN
      pcnt = 2.
    END IF
    IF ( Pfields(9) % RData(i,j) > .75) THEN
      pcnt = 3.
    END IF
    Pfields(9) % RData(i,j) = pcnt
  END DO
END DO



!---- Fill in field values ----
!
! If Nofill=1 use Fill_Pressure to check/fill in missing output field 
! values where Cb / embedded Cb mask indicates presence of Cb.
!
! Inputs/output for subroutine Fill_Pressure:
! In_field  - index of input field to be filled in. This will be adjusted  
!               within subroutine. On output will hold the adjusted field. 
! In_mask   - index of mask field used for comparison (will be either 3 or 7)
! Out_field - index of dummy field which will be overwritten and used within 
!               subroutine to compare field values. 
!
! On output, the adjusted field will be in both Pfields(In_field) and
! Pfields(Out_field), and either could be used. The code below uses 
! Pfields(In_field) for the adjusted field and Pfields(Out_field) as a 
! dummy field which is overwritten.
!

Nofill=1
IF ( Nofill == 1 ) THEN

 
!   Fill in Cloud Base Pressure values within area of Cb mask.

  In_field = 1     ! Field to be adjusted (will be replaced by adjusted field)
  In_mask = 3      ! Mask field
  Out_field = 4    ! Dummy field, will be overwritten in the subroutine
                   !  and replaced by adjusted field

! DEPENDS ON: fill_pressure
  CALL Fill_Pressure(In_field,In_mask,Out_field,Pfields,ErrorStatus)

 
!   Fill in Cloud Top Pressure values within area of Cb mask.

  In_field = 2     ! Field to be adjusted (will be replaced by adjusted field)
  In_mask = 3      ! Mask field
  Out_field = 4    ! Dummy field, will be overwritten in the subroutine
                   !  and replaced by adjusted field
                   
! DEPENDS ON: fill_pressure
  CALL Fill_Pressure(In_field,In_mask,Out_field,Pfields,ErrorStatus)

 
!   Fill in Cloud Base Pressure values within area of embedded Cb mask.

  In_field = 5     ! Field to be adjusted (will be replaced by adjusted field)
  In_mask = 7      ! Mask field
  Out_field = 8    ! Dummy field, will be overwritten in the subroutine
                   !  and replaced by adjusted field
                   
! DEPENDS ON: fill_pressure
  CALL Fill_Pressure(In_field,In_mask,Out_field,Pfields,ErrorStatus)

 
!   Fill in Cloud Top Pressure values within area of embedded Cb mask.

  In_field = 6     ! Field to be adjusted (will be replaced by adjusted field)
  In_mask = 7      ! Mask field
  Out_field = 8    ! Dummy field, will be overwritten in the subroutine
                   !  and replaced by adjusted field
                   
! DEPENDS ON: fill_pressure
  CALL Fill_Pressure(In_field,In_mask,Out_field,Pfields,ErrorStatus)

 
!   Fill in Cb horizontal extent values within area of Cb mask.

  In_field = 9      ! Field to be adjusted (will be replaced by adjusted field)
  In_mask = 3       ! Mask field
  Out_field = 8     ! Dummy field, will be overwritten in the subroutine
                    !  and replaced by adjusted field
                   
! DEPENDS ON: fill_pressure
  CALL Fill_Pressure(In_field,In_mask,Out_field,Pfields,ErrorStatus)


END IF


9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE ConvAct
#endif
