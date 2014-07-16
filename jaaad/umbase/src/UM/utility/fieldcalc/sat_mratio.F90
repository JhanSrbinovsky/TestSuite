#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Routine to calculate saturated mixing ratio


SUBROUTINE Sat_MRatio( NumFlds,       &  ! in
                       Fields1,       &  ! in 
                       Fields2,       &  ! in 
                       SatMR,         &  ! inout 
                       ErrorStatus)      ! inout

! Description:
!   Calculates the saturated mixing ratio (for use in relative humidity
!   calculation) given temperature and pressure as input.
!
! Method:
!   The headers for the output fields are copied from the first vector of input
!   fields (i.e. pressure). This should be of no consequence unless x_sat is 
!   being calculated in its own right, as in the relative humidity calculation
!   the headers are taken from the specific humidity field.
!
!   x_{sat} is calculated as:
!   x_{sat} = ( 0.622 * e_{s} ) / ( p - e_{s} )
!   where e_{s} is the saturated vapour pressure (SVP): 
!         e_{s} = 6.11 * 10^( 7.5*T / (237.7 + T) )
!   using the Murray formulation, with T in deg C and p in hPa
!   
!
! Owner: Dave Robinson
!
! History:
! Version Date     Comment
! ------- ----     -------
! 1.0     17/05/07 Original Code.  Sally Close
! 1.1     14/04/08 Adapted for input to UM.  Debi Turp
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumFlds
TYPE(PP_Field_type), INTENT(IN) :: Fields1(NumFlds)        ! pressure
TYPE(PP_Field_type), INTENT(IN) :: Fields2(NumFlds)        ! temperature

TYPE(PP_Field_type), INTENT(INOUT) :: SatMR(NumFlds)       ! saturated 
                                                           !   mixing ratio
INTEGER, INTENT(INOUT) :: ErrorStatus   

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Sat_MRatio"
#include "c_mdi.h"
#include "c_epslon.h"
#include "c_0_dg_c.h"

! Local Variables:
INTEGER :: i                                               ! Loop counter

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )


IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF


! Allocate output fields
DO i = 1, NumFlds

  IF ( (Fields1(i) % Hdr % NumCols /= Fields2(i) % Hdr % NumCols) .OR.     &
       (Fields1(i) % Hdr % NumRows /= Fields2(i) % Hdr % NumRows) ) THEN  
! DEPENDS ON: ereport
    CALL EReport( RoutineName, ErrorStatus,                                &
                  "Cannot multiply fields of different dimensions")
    ErrorStatus = StatusWarning
    GO TO 9999
  END IF
 
  IF ( ASSOCIATED( SatMR(i) % RData ) ) THEN
    DEALLOCATE( SatMR(i) % RData )
    NULLIFY( SatMR(i) % RData )
  END IF

  IF ( .NOT. ASSOCIATED( SatMR(i) % RData ) ) THEN
    ALLOCATE( SatMR(i) % RData(Fields1(i) % Hdr % NumCols,                 &
                               Fields1(i) % Hdr % NumRows) )
  END IF

END DO


! Set header
SatMR % Hdr = Fields1 % Hdr


! Calculate saturated mixing ratio using the two equations described
! in the header combined into one. Temperature is converted into 
! degrees C and pressure into hPa, hence in the equation 
!   T = Fields2(i)%RData - 273.15
!   p = Fields1(i)%RData / 100.0
! Note also the part of the equation (237.7+T) has been reduced to
!   Fields2(i)%RData - 35.45 
! which is equivalent to 237.7 + Fields2(i)%RData - 273.15

DO i = 1, NumFlds
  WHERE ( (Fields1(i) % RData /= Fields1(i) % Hdr % BMDI) .AND.  &
          (Fields2(i) % RData /= Fields2(i) % Hdr % BMDI) )

    SatMR(i) % RData = ( Epsilon * &
                       ( 6.11 * &
                       ( 10.0**( (7.5 * (Fields2(i) % RData - ZeroDegC)) &
                                      / (Fields2(i) % RData - 35.45) ) )  ) &
!                                              
                        /  &
!                       
                        ( (Fields1(i) % RData / 100.0) - &                      
                          ( 6.11 * &
                          ( 10.0**( (7.5 * (Fields2(i) % RData - ZeroDegC)) &
                                      / (Fields2(i) % RData - 35.45) ) ) ) ) )


  ELSEWHERE
    SatMR(i) % RData = SatMR(i) % Hdr % BMDI
  END WHERE

END DO


9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE Sat_MRatio
#endif
