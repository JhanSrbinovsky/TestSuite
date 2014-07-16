#if defined(FLDCALC)   
! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT 2004, Met Office, All Rights Reserved.            
! Please refer to file $UMDIR/vn$VN/copyright.txt for further details    
! *****************************COPYRIGHT*******************************  
!                                              
!+ Routine to calculate surface dust conc and ~2000-5000ft conc

SUBROUTINE DUST_CON (Numlevs,       & !in
                     Orog,          & !in
                     PFields,       & !in
                     Tfields,       & !in
                     Zfields,       & !in
                     dust1,         & !in
                     dust2,         & !in
                     dust3,         & !in
                     dust4,         & !in
                     dust5,         & !in
                     dust6,         & !in
                     dustc,         & !inout
                     dusty,         & !inout
                     ErrorStatus )    !inout

! Description: Calculates dust concentrations.
!
! Owner:       Glenn Greed
!
! Code Description:     
!   Language:           Fortran 90  
!   Software Standards: UMDP3 v6  
!

USE IO_Mod, ONLY:         & 
  PP_Header_type,         & 
  PP_Field_type            
USE Err_Mod, ONLY:        &
  StatusOK  
USE FldCodes_Mod, ONLY:             & 
  ST_Ptheta, ST_Ttheta,             &
  ST_Dust1,  ST_Dust2,  ST_Dust3,   &
  ST_Dust4,  ST_Dust5,  ST_Dust6,   &
  ST_DUSTC,  MO8_DUSTC,  PP_DUSTC,  & 
  ST_DUSTY,  MO8_DUSTY,  PP_DUSTY,  & 
  VC_Surface,                       &
  LV_Surface,  LV_Special
  
IMPLICIT None   
                              
! Subroutine Arguments:   
INTEGER, INTENT(IN) :: NumLevs 
TYPE(PP_Field_type), INTENT(IN) :: Orog              ! Model Orography
TYPE(PP_Field_type), INTENT(IN) :: ZFields(NumLevs)  !height on rho levels
TYPE(PP_Field_type), INTENT(IN) :: PFields(NumLevs)  !Pressure on theta,
TYPE(PP_Field_type), INTENT(IN) :: TFields(NumLevs)  !temperature on theta 
TYPE(PP_Field_type), INTENT(IN) :: dust1(NumLevs)    !and the  
TYPE(PP_Field_type), INTENT(IN) :: dust2(NumLevs)    !6 dust 
TYPE(PP_Field_type), INTENT(IN) :: dust3(NumLevs)    !particle size 
TYPE(PP_Field_type), INTENT(IN) :: dust4(NumLevs)    ! bins on theta
TYPE(PP_Field_type), INTENT(IN) :: dust5(NumLevs)
TYPE(PP_Field_type), INTENT(IN) :: dust6(NumLevs)

TYPE(PP_Field_type), INTENT(INOUT) :: dustc  ! surface dust concentration 
TYPE(PP_Field_type), INTENT(INOUT) :: dusty  ! 2000-5000ft dust concentration 

INTEGER, INTENT(INOUT) :: ErrorStatus 

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "DUST_CON" 

#include "c_r_cp.h"
#include "c_mdi.h"

! Layer dimensions in metres.
REAL, PARAMETER :: top_of_layer = 1524.0 
REAL, PARAMETER :: bottom_of_layer = 609.0 
REAL, PARAMETER :: layer_thickness = 915.0

! Local Variables:  
INTEGER :: i,j,k           ! Loop counters

! density of air on given model level  (P/RT) 
REAL :: rho_theta (dust1(1)%Hdr%NumCols, dust1(1)%Hdr%NumRows, Numlevs)
!  working array for total dust MMR
REAL :: tot_mmr (dust1(1)%Hdr%NumCols, dust1(1)%Hdr%NumRows, Numlevs)
!  thickness of dust layer
REAL :: thick (dust1(1)%Hdr%NumCols, dust1(1)%Hdr%NumRows)

! End of header -----------------------------------------------

CALL Timer( RoutineName, 3 ) 

IF ( ErrorStatus /= StatusOK ) THEN 
  ! Previous error - do not proceed  
  GO TO 9999         
END IF

IF ( ASSOCIATED( dustc % RData ) ) THEN   
  DEALLOCATE( dustc % RData )    
END IF  
IF ( ASSOCIATED( dusty % RData ) ) THEN   
  DEALLOCATE( dusty % RData )    
END IF  

dustc % Hdr = dust1(1) % Hdr
dustc % Hdr % LBVC     = VC_Surface
dustc % Hdr % MO8Level = LV_Surface
dustc % Hdr % BULEV    = 0.0 
dustc % Hdr % BHULEV   = 0.0 
dustc % Hdr % RLevel   = 0.0  
dustc % Hdr % RefLevel = 0.0 
dustc % Hdr % BHLEV    = 0.0  
dustc % Hdr % BHRLEV   = 0.0   
dustc % Hdr % BMDI     = RMDI
dustc % Hdr % PPCode   =  PP_DUSTC 
dustc % Hdr % MO8Type  = MO8_DUSTC 
dustc % Hdr % STCode   =  ST_DUSTC

dusty % Hdr = dustc % Hdr
dusty % Hdr % LBPROC   = 2048  !  weighted mean between two levels
dusty % Hdr % LBVC     = VC_Surface 
dusty % Hdr % MO8Level = LV_Special
dusty % Hdr % RLevel   = 5000.0 ! 5000ft (1524m) 
dusty % Hdr % RefLevel = 2000.0 ! 2000ft (609m)
dusty % Hdr % PPCode   =  PP_DUSTY 
dusty % Hdr % MO8Type  = MO8_DUSTY 
dusty % Hdr % STCode   =  ST_DUSTY

ALLOCATE( dustc % RData(dustc % Hdr % NumCols, &
                        dustc % Hdr % NumRows) ) 
ALLOCATE( dusty % RData(dusty % Hdr % NumCols, &
                        dusty % Hdr % NumRows) ) 

! loop over number of levels
! to calc total MMR at each point
! and the density at each point, in g per m cubed

DO k = 1,NumLevs 
   DO i = 1,dustc % Hdr % NumCols 
      DO j = 1,dustc % Hdr % NumRows 

        tot_mmr(i,j,k)= ( dust1(k)%RData(i,j) + dust2(k)%RData(i,j) + &
                          dust3(k)%RData(i,j) + dust4(k)%RData(i,j) + &
                          dust5(k)%RData(i,j) + dust6(k)%RData(i,j) )
              
        rho_theta(i,j,k)= ( ( PFields(k)%RData(i,j) * 1000.0)  / &
                            ( TFields(k)%RData(i,j) * R     ) )
                    
      END DO     
   END DO
END DO

! loop over model levels 
! here we consider conc between 2000 and 5000ft agl
! so approx 609 and 1524m and mask out data if outside this layer

dusty%RData(:,:) = 0.0

DO k = 2,Numlevs-1 

  thick=( ZFields(k+1)%RData - ZFields(k)%RData )
  
  WHERE ( (ZFields(k+1)%RData - Orog%RData) <  bottom_of_layer )
        thick(:,:)=0.0
  END WHERE
  
  WHERE ( ZFields(k)%RData - Orog%RData >  top_of_layer )
       thick(:,:)=0.0
  END WHERE     

  DO i = 1,dustc % Hdr % NumCols 
     DO j = 1,dustc % Hdr % NumRows 
       
       dusty%RData(i,j) = dusty%RData(i,j) +    &
                         ( thick(i,j) * tot_mmr(i,j,k) * rho_theta(i,j,k) )

     END DO     
  END DO

END DO

DO i = 1,dustc % Hdr % NumCols 
   DO j = 1,dustc % Hdr % NumRows 
  
     ! for final conc scale by ~thickness layer 
     dusty%RData(i,j) = dusty%RData(i,j) / layer_thickness
  
     ! surface conc much simpler here we take level 1
     dustc%RData(i,j) = tot_mmr(i,j,1) * rho_theta(i,j,1)

   END DO     
END DO


9999 CONTINUE 

CALL Timer( RoutineName, 4 )

END SUBROUTINE DUST_CON
#endif
