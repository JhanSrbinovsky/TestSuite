#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
SUBROUTINE Maximum( Field1,       &  ! in                                    
                    Field2,       &  ! in                                    
                    MaxField,     &  ! inout                                 
                    ErrorStatus )    ! inout                                 

! Description:                                                           
!   Finds the maximum of two fields at each point and outputs a          
!   field holding the maximum values.                                    
!                                                                        
! Method:                                                                
!   1) Allocates the maximum field structure                             
!   2) Where there are data values in at least one field, it takes       
!      the maximum of the two input fields using the Fortran             
!      MAX function.                                                     
!                                                                        
! Owner: Dave Robinson                                                   
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
TYPE(PP_Field_type), INTENT(IN)    :: Field1      ! Fields to find 
TYPE(PP_Field_type), INTENT(IN)    :: Field2      !    maximum of                 
TYPE(PP_Field_type), INTENT(INOUT) :: MaxField    ! Max of two fields      
INTEGER, INTENT(INOUT)             :: ErrorStatus                                    
INTEGER :: i,j                                                           

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Maximum"                   
#include "c_mdi.h"                                                 

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN                                      
  ! Previous error - do not proceed                                      
  GO TO 9999                                                             
END IF                                                                   
                                                                         
IF ( (Field1 % Hdr % NumCols /= Field2 % Hdr % NumCols) .OR. &           
     (Field1 % Hdr % NumRows /= Field2 % Hdr % NumRows) ) THEN           
! DEPENDS ON: ereport
  CALL EReport( RoutineName, StatusWarning, &                            
                "Cannot find maximum of fields of different dimensions" )
  ErrorStatus = StatusWarning                                            
  GO TO 9999                                                             
END IF                                                                   

! If there is data in MaxField, and it is neither of the two input       
! fields, get rid of it                                                  
IF ( ASSOCIATED( MaxField % RData )  ) THEN                              
  DEALLOCATE( MaxField % RData )                                         
END IF                                                                   
                                   
MaxField % Hdr = Field1 % Hdr                                            
MaxField % Hdr % BMDI = RMDI                                             
 
ALLOCATE( MaxField % RData(Field1 % Hdr % NumCols, &
                           Field1 % Hdr % NumRows) )                     


!Find the maximum of the two fields                                      
DO j = 1, Field1 % Hdr % NumRows                                               

  DO i = 1, Field1 % Hdr % NumCols                                       
  
    IF (Field1 % RData(i,j) >= Field2 % RData(i,j)) THEN               
      MaxField % RData(i,j) = Field1 % RData(i,j)                        
    ELSE                                                                 
      MaxField % RData(i,j) = Field2 % RData(i,j)                          
    END IF                                                               

  END DO                                                                       

END DO   
  

9999 CONTINUE                                                             

END SUBROUTINE Maximum                                                    
#endif
