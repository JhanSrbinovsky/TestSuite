#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Subroutine to extract data from D1 array. Converted from
!   REDIST_STOCHEM.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!  Use of an interface block means that a generic call can select
!  the appropriate one of these routines automatically from the
!  input variables.
!
!   Called from UKCA_MAIN1.
!
! Current code owner: Colin Johnson/Olaf Morgenstern
!                     Fiona O'Connor
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE UKCA_EXTRACT_D1_DATA3D(                             &                         
#include "argd1.h"                                                
         first,N,X)                                                     
                                                                        
        USE UKCA_D1_DEFS                                                
        IMPLICIT NONE                                                   
                                                                       
#include "parvars.h"                                            
#include "typsize.h"                                            
#include "typd1.h"                                                
                                                                        
        INTEGER, INTENT(IN) :: N         ! id of array                  
        LOGICAL, INTENT(IN) :: first                                    
        REAL, INTENT(INOUT) :: X(:,:,:)  ! extracted array              
                                                                        
        REAL, ALLOCATABLE   :: data1(:)                                 
        CHARACTER(LEN=72)   :: cmessage                                 
        INTEGER             :: errcode                                  
                                                                        
        IF(first .AND. SIZE(X) /= UkcaD1Codes(N)%length) THEN           
          errcode = N                                                   
          cmessage='Array sizes in local variable and D1 do not agree'  
          WRITE(6,*) cmessage, ' Error code: ',errcode,' PE: ',mype     
          WRITE(6,*) 'Expected size: ',SIZE(X),                        & 
                   ' Length in D1: ',UkcaD1Codes(N)%length              
! DEPENDS ON: ereport
          CALL EREPORT('extract_d1_data2d',ERRCODE,cmessage)            
        ELSE                                                            
          ALLOCATE(data1(UkcaD1Codes(N)%length))                        
          data1 = D1(UkcaD1Codes(N)%address:UkcaD1Codes(N)%address     & 
                  +UkcaD1Codes(N)%length-1)                             
          X=RESHAPE(data1,(/SIZE(X,DIM=1),SIZE(X,DIM=2),SIZE(X,DIM=3)/))
          DEALLOCATE(data1)                                             
        ENDIF                                                           
                                                                        
        RETURN                                                          
        END SUBROUTINE UKCA_EXTRACT_D1_DATA3D        
#endif
