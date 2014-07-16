#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT******************************* 
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT******************************* 
!                                                                       
! Description:
!  To set array dimensions using halo type.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Current Code Owner:       Colin Johnson/Olaf Morgenstern
!                           Fiona O'Connor
!
!
!  Code Description:
!   Language:  FORTRAN 90 (formatted)
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!                                                                       
      SUBROUTINE UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)                   
      USE UKCA_D1_DEFS                                                  
      IMPLICIT NONE                                                     
                                                                        
      INTEGER, INTENT(IN)     :: N  ! id of array                       
      INTEGER,  INTENT(OUT)   :: I1,I2,J1,J2 ! array bounds          
                                                                        
      CHARACTER(LEN=72)       :: cmessage                                     
      INTEGER :: ERRCODE                                                
      INTEGER :: halox      ! halo EW                                   
      INTEGER :: haloy      ! halo EW                                   
                                                                        
#include "parvars.h"                                            
                                                                        
!     set haloes

      halox=halosize(1,ukcaD1codes(N)%halo_type)                        
      haloy=halosize(2,ukcaD1codes(N)%halo_type)                        
                                                                        
      I1 = 1-halox                                                        
      I2 = ukcaD1codes(N)%len_dim1+halox                                  
      J1 = 1-haloy                                                        
      J2 = ukcaD1codes(N)%len_dim2+haloy                                  
                                                                        
      END SUBROUTINE UKCA_SET_ARRAY_BOUNDS    
#endif

