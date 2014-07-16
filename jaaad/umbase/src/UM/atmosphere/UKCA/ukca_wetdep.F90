#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT******************************* 
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT******************************* 
!                                                                       
! Purpose: Subroutine to to assign wet deposition rates in s-1 to       
!          array dpw in module asad_mod. This is passed for use in        
!          ASAD chemistry                                               
!          Based on routine from Cambridge TOMCAT model     
!                                                        
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD routine CDRIVE.                             
!                                                                       
! Current code owner: Colin Johnson/Olaf Morgenstern/Fiona O'Connor                                    
!                                                                       
!                                                                       
! Code description:                                                     
!   Language: FORTRAN 90                                                
!   This code is written to UMDP3 v6 programming standards.             
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
      SUBROUTINE UKCA_WETDEP(nlev,wetrt,n_points)                       

      USE ASAD_MOD,         ONLY: ndepw, nldepw, dpw
      IMPLICIT NONE                                                     

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "parvars.h"                                    
#include "typsize.h"                                    

        INTEGER, INTENT(IN) :: nlev                      ! Model level n
        INTEGER, INTENT(IN) :: n_points                  ! No of spatial
                                                                        
        REAL, INTENT(IN) :: wetrt(n_points,model_levels,jpdw) ! Wet dep rates
                                                                        
!       Local variables                                                 
                                                                        
        INTEGER :: i                         ! Loop variable            
        INTEGER :: js                        ! Loop variable            
        INTEGER :: nspec                     ! Pointer for species      
                                                                        
        DO js = 1,ndepw                                                 
          nspec = nldepw(js)                                            
          DO i = 1,n_points                                             
            dpw(i,nspec)= wetrt(i,nlev,js)                              
          ENDDO                                                         
        ENDDO                                                           
                                                                        
        RETURN                                                          
        END SUBROUTINE UKCA_WETDEP                                      
#endif

