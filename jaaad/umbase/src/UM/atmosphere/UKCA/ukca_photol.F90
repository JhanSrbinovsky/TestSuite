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
! Purpose: Subroutine to to assign photolysis rates in s-1 to           
!          array rk in module asad_mod. This is passed for use in         
!          asad chemistry.                                              
!          Based on photol.F from Cambridge TOMCAT model         
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!                                                                       
!          Called from ASAD routine CDRIVE.                             
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
        SUBROUTINE UKCA_PHOTOL(prt,n_points)                            

        USE ASAD_MOD,        ONLY:  rk, nprkx
        IMPLICIT NONE                                                   

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"

        INTEGER, INTENT(IN) :: n_points         ! No of spatial points  
                                                                        
        REAL, INTENT(IN) :: prt(n_points,jppj)  ! Photolysis rates      
                                                                        
! Local variables                                                       
                                                                        
        INTEGER :: jl                           ! Loop variable         
        INTEGER :: jr                           ! Loop variable         
                                                                        
        DO jl = 1, n_points                                             
          DO jr = 1, jppj                                               
            rk(jl,nprkx(jr)) = prt(jl,jr)                               
          ENDDO                                                         
        ENDDO                                                           
                                                                        
        RETURN                                                          
        END SUBROUTINE UKCA_PHOTOL 
#endif
