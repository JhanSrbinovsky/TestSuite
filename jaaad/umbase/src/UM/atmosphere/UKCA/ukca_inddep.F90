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
! Purpose: Routine to read in dry deposition velocities at 1 metre.      
!          Original version taken from the Cambridge TOMCAT model.      
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!                                                                       
!          Called from ASAD routine ASAD_CINIT.                         
!                                                                       
! Current code owner: Colin Johnson/Olaf Morgenstern/
!                     Fiona O'Connor                                    
!                                                                       
! Code description:                                                     
!   Language: FORTRAN 90                                                
!   This code is written to UMDP3 v6 programming standards.             
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
        SUBROUTINE UKCA_INDDEP                                          

        USE ASAD_MOD,           ONLY: depvel, jddepc, jddept
        USE UKCA_CHEM1_DAT
        IMPLICIT NONE                                                   
                                                                        
#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "parvars.h"

        INTEGER :: ns                   ! Loop variable                   
        INTEGER :: nt                   ! Loop variable                   
        INTEGER :: nc                   ! Loop variable                   

        CHARACTER(LEN=72)  :: cmessage  ! String for error message
                                                                        
!       Reading dry deposition velocities from module
!       Dataset is for 6 time points and 5 land categories.                            

        IF (SIZE(depvel_defs) /= SIZE(depvel)) THEN 
          cmessage='sizes of depvel_defs and depvel are inconsistent'
! DEPENDS ON: ereport
          CALL EREPORT('UKCA_INDDEP',1,cmessage)
        ENDIF 

        DO ns=1,jpdd
          DO nc=1,jddepc
            DO nt=1,jddept
              depvel(nt,nc,ns)=depvel_defs(nt,nc,ns)
            ENDDO
          ENDDO
        ENDDO

        RETURN                                                          
        END SUBROUTINE UKCA_INDDEP 
#endif
