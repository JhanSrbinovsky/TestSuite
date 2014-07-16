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
! Purpose: Subroutine to read in coefficients for calculating           
!          effective Henry's Law coefficients. Original version         
!          taken from the Cambridge TOMCAT model.                       
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!                                                                       
!          Called from ASAD routine ASAD_cinit.                         
!                                                                       
! Current code owner: Colin Johnson/Olaf Morgenstern                   
!                                                                       
!                                                                       
! Code description:                                                     
!   Language: FORTRAN 90                                                
!   This code is written to UMDP3 v6 programming standards.             
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
        SUBROUTINE UKCA_INWDEP                                          

        USE ASAD_MOD,         ONLY: kd298, k298, ddhr, dhr
        USE UKCA_CHEM1_DAT
        IMPLICIT NONE                                                   

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"

        INTEGER            :: ns       ! Loop variable                      

        CHARACTER (LEN=72) :: cmessage ! String for error message

!       Use module to define Henry constants

        IF (SIZE(henry_defs) /= jpdw*4) THEN
          cmessage='jpdw and henry_defs are inconsistent'
! DEPENDS ON: EREPORT
          CALL EREPORT('UKCA_INWDEP',1,cmessage)
        ENDIF

        DO ns=1,jpdw
          k298(ns)  = henry_defs(1,ns)
          dhr(ns)   = henry_defs(2,ns)
          kd298(ns) = henry_defs(3,ns)
          ddhr(ns)  = henry_defs(4,ns)
        ENDDO
                                                                        
        RETURN                                                          
        END SUBROUTINE UKCA_INWDEP  
#endif
