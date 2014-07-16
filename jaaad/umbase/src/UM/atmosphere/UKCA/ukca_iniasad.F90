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
! Purpose: Subroutine to initialise ASAD and fill ldepd and ldepw arrays
!          Adapted from original version written by Olaf Morgenstern.   
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from UKCA_MAIN1.                                      
!                                                                       
! Current code owner: Colin Johnson/Olaf Morgenstern
!                     Fiona O'Connor                   
!                                                                       
!                                                                       
! Code description:                                                     
!   Language: FORTRAN 90                                                
!   This code is written to UMDP3 v6 programming standards.             
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
        SUBROUTINE UKCA_INIASAD(npoints)                                

        USE UKCA_CHEM1_DAT
        USE ASAD_MOD
        IMPLICIT NONE                                                   
                                                                        
        INTEGER, INTENT(IN) :: npoints   ! no of spatial points         

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"

!       Local variables                                                 
                                                                        
        INTEGER :: k                    ! Loop variable                 
                                                                        
        CHARACTER (LEN=72) :: cmessage  ! Error message                 

        write(6,*) 'ASAD initialised from namelist using:'
        write(6,*) 'jpctr: ',jpctr,' jpspec: ',jpspec,' jpnr: ',jpnr
        write(6,*) 'jpbk: ',jpbk,' jptk: ',jptk, ' jppj: ',jppj
        write(6,*) 'jphk: ',jphk,' jpdd: ',jpdw, ' jpeq: ',jpeq

        CALL ASAD_MOD_INIT

!       Set up dry and wet deposition logicals using module info

        IF (size(chch_defs) /= jpspec) THEN 
          cmessage='size of chch_defs inconsistent with jpspec'
! DEPENDS ON: ereport
          CALL EREPORT('UKCA_INIASAD',1,cmessage)
        ENDIF

        DO k=1,jpspec
          ldepd(k) = (chch_defs(k)%switch1 == 1)
          ldepw(k) = (chch_defs(k)%switch2 == 1)
        ENDDO

! DEPENDS ON: asad_cinit
        CALL ASAD_CINIT(npoints)                                             
                                                                        
        RETURN                                                          
        END SUBROUTINE UKCA_INIASAD                                    
#endif
