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
! Purpose: Subroutine to calculate NOx lightning emissions for one      
!          vertical column based on cloud height above surface,         
!          latitude & land/ocean surface.                               
!                                                                       
!          Based on light.F from Cambridge TOMCAT model                 
!          written by Zoe Stockwell and adapted by Guang Zeng.          
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!                                                                       
! Method:  Cloud height above surface and surface type yield lightning  
!          flash frequency. Latitude yields CG/CC ratio (Price & Rind   
!          1993), assuming that all dH (zero degrees to cloud top) lie  
!          between 5.5-14km. Convert flashes to NOx production for the  
!          1/2 hr dynamical period.                                     
!                                                                       
!          Called from UKCA_LIGHT_CTL.                                  
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
        SUBROUTINE UKCA_LIGHT(ppress,niv,hkmb,hkmt,klc,klt,            &      
                              adlat,asfaera,asurf,anox)                 
        IMPLICIT NONE                                                   
                                                                        
#include "c_sulchm.h"                                          
                                                                        
        INTEGER, INTENT(IN) :: niv                ! No of vertical level
        INTEGER, INTENT(IN) :: klc                ! Level of cloud base 
        INTEGER, INTENT(IN) :: klt                ! Level of cloud top  
        INTEGER, INTENT(IN) :: asurf              ! Land (1) / sea (0) m
                                                  ! Assume less lightnin
                                                                        
        REAL, INTENT(IN) ::  hkmt                 ! Height of cloud top 
        REAL, INTENT(IN) ::  hkmb                 ! Height of cloud base
        REAL, INTENT(IN) ::  adlat                ! Latitude            
        REAL, INTENT(IN) ::  ppress(niv)          ! Pressures at model l
        REAL, INTENT(IN) ::  asfaera              ! Surf area * (radius 
                                                                        
        REAL, INTENT(OUT) :: anox(niv)            ! NOx lightning emissi
                                                                        
! Local variables                                                       
                                                                        
        INTEGER :: jniv                           ! Loop variable       
        INTEGER :: k                              ! Loop variable       
                                                                        
        REAL, PARAMETER ::  Min_clouddepth = 5.0  ! Minimum cloud depth 
        REAL, PARAMETER ::  convfac = 0.001158    ! Conversion factor to
        REAL, PARAMETER ::  N_nitrogen = 1.0e+26  ! no of molecules of N
                                                  ! (Noxon 1976, & proba
        REAL, PARAMETER ::  M_nitrogen = 14.01    ! Molecular mass of N 
                                                                        
        REAL :: aflash                 ! Flash frequency (flashes/min)  
        REAL :: adh                                                     
        REAL :: az                     ! Cloud-cloud/cloud-ground       
        REAL :: ap                     ! Cloud-ground flashes / total fl
        REAL :: acgfrac                ! Cloud-ground flash frequency (f
        REAL :: accfrac                ! Cloud-cloud flash frequency (fl
        REAL :: acgnox                                                  
        REAL :: accnox                                                  
        REAL :: dpcg                                                    
        REAL :: dpic                                                    
                                                                        
!       Initialise variables                                            
                                                                        
        aflash  = 0.0                                                   
        adh     = 0.0                                                   
        az      = 0.0                                                   
        ap      = 0.0                                                   
        acgfrac = 0.0                                                   
        accfrac = 0.0                                                   
        acgnox  = 0.0                                                   
        accnox  = 0.0                                                   
                                                                        
!       Set minimum cloud depth to 5 km                                 
                                                                        
        IF ((hkmt-hkmb) > Min_clouddepth) THEN                          
                                                                        
          IF (asurf == 0) THEN                  ! Ocean                 
            aflash = 6.40e-04*(hkmt**1.73)      ! Ocean flash frequency 
          ELSE                                  ! Land                  
            aflash = 3.44e-05*(hkmt**4.9)       ! Land flash frequency  
          ENDIF                                                         
                                                                        
!         Calculate flash rate in flashes/s/gridbox(km^2)               
          aflash = aflash*asfaera*1.0e-6/25.0/60.0                      
                                                                        
!         Work out proportion of flashes that are cloud to ground       
                                                                        
          adh = -6.64e-05*(adlat**2)-4.73e-03*adlat+7.34                
          az  = 0.021*(adh**4)-0.648*(adh**3)+7.493*(adh**2)           &     
                              -36.54*adh+63.09                          
          ap  = 1./(1.+az)                                              
          acgfrac = ap*aflash                                           
          accfrac = aflash-acgfrac                                      
                                                                        
!         Convert flash frequency to N (kg) per second.                 
                                                                        
          acgnox = acgfrac*N_nitrogen/AVOGADRO*M_nitrogen*1.0e-03       
          accnox = accfrac*N_nitrogen/AVOGADRO*M_nitrogen*1.0e-03*0.1   
                                                                        
!         Distribute over the column with each box having same vmr in an
!         multiply by CONVFAC, conversion factor that gives emissions if
!         were 100 flashes s**-1                                        
                                                                        
!         Work out which pressure is closest to 500 hPa                 
                                                                        
          LOOP: DO jniv = niv,1,-1                                      
            IF (ppress(jniv) >= 50000.0) EXIT LOOP                      
          ENDDO LOOP                                                    
                                                                        
!         3 from cloud base to 2 above top                              
!         KLT is the level above cloud top                              
                                                                        
          dpcg = ppress(1) - ppress(jniv-1)                             
          dpic = ppress(jniv) - ppress(klt)                             
                                                                        
          IF ((jniv-1) == 1) THEN                                       
            anox(1) = acgnox * convfac                                  
          ELSE                                                          
            DO k = 1,jniv-1                                             
              anox(k) = acgnox*(ppress(k)-ppress(k+1))                 &             
                        /dpcg * convfac                                 
            END DO                                                      
          ENDIF                                                         
                                                                        
          IF (ppress(jniv) <= ppress(klt)) THEN                         
            anox(klt-1) = anox(klt-1)+accnox*convfac                    
          ELSE                                                          
            DO k = jniv,klt-1                                           
              anox(k) = accnox*(ppress(k)-ppress(k+1))                 &     
                        /dpic * convfac                                 
            END DO                                                      
          ENDIF                                                         
                                                                        
        ENDIF                                                           
                                                                        
        RETURN                                                          
        END SUBROUTINE UKCA_LIGHT                                       
#endif                                                                  
