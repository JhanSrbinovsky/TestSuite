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
! Purpose: Subroutine to calculate wet deposition rates in s-1          
!          from convective and dynamic rainfall.                        
!
!          Based on routine from Cambridge TOMCAT model.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!                                                                       
!          Called from UKCA_CHEMISTRY_CTL.                              
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
        SUBROUTINE UKCA_WDEPRT(drain,crain,p_fieldda,q_levelsda,t,     &    
                              sinlat,tstep,                            &    
                              first_point,last_point,wetrt)             

        USE ASAD_MOD,        ONLY: ddhr, dhr, kd298, k298
        IMPLICIT NONE                                                   
                                                                        
#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "c_pi.h"
                                                                        
        INTEGER, INTENT(IN) :: p_fieldda            ! No of spatial pts
        INTEGER, INTENT(IN) :: q_levelsda           ! No of vertical levs
        INTEGER, INTENT(IN) :: first_point          ! First spatial pt
        INTEGER, INTENT(IN) :: last_point           ! Last spatial pt
                                                                        
        REAL, INTENT(IN) :: drain(p_fieldda,q_levelsda)                 
        REAL, INTENT(IN) :: crain(p_fieldda,q_levelsda)                 
        REAL, INTENT(IN) :: t(p_fieldda,q_levelsda)   ! Temperature     
        REAL, INTENT(IN) :: sinlat(p_fieldda)         ! Sine(latitude)  
        REAL, INTENT(IN) :: tstep                     ! Timestep        
                                                                        
        REAL, INTENT(OUT) :: wetrt(p_fieldda,q_levelsda,jpdw)           
                                                                        
!       Local variables                                                 
                                                                        
        INTEGER :: i                    ! Loop variable                 
        INTEGER :: k                    ! Loop variable                 
        INTEGER :: ns                   ! Loop variable                 
                                                                       
        REAL, PARAMETER :: fac   = 0.1          ! Factors to convert rainfall rate
        REAL, PARAMETER :: fcc   = 86400.E6     ! to cm/s               
        REAL, PARAMETER :: fc    = 0.3          ! Fraction of gridbox with conv rain
        REAL, PARAMETER :: csca  = 4.7          ! Conv rain scavenging rate
        REAL, PARAMETER :: dsca  = 2.4          ! Dyn rain scavenging rate
        REAL, PARAMETER :: clw   = 1.0E-6       ! Cloud liquid water concn
        REAL, PARAMETER :: Rgas  = 0.08314      ! Ideal gas constant  
        REAL, PARAMETER :: Nlat  = 65.0         ! Northern limit for scavenging
        REAL, PARAMETER :: Slat  = -65.0        ! Southern limit for scavenging
        REAL, PARAMETER :: Tmax  = 273.15       ! Max temp used to limit scav
        REAL, PARAMETER :: Tmin  = 253.15       ! Min temp used to limit scav
                                                                        
        REAL :: tmp1                              ! Temporary Variable 
        REAL :: tmp2                              ! Temporary Variable 
        REAL :: kaq                               ! Variable used to calculate faq 
        REAL :: rcw                               ! Conv scav rate 
        REAL :: rdw                               ! Dyn scav rate 
        REAL :: flim                              ! Limit to dyn scav 
        REAL :: rwcon(p_fieldda,q_levelsda,jpdw)  ! Gridbox conv scav rate
        REAL :: rwdyn(p_fieldda,q_levelsda,jpdw)  ! Gridbox dyn scav rate
        REAL :: hcoef(p_fieldda,q_levelsda,jpdw)  ! Effective Henry's coeff
        REAL :: faq(p_fieldda,q_levelsda,jpdw)    ! Fraction of species in aq phase
        REAL :: lat(p_fieldda)                    ! Latitudes           
        REAL :: ddrain(p_fieldda,q_levelsda)      ! Accumulated dynamic rain
        REAL :: ccrain(p_fieldda,q_levelsda)      ! Accumulated convective rain
                                                                        
!       Initialise arrays                                               
                                                                        
        DO ns = 1,jpdw                                                  
          DO k=1,q_levelsda                                             
            DO i=1,p_fieldda                                            
              wetrt(i,k,jpdw) = 0.0                                     
              rwdyn(i,k,jpdw) = 0.0                                     
              rwcon(i,k,jpdw) = 0.0                                     
            ENDDO                                                       
          ENDDO                                                         
        ENDDO                                                           
                                                                        
        ddrain = drain                                                  
        ccrain = crain                                                         
                                                                       
!       Add rain to the level below                                     
                                                                        
        DO i = first_point,last_point                                   
          DO k = q_levelsda,2                                           
            ccrain(i,k-1) = ccrain(i,k-1)+ccrain(i,k)                   
            ddrain(i,k-1) = ddrain(i,k-1)+ddrain(i,k)                   
          ENDDO                                                         
        ENDDO                                                           
                                                                        
        DO ns = 1,jpdw                                                  
          DO k = 1,q_levelsda                                           
            DO i = first_point,last_point                               
                                                                        
!             Calculate effective Henry's law coefficient taking into account
!             the effects of dissociation and complex formation upon 
!             solubility. pH of rain takes 5.0                          
!             see --- Christos thesis (1998) Eqs. 5.6-5.8 (pp.65)       
                                                                        
              tmp1 = (298.0 - t(i,k))/t(i,k)/298.0                      
              kaq  = kd298(ns)*exp(ddhr(ns)*tmp1)                       
              hcoef(i,k,ns) = k298(ns)*exp(dhr(ns)*tmp1)               &           
                                     *(1.0 + kaq/1.E-5)                 
                                                                        
!             Calculate the fraction of the tracer existing in the liquid phase
                                                                        
              tmp2 = clw*hcoef(i,k,ns)*Rgas*t(i,k)                      
              faq(i,k,ns) = tmp2/(1.0 + tmp2)                           
                                                                        
!             Calculate convective scavenging rate                      
                                                                        
              rcw = ccrain(i,k)*fac*csca*faq(i,k,ns)                    
                                                                        
!             Compute new scavenging rates                              
                                                                        
              rcw = rcw / fc                                            
                                                                        
!             Find effective scavenging rate for convective rain (sec-1)
                                                                        
              rwcon(i,k,ns) = (-alog(1-fc+fc*exp(-tstep*rcw)))/tstep    
                                                                        
            ENDDO                                                       
          ENDDO                                                         
        ENDDO                                                           
                                                                        
!       Latitude in degrees                                             
                                                                        
        DO i = first_point,last_point                                   
          lat(i) = asin(sinlat(i))*Recip_Pi_Over_180                    
        ENDDO                                                           
                                                                        
        DO ns = 1,jpdw                                                  
          DO k = 1,q_levelsda                                           
            DO i = first_point,last_point                               
                                                                        
!             For dynamic cloud assume cloud cover=1                    
!             Calculate dynamical scavenging rate                       
                                                                        
              rdw = ddrain(i,k)*fac*dsca*faq(i,k,ns)                    
                                                                        
!             Place limit to scavenging in polar regions                
                                                                        
              IF(lat(i) >= Nlat .OR. lat(i) <= Slat) THEN               
                                                                        
                IF (t(i,k) >= Tmax) THEN                                
                  flim = 1.0                                            
                ELSE IF(t(i,k) < Tmin) THEN                             
                  flim = 0.0  
                ELSE   
                  flim = 1.0 + 0.05*(t(i,k)-Tmax)                                                             
                END IF                                                  
                                                                        
                rdw = rdw * flim                                        
              ENDIF                                                     
                                                                        
!             Wet deposition rate for dynamical rain only(sec-1)        
                                                                        
              rwdyn(i,k,ns) = rdw                                       
                                                                        
            ENDDO                                                       
          ENDDO                                                         
        ENDDO                                                           
                                                                        
!       Wet deposition rate : first-order loss by                       
!       dynamic rain + convective rain                                  
                                                                        
        DO ns = 1,jpdw                                                  
          DO k = 1,q_levelsda                                           
            DO i = first_point,last_point                               
              wetrt(i,k,ns) = rwcon(i,k,ns) + rwdyn(i,k,ns)             
            ENDDO                                                       
          ENDDO                                                         
        ENDDO                                                           
                                                                        
        RETURN                                                          
        END SUBROUTINE UKCA_WDEPRT                                      
#endif                                                                  
                                                                        
                                                                        
