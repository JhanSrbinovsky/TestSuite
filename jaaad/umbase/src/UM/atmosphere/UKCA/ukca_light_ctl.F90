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
! Purpose: Subroutine to call the UKCA lightning emissions routine,     
!          which diagnoses NOx lightning emissions and returns values   
!          to UKCA_EMISSION_CTL.                                        
!          Based on routine provided by Olaf Morgenstern.               
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!                                                                       
!          Called from UKCA_EMISSION_CTL.                               
!                                                                       
! Current code owner: Colin Johnson/Olaf Morgenstern/Fiona O'Connor                   
!                                                                       
!                                                                       
! Code description:                                                     
!   Language: FORTRAN 90                                                
!   This code is written to UMDP3 v6 programming standards.             
!                                                                       
! --------------------------------------------------------------------- 
        SUBROUTINE UKCA_LIGHT_CTL(                                     &           
              rows, row_length, p_levels, conv_base_lev,               &           
              conv_top_lev, mask, lat, area,                           &           
              r_theta_levels, r_rho_levels,                            &           
              p_centres, p_boundaries,                                 &           
              an)                                                      
        IMPLICIT NONE                                                   
                                                                        
#include "c_a.h"                                                    
#include "c_g.h"                                                    
#include "c_pi.h"          
#include "c_v_m.h"
                                                                        
        INTEGER, INTENT(IN) :: rows                           ! No of la
        INTEGER, INTENT(IN) :: row_length                     ! No of lo
        INTEGER, INTENT(IN) :: p_levels                       ! No of ve
        INTEGER, INTENT(IN) :: conv_base_lev(row_length,rows) ! Levels o
        INTEGER, INTENT(IN) :: conv_top_lev(row_length,rows)  ! Levels o
        INTEGER, INTENT(IN) :: mask(row_length, rows)         ! land(1)/
                                                                        
        REAL, INTENT(IN) :: lat(row_length, rows)                       
        REAL, INTENT(IN) :: area(row_length, rows)
        REAL, INTENT(IN) :: r_theta_levels(row_length, rows, 0:p_levels)
        REAL, INTENT(IN) :: r_rho_levels(row_length, rows, p_levels)    
        REAL, INTENT(IN) :: p_centres(row_length, rows, p_levels)       
        REAL, INTENT(IN) :: p_boundaries(row_length, rows, 0:p_levels)  
                                                                        
        REAL, INTENT(OUT) :: an(row_length, rows, p_levels)             
                                                                        
! Local variables                                                       
                                                                        
        INTEGER :: i                                     ! Loop variable
        INTEGER :: j                                     ! Loop variable
        INTEGER :: k                                     ! Loop variable
        INTEGER :: klc                                   ! Cloud base
        INTEGER :: klt                                   ! Cloud top
        INTEGER :: asurf                                 ! Land/sea mask               
                                                                        
        REAL :: model_half_ht(row_length, rows, 0:p_levels)            
        REAL :: height(row_length, rows, p_levels)                      
        REAL :: cloud_base(row_length, rows)                            
        REAL :: cloud_top(row_length, rows)                             
        REAL :: delp_p(row_length, rows)                                
        REAL :: mass(row_length, rows)                                  
        REAL :: anox(p_levels)                                          
        REAL :: press(p_levels)                                         
        REAL :: hkmb                                                    
        REAL :: hkmt                                                    
        REAL :: asfaera                                                 
        REAL :: adlat                                                   
                                                                        
!       Find heights of model half levels                                     
                                                                        
        DO k = 0,p_levels                                               
          model_half_ht(:,:,k) = r_theta_levels(:,:,k)                 &          
                               - r_theta_levels(:,:,0)                  
        ENDDO                                                           
                                                                        
!       Find heights of model levels                                          
                                                                        
        DO k = 1,p_levels                                               
          height(:,:,k) = r_rho_levels(:,:,k) - r_theta_levels(:,:,0)   
        ENDDO                                                           
                                                                        
!       Convective cloud base and top height in km                            
                                                                        
        cloud_base = 0.0                                                
        cloud_top  = 0.0                                                
                                                                        
        do i = 1,rows                                                   
          do j = 1,row_length                                           
            IF (conv_base_lev(j,i) > 0) THEN                            
              cloud_base(j,i) = height(j,i,conv_base_lev(j,i))/1000.0   
            ENDIF                                                       
            IF (conv_top_lev(j,i) > 0) THEN                             
              cloud_top(j,i) = height(j,i,conv_top_lev(j,i))/1000.0     
            ENDIF                                                       
          ENDDO                                                         
        ENDDO                                                           
                                                                        
!       Initialize AN array                                                   
                                                                        
        an = 0.0                                                        
                                                                        
        DO i = 1,rows                                                   
          DO j = 1,row_length                                           
            hkmb      = cloud_base(j,i)                                 
            hkmt      = cloud_top(j,i)                                  
            adlat     = lat(j,i)                                        
            asurf     = mask(j,i)                                       
            klc       = conv_base_lev(j,i)                              
            klt       = conv_top_lev(j,i)                               
            asfaera   = Earth_Radius*Earth_Radius*area(j,i)             
            anox      = 0.0                                             
            press     = p_centres(j,i,:)                                

            IF (klc > 1 .AND. klt < p_levels) THEN   
!             Call UKCA_LIGHT to calculate NOx emissions (kg N/grid box/s)
! DEPENDS ON: ukca_light
              CALL UKCA_LIGHT(press,p_levels,hkmb,hkmt,klc,klt,adlat   & 
                   ,asfaera,asurf,anox)                                 
              an(j,i,:) = anox                                          
            ENDIF                                                       
          ENDDO                                                         
        ENDDO                                                           
                                                                        
!       Convert NOx emissions to kg(NO2)/kg(air)/s                            
                                                                        
        DO k = 1,p_levels                                               
          delp_p = p_boundaries(:,:,k-1) - p_boundaries(:,:,k)          
!         Calculate mass of each grid box                               
          mass = Earth_Radius*Earth_Radius*area*delp_p/g                
!         NOx emissions in kg(NO2)/kg(air)/s                            
          an(:,:,k) = an(:,:,k)*c_no2/c_n/mass                            
        ENDDO                                                           
                                                                        
        RETURN                                                          
        END SUBROUTINE UKCA_LIGHT_CTL                                   
#endif                                                                  
