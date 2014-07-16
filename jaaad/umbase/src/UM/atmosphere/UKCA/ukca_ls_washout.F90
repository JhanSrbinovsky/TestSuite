#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT******************************* 
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT******************************* 
!                                                                       
! Description:
!  To treat washout of UKCA dust tracers by LS precipitation.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Current Code Owner:       Colin Johnson
!
!
!  Code Description:
!   Language:  FORTRAN 90 (formatted)
!   This code is written to UMDP3 v6 programming standards.
!
! ######################################################################
!                                                                       
! Subroutine Interface:                                                 
      SUBROUTINE UKCA_LS_WASHOUT(                                      &       
        row_length, rows, model_levels, wet_model_levels,              &       
        r_rho_levels, r_theta_levels, tr_levels, halo_i, halo_j,       &       
        rho_r2, timestep, q, qcl, qcf, ls_rain3d, ls_snow3d,           &       
        dust_tracers,                                                  &       
#include "argsts.h"                                              
        STASHwork                                                      &       
        )                                                               
                                                                        
#include "c_dust_ndiv.h"                                    
#include "c_dustscav.h"                                      
                                                                        
      INTEGER, INTENT(IN) :: row_length                                 
      INTEGER, INTENT(IN) :: rows                                       
      INTEGER, INTENT(IN) :: model_levels                               
      INTEGER, INTENT(IN) :: wet_model_levels                           
      INTEGER, INTENT(IN) :: tr_levels                                  
      INTEGER, INTENT(IN) :: halo_i                                     
      INTEGER, INTENT(IN) :: halo_j                                     
                                                                        
      REAL, INTENT(IN) :: timestep  ! model timestep                    
      REAL, DIMENSION(row_length,rows,model_levels), INTENT(IN) ::     &   
            rho_r2                  ! density*r*r                       
      REAL, DIMENSION(row_length,rows,model_levels), INTENT(IN) ::     &   
            q                       ! specific humidity                 
      REAL, DIMENSION(row_length,rows,wet_model_levels), INTENT(IN) :: & 
            qcl                     ! cloud liquid water                
      REAL, DIMENSION(row_length,rows,wet_model_levels), INTENT(IN) :: & 
            qcf                     ! cloud frozen water                
      REAL, DIMENSION(row_length,rows,wet_model_levels), INTENT(IN) :: & 
            ls_rain3d               ! 3-D rain                          
      REAL, DIMENSION(row_length,rows,wet_model_levels), INTENT(IN) :: & 
            ls_snow3d               ! 3-D snow                          
      REAL, DIMENSION(1-halo_i:row_length+halo_i,                      & 
                         1-halo_j:rows+halo_j,0:model_levels),         & 
        INTENT(IN) :: r_theta_levels                                    
      REAL, DIMENSION(1-halo_i:row_length+halo_i,                      & 
                         1-halo_j:rows+halo_j,model_levels),           & 
        INTENT(IN) :: r_rho_levels                                      
      REAL, DIMENSION(row_length,rows,tr_levels,ndiv),                 & 
             INTENT(INOUT) :: dust_tracers   ! all dust tracers         
      REAL, DIMENSION(:), INTENT(INOUT) :: STASHwork                    
                                                                        
#include "typsts.h"                                              
#include "csubmodl.h"                                          
                                                                        
! Local variables                                                       
      REAL, DIMENSION(row_length,rows,wet_model_levels) ::             &        
                                                      lscav_dust_all    
! mass air per unit area in level                                       
      REAL, DIMENSION(row_length,rows,wet_model_levels) :: dm           
      REAL, DIMENSION(row_length,rows)                  :: rainrate     
      REAL, DIMENSION(row_length,rows)                  :: snowrate     
      REAL, DIMENSION(row_length,rows)                  :: delta_dust   
      REAL, DIMENSION(row_length,rows)                  :: rate         
      INTEGER :: IDIV,K                                                 
      INTEGER :: errorcode                                              
                                                                        
! Calculate mass of (dry) air per square metre   
! DEPENDS ON: mass_calc
      CALL MASS_CALC(                                                  &
        row_length, rows, model_levels, wet_model_levels,              &
        r_rho_levels(1:row_length,1:rows,1:model_levels),              &
        r_theta_levels(1:row_length,1:rows,0:model_levels),            &
        timestep, rho_r2(1:row_length,1:rows,1:model_levels),          &
        q, qcl, qcf,                                                   &
        dm )                                                            
                                                                        
      lscav_dust_all = 0.                                               
                                                                        
      DO K=wet_model_levels,1,-1                                        
! Deal with possible -ive precipn.                                      
        rainrate(:,:)=MAX(ls_rain3d(:,:,K),0.)                          
        snowrate(:,:)=MAX(ls_snow3d(:,:,K),0.)                          
! Calculate proportion of dust mixing ratio scavenged                   
        DO IDIV=1,ndiv                                                  
          krain=krain_dust(IDIV)                                        
          ksnow=ksnow_dust(IDIV)                                        
          rate=(rainrate*krain+snowrate*ksnow)*3600.0*timestep          
          delta_dust(:,:)=dust_tracers(:,:,K,IDIV)*                    &               
                    (1.0-1.0/(1.0+rate(:,:)))                           
! Calculate mass of dust removed                                        
          lscav_dust_all(:,:,IDIV)=lscav_dust_all(:,:,IDIV)+           &        
                    delta_dust(:,:)*dm(:,:,K)                           
! Decrement mixing ratio                                                
          dust_tracers(:,:,K,IDIV)=dust_tracers(:,:,K,IDIV)-           &         
                    delta_dust(:,:)                                     
        ENDDO ! IDIV                                                    
      ENDDO ! K                                                         
                                                                        
! Copy diagnostics into stashwork array                                 
                                                                        
        DO IDIV = 1,ndiv                                                
          IF (SF(230+K,4)) THEN                                         
! Convert "per timestep" diagnostic to "per second":                    
          lscav_dust_all(:,:,IDIV)=lscav_dust_all(:,:,IDIV)/timestep    

!! DEPENDS ON: copydiag
!            CALL COPYDIAG(stashwork(SI(230+IDIV,4,IM_INDEX)),          
!               lscav_dust_all,                                         
!               row_length,rows,0,0,0,0, at_extremity,                  
!               atmos_im,4,230+IDIV,                                    
!               errorcode,cmessage)                                     
                                                                        
!            IF (errorcode .GT. 0) THEN                                  
!              cmessage=": ERROR IN COPYDIAG(LSCAV_DUST_ALL)"            
!! DEPENDS ON: ereport
!              CALL EREPORT('UKCA_LS_WASH',errorcode,cmessage)           
!            ENDIF                                                       
          ENDIF                                                         
                                                                        
        ENDDO !NDIV                                                     
                                                                        
      END SUBROUTINE UKCA_LS_WASHOUT   
#endif
