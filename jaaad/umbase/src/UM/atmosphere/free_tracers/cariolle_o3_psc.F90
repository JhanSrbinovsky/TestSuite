#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Routine to calculate ozone using cariolle scheme      
      
      Subroutine cariolle_o3_psc( OZONE_TRACER,                      &
                                  O3_PROD_LOSS,   O3_P_L_VMR,         &
                                  O3_VMR,         O3_P_L_TEMP,        &
                                  O3_TEMP,        O3_P_L_COLO3,       &
                                  O3_COLO3,                           &
                                  theta,                              &
                                  p_theta_levels,                     &
                                  off_x,off_y, theta_off_size,        &
                                  rows,row_length,                    &
                                  dt,                                 &
                                  exner_theta_levels,                 &
                                  model_levels)       
                                                            
! Purpose:
!        This routine evaluates the cariolle ozone                  
! tracer   
!
! Method:    on input - the actual ozone distribution       
!            on output - the updated ozone distribution   
! p_theta_levels   : pressure on theta levels.                        
! row_length,rows  : number of elements on 1 pe (x and y dimensions)   
!                    (all on one pressure level)                            
! k                : level index, counting from 1 to model_levels            
!                    1        -> surface (997 hPa)                      
!                    model_levels->top(e.g. model_levels=64 -> 0.01hPa)
! dt               : local advection timestep                   
!                    DYN_TIMESTEP/NSWEEPS 
!                    (NSWEEPS=2 if max wind or divergence 
!                    check is switched on and any values greater than
!                    the specified values occur somewhere in the model)                             
                                                            
        Implicit none                                    
                                                     
#include "cprintst.h"                                          

        integer :: k , rows, row_length        
        integer :: off_x,off_y,theta_off_size                           
        integer :: model_levels                 

        real, parameter :: mrconvert =1.6551                         

! The ozone fields
        real :: OZONE_TRACER(1-off_x:row_length+off_x, 1-off_y:rows+off_y,&
                   model_levels)                                      

        real :: O3_PROD_LOSS(rows,model_levels)         
        real :: O3_P_L_VMR  (rows,model_levels)         
        real :: O3_VMR      (rows,model_levels)           
        real :: O3_P_L_TEMP (rows,model_levels)        
        real :: O3_TEMP     (rows,model_levels)         
        real :: O3_P_L_COLO3(rows,model_levels)        
        real :: O3_COLO3    (rows,model_levels)

        real :: O3_PROD_LOSS_3D(row_length,rows,model_levels)     
        real :: O3_P_L_VMR_3D  (row_length,rows,model_levels)       
        real :: O3_VMR_3D      (row_length,rows,model_levels)        
        real :: O3_P_L_TEMP_3D (row_length,rows,model_levels)      
        real :: O3_TEMP_3D     (row_length,rows,model_levels)     
        real :: O3_P_L_COLO3_3D(row_length,rows,model_levels)        
        real :: O3_COLO3_3D    (row_length,rows,model_levels)

        real :: theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,    &
                   model_levels)                                       
        real :: dt                                                 
        real :: p_theta_levels(1-off_x:row_length+off_x,               &
                   1-off_y:rows+off_y,model_levels)           
        real :: exner_theta_levels(1-off_x:row_length+off_x,           &
                   1-off_y:rows+off_y,model_levels)          
                                                        
#include "c_r_cp.h"                        
                                            
! defining auxiliary variables              
! temp : temperature (derived from theta and pstar)           
! lat  : latitudes for the cariolle parameters           
! dso3 : column integrated o3 differences             
! cosz : cos(solar zenith angle)                      
! dp   : weight for the dso3 integration   
                         
        integer :: i, j  
        real    :: temp                                         
        real    :: dp, level_dso3, to3     
        real    :: dso3(row_length,rows)

!       A Geer: need better way? Maybe not trapezium after all
        real :: dtracer_abv(row_length,rows)                   
                                              
        integer :: my_pe                                         

       If (PrintStatus  >=  PrStatus_Normal) THEN                       
          WRITE(6,*) 'Calculating Cariolle Ozone tracer'
       End If

! Expand the cariolle parameters to 3-dimensions

       Do k = 1, model_levels
         Do j = 1, rows
           Do i = 1, row_length
              
              O3_PROD_LOSS_3D(i,j,k)   = O3_PROD_LOSS(j,k)
              O3_P_L_VMR_3D(i,j,k)     = O3_P_L_VMR(j,k)
              O3_VMR_3D(i,j,k)         = O3_VMR(j,k)
              O3_P_L_TEMP_3D(i,j,k)    = O3_P_L_TEMP(j,k)
              O3_TEMP_3D(i,j,k)        = O3_TEMP(j,k)
              O3_P_L_COLO3_3D(i,j,k)   = O3_P_L_COLO3(j,k)
              O3_COLO3_3D(i,j,k)       = O3_COLO3(j,k)

           End Do
         End Do
       End Do


! loop over all field elements                          
       Do k=model_levels,1,-1   ! Loop going down from model levels to 1                  
         Do j=1,rows                      
           Do i=1, row_length      
! convert mixing ratio (in kgkg-1) to vmr               
            OZONE_TRACER(i,j,k) = OZONE_TRACER(i,j,k)/mrconvert   

! convert theta to temp                                  
            temp=theta(i,j,k) * exner_theta_levels(i,j,k)             
                                                              
! Calculate total column ozone differences above this level          
            If (k .eq. model_levels) then                            
! reset dso3 (vertical column o3 differences above)                
! for k=model_levels (top of model)                             
              dso3(i,j)=0                              
            Else                                       
! trapezium rule integration                        
              dp=p_theta_levels(i,j,k)-p_theta_levels(i,j,k+1)        
              to3=(OZONE_TRACER(i,j,k)-O3_VMR_3D(i,j,k))+dtracer_abv(i,j)
              level_dso3 = (dp*to3/2.)*2.152e25/101300.             
              dso3(i,j)=dso3(i,j)+level_dso3                      
            End If                                                
                                                   
! A. Geer hack to store tracer values for trapezium rule integration    
! because the layer above has been modified by this scheme already     
! and has been converted back to kgkg-1                                 

            dtracer_abv(i,j) = OZONE_TRACER(i,j,k)-O3_VMR_3D(i,j,k)     

            If (k==1) then

              OZONE_TRACER(i,j,k) = OZONE_TRACER(i,j,k) +              &
                                 (dt*O3_PROD_LOSS_3D(i,j,k)) +         &
                                 (dt*(O3_P_L_TEMP_3D(i,j,k)*           &
                                     (temp-O3_TEMP_3D(i,j,k)))) 
         
            Else
              OZONE_TRACER(i,j,k) = OZONE_TRACER(i,j,k) +              &
                                 (dt*O3_PROD_LOSS_3D(i,j,k)) -         &
                                 (dt*(O3_VMR_3D(i,j,k)*                &
                                          O3_P_L_VMR_3D(i,j,k))) +     &
                                 (dt*(O3_P_L_TEMP_3D(i,j,k)*           &
                                     (temp-O3_TEMP_3D(i,j,k)))) +      &
                                 (dt*(O3_P_L_COLO3_3D(i,j,k)*dso3(i,j)))
          
            End if

       OZONE_TRACER(i,j,k)=OZONE_TRACER(i,j,k)/(1.-(dt*O3_P_L_VMR_3D(i,j,k)))

! sometimes negative values over Antarctica near the ground               
            If (OZONE_TRACER(i,j,k) .lt. 0) OZONE_TRACER(i,j,k)=0        
                                                                      
! convert mixing ratio (in ppmv) back to ppm                          
             OZONE_TRACER(i,j,k) = OZONE_TRACER(i,j,k)*mrconvert        
           End Do                                                     
         End Do                                                    
       End Do                                                     

        return                                                
        end                                               
! END cariolle_o3_psc
!-------------------------------------------------------------------
#endif                                              
