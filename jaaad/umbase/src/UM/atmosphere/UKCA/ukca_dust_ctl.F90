#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT******************************* 
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT******************************* 
!                                                                       
! Description:
!   Subroutine to call Woodward dust scheme from UKCA
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!   Called from UKCA_MAIN1.
!
! Current code owner: Colin Johnson
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
! Subroutine Interface:                                                 
      SUBROUTINE UKCA_DUST_CTL(timestep,                               &                          
        land_pts, p_rho_levels, p_theta_levels, pstar, tstar,          &      
        t_modellevs, r_theta_levels, r_rho_levels,                     &      
        rho_r2, q, qcl, qcf, ls_rain3d, ls_snow3d,                     &      
        dust_tracers,                                                  &      
        tile_pts,tile_index,tile_frac,fland,                           &      
        pstar_land,tstar_tile,rhostar_land,soil_layer_moisture,        &      
        snow_tile,ustar_tile,mrel_land,clay_land,                      &      
        z_half, alpha_cd, ml_depth,                                    &
        rhokh_mix, rho_aresist, aresist, resist_b,                     &
        dtrdz_charney_grid, kent, we_lim,                              &
        t_frac, zrzi, kent_dsc,                                        &
        we_lim_dsc, t_frac_dsc,                                        &
        zrzi_dsc, zhsc, rb_dust_ndivs,                                 &
#include "arglndm.h"                                            
#include "argsts.h"                                              
        STASHwork                                                      &         
        )                                                               
! Data fields IN                                                        
!                                                                       
! Data fields IN/OUT                                                    
!                                                                       
! Data fields OUT                                                       
!                                                                       
!-----------------------------------------------------------------------
! Purpose: Interface for Dust Modelling, to include surface emissions   
!                                                                       
! Level 2 control routine                                               
!                                                                       
! Current owners of code:                C Johnson                      
!                                                                       
! History:                                                              
! Version     Date     Comment                                          
! -------     ----     -------                                          
!                                                                       
! Code Description:                                                     
!  Language: FORTRAN90 retaining F77 formatting                         
!  This code is written to UMDP3 v6 programming standards               
!                                                                       
! System components covered:                                            
!                                                                       
! System task:                                                          
!                                                                       
! Documentation: Not yet available                                      
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      USE UKCA_D1_DEFS                                                  
      IMPLICIT NONE                                                     
                                                                        
#include "parvars.h"  
#include "c_dust_ndiv.h"
#include "c_dustgen.h"                                        
#include "csubmodl.h"                                          
#include "typsts.h"                                              
#include "nstypes.h"                                            
#include "typsize.h"                                            
#include "typlndm.h"                                            
                                                                        
      INTEGER, INTENT(IN) :: land_pts       ! no of land points         
      INTEGER, INTENT(IN), DIMENSION(ntype)          :: tile_pts        
      INTEGER, INTENT(IN), DIMENSION(land_pts,ntype) :: tile_index      
      INTEGER, INTENT(IN) :: kent(1:row_length, 1:rows)
      INTEGER, INTENT(IN) :: kent_dsc(1:row_length, 1:rows)

      REAL, INTENT(IN)                             :: timestep          
      REAL, INTENT(IN), DIMENSION(land_pts)        :: pstar_land        
      REAL, INTENT(IN), DIMENSION(land_pts)        :: rhostar_land      
      REAL, INTENT(IN), DIMENSION(land_pts)        :: clay_land         
      REAL, INTENT(IN), DIMENSION(land_pts)        :: fland             
      REAL, INTENT(IN), DIMENSION(land_pts,ndiv)   :: mrel_land         
      REAL, INTENT(IN), DIMENSION(land_pts,ntiles) :: tile_frac         
      REAL, INTENT(IN), DIMENSION(land_pts,ntiles) :: tstar_tile        
      REAL, INTENT(IN), DIMENSION(land_pts,ntiles) :: snow_tile         
      REAL, INTENT(IN), DIMENSION(land_pts,ntiles) :: ustar_tile        
      REAL, INTENT(IN), DIMENSION(land_pts,sm_levels)                  &              
                                            :: soil_layer_moisture      
      REAL, INTENT(IN), DIMENSION(row_length,rows) :: pstar             
      REAL, INTENT(IN), DIMENSION(row_length,rows) :: tstar             
      REAL, INTENT(IN), DIMENSION(row_length,rows,model_levels)        &    
                                                   :: t_modellevs       
      REAL, DIMENSION(1-offx:row_length+offx, 1-offy:rows+offy,        &    
                     model_levels), INTENT(IN)    :: p_rho_levels       
      REAL, DIMENSION(1-offx:row_length+offx, 1-offy:rows+offy,        &    
                     model_levels), INTENT(IN)    :: p_theta_levels     
      REAL, DIMENSION(1-offx:row_length+offx, 1-offy:rows+offy,        &    
                   0:model_levels), INTENT(IN)    :: r_theta_levels     
      REAL, DIMENSION(1-offx:row_length+offx, 1-offy:rows+offy,        &    
                     model_levels), INTENT(IN)    :: r_rho_levels       
      REAL, DIMENSION(row_length,rows,model_levels), INTENT(IN) ::     &    
            rho_r2                  ! density*r*r                       
      REAL, DIMENSION(row_length,rows,wet_levels), INTENT(IN) ::       &    
            q                       ! specific humidity                 
      REAL, DIMENSION(row_length,rows,wet_levels), INTENT(IN) ::       &    
            qcl                     ! cloud liquid water                
      REAL, DIMENSION(row_length,rows,wet_levels), INTENT(IN) ::       &    
            qcf                     ! cloud frozen water                
      REAL, DIMENSION(row_length,rows,wet_levels), INTENT(IN) ::       &    
            ls_rain3d               ! 3-D rain                          
      REAL, DIMENSION(row_length,rows,wet_levels), INTENT(IN) ::       &    
            ls_snow3d               ! 3-D snow                          
      REAL, INTENT(IN) :: z_half(1:row_length,1:rows,1:bl_levels)
      REAL, INTENT(IN) :: alpha_cd(1:bl_levels)
      REAL, INTENT(IN) :: ml_depth(1:row_length,1:rows)
      REAL, INTENT(IN) :: rhokh_mix(1:row_length,1:rows,1:bl_levels)
      REAL, INTENT(IN) :: rho_aresist(1:row_length,1:rows)
      REAL, INTENT(IN) :: aresist(1:row_length,1:rows)
      REAL, INTENT(IN) :: resist_b(1:row_length,1:rows)
      REAL, INTENT(IN) :: dtrdz_charney_grid(1:row_length,1:rows,      &
                                               1:bl_levels)
      REAL, INTENT(IN) :: we_lim(1:row_length,1:rows,1:3)
      REAL, INTENT(IN) :: t_frac(1:row_length,1:rows,1:3)
      REAL, INTENT(IN) :: zrzi(1:row_length,1:rows,1:3)
      REAL, INTENT(IN) :: we_lim_dsc(1:row_length,1:rows,1:3)
      REAL, INTENT(IN) :: t_frac_dsc(1:row_length,1:rows,1:3)
      REAL, INTENT(IN) :: zrzi_dsc(1:row_length,1:rows,1:3)
      REAL, INTENT(IN) :: zhsc(1:row_length,1:rows)
      REAL, INTENT(IN) :: rb_dust_ndivs(1:row_length,1:rows,1:ndiv)

! Dust division 1-6:                                                    
      REAL, DIMENSION(1-offx:row_length+offx, 1-offy:rows+offy,        &    
            tr_levels,ndiv), INTENT(INOUT)        :: dust_tracers       
      REAL, DIMENSION(:), INTENT(INOUT) :: STASHwork                    
                                                                        
! Local variables                                                       
      INTEGER            :: I,J,K,L,M,N     ! loop counters             
      INTEGER            :: IDIV            ! loop counter              
      INTEGER            :: lev1=1          ! no. levs for vgrav calc   
      INTEGER            :: errorstatus     ! error code                
      INTEGER            :: sect=3          ! stash section             
      INTEGER            :: item            ! stash item                
                                                                        
      CHARACTER(LEN=72)  :: cmessage        ! for error reporting       
                                                                        
      REAL, DIMENSION(land_pts,ntiles,ndiv)                            &                      
              :: dust_flux_tile  ! dust flux (kg.m-2.s-2)               
      REAL, DIMENSION(land_pts,ntiles,ndiv)                            &  
              :: u_s_t_tile      ! threshold ustar on tiles             
      REAL, DIMENSION(land_pts,ntiles,ndiv)                            &  
              :: u_s_t_dry_tile  ! thresh. ustar on tiles excl. moisture
      REAL, DIMENSION(row_length,rows,0:model_levels)                  & 
              :: p_layer_boundaries
      REAL, DIMENSION(row_length,rows,0:model_levels)                  &
              :: p_layer_centres 
! Dust production flux (kg.m-2.s-2) 
      REAL, DIMENSION(row_length,rows,ndiv) :: dust_flux
! Threshold ustar 
      REAL, DIMENSION(row_length,rows,ndiv) :: u_s_t  
! Threshold ustar excluding moisture
      REAL, DIMENSION(row_length,rows,ndiv) :: u_s_t_dry
                                                                        
      REAL, DIMENSION(row_length,rows,model_levels) :: dustwork    
      REAL, DIMENSION(row_length,rows,bl_levels)    :: tr_flux 
      REAL, DIMENSION(row_length,rows)              :: dustwork2
      REAL, DIMENSION(row_length,rows)              :: dustwork3  
      REAL, DIMENSION(row_length,rows)              :: res_factor 
      REAL, DIMENSION(row_length,rows)              :: vstokes1    
      REAL, DIMENSION(row_length,rows)              :: drydep_str 
                                                                        
! set p at layer boundaries.   
! NB: pstar has no halo, p_theta and p_rho arrays have a halo    
      DO j = 1, rows  
        DO i = 1, row_length 
          p_layer_boundaries(i,j,0) = pstar(i,j)     
          p_layer_centres(i,j,0) = pstar(i,j)   
        ENDDO   
      ENDDO
      DO k = 1, model_levels - 1   
        DO j = 1, rows   
          DO i = 1, row_length 
            p_layer_boundaries(i,j,k) = p_rho_levels(i,j,k+1)  
            p_layer_centres(i,j,k) = p_theta_levels(i,j,k)   
          ENDDO   
        ENDDO     
      ENDDO
      DO j = 1, rows 
        DO i = 1, row_length  
          p_layer_boundaries(i,j,model_levels) = 0.0  
          p_layer_centres(i,j,model_levels) =                          & 
                                 p_theta_levels(i,j,model_levels)  
        ENDDO  
      ENDDO  
      
! Initialise dust flux and friction velocity    
      dust_flux=0.0   
      u_s_t=0.0    
      u_s_t_dry=0.0      
      u_s_t_tile=0.0      
      u_s_t_dry_tile=0.0   
      dust_flux_tile=0.0
      
                  
! Find dust emission fluxes on land tiles
! DEPENDS ON: dust_srce
       CALL DUST_SRCE(                                                 &   
        land_pts,ntiles,sm_levels,tile_pts,tile_index,tile_frac,fland, &  
        pstar_land,tstar_tile,rhostar_land,soil_layer_moisture,        &  
        snow_tile,ustar_tile,mrel_land,clay_land,                      &  
        dust_flux_tile,u_s_t_tile,u_s_t_dry_tile                       &  
        )                                                               
                                                                        
! Calculate gridbox mean values of dust flux and ustar diagnostics      
      DO M = 1,ntiles                                                   
        DO N = 1,tile_pts(M)                                            
          L = tile_index(N,M)                                           
          J = (land_index(L)-1)/row_length + 1                          
          I = land_index(L) - (J-1)*row_length                          
          DO IDIV = 1,ndiv                                              
            dust_flux(I,J,IDIV) = dust_flux(I,J,IDIV) +                &  
              dust_flux_tile(L,M,IDIV)*tile_frac(L,M)                   
            u_s_t(I,J,IDIV) = u_s_t(I,J,IDIV) +                        &  
              u_s_t_tile(L,M,IDIV)*tile_frac(L,M)                       
            u_s_t_dry(I,J,IDIV) = u_s_t_dry(I,J,IDIV) +                &            
              u_s_t_dry_tile(L,M,IDIV)*tile_frac(L,M)                   
          ENDDO  !ndiv                                                  
        ENDDO !tile_pts                                                 
      ENDDO !ntiles                                                     
                                                                        
      DO IDIV = 1,NDIV                                                  
                                                                        
! Write the dust emission flux diagnostics                              
!        ITEM = 400+IDIV  ! dust emission flux                          
!        IF (sf(item,sect)) THEN 
! DEPENDS ON: copydiag
!          CALL COPYDIAG (stashwork(si(item,sect,IM_INDEX)),            
!     &      dust_flux(1:row_length,1:rows,IDIV),                       
!     &      row_length,rows,0,0,0,0, at_extremity,                     
!     &      atmos_im,sect,item,                                        
!     &      errorstatus,cmessage)                                      
!          IF (errorstatus.GT.0) THEN                                   
!            cmessage=": ERROR IN COPYDIAG (ITEM 401-6)"                
! DEPENDS ON: ereport
!            CALL EREPORT( 'UKCA_DUST_CTL', errorstatus, cmessage )     
!          ENDIF                                                        
!        ENDIF                                                          
                                                                        
! Do gravitational settling, dry deposition and add emissions           
        dustwork=dust_tracers(1:row_length,1:rows,:,IDIV)               

! DEPENDS ON: vgrav
        CALL VGRAV(                                                    &                          
          row_length,rows,lev1,drep(IDIV),                             &                          
          rhop,pstar,tstar,                                            &                          
          vstokes1,dustwork2,dustwork3                                 &                          
          )                                                             
                                                                        
        res_factor(:,:)=aresist(:,:)*(vstokes1(:,:)+                   &             
          1./(aresist(:,:)+rb_dust_ndivs(:,:,IDIV)+                    &             
          aresist(:,:)*rb_dust_ndivs(:,:,IDIV)*vstokes1(:,:)) )      
        
! DEPENDS ON: tr_mix
        CALL TR_MIX(                                                   &                                     
             halo_i, halo_j, row_length, rows, bl_levels               &           
            ,offx, offy, alpha_cd                                      &
            ,rhokh_mix(1,1,2), rho_aresist                             &
            ,dtrdz_charney_grid, r_rho_levels(:,:,1:bl_levels)         &
            ,r_theta_levels(:,:,0:bl_levels), timestep                 &
            ,tr_flux, dustwork( :,:,1:bl_levels)                       &
            ,dust_flux(:,:,IDIV), res_factor                           &
            ,drydep_str                                                &
            ,kent, we_lim, t_frac, zrzi                                &
            ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                &
            ,ml_depth ,zhsc, z_half                                    &
            ,errorstatus, .false.                                      &
             )                                                          
                                                                        
        IF (errorstatus /=0) THEN                                       
          cmessage='Error from TR_MIX'    
! DEPENDS ON: ereport
          CALL EREPORT( 'UKCA_DUST_CTL', errorstatus, CMESSAGE )        
        ENDIF                                                           
                                                                        
        item = 440 + IDIV      !dust dry dep flux, from lowest layer    
        IF ( SF(item,sect) ) THEN                                       
! Change sign of dry dep flux (otherwise negative)                      
          dustwork2=-drydep_str                                         
          write(6,*) 'UKCA dustcntl: ',item,' selected for o/p'         
          write(6,*) 'COPYDIAG not called!'                             
           
! DEPENDS ON: copydiag
!          CALL COPYDIAG(STASHWORK(SI(ITEM,SECT,IM_INDEX)),dustwork2    
!            ROW_LENGTH,ROWS,0,0,0,0,AT_EXTREMITY,                      
!            ATMOS_IM,sect,item,                                        
!            errorstatus,cmessage)                                      
!                                                                       
!          IF (errorstatus /=0) THEN  
! DEPENDS ON: ereport
!            CALL EREPORT( ROUTINENAME, ERRORSTATUS, CMESSAGE )         
!          ENDIF                                                        
                                                                        
        ENDIF !stashflag                                                
        
! DEPENDS ON: gravsett
        CALL GRAVSETT(                                                 &                                              
          row_length,rows,model_levels,dustwork,                       &            
          drep(IDIV),rhop,p_layer_centres,p_layer_boundaries,          &           
          t_modellevs,timestep,dustwork2,.false.                       &           
          )                                                             
                                                                        
        item = 450 + IDIV      !dust dry dep flux, from second layer    
        IF ( SF(item,sect) ) THEN                                       
          write(6,*) 'UKCA dustcntl: ',item,' selected for o/p'         
          write(6,*) 'COPYDIAG not called!'                             
           
! DEPENDS ON: copydiag
!          CALL COPYDIAG(STASHWORK(SI(item,sect,IM_INDEX)),dustwork2,   
!            row_length,rows,0,0,0,0,at_extremity,                      
!            atmos_im,sect,item,                                        
!            errorstatus,cmessage)                                      
!                                                                       
!          IF (errorstatus /=0) THEN                                    
!            cmessage='Error from COPYDIAG' 
! DEPENDS ON: ereport
!            CALL EREPORT( 'UKCA_DUST_CTL', ERRORSTATUS, CMESSAGE )     
!          ENDIF                                                        
                                                                        
        ENDIF  !Stash flag                                              
                                                                        
! write temporary array back into appropriate dust division             
          dust_tracers(1:row_length,1:rows,:,IDIV)=dustwork             
      ENDDO !ndiv                                                       
      
! DEPENDS ON: ukca_ls_washout
      CALL UKCA_LS_WASHOUT(                                            &                                        
        row_length, rows, model_levels, wet_levels,                    &      
        r_rho_levels, r_theta_levels, tr_levels, halo_i, halo_j,       &      
        rho_r2, timestep, q, qcl, qcf, ls_rain3d, ls_snow3d,           &      
        dust_tracers,                                                  &      
#include "argsts.h"                                              
        STASHwork                                                      &      
        )                                                               
                                                                        
      RETURN                                                                        
      END SUBROUTINE UKCA_DUST_CTL                                      
! ######################################################################
#endif
