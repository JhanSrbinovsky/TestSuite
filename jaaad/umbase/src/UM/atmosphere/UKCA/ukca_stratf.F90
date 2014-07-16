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
! Purpose: Subroutine to overwrite values at top of model            
!          using interpolated 5-day fields from the 2-d model.          
!          Based on STRATF.F from Cambridge TOMCAT model and            
!          modified by Olaf Morgenstern to allow for flexible           
!          positioning of the NOy species.                              
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
! Code description:                                                     
!   Language: FORTRAN 90                                                
!   This code is written to UMDP3 v6 programming standards.             
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
      SUBROUTINE UKCA_STRATF(i_day_number, lon, lat, lev, p_field,     &
                             first_row, glon, glat, ntracer, sinlat,   & 
                             pl, um_ozone3d, p_above, tracer)         

      USE ASAD_MOD,          ONLY: sc, speci
      USE UKCA_TROPOPAUSE
      USE UKCA_D1_DEFS
      USE UKCA_CSPECIES
      IMPLICIT NONE  

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "param2d.h"
#include "c_v_m.h"
#include "parvars.h"                                            
#include "typsize.h"

      INTEGER, INTENT(IN) :: lon        ! No of longitudes            
      INTEGER, INTENT(IN) :: lat        ! No of latitudes             
      INTEGER, INTENT(IN) :: lev        ! No of levels                
      INTEGER, INTENT(IN) :: p_field    ! No of spatial points        
      INTEGER, INTENT(IN) :: first_row  ! First global latitude       
      INTEGER, INTENT(IN) :: glat       ! No of global latitudes      
      INTEGER, INTENT(IN) :: glon       ! No of global longitudes     
      INTEGER, INTENT(IN) :: ntracer    ! No of chemical tracers
                                                                              
      REAL, INTENT(IN) :: sinlat(p_field) ! Sine(latitude)            
      REAL, INTENT(IN) :: pl(lon,lat,lev) ! Model pressures                  
      REAL, INTENT(IN) :: um_ozone3d(lon,lat,lev) ! UM O3
      REAL, INTENT(IN) :: p_above  ! above which tbc are applied            

      REAL, INTENT(INOUT) :: tracer(lon,lat,lev,ntracer)
                                                                              
!     Local variables                                                       
                                                                               
      INTEGER, PARAMETER :: n_2d = 74  ! no of 2d fields in file 

      INTEGER :: idofy          ! Day number                          
      INTEGER :: ipos           ! Position in 2D file                        
      INTEGER :: i_day_number   ! Day number                                 
      INTEGER :: info           ! Tag used in communication                  
      INTEGER :: i,ij,j,k,l     ! Loop variables                                                 
      INTEGER :: kk             ! Pointer/index
      INTEGER :: m              ! Loop variable

      LOGICAL :: MASK(lon,lat,lev)  ! mask for overwriting stratosphere
      LOGICAL, SAVE :: firstcall = .true. 
      
      REAL, PARAMETER :: o3_hno3_ratio = 1.0/1000.0 ! kg[N]/kg[O3] from 
                                                    ! Murphy and Fahey 1994  
 
      REAL, SAVE :: ch4_2d(nolat,nolev, n_2d)! 2D data for ch4 
      REAL, SAVE :: o3_2d(nolat,nolev, n_2d) ! 2D data for o3 
      REAL, SAVE :: noy_2d(nolat,nolev, n_2d)! 2D data for noy 
                                                                                
      REAL :: fpos                        ! Decides which 2D field to 
      REAL :: delpos                                                         
      REAL :: fac                         ! Factor used to calculate zmean   
      REAL :: noy(nolat,nolev)            ! 2D field interpolated in time
      REAL :: noya(nolat,nolev)           ! 2D field stradding day no        
      REAL :: noyb(nolat,nolev)           ! 2D field stradding day no        
      REAL :: o3(nolat,nolev)             ! 2D field interpolated in time
      REAL :: o3a(nolat,nolev)            ! 2D field stradding day no        
      REAL :: o3b(nolat,nolev)            ! 2D field stradding day no        
      REAL :: ch4(nolat,nolev)            ! 2D field interpolated in time
      REAL :: ch4a(nolat,nolev)           ! 2D field stradding day no        
      REAL :: ch4b(nolat,nolev)           ! 2D field stradding day no        
      REAL :: o33d(lon,lat,lev)           ! 2D field interpolated onto 3D
      REAL :: noy3d(lon,lat,lev)          ! 2D field interpolated onto 3D    
      REAL :: ch43d(lon,lat,lev)          ! 2D field interpolated onto 3D
      REAL :: hno33d(lon,lat,lev)         ! 3D field from fixed o3:hno3 ratio      
      REAL :: noyzm(lat,lev)              ! NOy zonal mean            
      REAL :: hno3zm(lat,lev)             ! HNO3 zonal mean                  
      REAL :: hno4zm(lat,lev)             ! HNO4 zonal mean                  
      REAL :: nozm(lat,lev)               ! NOz zonal mean            
      REAL :: no2zm(lat,lev)              ! NO2 zonal mean                   
      REAL :: no3zm(lat,lev)              ! NO3 zonal mean                   
      REAL :: n2o5zm(lat,lev)             ! N2O5 zonal mean           
      REAL :: st(lon,lat,lev,jpspec)      ! Tracer MMRs               

!     Parameter to use UM ancil O3 instead of 2-D O3

      LOGICAL, PARAMETER :: L_use_UMO3=.true.
                                                
!     Parameter to use fixed o3 to hno3 ratio rather than 2d NOy

      LOGICAL, PARAMETER :: L_use_O3HNO3ratio =.true.

!     Parameter to overwrite stratosphere (fixed no of levels above tropopause)

      LOGICAL, PARAMETER :: L_all_strat   = .true.
      INTEGER, PARAMETER :: no_above_trop = 3                                
                                                                        
!     Check input unit numbers if change these                        
                                                                        
      IF (mype == 0 .AND. firstcall) THEN                                             

        OPEN(70,FILE='/data/cr/cce/hadcj/tropdata/ch4_topbound.dat')
        OPEN(71,FILE='/data/cr/cce/hadcj/tropdata/o3_topbound.dat')
        OPEN(72,FILE='/data/cr/cce/hadcj/tropdata/noy_topbound.dat')                                                    
                                                                        
!       Add up NOy - includes ClONO2 (~10-15%) in addition to         
!       NOy in 3-D. Also, NO3 is included implicitly in 2-D NOx + N2O5
                                                                        
        DO k = 1,n_2d                                                 
          DO j = 1,nolev                                              
            read(70,*) (ch4a (i,j),i=nolat,1,-1)                      
            read(71,*) (o3a  (i,j),i=nolat,1,-1)                      
            read(72,*) (noya (i,j),i=nolat,1,-1)                      
          ENDDO                                                       
          ch4_2d(:,:,k) = ch4a(:,:) 
          o3_2d (:,:,k) = o3a (:,:) 
          noy_2d(:,:,k) = noya(:,:) 
        ENDDO                                                         

        CLOSE(70)
        CLOSE(71)
        CLOSE(72)                                                               

      END IF         ! end if (mype.eq.0 .AND. firstcall)                             
         
      IF (firstcall) THEN 
        CALL GC_RBCAST(1,nolat*nolev*n_2d,0,nproc,info,ch4_2d)                  
        CALL GC_RBCAST(2,nolat*nolev*n_2d,0,nproc,info,o3_2d)                   
        CALL GC_RBCAST(3,nolat*nolev*n_2d,0,nproc,info,noy_2d)                  
        CALL GC_SSYNC (nproc,info)                                       
      ENDIF                                

!     Work out position for current timestep 
         
      idofy = i_day_number 
      fpos  = idofy/5.0 + 1.0 
      ipos  = nint(fpos) 
      IF ((fpos-ipos*1.0) < 0.0) ipos = ipos-1 
      delpos = fpos - ipos*1.0 
         
!     Reset IPOS if nearest day=365 or zero 
         
      IF (ipos == n_2d) ipos = 1 
         
      ch4a(:,:) = ch4_2d(:,:,ipos) 
      o3a (:,:) = o3_2d (:,:,ipos) 
      noya(:,:) = noy_2d(:,:,ipos) 
  
      ch4b(:,:) = ch4_2d(:,:,ipos+1) 
      o3b (:,:) = o3_2d (:,:,ipos+1) 
      noyb(:,:) = noy_2d(:,:,ipos+1)                                                                         

!     Interpolate in time                                             
                                                                        
      DO j = 1, nolev                                                 
        DO i = 1, nolat                                               
          ch4(i,j) = (ch4b(i,j)-ch4a(i,j))*delpos + ch4a(i,j)         
          o3 (i,j) = (o3b (i,j)-o3a (i,j))*delpos + o3a (i,j)         
          noy(i,j) = (noyb(i,j)-noya(i,j))*delpos + noya(i,j)         
        ENDDO                                                         
      ENDDO                                                           
                                                                        
!     Interpolate onto 3-D lats and levs                              
                                                                        
      IF (L_use_UMO3) THEN
        o33d(:,:,:)=um_ozone3d(:,:,:)/c_o3  ! convert to vmr
      ELSE
! DEPENDS ON: ukca_interp
        CALL UKCA_INTERP(o3, pl,lon,lat,lev,p_field,sinlat,o33d)
      ENDIF

      IF (L_use_O3HNO3ratio) THEN
        hno33d(:,:,:) = o33d(:,:,:)*o3_hno3_ratio*c_o3/c_n
      ELSE
! DEPENDS ON: ukca_interp
        CALL UKCA_INTERP(noy,pl,lon,lat,lev,p_field,sinlat,noy3d)       
      ENDIF  
! DEPENDS ON: ukca_interp
      CALL UKCA_INTERP(ch4,pl,lon,lat,lev,p_field,sinlat,ch43d)       
                                                                        
!     Set up 3D species array                                         
                                                                        
      st = 0.0                                                        
      DO j = 1,lat                                                    
        DO i = 1,lon                                                  
          k           = i+(j-1)*lon                                   
          st(i,j,:,:) = sc(k,:,:)                                     
        ENDDO                                                         
      ENDDO                                                           

      IF (.NOT. L_use_O3HNO3ratio) THEN  
                                                                            
!       Calculate zonal means for NOy species                            
! DEPENDS ON: ukca_calc_noy_zmeans 
         CALL UKCA_CALC_NOY_ZMEANS(lon, lat, lev,                   & 
     &                            glon, glat, first_row,            & 
     &                            st, noyzm, nozm, no2zm, no3zm,    & 
     &                            n2o5zm, hno3zm, hno4zm)                  

      ENDIF
                                                                              
!     Overwrite values in species array st
     
!     New N fields = zonal 3D field * (2d NOy/3d NOy) - at present onl
!     scale to total NOy - the partitioning between NOy species is    
!     not really correct because the species are treated differently  
!     in the two models. This can be improved. --- Comments from TOMCA
                                                                        
      IF (L_use_O3HNO3ratio .AND. (.NOT. L_all_strat)) THEN

!       Overwrite o3, hno3, and ch4 vmr at all gridboxes above p_above
!       O3   - either UM ancillary or Cambridge 2d
!       HNO3 - using fixed o3:hno3 ratio
!       CH4  - from Cambridge 2D model

        MASK(:,:,:) = .FALSE.
        DO l = 1, lev  
          DO k = 1, lat
            DO j = 1, lon                                                
             IF (pl(j,k,l) <= p_above) MASK(j,k,l) = .TRUE.
            ENDDO
          ENDDO
        ENDDO

        WHERE(MASK(:,:,:))                     
          st(:,:,:,nn_o3)    = o33d(:,:,:)    
          st(:,:,:,nn_o3s)   = o33d(:,:,:) 
          st(:,:,:,nn_hono2) = hno33d(:,:,:) 
          st(:,:,:,nn_ch4)   = ch43d(:,:,:) 
        ENDWHERE            
         
      ELSE IF (L_use_O3HNO3ratio .AND. L_all_strat) THEN 

!       Overwrite o3, hno3, and ch4 at all gridboxes a fixed
!       number of model levels above the tropopause
!       O3   - either UM ancillary or Cambridge 2d
!       HNO3 - using fixed o3:hno3 ratio
!       CH4  - from Cambridge 2D model

        MASK(:,:,:) = .FALSE.
        DO l = lev,1,-1                                                                                                                               
          MASK(:,:,l) = tropopause_level(:,:)+no_above_trop <= l
        ENDDO
                                  
        WHERE (MASK(:,:,:))                         
          st(:,:,:,nn_o3)    = o33d(:,:,:)   
          st(:,:,:,nn_o3s)   = o33d(:,:,:) 
          st(:,:,:,nn_hono2) = hno33d(:,:,:) 
          st(:,:,:,nn_ch4)   = ch43d(:,:,:) 
        ENDWHERE            

      ELSE

!       Overwrite o3, NOy and CH4 at all gridboxes above p_above
!       O3  - use either UM ancillary or Cambridge 2D
!       NOy - use Cambridge 2D
!       CH4 - use Cambridge 2D
        
        DO l = lev,1,-1                                                 
        IF (ANY(pl(:,:,l) <= p_above)) THEN                            
          DO j = 1,lat                                                
            DO i = 1,lon                                              
              IF (pl(i,j,l) <= p_above) THEN                           
                DO k = 1,jpspec                                       
                  SELECT CASE (speci(k))                              
                  CASE ('O3        ')                           ! O3                      
                    st(i,j,l,k)        = o33d(i,j,l)            ! vmr                       
                  CASE ('O3S       ')                           ! Strat O3        
                    st(i,j,l,k)        = o33d(i,j,l)                   
                  CASE ('N         ')                           ! N                          
                    st(i,j,l,k) = 0.0                                 
                  CASE ('NO        ')                           ! NO                         
                    st(i,j,l,k)=nozm(j,l)*noy3d(i,j,l)/noyzm(j,l)     
                  CASE ('NO3       ')                           ! NO3                        
                    st(i,j,l,k)=no3zm(j,l)*noy3d(i,j,l)/noyzm(j,l)    
                  CASE ('NO2       ')                           ! NO2                        
                    st(i,j,l,k)=no2zm(j,l)*noy3d(i,j,l)/noyzm(j,l)    
                  CASE ('N2O5      ')                           ! N2O5                       
                    st(i,j,l,k)=n2o5zm(j,l)*noy3d(i,j,l)/noyzm(j,l)   
                  CASE ('HO2NO2    ')                           ! HNO4                       
                    st(i,j,l,k)=hno4zm(j,l)*noy3d(i,j,l)/noyzm(j,l)  
                  CASE ('HONO2     ')                           ! HNO3                       
                    st(i,j,l,k)=hno3zm(j,l)*noy3d(i,j,l)/noyzm(j,l) 
                  CASE ('CH4       ')                           ! CH4                        
                    st(i,j,l,k)=ch43d(i,j,l)                         
                  END SELECT                                          
                ENDDO                                                 
              ENDIF                                                   
                                                                        
            ENDDO                                                     
          ENDDO                                                       
        ENDIF                                                         
        ENDDO                                                           
       
      ENDIF    ! End of IF(L_use_O3HNO3ratio and L_all_strat) statement                                                                              

      DO j = 1,lat                                                    
        DO i = 1,lon                                                  
          k = i+(j-1)*lon                                             
          sc(k,:,:) = st(i,j,:,:)                                     
        ENDDO                                                         
      ENDDO                                                           

!     Overwrite tracers using st array

      IF (L_ukca_family) THEN   ! family chemistry
        tracer(:,:,:,n_nox) = 0.0
        DO m=1,jpspec
          SELECT CASE (speci(m))
            CASE ('O3        ')
              tracer(:,:,:,n_ox)     = st(:,:,:,m)*c_o3        ! mmr
            CASE ('O3S       ')
              tracer(:,:,:,n_sx)     = st(:,:,:,m)*c_o3        ! mmr
            CASE ('N         ','NO        ',                           &
                  'NO2       ','NO3       ')
              tracer(:,:,:,n_nox)    = tracer(:,:,:,n_nox)             &
                                     + st(:,:,:,m)*c_species(m)
            CASE ('N2O5      ')
              tracer(:,:,:,n_n2o5)   = st(:,:,:,m)*c_n2o5
            CASE ('HO2NO2    ')
              tracer(:,:,:,n_ho2no2) = st(:,:,:,m)*c_ho2no2
            CASE ('HONO2     ')
              tracer(:,:,:,n_hono2)  = st(:,:,:,m)*c_hono2
            CASE ('CH4       ')
              tracer(:,:,:,n_ch4)    = st(:,:,:,m)*c_ch4
          END SELECT
        ENDDO

      ELSE                      ! non-family chemistry

        tracer(:,:,:,n_o3)     = st(:,:,:,nn_o3)*c_o3        ! mmr 
        tracer(:,:,:,n_o3s)    = st(:,:,:,nn_o3s)*c_o3       ! mmr 
        tracer(:,:,:,n_no)     = st(:,:,:,nn_no)*c_no 
        tracer(:,:,:,n_no2)    = st(:,:,:,nn_no2)*c_no2 
        tracer(:,:,:,n_no3)    = st(:,:,:,nn_no3)*c_no3 
        tracer(:,:,:,n_n2o5)   = st(:,:,:,nn_n2o5)*c_n2o5 
        tracer(:,:,:,n_ho2no2) = st(:,:,:,nn_ho2no2)*c_ho2no2 
        tracer(:,:,:,n_hono2)  = st(:,:,:,nn_hono2)*c_hono2 
        tracer(:,:,:,n_ch4)    = st(:,:,:,nn_ch4)*c_ch4 

      END IF           ! end of L_ukca_family statement                                                                        

      firstcall = .false. 
      
      RETURN                                                          
      END SUBROUTINE UKCA_STRATF                                      
#endif                                                                  
