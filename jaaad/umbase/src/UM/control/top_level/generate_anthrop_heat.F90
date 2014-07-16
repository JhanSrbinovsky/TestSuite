! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the anthropogenic contribution to surface heat fluxes for 
! urban tiles by linear interpolation of monthly values. The value is 
! then passed in anthrop_heat(n), that has a value of 0.0 except when
! n=6 (urban) and l_anthrop_heat=.true., and added to the surface 
! heat fluxes in sf_expl, sf_impl and sf_impl_2.
!
! Original code from Martin Best and Peter Clark (December 2005).
! Updated for UM7.1 by Jorge Bornemann (May 2008)
! 

      Subroutine Generate_Anthropogenic_Heat(                           &
     & val_year, val_day_number, val_hour, val_minute, val_second       &
     &, ntiles, anthrop_heat, l_anthrop_heat_src                        &
     & ) 

      IMPLICIT NONE

! PARAMETER definition of urban tile
      INTEGER , PARAMETER :: urban_tile=6
          
! IN time information for current timestep  
      INTEGER, INTENT(IN) ::                                            &
     &  val_year                                                        & 
     &, val_day_number                                                  &
     &, val_hour                                                        &
     &, val_minute                                                      &
     &, val_second
! IN Number of tiles   
      INTEGER, INTENT(IN) :: ntiles
! IN Switch for  anthropogenic heat source 
      LOGICAL, INTENT(IN) :: l_anthrop_heat_src 
! OUT   
      REAL, INTENT(OUT) :: anthrop_heat(ntiles) 
!                          ! OUT Additional heat source on tiles 
!                          ! For anthropogenic urban heat source  W/m2
! Local Variables
      REAL                                                              &
     & urban_month(12)                                                  &
                           ! Monthly values of anthropogenic
                           ! contribution to surface heat fluxes W/m2
     & ,mm, dpm
      INTEGER i,im,im1

!
! urban anthropogenic heat source is taken from the digest of energy 
! statistics (1995 - 2003) monthly averaged to 9 years, converted
! to w/m2 and adjusted to fraction dissipated in urban areas
!

      DATA urban_month /                                                & 
     &  25.1447                                                         &
     &, 23.4529                                                         &
     &, 24.6175                                                         &
     &, 20.7802                                                         &
     &, 18.8677                                                         &
     &, 18.1567                                                         &
     &, 17.1881                                                         &
     &, 16.6732                                                         &
     &, 18.5490                                                         &
     &, 20.7435                                                         &
     &, 22.9870                                                         &
     &, 26.1500 / 
                  
      dpm=365.0/12.0    
      DO i=1,ntiles
        anthrop_heat(i)=0.0
      END DO  
      
      mm=val_day_number/dpm-0.5
      
      IF (mm .lt. 0.0) THEN
      
        im=12
        im1=1
        mm=mm+1.0
        
      ELSE
      
        im=int(mm)
        mm=mm-im
        im=im+1
        im1=im+1
        
      END IF
      
      IF (l_anthrop_heat_src) THEN
        anthrop_heat(urban_tile)= mm * urban_month(im1)                  &
     &                           + (1.0-mm) * urban_month(im)
      END IF
      
      RETURN
      END
