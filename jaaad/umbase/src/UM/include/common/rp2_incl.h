! ------------------------------------------------------------------
! Description
!   This comdeck includes the variables in STPH_RP that need to be
!   saved from timestep to timestep.
!
! ------------------------------------------------------------------
!
!  Variables read from already existing namelists
!
 REAL :: cape_mean      ! cape closure timescale mean value
 REAL :: RHCRITmean     ! critical relative humidity mean  
                         ! value - eq. rhcrit(3)
 REAL :: GWD_FRC_mean   ! gravity wave constant  mean value
 REAL :: KAY_GWAVE_mean ! critical froude number mean value
 REAL :: Charnock_mean  ! Charnock parameter mean value
 REAL :: M_CI_mean      ! multiplication parameter for CI and CIC
                        ! in 3C microphysics
 REAL :: G0_RP_mean     ! Stability function modification
                        ! parameter
 REAL :: par_mean       ! Neutral mixing length modification
                        ! parameter at t-1
                        
!
! The other variables are hardcoded in the UM code, and modified by
! STPH_RP where needed
!
 REAL    :: entcoef0    ! Entrainment rate coef. at t-1
 REAL    :: RHCRIT0     ! critical relative humidity at t-1
 REAL    :: par_mezcla0 ! Neutral mixing length modification
                        ! parameter at t-1
 REAL    :: G0_RP_0     ! Stability function modification
                        ! parameter at t-1
!              
! Only define if using 3B Microphysics scheme
!       
#if defined(A04_3B)
 REAL    :: CI_0             ! Ice fall speed variable at t-1
#endif
#if defined(A04_3C) || defined(A04_3D)
 REAL    :: M_CI_0           ! variable to modify ice fall speed at t-1
#endif
 REAL    :: cape_timescale_0 ! cape closure timescale at t-1
 REAL    :: GWD_FRC_0        ! gravity wave constant at t-1
 REAL    :: KAY_GWAVE_0      ! critical froude number at t-1
 REAL    :: Charnock_0      ! critical froude number at t-1
!
! Only include ice fall speed (M_CI_0) if using 3C/3D Microphysics scheme  
! Only include ice fall speed (CI_0) if using 3B Microphysics scheme. 
!   
#if defined(A04_3C)|| defined(A04_3D)
 COMMON/RP_INCL/cape_mean,RHCRITmean,Charnock_mean,M_CI_mean,     &
        G0_RP_mean, par_mean, GWD_FRC_mean,KAY_GWAVE_mean,        &
        entcoef0,RHCRIT0,G0_RP_0,M_CI_0,cape_timescale_0,         &
        GWD_FRC_0,KAY_GWAVE_0,par_mezcla0,                        &
        Charnock_0
#else
 COMMON/RP_INCL/cape_mean,RHCRITmean,Charnock_mean,G0_RP_mean,    &
        par_mean, GWD_FRC_mean,KAY_GWAVE_mean,entcoef0,           &
        RHCRIT0,G0_RP_0,cape_timescale_0,GWD_FRC_0,KAY_GWAVE_0,   &
        par_mezcla0,Charnock_0
#endif
#if defined(A04_3B)
 COMMON/RP_INCL/cape_mean,RHCRITmean,Charnock_mean,G0_RP_mean,    &
        par_mean,GWD_FRC_mean,KAY_GWAVE_mean,entcoef0,            &
        RHCRIT0,G0_RP_0,CI_0,cape_timescale_0,GWD_FRC_0,          &
        KAY_GWAVE_0,par_mezcla0,Charnock_0
#endif
