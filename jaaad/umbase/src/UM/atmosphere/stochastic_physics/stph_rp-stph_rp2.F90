#if defined(A35_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
       SUBROUTINE stph_rp2(max_model_levels,                            &
     &            RHCRIT,RHCRIT_max,RHCRIT_min,                         &
     &            CI_max,CI_min,                                        &
     &            GWD_FRC,GWD_FRC_max,GWD_FRC_min,                      &
     &            KAY_GWAVE,KAY_GWAVE_max,KAY_GWAVE_min,                &
     &            par_mezcla_max,par_mezcla_min,                        &
     &            G0_max,G0_min,G0_RP,par_mezcla,                       &    
     &            M_CI,M_CI_max,M_CI_min,                               &
     &            Charnock,Charnock_max,Charnock_min)
     
     Use cv_run_mod, Only:                                              &
           cape_timescale, cape_timescale_min, cape_timescale_max,      &
           entcoef_min, entcoef_max
           
       IMPLICIT NONE
!
! Description:  Version 2 of the random parameters scheme - based on stph_rp
!               add charnock parameter and code to modify ice fall speed in
!               3C microphysics
!
! Method:       The value of a given PARAMETER is chosen randomly
!               between a maximum and minimum values. PARAMETER's values
!               are temporally correlated (1st order markov process)
!
! Code author:  Sarah Beare
!               sarah.beare@metoffice.gov.uk
!
! Language:     FORTRAN 90 + common extensions
!
!-------------------------------------------------------------
! Variable definition
!-------------------------------------------------------------
!
! IN variables
!
 INTEGER, INTENT(in) :: max_model_levels
!
! Maximum and minimum values for the STPH_RP scheme
! Large Scale Precipitation
!
 REAL, INTENT(in) :: RHCRIT_max  ! Max value critical relative humidity
 REAL, INTENT(in) :: RHCRIT_min  ! Min value critical relative humidity
 REAL, INTENT(in) :: CI_max      ! Max value ice fall speed parameter
 REAL, INTENT(in) :: CI_min      ! Min value ice fall speed parameter
!
! Maximum and minimum values for the STPH_RP scheme
! Gravity Wave drag
!
 REAL, INTENT(in) :: GWD_FRC_max   ! Max value critical froude number
 REAL, INTENT(in) :: GWD_FRC_min   ! Min value critical froude number                 
 REAL, INTENT(in) :: KAY_GWAVE_max ! Max value gravity wave parameter
 REAL, INTENT(in) :: KAY_GWAVE_min ! Min value gravity wave parameter
!
! Maximum and minimum values for the STPH_RP scheme
! Boundary Layer
!
 REAL, INTENT(in) :: par_mezcla_max ! Max value neut. mix length param.
 REAL, INTENT(in) :: par_mezcla_min ! Min value neut. mix length param.
 REAL, INTENT(in) :: G0_max         ! Max value of parameter to modify 
                                    ! stability functions
 REAL, INTENT(in) :: G0_min         ! Min value of parameter to modify 
                                    ! stability functions
 REAL, INTENT(in) :: Charnock_max   ! Max value Charnock param.
 REAL, INTENT(in) :: Charnock_min   ! Min value Charnock param.                             
 REAL, INTENT(in) :: M_CI_max       ! Max value of multiplication factor for CI  
 REAL, INTENT(in) :: M_CI_min       ! Min value of multiplication factor for CI                             
!                                  
!                                  
! IN/OUT variables
!
 REAL,INTENT(InOut) :: KAY_GWAVE      ! Surface stress constant GWD
 REAL,INTENT(InOut) :: RHCRIT(max_model_levels)! Crit RH for layer 
                                               ! cloud formation
 REAL,INTENT(InOut) :: GWD_FRC        ! Critical Froude Number 
                                      ! for 4A Scheme
 REAL,INTENT(InOut) :: par_mezcla     ! Used to vary LAMBDAH,LAMBDAM 
                                      ! (neutral mixing length)
 REAL,INTENT(InOut) :: G0_RP          ! Used to modify stability 
                                      ! functions in EXCOEF
 REAL,INTENT(InOut) :: M_CI           ! used to vary CI in 3c microphysics
 REAL,INTENT(InOut) :: Charnock       ! Charnock Parameter
!
! GLOBAL VARIABLES
!
#include "entcnst.h"
#include "c_lspdrp.h"
#include "cprintst.h"
#include "parvars.h"
#include "cntlatm.h"
!
! LOCAL VARIABLES
#include "rp2_incl.h"
!
! Local variables to assign values to RHCRIT, which is level
! dependant

 INTEGER,ALLOCATABLE,SAVE :: interval(:)
 REAL,   ALLOCATABLE,SAVE :: drhcrit(:)     ! used to store differences
                                            ! in rhcrit between model levels
 INTEGER,SAVE             :: max_interval
!
 LOGICAL,SAVE   :: primero=.true. ! Logical to calculate the random
                                  ! seed only the first time
 LOGICAL        :: local_switch   ! Logical to calculate the first new 
                                  ! rhcrit value.
 INTEGER        :: i,j            ! loop variables
!
! Variables associated with random number generator
!
 REAL    :: k                     ! random number
 REAL    :: shock                 ! random value used in first-order
                                  ! autoregression the random seed is 
                                  ! based on actual time (milisecond)
 INTEGER :: DT(8)                 ! DT values to keep the time info for
                                  ! the random seed
 INTEGER :: istat                 ! To record call in the synch.
 INTEGER :: tam                   ! size of the random seed
 INTEGER,ALLOCATABLE :: kk(:)     ! New seed for random num generator
 INTEGER,ALLOCATABLE :: prev(:)   ! Old seed for random num generator
 INTEGER             :: icode         ! Error return strings
!
!Variables associated with autoregression model
!
 REAL,PARAMETER :: corr=0.95      ! Correlation value for the
                                  ! first-order autoregression
 REAL           :: maxran, minran ! max and min values for the random 
                                  ! value in first-order autoregression
                                  ! Same variable (different values)
                                  ! for all params
!
! Define mean values
!
#if defined(A05_4A)
 REAL,PARAMETER :: entcoefmean=3.0 ! entrainment rate coefficient mean
#endif
! REAL,PARAMETER :: par_mean=0.15   ! neutral mixing length param mean
! REAL,PARAMETER :: G0_RP_mean=10.0    ! stability function parameter mean
#if defined(A04_3B)
 REAL,PARAMETER :: CI_mean=25.2    ! Ice fall speed mean value
#endif
!-------------------------------------------------------------
! Print old parameter values if requested
!-------------------------------------------------------------
IF (PrintStatus  >  PrStatus_Normal) THEN
    WRITE(6,*) ' STPH_RP2 '
    WRITE(6,*) ' STPH_RP2 OLD PARAM VALUES'
#if defined(A05_4A)
    WRITE(6,*) 'entcoef .........', entcoef
#endif
    WRITE(6,*) 'rhcrit(4,7,10) ..', rhcrit(4),rhcrit(7),rhcrit(10)
#if defined(A04_3B)
    WRITE(6,*) 'CI (ice speed)...', CI
#endif
#if defined(A04_3C)|| defined(A04_3D)
    WRITE(6,*) 'M_CI (ice speed)...', M_CI
#endif
    WRITE(6,*) 'coef for LAMBDA .', par_mezcla
    WRITE(6,*) 'G0_RP ..............', G0_RP
    WRITE(6,*) 'CAPE ............', cape_timescale
    WRITE(6,*) 'GWD_FRC..........', gwd_frc
    WRITE(6,*) 'KAY_GWAVE........', kay_gwave
    WRITE(6,*) 'Charnock..', Charnock
ENDIF
!-------------------------------------------------------------
!                  STOCH PHYS
! If this is a new simulation a random seed is generated the first 
! time random parameters is called.  Original variable values are 
! passed as well .
!-----------------------------------------------------------------
 IF (primero) THEN
    IF (PrintStatus  >=  PrStatus_Normal) THEN
        WRITE(6,*) 'L_RPSEED_READ', L_RPSEED_READ
        WRITE(6,*) 'L_RPSEED_WRITE', L_RPSEED_WRITE
    ENDIF
    IF (L_RPSEED_READ) THEN
       ! ---------------------------------------
       ! read in random seed from file.
       ! ---------------------------------------
        CALL date_and_time(VALUES=DT)
        CALL random_seed(SIZE=tam)
        IF (.not.allocated(prev)) ALLOCATE(prev(tam))
        IF (.not.allocated(kk)) ALLOCATE(kk(tam))
        
        CALL random_seed(GET=prev(1:tam))
        
        IF (PrintStatus  >=  PrStatus_Normal) THEN
           WRITE(6,*)  'OPENING RANDOM SEED FILE'
        ENDIF
        ! DEPENDS ON: RPOpeninput
        CALL RPOpenInput()
     
        IF (PrintStatus  >=  PrStatus_Normal) THEN
            WRITE(6,*)  'READING RANDOM SEED FILE'
        ENDIF
        ! DEPENDS ON: RPReadEntry
        CALL RPReadEntry(kk,tam)

        IF (PrintStatus  >=  PrStatus_Normal) THEN
            WRITE(6,*)  'Closing RANDOM SEED FILE'
        ENDIF
        ! DEPENDS ON: RPCloseinput
        CALL RPCloseInput()

        CALL random_seed(PUT=kk(1:tam))
        IF (PrintStatus  >=  PrStatus_Normal) THEN
          WRITE(6,*)  'READING RANDOM SEED FROM FILE'
          WRITE(6,*) ' Size .......', tam
          WRITE(6,*) ' Prev. seed .', prev
          WRITE(6,*) ' seed ...', kk
       ENDIF
       
       DEALLOCATE(prev)
       DEALLOCATE(kk)
    ELSE
       ! ---------------------------------------
       ! generate new random seed and write to file .
       ! ---------------------------------------
       CALL date_and_time(VALUES=DT)
       CALL random_seed(SIZE=tam)
       
       IF (.not.allocated(prev)) ALLOCATE(prev(tam))
       IF (.not.allocated(kk)) ALLOCATE(kk(tam))
       
       CALL random_seed(GET=prev(1:tam))
       IF (DT(8) <  100) THEN
           DT(8)=DT(8)+100
       ENDIF
       
       kk(:)=abs(prev(:)-(DT(8)*DT(8)*100000000000))
       CALL random_seed(PUT=kk(1:tam))
       
       ! DEPENDS ON: RPOpenOutput
       CALL RPOpenoutput()
       
        ! DEPENDS ON: RPwriteEntry
       CALL RPwriteEntry(kk,tam)
       
        ! DEPENDS ON: RPCloseoutput
       CALL RPCloseoutput()
       
       IF (PrintStatus  >=  PrStatus_Normal) THEN
          WRITE(6,*) ' GENERATION OF SEED '
          WRITE(6,*) ' Size .......', tam
          WRITE(6,*) ' Prev. seed .', prev
          WRITE(6,*) ' New seed ...', kk
       ENDIF
       
       DEALLOCATE(prev)
       DEALLOCATE(kk)
    ENDIF
    ! ---------------------------------------
    ! Allocate interval the first time STPH_RP is called only
    ! This will be used when perturbing RHCRIT.
    ! ---------------------------------------
    IF (.not.allocated(interval)) THEN
      ALLOCATE (interval(max_model_levels))
      interval(:)=0
    ENDIF
    IF (.not.allocated(drhcrit)) THEN
      ALLOCATE (drhcrit(max_model_levels))
      drhcrit(:)=0
    ENDIF
    j=1
    DO i=1,max_model_levels
     IF(i >  2) THEN
      IF(rhcrit(i) <  rhcrit(i-1)) THEN          
         drhcrit(i)=rhcrit(i-1)-rhcrit(i)  
         interval(i)=j
         j=j+1
      ENDIF
     ENDIF
    ENDDO
    max_interval=j-1
!
    DO i=1,max_model_levels
     IF(i >  3) THEN
      IF(interval(i) <  interval(i-1)) THEN
            interval(i)=interval(i-1)
      ENDIF
      IF(drhcrit(i) <  drhcrit(i-1)) THEN
            drhcrit(i)=drhcrit(i-1)
      ENDIF
     ENDIF
    ENDDO
 primero=.false.
!
! Pass original values (as read from UMUI) of some variables
!
    cape_mean  = cape_timescale ! cape closure timescale mean value
    RHCRITmean = rhcrit(3)      ! critical RH mean value 
                                 ! - eq. rhcrit(3)
    GWD_FRC_mean   = GWD_FRC    ! Gravity wave stress mean value
    KAY_GWAVE_mean = KAY_GWAVE  ! Critical froude number mean value
    Charnock_mean  = Charnock   ! Charnock parameter mean value
    G0_RP_mean     = G0_RP      ! parameter to modify 
                                ! stability functions mean value
    par_mean       = par_mezcla ! neut. mix length param mean value
#if defined(A04_3C)|| defined(A04_3D)
    M_CI_mean      = M_CI       ! multiplication factor for CI and CIC
                                ! in 3C microphysics
#endif
!                                         
! Timestep 1: set values at t-1 (rp_incl.h) to be equal to mean value
! 
    cape_timescale_0 = cape_mean
    RHCRIT0          = RHCRITmean
    GWD_FRC_0        = GWD_FRC_mean
    KAY_GWAVE_0      = KAY_GWAVE_mean
    Charnock_0       = Charnock_mean
!
#if defined(A05_4A)
    entcoef0         = entcoefmean
#endif
    par_mezcla0      = par_mean
    G0_RP_0          = G0_RP_mean
#if defined(A04_3B)
    CI_0             = CI_mean
#endif
#if defined(A04_3C)|| defined(A04_3D)
    M_CI_0             = M_CI_mean
#endif
 ENDIF
!----------------------------------------------------------------------
! Each time random parameters is called broadcast the random number
! to all processors
!----------------------------------------------------------------------
 IF(mype == 0) then
    CALL random_number(k)
 END IF
 CALL gc_rbcast(3145,1,0,nproc,istat,k)
!
#if defined(A05_4A)
!----------------------------------------------------------------------
!Bit dedicated to the entrainment rate coefficient
!(Convection)
!----------------------------------------------------------------------
 maxran=(entcoef_max-entcoef_min)/3
 minran=maxran*(-1)
 shock=k*(maxran-minran)+minran
!
 entcoef=entcoefmean+((entcoef0-entcoefmean)*corr)+shock
 IF (entcoef <  entcoef_min) entcoef=entcoef_min
 IF (entcoef >  entcoef_max) entcoef=entcoef_max
 entcoef0=entcoef
#endif
!----------------------------------------------------------------------
!Bit dedicated to CAPE_TIMESCALE (CONVECTION)
!----------------------------------------------------------------------
 maxran=(cape_timescale_max-cape_timescale_min)/3
 minran=maxran*(-1)
 shock=k*(maxran-minran)+minran
!
 cape_timescale=cape_mean+((cape_timescale_0-cape_mean)*corr)+shock
! 
 IF (cape_timescale < cape_timescale_min)                               &
     cape_timescale=cape_timescale_min
 IF (cape_timescale >  cape_timescale_max)                              &
     cape_timescale=cape_timescale_max
 cape_timescale_0=cape_timescale
!
!----------------------------------------------------------------------
!Bit dedicated to the RHcrit (LSP)
!----------------------------------------------------------------------
 maxran=(RHCRIT_max-RHCRIT_min)/3
 minran=maxran*(-1)
 shock=k*(maxran-minran)+minran
!
 local_switch = .TRUE. 
 DO j=1,max_interval
  DO i=1,max_model_levels
   IF (j == 1) THEN
    IF (interval(i) == j) THEN
     IF (local_switch) THEN 
       RHCRIT(i)=RHCRITmean+((RHCRIT0-RHCRITmean)*corr)+shock
       IF (RHCRIT(i) <  RHCRIT_min) RHCRIT(i)=RHCRIT_min
       IF (RHCRIT(i) >  RHCRIT_max) RHCRIT(i)=RHCRIT_max
       RHCRIT0=RHCRIT(i)
       local_switch=.FALSE.
     ELSE
       RHCRIT(i)=RHCRIT(i-1)
     ENDIF  
    ENDIF
   ELSE
    IF (interval(i) == j) THEN
      RHCRIT(i)=RHCRIT0-(drhcrit(i)*(j-1))
    ENDIF
   ENDIF
  ENDDO
 ENDDO
! RHCRIT0=RHCRIT(3)
!
#if defined(A04_3B)
!----------------------------------------------------------------------
!Bit dedicated to CI  (Ice fall speed) - 3B Microphysics 
!----------------------------------------------------------------------
 maxran=(CI_max-CI_min)/3
 minran=maxran*(-1)
 shock=k*(maxran-minran)+minran
!
 CI=CI_mean+((CI_0-CI_mean)*corr)+shock
 IF (CI <  CI_min) CI=CI_min
 IF (CI >  CI_max) CI=CI_max
 CI_0=CI
!
#endif
#if defined(A04_3C)|| defined(A04_3D)
!----------------------------------------------------------------------
!Bit dedicated to M_CI  (Ice fall speed) - 3C Microphysics LSPCON3C
!----------------------------------------------------------------------
       maxran=(M_CI_max-M_CI_min)/3
       minran=maxran*(-1)
       shock=k*(maxran-minran)+minran
!
       M_CI=M_CI_mean+((M_CI_0-M_CI_mean)*corr)+shock
       IF (M_CI.lt.M_CI_min) M_CI=M_CI_min
       IF (M_CI.gt.M_CI_max) M_CI=M_CI_max
       M_CI_0=M_CI
#endif
!----------------------------------------------------------------------
! Bit dedicated to neutral mixing length - EXCOEF (PBL)
!----------------------------------------------------------------------
 maxran=(par_mezcla_max-par_mezcla_min)/3
 minran=maxran*(-1)
 shock=k*(maxran-minran)+minran
!
 par_mezcla=par_mean+((par_mezcla0-par_mean)*corr)+shock
 IF (par_mezcla <  par_mezcla_min) par_mezcla=par_mezcla_min
 IF (par_mezcla >  par_mezcla_max) par_mezcla=par_mezcla_max
 par_mezcla0=par_mezcla
!----------------------------------------------------------------------
!Bit dedicated to Stability function parameter in EXCOEF (PBL)
!----------------------------------------------------------------------
 maxran=(G0_max-G0_min)/3
 minran=maxran*(-1)
 shock=k*(maxran-minran)+minran
!
 G0_RP=G0_RP_mean+((G0_RP_0-G0_RP_mean)*corr)+shock
 IF (G0_RP <  G0_min) G0_RP=G0_min
 IF (G0_RP >  G0_max) G0_RP=G0_max
 G0_RP_0=G0_RP
!----------------------------------------------------------------------
!Bit dedicated to Charnock parameter (PBL)
!----------------------------------------------------------------------
 maxran=(Charnock_max-Charnock_min)/3
 minran=maxran*(-1)
 shock=k*(maxran-minran)+minran
!
 Charnock=Charnock_mean+((Charnock_0-Charnock_mean)*corr)+shock
 IF (Charnock <  Charnock_min) Charnock=Charnock_min
 IF (Charnock >  Charnock_max) Charnock=Charnock_max
 Charnock_0=Charnock
!----------------------------------------------------------------------
! Bit dedicated to GWD_FRC - surface stress constant (GWD)
!----------------------------------------------------------------------
 maxran=(GWD_FRC_max-GWD_FRC_min)/3
 minran=maxran*(-1)
 shock=k*(maxran-minran)+minran
!
 GWD_FRC=GWD_FRC_mean+((GWD_FRC_0-GWD_FRC_mean)*corr)+shock
 IF (GWD_FRC <  GWD_FRC_min) GWD_FRC=GWD_FRC_min
 IF (GWD_FRC >  GWD_FRC_max) GWD_FRC=GWD_FRC_max
 GWD_FRC_0=GWD_FRC
!----------------------------------------------------------------------
! Bit dedicated to KAY_GWAVE - critical froude number (GWD)
!----------------------------------------------------------------------
 maxran=(KAY_GWAVE_max-KAY_GWAVE_min)/3
 minran=maxran*(-1)
 shock=k*(maxran-minran)+minran
!
 KAY_GWAVE=KAY_GWAVE_mean+((KAY_GWAVE_0-KAY_GWAVE_mean)*corr)+shock
 IF (KAY_GWAVE <  KAY_GWAVE_min) KAY_GWAVE=KAY_GWAVE_min
 IF (KAY_GWAVE >  KAY_GWAVE_max) KAY_GWAVE=KAY_GWAVE_max
 KAY_GWAVE_0=KAY_GWAVE
!----------------------------------------------------------------------
! Print new values of stochastic parameters if requested.
!----------------------------------------------------------------------
 IF (PrintStatus  >=   PrStatus_Normal) THEN
    WRITE(6,*) ' STPH_RP2 NEW PARAM VALUES'
#if defined(A05_4A)
    WRITE(6,*) 'entcoef .........', entcoef
#endif
    WRITE(6,*) 'rhcrit(4,7,10) ..', rhcrit(4),rhcrit(7),rhcrit(10)
#if defined(A04_3B)
    WRITE(6,*) 'CI (ice speed)...', CI
#endif
#if defined(A04_3C)|| defined(A04_3D)
       WRITE(6,*) 'M_CI........', M_CI
#endif
    WRITE(6,*) 'coef for LAMBDA .', par_mezcla
    WRITE(6,*) 'G0_RP ..............', G0_RP
    WRITE(6,*) 'CAPE ............', cape_timescale
    WRITE(6,*) 'GWD_FRC..........', gwd_frc
    WRITE(6,*) 'KAY_GWAVE........', kay_gwave
    WRITE(6,*) 'RHCRIT all..', rhcrit
    WRITE(6,*) 'Charnock..', Charnock
 ENDIF
!
RETURN
END SUBROUTINE stph_rp2
#endif
