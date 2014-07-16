#if defined(A04_3B) || defined(A04_3C) || defined(A04_3D)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!*LL  SUBROUTINE LSP_SCAV-----------------------------------------------
!LL
!LL  Purpose: Scavenge aerosol by large scale precipitation.
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.4  15/08/94  New routine. Pete Clark.
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3,
!LL                        Version 7, dated 11/3/93.
!LL
!LL  Logical component covered: Part of P26.
!LL
!LL  System task:
!LL
!LL  Documentation: Unified Model Documentation Paper No 26.
!*
!*L  Arguments:---------------------------------------------------------
      SUBROUTINE LSP_SCAV(                                              &
     & TIMESTEP,POINTS,RAIN,SNOW,AEROSOL                                &
     &)
      IMPLICIT NONE
      INTEGER                                                           &
                      ! Input integer scalar :-
     & POINTS         ! IN Number of points to be processed.
      REAL                                                              &
                      ! Input real scalar :-
     & TIMESTEP       ! IN Timestep (s).
      REAL                                                              &
                      ! Input real arrays :-
     & RAIN(POINTS)                                                     &
                      ! IN Rate of rainfall in this layer from
!                     !       above
!*                    !       (kg per sq m per s).
     &,SNOW(POINTS)   ! IN Rate of snowfall in this layer from
!                     !       above
!*                    !       (kg per sq m per s).
      REAL                                                              &
                      ! Updated real arrays :-
     & AEROSOL(POINTS) ! INOUT Aerosol mixing ratio
!*L   External subprogram called :-
!     EXTERNAL None
!-----------------------------------------------------------------------
!  Define local scalars.
!-----------------------------------------------------------------------
!  (a) Reals effectively expanded to workspace by the Cray (using
!      vector registers).
      REAL                                                              &
                      ! Real workspace.
     & KRAIN,KSNOW
      PARAMETER(KRAIN=1.0E-4,KSNOW=1.0E-4)
      REAL                                                              &
                      ! Real workspace.
     & RRAIN,RSNOW
!  (b) Others.
      INTEGER I       ! Loop counter (horizontal field index).
!
! Overall rate = KRAIN*(R) where R is in mm/hr=kg/m2/s*3600.0
      RRAIN=KRAIN*TIMESTEP*3600.0
      RSNOW=KSNOW*TIMESTEP*3600.0
      DO I=1,POINTS
        AEROSOL(I)=AEROSOL(I)/(1.0+RRAIN*RAIN(I)+RSNOW*SNOW(I))
      END DO
      RETURN
      END SUBROUTINE LSP_SCAV
#endif
