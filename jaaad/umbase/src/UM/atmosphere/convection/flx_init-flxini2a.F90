#if defined(A05_4A) || defined(A05_5A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  SUBROUTINE FLX_INIT--------------------------------------------------
!
!  PURPOSE : CALCULATE INITIAL DOWNDRAUGHT MASSFLUX
!
!  SUITABLE FOR SINGLE COLUMN MODEL USE
!
!  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY SUMMER 1992
!
!  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!  VERSION NO. 4  DATED 5/2/92
!
!  SYSTEM TASK : P27
!
!  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
!
! ----------------------------------------------------------------------
!
! ARGUMENTS-------------------------------------------------------------
!
      SUBROUTINE FLX_INIT(NPNTS,KCT,ICCB,ICCT,FLX,FLX_DD_K,BDDI         &
                         ,FLX_STRT)

      Use cv_run_mod,                                                   &
          Only: dd_opt

      IMPLICIT NONE
!
!----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!----------------------------------------------------------------------
!
      INTEGER I                 ! LOOP COUNTER
!
      INTEGER NPNTS             ! IN NUMBER OF POINTS
!
      INTEGER KCT               ! IN CONVECTIVE CLOUD TOP
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!----------------------------------------------------------------------
!
      INTEGER ICCB(NPNTS)       ! IN CONVECTIVE CLOUD BASE
!
      INTEGER ICCT(NPNTS)       ! IN CONVECTIVE CLOUD TOP
!
      REAL FLX(NPNTS,KCT+1)     ! IN CONVECTIVE MASSFLUX (PA/S)
!
      LOGICAL BDDI(NPNTS)       ! IN MASK FOR THOSE POINTS WHERE
                                !    DOWNDRAUGHT MAY INITIATE
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE OUTPUT
!----------------------------------------------------------------------
!
      REAL FLX_DD_K(NPNTS)      ! OUT DOWNDRAUGHT MASSFLUX OF LAYER K
                                !     (PA/S)
!
      REAL FLX_STRT(NPNTS)      ! OUT UPDRAUGHT MASSFLUX AT LEVEL
                                !     DOWNDRAUGHT STARTS (PA/S)
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!----------------------------------------------------------------------
!
      INTEGER KDDREF            ! REFERENCE LEVEL FOR DOWNDRAUGHT
                                ! MASSFLUX
      REAL FLXSCALE             ! THE SCALING FACTOR FOR THE INITIAL
                                ! DOWNDRAUGHT MASSFLUX  
!
!----------------------------------------------------------------------
! CALCULATE DOWNDRAUGHT MASSFLUX BASED ON A REFERENCE LEVEL WHICH IS
! 3/4 CLOUD DEPTH
!----------------------------------------------------------------------
!
      IF (DD_OPT == 1) THEN
       FLXSCALE=0.1
      ELSE
       FLXSCALE=0.05
      END IF

      DO I=1,NPNTS
       IF (BDDI(I)) THEN
          KDDREF = ICCB(I) + 0.75*(ICCT(I) - ICCB(I))
          IF (KDDREF  >=  ICCT(I)-1) KDDREF=ICCT(I)-1
          FLX_STRT(I) = FLX(I,KDDREF)
          FLX_DD_K(I) = FLX_STRT(I) * FLXSCALE
       END IF
      END DO
!
      RETURN
      END SUBROUTINE FLX_INIT
!
#endif
