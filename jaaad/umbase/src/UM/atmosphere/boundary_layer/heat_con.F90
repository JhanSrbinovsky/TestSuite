#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)    \
 || defined(A08_7A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE HEAT_CON----------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE HEAT_CON(NPNTS,HCON,STHU,STHF,                         &
     &                    V_SAT,HCONS                                   &
! LOGICAL LTIMER
     &,LTIMER                                                           &
     &)

! module for land-surface namelist
      USE LAND_SURF_MOD, ONLY :                                         &
     & SOILHC_METHOD
!        SOILHC_METHOD=1: Method of Cox et al (1999).
!        SOILHC_METHOD=2: Simplified Johansen (1975). 

      IMPLICIT NONE
!
! Description:
!    Calculates the soil thermal conductivity including the
!    effects of water and ice.
!
! Method == 1
!      Described in Cox et al (1999), Appendix B.
!      http://www.springerlink.com/content/9b459pyfhyjwk1ln/
!
! Method == 2
!      Simplified Johansen (1975). 
!      See http://www-nwp/~frid/thermal_conductivity.pdf
!                                             (Dharssi, Jan 2008)  
!
! Documentation : UM Documentation Paper 25
!
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.1               New deck.    Peter Cox
!LL   4.5   18/06/98  Changed Timer calls to indicate non-barrier
!LL                                                   P.Burton
!  5.5   17/04/03 Remove references to obsolete sections
!                 A03_5A,7A. T.White
!
! Code Description:
! Fortran 90 
!
! System component covered: P25
! System Task: P25
!


! Subroutine arguments:
!   Scalar arguments with intent(IN) :
      INTEGER                                                           &
     & NPNTS              ! IN Number of gridpoints

      REAL                                                              &
     & HCON(NPNTS)                                                      &
                          ! IN Dry soil thermal conductivity (W/m/K).
     &,STHU(NPNTS)                                                      &
                          ! IN Fractional saturation of unfrozen water
!                         !    at layer boundaries.
     &,STHF(NPNTS)                                                      &
                          ! IN Fractional saturation of frozen water
!                         !    at layer boundaries.
     &,V_SAT(NPNTS)       ! IN Volumetric soil moisture concentration
!                         !    at saturation (m3/m3 soil).
!
      LOGICAL LTIMER      ! Logical switch for TIMER diags

!   Array arguments with intent(OUT) :
      REAL                                                              &
     & HCONS(NPNTS)       ! OUT The thermal conductivity between adjacen
!                         !     layers including effects of water and ic
!                         !     (W/m/K).
!
#include "c_soilh.h"
! Local scalars:
      INTEGER                                                           &
     & I,J                ! WORK Loop counter.

! Local arrays:
      REAL                                                              &
     & HCSAT(NPNTS)                                                     &
                          ! WORK The thermal conductivity of the
!                         !  saturated  soil at current ratio of ice to
!                         !      liquid water (W/m/K).
     &,STH(NPNTS)                                                       &
                          ! WORK Fractional saturation of water
!                         !     (liquid+ice) at layer boundaries.
     &,THICE(NPNTS)                                                     &
                          ! WORK The concentration of ice at saturation
!                         !      for the current mass fraction of liquid
!                         !      water (m3 H2O/m3 soil).
     &,THWAT(NPNTS)                                                     &
                          ! WORK The concentration of liquid water at
!                         !      saturation for the current mass
!                         !  fraction of liquid water  (m3 H2O/m3 soil).
     &,KE(NPNTS)          ! WORK Kersten number
!-----------------------------------------------------------------------
! Local Parameters (Source: "The Frozen Earth" p.90)
!-----------------------------------------------------------------------
      REAL                                                              &
     & HCAIR                                                            &
                          ! Thermal conductivity of air (W/m/K).
     &,HCICE                                                            &
                          ! Thermal conductivity of ice (W/m/K).
     &,HCWAT              ! Thermal conductivity of liquid water (W/m/K)
      PARAMETER (HCAIR=0.025,HCICE=2.24,HCWAT=0.56)

!-----------------------------------------------------------------------
! Following local parameter values determined by linear regression,
! see Dharssi(2008). 
!-----------------------------------------------------------------------
      REAL                                                              &
     & HCSAT_MAX_ALLOWED                                                &
     &,HCSAT_MIN_ALLOWED                                                &
     &,HCSAT_GRADIENT                                                   &
     &,HCSAT_INTERCEPT
      PARAMETER ( HCSAT_MAX_ALLOWED = 2.20, HCSAT_MIN_ALLOWED = 1.58    &
     &           ,HCSAT_GRADIENT    = 12.4, HCSAT_INTERCEPT   = 0.25 )

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('HEATCON ',103)
      ENDIF

!----------------------------------------------------------------------
! Initialise all points
!----------------------------------------------------------------------
      DO I=1,NPNTS
        IF (V_SAT(I) >  0.0) THEN ! Soil points
          HCONS(I)=HCON(I)
        ELSE ! Ice points
          HCONS(I)=SNOW_HCON
        ENDIF
      ENDDO

      IF (SOILHC_METHOD==1) THEN
      DO I=1,NPNTS
!---------------------------------------------------------------
! Only do calculation for non land-ice pts
! V_SAT is set to zero for land-ice points
!---------------------------------------------------------------
        IF (V_SAT(I) >  0.0) THEN

          IF (STHU(I) >  0.0) THEN
            THWAT(I)=V_SAT(I)*STHU(I)/(STHU(I)+STHF(I))
          ELSE
            THWAT(I)=0.0
          ENDIF

          IF (STHF(I) >  0.0) THEN
            THICE(I)=V_SAT(I)*STHF(I)/(STHU(I)+STHF(I))
          ELSE
            THICE(I)=0.0
          ENDIF

          STH(I)=STHU(I)+STHF(I)
          HCSAT(I)=HCON(I)*(HCWAT**THWAT(I))*(HCICE**THICE(I))          &
     &                   /(HCAIR**V_SAT(I))
          HCONS(I)=(HCSAT(I)-HCON(I))*STH(I)+HCON(I)
        ENDIF

      ENDDO
      ENDIF ! (SOILHC_METHOD==1)
      
      IF (SOILHC_METHOD==2) THEN
      DO I=1,NPNTS
        IF (V_SAT(I) >  0.0) THEN

          IF (STHF(I) >  0.0) THEN
            THICE(I)=V_SAT(I)*STHF(I)/(STHU(I)+STHF(I))
          ELSE
            THICE(I)=0.0
          ENDIF
  
          THWAT(I)=V_SAT(I)-THICE(I)

          STH(I)=STHU(I)+STHF(I)
  
          HCSAT(I)=HCSAT_MIN_ALLOWED                                    &
     &            +HCSAT_GRADIENT*(HCON(I)-HCSAT_INTERCEPT)

          IF(HCSAT(I) > HCSAT_MAX_ALLOWED) HCSAT(I)=HCSAT_MAX_ALLOWED
          IF(HCSAT(I) < HCSAT_MIN_ALLOWED) HCSAT(I)=HCSAT_MIN_ALLOWED
  
          ! Adjust HCSAT for frozen soil water  
          HCSAT(I)=HCSAT(I)*(HCWAT**THWAT(I))*(HCICE**THICE(I))         &
     &            /(HCWAT**V_SAT(I))

          IF(STH(I) <= 0.1) THEN
            KE(I)=0.0
          ELSE
            KE(I)=1.0 + LOG10(STH(I))
          ENDIF

          HCONS(I)=(HCSAT(I)-HCON(I))*KE(I)+HCON(I)
        ENDIF

      ENDDO      
      ENDIF ! (SOILHC_METHOD==2)

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('HEATCON ',104)
      ENDIF
      RETURN
      END SUBROUTINE HEAT_CON
#endif
 
