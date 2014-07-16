! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Subroutine CALC_ZW-------------------------------------------------
!
!   Purpose: To calculate the mean water table depth from the soil
!            moisture deficit as described in Koster et al., 2000.,
!            using the Newton-Raphson method.
!
! Documentation: UNIFIED MODEL DOCUMENTATION PAPER NO 25
!
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.


      SUBROUTINE CALC_ZW(NPNTS,NSHYD,SOIL_PTS,SOIL_INDEX                &
     &  ,BEXP,SATHH,SMCL,SMCLZW,SMCLSAT,SMCLSATZW,V_SAT,ZW)

      IMPLICIT NONE

! Global variables:
#include "c_densty.h"
#include "c_topog.h"

! Subroutine arguments
!   Scalar arguments with intent(IN) :
      INTEGER                                                           &
     & NPNTS                                                            &
                            ! IN Number of gridpoints.
     &,NSHYD                                                            &
                            ! IN Number of soil moisture levels.
     &,SOIL_PTS             ! IN Number of soil points.

!   Array arguments with intent(IN) :
      INTEGER                                                           &
     & SOIL_INDEX(NPNTS)    ! IN Array of soil points.

      REAL                                                              &
     & BEXP(NPNTS)                                                      &
                            ! IN Clapp-Hornberger exponent.
     &,SATHH(NPNTS)                                                     &
                            ! IN Saturated soil water pressure (m).
     &,SMCL(NPNTS,NSHYD)                                                &
                            ! IN Total soil moisture contents
!                           !      of each layer (kg/m2).
     &,SMCLSAT(NPNTS,NSHYD)                                             &
                            ! IN Soil moisture contents of each
!                           !      layer at saturation (kg/m2).
     &,SMCLZW(NPNTS)                                                    &
                            ! IN moisture content in deep layer.
     &,SMCLSATZW(NPNTS)                                                 &
                            ! IN moisture content in deep layer
!                           !      at saturation.
     &,V_SAT(NPNTS)         ! IN Volumetric soil moisture.

!
!   Array arguments with intent(INOUT) :
      REAL                                                              &
     & ZW(NPNTS)            ! INOUT Water table depth (m).

      REAL                                                              &
     & ZW_OLD                                                           &
                  ! WORK zw from last timestep
     &,FN                                                               &
                  ! WORK
     &,DFN                                                              &
                  ! WORK
     &,ZWEST                                                            &
                  ! WORK
     &,SMD                                                              &
                  ! WORK soil moisture deficit
     &,PSISAT                                                           &
                  ! WORK
     &,DIFFINTER  ! WORK

      INTEGER                                                           &
     & I,J,N                                                            &
                  ! Loop counters
     &,IT                                                               &
                  ! Loop counters
     &,NITER                                                            &
                  ! Number of iterations
     &,MAXINTER
      PARAMETER(NITER=3)

      DO J=1,SOIL_PTS
      I=SOIL_INDEX(J)
        ZW_OLD=ZW(I)

! Calculate soil moisture deficit:
        SMD=0.0
        DO N=1,NSHYD
           SMD=SMD+(SMCLSAT(I,N)-SMCL(I,N))/RHO_WATER
        ENDDO
        SMD=SMD+(SMCLSATZW(I)-SMCLZW(I))/RHO_WATER
        PSISAT=-SATHH(I)
        ZWEST=ZW_OLD

        MAXINTER=0
        DO IT=1,NITER
           IF(ZWEST <  0)ZWEST=0.0
           IF(ZWEST >  ZW_MAX)ZWEST=ZW_MAX
           ZW_OLD=ZWEST
           ZWEST=ZW_OLD

! Newton-Raphson. zw(next)=zw-f(zw)/f'(zw).

           FN = ZWEST*V_SAT(I) - SMD -                                  &
     &        V_SAT(I) * BEXP(I)/(BEXP(I)-1.0) * PSISAT *               &
     &        ( 1.0 - ( 1.0-ZWEST/PSISAT ) **                           &
     &                ( 1.0 - 1.0/BEXP(I) ) )       !   f(zw)
           DFN =  V_SAT(I) +V_SAT(I) *                                  &
     &        ( 1.0-ZWEST/PSISAT ) ** (-1.0/BEXP(I))

           IF(ABS(DFN) >  1.E-10)ZWEST=ZWEST-FN/DFN
           IF(ZWEST <  0)ZWEST=0.0
           IF(ZWEST >  ZW_MAX)ZWEST=ZW_MAX
           DIFFINTER=ABS(ZW_OLD-ZWEST)
        ENDDO

        IF(DIFFINTER >  0.01)MAXINTER=MAXINTER+1
        ZW(I)=ZWEST
        IF(ZW(I) <  0)ZW(I)=0.0
        IF(ZW(I) >  ZW_MAX)ZW(I)=ZW_MAX
      ENDDO

      END SUBROUTINE CALC_ZW
