#if defined(SCMA) || defined(RECON) || defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE FREEZE_SOIL -----------------------------------------

!
! Subroutine Interface:
      SUBROUTINE FREEZE_SOIL (NPNTS,NSHYD,B,DZ                          &
     &,                       SATHH,SMCL,TSOIL,V_SAT,STHU,STHF)


      IMPLICIT NONE
!
! Description:
!     Calculates the unfrozen and frozen water within a soil layer
!     as a fraction of saturation.                          (Cox, 6/95)
!
! Documentation : UM Documentation Paper 25
!
! Current Code Owner : David Gregory
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  3.4      6/95     Original code    Peter Cox
!  4.4  18/09/97     Rename from FREEZE to FREEZE_SOIL. D. Robinson.
!  6.2  14/06/06     Add CONTROL def block. Clive Jones
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! System component covered: P25
! System Task: P25
!


! Subroutine arguments
!   Scalar arguments with intent(IN) :
      INTEGER                                                           &
     & NPNTS                                                            &
                            ! IN Number of gridpoints.
     &,NSHYD                ! IN Number of soil layers.


!   Array arguments with intent(IN) :

      REAL                                                              &
     & B(NPNTS)                                                         &
                            ! IN Clapp-Hornberger exponent.
     &,DZ(NSHYD)                                                        &
                            ! IN Thicknesses of the soil layers (m).
     &,SATHH(NPNTS)                                                     &
                            ! IN Saturated soil water pressure (m).
     &,SMCL(NPNTS,NSHYD)                                                &
                            ! IN Soil moisture content of
                            !    layers (kg/m2).
     &,TSOIL(NPNTS,NSHYD)                                               &
                            ! IN Sub-surface temperatures (K).
     &,V_SAT(NPNTS)         ! IN Volumetric soil moisture
                            !    concentration at saturation
                            !    (m3 H2O/m3 soil).

!   Array arguments with intent(OUT) :
      REAL                                                              &
     & STHF(NPNTS,NSHYD)                                                &
                            ! OUT Frozen soil moisture content of
                            !     the layers as a fraction of
                            !     saturation.
     &,STHU(NPNTS,NSHYD)    ! OUT Unfrozen soil moisture content of
                            !     the layers as a fraction of
                            !     saturation.

! Local scalars:
      INTEGER                                                           &
     & I,J,N                ! WORK Loop counters.

! Local arrays:
      REAL                                                              &
     & SMCLF(NPNTS,NSHYD)                                               &
                            ! WORK Frozen moisture content of the
                            !      soil layers (kg/m2).
     &,SMCLU(NPNTS,NSHYD)                                               &
                            ! WORK Unfrozen moisture content of the
                            !      soil layers (kg/m2).
     &,SMCLSAT(NPNTS,NSHYD)                                             &
                            ! WORK The saturation moisture content of
                            !      the layers (kg/m2).
     &,TMAX(NPNTS)                                                      &
                            ! WORK Temperature above which all water is
                            !      unfrozen (Celsius)
     &,TSL(NPNTS,NSHYD)     ! WORK Soil layer temperatures (Celsius).

#include "c_densty.h"
#include "c_perma.h"
#include "c_0_dg_c.h"

      DO N=1,NSHYD

        DO I=1,NPNTS
!-----------------------------------------------------------------------
! Calculate TMAX, the temperature above which all soil water is
! unfrozen
!-----------------------------------------------------------------------
          SMCLSAT(I,N)=RHO_WATER*DZ(N)*V_SAT(I)
          TSL(I,N)=TSOIL(I,N)-ZERODEGC
          IF (SMCL(I,N) >  0.0) THEN
            TMAX(I)=-SATHH(I)/DPSIDT*(SMCLSAT(I,N)/SMCL(I,N))**(B(I))
          ELSE
            TMAX(I)=-273.15
          ENDIF

!--------------------------------------------------------------------
! Diagnose unfrozen and frozen water contents
!--------------------------------------------------------------------
          IF (TSL(I,N) >= TMAX(I)) THEN
            SMCLU(I,N)=SMCL(I,N)
            SMCLF(I,N)=0.0
          ELSE
!-----------------------------------------------------------------
! For ice points (V_SAT=0) set SMCLU=0.0 and SMCLF=0.0
!-----------------------------------------------------------------
            IF (V_SAT(I) == 0.0) THEN
              SMCLU(I,N)=0.0
              SMCLF(I,N)=0.0
            ELSE
              SMCLU(I,N)=SMCLSAT(I,N)                                   &
     &                    *(-DPSIDT*TSL(I,N)/SATHH(I))**(-1.0/B(I))
              SMCLF(I,N)=SMCL(I,N)-SMCLU(I,N)
            ENDIF
          ENDIF
          IF (SMCLSAT(I,N) >  0.0) THEN
            STHF(I,N)=SMCLF(I,N)/SMCLSAT(I,N)
            STHU(I,N)=SMCLU(I,N)/SMCLSAT(I,N)
          ELSE
            STHF(I,N)=0.0
            STHU(I,N)=0.0
          ENDIF
        ENDDO

      ENDDO

      RETURN
      END SUBROUTINE FREEZE_SOIL
#endif
