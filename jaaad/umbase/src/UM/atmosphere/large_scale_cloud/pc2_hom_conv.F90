#if defined(A09_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Cloud Scheme: Homogenous forcing and Turbulence (non-updating)
! Subroutine Interface:
      SUBROUTINE PC2_HOM_CONV(                                          &
!      Pressure related fields
     & p_theta_levels                                                   &
!      Array dimensions
     &,levels, row_length, rows                                         &
!      Timestep
     &,timestep                                                         &
!      Prognostic Fields
     &,T, Q, QCL, CF, CFL, CFF                                          &
!      Forcing quantities for driving the homogenous forcing
     &,DTIN, DQIN, DQCLIN, DPDT, DCFLIN                                 &
!      Output increments to the prognostic fields
     &,DTPC2, DQPC2, DQCLPC2, DCFPC2, DCFLPC2                           &
!      Other quantities for the turbulence
     &,DBSDTBS0,DBSDTBS1,l_mixing_ratio)
!
      IMPLICIT NONE
!
! Description:
!   This subroutine calculates the change in liquid content, liquid
!   cloud fraction and total cloud fraction as a result of homogenous
!   forcing of the gridbox with temperature, pressue, vapour and liquid
!   increments and turbulent mixing effects.
!
! Method:
!   Uses the method in Gregory et al (2002, QJRMS 128 1485-1504) and
!   Wilson and Gregory (2003, QJRMS 129 967-986)
!   which considers a probability density distribution whose
!   properties are only influenced by a change of width due to
!   turbulence. There is NO checking that the condensate output
!   values from this routine are sensible.
!
! Current Owner of Code: Section 09 Code Owner
!
! History:
! Version   Date     Comment
!  6.1    25-05-04   Original code based upon PC2_HOMOG_PLUS_TURB.
!                                                     Damian Wilson
!  6.2    15/08/05   Free format fixes. P.Selwood
!  6.4    29-01-07   Link erosion rate to RH (Damian Wilson)
!  6.4    18/08/06   Use mixing ratios. Damian Wilson
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: PC2 cloud scheme documentation
!
!  Global Variables:----------------------------------------------------
#include "c_pc2pdf.h"
#include "c_lheat.h"
#include "c_r_cp.h"
#include "c_epslon.h"
!
!  Subroutine Arguments:------------------------------------------------
!
! arguments with intent in. ie: input variables.
!
      INTEGER                                                           &
                        !, INTENT(IN)
     & LEVELS                                                           &
!       No. of levels being processed.
     &, row_length,rows
!       Row length and number of rows being processed.
!
      REAL                                                              &
                        !, INTENT(IN)
     & timestep                                                         &
!       Model timestep (s)
     &,DBSDTBS0                                                         &
!       Value of dbs/dt / bs which is independent of forcing (no units)
     &,DBSDTBS1
!       Value of dbs/dt / bs which is proportional to the forcing
!       ( (kg kg-1 s-1)-1 )
!
      REAL                                                              &
                        !, INTENT(IN)
     & p_theta_levels(row_length,rows,levels)                           &
!       pressure at all points (Pa)
     &,T(row_length,rows,LEVELS)                                        &
!       Temperature (K)
     &,Q(row_length,rows,LEVELS)                                        &
!       Vapour content (kg water per kg air)
     &,QCL(row_length,rows,LEVELS)                                      &
!       Liquid content (kg water per kg air)
     &,CF(row_length,rows,LEVELS)                                       &
!       Total cloud fraction (no units)
     &,CFL(row_length,rows,LEVELS)                                      &
!       Liquid cloud fraction (no units)
     &,CFF(row_length,rows,LEVELS)                                      &
!       Ice cloud fraction (no units)
     &,DTIN(row_length,rows,LEVELS)                                     &
!       Increment of temperature from forcing mechanism (K)
     &,DQIN(row_length,rows,LEVELS)                                     &
!       Increment of vapour from forcing mechanism (kg kg-1)
     &,DQCLIN(row_length,rows,LEVELS)                                   &
!       Increment of liquid from forcing mechanism (kg kg-1)
     &,DPDT(row_length,rows,LEVELS)                                     &
!       Increment in pressure from forcing mechanism (Pa)
     &,DCFLIN(row_length,rows,LEVELS)
!       Increment in liquid cloud fraction (no units)
!
      Logical                                                           &
                        !, INTENT(IN)
     & l_mixing_ratio   ! Use mixing ratio formulation
!
! arguments with intent out. ie: output variables.
!
      REAL                                                              &
                        !, INTENT(OUT)
     & DTPC2(row_length,rows,LEVELS)                                    &
!       PC2 Increment to Temperature (K)
     &,DQPC2(row_length,rows,LEVELS)                                    &
!       PC2 Increment to Vapour content (kg water per kg air)
     &,DQCLPC2(row_length,rows,LEVELS)                                  &
!       PC2 Increment to Liquid content (kg water per kg air)
     &,DCFPC2(row_length,rows,LEVELS)                                   &
!       PC2 Increment to Total cloud fraction (no units)
     &,DCFLPC2(row_length,rows,LEVELS)
!       PC2 Increment to Liquid cloud fraction (no units)
!
!  External subroutine calls: ------------------------------------------
!
!  Local parameters and other physical constants------------------------
      REAL                                                              &
     & LCRCP                                                            &
!       Latent heat of condensation divided by heat capacity of air.
     &,B_FACTOR
!       Premultiplier to calculate the amplitude of the probability
!       density function at the saturation boundary (G_MQC).
!
      PARAMETER(                                                        &
     &          LCRCP=LC/CP                                             &
     &,         B_FACTOR=(PDF_POWER+1.0)/(PDF_POWER+2.0)                &
     &          )
!
!  Local scalars--------------------------------------------------------
!
      REAL                                                              &
     & ALPHA                                                            &
                  ! Rate of change of saturation specific humidity with
!                   temperature calculated at dry-bulb temperature
!                   (kg kg-1 K-1)
     &,ALPHA_P                                                          &
                  ! Rate of change of saturation specific humidity with
!                   pressure calculated at dry-bulb temperature (Pa K-1)
     &,AL                                                               &
                  ! 1 / (1 + alpha L/cp)  (no units)
     &,C_1                                                              &
                  ! Mid-timestep liquid cloud fraction (no units)
     &,DBSDTBS                                                          &
                  ! Relative rate of change of distribution width (s-1)
     &,DQCDT                                                            &
                  ! Forcing of QC (kg kg-1 s-1)
     &,DELTAL                                                           &
                  ! Change in liquid content (kg kg-1)
     &,CFL_TO_M                                                         &
                  ! CFL(i,j,k)**PDF_MERGE_POWER
     &,SKY_TO_M                                                         &
                  ! (1-CFL(i,j,k))**PDF_MERGE_POWER
     &,G_MQC                                                            &
                  ! Amplitude of the probability density function at
!                   the saturation boundary (kg kg-1)-1
     &,QC                                                               &
                  ! aL (q + l - qsat(TL) )  (kg kg-1)
     &,SD                                                               &
                  ! Saturation deficit (= aL (q - qsat(T)) )  (kg kg-1)
     &,DCS        ! Injected cloud fraction
!
!  (b) Others.
      INTEGER K,I,J       ! Loop counters: K - vertical level index
!                           I,J - horizontal position index
!
!  Local arrays---------------------------------------------------------
!
      REAL                                                              &
     & QSL_T(row_length,rows)                                           &
!       Saturated specific humidity for dry bulb temperature T
     &,QSL_TL(row_length,rows)                                          &
!       Saturated specific humidity for liquid temperature TL
     &,TL(row_length,rows)
!       Liquid temperature (= T - L/cp QCL)  (K)
!
!- End of Header
!
! ==Main Block==--------------------------------------------------------
!
! Loop round levels to be processed
! Levels_do1:
      DO k=LEVELS,1,-1
!
! ----------------------------------------------------------------------
! 1. Calculate liquid water temperature TL
! ----------------------------------------------------------------------
!
! Rows_do1:
        DO j=1,rows
! Row_length_do1:
          DO I=1,row_length
            TL(i,j) = T(i,j,k) - LCRCP*QCL(i,j,k)
          END DO  ! Row_length_do1
        END DO  ! Rows_do1
!
! ----------------------------------------------------------------------
! 2. Calculate Saturated Specific Humidity with respect to liquid water
!    for both dry bulb and wet bulb temperatures.
! ----------------------------------------------------------------------
!
! DEPENDS ON: qsat_wat_mix
        call qsat_wat_mix(qsl_T,T(1,1,k),p_theta_levels(1,1,k),         &
     &                row_length*rows,l_mixing_ratio)
! DEPENDS ON: qsat_wat_mix
        call qsat_wat_mix(qsl_TL,TL,p_theta_levels(1,1,k),              &
     &                row_length*rows,l_mixing_ratio)
!
! Rows_do2:
        DO j=1,rows
! Row_length_do2:
          DO I=1,row_length
!
! There is no need to perform the total cloud fraction calculation in
! this subroutine if there is no, or full, liquid cloud cover.
!
            IF (CFL(i,j,k)  >   0.0 .AND. CFL(i,j,k)  <   1.0) THEN
!
! ----------------------------------------------------------------------
! 3. Calculate the parameters relating to the probability density func.
! ----------------------------------------------------------------------
!
! Need to estimate the rate of change of saturated specific humidity
! with respect to temperature (alpha) first, then use this to calculate
! factor aL. Also estimate the rate of change of qsat with pressure.
              ALPHA = EPSILON*LC*QSL_T(i,j) / (R*T(i,j,k)**2)
              AL = 1.0 / (1.0 + LCRCP*ALPHA)
              ALPHA_P = -QSL_T(i,j) / P_theta_levels(i,j,k)
!
! Calculate the saturation deficit SD
!
              SD = AL*(QSL_T(i,j)-Q(i,j,k))
!
! Calculate the amplitude of the probability density function at the
! saturation boundary.
!
              IF (QCL(i,j,k)  >   0.0 .AND. SD  >   0.0) THEN
!
                CFL_TO_M = CFL(i,j,k)**PDF_MERGE_POWER
                SKY_TO_M = (1.0 - CFL(i,j,k))**PDF_MERGE_POWER
!
                G_MQC = B_FACTOR * ( (1.0-CFL(i,j,k))**2 *              &
     &           CFL_TO_M / (SD * (CFL_TO_M + SKY_TO_M))                &
     &                +                   CFL(i,j,k)**2 *               &
     &           SKY_TO_M / (QCL(i,j,k) * (CFL_TO_M + SKY_TO_M)) )
!
              ELSE
                G_MQC = 0.0
              END IF
!
! Calculate the rate of change of Qc due to the forcing
!
             DQCDT = AL * ( DQIN(i,j,k) - ALPHA*DTIN(i,j,k)             &
     &              -ALPHA_P*DPDT(i,j,k) ) + DQCLIN(i,j,k)
!
! For the background homogeneous forcing from the convection there is
! an additional term because the detrained plume must be saturated.
! This can also be written as a forcing.
!
!              IF (CFL(i,j,k)  >   0.0 .AND. CFL(i,j,k)  <   1.0) THEN
! This if test is already guaranteed
              IF (DCFLIN(i,j,k)  >   0.0) THEN
                DCS = DCFLIN(i,j,k) / (1.0 - CFL(i,j,k))
              ELSE IF (DCFLIN(i,j,k)  <   0.0) THEN
                DCS = DCFLIN(i,j,k) / (    - CFL(i,j,k))
              ELSE
                DCS = 0.0
              ENDIF
! Limit DCS to 0 and 1
              DCS = MAX(  MIN(DCS,1.0)  ,0.0)
!
              DQCDT = DQCDT - AL * DCS * (QSL_T(i,j)-Q(i,j,k))
!
! Calculate Qc
!
              QC = AL * (Q(i,j,k) + QCL(i,j,k) - QSL_TL(i,j))
!
! Calculate the relative rate of change of width of the distribution
! dbsdtbs from the forcing rate
!
              DBSDTBS = (DBSDTBS0 * TIMESTEP + DQCDT * DBSDTBS1) *      &
     &           exp(-dbsdtbs_exp * QC / (aL * qsl_tl(i,j)))
!
! ----------------------------------------------------------------------
! 4. Calculate the change of liquid cloud fraction. This uses the
! arrival value of QC for better behaved numerics.
! ----------------------------------------------------------------------
!
! DQCDT is the homogeneous forcing part, (QC+DQCDT)*DBSDTBS is the
! width narrowing part
!
              DCFLPC2(i,j,k) = G_MQC * ( DQCDT - (QC + DQCDT)*DBSDTBS)
!
! Calculate the condensation amount DELTAL. This uses a mid value
! of cloud fraction for better numerical behaviour.
!
              C_1 = MAX( 0., MIN( (CFL(i,j,k) + DCFLPC2(i,j,k)), 1.) )
!
              DCFLPC2(i,j,k) = C_1 - CFL(i,j,k)
!
              C_1 = 0.5 * (C_1 + CFL(i,j,k))
              DELTAL = C_1 * DQCDT +                                    &
     &        MAX( (QCL(i,j,k) - QC * C_1) * DBSDTBS , (-QCL(i,j,k)) )
!
! ----------------------------------------------------------------------
! 5. Calculate change in total cloud fraction.
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! The following If test is a copy of the PC2_TOTAL_CF subroutine.
! ----------------------------------------------------------------------
! PC2_TCF_if1:
              IF (DCFLPC2(i,j,k)  >   0.0) THEN
! ...  .AND. CFL(i,j,k)  <   1.0 already assured.
                DCFPC2(i,j,k) = DCFLPC2(i,j,k) * (1.0 - CF(i,j,k)) /    &
     &                                           (1.0 - CFL(i,j,k))
              ELSE IF (DCFLPC2(i,j,k)  <   0.0) THEN
! ...  .AND. CFL(i,j,k)  >   0.0 already assured.
                DCFPC2(i,j,k) = DCFLPC2(i,j,k) *                        &
     &                  (CF(i,j,k)- CFF(i,j,k)) / CFL(i,j,k)
              ELSE
                DCFPC2(i,j,k) = 0.0
              ENDIF ! PC2_TCF_if1
! ----------------------------------------------------------------------
!
            ELSE IF (CFL(i,j,k)  ==  1.0) THEN
!
! Cloud fraction is 1
!
              DCFPC2(i,j,k)  = 0.0
              DCFLPC2(i,j,k) = 0.0
!
              ALPHA = EPSILON * LC * QSL_T(i,j) / (R * T(i,j,k)**2)
              AL = 1.0 / (1.0 + LCRCP*ALPHA)
              ALPHA_P = -QSL_T(i,j) / P_theta_levels(i,j,k)
              DELTAL = AL * (DQIN(i,j,k) - ALPHA*DTIN(i,j,k)            &
     &              -ALPHA_P*DPDT(i,j,k)) + DQCLIN(i,j,k)
!
            ELSE
!
! Cloud fraction is 0
!
              DCFPC2(i,j,k)  = 0.0
              DCFLPC2(i,j,k) = 0.0
!
              DELTAL = 0.0
!
            END IF
!
! Increment water contents and temperature due to latent heating.
! This subroutine will output only the condensation increments
! hence we comment out updates to qcl, q and t
!           QCL(i,j,k) = QCL(i,j,k) + DELTAL
! Q = input Q + Forcing - Condensation
!           Q(i,j,k)   = Q(i,j,k) + DQIN(i,j,k) - (DELTAL - DLIN(i,j,k))
!           T(i,j,k)   = T(i,j,k) + DTIN(i,j,k)
!    &                            + LCRCP * (DELTAL - DLIN(i,j,k))
!
! These are the condensation increments
            DQCLPC2(i,j,k) = DELTAL - DQCLIN(i,j,k)
            DQPC2(i,j,k)   = - DQCLPC2(i,j,k)
            DTPC2(i,j,k)   = LCRCP * DQCLPC2(i,j,k)
!
          END DO  ! Row_length_do2
        END DO  ! Rows_do2
!
      END DO  ! Levels_do1
!
! End of the subroutine
!
      RETURN
      END SUBROUTINE PC2_HOM_CONV
#endif
