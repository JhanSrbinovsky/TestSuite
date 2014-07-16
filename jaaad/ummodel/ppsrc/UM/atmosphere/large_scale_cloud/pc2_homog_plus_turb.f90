
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ PC2 Cloud Scheme: Homogeneous forcing and Turbulence
! Subroutine Interface:
      SUBROUTINE PC2_HOMOG_PLUS_TURB(                                   &
!      Pressure related fields
     & p_theta_levels                                                   &
!      Array dimensions
     &,levels, row_length,rows                                          &
!      Timestep
     &,timestep                                                         &
!      Prognostic Fields
     &,T, CF, CFL, CFF, Q, QCL                                          &
!      Forcing quantities for driving the homogeneous forcing
     &,DTDT,DQDT,DLDT,DPDT                                              &
!      Other quantities for the turbulence
     &,DBSDTBS0,DBSDTBS1,l_mixing_ratio)
!
      IMPLICIT NONE
!
! Purpose:
!   This subroutine calculates the change in liquid content, liquid
!   cloud fraction and total cloud fraction as a result of homogeneous
!   forcing of the gridbox with temperature, pressure, vapour and liquid
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
!  5.4    22-07-02   Original Code (Damian Wilson)
!  6.1    09-07-02   Update comment lines (Damian Wilson)
!  6.4    08-01-07   Link erosion rate to RH (Damian Wilson)
!  6.4    18-08-06   Use mixing ratio formulation. Damian Wilson
!
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: PC2 Cloud Scheme Documentation
!
!  Global Variables:----------------------------------------------------
!
! Parameters required for the PC2 cloud scheme.
!
! Version    Date      Modification
!   5.4    16/08/02    Original Code.            Damian Wilson
!   5.5    13/03/03    Minor fix to defs. P.Selwood.
!   6.1    25/05/04    Change pdf_merge_power and pdf_power and
!                      turbulence in convection parameters.
!                                                Damian Wilson
!   6.2    28/2/05     Adjust cloud fraction tolerances
!                                                Damian Wilson
!   6.2    26/07/06    Change defs to include SCMA. R Barnes
!   6.4    02/02/07    Include parameter to control RH dependent
!                      PC2 erosion of cloud
!
      ! Number of iterations in the initiation
      INTEGER,PARAMETER:: INIT_ITERATIONS=10

      ! Tolerance of cloud fraction for resetting of cloud
      ! Cloud_pc2_tol_2 should always be less than cloud_pc2_tol
      REAL,PARAMETER:: CLOUD_PC2_TOL   = 0.005
      REAL,PARAMETER:: CLOUD_PC2_TOL_2 = 0.001

      ! Tolerance of critical relative humidity for initiation of cloud
      REAL,PARAMETER:: RHCRIT_TOL=0.01

      ! Power that is used to weight the two values of G when they are
      ! merged together
      REAL,PARAMETER:: PDF_MERGE_POWER=0.5

      ! Power that describes the way G varies with s near the boundary
      ! of the cloud probability density function. For a "top-hat" it
      ! is equal to 0, for a "triangular" distribution it is equal to 1.
      REAL,PARAMETER:: PDF_POWER=0.0

      ! Parameters that govern the turbulent decrease in width of the
      ! PDF.  (dbs/dt)/bs = (DBSDTBS_TURB_0 + dQc/dt DBSDTBS_TURB_1)
      !                 * exp( - dbsdtbs_exp Q_c / (a_L qsat(T_L)))
      ! dbsdtbs_turb_0 is now set in the UMUI
!     REAL,PARAMETER:: DBSDTBS_TURB_0 = -2.25E-5
      REAL,PARAMETER:: DBSDTBS_TURB_1 = 0.0
      REAL,PARAMETER:: DBSDTBS_CONV   = 0.0
      Real,Parameter:: dbsdtbs_exp    = 10.05

! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------
!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

!*----------------------------------------------------------------------
!
!  Subroutine Arguments:------------------------------------------------
      INTEGER                                                           &
                        !, INTENT(IN)
     & LEVELS                                                           &
!       No. of levels being processed.
     &, row_length,rows
!       Row length and number of rows being processed.
!
      REAL                                                              &
                        !, INTENT(IN)
     & p_theta_levels(row_length,rows,levels)                           &
!       pressure at all points (Pa)
     &,timestep                                                         &
!       Model timestep (s)
     &,CFF(row_length,rows,LEVELS)                                      &
!       Ice cloud fraction (no units)
     &,DTDT(row_length,rows,LEVELS)                                     &
!       Increment of temperature from forcing mechanism (K)
     &,DQDT(row_length,rows,LEVELS)                                     &
!       Increment of vapour from forcing mechanism (kg kg-1)
     &,DLDT(row_length,rows,LEVELS)                                     &
!       Increment of liquid from forcing mechanism (kg kg-1)
     &,DPDT(row_length,rows,LEVELS)                                     &
!       Increment in pressure from forcing mechanism (Pa)
     &,DBSDTBS0                                                         &
!       Value of dbs/dt / bs which is independent of forcing (no units)
     &,DBSDTBS1
!       Value of dbs/dt / bs which is proportional to the forcing
!       ( (kg kg-1 s-1)-1 )
!
      Logical                                                           &
                        !, INTENT(IN)
     & l_mixing_ratio   ! Use mixing ratio formulation
!
      REAL                                                              &
                        !, INTENT(INOUT)
     & T(row_length,rows,LEVELS)                                        &
!       Temperature (K)
     &,CF(row_length,rows,LEVELS)                                       &
!       Total cloud fraction (no units)
     &,CFL(row_length,rows,LEVELS)                                      &
!       Liquid cloud fraction (no units)
     &,Q(row_length,rows,LEVELS)                                        &
!       Vapour content (kg water per kg air)
     &,QCL(row_length,rows,LEVELS)
!       Liquid content (kg water per kg air)
!
!  External functions:
!
!  Local parameters and other physical constants------------------------
      REAL                                                              &
     & LCRCP                                                            &
!       Latent heat of condensation divided by heat capacity of air.
     &,B_FACTOR
!       Premultiplier to calculate the amplitude of the probability
!       density function at the saturation boundary (G).
!
      PARAMETER(                                                        &
     &          LCRCP=LC/CP                                             &
     &,         B_FACTOR=(PDF_POWER+1.0)/(PDF_POWER+2.0)                &
     &          )
!
!  Local scalars--------------------------------------------------------
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
     &,G                                                                &
                  ! Amplitude of the probability density function at
!                   the saturation boundary (kg kg-1)-1
     &,QC                                                               &
                  ! aL (q + l - qsat(TL) )  (kg kg-1)
     &,SD         ! Saturation deficit (= aL (q - qsat(T)) )  (kg kg-1)
!
!  (b) Others.
      INTEGER K,I,J                                                     &
                          ! Loop counters: K - vertical level index
!                           I,J - horizontal position index
     &,       NPT         ! Number of point on which to perform
!                           calculations
!
!  Local arrays---------------------------------------------------------
      REAL                                                              &
     & QSL_T(row_length,rows)                                           &
!       Saturated specific humidity for dry bulb temperature T
     &,QSL_TL(row_length,rows)                                          &
!       Saturated specific humidity for liquid temperature TL
     &,TL(row_length,rows)                                              &
!       Liquid temperature (= T - L/cp QCL)  (K)
     &,CF_C(row_length*rows)                                            &
!       Total cloud fraction on condensed points
     &,CFL_C(row_length*rows)                                           &
!       Liquid cloud fraction on condensed points
     &,CFF_C(row_length*rows)                                           &
!       Ice cloud fraction on condensed points
     &,DELTAC_C(row_length*rows)                                        &
!       Change in total cloud fraction (no units)
     &,DELTACL_C(row_length*rows)                                       &
!       Change in liquid cloud fraction (no units)
     &,DELTACF_C(row_length*rows)
!       Change in ice cloud fraction (no units)
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
            TL(i,j)=T(i,j,k)-LCRCP*QCL(i,j,k)
! Row_length_do1:
          END DO
! Rows_do1:
        END DO
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
        NPT=0
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
! Increment counter number
!
              NPT=NPT+1
!
! Write condensed cloud fraction variables
!
              CFL_C(NPT) = CFL(i,j,k)
              CF_C(NPT)  = CF(i,j,k)
              CFF_C(NPT) = CFF(i,j,k)
!
! ----------------------------------------------------------------------
! 3. Calculate the parameters relating to the probability density func.
! ----------------------------------------------------------------------
!
! Need to estimate the rate of change of saturated specific humidity
! with respect to temperature (alpha) first, then use this to calculate
! factor aL. Also estimate the rate of change of qsat with pressure.
              ALPHA=EPSILON*LC*QSL_T(i,j)/(R*T(i,j,k)**2)
              AL=1.0/(1.0+LCRCP*ALPHA)
              ALPHA_P = -QSL_T(i,j)/P_theta_levels(i,j,k)
!
! Calculate the saturation deficit SD
!
              SD=AL*(QSL_T(i,j)-Q(i,j,k))
!
! Calculate the amplitude of the probability density function at the
! saturation boundary.
!
              IF (QCL(i,j,k)  >   0.0 .AND. SD  >   0.0) THEN
!
                G=B_FACTOR*(   CFL(i,j,k)**PDF_MERGE_POWER              &
     &           *(1.0-CFL(i,j,k))**2/SD                                &
     &           +         (1.0-CFL(i,j,k))**PDF_MERGE_POWER            &
     &           *CFL(i,j,k)**2/QCL(i,j,k)  )                           &
     &           /(   CFL(i,j,k)**PDF_MERGE_POWER+(1.0-CFL(i,j,k))      &
     &           **PDF_MERGE_POWER   )
!
              ELSE
                G=0.0
              END IF
!
! Calculate the rate of change of Qc due to the forcing
!
              DQCDT=AL * ( DQDT(i,j,k)-ALPHA*DTDT(i,j,k)                &
     &                    -ALPHA_P*DPDT(i,j,k) ) + DLDT(i,j,k)
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
              DELTACL_C(NPT) = G * ( DQCDT - (QC + DQCDT) * DBSDTBS)
              DELTACF_C(NPT) = 0.0
!
! Calculate the condensation amount DELTAL. This uses a mid value
! of cloud fraction for better numerical behaviour.
!
              C_1 = CFL(i,j,k) + DELTACL_C(NPT)
              IF (C_1  >   1.0) THEN
                DELTACL_C(NPT) = 1.0 - CFL(i,j,k)
                C_1=1.0
              ELSE IF (C_1  <   0.0) THEN
                C_1=0.0
                DELTACL_C(NPT) = (- CFL(i,j,k) )
              END IF
                C_1 = 0.5 * (C_1 + CFL(i,j,k))
              DELTAL = C_1 * DQCDT + (QCL(i,j,k) - QC * C_1) * DBSDTBS
!
            ELSE IF (CFL(i,j,k)  ==  1.0) THEN
!
! Cloud fraction is 1
!
              ALPHA=EPSILON*LC*QSL_T(i,j)/(R*T(i,j,k)**2)
              AL=1.0/(1.0+LCRCP*ALPHA)
              ALPHA_P = -QSL_T(i,j)/P_theta_levels(i,j,k)
              DELTAL=AL * (DQDT(i,j,k)-ALPHA*DTDT(i,j,k)                &
     &                     -ALPHA_P*DPDT(i,j,k)) + DLDT(i,j,k)
!
            ELSE
!
! Cloud fraction is 0
!
              DELTAL = 0.0
!
            END IF
!
! Update water contents and temperature due to latent heating
!
            QCL(i,j,k) = QCL(i,j,k) + DELTAL
! Q = input Q + Forcing - Condensation
            Q(i,j,k)   = Q(i,j,k)   + DQDT(i,j,k)                       &
     &                  - (DELTAL - DLDT(i,j,k))
            T(i,j,k)   = T(i,j,k)   + DTDT(i,j,k)                       &
     &                    + LCRCP * (DELTAL - DLDT(i,j,k))
!
! Row_length_do2:
          END DO
! Rows_do2:
        END DO
!
! ----------------------------------------------------------------------
! 5. Now update cloud fractions.
! ----------------------------------------------------------------------
!
! Calculate change in total cloud fraction.
!
        IF (NPT  >   0) THEN
! DEPENDS ON: pc2_total_cf
          CALL PC2_TOTAL_CF(                                            &
     &          NPT,CFL_C,CFF_C,DELTACL_C,DELTACF_C,CF_C)
        END IF
        NPT=0
!
! Rows_do3:
        DO j=1,rows
! Row_length_do3:
          DO I=1,row_length
!
            IF (CFL(i,j,k)  >   0.0 .AND. CFL(i,j,k)  <   1.0) THEN
              NPT=NPT+1
!
! Update cloud fractions
!
              CF(i,j,k)  = CF_C(NPT)
              CFL(i,j,k) = CFL(i,j,k) + DELTACL_C(NPT)
!
! End if for CFL gt 0 and CFL lt 1
            END IF
! Row_length_do3:
          END DO
! Rows_do3:
        END DO
! Levels_do1:
      END DO
!
! End of the subroutine
!
      RETURN
      END SUBROUTINE PC2_HOMOG_PLUS_TURB
