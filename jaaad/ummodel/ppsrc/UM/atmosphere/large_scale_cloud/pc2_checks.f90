
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ PC2 Cloud Scheme: Checking cloud parameters
! Subroutine Interface:
      SUBROUTINE PC2_CHECKS(                                            &
!      Pressure related fields
     & p_theta_levels                                                   &
!      Array dimensions
     &,levels, row_length,rows                                          &
!      Prognostic Fields
     &,T, CF, CFL, CFF, Q, QCL, QCF                                     &
!      Logical control
     &,l_mixing_ratio)
!
      IMPLICIT NONE
!
! Purpose:
!   This subroutine checks that cloud fractions, liquid water and
!   vapour contents take on physically reasonable values.
!
! Method:
!   Apply checks sequentially to the input values. It is more important
!   to ensure that liquid does not go negative than to ensure that
!   saturation deficit is correct.
!
! Current Owner of Code: D. R. Wilson
!
! History:
! Version   Date     Comment
!  5.4    22-07-02   Original Code (Damian Wilson)
!  5.5    06-01-03   Adjust the treatment of ice when ice-cloud-
!                    fraction is zero. Damian Wilson
!  6.1    08-07-04   Apply a small minimum limit to the qcl and qcf
!                    in a cloud rather than a zero limit. Damian Wilson
!  6.4    18-08-06   Use mixing ratio formulation.  Damian Wilson
!
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: PC2 Cloud Scheme Documentation
!
!  Global Variables:----------------------------------------------------
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
! Description:
!
!  Contains various cloud droplet parameters, defined for
!  land and sea areas.
!
!  NTOT_* is the total number concentration (m-3) of cloud droplets;
!  KPARAM_* is the ratio of the cubes of the volume-mean radius
!                                           and the effective radius;
!  DCONRE_* is the effective radius (m) for deep convective clouds;
!  DEEP_CONVECTION_LIMIT_* is the threshold depth (m) bewteen shallow
!                                          and deep convective cloud.
!
! Current Code Owner: Andy Jones
!
! History:
!
! Version   Date     Comment
! -------   ----     -------
!    1     040894   Original code.    Andy Jones
!  5.2     111000   Updated in line with Bower et al. 1994 (J. Atmos.
!                   Sci., 51, 2722-2732) and subsequent pers. comms.
!                   Droplet concentrations now as used in HadAM4.
!                                     Andy Jones
!  5.4     02/09/02 Moved THOMO here from C_LSPMIC.      Damian Wilson
!  6.2     17/11/05 Remove variables that are now in UMUI. D. Wilson
!
!     REAL,PARAMETER:: NTOT_LAND is set in UMUI
!     REAL,PARAMETER:: NTOT_SEA is set in UMUI
      REAL,PARAMETER:: KPARAM_LAND = 0.67
      REAL,PARAMETER:: KPARAM_SEA = 0.80
      REAL,PARAMETER:: DCONRE_LAND = 9.5E-06
      REAL,PARAMETER:: DCONRE_SEA = 16.0E-06
      REAL,PARAMETER:: DEEP_CONVECTION_LIMIT_LAND = 500.0
      REAL,PARAMETER:: DEEP_CONVECTION_LIMIT_SEA = 1500.0
!
! Maximum Temp for homogenous nucleation (deg C)
      REAL,PARAMETER:: THOMO = -40.0
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

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
     & p_theta_levels(row_length,rows,levels)
!       pressure at all points (Pa)
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
     &,CFF(row_length,rows,LEVELS)                                      &
!       Ice cloud fraction (no units)
     &,Q(row_length,rows,LEVELS)                                        &
!       Vapour content (kg water per kg air)
     &,QCL(row_length,rows,LEVELS)                                      &
!       Liquid content (kg water per kg air)
     &,QCF(row_length,rows,LEVELS)
!       Ice content (kg water per kg air)
!
!  External functions:
!
!  Local parameters and other physical constants------------------------
      REAL                                                              &
     & LCRCP                                                            &
!       Latent heat of condensation divided by heat capacity of air.
     &,LSRCP                                                            &
!       Latent heat of sublimation divided by heat capacity of air.
     &,LFRCP                                                            &
!       Latent heat of freezing divided by heat capacity of air.
     &,ONE_OVER_QCF0                                                    &
!       One_over_qcf0 is the reciprocal of a reasonable ice content.
     &,CONDENSATE_LIMIT   ! Minimum value of condensate
!
      PARAMETER(                                                        &
     &          LCRCP = LC/CP                                           &
     &,         LSRCP =(LC+LF)/CP                                       &
     &,         LFRCP = LF/CP                                           &
     &,         ONE_OVER_QCF0 = 1.0E4                                   &
     &,         CONDENSATE_LIMIT = 1.0E-10                              &
     &          )
!
!  Local scalars--------------------------------------------------------
!
!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
      REAL                                                              &
     & ALPHA                                                            &
                  ! Rate of change of saturation specific humidity with
!                   temperature calculated at dry-bulb temperature
!                   (kg kg-1 K-1)
     &,AL                                                               &
                  ! 1 / (1 + alpha L/cp)  (no units)
     &,SD         ! Saturation deficit (= aL (q - qsat(T)) )  (kg kg-1)
!
!  (b) Others.
      INTEGER K,I,J       ! Loop counters: K - vertical level index
!                           I,J - horizontal position index
!
!  Local dynamic arrays-------------------------------------------------
!    1 block of real workspace is required.
      REAL                                                              &
     & QSL_T(row_length,rows)
!       Saturated specific humidity for dry bulb temperature T
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
! 1. Calculate Saturated Specific Humidity with respect to liquid water
!    for dry bulb temperatures.
! ----------------------------------------------------------------------
!
! DEPENDS ON: qsat_wat_mix
        call qsat_wat_mix(qsl_T,T(1,1,k),p_theta_levels(1,1,k),         &
     &                row_length*rows,l_mixing_ratio)
!
! Rows_do1:
        DO j=1,rows
! Row_length_do1:
          DO I=1,row_length
!
!----------------------------------------------------------------------
! 2. Calculate the saturation deficit.
! ----------------------------------------------------------------------
!
! Need to estimate the rate of change of saturated specific humidity
! with respect to temperature (alpha) first, then use this to calculate
! factor aL.
            ALPHA=EPSILON*LC*QSL_T(i,j)/(R*T(i,j,k)**2)
            AL=1.0/(1.0+LCRCP*ALPHA)
!
! Calculate the saturation deficit SD
!
            SD=AL*(QSL_T(i,j)-Q(i,j,k))
!
! ----------------------------------------------------------------------
!  3. Checks are applied here for liquid cloud
! ----------------------------------------------------------------------
!
! Earlier versions checked whether saturation deficit is zero (or less
! than zero). If so, then the liquid cloud fraction was forced to one.
! This check has been suspended for numerical reasons.
!            IF (SD  <=  0.0 .OR. CFL(i,j,k)  >   1.0) THEN
!
! Instead, check simply whether input values of liquid cloud fraction
! are, or are between, zero and one. If not, adjust them to zero or one.
! Adjust the total cloud fractions accordingly.
!
            IF (CFL(i,j,k)  >   1.0) THEN
              CFL(i,j,k)=1.0
              CF(i,j,k) =1.0
            END IF
!
! Check also whether the liquid water content is less than zero, and
! set liquid cloud fraction to zero if it is.
!
            IF (QCL(i,j,k)  <   CONDENSATE_LIMIT                        &
     &          .OR. CFL(i,j,k)  <   0.0) THEN
              CFL(i,j,k)=0.0
              CF(i,j,k) =CFF(i,j,k)
            END IF
!
! Check whether the saturation deficit is less than zero, or whether it
! is greater than zero but the liquid cloud fraction is one. If it is
! then condense or evaporate some liquid (provided there is enough
! liquid) to bring the saturation deficit to zero. Adjust the
! temperature for the latent heating.
!
            IF (  SD  <   0.0 .OR.                                      &
     &         (SD  >   0.0 .AND. CFL(i,j,k)  ==  1.0                   &
     &          .AND. QCL(i,j,k)  >   SD)  ) THEN
              Q(i,j,k)   = Q(i,j,k)   + SD
              QCL(i,j,k) = QCL(i,j,k) - SD
              T(i,j,k)   = T(i,j,k)   - SD * LCRCP
            END IF
!
! Check whether the liquid content is less than zero, or whether it is
! greater than zero but the liquid cloud fraction is zero. If so then
! condense or evaporate liquid to bring the liquid water to zero. Adjust
! the temperature for latent heating.
!
            IF (  QCL(i,j,k)  <   CONDENSATE_LIMIT .OR.                 &
     &         (QCL(i,j,k)  >   0.0 .AND. CFL(i,j,k)  ==  0.0)  ) THEN
              Q(i,j,k)   = Q(i,j,k) + QCL(i,j,k)
              T(i,j,k)   = T(i,j,k) - QCL(i,j,k) * LCRCP
              QCL(i,j,k) = 0.0
            END IF
!
! ----------------------------------------------------------------------
!  4. Check that ice content and ice cloud fraction are sensible.
! ----------------------------------------------------------------------
!
! Check whether ice content is zero (or less than zero). If so then
! force the ice cloud fraction to zero. Also check whether input values
! of ice cloud fraction are, or are between, zero and one. If not,
! adjust them to zero or one. Adjust the total cloud fractions
! accordingly.
!
            IF (CFF(i,j,k)  >   1.0) THEN
              CFF(i,j,k)=1.0
              CF(i,j,k) =1.0
            END IF
!
            IF (QCF(i,j,k)  <   CONDENSATE_LIMIT                        &
     &          .OR. CFF(i,j,k)  <   0.0) THEN
              CFF(i,j,k)=0.0
              CF(i,j,k) =CFL(i,j,k)
            END IF
!
! If ice content is negative then condense some vapour to remove the
! negative part. Adjust the temperature for the latent heat.
!
            IF (QCF(i,j,k)  <   CONDENSATE_LIMIT) THEN
              Q(i,j,k)   = Q(i,j,k) + QCF(i,j,k)
              T(i,j,k)   = T(i,j,k) - QCF(i,j,k) * LSRCP
              QCF(i,j,k) = 0.0
            END IF
!
! If ice content is positive but ice cloud fraction negative, create
! some ice cloud fraction.
!
            IF (QCF(i,j,k)  >   0.0 .AND. CFF(i,j,k)  ==  0.0) THEN
              CFF(i,j,k) = QCF(i,j,k) * ONE_OVER_QCF0
!             CF is not adjusted here but is checked below
            END IF
!
! ----------------------------------------------------------------------
!  5. Check that total cloud fraction is sensible.
! ----------------------------------------------------------------------
!
! Total cloud fraction must be bounded by
! i) The maximum of the ice and liquid cloud fractions (maximum overlap)
!
            IF (CF(i,j,k)  <   MAX(CFL(i,j,k),CFF(i,j,k))) THEN
              CF(i,j,k) = MAX(CFL(i,j,k),CFF(i,j,k))
            END IF
!
! ii) The sum of the ice and liquid cloud fractions or one, whichever
! is the lower (minimum overlap)
!
            IF (CF(i,j,k)  >   MIN( (CFL(i,j,k)+CFF(i,j,k)),1.0 ) ) THEN
              CF(i,j,k) = MIN( (CFL(i,j,k)+CFF(i,j,k)),1.0 )
            END IF
!
! ----------------------------------------------------------------------
!  6. Homogeneous nucleation.
! ----------------------------------------------------------------------
!
            IF (T(i,j,k)  <   (ZERODEGC+THOMO)) THEN
! Turn all liquid to ice
              IF (QCL(i,j,k)  >   0.0) THEN
                CFF(i,j,k) = CF(i,j,k)
                CFL(i,j,k) = 0.0
              END IF
              QCF(i,j,k) = QCF(i,j,k) + QCL(i,j,k)
              T(i,j,k)   = T(i,j,k)   + QCL(i,j,k) * LFRCP
              QCL(i,j,k) = 0.0
            END IF
!
! Row_length_do1:
          END DO
! Rows_do1:
        END DO
! Levels_do1:
      END DO
!
! End of the subroutine
!
      RETURN
      END SUBROUTINE PC2_CHECKS
