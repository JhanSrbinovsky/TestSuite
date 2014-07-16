
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ PC2 Cloud Scheme: Reset clouds for extreme relative total humidity.
! Subroutine Interface:
      SUBROUTINE PC2_CHECKS2(                                           &
!      Pressure related fields
     & p_theta_levels, RHCRIT                                           &
!      Array dimensions
     &,levels, row_length,rows                                          &
     &,rhc_row_length,rhc_rows                                          &
!      Prognostic Fields
     &,T, CF, CFL, CFF, Q, QCL                                          &
!      Logical control
     &,l_mixing_ratio)
!
      IMPLICIT NONE
!
! Purpose:
!   Check that cloud fraction is either zero or one when relative
!   total humidity is small or large.
!
! Method:
!   Calculate relative total humidity, compare to RHcrit and adjust
!   the cloud variables appropriately.
!
! Current Owner of Code: D. R. Wilson
!
! History:
! Version   Date     Comment
!  5.4    22-07-02   Original Code (Damian Wilson)
!  6.1    17-11-04   Correct indexing and sign errors (Damian Wilson)
!  6.2    28-02-05   Include second set of thresholds (Damian Wilson)
!  6.4    18-08-06   Use mixing ratio formulation.  Damian Wilson
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
     &, row_length,rows                                                 &
!       Row length and number of rows being processed.
     &, rhc_row_length,rhc_rows
!       Dimensions of the RHCRIT variable.
!
      REAL                                                              &
                        !, INTENT(IN)
     & p_theta_levels(row_length,rows,levels)                           &
!       pressure at all points (Pa)
     &,RHCRIT(rhc_row_length,rhc_rows,LEVELS)                           &
!       Critical relative humidity.  See the the paragraph incorporating
!       eqs P292.11 to P292.14; the values need to be tuned for the give
!       set of levels.
     &,CFF(row_length,rows,LEVELS)
!       Ice cloud fraction (no units)
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
     &,C_THRESH_LOW,C_THRESH_LOW_2                                      &
!       Low cloud fraction thresholds
     &,C_THRESH_HIGH,C_THRESH_HIGH_2
!       High cloud fraction thresholds
!
      PARAMETER(                                                        &
     &          LCRCP=LC/CP                                             &
     &,         C_THRESH_LOW    =       CLOUD_PC2_TOL                   &
     &,         C_THRESH_HIGH   = 1.0 - CLOUD_PC2_TOL                   &
     &,         C_THRESH_LOW_2  =       CLOUD_PC2_TOL_2                 &
     &,         C_THRESH_HIGH_2 = 1.0 - CLOUD_PC2_TOL_2                 &
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
     &,RHT                                                              &
                  ! Relative total humidity
     &,SD         ! Saturation deficit
!
!  (b) Others.
      INTEGER K,I,J                                                     &
                        ! Loop counters: K - vertical level index
!                         I,J - horizontal position index
     &,       IRHI,IRHJ                                                 &
                        ! Indices for RHcrit array
     &,       MULTRHC   ! Zero if (rhc_row_length*rhc_rows) le 1, else 1
!
!  Local dynamic arrays-------------------------------------------------
!    3 blocks of real workspace are required.
      REAL                                                              &
     & QSL_T(row_length,rows)                                           &
!       Saturated specific humidity for temperature T
     &,QSL_TL(row_length,rows)                                          &
!       Saturated specific humidity for liquid temperature TL
     &,TL(row_length,rows)
!       Liquid temperature (= T - L/cp QCL)  (K)
!
!- End of Header
!
! Set up a flag to state whether RHcrit is a single parameter or defined
! on all points.
!
      IF (rhc_row_length*rhc_rows  >   1) THEN
        MULTRHC=1
      ELSE
        MULTRHC=0
      ENDIF
!
! ==Main Block==--------------------------------------------------------
!
! Loop round levels to be processed
! Levels_do1:
      DO k=LEVELS,1,-1
!
! ----------------------------------------------------------------------
! 1. Calculate TL
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
!    for dry bulb temperature and liquid temperature.
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
          DO i=1,row_length
!

! Set up index pointers to critical relative humidity value
!
            IRHI = (MULTRHC * (i - 1)) + 1
            IRHJ = (MULTRHC * (j - 1)) + 1
!
! Calculate relative total humidity with respect to the liquid
! temperature and threshold relative humidity
!
            RHT=(Q(i,j,k)+QCL(i,j,k))/QSL_TL(i,j)
!
! ----------------------------------------------------------------------
! 3. Determine whether resetting is required and calculate appropriate
!    increments if this is the case.
! ----------------------------------------------------------------------
!
! Is the relative total humidity greater than the threshold amount?
! If so, evaporate some liquid to set the saturation deficit to zero.
!
            IF ( (RHT  >   (2.0-RHCRIT(IRHI,IRHJ,k)) .AND.              &
     &        (CFL(i,j,k)  >=  C_THRESH_HIGH) ) .OR.                    &
     &        (CFL(i,j,k)  >=  C_THRESH_HIGH_2) ) THEN
              CFL(i,j,k) = 1.0
              CF(i,j,k)  = 1.0
!
! Calculate the saturation deficit
!
              ALPHA=EPSILON*LC*QSL_T(i,j)/(R*T(i,j,k)**2)
              AL=1.0/(1.0+LCRCP*ALPHA)
              SD=AL*(QSL_T(i,j)-Q(i,j,k))
!
! Update the water contents
!
              QCL(i,j,k) = QCL(i,j,k) - SD
              Q(i,j,k)   = Q(i,j,k)   + SD
              T(i,j,k)   = T(i,j,k)   - SD * LCRCP
!
            END IF
!
! Is the relative total humidity less than the threshold amount?
! If so, evaporate all the liquid.
!
            IF ( (RHT  <   (RHCRIT(IRHI,IRHJ,k)) .AND.                  &
     &        (CFL(i,j,k)  <=  C_THRESH_LOW) )  .OR.                    &
     &        (CFL(i,j,k)  <=  C_THRESH_LOW_2)  ) THEN
              CFL(i,j,k) = 0.0
              CF(i,j,k)  = CFF(i,j,k)
              Q(i,j,k)   = Q(i,j,k) + QCL(i,j,k)
              T(i,j,k)   = T(i,j,k) - QCL(i,j,k) * LCRCP
              QCL(i,j,k) = 0.0
            END IF
!
! Row_length_do2:
          END DO
! Rows_do2:
        END DO
! Levels_do1:
      END DO
!
! End of the subroutine
!
      RETURN
      END SUBROUTINE PC2_CHECKS2
