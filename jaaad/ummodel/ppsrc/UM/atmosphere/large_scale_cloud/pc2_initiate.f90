
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ PC2 Cloud Scheme: Initiation
! Subroutine Interface:
      SUBROUTINE PC2_INITIATE(                                          &
!      Pressure related fields
     & p_theta_levels, ccb, cumulus, rhcrit                             &
!      Array dimensions
     &,levels, row_length,rows                                          &
     &,rhc_row_length,rhc_rows                                          &
!      Prognostic Fields
     &,  T, CF, CFL, CFF, Q, QCL, RHTS                                  &
!      Logical control
     &,  l_mixing_ratio)
!
      IMPLICIT NONE
!
! Purpose:
!   Initiate liquid and total cloud fraction and liquid water content
!
! Method:
!   Uses the method proposed in Annex C of the PC2 cloud scheme project
!   report, which considers a Smith-like probability density
!   distribution whose width is given by a prescribed value of RHcrit.
!
! Current Owner of Code: D. R. Wilson
!
! History:
! Version   Date     Comment
!  5.4    22-07-02   Original Code (Damian Wilson)
!  6.1    25-05-04   Restrict initiation if boundary layer convection
!                    is present. Damian Wilson
!  6.2    15-08-05   Free format fixes. P.Selwood
!  6.2    31-03-05   Extra comparison of RH_T with RH(TL)  Damian Wilson
!  6.4    18-08-06   Use mixing ratio formulation. Damian Wilson
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
     &, row_length,rows                                                 &
!       Row length and number of rows being processed.
     &, rhc_row_length,rhc_rows                                         &
!       Dimensions of the RHCRIT variable.
     &,ccb(row_length,rows)
!       convective cloud base
!
      LOGICAL                                                           &
                        !, INTENT(IN)
     & cumulus(row_length,rows)                                         &
!       Is this a boundary layer cumulus point
     &,l_mixing_ratio  ! Use mixing ratio formulation
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
     &,QCL(row_length,rows,LEVELS)                                      &
!       Liquid content (kg water per kg air)
     &,RHTS(row_length,rows,LEVELS)
!       Variable carrying initial RHT wrt TL from start of timestep
!
!  External functions:
!
!  Local parameters and other physical constants------------------------
      REAL                                                              &
     & LCRCP
!       Latent heat of condensation divided by heat capacity of air.
!
      PARAMETER(                                                        &
     &          LCRCP=LC/CP                                             &
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
     &,BS                                                               &
                  ! Width of distribution (kg kg-1)
     &,C_1                                                              &
                  ! New liquid cloud fraction
     &,DELTAL                                                           &
                  ! Change in liquid (kg kg-1)
     &,DESCENT_FACTOR                                                   &
                      ! Rate at which to relax to new liquid water
!                      content
     &,L_BS                                                             &
                  ! New liquid water content divided by PDF width
!                   (no units)
     &,L_OUT                                                            &
                  ! New liquid water content (kg kg-1)
     &,Q_OUT                                                            &
                  ! New vapour content (kg kg-1)
     &,QC                                                               &
                  ! aL (q + l - qsat(TL) )  (kg kg-1)
     &,QN                                                               &
                  ! Normalized value of QC for the probability density
!                   function.
     &,RHT                                                              &
                  ! Total relative humidity (liquid+vapour)/qsat
!                  (kg kg-1)
     &,RH0                                                              &
                  ! Equivalent critical relative humidity
     &,SD         ! Saturation deficit (= aL (q - qsat(T)) )  (kg kg-1)
!
!  (b) Others.
      INTEGER K,I,J,L                                                   &
                        ! Loop counters: K - vertical level index
!                         I,J - horizontal position index
!                         L   - counter for iterations
     &,       IRHI,IRHJ                                                 &
                        ! Indices for RHcrit array
     &,       MULTRHC                                                   &
                        ! Zero if (rhc_row_length*rhc_rows) le 1, else 1
     &,       NPTI      ! Number of points to iterate over
!
!  Local dynamic arrays-------------------------------------------------
!    18 blocks of real workspace are required.
      REAL                                                              &
     & QSL_T(row_length,rows)                                           &
!       Saturated specific humidity for dry bulb temperature T
     &,QSL_TL(row_length,rows)                                          &
!       Saturated specific humidity for liquid temperature TL
     &,TL(row_length,rows)                                              &
!       Liquid temperature (= T - L/cp QCL)  (K)
     &,CF_C(row_length*rows)                                            &
!       Total cloud fraction on compressed points (no units)
     &,CFL_C(row_length*rows)                                           &
!       Liquid cloud fraction on compressed points (no units)
     &,CFF_C(row_length*rows)                                           &
!       Ice cloud fraction on compressed points (no units)
     &,DELTACL_C(row_length*rows)                                       &
!       Change in liquid cloud fraction (no units)
     &,DELTACF_C(row_length*rows)                                       &
!       Change in ice cloud fraction (no units)
     &,p_theta_levels_c(row_length*rows)                                &
!       Pressure on condensed points (Pa)
     &,Q_C(row_length*rows)                                             &
!       Vapour content on compressed points (kg kg-1)
     &,QCL_C(row_length*rows)                                           &
!       Liquid water content on compressed points (kg kg-1)
     &,QN_C(row_length*rows)                                            &
!      QN on compressed points (no units)
     &,QSL_T_C(row_length*rows)                                         &
!       Saturated specific humidity for dry-bulb temperature T on
!       compressed points (kg kg-1)
     &,QSL_TL_C(row_length*rows)                                        &
!       Saturated specific humidity for liquid temperature TL on
!       compressed points (kg kg-1)
     &,RH0_C(row_length*rows)                                           &
!       Equivalent critical relative humidity on compressed points
!       (no units)
     &,T_C(row_length*rows)
!       Temperature on compressed points (K)
!
      INTEGER                                                           &
     & NI(row_length*rows)                                              &
!       Condensed point counter
     &,NJ(row_length*rows)
!       Condensed point counter
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
          DO i=1,row_length
            TL(i,j)=T(i,j,k)-LCRCP*QCL(i,j,k)
! Row_length_do1:
          END DO
! Rows_do1:
        END DO
!
! ----------------------------------------------------------------------
! 2. Calculate Saturated Specific Humidity with respect to liquid water
!    for liquid temperatures.
! ----------------------------------------------------------------------
!
! DEPENDS ON: qsat_wat_mix
        call qsat_wat_mix(qsl_TL,TL,p_theta_levels(1,1,k),              &
     &                row_length*rows,l_mixing_ratio)
!
! Set number of points to perform the compressed calculation over
! to zero
!
        NPTI=0.0
! Rows_do2:
        DO j=1,rows
! Row_length_do2:
          DO i=1,row_length
!
! ----------------------------------------------------------------------
! 3. Calculate equivalent critical relative humidity values (RH0) for
!    initiation.
! ----------------------------------------------------------------------
!
! Set up index pointers to critical relative humidity value
!
            IRHI = (MULTRHC * (i - 1)) + 1
            IRHJ = (MULTRHC * (j - 1)) + 1
!
! Calculate relative total humidity with respect to the liquid
! temperature
!
            RHT=(Q(i,j,k)+QCL(i,j,k))/QSL_TL(i,j)
!
! Initialize the equivalent critical relative humidity to a dummy
! negative value
!
            RH0 = -1.0
!
! Do we need to force some initiation from zero liquid cloud amount?
! If so, set RH0 to equal the critical relative humidity
!
             IF ( (CFL(i,j,k)  ==  0.0                                  &
     &          .OR.(CFL(i,j,k)  <   0.05 .AND. T(i,j,k)  <   ZERODEGC))&
     &          .AND. (.not. cumulus(i,j))                              &
     &          .AND.  RHT  >   RHTS(i,j,k)                             &
     &          .AND. RHT  >   (RHCRIT(IRHI,IRHJ,K)+RHCRIT_TOL)  ) THEN
              RH0 = RHCRIT(IRHI,IRHJ,K)
            END IF
!
! Do we need to force some initiation from total liquid cloud cover?
! If so, set RH0 to equal the critical relative humidity
!
            IF (  CFL(i,j,k)  ==  1.0                                   &
     &      .AND.  RHT  <   RHTS(i,j,k)                                 &
     &          .AND. RHT  <   (2.0-RHCRIT(IRHI,IRHJ,K)                 &
     &          - RHCRIT_TOL)  )  THEN
              RH0 = RHCRIT(IRHI,IRHJ,K)
            END IF
!
! ----------------------------------------------------------------------
! 4. Calculate the cloud fraction to initiate to
! ----------------------------------------------------------------------
!
            IF (RH0  >   0.0) THEN
              NPTI=NPTI+1
!
! Is the liquid cloud cover coming down from one or up from zero?
! If it is coming down then reverse the value of RHT and work with
! saturation deficit taking the place of liquid water content.
!
              IF (CFL(i,j,k)  >   0.5) THEN
                RHT=2.0-RHT
              END IF
!
! Calculate the new liquid cloud fraction. Start by calculating QN then
! use the Smith scheme relationships to convert this to a cloud fraction
!
              QN=(RHT-1.0)/(1.0-RH0)
              IF (QN  <=  -1.0) THEN
                C_1 = 0.0
              ELSE IF (QN  <   0.0) THEN
                C_1 = 0.5*(1.0+QN)**(PDF_POWER+1.0)
              ELSE IF (QN  <   1.0) THEN
                C_1 = 1.0 - 0.5*(1.0-QN)**(PDF_POWER+1.0)
              ELSE
                C_1 = 1.0
              ENDIF
!
! Reverse the new cloud fraction back to its correct value if the cloud
! fraction is being decreased from a high value.
!
              IF (CFL(i,j,k)  >   0.5) THEN
                C_1=1.0-C_1
              END IF
!
! Calculate change in total cloud fraction. This depends upon the sign
! of the change of liquid cloud fraction. Change in ice cloud fraction
! is zero.
!
              DELTACL_C(NPTI) = C_1 - CFL(i,j,k)
              DELTACF_C(NPTI) = 0.0
!
! ----------------------------------------------------------------------
! 5. Calculate the liquid water content to initiate to.
! ----------------------------------------------------------------------
!
! We need to iterate to obtain the liquid content because the
! rate of change of gradient of the saturation specific humidity is
! specified as a function of the temperature, not the liquid
! temperature. We will only iterate over the necessary points where
! the initiation works out that a change in cloud fraction is required.
!
! Gather variables
!
              NI(NPTI)=i
              NJ(NPTI)=j
              T_C(NPTI)=T(i,j,k)
              QN_C(NPTI)=QN
              RH0_C(NPTI)=RH0
              QSL_TL_C(NPTI)=QSL_TL(i,j)
              CF_C(NPTI) =CF(i,j,k)
              CFL_C(NPTI)=CFL(i,j,k)
              CFF_C(NPTI)=CFF(i,j,k)
              QCL_C(NPTI)=QCL(i,j,k)
              Q_C(NPTI)=Q(i,j,k)
              p_theta_levels_c(NPTI)=p_theta_levels(i,j,k)
!
! End if for RH0 gt 0
            END IF
! Row_length_do1:
          END DO
! Rows_do1:
        END DO
!
! Calculate change in total cloud fraction.
!
        IF (NPTI  >   0) THEN
! DEPENDS ON: pc2_total_cf
          CALL PC2_TOTAL_CF(                                            &
     &          NPTI,CFL_C,CFF_C,DELTACL_C,DELTACF_C,CF_C)
        END IF
!
! Iterations_do1:
        DO L=1,INIT_ITERATIONS
!
! Calculate saturated specific humidity with respect to the temperature
! (not the liquid temperature).
!
! DEPENDS ON: qsat_wat_mix
          call qsat_wat_mix(qsl_T_c,T_c,p_theta_levels_c,npti,          &
     &                  l_mixing_ratio)
!
! Now loop over the required points
!
! Points_do1
          DO i=1,NPTI
            ALPHA=EPSILON*LC*QSL_T_C(i)/(R*T_C(i)**2)
            AL=1.0/(1.0+LCRCP*ALPHA)
            BS=AL*(1.0-RH0_C(i))*QSL_TL_C(i)
!
            IF (QN_C(i)  <=  -1.0) THEN
              L_BS = 0.0
            ELSE IF (QN_C(i)  <   0.0) THEN
              L_BS = 0.5/(PDF_POWER+2.0)*(1.0+QN_C(i))**(PDF_POWER+2.0)
            ELSE IF (QN_C(i)  <   1.0) THEN
              L_BS = (QN_C(i)+0.5/(PDF_POWER+2.0)                       &
     &               *(1.0-QN_C(i))**(PDF_POWER+2.0))
            ELSE
              L_BS = QN_C(i)
            ENDIF
!
            L_OUT = L_BS * BS
!
! Are we working with saturation deficit instead of liquid water?
! If so, convert back to be the correct way around.
!
            IF (CFL_C(i)  >   0.5) THEN
              Q_OUT = QSL_T_C(i) - L_OUT / AL
              L_OUT = QCL_C(i) + Q_C(i) - Q_OUT
            END IF

!
! Calculate amount of condensation. To ensure convergence move
! slowly towards the calculated liquid water content L_OUT
!
            DESCENT_FACTOR=AL
            DELTAL=DESCENT_FACTOR*(L_OUT-QCL_C(i))
            Q_C(i)   = Q_C(i)   - DELTAL
            QCL_C(i) = QCL_C(i) + DELTAL
            T_C(i)   = T_C(i)   + DELTAL*LCRCP
!
! Points_do1
          END DO
! Iterations_do1:
        END DO
!
! Now scatter back values which have been changed
!
! Points_do2:
        DO i=1,NPTI
!
          Q(NI(i),NJ(i),k)   = Q_C(i)
          QCL(NI(i),NJ(i),k) = QCL_C(i)
          T(NI(i),NJ(i),k)   = T_C(i)
          CF(NI(i),NJ(i),k)  = CF_C(i)
          CFL(NI(i),NJ(i),k)  = CFL_C(i)  + DELTACL_C(i)
!
! Points_do2:
        END DO
! Levels_do1:
      END DO
!
! End of the subroutine
!
      RETURN
      END SUBROUTINE PC2_INITIATE
