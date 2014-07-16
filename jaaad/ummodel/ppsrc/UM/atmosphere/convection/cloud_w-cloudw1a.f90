
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE CLOUD_W------------------------------------------------
!LL
!LL  PURPOSE : CLOUD MICROPHYSICS ROUTINE
!LL
!LL            CALCULATES PRECIPITATION PRODUCED IN LIFTING PARCEL
!LL            FROM LAYER K TO K+1
!LL
!LL            CALL CON_RAD TO CALCULATE PARAMETERS FOR RADIATION
!LL            CALCULATION
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
!LL
!LLEND-----------------------------------------------------------------
!
!
!+ Convection Cloud Microphysics Scheme.
! Subroutine Interface:
      SUBROUTINE CLOUD_W (K, NPNTS, XPKP1, QCLPKP1, QCFPKP1, PREKP1, XSQKP1   &
               , BLOWST, FLXKP1, XPK, QCLPK, QCFPK, THEKP1,QEKP1, BWKP1       &
               , BLAND, QSEKP1, BGMKP1, BTERM, CCA, ICCB, ICCT, TCW, DEPTH    &
               , EKP14, EKP34, DELEXKP1, CCLWP, DELPKP1, CCW, LCCA, LCBASE    &
               , LCTOP, LCCLWP, L_SHALLOW, L_Q_INTERACT, START_LEV)

      Use cv_cntl_mod, Only:                                                  &
          lcv_ccrad

      Use cv_run_mod, Only:                                                   &
          l_ccw, l_no_dcrit, l_fix_udfactor,                                  &
          ud_factor, mparwtr, ccw_for_precip_opt

      IMPLICIT NONE
!
!----------------------------------------------------------------------
! MODEL CONSTANTS
!----------------------------------------------------------------------
!
! CRITDEP start

      ! critical depth of cloud for the formation of
      ! convective precipitation over sea (m)
      REAL,PARAMETER:: CRITDSEA = 1.5E3

      ! critical depth of cloud for the formation of convective
      ! precipitation over land (m)
      REAL,PARAMETER:: CRITDLND = 4.0E3

      ! critical depth of a glaciated cloud for the formation of
      ! convective precipitation (m)
      REAL,PARAMETER:: CRITDICE = 1.0E3

! CRITDEP end
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
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
!----------------------------------------------------------------------
! VECTOR LENGTH AND LOOP COUNTERS
!----------------------------------------------------------------------
!
      INTEGER NPNTS          ! IN VECTOR LENGTH
!
      INTEGER K              ! IN PRESENT MODEL LAYER
!
      INTEGER I              ! LOOP COUNTER
!
! ----------------------------------------------------------------------
! Arguments with Intent IN.  ie:  Input variables.
! ----------------------------------------------------------------------
!
      REAL THEKP1(NPNTS)     ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (K)
!
      REAL QEKP1(NPNTS)      ! IN MIXING RATIO OF CLOUD ENVIRONMENT
                             !    IN LAYER K+1 (KG/KG)
!
      REAL QSEKP1(NPNTS)     ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
!
      REAL XPK(NPNTS)        ! IN PARCEL CLOUD WATER IN LAYER K (KG/KG)
!
      REAL QCLPK(NPNTS)      ! IN PARCEL LIQUID CONDENSATE MIXING RATIO
!                                 IN LAYER K (KG/KG)
!
      REAL QCFPK(NPNTS)      ! IN PARCEL FROZEN CONDENSATE MIXING RATIO
!                                 IN LAYER K (KG/KG)
!
      LOGICAL BWKP1(NPNTS)   ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
!
      LOGICAL BGMKP1(NPNTS)  ! IN MASK FOR PARCELS WHICH ARE
                             !    SATURATED IN LAYER K+1
!
      LOGICAL BLAND(NPNTS)   ! IN LAND/SEA MASK
!
      LOGICAL BTERM(NPNTS)   ! IN MASK FOR PARCELS WHICH TERMINATE IN
                             !    LAYER K+1
!
      LOGICAL BLOWST(NPNTS)  ! IN MASK FOR THOSE POINTS AT WHICH
                             !    STABILITY IS LOW ENOUGH FOR
                             !    CONVECTION TO OCCUR
!
      LOGICAL L_SHALLOW(NPNTS) ! IN MASK FOR POINTS WHERE SHALLOW
                               !    CONVECTION IS LIKELY

      LOGICAL L_Q_INTERACT   ! IN Switch allows overwriting of parcel
!                                 variables (will alter results).
!
      REAL FLXKP1(NPNTS)     ! IN PARCEL MASSFLUX IN LAYER K+1 (PA/S)
!
      REAL XSQKP1(NPNTS)     ! IN EXCESS PARCEL MIXING RATIO IN
                             !    LAYER K+1 (KG/KG)
!
      REAL EKP14(NPNTS)      ! IN ENTRAINMENT RATE AT LEVEL K+1/4
                             !    MULTIPLIED BY APPROPRIATE LAYER
                             !    THICKNESS
!
      REAL EKP34(NPNTS)      ! IN ENTRAINEMNT RATE AT LEVEL K+3/4
                             !    MULTIPLIED BY APPROPRIATE LAYER
                             !    THICKNESS
!
      REAL DELEXKP1(NPNTS)   ! IN DIFFERENCE IN EXNER RATIO ACROSS
                             !    LAYER K+1 (PA)
!
      REAL DELPKP1(NPNTS)    ! IN PRESSURE DIFFERENCE ACROSS LAYER K+1
!
      INTEGER START_LEV(NPNTS)  ! IN LEVEL AT WHICH CONVECTION INITIATED
!
! ----------------------------------------------------------------------
! Arguments with Intent IN/OUT. ie: input variables changed on output.
! ----------------------------------------------------------------------
!
      REAL TCW(NPNTS)        ! INOUT
                             ! IN  TOTAL CONDENSED WATER SUMMED UPTO
                             !     LAYER K (KG/M**2/S)
                             ! OUT TOTAL CONDENSED WATER SUMMED UPTO
                             !     LAYER K+1 (KG/M**2/S)
!
      REAL DEPTH(NPNTS)      ! INOUT
                             ! IN  DEPTH OF CONVECTIVE CLOUD TO
                             !     LAYER K (M)
                             ! OUT DEPTH OF CONVECTIVE CLOUD TO
                             !     LAYER K+1 (M)
!
      REAL CCLWP(NPNTS)      ! INOUT
!                              IN  CONDENSED WATER PATH SUMMED UP TO
!                                  LAYER K (KG/M**2)
!                              OUT CONDENSED WATER PATH SUMMED UP TO
!                                  LAYER K+1 (KG/M**2)
!
      REAL QCLPKP1(NPNTS)    ! INOUT
!                              IN  PARCEL LIQUID CONDENSATE MIXING RATIO
!                                  IN LAYER K+1 BEFORE EXCESS MIXING
!                                  RATIO WATER IS ADDED (KG/KG)
!                              OUT PARCEL LIQUID CONDENSATE MIXING RATIO
!                                  IN LAYER K+1 (KG/KG)
!
      REAL QCFPKP1(NPNTS)    ! INOUT
!                              IN  PARCEL FROZEN CONDENSATE MIXING RATIO
!                                  IN LAYER K+1 BEFORE EXCESS MIXING
!                                  RATIO WATER IS ADDED (KG/KG)
!                              OUT PARCEL FROZEN CONDENSATE MIXING RATIO
!                                  IN LAYER K+1 (KG/KG)
!
! ----------------------------------------------------------------------
! Arguments with Intent OUT. ie: Output variables.
! ----------------------------------------------------------------------
!
      REAL PREKP1(NPNTS)     ! OUT PRECIPITATION FROM PARCEL AS IT RISES
!                                  FROM LAYER K TO K+1 (KG/M**2/S)
!
      REAL XPKP1(NPNTS)      ! OUT PARCEL CLOUD WATER IN LAYER K+1
                             !     (KG/KG)
!
      REAL CCA(NPNTS)        ! OUT CONVECTIVE CLOUD AMOUNT (%)
!
      INTEGER ICCB(NPNTS)    ! OUT CONVECTIVE CLOUD BASE LEVEL
!
      INTEGER ICCT(NPNTS)    ! OUT CONVECTIVE CLOUD TOP LEVEL
!
      REAL CCW(NPNTS)        ! OUT CONVECTIVE CLOUD LIQUID WATER
                             ! (G/KG) ON MODEL LEVELS
!
      REAL LCCA(NPNTS)       ! OUT LOWEST CONV.CLOUD AMOUNT (%)
!
      INTEGER LCBASE(NPNTS)  ! OUT LOWEST CONV.CLOUD BASE LEVEL
!
      INTEGER LCTOP(NPNTS)   ! OUT LOWEST CONV.CLOUD TOP LEVEL
!
      REAL LCCLWP(NPNTS)     ! OUT LOWEST CONV.CLOUD LIQ.WATER PATH
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE LOCALLY DEFINED
!----------------------------------------------------------------------
!
      REAL DCRIT             ! CRITICAL DEPTH AT WHICH PRECIPITATION
                             ! MAY FORM (M)

      Real ::              &
      mparmult               ! Factor used to multiply mparwtr, value
                             ! between 1. and ~3. being 1.0 for deeper
                             ! clouds.   
!
      REAL XMIN              ! AMOUNT OF CLOUD WATER RETAINED BY THE
                             ! PARCEL ON PRECIPITATION (KG/KG)
!
      REAL EPSS              ! (1.0+EKP14)*(1.0+EKP34)
!
      REAL CCW_UD(NPNTS)     ! Cloud water for radiation
!                            !
      LOGICAL L_FALLOUT      ! .true.  is DEFAULT. User DO NOT TOUCH!
!                              .false. KILLs precipitation
!                                      DIAGNOSTIC TEST CASE ONLY:
!                                     (prevent convection updating qT).
      PARAMETER(L_FALLOUT=.TRUE.)

!----------------------------------------------------------------------
!
!- End of Header
!
! ==Main Block==--------------------------------------------------------
!
! ----------------------------------------------------------------------
! CALCULATE PARCEL CLOUD CONDENSATE BEFORE PRECIPITATION
!
! UM DOCUMENTATION PAPER 27
! SECTION (2B), EQUATION (13A or 13C)
! ----------------------------------------------------------------------
! Cloudw_do1:
      DO I=1, NPNTS
        EPSS = (1.+EKP14(I)) * (1.+EKP34(I))
        XPKP1(I) = ( XPK(I)/EPSS ) + XSQKP1(I)
!
        IF (L_Q_INTERACT) THEN
          IF (BWKP1(I)) THEN
            QCLPKP1(I) = QCLPKP1(I) + XSQKP1(I)
          ELSE
            QCFPKP1(I) = QCFPKP1(I) + XSQKP1(I)
          ENDIF
        ENDIF
      END DO ! Cloudw_do1
!
      IF (L_Q_INTERACT) THEN
! Cloudw_do2:
        DO I=1, NPNTS
!         Xpk(p1) are used for inputs only: need to work out how to use
!          Qclpk and Qcfpk instead. Meanwhile, overwrite Xpk(p1).
          XPK  (I) = QCLPK  (I) + QCFPK  (I)
          XPKP1(I) = QCLPKP1(I) + QCFPKP1(I)
        END DO ! Cloudw_do2
!
      ENDIF
!
!L----------------------------------------------------------------------
!L STORE CONVECTIVE CLOUD LIQUID WATER BEFORE PRECIPITATION
!L----------------------------------------------------------------------
!L
      DO I=1,NPNTS
        CCW(I) = XPKP1(I)
      END DO
      IF (.NOT. L_CCW) THEN
!L
!L----------------------------------------------------------------------
!L CALCULATE CONVECTIVE CLOUD BASE, CONVECTIVE CLOUD TOP , TOTAL
!L CONDENSED WATER/ICE AND CONVECTIVE CLOUD AMOUNT
!L
!L SUBROUTINE CON_RAD
!L
!L UM DOCUMENTATION PAPER 27
!L SECTION (9)
!L----------------------------------------------------------------------

! DEPENDS ON: con_rad
        CALL CON_RAD(K,XPK,XPKP1,FLXKP1,BTERM,CCA,ICCB,ICCT,START_LEV,  &
             TCW,CCW,CCLWP,DELPKP1,LCCA,LCBASE,LCTOP,LCCLWP,NPNTS,      &
             L_Q_INTERACT)
      ENDIF

!L----------------------------------------------------------------------
!L CALCULATE CLOUD DEPTH AND ASSIGN CRITICAL CLOUD DEPTHS
!L
!L UM DOCUMENTATION PAPER 27
!L SECTION (8), EQUATION (34), (35)
!L----------------------------------------------------------------------

      Do I=1,NPNTS
        IF ( BLOWST(I) ) DEPTH(I) = 0.

! This could be improved and simplified by using actual model heights

        IF ( BGMKP1(I) )                                                &
          DEPTH(I) = DEPTH(I) + ( CP * THEKP1(I) *                      &
                                   (1.0+C_VIRTUAL*QEKP1(I)) *           &
                                               DELEXKP1(I)/G )
      End Do  

!L----------------------------------------------------------------------
!L CALCULATE PRECIPITATION FROM LAYER K+1 AND ADJUST CLOUD WATER
!L
!L UM DOCUMENTATION PAPER 27
!L SECTION (8), EQUATION (36)  (ccw_for_precip_opt=0)
!L----------------------------------------------------------------------

      Select Case (ccw_for_precip_opt)

      Case (3)            ! Manoj's function for congestus
                          ! 
                          
        Do i=1,npnts

          mparmult = 1.5 + 0.5*tanh((3000.0 -depth(i))/800.0)
          xmin = MIN(mparwtr*mparmult, 0.5*qsekp1(i)) 

          ! If a land point and liquid water in the layer 
          ! increase the minimum cloud water for precipitation
          ! The reasons for this are an attempt to take some account of 
          ! more aerosols over land leading to more small cloud drops and
          ! therefore more cloud water before precipitation.  

          If (bwkp1(i) .and. bland(i)) xmin =xmin*2.0

          ! limit max value to 0.003 kg/kg
          xmin = MIN(xmin,0.003)          

          If (l_q_interact) Then   ! PC2
          ! Limit xmin
             xmin = max(xmin, 2.0E-4)
          End If

          ! Precipitate if cloud water in the layer > xmin

          If ( xpkp1(i) > xmin ) Then

            prekp1(i) = (xpkp1(i) -xmin) * flxkp1(i) /g 

            If (l_q_interact) Then        ! PC2
            ! Dodgy flux calculation : hydrostatic approximation?
            ! This calculation needs sorting properly. TEMPORARY getaround.
               qclpkp1(i) = qclpkp1(i) * xmin / xpkp1(i)
               qcfpkp1(i) = qcfpkp1(i) * xmin / xpkp1(i)
            End If 

            xpkp1(i) = xmin      ! cloud water in layer after precip
            ccw_ud(i)= xpkp1(i)*ud_factor 

          Else     ! no precipitation
            prekp1(i) = 0.0

          ! UD_FACTOR ought to be applied here too but, in order to 
          ! simplify the code, this is corrected by moving the point
          ! at which UD_FACTOR is applied to after convection 
          ! in CONV_CTL and passing UD_FACTOR of 1 to here
            ccw_ud(i) = xpkp1(i)  

          End If     ! test on whether to precipitate

        End Do

      Case (2)            ! Manoj's function dependent on depth of cloud
                          ! Also removed extra tests on l_ccw (usually true)
                          ! and l_fallout

        Do i=1,npnts

          mparmult = 2.0 + 1.0*tanh((1500.0 -depth(i))/1000.0)
          xmin = MIN(mparwtr*mparmult, 0.5*qsekp1(i)) 

          ! If a land point and liquid water in the layer 
          ! increase the minimum cloud water for precipitation
          ! The reasons for this are an attempt to take some account of 
          ! more aerosols over land leading to more small cloud drops and
          ! therefore more cloud water before precipitation.  

          If (bwkp1(i) .and. bland(i)) xmin =xmin*2.0

          If (l_q_interact) Then   ! PC2
          ! Limit xmin
             xmin = max(xmin, 2.0E-4)
          End If

          ! Precipitate if cloud water in the layer > xmin

          If ( xpkp1(i) > xmin ) Then

            prekp1(i) = (xpkp1(i) -xmin) * flxkp1(i) /g 

            If (l_q_interact) Then        ! PC2
            ! Dodgy flux calculation : hydrostatic approximation?
            ! This calculation needs sorting properly. TEMPORARY getaround.
               qclpkp1(i) = qclpkp1(i) * xmin / xpkp1(i)
               qcfpkp1(i) = qcfpkp1(i) * xmin / xpkp1(i)
            End If 

            xpkp1(i) = xmin      ! cloud water in layer after precip
            ccw_ud(i)= xpkp1(i)*ud_factor 

          Else     ! no precipitation
            prekp1(i) = 0.0

          ! UD_FACTOR ought to be applied here too but, in order to 
          ! simplify the code, this is corrected by moving the point
          ! at which UD_FACTOR is applied to after convection 
          ! in CONV_CTL and passing UD_FACTOR of 1 to here
            ccw_ud(i) = xpkp1(i)  

          End If     ! test on whether to precipitate

        End Do

      Case (1)            ! As L_no_dcrit = .true. option
                          ! No test on a critical depth for precipitation
                          ! Option in use for HadGEM2

        Do i=1,npnts
          xmin = MIN(mparwtr, 0.5*qsekp1(i)) 

          ! If a land point and liquid water in the layer 
          If (bwkp1(i) .and. bland(i)) xmin =xmin*2.0

          If (l_q_interact) Then
          ! Limit xmin
             xmin = max(xmin, 2.0E-4)
          End If

          ! Precipitate if 
          ! either  shallow  or (not shallow and l_ccw)
          ! and cloud water in the layer > xmin

          If ( ( (l_shallow(i)).or.((.not.l_shallow(i)) .and. l_ccw ) )  &
               .and. (xpkp1(i) > xmin) .and. l_fallout) Then

            prekp1(i) = (xpkp1(i) -xmin) * flxkp1(i) /g 

            If (l_q_interact) Then        ! PC2
            ! Dodgy flux calculation : hydrostatic approximation?
            ! This calculation needs sorting properly. TEMPORARY getaround.
               qclpkp1(i) = qclpkp1(i) * xmin / xpkp1(i)
               qcfpkp1(i) = qcfpkp1(i) * xmin / xpkp1(i)
            End If 

            xpkp1(i) = xmin
            ccw_ud(i)= xpkp1(i)*ud_factor 

          Else     ! no precipitation
            prekp1(i) = 0.0

           ! UD_FACTOR ought to be applied here too but, in order to 
           ! simplify the code, this is corrected by moving the point
           ! at which UD_FACTOR is applied to after convection 
           ! in CONV_CTL and passing UD_FACTOR of 1 to here
            ccw_ud(i) = xpkp1(i)  

          End If     ! test on whether to precipitate

        End Do

      Case Default       ! 0 - original convection code using a critical depth

        Do i=1,npnts

          ! Assign critical cloud depths
  
          If (.not.bwkp1(i)) Then    ! Ice water cloud in this layer 
            dcrit = critdice
          Else If (bland(i)) Then    ! liquid water over land
            dcrit = critdlnd
          Else                       ! liquid water over sea 
            dcrit = critdsea
          End If
          
          xmin = MIN(mparwtr, 0.5*qsekp1(i)) 

          If (l_q_interact) Then
          ! Limit xmin
             xmin = max(xmin, 2.0E-4)
          End If

          ! precipitate if 
          ! Either  depth of cloud > dcrit OR convection is not shallow
          ! and    cloud water in the layer > xmin 

          If ( ( (depth(i) > dcrit).or.((.not.l_shallow(i)) .and. l_ccw ) )  &
               .and. (xpkp1(i) > xmin) .and. l_fallout) Then

            prekp1(i) = (xpkp1(i) -xmin) * flxkp1(i) /g 

            If (l_q_interact) Then        ! PC2
            ! Dodgy flux calculation : hydrostatic approximation?
            ! This calculation needs sorting properly. TEMPORARY getaround.
               qclpkp1(i) = qclpkp1(i) * xmin / xpkp1(i)
               qcfpkp1(i) = qcfpkp1(i) * xmin / xpkp1(i)
            End If 
            xpkp1(i) = xmin
            ccw_ud(i)= xpkp1(i)*ud_factor 

          Else     ! no precipitation
            prekp1(i) = 0.0

           ! UD_FACTOR ought to be applied here too but, in order to 
           ! simplify the code, this is corrected by moving the point
           ! at which UD_FACTOR is applied to after convection 
           ! in CONV_CTL and passing UD_FACTOR of 1 to here
            ccw_ud(i) = xpkp1(i)  

          End If     ! test on whether to precipitate

        End Do

      End Select       ! test on value of ccw_for_precip_opt 

      If (L_CCW) Then

!L----------------------------------------------------------------------
!L CALCULATE CONVECTIVE CLOUD BASE, CONVECTIVE CLOUD TOP , TOTAL
!L CONDENSED WATER/ICE AND CONVECTIVE CLOUD AMOUNT
!L
!L SUBROUTINE CON_RAD - MOVED TO AFTER RAIN OUT HAS OCCURRED IF L_CCW
!L IS TRUE (SET IN UMUI).
!L UM DOCUMENTATION PAPER 27
!L SECTION (9)
!L----------------------------------------------------------------------

! DEPENDS ON: con_rad
        CALL CON_RAD(K,XPK,CCW_UD,FLXKP1,BTERM,CCA,ICCB,ICCT,START_LEV, &
             TCW,CCW,CCLWP,DELPKP1,LCCA,LCBASE,LCTOP,LCCLWP,NPNTS,      &
             L_Q_INTERACT)

!-----------------------------------------------------------------------
! STORE CONVECTIVE CLOUD LIQUID WATER AFTER PRECIPITATION
!-----------------------------------------------------------------------

        Do i=1,npnts
          ccw(i) = ccw_ud(i)
        End Do
      End If

      RETURN
      END SUBROUTINE CLOUD_W
