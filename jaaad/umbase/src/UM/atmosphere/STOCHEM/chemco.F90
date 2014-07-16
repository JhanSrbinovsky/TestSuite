#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE CHEMCO(rc,tc,m,h2o,asize)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : CALCULATES RATE COEFFICIENTS
!-
!-   Inputs  : TC,M,H2O
!-   Outputs : RC
!-   Controls:
!-
!
! Current Owner of Code: M.G. Sanderson
!
! History:
! Version   Date                    Comment
!  3.4    09/12/93  Created.  W.J. Collins
!  4.4    26/11/96  JPL 1994 rate coefficients, Atkinson 1994 and
!                   MCM 1996. R.G. Derwent
!  4.4    07/01/97  Added aqueous phase reactions RC260-262.
!                   W.J. Collins.
!  5.5    21/10/03  Vectorised code. M.G. Sanderson
!  6.1    20/07/04  Minor rate constant corrections. M.G. Sanderson
!  6.2    28/03/06  Minor changes for vn6.2  M.G. Sanderson
!-
!VVV  V3.0  CHEMCO 20/10/04 - 76 species chemistry
!----------------------------------------------------------------------
      USE IN_STOCHEM_CHM
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: asize                 ! Array size
      REAL, DIMENSION(chunk), INTENT(IN) :: tc     ! Temperature (K)
      REAL, DIMENSION(chunk), INTENT(IN) :: m      ! Air density (molec
      REAL, DIMENSION(chunk), INTENT(IN) :: h2o    ! Water concn (vmr)
      REAL, DIMENSION(chunk,nr), INTENT(OUT) :: rc ! Rate Constants
!
      INTEGER, PARAMETER :: npdep = 10             ! No. press-dep react
      INTEGER :: i, j, k                           ! Loop counts
      INTEGER :: r                                 ! Reaction type
      INTEGER :: nklo, nkhi, nres ! Reaction nos of klo, khi and locatio
!                              of result for pressure-dependent reaction
!
      REAL :: one_over_e = 0.3678794  ! 1/exp(1)
      REAL :: tfc = 1.0 / 300.0  ! Used in k's above
      REAL :: taq = 1.0 / 298.0  ! For aqueous chemistry
      REAL :: ar                 ! Arrhenius pre-exponential factor
      REAL :: n                  ! Temperature power (K)
      REAL :: ea                 ! Activation energy (K)
      REAL :: fac1               ! Temporary store
      REAL :: fac2               ! Temporary store
      REAL :: fac3               ! Temporary store
      REAL :: brn                ! Temporary store
      REAL :: f                  ! Appropriate value of Fc factor
      REAL :: fc                 ! Fc factor for R29
!
! Set up rate coefficient data type
! Rate constant k = A.(T/300)**n.exp(E/T)
! E is -Ea/R, where Ea is the activation energy and R is the ideal gas c
! Hence E has uniits of K.
! r is the reaction type:
!  0 - no reaction
!  1 - k = A  (temperature-independent)
!  2 - k = A.(T/300)**n
!  3 - k = A.exp(E/T)
!  4 - k = A.(T/300)**n.exp(E/T)
!  5 - k = aqueous phase
!
      TYPE RR
        REAL :: a    ! Pre-exponential factor
        REAL :: n    ! Power for temperature
        REAL :: e    ! Activation Energy (actually -Ea/R)
        INTEGER :: r ! Reaction type (see above)
      ENDTYPE RR
!
! Set up type to hold information needed to calculate pressure-dependent
! rate constants
      TYPE PD
        INTEGER :: klo  ! Reaction number holding low-pressure rate cons
        INTEGER :: khi  ! Reaction number holding high-pressure rate con
        INTEGER :: res  ! Rate constant to store final result in
        REAL    :: fc   ! Fc factor
      ENDTYPE PD
!
      TYPE(RR), DIMENSION(nr), SAVE :: rk  ! Holds A, n, Ea for each rat
      TYPE(RR), SAVE :: ZERO_RATE          ! Null data
      TYPE(PD), DIMENSION(npdep), SAVE :: pdep ! Pressure-dependent reac
!
      LOGICAL, SAVE :: first = .TRUE.
!
      IF (first) THEN   ! Set up arrays zero, rk and pdep on first call
      zero_rate = RR(0.0, 0.0, 0.0, 0)
!
! Set up rk array 19 lines at a time. 19 continuation lines are Fortran-
! standard limit (according to warning messages when more are specified)
!
!              A       n     -Ea/R  rtype   Num   Reaction Details
      rk(1:19) = (/                                                     &
     &  RR( 6.0e-34, -2.3,     0.0,  2),                                &
                                          !   1 : O + O2 + M = O3 + M
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  RR( 9.0e-32, -1.5,     0.0,  2),                                &
                                          !   4 : O + NO + M = NO2 + M
     &  RR( 3.0e-11,  0.0,     0.0,  1),                                &
                                          !   5 : O + NO + M = NO2 + M
     &  RR( 3.2e-11,  0.0,    70.0,  3),                                &
                                          !   6 : O(1D) + M = O(3P) + M
     &  RR( 1.8e-11,  0.0,   110.0,  3),                                &
                                          !   7 : O(1D) + M 2nd expressi
     &  RR( 2.2e-10,  0.0,     0.0,  1),                                &
                                          !   8 : O(1D) + H2O = 2 OH
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  RR( 2.0e-12,  0.0, -1400.0,  3),                                &
                                          !  11 : NO + O3 = NO2 + O2
     &  RR( 1.2e-13,  0.0, -2450.0,  3),                                &
                                          !  12 : NO2 + O3 = NO3 + O2
     &  RR( 1.6e-12,  0.0,  -940.0,  3),                                &
                                          !  13 : OH + O3 = HO2 + O2
     &  RR( 1.1e-14,  0.0,  -500.0,  3),                                &
                                          !  14 : HO2 + O3 = OH + O2 + O
     &  RR( 1.5e-11,  0.0,   170.0,  3),                                &
                                          !  15 : NO + NO3 = NO2 + NO2
     &  RR( 6.5e-12,  0.0,   120.0,  3),                                &
                                          !  16 : NO2 + O = NO + O2 - ad
     &  RR( 3.7e-12,  0.0,   250.0,  3),                                &
                                          !  17 : NO + HO2 = OH + NO2
     &  zero_rate,                                                      &
     &  RR( 4.5e-14,  0.0, -1260.0,  3)/) !  19 : NO2 + NO3 = NO + NO2 +
      rk(20:38) = (/                                                    &
     &  RR( 2.2e-30, -3.9,     0.0,  2),                                &
                                          !  20 : NO2 + NO3 + M = N2O5 +
     &  RR( 2.6e-30, -3.2,     0.0,  2),                                &
                                          !  21 : NO2 + OH + M = HNO3 +
     &  RR( 1.8e-31, -3.2,     0.0,  2),                                &
                                          !  22 : NO2 + HO2 + M = HO2NO2
     &  RR( 2.1e-27,  0.0, 10900.0,  3),                                &
                                          !  23 : HO2NO2 + M = HO2 + NO2
     &  RR( 1.3e-12,  0.0,   380.0,  3),                                &
                                          !  24 : OH + HO2NO2 = H2O + NO
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  RR( 8.5e-13,  0.0, -2450.0,  3),                                &
                                          !  27 : NO3 + NO3 = NO2 + NO2
     &  RR( 3.5e-12,  0.0,  -925.0,  3),                                &
                                          !  28 : OH + NH3  = NH2 + H2O
     &  RR( 2.2e-03, -4.4,-11080.0,  4),                                &
                                          !  29 : N2O5 + M = NO2 + NO3 +
     &  RR( 4.8e-11,  0.0,   250.0,  3),                                &
                                          !  30 : HO2 + OH = H2O + O2
     &  RR( 2.9e-12,  0.0,  -160.0,  3),                                &
                                          !  31 : OH + H2O2 = H2O + HO2
     &  RR( 9.2e-13,  0.0,     0.0,  1),                                &
                                          !  32 : NO3 + HO2 = HNO3 + O2
     &  RR( 5.5e-12,  0.0, -2000.0,  3),                                &
                                          !  33 : OH + H2 (+O2) = HO2 +
     &  RR( 3.6e-12,  0.0,     0.0,  1),                                &
                                          !  34 : NO3 + HO2 = OH + NO2 +
     &  RR( 7.2e-15,  0.0,   785.0,  3),                                &
                                          !  35 : OH + HNO3 = NO3 + H2O
     &  RR( 2.3e-13,  0.0,   600.0,  3),                                &
                                          !  36 : HO2 + HO2 (+M,H2O) = H
     &  RR( 1.9e-33,  0.0,   890.0,  3),                                &
                                          !  37 : HO2 + HO2 (+M,H2O) = H
     &  RR( 1.4e-21,  0.0,  2200.0,  3)/) !  38 : HO2 + HO2 (+M,H2O) = H
      rk(39:57) = (/                                                    &
     &  RR( 3.0e-31, -3.3,     0.0,  2),                                &
                                          !  39 : OH + SO2 + M = HOSO2 +
     &  RR( 4.0e-17,  0.0,     0.0,  1),                                &
                                          !  40 : SO2 + CH3O2 = PRODUCTS
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  RR( 5.0e-06,  0.0,     0.0,  1),                                &
                                          !  43 : UPTAKE OF NOy BY AEROS
     &  RR( 1.9e-12,  0.0,   190.0,  3),                                &
                                          !  44 : OH + CH3OOH = CH3O2 +
     &  RR( 1.0e-12,  0.0,   190.0,  3),                                &
                                          !  45 : OH + CH3OOH = OH + HCH
     &  RR( 1.5e-12, -0.7,     0.0,  2),                                &
                                          !  46 : NO2 + NO3 + M = N2O5 +
     &  RR( 2.4e-11, -1.3,     0.0,  2),                                &
                                          !  47 : NO2 + OH + M = HNO3 +
     &  RR( 4.7e-12, -1.4,     0.0,  2),                                &
                                          !  48 : NO2 + HO2 + M = HO2NO2
     &  RR( 9.7e+14,  0.1,-11080.0,  4),                                &
                                          !  49 : N2O5 + M = NO2 + NO3 +
     &  RR( 4.1e-16,  0.0,  1440.0,  3),                                &
                                          !  50 : OH + HNO3 = NO3 + H2O
     &  RR( 1.9e-33,  0.0,   725.0,  3),                                &
                                          !  51 : OH + HNO3 = NO3 + H2O
     &  RR( 1.5e-12,  0.0,     0.0,  1),                                &
                                          !  52 : OH + SO2 + M = HOSO2 +
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate /)
      rk(58:76) = (/                                                    &
     &  RR( 0.08163,  0.0,     0.0,  1),                                &
                                          !  58 : Stratospheric methane
     &  RR(2.374e-13,2.58, -1082.0,  4),                                &
                                          !  59 : OH + CH4 = CH3O2 + H2O
     &  RR( 4.2e-12,  0.0,   180.0,  3),                                &
                                          !  60 : NO + CH3O2 = HCHO + HO
     &  RR(9.14e-14,  0.0,   416.0,  3),                                &
                                          !  61 : CH3O2 + CH3O2 = 2 HCHO
     &  RR( 25.0,     0.0, -1165.0,  3),                                &
                                          !  62 : branching ratio for R6
!                                                               = CH3OH
     &  RR( 6.7e-12,  0.0,  -600.0,  3),                                &
                                          !  63 : CH3OH + OH = HO2 + HCH
     &  zero_rate,                                                      &
     &  RR( 3.8e-13,  0.0,   800.0,  3),                                &
                                          !  65 : CH3O2 + HO2 = CH3OOH +
     &  RR( 1.0e-11,  0.0,     0.0,  1),                                &
                                          !  66 : OH + HCHO = HO2 + CO +
     &  RR( 5.8e-16,  0.0,     0.0,  1),                                &
                                          !  67 : NO3 + HCHO = HO2 + CO
     &  zero_rate,                                                      &
     &  RR(3.54e-33,  0.0,     0.0,  1),                                &
                                          !  69 : OH + CO Press dependen
     &  RR( 1.5e-13,  0.0,     0.0,  1),                                &
                                          !  70 : OH + CO = HO2 + CO2
     &  RR(1.359e-12, 2.0,  -492.0,  4),                                &
                                          !  71 : OH + C2H6 = C2H5O2 - A
     &  RR( 8.7e-12,  0.0,     0.0,  1),                                &
                                          !  72 : C2H5O2 + NO = CH3CHO +
     &  RR( 2.0e-13,  0.0,     0.0,  1),                                &
                                          !  73 : C2H5O2 + CH3O2 = CH3CH
     &  RR( 4.4e5,    0.0, -3910.0,  3),                                &
                                          !  74 : branching ratio for R8
     &  RR( 6.0e-12,  0.0,   250.0,  3),                                &
                                          !  75 : OH + CH3CHO = CH3COO2
     &  RR( 2.7e-28, -7.1,     0.0,  2)/) !  76 : CH3COO2 + NO2 + M = CH
      rk(77:95) = (/                                                    &
     &  RR( 1.2e-11, -0.9,     0.0,  2),                                &
                                          !  77 : CH3COO2 + NO2 + M = CH
     &  RR( 5.5e-03,  0.0,-12064.0,  3),                                &
                                          !  78 : CH3COO2NO2 + M = CH3CO
     &  RR( 2.0e-11,  0.0,     0.0,  1),                                &
                                          !  79 : CH3COO2 + NO = CH3O2 +
     &  RR( 5.1e-12,  0.0,   272.0,  3),                                &
                                          !  80 : CH3O2 + CH3COO2 = HCHO
     &  RR(1.359e-12, 2.0,   190.0,  4),                                &
                                          !  81 : OH + n-C4H10 (+O2) = s
     &  RR( 3.9e+16,  0.0,-13628.0,  3),                                &
                                          !  82 : CH3COO2NO2 + M = CH3CO
     &  RR( 4.1e-12,  0.0,     0.0,  1),                                &
                                          !  83 : s-C4H9O2 + NO = CH3COC
     &  RR( 1.0e-14,  0.0,     0.0,  1),                                &
                                          !  84 : s-C4H9O2 + CH3O2 = CH3
     &  zero_rate,                                                      &
     &  RR(2.916e-13, 2.0,   414.0,  4),                                &
                                          !  86 : OH + CH3COC2H5 (+O2) =
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  RR( 9.8e-14,  0.0,  -100.0,  3),                                &
                                          !  90 : C2H5O2 + C2H5O2 = C2H5
     &  RR( 2.8e-12,  0.0,   530.0,  3),                                &
                                          !  91 : CH3COO2 + CH3COO2 = 2
     &  RR( 1.35e-12, 2.0,    44.0,  4),                                &
                                          !  92 : OH + C3H8 = i-C3H7O2 +
     &  RR( 4.8e-12,  0.0,     0.0,  1),                                &
                                          !  93 : i-C3H7O2 + NO = NO2 +
     &  RR(4.806e-13, 2.0,  -230.0,  4),                                &
                                          !  94 : CH3COCH3 + OH = CH3COC
     &  RR( 5.3e-12,  0.0,     0.0,  1)/) !  95 : CH3COCH2O2 + NO = NO2
      rk(96:114) = (/                                                   &
     &  RR( 1.2e-12,  0.0,     0.0,  1),                                &
                                          !  96 : CH3COCH2O2 + CH3O2 = H
     &  RR( 4.0e-14,  0.0,     0.0,  1),                                &
                                          !  97 : C3H7O2 + CH3O2 = HCHO
     &  RR( 9.5e-13,  0.0,  -650.0,  3),                                &
                                          !  98 : OH + PAN = CH3COO2 + H
     &  RR( 7.5e-13,  0.0,   700.0,  3),                                &
                                          !  99 : HO2 + C2H5O2 = C2H5OOH
     &  RR( 1.0e-11,  0.0,     0.0,  1),                                &
                                          ! 100 : OH + C2H5OOH = CH3CHO
     &  RR(1.92e-13,  0.0,  1250.0,  3),                                &
                                          ! 101 : HO2 + C3H7O2 = C3H7OOH
     &  RR(2.42e-11,  0.0,     0.0,  1),                                &
                                          ! 102 : OH + C3H7OOH = CH3COCH
     &  RR(2.24e-13,  0.0,  1250.0,  3),                                &
                                          ! 103 : HO2 + SC4H9O2 = SC4H9O
     &  RR( 3.2e-11,  0.0,     0.0,  1),                                &
                                          ! 104 : OH + SC4H9OOH = CH3COC
     &  RR( 5.0e-12,  0.0,     0.0,  1),                                &
                                          ! 105 : NO + CH3COCH(O2)CH3 =
     &  RR( 1.0e-14,  0.0,     0.0,  1),                                &
                                          ! 106 : CH3O2 + CH3COCH(O2)CH3
!                                                      HCHO + HO2 + CH3C
     &  zero_rate,                                                      &
     &  RR( 7.0e-29, -3.1,     0.0,  2),                                &
                                          ! 108 : OH + C2H4 (+O2) = HOC2
     &  RR( 9.0e-12,  0.0,     0.0,  1),                                &
                                          ! 109 : OH + C2H4 (+O2) = HOC2
     &  RR( 9.0e-12,  0.0,     0.0,  1),                                &
                                          ! 110 : HOC2H4O2 + NO = 2 HCHO
     &  RR( 1.0e-12,  0.0,     0.0,  1),                                &
                                          ! 111 : CH3O2 + HOC2H4O2 = 3 H
     &  RR( 1.2e-14,  0.0, -2630.0,  3),                                &
                                          ! 112 : O3 + C2H4 = HCHO + 0.4
!                                                 + 0.22 CO2 + 0.31 H2O
     &  zero_rate,                                                      &
     &  zero_rate /)
      rk(115:133) = (/                                                  &
     &  zero_rate,                                                      &
                                          ! 115
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
                                          ! 120
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  RR( 4.0e-15,  0.0, -1900.0,  3),                                &
                                          ! 123 : O3 + C3H6 = HCHO + 0.3
!                                   + 0.60 CO2 + 0.28 OH + 0.12 CH3OH +
     &  RR( 2.6e-15,  0.0, -1900.0,  3),                                &
                                          ! 124 : O3 + C3H6 = CH3CHO + 0
!                                                             + 0.42 CO2
     &  RR( 8.0e-27, -3.5,     0.0,  2),                                &
                                          ! 125 : OH + C3H6 (+O2) = HOC3
     &  RR( 9.0e-12,  0.0,     0.0,  1),                                &
                                          ! 126 : HOC3H6O2 + NO = HCHO +
     &  RR( 1.0e-14,  0.0,     0.0,  1),                                &
                                          ! 127 : CH3O2 + HOC3H6O2 = 2 H
!
     &  RR(7.86e-15,  0.0, -1913.0,  3),                                &
                                          ! 128 : O3 + C5H8 = MVK + 0.78
!                                                                      +
     &  RR(7.56e-16,  0.0, -1521.0,  3),                                &
                                          ! 129 : O3 + MVK = MGLYOX + 0.
!                                                                      +
     &  RR( 1.0e-14,  0.0,     0.0,  1),                                &
                                          ! 130 : CH2OO + NO = NO2 + HCH
     &  RR( 1.0e-15,  0.0,     0.0,  1),                                &
                                          ! 131 : CH2OO + NO2 = NO3 + HC
     &  RR( 5.8e-17,  0.0,     0.0,  1),                                &
                                          ! 132 : CH2OO + H2O = HCOOH +
     &  zero_rate /)
      rk(134:152) = (/                                                  &
     &  RR( 3.0e-11,  0.0,     0.0,  1),                                &
                                          ! 134 : OH + C3H6 (+O2) = HOC3
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  RR( 37.0,     0.0,  -660.0,  3),                                &
                                          ! 139 : Branching ratio for R1
     &  RR( 4.3e-13,  0.0,  1040.0,  3),                                &
                                          ! 140 : CH3CO3 + HO2 = y O3 +
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
                                          ! 145
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
                                          ! 150
     &  zero_rate,                                                      &
     &  zero_rate /)                      ! 152
      rk(153:171) = (/                                                  &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
                                          ! 155
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
                                          ! 160
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
                                          ! 165
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
                                          ! 170
     &  zero_rate /)                      ! 171
      rk(172:190) = (/                                                  &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
                                          ! 175
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
                                          ! 180
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
                                          ! 185
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate /)                      ! 190
      rk(191:209) = (/                                                  &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
                                          ! 195
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
                                          ! 200
     &  RR( 5.7e-12,  0.0, -4426.0,  3),                                &
                                          ! 201 : NO3 + C2H6 = C2H5O2 +
     &  RR(2.76e-12,  0.0, -3279.0,  3),                                &
                                          ! 202 : NO3 + n-C4H10 = s-C4H9
     &  RR(4.392e-13, 2.0, -2282.0,  4),                                &
                                          ! 203 : NO3 + C2H4 = C2H4NO3 =
     &  RR(4.95e-12,  0.0,     0.0,  1),                                &
                                          ! 204 : CH2(NO3)CHO + OH = HCH
     &  RR(4.59e-13,  0.0, -1156.0,  3),                                &
                                          ! 205 : NO3 + C3H6 = C3H6NO3 =
     &  RR(5.25e-12,  0.0,     0.0,  1),                                &
                                          ! 206 : CH3CH(NO3)CHO + OH = C
     &  zero_rate,                                                      &
     &  RR( 1.4e-12,  0.0, -1900.0,  3),                                &
                                          ! 208 : NO3 + CH3CHO = CH3COO2
     &  RR(3.03e-12,  0.0,  -446.0,  3)/) ! 209 : NO3 + C5H8 = (NO3)C4H6
      rk(210:228) = (/                                                  &
     &  RR(4.16e-11,  0.0,     0.0,  1),                                &
                                          ! 210 : (NO3)C4H6CHO + OH = MV
     &  RR( 1.7e-42,  0.0,  7810.0,  3),                                &
                                          ! 211 : OH + DMS = DMSO + HO2
     &  RR( 1.0e-11,  0.0,     0.0,  1),                                &
                                          ! 212 : HO2 + TOLP1 = MEMALD +
     &  RR( 1.0e-11,  0.0,     0.0,  1),                                &
                                          ! 213 : NO2 + TOLP1 = ORGNIT
     &  RR( 9.6e-12,  0.0,  -234.0,  3),                                &
                                          ! 214 : OH + DMS = CH3SO + HCH
     &  RR( 5.5e-31,  0.0,  7460.0,  3),                                &
                                          ! 215 : OH + DMS = DMSO + HO2
     &  RR( 6.0e-13,  0.0,     0.0,  1),                                &
                                          ! 216 : CH3SO + O3 = CH3SO2
     &  RR( 8.0e-12,  0.0,     0.0,  1),                                &
                                          ! 217 : CH3SO + NO2 = CH3SO2 +
     &  RR( 3.0e-13,  0.0,     0.0,  1),                                &
                                          ! 218 : CH3SO2 + O3 = CH3SO3
     &  RR( 4.0e-12,  0.0,     0.0,  1),                                &
                                          ! 219 : CH3SO2 + NO2 = CH3SO3
     &  RR( 5.0e+13,  0.0, -8656.0,  3),                                &
                                          ! 220 : CH3SO2 + O2 = CH3O2 +
     &  RR( 5.0e-11,  0.0,     0.0,  1),                                &
                                          ! 221 : CH3SO3 + HO2 = MSA
     &  RR( 5.0e+13,  0.0,-11071.0,  3),                                &
                                          ! 222 : CH3SO3 + O2 = CH3O2 +
     &  RR( 1.6e-15,  0.0,     0.0,  1),                                &
                                          ! 223 : CH3SO3 + HCHO = MSA +
     &  RR( 5.8e-11,  0.0,     0.0,  1),                                &
                                          ! 224 : OH + DMSO =  DMSO2 + H
     &  RR( 1.0e-12,  0.0,     0.0,  1),                                &
                                          ! 225 : OH + DMSO2 = CH3SO2CH2
     &  RR( 4.1e-12,  0.0,   180.0,  3),                                &
                                          ! 226 : CH3SO2CH2O2 + NO = NO2
     &  RR( 3.0e-13,  0.0,     0.0,  1),                                &
                                          ! 227 : CH3SO2CH2O2 + CH3O2 =
     &  RR( 1.9e-13,  0.0,   520.0,  3)/) ! 228 : NO3 + DMS = CH3SO + HC
      rk(229:247) = (/                                                  &
     &  zero_rate,                                                      &
     &  RR(1.37e-11,  0.0,     0.0,  1),                                &
                                          ! 230 : OH + (C6H4)(CH3)2 = HO
     &  RR( 1.0e-11,  0.0,     0.0,  1),                                &
                                          ! 231 : OXYL1 + NO2 = ORGNIT
     &  RR( 5.6e-11,  0.0,     0.0,  1),                                &
                                          ! 232 : OH + CH3COCH=CHCHO = C
     &  RR( 9.0e-12,  0.0,     0.0,  1),                                &
                                          ! 233 : CH3COCH(OH)CH(O2)CHO +
!                                                              HO2 + CH3
     &  RR(5.60e-12,  0.0,     0.0,  1),                                &
                                          ! 234 : OH + TOLUENE = MEMALD
     &  RR( 3.6e-13,  0.0,     0.0,  1),                                &
                                          ! 235 : OH + TOLUENE = TOLP1,
     &  RR(1.37e-11,  0.0,     0.0,  1),                                &
                                          ! 236 : OH + O-XYLENE = ORGNIT
     &  RR( 2.7e-12,  0.0,     0.0,  1),                                &
                                          ! 237 : OH + ORGNIT = MEMALD +
     &  RR( 7.0e-14,  0.0,     0.0,  1),                                &
                                          ! 238 : NO3 + ORGNIT = MEMALD
     &  zero_rate,                                                      &
     &  RR( 1.0e-11,  0.0,     0.0,  1),                                &
                                          ! 240 : OXYL1 + HO2 = MGLYOX +
     &  RR( 1.0e-13,  0.0,     0.0,  1),                                &
                                          ! 241 : MEMALDIAL1 + CH3O2 =
!                                                               2HO2 + H
     &  RR( 5.0e-13,  0.0,     0.0,  1),                                &
                                          ! 242 : RO2IP1 + CH3O2 = 2HO2
     &  RR( 2.0e-12,  0.0,     0.0,  1),                                &
                                          ! 243 : RO2IP2 + CH3O2 = 2HO2
     &  RR(2.45e-13,  0.0,  1250.0,  3),                                &
                                          ! 244 : RO2IP1 + HO2 = ISOPOOH
     &  RR( 4.2e-11,  0.0,     0.0,  1),                                &
                                          ! 245 : ISOPOOH + OH = MVK + H
     &  RR(2.23e-13,  0.0,  1250.0,  3),                                &
                                          ! 246 : RO2IP2 + HO2 = MVKOOH
     &  RR(5.77e-11,  0.0,     0.0,  1)/) ! 247 : MVKOOH + OH = MGLYOX +
      rk(248:266) = (/                                                  &
     &  RR(1.72e-11,  0.0,     0.0,  1),                                &
                                          ! 248 : MGLYOX + OH = CH3COO2
     &  RR(1.14e-11,  0.0,     0.0,  1),                                &
                                          ! 249 : GLYOX + OH = HO2 + 2 C
     &  zero_rate,                                                      &
     &  RR(2.54e-11,  0.0,   410.0,  3),                                &
                                          ! 251 : OH + C5H8 = (HO)C5H8O2
     &  RR(2.08e-12,  0.0,   180.0,  3),                                &
                                          ! 252 : (HO)C5H8O2 + NO = MVK
     &  RR(4.13e-12,  0.0,   452.0,  3),                                &
                                          ! 253 : OH + MVK = (HO)MVKO2 +
     &  RR(2.46e-12,  0.0,   180.0,  3),                                &
                                          ! 254 : (HO)MVKO2 + NO = CH3CO
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  RR( 5.2e6,    0.0, -3650.0,  5),                                &
                                          ! 260 : HSO3(-) + H2O2(aq) = H
     &  RR( 4.2e5,    0.0, -4131.0,  5),                                &
                                          ! 261 : HSO3(-) + O3 = H(+) +
     &  RR( 1.5e9,    0.0,  -996.0,  5),                                &
                                          ! 262 : SO3(2-) + O3 = SO4(2-)
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate /)                      ! 266
      rk(267:285) = (/                                                  &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
                                          ! 270
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
                                          ! 275
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
                                          ! 280
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate /)                      ! 285
      rk(286:295) = (/                                                  &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  RR(1.51368e-07,0.0,    0.0,  1),                                &
                                          ! 290 : Be-7 radioactive decay
     &  RR(2.09465e-06,0.0,    0.0,  1),                                &
                                          ! 291 : Rn-222 = Pb-210 rad dk
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate,                                                      &
     &  zero_rate /)                      ! 295
!
! Pressure-dependent reactions
!
!                  klo khi res  Fc
      pdep = (/ PD(  4,  5,  5, 0.6),                                   &
                                          ! O + NO2 + M = NO2 + M
     &          PD( 20, 46, 20, 0.6),                                   &
                                          ! NO2 + NO3 + M = N2O5 + M
     &          PD( 21, 47, 21, 0.6),                                   &
                                          ! NO2 + OH + M = HNO3 + M
     &          PD( 22, 48, 22, 0.6),                                   &
                                          ! NO2 + HO2 + M = HO2NO2 + M
     &          PD( 29, 49, 29, 0.0),                                   &
                                          ! N2O5 + M = NO2 + NO3 + M
     &          PD( 39, 52, 39, 0.6),                                   &
                                          ! OH + SO2 + M = HOSO2 + M
     &          PD( 76, 77, 77, 0.3),                                   &
                                          ! CH3COO2 + NO2 + M = CH3COO2N
     &          PD( 78, 82, 78, 0.3),                                   &
                                          ! CH3COO2NO2 + M = CH3COO2 + N
     &          PD(108,109,109, 0.7),                                   &
                                          ! OH + C2H4 + M (+O2) = HOC2H4
     &          PD(125,134,125, 0.5) /)   ! OH + C3H6 + M (+O2) = HOC3H6
!
      first = .FALSE.
!
      END IF  ! If first = .TRUE.
!
      rc = 0.0
!
      DO i = 1, nr
        r = rk(i)%r
        ar = rk(i)%a
        n = rk(i)%n
        ea = rk(i)%e
        IF (r == 1) THEN
          DO j = 1, chunk
            rc(j,i) = ar
          END DO
        ELSE IF (r == 2) THEN
          DO j = 1, chunk
            rc(j,i) = ar * ((tc(j)*tfc) ** n)
          END DO
        ELSE IF (r == 3) THEN
          DO j = 1, chunk
            rc(j,i) = ar * exp(ea / tc(j))
          END DO
        ELSE IF (r == 4) THEN
          DO j = 1, chunk
            rc(j,i) = ar * ((tc(j)*tfc) ** n) * exp(ea / tc(j))
          END DO
        ELSE IF (r == 5) THEN
          DO j = 1, chunk
            fac1 = (1.0/tc(j)) - taq
            rc(j,i) = ar * exp(ea * fac1)
          END DO
        END IF
!
      END DO
!
! Following lines may be uncommented to check reaction types are correct
!
!     write(6,*) 'Checking CHEMCO for faulty rate constants'
!     DO i = 1, nr
!       IF (r == 1) THEN
!         IF (n /= 0.0 .OR. ea /= 0.0) WRITE(6,*) i
!       END IF
!       IF (r == 2) THEN
!         IF (n == 0.0 .OR. ea /= 0.0) WRITE(6,*) i
!       END IF
!       IF (r == 3 .OR. r == 5) THEN
!         IF (n /= 0.0 .OR. ea == 0.0) WRITE(6,*) i
!       END IF
!       IF (r == 4) THEN
!         IF (n == 0.0 .OR. ea == 0.0) WRITE(6,*) i
!       END IF
!     END DO
!
!     write(6,*) 'Finished checking CHEMCO'
!
! Pressure-dependent reactions
!
      DO k = 1, npdep
        nklo = pdep(k)%klo                          ! Reaction no. for k
        nkhi = pdep(k)%khi                          ! Reaction no. for k
        nres = pdep(k)%res                          ! Reaction no. for r
        f = pdep(k)%fc
        IF (f > 0.01) THEN ! f is zero for R29 to flag this reaction
          fc = f
          DO j = 1, chunk
            rc(j,nklo) = rc(j,nklo) * m(j)          ! klo * M
            brn = 0.75 - 1.27 * LOG10(fc)
            fac1 = rc(j,nklo) / rc(j,nkhi)          ! khi
            fac2 = rc(j,nklo) / (1.0 + fac1)
            fac3 = 1.0 + ((LOG10(fac1) / brn) ** 2)
            rc(j,nres) = fac2 * (fc ** (1.0 / fac3))
          END DO
        ELSE
          DO j = 1, chunk
            fc = EXP(-tc(j)/250.0) + EXP(-1050.0/tc(j))      ! For R29
            rc(j,nklo) = rc(j,nklo) * m(j)          ! klo * M
            brn = 0.75 - 1.27 * LOG10(fc)
            fac1 = rc(j,nklo) / rc(j,nkhi)          ! khi
            fac2 = rc(j,nklo) / (1.0 + fac1)
            fac3 = 1.0 + ((LOG10(fac1) / brn) ** 2)
            rc(j,nres) = fac2 * (fc ** (1.0 / fac3))
          END DO
        END IF
      END DO
!
! Do remaining calculations for certain reactions. This section
! needs to be edited by hand
!
      DO j = 1, chunk
!
! Reaction   7 : O(1D) + M
        rc(j,7) = (o2frac*rc(j,6) + n2frac*rc(j,7)) * m(j)
!
! Reaction  23 : HO2NO2 (+ M) -> HO2 + NO2
!   needs reverse pressure-dependent rate constant
        rc(j,23) = rc(j,22) / rc(j,23)
!
! Reaction  35 : OH + HNO3
        fac1 = rc(j,51) * m(j)
        rc(j,35) = rc(j,35) + fac1 / (1.0 + fac1 / rc(j,50))
!
! Reaction  36 : HO2 + HO2 (+ M, H2O)
        rc(j,36) = (rc(j,36) + rc(j,37) * m(j)) *                       &
     &    (1.0 + rc(j,38) * h2o(j))
!
! Reactions  61 and 62 : CH3O2 + CH3O2
        fac1 = rc(j,62) / (1.0 + rc(j,62))
        rc(j,62) = rc(j,61) * (1.0 - fac1)  ! -> CH3OH + HCHO + O2
        rc(j,61) = rc(j,61) * fac1          ! -> 2 HCHO + 2 HO2
!
! Reaction  70 : OH + CO
        rc(j,70) = rc(j,70) + (rc(j,69) * m(j))
!
! Reactions  74 and 80 : CH3O2 + CH3COO2
        fac1 = rc(j,74) / (1.0 + rc(j,74))
        rc(j,74) = rc(j,80) * (1.0-fac1)    ! -> 2 HCHO + O2
        rc(j,80) = rc(j,80) * fac1          ! -> HCHO + HO2 + CH3O2 + CO
!
! Reaction 132 : CH2OO + H2O
        rc(j,132) = rc(j,132) * h2o(j)
!
! Reaction 220 : CH3SO2 + O2 -> CH3O2 + SO2
        rc(j,220) = rc(j,220) * one_over_e
!
! Reaction 222 : CH3SO3 + O2 -> CH3O2 + SA
        rc(j,222) = rc(j,222) * one_over_e
!
! Set up O2 concentration
        fac1 = m(j) * o2frac
!
! Reaction   1 : O + O2 (+ M) = O3 + M
        rc(j,1) = rc(j,1) * m(j) * fac1
!
! Reaction 215 : OH + DMS = DMSO + HO2
        rc(j,215) = (rc(j,211) * fac1) / (1.0 + rc(j,215) * fac1)
!
      END DO
!
      END SUBROUTINE CHEMCO
#endif
