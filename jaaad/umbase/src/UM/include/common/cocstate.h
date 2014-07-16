!==================================================================
!    COMDECK COCSTATE
!    ----------------
!LL 4.4 Introduce coeff_T and Temp0 to enable temperature to be
!LL     calculated from theta and salinity (M. J. Bell)
!     5.5  25/02/03  Added coefficients for McDougall EOS. I. Stevens
      INTEGER,PARAMETER :: MAXLEV=99

      REAL::C(MAXLEV,9)      ! Polynomial coeffs (middle of box)
      REAL::CI(MAXLEV,9,2)   ! ------"---------- (top & bottom of box)
      REAL::coeff_T(MAXLEV,9)! Polynomial coeffs (middle of box) for T
      REAL::TO(MAXLEV)       ! Reference pot. temp. (middle of box)
      REAL::TOI(MAXLEV,2)    ! ----------"--------  (top & bottom of box
      REAL::SO(MAXLEV)       ! Reference salinities (middle of box)
      REAL::SOI(MAXLEV,2)    ! ----------"--------- (top & bottom of box
      REAL::SIGO(MAXLEV)     ! Reference density (middle of box)
      REAL::TempO(MAXLEV)    ! Reference temperatures  (middle of box)
      REAL p_dbar(maxlev)                                               &
     &,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11                            &
     &,two_a2,three_a3,two_a6,two_a8,two_a11                            &
     &,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12                        &
     &,two_b2,three_b3,four_b4,three_b7,onep5_b8                        &
     &,onep5_b9,two_b9,three_b11

      COMMON/CSTATE/ C,TO,SO,CI,TOI,SOI,sigo,coeff_T,TempO              &
     &,p_dbar                                                           &
     &,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11                            &
     &,two_a2,three_a3,two_a6,two_a8,two_a11                            &
     &,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12                        &
     &,two_b2,three_b3,four_b4,three_b7,onep5_b8                        &
     &,onep5_b9,two_b9,three_b11

! COCSTATE end
