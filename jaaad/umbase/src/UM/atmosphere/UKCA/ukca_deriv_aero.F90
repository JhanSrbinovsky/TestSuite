#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    To evaluate concentrations from rate coefficients, J values,
!    and dry and wet deposition rates.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Current Code Owner:       Colin Johnson/Jamie Rae
!                            Fiona O'Connor
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
      SUBROUTINE UKCA_DERIV_AERO(nr_therm, nr_phot, n_be_calls,         &
                                 n_pnts, ndiags, rc, wdep, ddep,        &
                                 dj, k_dms, h2o, m, o2, vol,            &
                                 dts, L_ddep, L_wdep, L_flux,           &
                                 y, d,                                  &
                                 delta_SO2_wetox_H2O2,                  &
                                 delta_SO2_wetox_O3,                    &
                                 delta_SO2_dryox_OH)

      USE UKCA_CONSTANTS,   ONLY: avogadro
      IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "parvars.h"
#include "typsize.h"

      INTEGER, INTENT(IN) :: nr_therm      ! No of thermal reactions
      INTEGER, INTENT(IN) :: nr_phot       ! No of photolytic reactions
      INTEGER, INTENT(IN) :: n_be_calls    ! No diagnostics
      INTEGER, INTENT(IN) :: n_pnts        ! No of points
      INTEGER, INTENT(IN) :: ndiags        ! No diagnostics

      LOGICAL, INTENT(IN) :: L_ddep(jpspec)  ! T if dry deposits
      LOGICAL, INTENT(IN) :: L_wdep(jpspec)  ! T if wet deposits
      LOGICAL, INTENT(IN) :: L_flux          ! T for reaction diagnostics

      REAL, INTENT(IN) :: rc(theta_field_size,nr_therm) ! Rate coeffs
      REAL, INTENT(IN) :: wdep(theta_field_size,jpspec) ! Wet dep rates
      REAL, INTENT(IN) :: ddep(theta_field_size,jpspec) ! Dry dep rates
      REAL, INTENT(IN) :: dj(theta_field_size,jppj)     ! Phot rates
      REAL, INTENT(IN) :: k_dms(theta_field_size,5)     ! DMS rate coeffs
      REAL, INTENT(IN) :: h2o(theta_field_size)         ! H2O
      REAL, INTENT(IN) :: m(theta_field_size)           ! Air density
      REAL, INTENT(IN) :: o2(theta_field_size)          ! Oxygen
      REAL, INTENT(IN) :: vol(theta_field_size)         ! Cell volume
      REAL, INTENT(IN) :: dts                           ! timestep

! Species concentrations
      REAL, INTENT(INOUT) :: y(theta_field_size,jpspec)

! Fluxes of SO2 by oxidation pathway
      REAL, INTENT(OUT) :: delta_SO2_wetox_H2O2(theta_field_size)
      REAL, INTENT(OUT) :: delta_SO2_wetox_O3(theta_field_size)
      REAL, INTENT(OUT) :: delta_SO2_dryox_OH(theta_field_size)
! Chemical fluxes
      REAL, INTENT(OUT) :: d(theta_field_size,ndiags)

! Local variables
      REAL :: p1(theta_field_size)        ! Production and loss terms
      REAL :: p2(theta_field_size)
      REAL :: p(theta_field_size)
      REAL :: l(theta_field_size)
      REAL :: l1(theta_field_size)
      REAL :: l2(theta_field_size)
      REAL :: l3(theta_field_size)
      REAL :: r1(theta_field_size)
      REAL :: r2(theta_field_size)
      REAL :: f_so2(theta_field_size)     ! fraction of dms to SO2
      REAL :: f_so4(theta_field_size)     ! fraction of dms to SO4
      REAL :: f_msa(theta_field_size)     ! fraction of dms to MSA
      REAL :: yp(theta_field_size,jpspec) ! initial concentrations


      INTEGER, PARAMETER :: nit=8   ! No of iterations
      INTEGER :: i,j,n              ! loop counters
      INTEGER :: icnt_dd            ! counter
      INTEGER :: icnt_wd            ! counter

      DO n = 1, n_be_calls

        DO j = 1,jpspec
          DO i = 1, n_pnts
            yp(i,j) = y(i,j)
          END DO
        END DO

!      Set constant fields

       YP(:,12) = 5.0e-7*M(:)     ! H2
       YP(:,16) = 350.0e-6*M(:)   ! CO2
       YP(:,19) = H2O(:)          ! H2O
       YP(:,22) = O2(:)           ! O2
       YP(:,23) = M(:) - O2(:)    ! N2

       Y (:,12) = YP(:,12)        ! H2
       Y (:,16) = YP(:,16)        ! CO2
       Y (:,19) = YP(:,19)        ! H2O
       Y (:,22) = YP(:,22)        ! O2
       Y (:,23) = YP(:,23)        ! N2

! Initialise sulphur oxidation fluxes
       delta_SO2_wetox_H2O2(:) = 0.0
       delta_SO2_wetox_O3(:) = 0.0
       delta_SO2_dryox_OH(:) = 0.0

! Iteration start
       DO I=1,NIT

!
!   This section written automatically by MECH9GEN from the file ukca_tropch3
!   with 129 equations and 46 species.
!
!          O(3P)        Y( 1)
      P(:) =                                                            &
      +(DJ(:,15) *Y(:,3 ))                                              &
      +(DJ(:,12) *Y(:,5 ))        +(DJ(:,13) *Y(:,22)*2.00)             &
      +(RC(:,76) *Y(:,10)*Y(:,10))+(DJ(:,10) *Y(:,6 ))                  &
      +(RC(:,49) *Y(:,2 )*Y(:,23))+(RC(:,50) *Y(:,2 )*Y(:,22))
      L(:) = DDEP(:,1) + WDEP(:,1)                                      &
      +(RC(:,44) *Y(:,3 ))+(RC(:,83) *Y(:,6 ))+(RC(:,96) *Y(:,22))
      Y(:, 1) = P(:)/L(:)
!
!          O(1D)        Y( 2)
      P(:) =                                                            &
      +(DJ(:,14) *Y(:,3 ))
      L(:) = DDEP(:,2) + WDEP(:,2)                                      &
      +(RC(:,48) *Y(:,19))+(RC(:,49) *Y(:,23))+(RC(:,50) *Y(:,22))      &
      +(RC(:,45) *Y(:,14))+(RC(:,46) *Y(:,14))+(RC(:,47) *Y(:,14))
      Y(:, 2) = P(:)/L(:)
!
!          OH           Y(10)
      P(:) =                                                            &
      +(DJ(:,1)  *Y(:,34))        +(DJ(:,1)  *Y(:,39))                  &
      +(DJ(:,17) *Y(:,21))        +(DJ(:,1)  *Y(:,33))                  &
      +(DJ(:,6)  *Y(:,9 ))        +(DJ(:,1)  *Y(:,20))                  &
      +(DJ(:,1)  *Y(:,26))        +(DJ(:,2)  *Y(:,13)*2.00)             &
      +(RC(:,82) *Y(:,10)*Y(:,34))+(RC(:,85) *Y(:,45)*Y(:,44))          &
      +(RC(:,72) *Y(:,10)*Y(:,39))+(RC(:,80) *Y(:,10)*Y(:,33))          &
      +(RC(:,57) *Y(:,10)*Y(:,26))+(RC(:,66) *Y(:,10)*Y(:,20))          &
      +(RC(:,45) *Y(:,2 )*Y(:,14))+(RC(:,48) *Y(:,2 )*Y(:,19)*2.00)     &
      +(RC(:,3)  *Y(:,11)*Y(:,3 ))+(RC(:,10) *Y(:,11)*Y(:,28))          &
      +(RC(:,1)  *Y(:,11)*Y(:,4 ))+(RC(:,2)  *Y(:,11)*Y(:,5 ))
      L(:) = DDEP(:,10) + WDEP(:,10)                                    &
      +(RC(:,99) *Y(:,10))                                              &
      +(RC(:,97) *Y(:,4 ))+(RC(:,98) *Y(:,6 ))+(RC(:,99) *Y(:,10))      &
      +(RC(:,80) *Y(:,33))+(RC(:,81) *Y(:,34))+(RC(:,82) *Y(:,34))      &
      +(RC(:,77) *Y(:,29))+(RC(:,78) *Y(:,40))+(RC(:,79) *Y(:,33))      &
      +(RC(:,75) *Y(:,3 ))+(RC(:,76) *Y(:,10))+(RC(:,76) *Y(:,10))      &
      +(RC(:,72) *Y(:,39))+(RC(:,73) *Y(:,27))+(RC(:,74) *Y(:,5 ))      &
      +(RC(:,69) *Y(:,37))+(RC(:,70) *Y(:,37))+(RC(:,71) *Y(:,39))      &
      +(RC(:,66) *Y(:,20))+(RC(:,67) *Y(:,20))+(RC(:,68) *Y(:,41))      &
      +(RC(:,63) *Y(:,8 ))+(RC(:,64) *Y(:,9 ))+(RC(:,65) *Y(:,21))      &
      +(RC(:,60) *Y(:,13))+(RC(:,61) *Y(:,17))+(RC(:,62) *Y(:,11))      &
      +(RC(:,57) *Y(:,26))+(RC(:,58) *Y(:,26))+(RC(:,59) *Y(:,12))      &
      +(RC(:,54) *Y(:,30))+(RC(:,55) *Y(:,15))+(RC(:,56) *Y(:,35))      &
      +(RC(:,51) *Y(:,14))+(RC(:,52) *Y(:,24))+(RC(:,53) *Y(:,30))      &
      +(RC(:,103)*Y(:,48))                                              &
      +(RC(:,105)*Y(:,47))+(RC(:,106)*Y(:,47))
      Y(:,10) = P(:)/L(:)
!
!          O3           Y( 3)
      P(:) =                                                            &
      +(RC(:,96) *Y(:,1 )*Y(:,22))                                      &
      +(RC(:,9)  *Y(:,11)*Y(:,28))+(RC(:,14) *Y(:,11)*Y(:,36))
      L(:) = DDEP(:,3) + WDEP(:,3)                                      &
      +(DJ(:,15) )                                                      &
      +(RC(:,44) *Y(:,1 ))+(RC(:,75) *Y(:,10))+(DJ(:,14) )              &
      +(RC(:,3)  *Y(:,11))+(RC(:,37) *Y(:,4 ))+(RC(:,38) *Y(:,6 ))
      Y(:, 3) = (YP(:, 3)+DTS*P(:))/(1.0+DTS*L(:))
!
!          NO           Y( 4)
      P(:) =                                                            &
      +(DJ(:,11) *Y(:,5 ))        +(DJ(:,17) *Y(:,21))                  &
      +(RC(:,83) *Y(:,1 )*Y(:,6 ))+(DJ(:,10) *Y(:,6 ))
      L(:) = DDEP(:,4) + WDEP(:,4)                                      &
      +(RC(:,36) *Y(:,5 ))+(RC(:,37) *Y(:,3 ))+(RC(:,97) *Y(:,10))      &
      +(RC(:,30) *Y(:,32))+(RC(:,32) *Y(:,36))+(RC(:,34) *Y(:,38))      &
      +(RC(:,23) *Y(:,25))+(RC(:,26) *Y(:,28))+(RC(:,28) *Y(:,31))      &
      +(RC(:,1)  *Y(:,11))+(RC(:,16) *Y(:,18))+(RC(:,17) *Y(:,18))
      Y(:, 4) = (YP(:, 4)+DTS*P(:))/(1.0+DTS*L(:))
!
!          NO2          Y( 6)
      P(:) =                                                            &
      +(DJ(:,20) *Y(:,41))                                              &
      +(DJ(:,16) *Y(:,29))        +(DJ(:,16) *Y(:,40))                  &
      +(DJ(:,9)  *Y(:,7 ))        +(DJ(:,12) *Y(:,5 ))                  &
      +(DJ(:,5)  *Y(:,8 ))        +(DJ(:,6)  *Y(:,9 ))                  &
      +(RC(:,94) *Y(:,7 ))        +(RC(:,101)*Y(:,40))                  &
      +(RC(:,91) *Y(:,8 ))        +(RC(:,93) *Y(:,29))                  &
      +(RC(:,77) *Y(:,10)*Y(:,29))+(RC(:,78) *Y(:,10)*Y(:,40))          &
      +(RC(:,68) *Y(:,10)*Y(:,41))+(RC(:,74) *Y(:,10)*Y(:,5 ))          &
      +(RC(:,63) *Y(:,10)*Y(:,8 ))+(RC(:,65) *Y(:,10)*Y(:,21))          &
      +(RC(:,36) *Y(:,4 )*Y(:,5 )*2.00)+(RC(:,37) *Y(:,4 )*Y(:,3 ))     &
      +(RC(:,34) *Y(:,38)*Y(:,4 ))+(RC(:,35) *Y(:,38)*Y(:,5 ))          &
      +(RC(:,32) *Y(:,36)*Y(:,4 ))+(RC(:,33) *Y(:,36)*Y(:,5 ))          &
      +(RC(:,30) *Y(:,32)*Y(:,4 ))+(RC(:,31) *Y(:,32)*Y(:,5 ))          &
      +(RC(:,28) *Y(:,31)*Y(:,4 ))+(RC(:,29) *Y(:,31)*Y(:,5 ))          &
      +(RC(:,26) *Y(:,28)*Y(:,4 ))+(RC(:,27) *Y(:,28)*Y(:,5 ))          &
      +(RC(:,23) *Y(:,25)*Y(:,4 ))+(RC(:,24) *Y(:,25)*Y(:,5 ))          &
      +(RC(:,16) *Y(:,18)*Y(:,4 ))+(RC(:,18) *Y(:,18)*Y(:,5 ))          &
      +(RC(:,1)  *Y(:,11)*Y(:,4 ))+(RC(:,2)  *Y(:,11)*Y(:,5 ))
      L(:) = DDEP(:,6) + WDEP(:,6)                                      &
      +(RC(:,100)*Y(:,36))+(DJ(:,10) )                                  &
      +(RC(:,92) *Y(:,28))+(RC(:,95) *Y(:,5 ))+(RC(:,98) *Y(:,10))      &
      +(RC(:,38) *Y(:,3 ))+(RC(:,83) *Y(:,1 ))+(RC(:,90) *Y(:,11))
      Y(:, 6) = (YP(:, 6)+DTS*P(:))/(1.0+DTS*L(:))
!
!          NO3          Y( 5)
!      P(:) =
!     &+(RC(:,94) *Y(:,7 ))        +(DJ(:,9)  *Y(:,7 ))
!     &+(RC(:,38) *Y(:,6 )*Y(:,3 ))+(RC(:,64) *Y(:,10)*Y(:,9 ))
!      L(:) = DDEP(:,5) + WDEP(:,5)
!     &+(DJ(:,11) )        +(DJ(:,12) )
!     &+(RC(:,42) *Y(:,37))+(RC(:,74) *Y(:,10))+(RC(:,95) *Y(:,6 ))
!     &+(RC(:,39) *Y(:,17))+(RC(:,40) *Y(:,27))+(RC(:,41) *Y(:,35))
!     &+(RC(:,33) *Y(:,36))+(RC(:,35) *Y(:,38))+(RC(:,36) *Y(:,4 ))
!     &+(RC(:,27) *Y(:,28))+(RC(:,29) *Y(:,31))+(RC(:,31) *Y(:,32))
!     &+(RC(:,2)  *Y(:,11))+(RC(:,18) *Y(:,18))+(RC(:,24) *Y(:,25))
!     &+(RC(:,107) *Y(:,47))
!      Y(:, 5) = (YP(:, 5)+DTS*P(:))/(1.0+DTS*L(:))
!
!          N2O5         Y( 7)
!      P(:) =
!     &+(RC(:,95) *Y(:,6 )*Y(:,5 ))
!      L(:) = DDEP(:,7) + WDEP(:,7)
!     &+(RC(:,43) *Y(:,19))+(RC(:,94) )        +(DJ(:,9)  )
!      Y(:, 7) = (YP(:, 7)+DTS*P(:))/(1.0+DTS*L(:))
!
!        NO3/N2O5  Y(:, 5)/Y(:, 7)
       P1(:) =                                                          &
      +(RC(:,38) *Y(:,6 )*Y(:,3 ))+(RC(:,64) *Y(:,10)*Y(:,9 ))
       L(:) = DDEP(:,5) + WDEP(:,5)                                     &
      +(DJ(:,11) )        +(DJ(:,12) )                                  &
      +(RC(:,42) *Y(:,37))+(RC(:,74) *Y(:,10))+(RC(:,95) *Y(:,6 ))      &
      +(RC(:,39) *Y(:,17))+(RC(:,40) *Y(:,27))+(RC(:,41) *Y(:,35))      &
      +(RC(:,33) *Y(:,36))+(RC(:,35) *Y(:,38))+(RC(:,36) *Y(:,4 ))      &
      +(RC(:,27) *Y(:,28))+(RC(:,29) *Y(:,31))+(RC(:,31) *Y(:,32))      &
      +(RC(:,2)  *Y(:,11))+(RC(:,18) *Y(:,18))+(RC(:,24) *Y(:,25))      &
      +(RC(:,107)*Y(:,47))
      P2(:) = 0.0
      R1(:) = RC(:,94) + DJ(:,9)
      R2(:) = RC(:,95) *Y(:,6 )
      L1(:) = DDEP(:,7) + WDEP(:,7)                                     &
      +(RC(:,43) *Y(:,19))+(RC(:,94) )        +(DJ(:,9)  )
      L2(:) = 1.0+L(:)*DTS
      L3(:) = 1.0+L1(:)*DTS
      Y(:,5) = (L3(:)*(YP(:,5)+P1(:)*DTS)+R1(:)*DTS*YP(:,7))/           &
             ((L3(:)*L2(:))-R1(:)*R2(:)*DTS**2.0)
      Y(:,7) = (YP(:,7) + P2(:)*DTS + R2(:)*DTS*Y(:,5))/L3(:)
!
!          HO2NO2       Y( 8)
      P(:) =                                                            &
      +(RC(:,90) *Y(:,11)*Y(:,6 ))
      L(:) = DDEP(:,8) + WDEP(:,8)                                      &
      +(RC(:,63) *Y(:,10))+(RC(:,91) )        +(DJ(:,5)  )
      Y(:, 8) = (YP(:, 8)+DTS*P(:))/(1.0+DTS*L(:))
!
!          HONO2        Y( 9)
      P(:) =                                                            &
      +(RC(:,43) *Y(:,7 )*Y(:,19)*2.00)+(RC(:,98) *Y(:,10)*Y(:,6 ))     &
      +(RC(:,41) *Y(:,5 )*Y(:,35))+(RC(:,42) *Y(:,5 )*Y(:,37))          &
      +(RC(:,39) *Y(:,5 )*Y(:,17))+(RC(:,40) *Y(:,5 )*Y(:,27))          &
      +(RC(:,107)*Y(:,5 )*Y(:,47))
      L(:) = DDEP(:,9) + WDEP(:,9)                                      &
      +(RC(:,64) *Y(:,10))+(DJ(:,6)  )
      Y(:, 9) = (YP(:, 9)+DTS*P(:))/(1.0+DTS*L(:))
!
!          HO2          Y(11)
      P(:) =                                                            &
      +(DJ(:,1)  *Y(:,34))        +(DJ(:,20) *Y(:,41))                  &
      +(DJ(:,18) *Y(:,35))        +(DJ(:,1)  *Y(:,33))                  &
      +(DJ(:,7)  *Y(:,27))        +(DJ(:,1)  *Y(:,20))                  &
      +(DJ(:,3)  *Y(:,17)*2.00)        +(DJ(:,5)  *Y(:,8 ))             &
      +(RC(:,91) *Y(:,8 ))        +(DJ(:,1)  *Y(:,26))                  &
      +(RC(:,75) *Y(:,10)*Y(:,3 ))+(RC(:,84) *Y(:,46)*Y(:,44))          &
      +(RC(:,61) *Y(:,10)*Y(:,17))+(RC(:,74) *Y(:,10)*Y(:,5 ))          &
      +(RC(:,59) *Y(:,10)*Y(:,12))+(RC(:,60) *Y(:,10)*Y(:,13))          &
      +(RC(:,47) *Y(:,2 )*Y(:,14)*2.00)+(RC(:,55) *Y(:,10)*Y(:,15))     &
      +(RC(:,31) *Y(:,32)*Y(:,5 ))+(RC(:,39) *Y(:,5 )*Y(:,17))          &
      +(RC(:,29) *Y(:,31)*Y(:,5 ))+(RC(:,30) *Y(:,32)*Y(:,4 ))          &
      +(RC(:,25) *Y(:,25)*Y(:,28))+(RC(:,28) *Y(:,31)*Y(:,4 ))          &
      +(RC(:,23) *Y(:,25)*Y(:,4 ))+(RC(:,24) *Y(:,25)*Y(:,5 ))          &
      +(RC(:,20) *Y(:,18)*Y(:,18)*2.00)+(RC(:,21) *Y(:,18)*Y(:,28))     &
      +(RC(:,16) *Y(:,18)*Y(:,4 ))+(RC(:,18) *Y(:,18)*Y(:,5 ))
      L(:) = DDEP(:,11) + WDEP(:,11)                                    &
      +(RC(:,89) *Y(:,11))+(RC(:,90) *Y(:,6 ))                          &
      +(RC(:,15) *Y(:,38))+(RC(:,62) *Y(:,10))+(RC(:,89) *Y(:,11))      &
      +(RC(:,12) *Y(:,32))+(RC(:,13) *Y(:,36))+(RC(:,14) *Y(:,36))      &
      +(RC(:,9)  *Y(:,28))+(RC(:,10) *Y(:,28))+(RC(:,11) *Y(:,31))      &
      +(RC(:,6)  *Y(:,18))+(RC(:,7)  *Y(:,25))+(RC(:,8)  *Y(:,28))      &
      +(RC(:,4)  *Y(:,11))+(RC(:,4)  *Y(:,11))+(RC(:,5)  *Y(:,18))      &
      +(RC(:,1)  *Y(:,4 ))+(RC(:,2)  *Y(:,5 ))+(RC(:,3)  *Y(:,3 ))
      Y(:,11) = (YP(:,11)+DTS*P(:))/(1.0+DTS*L(:))
!
!          H2           Y(12)
!      P(:) =                                                            &
!     &+(RC(:,46) *Y(:,2 )*Y(:,14))+(DJ(:,4)  *Y(:,17))
!      L(:) = DDEP(:,12) + WDEP(:,12)
!     &+(RC(:,59) *Y(:,10))
!      Y(:,12) = (YP(:,12)+DTS*P(:))/(1.0+DTS*L(:))
       Y(:,12) = YP(:,12)
!
!          H2O2         Y(13)
      P(:) =                                                            &
      +(RC(:,99) *Y(:,10)*Y(:,10))                                      &
      +(RC(:,4)  *Y(:,11)*Y(:,11))+(RC(:,89) *Y(:,11)*Y(:,11))
      L(:) = DDEP(:,13) + WDEP(:,13)                                    &
      +(RC(:,60) *Y(:,10))+(DJ(:,2)  )                                  &
      +(RC(:,104)*Y(:,48))
      Y(:,13) = (YP(:,13)+DTS*P(:))/(1.0+DTS*L(:))
!
!          CH4          Y(14)
      P(:) =                                                            &
      +(DJ(:,8)  *Y(:,27))
      L(:) = DDEP(:,14) + WDEP(:,14)                                    &
      +(RC(:,51) *Y(:,10))                                              &
      +(RC(:,45) *Y(:,2 ))+(RC(:,46) *Y(:,2 ))+(RC(:,47) *Y(:,2 ))
      Y(:,14) = (YP(:,14)+DTS*P(:))/(1.0+DTS*L(:))
!
!          CO           Y(15)
      P(:) =                                                            &
      +(DJ(:,18) *Y(:,35))                                              &
      +(DJ(:,7)  *Y(:,27))        +(DJ(:,8)  *Y(:,27))                  &
      +(DJ(:,3)  *Y(:,17))        +(DJ(:,4)  *Y(:,17))                  &
      +(RC(:,39) *Y(:,5 )*Y(:,17))+(RC(:,61) *Y(:,10)*Y(:,17))
      L(:) = DDEP(:,15) + WDEP(:,15)                                    &
      +(RC(:,55) *Y(:,10))
      Y(:,15) = (YP(:,15)+DTS*P(:))/(1.0+DTS*L(:))
!
!          CO2          Y(16)
!      P(:) =                                                            &
!     &+(RC(:,32) *Y(:,36)*Y(:,4 ))+(RC(:,33) *Y(:,36)*Y(:,5 ))
!     &+(RC(:,26) *Y(:,28)*Y(:,4 ))+(RC(:,27) *Y(:,28)*Y(:,5 ))
!      L(:) = DDEP(:,16) + WDEP(:,16)
      Y(:,16) = YP(:,16)
!
!          HCHO         Y(17)
      P(:) =                                                            &
      +(DJ(:,20) *Y(:,41))                                              &
      +(DJ(:,1)  *Y(:,20))        +(DJ(:,1)  *Y(:,39))                  &
      +(RC(:,68) *Y(:,10)*Y(:,41))+(RC(:,77) *Y(:,10)*Y(:,29))          &
      +(RC(:,47) *Y(:,2 )*Y(:,14))+(RC(:,66) *Y(:,10)*Y(:,20))          &
      +(RC(:,35) *Y(:,38)*Y(:,5 ))+(RC(:,46) *Y(:,2 )*Y(:,14))          &
      +(RC(:,22) *Y(:,18)*Y(:,28))+(RC(:,34) *Y(:,38)*Y(:,4 ))          &
      +(RC(:,20) *Y(:,18)*Y(:,18)*2.00)+(RC(:,21) *Y(:,18)*Y(:,28))     &
      +(RC(:,18) *Y(:,18)*Y(:,5 ))+(RC(:,19) *Y(:,18)*Y(:,18))          &
      +(RC(:,6)  *Y(:,11)*Y(:,18))+(RC(:,16) *Y(:,18)*Y(:,4 ))
      L(:) = DDEP(:,17) + WDEP(:,17)                                    &
      +(DJ(:,4)  )                                                      &
      +(RC(:,39) *Y(:,5 ))+(RC(:,61) *Y(:,10))+(DJ(:,3)  )
      Y(:,17) = (YP(:,17)+DTS*P(:))/(1.0+DTS*L(:))
!
!          MeOO         Y(18)
      P(:) =                                                            &
      +(DJ(:,7)  *Y(:,27))        +(DJ(:,19) *Y(:,37))                  &
      +(RC(:,51) *Y(:,10)*Y(:,14))+(RC(:,67) *Y(:,10)*Y(:,20))          &
      +(RC(:,27) *Y(:,28)*Y(:,5 ))+(RC(:,45) *Y(:,2 )*Y(:,14))          &
      +(RC(:,25) *Y(:,25)*Y(:,28))+(RC(:,26) *Y(:,28)*Y(:,4 ))          &
      +(RC(:,10) *Y(:,11)*Y(:,28))+(RC(:,21) *Y(:,18)*Y(:,28))
      L(:) = DDEP(:,18) + WDEP(:,18)                                    &
      +(RC(:,21) *Y(:,28))+(RC(:,22) *Y(:,28))                          &
      +(RC(:,19) *Y(:,18))+(RC(:,20) *Y(:,18))+(RC(:,20) *Y(:,18))      &
      +(RC(:,17) *Y(:,4 ))+(RC(:,18) *Y(:,5 ))+(RC(:,19) *Y(:,18))      &
      +(RC(:,5)  *Y(:,11))+(RC(:,6)  *Y(:,11))+(RC(:,16) *Y(:,4 ))
      Y(:,18) = (YP(:,18)+DTS*P(:))/(1.0+DTS*L(:))
!
!          H2O          Y(19)
!      P(:) =                                                            &
!     &+(RC(:,81) *Y(:,10)*Y(:,34))+(RC(:,86) *Y(:,43)*Y(:,19))          &
!     &+(RC(:,79) *Y(:,10)*Y(:,33))+(RC(:,80) *Y(:,10)*Y(:,33))          &
!     &+(RC(:,77) *Y(:,10)*Y(:,29))+(RC(:,78) *Y(:,10)*Y(:,40))          &
!     &+(RC(:,73) *Y(:,10)*Y(:,27))+(RC(:,76) *Y(:,10)*Y(:,10))          &
!     &+(RC(:,70) *Y(:,10)*Y(:,37))+(RC(:,71) *Y(:,10)*Y(:,39))          &
!     &+(RC(:,68) *Y(:,10)*Y(:,41))+(RC(:,69) *Y(:,10)*Y(:,37))          &
!     &+(RC(:,66) *Y(:,10)*Y(:,20))+(RC(:,67) *Y(:,10)*Y(:,20))          &
!     &+(RC(:,64) *Y(:,10)*Y(:,9 ))+(RC(:,65) *Y(:,10)*Y(:,21))          &
!     &+(RC(:,62) *Y(:,10)*Y(:,11))+(RC(:,63) *Y(:,10)*Y(:,8 ))          &
!     &+(RC(:,60) *Y(:,10)*Y(:,13))+(RC(:,61) *Y(:,10)*Y(:,17))          &
!     &+(RC(:,58) *Y(:,10)*Y(:,26))+(RC(:,59) *Y(:,10)*Y(:,12))          &
!     &+(RC(:,56) *Y(:,10)*Y(:,35))+(RC(:,57) *Y(:,10)*Y(:,26))          &
!     &+(RC(:,53) *Y(:,10)*Y(:,30))+(RC(:,54) *Y(:,10)*Y(:,30))          &
!     &+(RC(:,51) *Y(:,10)*Y(:,14))+(RC(:,52) *Y(:,10)*Y(:,24))
!      L(:) = DDEP(:,19) + WDEP(:,19)                                    &
!     &+(RC(:,43) *Y(:,7 ))+(RC(:,48) *Y(:,2 ))+(RC(:,86) *Y(:,43))
      Y(:,19) = YP(:,19)
!
!          MeOOH        Y(20)
      P(:) =                                                            &
      +(RC(:,5)  *Y(:,11)*Y(:,18))
      L(:) = DDEP(:,20) + WDEP(:,20)                                    &
      +(RC(:,66) *Y(:,10))+(RC(:,67) *Y(:,10))+(DJ(:,1)  )
      Y(:,20) = (YP(:,20)+DTS*P(:))/(1.0+DTS*L(:))
!
!          HONO         Y(21)
      P(:) =                                                            &
      +(RC(:,97) *Y(:,10)*Y(:,4 ))
      L(:) = DDEP(:,21) + WDEP(:,21)                                    &
      +(RC(:,65) *Y(:,10))+(DJ(:,17) )
      Y(:,21) = (YP(:,21)+DTS*P(:))/(1.0+DTS*L(:))
!
!          O2           Y(22)
!      P(:) =                                                            &
!     &+(DJ(:,15) *Y(:,44))
!     &+(DJ(:,15) *Y(:,3 ))        +(DJ(:,14) *Y(:,44))
!     &+(DJ(:,11) *Y(:,5 ))        +(DJ(:,14) *Y(:,3 ))
!     &+(RC(:,88) *Y(:,43)*Y(:,22))+(RC(:,89) *Y(:,11)*Y(:,11))
!     &+(RC(:,84) *Y(:,46)*Y(:,44))+(RC(:,85) *Y(:,45)*Y(:,44))
!     &+(RC(:,75) *Y(:,10)*Y(:,3 ))+(RC(:,83) *Y(:,1 )*Y(:,6 ))
!     &+(RC(:,44) *Y(:,1 )*Y(:,3 )*2.00)+(RC(:,50) *Y(:,2 )*Y(:,22))
!     &+(RC(:,3)  *Y(:,11)*Y(:,3 ))+(RC(:,13) *Y(:,11)*Y(:,36))
!      L(:) = DDEP(:,22) + WDEP(:,22)
!     &+(RC(:,102)*Y(:,42))+(DJ(:,13) )
!     &+(RC(:,50) *Y(:,2 ))+(RC(:,88) *Y(:,43))+(RC(:,96) *Y(:,1 ))
      Y(:,22) = YP(:,22)
!
!          N2           Y(23)
!      P(:) =                                                            &
!     &+(RC(:,49) *Y(:,2 )*Y(:,23))+(RC(:,87) *Y(:,43)*Y(:,23))
!      L(:) = DDEP(:,23) + WDEP(:,23)
!     &+(RC(:,49) *Y(:,2 ))+(RC(:,87) *Y(:,43))
      Y(:,23) = YP(:,23)
!
!          C2H6         Y(24)
      P(:) = 0.0
      L(:) = DDEP(:,24) + WDEP(:,24)                                    &
      +(RC(:,52) *Y(:,10))
      Y(:,24) = (YP(:,24)+DTS*P(:))/(1.0+DTS*L(:))
!
!          EtOO         Y(25)
      P(:) =                                                            &
      +(DJ(:,18) *Y(:,35))                                              &
      +(RC(:,52) *Y(:,10)*Y(:,24))+(RC(:,58) *Y(:,10)*Y(:,26))          &
      +(RC(:,32) *Y(:,36)*Y(:,4 ))+(RC(:,33) *Y(:,36)*Y(:,5 ))
      L(:) = DDEP(:,25) + WDEP(:,25)                                    &
      +(RC(:,25) *Y(:,28))                                              &
      +(RC(:,7)  *Y(:,11))+(RC(:,23) *Y(:,4 ))+(RC(:,24) *Y(:,5 ))
      Y(:,25) = (YP(:,25)+DTS*P(:))/(1.0+DTS*L(:))
!
!          EtOOH        Y(26)
      P(:) =                                                            &
      +(RC(:,7)  *Y(:,11)*Y(:,25))
      L(:) = DDEP(:,26) + WDEP(:,26)                                    &
      +(RC(:,57) *Y(:,10))+(RC(:,58) *Y(:,10))+(DJ(:,1)  )
      Y(:,26) = (YP(:,26)+DTS*P(:))/(1.0+DTS*L(:))
!
!          MeCHO        Y(27)
      P(:) =                                                            &
      +(RC(:,78) *Y(:,10)*Y(:,40))+(DJ(:,1)  *Y(:,26))                  &
      +(RC(:,25) *Y(:,25)*Y(:,28))+(RC(:,57) *Y(:,10)*Y(:,26))          &
      +(RC(:,23) *Y(:,25)*Y(:,4 ))+(RC(:,24) *Y(:,25)*Y(:,5 ))
      L(:) = DDEP(:,27) + WDEP(:,27)                                    &
      +(DJ(:,8)  )                                                      &
      +(RC(:,40) *Y(:,5 ))+(RC(:,73) *Y(:,10))+(DJ(:,7)  )
      Y(:,27) = (YP(:,27)+DTS*P(:))/(1.0+DTS*L(:))
!
!          MeCO3        Y(28)
      P(:) =                                                            &
      +(DJ(:,19) *Y(:,37))        +(DJ(:,1)  *Y(:,39))                  &
      +(RC(:,93) *Y(:,29))        +(DJ(:,16) *Y(:,29))                  &
      +(RC(:,40) *Y(:,5 )*Y(:,27))+(RC(:,73) *Y(:,10)*Y(:,27))          &
      +(RC(:,34) *Y(:,38)*Y(:,4 ))+(RC(:,35) *Y(:,38)*Y(:,5 ))
      L(:) = DDEP(:,28) + WDEP(:,28)                                    &
      +(RC(:,26) *Y(:,4 ))+(RC(:,27) *Y(:,5 ))+(RC(:,92) *Y(:,6 ))      &
      +(RC(:,21) *Y(:,18))+(RC(:,22) *Y(:,18))+(RC(:,25) *Y(:,25))      &
      +(RC(:,8)  *Y(:,11))+(RC(:,9)  *Y(:,11))+(RC(:,10) *Y(:,11))
      Y(:,28) = (YP(:,28)+DTS*P(:))/(1.0+DTS*L(:))
!
!          PAN          Y(29)
      P(:) =                                                            &
      +(RC(:,92) *Y(:,28)*Y(:,6 ))
      L(:) = DDEP(:,29) + WDEP(:,29)                                    &
      +(RC(:,77) *Y(:,10))+(RC(:,93) )        +(DJ(:,16) )
      Y(:,29) = (YP(:,29)+DTS*P(:))/(1.0+DTS*L(:))
!
!          C3H8         Y(30)
      P(:) = 0.0
      L(:) = DDEP(:,30) + WDEP(:,30)                                    &
      +(RC(:,53) *Y(:,10))+(RC(:,54) *Y(:,10))
      Y(:,30) = (YP(:,30)+DTS*P(:))/(1.0+DTS*L(:))
!
!          n-PrOO       Y(31)
      P(:) =                                                            &
      +(RC(:,53) *Y(:,10)*Y(:,30))+(RC(:,79) *Y(:,10)*Y(:,33))
      L(:) = DDEP(:,31) + WDEP(:,31)                                    &
      +(RC(:,11) *Y(:,11))+(RC(:,28) *Y(:,4 ))+(RC(:,29) *Y(:,5 ))
      Y(:,31) = (YP(:,31)+DTS*P(:))/(1.0+DTS*L(:))
!
!          i-PrOO       Y(32)
      P(:) =                                                            &
      +(RC(:,54) *Y(:,10)*Y(:,30))+(RC(:,81) *Y(:,10)*Y(:,34))
      L(:) = DDEP(:,32) + WDEP(:,32)                                    &
      +(RC(:,12) *Y(:,11))+(RC(:,30) *Y(:,4 ))+(RC(:,31) *Y(:,5 ))
      Y(:,32) = (YP(:,32)+DTS*P(:))/(1.0+DTS*L(:))
!
!          n-PrOOH      Y(33)
      P(:) =                                                            &
      +(RC(:,11) *Y(:,11)*Y(:,31))
      L(:) = DDEP(:,33) + WDEP(:,33)                                    &
      +(RC(:,79) *Y(:,10))+(RC(:,80) *Y(:,10))+(DJ(:,1)  )
      Y(:,33) = (YP(:,33)+DTS*P(:))/(1.0+DTS*L(:))
!
!          i-PrOOH      Y(34)
      P(:) =                                                            &
      +(RC(:,12) *Y(:,11)*Y(:,32))
      L(:) = DDEP(:,34) + WDEP(:,34)                                    &
      +(RC(:,81) *Y(:,10))+(RC(:,82) *Y(:,10))+(DJ(:,1)  )
      Y(:,34) = (YP(:,34)+DTS*P(:))/(1.0+DTS*L(:))
!
!          EtCHO        Y(35)
      P(:) =                                                            &
      +(RC(:,80) *Y(:,10)*Y(:,33))+(DJ(:,1)  *Y(:,33))                  &
      +(RC(:,28) *Y(:,31)*Y(:,4 ))+(RC(:,29) *Y(:,31)*Y(:,5 ))
      L(:) = DDEP(:,35) + WDEP(:,35)                                    &
      +(RC(:,41) *Y(:,5 ))+(RC(:,56) *Y(:,10))+(DJ(:,18) )
      Y(:,35) = (YP(:,35)+DTS*P(:))/(1.0+DTS*L(:))
!
!          EtCO3        Y(36)
      P(:) =                                                            &
      +(RC(:,101)*Y(:,40))        +(DJ(:,16) *Y(:,40))                  &
      +(RC(:,41) *Y(:,5 )*Y(:,35))+(RC(:,56) *Y(:,10)*Y(:,35))
      L(:) = DDEP(:,36) + WDEP(:,36)                                    &
      +(RC(:,33) *Y(:,5 ))+(RC(:,100)*Y(:,6 ))                          &
      +(RC(:,13) *Y(:,11))+(RC(:,14) *Y(:,11))+(RC(:,32) *Y(:,4 ))
      Y(:,36) = (YP(:,36)+DTS*P(:))/(1.0+DTS*L(:))
!
!          Me2CO        Y(37)
      P(:) =                                                            &
      +(RC(:,82) *Y(:,10)*Y(:,34))+(DJ(:,1)  *Y(:,34))                  &
      +(RC(:,30) *Y(:,32)*Y(:,4 ))+(RC(:,31) *Y(:,32)*Y(:,5 ))
      L(:) = DDEP(:,37) + WDEP(:,37)                                    &
      +(DJ(:,19) )                                                      &
      +(RC(:,42) *Y(:,5 ))+(RC(:,69) *Y(:,10))+(RC(:,70) *Y(:,10))
      Y(:,37) = (YP(:,37)+DTS*P(:))/(1.0+DTS*L(:))
!
!          MeCOCH2O     Y(38)
      P(:) =                                                            &
      +(RC(:,70) *Y(:,10)*Y(:,37))+(RC(:,71) *Y(:,10)*Y(:,39))          &
      +(RC(:,42) *Y(:,5 )*Y(:,37))+(RC(:,69) *Y(:,10)*Y(:,37))
      L(:) = DDEP(:,38) + WDEP(:,38)                                    &
      +(RC(:,15) *Y(:,11))+(RC(:,34) *Y(:,4 ))+(RC(:,35) *Y(:,5 ))
      Y(:,38) = (YP(:,38)+DTS*P(:))/(1.0+DTS*L(:))
!
!          MeCOCH2O     Y(39)
      P(:) =                                                            &
      +(RC(:,15) *Y(:,11)*Y(:,38))
      L(:) = DDEP(:,39) + WDEP(:,39)                                    &
      +(RC(:,71) *Y(:,10))+(RC(:,72) *Y(:,10))+(DJ(:,1)  )
      Y(:,39) = (YP(:,39)+DTS*P(:))/(1.0+DTS*L(:))
!
!          PPAN         Y(40)
      P(:) =                                                            &
      +(RC(:,100)*Y(:,36)*Y(:,6 ))
      L(:) = DDEP(:,40) + WDEP(:,40)                                    &
      +(RC(:,78) *Y(:,10))+(RC(:,101))        +(DJ(:,16) )
      Y(:,40) = (YP(:,40)+DTS*P(:))/(1.0+DTS*L(:))
!
!          MeONO2       Y(41)
      P(:) =                                                            &
      +(RC(:,17) *Y(:,18)*Y(:,4 ))
      L(:) = DDEP(:,41) + WDEP(:,41)                                    &
      +(RC(:,68) *Y(:,10))+(DJ(:,20) )
      Y(:,41) = (YP(:,41)+DTS*P(:))/(1.0+DTS*L(:))
!
!          O(3P)S       Y(42)
      P(:) =                                                            &
      +(DJ(:,15) *Y(:,44))                                              &
      +(RC(:,87) *Y(:,43)*Y(:,23))+(RC(:,88) *Y(:,43)*Y(:,22))
      L(:) = DDEP(:,42) + WDEP(:,42)                                    &
      +(RC(:,102)*Y(:,22))
      Y(:,42) = (YP(:,42)+DTS*P(:))/(1.0+DTS*L(:))
!
!          O(1D)S       Y(43)
      P(:) =                                                            &
      +(DJ(:,14) *Y(:,44))
      L(:) = DDEP(:,43) + WDEP(:,43)                                    &
      +(RC(:,86) *Y(:,19))+(RC(:,87) *Y(:,23))+(RC(:,88) *Y(:,22))
      Y(:,43) = (YP(:,43)+DTS*P(:))/(1.0+DTS*L(:))
!
!          O3S          Y(44)
      P(:) =                                                            &
      +(RC(:,102)*Y(:,42)*Y(:,22))
      L(:) = DDEP(:,44) + WDEP(:,44)                                    &
      +(DJ(:,15) )                                                      &
      +(RC(:,84) *Y(:,46))+(RC(:,85) *Y(:,45))+(DJ(:,14) )
      Y(:,44) = (YP(:,44)+DTS*P(:))/(1.0+DTS*L(:))
!
!          OHS          Y(45)
      P(:) = 0.0
      L(:) = DDEP(:,45) + WDEP(:,45)                                    &
      +(RC(:,85) *Y(:,44))
      Y(:,45) = (YP(:,45)+DTS*P(:))/(1.0+DTS*L(:))
!
!          HO2S         Y(46)
      P(:) = 0.0
      L(:) = DDEP(:,46) + WDEP(:,46)                                    &
      +(RC(:,84) *Y(:,44))
      Y(:,46) = (YP(:,46)+DTS*P(:))/(1.0+DTS*L(:))
!
! Fraction of DMS oxidation to products
      f_so2(:) = k_dms(:,1)/                                            &
                 (k_dms(:,1)+k_dms(:,2)*y(:,3)+k_dms(:,3)*y(:,6))
      f_so4(:) = (1.0-f_so2)*k_dms(:,5)/(k_dms(:,5)+k_dms(:,4)*y(:,11))
      f_msa(:) = (1.0-f_so2)*k_dms(:,4)*y(:,11)/                        &
                 (k_dms(:,5)+k_dms(:,4)*y(:,11))
!
!          DMS          Y(47)
      P(:) = 0.0
      L(:) = DDEP(:,47) + WDEP(:,47)                                    &
      +(RC(:,105) *Y(:,10))+(RC(:,106) *Y(:,10))                        &
      +(RC(:,107) *Y(:,5))
      Y(:,47) = (YP(:,47)+DTS*P(:))/(1.0+DTS*L(:))
!
!          SO2          Y(48)
      P(:) = 0.0                                                        &
      +(RC(:,105)*Y(:,10)+RC(:,106)*Y(:,10)                             &
      + RC(:,107)*Y(:,5))*Y(:,47)*f_so2(:)
      L(:) = 0.0 + DDEP(:,48) + WDEP(:,48)                              &
      +(RC(:,103) * Y(:,10))                                            &
      +(RC(:,104) * Y(:,13))
      Y(:,48) = (YP(:,48)+DTS*P(:))/(1.0+DTS*L(:))
!
!          H2SO4        Y(49)
      P(:) = 0.0                                                        &
      +(RC(:,103) * Y(:,48) * Y(:,10))                                  &
      +(RC(:,105)*Y(:,10)+RC(:,106)*Y(:,10)                             &
      + RC(:,107)*Y(:,5))*Y(:,47)*f_so4(:)
      L(:) = 0.0
      Y(:,49) = (YP(:,49)+DTS*P(:))/(1.0+DTS*L(:))
!
!          MSA          Y(50)
      P(:) = 0.0                                                        &
      +(RC(:,105)*Y(:,10)+RC(:,106)*Y(:,10)                             &
      + RC(:,107)*Y(:,5))*Y(:,47)*f_msa(:)
      L(:) = 0.0
      Y(:,50) = (YP(:,50)+DTS*P(:))/(1.0+DTS*L(:))

!      iteration loop stop
      ENDDO

! Calculate flux terms at end of iteration.

! Fluxes to aqueous sulphate (always required for MODE):
      delta_SO2_wetox_H2O2(:)= delta_SO2_wetox_H2O2(:) +                &
                               rc(:,104)*y(:,13)*y(:,48)*dts
      delta_SO2_wetox_O3(:) = 0.0
      delta_SO2_dryox_OH(:) = delta_SO2_dryox_OH(:) +                   &
                                rc(:,103)*y(:,10)*y(:,48)*dts

      D(:,:) = 0.0
      IF (L_flux) THEN

! R1    HO2       +NO        =OH        +NO2
        D(:,1) = D(:,1) + RC(:,1)*Y(:,11)*Y(:,4)*DTS*Vol(:)

! R2    HO2       +NO3       =OH        +NO2
        D(:,2) = D(:,2) + RC(:,2)*Y(:,11)*Y(:,5)*DTS*Vol(:)

! R3    HO2       +O3        =OH        +O2
        D(:,3) = D(:,3) + RC(:,3)*Y(:,11)*Y(:,3)*DTS*Vol(:)

! R4    HO2       +HO2       =H2O2
        D(:,4) = D(:,4) + RC(:,4)*Y(:,11)*Y(:,11)*DTS*Vol(:)

! R5    HO2       +MeOO      =MeOOH
        D(:,5) = D(:,5) + RC(:,5)*Y(:,11)*Y(:,18)*DTS*Vol(:)

! R6    HO2       +MeOO      =HCHO
        D(:,6) = D(:,6) + RC(:,6)*Y(:,11)*Y(:,18)*DTS*Vol(:)

! R7    HO2       +EtOO      =EtOOH
        D(:,7) = D(:,7) + RC(:,7)*Y(:,11)*Y(:,25)*DTS*Vol(:)

! R8    HO2       +MeCO3     =NULL
        D(:,8) = D(:,8) + RC(:,8)*Y(:,11)*Y(:,28)*DTS*Vol(:)

! R9    HO2       +MeCO3     =O3
        D(:,9) = D(:,9) + RC(:,9)*Y(:,11)*Y(:,28)*DTS*Vol(:)

! R10   HO2       +MeCO3     =OH + MeOO
        D(:,10) = D(:,10) + RC(:,10)*Y(:,11)*Y(:,28)*DTS*Vol(:)

! R11   HO2       +n-PrOO    =n-PrOOH
        D(:,11) = D(:,11) + RC(:,11)*Y(:,11)*Y(:,31)*DTS*Vol(:)

! R12   HO2       +i-PrOO    =i-PrOOH
        D(:,12) = D(:,12) + RC(:,12)*Y(:,11)*Y(:,32)*DTS*Vol(:)

! R13   HO2       +EtCO3     =O2
        D(:,13) = D(:,13) + RC(:,13)*Y(:,11)*Y(:,36)*DTS*Vol(:)

! R14   HO2       +EtCO3     =O3
        D(:,14) = D(:,14) + RC(:,14)*Y(:,11)*Y(:,36)*DTS*Vol(:)

! R15   HO2       +MeCOCH2OO =MeCOCH2OOH
        D(:,15) = D(:,15) + RC(:,15)*Y(:,11)*Y(:,38)*DTS*Vol(:)

! R16   MeOO      +NO        =HO2       +HCHO      +NO2
        D(:,16) = D(:,16) + RC(:,16)*Y(:,18)*Y(:,4)*DTS*Vol(:)

! R17   MeOO      +NO        =MeONO2
        D(:,17) = D(:,17) + RC(:,17)*Y(:,18)*Y(:,4)*DTS*Vol(:)

! R18   MeOO      +NO3       =HO2       +HCHO      +NO2
        D(:,18) = D(:,18) + RC(:,18)*Y(:,18)*y(:,5)*DTS*Vol(:)

! R19   MeOO      +MeOO      =HCHO
        D(:,19) = D(:,19) + RC(:,19)*Y(:,18)*Y(:,18)*DTS*Vol(:)

! R20   MeOO      +MeOO      =HO2       +HCHO
        D(:,20) = D(:,20) + RC(:,20)*Y(:,18)*Y(:,18)*DTS*Vol(:)

! R21   MeOO      +MeCO3     =HO2       +HCHO      +MeOO
        D(:,21) = D(:,21) + RC(:,21)*Y(:,18)*Y(:,28)*DTS*Vol(:)

! R22   MeOO      +MeCO3     =HCHO
        D(:,22) = D(:,22) + RC(:,22)*Y(:,18)*Y(:,28)*DTS*Vol(:)

! R23   EtOO      +NO        =MeCHO     +HO2       +NO2
        D(:,23) = D(:,23) + RC(:,23)*Y(:,25)*Y(:,4)*DTS*Vol(:)

! R24   EtOO      +NO3       =MeCHO     +HO2       +NO2
        D(:,24) = D(:,24) + RC(:,24)*Y(:,25)*y(:,5)*DTS*Vol(:)

! R25   EtOO      +MeCO3     =MeCHO     +HO2       +MeOO
        D(:,25) = D(:,25) + RC(:,25)*Y(:,25)*Y(:,28)*DTS*Vol(:)

! R26   MeCO3     +NO        =MeOO      +CO2       +NO2
        D(:,26) = D(:,26) + RC(:,26)*Y(:,28)*Y(:,4)*DTS*Vol(:)

! R27   MeCO3     +NO3       =MeOO      +CO2       +NO2
        D(:,27) = D(:,27) + RC(:,27)*Y(:,28)*Y(:,5)*DTS*Vol(:)

! R28   n-PrOO    +NO        =EtCHO     +HO2       +NO2
        D(:,28) = D(:,28) + RC(:,28)*Y(:,31)*Y(:,4)*DTS*Vol(:)

! R29   n-PrOO    +NO3       =EtCHO     +HO2       +NO2
        D(:,29) = D(:,29) + RC(:,29)*Y(:,31)*y(:,5)*DTS*Vol(:)

! R30   i-PrOO    +NO        =Me2CO     +HO2       +NO2
        D(:,30) = D(:,30) + RC(:,30)*Y(:,32)*Y(:,4)*DTS*Vol(:)

! R31   i-PrOO    +NO3       =Me2CO     +HO2       +NO2
        D(:,31) = D(:,31) + RC(:,31)*Y(:,32)*Y(:,5)*DTS*Vol(:)

! R32   EtCO3     +NO        =EtOO      +CO2       +NO2*Vol(:)
        D(:,32) = D(:,32) + RC(:,32)*Y(:,36)*Y(:,4)*DTS*Vol(:)

! R33   EtCO3     +NO3       =EtOO      +CO2       +NO2
        D(:,33) = D(:,33) + RC(:,33)*Y(:,36)*Y(:,5)*DTS*Vol(:)

! R34   MeCOCH2OO +NO        =MeCO3     +HCHO      +NO2
        D(:,34) = D(:,34) + RC(:,34)*Y(:,38)*Y(:,4)*DTS*Vol(:)

! R35   MeCOCH2OO +NO3       =MeCO3     +HCHO      +NO2
        D(:,35) = D(:,35) + RC(:,35)*Y(:,38)*Y(:,5)*DTS*Vol(:)

! R36   NO        +NO3       =NO2
        D(:,36) = D(:,36) + RC(:,36)*Y(:,4)*Y(:,5)*DTS*Vol(:)

! R37   NO        +O3        =NO2
        D(:,37) = D(:,37) + RC(:,37)*Y(:,4)*Y(:,3)*DTS*Vol(:)

! R38   NO2       +O3        =NO3
        D(:,38) = D(:,38) + RC(:,38)*Y(:,6)*Y(:,3)*DTS*Vol(:)

! R39   NO3       +HCHO      =HONO2     +HO2       +CO
        D(:,39) = D(:,39) + RC(:,39)*Y(:,5)*Y(:,17)*DTS*Vol(:)

! R40   NO3       +MeCHO     =HONO2     +MeCO3
        D(:,40) = D(:,40) + RC(:,40)*Y(:,5)*Y(:,27)*DTS*Vol(:)

! R41   NO3       +EtCHO     =HONO2     +EtCO3
        D(:,41) = D(:,41) + RC(:,41)*Y(:,5)*Y(:,35)*DTS*Vol(:)

! R42   NO3       +Me2CO     =HONO2     +MeCOCH2OO
        D(:,42) = D(:,42) + RC(:,42)*Y(:,5)*Y(:,37)*DTS*Vol(:)

! R43   N2O5      +H2O       =HONO2
        D(:,43) = D(:,43) + RC(:,43)*Y(:,7)*H2O(:)*DTS*Vol(:)

! R44   O(3P)     +O3        =O2
        D(:,44) = D(:,44) + RC(:,44)*Y(:,1)*Y(:,3)*DTS*Vol(:)

! R45   O(1D)     +CH4       =OH        +MeOO
        D(:,45) = D(:,45) + RC(:,45)*Y(:,2)*Y(:,14)*DTS*Vol(:)

! R46   O(1D)     +CH4       =HCHO      +H2
        D(:,46) = D(:,46) + RC(:,46)*Y(:,2)*Y(:,14)*DTS*Vol(:)

! R47   O(1D)     +CH4       =HCHO      +HO2
        D(:,47) = D(:,47) + RC(:,47)*Y(:,2)*Y(:,14)*DTS*Vol(:)

! R48   O(1D)     +H2O       =OH
        D(:,48) = D(:,48) + RC(:,48)*Y(:,2)*H2O(:)*DTS*Vol(:)

! R49   O(1D)     +N2        =O(3P)     +N2
        D(:,49) = D(:,49) + RC(:,49)*Y(:,2)*Y(:,23)*DTS*Vol(:)

! R50   O(1D)     +O2        =O(3P)     +O2
        D(:,50) = D(:,50) + RC(:,50)*Y(:,2)*Y(:,22)*DTS*Vol(:)

! R51   OH        +CH4       =H2O       +MeOO
        D(:,51) = D(:,51) + RC(:,51)*Y(:,10)*Y(:,14)*DTS*Vol(:)

! R52   OH        +C2H6      =H2O       +EtOO
        D(:,52) = D(:,52) + RC(:,52)*Y(:,10)*Y(:,24)*DTS*Vol(:)

! R53   OH        +C3H8      =n-PrOO    +H2O
        D(:,53) = D(:,53) + RC(:,53)*Y(:,10)*Y(:,30)*DTS*Vol(:)

! R54   OH        +C3H8      =i-PrOO    +H2O
        D(:,54) = D(:,54) + RC(:,54)*Y(:,10)*Y(:,30)*DTS*Vol(:)

! R55   OH        +CO        =HO2
        D(:,55) = D(:,55) + RC(:,55)*Y(:,10)*Y(:,15)*DTS*Vol(:)

! R56   OH        +EtCHO     =H2O       +EtCO3
        D(:,56) = D(:,56) + RC(:,56)*Y(:,10)*Y(:,35)*DTS*Vol(:)

! R57   OH        +EtOOH     =H2O       +MeCHO     +OH
        D(:,57) = D(:,57) + RC(:,57)*Y(:,10)*Y(:,26)*DTS*Vol(:)

! R58   OH        +EtOOH     =H2O       +EtOO
        D(:,58) = D(:,58) + RC(:,58)*Y(:,10)*Y(:,26)*DTS*Vol(:)

! R59   OH        +H2        =H2O       +HO2
        D(:,59) = D(:,59) + RC(:,59)*Y(:,10)*Y(:,12)*DTS*Vol(:)

! R60   OH        +H2O2      =H2O       +HO2
        D(:,60) = D(:,60) + RC(:,60)*Y(:,10)*Y(:,13)*DTS*Vol(:)

! R61   OH        +HCHO      =H2O       +HO2       +CO
        D(:,61) = D(:,61) + RC(:,61)*Y(:,10)*Y(:,17)*DTS*Vol(:)

! R62   OH        +HO2       =H2O
        D(:,62) = D(:,62) + RC(:,62)*Y(:,10)*Y(:,11)*DTS*Vol(:)

! R63   OH        +HO2NO2    =H2O       +NO2
        D(:,63) = D(:,63) + RC(:,63)*Y(:,10)*Y(:,8)*DTS*Vol(:)

! R64   OH        +HONO2     =H2O       +NO3
        D(:,64) = D(:,64) + RC(:,64)*Y(:,10)*Y(:,9)*DTS*Vol(:)

! R65   OH        +HONO      =H2O       +NO2
        D(:,65) = D(:,65) + RC(:,65)*Y(:,10)*Y(:,21)*DTS*Vol(:)

! R66   OH        +MeOOH     =H2O       +HCHO      +OH
        D(:,66) = D(:,66) + RC(:,66)*Y(:,10)*Y(:,20)*DTS*Vol(:)

! R67   OH        +MeOOH     =H2O       +MeOO
        D(:,67) = D(:,67) + RC(:,67)*Y(:,10)*Y(:,20)*DTS*Vol(:)

! R68   OH        +MeONO2    =HCHO      +NO2       +H2O
        D(:,68) = D(:,68) + RC(:,68)*Y(:,10)*Y(:,41)*DTS*Vol(:)

! R69   OH        +Me2CO     =H2O       +MeCOCH2OO
        D(:,69) = D(:,69) + RC(:,69)*Y(:,10)*Y(:,37)*DTS*Vol(:)

! R70   OH        +Me2CO     =H2O       +MeCOCH2OO
        D(:,70) = D(:,70) + RC(:,70)*Y(:,10)*Y(:,37)*DTS*Vol(:)

! R71   OH        +MeCOCH2OOH=H2O       +MeCOCH2OO
        D(:,71) = D(:,71) + RC(:,71)*Y(:,10)*Y(:,39)*DTS*Vol(:)

! R72   OH        +MeCOCH2OOH=H2O
        D(:,72) = D(:,72) + RC(:,72)*Y(:,10)*Y(:,39)*DTS*Vol(:)

! R73   OH        +MeCHO     =H2O       +MeCO3
        D(:,73) = D(:,73) + RC(:,73)*Y(:,10)*Y(:,27)*DTS*Vol(:)

! R74   OH        +NO3       =HO2       +NO2
        D(:,74) = D(:,74) + RC(:,74)*Y(:,10)*Y(:,5)*DTS*Vol(:)

! R75   OH        +O3        =HO2       +O2
        D(:,75) = D(:,75) + RC(:,75)*Y(:,10)*Y(:,3)*DTS*Vol(:)

! R76   OH        +OH        =H2O       +O(3P)
        D(:,76) = D(:,76) + RC(:,76)*Y(:,10)*Y(:,10)*DTS*Vol(:)

! R77   OH        +PAN       =HCHO      +NO2       +H2O
        D(:,77) = D(:,77) + RC(:,77)*Y(:,10)*Y(:,29)*DTS*Vol(:)

! R78   OH        +PPAN      =MeCHO     +NO2       +H2O
        D(:,78) = D(:,78) + RC(:,78)*Y(:,10)*Y(:,40)*DTS*Vol(:)

! R79   OH        +n-PrOOH   =n-PrOO    +H2O
        D(:,79) = D(:,79) + RC(:,79)*Y(:,10)*Y(:,33)*DTS*Vol(:)

! R80   OH        +n-PrOOH   =EtCHO     +H2O       +OH
        D(:,80) = D(:,80) + RC(:,80)*Y(:,10)*Y(:,33)*DTS*Vol(:)

! R81   OH        +i-PrOOH   =i-PrOO    +H2O
        D(:,81) = D(:,81) + RC(:,81)*Y(:,10)*Y(:,34)*DTS*Vol(:)

! R82   OH        +i-PrOOH   =Me2CO     +OH
        D(:,82) = D(:,82) + RC(:,82)*Y(:,10)*Y(:,34)*DTS*Vol(:)

! R83   O3P       +NO2       =NO        +O2
        D(:,83) = D(:,83) + RC(:,83)*Y(:,1)*Y(:,6)*DTS*Vol(:)

! R84   HO2S      +O3S       =HO2S      +O2
        D(:,84) = D(:,84) + RC(:,84)*Y(:,46)*Y(:,44)*DTS*Vol(:)

! R85   OHS       +O3S       =HO2S      +O2
        D(:,85) = D(:,85) + RC(:,85)*Y(:,45)*Y(:,44)*DTS*Vol(:)

! R86   O(1D)S    +H2O       =H2O
        D(:,86) = D(:,86) + RC(:,86)*Y(:,43)*H2O(:)*DTS*Vol(:)

! R87   O(1D)S    +N2        =O(3P)S    +N2
        D(:,87) = D(:,87) + RC(:,87)*Y(:,43)*Y(:,23)*DTS*Vol(:)

! R88   O(1D)S    +O2        =O(3P)S    +O2
        D(:,88) = D(:,88) + RC(:,88)*Y(:,43)*Y(:,22)*DTS*Vol(:)

! R89   HO2       +HO2       =H2O2      +O2
        D(:,89) = D(:,89) + RC(:,89)*Y(:,11)*Y(:,11)*DTS*Vol(:)

! R90   HO2       +NO2       =HO2NO2
        D(:,90) = D(:,90) + RC(:,90)*Y(:,11)*y(:,6)*DTS*Vol(:)

! R91   HO2NO2    +M         =HO2       +NO2
        D(:,91) = D(:,91) + RC(:,91)*Y(:,8)*M(:)*DTS*Vol(:)

! R92   MeCO3     +NO2       =PAN
        D(:,92) = D(:,92) + RC(:,92)*Y(:,28)*Y(:,6)*DTS*Vol(:)

! R93   PAN       +M         =MeCO3     +NO2
        D(:,93) = D(:,93) + RC(:,93)*Y(:,29)*M(:)*DTS*Vol(:)

! R94   N2O5      +M         =NO2       +NO3
        D(:,94) = D(:,94) + RC(:,94)*Y(:,7)*M(:)*DTS*Vol(:)

! R95   NO2       +NO3       =N2O5
        D(:,95) = D(:,95) + RC(:,95)*Y(:,6)*Y(:,5)*DTS*Vol(:)

! R96   O(3P)     +O2        =O3
        D(:,96) = D(:,96) + RC(:,96)*Y(:,1)*Y(:,22)*DTS*Vol(:)

! R97   OH        +NO        =HONO
        D(:,97) = D(:,97) + RC(:,97)*Y(:,10)*Y(:,4)*DTS*Vol(:)

! R98   OH        +NO2       =HONO2
        D(:,98) = D(:,98) + RC(:,98)*Y(:,10)*Y(:,6)*DTS*Vol(:)

! R99   OH        +OH        =H2O2
        D(:,99) = D(:,99) + RC(:,99)*Y(:,10)*Y(:,10)*DTS*Vol(:)

! R100  EtCO3     +NO2       =PPAN
        D(:,100) = D(:,100) + RC(:,100)*Y(:,36)*Y(:,6)*DTS*Vol(:)

! R101  PPAN      +M         =EtCO3     +NO2
        D(:,101) = D(:,101) + RC(:,101)*Y(:,40)*M(:)*DTS*Vol(:)

! R102  O(3P)S    +O2        =O3S
        D(:,102) = D(:,102) + RC(:,102)*Y(:,42)*Y(:,22)*DTS*Vol(:)

! R103  SO2       +OH        =H2SO4
        D(:,103) = D(:,103)+(RC(:,103)*Y(:,48)*Y(:,10))*vol(:)

! R104  HSO3      +H2O2      =
        D(:,104) = D(:,104)+(RC(:,104)*Y(:,48)*Y(:,13))*vol(:)

! R105  DMS       +OH        =SO2
        D(:,105) = D(:,105)+(RC(:,105)*Y(:,10)+RC(:,106)*Y(:,10))*      &
                   Y(:,47)*f_so2(:)*vol(:)

! R106  DMS       +OH        =H2SO4
        D(:,106) = D(:,106)+(RC(:,105)*Y(:,10)+RC(:,106)*Y(:,10))*      &
                   Y(:,47)*f_so4(:)*vol(:)

! R107  DMS       +OH        =MSA
        D(:,107) = D(:,107)+(RC(:,105)*Y(:,10)+RC(:,106)*Y(:,10))*      &
                   Y(:,47)*f_msa(:)*vol(:)

! R108  DMS       +NO3       =SO2
        D(:,108) = D(:,108)+(RC(:,107)*Y(:,5)*Y(:,47)*f_so2(:))*vol(:)

! R109  DMS       +NO3       =H2SO4
        D(:,109) = D(:,109)+(RC(:,107)*Y(:,47)*Y(:,5)*f_so4(:))*vol(:)

! R110  DMS       +NO3       =MSA
        D(:,110) = D(:,110)+(RC(:,107)*Y(:,47)*Y(:,5)*f_msa(:))*vol(:)

!      --------------------
!      PHOTOLYTIC REACTIONS
!      --------------------
!
!       EtOOH     =MeCHO     +HO2       +OH
        D(:,nr_therm+1)=D(:,nr_therm+1)+DJ(:,1)*Y(:,26)*DTS*Vol(:)

!       H2O2      =OH
        D(:,nr_therm+2)=D(:,nr_therm+2)+DJ(:,2)*Y(:,13)*DTS*Vol(:)

!       HCHO      =HO2       +CO
        D(:,nr_therm+3)=D(:,nr_therm+3)+DJ(:,3)*Y(:,17)*DTS*Vol(:)

!       HCHO      =H2        +CO
        D(:,nr_therm+4)=D(:,nr_therm+4)+DJ(:,4)*Y(:,17)*DTS*Vol(:)

!       HO2NO2    =HO2       +NO2
        D(:,nr_therm+5)=D(:,nr_therm+5)+DJ(:,5)*Y(:,8)*DTS*Vol(:)

!       HONO2     =OH        +NO2
        D(:,nr_therm+6)=D(:,nr_therm+6)+DJ(:,6)*Y(:,9)*DTS*Vol(:)

!       MeCHO     =MeOO      +HO2       +CO
        D(:,nr_therm+7)=D(:,nr_therm+7)+DJ(:,7)*Y(:,27)*DTS*Vol(:)

!       MeCHO     =CH4       +CO
        D(:,nr_therm+8)=D(:,nr_therm+8)+DJ(:,8)*Y(:,27)*DTS*Vol(:)

!       MeOOH     =HO2       +HCHO      +OH
        D(:,nr_therm+9)=D(:,nr_therm+9)+DJ(:,1)*Y(:,20)*DTS*Vol(:)

!       N2O5      =NO3       +NO2
        D(:,nr_therm+10)=D(:,nr_therm+10)+DJ(:,9)*Y(:,7)*DTS*Vol(:)

!       NO2       =NO        +O(3P)
        D(:,nr_therm+11)=D(:,nr_therm+11)+DJ(:,10)*y(:,6)*Vol(:)

!       NO3       =NO        +O2
        D(:,nr_therm+12)=D(:,nr_therm+12)+DJ(:,11)*y(:,5)*DTS*Vol(:)

!       NO3       =NO2       +O(3P)
        D(:,nr_therm+13)=D(:,nr_therm+13)+DJ(:,12)*y(:,5)*DTS*Vol(:)

!       O2        =O(3P)
        D(:,nr_therm+14)=D(:,nr_therm+14)+DJ(:,13)*Y(:,22)*DTS*Vol(:)

!       O3        =O2        +O(1D)
        D(:,nr_therm+15)=D(:,nr_therm+15)+DJ(:,14)*Y(:,3)*DTS*Vol(:)

!       O3        =O2        +O(3P)
        D(:,nr_therm+16)=D(:,nr_therm+16)+DJ(:,15)*Y(:,3)*DTS*Vol(:)

!       PAN       =MeCO3     +NO2
        D(:,nr_therm+17)=D(:,nr_therm+17)+DJ(:,16)*Y(:,29)*DTS*Vol(:)

!       HONO      =OH        +NO
        D(:,nr_therm+18)=D(:,nr_therm+18)+DJ(:,17)*Y(:,21)*DTS*Vol(:)

!       EtCHO     =EtOO      +HO2       +CO
        D(:,nr_therm+19)=D(:,nr_therm+19)+DJ(:,18)*Y(:,35)*DTS*Vol(:)

!       Me2CO     =MeCO3     +MeOO
        D(:,nr_therm+20)=D(:,nr_therm+20)+DJ(:,19)*Y(:,37)*DTS*Vol(:)

!       n-PrOOH   =EtCHO     +HO2       +OH
        D(:,nr_therm+21)=D(:,nr_therm+21)+DJ(:,1)*Y(:,33)*DTS*Vol(:)

!       i-PrOOH   =Me2CO     +HO2       +OH
        D(:,nr_therm+22)=D(:,nr_therm+22)+DJ(:,1)*Y(:,34)*DTS*Vol(:)

!       MeCOCH2OOH=MeCO3     +HCHO      +OH
        D(:,nr_therm+23)=D(:,nr_therm+23)+DJ(:,1)*Y(:,39)*DTS*Vol(:)

!       PPAN      =EtCO3     +NO2
        D(:,nr_therm+24)=D(:,nr_therm+24)+DJ(:,16)*Y(:,40)*DTS*Vol(:)

!       MeONO2    =HO2       +HCHO      +NO2
        D(:,nr_therm+25)=D(:,nr_therm+25)+DJ(:,20)*Y(:,41)*DTS*Vol(:)

!       O3S       =O2        +O(1D)S
        D(:,nr_therm+26)=D(:,nr_therm+26)+DJ(:,14)*Y(:,44)*DTS*Vol(:)

!       O3S       =O2        +O(3P)S
        D(:,nr_therm+27)=D(:,nr_therm+27)+DJ(:,15)*Y(:,44)*DTS*Vol(:)
!
!      --------------------
!      DRY DEPOSITION FLUXES
!      --------------------
!
        icnt_dd = nr_therm+nr_phot
        DO j=1,jpspec
          IF (L_ddep(j)) THEN
           icnt_dd = icnt_dd + 1
           D(:,icnt_dd) = D(:,icnt_dd) + DDEP(:,j)*Y(:,j)*Vol(:)
         ENDIF
        ENDDO
!
!      --------------------
!      WET DEPOSITION FLUXES
!      --------------------
!
        icnt_wd = icnt_dd
        DO j=1,jpspec
          IF (L_wdep(j)) THEN
           icnt_wd = icnt_wd + 1
           D(:,icnt_wd) = D(:,icnt_wd) + WDEP(:,j)*Y(:,j)*Vol(:)
         ENDIF
        ENDDO
      ENDIF     ! L_flux

      END DO  ! n_be_calls

! Convert fluxes from molecules -> moles
      IF (L_flux) D(:,:) = D(:,:)/avogadro


      END SUBROUTINE UKCA_DERIV_AERO
#endif
