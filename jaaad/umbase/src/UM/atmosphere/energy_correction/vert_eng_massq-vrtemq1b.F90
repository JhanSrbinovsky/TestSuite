#if defined(A14_1B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE VERT_ENG_MASSQ
!LL
!LL  PURPOSE : PART OF ENERGY CORRECTION SUITE OF ROUTINES
!LL            Given the prognostic input fields it integrates
!LL            vertically various fields required by the energy
!LL            correction.
!LL            Also if call from section 30  (climate diagnostics)
!LL            integrates qu fluxes etc
!LL
!LL
!LL  MODEL            MODIFICATION HISTORY:
!LL VERSION  DATE
!LL   5.1  1/03/00  New subroutine added for use by the energy
!LL                 correction and section 30 diagnostics to ensure
!LL                 calculations done in same way. R A Stratton
!LL   5.2 14/09/00  Correct potential energy calculation,
!LL                 reproducible results.
!LL   5.3 23/04/01  Take account of wet_model_level+1 in same way
!LL                 as flux_rho now does. R A Stratton
!LL   5.3 08/06/01  Remove duplicate type declaration.  A van der Wal
!     5.4 03/04/02  correction to rho_dry                  A. Malcolm
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL
!LL  DOCUMENTATION :  Energy correction documentation
!LL
!LLEND-----------------------------------------------------------------
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE VERT_ENG_MASSQ(                                        &
     &                      halo_i, halo_j, off_x, off_y                &
     &,                     global_row_length, proc_row_group           &
     &,                     at_extremity, n_proc, n_procx, n_procy      &
     &,                     neighbour                                   &
     &,                     mype                                        &
     &,                     row_length, rows, n_rows                    &
     &,                     model_domain                                &
     &,                     model_levels,wet_model_levels               &
     &,                     r_theta_levels,r_rho_levels                 &
     &,                     cos_theta_longitude,sin_theta_longitude     &
     &,                     theta ,u,v,w, rho_r2, q,qcl,qcf             &
     &,                     exner_theta_levels                          &
     &,                     Lqflux                                      &
! output fields
     &,                     rho_dry,dry_to_wet                          &
     &,                     vert_int_array,vert_qflux)


      IMPLICIT NONE

! constants required for calculations
#include "c_r_cp.h"
#include "c_g.h"
#include "c_lheat.h"
#include "c_a.h"
! Required for interpolation
#include "fldtype.h"

!----------------------------------------------------------------------
! INPUT variables
!----------------------------------------------------------------------
! MPP & halo info

      Integer off_x, off_y, halo_i, halo_j   ! halo info
      Integer                                                           &
     &  global_row_length                                               &
                           ! global row length
     &, n_proc                                                          &
                           ! Total number of processors
     &, n_procx                                                         &
                           ! Number of processors in longitude
     &, n_procy                                                         &
                           ! Number of processors in latitude
     &, proc_row_group                                                  &
                           ! Group id for processors on the same row
     &, neighbour(4)                                                    &
                       ! Array with the Ids of the four neighbours in
                       ! in the horizontal plane
     &, mype               ! processor number


      Logical                                                           &
     &  at_extremity(4)! Indicates if this processor is at north, south,
                       ! east or west of the processor grid

! Parameters
      Integer                                                           &
     &   PNorth,                                                        &
                      ! North processor address in the neighbor array
     &   PEast,                                                         &
                      ! East processor address in the neighbor array
     &   PSouth,                                                        &
                      ! South processor address in the neighbor array
     &   PWest,                                                         &
                      ! West processor address in the neighbor array
     &   NoDomain     ! Value in neighbor array if the domain has
                      !  no neighbor in this direction. Otherwise
                      !  the value will be the tid of the neighbor
      Parameter (                                                       &
     &   PNorth   = 1,                                                  &
     &   PEast    = 2,                                                  &
     &   PSouth   = 3,                                                  &
     &   PWest    = 4,                                                  &
     &   NoDomain = -1)

! input model dimensional info
      Integer                                                           &
     &   row_length                                                     &
                             ! number of points per row
     &,  rows                                                           &
                             ! number of rows on p grid
     &,  n_rows                                                         &
                             ! number of rows on v grid
     &,  model_domain                                                   &
                             ! global or limited area etc.
     &,  model_levels                                                   &
                             ! number of model levels
     &,  wet_model_levels    ! number of wet model levels

! Input radius of model levels
      REAL                                                              &
     &  R_THETA_LEVELS(1-halo_i:row_length+halo_i,                      &
     &            1-halo_j:rows+halo_j,0:model_levels)                  &
     &, R_RHO_LEVELS(1-halo_i:row_length+halo_i,                        &
     &            1-halo_j:rows+halo_j,model_levels)

! Input trigonometric functions (required by polar_wind)
      Real                                                              &
     &  cos_theta_longitude (row_length, rows)                          &
     &, sin_theta_longitude (row_length, rows)

!Input model data
      REAL                                                              &
     &  theta(1-off_x:row_length+off_x,                                 &
     &       1-off_y:rows+off_y, model_levels)                          &
                                                 ! theta
     &, U(1-off_x:row_length+off_x,                                     &
     &       1-off_y:rows+off_y, model_levels)                          &
                                                 !U COMPONENT OF WIND
     &, V(1-off_x:row_length+off_x,                                     &
     &       1-off_y:n_rows+off_y, model_levels)                        &
                                                 !V COMPONENT OF WIND
     &, w(1-off_x:row_length+off_x,                                     &
     &       1-off_y:rows+off_y,0:model_levels)                         &
                                                 ! w
     &, rho_r2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &        model_levels)                                             &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &      wet_model_levels)                                           &
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_model_levels)                                         &
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_model_levels)                                         &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)

!In
      logical                                                           &
     &  lqflux            ! true if q flux calculations required.

!In/OUT  Vertical integrals plus info on rho_dry

      REAL                                                              &
     &  vert_int_array(row_length,rows,10)                              &
                                              ! vertical integrals
     &, vert_qflux(row_length,rows,3)                                   &
                                              ! q flux integrals
     &, rho_dry(row_length,rows,model_levels)                           &
                                                ! rho dry x r^2
     &, dry_to_wet(row_length,rows,model_levels)

!----------------------------------------------------------------------
! local variables
!----------------------------------------------------------------------
! pointers for sum_array
      Integer                                                           &
     & ip_dry_mass                                                      &
                        ! dry mass
     &,ip_wet_mass                                                      &
                        ! wet mass
     &,ip_cvT                                                           &
                        ! cvT
     &,ip_gr                                                            &
                        ! gr
     &,ip_keu                                                           &
                        ! keu
     &,ip_kev                                                           &
                        ! kev
     &,ip_kew                                                           &
                        ! kew
     &,ip_q                                                             &
                        ! q
     &,ip_qcl                                                           &
                        ! qcl
     &,ip_qcf                                                           &
                        ! qcf
     &,ip_qu                                                            &
                        ! qu
     &,ip_qv                                                            &
                        ! qv
     &,ip_qw            ! qw

      Parameter (ip_dry_mass=1, ip_wet_mass=2, ip_cvT=3, ip_gr=4,       &
     &           ip_keu=5, ip_kev=6, ip_kew=7, ip_q=8, ip_qcl=9,        &
     &           ip_qcf=10)
      Parameter (ip_qu=1, ip_qv=2, ip_qw=3)

      Real                                                              &
     &  weight1, weight2, weight3                                       &
                                  ! weights
     &, tempd, tempw, ww2, ww1                                          &
     &, cv                        ! specific heat at constant vol

      Real                                                              &
     &  rho(row_length,rows,model_levels)                               &
                                                ! rho only
     &, delr_rho(row_length,rows,model_levels)                          &
                                                ! dr
     &, T(row_length,rows,model_levels)                                 &
                                                ! TEMPERATURE
     &, rho_theta(row_length,rows,model_levels)                         &
                                                 ! not used
     &, t_rho(row_length,rows,model_levels)                             &
     &, u_rho(row_length,rows,model_levels)                             &
     &, v_rho(row_length,rows,model_levels)                             &
     &, w_rho(row_length,rows,model_levels)                             &
     &, q_rho(row_length,rows,wet_model_levels+1)                       &
     &, qcl_rho(row_length,rows,wet_model_levels+1)                     &
     &, qcf_rho(row_length,rows,wet_model_levels+1)

      Real                                                              &
     &  mag_vector_np(model_levels)                                     &
                                     ! magnitude of the vector wind NP
     &, dir_vector_np(model_levels)                                     &
                                     ! direction of the vector wind NP
     &, mag_vector_sp(model_levels)                                     &
                                     ! magnitude of the vector wind SP
     &, dir_vector_sp(model_levels)  ! direction of the vector wind SP

      INTEGER I, J, K                ! LOOP COUNTER

!----------------------------------------------------------------------
! EXTERNAL SUBROUTINE CALLS  -
!----------------------------------------------------------------------
      External                                                          &
     &  u_to_p,v_to_p, polar_vector_wind_n,swap_bounds
!----------------------------------------------------------------------
! constant need later

      cv=cp-R       ! value for dry air

!----------------------------------------------------------------------
! zero output arrays
!----------------------------------------------------------------------
      Do k=1,10         ! energy integrals
        Do j=1,rows
          Do i=1,row_length
            vert_int_array(i,j,k) = 0.0
          End Do
        End Do
      End Do
      IF (lqflux) then
        Do k=1,3          ! q flux integrals
          Do j=1,rows
            Do i=1,row_length
              vert_qflux(i,j,k) = 0.0
            End Do
          End Do
        End Do
      Endif

!----------------------------------------------------------------------
! First convert theta to temperature
!----------------------------------------------------------------------
      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            T(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
          End Do
        End Do
      End Do

!----------------------------------------------------------------------
! CALCULATE layer thickness for rho layers and calculate rho
!----------------------------------------------------------------------
      DO K=1,model_levels
        Do j = 1, rows
          DO I=1,row_length

            DELR_RHO(I,J,K) = R_THETA_LEVELS(I,J,K) -                   &
     &                        R_THETA_LEVELS(I,J,K-1)
            rho(i,j,k) = rho_r2(i,j,k)/                                 &
     &           (r_rho_levels(i,j,k)*r_rho_levels(i,j,k))

          END DO
        END DO
      END DO

!----------------------------------------------------------------------
! convert rho to rho dry in the same way as done in dynamics - Flux_rho
! Note now uses linear vertical intepolation to be consistent with
! new dynamics code at UM 5.1.
!----------------------------------------------------------------------
      k = 1
        Do j = 1, rows
          Do i = 1, row_length
            tempd = 1. - q(i,j,k)- qcl(i,j,k) - qcf(i,j,k)
            rho_dry(i,j,k) = rho_r2(i,j,k) * tempd
            dry_to_wet(i,j,k)= 1. / tempd
          End Do
        End Do

      Do k = 2, wet_model_levels
        Do j = 1, rows
          Do i = 1, row_length
! linear interpolation weights
            weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
            weight2 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
            weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)

            tempd = ( weight2 *                                         &
     &              (1. - q(i,j,k)- qcl(i,j,k) - qcf(i,j,k) ) +         &
     &                weight1 *                                         &
     &              (1. - q(i,j,k-1)- qcl(i,j,k-1) - qcf(i,j,k-1) ) )   &
     &               / weight3

            dry_to_wet(i,j,k)= 1./tempd
            rho_dry(i,j,k) = rho_r2(i,j,k) * tempd
          End Do
        End Do
      End Do

! Special case of wet_model_levels < model_levels

      If (wet_model_levels <  model_levels) then
        k=wet_model_levels+1
        Do j = 1, rows
          Do i = 1, row_length
            weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
            weight2 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
            weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
            tempd = ( weight2 +                                         &
     &                   weight1 *                                      &
     &      (1. - q(i,j,k-1)- qcl(i,j,k-1) - qcf(i,j,k-1)) ) / weight3

            dry_to_wet(i,j,k)= 1./tempd
            rho_dry(i,j,k) = rho_r2(i,j,k) * tempd
          End Do
        End Do

        Do k = wet_model_levels+2,model_levels
          Do j = 1, rows
            Do i = 1, row_length
              rho_dry(i,j,k) = rho_r2(i,j,k)
              dry_to_wet(i,j,k)= 1.
            End Do
          End Do
        End Do
      Endif      ! test on wet levels < model_levels

!-------------------------------------------------------------------
! Intepolate T, w & q to rho points
! Note needs to remain consistent with intepolation used in dynamics.
! Using linear interpolation
!
!                      K               for rho
!      K                          K-1  for theta
!      X<--- w1------->X<-- w2--->X
!       <----------w3------------>
!-------------------------------------------------------------------
      k=1
        Do j = 1, rows
          DO I=1,row_length
! assume bottom rho level value equal to bottom theta level value
            q_rho(i,j,k)   = q(i,j,k)
            qcl_rho(i,j,k) = qcl(i,j,k)
            qcf_rho(i,j,k) = qcf(i,j,k)
            T_rho(i,j,k)   = T(i,j,k)
! only w has a value at surface
            weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
            weight2 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
            weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
            ww1 = weight1/weight3
            ww2 = weight2/weight3
            w_rho(i,j, K) =  ww2 * w(i,j,K)  + ww1 * w(i,j,K-1)
          END DO
        END DO

      DO K=2,wet_model_levels
        Do j = 1, rows
          DO I=1,row_length
            weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
            weight2 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
            weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
            ww1 = weight1/weight3
            ww2 = weight2/weight3
            t_rho(i,j, K) = ww2 * t(i,j,K) + ww1 * t(i,j,K-1)
            w_rho(i,j, K) = ww2 * w(i,j,K) + ww1 * w(i,j,K-1)
            q_rho (i,j, K)  = ww2 * q(i,j,K)   + ww1 * q(i,j,K-1)
            qcl_rho(i,j, K) = ww2 * qcl(i,j,K) + ww1 * qcl(i,j,K-1)
            qcf_rho(i,j, K) = ww2 * qcf(i,j,K) + ww1 * qcf(i,j,K-1)
          END DO
        END DO
      END DO

! Special case of wet_model_levels < model_levels
! wet to dry conversions imply q_rho has a value on wet_model_levels+1
! Assume also true for qcl and qcf ?

      if (wet_model_levels <  model_levels) then
        k=wet_model_levels+1
        Do j = 1, rows
          DO I=1,row_length
            weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
            weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
            ww1 = weight1/weight3
            q_rho (i,j, K)  = ww1 * q(i,j,K-1)
            qcl_rho(i,j, K) = ww1 * qcl(i,j,K-1)
            qcf_rho(i,j, K) = ww1 * qcf(i,j,K-1)
          END DO
        END DO
      endif

! Any dry levels at top of model
      DO K=wet_model_levels+1,model_levels
        Do j = 1, rows
          DO I=1,row_length
            weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
            weight2 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
            weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
            ww1 = weight1/weight3
            ww2 = weight2/weight3
            t_rho(i,j, K) = ww2 * t(i,j,K) + ww1 * t(i,j,K-1)
            w_rho(i,j, K) = ww2 * w(i,j,K) + ww1 * w(i,j,K-1)
          END DO
        END DO
      END DO
!-------------------------------------------------------------------
! Interpolate winds to rho points (already on same vertical level)
!-------------------------------------------------------------------
! Need to call swap bounds as halo points not set
! DEPENDS ON: swap_bounds
      CALL swap_bounds(u,row_length,rows,model_levels,                  &
     &                       off_x,off_y,fld_type_u,.true.)
! DEPENDS ON: swap_bounds
      CALL swap_bounds(v,row_length,n_rows,model_levels,                &
     &                       off_x,off_y,fld_type_v,.true.)

! DEPENDS ON: v_to_p
      CALL v_to_p(v,row_length,rows,n_rows,model_levels,                &
     &                  off_x,off_y,model_domain,at_extremity,v_rho)

! DEPENDS ON: u_to_p
      CALL u_to_p(u,row_length,rows,model_levels,                       &
     &                  off_x,off_y,model_domain,at_extremity,u_rho)

! Problem of u & v values at poles
! Uses v on row next to poles to calculate wind at the polar point

        If (model_domain  ==  1 ) Then    ! global model

! DEPENDS ON: polar_vector_wind_n
          Call Polar_vector_wind_n                                      &
     &                      (v,                                         &
     &                       sin_theta_longitude,cos_theta_longitude,   &
     &                       row_length,n_rows, model_levels,           &
     &                       mag_vector_np,dir_vector_np,               &
     &                       mag_vector_sp,dir_vector_sp,               &
     &                       off_x, off_y, global_row_length,           &
     &                       proc_row_group, at_extremity)

          If (at_extremity(PSouth) ) Then
            Do k = 1,model_levels
              Do i = 1, row_length
                u_rho(i,1,k) = 0.0
                v_rho(i,1,k) = mag_vector_sp(k)
              End Do
            End Do
          End If
         If (at_extremity(PNorth) ) Then
            Do k = 1,model_levels
              Do i = 1, row_length
                u_rho(i,rows,k) = 0.0
                v_rho(i,rows,k) = mag_vector_np(k)
              End Do
            End Do
          Endif
        Endif              ! end test on global model
!----------------------------------------------------------------------
! Integrals over fields using values interpolated to rho grid
! At present all integrals done after interpolation to rho points.
! Integrals for energy involve using dry mass. This may need to be
! altered.
!----------------------------------------------------------------------
! (a) fields on all model levels

      DO K=1,model_levels
        DO j = 1, rows
          DO I=1,row_length

            tempw = RHO_r2(I,j,k)*DELR_RHO(I,j,k)   ! wet mass

            tempd = RHO_dry(I,j,k)*DELR_RHO(I,j,k)  ! dry mass

! dry mass
            vert_int_array(I,J,ip_dry_mass) =                           &
     &                    vert_int_array(i,j,ip_dry_mass) + tempd

! total mass wet + dry
            vert_int_array(I,J,ip_wet_mass) =                           &
     &                    vert_int_array(i,j,ip_wet_mass) + tempw

! gr term
            vert_int_array(I,J,ip_gr) = vert_int_array(i,j,ip_gr) +     &
     &         g*(r_rho_levels(i,j,k)-earth_radius)*tempd

! KE terms for u, v & w
            vert_int_array(I,J,ip_keu) = vert_int_array(i,j,ip_keu) +   &
     &                 0.5*u_rho(i,j,k)*u_rho(i,j,k)*tempd
            vert_int_array(I,J,ip_kev) = vert_int_array(i,j,ip_kev) +   &
     &                 0.5*v_rho(i,j,k)*v_rho(i,j,k)*tempd
            vert_int_array(I,J,ip_kew) = vert_int_array(i,j,ip_kew) +   &
     &                 0.5*w_rho(i,j,k)*w_rho(i,j,k)*tempd

! cvT term
            vert_int_array(I,J,ip_cvt) = vert_int_array(i,j,ip_cvt) +   &
     &                            cv*T_rho(i,j,k)*tempd
          END DO
        END DO
      END DO

! (b) fields on wet levels only

      DO K=1,wet_model_levels
        DO j = 1, rows
          DO I=1,row_length
            tempw = RHO_r2(I,j,k)*DELR_RHO(I,j,k)  ! wet mass
! q*rho
            vert_int_array(I,J,ip_q) = vert_int_array(i,j,ip_q) +       &
     &                                      q_rho(i,j,k)*tempw
! qcl*rho
            vert_int_array(I,J,ip_qcl) = vert_int_array(i,j,ip_qcl) +   &
     &                                      qcl_rho(i,j,k)*tempw
! qcf*rho
            vert_int_array(I,J,ip_qcf) = vert_int_array(i,j,ip_qcf) +   &
     &                                      qcf_rho(i,j,k)*tempw

          END DO
        END DO
      END DO

! Special case of wet_model_levels < model_levels
      If (wet_model_levels <  model_levels) then
        k=wet_model_levels+1
        DO j = 1, rows
          DO I=1,row_length
            tempw = RHO_r2(I,j,k)*DELR_RHO(I,j,k)  ! wet mass
            vert_int_array(I,J,ip_q) = vert_int_array(i,j,ip_q) +       &
     &                                      q_rho(i,j,k)*tempw
            vert_int_array(I,J,ip_qcl) = vert_int_array(i,j,ip_qcl) +   &
     &                                      qcl_rho(i,j,k)*tempw
            vert_int_array(I,J,ip_qcf) = vert_int_array(i,j,ip_qcf) +   &
     &                                      qcf_rho(i,j,k)*tempw
          Enddo
        Enddo
      Endif   !         wet_model_levels < model_levels

!----------------------------------------------------------------------
! Additional q flux integrals not required by energy correction
! At present all integrals done after interpolation to rho points.
!----------------------------------------------------------------------
! (b) fields on wet levels only

      If (Lqflux) then
        DO K = 1,wet_model_levels
          Do j = 1, rows
            DO I = 1,row_length
              tempw = RHO_r2(I,j,k)*DELR_RHO(I,j,k)
              vert_qflux(I,J,ip_qu) = vert_qflux(i,j,ip_qu) +           &
     &                               u_rho(i,j,k)*q_rho(i,j,k)*tempw
              vert_qflux(I,J,ip_qv) = vert_qflux(i,j,ip_qv) +           &
     &                               v_rho(i,j,k)*q_rho(i,j,k)*tempw
              vert_qflux(I,J,ip_qw) = vert_qflux(i,j,ip_qw) +           &
     &                               w_rho(i,j,k)*q_rho(i,j,k)*tempw

            END DO
          END DO
        END DO

! Special case of wet_model_levels < model_levels

        If (wet_model_levels <  model_levels) then
          k=wet_model_levels+1
          DO j = 1, rows
            DO I=1,row_length
              vert_qflux(I,J,ip_qu) = vert_qflux(i,j,ip_qu) +           &
     &                               u_rho(i,j,k)*q_rho(i,j,k)*tempw
              vert_qflux(I,J,ip_qv) = vert_qflux(i,j,ip_qv) +           &
     &                               v_rho(i,j,k)*q_rho(i,j,k)*tempw
              vert_qflux(I,J,ip_qw) = vert_qflux(i,j,ip_qw) +           &
     &                               w_rho(i,j,k)*q_rho(i,j,k)*tempw
            Enddo
          Enddo
        Endif  !           wet_model_levels < model_levels
      Endif
!----------------------------------------------------------------------

      RETURN
      END SUBROUTINE VERT_ENG_MASSQ
#endif
