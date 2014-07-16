
! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
SUBROUTINE DEEP_TURB_GRAD_STRESS(n_dp, nlev, max_cldlev                &
     ,                       nclev                                     &
     ,                       timestep                                  &
     ,                       zcld                                      &
     ,                       mf_cld, w_up, k_func                      &
     ,                       u,v                                       &
     ,                       r2rho,r2rho_th,dr_across_rh,dr_across_th  &
        ! output arguements
     ,                       uw_cld,vw_cld)
! Description:
!   This routine calculates the in-cloud gradient part of the deep CMT stress
!   using turbulence ideas. This version is designed for use with the mass flux
!   convection scheme. Note at this stage this version is different from the 
!   shallow turbulence version in its assumptions at cloud top.
!
!  Solves the equation
!      du/dt  = -d(uw)/dz for the gradient part of uw.
!
!  The above equation is solved implicitly
!
!  The gradient part of uw = - mf(wup/wcld)*K(eta)du/dz        
!
! ------------------------------------------------------------------------------

  IMPLICIT NONE

! ------------------------------------------------------------------------------
! Arguments with intent IN:
!

    Integer, intent(in) :: &
      nlev                 & ! No. of model layers
    , max_cldlev           & ! Maximum number of cloud levels
    , n_dp                 & ! No. of deep convection points
    , nclev(n_dp)            ! number of cloud levels

    Real, intent(in) ::  &
      timestep                ! model timestep (s)  

    Real, intent(in) ::       &
      zcld(n_dp)              & ! cloud layer depth (m) for CMT
    , w_up(n_dp,nlev)         & ! wup/wcld 
    , mf_cld(n_dp,nlev)       & ! mass flux in (th lev) (m/s)
    , k_func(n_dp,nlev)       & ! k function (dependent on eta)
    , u(n_dp,nlev)            & ! U-component of mean wind (m/s)
    , v(n_dp,nlev)            & ! V-component of mean wind (m/s)
    , r2rho(n_dp,nlev)        & ! r2*rho on rho levels (kg/m)
    , r2rho_th(n_dp,nlev)     & ! r2*rho on theta levels (kg/m) 
    , dr_across_th(n_dp,nlev) & ! thickness of theta layers (m)
    , dr_across_rh(n_dp,nlev)   ! thickness of rho layers (m)

!
! Arguments with intent out:
!
    Real, intent(out) ::    &
      uw_cld(n_dp,nlev)     & ! U-component of stress (m2/s2)
    , vw_cld(n_dp,nlev)       ! V-component of stress (m2/s2)

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
!
    Integer :: i,k    ! counters

    Real ::                       &
      visc(n_dp,max_cldlev+1)     & ! viscosity profile (M2S-1)
    , A(n_dp,max_cldlev+1)        & ! tridiagonal matrix elements
    , B(n_dp,max_cldlev+1)        &
    , C(n_dp,max_cldlev+1)        &
    , ue_tp1(n_dp,max_cldlev+1)   & ! after timestep velocity vectors
    , ve_tp1(n_dp,max_cldlev+1)   &
    , dz,dz12                 

    Integer ::           &
      max_cldlev1          ! max cloud levels plus 1
!
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
! Calculate the eddy viscosity profile
!
!        K'= K(z/zcld)*mf*zcld*wup/wcld 
!
!   on theta levels - ie uw stress levels
! Note no specail conditions for values at cloud top unlike shallow
! code.
!-----------------------------------------------------------------------

      max_cldlev1 = max_cldlev+1

!CDIR NOUNROLL
      Do k=1,max_cldlev1  ! from level above cloud base 
        Do i=1,n_dp

          If (k <=  nclev(i)) then         ! in cloud

            visc(i,k)=k_func(i,k)*mf_cld(i,k)*W_up(i,k)*zcld(i) 

          Else
            visc(i,k)=0.0
          Endif
        End Do
      End Do

!
! Calculate GRADIENT component of stress. Use implicit timestepping
!

      k=1
        Do i=1,n_dp
          dz   = dr_across_rh(i,k)*r2rho(i,k)
          dz12 = dr_across_th(i,k)
          A(i,k) = 0.0
          C(i,k) = -visc(i,k)*r2rho_th(i,k)*timestep/(dz*dz12)  
          B(i,k) = 1.0 - A(i,k) - C(i,k)
        End do

!CDIR NOUNROLL
      Do k=2,max_cldlev1
        Do i=1,n_dp
          dz = dr_across_rh(i,k)*r2rho(i,k)

          If (k < nclev(i)) Then

            dz12 = dr_across_th(i,k-1)
            A(i,k) = -visc(i,k-1)*r2rho_th(i,k-1)*timestep/(dz*dz12) 

            dz12 = dr_across_th(i,k)
            C(i,k) = -visc(i,k)*r2rho_th(i,k)*timestep/(dz*dz12)

          Else if (k == nclev(i)) Then

            dz12 = dr_across_th(i,k-1)
            A(i,k) = -visc(i,k-1)*r2rho_th(i,k-1)*timestep/(dz*dz12)
            C(i,k) = 0.0

          Else    ! elements not required in calculation (zero)

            C(i,k) = 0.0
            A(i,k) = 0.0

          End if

          B(i,k) = 1.0 - A(i,k) - C(i,k)

        End do
      End do

!
! Calculate NEW timestep wind conponents using tridiagonal matrix solver
!
! DEPENDS ON: tridiag_all
       CALL TRIDIAG_all(max_cldlev1,n_dp,nclev,A,B,C,u,ue_tp1)
! DEPENDS ON: tridiag_all
       CALL TRIDIAG_all(max_cldlev1,n_dp,nclev,A,B,C,v,ve_tp1)


!
! Calculate stress profile -Kdu/dz from latest u/v values
!

!CDIR NOUNROLL
      Do k=1,max_cldlev1
        Do i=1,n_dp

          If (k < nclev(i)) then   ! in cloud

            dz = dr_across_th(i,k)
            uw_cld(i,k)=-visc(i,k)*(ue_tp1(i,k+1)-ue_tp1(i,k))/dz
            vw_cld(i,k)=-visc(i,k)*(ve_tp1(i,k+1)-ve_tp1(i,k))/dz

          Else        ! not in cloud  set to zero

            uw_cld(i,k) = 0.0
            vw_cld(i,k) = 0.0

          End if

        End Do
      End Do

!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE DEEP_TURB_GRAD_STRESS
