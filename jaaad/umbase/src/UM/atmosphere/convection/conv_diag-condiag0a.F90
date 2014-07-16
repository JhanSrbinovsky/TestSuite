#if defined(A05_0A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Diagnosis of convective occurance and type
!
      SUBROUTINE CONV_DIAG(                                             &

! IN values defining field dimensions and subset to be processed :
     & row_length, rows, n_rows, halo_i, halo_j, off_x, off_y           &
     &, global_row_length, proc_row_group, at_extremity                 &

! IN values defining vertical grid of model atmosphere :
     &, model_domain, bl_levels, model_levels, wet_model_levels         &
     &, p, P_theta_lev,exner_rho, r_rho_levels,r_theta_levels           &
     &, rho_only, rho_theta, z_full, z_half                             &
     &, SIN_THETA_LONGITUDE,COS_THETA_LONGITUDE                         &

! IN Model switches
     &, l_mixing_ratio                                                  &

! IN Cloud data :
     &, qcf, qcl, cloud_fraction                                        &

! IN everything not covered so far :
     &, pstar, q, theta, exner_theta_levels, u, v, u_0_p, v_0_p         &
     &, L_flux_bc, flux_e, flux_h, L_spec_z0, z0h_scm                   &
     &, tstar, land_mask, timestep                                      &
     &, w_copy                                                          &

! SCM Diagnostics (dummy values in full UM)
     &, nSCMDpkgs,L_SCMDiags                                            &

! OUT data required elsewhere in UM system :
     &, zh,zhpar,z_lcl,z_lcl_uv,delthvu,ql_ad,NTML,ntpar,NLCL,cumulus   &
     &, L_shallow,l_congestus, l_congestus2, CIN, CAPE                  &

     &, ERROR                                                           &
     & )

! Purpose:
!   Diagnosis of convective occurance - version for use with no
!     convection scheme.
!
!   Called by Atmos_physics2
!
! Current owners of code: Convection code owner
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
!
      IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!
!
! (a) Defining horizontal grid and subset thereof to be processed.

      integer, intent(in) ::                                            &
     &  row_length                                                      &
                           ! Local number of points on a row
     &, rows                                                            &
                           ! Local number of rows in a theta field
     &, n_rows                                                          &
                           ! Local number of rows in a v field
     &, halo_i                                                          &
                           ! Size of halo in i direction.
     &, halo_j                                                          &
                           ! Size of halo in j direction.
     &, off_x                                                           &
                           ! Size of small halo in i
     &, off_y                                                           &
                           ! Size of small halo in j.
     &, global_row_length                                               &
                           ! number of points on a row
     &, proc_row_group     ! Group id for processors on the same row

      Logical,  intent(in) ::                                           &
     &  at_extremity(4)    ! Indicates If this processor is at north,
                           !  south,east or west of the processor grid

! (b) Defining vertical grid of model atmosphere.

      integer, intent(in) ::                                            &
     &  model_domain                                                    &
     &, bl_levels                                                       &
                           ! Max. no. of "boundary" levels allowed.
     &, model_levels                                                    &
                           ! number of model levels
     &, wet_model_levels
                           ! number of wet model levels

      real, intent(in) ::                                               &
     &  P(1-off_x:row_length+off_x,1-off_y:rows+off_y, model_levels)    &
                   ! pressure on rho levels  (Pa)
     &, P_theta_lev(row_length, rows, model_levels)                     &
                   ! pressure on theta levels  (Pa)
     &, exner_rho(1-off_x:row_length+off_x,1-off_y:rows+off_y,          &
                                              model_levels)             &
                   ! Exner pressure on rho levels 
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j,model_levels)                 &
                   ! radius of rho levels (m)
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:model_levels)             &
                   ! radius of theta levels (m)
     &, rho_only(row_length,rows,model_levels)                          &
                                               ! density (kg/m3)
     &, rho_theta(row_length,rows,model_levels-1)                       &
                                                  ! rho th lev (kg/m3)
     &, z_full(row_length,rows,model_levels)                            &
                                               ! height th lev (m)
     &, z_half(row_length,rows,model_levels)                            &
                                               ! height rho lev (m)
     &, sin_theta_longitude (row_length,rows)                           &
                                                 ! sin(long)
     &, cos_theta_longitude (row_length,rows)    ! cos(long)


      Logical,  intent(in) ::                                           &
     &  L_mixing_ratio                                                  &
                           ! true if input moisture fields mixing ratio
                           ! false input moisture fields specific
     &, L_flux_bc                                                       &
                           ! true if using specified surface fluxes
     &, L_spec_z0          ! true if using specified roughness lengths

! (c) Cloud data.

      real, intent(in) ::                                               &
     & QCF(row_length,rows,wet_model_levels)                            &
!                                 ! Cloud ice (kg per kg air)
     &,QCL(row_length,rows,wet_model_levels)                            &
!                                 ! Cloud liquid water (kg/kg air).
     &,cloud_fraction(row_length, rows, wet_model_levels)
!                                 ! Cloud fraction

! (d) Atmospheric + any other data not covered so far, incl control.

      real, intent(in) ::                                               &
     &  PSTAR(row_length, rows)                                         &
                                       ! Surface pressure (Pascals).
     &, Q(row_length,rows,wet_model_levels)                             &
                                             ! water vapour (kg/kg)
     &, theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &          model_levels)                                           &
                                       ! Theta (Kelvin)
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)            &
                                       ! exner pressure theta lev (Pa)
     &, U(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               & 
                                       !  W'ly wind component (m/s).
     &, V(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)                                               & 
                                       !  S'ly wind component (m/s).
     &, U_0_P(row_length,rows)                                          & 
                                  ! W'ly component of surface current
                                  !    (metres per second) on P-grid.
     &, V_0_P(row_length,rows)                                          & 
                                  ! S'ly component of surface current
                                  !    (metres per second) on P-grid.
     &,flux_e(row_length,rows)                                          &
                                  !  Specified surface
                                  !    latent heat flux (W/m^2)
     &,flux_h(row_length,rows)                                          &
                                  !  Specified surface
                                  !    sensible heat fluxes (in W/m2)
     &,z0h_scm(row_length,rows)   !  Namelist input z0h (if >0)
                                  !    (if <=0 use Z0HSEA)

      real, intent(inout) ::                                            &
     & TSTAR(row_length,rows)    ! Surface temperature (= top soil
                                 !  layer temperature) (K).
                                 ! Note only INOUT if SCM
      logical, intent(in) ::                                            &
     & Land_MASK(row_length, rows)  ! T If land, F Elsewhere.


      real, intent(in) ::                                               &
     &  TIMESTEP                                                        &
                                     ! Timestep (seconds).
     &, w_copy(row_length,rows,                                         &
     &          0:model_levels)      ! vertical velocity (m/s)


! Additional variables for SCM diagnostics which are dummy in full UM
      integer                                                           &
     &  nSCMDpkgs          ! No of diagnostics packages
!
      logical                                                           &
     &  L_SCMDiags(nSCMDpkgs) ! Logicals for diagnostics packages


      real, intent(inout) ::                                            &
     & zh(row_length,rows)    ! Height above surface of top
                              !  of boundary layer (metres).
      real, intent(out) ::                                              &
     & zhpar(row_length,rows)                                           &
                              ! Height of max parcel ascent (m)
     &,z_lcl(row_length,rows)                                           &
                              ! Height of lIfting condensation level (m)
     &,z_lcl_uv(row_length,rows)                                        &
                                 ! Height of lIfting condensation
                                 !     level on uv grid (m)
                                 ! ONLY REQUIRED BY CONVECTION SCHEME
     &,delthvu(row_length,rows)                                         &
                                 ! Integral of undilute parcel buoyancy
                                 ! over convective cloud layer
                                 ! (for convection scheme)
     &,ql_ad(row_length,rows)                                           &
                                 ! Adiabatic liquid water content at 
                                 ! inversion (kg/kg) (5A code only)
     &,CAPE(row_length,rows)                                            &
                                 ! CAPE from parcel ascent (m2/s2)
     &,CIN(row_length,rows)      ! CIN from parcel ascent (m2/s2)

      integer, intent(out)  ::                                          &
     & NTML(row_length,rows)                                            &
                               ! Number of model levels in the
                               ! turbulently mixed layer.
     &,ntpar(row_length,rows)                                           &
                               ! Max levels for parcel ascent
     &,nlcl(row_length,rows)   ! No. of model layers below the
                               ! lIfting condensation level.


      logical, intent(out) ::                                           &
     & cumulus(row_length,rows)                                         &
                                  ! Logical indicator for convection
     &,L_shallow(row_length,rows)                                       &
                                  ! Logical indicator for shallow Cu
     &,L_congestus(row_length,rows)                                     &
                                    ! Logical indicator for congestus Cu
     &,L_congestus2(row_length,rows) ! Logical ind 2 for congestus Cu

      integer, intent(out)  ::                                          &
     & ERROR          ! 0 - no error in this routine

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

      Character*(*), Parameter ::  RoutineName = 'conv_diag'

      integer ::                                                        &
     & I,J                                                              &
                  ! LOCAL Loop counter (horizontal field index).
     &,K                                                                &
                  ! LOCAL Loop counter (vertical level index).
     &,MBL        ! Maximum number of model layers allowed in the
!                 ! mixing layer; set to bl_levels-1.

      real ::                                                           &
     &  MAG_VECTOR_NP(1)                                                &
     &, DIR_VECTOR_NP(1)                                                &
     &, MAG_VECTOR_SP(1)                                                &
     &, DIR_VECTOR_SP(1)

      logical ::                                                        &
     & topbl(row_length, rows)                                          &
                                 ! Flag set when top of boundary layer
!                                ! is reached.
     &,above_lcl(row_length,rows)                                       &
                                 ! Flag set when parcel above LCL.
     &,UNSTABLE(row_length, rows)  ! Flag to indicate an unstable
!                                  ! boundary layer forced from
!                                  ! the surface.

      real ::                                                           &
     & T(row_length, rows, model_levels)                                &
                                         ! temperature (from theta)
     &,QS(row_length, rows, wet_model_levels)                           &
!                              ! Saturated sp humidity at pressure and
!                              ! temperature of sucessive levels.
     &,QS_STAR(row_length, rows)                                        &
                                 ! Saturated sp humidity at sea surface
     &,TL(row_length, rows, model_levels)                               &
!                              ! Ice/liquid water temperature,
!                                 but replaced by T in LS_CLD.
     &,QW(row_length, rows, wet_model_levels)                           &
                                               ! Total water content
     &,SVL(row_length, rows, wet_model_levels)                          &
!                              ! Liquid/frozen water virtual
!                              ! static energy over CP.
     &,T_LCL(row_length, rows)                                          &
                                      ! Temperature at lIftng condensati
!                                     ! level.
     &,P_LCL(row_length, rows)                                          &
                                      ! Pressure at lIfting condensation
!                                     ! level.
     &,sl_plume(row_length, rows)                                       &
                                      ! Liquid/frozen water static
!                                     ! energy over CP for a plume
!                                     ! rising without dilution
!                                     ! from level 1.
     &,QW_PLUME(row_length, rows)                                       &
                                      ! QW for a plume rising without
!                                     ! dilution from level 1.
     &,th_ref(row_length, rows)                                         &
                                      ! theta - reference value
     &,t_ref(row_length, rows)                                          &
                                      ! t - reference value
     &,qsat_lev(row_length, rows)                                       &
                                      ! qsat for refreence temperature
!
     &,env_svl(row_length, rows,wet_model_levels)                       &
                                      ! Density (virtual) static energy
!                                     ! over CP for layer.
     &,par_svl(row_length, rows,wet_model_levels)                       &
                                      ! Density (virtual) static energy
                                      ! over CP of parcel for level.
     &,thv_pert(row_length, rows)                                       &
                                      ! threshold thv of parcel
     &,dtv_min(row_length, rows)                                        &
                                      ! min TV of parcel in cld layer
     &,TV1_SD(row_length, rows)                                         &
                                      ! Approx to standard dev of level
!                                     ! 1 virtual temperature (K).
     &,fb_surf(row_length, rows)                                        &
                                      ! Change in theta_v from surface
!                                     ! to layer 1 (note dIff from BL)
     &,U_P(row_length, rows)                                            &
                                      ! U(1) on P-grid.
     &,V_P(row_length, rows)                                            &
                                      ! V(1) on P-grid.
     &,DQSDT(row_length,rows)                                           &
                                      ! d(QSAT)/dT
     &,Z0H(row_length,rows)           ! local roughness length (m)

      integer ::                                                        &
     & kcucheck(row_length, rows)                                       &
                                  ! Position of level just below 2.5km
!                                ! (used for gradient check to
!                                !  diagnose convection)
     &,k_plume(row_length, rows) ! start level for surface-driven plume


!
! parcel calculation
!
      real ::                                                           &
     & q_vap_parc                                                       &
                         ! Vapour content of parcel.
     &,q_liq_parc                                                       &
                         ! Condensed water content of parcel.
     &,T_PARc                                                           &
                         ! Temperature of parcel.
     &,t_dens_parc(row_length, rows, wet_model_levels)                  &
                         ! Density potential temperature of parcel.
     &,t_dens_env(row_length, rows, wet_model_levels)                   &
                         ! Density potential temperature of environment.
     &,denv_bydz(row_length, rows, wet_model_levels)                    &
                         ! Gradient of density potential
                         ! temperature in the environment.
     &,dpar_bydz(row_length, rows, wet_model_levels)                    &
                         ! Gradient of density potential
                         ! temperature of the parcel.
     &,th_par_km1(row_length, rows)                                     &
                                     ! parcel theta at level below
     &,Z_SURF                                                           &
                         ! approx height of top of surface layer
     &,vap_press                                                        &
                         ! Vapour pressure.
     &,grad_cld                                                         &
                         ! SVL gradient in layer above LCL.
     &,grad_sub                                                         &
                         ! SVL gradient in layer below LCL.
     &,q_liq_env                                                        &
                         ! Condensed water content of environment.
     &,dq_sat_env                                                       &
                         ! DQSAT/DT for environment
     &,theta1                                                           &
                         ! Potential temperature in layer 1
     &,USHEAR                                                           &
                         ! U wind shear from level 1 to surface
     &,VSHEAR                                                           &
                         ! V wind shear from level 1 to surface
     &,WSHR1                                                            &
                         ! magnitude of surface wind shear
     &,WSHR2                                                            &
                         ! (WSHR1)**2
     &,rhostar                                                          &
                         ! surface air density
     &,theta_star                                                       &
                         ! theta at surface
     &,wthvbar                                                          &
                         ! surface buoyancy flux
     &,CD                                                               &
                         ! bulk transfer coefficient for momentum
!                        ! (Approximation. Also assume CH(heat)=CD.)
     &,ustar                                                            &
                         ! surface friction velocity
     &, Lrcp_const                                                      &
                          ! lc or lc+lf over cp
     &, Lrcp_const_parc                                                 &
                               ! lc or lc+lf over cp
     &, L_const                                                         &
                         ! lc or lc+lf
     &, dz                                                              &
                         ! layer depth
     &, z_pr, th_par                                                    &
     &, dmass, inc

      real ::                                                           &
     &  r2rho_theta(row_length,rows,model_levels)                       &
                                                   ! r**2rho on theta
     &, w_avg(row_length,rows)                                          &
                                      ! mean w over layer (m/s)
     &, mass(row_length,rows)         ! mass for layer (kg)
!
! Model constants:
!

#include "c_lheat.h"
#include "c_r_cp.h"
#include "c_g.h"
#include "c_a.h"
#include "c_0_dg_c.h"
#include "c_epslon.h"
#include "c_rough.h"
#include "c_vkman.h"
#include "c_bolton1980.h"
#include "fldtype.h"
#include "domtyp.h"
#if defined(SCMA)
! Include parameters necessary for calls to SCMoutput...
#include "s_scmop.h"
#endif


      real, parameter ::                                                &
     & GRCP=G/CP                                                        &
                                ! Adiabatic lapse rate.
     &,LCRCP=LC/CP                                                      &
                                ! Latent heat of condensation / CP.
     &,LFRCP=LF/CP                                                      &
                                ! Latent heat of fusion / CP.
     &,LS=LC+LF                                                         &
                                ! Latent heat of sublimation.
     &,LSRCP=LS/CP                                                      &
                                ! Latent heat of sublimation / CP.
     &,ra2=1./(Earth_Radius*Earth_Radius)   ! 1./(a*a)

      real, parameter ::                                                &
     & A_PLUME = 0.2                                                    &
     &,B_PLUME = 3.26                                                   &
     &,MAX_T_GRAD = 1.0E-3                                              &
     &,SC_CFTOL   = 0.1         ! Cloud fraction required for a
                                ! cloud layer to be diagnosed.
!
! Parameters
!
      Integer, parameter ::                                             &
     &  PNorth = 1                                                      &
                     ! North processor address in the neighbor array
     &, PEast  = 2                                                      &
                     ! East processor address in the neighbor array
     &, PSouth = 3                                                      &
                     ! South processor address in the neighbor array
     &, PWest  = 4                                                      &
                     ! West processor address in the neighbor array
     &, NoDomain = -1 ! Value in neighbor array If the domain has
                      !  no neighbor in this direction. Otherwise
                      !  the value will be the tid of the neighbor



!-----------------------------------------------------------------------
! Mixing ratio, r,  versus specific humidity, q
!
! In most cases the expression to first order are the same
!
!  Tl = T - (lc/cp)qcl - [(lc+lf)/cp]qcf
!  Tl = T - (lc/cp)rcl - [(lc+lf)/cp]rcf  - equally correct definition
!
! thetav = theta(1+cvq)         accurate
!        = theta(1+r/epsilon)/(1+r) ~ theta(1+cvr) approximate
! 
! svl = (Tl+gz/cp)*(1+(1/epsilon-1)qt)   
!     ~ (Tl+gz/cp)*(1+(1/epsilon-1)rt)
!
! dqsat/dT = epsilon*Lc*qsat/(R*T*T) 
! drsat/dT = epsilon*Lc*rsat/(R*T*T)  equally approximate
!
! Only altering the expression for vapour pressure
!
!  e = qp/epsilon       - approximation 
!  e = rp/(epsilon+r)   - accurate    
!-----------------------------------------------------------------------
! 1.0 Initialisation
!-----------------------------------------------------------------------

      ERROR = 0

!-----------------------------------------------------------------------
! 1.1 Verify grid/subset definitions.
!-----------------------------------------------------------------------

      If ( bl_levels <  1 .OR. ROWS <  1 .OR. MODEL_levels <  1 .OR.    &
     &      wet_model_levels <  1 ) then
        ERROR = 1
        GOTO 9999

      End If


!-----------------------------------------------------------------------
! initialisation of output arrays
!-----------------------------------------------------------------------
      do j=1,rows
        do i=1,row_length

          cumulus(i,j)     = .false.

! set false as not used if no convection

          L_shallow(i,j)   = .false.
          L_congestus(i,j) = .false.
          L_congestus2(i,j) = .false.
!
! Requried by Boundatry layer scheme
!
          ntml(i,j) = 1
          nlcl(i,j) = 1
!
! Only used in this routine as no convection call
!
          delthvu(i,j) = 0.0
!
! required by Boundary layer scheme
!
          z_lcl(i,j) = z_half(i,j,nlcl(i,j)+1)

!
! set to zero as not used
!
          z_lcl_uv(i,j)= 0.0

        End Do
      End Do
!-----------------------------------------------------------------------
! initialisation of local arrays
!-----------------------------------------------------------------------
      do j=1,rows
        do i=1,row_length
          topbl(i,j) = .false.
          kcucheck(i,j) = 1
          dtv_min(i,j)  = 0.0
          cin(i,j)      = 0.0
          cape(i,j)     = 0.0
        End Do
      End Do

!-----------------------------------------------------------------------
! set variables
!-----------------------------------------------------------------------

!  Set MBL, "maximum number of boundary levels" for the purposes of
!  boundary layer height calculation.

      MBL = bl_levels - 1

!-----------------------------------------------------------------------
! 1.2 Put level 1 wind components on P grid
!-----------------------------------------------------------------------

! DEPENDS ON: u_to_p
        CALL U_TO_P(U,row_length,rows,1,                                &
     &            off_x, off_y, model_domain,                           &
     &            at_extremity, u_p)

! DEPENDS ON: v_to_p
        CALL V_TO_P(V,row_length,rows,n_rows,1,                         &
     &            off_x, off_y, model_domain,                           &
     &            at_extremity, v_p)

      IF(MODEL_DOMAIN  ==  mt_global) THEN

! Overwrite values of U_P, V_P at the poles with the magnitude of
! the vector wind.

! DEPENDS ON: polar_vector_wind_n
        Call Polar_vector_wind_n(                                       &
     &                       v,                                         &
     &                       sin_theta_longitude,                       &
     &                       cos_theta_longitude, row_length,           &
     &                       n_rows, 1, mag_vector_np,                  &
     &                       dir_vector_np, mag_vector_sp,              &
     &                       dir_vector_sp,                             &
     &                       off_x, off_y, global_row_length,           &
     &                       proc_row_group, at_extremity)
        If (at_extremity(PSouth) ) Then
            DO I=1,ROW_LENGTH
              V_P(I,1) = MAG_VECTOR_SP(1)
              U_P(I,1) = 0.0
            END DO
        End If
        If (at_extremity(PNorth) ) Then
            DO I=1,ROW_LENGTH
              V_P(I,rows) = MAG_VECTOR_NP(1)
              U_P(I,rows) = 0.0
            END DO
        End If

      End If

!-----------------------------------------------------------------------
! 1.3 Calculate layer depths and heights.
!-----------------------------------------------------------------------
      Do K=1, model_levels-1
        Do J=1, rows
          Do I=1, row_length
            r2rho_theta(i,j,k) = rho_theta(i,j,k)*                      &
     &                            (r_theta_levels(i,j,k)/earth_radius)  &
     &                           *(r_theta_levels(i,j,k)/earth_radius)
          End Do
        End Do
      End Do

!-----------------------------------------------------------------------
! 2.0 Calculate various quantities required by the parcel calculations
!-----------------------------------------------------------------------

! ----------------------------------------------------------------------
! 2.1 Calculate T at old time level. (Strictly partway through timestep
!     value).
! ---------------------------------------------------------------------

      Do k = 1, model_levels
         Do j = 1, rows
            Do i = 1, row_length
               T(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
            End Do
         End Do
      End Do

!-----------------------------------------------------------------------
! 2.2 Calculate total water content, QW and Liquid water temperature, TL
!-----------------------------------------------------------------------

      DO K=1,wet_model_levels
        Do j = 1, rows
          DO I=1,row_length
            QW(i,J,K) = Q(i,J,K) + qcl(i,j,k) + qcf(i,j,k)
                                ! BL doc P243.10
            TL(i,J,K) = T(i,J,K) - LCRCP*qcl(i,j,k) - LSRCP*qcf(i,j,k)
                                ! BL doc P243.9

!
! Calculate SVL: conserved variable  svl = Tl+gz(1+(1/epsilon-1)qt)/cp
!   static energy /cp
!
            SVL(i,j,K) = ( TL(i,j,K) + GRCP * z_full(i,j,K) )           &
     &                               * ( 1.0 + c_virtual*QW(i,j,K) )
          End Do
        End Do
      End Do

!
! qsat at surface
!
! DEPENDS ON: qsat_mix
      call qsat_mix(QS_STAR,TSTAR,PSTAR,row_length*rows,l_mixing_ratio)

!-----------------------------------------------------------------------
! 2.3 Set appropriate roughness length
!-----------------------------------------------------------------------
      ! Default z0h for standard TSTAR boundary condition
      DO J=1, rows
      DO I=1, row_length
        IF (LAND_MASK(I,j)) THEN
!         ! Approximate Z0H for land as 0.1.
          Z0H(i,j) = 0.1
        ELSE
          Z0H(i,j) = Z0HSEA
        End If
      End Do ! I
      End Do ! J

      IF ( L_spec_z0 ) THEN
        ! Code to use z0h_scm if Namelist specifies it
        DO J=1, rows
        DO I=1, row_length
          IF ( Z0H_SCM(i,j)  >   0.0 ) THEN
            Z0H(i,j) = Z0H_SCM(i,j)
          End If ! z0h_scm
        End Do ! I
        End Do ! J
      End If
!-----------------------------------------------------------------------
! 3.0 Set topbl to .false. and calculate boundary layer top using
!      a non-local, plume method.
! ----------------------------------------------------------------------
      IF ( .NOT. L_flux_bc) THEN
!      !----------------------------------------------------------------
!      ! Standard UM code for surface T boundary condition:
!      ! Calculate the surface buoyancy flux using
!      ! approximation for unstable CD as 1.5*neutral value (defined
!      ! Garratt p54) and that CH=CD.
!      !----------------------------------------------------------------
       K=1
       do j=1,rows
       do i=1,row_length

         theta1 = theta(i,j,1)
         theta_star = TSTAR(i,j)*((100000.0/PSTAR(i,j))**KAPPA)
         rhostar = PSTAR(i,j) / ( R*TSTAR(i,j) )

         USHEAR = U_P(i,j) - U_0_P(i,j)
         VSHEAR = V_P(i,j) - V_0_P(i,j)
         WSHR2 = MAX (1.0E-6 , USHEAR*USHEAR + VSHEAR*VSHEAR)
         WSHR1 = SQRT(WSHR2)
         CD = 1.5 * ( VKMAN/ALOG(z_full(i,j,1)/Z0H(I,j)) )**2

         If (Land_MASK(i,j)) then         ! land
           wthvbar = WSHR1 * CD * ( theta_star - theta1 )
         Else                             ! sea
           wthvbar = WSHR1 * CD *                                       &
     &     ( theta_star - theta1 + 0.61*theta1*(QS_STAR(i,j)-Q(i,j,1)) )
         End If

         ustar = SQRT(CD * WSHR2)
         fb_surf(i,j) = G * wthvbar /                                   &
     &                       ( rhostar * theta1*(1.0+0.61*Q(i,j,1)) )

         If (fb_surf(i,j)  >   0.0) then
           TV1_SD(i,j) = 1.93 * wthvbar / ( rhostar * ( 75.0 *          &
     &            fb_surf(i,j) + ustar*ustar*ustar)**(1.0/3.0) )
         Else
           TV1_SD(i,j) = 0.0
         End If

       End Do ! I
       End Do ! J

      ELSE ! if L_flux_bc
!       !-----------------------------------------------------------
!       ! Code for specified surface flux boundary condition.
!       ! Assumes a saturated surface (ie. only appropriate for sea)
!       !-----------------------------------------------------------
        DO J=1, rows
        DO I=1, row_length
          ! For taylor expansion about T0=SL(K=1)
          TSTAR(I,J) = T(I,J,1)+GRCP*Z_FULL(I,J,1)
        End Do
        End Do

! DEPENDS ON: qsat_mix
        call qsat_mix(QS_STAR,TSTAR,PSTAR,row_length*rows,l_mixing_ratio)

        DO J=1, rows
        DO I=1, row_length
          DQSDT(I,J) = (EPSILON * LC * QS_STAR(I,J))                    &
     &                 / ( R * TSTAR(I,J) * TSTAR(I,J) )
        End Do
        End Do

        do j=1,rows
        do i=1,row_length

         USHEAR = U_P(I,j) - U_0_P(I,j)
         VSHEAR = V_P(I,j) - V_0_P(I,j)
!        ! Need to have a higher minimum wind speed limit with
!        ! specified fluxes in order not to generate huge TSTAR
         WSHR2 = MAX (0.1, USHEAR*USHEAR + VSHEAR*VSHEAR)
         WSHR1 = SQRT(WSHR2)

         ! Calculate WTHV from namelist flux_h and flux_e (in W/m2)
         WTHVBAR = ((flux_h(I,j)/CP)+0.61*(flux_e(I,j)/LC))             &
     &           * ((100000.0/PSTAR(I,j))**KAPPA)

         CD = 1.5 * ( VKMAN/ALOG(Z_FULL(I,j,1)/Z0H(i,j)) )**2

         THETA1 = theta(I,j,1)

         ! Taylor expansion for qsat(T*) about SL(k=1)
         TSTAR(I,j) = ( THETA1 + (WTHVBAR/(WSHR1*CD))                   &
     &              -   0.61*THETA1                                     &
     &              *   (QS_STAR(I,j)-Q(I,j,1)-DQSDT(I,j)*TSTAR(I,j)))  &
     &              /   ( (100000.0/PSTAR(I,j))**KAPPA +                &
     &                       0.61*THETA1*DQSDT(I,j) )

         RHOSTAR = PSTAR(I,j) / ( R*TSTAR(I,j) )

         USTAR = SQRT(CD * WSHR2)
         FB_SURF(I,j) = G * WTHVBAR /                                   &
     &                       ( RHOSTAR * THETA1*(1.0+0.61*Q(I,j,1)) )

         IF (FB_SURF(I,j)  >   0.0) THEN
           TV1_SD(I,j) = 1.93 * WTHVBAR / ( RHOSTAR * ( 75.0 *          &
     &                   FB_SURF(I,j) + USTAR*USTAR*USTAR)**(1.0/3.0) )
         ELSE
           TV1_SD(I,j) = 0.0
         End If

        End Do ! I
        End Do ! J

      End If !  L_flux_bc

      do j=1,rows
      do i=1,row_length

        k_plume(i,j) = 1

        UNSTABLE(i,j) = (fb_surf(i,j)  >   0.0)

        If ( UNSTABLE(i,j) ) then
!-----------------------------------------------------------------------
! Only perform parcel ascent If UNSTABLE
! Start plume ascent from grid-level above top of surface layer, taken
! to be at a height, z_surf, given by 0.1*zh
!-----------------------------------------------------------------------
         Z_SURF = 0.1 * zh(i,j)

         Do While( z_full(i,j,k_plume(i,j))  <   Z_SURF .and.           &
!                   ! not reached Z_SURF
     &              SVL(i,j,k_plume(i,j)+1)  <   SVL(i,j,k_plume(i,j)) )
!                   ! not reached inversion

            k_plume(i,j) = k_plume(i,j) + 1

         End Do

         sl_plume(i,j) = TL(i,j,k_plume(i,j))                           &
     &                        + GRCP * z_full(i,j,k_plume(i,j))
         thv_pert(i,j) = MAX( A_PLUME,                                  &
     &         MIN( MAX_T_GRAD*zh(i,j), B_PLUME*TV1_SD(i,j) ) )
         qw_plume(i,j) = QW(i,j,k_plume(i,j))

! Added for more acturate parcel cal later

         th_ref(i,j) = tl(i,j,k_plume(i,j))                             &
     &                        /exner_theta_levels(i,j,k_plume(i,j))
         th_par_km1(i,j) = th_ref(i,j)

!-----------------------------------------------------------------------
!! 0.2 Calculate temperature and pressure of lifting condensation level
!!     using approximations from Bolton (1980)
!-----------------------------------------------------------------------
         If (l_mixing_ratio) then
! expression for mixing ratio
           vap_press = 0.01*Q(i,j,k_plume(i,j)) *                       &
     &                                 P_theta_lev(i,j,k_plume(i,j))    &
     &                      / (EPSILON+Q(i,j,k_plume(i,j)) )
         Else
! expression for specific humidity
           vap_press = Q(i,j,k_plume(i,j)) *                            &
     &         P_theta_lev(i,j,k_plume(i,j)) / ( 100.0*EPSILON )
         End If
         If (vap_press  >   0.0) then
           T_LCL(i,j) = a_bolton + b_bolton/                            &
     &                      (c_bolton*LOG(T(i,j,k_plume(i,j)))          &
     &                                    - LOG(vap_press) - d_bolton )
           P_LCL(i,j) = P_theta_lev(i,j,k_plume(i,j)) *                 &
     &                ( T_LCL(i,j) / T(i,j,k_plume(i,j)) )**(1.0/KAPPA)
         Else
           P_LCL(i,j) = pstar(i,j)
         End If

        End If   ! test on UNSTABLE
      End Do
      End Do
!
! Reset zh  (at this point in the code ntml is initialised as =1)
!
      do j=1,rows
        do i=1,row_length
          zh(i,j) = z_half(i,j,ntml(i,j)+1)
        End Do
      End Do
!
!-----------------------------------------------------------------------
! Find NLCL
!-----------------------------------------------------------------------
!
!    ---------------   p      nlcl+1  , p_theta(nlcl+2)
!
!    - - - - - - - -   uv     nlcl+1,  z_lcl , p(nlcl+1)    either
!     + + + + + + + +   lcl, Plcl, not a model level        lower part
!    ---------------   p      nlcl , p_theta(nlcl+1)         of layer
!
!    - - - - - - - -   uv     nlcl   p(nlcl)
!
!-----------------------------------------------------------------------
!
!    ---------------   p      nlcl+1  , p_theta(nlcl+2)
!     + + + + + + + +   lcl, Plcl, not a model level
!
!    - - - - - - - -   uv     nlcl+1,  z_lcl , p(nlcl+1)     or
!                                                          upper part
!    ---------------   p      nlcl , p_theta_lev(nlcl+1)        of layer
!
!    - - - - - - - -   uv     nlcl   p(nlcl)
!
!-----------------------------------------------------------------------
      DO  K = 2,wet_model_levels
        do j=1,rows
         do i=1,row_length
          If ( UNSTABLE(i,j) ) then
           If ( P_LCL(i,j)  <   P(i,j,K) ) then
             nlcl(i,j) = K-1
             z_lcl(i,j) = z_half(i,j,nlcl(i,j)+1)

           End If
          End If
         End Do
        End Do
      End Do

!-----------------------------------------------------------------------
! 4.0 Parcel ascent - only perform parcel ascent If UNSTABLE
!-----------------------------------------------------------------------
! Note initial parcel conditions different from those used in the G-R
! mass flux convection scheme.
!
! ---------------------------------------------------------------------
! Calculate parcel water by linearising qsat about the Parcel's
! temperature extrapolated up to the next grid_level.
!----------------------------------------------------------------------

       DO  K = 1,wet_model_levels


        do j=1,rows
         do i=1,row_length
           If ( UNSTABLE(i,j) ) then
             t_ref(i,j) = th_ref(i,j)*exner_theta_levels(i,j,k)
           End If
         End Do
        End Do
! DEPENDS ON: qsat_mix
        call qsat_mix(qsat_lev,t_ref,p_theta_lev(1,1,k)                 &
     &                                 ,row_length*rows,l_mixing_ratio)

        do j=1,rows
        do i=1,row_length
         If ( UNSTABLE(i,j) ) then
          If(T_ref(i,j) >  TM) then
             lrcp_const = lcrcp
             l_const    = lc
          Else
             lrcp_const = lsrcp
             l_const    = ls
          End If


          dq_sat_env = epsilon*l_const*qsat_lev(i,j)/(R*T_ref(i,j)**2.0)

          q_liq_parc = MAX( 0.0, ( qw_plume(i,j) - QSat_lev(i,j)        &
     &       -dq_sat_env*( sl_plume(i,j)-grcp*z_full(i,j,K)-T_ref(i,j) )&
     &                             ) / (1.0+lrcp_const*dq_sat_env) )
          q_liq_env  = MAX( 0.0, ( qw(i,j,K) - QSat_lev(i,j)            &
     &        -dq_sat_env*( TL(i,j,K)               - T_ref(i,j) )      &
     &                             ) / (1.0+Lrcp_const*dq_sat_env) )
!
! add on the dIfference in the environment's ql as calculated by the
! UM cloud scheme (using some RH_CRIT value) and what it
! would be If RH_CRIT=1. This then imitates partial condensation
! in the parcel.
!
          q_liq_parc = q_liq_parc + qcl(i,j,k)                          &
     &                            + qcf(i,j,k)- q_liq_env
          T_PARC=sl_plume(i,j)-GRCP*z_full(i,j,K)                       &
     &                                  +lrcp_const*q_liq_parc

! May need to recalculate if T_parc is > Tm and T_ref < Tm

          If (T_ref(i,j) <= TM.and.T_parc >  TM) then

! recalculate using corrected latent heats
            lrcp_const_parc = lcrcp
            q_liq_parc = MAX( 0.0, ( qw_plume(i,j) - Qsat_lev(i,j)      &
     &      -dq_sat_env*( sl_plume(i,j)-grcp*z_full(i,j,K)-T_ref(i,j) ) &
     &                     ) / (1.0+lrcp_const_parc*dq_sat_env) )
            q_liq_parc = q_liq_parc + qcl(i,j,k)                        &
     &                               + qcf(i,j,k)- q_liq_env

! revised at parcel calculation

            T_PARC=sl_plume(i,j)-GRCP*z_full(i,j,K)                     &
     &                                  +lrcp_const_parc*q_liq_parc

          End If

          q_vap_parc=qw_plume(i,j)-q_liq_parc
!
          t_dens_parc(i,j,k)=T_PARC*                                    &
     &                    (1.0+c_virtual*q_vap_parc-q_liq_parc)


          t_dens_env(i,j,k)=T(i,j,K)*                                   &
     &                    (1.0+c_virtual*Q(i,j,K)-qcl(i,j,k)-qcf(i,j,k))

          env_svl(i,j,k) = t_dens_env(i,j,k)  + GRCP*z_full(i,j,K)
          par_svl(i,j,k) = t_dens_parc(i,j,k) + GRCP*z_full(i,j,K)

          if (k >= 2) then

!         !-------------------------------------------------------------
!         ! Find vertical gradients in parcel and environment SVL
!         ! (using values from level below (i.e. K-1)).
!         !-------------------------------------------------------------

            dz = z_full(i,j,K) - z_full(i,j,K-1)

            dpar_bydz(i,j,k) = (par_svl(i,j,k) - par_svl(i,j,k-1))/dz
            denv_bydz(i,j,k) = (env_svl(i,j,k) - env_svl(i,j,k-1))/dz


          End If   ! test on k

! calculate t_ref for next level
          If (k >  1 .and. k <   wet_model_levels-1) then
            z_pr = (z_full(i,j,k+1)-z_full(i,j,k))                      &
     &                         /(z_full(i,j,k)-z_full(i,j,K-1))
            th_par = t_parc/exner_theta_levels(i,j,k)
            th_ref(i,j) = th_par*(1.+z_pr) - th_par_km1(i,j)*z_pr
            th_par_km1(i,j) = th_par
          End If

         End If  ! test on UNSTABLE
        End Do   ! i loop
        End Do   ! j loop
      End Do      ! level loop

!-----------------------------------------------------------------------
! tests on parcel ascent
!-----------------------------------------------------------------------
!   Now compare plume s_VL with each model layer s_VL in turn to
!     find the first time that plume has negative buoyancy.
!-----------------------------------------------------------------------

      DO  K = 2,wet_model_levels

        do j=1,rows
        do i=1,row_length

         If ( UNSTABLE(i,j) ) then
!-----------------------------------------------------------------------
! Only perform tests if parcel ascent If UNSTABLE
!-----------------------------------------------------------------------

!  Find level just below 2.5km (for use in Cu diagnosis)

          If ( z_full(i,j,K)  >   2500.0                                &
     &         .and. kcucheck(i,j)  ==  1 ) kcucheck(i,j) = K-1

!-----------------------------------------------------------------------
! Set flag to true when level below is at least one level above the lcl
! and above the lcl transition zone
! Code implies ABOVE_LCL at NLCL+3 or greater.

          If (k-1 >  nlcl(i,j)+1                                        &
     &                   .and. z_full(i,j,k-1) >  1.1*z_lcl(i,j)) then
            above_lcl(i,j)=.true.
          Else
            above_lcl(i,j)=.false.
          End If

!-----------------------------------------------------------------------
! Tests applied to find parcel top
!
!-----------------------------------------------------------------------

          If ( .not.topbl(i,j) .and. K  >   k_plume(i,j) .and.          &
     & (  (t_dens_parc(i,j,k)-t_dens_env(i,j,k) <=  - thv_pert(i,j)).OR.&
!
!                      plume non buoyant
!
     & (above_lcl(i,j).and.(denv_bydz(i,j,k) >  1.25*dpar_bydz(i,j,k))) &
!
!                      or environmental virtual temperature gradient
!                      signIficantly larger than parcel gradient
!                      above lIfting condensation level
!
     &           .OR. (K  >   wet_model_levels-1)                       &
!                      or reached top of model
     &         )                                                        &
     &         ) then
!
            topbl(i,j) = .TRUE.
            zh(i,j) = z_half(i,j,K)
            ntml(i,j) = K-1
            If ( delthvu(i,j)  >   0.0) then
! compensate for any negative buoyancy of parcel in cloud layer
              delthvu(i,j) = delthvu(i,j) - dtv_min(i,j) *              &
     &                                    ( zh(i,j) - z_lcl(i,j) )
            End If
          End If

!-----------------------------------------------------------------------
! While doing parcel ascent
! (a) find minimum buoyancy
! (b) integrate CAPE over the ascent
!-----------------------------------------------------------------------

          If (K  >   nlcl(i,j) .and. .not. topbl(i,j)) then
            dtv_min(i,j) = MIN( dtv_min(i,j),                           &
     &                         (t_dens_parc(i,j,k)-t_dens_env(i,j,k))   &
     &                    /exner_theta_levels(i,j,k)  )

            delthvu(i,j) = delthvu(i,j) +                               &
     &                      (t_dens_parc(i,j,k)-t_dens_env(i,j,k)) *    &
     &                     ( z_half(i,j,K+1) - z_half(i,j,K) )          &
     &                        /exner_theta_levels(i,j,k)

! calculation of CIN and CAPE from profiles

            inc =g*(t_dens_parc(i,j,k)-t_dens_env(i,j,k))               &
     &             *(z_half(i,j,k+1) - z_half(i,j,k))/t_dens_env(i,j,k)
            If (inc <  0.0) then
              CIN(i,j) = CIN(i,j) + inc
            End If
            CAPE(i,j) = CAPE(i,j) + inc
          End If
!
!-----------------------------------------------------------------------
         End If  ! test on UNSTABLE
        End Do   ! i loop
        End Do   ! j loop
      End Do

!-----------------------------------------------------------------------
!! 0.4 Save parcel ascent top: this will be used to allow mixing and
!!     entrainment into decoupled Sc of single layer thickness when it
!!     occurs above Cu.
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        zhpar(i,j) = zh(i,j)
        ntpar(i,j) = ntml(i,j)
      End Do
      End Do

!-----------------------------------------------------------------------
!     Test height derived above against lifting condensation level
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
!-----------------------------------------------------------------------
!     Check lifting condensation levels against height of parcel ascent,
!     if lifting condensation level lower than parcel ascent, and is
!     within BL_LEVELS, then decide
!     on type of cloudy layer. If lifting condensation level at or below
!     low grid point, assume fog layer and turbulent mixing. For
!     gradient tests assume any if LCL and top of parcel ascent is less
!     than two levels then stratocumulus.
!-----------------------------------------------------------------------
! altered for 5A to ensure nptar -ntml >=2
! (at this point ntml holds ntpar)

        IF ( ntml(i,j)-nlcl(i,j)  >=  3                                 &
     &                        .and. nlcl(i,j)  >   k_plume(i,j)         &
     &                              .and. nlcl(i,j)  <   MBL-1 ) then
!-----------------------------------------------------------------------
!     Cloudy boundary layer, diagnose whether stratocumulus or cumulus.
!     For stratocumulus top of mixed layer = zh
!     For cumulus top of mixed layer = ZLCL
!     NTML >= MBL indicates convection. (ntml hold parcel top here.)
!     Diagnosis is done by comparing gradients
!-----------------------------------------------------------------------

          If (ntml(i,j)  >=  MBL) then
            cumulus(i,j) = .TRUE.
          Else


! Current test is against a height of ~<2.5km
! This could be replaced by a scale height if a suitable method
! for determining a sensible height was possible from profile/cumulus
! depth information available in this routine

            If (ntml(i,j)  >   kcucheck(i,j)                            &
     &             .and. nlcl(i,j)  <=  kcucheck(i,j)-2) then

              grad_cld =  ABS( QW(i,j,kcucheck(i,j)) -                  &
     &                                  QW(i,j,nlcl(i,j)) ) /           &
     &          ( z_full(i,j,kcucheck(i,j)) - z_full(i,j,nlcl(i,j)) )
            Else
              grad_cld =  ABS( QW(i,j,ntml(i,j)) -                      &
     &                                  QW(i,j,nlcl(i,j)) ) /           &
     &              ( z_full(i,j,ntml(i,j)) - z_full(i,j,nlcl(i,j)) )
            End If

            grad_sub =  ABS( QW(i,j,nlcl(i,j)) -                        &
     &                                QW(i,j,k_plume(i,j)) ) /          &
     &           ( z_full(i,j,nlcl(i,j)) - z_full(i,j,k_plume(i,j)) )

            If (grad_cld  >   1.10*grad_sub) then
!-----------------------------------------------------------------------
!     Not well mixed, however it is possible that the depth of a well
!     mixed boundary layer has increased but not yet been mixed yet so
!     test gradient from next level down.
!     Note typical cumulus profiles are expected to have a fairly
!     uniform q profile from the surface to the cloud base and then a
!     decreasing profile of q above this in the cloud. Typical the
!     decreasing gradient from the cloud base to 2.5km will be the
!     order of > 1.10 the below cloud value.
!-----------------------------------------------------------------------

! test against a height ~ 2.5km

              If (ntml(i,j)  <=  kcucheck(i,j)) then
              grad_cld =  ABS( QW(i,j,ntml(i,j)-1) -                    &
     &                                  QW(i,j,nlcl(i,j)) ) /           &
     &           ( z_full(i,j,ntml(i,j)-1) - z_full(i,j,nlcl(i,j)) )
              End If

              If ( grad_cld  >   1.10*grad_sub) then
!-----------------------------------------------------------------------
!      Diagnose a cumulus layer
!-----------------------------------------------------------------------
                cumulus(i,j) = .TRUE.
              End If
            Else

! Diagnosed well-mixed, but now check that LCL hasn't risen or fallen
! and not yet been mixed (so could have been erroneously identIfied as
! well-mixed)

! First check using level below (recalculate grad_sub)

              If (nlcl(i,j) - k_plume(i,j)  >=  2) then
                 grad_sub =  ABS( QW(i,j,nlcl(i,j)-1) -                 &
     &                                QW(i,j,k_plume(i,j)) ) /          &
     &           ( z_full(i,j,nlcl(i,j)-1) - z_full(i,j,k_plume(i,j)) )

                 If ( grad_cld  >   1.10*grad_sub) then
                   cumulus(i,j) =.TRUE.
                 End If

              End If

! If still diagnosing well-mixed, check using level above
! (recalculate grad_cld)

              If (.not. cumulus(i,j) ) then

               If (ntml(i,j)  >   kcucheck(i,j)                         &
     &             .and. nlcl(i,j)  <=  kcucheck(i,j)-2) then

                grad_cld =  ABS( QW(i,j,kcucheck(i,j)) -                &
     &                                  QW(i,j,nlcl(i,j)+1) ) /         &
     &         ( z_full(i,j,kcucheck(i,j)) - z_full(i,j,nlcl(i,j)+1) )
               Else
                grad_cld =  ABS( QW(i,j,ntml(i,j)) -                    &
     &                                  QW(i,j,nlcl(i,j)+1) ) /         &
     &            ( z_full(i,j,ntml(i,j)) - z_full(i,j,nlcl(i,j)+1) )
               End If

               If ( grad_cld  >   1.10*grad_sub) then
                 cumulus(i,j) =.TRUE.
               End If

              End If
            End If
          End If
        End If

!-----------------------------------------------------------------------
!      Check that a cumulus layer has not been erroneously diagnosed in
!      a deep cloudy region
!      As the above checks are done on the total water rather than q it
!      is possible the conditions can be met in areas where the level of
!      prognostic qcl or qcf is high. The type of mistake is only
!      thought to occur over land.
!-----------------------------------------------------------------------
        K=nlcl(i,j)

        If ( Land_MASK(i,j) .and. cumulus(i,j) .and.                    &
     &                                 ntpar(i,j)  <   MBL ) then
          Do While( K  <=  ntpar(i,j) .and. cloud_fraction(i,j,K)       &
     &                                            >=  SC_CFTOL )
            K = K + 1
          End Do
          If (K  ==  ntpar(i,j)+1) cumulus(i,j) = .false.
        End If


!-----------------------------------------------------------------------
! No separation into shallow and congestus as no convection scheme
! will be called alter.
!-----------------------------------------------------------------------


        IF ( CUMULUS(I,j) ) THEN

! No separation into shallow deep or congestus required if no convection
! scheme in use.

!-----------------------------------------------------------------------
!      Set mixed layer depth to z_lcl
!-----------------------------------------------------------------------
          If (P_LCL(i,j)  <   (P_theta_lev(i,j,nlcl(i,j)+1))) then
!-----------------------------------------------------------------------
!      If LCL is diagnosed in the upper half of the layer set z_lcl to
!      the height of the upper layer interface
!      (in code above LCL is always set to the lower interface).
!-----------------------------------------------------------------------
             nlcl(i,j) = nlcl(i,j)+1
             z_lcl(i,j) = z_half(i,j,nlcl(i,j)+1)

          End If
          zh(i,j) = z_lcl(i,j)
          ntml(i,j) = nlcl(i,j)

!      If cumulus has been diagnosed but delthvu is negative, reset
!      cumulus to FALSE but leave zh and NTML at LCL

          If (delthvu(i,j)  <=  0.0) then
            cumulus(i,j) = .false.
          End If

        Else

!-----------------------------------------------------------------------
!      If not cumulus, reset parameters to within bl_levels
!-----------------------------------------------------------------------
          If (ntml(i,j)  >   MBL) then
            ntml(i,j)  = MBL
            ntpar(i,j) = MBL
            zh(i,j)    = z_half(i,j,MBL+1)
            zhpar(i,j) = zh(i,j)
          End If
          If (nlcl(i,j)  >   MBL) then
            nlcl(i,j)    = MBL
            z_lcl(i,j)   = zh(i,j)
          End If

        End If

      End Do
      End Do

#if defined(SCMA)
        If (L_SCMDiags(SCMDiag_bl)) Then
! DEPENDS ON: scmoutput
          call SCMoutput(env_svl,                                       &
               'thv_env','Environment virtual temperature + gz/cp','K', &
               t_avg,d_wet,default_streams,'',RoutineName)
! DEPENDS ON: scmoutput
          call SCMoutput(par_svl,                                       &
               'thv_par','Parcel virtual temperature + gz/cp','K',      &
               t_avg,d_wet,default_streams,'',RoutineName)

        End If
#endif

9999   CONTINUE  ! Branch for error exit.

      RETURN
      END SUBROUTINE CONV_DIAG
#endif
