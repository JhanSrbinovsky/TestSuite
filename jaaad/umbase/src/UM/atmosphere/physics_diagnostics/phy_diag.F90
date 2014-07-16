#if defined(A16_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculate various diagnostics related to dynamics variables.
!
! Subroutine Interface:
      SUBROUTINE Phy_diag(                                              &
! Primary data: in
     & p_star,p,rho,u,v,w,theta,q                                       &
     &,p_theta_levels,exner_rho_levels,exner_theta_levels               &
! Grid sizes and definition: in
     &,rows,n_rows,row_length,model_levels,wet_levels,bl_levels         &
     &,theta_field_size,u_field_size,v_field_size                       &
     &,global_row_length, global_rows                                   &
     &,eta_theta_levels,eta_rho_levels                                  &
     &,r_theta_levels,r_rho_levels                                      &
     &,delta_lambda,delta_phi,sec_theta_latitude                        &
! Control information: in
     &,Model_domain,npmsl_height                                        &
! Pressure levels for output arrays: in
     &,hts_p_press,t_p_press                                            &
     &,RHice_p_press,RHwat_p_press,wbpt_p_press                         &
! Flags to request each diagnostic output field: in
     &,qT_m                                                             &
     &,qhts_theta                                                       &
     &,qhts_p,qt_p,qRHice_p,qwbpt_p                                     &
     &,qp_MSL                                                           &
     &,qhts_rho                                                         &
     &,qRHwat_p                                                         &
! Diagnostics lengths: in
     &,hts_p_levs,t_p_levs                                              &
     &,RHice_p_levs,RHwat_p_levs,wbpt_p_levs                            &
! Diagnostic arrays: out
     &,T                                                                &
     &,hts_theta                                                        &
     &,hts_p,t_p,RHice_p,wbpt_p                                         &
     &,p_MSL                                                            &
     &,hts_rho                                                          &
     &,RHwat_p                                                          &
     &     )

#if defined(FLUME)
! FLUME-STASH 
USE MATMFlumeModel, only:flumeSendDiag
USE SharedID
USE flumerun
#endif

      IMPLICIT NONE
!
! Description:
!   Calculate physics-related diagnostics - held in STASH section 16 -
!   which may include interpolation onto pressure surfaces. Diagnostics
!   currently supported:
!   STASH item
!     4 temperature on model levels
!   201 geopotential height on theta levels
!   202 geopotential height on pressure surfaces
!   203 temperature         on pressure surfaces
!   204 relative humidity wrt ice on pressure surfaces
!   205 wet bulb potential temperature on pressure surfaces
!   222 mean sea level pressure
!   255 geopotential height on rho levels
!   256 relative humidity wrt water on pressure surface
!
! Method:
!   Required level lists and logical switches are determined by the
!   calling routine from STASH requests and STASHflags.
!   Primary model data is input, and each diagnostic is calculated - if
!   its flag is set - in simple sequential order. Where the extraction
!   of the diagnostic quantity requires further calculation, a lower
!   level diagnostic-specific routine is called.
!
! Current Code Owner: <Rick Rawlins>
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.0  02/06/99   Extensive revision of PHYDIA1A deck for 'C-P C
!                  dynamics upgrade' project. Rick Rawlins.
!  5.1  14/02/00   Correct interface to routine vert_interp_mdi2 for
!                  RH diagnostic (16,204).
!                  Replace rmdi_pp by rmdi. R Rawlins
!  5.1  17/02/00   Add diagnostic for T on model levels. R Rawlins
!  5.2  19/01/01   Add diagnostic for wet bulb potential temperature
!                  (16,205). R Rawlins
!  5.2  19/03/01   Linear interpolate in exner instead of height.
!                  (T on p levels). C. Wilson
!  5.2  19/03/01   Linear interpolate in exner instead of height.
!                  C. Wilson
!  5.2  19/03/01   Change to use exner pressures instead of p in
!                  call to vert_interp_mdi2 (now vert_interp2)
! 6/6/01   5.3         Pass gc_all_proc_group to CALC_PMSL  P.Burton
!  5.3  07/11/01      Add Clive Wilson change to PMSL calc    A.Malcolm
!  5.3  19/11/01   Subroutine Thetaw has been generalized.
!                  Population of pressure array now done in PHYDIAG.
!                  D.M. Goddard
!  5.3  05/12/01   Add geopotential height on theta levels (16,201)
!                  and geopotential height on rho levels (16,255)
!                  D.M. Goddard
!  5.4  10/04/02   Add Relative Humidity wrt water on presure levels
!                  (16,256). D.M. Goddard
!  6.0  26/09/03   Fix reference to uninitialised height array
!                  when (16,255) is requested without (16,202).
!                  Restrict (16,201) to 1,model_levels.
!                  P. Selwood
!  6.1  08/07/04   send P_theta_levels into vert_h_onto_p
!                                                    Michael Hughes
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!   Documentation: UMDP 80
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Global variables (*CALLed COMDECKs etc...):
#include "parvars.h"
#include "cphyscon.h"
#include "interpor.h"

! Subroutine arguments
!   Scalar arguments with intent(in):
! Grid sizes:
      INTEGER                                                           &
     & rows,n_rows,row_length,model_levels,wet_levels,bl_levels         &
     &,theta_field_size,u_field_size,v_field_size                       &
     &,global_row_length, global_rows                                   &
! Control information: in
     &,Model_domain                                                     &
                       ! Domain of atmosphere model:
!                                   global,LAM,cyclic LAM,single column

! Diagnostics lengths: IN
     &,hts_p_levs                                                       &
                      ! NO OF LEVS ON WHICH TO INTERP hts_p
     &,t_p_levs                                                         &
                      ! NO OF LEVS ON WHICH TO INTERP t_p
     &,RHice_p_levs                                                     &
                      ! NO OF LEVS ON WHICH TO INTERP RHice_p
     &,RHwat_p_levs                                                     &
                      ! NO OF LEVS ON WHICH TO INTERP RHwat_p
     &,wbpt_p_levs    ! NO OF LEVS ON WHICH TO INTERP wbpt_p

! Grid definition: IN
      REAL                                                              &
     & delta_lambda,delta_phi                                           &
! Orographic height threshold for new pmsl calculation: IN
     &,npmsl_height

! Flags to request each diagnostic output field: IN
      LOGICAL                                                           &
     & qhts_p                                                           &
                   ! Flag for geopotential heights   on pressure levels
     &,qt_p                                                             &
                   ! Flag for temperature            on pressure levels
     &,qRHice_p                                                         &
                   ! Flag for relatice humidity wrt ice
                   !   on pressure surfaces
     &,qRHwat_p                                                         &
                   ! Flag for relative humidity wrt water
                   !   on pressure surfaces
     &,qwbpt_p                                                          &
                   ! Flag for wet bulb potential T   on pressure levels
     &,qp_MSL                                                           &
                   ! Flag for mean sea level pressure
     &,qT_m                                                             &
                   ! Flag for temperature on model levels
     &,qhts_theta                                                       &
                   ! Flag for geopotential heights   on theta levels
     &,qhts_rho    ! Flag for geopotential heights   on rho levels
!   Array  arguments with intent(in):
! Primary data: IN
      REAL                                                              &
     & p_star(row_length, rows)                                         &
     &,p    (1-offx:row_length+offx, 1-offy:rows+offy,    model_levels) &
     &,rho  (1-offx:row_length+offx, 1-offy:rows+offy,    model_levels) &
     &,u    (1-offx:row_length+offx, 1-offy:rows+offy,    model_levels) &
     &,v    (1-offx:row_length+offx, 1-offy:n_rows+offy,  model_levels) &
     &,w    (1-offx:row_length+offx, 1-offy:rows+offy,  0:model_levels) &
     &,theta(1-offx:row_length+offx, 1-offy:rows+offy,    model_levels) &
     &,q (1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,  wet_levels) &
     &,p_theta_levels                                                   &
     &      (1-offx:row_length+offx, 1-offy:rows+offy,    model_levels) &
     &,exner_rho_levels                                                 &
     &      (1-offx:row_length+offx, 1-offy:rows+offy,    model_levels) &
     &,exner_theta_levels                                               &
     &      (1-offx:row_length+offx, 1-offy:rows+offy,    model_levels) &


! Vertical grid definition: IN
     &,eta_theta_levels(0:model_levels)                                 &
                                        ! vertical grid for theta vars
     &,eta_rho_levels    (model_levels)                                 &
                                        ! vertical grid for rho   vars
     &,r_theta_levels  (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j,       0:model_levels)     &
     &,r_rho_levels    (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j,         model_levels)     &
! Derived from horizontal grid: IN
     &,sec_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy)      &

! Pressure levels for output arrays: IN
     &,hts_p_press(hts_p_levs)                                          &
                                      ! for heights      on p surfaces
     &,t_p_press  (  t_p_levs)                                          &
                                      ! for temperature  on p surfaces
     &,RHice_p_press ( RHice_p_levs)                                    &
                                      ! for r.h. wrt ice on p surf
     &,RHwat_p_press ( RHwat_p_levs)                                    &
                                      ! for r.h. wrt water on p surf
     &,wbpt_p_press(wbpt_p_levs)      ! for wet bulb T   on p surfaces

!   Scalar arguments with intent(InOut):

!   Array  arguments with intent(InOut):

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
! Diagnostic arrays: OUT
      REAL                                                              &
     & hts_p(row_length,rows,hts_p_levs)                                &
                                         ! heights           at p levels
     &,t_p  (row_length,rows,t_p_levs)                                  &
                                         ! temperature       at p levels
     &,RHice_p (row_length,rows,RHice_p_levs)                           &
                                              ! rh wrt ice at p levels
     &,RHwat_p (row_length,rows,RHwat_p_levs)                           &
                                              ! rh wrt water at p-levels
     &,wbpt_p(row_length,rows,wbpt_p_levs)                              &
                                          !wet bulb pot temp at p levels
     &,p_MSL(row_length,rows)                                           &
                                         ! pressure at mean sea level
     &,T(row_length,rows,model_levels)                                  &
                                         ! temperature on model levels
     &,hts_theta(row_length,rows,model_levels)                          &
                                         ! heights on theta levels
     &,hts_rho(row_length,rows,model_levels)
                                         ! heights on rho levels
! Local parameters:
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='Phy_diag')
      LOGICAL, PARAMETER ::                                             &
     & l_potential=.true.              ! Wet bulb potential temperature
                                       ! required from subroutine ThetaW
! Local scalars:
!   ErrorStatus
      INTEGER      ErrorStatus                                          &
                                        ! Error flag (0 = OK)
     &,i,j,k                                                            &
                                        ! loop counters
     &,interp_order    !  order of vertical interpolation

      CHARACTER*256                                                     &
     & CMessage         ! Error message if return code >0

      REAL                                                              &
     & dummy                                                            &
                                ! dummy argument - not referenced
     &,pressure_pa                                                      &
                                ! pressure in pascals
     &,pressure_ex              ! exner pressure

      ! Arguments for FlumeSendDiag
      INTEGER im ! model
      INTEGER is ! section
      INTEGER ie ! item
      INTEGER levels ! number of vertical levels

! Local dynamic arrays:
      REAL                                                              &
     & p_at_theta                                                       &
     &       (row_length,rows,model_levels)                             &
                                            ! pressure at theta points
     &,rh    (row_length,rows,model_levels)                             &
                                            ! workspace for RH
     &,height(row_length,rows,model_levels)                             &
                                            ! height at rho levels
     &,work1  (row_length,rows)                                         &
                                            ! temporary space
     &,work2  (row_length,rows)                                         &
                                            ! temporary space
     &,work3 (row_length,rows)              ! temporary space

! Function & Subroutine calls:
      External Calc_PMSL,Ereport,Qsat,Thetaw,T_vert_interp_to_p         &
     &,vert_interp2,vert_h_onto_p

!- End of header


! ----------------------------------------------------------------------
! Section 0.  Initialisation.
! ----------------------------------------------------------------------

! set Error code to zero
      ErrorStatus = 0

! Set order of vertical interpolation
      interp_order = interp_order_linear

!  Calculate p at theta points. Store in p_at_theta

      DO k = 1, model_levels
        DO j = 1, rows
        DO i = 1, row_length
             p_at_theta(i,j,k) = p_theta_levels(i,j,k)
        ENDDO ! i
        ENDDO ! j
      ENDDO ! k

!   Remove radius of earth from rho levels to create geopotential height
!   of rho level.
!   Required for either 202 (geopotential height on pressure levels) or
!                       255 (geopotential height on rho levels)
      IF (qhts_p .OR. qhts_rho) THEN

      DO k = 1, model_levels
        DO j = 1, rows
        DO i = 1, row_length
             height(i,j,k) = r_rho_levels(i,j,k) - earth_radius
        ENDDO ! i
        ENDDO ! j
      ENDDO ! k

      ENDIF ! on relevant STASHflags

!   Calculate temperature at theta points
      IF(qt_p .OR. qRHice_p .OR. qRHwat_p .OR. qT_m .OR. qwbpt_p) THEN

      DO k = 1, model_levels
        DO j = 1, rows
        DO i = 1, row_length
             T(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
        ENDDO ! i
        ENDDO ! j
      ENDDO ! k

      ENDIF ! on relevant STASHflags

! ----------------------------------------------------------------------
! STASH item 201: geopotential height      on  theta levels
! ----------------------------------------------------------------------
      IF(qhts_theta) THEN

!   Remove radius of earth from height field

        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              hts_theta(i,j,k) = r_theta_levels(i,j,k) - earth_radius
            ENDDO ! i
          ENDDO ! j
        ENDDO ! k

#if defined(FLUME)
!  This item is output via copydiag3d so call flumeSendDiag not required
#endif
      ENDIF ! on relevant STASHflags

      IF(qhts_rho) THEN
! ----------------------------------------------------------------------
! STASH item 255: geopotential height on  rho levels
! ----------------------------------------------------------------------

        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              hts_rho(i,j,k) = height(i,j,k)
            ENDDO ! i
          ENDDO ! j
        ENDDO ! k

#if defined(FLUME)
!  This item is output via copydiag3d so call flumeSendDiag not required
#endif
      ENDIF ! on relevant STASHflags
! ----------------------------------------------------------------------
! STASH item 202: geopotential height      on pressure surfaces
! ----------------------------------------------------------------------
      IF(qhts_p) THEN

        DO k = 1, hts_p_levs

         pressure_pa = hts_p_press(k)*100. ! convert to Pascals
! DEPENDS ON: vert_h_onto_p
         CALL vert_h_onto_p(                                            &
     &        height(1,1,1), row_length,rows, model_levels              &
     &       ,pressure_pa,r_rho_levels, r_theta_levels                  &
     &, p_theta_levels                                                  &
     &       ,theta, exner_theta_levels, exner_rho_levels               &
     &       ,R, g, Lapse,bl_levels                                     &
     &       ,offx, offy, halo_i, halo_j                                &
     &       ,p, interp_order,kappa, p_zero,cp, hts_p(1,1,k) )

        ENDDO ! over output STASH pressure levels

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 16
          ie = 202
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = hts_p_levs
          CALL flumeSendDiag (hts_p,im,is,ie,row_length,rows,levels)
        END IF
#endif

      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 203: temperature              on pressure surfaces
! ----------------------------------------------------------------------
      IF(qt_p) THEN

        DO k = 1, t_p_levs

         pressure_pa = t_p_press(k)*100. ! convert to Pascals
! DEPENDS ON: t_vert_interp_to_p
         CALL T_vert_interp_to_p(                                       &
     &        T, theta, row_length, rows                                &
     &        ,model_levels, pressure_pa, offx, offy,halo_i,halo_j      &
     &        ,p_theta_levels, Lapse, R, g                              &
     &        ,bl_levels                                                &
     &        ,exner_theta_levels                                       &
     &        ,r_theta_levels                                           &
     &        ,kappa, p_zero,  T_p(1,1,k) )

        ENDDO ! over output STASH pressure levels

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 16
          ie = 203
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = t_p_levs
          CALL flumeSendDiag (T_p,im,is,ie,row_length,rows,levels)
        END IF
#endif

      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 204, 256 : RH: relative humidity on pressure surfaces
! ----------------------------------------------------------------------

! STASH item 204: Relative humidity wrt ice
      IF(qRHice_p) THEN

        DO k = 1, wet_levels
!  Find humidity saturation at theta points - store in rh
! DEPENDS ON: qsat
          CALL QSAT(rh(1,1,k),T(1,1,k),p_at_theta(1,1,k),               &
     &              theta_field_size)
!  And convert to relative humidity
          DO j = 1, rows
            DO i = 1, row_length

              rh(i,j,k) = q(i,j,k)/rh(i,j,k)*100.
!  Supersaturation (>100%) can occur with mixed phase scheme but
!  negative humidity is removed from the diagnostic:

              IF(rh(i,j,k) <  0.0) THEN
                rh(i,j,k) = 0.
              END IF
            END DO ! i
          END DO ! j

        END DO ! k wet_levels

!  Interpolate
        DO k = 1, RHice_p_levs

          pressure_pa = RHice_p_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          CALL vert_interp2(                                            &
     &          rh, row_length, rows, wet_levels                        &
     &        , pressure_ex                                             &
     &        , 0, 0, offx, offy                                        &
     &        , exner_theta_levels, interp_order                        &
     &        , RHice_p(1,1,k) )

        ENDDO ! k over output STASH pressure levels

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 16
          ie = 204
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = RHice_p_levs
          CALL flumeSendDiag (RHice_p,im,is,ie,row_length,rows,levels)
        END IF
#endif

      END IF  ! on STASHflag

! STASH item 256: Relative humidity wrt water
      IF(qRHwat_p) THEN

        DO k = 1, wet_levels
!  Find humidity saturation at theta points - store in rh
! DEPENDS ON: qsat_wat
          CALL QSAT_WAT(rh(1,1,k),T(1,1,k),p_at_theta(1,1,k),           &
     &              theta_field_size)
!  And convert to relative humidity
          DO j = 1, rows
            DO i = 1, row_length

              rh(i,j,k) = q(i,j,k)/rh(i,j,k)*100.
!  Supersaturation (>100%) can occur with mixed phase scheme but
!  negative humidity is removed from the diagnostic:

              IF(rh(i,j,k) <  0.0) THEN
                rh(i,j,k) = 0.
              END IF
            END DO ! i
          END DO ! j

        END DO ! k wet_levels

!  Interpolate
        DO k = 1, RHwat_p_levs

          pressure_pa = RHwat_p_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          CALL vert_interp2(                                            &
     &          rh, row_length, rows, wet_levels                        &
     &        , pressure_ex                                             &
     &        , 0, 0, offx, offy                                        &
     &        , exner_theta_levels, interp_order                        &
     &        , RHwat_p(1,1,k) )

        ENDDO ! k over output STASH pressure levels

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 16
          ie = 256
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = RHwat_p_levs
          CALL flumeSendDiag (RHwat_p,im,is,ie,row_length,rows,levels)
        END IF
#endif

      END IF  ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 205: wet bulb potential temperature on pressure surfaces
! ----------------------------------------------------------------------
      IF(qwbpt_p) THEN

        DO k = 1, wbpt_p_levs

         pressure_pa = wbpt_p_press(k)*100. ! convert to Pascals

! Interpolate T onto required pressure level (work1)
! DEPENDS ON: t_vert_interp_to_p
         CALL T_vert_interp_to_p(                                       &
     &        T, theta, row_length, rows                                &
     &        ,model_levels, pressure_pa, offx, offy, halo_i, halo_j    &
     &        ,p_theta_levels, Lapse, R, g                              &
     &        ,bl_levels                                                &
     &        ,exner_theta_levels, r_theta_levels                       &
     &        ,kappa, p_zero, work1 )

! Interpolate q onto required pressure level (work2)
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          CALL vert_interp2(                                            &
     &         q, row_length, rows, wet_levels                          &
     &        ,pressure_ex                                              &
     &        ,halo_i, halo_j, offx, offy                               &
     &        ,exner_theta_levels, interp_order                         &
     &        ,work2 )

! Generate pressure array for required pressure level (work3)
         DO j = 1, rows
           DO i = 1, row_length
             work3(i,j) = pressure_pa
           END DO
         END DO
! DEPENDS ON: thetaw
         CALL Thetaw(                                                   &
     &        theta_field_size,work1,work2,work3,l_potential,           &
                                                               ! in
     &        wbpt_p(1,1,k))                                  ! out

        ENDDO ! over output STASH pressure levels

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 16
          ie = 205
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          levels = wbpt_p_levs
          CALL flumeSendDiag (wbpt_p,im,is,ie,row_length,rows,levels)
        END IF
#endif


      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 222: mean sea level pressure
! ----------------------------------------------------------------------
      IF(qp_MSL) THEN

! DEPENDS ON: calc_pmsl
         CALL Calc_PMSL(theta, exner_theta_levels, p                    &
     &                 ,r_theta_levels, r_rho_levels, eta_theta_levels  &
     &                 ,g, R, Lapse, earth_radius                       &
     &                 ,row_length, rows, model_levels, bl_levels       &
     &                 ,offx, offy, halo_i, halo_j                      &
     &                 , p_MSL, p_star                                  &
     &                 ,mype, nproc, nproc_x, nproc_y                   &
     &                 ,g_datastart                                     &
     &                 ,neighbour, at_extremity                         &
     &                 ,gc_all_proc_group, model_domain                 &
     &                 ,delta_lambda, delta_phi                         &
     &                 ,Cp, npmsl_height, sec_theta_latitude            &
     &                 ,global_row_length, global_rows)

#if defined(FLUME)
!       FLUME-STASH 
        IF (Flume_run) THEN
          im = 1
          is = 16
          ie = 222
        ! Keep the current active model, section, item ID 
        ! in case they are required by the getvarid function
          model=im
          section=is
          item=ie
          CALL flumeSendDiag (p_msl,im,is,ie,row_length,rows)
        END IF
#endif

      ENDIF ! on STASHflag

! ---------------------------------------------------------------------
! Check error condition
      IF(ErrorStatus >  0) THEN
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

      RETURN
      END SUBROUTINE Phy_diag
#endif
