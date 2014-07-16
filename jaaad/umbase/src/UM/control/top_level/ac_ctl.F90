#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine AC_CTL   -----------------------------------------------
!LL
!LL Level 2 control routine
!LL version for CRAY YMP
!LL
!LL programming standard : unified model documentation paper No 3
!LL
!LL Logical components covered : P1
!LL
!LL Project task : P0
!LLEND -----------------------------------------------------------------
!*L Arguments
      SUBROUTINE AC_CTL(INT18,P_FIELDDA,Q_LEVELSDA,                     &
     & model_levelsda, theta_star_inc_halo,                             &
     & q_star, qcl_star, qcf_star,                                      &
     & OBS_FLAG,OBS,obs_flag_len,obs_len,                               &
     & p, p_theta_levels, exner_theta_levels, r_theta_levels,           &
     & fv_cos_theta_latitude,                                           &
     & cf_area, cf_bulk, cf_liquid, cf_frozen,                          &
     & pstar, ntml, cumulus,                                            &
     & STASHwork,                                                       &
#include "argduma.h"
#include "argsts.h"
#include "argppx.h"
     &                  l_mixing_ratio, ICODE, CMESSAGE)

Use ac_diagnostics_mod

      IMPLICIT NONE

#include "cmaxsize.h"
#include "csubmodl.h"
#include "parvars.h"
#include "typsize.h"
#include "typduma.h"
#include "typsts.h"
#include "typcona.h"
#include "ppxlook.h"
#include "nstypes.h"

      INTEGER       INT18        ! Dummy variable for STASH_MAXLEN(18)
      INTEGER       P_FIELDDA,Q_LEVELSDA
      INTEGER model_levelsda             ! copy of model_levels
                                         !    for dynamic allocation
      INTEGER       ICODE        ! Return code : 0 Normal Exit
!                                !             : > 0 Error

      REAL theta_star_inc_halo(1-offx:row_length+offx,                  &
                                                       ! latest
     &               1-offy:rows+offy, model_levelsda) ! theta
      REAL q_star(P_FIELDDA,Q_LEVELSDA)            ! specific humidity,
      REAL qcl_star(P_FIELDDA,Q_LEVELSDA)          ! cloud liquid,
      REAL qcf_star(P_FIELDDA,Q_LEVELSDA)          ! cloud ice content
      CHARACTER*(*) CMESSAGE     ! Error message if return code >0
      INTEGER :: obs_flag_len,obs_len
      INTEGER :: OBS_FLAG(obs_flag_len)
      REAL    :: OBS(obs_len)
      Logical, intent(in)::                                             &
     & l_mixing_ratio            ! Use mixing ratio (if code available)

      Real, Intent (InOut) ::                                           &
     &  p(1-offx:row_length+offx, 1-offy:rows+offy,                     &
     &      model_levels+1)                                             &
     &, p_theta_levels(1-offx:row_length+offx,                          &
     &                   1-offy:rows+offy, model_levels)                &
     &, exner_theta_levels(1-offx:row_length+offx,                      &
     &                     1-offy:rows+offy, model_levels)              &
     &, r_theta_levels(1-offx:row_length+offx,                          &
     &                 1-offy:rows+offy, model_levels)                  &
     &, fv_cos_theta_latitude(1-offx:row_length+offx,                   &
     &                        1-offy:rows+offy)                         &
     &, cf_area(row_length, rows, wet_levels)                           &
     &, cf_bulk(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &           wet_levels)                                            &
     &, cf_liquid(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,     &
     &             wet_levels)                                          &
     &, cf_frozen(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,     &
     &             wet_levels)                                          &
     &, pstar(row_length,rows)                                          &
     &, ntml(row_length,rows)                                           &
     &, cumulus(row_length,rows) 

#include "c_mdi.h"
#include "chsunits.h"
#include "ccontrol.h"
#include "ctime.h"
#include "csizeobs.h"
#include "chistory.h"
#include "ctracera.h"
#include "cruntimc.h"

#include "c_lheat.h"
#include "c_r_cp.h"
#include "acparm.h"
#include "comacp.h"

      REAL lcrcp
      PARAMETER (lcrcp=lc/cp)

! PC2 options are only important if L_pc2 = .true.
! Seek advice from the PC2 team before altering these parameters
! from .true., you will need to have put in place large amounts
! of extra code first.
      LOGICAL,PARAMETER:: L_pc2_cond=.true.  ! Do condensation and
! liquid cloud fraction changes.
      LOGICAL,PARAMETER:: L_pc2_cfl =.true.  ! Do liquid cloud
! fraction changes. This requires that the condensation as a
! result of assimilation is calculated directly.
! Note: One must not try to run with condensation on but liquid cloud
! fraction off with this code.
      LOGICAL,PARAMETER:: L_pc2_cff =.true.  ! Do ice cloud fraction
! changes. This requires that the ice increment from assimilation is
! calculated directly.

!L Dynamically allocated area for stash processing
      REAL STASHWORK(*)


      INTEGER                                                           &
     &       STASHMACRO_TAG,                                            &
                                       ! STASHmacro tag number
     &       MDI,                                                       &
                                       ! Missing data indicator
     &       K, ERROR                  ! do loop variable/ error

      INTEGER J                        ! DO Loop Variable.
      INTEGER I, IIND                  ! temporary scalars
      INTEGER im_index                 ! Internal model index
      INTEGER Q_LEVELS, P_FIELD          ! copies of Q_LEVELSDA,
                                         ! P_FIELDDA to avoid edits
      INTEGER rhc_row_length, rhc_rows   ! rhcrit dimensions
      INTEGER int_p                      ! field pointers
      INTEGER j_start, j_end             ! start and end indices
      INTEGER i_start, i_end             ! for processors
      INTEGER ji                         ! 2D array index for halo
                                         ! i/j variables
      INTEGER levels_per_level, large_levels ! needed by ls_arcld
      REAL RHCPT (1,1,Q_LEVELSDA)        ! rhcrit array
      REAL p_layer_centres(row_length, rows, 0:model_levelsda)
!          pressure at layer centres. Same as p_theta_levels
!          except bottom level = pstar, and at top = 0
      REAL p_layer_boundaries(row_length, rows, 0:model_levelsda)
!                pressure at layer boundaries. Same as p except at
!                bottom level = pstar, and at top = 0.
      REAL exner_layer_centres(row_length, rows, model_levelsda)
!          exner at layer centres. Same as exner_theta_levels
      REAL WORK(P_FIELDDA,Q_LEVELSDA)    ! array for large-scale
                                         ! latent heating
                                         ! or ls_cld dummy output
      REAL WORK2(P_FIELDDA,Q_LEVELSDA)   ! convective heating
                                         ! or ls_cld dummy output
      REAL theta_star(P_FIELDDA,model_levelsda) ! non-halo theta
      REAL bulk_cloud_nohalo(row_length,rows,q_levelsda)

      REAL DUMMY_FIELD(P_FIELDDA)

      ! These variables are for PC2 calculations
      REAL, DIMENSION(:,:,:), ALLOCATABLE::                             &
     &           q_work                                                 &
                                     ! Vapour work array (kg kg-1)
     &,          qcl_work                                               &
                                     ! Liquid work array (kg kg-1)
     &,          qcf_work                                               &
                                     ! Ice work array    (kg kg-1)
     &,          t_work                                                 &
                                     ! Temperature/theta work array
!                                    ! please see inline comments (K)
     &,          bulk_cloud_fraction                                    &
                                       ! Bulk cloud fraction (no units)
     &,          cloud_fraction_liquid                                  &
                                       ! Liquid cloud fraction  "
     &,          cloud_fraction_frozen                                  &
                                       ! Ice cloud fraction     "
     &,          delta_q                                                &
                                 ! Change in q to force condensation
     &,          delta_qcl                                              &
                                 ! Change in qcl to force condensation
     &,          delta_qcf                                              &
                                 ! Change in qcf
     &,          delta_t                                                &
                                 ! Change in t to force condensation
     &,          pc2_work        ! A work array for PC2.

! Declare allocatable arrays for passing cloud fractions
! to LS_ACF_Brooks
      Real, DIMENSION (:,:,:), ALLOCATABLE::                            &
     & cf_bulk_nohalo, cf_liquid_nohalo, cf_frozen_nohalo
!
! Start Routine
!
      Q_LEVELS = Q_LEVELSDA
      P_FIELD  = P_FIELDDA
      DO J=1,P_FIELDDA
         DUMMY_FIELD(J) =0.0
      END DO
!L
!L 1.0 Get address for each field from its STASH section/item code
!L     and STASHmacro tag  (searching only on STASHmacro tag)
      MDI            = IMDI
      STASHMACRO_TAG = 30

! Initialise STASHWORK for section 18.
      DO J = 1, INT18
        STASHWORK(J) = RMDI

      END DO

! Create local versions of exner and p on layer centres/boundaries.
! These local arrays have NO halo.
      Do j=1,rows
        Do i=1,row_length
          p_layer_centres(i,j,0) = pstar(i,j)
          p_layer_boundaries(i,j,0) = pstar(i,j)
        End Do
      End Do

      Do k = 1, model_levels - 1
        Do j = 1, rows
          Do i=1,  row_length
            p_layer_boundaries(i,j,k) = p(i,j,k+1) 
            p_layer_centres(i,j,k) = p_theta_levels(i,j,k) 
            exner_layer_centres(i,j,k)=exner_theta_levels(i,j,k) 
          End Do
        End Do
      End Do
      k=model_levels
      Do j = 1, rows
        Do i=1,  row_length
          p_layer_boundaries(i,j,k) = 0.0
          p_layer_centres(i,j,k) = p_theta_levels(i,j,k) 
          exner_layer_centres(i,j,k)=exner_theta_levels(i,j,k) 
        End Do
      End Do

!L 1.5 large scale rainfall rate LSRR stored in ac_diagnostics module

      IF (LSRR(1)  ==  rmdi) THEN
        ICODE    = 4203
        CMESSAGE = "AC_CTL: large scale rainfall rate not available"
        write(6,*)'AC_CTL 4203 ',ICODE,CMESSAGE
      END IF

      IF (ICODE  >   0) GOTO 999

!L 1.6 large scale snowfall rate LSSR stored in ac_diagnostics module

      IF (LSSR(1)  ==  rmdi) THEN
        ICODE    = 4204
        CMESSAGE = "AC_CTL: large scale snowfall rate not available"
        write(6,*)'AC_CTL 4204 ',ICODE,CMESSAGE
      END IF

      IF (ICODE  >   0) GOTO 999
      
      IF (USE_CONV_IN_MOPS) THEN 

!L 1.7 convective rainfall rate CVRR stored in ac_diagnostics module

       IF (CVRR(1)  ==  rmdi) THEN
         ICODE    = 5205
         CMESSAGE = "AC_CTL: convective rainfall rate not available"

       END IF

       IF (ICODE  >   0) GOTO 999

!L 1.8 convective snowfall rate CVSR stored in ac_diagnostics module

       IF (CVSR(1)  ==  rmdi) THEN
         ICODE    = 5206
         CMESSAGE = "AC_CTL: convective snowfall rate not available"

       END IF

       IF (ICODE  >   0) GOTO 999

!L 1.10 convective cloud cover on each model level CONVCC stored in module

       IF (CONVCC(1,1)  ==  rmdi) THEN
         ICODE    = 5212
         CMESSAGE = "AC_CTL: convective cloud amount not available"

       END IF

       IF (ICODE  >   0) GOTO 999

!L 1.11 bulk cloud fraction (liquid+ice) after large scale cloud CF_LSC 

       IF (CF_LSC(1,1)  ==  rmdi) THEN
         ICODE    = 9201
         CMESSAGE = "AC_CTL: cf after large scale cloud not available"
       ENDIF

       IF (ICODE  >   0) GOTO 999

      ENDIF  !(IF USE_CONV_IN_MOPS)
     
      IF( L_LHN ) THEN
!  seek convective heating rate and
!  diagnostics for calculating large-scale latent heating rate

      IF (USE_CONV_IN_MOPS) THEN 
!L 1.13 temperature increment across convection TINC_CVN stored in module

       IF (TINC_CVN(1,1)  ==  rmdi) THEN
         ICODE    = 5181
         CMESSAGE = "AC_CTL: temp incrs across conv'n not available"
       ENDIF

       IF (ICODE  >   0) GOTO 999

      ENDIF  

!L 1.14 temperature increment across large scale precipitation TINC_PPN

      IF (TINC_PPN(1,1)  ==  rmdi) THEN
        ICODE    = 4181
        CMESSAGE = "AC_CTL: temp incrs across ls_ppn not available"
      ENDIF

      IF (ICODE  >   0) GOTO 999

!L 1.15 cloud liquid water after large scale cloud QCL_LSC stored in module

      IF (QCL_LSC(1,1)  ==  rmdi) THEN
        ICODE    = 9206
        CMESSAGE = "AC_CTL: qcl after large scale cloud not available"
      ENDIF

      IF (ICODE  >   0) GOTO 999

!L 1.16 cloud liquid water after advection QCL_ADV stored in ac_diagnostics module

      IF (QCL_ADV(1,1)  ==  rmdi) THEN
        ICODE    = 12254
        CMESSAGE = "AC_CTL: qcl after advection not available"
      ENDIF

      IF (ICODE  >   0) GOTO 999


! 1.25  Calculate 'large-scale' latent heating contributions
!       ----------------------------------------------------
       DO K=1,Q_LEVELS
         DO J=1,P_FIELD
          WORK(J,K) = TINC_PPN(J,K) +                                   &
     &              lcrcp*( QCL_LSC(J,K)- QCL_ADV(J,K)  )
         END DO
       END DO
!  large scale latent heating currently dT/dt in K/timestep
!  as is convective heating - convert both to dtheta/dt in K/s
       IF (USE_CONV_IN_MOPS) THEN
        Do k = 1, q_levels
          Do j = 1, rows
            Do i=1,  row_length
              int_p = (j-1) * row_length + i
              WORK(int_p,k)  =  WORK(int_p,k) /                          &
     &           (exner_layer_centres(i,j,k) * SECS_PER_STEPim(atmos_im))
              WORK2(int_p,k) =  TINC_CVN(int_p,k) /                      &
     &           (exner_layer_centres(i,j,k) * SECS_PER_STEPim(atmos_im))
            End Do
          End Do
        End Do

       ELSE     ! IF USE_CONV_IN_MOP is false

        Do k = 1, q_levels
          Do j = 1, rows
            Do i=1,  row_length
              int_p = (j-1) * row_length + i
              WORK(int_p,k)  =  WORK(int_p,k) /                          &
     &           (exner_layer_centres(i,j,k) * SECS_PER_STEPim(atmos_im))
              WORK2(int_p,k) = 0.0
            End Do
          End Do
        End Do
       ENDIF

      ELSE     !  if LHN not selected
!                 initialise dummy heating rate array to pass to AC
       DO K=1,Q_LEVELS
         DO J=1,P_FIELD
           WORK(J,K) =  0.0
           WORK2(J,K) =  0.0
         END DO
       END DO

      END IF   !  L_LHN

!L----------------------------------------------------------------------
!L 2. --- Section 18 Data Assimilation ------

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('AC      ', 3)

      END IF

      im_index=internal_model_index(A_IM)


!  copy non halo values from theta_star_inc_halo to
!  2-d array theta_star
      Do K = 1, model_levels
        Do J = 1, rows
          DO I=1,  row_length
            int_p = (J-1) * row_length + I
            theta_star(int_p,K) = theta_star_inc_halo(I,J,K)
          End Do
        End Do
      End Do

        If (L_pc2 .and.                                                 &
     &     (L_pc2_cond .or. L_pc2_cfl .or. L_pc2_cff)   ) then
! Input values needed for the PC2 scheme
        Allocate ( q_work  (row_length,rows,q_levelsda) )
        Allocate ( qcl_work(row_length,rows,q_levelsda) )
        Allocate ( qcf_work(row_length,rows,q_levelsda) )
        Allocate ( t_work  (row_length,rows,q_levelsda) )
        Allocate ( pc2_work(row_length,rows,q_levelsda) )
        Do k = 1, q_levelsda
          Do j = 1, rows
            Do i= 1, row_length
              int_p = (j-1) * row_length + i
              ! Copy input values into the work arrays
              q_work  (i,j,k) = q_star  (int_p,k)
              qcl_work(i,j,k) = qcl_star(int_p,k)
              qcf_work(i,j,k) = qcf_star(int_p,k)
              t_work  (i,j,k) = theta_star(int_p,k)
              ! t_work is now theta before AC assimilation
              pc2_work(i,j,k) = exner_layer_centres(i,j,k)
              ! pc2_work is now exner before AC assimilation
            End Do
          End Do
        End Do
      End If  ! L_pc2 .and. (L_pc2_cond.or.L_pc2_cfl.or.L_pc2_cff)

      IF (USE_CONV_IN_MOPS) THEN

! DEPENDS ON: ac
       CALL AC (                                                         &
     &   model_levels, wet_levels, row_length, rows, BL_LEVELS,          &
     &   A_MAX_OBS_SIZE, A_MAX_NO_OBS, theta_field_size,                 &
     &   STEPim(atmos_im), SECS_PER_STEPim(atmos_im),                    &
     &   exner_layer_centres, pstar,                                     &
     &   p_layer_centres(1,1,1),                                         &
     &   theta_star, q_star, qcl_star, qcf_star,                         &
     &   CONVCC, LSRR, LSSR, CVRR, CVSR,                                 &
     &   CF_LSC,WORK2,WORK, RHCRIT, L_eacf,                              &
     &   OBS_FLAG,OBS,                                                   &
     &   STINDEX(1,1,18,im_index),                                       &
     &   STLIST, LEN_STLIST, SI(1,18,im_index), SF(1,18),                &
     &   STASHWORK, STASH_LEVELS,                                        &
     &   NUM_STASH_LEVELS, STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,        &
#include "argppx.h"
     &   ICODE, CMESSAGE)
      ELSE
       CALL AC (                                                         &
     &   model_levels, wet_levels, row_length, rows, BL_LEVELS,          &
     &   A_MAX_OBS_SIZE, A_MAX_NO_OBS, theta_field_size,                 &
     &   STEPim(atmos_im), SECS_PER_STEPim(atmos_im),                    &
     &   exner_layer_centres, pstar,                                     &
     &   p_layer_centres(1,1,1),                                         &
     &   theta_star, q_star, qcl_star, qcf_star,                         &
     &   WORK2, LSRR, LSSR, DUMMY_FIELD, DUMMY_FIELD,                                 &
     &   CF_LSC,WORK2,WORK, RHCRIT, L_eacf,                              &
     &   OBS_FLAG,OBS,                                                   &
     &   STINDEX(1,1,18,im_index),                                       &
     &   STLIST, LEN_STLIST, SI(1,18,im_index), SF(1,18),                &
     &   STASHWORK, STASH_LEVELS,                                        &
     &   NUM_STASH_LEVELS, STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,        &
#include "argppx.h"
     &   ICODE, CMESSAGE)
      ENDIF

      ! For PC2 we need to calculate the forcing values of qT and TL
        If (L_pc2) then

        If (L_pc2_cond .or. L_pc2_cfl .or. L_pc2_cff) then

        Allocate ( delta_q  (row_length,rows,q_levelsda) )
        Allocate ( delta_qcl(row_length,rows,q_levelsda) )
        Allocate ( delta_qcf(row_length,rows,q_levelsda) )
        Allocate ( delta_t  (row_length,rows,q_levelsda) )
        Allocate ( bulk_cloud_fraction(row_length,rows,q_levelsda)  )
        Allocate ( cloud_fraction_liquid(row_length,rows,q_levelsda))
        Allocate ( cloud_fraction_frozen(row_length,rows,q_levelsda))
        Do k = 1, q_levelsda
          Do j = 1, rows
            Do i= 1, row_length

              int_p = (j-1) * row_length + i
              ! Calculate the increments across the AC scheme
              delta_q  (i,j,k) =  q_star  (int_p,k) - q_work  (i,j,k)
              delta_qcl(i,j,k) =  qcl_star(int_p,k) - qcl_work(i,j,k)
              delta_qcf(i,j,k) =  qcf_star(int_p,k) - qcf_work(i,j,k)

              ! t_work is currently theta before AC assimilation (K)
              delta_t  (i,j,k) =  exner_layer_centres(i,j,k)            &
     &                         * (theta_star(int_p,k) - t_work(i,j,k) )
              ! deltat is now change in temperature associated with
              ! AC assimilation minus the adiabatic contribution
              ! associated with the change in pressure. This last
              ! part is dealt with in pc2_pressure_ctl.
              ! Pc2_work currently contains exner before AC assim
              t_work  (i,j,k) =  pc2_work(i,j,k) * t_work(i,j,k)
              ! t_work now contains temperature before AC assimilation

              pc2_work(i,j,k) =  0.0
              ! pc2_work now contains an array of zeros

              ! Now copy cloud fraction information from the prognostics
              ! arrays
              bulk_cloud_fraction  (i,j,k) = cf_bulk(i,j,k) 
              cloud_fraction_liquid(i,j,k) = cf_liquid(i,j,k) 
              cloud_fraction_frozen(i,j,k) = cf_frozen(i,j,k) 

            End Do
          End Do
        End Do

        ! Now force condensation and cloud fraction updates.
        ! PC2 is not set up for prognostic area_cloud_fraction.

! DEPENDS ON: pc2_assim
        Call pc2_assim( q_levelsda, row_length, rows                    &
     &,                 SECS_PER_STEPim(atmos_im), l_pc2_cfl, l_pc2_cff &
     &,                 l_mixing_ratio                                  &
     &,                 t_work, bulk_cloud_fraction                     &
     &,                 cloud_fraction_liquid, cloud_fraction_frozen    &
     &,                 q_work, qcl_work, qcf_work                      &
     &,              p_layer_centres(1:row_length,1:rows,1:q_levelsda)  &
     &,                 delta_t, delta_q                                &
     &,                 delta_qcl, delta_qcf, pc2_work                  &
     &                  )

        ! q_work, qcl_work and t_work have now been updated by the
        ! forcing and condensation terms together.
        ! cloud_fraction_liquid (and _frozen and _bulk) has been
        ! updated by the condensation. qcf_work and p_work are not
        ! updated.

        ! Now copy q_work, qcl_work and t_work variables back
        ! into the inout variables if we haven't a confident
        ! direct estimate of condensation from another source.
        If (L_pc2_cond) then
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i= 1, row_length
                int_p = (j-1) * row_length + i
                q_star  (int_p,k) = q_work  (i,j,k)
                qcl_star(int_p,k) = qcl_work(i,j,k)
                ! Remember the inout variable is theta, not temp.
                theta_star(int_p,k) = t_work(i,j,k)                     &
     &                              / exner_layer_centres(i,j,k)
              End Do
            End Do
          End Do
        End If  ! L_pc2_cond

        ! Update the D1 cloud fractions
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i= 1, row_length
              cf_bulk(i,j,k) = bulk_cloud_fraction(i,j,k) 
              cf_liquid(i,j,k) = cloud_fraction_liquid(i,j,k) 
              cf_frozen(i,j,k) = cloud_fraction_frozen(i,j,k) 
            End Do
          End Do
        End Do

          If (L_CLD_AREA) Then
            If (L_ACF_Brooks) Then
              Allocate ( cf_bulk_nohalo(row_length,rows,wet_levels) )
              Allocate ( cf_liquid_nohalo(row_length,rows,wet_levels) )
              Allocate ( cf_frozen_nohalo(row_length,rows,wet_levels) )

              i_start = 1
              i_end = row_length
              j_start = 1
              j_end = rows
              If (model_domain  ==  mt_lam) Then
                If(at_extremity(PSouth)) j_start = 2
                If(at_extremity(PNorth)) j_end = rows-1
                If(at_extremity(PWest)) i_start = 2
                If(at_extremity(PEast)) i_end = row_length-1
              End If
              If (model_domain  ==  mt_cyclic_lam) Then
                If (at_extremity(PSouth)) j_start = 2
                If (at_extremity(PNorth)) j_end = rows-1
              End If

! Place bulk, liquid and frozen cloud fractions in halo-free arrays
              Do k = 1, wet_levels
                Do j = j_start, j_end
                  Do i = i_start, i_end
                    ji = i+halo_i + (j+halo_j-1) * (row_length+2*halo_i)
                    cf_bulk_nohalo(i,j,k)  = cf_bulk(i,j,k)
                    cf_liquid_nohalo(i,j,k)= cf_liquid(i,j,k)
                    cf_frozen_nohalo(i,j,k)= cf_frozen(i,j,k)
                  End Do
                End Do
              End Do

! DEPENDS ON: ls_acf_brooks
              Call LS_ACF_Brooks (                                      &
     &             halo_i, halo_j, Offx, Offy                           &
     &,            row_length, rows, model_levels, wet_levels           &
     &,            r_theta_levels, delta_lambda, delta_phi              &
     &,            FV_cos_theta_latitude                                &
     &,            cf_bulk_nohalo, cf_liquid_nohalo                     &
     &,            cf_frozen_nohalo, cumulus                            &
     &,            CF_AREA )

              Deallocate ( cf_bulk_nohalo )
              Deallocate ( cf_liquid_nohalo )
              Deallocate ( cf_frozen_nohalo )

            End If ! L_ACF_Brooks
          End If ! L_cld_area

        ! Now deallocate the arrays
        Deallocate(q_work  )
        Deallocate(qcl_work)
        Deallocate(qcf_work)
        Deallocate(t_work  )
        Deallocate(pc2_work)
        Deallocate(delta_q  )
        Deallocate(delta_qcl)
        Deallocate(delta_qcf)
        Deallocate(delta_t  )
        Deallocate(bulk_cloud_fraction  )
        Deallocate(cloud_fraction_liquid)
        Deallocate(cloud_fraction_frozen)

        End If !  L_pc2_cond.or.L_pc2_cfl.or.L_pc2_cff

      Else  ! L_pc2

!  2.1 call cloud scheme to 'rebalance' thermodynamic fields
!  calculate Tl and qt
         Do k=1,q_levels
           Do j = 1, rows
             Do i=1,  row_length
               int_p = (j-1) * row_length + i
               theta_star(int_p,k) = theta_star(int_p,k)*               &
     &             exner_layer_centres(i,j,k) - lcrcp*qcl_star(int_p,k)
               q_star(int_p,k) = q_star(int_p,k) + qcl_star(int_p,k)
             End Do
           End Do
         End Do
! set up some arguments for ls_cld

      rhc_row_length = 1
      rhc_rows       = 1

        Do k = 1, q_levelsda
          Do j = 1, rows
            Do i= 1, row_length
              bulk_cloud_nohalo(i,j,k) = cf_bulk(i,j,k) 
            End Do
          End Do
        End Do
      DO K = 1, Q_LEVELS
        RHCPT (1,1,K) = RHCRIT(K)
      END DO
! Determine number of sublevels for vertical gradient area cloud
! Want an odd number of sublevels per level: 3 is hardwired in do loops
        levels_per_level = 3
        large_levels = ((q_levels - 2)*levels_per_level) + 2

! DEPENDS ON: ls_arcld
        CALL ls_arcld( p_layer_centres, RHCPT, p_layer_boundaries,      &
     &                 model_levels, q_levels, row_length, rows,        &
     &                 rhc_row_length, rhc_rows, bl_levels,             &
     &                 cloud_fraction_method,overlap_ice_liquid,        &
     &                 ice_fraction_method,ctt_weight,t_weight,         &
     &                 qsat_fixed,sub_cld,                              &
     &                 levels_per_level, large_levels,                  &
     &                 L_cld_area,L_ACF_Cusack,L_ACF_Brooks,L_eacf,     &
     &                 halo_i, halo_j, offx, offy,                      &
     &                 delta_lambda, delta_phi,                         &
     &                 r_theta_levels, FV_cos_theta_latitude,           &
     &                 ntml, cumulus, l_mixing_ratio,                   &
     &                 qcf_star, theta_star, q_star, qcl_star,          &
     &                 CF_AREA(:,:,1), bulk_cloud_nohalo,               &
     &                 WORK2, WORK,                                     &
     &                 ERROR,mype)

        Do k = 1, wet_levels
          Do j = 1, rows
            Do i= 1, row_length
              cf_bulk(i,j,k) = bulk_cloud_nohalo(i,j,k) 
            End Do
          End Do
        End Do

! "1.24.4"  convert t back to theta
         Do k=1,q_levels
           Do j = 1, rows
             Do i=1,  row_length
               int_p = (j-1) * row_length + i
               theta_star(int_p,k) = theta_star(int_p,k)/               &
     &                         exner_layer_centres(i,j,k)
             End Do
           End Do
         End Do

      End If  ! L_pc2

!  copy back theta_star into theta_star_inc_halo
      Do K = 1, model_levels
        Do J = 1, rows
          DO I=1,  row_length
            int_p = (J-1) * row_length + I
            theta_star_inc_halo(I,J,K) = theta_star(int_p,K)
          End Do
        End Do
      End Do


        IF (LTIMER) THEN
! DEPENDS ON: timer
          CALL TIMER('AC      ', 4)

        END IF

!! DEPENDS ON: stash
!        CALL STASH(a_sm, a_im, 18, STASHWORK,                           &
!#include "argdumga.h"
!#include "argsts.h"
!#include "argppx.h"
!     &                                      ICODE, CMESSAGE)
!
!        IF (ICODE >  0) GOTO 999

 999  CONTINUE
      RETURN
      END SUBROUTINE AC_CTL
#endif
