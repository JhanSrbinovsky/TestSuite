#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Subroutine to do boundary layer mixing of UKCA tracers,
!          simultaneously add surface emissions fields onto tracer
!          fields. Also, adds aircraft and lightning emissions.
!          Adapted from version supplied by Olaf Morgenstern.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from UKCA_MAIN1.
!
! Current code owner: Colin Johnson/Olaf Morgenstern
!                     Fiona O'Connor
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
      SUBROUTINE UKCA_EMISSION_CTL(                                    &
          n_chem_tracers, n_emissions,                                 &
          timestep, em_spec,                                           &
          f3_at_u, r_rho_levels, r_theta_levels, sin_theta_latitude,   &
          FV_cos_theta_latitude,                                       &
          true_longitude, delta_lambda, delta_phi,                     &
          iyear, imonth, iday,                                         &
          tropopause_height,                                           &
          ls_mask, conv_cloud_base, conv_cloud_top,                    &
          theta, q, qcl, qcf,                                          &
          exner_rho_levels, rho_r2,                                    &
          p_layer_boundaries, p_theta_levels,                          &
          emissions, aircraftems,                                      &
          SO2emiss_3D,                                                 &
          z_half, alpha_cdx, ml_depth,                                 &
          rhokh_mix, rho_aresist, aresist, resist_b,                   &
          dtrdz_charney_grid, kent, we_lim,                            &
          t_frac, zrzi, kent_dsc,                                      &
          we_lim_dsc, t_frac_dsc,                                      &
          zrzi_dsc, zhsc, rb_dust_ndivs,                               &
          ch4_wetl_emiss, tracers,                                     &
          semiss_budget, light_budget,                                 &
          n_boundary_vals, lbc_spec, lbc_mmr,                          &
          mass)

      USE ASAD_MOD,          ONLY: advt
      USE UKCA_CONSTANTS
      IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "parvars.h"
#include "typsize.h"
#include "c_dust_ndiv.h"

      INTEGER, INTENT(IN) :: n_chem_tracers
      INTEGER, INTENT(IN) :: n_emissions
      INTEGER, INTENT(IN) :: n_boundary_vals   ! No species with  b.c's
      INTEGER, INTENT(IN) :: conv_cloud_base(1:row_length, 1:rows)
      INTEGER, INTENT(IN) :: conv_cloud_top(1:row_length, 1:rows)
      INTEGER, INTENT(IN) :: kent(1:row_length, 1:rows)
      INTEGER, INTENT(IN) :: kent_dsc(1:row_length, 1:rows)
      INTEGER, INTENT(IN) :: iyear
      INTEGER, INTENT(IN) :: imonth
      INTEGER, INTENT(IN) :: iday

      LOGICAL, INTENT(IN) :: ls_mask(1:row_length, 1:rows)

      REAL, INTENT(IN) :: timestep
      REAL, INTENT(IN) :: emissions(row_length,rows,n_emissions)
      REAL, INTENT(IN) :: aircraftems(row_length,rows,model_levels)
      REAL, INTENT(IN) :: SO2emiss_3D(row_length,rows,model_levels)
      REAL, INTENT(IN) :: lbc_mmr(n_boundary_vals)
      REAL, INTENT(IN) :: f3_at_u(1:row_length,1:rows)
      REAL, INTENT(IN) :: r_rho_levels(1:row_length, 1:rows,            &
                             1:model_levels)        ! ht of rho levs
      REAL, INTENT(IN) :: r_theta_levels(1:row_length, 1:rows,          &
                             0:model_levels)        ! ht of theta levs
      REAL, INTENT(IN) :: sin_theta_latitude(1:row_length,1:rows)
      REAL, INTENT(IN) :: FV_cos_theta_latitude(1:row_length,1:rows)
      REAL, INTENT(IN) :: true_longitude(row_length, rows)
      REAL, INTENT(IN) :: delta_lambda
      REAL, INTENT(IN) :: delta_phi
      REAL, INTENT(IN) :: theta(1:row_length,1:rows,model_levels)
      REAL, INTENT(IN) :: q(1:row_length,1:rows,1:model_levels)
      REAL, INTENT(IN) :: qcl(1:row_length,1:rows,1:model_levels)
      REAL, INTENT(IN) :: qcf(1:row_length,1:rows,1:model_levels)
      REAL, INTENT(IN) :: exner_rho_levels(1:row_length,1:rows,         &
                                             1:model_levels+1)
      REAL, INTENT(IN) :: rho_r2(1:row_length,1:rows,model_levels)
      REAL, INTENT(IN) :: p_layer_boundaries(1:row_length,1:rows,       &
                                               0:model_levels)
      REAL, INTENT(IN) :: p_theta_levels(1:row_length,1:rows,           &
                                           1:model_levels)
      REAL, INTENT(IN) :: z_half(1:row_length,1:rows,1:bl_levels)
      REAL, INTENT(IN) :: alpha_cdx(1:bl_levels)
      REAL, INTENT(IN) :: ml_depth(1:row_length,1:rows)
      REAL, INTENT(IN) :: tropopause_height(row_length, rows)
      REAL, INTENT(IN) :: rhokh_mix(1:row_length,1:rows,1:bl_levels)
      REAL, INTENT(IN) :: rho_aresist(1:row_length,1:rows)
      REAL, INTENT(IN) :: aresist(1:row_length,1:rows)
      REAL, INTENT(IN) :: resist_b(1:row_length,1:rows)
      REAL, INTENT(IN) :: dtrdz_charney_grid(1:row_length,1:rows,       &
                                               1:bl_levels)
      REAL, INTENT(IN) :: we_lim(1:row_length,1:rows,1:3)
      REAL, INTENT(IN) :: t_frac(1:row_length,1:rows,1:3)
      REAL, INTENT(IN) :: zrzi(1:row_length,1:rows,1:3)
      REAL, INTENT(IN) :: we_lim_dsc(1:row_length,1:rows,1:3)
      REAL, INTENT(IN) :: t_frac_dsc(1:row_length,1:rows,1:3)
      REAL, INTENT(IN) :: zrzi_dsc(1:row_length,1:rows,1:3)
      REAL, INTENT(IN) :: zhsc(1:row_length,1:rows)
      REAL, INTENT(IN) :: rb_dust_ndivs(1:row_length,1:rows,1:ndiv)
      REAL, INTENT(IN) :: mass(row_length, rows, model_levels)

      CHARACTER*10, INTENT(IN) :: em_spec(n_emissions)
      CHARACTER*10, INTENT(IN) :: lbc_spec(n_boundary_vals)

! Budget of surface emissions
      REAL, INTENT(OUT) :: semiss_budget(row_length,rows,model_levels)
! Budget of lightning emissions
      REAL, INTENT(OUT) :: light_budget(row_length,rows,model_levels)
! Methane wetland emissions
      REAL,INTENT(INOUT) :: ch4_wetl_emiss(1:row_length,1:rows)
! Tracer MMRs
      REAL,INTENT(INOUT) :: tracers(row_length,rows,model_levels,       &
                                      n_chem_tracers)

! Local variables

      INTEGER, SAVE :: inox                 ! Index for NO/NOx tracer
      INTEGER, SAVE :: iso2                 ! Index for SO2 tracer
      INTEGER :: j                          ! Loop variable
      INTEGER :: k                          ! Loop variable
      INTEGER :: l                          ! Loop variable
      INTEGER :: errorstatus                ! Error code from tr_mix
      INTEGER :: jso2_high                  ! HL SO2 emissn index
      INTEGER :: ils_mask(row_length,rows)  ! Land/sea mask (1/0)
      INTEGER :: em_count                   ! counter

      REAL :: res_factor(row_length,rows)                      ! dry
      REAL :: tr_flux(row_length,rows,bl_levels,n_chem_tracers)! trac
      REAL :: surf_dep_flux(row_length,rows,n_chem_tracers)    ! surf
      REAL :: em_field(row_length,rows,n_chem_tracers)         ! surf
      REAL :: conv_aircraftems(row_length,rows,model_levels)   ! airc
      REAL :: lightningems(row_length,rows,model_levels)       ! nox ligh  
      REAL :: conv_SO2emiss_3D(row_length,rows,model_levels)   ! 3D SO2
    
      REAL, SAVE, ALLOCATABLE :: surf_area(:,:)                 ! gridbox area
      REAL, SAVE, ALLOCATABLE :: theta_latitude(:,:)
      REAL, SAVE, ALLOCATABLE :: molmass(:)
      REAL, SAVE, ALLOCATABLE :: lbc_molmass(:)
      LOGICAL, SAVE :: firstcall = .TRUE.

      CHARACTER*72  :: cmessage                          ! Error message

!     Initialise variables

      res_factor(:,:)  = 0.0      ! Set dry deposition to zero
      em_field (:,:,:) = 0.0      ! Initial ems field for each tracer

      IF (firstcall) THEN
        inox             = -99      ! Initial index for nox tracer
        iso2             = -99      ! Initial index for so2 tracer
        jso2_high        = -99      ! Initial index for SO2 HL emissions

        ALLOCATE(theta_latitude(row_length,rows))
        theta_latitude = asin(sin_theta_latitude)

!       Calculate the gridbox surface area

        ALLOCATE(surf_area(row_length,rows))
        DO k = 1, rows 
          DO j = 1, row_length 
            surf_area(j,k) = r_theta_levels(j,k,0)                   & 
                           * r_theta_levels(j,k,0)                   & 
                           * delta_lambda * delta_phi                & 
                           * FV_cos_theta_latitude(j,k) 
          ENDDO 
        ENDDO 
        
!       Find index for SO2 and NOx tracer

        DO k = 1,jpctr
          SELECT CASE (advt(k))
            CASE('NOx       ')
              inox = k
           CASE('NO        ')
              inox = k
           CASE('SO2       ')
              iso2 = k
          END SELECT
        ENDDO

        IF (inox == -99) then
          cmessage = 'Did not find NO or NOx tracer'
! DEPENDS ON: ereport
          CALL EREPORT('UKCA_EMISSION_CTL',1,cmessage)
        ENDIF

        IF (L_ukca_aerchem .AND. iso2 == -99) THEN
          cmessage = 'Did not find SO2 tracer'
! DEPENDS ON: ereport
          CALL EREPORT('UKCA_EMISSION_CTL',1,cmessage)
        ENDIF
         
        IF (L_ukca_strat .OR. L_ukca_stratcfc .OR.                     &
          L_ukca_strattrop) THEN
          ALLOCATE(molmass(n_emissions))
          ALLOCATE(lbc_molmass(n_boundary_vals))

          molmass = 0.
          WHERE (em_spec == 'NOx       ') molmass = m_no2
          WHERE (em_spec == 'NO2       ') molmass = m_no2
          WHERE (em_spec == 'NO        ') molmass = m_no
          WHERE (em_spec == 'CH4       ') molmass = m_ch4
          WHERE (em_spec == 'CO        ') molmass = m_co
          WHERE (em_spec == 'HCHO      ') molmass = m_hcho
          WHERE (em_spec == 'C2H6      ') molmass = m_c2h6
          WHERE (em_spec == 'C3H8      ') molmass = m_c3h8
          WHERE (em_spec == 'Me2CO     ') molmass = m_me2co
          WHERE (em_spec == 'MeCHO     ') molmass = m_mecho
          WHERE (em_spec == 'SO2       ') molmass = m_so2
          WHERE (em_spec == 'Me2S      ') molmass = m_me2s
          WHERE (em_spec == 'COS       ') molmass = m_ocs
          WHERE (em_spec == 'NH3       ') molmass = m_nh3
          WHERE (em_spec == 'N2O       ') molmass = m_n2o
          WHERE (em_spec == 'CF2Cl2    ') molmass = m_cf2cl2
          WHERE (em_spec == 'CFCl3     ') molmass = m_cfcl3
          WHERE (em_spec == 'MeBr      ') molmass = m_mebr
          WHERE (em_spec == 'CF2ClCFCl2') molmass = m_cf2clcfcl2
          WHERE (em_spec == 'MeCl      ') molmass = m_mecl
          WHERE (em_spec == 'MeCCl3    ') molmass = m_meccl3
          WHERE (em_spec == 'CHF2Cl    ') molmass = m_chf2cl
          WHERE (em_spec == 'CFCl3     ') molmass = m_cfcl3
          WHERE (em_spec == 'MeBr      ') molmass = m_mebr
          WHERE (em_spec == 'CF2ClCFCl2') molmass = m_cf2clcfcl2
          WHERE (em_spec == 'MeCl      ') molmass = m_mecl
          WHERE (em_spec == 'MeCCl3    ') molmass = m_meccl3
          WHERE (em_spec == 'CHF2Cl    ') molmass = m_chf2cl
          WHERE (em_spec == 'CCl4      ') molmass = m_ccl4
          WHERE (em_spec == 'CF2ClBr   ') molmass = m_cf2clbr
          WHERE (em_spec == 'CF3Br     ') molmass = m_cf3br
          WHERE (em_spec == 'C5H8      ') molmass = m_isop
          WHERE (em_spec == 'SO4       ') molmass = m_so4
          WHERE (em_spec == 'H2        ') molmass = m_h2
!          WHERE (em_spec == 'C4H10     ') molmass = m_c4h10
!          WHERE (em_spec == 'C10H16    ') molmass = m_c10h16
!          WHERE (em_spec == 'N2O       ') molmass = m_n2o
!          WHERE (em_spec == 'MeOH      ') molmass = m_meoh
!          WHERE (em_spec == 'Rn        ') molmass = m_rn
!          WHERE (em_spec == 'C2H4      ') molmass = m_c2h4
!          WHERE (em_spec == 'C3H6      ') molmass = m_c3h6

          IF (minval(molmass) == 0.) THEN
            cmessage=' Species missing from molmass list.'
! DEPENDS ON: ereport
            CALL EREPORT('UKCA_EMISSION_CTL',1,cmessage)
          ENDIF

! Set MMRs of N2O, CFCs and halons to specified global constants,
! which may follow a time dependence. Set values of inorganic
! Cl and Br compounds to 0 at lower boundary.
! Adjust this if new stratospheric species are introduced.
! The diagnsotic (surface emission) contains the additional
! amount of tracer in mols added globally, which is negative
! in case of sink gases.

          lbc_molmass = 0.
          WHERE (lbc_spec == 'no_N2O    ') lbc_molmass = m_n2o
          WHERE (lbc_spec == 'no_H2     ') lbc_molmass = m_h2
          WHERE (lbc_spec == 'N2O       ') lbc_molmass = m_n2o
          WHERE (lbc_spec == 'CFCl3     ') lbc_molmass = m_cfcl3
          WHERE (lbc_spec == 'CF2Cl2    ') lbc_molmass = m_cf2cl2
          WHERE (lbc_spec == 'HCl       ') lbc_molmass = m_hcl
          WHERE (lbc_spec == 'HOCl      ') lbc_molmass = m_hocl
          WHERE (lbc_spec == 'ClONO2    ') lbc_molmass = m_clono2
          WHERE (lbc_spec == 'Clx       ') lbc_molmass = m_clo
          WHERE (lbc_spec == 'OClO      ') lbc_molmass = m_oclo
          WHERE (lbc_spec == 'MeBr      ') lbc_molmass = m_mebr
          WHERE (lbc_spec == 'HBr       ') lbc_molmass = m_hbr
          WHERE (lbc_spec == 'HOBr      ') lbc_molmass = m_hobr
          WHERE (lbc_spec == 'BrONO2    ') lbc_molmass = m_brono2
          WHERE (lbc_spec == 'Brx       ') lbc_molmass = m_bro
          WHERE (lbc_spec == 'BrCl      ') lbc_molmass = m_brcl
          WHERE (lbc_spec == 'CF2ClCFCl2') lbc_molmass = m_cf2clcfcl2
          WHERE (lbc_spec == 'MeCl      ') lbc_molmass = m_mecl
          WHERE (lbc_spec == 'MeCCl3    ') lbc_molmass = m_meccl3
          WHERE (lbc_spec == 'CHF2Cl    ') lbc_molmass = m_chf2cl
          WHERE (lbc_spec == 'CCl4      ') lbc_molmass = m_ccl4
          WHERE (lbc_spec == 'CF2ClBr   ') lbc_molmass = m_cf2clbr
          WHERE (lbc_spec == 'CF3Br     ') lbc_molmass = m_cf3br
          WHERE (lbc_spec == 'H2        ') lbc_molmass = m_h2
          WHERE (lbc_spec == 'COS       ') lbc_molmass = m_ocs
          WHERE (lbc_spec == 'Cl        ') lbc_molmass = m_cl
          WHERE (lbc_spec == 'ClO       ') lbc_molmass = m_clo
          WHERE (lbc_spec == 'Cl2O2     ') lbc_molmass = m_cl2o2
          WHERE (lbc_spec == 'Br        ') lbc_molmass = m_br
          WHERE (lbc_spec == 'TOT_Cl    ') lbc_molmass = m_cl
          WHERE (lbc_spec == 'TOT_Br    ') lbc_molmass = m_br
          WHERE (lbc_spec == 'BrO       ') lbc_molmass = m_bro
          WHERE (lbc_spec == 'AGE       ') lbc_molmass = 1.

          IF (minval(lbc_molmass) == 0.) THEN
            cmessage=' Species missing from LBC molmass list.'
! DEPENDS ON: ereport
            CALL EREPORT('UKCA_EMISSION_CTL',1,cmessage)
          ENDIF
        ENDIF      ! L_ukca_strat etc

        firstcall = .FALSE.
      ENDIF      ! firstcall etc

!     Set methane emissions to zero over non-land surfaces

      WHERE(ch4_wetl_emiss < 0.0) ch4_wetl_emiss = 0.0

      em_count = 0   ! set emission counter to 0

      DO k = 1,jpctr              ! loop over tracers

!       Check if tracer has surface emissions and set emission.
!       Otherwise emission field is zero from initialisation.

        DO l=1,n_emissions
          IF (advt(k) == em_spec(l) .AND.                              &
            em_spec(l) == 'NO      ' ) THEN
!          Convert from kg NO2/m2/s to kg NO/m2/s
            em_field(:,:,k) = emissions(:,:,l)*m_no/m_no2
          ELSE IF (advt(k) == em_spec(l)(1:3) .AND.                    &
            em_spec(l) == 'SO2_low ' ) THEN
!          Convert from kg S/m2/s to kg SO2/m2/s
            em_field(:,:,k) = emissions(:,:,l)*m_so2/m_s
          ELSE IF (advt(k) == em_spec(l) .AND.                         &
            em_spec(l) == 'DMS     ' ) THEN
!          Convert from kg S/m2/s to kg DMS/m2/s
            em_field(:,:,k) = emissions(:,:,l)*m_dms/m_s
          ELSE IF (advt(k) == em_spec(l) ) THEN
            em_field(:,:,k) = emissions(:,:,l)
          ENDIF

          IF (L_ukca_budget2 .AND. (advt(k) == em_spec(l)) ) THEN
! calculate emission budget in unit of mols.
            em_count = em_count + 1
            IF (em_count > model_levels) THEN
              cmessage=' Increase size of semiss_budget'
! DEPENDS ON: ereport
              CALL ereport('UKCA_EMISSION_CTL',1,cmessage)
            ENDIF
            semiss_budget(:,:,em_count) = emissions(:,:,l) *           &
                      surf_area * timestep * 1000. / molmass(l)
          ENDIF

        ENDDO       ! l=1,n_emissions

!       Add on wetland methane emissions,
!       converting from ugC/m2/s to kgCH4/m2/s

        IF (advt(k) == 'CH4     ') THEN
          em_field(:,:,k) = em_field(:,:,k) +                          &
                            ch4_wetl_emiss*m_ch4*1.0e-9/m_c
        END IF

        IF (L_ukca_strat .OR. L_ukca_stratcfc .OR.                   &
            L_ukca_strattrop)  THEN
          DO l=1,n_boundary_vals
            IF (advt(k) == lbc_spec(l)) THEN
! Set emissions equal to difference of tracer at surface and intented value,
! scaled with mass in gridbox, area and timestep, to turn it into a
! surface emission rate.
              em_field(:,:,k) = (lbc_mmr(l) - tracers(:,:,1,k))      &
                  * mass(:,:,1) / surf_area / timestep

              em_count = em_count + 1
              IF (em_count > model_levels) THEN
                  cmessage=' Increase size of semiss_budget'
! DEPENDS ON: ereport
                CALL ereport('UKCA_EMISSION_CTL',1,cmessage)
              ENDIF

              semiss_budget(:,:,em_count) = em_field(:,:,k) *          &
                    surf_area * timestep * 1000. / lbc_molmass(l)
            ENDIF
          ENDDO
        ENDIF   ! L_ukca_strat etc

!      Call boundary layer mixing and add surface emissions. 
!      Exclude H2O tracer here if an advected tracer. Note 
!      that in climate model q is not mixed either.

        IF (advt(k) .NE. 'H2O       ') THEN
! DEPENDS ON: tr_mix
            CALL TR_MIX(                                               &
            0, 0, row_length, rows, bl_levels,                         &
            offx, offy, alpha_cd,                                      &
            rhokh_mix(1,1,2), rho_aresist,                             &
            dtrdz_charney_grid, r_rho_levels(:,:,1:bl_levels),         &
            r_theta_levels(:,:,0:bl_levels), timestep,                 &
            tr_flux(:,:,1:bl_levels,k), tracers(:,:,1:bl_levels,k),    &
            em_field(:,:,k), res_factor(:,:),                          &
            surf_dep_flux(:,:,k),                                      &
            kent, we_lim, t_frac, zrzi,                                &
            kent_dsc, we_lim_dsc, t_frac_dsc,                          &
            zrzi_dsc, ml_depth, zhsc, z_half,                          &
            errorstatus, .false.                                       &
            )
            IF (errorstatus > 0) THEN
              cmessage=' Error in TR_MIX'
! DEPENDS ON: ereport
            CALL EREPORT('UKCA_EMISSION_CTL',K,cmessage)
          ENDIF
        ENDIF

      ENDDO                       ! end of loop over tracers

!     Set up integer land/sea mask 
         
      DO l = 1,rows 
        DO k = 1,row_length 
          IF (ls_mask(k,l)) THEN 
            ils_mask(k,l) = 1.0 
          ELSE 
            ils_mask(k,l) = 0.0 
          ENDIF 
        ENDDO 
      ENDDO 
         
!     Diagnose NO2 lightning emissions 
         
! DEPENDS ON: ukca_light_ctl 
      CALL UKCA_LIGHT_CTL(                                             & 
          rows,row_length,model_levels,                                & 
          conv_cloud_base(1:row_length,1:rows),                        & 
          conv_cloud_top(1:row_length,1:rows),                         & 
          ils_mask(1:row_length,1:rows),                               & 
          Recip_Pi_Over_180                                            & 
          *asin(f3_at_u(1:row_length,1:rows)/two_omega),               & 
          surf_area(1:row_length,1:rows),                              &
          r_theta_levels(1:row_length,1:rows,0:model_levels),          & 
          r_rho_levels(1:row_length, 1:rows,1:model_levels),           & 
          p_theta_levels(1:row_length,1:rows,1:model_levels),          & 
          p_layer_boundaries(1:row_length,1:rows,0:model_levels),      & 
          lightningems(1:row_length,1:rows,1:model_levels)) 

!       Convert aircraft emissions from kg NO2/gridbox/s to
!       kg NO2/m2/s (for family chemistry) or kg NO/m2/s (for
!       non-family chemistry)

        IF (L_ukca_family) THEN

!         Convert aircraft emissions from kg NO2/gridbox/s
!         to kg NO2/m2/s

          DO l=1,model_levels
            DO k=1,rows
              DO j=1,row_length
                conv_aircraftems(j,k,l)  = aircraftems(j,k,l)          & 
                                         / surf_area(j,k) 
              ENDDO
            ENDDO
          ENDDO

        ELSE

!         Convert aircraft emissions from kg NO2/gridbox/s
!         to kg NO/m2/s

          DO l=1,model_levels
            DO k=1,rows
              DO j=1,row_length
                conv_aircraftems(j,k,l)  = aircraftems(j,k,l)          & 
                                         * m_no/(surf_area(j,k)*m_no2) 
              ENDDO
            ENDDO
          ENDDO
       ENDIF

!      Setup SO2 emissions

       IF (L_ukca_aerchem) THEN

!        Convert 3-D natural SO2 from kg S/m2/s to kg SO2/m2/s.

         conv_SO2emiss_3D(:,:,:)=SO2emiss_3D(:,:,:)*m_so2/m_s

!        Find index for SO2 emissions, and convert to kg SO2/m2/s.
         DO k=1,n_emissions
           IF (em_spec(k)(1:8) == 'SO2_high') THEN
             em_field(:,:,k) = emissions(:,:,k)*m_so2/m_s
             jso2_high       = k
             EXIT
           ENDIF
         END DO

!        Check if emissions are present

         IF (jso2_high == -99) THEN
            cmessage = 'Did not find High Level SO2 emissions'
            CALL EREPORT('UKCA_EMISSION_CTL',1,cmessage)
         ENDIF

       ENDIF    ! L_ukca_aerchem

!       Add aircraft emissions to NO or NOx tracer

        DO k = 1,model_levels   ! Loop over levels

! DEPENDS ON: trsrce
          CALL TRSRCE(                                                 &
          rows, row_length, 0, 0, 0, 0,                                &
          model_levels, wet_levels,                                    &
          0, 0,                                                        &
          r_rho_levels, r_theta_levels,                                &
          theta, q , qcl , qcf , exner_rho_levels, rho_r2,             &
          tracers(:,:,k,inox), conv_aircraftems(:,:,k), k,             &
          timestep, 1, 1, 0.0)

        ENDDO                  ! End of looping over levels

!       Update tracer fields with NO/NOx lightning emissions 
                 
        IF (L_ukca_family) THEN 
          tracers(1:row_length,1:rows,:,inox) =                        & 
          tracers(1:row_length,1:rows,:,inox) + timestep * lightningems 
        ELSE 
          tracers(1:row_length,1:rows,:,inox) =                        & 
          tracers(1:row_length,1:rows,:,inox) +                        & 
          timestep*lightningems*m_no/m_no2 
        ENDIF 

!       Add 3-D volcanic and high-level anthropogenic emissions 
!       to SO2 tracer

        IF(L_ukca_aerchem) THEN
          DO k = 1,model_levels   ! Loop over levels
            IF (k == SO2_high_level) THEN
              conv_SO2emiss_3D(:,:,k) = conv_SO2emiss_3D(:,:,k) +      &
                                         em_field(:,:,jso2_high)
              EXIT
            ENDIF
          ENDDO
! DEPENDS ON: trsrce
          CALL TRSRCE( rows, row_length, 0, 0, 0, 0,                   &
             model_levels, wet_levels, 0, 0,                           &
             r_rho_levels, r_theta_levels,                             &
             theta, q , qcl , qcf , exner_rho_levels, rho_r2,          &
             tracers(:,:,k,iso2), conv_SO2emiss_3D(:,:,k),             &
             SO2_high_level, timestep, 1, 1, 0.0)
        ENDIF

        IF (L_ukca_budget2) THEN
! Lightning production in mols.
          light_budget = timestep * lightningems * mass * 1000./m_no2
        ENDIF

      IF ((L_ukca_strat .OR. L_ukca_stratcfc .OR.                   &
          L_ukca_strattrop) .AND. iso2 > 0)  THEN
! Perform emission of volcanic SO2 from explosive volcanic eruptions into
!  Stratosphere
! DEPENDS ON: ukca_volcanic_so2
          CALL UKCA_VOLCANIC_SO2(                                       &
                 tracers(1:row_length, 1:rows, :, iso2),                &
                 mass, theta_latitude, true_longitude,                  &
                 delta_phi, delta_lambda,                               &
                 row_length, rows, model_levels,                        &
                 iyear, imonth, iday, timestep,                         &
                 tropopause_height,                                     &
                 r_theta_levels(1:row_length, 1:rows, 1:model_levels))
      ENDIF



      RETURN
      END SUBROUTINE UKCA_EMISSION_CTL
#endif
