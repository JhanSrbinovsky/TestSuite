#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to diagnose ISCCP cloud types
!
! Purpose:
!   This subroutine diagnoses ISCCP cloud types
!
! Method:
!   Argument list mirrors that of R2_SWRAD
!
! Current Owner of Code: M. J. Webb
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- --------------------------------------------------------------------
      SUBROUTINE ISCCP(IERR, I_SEGMENT                                  &
!
     &   , H2O                                                          &
!                       Pressure Fields
     &   , TSTAR, PSTAR, P_LAYER_BOUNDARIES, P_LAYER_CENTRES            &
!                       Temperatures
     &   , TAC                                                          &
!                       Options for treating clouds
     &   , LCA_AREA, LCA_BULK                                           &
!                       Convective Cloud Fields
     &   , CCA, CCCWP, CCB, CCT                                         &
!                       Surface Fields
!                       Solar Fields
     &   , LIT, LIST, SCS                                               &
     &   , TRINDX_FLD                                                   &
!
!                       General Diagnostics
     &   , LW_diag, SW_diag, row_list, col_list                         &
!
!                       Physical Dimensions
     &   , NLIT                                                         &
     &   , N_POINTS, NLEVS, N_LAYER, NCLDS                              &
!          Number of points in segment, number of levels,
!          cloud levels
     &   , NWET, NOZONE                                                 &
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_COLUMN                &
     &   , N_CCA_LEV                                                    &
           ! Variables needed to calculate layer masses
     &   , rho_r2, r_rho_levels, r_theta_levels                         &
     &   , q, qcl, qcf, qcf2, qrain, qgraup                             &
     &   , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio        &
     &   )
!
!
!
      Use swrdiag_mod, Only:                                            &
          StrSWDiag

      Use lwrdiag_mod, Only:                                            &
          StrLWDiag

      Use cv_cntl_mod, Only:                                            &
          lcv_3d_cca

      IMPLICIT NONE
!
!
!
!     COMDECKS INCLUDED
!     UNIT NUMBERS FOR PRINTED OUTPUT
#include "stdio3a.h"
!     ERROR FLAGS
#include "error3a.h"
!     Radiation user interface data (contains IS_NCOL)
#include "rad_com.h"
!
!     DIMENSIONS OF ARRAYS:
      INTEGER IERR   ! error code

      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_FIELD                                                    &
!             FIELD SIZE IN CALLING PROGRAM
     &   , NPD_PROFILE                                                  &
!             SIZE OF ARRAY OF PROFILES
     &   , NPD_LAYER                                                    &
!             ARRAY SIZES FOR LAYERS
     &   , NPD_COLUMN
!             NUMBER OF COLUMNS PER POINT

      INTEGER                                                           &
     &     TOP_HEIGHT                                                   &
!             FLAG TO SPECIFY METHOD FOR ISCCP CLOUD TOP HEIGHT
     &   , OVERLAP
!             FLAG TO SPECIFY METHOD FOR ISCCP CLOUD OVERLAP

!
!
!     ACTUAL SIZES USED:
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_POINTS                                                     &
!             Total number of points, including unlit ones
     &   , NWET                                                         &
!             NUMBER OF WET LEVELS
     &   , NOZONE                                                       &
!             NUMBER OF LEVELS WITH OZONE
     &   , NLEVS                                                        &
!             NUMBER OF ATMOSPHERIC LAYERS
     &   , N_LAYER                                                      &
!             Number of layers seen in the radiation scheme
     &   , NCLDS                                                        &
!             NUMBER OF CLOUDY LEVELS
     &   , N_CCA_LEV
!             NUMBER OF CONVECTIVE CLOUD LEVELS
!
!
!
!
!     GASEOUS MIXING RATIOS
      REAL                                                              &
                !, INTENT(IN)
     &     H2O(NPD_FIELD, NWET)
!             MASS MIXING RATIO OF WATER
!
!     GENERAL ATMOSPHERIC PROPERTIES:
      REAL                                                              &
                !, INTENT(IN)
     &     TSTAR(NPD_FIELD)                                             &
!             SURFACE TEMPERATURES
     &   , PSTAR(NPD_FIELD)                                             &
!             SURFACE PRESSURES
     &   , P_LAYER_BOUNDARIES(NPD_FIELD,0:NLEVS)                        &
!            PRESSURE AT BOUNDARIES OF LAYERS
     &   , P_LAYER_CENTRES(NPD_FIELD,0:NLEVS)                           &
!            PRESSURE AT CENTRES OF LAYERS
     &   , TAC(NPD_FIELD, NLEVS)
!             TEMPERATURES AT CENTRES OF LAYERS
!
!     INCIDENT SOLAR RADIATION:
      INTEGER                                                           &
                !, INTENT(IN)
     &     NLIT                                                         &
!             NUMBER OF LIT POINTS
     &   , LIST(NPD_FIELD)
!             LIST OF LIT POINTS
      REAL                                                              &
                !, INTENT(IN)
     &     SCS                                                          &
!             SCALING OF SOLAR INCIDENT FIELD
     &   , LIT(NPD_FIELD)
!             FRACTION OF TIME POINT IS LIT
!
      REAL                                                              &
                !, INTENT(IN)
     &     LCA_AREA(NPD_FIELD, NCLDS+1/(NCLDS+1))                       &
!             AREA FRACTIONS OF LAYER CLOUDS OUTSIDE CONVECTIVE TOWERS
     &   , LCA_BULK(NPD_FIELD, NCLDS+1/(NCLDS+1))
!             BULK FRACTIONS OF LAYER CLOUDS OUTSIDE CONVECTIVE TOWERS
!
!     PROPERTIES OF CONVECTIVE CLOUDS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     CCB(NPD_FIELD)                                               &
!             BASE OF CONVECTIVE CLOUD
     &   , CCT(NPD_FIELD)
!             TOP OF CONVECTIVE CLOUD
      REAL                                                              &
                !, INTENT(IN)
     &     CCCWP(NPD_FIELD)                                             &
!             WATER PATH OF CONVECTIVE CLOUD
     &   , CCA(NPD_FIELD,N_CCA_LEV)
!             FRACTION OF CONVECTIVE CLOUD
      LOGICAL                                                           &
                !, INTENT(IN)
     &     l_mcr_qcf2                                                   &
                          ! Use second ice category
     &,    l_mcr_qrain                                                  &
                          ! Use prognostic rain
     &,    l_mcr_qgraup                                                 &
                          ! Use graupel
     &,    l_mixing_ratio ! Use mixing ratios in layer mass calculation
!
!                       Level of tropopause
      INTEGER                                                           &
     &     TRINDX_FLD(NPD_FIELD)
!             THE LAYER BOUNDARY OF THE TROPOPAUSE
!     Information for the calculation of layer masses
      Real, intent(in)::                                                &
     &  rho_r2(npd_field,nlevs)                                         &
                                ! Air density*radius of earth**2 / kg m-1
     &, r_rho_levels(npd_field,nlevs)                                   &
                                      ! Height of rho levels / m
     &, r_theta_levels(npd_field,0:nlevs)                               &
                                           ! Height of theta levels / m
     &, q(npd_field,nwet)                                               &
                                ! Water vapour mixing ratio / kg kg-1
     &, qcl(npd_field,nwet)                                             &
                                ! Liquid water mixing ratio / kg kg-1
     &, qcf(npd_field,nwet)                                             &
                                ! Ice mixing ratio / kg kg-1
     &, qcf2(npd_field,nwet)                                            &
                                ! Second ice category mr / kg kg-1
     &, qrain(npd_field,nwet)                                           &
                                ! Rain mixing ratio / kg kg-1
     &, qgraup(npd_field,nwet)  ! Graupel mixing ratio / kg kg-1
!
!     CALCULATED FLUXES:
      REAL                                                              &
                !, INTENT(OUT)
     &     SWSEA(NPD_FIELD)
!             SEA-SURFACE COMPONENTS OF FLUX
!
!     DIAGNOSTICS:
!
!     Definition of the diagnostic structure
!
      Type (StrLWDiag) :: LW_diag
      Type (StrSWDiag) :: SW_diag
!
      Integer, Intent(IN) :: row_list(npd_field)
!                              List of row indices of lit points
      Integer, Intent(IN) :: col_list(npd_field)
!                              List of column indices of lit points
!
!
!     CALCULATED LOCAL DIAGNOSTICS:
      REAL                                                              &
                !, INTENT(OUT)
     &     PFULL_FLD(NPD_FIELD, NCLDS)                                  &
     &   , PHALF_FLD(NPD_FIELD, NCLDS+1)                                &
     &   , QV_FLD(NPD_FIELD, NCLDS)                                     &
     &   , CC_FLD(NPD_FIELD, NCLDS)                                     &
     &   , CONV_FLD(NPD_FIELD, NCLDS)                                   &
     &   , DTAU_S_FLD(NPD_FIELD, NCLDS)                                 &
     &   , DTAU_C_FLD(NPD_FIELD, NCLDS)                                 &
     &   , SKT_FLD(NPD_FIELD)                                           &
     &   , EMSFC_LW_FLD(NPD_FIELD)                                      &
     &   , AT_FLD(NPD_FIELD, NCLDS)                                     &
     &   , DEM_S_FLD(NPD_FIELD, NCLDS)                                  &
     &   , DEM_C_FLD(NPD_FIELD, NCLDS)                                  &
     &   , FQ_ISCCP_FLD(NPD_FIELD, 7, 7)                                &
     &   , MEANALBEDOCLD(NPD_FIELD)                                     &
     &   , MEANTAUCLD(NPD_FIELD)                                        &
     &   , MEANPTOP(NPD_FIELD)                                          &
     &   , TOTALCLDAREA(NPD_FIELD)
!
!     LOCAL VARIABLES.
!
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , J                                                            &
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE

!     GENERAL ATMOSPHERIC PROPERTIES:
      REAL                                                              &
     &     D_MASS(NPD_PROFILE, NPD_LAYER)                               &
!             MASS THICKNESSES OF LAYERS
     &   , P(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             PRESSURE FIELD
     &   , T(NPD_PROFILE, 0: NPD_LAYER)
!             TEMPERATURE FIELD
!
      REAL :: layer_heat_capacity(npd_profile, npd_layer)
!             Specific heat capacity of layer * d_mass
!
      REAL                                                              &
     &     DUMMY,DUMMY1D(1),DUMMY2D(1,1)
!     ADDITIONAL ARGUMENTS
      INTEGER                                                           &
     &     I_SEGMENT


      IF (IS_NCOL  <   10) THEN
       WRITE(0,*) 'WARNING: No. of ISCCP sampling columns less than 10'
       WRITE(0,*) 'This may introduce _systematic_ biases'
        IF (IS_NCOL  <   1) THEN
          WRITE(IU_ERR, '(/A, /A)')                                     &
     &      'ERROR: Number of ISCCP sampling columns not specified'
          IERR=I_ERR_FATAL
          RETURN
        ENDIF
      ENDIF

      OVERLAP=3
      TOP_HEIGHT=1


!
!     SET THE THERMODYNAMIC PROPERTIES OF THE ATMOSPHERE.
! DEPENDS ON: r2_set_thermodynamic
      CALL R2_SET_THERMODYNAMIC(NLIT, NLEVS, N_LAYER, nwet, LIST        &
     &   , .FALSE., .FALSE.                                             &
     &   , PSTAR, DUMMY1D, DUMMY1D, DUMMY1D                             &
     &   , P_LAYER_BOUNDARIES                                           &
     &   , P_LAYER_CENTRES                                              &
     &   , DUMMY2D                                                      &
     &   , DUMMY2D, TAC                                                 &
     &   , rho_r2, r_rho_levels, r_theta_levels                         &
     &   , q, qcl, qcf, qcf2, qrain, qgraup                             &
     &   , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio        &
     &   , P(:,1:), T(:,1:), DUMMY2D, DUMMY1D, DUMMY1D, DUMMY1D, D_MASS &
     &   , layer_heat_capacity                                          &
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER                            &
     &   )
!
      IF (LW_diag%L_isccp_cf .OR.                                       &
     &    LW_diag%L_isccp_cf_tau_0_to_p3 .OR.                           &
     &    LW_diag%L_isccp_cf_tau_p3_to_1p3 .OR.                         &
     &    LW_diag%L_isccp_cf_tau_1p3_to_3p6 .OR.                        &
     &    LW_diag%L_isccp_cf_tau_3p6_to_9p4 .OR.                        &
     &    LW_diag%L_isccp_cf_tau_9p4_to_23 .OR.                         &
     &    LW_diag%L_isccp_cf_tau_23_to_60 .OR.                          &
     &    LW_diag%L_isccp_cf_tau_ge_60) THEN
          IF (.NOT.LW_diag%L_isccp_weights) THEN
            WRITE(IU_ERR, '(/A, /A)')                                   &
     &         '*** ERROR: ISCCP CLOUD FRACTIONS MAY BE DIAGNOSED'      &
     &         , 'ONLY IN CONJUNCTION WITH THE CORREPONDING WEIGHT.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF


!     Set up arrays for ISCCP_CLOUDTYPES_FLD

      DO I=1,NCLDS+1
        DO L=1,N_POINTS
          PHALF_FLD(L,I)=P_LAYER_BOUNDARIES(L,(I-1))
        ENDDO
      ENDDO

       IF (LCV_3D_CCA) THEN
         DO I=1, NCLDS
            DO L=1, N_POINTS
               CONV_FLD(L,I)=CCA(L,I)
            ENDDO
         ENDDO
       ELSE
         DO I=1, NCLDS
            DO L=1, N_POINTS
               IF ( (I <= CCT(L)-1).AND.(I >= CCB(L)) ) THEN
                  CONV_FLD(L,I)=CCA(L,1)
               ELSE
                  CONV_FLD(L,I)=0.0E+00
               ENDIF
            ENDDO
         ENDDO
       ENDIF

      DO I=1,NCLDS
        DO L=1,N_POINTS
          PFULL_FLD(L,I)=P_LAYER_BOUNDARIES(L,(I-1))
          QV_FLD(L,I)=H2O(L,I)
          CC_FLD(L,I)=CONV_FLD(L,I)+(1-CONV_FLD(L,I))*LCA_AREA(L,I)
          DTAU_S_FLD(L,I)=0        ! lit points filled in below
          DTAU_C_FLD(L,I)=0        ! lit points filled in below
          SKT_FLD(L)=TSTAR(L)
          EMSFC_LW_FLD(L)=0.99
          AT_FLD(L,I)=TAC(L,I)
          DEM_S_FLD(L,I)=0         ! lit points filled in below
          DEM_C_FLD(L,I)=0         ! lit points filled in below
        ENDDO
      ENDDO

      DO I=1,7
        DO J=1,7
          DO L=1,N_POINTS
            FQ_ISCCP_FLD(L,I,J)=0
          ENDDO
        ENDDO
      ENDDO
      DO L=1,N_POINTS
        MEANALBEDOCLD(L)=0.0
        MEANTAUCLD(L)=0.0
        MEANPTOP(L)=0.0
        TOTALCLDAREA(L)=0.0
      ENDDO

!     Fill in optical thickness and emissivity for lit points.
!     1.0E-300 is used rather than 0.0 on the IF test to stop very
!     small weights from causing an overflow failure

      DO I=1,NCLDS
        DO L=1,NLIT

          IF (SW_diag%ls_cloud_weight_extinction                        &
     &    (col_list(L), row_list(L),I) >  1.0E-300) THEN
            DTAU_S_FLD(LIST(L),I)=D_MASS(L,NLEVS+1-I)                   &
     &       *SW_diag%ls_cloud_extinction                               &
     &       (col_list(L), row_list(L),I)                               &
     &       /SW_diag%ls_cloud_weight_extinction                        &
     &       (col_list(L), row_list(L),I)
          ENDIF


          IF (SW_diag%cnv_cloud_weight_extinction                       &
     &    (col_list(L), row_list(L),I) >  1.0E-300) THEN
            DTAU_C_FLD(LIST(L),I)=D_MASS(L,NLEVS+1-I)                   &
     &       *SW_diag%cnv_cloud_extinction                              &
     &       (col_list(L), row_list(L),I)                               &
     &       /SW_diag%cnv_cloud_weight_extinction                       &
     &       (col_list(L), row_list(L),I)
          ENDIF

          IF (LW_diag%ls_cloud_weight_absorptivity                      &
     &    (col_list(L), row_list(L),I) >  1.0E-300) THEN
            DEM_S_FLD(LIST(L),I)=1.0-EXP(-1.666*D_MASS(L,NLEVS+1-I)     &
     &       *LW_diag%ls_cloud_absorptivity                             &
     &       (col_list(L), row_list(L),I)                               &
     &       /LW_diag%ls_cloud_weight_absorptivity                      &
     &       (col_list(L), row_list(L),I))
          ENDIF

          IF (LW_diag%cnv_cloud_weight_absorptivity                     &
     &    (col_list(L), row_list(L),I) >  1.0E-300) THEN
            DEM_C_FLD(LIST(L),I)=1.0-EXP(-1.666*D_MASS(L,NLEVS+1-I)     &
     &       *LW_diag%cnv_cloud_absorptivity                            &
     &       (col_list(L), row_list(L),I)                               &
     &       /LW_diag%cnv_cloud_weight_absorptivity                     &
     &       (col_list(L), row_list(L),I))
          ENDIF
        ENDDO
      ENDDO



! DEPENDS ON: isccp_cloudtypes_fld
      CALL ISCCP_CLOUDTYPES_FLD(I_SEGMENT                               &
!          Input
     &   , NPD_FIELD                                                    &
                        ! Number of points in whole field
     &   , N_POINTS                                                     &
                        ! Number of points in segment
     &   , NLIT                                                         &
                        ! Number of lit points
     &   , LIST                                                         &
                        ! Indexes of lit points
     &   , NCLDS                                                        &
                        ! Number of levels
     &   , IS_NCOL                                                      &
                        ! Number of columns for gridbox decomposition
     &   , PFULL_FLD                                                    &
                        ! Pressure on full levels
     &   , PHALF_FLD                                                    &
                        ! Pressure on full levels
     &   , QV_FLD                                                       &
                        ! Water vapour
     &   , CC_FLD                                                       &
                        ! Total cloud amount (cc+(1-cc)*ls)
     &   , CONV_FLD                                                     &
                        ! Convective cloud amount
     &   , DTAU_S_FLD                                                   &
                        ! Large-scale optical thickness
     &   , DTAU_C_FLD                                                   &
                        ! Convective optical thickness
     &   , TOP_HEIGHT                                                   &
                        ! Flag to specify cloud top method
     &   , OVERLAP                                                      &
                        ! Flag to specify overlap method
     &   , SKT_FLD                                                      &
                        ! Surface temperature
     &   , EMSFC_LW_FLD                                                 &
                        ! Surface emissivity
     &   , AT_FLD                                                       &
                        ! Atmospheric temperature
     &   , DEM_S_FLD                                                    &
                        ! Large-scale cloud emissivities
     &   , DEM_C_FLD                                                    &
                        ! Convective cloud emissivities
     &   , TRINDX_FLD                                                   &
                        ! tropopause index
!          Output
     &   , FQ_ISCCP_FLD                                                 &
                        ! Cloud amounts in various ISCCP classes
     &   , MEANALBEDOCLD                                                &
                        ! Weighted mean cloud albedo
     &   , MEANTAUCLD                                                   &
                        ! Weighted mean cloud optical depth
     &   , MEANPTOP                                                     &
                        ! Weighted mean cloud top pressure
     &   , TOTALCLDAREA                                                 &
                        ! Total cloud area
     &   )


!     ISCCP diagnostics

        DO I=1,7
          DO L=1,NLIT

            IF (LW_diag%L_isccp_cf_tau_0_to_p3)                         &
     &           LW_diag%isccp_cf_tau_0_to_p3                           &
     &             (col_list(L),row_list(L),I)                          &
     &             =FQ_ISCCP_FLD(LIST(L),1,I)

            IF (LW_diag%L_isccp_cf_tau_p3_to_1p3)                       &
     &           LW_diag%isccp_cf_tau_p3_to_1p3                         &
     &             (col_list(L),row_list(L),I)                          &
     &             =FQ_ISCCP_FLD(LIST(L),2,I)

            IF (LW_diag%L_isccp_cf_tau_1p3_to_3p6)                      &
     &           LW_diag%isccp_cf_tau_1p3_to_3p6                        &
     &             (col_list(L),row_list(L),I)                          &
     &             =FQ_ISCCP_FLD(LIST(L),3,I)

            IF (LW_diag%L_isccp_cf_tau_3p6_to_9p4)                      &
     &           LW_diag%isccp_cf_tau_3p6_to_9p4                        &
     &             (col_list(L),row_list(L),I)                          &
     &             =FQ_ISCCP_FLD(LIST(L),4,I)

            IF (LW_diag%L_isccp_cf_tau_9p4_to_23)                       &
     &           LW_diag%isccp_cf_tau_9p4_to_23                         &
     &             (col_list(L),row_list(L),I)                          &
     &             =FQ_ISCCP_FLD(LIST(L),5,I)

            IF (LW_diag%L_isccp_cf_tau_23_to_60)                        &
     &           LW_diag%isccp_cf_tau_23_to_60                          &
     &             (col_list(L),row_list(L),I)                          &
     &             =FQ_ISCCP_FLD(LIST(L),6,I)

            IF (LW_diag%L_isccp_cf_tau_ge_60)                           &
     &           LW_diag%isccp_cf_tau_ge_60                             &
     &             (col_list(L),row_list(L),I)                          &
     &             =FQ_ISCCP_FLD(LIST(L),7,I)

            IF (LW_diag%L_isccp_cf) THEN
              LW_diag%isccp_cf(col_list(L),row_list(L),I)=0.0
              DO J=1,7
                LW_diag%isccp_cf(col_list(L),row_list(L),I)=            &
     &          LW_diag%isccp_cf(col_list(L),row_list(L),I)+            &
     &             FQ_ISCCP_FLD(LIST(L),J,I)
              ENDDO
            ENDIF
          ENDDO
        ENDDO

        IF (LW_diag%L_meanalbedocld .AND. LW_diag%L_totalcldarea) THEN
          DO L=1,NLIT
            LW_diag%meanalbedocld(col_list(L), row_list(L))             &
     &             =MEANALBEDOCLD(LIST(L)) * TOTALCLDAREA(LIST(L))
          ENDDO
        ENDIF

        IF (LW_diag%L_meantaucld .AND. LW_diag%L_totalcldarea) THEN
          DO L=1,NLIT
            LW_diag%meantaucld(col_list(L), row_list(L))                &
     &             =MEANTAUCLD(LIST(L)) * TOTALCLDAREA(LIST(L))
          ENDDO
        ENDIF

        IF (LW_diag%L_meanptop .AND. LW_diag%L_totalcldarea) THEN
          DO L=1,NLIT
            LW_diag%meanptop(col_list(L), row_list(L))                  &
     &             =MEANPTOP(LIST(L)) * TOTALCLDAREA(LIST(L))
          ENDDO
        ENDIF

        IF (LW_diag%L_totalcldarea) THEN
          DO L=1,NLIT
            LW_diag%totalcldarea(col_list(L), row_list(L))              &
     &             =TOTALCLDAREA(LIST(L))
          ENDDO
        ENDIF

        IF (LW_diag%L_isccp_weights) THEN
          DO L=1,NLIT
            LW_diag%isccp_weights(col_list(L), row_list(L))=1.0
          ENDDO
        ENDIF

!
!     sychronise PEs
!     call barrier()
!     Force core dump
!      if (node == -1) DUMMY=sqrt(real((-nlit-1)))
!     Force other PE's stop here
!      call barrier()

      RETURN
      END SUBROUTINE ISCCP
#endif
