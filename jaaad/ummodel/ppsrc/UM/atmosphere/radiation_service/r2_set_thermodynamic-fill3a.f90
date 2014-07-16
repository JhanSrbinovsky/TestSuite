


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set thermodynamic properties
!
! Purpose:
!   Pressures, temperatures at the centres and edges of layers
!   and the masses in layers are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: James Manners
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_THERMODYNAMIC(N_PROFILE, NLEVS, N_LAYER,q_levels&
     &   , I_GATHER, L_EXTRA_TOP, L_BOUNDARY_TEMPERATURE                &
     &   , PSTAR, TSTAR, TSTAR_SOLID, TSTAR_SEA                         &
     &   , P_LAYER_BOUNDARIES                                           &
     &   , P_LAYER_CENTRES                                              &
     &   , HEIGHT_THETA                                                 &
     &   , HEIGHT_RHO                                                   &
     &   , TAC                                                          &
           ! IN. Height and moisture information for the calculation
           ! of layer masses
     &   , rho_r2, r_rho_levels, r_theta_levels                         &
     &   , q, qcl, qcf, qcf2, qrain, qgraup                             &
     &   , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio        &
           ! OUT
     &   , P, T, T_BDY, T_SURFACE, T_SOLID, T_SEA, D_MASS               &
     &   , layer_heat_capacity                                          &
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER                            &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     INCLUDED COMDECKS
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------
! C_PERMA start

      ! Specific heat capacity of water vapour (J/kg/K)
      REAL,PARAMETER:: HCAPV=1850.0

      ! Specific heat capacity of water (J/kg/K)
      REAL,PARAMETER:: HCAPW=4180.0

      ! Specific heat capacity of ice (J/kg/K)
      REAL,PARAMETER:: HCAPI=2100.0

      ! Density of ice (kg/m3)
      REAL,PARAMETER:: RHO_ICE=917

      ! Rate of change of ice potential with temperature
      ! RHO_ICE*LF/ZERODEGC*1/(RHO_WATER*G) (m/K)
      REAL,PARAMETER:: DPSIDT=114.3

! C_PERMA end
!
!     DUMMY ARGUMENTS.
!     SIZES OF ARRAYS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_FIELD                                                    &
!             SIZE OF ARRAY FROM UM
     &   , NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     SIZES USED:
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , NLEVS                                                        &
!             Number of levels in the main model
     &   , q_levels                                                     &
                     ! Number of moist model levels
     &   , N_LAYER                                                      &
!             Number of layers seen by radiation
!             NUMBER OF LEVELS
     &   , I_GATHER(NPD_FIELD)
!             LIST OF POINTS TO GATHER
!
      REAL                                                              &
                !, INTENT(IN)
     &     PSTAR(NPD_FIELD)                                             &
!             SURFACE PRESSURES
     &   , TSTAR(NPD_FIELD)                                             &
!             SURFACE TEMPERSTURES
     &   , TSTAR_SOLID(NPD_FIELD)                                       &
!             SOLID SURFACE TEMPERATURES
     &   , TSTAR_SEA(NPD_FIELD)                                         &
!             OPEN SEA SURFACE TEMPERATURES
     &   ,P_LAYER_BOUNDARIES(NPD_FIELD,0:NLEVS)                         &
!             PRESSURE AT EDGES OF LAYERS
     &   ,P_LAYER_CENTRES(NPD_FIELD,0:NLEVS)                            &
!             PRESSURE AT CENTRES OF LAYERS
     &   ,HEIGHT_THETA(NPD_FIELD,0:NLEVS)                               &
!             HEIGHT AT CENTRES OF LAYERS
     &   ,HEIGHT_RHO(NPD_FIELD,NLEVS)                                   &
!             HEIGHT AT EDGES OF LAYERS
     &   , TAC(NPD_FIELD, NLEVS)
!             TEMPERATURES AT CENTRES OF LAYERS
       Real, intent(in) ::                                              &
     &  rho_r2(npd_field,nlevs)                                         &
                                ! Air density*radius of earth**2 / kg m-1
     &, r_rho_levels(npd_field,nlevs)                                   &
                                      ! Height of rho levels / m
     &, r_theta_levels(npd_field,0:nlevs)                               &
                                           ! Height of theta levels / m
     &, q(npd_field,q_levels)                                           &
                                      ! Vapour content / kg kg-1
     &, qcl(npd_field,q_levels)                                         &
                                      ! Liquid water content / kg kg-1
     &, qcf(npd_field,q_levels)                                         &
                                      ! Ice content / kg kg-1
     &, qcf2(npd_field,q_levels)                                        &
                                      ! Second ice content / kg kg-1
     &, qrain(npd_field,q_levels)                                       &
                                      ! Rain water content / kg kg-1
     &, qgraup(npd_field,q_levels)    ! Graupel content / kg kg-1

       Logical, intent(in) ::                                           &
     &  l_mcr_qcf2                                                      &
                                      ! Use second prognostic ice
     &, l_mcr_qrain                                                     &
                                      ! Use prognostic rain
     &, l_mcr_qgraup                                                    &
                                      ! Use prognostic graupel
     &, l_mixing_ratio                ! Use mixing ratios
!
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_EXTRA_TOP                                                  &
!             Flag to create an extra top layer for radiation
     &   , L_BOUNDARY_TEMPERATURE
!             FLAG TO CALCULATE TEMPERATURES AT BOUNADRIES OF LAYERS.
!
!
      REAL                                                              &
                !, INTENT(OUT)
     &     D_MASS(NPD_PROFILE, NPD_LAYER)                               &
!             MASS THICKNESSES OF LAYERS
     &   , P(NPD_PROFILE, NPD_LAYER)                                    &
!             PRESSURE FIELD
     &   , T(NPD_PROFILE, NPD_LAYER)                                    &
     &   , T_BDY(NPD_PROFILE, 0: NPD_LAYER)                             &
!             TEMPERATURES AT EDGES OF LAYERS
     &   , T_SURFACE(NPD_PROFILE)                                       &
!             GATHERED TEMPERATURE OF SURFACE
     &   , T_SOLID(NPD_PROFILE)                                         &
!             GATHERED TEMPERATURE OF SOLID SURFACE
     &   , T_SEA(NPD_PROFILE)
!             GATHERED OPEN SEA TEMPERATURE
!
      REAL, INTENT(OUT) :: layer_heat_capacity(npd_profile, npd_layer)
!             Specific heat capacity of layer * d_mass
!
!     LOCAL VARIABLES.
!
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , II                                                           &
!             LOOP VARIABLE
     &   , L                                                            &
!             LOOP VARIABLE
     &   , LG                                                           &
!             INDEX TO GATHER
     &   , I_TOP_COPY                                                   &
!             Topmost layer where properties are set by copying the
!             input fields.
     &   , level_number  ! Model level in the upward counting
                         ! coordinate system
!
      REAL                                                              &
     &     PU                                                           &
!             PRESSURE FOR UPPER LAYER
     &   , PL                                                           &
!             PRESSURE FOR LOWER LAYER
     &   , PML1                                                         &
!             PRESSURE FOR INTERPOLATION
     &   , WTL                                                          &
!             WEIGHT FOR LOWER LAYER
     &   , WTU                                                          &
!             WEIGHT FOR UPPER LAYER
     &   , rho1(npd_profile)                                            &
                                         ! Air density / kg m-3
     &   , rho2(npd_profile)                                            &
                                         ! Air density / kg m-3
     &   , deltaz(npd_profile, npd_layer)                               &
                                         ! Layer thickness / m
     &   , rhodz_moist(npd_profile,npd_layer)                           &
                                              ! Thick*density / kg m-2
     &   , q_total(npd_profile)          ! Total moisture mr / kg kg-1
!
!
!
!
!     Calculate properties at the centres of layers: a top layer is
!     artificially created, if required, reaching from the top of
!     the model's actual atmosphere to zero pressure, assuming an
!     isothermal profile.
      IF (L_EXTRA_TOP) THEN
!       Note that we will continue to use this calculation for the
!       artificial top layer when the more accurate mixing ratio code
!       is used for the lower levels because we do not have the model
!       information to make a better estimate.
        DO L=1, N_PROFILE
          LG=I_GATHER(L)
          P(L, 1)=0.5E+00*P_LAYER_BOUNDARIES(LG, NLEVS)
          T(L, 1)=TAC(LG, NLEVS)
          D_MASS(L, 1)=P_LAYER_BOUNDARIES(LG, NLEVS)/G
          layer_heat_capacity(l, 1)=cp
        ENDDO
!       The second radiative layer will be the first to have properties
!       set by copying input fields.
        I_TOP_COPY=2
      ELSE
!       The second radiative layer will be the first to have properties
!       set by copying input fields.
        I_TOP_COPY=1
      ENDIF
      DO I=I_TOP_COPY, N_LAYER
        DO L=1, N_PROFILE
          LG=I_GATHER(L)
          P(L, I)=P_LAYER_CENTRES(LG, N_LAYER+1-I)
          T(L, I)=TAC(LG, N_LAYER+1-I)
          if (.not. l_mixing_ratio) then
            D_MASS(L, I)                                                &
     &        =ABS(P_LAYER_BOUNDARIES(LG, N_LAYER-I)                    &
     &        -P_LAYER_BOUNDARIES(LG, N_LAYER+1-I))/G
          End if
        ENDDO
      ENDDO
      if (l_mixing_ratio) then
        ! Perform more accurate calculation of the mass of the
        ! layer

        do i=i_top_copy, n_layer      ! This is the level in the
                                      ! downward counting coord. system

          level_number = n_layer+1-i  ! This is the level in the
                                      ! upward counting coord. system

          ! We note that rho_r2, r_rho_levels, q, qcl, qcf, qrain
          ! qcf2 and qgraup are in the upward
          ! counting label system and defined on all points. rho1,
          ! rho2, deltaz, rhodz_moist, rhodz_dry and q_total are in the
          ! downward counting label system and defined on gathered
          ! points.

          do l=1, n_profile
            lg=i_gather(l)

            ! Calculate densities at the boundaries of the layer
            ! by removing the r**2 term from rho_r2
            ! Rho1 is the density at the lower boundary
            rho1(l)= rho_r2(lg,level_number)                            &
     &                       /( r_rho_levels(lg,level_number) *         &
     &                          r_rho_levels(lg,level_number) )


            ! Check whether there is a rho level above the current
            ! moist level.
            if (level_number  <   nlevs) then
              ! Rho2 is the density at the upper boundary.
              rho2(l)= rho_r2(lg,level_number+1)                        &
     &                         /( r_rho_levels(lg,level_number+1) *     &
     &                            r_rho_levels(lg,level_number+1) )

              ! Calculate the average value of rho across the layer
              ! multiplied by the layer thickness and the layer
              ! thickness.
              rhodz_moist(l,i) =                                        &
     &                    rho2(l) * ( r_theta_levels(lg,level_number)   &
     &                              - r_rho_levels(lg,level_number) )   &
     &                 +  rho1(l) * ( r_rho_levels(lg,level_number+1)   &
     &                              - r_theta_levels(lg,level_number))
              deltaz(l,i) = r_rho_levels(lg,level_number+1)             &
     &                          - r_rho_levels(lg,level_number)

              If (level_number  ==  1) Then
                ! For the lowest layer we need to extend the lower
                ! boundary from the lowest (in height) rho level
                ! to the surface. The surface is the 0'th theta level
                ! when counting upwards.
                deltaz(l,i) = r_rho_levels(lg,2)                        &
     &                          - r_theta_levels(lg,0)
                rhodz_moist(l,i) = rhodz_moist(l,i)*deltaz(l,i)         &
     &                 / (r_rho_levels(lg,2)-r_rho_levels(lg,1))
              End if  ! level_number  ==  1

            Else
              ! For a top layer higher than the highest rho level
              ! we can calculate a pseudo rho level. We will assume
              ! it has a similar density to the rho level below
              ! and that the intervening theta level is in the centre
              ! of the layer.
              deltaz(l,i) = 2.0*(r_theta_levels(lg,level_number)        &
     &                            -r_rho_levels(lg,level_number))
              rhodz_moist(l,i) = rho1(l) * deltaz(l,i)

            End if  ! level_number  <   nlevs

            If (level_number  <=  q_levels) then
              ! Calculate total moisture
              q_total(l) = q(lg,level_number) + qcl(lg,level_number)    &
     &                                        + qcf(lg,level_number)

              ! Calculate the specific heat capacity of the moist air
              layer_heat_capacity(l,i)=cp + q(lg,level_number)*hcapv +  &
     &          qcl(lg,level_number)*hcapw + qcf(lg,level_number)*hcapi

              ! Add on contributions from optional prognostics
              If (l_mcr_qcf2) Then
                q_total(l) = q_total(l) + qcf2(lg,level_number)
                layer_heat_capacity(l,i) = layer_heat_capacity(l,i) +   &
     &            qcf2(lg,level_number)*hcapi
              End if  ! l_mcr_qcf2

              If (l_mcr_qrain) Then
                q_total(l) = q_total(l) + qrain(lg,level_number)
                layer_heat_capacity(l,i) = layer_heat_capacity(l,i) +   &
     &            qrain(lg,level_number)*hcapw
              End if  ! l_mcr_qrain

              If (l_mcr_qgraup) Then
                q_total(l) = q_total(l) + qgraup(lg,level_number)
                layer_heat_capacity(l,i) = layer_heat_capacity(l,i) +   &
     &            qgraup(lg,level_number)*hcapi
              End if  ! l_mcr_qgraup

            Else

              ! There is no moisture in the model
              q_total(l) = 0.0
              layer_heat_capacity(l,i) = cp

            End if  ! level_number le q_levels

            ! Convert from the moist density of air to the dry
            ! density of air
            d_mass(l,i) = rhodz_moist(l,i) / (1.0 + q_total(l))

          End Do  ! l
        End Do  ! i

        ! If l_mixing ratio is true, layer_heat_capacity is the true
        ! heat capacity of the moist air including condensate.
        Do i=1,n_layer
          Do l=1,n_profile
            layer_heat_capacity(l,i)=d_mass(l,i)*                       &
     &                                 layer_heat_capacity(l,i)
          End Do
        End Do

      Else

        ! If l_mixing_ratio is false, layer_heat_capacity is based
        ! on the specific heat of dry air and a mass derived from the 
        ! hydrostatic approximation.
        Do i=1,n_layer
          Do l=1,n_profile
            layer_heat_capacity(l,i)=d_mass(l,i)*cp
          End Do
        End Do

      End if  ! L_mixing_ratio
!
!
      IF (L_BOUNDARY_TEMPERATURE) THEN
!
!        GATHER THE SURFACE TEMPERATURE.
         DO L=1, N_PROFILE
            LG=I_GATHER(L)
            T_SURFACE(L)=TSTAR(LG)
            T_SOLID(L)=TSTAR_SOLID(LG)
            T_SEA(L)=TSTAR_SEA(LG)
         ENDDO
!
!        INTERPOLATE TEMPERATURES AT THE BOUNDARIES OF LAYERS
!        FROM THE EXNER FUNCTION.
         DO L=1, N_PROFILE
            LG=I_GATHER(L)
!
!           TAKE THE TEMPERATURE OF THE AIR JUST ABOVE THE SURFACE AS
!           THE TEMPERATURE AT THE MIDDLE OF THE BOTTOM LAYER.
            T_BDY(L, N_LAYER)=TAC(LG, 1)
!        Set the temperature at the top of the model by isothermal
!        extrapolation of the temperature in the middle of the topmost
!        layer used in the full atmospheric model.
!        otherwise.
         T_BDY(L, 0)=TAC(LG, NLEVS)
!
         ENDDO
!
!        If an extra layer is to be added, the temperature at the
!        bottom of this layer must be set by isothermal extrapolation.
         IF (L_EXTRA_TOP) THEN
           DO L=1, N_PROFILE
             T_BDY(L, 1)=TAC(I_GATHER(L), NLEVS)
           ENDDO
         ENDIF
!
!
!        Now set the temperatures at interior boundaries of the
!        profile supplied using interpolation with the EXNER function.
         DO I=I_TOP_COPY, N_LAYER-1
            II=N_LAYER-I
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               WTU=HEIGHT_RHO(LG,II+1)-HEIGHT_THETA(LG,II)
               WTL=HEIGHT_THETA(LG,II+1)-HEIGHT_RHO(LG,II+1)
               T_BDY(L, I)=(WTU*TAC(LG, N_LAYER+1-I)                    &
     &            +WTL*TAC(LG, N_LAYER-I))/(WTL+WTU)
            ENDDO
         ENDDO
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE R2_SET_THERMODYNAMIC
