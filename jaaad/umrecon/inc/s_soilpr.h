! Description:
!   This deck defines the soil parameters used in the SCM.
!
! Current Code Owner: R Wong
!
!
! Declarations:
      Integer nsoilp            ! Number of possible soil parameters
      Parameter(nsoilp=4)       !

      Real :: b_exp_typ(nsoilp)                  
                                ! Single Layer :
                                !     (C_EAG in code) Eagleson's exponent for
                                !     calc. sub surf. runoff, P253.4
                                ! Multilayer Hydrology:
                                !     Exponent used in calculation of soil
                                !     water suction and hydraulic conductivity
                                !     (known as B_WAG in HYDROL2A)
                                ! MOSES:
                                !     Exponent used in calculation of soil
                                !     water suction and hydraulic conductivity
                                !     (known as B Clapp-Hornberger exponent)

     Real :: sathh_typ(nsoilp)  ! Single layer :
                                !     Dummy
                                ! MOSES/Multilayer hydrology :
                                !     Saturated soil water suction

     Real :: satcon_typ(nsoilp) ! Saturated hydrological conductivity of
                                ! the soil (kg/m2/s)

     Real :: hcap_typ(nsoilp)   ! Soil heat capacity (J/K/m**3)

     Real :: hcon_typ(nsoilp)   ! Soil thermal conductivity (W/m/K)



     Real :: v_crit_typ(nsoilp) ! Volumetric soil moisture content the
                                ! critical point; below this value
                                ! evaporation falls below its max
                                ! (m**3 per m**3 soil)

     Real :: v_wilt_typ(nsoilp) ! Volumetric soil moisture content at wilting
                                ! point (m**3/m**3)

     Real :: v_sat_typ(nsoilp)  ! Volumetric soil moisture content at
                                ! saturation(m**3/m**3 soil)


!    (i)  Values of SATCON, V_SAT, B_WAG (B_EXP), SATHH are given by
!         Wageningen for the Maulem/Van Genuchten curves:
!          SATHH = 1 / ALPHA
!          B_WAG = (1-M) / M = 1 / (N-1)
!    (ii) Values of V_CRIT , V_SAT, V_WILT are used which apply
!         for the soil moisture variable defined as the toal soil
!         moisture minus the residuals given by Wageningen.
!    (iii)Soil types ICE, CLAY, LOAM, LOAMY SAND

      COMMON /SOIL_DATA/ b_exp_typ,  hcap_typ, hcon_typ, satcon_typ           &
                       , sathh_typ, v_crit_typ, v_sat_typ, v_wilt_typ

!-----------------------------------------------------------------------------
