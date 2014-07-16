#if defined(A71_1A) && !defined(A17_2A) && defined(A17_2B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Function to calculate cloud droplet number concentration.
!
! Purpose:
!   Cloud droplet number concentration is calculated from aerosol
!   concentration or else fixed values are assigned.
!
! Method:
!   Sulphate aerosol mass concentration is converted to number
!   concentration by assuming a log-normal size distribution.
!   Sea-salt and/or biomass-burning aerosols may then
!   be added if required. The total is then converted to cloud
!   droplet concentration following the parametrization of
!   Jones et al. (1994) and lower limits are imposed.
!   Alternatively, fixed droplet values are assigned if the
!   parametrization is not required.
!
! Current Owner of Code: A. Jones
!
! History:
!       Version         Date                    Comment
!       6.2             20-10-05                Based on NDROP1.
!                                               Aitken-mode SO4 no
!                                               longer used & aerosol
!                                               parameters revised.
!                                               (A. Jones)
!
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      FUNCTION NUMBER_DROPLET(L_AEROSOL_DROPLET, L_NH42SO4              &
     &   , AITKEN_SULPHATE, ACCUM_SULPHATE, DISS_SULPHATE               &
     &   , L_SEASALT_CCN, SEA_SALT_FILM, SEA_SALT_JET                   &
     &   , L_BIOGENIC_CCN, BIOGENIC                                     &
     &   , L_BIOMASS_CCN, BIOMASS_AGED, BIOMASS_CLOUD                   &
     &   , L_OCFF_CCN, OCFF_AGED, OCFF_CLOUD                            &
     &   , DENSITY_AIR                                                  &
     &   , SNOW_DEPTH                                                   &
     &   , LAND_FRACT                                                   &
     &   , NTOT_LAND, NTOT_SEA                                          &
     &   )



      IMPLICIT NONE


!     Comdecks included:
#include "c_pi.h"


      LOGICAL                                                           &
     &     L_AEROSOL_DROPLET                                            &
!             Flag to use aerosols to find droplet number
     &   , L_NH42SO4                                                    &
!             Is the input "sulphate" aerosol in the form of
!             ammonium sulphate (T) or just sulphur (F)?
     &   , L_SEASALT_CCN                                                &
!             Is sea-salt aerosol to be used?
     &   , L_BIOMASS_CCN                                                &
!             Is biomass smoke aerosol to be used?
     &   , L_BIOGENIC_CCN                                               &
!             Is biogenic aerosol to be used?
     &   , L_OCFF_CCN
!             Is fossil-fuel organic carbon aerosol to be used?

      REAL                                                              &
     &     AITKEN_SULPHATE                                              &
!             Dummy in version 2 of the routine
     &   , ACCUM_SULPHATE                                               &
!             Mixing ratio of accumulation-mode sulphate aerosol
     &   , DISS_SULPHATE                                                &
!             Mixing ratio of dissolved sulphate aerosol
     &   , BIOMASS_AGED                                                 &
!             Mixing ratio of aged biomass smoke
     &   , BIOMASS_CLOUD                                                &
!             Mixing ratio of in-cloud biomass smoke
     &   , SEA_SALT_FILM                                                &
!             Number concentration of film-mode sea salt aerosol (m-3)
     &   , SEA_SALT_JET                                                 &
!             Number concentration of jet-mode sea salt aerosol (m-3)
     &   , BIOGENIC                                                     &
!             Mixing ratio of biogenic aerosol
     &   , OCFF_AGED                                                    &
!             Mixing ratio of aged fossil-fuel organic carbon
     &   , OCFF_CLOUD                                                   &
!             Mixing ratio of in-cloud fossil-fuel organic carbon
     &   , DENSITY_AIR                                                  &
!             Density of air (kg m-3)
     &   , SNOW_DEPTH                                                   &
!             Snow depth (m; >5000 is flag for ice-sheets)
     &   , LAND_FRACT                                                   &
!             Land fraction
     &   , NTOT_LAND                                                    &
!             Droplet number over land if parameterization is off (m-3)
     &   , NTOT_SEA
!             Droplet number over sea if parameterization is off (m-3)

      REAL                                                              &
     &     NUMBER_DROPLET
!             Returned number concentration of cloud droplets (m-3)



!     Local variables:

      REAL                                                              &
     &     PARTICLE_VOLUME                                              &
!             Mean volume of aerosol particle
     &   , N_CCN
!             Number density of CCN

      REAL                                                              &
     &     RADIUS_0_NH42SO4                                             &
!             Median radius of log-normal distribution for (NH4)2SO4
     &   , SIGMA_0_NH42SO4                                              &
!             Geometric standard deviation of same
     &   , DENSITY_NH42SO4                                              &
!             Density of ammonium sulphate aerosol
     &   , RADIUS_0_BIOMASS                                             &
!             Median radius of log-normal distribution for biomass smoke
     &   , SIGMA_0_BIOMASS                                              &
!             Geometric standard deviation of same
     &   , DENSITY_BIOMASS                                              &
!             Density of biomass smoke aerosol
     &   , RADIUS_0_BIOGENIC                                            &
!             Median radius of log-normal dist. for biogenic aerosol
     &   , SIGMA_0_BIOGENIC                                             &
!             Geometric standard deviation of same
     &   , DENSITY_BIOGENIC                                             &
!             Density of biogenic aerosol
     &   , RADIUS_0_OCFF                                                &
!             Median radius of log-normal dist. for OCFF aerosol
     &   , SIGMA_0_OCFF                                                 &
!             Geometric standard deviation of same
     &   , DENSITY_OCFF
!             Density of OCFF aerosol

      PARAMETER(                                                        &
     &     RADIUS_0_NH42SO4=9.5E-8                                      &
     &   , SIGMA_0_NH42SO4=1.4                                          &
     &   , DENSITY_NH42SO4=1.769E+03                                    &
     &   , RADIUS_0_BIOMASS=1.2E-07                                     &
     &   , SIGMA_0_BIOMASS=1.30                                         &
     &   , RADIUS_0_BIOGENIC=9.5E-08                                    &
     &   , SIGMA_0_BIOGENIC=1.50                                        &
     &   , DENSITY_BIOGENIC=1.3E+03                                     &
     &   , DENSITY_BIOMASS=1.35E+03                                     &
     &   , RADIUS_0_OCFF=0.12E-06                                       &
     &   , SIGMA_0_OCFF=1.30                                            &
     &   , DENSITY_OCFF=1350.0                                          &
     &   )



      IF (L_AEROSOL_DROPLET) THEN

!        If active, aerosol concentrations are used to calculate the
!        number of CCN, which is then used to determine the number
!        concentration of cloud droplets (m-3).

         PARTICLE_VOLUME=(4.0E+00*PI/3.0E+00)*RADIUS_0_NH42SO4**3       &
     &                         *EXP(4.5E+00*(LOG(SIGMA_0_NH42SO4))**2)

         IF (L_NH42SO4) THEN
!           Input data have already been converted to ammonium sulphate.
            N_CCN=(ACCUM_SULPHATE+DISS_SULPHATE)                        &
     &                  *DENSITY_AIR/(DENSITY_NH42SO4*PARTICLE_VOLUME)
         ELSE
!           Convert m.m.r. of sulphur to ammonium sulphate by
!           multiplying by ratio of molecular weights:
            N_CCN=(ACCUM_SULPHATE+DISS_SULPHATE)*4.125                  &
     &                  *DENSITY_AIR/(DENSITY_NH42SO4*PARTICLE_VOLUME)
         ENDIF

         IF (L_SEASALT_CCN) THEN
            N_CCN=N_CCN+SEA_SALT_FILM+SEA_SALT_JET
         ENDIF

         IF (L_BIOMASS_CCN) THEN
            PARTICLE_VOLUME=(4.0E+00*PI/3.0E+00)*RADIUS_0_BIOMASS**3    &
     &                         *EXP(4.5E+00*(LOG(SIGMA_0_BIOMASS))**2)
            N_CCN=N_CCN+((BIOMASS_AGED+BIOMASS_CLOUD)*DENSITY_AIR       &
     &                             /(DENSITY_BIOMASS*PARTICLE_VOLUME))
         ENDIF

         IF (L_BIOGENIC_CCN) THEN
            PARTICLE_VOLUME=(4.0E+00*PI/3.0E+00)*RADIUS_0_BIOGENIC**3   &
     &                         *EXP(4.5E+00*(LOG(SIGMA_0_BIOGENIC))**2)
            N_CCN=N_CCN+(BIOGENIC*DENSITY_AIR                           &
     &                             /(DENSITY_BIOGENIC*PARTICLE_VOLUME))
         ENDIF

         IF (L_OCFF_CCN) THEN
            PARTICLE_VOLUME=(4.0E+00*PI/3.0E+00)*RADIUS_0_OCFF**3       &
     &                         *EXP(4.5E+00*(LOG(SIGMA_0_OCFF))**2)
            N_CCN=N_CCN+((OCFF_AGED+OCFF_CLOUD)*DENSITY_AIR             &
     &                             /(DENSITY_OCFF*PARTICLE_VOLUME))
         ENDIF

!        Apply relation of Jones et al. (1994) to get droplet number
!        and apply minimum value (equivalent to 5 cm-3):

         NUMBER_DROPLET=3.75E+08*(1.0E+00-EXP(-2.5E-9*N_CCN))

         IF (NUMBER_DROPLET  <   5.0E+06) THEN
            NUMBER_DROPLET=5.0E+06
         ENDIF

!        If gridbox is more than 20% land AND this land is not covered
!        by an ice-sheet, use larger minimum droplet number (=35 cm-3):

         IF (LAND_FRACT  >   0.2 .AND. SNOW_DEPTH  <   5000.0           &
     &                        .AND. NUMBER_DROPLET  <   35.0E+06) THEN
            NUMBER_DROPLET=35.0E+06
         ENDIF

      ELSE

!        Without aerosols, the number of droplets is fixed; a simple
!        50% criterion is used for land or sea in this case.

         IF (LAND_FRACT  >=  0.5) THEN
            NUMBER_DROPLET=NTOT_LAND
         ELSE
            NUMBER_DROPLET=NTOT_SEA
         ENDIF
!
      ENDIF
!
!
!
      RETURN
      END FUNCTION NUMBER_DROPLET
#endif
