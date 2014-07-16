#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the parametrization schemes for clouds.
!
! Purpose:
!   The parametrization schemes for each component within a cloud
!   are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: James Manners
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_CLOUD_PARAMETRIZATION(IERR, N_BAND              &
     &   , I_ST_WATER, I_CNV_WATER, I_ST_ICE, I_CNV_ICE                 &
     &   , L_DROP_TYPE, I_DROP_PARAMETRIZATION, DROP_PARAMETER_LIST     &
     &   , DROP_PARM_MIN_DIM, DROP_PARM_MAX_DIM                         &
     &   , L_ICE_TYPE, I_ICE_PARAMETRIZATION, ICE_PARAMETER_LIST        &
     &   , ICE_PARM_MIN_DIM, ICE_PARM_MAX_DIM                           &
     &   , I_CONDENSED_PARAM, CONDENSED_PARAM_LIST                      &
     &   , CONDENSED_MIN_DIM, CONDENSED_MAX_DIM                         &
     &   , NPD_BAND, NPD_DROP_TYPE, NPD_ICE_TYPE, NPD_CLOUD_PARAMETER   &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED.
#include "dimfix3a.h"
#include "cldcmp3a.h"
#include "stdio3a.h"
#include "error3a.h"
!
!
!     DUMMY ARGUMENTS:
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
!
!     SIZES OF ARRAYS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_BAND                                                     &
!             MAXIMUM NUMBER OF SPECTRAL BANDS
     &   , NPD_DROP_TYPE                                                &
!             MAXIMUM NUMBER OF TYPES OF DROPLETS
     &   , NPD_ICE_TYPE                                                 &
!             MAXIMUM NUMBER OF TYPES OF ICE CRYSTALS
     &   , NPD_CLOUD_PARAMETER
!             MAXIMUM NUMBER OF PARAMETERS FOR CLOUDS
!
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_BAND
!             NUMBER OF SPECTRAL BANDS
!
!     TYPES OF DROPLETS AND CRYSTALS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_ST_WATER                                                   &
!             TYPE OF WATER DROPLETS IN STRATIFORM CLOUDS
     &   , I_CNV_WATER                                                  &
!             TYPE OF WATER DROPLETS IN CONVECTIVE CLOUDS
     &   , I_ST_ICE                                                     &
!             TYPE OF ICE CRYSTALS IN STRATIFORM CLOUDS
     &   , I_CNV_ICE
!             TYPE OF ICE CRYSTALS IN CONVECTIVE CLOUDS
!
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_DROP_TYPE(NPD_DROP_TYPE)                                   &
!             FLAGS FOR TYPES OF DROPLET PRESENT
     &   , L_ICE_TYPE(NPD_ICE_TYPE)
!             FLAGS FOR TYPES OF ICE CRYSTAL PRESENT
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_DROP_PARAMETRIZATION(NPD_DROP_TYPE)                        &
!             PARAMETRIZATIONS OF TYPES OF DROPLETS
     &   , I_ICE_PARAMETRIZATION(NPD_ICE_TYPE)
!             PARAMETRIZATIONS OF TYPES OF ICE CRYSTALS
      REAL                                                              &
                !, INTENT(IN)
     &     DROP_PARAMETER_LIST(NPD_CLOUD_PARAMETER                      &
     &        , NPD_BAND, NPD_DROP_TYPE)                                &
!             PARAMETERS FOR OPTICAL PARAMETRIZATIONS OF DROPLETS
     &   , DROP_PARM_MIN_DIM(NPD_DROP_TYPE)                             &
!             MINIMUM SIZE OF DROPLETS PERMITTED IN PARAMETRIZATIONS
     &   , DROP_PARM_MAX_DIM(NPD_DROP_TYPE)                             &
!             MAXIMUM SIZE OF DROPLETS PERMITTED IN PARAMETRIZATIONS
     &   , ICE_PARAMETER_LIST(NPD_CLOUD_PARAMETER                       &
     &        , NPD_BAND, NPD_ICE_TYPE)                                 &
!             PARAMETERS FOR OPTICAL PARAMETRIZATIONS OF ICE CRYSTALS
     &   , ICE_PARM_MIN_DIM(NPD_ICE_TYPE)                               &
!             MINIMUM SIZE OF ICE CRYSTALS PERMITTED IN PARAMETRIZATIONS
     &   , ICE_PARM_MAX_DIM(NPD_ICE_TYPE)
!             MAXIMUM SIZE OF ICE CRYSTALS PERMITTED IN PARAMETRIZATIONS
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     I_CONDENSED_PARAM(NPD_CLOUD_COMPONENT)
!             TYPES OF PARAMETRIZATION USED FOR CONDENSED
!             COMPONENTS IN CLOUDS
      REAL                                                              &
                !, INTENT(OUT)
     &     CONDENSED_PARAM_LIST(NPD_CLOUD_PARAMETER                     &
     &        , NPD_CLOUD_COMPONENT, NPD_BAND)                          &
!             COEFFICIENTS FOR PARAMETRIZATION OF CONDENSED PHASES
     &   , CONDENSED_MIN_DIM(NPD_CLOUD_COMPONENT)                       &
!             MINIMUM DIMENSION OF EACH CONDENSED COMPONENT
     &   , CONDENSED_MAX_DIM(NPD_CLOUD_COMPONENT)
!             MAXIMUM DIMENSION OF EACH CONDENSED COMPONENT
!
!
!     LOCAL VARIABLES:
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , J                                                            &
!             LOOP VARIABLE
     &   , I_SCHEME
!             PARAMETRIZATION SCHEME
!
!     FUNCTIONS CALLED:
      INTEGER                                                           &
     &     SET_N_CLOUD_PARAMETER
!             FUNCTION TO FIND NUMBER OF PARAMETERS FOR CLOUDS
      EXTERNAL                                                          &
     &     SET_N_CLOUD_PARAMETER
!
!
!
!     SELECT PARAMETRIZATION FOR WATER IN STRATIFORM CLOUDS:
!
      IF ( (I_ST_WATER <= NPD_DROP_TYPE).AND.                           &
     &     (L_DROP_TYPE(I_ST_WATER)) ) THEN
         I_SCHEME=I_DROP_PARAMETRIZATION(I_ST_WATER)
         I_CONDENSED_PARAM(IP_CLCMP_ST_WATER)=I_SCHEME
         CONDENSED_MIN_DIM(IP_CLCMP_ST_WATER)                           &
     &      =DROP_PARM_MIN_DIM(I_ST_WATER)
         CONDENSED_MAX_DIM(IP_CLCMP_ST_WATER)                           &
     &      =DROP_PARM_MAX_DIM(I_ST_WATER)
      ELSE
         WRITE(IU_ERR, '(/A, /A)') '*** ERROR: NO DATA EXIST FOR TYPE ' &
     &     , 'OF DROPLET SELECTED IN STRATIFORM WATER CLOUDS.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      DO I=1, N_BAND
! DEPENDS ON: set_n_cloud_parameter
         DO J=1, SET_N_CLOUD_PARAMETER(I_SCHEME, IP_CLCMP_ST_WATER)
            CONDENSED_PARAM_LIST(J, IP_CLCMP_ST_WATER, I)               &
     &         =DROP_PARAMETER_LIST(J, I, I_ST_WATER)
         ENDDO
      ENDDO
!
!
!     SELECT PARAMETRIZATION FOR WATER IN CONVECTIVE CLOUDS:
!
      IF ( (I_CNV_WATER <= NPD_DROP_TYPE).AND.                          &
     &     (L_DROP_TYPE(I_CNV_WATER)) ) THEN
         I_SCHEME=I_DROP_PARAMETRIZATION(I_CNV_WATER)
         I_CONDENSED_PARAM(IP_CLCMP_CNV_WATER)=I_SCHEME
         CONDENSED_MIN_DIM(IP_CLCMP_CNV_WATER)                          &
     &      =DROP_PARM_MIN_DIM(I_CNV_WATER)
         CONDENSED_MAX_DIM(IP_CLCMP_CNV_WATER)                          &
     &      =DROP_PARM_MAX_DIM(I_CNV_WATER)
      ELSE
         WRITE(IU_ERR, '(/A, /A)') '*** ERROR: NO DATA EXIST FOR TYPE ' &
     &     , 'OF CRYSTAL SELECTED IN CONVECTIVE WATER CLOUDS.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      DO I=1, N_BAND
! DEPENDS ON: set_n_cloud_parameter
         DO J=1, SET_N_CLOUD_PARAMETER(I_SCHEME, IP_CLCMP_CNV_WATER)
            CONDENSED_PARAM_LIST(J, IP_CLCMP_CNV_WATER, I)              &
     &         =DROP_PARAMETER_LIST(J, I, I_CNV_WATER)
         ENDDO
      ENDDO
!
!
!     SELECT PARAMETRIZATION FOR ICE IN STRATIFORM CLOUDS:
!
      IF ( (I_ST_ICE <= NPD_ICE_TYPE).AND.                              &
     &     (L_ICE_TYPE(I_ST_ICE)) ) THEN
         I_SCHEME=I_ICE_PARAMETRIZATION(I_ST_ICE)
         I_CONDENSED_PARAM(IP_CLCMP_ST_ICE)=I_SCHEME
         CONDENSED_MIN_DIM(IP_CLCMP_ST_ICE)                             &
     &      =ICE_PARM_MIN_DIM(I_ST_ICE)
         CONDENSED_MAX_DIM(IP_CLCMP_ST_ICE)                             &
     &      =ICE_PARM_MAX_DIM(I_ST_ICE)
      ELSE
         WRITE(IU_ERR, '(/A, /A)') '*** ERROR: NO DATA EXIST FOR TYPE ' &
     &      , 'OF CRYSTAL SELECTED IN STRATIFORM ICE CLOUDS.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      DO I=1, N_BAND
! DEPENDS ON: set_n_cloud_parameter
         DO J=1, SET_N_CLOUD_PARAMETER(I_SCHEME, IP_CLCMP_ST_ICE)
            CONDENSED_PARAM_LIST(J, IP_CLCMP_ST_ICE, I)                 &
     &         =ICE_PARAMETER_LIST(J, I, I_ST_ICE)
         ENDDO
      ENDDO
!
!
!     SELECT PARAMETRIZATION FOR ICE IN CONVECTIVE CLOUDS:
!
      IF ( (I_CNV_ICE <= NPD_ICE_TYPE).AND.                             &
     &     (L_ICE_TYPE(I_CNV_ICE)) ) THEN
         I_SCHEME=I_ICE_PARAMETRIZATION(I_CNV_ICE)
         I_CONDENSED_PARAM(IP_CLCMP_CNV_ICE)=I_SCHEME
         CONDENSED_MIN_DIM(IP_CLCMP_CNV_ICE)                            &
     &      =ICE_PARM_MIN_DIM(I_CNV_ICE)
         CONDENSED_MAX_DIM(IP_CLCMP_CNV_ICE)                            &
     &      =ICE_PARM_MAX_DIM(I_CNV_ICE)
      ELSE
         WRITE(IU_ERR, '(/A, /A)') '*** ERROR: NO DATA EXIST FOR TYPE ' &
     &      , 'OF CRYSTAL SELECTED IN CONVECTIVE ICE CLOUDS.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
      DO I=1, N_BAND
! DEPENDS ON: set_n_cloud_parameter
         DO J=1, SET_N_CLOUD_PARAMETER(I_SCHEME, IP_CLCMP_CNV_ICE)
            CONDENSED_PARAM_LIST(J, IP_CLCMP_CNV_ICE, I)                &
     &         =ICE_PARAMETER_LIST(J, I, I_CNV_ICE)
         ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE R2_SET_CLOUD_PARAMETRIZATION
#endif
