#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set pointers to types of clouds
!
! Method:
!       The types of condensate included are examined. Their phases
!       are set and depending on the representation of clouds adopted
!       it is determined to which type of cloud they contribute.
!
! Current owner of code: James Manners
!
! Description of code:
!   Fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SET_CLOUD_POINTER(IERR                                 &
     &   , N_CONDENSED, TYPE_CONDENSED, I_CLOUD_REPRESENTATION          &
     &   , L_DROP, L_ICE                                                &
     &   , I_PHASE_CMP, I_CLOUD_TYPE, L_CLOUD_CMP                       &
     &   , ND_CLOUD_COMPONENT                                           &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     Dimensions of arrays
      INTEGER, INTENT(IN) ::                                            &
     &  ND_CLOUD_COMPONENT
!         Maximum number of condensed components allowed in clouds
!
!     Inclusion of header files.
#include "def_std_io_icf3z.h"
#include "error_pcf3z.h"
#include "cloud_component_pcf3z.h"
#include "cloud_representation_pcf3z.h"
#include "cloud_type_pcf3z.h"
#include "phase_pcf3z.h"
!
!     Dummy variables.
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
      INTEGER, INTENT(IN) ::                                            &
     &    N_CONDENSED                                                   &
!           Number of condensed components
     &  , TYPE_CONDENSED(ND_CLOUD_COMPONENT)                            &
!           Types of components
     &  , I_CLOUD_REPRESENTATION
!           Representation of clouds used
      LOGICAL, INTENT(IN) ::                                            &
     &    L_DROP                                                        &
!           Flag for inclusion of droplets
     &  , L_ICE
!           Flag for inclusion of ice crystals
!
      INTEGER, INTENT(OUT) ::                                           &
     &    I_PHASE_CMP(ND_CLOUD_COMPONENT)                               &
!           Phases of components
     &  , I_CLOUD_TYPE(ND_CLOUD_COMPONENT)
!           Types of cloud to which each component contributes
      LOGICAL, INTENT(OUT) ::                                           &
     &    L_CLOUD_CMP(ND_CLOUD_COMPONENT)
!           Logical switches to `include' components
!
!
!     Local variables
      INTEGER                                                           &
     &    K
!           Loop variable
!
!
!
      DO K=1, N_CONDENSED
!
        I_CLOUD_TYPE(K)=0
!       Set pointers for valid condensed components.
        IF (I_CLOUD_REPRESENTATION == IP_CLOUD_HOMOGEN) THEN
!
          IF (TYPE_CONDENSED(K) == IP_CLCMP_ST_WATER) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_HOMOGEN
          ELSEIF (TYPE_CONDENSED(K) == IP_CLCMP_ST_ICE) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_HOMOGEN
          ENDIF
!
        ELSE IF (I_CLOUD_REPRESENTATION == IP_CLOUD_ICE_WATER) THEN
!
          IF (TYPE_CONDENSED(K) == IP_CLCMP_ST_WATER) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_WATER
          ELSEIF (TYPE_CONDENSED(K) == IP_CLCMP_ST_ICE) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_ICE
          ENDIF
!
        ELSE IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CONV_STRAT) THEN
!
          IF (TYPE_CONDENSED(K) == IP_CLCMP_ST_WATER) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_STRAT
          ELSE IF (TYPE_CONDENSED(K) == IP_CLCMP_ST_ICE) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_STRAT
          ELSEIF (TYPE_CONDENSED(K) == IP_CLCMP_CNV_WATER) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_CONV
          ELSE IF (TYPE_CONDENSED(K) == IP_CLCMP_CNV_ICE) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_CONV
          ENDIF
!
        ELSE IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CSIW) THEN
!
          IF (TYPE_CONDENSED(K) == IP_CLCMP_ST_WATER) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_SW
          ELSEIF (TYPE_CONDENSED(K) == IP_CLCMP_ST_ICE) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_SI
          ELSEIF (TYPE_CONDENSED(K) == IP_CLCMP_CNV_WATER) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_CW
          ELSEIF (TYPE_CONDENSED(K) == IP_CLCMP_CNV_ICE) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_CI
          ENDIF
#if defined(XICE)
!
        ELSE IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CSIW_CRYS) THEN
!
          IF (TYPE_CONDENSED(K) == IP_CLCMP_ST_WATER) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_SW
          ELSE IF (TYPE_CONDENSED(K) == IP_CLCMP_CNV_WATER) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_CW
          ELSE IF (TYPE_CONDENSED(K) == IP_CLCMP_ST_ICE_COL) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_SI
          ELSEIF (TYPE_CONDENSED(K) == IP_CLCMP_CNV_ICE_COL) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_CI
          ELSEIF (TYPE_CONDENSED(K) == IP_CLCMP_ST_ICE_ROS) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_SI
          ELSEIF (TYPE_CONDENSED(K) == IP_CLCMP_CNV_ICE_ROS) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_CI
          ELSEIF (TYPE_CONDENSED(K) == IP_CLCMP_ST_ICE_PLT) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_SI
          ELSEIF (TYPE_CONDENSED(K) == IP_CLCMP_CNV_ICE_PLT) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_CI
          ELSEIF (TYPE_CONDENSED(K) == IP_CLCMP_ST_ICE_PYC) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_SI
          ELSEIF (TYPE_CONDENSED(K) == IP_CLCMP_CNV_ICE_PYC) THEN
            I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_CI
          ENDIF
#endif
!
        ENDIF
!
!       Check for 0 flagging illegal types.
        IF (I_CLOUD_TYPE(K) == 0) THEN
          WRITE(IU_ERR, '(/A)')                                         &
     &      '*** Error: A component is not compatible with the '        &
     &      //'representation of clouds selected.'
          IERR=I_ERR_FATAL
          RETURN
        ENDIF
!
        IF (TYPE_CONDENSED(K) == IP_CLCMP_ST_WATER) THEN
!
          I_PHASE_CMP(K)=IP_PHASE_WATER
          L_CLOUD_CMP(K)=L_DROP
!
        ELSE IF (TYPE_CONDENSED(K) == IP_CLCMP_ST_ICE) THEN
!
          I_PHASE_CMP(K)=IP_PHASE_ICE
          L_CLOUD_CMP(K)=L_ICE
!
        ELSE IF (TYPE_CONDENSED(K) == IP_CLCMP_CNV_WATER) THEN
!
          I_PHASE_CMP(K)=IP_PHASE_WATER
          L_CLOUD_CMP(K)=L_DROP
!
        ELSE IF (TYPE_CONDENSED(K) == IP_CLCMP_CNV_ICE) THEN
!
          I_PHASE_CMP(K)=IP_PHASE_ICE
          L_CLOUD_CMP(K)=L_ICE
!
        ENDIF
!
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SET_CLOUD_POINTER
#endif
