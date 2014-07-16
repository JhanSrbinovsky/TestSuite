#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to aggregate clouds into regions.
!
! Method:
!       The clouds in a layer are combined in groups to form regions
!       which will be considered as bulk entities in the solution of the
!       equation of transfer. The extents of these regions are also
!       determined.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE AGGREGATE_CLOUD(IERR                                   &
     &   , N_PROFILE, N_LAYER, N_CLOUD_TOP                              &
     &   , I_CLOUD, I_CLOUD_REPRESENTATION, N_CLOUD_TYPE                &
     &   , FRAC_CLOUD                                                   &
     &   , I_REGION_CLOUD, FRAC_REGION                                  &
     &   , ND_PROFILE, ND_LAYER, ND_CLOUD_TYPE, ND_REGION               &
     &   , ID_CT                                                        &
     &   )
!
!
!
!
      IMPLICIT NONE
!
!
!     Dummy array sizes
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Maximum number of profiles
     &  , ND_LAYER                                                      &
!           Maximum number of layers
     &  , ND_CLOUD_TYPE                                                 &
!           Maximum number of types of cloud
     &  , ND_REGION                                                     &
!           Maximum number of cloudy regions
     &  , ID_CT
!           Topmost declared cloudy layer
!
!     Inclusion of header files.
#include "def_std_io_icf3z.h"
#include "error_pcf3z.h"
#include "cloud_representation_pcf3z.h"
#include "cloud_type_pcf3z.h"
#include "cloud_region_pcf3z.h"
#include "cloud_scheme_pcf3z.h"
#include "c_kinds.h"
!
!     Dummy variables.
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , N_CLOUD_TOP
!           Topmost cloudy layer
      INTEGER, INTENT(IN) ::                                            &
     &    I_CLOUD                                                       &
!           Cloud scheme used
     &  , I_CLOUD_REPRESENTATION                                        &
!           Representation of clouds used
     &  , N_CLOUD_TYPE
!           Number of types of cloud
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FRAC_CLOUD(ND_PROFILE, ID_CT: ND_LAYER, ND_CLOUD_TYPE)
!           Fractions of each type of cloud
!
      INTEGER, INTENT(OUT) ::                                           &
     &    I_REGION_CLOUD(ND_CLOUD_TYPE)
!           Regions in which particular types of cloud fall
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FRAC_REGION(ND_PROFILE, ID_CT: ND_LAYER, ND_REGION)
!           Fractions of total cloud occupied by each region
!
!
!     Local variables
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L                                                             &
!           Loop variable
     &  , K
!           Loop variable
!
!
!
      IF ( (I_CLOUD == IP_CLOUD_TRIPLE).OR.                             &
     &     (I_CLOUD == IP_CLOUD_PART_CORR_CNV) ) THEN
!
        IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CSIW) THEN
!
          DO K=1, N_CLOUD_TYPE
            IF (K == IP_CLOUD_TYPE_SW) THEN
              I_REGION_CLOUD(K)=IP_REGION_STRAT
            ELSE IF (K == IP_CLOUD_TYPE_SI) THEN
              I_REGION_CLOUD(K)=IP_REGION_STRAT
            ELSE IF (K == IP_CLOUD_TYPE_CW) THEN
              I_REGION_CLOUD(K)=IP_REGION_CONV
            ELSE IF (K == IP_CLOUD_TYPE_CI) THEN
              I_REGION_CLOUD(K)=IP_REGION_CONV
            ENDIF
          ENDDO
!
          DO I=N_CLOUD_TOP, N_LAYER
            DO L=1, N_PROFILE
              FRAC_REGION(L, I, IP_REGION_STRAT)                        &
     &          =FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SW)                     &
     &          +FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SI)
              FRAC_REGION(L, I, IP_REGION_CONV)                         &
     &          =FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CW)                     &
     &          +FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CI)
            ENDDO
          ENDDO
!
        ELSE IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CONV_STRAT) THEN
!
          DO K=1, N_CLOUD_TYPE
            IF (K == IP_CLOUD_TYPE_STRAT) THEN
              I_REGION_CLOUD(K)=IP_REGION_STRAT
            ELSE IF (K == IP_CLOUD_TYPE_CONV) THEN
              I_REGION_CLOUD(K)=IP_REGION_CONV
            ENDIF
          ENDDO
!
          DO I=N_CLOUD_TOP, N_LAYER
            DO L=1, N_PROFILE
              FRAC_REGION(L, I, IP_REGION_STRAT)                        &
     &          =FRAC_CLOUD(L, I, IP_CLOUD_TYPE_STRAT)
              FRAC_REGION(L, I, IP_REGION_CONV)                         &
     &          =FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CONV)
             ENDDO
          ENDDO
!
!
        ELSE
          WRITE(IU_ERR, '(/A, /A)')                                     &
     &      '*** Error: This representation of clouds is not '          &
     &      //'compatible with separate '                               &
     &      , 'convective and stratiform and overlap.'
          IERR=I_ERR_FATAL
          RETURN
        ENDIF
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE AGGREGATE_CLOUD
#endif
