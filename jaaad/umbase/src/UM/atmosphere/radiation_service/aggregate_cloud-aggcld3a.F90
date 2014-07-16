#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
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
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       HADAM3          05-06-96                Original Code
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
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
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     DUMMY ARRAY SIZES
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     INCLUDE COMDECKS
#include "stdio3a.h"
#include "error3a.h"
#include "dimfix3a.h"
#include "clrepp3a.h"
#include "cldtyp3a.h"
#include "clschm3a.h"
#include "cldreg3a.h"
!
!     DUMMY VARIABLES.
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , N_CLOUD_TOP
!             TOPMOST CLOUDY LAYER
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_CLOUD                                                      &
!             CLOUD SCHEME USED
     &   , I_CLOUD_REPRESENTATION                                       &
!             REPRESENTATION OF CLOUDS USED
     &   , N_CLOUD_TYPE
!             NUMBER OF TYPES OF CLOUD
!
      REAL                                                              &
                !, INTENT(OUT)
     &     FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             FRACTIONS OF EACH TYPE OF CLOUD
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     I_REGION_CLOUD(NPD_CLOUD_TYPE)
!             REGIONS IN WHICH PARTICULAR TYPES OF CLOUD FALL
      REAL                                                              &
                !, INTENT(OUT)
     &     FRAC_REGION(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             FRACTIONS OF TOTAL CLOUD OCCUPIED BY EACH REGION
!
!
!     LOCAL VARIABLES
      INTEGER                                                           &
     &     I                                                            &
!            LOOP VARIABLE
     &   , L                                                            &
!            LOOP VARIABLE
     &   , K
!            LOOP VARIABLE
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
                  FRAC_REGION(L, I, IP_REGION_STRAT)                    &
     &               =FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SW)                &
     &               +FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SI)
                  FRAC_REGION(L, I, IP_REGION_CONV)                     &
     &               =FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CW)                &
     &               +FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CI)
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
                  FRAC_REGION(L, I, IP_REGION_STRAT)                    &
     &               =FRAC_CLOUD(L, I, IP_CLOUD_TYPE_STRAT)
                  FRAC_REGION(L, I, IP_REGION_CONV)                     &
     &               =FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CONV)
               ENDDO
            ENDDO
!
!
         ELSE
            WRITE(IU_ERR, '(/A)')                                       &
     &         '*** ERROR: THIS REPRESENTATION OF CLOUDS IS NOT '       &
     &         //'COMPATIBLE WITH THE TRIPLE OVERLAP.'
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
#endif
