#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set arrays describing the spherical truncation.
!
! Purpose:
!   This routine sets an arrays of pointers to control the allocation
!   of memory.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SET_TRUNCATION(IERR                                    &
     &  , I_TRUNCATION, LS_GLOBAL_TRUNC                                 &
     &  , LS_MAX_ORDER, LS_LOCAL_TRUNC                                  &
     &  , MS_MIN, MS_MAX, MS_TRUNC                                      &
     &  , IA_SPH_MM, N_ORDER_PHASE                                      &
     &  , ND_MAX_ORDER                                                  &
     &  )
!
!
!
      IMPLICIT NONE
!
!
      INTEGER, INTENT(IN) ::                                            &
     &    ND_MAX_ORDER
!           Size allocated for orders of spherical harmonics
!
!     Header files
#include "c_kinds.h"
#include "sph_truncation_pcf3z.h"
#include "error_pcf3z.h"
#include "def_std_io_icf3z.h"
!
!
!     Dummy arguments
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
      INTEGER, INTENT(IN) ::                                            &
     &    I_TRUNCATION                                                  &
!           Type of spherical truncation
     &  , LS_GLOBAL_TRUNC                                               &
!           Global order of truncation
     &  , MS_MIN                                                        &
!           Lowest value of of the azimuthal order calculated
     &  , MS_MAX
!           Highest value of of the azimuthal order calculated
      INTEGER, INTENT(OUT) ::                                           &
     &    LS_MAX_ORDER                                                  &
!           Maximum order of spherical harmonic terms required
     &  , LS_LOCAL_TRUNC(0: ND_MAX_ORDER)                               &
!           Truncating order for individual azimuthal quantum numbers
     &  , MS_TRUNC(0: ND_MAX_ORDER)                                     &
!           Maximum azimuthal quantum number for each order
     &  , IA_SPH_MM(0: ND_MAX_ORDER)                                    &
!           Address of spherical coefficients of (m, m) for each m
     &  , N_ORDER_PHASE
!           Order of terms in the phase function to be used in
!           direct calculation of spherical harmonics
!
!
!     Local Variables
      INTEGER                                                           &
     &    LS                                                            &
!           Order of spherical harmonic
     &  , MS
!           Azimuthal order of spherical harmonic
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    EVAL_UPLM
!
!
!
!     Carry out a preliminary check that the truncation is appropriate.
      IF ( (I_TRUNCATION == IP_TRUNC_AZIM_SYM).AND.                     &
     &   (MS_MAX >  0) ) THEN
!
        WRITE(IU_ERR, '(/A)')                                           &
     &    '*** Error: An azimuthally symmetric truncation is not '      &
     &    //'appropriate if MS_MAX > 0.'
        IERR=I_ERR_FATAL
        RETURN
!
      ENDIF
!
      IF ( (I_TRUNCATION == IP_TRUNC_TRIANGULAR).OR.                    &
     &   (I_TRUNCATION == IP_TRUNC_ADAPTIVE) ) THEN
!
!       Set the characteristics for a triangular truncation:
!       azimuthal orders from MS_MIN to MS_MAX are needed.
!       In order to keep an even number of harmonics for all
!       azimuthal orders, the maximum order must be set
!       1 greater than the (odd) order of the truncation. In
!       addition, an extra order beyond the truncation for the
!       particular m is required for the case of solar radiation.
!       ("+4" appears in the loop below because MS is one greater
!       then the order for which space is being reserved.) Note
!       finally that space is allocated for the unused harmonic
!       LS_GLOBAL_TRUNC+2 for even values of MS just to keep the
!       programming simple.
!
!       The adaptive truncation comes here as well since the maximum
!       conceivable number of harmonics might be required for each
!       azimuthal order.
        LS_MAX_ORDER=LS_GLOBAL_TRUNC+1
        MS_TRUNC(MS_MIN)=MS_MIN
        DO LS=MS_MIN+1, LS_MAX_ORDER
          MS_TRUNC(LS)=MIN(MIN(MS_MAX, LS), LS_GLOBAL_TRUNC)
        ENDDO
        IA_SPH_MM(MS_MIN)=1
        DO MS=MS_MIN+1, MS_MAX
          IA_SPH_MM(MS)=IA_SPH_MM(MS-1)+LS_GLOBAL_TRUNC+4-MS
        ENDDO
!       For each MS an even number of terms must be calculated. The
!       global truncation will be odd.
        DO MS=MS_MIN, MS_MAX
          LS_LOCAL_TRUNC(MS)=LS_GLOBAL_TRUNC+MOD(MS, 2)
        ENDDO
!
      ELSE IF (I_TRUNCATION == IP_TRUNC_RHOMBOHEDRAL) THEN
!
!       Set the characteristics for a rhombohedral truncation.
!       If calculation begins with an odd azimuthal order, one
!       extra order will be required to ensure that even of polar
!       orders are calculated.
        LS_MAX_ORDER=LS_GLOBAL_TRUNC+MOD(MS_MIN, 2)
!       N.B. If LS_MAX_ORDER is ever changed, make sure that the
!       following code is still valid. Here LS_MAX_ORDER logically
!       means LS_GLOBAL_TRUNC+MOD(MS_MIN, 2).
!
!       Reset the maximum azimuthal order if it has been set too high
        MS_TRUNC(MS_MIN)=MS_MIN
        DO LS=MS_MIN+1, LS_MAX_ORDER
          MS_TRUNC(LS)=MIN(MS_MAX, MIN(LS_MAX_ORDER-LS+MS_MIN, LS))
        ENDDO
        IA_SPH_MM(MS_MIN)=1
!       The "+4" rather than "+3" below allows for the fact that the
!       requisite number of harmonics does not fall off exactly
!       linearly with the azimuthal order, but does so in steps.
        DO MS=MS_MIN+1, MS_MAX
          IA_SPH_MM(MS)=IA_SPH_MM(MS-1)                                 &
     &      +LS_GLOBAL_TRUNC+4+MS_MIN-2*MS+1
        ENDDO
!       For each MS an even number of terms must be calculated. The
!       global truncation will be odd.
        DO MS=MS_MIN, MS_MAX
          LS_LOCAL_TRUNC(MS)=LS_GLOBAL_TRUNC                            &
     &      +MOD(MS_MIN, 2)-(MS-MS_MIN)
        ENDDO
!
      ELSE IF (I_TRUNCATION == IP_TRUNC_AZIM_SYM) THEN
!
!       Set the characteristics for an azimuthally symmetric truncation.
!       This will be the normal case in the infra-red region.
        LS_MAX_ORDER=LS_GLOBAL_TRUNC
        MS_TRUNC(0)=0
        IA_SPH_MM(0)=1
        DO LS=1, LS_MAX_ORDER
          MS_TRUNC(LS)=0
        ENDDO
!       Set the address of one extra order for use with solar sources.
        LS_LOCAL_TRUNC(0)=LS_GLOBAL_TRUNC
!
      ELSE
!
        WRITE(IU_ERR, '(/A)')                                           &
     &    '***Error: An illegal truncation has been requested.'
        IERR=I_ERR_FATAL
        RETURN
!
      ENDIF
!
!     Calculate enough terms of the phase function to satsify the
!     truncation at every azimuthal order.
      N_ORDER_PHASE=1
      DO MS=MS_MIN, MS_MAX
        N_ORDER_PHASE=MAX(N_ORDER_PHASE, LS_LOCAL_TRUNC(MS))
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SET_TRUNCATION
#endif
