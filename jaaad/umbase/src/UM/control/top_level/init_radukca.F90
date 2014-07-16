#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+  initialization of radiative feedback
!
! Subroutine Interface:
      SUBROUTINE INIT_RADUKCA(                                          &
#include "argd1.h"
#include "argptra.h"
     &      ngrgas,grgas_addr)

      IMPLICIT NONE
!
! Description:
! initialization of radiative feedback
!  Initialise and check address array to feed chemical tracers from
!  UKCA into radiation scheme.

! Current Code Owner: Bill Collins
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!
!
#include "cmaxsize.h"
#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
#include "typptra.h"
#include "nstypes.h"
#include "cruntimc.h"
#include "cntlatm.h"
#include "csubmodl.h"
#include "cprintst.h"

! include switches for minor gases
#include "lwopt3a.h"
#include "lwcopt3a.h"
#include "feedback.h"

! Subroutine arguments
!   Scalar arguments with Intent(In):
     INTEGER, INTENT(IN) :: ngrgas

!   Array  arguments with Intent(In):
!   Scalar arguments with Intent(InOut):
!   Array  arguments with Intent(InOut):
!   Scalar arguments with Intent(Out):
!   Array  arguments with Intent(Out):
      INTEGER, INTENT(OUT) :: grgas_addr(ngrgas)


! Local parameters:
      CHARACTER (Len=* ), Parameter :: RoutineName='Init_RadUKCA'

! STASH codes of tracers
      INTEGER, PARAMETER :: i_o3  =  1, i_ch4 = 9 , i_n2o = 49,         &
     &                      i_f11 = 55, i_f12 = 56,                     &
     &                      i_f113= 64, i_h22 = 65, i_cf3chf2 = -1,     &
     &                      i_chf2chf2 = -1, i_ro3 = 59

! Local scalars:
      INTEGER ::  I   ! Loop counter
      INTEGER :: section, m_atm_modl, address, item
      INTEGER, PARAMETER :: ukca_sect=34
      CHARACTER (Len=80) :: cmessage

!- End of header
      m_atm_modl = submodel_for_sm(a_im)
      DO i=1,no_obj_d1(m_atm_modl)
        section = d1_addr(d1_section, i, m_atm_modl)
        IF((d1_addr(d1_object_type,i,m_atm_modl)==prognostic).AND.  &
     &    (section == ukca_sect)) THEN
          item    = d1_addr(d1_item,    i, m_atm_modl)
          address = d1_addr(d1_address, i, m_atm_modl)
          SELECT CASE (item)
            CASE(i_o3)
              IF (L_ukca_rado3  ) grgas_addr(p_o3  )= address
            CASE(i_ch4)
              IF (L_ukca_radch4 ) grgas_addr(p_ch4 )= address
            CASE(i_n2o)
              IF (L_ukca_radn2o ) grgas_addr(p_n2o )= address
            CASE(i_f11)
              IF (L_ukca_radf11 ) grgas_addr(p_f11 )= address
            CASE(i_f12)
              IF (L_ukca_radf12 ) grgas_addr(p_f12 )= address
            CASE(i_f113)
              IF (L_ukca_radf113) grgas_addr(p_f113)= address
            CASE(i_h22)  
              IF (L_ukca_radf22 ) grgas_addr(p_f22 )= address
            CASE(i_ro3)
              IF (L_ukca_family .AND. L_ukca_rado3  )               &
     &          grgas_addr(p_o3  )= address
          END SELECT
        END IF
      END DO

! Turn grgas_addr into number of tracer in free_tracers array
      IF (ukca_sect == 33) THEN
        grgas_addr = (grgas_addr - jtracer(1,1)) /                  &
     &    (theta_off_size * tr_levels) + 1
      ELSEIF (ukca_sect == 34) THEN
        grgas_addr = (grgas_addr - jtr_ukca(1,1)) /                 &
     &    (theta_off_size * tr_levels) + 1
      ELSE
        cmessage=' Tracer section not identified'
! DEPENDS ON: ereport
        CALL ereport(RoutineName,1,cmessage)
      ENDIF

      IF(mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
        IF ((grgas_addr(p_o3) < 0) .AND. L_ukca_rado3) THEN
          cmessage='WARNING: Ox not found among chemical species.'
! DEPENDS ON: ereport
          CALL ereport(RoutineName,-1,cmessage)
        END IF

        IF ((grgas_addr(p_ch4) < 0) .AND. L_ukca_radch4) THEN
          cmessage='WARNING: CH4 not found among chemical species.'
! DEPENDS ON: ereport
          CALL ereport(RoutineName,-2,cmessage)
        END IF

        IF (ozone_levels /= model_levels) THEN
          cmessage = 'Ozone levels must equal model levels.'
! DEPENDS ON: ereport
          CALL ereport(RoutineName,1,cmessage)
        ENDIF

        IF (ozone_levels /= tr_levels) THEN
          cmessage='Tracer levels must equal model levels.'
! DEPENDS ON: ereport
          CALL ereport(RoutineName,2,cmessage)
        ENDIF

        IF ((grgas_addr(p_n2o) > 0) .AND. (.NOT. l_n2o_lw)) THEN 
          cmessage='N2O found but absorption by N2O not selected.'
! DEPENDS ON: ereport
          CALL ereport(RoutineName,-3,cmessage)
        END IF
        IF ((grgas_addr(p_f11) > 0) .AND. (.NOT. l_cfc11_lw)) THEN 
          cmessage='CFCl1 found but absorption by CFCl1 not selected.'
! DEPENDS ON: ereport
          CALL ereport(RoutineName,-4,cmessage)
        END IF

        IF ((grgas_addr(p_f12) > 0) .AND. (.NOT. l_cfc12_lw)) THEN 
          cmessage='CF2Cl2 found but absorption by CF2Cl2 not selected.'
! DEPENDS ON: ereport
          CALL ereport(RoutineName,-5,cmessage)
        END IF

        IF ((grgas_addr(p_f113) > 0) .AND. (.NOT. l_cfc113_lw)) THEN
          cmessage=                                                   &
     &    'CFC-113 found but absorption by CFC-113 not selected.'
! DEPENDS ON: ereport
          CALL ereport(RoutineName,-6,cmessage)
        END IF

        IF ((grgas_addr(p_f22) > 0) .AND. (.NOT. l_hcfc22_lw)) THEN
          cmessage=                                                   &
     &    'CHF2Cl found but absorption by CHF2Cl not selected.'
! DEPENDS ON: ereport
          CALL ereport(RoutineName,-7,cmessage)
        END IF

        IF ((grgas_addr(p_f11) > 0) .AND. (l_cfc11_lw) .AND.          &
     &    .NOT.(L_ukca_stratcfc))  THEN
          cmessage=                                                   &
     &    'L_ukca_stratcfc not selected, CFC-11 not treated properly.'               
! DEPENDS ON: ereport
          CALL ereport(RoutineName,-8,cmessage)
        END IF

        IF ((grgas_addr(p_f12) > 0) .AND. (l_cfc12_lw) .AND.          &
     &    .NOT.(L_ukca_stratcfc)) THEN
          cmessage=                                                   &
     &    'L_ukca_stratcfc not selected, CFC-12 not treated properly.'               
! DEPENDS ON: ereport
          CALL ereport(RoutineName,-9,cmessage)
        END IF
      END IF
      RETURN
      END SUBROUTINE INIT_RADUKCA
#endif
