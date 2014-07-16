#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine UP_ANCIL
!LL
!LL CW, SI      <- programmer of some or all of previous code or changes
!LL
!LL Programing standard : UM documentation paper no3,
!LL                       version no1, dated 15/01/90
!LL
!LL System components covered : C71
!LL
!LL System task C7
!LL
!LL Purpose: The routine is entered when any of the ancillary
!LL         fields have to be updated. The list of fields is
!LL         searched for update requirements. The position of
!LL         the data required is found from the header information
!LL         read in from subroutine IN_ANCIL. The data is read in
!LL         and updates the existing information.
!LL
!LL Documentation: Unified Model documentation paper no C7.
!LL                version no 4, dated 15/06/90
!LL
!LLEND

      SUBROUTINE UP_ANCIL(                                              &
#include "argd1.h"
#include "argduma.h"
#include "argptra.h"
#include "arganc.h"
     &                    I_AO,                                         &
#include "argppx.h"
     &                    ICODE,CMESSAGE)

      IMPLICIT NONE

!*L Arguments
!L
#include "parparm.h"
#include "typsize.h"
#include "typd1.h"
#include "typduma.h"
#include "typptra.h"
#include "typanc.h"
      INTEGER                                                           &
     &       I_AO,                                                      &
                       ! Sub-model indicator = 1 Atmosphere
     &       ICODE     ! Return code

      CHARACTER*80                                                      &
     &       CMESSAGE  ! Error message
      CHARACTER (Len=*), Parameter :: RoutineName = 'UP_ANCIL'

! Local Storage
      INTEGER   ANCIL_REF_DAYS,ANCIL_REF_SECS
      INTEGER   ANCIL_OFFSET_STEPS ! offset of ref. from basis time
!*
! Include COMDECKS
#include "chsunits.h"
#include "cmaxsize.h"
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "ccontrol.h"
#include "ctime.h"

      LOGICAL SMC_UPDATED  ! T if soil moisture updated

! Local Storage
      INTEGER   A_STEPS_PER_HR

      REAL DZ_SOIL(SM_LEVELS)  ! soil level thicknesses

      ICODE=0
      CMESSAGE=' '

!L  Convert ancillary reference time to days & secs
! DEPENDS ON: time2sec
        CALL TIME2SEC(ANCIL_REFTIME(1),ANCIL_REFTIME(2),                &
     &                ANCIL_REFTIME(3),ANCIL_REFTIME(4),                &
     &                ANCIL_REFTIME(5),ANCIL_REFTIME(6),                &
     &                    0,0,ANCIL_REF_DAYS,ANCIL_REF_SECS,LCAL360)

!L  Compute offset in timesteps of basis time from ancillary ref.time
! DEPENDS ON: tim2step
        CALL TIM2STEP(BASIS_TIME_DAYS-ANCIL_REF_DAYS,                   &
     &                BASIS_TIME_SECS-ANCIL_REF_SECS,                   &
     &                STEPS_PER_PERIODim(I_AO),SECS_PER_PERIODim(I_AO), &
     &                ANCIL_OFFSET_STEPS)

      IF(I_AO == 1) THEN
! Compute A_STEPS_PER_HR for use in REPLANCA
      A_STEPS_PER_HR = 3600*STEPS_PER_PERIODim(a_im)/                   &
     &                       SECS_PER_PERIODim(a_im)

#if defined(ATMOS)

! DEPENDS ON: replanca
        CALL REPLANCA(I_YEAR,I_MONTH,I_DAY,I_HOUR,I_MINUTE,I_SECOND,    &
     &                I_DAY_NUMBER,ANCIL_REFTIME,ANCIL_OFFSET_STEPS,    &
     &                THETA_FIELD_SIZE,ROWS,U_FIELD_SIZE,               &
     &                V_FIELD_SIZE,                                     &
     &                D1,D1(JLAND),                                     &
     &                STEPim(I_AO),LAND_FIELD,A_STEPS_PER_HR,           &
     &                D1(JFRAC_LAND),                                   &
     &                D1(JTSTAR_LAND),D1(JTSTAR_SEA),                   &
     &                D1(JTSTAR_SICE),                                  &
     &                D1(JICE_FRACTION),D1(JTSTAR),                     &
     &                D1(JTSTAR_ANOM),SM_LEVELS,                        &
     &                DZ_SOIL,SMC_UPDATED,                              &
     &                A_REALHD(2),A_REALHD(3),                          &
     &                LEN1_LOOKUP,LEN_FIXHD,                            &
     &                PP_LEN_INTHD,PP_LEN_REALHD,LEN_TOT,               &
     &                FIXHD_ANCILA,INTHD_ANCILA,REALHD_ANCILA,          &
     &                LOOKUP_ANCILA,LOOKUP_ANCILA,                      &
     &                FTNANCILA,LOOKUP_START_ANCILA,                    &
     &                NANCIL_DATASETSA,NANCIL_LOOKUPSA,                 &
#include "argppx.h"
     &                ICODE,CMESSAGE,LCAL360)

        IF (ICODE  /=  0) THEN
          Write (6,*) ' UP_ANCIL : Error in REPLANCA.'
! DEPENDS ON: ereport
          CALL Ereport (RoutineName,ICODE,Cmessage)
        ENDIF

! Update vegetation parameters
! DEPENDS ON: update_veg
         CALL UPDATE_VEG(                                               &
#include "argd1.h"
     &                   I_AO                                           &
     &                   )

! Update partitioning of unfrozen and frozen soil moisture
! if soil moisture updated
      IF(SMC_UPDATED) THEN
! DEPENDS ON: update_smc
         CALL UPDATE_SMC(                                               &
#include "argd1.h"
#include "argptra.h"
     &                   DZ_SOIL)
      END IF
#endif
      END IF

      RETURN
      END SUBROUTINE UP_ANCIL

!-----------------------------------------------------------------------

#endif
