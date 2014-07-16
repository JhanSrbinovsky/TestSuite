#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine INANCCTL
!LL
!LL Control routine for CRAY YMP
!LL
!LL Programming standard; Unified Model Documentation Paper No. 3
!LL draft version no. 3, dated 12/7/89
!LL
!LL System components covered : C710
!LL
!LL System task : C7
!LL
!LL Purpose : Takes as input,the code defining the frequency of update
!LL           of ancillary fields as set by the user interface.
!LL           Converts them into a list of numbers of timesteps after
!LL           which each field must be updated, and calculates the
!LL           frequency with which this list must be interogated.
!LL           Where the update interval is in months or years,
!LL           the check will be carried out each day. The physical
!LL           files required are also determined by input code,
!LL           and the headers and lookup tables are read into
!LL           COMMON/ANCILHDA/ or COMMON/ANCILHDO/ (*COMDECK CANCIL)
!LL
!LL Documentation : Unified Model Documentation Paper No C7
!LL                 Version No.4  dated 15/06/90
!LLEND
      SUBROUTINE INANCCTL(                                              &
#include "argduma.h"
#include "argsts.h"
#include "argptra.h"
#include "arganc.h"
#include "argppx.h"
     &                   ICODE,CMESSAGE)

!*
      IMPLICIT NONE
!*L Arguments
!L
#include "csubmodl.h"
#include "parvars.h"
#include "typsize.h"
#include "typduma.h"
#include "typsts.h"
#include "typptra.h"
#include "typanc.h"
#include "ppxlook.h"

!*L
      INTEGER                                                           &
     &        ICODE            ! Out return code :0 Nor al Exit
!                              !                 :>0 Error

      CHARACTER*80                                                      &
     &        CMESSAGE         ! Out error message if ICODE >0

!*

#include "chsunits.h"
#include "cmaxsize.h"
#include "ccontrol.h"
#include "clookadd.h"
#include "ctime.h"
#if defined(MPP)
#include "decomptp.h"
#endif


! Comdecks for ancillary files/fields.
#include "cancftna.h"
#include "canctita.h"

! Local variables
      INTEGER IM_IDENT       ! internal model identifier
      INTEGER IM_INDEX       ! internal model index for STASH arrays
      INTEGER IM_INDEX_S     ! as IM_INDEX for Slab Model
      INTEGER I
      INTEGER STEPS_PER_HR   ! steps per hour for atmos/wave sub_models
#if defined(MPP)
      INTEGER decomp_type   ! domain decomposition type
#endif


!L initialise reference time for time interpolation of ancillaries.
      IF ( ANCIL_REFTIME(1) == 0 .and. ANCIL_REFTIME(2) == 0            &
     &.and.ANCIL_REFTIME(3) == 0 .and. ANCIL_REFTIME(4) == 0            &
     &.and.ANCIL_REFTIME(5) == 0 .and. ANCIL_REFTIME(5) == 0 ) THEN
      WRITE(6,*)' ANCIL_REFTIME set same as MODEL_BASIS_TIME'
        ANCIL_REFTIME(1) = MODEL_BASIS_TIME(1)
        ANCIL_REFTIME(2) = MODEL_BASIS_TIME(2)
        ANCIL_REFTIME(3) = MODEL_BASIS_TIME(3)
        ANCIL_REFTIME(4) = MODEL_BASIS_TIME(4)
        ANCIL_REFTIME(5) = MODEL_BASIS_TIME(5)
        ANCIL_REFTIME(6) = MODEL_BASIS_TIME(6)
      WRITE(6,*)' ANCIL_REFTIME = MODEL_BASIS_TIME = ',ANCIL_REFTIME
      ELSE
      WRITE(6,*)' ANCIL_REFTIME set by User Interface = ',ANCIL_REFTIME
      END IF

#if defined(ATMOS)
!  Set up internal model identifier and STASH index
      im_ident = atmos_im
      im_index = internal_model_index(im_ident)
#if defined(MPP)
! Check that current decomposition is correct for ancillaries
! for this sub-model
      decomp_type=decomp_standard_atmos
      IF (current_decomp_type  /=  decomp_type) THEN
! DEPENDS ON: change_decomposition
        CALL CHANGE_DECOMPOSITION(decomp_type,icode)
        IF (icode  >   0) THEN
          WRITE(6,*) 'INANCCTL : Error'
          WRITE(6,*) 'Call to CHANGE_DECOMPOSITION failed with ',       &
     &               'decomposition type ',decomp_type
        CMESSAGE='INANCCTL;Unsupported decomposition for MPP code'
          GO TO 9999
        END IF
      END IF
#endif

!L Initialise Fortran file numbers
      DO I=1,NANCIL_DATASETSA
        FTNANCILA(I)=FTN_ANCIL_A(I)
      ENDDO
      STEPS_PER_HR = 3600*STEPS_PER_PERIODim(a_im)/                     &
     &                      SECS_PER_PERIODim(a_im)

! DEPENDS ON: inancila
      CALL INANCILA (LEN_FIXHD,PP_LEN_INTHD,PP_LEN_REALHD,A_LEN1_LEVDEPC&
     &              ,A_LEN2_LEVDEPC,FIXHD_ANCILA,INTHD_ANCILA,          &
     &              REALHD_ANCILA,LOOKUP_ANCILA,                        &
     &              A_FIXHD,A_REALHD,A_LEVDEPC,                         &
     &              NANCIL_DATASETSA,NANCIL_LOOKUPSA,FTNANCILA,         &
     &              LOOKUP_START_ANCILA,LEN1_LOOKUP,                    &
!  Following lines just to get valid compilation: entire routine
!  needs further development for new dynamics.
     &              glsize(1,fld_type_p),glsize(2,fld_type_p),          &
     &              glsize(2,fld_type_v),                               &
     &              MODEL_LEVELS,TR_LEVELS,                             &
     &              ST_LEVELS,SM_LEVELS,OZONE_LEVELS,tpps_ozone_levels, &
     &              TITLE_ANCIL_A,                                      &
!added entire SI array - kdcorbin, 05/10
     &              SI,nitems,nsects,n_internal_model,                  &
!                   all prognostics assumed to be in section 0
     &              SI(1,0,im_index),                                   &
                                          ! SI for atmos_im, sect 0
     &              SI(1,0,im_index),                                   &
                                          ! SI for atmos_im, sect 0
     &              NITEMS,                                             &
     &              ANCILLARY_STEPSim(a_im),STEPS_PER_HR,               &
#include "argppx.h"
     &              ICODE,CMESSAGE,LCAL360)

      IF(ICODE >  0) THEN
         write(6,*) 'INANCCTL: Error return from INANCILA ',ICODE
         write(6,*) CMESSAGE
         GO TO 9999            ! Jump to end of routine
      ENDIF

#endif

 9999 CONTINUE
      RETURN
      END SUBROUTINE INANCCTL
#endif
