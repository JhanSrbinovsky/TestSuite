#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine UP_BOUND
!
! Purpose: Updates the boundary conditions in D1 with data from disk.
!          Atmosphere: On the first call, the orography LBC is read
!                      in, together with the LBCs at the start
!                      ( D1(JLBC) ) and end ( D1(JLBC_TEND) ) of the
!                      update period. (The only exception to this is
!                      when the LBC_TEND data is to come from a
!                      different boundary file. In this case, the data
!                      from the first file that would normally be read
!                      intoLBC is read into LBC_TEND instead. An extra
!                      call is then made to this routine in which the
!                      LBC_TEND data is copied into LBC, and data from
!                      the second boundary file is read into LBC_TEND.)
!                        At subsequent boundary updating steps, the
!                      value that was in D1(JLBC_TEND) is copied
!                      to D1(J_LBC) (this is now the LBC at the start
!                      of the new period) and the next record is read
!                      from disk into D1(J_LBC_TEND).
!          Other     : At the first step, two records are read from
!                      disk. The first record is stored in the dump,
!                      as is the tendency (the difference between the
!                      first and second records).
!                      At subsequent boundary updating sets, the next
!                      record is read, and the difference between this
!                      and the current LBC value is stored in the dump
!                      as the new tendency.
!
      SUBROUTINE UP_BOUND(I_AO,                                         &
#include "argd1.h"
#include "argduma.h"
#include "argptra.h"
#include "argbnd.h"
#include "argppx.h"
     &                  ICODE,CMESSAGE)

      IMPLICIT NONE

#include "cmaxsize.h"
#include "parparm.h"
#include "typsize.h"
#include "typd1.h"
#include "typduma.h"
#include "typptra.h"

#include "typbnd.h"

      INTEGER                                                           &
     &       I_AO,                                                      &
                           !  atmosphere/Ocean indicator
     &       ICODE         ! Error code = 0 Normal Exit
!                          !            > 0 Error Condition

      CHARACTER*(80)                                                    &
     &       CMESSAGE      ! Error message

!*
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "clookadd.h"
#include "ctime.h"
#include "cprintst.h"
#include "chsunits.h"
#include "cntlall.h"
#include "cntlatm.h"

! Local variables

      INTEGER                                                           &
     &       I,J                                                        &
     &,      NFTIN                                                      &
     &,      pretend_ts                                                 &
     &,      first_ts                                                   &
                                  ! first timestep
     &,      last_ts                                                    &
                                  ! last  timestep
     &,      steps_to_next_update                                       &
     &,      len_buf                                                    &
                                  ! length of buffer for readflds
     &,      item_in_file         ! stash item code
      LOGICAL                                                           &
     &       PERIODIC,                                                  &
                                ! True if periodic lateral boundary data
     &       Between_ALBC_files ! True if about to swap atmos bndy files

#if defined(ATMOS) && !defined(GLOBAL)
      Integer           :: steps_from_bdi_start ! Timesteps between
                                                ! start of bdi and
                                                ! start of run
      Logical, Save     :: first_atm_call = .true.
      Character (Len=4) :: ch_lbc_time
#endif

! --------------------------------------------------------------------

      ICODE=0
      CMESSAGE=' '

#if defined(ATMOS) && !defined(GLOBAL)

! 1.1   See whether the LBC and LBC_TEND data is to come from separate
!       boundary files, a situation that could arise if we're using two
!       boundary files. If this is the case, the LBC data from the first
!       file will initially be read into LBC_TEND. An extra call to this
!       routine will then be made after the second boundary file has
!       been opened. On this extra call, the LBC_TEND data will be
!       copied to LBC, and the LBC_TEND data will be replaced with data
!       from the second boundary file.

        ! Initialise:
        Between_ALBC_files = .FALSE.
        Steps_from_bdi_start = BNDARY_OFFSETim(a_im)

        IF (Num_ALBCS == 2)                                             &
     &    Between_ALBC_files = ALBC_num == 1 .AND.                      &
     &                         STEPim(a_im) >= ALBC_SwapStep

! 1.2 Read atmosphere lateral boundary field, first step.

      IF (first_atm_call .AND. I_AO == 1) THEN


        IF(BOUND_FIELDCODE(1) <= 0) THEN
          CMESSAGE= 'UP_BOUND: NO LBCs for atmosphere LAM'
          ICODE=1
          GOTO 9999
        END IF

        NFTIN=125

! 1.2.1 Read Orography data

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,1,LOOKUP_BOUNDA(1,1),LEN1_LOOKUP,         &
     &                D1(JOROG_LBC),                                    &
     &                global_LENRIMA(fld_type_p,halo_type_extended,     &
     &                               rima_type_orog),                   &
     &                FIXHD_BOUNDA(1,1),                                &
#include "argppx.h"
     &                ICODE,CMESSAGE)

         IF (ICODE  >   0) THEN
           WRITE(6,*) 'UP_BOUND : Problem in READFLDS reading ',        &
     &                'atmosphere orography'
           WRITE(6,*) 'ICODE: ',ICODE
           WRITE(6,*) 'CMESSAGE ',CMESSAGE
           GOTO 9999
         ENDIF

! 1.2.2 Update LOOKUP_BOUNDA with the correct information for the
!       current set of LOOKUP headers

        DO i=2,RIM_LOOKUPSA     ! Loop over lookup headers for the
                                ! first record (ignoring orog.)

          j=NBOUND_LOOKUP(1)+i-2  ! The "real" lookup header number
                                  ! from the LBC file

          LOOKUP_BOUNDA(LBYR,i)=LOOKUP_COMP_BOUNDA(LBCC_LBYR,j)
          LOOKUP_BOUNDA(LBMON,i)=LOOKUP_COMP_BOUNDA(LBCC_LBMON,j)
          LOOKUP_BOUNDA(LBDAT,i)=LOOKUP_COMP_BOUNDA(LBCC_LBDAT,j)
          LOOKUP_BOUNDA(LBHR,i)=LOOKUP_COMP_BOUNDA(LBCC_LBHR,j)
          LOOKUP_BOUNDA(LBMIN,i)=LOOKUP_COMP_BOUNDA(LBCC_LBMIN,j)
          LOOKUP_BOUNDA(LBDAY,i)=LOOKUP_COMP_BOUNDA(LBCC_LBDAY,j)
          LOOKUP_BOUNDA(LBEGIN,i)=LOOKUP_COMP_BOUNDA(LBCC_LBEGIN,j)
          LOOKUP_BOUNDA(NADDR,i)=LOOKUP_COMP_BOUNDA(LBCC_NADDR,j)

        ENDDO ! i

        IF (.NOT.Between_ALBC_files) THEN

! 1.2.3 Read in the LBCs for timestep 0

! DEPENDS ON: read_atmos_lbcs
          CALL READ_ATMOS_LBCS(                                         &
     &      LENRIMA(1,1,rima_type_norm),                                &
     &      global_LENRIMA(1,1,rima_type_norm),                         &
     &      MODEL_LEVELS,WET_LEVELS,TR_VARS,TR_LEVELS,                  &
     &      L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_pc2,               &
     &      L_mcr_qcf2_lbc, L_mcr_qrain_lbc,                            &
     &      L_mcr_qgraup_lbc, L_pc2_lbc, L_murk, L_murk_lbc,            &
     &      LEN1_LOOKUP, LEN_FIXHD, NFTIN,                              &
     &      RIM_LOOKUPSA-1, LOOKUP_BOUNDA(1,2), FIXHD_BOUNDA,           &
     &      D1(JU_LBC), D1(JV_LBC), D1(JW_LBC), D1(JRHO_LBC),           &
     &      D1(JTHETA_LBC), D1(JQ_LBC), D1(JQCL_LBC), D1(JQCF_LBC),     &
     &      D1(JQCF2_LBC), D1(JQRAIN_LBC), D1(JQGRAUP_LBC),             &
     &      D1(JCF_BULK_LBC), D1(JCF_LIQUID_LBC), D1(JCF_FROZEN_LBC),   &
     &      D1(JEXNER_LBC), D1(JU_ADV_LBC), D1(JV_ADV_LBC),             &
     &      D1(JW_ADV_LBC), D1(JMURK_LBC), D1(JTRACER_LBC(1)),          &
#include "argppx.h"
     &      ICODE,CMESSAGE)

          IF (ICODE  >   0) THEN
            WRITE(6,*) 'Problem with READ_ATMOS_LBCS reading initial ', &
     &                 'boundary data for timestep 0'
            WRITE(6,*) 'ICODE= ',ICODE
            WRITE(6,*) 'CMESSAGE= ',CMESSAGE
            GOTO 9999
          ENDIF

          If (PrintStatus >= PrStatus_Normal) Then
            write (ch_lbc_time(1:2),'(I2.2)') LOOKUP_BOUNDA(LBHR,2)
            write (ch_lbc_time(3:4),'(I2.2)') LOOKUP_BOUNDA(LBMIN,2)
            write (6,*)                                                 &
     &      'Up_Bound: Timestep ',STEPim(a_im),' : LBCs read in for ',  &
     &      ch_lbc_time,'Z ',LOOKUP_BOUNDA(LBDAT,2),'/',                &
     &      LOOKUP_BOUNDA(LBMON,2),'/',LOOKUP_BOUNDA(LBYR,2)
          End If

! 1.2.4  Increment lookup header to start of next update period

          NBOUND_LOOKUP(1)=NBOUND_LOOKUP(1)+(RIM_LOOKUPSA-1)

          DO i=2,RIM_LOOKUPSA     ! Loop over lookup headers for the
                                  ! first record (ignoring orog.)

            j=NBOUND_LOOKUP(1)+i-2  ! The "real" lookup header number
                                    ! from the LBC file

            LOOKUP_BOUNDA(LBYR,i)=LOOKUP_COMP_BOUNDA(LBCC_LBYR,j)
            LOOKUP_BOUNDA(LBMON,i)=LOOKUP_COMP_BOUNDA(LBCC_LBMON,j)
            LOOKUP_BOUNDA(LBDAT,i)=LOOKUP_COMP_BOUNDA(LBCC_LBDAT,j)
            LOOKUP_BOUNDA(LBHR,i)=LOOKUP_COMP_BOUNDA(LBCC_LBHR,j)
            LOOKUP_BOUNDA(LBMIN,i)=LOOKUP_COMP_BOUNDA(LBCC_LBMIN,j)
            LOOKUP_BOUNDA(LBDAY,i)=LOOKUP_COMP_BOUNDA(LBCC_LBDAY,j)
            LOOKUP_BOUNDA(LBEGIN,i)=LOOKUP_COMP_BOUNDA(LBCC_LBEGIN,j)
            LOOKUP_BOUNDA(NADDR,i)=LOOKUP_COMP_BOUNDA(LBCC_NADDR,j)

          ENDDO ! i

        END IF ! (.NOT.Between_ALBC_files)

! 1.2.5 Read in the LBC_TEND for the end of the update period

! DEPENDS ON: read_atmos_lbcs
        CALL READ_ATMOS_LBCS(                                           &
     &    LENRIMA(1,1,rima_type_norm),                                  &
     &    global_LENRIMA(1,1,rima_type_norm),                           &
     &    MODEL_LEVELS,WET_LEVELS,TR_VARS,TR_LEVELS,                    &
     &    L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_pc2,                 &
     &    L_mcr_qcf2_lbc, L_mcr_qrain_lbc,                              &
     &    L_mcr_qgraup_lbc, L_pc2_lbc, L_murk, L_murk_lbc,              &
     &    LEN1_LOOKUP, LEN_FIXHD, NFTIN,                                &
     &    RIM_LOOKUPSA-1, LOOKUP_BOUNDA(1,2), FIXHD_BOUNDA,             &
     &    D1(JU_LBC_TEND), D1(JV_LBC_TEND), D1(JW_LBC_TEND),            &
     &    D1(JRHO_LBC_TEND), D1(JTHETA_LBC_TEND), D1(JQ_LBC_TEND),      &
     &    D1(JQCL_LBC_TEND), D1(JQCF_LBC_TEND),                         &
     &    D1(JQCF2_LBC_TEND), D1(JQRAIN_LBC_TEND),                      &
     &    D1(JQGRAUP_LBC_TEND), D1(JCF_BULK_LBC_TEND),                  &
     &    D1(JCF_LIQUID_LBC_TEND), D1(JCF_FROZEN_LBC_TEND),             &
     &    D1(JEXNER_LBC_TEND),                                          &
     &    D1(JU_ADV_LBC_TEND), D1(JV_ADV_LBC_TEND),                     &
     &    D1(JW_ADV_LBC_TEND), D1(JMURK_LBC_TEND),                      &
     &    D1(JTRACER_LBC_TEND(1)),                                      &
#include "argppx.h"
     &    ICODE,CMESSAGE)

        IF (ICODE  >   0) THEN
          WRITE(6,*) 'Problem with READ_ATMOS_LBCS reading tendency ',  &
     &               'boundary data for timestep 0'
          WRITE(6,*) 'ICODE= ',ICODE
          WRITE(6,*) 'CMESSAGE= ',CMESSAGE
          GOTO 9999
        ENDIF

        If (PrintStatus >= PrStatus_Normal) Then
          write (ch_lbc_time(1:2),'(I2.2)') LOOKUP_BOUNDA(LBHR,2)
          write (ch_lbc_time(3:4),'(I2.2)') LOOKUP_BOUNDA(LBMIN,2)
          write (6,*)                                                   &
     &    'Up_Bound: Timestep ',STEPim(a_im),' : LBCs read in for ',    &
     &    ch_lbc_time,'Z ',LOOKUP_BOUNDA(LBDAT,2),'/',                  &
     &    LOOKUP_BOUNDA(LBMON,2),'/',LOOKUP_BOUNDA(LBYR,2)
        End If

! Increment the lookup header, ready for the next data

        NBOUND_LOOKUP(1)=NBOUND_LOOKUP(1)+(RIM_LOOKUPSA-1)

! Check to see if we need to increment the LBCs to the correct timestep
! (If between boundary files, this step will be done in an additional
! call to this routine.)

        IF (MOD(Steps_from_bdi_start+STEPim(a_im),RIM_STEPSA) /= 0      &
     &      .and. .NOT.Between_ALBC_files) THEN

          IF ((STEPim(a_im) < RIM_STEPSA) .AND.                         &
     &        (STEPim(a_im) /= 0))  THEN

            ! First timestep between two LBC update intervals
            ! We need to increment the LBC to the correct timestep

            first_ts = 0
            last_ts  = STEPim(a_im)-1

          ELSE

            ! The current timestep falls between two LBC update
            ! intervals (this is probably a continuation run).
            ! We need to increment the LBC to the correct timestep

            first_ts = 0
            last_ts  = MOD(Steps_from_bdi_start + STEPim(a_im),         &
     &                     RIM_STEPSA) - 1

          END IF

          DO pretend_ts = first_ts, last_ts

            steps_to_next_update=RIM_STEPSA-pretend_ts

! DEPENDS ON: increment_atmos_lbcs
            CALL INCREMENT_ATMOS_LBCS(                                  &
     &        steps_to_next_update,                                     &
     &        LENRIMA(1,1,rima_type_norm),                              &
     &        MODEL_LEVELS,WET_LEVELS,TR_VARS,TR_LEVELS,                &
     &        L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc,        &
     &        L_pc2_lbc, L_murk_lbc, L_int_uvw_lbc,                     &
     &        D1(JU_LBC),D1(JU_LBC_TEND),                               &
     &        D1(JV_LBC),D1(JV_LBC_TEND),                               &
     &        D1(JW_LBC),D1(JW_LBC_TEND),                               &
     &        D1(JRHO_LBC),D1(JRHO_LBC_TEND),                           &
     &        D1(JTHETA_LBC),D1(JTHETA_LBC_TEND),                       &
     &        D1(JQ_LBC),D1(JQ_LBC_TEND),                               &
     &        D1(JQCL_LBC),D1(JQCL_LBC_TEND),                           &
     &        D1(JQCF_LBC),D1(JQCF_LBC_TEND),                           &
     &        D1(JQCF2_LBC),D1(JQCF2_LBC_TEND),                         &
     &        D1(JQRAIN_LBC),D1(JQRAIN_LBC_TEND),                       &
     &        D1(JQGRAUP_LBC),D1(JQGRAUP_LBC_TEND),                     &
     &        D1(JCF_BULK_LBC  ),D1(JCF_BULK_LBC_TEND  ),               &
     &        D1(JCF_LIQUID_LBC),D1(JCF_LIQUID_LBC_TEND),               &
     &        D1(JCF_FROZEN_LBC),D1(JCF_FROZEN_LBC_TEND),               &
     &        D1(JEXNER_LBC),D1(JEXNER_LBC_TEND),                       &
     &        D1(JU_ADV_LBC),D1(JU_ADV_LBC_TEND),                       &
     &        D1(JV_ADV_LBC),D1(JV_ADV_LBC_TEND),                       &
     &        D1(JW_ADV_LBC),D1(JW_ADV_LBC_TEND),                       &
     &        D1(JMURK_LBC),D1(JMURK_LBC_TEND),                         &
     &        D1(JTRACER_LBC(1)),D1(JTRACER_LBC_TEND(1)),               &
     &        RIM_STEPSA,                                               &
     &        ICODE,CMESSAGE)

            IF (ICODE  >   0) THEN
              WRITE(6,*) 'Failure in INCREMENT_ATMOS_LBCS while ',      &
     &                   'attempting to set LBCS for first timestep'
              GOTO 9999
            ENDIF

          ENDDO  ! pretend_ts

        ENDIF  ! IF (MOD(STEPim(a_im),RIM_STEPSA)  /=  0 .AND.
               !     .NOT.Between_ALBC_files)


      ENDIF  ! (first_atm_call .AND. I_AO == 1)

#endif


#if defined(ATMOS) && !defined(GLOBAL)


! 2.1 Read atmosphere lateral boundary fields, general update step


      IF (.NOT.first_atm_call .AND. I_AO == 1) THEN

        NFTIN=125

        IF (BOUND_FIELDCODE(1) >  0 .AND.                               &
                                                     ! If this is an
     &      MOD(steps_from_bdi_start + STEPim(a_im),                    &
                                                     ! LBC update step
     &          RIM_STEPSA) == 0 ) THEN

          IF (NBOUND_LOOKUP(1)  >=  FIXHD_BOUNDA(152,1)) THEN
            ICODE=11
            CMESSAGE='UP_BOUND : Reached end of atmosphere LBC file'
            GOTO 9999
          ENDIF

! 2.1.1 Copy data from LBC_TEND to LBC - this is the new LBC data
!       starting point for the new LBC update period

! DEPENDS ON: copy_atmos_lbcs
          CALL COPY_ATMOS_LBCS(                                         &
     &      LENRIMA(1,1,rima_type_norm),                                &
     &      MODEL_LEVELS,WET_LEVELS,TR_VARS,TR_LEVELS,                  &
     &      L_mcr_qcf2_lbc, L_mcr_qrain_lbc,                            &
     &      L_mcr_qgraup_lbc, L_pc2_lbc, L_murk_lbc,                    &
     &      D1(JU_LBC), D1(JU_LBC_TEND),                                &
     &      D1(JV_LBC), D1(JV_LBC_TEND),                                &
     &      D1(JW_LBC), D1(JW_LBC_TEND),                                &
     &      D1(JRHO_LBC), D1(JRHO_LBC_TEND),                            &
     &      D1(JTHETA_LBC), D1(JTHETA_LBC_TEND),                        &
     &      D1(JQ_LBC), D1(JQ_LBC_TEND),                                &
     &      D1(JQCL_LBC), D1(JQCL_LBC_TEND),                            &
     &      D1(JQCF_LBC), D1(JQCF_LBC_TEND),                            &
     &      D1(JQCF2_LBC), D1(JQCF2_LBC_TEND),                          &
     &      D1(JQRAIN_LBC), D1(JQRAIN_LBC_TEND),                        &
     &      D1(JQGRAUP_LBC), D1(JQGRAUP_LBC_TEND),                      &
     &      D1(JCF_BULK_LBC  ), D1(JCF_BULK_LBC_TEND  ),                &
     &      D1(JCF_LIQUID_LBC), D1(JCF_LIQUID_LBC_TEND),                &
     &      D1(JCF_FROZEN_LBC), D1(JCF_FROZEN_LBC_TEND),                &
     &      D1(JEXNER_LBC), D1(JEXNER_LBC_TEND),                        &
     &      D1(JU_ADV_LBC), D1(JU_ADV_LBC_TEND),                        &
     &      D1(JV_ADV_LBC), D1(JV_ADV_LBC_TEND),                        &
     &      D1(JW_ADV_LBC), D1(JW_ADV_LBC_TEND),                        &
     &      D1(JMURK_LBC), D1(JMURK_LBC_TEND),                          &
     &      D1(JTRACER_LBC(1)), D1(JTRACER_LBC_TEND(1))                 &
     &      )

! 2.1.2 Update LOOKUP_BOUNDA with the correct information for the
!       current set of LOOKUP headers

          DO i=2,RIM_LOOKUPSA     ! Loop over lookup headers for the
                                  ! first record (ignoring orog.)

            j=NBOUND_LOOKUP(1)+i-2  ! The "real" lookup header number
                                    ! from the LBC file

            LOOKUP_BOUNDA(LBYR,i)=LOOKUP_COMP_BOUNDA(LBCC_LBYR,j)
            LOOKUP_BOUNDA(LBMON,i)=LOOKUP_COMP_BOUNDA(LBCC_LBMON,j)
            LOOKUP_BOUNDA(LBDAT,i)=LOOKUP_COMP_BOUNDA(LBCC_LBDAT,j)
            LOOKUP_BOUNDA(LBHR,i)=LOOKUP_COMP_BOUNDA(LBCC_LBHR,j)
            LOOKUP_BOUNDA(LBMIN,i)=LOOKUP_COMP_BOUNDA(LBCC_LBMIN,j)
            LOOKUP_BOUNDA(LBDAY,i)=LOOKUP_COMP_BOUNDA(LBCC_LBDAY,j)
            LOOKUP_BOUNDA(LBEGIN,i)=LOOKUP_COMP_BOUNDA(LBCC_LBEGIN,j)
            LOOKUP_BOUNDA(NADDR,i)=LOOKUP_COMP_BOUNDA(LBCC_NADDR,j)

          ENDDO ! i

! 2.1.3 Read in the new set of LBC_TENDs

! DEPENDS ON: read_atmos_lbcs
          CALL READ_ATMOS_LBCS(                                         &
     &      LENRIMA(1,1,rima_type_norm),                                &
     &      global_LENRIMA(1,1,rima_type_norm),                         &
     &      MODEL_LEVELS,WET_LEVELS,TR_VARS,TR_LEVELS,                  &
     &      L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_pc2,               &
     &      L_mcr_qcf2_lbc, L_mcr_qrain_lbc,                            &
     &      L_mcr_qgraup_lbc, L_pc2_lbc, L_murk, L_murk_lbc,            &
     &      LEN1_LOOKUP, LEN_FIXHD, NFTIN,                              &
     &      RIM_LOOKUPSA-1, LOOKUP_BOUNDA(1,2), FIXHD_BOUNDA,           &
     &      D1(JU_LBC_TEND), D1(JV_LBC_TEND), D1(JW_LBC_TEND),          &
     &      D1(JRHO_LBC_TEND), D1(JTHETA_LBC_TEND), D1(JQ_LBC_TEND),    &
     &      D1(JQCL_LBC_TEND), D1(JQCF_LBC_TEND),                       &
     &      D1(JQCF2_LBC_TEND), D1(JQRAIN_LBC_TEND),                    &
     &      D1(JQGRAUP_LBC_TEND),D1(JCF_BULK_LBC_TEND),                 &
     &      D1(JCF_LIQUID_LBC_TEND), D1(JCF_FROZEN_LBC_TEND),           &
     &      D1(JEXNER_LBC_TEND),                                        &
     &      D1(JU_ADV_LBC_TEND), D1(JV_ADV_LBC_TEND),                   &
     &      D1(JW_ADV_LBC_TEND), D1(JMURK_LBC_TEND),                    &
     &      D1(JTRACER_LBC_TEND(1)),                                    &
#include "argppx.h"
     &      ICODE,CMESSAGE)

          IF (ICODE  >   0) THEN
            WRITE(6,*) 'Problem with READ_ATMOS_LBCS reading ',         &
     &                 'LBC_TENDs for general update step.'
            WRITE(6,*) 'ICODE= ',ICODE
            WRITE(6,*) 'CMESSAGE= ',CMESSAGE
            GOTO 9999
          ENDIF

          If (PrintStatus >= PrStatus_Normal) Then
            write (ch_lbc_time(1:2),'(I2.2)') LOOKUP_BOUNDA(LBHR,2)
            write (ch_lbc_time(3:4),'(I2.2)') LOOKUP_BOUNDA(LBMIN,2)
            write (6,*)                                                 &
     &      'Up_Bound: Timestep ',STEPim(a_im),' : LBCs read in for ',  &
     &      ch_lbc_time,'Z ',LOOKUP_BOUNDA(LBDAT,2),'/',                &
     &      LOOKUP_BOUNDA(LBMON,2),'/',LOOKUP_BOUNDA(LBYR,2)
          End If

          IF (MOD(STEPim(a_im),RIM_STEPSA)  /=  0 .AND.                 &
     &        Between_ALBC_files) THEN

            ! The current timestep falls between two LBC update
            ! intervals (this is probably a continuation run).
            ! We need to increment the LBC to the correct timestep

            DO pretend_ts=0,MOD(STEPim(a_im),RIM_STEPSA)-1

              steps_to_next_update=RIM_STEPSA-pretend_ts

! DEPENDS ON: increment_atmos_lbcs
              CALL INCREMENT_ATMOS_LBCS(                                &
     &          steps_to_next_update,                                   &
     &          LENRIMA(1,1,rima_type_norm),                            &
     &          MODEL_LEVELS,WET_LEVELS,TR_VARS,TR_LEVELS,              &
     &          L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc,      &
     &          L_pc2_lbc, L_murk_lbc, L_int_uvw_lbc,                   &
     &          D1(JU_LBC),D1(JU_LBC_TEND),                             &
     &          D1(JV_LBC),D1(JV_LBC_TEND),                             &
     &          D1(JW_LBC),D1(JW_LBC_TEND),                             &
     &          D1(JRHO_LBC),D1(JRHO_LBC_TEND),                         &
     &          D1(JTHETA_LBC),D1(JTHETA_LBC_TEND),                     &
     &          D1(JQ_LBC),D1(JQ_LBC_TEND),                             &
     &          D1(JQCL_LBC),D1(JQCL_LBC_TEND),                         &
     &          D1(JQCF_LBC),D1(JQCF_LBC_TEND),                         &
     &          D1(JQCF2_LBC),D1(JQCF2_LBC_TEND),                       &
     &          D1(JQRAIN_LBC),D1(JQRAIN_LBC_TEND),                     &
     &          D1(JQGRAUP_LBC),D1(JQGRAUP_LBC_TEND),                   &
     &          D1(JCF_BULK_LBC  ),D1(JCF_BULK_LBC_TEND  ),             &
     &          D1(JCF_LIQUID_LBC),D1(JCF_LIQUID_LBC_TEND),             &
     &          D1(JCF_FROZEN_LBC),D1(JCF_FROZEN_LBC_TEND),             &
     &          D1(JEXNER_LBC),D1(JEXNER_LBC_TEND),                     &
     &          D1(JU_ADV_LBC),D1(JU_ADV_LBC_TEND),                     &
     &          D1(JV_ADV_LBC),D1(JV_ADV_LBC_TEND),                     &
     &          D1(JW_ADV_LBC),D1(JW_ADV_LBC_TEND),                     &
     &          D1(JMURK_LBC), D1(JMURK_LBC_TEND),                      &
     &          D1(JTRACER_LBC(1)),D1(JTRACER_LBC_TEND(1)),             &
     &          RIM_STEPSA,                                             &
     &          ICODE,CMESSAGE)

              IF (ICODE  >   0) THEN
                WRITE(6,*) 'Failure in INCREMENT_ATMOS_LBCS while ',    &
     &                     'attempting to set LBCS for first timestep'
                GOTO 9999
              ENDIF

            ENDDO  ! pretend_ts

          ENDIF  ! IF (MOD(STEPim(a_im),RIM_STEPSA)  /=  0 .AND.
                 !     Between_ALBC_files)

! 2.1.4  Increment lookup header to start of next update period

          NBOUND_LOOKUP(1)=NBOUND_LOOKUP(1)+(RIM_LOOKUPSA-1)

        ENDIF ! IF it's an LBC update step

      ELSE

        first_atm_call = .false.

      ENDIF !  (.NOT.first_atm_call .AND. I_AO == 1)

#endif

 9999 CONTINUE
      RETURN
      END SUBROUTINE UP_BOUND

#endif
