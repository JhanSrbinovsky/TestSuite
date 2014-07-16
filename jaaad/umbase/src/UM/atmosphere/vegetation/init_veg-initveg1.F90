#if defined(A19_1A) || defined(A19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Calls routines to initialize veg parameters and accumulated C fluxes
!
! Subroutine Interface:
      SUBROUTINE INIT_VEG(A_STEP,                                       &
#include "argd1.h"
#include "argduma.h"
#include "argptra.h"
#include "arglndm.h" 
     &              ICODE,CMESSAGE)

      IMPLICIT NONE
!
! Description:
!   Initializes vegetation parameters from fractions of surface types
!   and initializes accumulated carbon fluxes to zero if a new TRIFFID
!   calling period is starting.
!
! Method:
!   Calls routine SPARM to initialize vegetation parameters.
!   Calls routine INIT_ACC to initialize accumulated carbon fluxes.
!
! Current Code Owner:  Richard Betts
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   4.4   10/10/97   Original code.  Richard Betts
!   5.2   15/11/00   Re-Written for New Dynamics    M. Best
!   5.3   08/10/01   Bug fix for IF block to initialise G_LEAF_ACC
!                    to zero whenerver phenology is used     M. Best
!   5.4   28/08/02   Add arguements for canopy snow to call of
!                    SPARM.   R. Essery.
!   6.2   21/02/06   Bug fix for initialisation of G_LEAF_ACC, extra
!                    loop added to set all 5 pft's to zero  R. Harrison

!  6.2  01/03/06  replace JSOIL_CARB, JRSP_S_ACC with RothC pointers.
!                                                           C.D. Jones
!   6.2   21/03/06   Included nstypes.h J Ridley.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!

#include "nstypes.h"
#include "cmaxsize.h"
#include "csubmodl.h"
#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
#include "typduma.h"
#include "typlndm.h" 
#include "typptra.h"
#include "chsunits.h"
#include "ccontrol.h"
#include "cruntimc.h"
#include "ctime.h"

      INTEGER                                                           &
     & A_STEP             ! IN Current timestep in atmosphere model

      INTEGER                                                           &
     & TILE_PTS(NTYPE)                                                  &
                                    ! LOCAL Number of land points which
!                                   !       include the nth surface type
     &,TILE_INDEX(LAND_POINTS,NTYPE)                                    &
                                    ! LOCAL Indices of land points which
!                                   !       include the nth surface type
     &,NSTEP_TRIF                                                       &
                                    ! LOCAL Number of atmospheric
!                                   !       timesteps between calls to
!                                   !       TRIFFID.
     &,I                                                                &
                                    ! LOCAL Loop counter for all points
     &,L                                                                &
                                    ! LOCAL Loop counter for land points
     &,N          ! ** TEMPORARY ** loop counter for types

      INTEGER ICODE                 ! LOCAL Internal return code
      CHARACTER*80 CMESSAGE         ! LOCAL Internal error message

      EXTERNAL                                                          &
     & INIT_MIN                                                         &
     &,TILEPTS                                                          &
     &,SPARM                                                            &
     &,INIT_ACC

!-----------------------------------------------------------------------
! If TRIFFID on, call INIT_MIN to ensure PFT fractions are GE minimum
! fraction except where vegetation excluded by ice, water or urban
!-----------------------------------------------------------------------
      IF (L_TRIFFID) THEN
! DEPENDS ON: init_min
        CALL INIT_MIN(LAND_POINTS,D1(JFRAC_TYP),D1(JSOIL_CARB1))
      ENDIF

!-----------------------------------------------------------------------
! Call TILEPTS to initialise TILE_PTS and TILE_INDEX
!-----------------------------------------------------------------------
! DEPENDS ON: tilepts
      CALL TILEPTS(LAND_POINTS,D1(JFRAC_TYP),TILE_PTS,TILE_INDEX)

!-----------------------------------------------------------------------
! Initialise tiled and gridbox mean vegetation parameters
!-----------------------------------------------------------------------
      WRITE(6,*) 'INITVEG: CALLING SPARM'
! DEPENDS ON: sparm
      CALL SPARM (LAND_POINTS,NTILES,CAN_MODEL,TILE_PTS,TILE_INDEX,     &
     &            D1(JFRAC_TYP),D1(JCANHT_PFT),                         &
     &            D1(JLAI_PFT),D1(JSAT_SOIL_COND),D1(JCATCH_SNOW),      &
     &            D1(JCATCH_TILE),D1(JINFIL_TILE),D1(JZ0_TILE))

      IF (L_TRIFFID) THEN
!-----------------------------------------------------------------------
! If this is an NRUN and re-start from mid-way through a TRIFFID calling
! period has not been requested: (i) initialise accumulation prognostics
! to zero, (ii) set TRIFFID_PERIOD in integer header, and
! (iii) initialise ASTEPS_SINCE_TRIFFID integer header to zero.
! If mid-period restart is requested then leave the accumulated fields
! unchanged, and if a new calling period is specified then reset
! calling period header the new value provided that the number of
! atmosphere timesteps since the last call to TRIFFID does not exceed
! the new calling period .
! A_INTHD(22) holds TRIFFID_PERIOD in days.
! A_INTHD(23) holds the number of atmosphere timesteps since the last
! call to TRIFFID.
!-----------------------------------------------------------------------
        IF (A_STEP == 0) THEN
          IF (L_NRUN_MID_TRIF) THEN

            IF (TRIFFID_PERIOD /= A_INTHD(22)) THEN
              NSTEP_TRIF=INT(86400.0*TRIFFID_PERIOD/                    &
     &        SECS_PER_STEPim(atmos_im))

              IF (A_INTHD(23) >  NSTEP_TRIF) THEN
                WRITE(6,*) '**ERROR IN TRIFFID** YOU HAVE SELECTED TO'
                WRITE(6,*) 'START MID-WAY THROUGH A TRIFFID CALLING'
                WRITE(6,*) 'PERIOD BUT YOUR INITIAL DUMP CONTAINS'
                WRITE(6,*) 'PROGNOSTICS ACCUMULATED OVER A PERIOD'
                WRITE(6,*) 'LONGER THAN THE NEW CALLING PERIOD'

              ELSE

                A_INTHD(22)=TRIFFID_PERIOD

              ENDIF
            ENDIF

          ELSE

! DEPENDS ON: init_acc
            CALL INIT_ACC(LAND_POINTS,                                  &
     &                    D1(JNPP_PFT_ACC),                             &
     &                    D1(JG_PHLF_PFT_ACC),D1(JRSP_W_PFT_ACC),       &
     &                    D1(JRSP_S_ACC1),ICODE,CMESSAGE)

            A_INTHD(22)=TRIFFID_PERIOD
            A_INTHD(23)=0

          ENDIF
        ENDIF

      ENDIF

      IF (L_PHENOL) THEN

!-----------------------------------------------------------------------
! Initialise accumulated leaf turnover rate to zero
!-----------------------------------------------------------------------
       DO N=1,5

      DO L = 1,LAND_POINTS
       D1((JG_LF_PFT_ACC+L-1)+((N-1)*LAND_POINTS)) = 0.0
       ENDDO
       ENDDO


      ENDIF

      RETURN
      END SUBROUTINE INIT_VEG
#endif
