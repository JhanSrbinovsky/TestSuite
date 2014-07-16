#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: SETGRCTL -------------------------------------------------
!LL
!LL  Purpose: Sets timestep group control switches.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Author:   T.C.Johns
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.1   8/02/93 : Changed order of comdecks to define NUNITS for
!LL                   comdeck CCONTROL.
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL   3.3  02/12/93  Generalise code for submodels and shared data
!LL                  partitions in coupled models (eg. SLAB).   (TCJ)
!LL
!LL   3.5  18/04/95  Stage 1 of submodel project: partial generalise
!LL                  to arbitrary submodels. R. Rawlins
!LL  4.1  17/04/96  Introduce wave sub-model.  RTHBarnes.
!LL  6.2  23/11/05  Removed all references to the wavemodel.
!LL                 T.Edwards
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: C0
!LL
!LL  Project task: C0
!LL
!LL  External documentation: On-line UM document C0 - The top-level
!LL                          control system
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE SETGRCTL (internal_model,submodel,NGROUP,              &
     &                     ICODE,CMESSAGE)
!
      IMPLICIT NONE
!
      INTEGER MODEL_DUMP_NUMBER(4)
      INTEGER internal_model  ! OUT - internal model id to run next
      INTEGER submodel        ! OUT - submodel id for dump partition
      INTEGER NGROUP          ! OUT - Number of steps in "group"
      INTEGER ICODE           ! Out - Return code
      CHARACTER*(80) CMESSAGE ! Out - Error message
!
!*----------------------------------------------------------------------
!  Common blocks
!
#include "cmaxsize.h"
#include "csubmodl.h"
#include "chsunits.h"
#include "chistory.h"
#include "ccontrol.h"
#include "ctime.h"
!
!  Local variables
! Temporary assignment to be replaced by node navigate at later stage
      INTEGER im      ! temporary internal model id for ocean or slab

!
!L----------------------------------------------------------------------
!L 1. Set timestep group control data using history file information,
!L    and model step numbers accumulated in CTIME
!L
!
! Hardwire settings follow, awaiting replacement by node navigation code
!
#if defined(ATMOS)
#if defined(OCEAN)
      im=ocean_im
#endif
#if defined(SLAB)
      im= slab_im
#endif
#if defined(OCEAN) || defined(SLAB)

! Check if ocean/slab has completed the same number of groups as atmos
      IF( (STEPim(atmos_im)/GROUPim(atmos_im) ) ==                      &
     &    (STEPim(      im)/GROUPim(      im) ) ) THEN
        internal_model=atmos_im
      ELSE
        internal_model=      im    ! either slab or ocean
      ENDIF
#else
      internal_model=atmos_im
#endif

#else
#if defined(OCEAN)
      internal_model=ocean_im
#else
      ICODE=1
      CMESSAGE="SETGRCTL : Illegal sub-model type, not ATMOS or OCEAN"
#endif
#endif

!!
!! 1.1 Determine if a new internal model or submodel for next group.
!!     {More generalisation needed later to cater for more complex
!!     coupling cases.}
!!

      IF(STEPim(internal_model) == 0 ) THEN ! Initial time, must be new
         new_im=.true.
         new_sm=.true.
      ELSE
         new_im=.false.
         new_sm=.false.
      ENDIF  ! Test for initial step

      IF( N_INTERNAL_MODEL >  1) THEN ! ie coupled model
! Check if end of group reached
         IF(mod(                                                        &
     &       STEPim(internal_model),GROUPim(internal_model)) == 0) THEN
             new_im=.true.            ! New internal model next

             IF( N_SUBMODEL_PARTITION >  1) THEN ! ie coupled submodels
                new_sm=.true.
             ENDIF                        ! Coupled submodel

         ENDIF                    ! Timestep at end of group
      ENDIF              ! Coupled model

!!
!! 1.2 Find submodel partition (ie D1/dump) identifier.
!!
      submodel=SUBMODEL_PARTITION_INDEX(internal_model)

!L   Find group of timesteps for next internal model
      NGROUP  = GROUPim(internal_model)

!L   Set switches as necessary for control variables held in
!L   CHISTORY and CCONTROL. {These are held over from 3.4 and
!L   could probably be removed, to be replaced by more generic items.}
      LATMOSNEXT=.FALSE.
      LOCEANNEXT=.FALSE.
      RUN_OCEAN_FIRST="N"
      IF(internal_model == atmos_im) THEN
          LATMOSNEXT=.TRUE.
      ELSEIF(internal_model == ocean_im) THEN
          LOCEANNEXT=.TRUE.
          RUN_OCEAN_FIRST="Y"
      ELSEIF(internal_model /= slab_im .and.                            &
     &       internal_model /= wave_im) THEN
          ICODE=1
          CMESSAGE=                                                     &
     &  "SETGRCTL : Illegal sub-model type, not ATMOS, OCEAN or WAVE"
          write(6,*) CMESSAGE
          write(6,*) 'illegal internal_model=',internal_model
      ENDIF

! diagnostic write start
      write(6,*) 'im,sm,ngroup,new_im,new_sm',                          &
     &  internal_model,submodel,ngroup,new_im,new_sm
! diagnostic write end
!
      RETURN
!L----------------------------------------------------------------------
      END SUBROUTINE SETGRCTL
#endif
