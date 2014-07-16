! ----------------------- Comdeck: CHISTG   ----------------------------
! Description: COMDECK defining Character History variables for
!              generic aspects of internal models
!              Generic means values likely to be common to the control
!              of any sub-model/internal model.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.4  30/05/97  Added vars LASTATMim, CURRATMim, LASTDMPim.  K Rogers
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      CHARACTER*14 END_DUMPim(N_INTERNAL_MODEL_MAX)!most recent dumpname
      CHARACTER*80 RESTARTim(N_INTERNAL_MODEL_MAX) !current restart dump
      CHARACTER*14 SAFEDMPim(N_INTERNAL_MODEL_MAX)
! Name of old safe restart dump
      CHARACTER*14 NEWSAFEim(N_INTERNAL_MODEL_MAX)
! Name of new safe restart dump
      CHARACTER*14 LASTATMim(N_INTERNAL_MODEL_MAX) ! Keep name of last
!                                                  ! atmos restart dump
!                                                  ! until ocean dump
      CHARACTER*14 CURRATMim(N_INTERNAL_MODEL_MAX) ! Keep name of
!                                                  ! current atmos
!                                                  ! restart dump
      CHARACTER*14 LASTDMPim(N_INTERNAL_MODEL_MAX) ! Keep name of last
!                                                  ! atmos/ocean dumps
!                                                  ! until meaning done

!
!
! History Common Block for generic model characters variables.
!
      COMMON /CHISTG/                                                   &
     &  END_DUMPim, RESTARTim,                                          &
     &  SAFEDMPim, NEWSAFEim, LASTATMim, CURRATMim, LASTDMPim

      NAMELIST /NLCHISTG/                                               &
     &  END_DUMPim, RESTARTim,                                          &
     &  SAFEDMPim, NEWSAFEim, LASTATMim, CURRATMim, LASTDMPim

! CHISTG end
