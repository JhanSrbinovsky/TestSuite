
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      Subroutine READWRITD(UNIT,ICODE,CMESSAGE)
      IMPLICIT NONE
!
! Purpose: Reads information for time-step control of WRITEDATA from
!          NAMLST file
!
! Current code owner: S.J.Swarbrick
!
!  Model               Modification history:
! version  Date
! -------  ----        --------------------
!   3.4    28/07/94    Original code - S.J.Swarbrick
!  3.5  08/06/95  Add UNIT - moved from INITIAL to READCNTL. RTHBarnes
!  4.1  31/05/96  Add L_WRIT_WAVSTEP for wave sub-model. RTHBarnes.
!  5.3  23/08/01  Namelists now in C_WRITD. Use "write data" rather
!                 than writd1 as the d1 reference is no longer relevant
!                 S.D.Mullerworth
!
! Code description:
!   FORTRAN 77 + common FORTRAN 90 extensions. Written to UM
!   programming standards vn. 7.
!
!
!
! This Comdeck declares and stores the logical and integer
! variables used in the time-step control for writing general
! data.
!
! Code owner: S.J.Swarbrick
!
! Switch which activates output of arrays for Peer utility
      LOGICAL L_PEER

! Switches which activate dump writing
      LOGICAL                                                           &
     &  L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                  &
     &  ,L_WRIT_INIT

! Timesteps for dump writing
      INTEGER                                                           &
     &  T_WRITD1_START                                                  &
                              ! First timestep
     &  ,T_WRITD1_END                                                   &
                              ! Last timestep
     &  ,T_WRITD1_INT         ! Timestep interval between dumps

      INTEGER                                                           &
     &  PEER_VN                  ! Version of PEER utility

      NAMELIST/NLSTWRITEDATA/                                           &
     &  L_PEER,PEER_VN                                                  &
     &  ,T_WRITD1_START  ,T_WRITD1_END  ,T_WRITD1_INT                   &
     &  ,L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                 &
     &  ,L_WRIT_INIT

      COMMON/WRITEDATA/                                                 &
     &  L_PEER,PEER_VN                                                  &
     &  ,T_WRITD1_START  ,T_WRITD1_END  ,T_WRITD1_INT                   &
     &  ,L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                 &
     &  ,L_WRIT_INIT
!
!
      CHARACTER*80 CMESSAGE  ! Error return message
!
      INTEGER ICODE          ! Indicator for normal or error return
      INTEGER UNIT           ! Namelist input unit no.
!
! --------------------------------------------------------------------
!
!
!  Read control variables for WRITE_DATA
!
      READ(UNIT,NLSTWRITEDATA,ERR=99)
!
!  Normal return
!
      ICODE=0
      RETURN
!
!  Error return
!
  99  ICODE=1
      CMESSAGE='READWRITD: error reading namelist NLSTWRITDATA'
      RETURN
      END SUBROUTINE READWRITD
