
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine GET_DIM
!LL
!LL  Written by D Robinson 28/7/92
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!    6.0     05/09/03  Add def for use with makebc. R. Sempers
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LL  Purpose : Get dimensions of data set components
!LL            from fixed length header.
!LL
!LLEND------------------------------------------------------------------
      SUBROUTINE GET_DIM (FIXHD,LEN_FIXHD,                              &
     &                    LEN_INTHD,LEN_REALHD,                         &
     &                    LEN1_LEVDEPC,LEN2_LEVDEPC,                    &
     &                    LEN1_ROWDEPC,LEN2_ROWDEPC,                    &
     &                    LEN1_COLDEPC,LEN2_COLDEPC,                    &
     &                    LEN1_FLDDEPC,LEN2_FLDDEPC,                    &
     &                    LEN_EXTCNST,LEN_DUMPHIST,                     &
     &                    LEN_CFI1,LEN_CFI2,LEN_CFI3,                   &
     &                    LEN1_LOOKUP,LEN2_LOOKUP,                      &
     &                    LEN_DATA)

!L----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER                                                           &
                       !  Dimension of :-
     &  LEN_FIXHD                                                       &
                       !  Fixed length header
     & ,LEN_INTHD                                                       &
                       !  Integer header
     & ,LEN_REALHD                                                      &
                       !  Real header
     & ,LEN1_LEVDEPC                                                    &
                       !  Level dependent constants (1st)
     & ,LEN2_LEVDEPC                                                    &
                       !  Level dependent constants (2nd)
     & ,LEN1_ROWDEPC                                                    &
                       !  Rows  dependent constants (1st)
     & ,LEN2_ROWDEPC                                                    &
                       !  Rows  dependent constants (2nd)
     & ,LEN1_COLDEPC                                                    &
                       !  Col   dependent constants (1st)
     & ,LEN2_COLDEPC                                                    &
                       !  Col   dependent constants (2nd)
     & ,LEN1_FLDDEPC                                                    &
                       !  Field dependent constants (1st)
     & ,LEN2_FLDDEPC                                                    &
                       !  Field dependent constants (2nd)
     & ,LEN_EXTCNST                                                     &
                       !  Extra constants
     & ,LEN_DUMPHIST                                                    &
                       !  Dump history
     & ,LEN_CFI1                                                        &
                       !  Compressed field index 1
     & ,LEN_CFI2                                                        &
                       !  Compressed field index 2
     & ,LEN_CFI3                                                        &
                       !  Compressed field index 3
     & ,LEN1_LOOKUP                                                     &
                       !  Look up table (1st)
     & ,LEN2_LOOKUP                                                     &
                       !  Look up table (2nd)
     & ,LEN_DATA       !  Data section

      INTEGER  FIXHD(LEN_FIXHD)   ! IN  Fixed length header

      LEN_INTHD    = FIXHD(101)
      LEN_REALHD   = FIXHD(106)
      LEN1_LEVDEPC = FIXHD(111)
      LEN2_LEVDEPC = FIXHD(112)
      LEN1_ROWDEPC = FIXHD(116)
      LEN2_ROWDEPC = FIXHD(117)
      LEN1_COLDEPC = FIXHD(121)
      LEN2_COLDEPC = FIXHD(122)
      LEN1_FLDDEPC = FIXHD(126)
      LEN2_FLDDEPC = FIXHD(127)
      LEN_EXTCNST  = FIXHD(131)
      LEN_DUMPHIST = FIXHD(136)
      LEN_CFI1     = FIXHD(141)
      LEN_CFI2     = FIXHD(143)
      LEN_CFI3     = FIXHD(145)
      LEN1_LOOKUP  = FIXHD(151)
      LEN2_LOOKUP  = FIXHD(152)
      LEN_DATA     = FIXHD(161)

      RETURN
      END SUBROUTINE GET_DIM
