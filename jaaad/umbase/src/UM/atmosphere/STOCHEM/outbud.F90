#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE OUTBUD(nm,mass,m0,tflux,clist,nchem,fnames,            &
     &  nflux,totm0,daym0,totmas,daymas,totavg,totflu,navg,flist,       &
     &  stash_code,month,year,out,fixhd12,period,umstepno,              &
     &  z_top_of_model,first_constant_r_rho_level)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Output Budgets, Inventories
!-                         and Concentration fields
!-
!-   Inputs  : NM,MASS,M0,TFLUX,CLIST,NCHEM,FNAMES,NFLUX,FLIST,TOTM0,
!-             TOTMAS, TOTAVG,NAVG,MCONC,DAYM0,DAYMAS
!-   Outputs : NONE
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  3.4    09/12/93  Created.  W.J. Collins
!  4.5    14/05/99  Now outputs 3D fluxes in PP format. C.E. Johnson
!  6.0    12/09/03  Fix problems with incorrect WRITE syntax. Introduce
!                   standard UM modification  history.  P.Dando
!  6.1    21/10/04  No change.
!  6.2    02/03/06  Now outputs stash code read in from variables
!                   in module.
!
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_OUT
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER,                             INTENT(IN) :: navg
      INTEGER,                             INTENT(IN) :: nflux
      INTEGER,                             INTENT(IN) :: nchem
      INTEGER,                             INTENT(IN) :: month
      INTEGER,                             INTENT(IN) :: year
      INTEGER,                             INTENT(IN) :: fixhd12
      INTEGER, DIMENSION(nlnpe,nlpe,nlev), INTENT(IN) :: nm
      INTEGER, DIMENSION(2,numflux),       INTENT(IN) :: flist
      INTEGER, DIMENSION(numchem),         INTENT(IN) :: clist
      INTEGER, DIMENSION(numflux),         INTENT(IN) :: stash_code
      INTEGER,                             INTENT(IN) ::                &
     &                                      first_constant_r_rho_level
      INTEGER,                             INTENT(IN) :: umstepno

      REAL,                                INTENT(IN) :: daym0
      REAL,                                INTENT(IN) :: daymas
      REAL,                                INTENT(IN) :: period
      REAL,                                INTENT(IN) :: z_top_of_model
      REAL, DIMENSION(2,nc),               INTENT(IN) :: totm0
      REAL, DIMENSION(2,nc),               INTENT(IN) :: totmas
      REAL, DIMENSION(2,nc),               INTENT(IN) :: totavg
      REAL, DIMENSION(numflux),            INTENT(IN) :: totflu
      REAL, DIMENSION(numchem,nlnpe,nlpe,nlev), INTENT(IN) :: mass
      REAL, DIMENSION(numchem,nlnpe,nlpe,nlev), INTENT(IN) :: m0
      REAL, DIMENSION(num3dflux_dim,nlnpe,nlpe,nlev),INTENT(IN) :: tflux
      REAL, DIMENSION(numchem,nlnpe,nlpe,nlev), INTENT(IN) :: out

      CHARACTER(LEN=40), DIMENSION(numflux), INTENT(IN) :: fnames

      INTEGER :: i
      INTEGER :: j
      INTEGER :: l
      INTEGER :: ll
      INTEGER :: len_io
      INTEGER :: len_in
      INTEGER :: Word_Address
      INTEGER :: nhdr
      INTEGER :: Nfield
      INTEGER :: icode
      INTEGER :: ihour
      INTEGER :: iday
      INTEGER, DIMENSION(64,nlev*(num3dflux+1+numchem)) :: pplook
      INTEGER, DIMENSION(256)                           :: flheader

      REAL    :: a_io
      REAL, DIMENSION(nlong,mnlat,nlev)                 :: global

      CHARACTER(LEN=72)                                 :: cmessage
      CHARACTER(LEN=14)                                 :: filename
      CHARACTER(LEN=46)                                 :: cstring

#include "cntl_io.h"

      COMMON /SUM3_D/ global

      iday=INT(daymas)
      ihour=INT((daymas-iday)*24)
      IF(mype==0) THEN
! Output inventories.
! DEPENDS ON: findname
        CALL FINDNAME('t',filetype2,'c',0,0,filename,                   &
     &     month,year)
        OPEN(18,FILE=filename,STATUS='UNKNOWN')
        WRITE(18,*) ' Inventory and Budget file from day ',INT(daym0),  &
     &    ' to ',INT(daymas),' for ',expt_id,job_id
        WRITE(18,'(A,I3,A,I5)') ' MONTH: ',month,' YEAR: ',year
        WRITE(18,*) ' Accumulated Time was: ',period,' days.'
        WRITE(18,*) ' NCHEM: ',nc,' NAVG:  ',navg
        WRITE(18,'(A,A,F7.3,/)') ' TOTAL AND TROPOSPHERIC INVENTORIES', &
     &    ' IN MOLES AT DAY:', daym0
        WRITE(18,'(A,'': '',2E20.8)') (cnames(l),                       &
     &    totm0(1,l)*lmolec/na, totm0(2,l)*lmolec/na, l=1,nc)
        WRITE(18,*)
        WRITE(18,'(A,A,F7.3,/)') ' TOTAL AND TROPOSPHERIC INVENTORIES', &
     &    ' IN MOLES AT DAY:', DAYMAS
        WRITE(18,'(A,'': '',2E20.8)') (cnames(l),                       &
     &    totmas(1,l)*lmolec/na, totmas(2,l)*lmolec/na, l=1,nc)
! Output average inventories.
        WRITE(18,*)
        WRITE(18,'(A/)') ' AVERAGE INVENTORIES OVER PERIOD IN MOLES'
        WRITE(18,'(A,'': '',2E20.8)') (cnames(l),                       &
     &    totavg(1,l)*lmolec/na, totavg(2,l)*lmolec/na, l=1,nc)
! Output fluxes.
        WRITE(18,*) 'NFLUX: ',nflux
        WRITE(18,'(/A,F6.3,A/A,F6.3)')                                  &
     &    ' FLUX ANALYSIS IN MOLES PER ',daymas-daym0,                  &
     &    ' DAY(S)',' DATA COLLECTED FROM DAY ',daymas
        DO l = 1, nflux
          IF (flist(1,l) < 1700) THEN
            WRITE(18,'(1X,I3,1X,A,'': '',1PE20.8)')                     &
     &        l,fnames(l),totflu(l)*lmolec/na
          ELSE
            WRITE(18,'(1X,I3,1X,A,'': '',1PE20.8)')                     &
     &        l,fnames(l),totflu(l)
          END IF
        END DO
! Output 3-D flux names.
        l = len_flx_str-1
        WRITE(18,*) ' '
        WRITE(18,*) '3-D fluxes requested:'
        DO j=1,numflux
          cstring=fluxnames(j)
          IF (cstring(l:len_flx_str) == ' 1') WRITE(18,'(A46)') cstring
        END DO
! Write out name of emission scenario used.
        WRITE(18,*) ' '
        SELECT CASE(scenario)
        CASE('fi')
          WRITE(18,*) 'SRES A1FI Emission Scenario Used'
        CASE('a2')
          WRITE(18,*) 'SRES A2 Emission Scenario Used'
        CASE('b1')
          WRITE(18,*) 'SRES B1 Emission Scenario Used'
        CASE('b2')
          WRITE(18,*) 'SRES B2 Emission Scenario Used'
        CASE('ab')
          WRITE(18,*) 'SRES A1B Emission Scenario Used'
        CASE('pi')
          WRITE(18,*) 'PRE-INDUSTRIAL Emission Scenario Used'
        CASE('A2')
          WRITE(18,*) 'SRES A2 2030 for IPCC4AR Emission Scenario Used'
        CASE('bu')
          WRITE(18,*) 'IIASA CLE Emission Scenario Used'
        CASE('mf')
          WRITE(18,*) 'IIASA MFR Emission Scenario Used'
        CASE DEFAULT
          WRITE(18,*) 'Unknown Emission Scenario Used'
        END SELECT
        CLOSE(18)
      END IF

!  Set and write out the fixed header
! DEPENDS ON: writehdr
      CALL WRITEHDR(month,year,daym0,daymas,filename,                   &
     &  flheader,fixhd12,umstepno,z_top_of_model,                       &
     &  first_constant_r_rho_level)

! Set PP lookup headers for all output fields
      nhdr=nlev*(num3dflux+1+nchem)
! DEPENDS ON: setlook
      CALL SETLOOK(month,year,daym0,daymas,flist,clist,stash_code,      &
     &  nchem,nhdr,flheader,pplook)
! Buffer lookup headers to file
      icode=0
      Word_Address=flheader(150)-1
! DEPENDS ON: setpos
      CALL SETPOS(20,Word_Address,icode)
      IF (icode /= 0) THEN
        cmessage='Error in SETPOS'
        WRITE(6,*)  cmessage
! DEPENDS ON: ereport
        CALL EREPORT('OUTBUD',1,cmessage)
      END IF
! DEPENDS ON: buffout
      CALL BUFFOUT(20,pplook,64*nhdr,len_io,a_io)
      IF (a_io /= -1.0.OR.len_io /= 64*nhdr) THEN
        cmessage='Error writing PPfile data:1 Buffer out of data'
        WRITE(6,*)  cmessage,a_io,len_io,64*nhdr
! DEPENDS ON: ereport
        CALL EREPORT('OUTBUD',1,cmessage)
      END IF

! Write out the cell distribution array
      len_in=nlong*mnlat
      Nfield=0
      DO i=1,nlev
! DEPENDS ON: sum3d
        CALL SUM3D(REAL(nm))
        Nfield = Nfield + 1
        Word_Address = PPLOOK(29,Nfield)
! DEPENDS ON: setpos
        CALL SETPOS(20,Word_Address,icode)
        IF (icode /= 0) THEN
          cmessage = 'Error in SETPOS'
          WRITE(6,*)  cmessage
! DEPENDS ON: ereport
          CALL EREPORT('OUTBUD',1,cmessage)
        END IF
! DEPENDS ON: buffout
        CALL BUFFOUT(20,global(:,:,i),len_in,len_io,a_io)
        IF (a_io /= -1.0.OR.len_io /= len_in) THEN
          cmessage='Error writing PPfile data:2 buffer out of data'
          WRITE(6,*)  cmessage,a_io,len_io,len_in
! DEPENDS ON: ereport
          CALL EREPORT('OUTBUD',1,cmessage)
        END IF
      END DO

! Write out the 3-D fluxes
      DO l=1,nflux
        IF (flist(2,l) > 0) THEN
! Range for lightning flashes, surf dep resistances, AOT40
          IF (flist(1,l) > 1700) THEN
! DEPENDS ON: sum3d
            CALL SUM3D(tflux(flist(2,l),:,:,:)*1.0/REAL(navg))
          ELSE
! DEPENDS ON: sum3d
            CALL SUM3D(tflux(flist(2,l),:,:,:)*lmolec/na)
          END IF
          DO i=1,nlev
            Nfield = Nfield + 1
            Word_Address = PPLOOK(29,Nfield)
! DEPENDS ON: setpos
            CALL SETPOS(20,word_address,icode)
            IF (icode /= 0) THEN
              cmessage='Error in SETPOS'
              WRITE(6,*)  cmessage
! DEPENDS ON: ereport
              CALL EREPORT('OUTBUD',1,cmessage)
            END IF
! DEPENDS ON: buffout
            CALL BUFFOUT(20,global(:,:,i),len_in,len_io,a_io)
            IF (a_io /= -1.0.OR.len_io /= len_in) THEN
              cmessage='Error writing PPfile data:3 buffer out of data'
              WRITE(6,*)  cmessage,a_io,len_io,len_in
! DEPENDS ON: ereport
              CALL EREPORT('OUTBUD',1,cmessage)
            END IF
          END DO
        END IF
      END DO

! For each output species in turn, create a global 3-D field,
! then write out the field.
      DO ll=1,nchem
        l=clist(ll)
! DEPENDS ON: sum3d
        CALL SUM3D(out(ll,:,:,:))
        DO i=1,nlev
          Nfield=Nfield+1
          Word_Address=PPLOOK(29,Nfield)
! DEPENDS ON: setpos
          CALL SETPOS(20,Word_Address,icode)
          IF (icode /= 0) THEN
            cmessage='Error in SETPOS'
            WRITE(6,*)  cmessage
! DEPENDS ON: ereport
            CALL EREPORT('OUTBUD',1,cmessage)
          END IF
! DEPENDS ON: buffout
          CALL BUFFOUT(20,global(:,:,i),len_in,len_io,a_io)
          IF (a_io /= -1.0.OR.len_io /= len_in) THEN
            cmessage='Error writing PPfile data:4 buffer out of data'
            WRITE(6,*)  cmessage,a_io,len_io,len_in
! DEPENDS ON: ereport
            CALL EREPORT('OUTBUD',1,cmessage)
          END IF
        END DO
      END DO

! DEPENDS ON: file_close
      CALL FILE_CLOSE(20,filename,14,1,0,icode)
      IF (icode /= 0) THEN
        cmessage='Error in FILE_CLOSE'
        WRITE(6,*)  cmessage
! DEPENDS ON: ereport
        CALL EREPORT('OUTBUD',1,cmessage)
      END IF
      WRITE(6,*) 'OUTBUD: Wrote out ',Nfield,' 3-D fields to file: ',   &
     &  FILENAME

      END SUBROUTINE OUTBUD
#endif
