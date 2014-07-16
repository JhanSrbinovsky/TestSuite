#if defined(C80_1A) || defined(UTILIO) || defined(RECON)               \
 || defined(VAROPSVER)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE PR_FIXHD---------------------------------------
!LL
!LL  Purpose: Prints out fixed length header record and checks
!LL           validity of information.
!LL
!LL AD, DR, SI  <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  3.1   22/12/92     Allow use by ancillary field headers
!LL                     Author A. Dickinson    Reviewer C. Wilson
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL  3.2   19/04/93     Code for new real missing data indicator.
!LL                     Author T.Johns         Reviewer A.Dickinson
!LL  3.4   5/07/94  Allowed for new PF model vert coord (density*r*r)
!LL                 and for Varobs file.    Author Colin Parrett.
!LL  4.0  21/11/95   Allowed for Covariance files.   C.Parrett
!LL  4.1  09/04/96  Introduce wave sub-model.  RTHBarnes.
!LL  4.1   16/04/96  Allowed for OPS Obstore files.  Author Colin Parret
!LL  4.1  23/05/96   Removed check of data length for MPP code  P.Burton
!LL  4.1  22/05/96   Print out UM Version Number   D. Robinson
!LL  5.1  10/04/00  New reconfiguration and control of stdout with
!LL                 PrintStatus. P.Selwood.
!LL  5.1  02/05/00   Increase format sizes. D Robinson.
!LL  5.5  17/02/03   Upgrade Wave model from 4.1 to 5.5 D.Holmes-Bell
!LL  6.2  23/11/05   Removed all references to the wavemodel.
!LL                  T.Edwards
!LL  6.2  19/01/06   Change labelling for Validity Times printed from
!LL                  boundary datasets. D Robinson
!LL
!LL  Programming standard:
!LL           Unified Model Documentation Paper No 3
!LL           Version No 2  09/09/91
!LL
!LL  System component: C25
!LL
!LL  System task: F3
!LL
!LL  Documentation:
!LL           Unified Model Documentation Paper No F3
!LL           Version No 5 9/2/90
!LL
!LL------------------------------------------------------------
!*L Arguments:-------------------------------------------------
      SUBROUTINE PR_FIXHD                                               &
     &(FIXHD,LEN_FIXHD,LEN_INTHD,LEN_REALHD,LEN1_LEVDEPC                &
     &,LEN2_LEVDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC,LEN1_COLDEPC,LEN2_COLDEPC &
     &,LEN1_FLDDEPC,LEN2_FLDDEPC,LEN_EXTCNST,LEN_DUMPHIST,LEN_CFI1      &
     &,LEN_CFI2,LEN_CFI3,LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA               &
     &,ICODE,CMESSAGE)

#if defined(RECON) || defined(VAROPSVER)
      USE Rcf_PrintStatus_Mod, Only : PrintStatus, PrStatus_Min,        &
     &                            PrStatus_Oper
#endif
      IMPLICIT NONE

      INTEGER                                                           &
     & LEN_FIXHD                                                        &
                     !IN Length of fixed length header
     &,LEN_INTHD                                                        &
                     !IN Length of integer header
     &,LEN_REALHD                                                       &
                     !IN Length of real header
     &,LEN1_LEVDEPC                                                     &
                     !IN 1st dim of level dep consts
     &,LEN2_LEVDEPC                                                     &
                     !IN 2nd dim of level dep consts
     &,LEN1_ROWDEPC                                                     &
                     !IN 1st dim of row dep consts
     &,LEN2_ROWDEPC                                                     &
                     !IN 2nd dim of row dep consts
     &,LEN1_COLDEPC                                                     &
                     !IN 1st dim of column dep consts
     &,LEN2_COLDEPC                                                     &
                     !IN 2nd dim of column dep consts
     &,LEN1_FLDDEPC                                                     &
                     !IN 1st dim of field dep consts
     &,LEN2_FLDDEPC                                                     &
                     !IN 2nd dim of field dep consts
     &,LEN_EXTCNST                                                      &
                     !IN Length of extra constants
     &,LEN_DUMPHIST                                                     &
                     !IN Length of history block
     &,LEN_CFI1                                                         &
                     !IN Length of comp field index 1
     &,LEN_CFI2                                                         &
                     !IN Length of comp field index 2
     &,LEN_CFI3                                                         &
                     !IN Length of comp field index 3
     &,LEN1_LOOKUP                                                      &
                     !IN 1st dim of lookup
     &,LEN2_LOOKUP   !IN 2nd dim of lookup

      INTEGER                                                           &
     & FIXHD(LEN_FIXHD)                                                 &
                        !IN Fixed length header
     &,LEN_DATA                                                         &
                        !IN Length of real data
     &,ICODE          !OUT Return code; successful=0
                      !                 error > 0

      CHARACTER*(80)                                                    &
     & CMESSAGE       !OUT Error message if ICODE > 0

! -------------------------------------------------------------
! Workspace usage:---------------------------------------------
! None
! -------------------------------------------------------------
!*L External subroutines called:-------------------------------
! None
!*-------------------------------------------------------------
#include "c_mdi.h"
#if !defined(RECON) && !defined(VAROPSVER)
#include "cprintst.h"
#endif
! Local variables:---------------------------------------------
      INTEGER I
!--------------------------------------------------------------

      ICODE=0
      CMESSAGE=' '

      If ( PrintStatus >= PrStatus_Oper ) Then
      WRITE(6,'('' '')')
      WRITE(6,'('' FIXED LENGTH HEADER'')')
      WRITE(6,'('' -------------------'')')

      WRITE(6,'('' Dump format version'',I6)')FIXHD(1)
      WRITE(6,'('' UM Version No      '',I6)')FIXHD(12)
      End If

      IF(FIXHD(2) == 1)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Atmospheric data'')')
        End If
      ELSE IF(FIXHD(2) == 2)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Oceanic data'')')
        End If
      ELSE IF (FIXHD(2) == 4) THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Wave sub-model data'')')
        End If
      ELSE
        If ( PrintStatus >= PrStatus_Min ) Then
        WRITE(6,'('' ***FATAL ERROR*** Invalid data type: FIXHD(2)='',  &
     &I9)')FIXHD(2)
        End If
       ICODE=4
       CMESSAGE='PR_FIXHD: Consistency check'
       RETURN
      ENDIF

      IF(FIXHD(3) == 1)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' On hybrid levels'')')
        End If
      ELSEIF(FIXHD(3) == 2)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' On sigma levels'')')
        End If
      ELSEIF(FIXHD(3) == 3)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' On pressure levels'')')
        End If
      ELSEIF(FIXHD(3) == 4)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' On depth levels'')')
        End If
      ELSEIF(FIXHD(3) == 5)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Charney-Phillips on radius levels'')')
        End If
      ELSEIF (FIXHD(3) == 6) THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
          WRITE(6,'('' Wave model direction levels and frequency pseudo &
     &-levels'')')
        End If
      ELSEIF(FIXHD(3) == IMDI.AND.FIXHD(5) == 4)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Missing data indicator used for vert coord type    &
     &'')')
        End If
      ELSE
        If ( PrintStatus >= PrStatus_Min ) Then
          WRITE(6,'('' ***FATAL ERROR*** Invalid level type: FIXHD(3)=  &
     &'', I9)')FIXHD(3)
        End If
       ICODE=4
       CMESSAGE='PR_FIXHD: Consistency check'
       RETURN
      ENDIF

      IF(FIXHD(4) == 0)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Over global domain'')')
        End If
      ELSEIF(FIXHD(4) == 1)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Over N. Hemispheric domain'')')
        End If
      ELSEIF(FIXHD(4) == 2)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Over S. Hemispheric domain'')')
        End If
      ELSEIF(FIXHD(4) == 3)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Over LA domain with no wrap around'')')
        End If
      ELSEIF(FIXHD(4) == 4)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Over LA domain with wrap around'')')
        End If
      ELSEIF(FIXHD(4) == 103)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Over rotated LA domain'')')
        End If
      ELSE
        If ( PrintStatus >= PrStatus_Min ) Then
        WRITE(6,'('' ***FATAL ERROR*** Invalid domain: FIXHD(4)='',     &
     &I9)')FIXHD(4)
        End If
       ICODE=4
       CMESSAGE='PR_FIXHD: Consistency check'
       RETURN
      ENDIF

      IF(FIXHD(5) == 1)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Instantaneous dump'')')
        End If
      ELSEIF(FIXHD(5) == 2)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Meaned dump'')')
        End If
      ELSEIF(FIXHD(5) == 3)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' FIELDS file'')')
        End If
      ELSEIF(FIXHD(5) == 4)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Ancillary dataset'')')
        End If
      ELSEIF(FIXHD(5) == 5)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Boundary dataset'')')
        End If
      ELSEIF(FIXHD(5) == 6)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' AC Observation File'')')
        End If
      ELSEIF(FIXHD(5) == 7)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Var Observation File'')')
        End If
      ELSEIF(FIXHD(5) == 8)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Cx file (model columns at ob locations)'')')
        End If
       ELSEIF(FIXHD(5) == 9)THEN
         If ( PrintStatus >= PrStatus_Oper ) Then
         WRITE(6,'('' Covariance File'')')
         End If
      ELSE IF (FIXHD(5) == 10) THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE (6, '('' OPS Obstore file '')')
        End If
      ELSE
        If ( PrintStatus >= PrStatus_Min ) Then
        WRITE(6,'('' ***FATAL ERROR*** Invalid dump type: FIXHD(5)='',  &
     &I9)')FIXHD(5)
        End If
       ICODE=4
       CMESSAGE='PR_FIXHD: Consistency check'
       RETURN
      ENDIF

      If ( PrintStatus >= PrStatus_Oper ) Then
      WRITE(6,'('' Exp No ='',I6,'' Run Id ='',I6)') FIXHD(7),FIXHD(6)
      End If

      IF(FIXHD(8) == 1)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Gregorian calendar'')')
        End If
      ELSEIF(FIXHD(8) == 2)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' 360-day calendar'')')
        End If
      ELSEIF(FIXHD(8) == IMDI.AND.FIXHD(5) == 4)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Missing data indcator used as calendar'')')
        End If
      ELSE
        If ( PrintStatus >= PrStatus_Min ) Then
          WRITE(6,'('' *** FATAL ERROR *** Invalid calendar type:       &
     &    FIXHD(8) ='',I9)')FIXHD(8)
        End If
       ICODE=4
       CMESSAGE='PR_FIXHD: Consistency check'
       RETURN
      ENDIF

      IF(FIXHD(9) == 1)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Arakawa A grid'')')
        End If
      ELSEIF(FIXHD(9) == 2)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Arakawa B grid'')')
        End If
      ELSEIF(FIXHD(9) == 3)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Arakawa C grid'')')
        End If
      ELSEIF(FIXHD(9) == 4)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Arakawa D grid'')')
        End If
      ELSEIF(FIXHD(9) == 5)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Arakawa E grid'')')
        End If
      ELSEIF(FIXHD(9) == IMDI.AND.FIXHD(5) == 4)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Missing data indicator used for grid type'')')
        End If
      ELSEIF(FIXHD(9) == IMDI.AND.FIXHD(5) == 5)THEN
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Missing data indicator used for grid type'')')
        End If
      ELSE
        If ( PrintStatus >= PrStatus_Min ) Then
          WRITE(6,'('' *** FATAL ERROR *** Invalid grid type:           &
     &    FIXHD(9)=''   ,I9)')FIXHD(9)
        End If
       ICODE=4
       CMESSAGE='PR_FIXHD: Consistency check'
       RETURN
      ENDIF

      If ( PrintStatus >= PrStatus_Oper ) Then

      If ( FIXHD(5) == 5 ) THEN    !  Boundary dataset

      WRITE(6,'(''                       Year  Month Day Hour Min  Sec  &
     & DayNo  '')')
      WRITE(6,'('' First Validity time ='',7I5)')(FIXHD(I),I=21,27)
      WRITE(6,'('' Last  Validity time ='',7I5)')(FIXHD(I),I=28,34)
      WRITE(6,'('' Interval            ='',7I5)')(FIXHD(I),I=35,41)

      Else

      WRITE(6,'(''                 Year  Month Day Hour Min  Sec DayNo  &
     &  '')')
      WRITE(6,'('' Data time     ='',7I5)')(FIXHD(I),I=21,27)
      WRITE(6,'('' Validity time ='',7I5)')(FIXHD(I),I=28,34)
      WRITE(6,'('' Creation time ='',7I5)')(FIXHD(I),I=35,41)

      End If

      WRITE(6,'(''                        Start     1st dim    2nd dim  &
     &  1st parm    2nd parm'')')
      WRITE(6,'('' Integer Consts   '',2I11,11X,I11)')FIXHD(100),       &
     &  FIXHD(101),LEN_INTHD
      WRITE(6,'('' Real Consts      '',2I11,11X,I11)')FIXHD(105),       &
     &  FIXHD(106),LEN_REALHD
      WRITE(6,'('' Level Dep Consts '',5I11)')FIXHD(110),               &
     &  FIXHD(111),FIXHD(112),LEN1_LEVDEPC,LEN2_LEVDEPC
      WRITE(6,'('' Row Dep Consts   '',5I11)')FIXHD(115),               &
     &  FIXHD(116),FIXHD(117),LEN1_ROWDEPC,LEN2_ROWDEPC
      WRITE(6,'('' Column Dep Consts'',5I11)')FIXHD(120),               &
     &  FIXHD(121),FIXHD(122),LEN1_COLDEPC,LEN2_COLDEPC
      WRITE(6,'('' Fields of Consts '',5I11)')FIXHD(125),               &
     &  FIXHD(126),FIXHD(127),LEN1_FLDDEPC,LEN2_FLDDEPC
      WRITE(6,'('' Extra Consts     '',2I11,11X,I11)')FIXHD(130),       &
     &  FIXHD(131),LEN_EXTCNST
      WRITE(6,'('' History Block    '',2I11,11X,I11)')FIXHD(135),       &
     &  FIXHD(136),LEN_DUMPHIST
      WRITE(6,'('' CFI No 1         '',2I11,11X,I11)')FIXHD(140),       &
     &  FIXHD(141),LEN_CFI1
      WRITE(6,'('' CFI No 2         '',2I11,11X,I11)')FIXHD(142),       &
     &  FIXHD(143),LEN_CFI2
      WRITE(6,'('' CFI No 3         '',2I11,11X,I11)')FIXHD(144),       &
     &  FIXHD(145),LEN_CFI3
      WRITE(6,'('' Lookup Tables    '',5I11)')FIXHD(150),               &
     &  FIXHD(151),FIXHD(152),LEN1_LOOKUP,LEN2_LOOKUP
      WRITE(6,'('' Model Data       '',2I11,11X,I11)')FIXHD(160),       &
     &  FIXHD(161),LEN_DATA
      End If

! Check model parameters against header record entries

      IF(FIXHD(101) >  0)THEN
        IF(LEN_INTHD /= FIXHD(101))THEN
        If ( PrintStatus >= PrStatus_Min ) Then
        WRITE(6,'('' *ERROR* Integer Consts'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        End If
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
      ENDIF

      IF(FIXHD(105) >  0)THEN
        IF(LEN_REALHD /= FIXHD(106))THEN
        If ( PrintStatus >= PrStatus_Min ) Then
        WRITE(6,'('' *ERROR* Real Consts'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        End If
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
      ENDIF

      IF(FIXHD(110) >  0)THEN
       IF(LEN1_LEVDEPC /= 0)THEN
        IF(LEN1_LEVDEPC /= FIXHD(111).OR.LEN2_LEVDEPC /= FIXHD(112))THEN
        If ( PrintStatus >= PrStatus_Min ) Then
        WRITE(6,'('' *ERROR* Level Dep Consts'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        End If
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
       ENDIF
      ENDIF

      IF(FIXHD(115) >  0)THEN
       IF(LEN1_ROWDEPC /= 0)THEN
        IF(LEN1_ROWDEPC /= FIXHD(116).OR.LEN2_ROWDEPC /= FIXHD(117))THEN
        If ( PrintStatus >= PrStatus_Min ) Then
        WRITE(6,'('' *ERROR* Row Dep Consts'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        End If
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
       ENDIF
      ENDIF

      IF(FIXHD(120) >  0)THEN
       IF(LEN1_COLDEPC /= 0)THEN
        IF(LEN1_COLDEPC /= FIXHD(121).OR.LEN2_COLDEPC /= FIXHD(122))THEN
        If ( PrintStatus >= PrStatus_Min ) Then
        WRITE(6,'('' *ERROR* Column Dep Consts'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        End If
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
       ENDIF
      ENDIF

      IF(FIXHD(125) >  0)THEN
       IF(LEN1_FLDDEPC /= 0)THEN
        IF(LEN1_FLDDEPC /= FIXHD(126).OR.LEN2_FLDDEPC /= FIXHD(127))THEN
        If ( PrintStatus >= PrStatus_Min ) Then
        WRITE(6,'('' *ERROR* Fields of Consts'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        End If
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
       ENDIF
      ENDIF

      IF(FIXHD(130) >  0)THEN
       IF(LEN_EXTCNST /= 0)THEN
        IF(LEN_EXTCNST /= FIXHD(131))THEN
        If ( PrintStatus >= PrStatus_Min ) Then
        WRITE(6,'('' *ERROR* Extra Consts'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        End If
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
       ENDIF
      ENDIF

      IF(FIXHD(135) >  0)THEN
       IF(LEN_DUMPHIST /= 0)THEN
        IF(LEN_DUMPHIST /= FIXHD(136))THEN
        If ( PrintStatus >= PrStatus_Min ) Then
        WRITE(6,'('' *ERROR* History File'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        End If
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
       ENDIF
      ENDIF

      IF(FIXHD(140) >  0)THEN
       IF(LEN_CFI1 /= 0)THEN
        IF(LEN_CFI1 /= FIXHD(141))THEN
        If ( PrintStatus >= PrStatus_Min ) Then
        WRITE(6,'('' *ERROR* CFI No 1'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        End If
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
       ENDIF
      ENDIF

      IF(FIXHD(142) >  0)THEN
       IF(LEN_CFI2 /= 0)THEN
        IF(LEN_CFI2 /= FIXHD(143))THEN
        If ( PrintStatus >= PrStatus_Min ) Then
        WRITE(6,'('' *ERROR* CFI No 2'')')
        End If
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Parameter and header values dont match'')')
        End If
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
       ENDIF
      ENDIF

      IF(FIXHD(144) >  0)THEN
       IF(LEN_CFI3 /= 0)THEN
        IF(LEN_CFI3 /= FIXHD(145))THEN
        If ( PrintStatus >= PrStatus_Min ) Then
        WRITE(6,'('' *ERROR* CFI No 3'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        End If
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
       ENDIF
      ENDIF

      IF(FIXHD(150) >  0)THEN
        IF(LEN2_LOOKUP /= IMDI)THEN
        IF(LEN1_LOOKUP /= FIXHD(151).OR.LEN2_LOOKUP /= FIXHD(152))THEN
        If ( PrintStatus >= PrStatus_Min ) Then
        WRITE(6,'('' *ERROR* Lookup Table'')')
        End If
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Parameter and header values dont match'')')
        End If
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
        ENDIF
      ENDIF

#if !defined(MPP)
      IF(FIXHD(160) >  0)THEN
        IF(LEN_DATA /= IMDI)THEN
        IF(LEN_DATA /= FIXHD(161))THEN
        If ( PrintStatus >= PrStatus_Min ) Then
        WRITE(6,'('' *ERROR* Model Data'')')
        End If
        If ( PrintStatus >= PrStatus_Oper ) Then
        WRITE(6,'('' Parameter and header values dont match'')')
        End If
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
        ENDIF
      ENDIF
#endif

      RETURN
      END SUBROUTINE PR_FIXHD
#endif
