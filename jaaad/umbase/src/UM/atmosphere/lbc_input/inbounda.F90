#if defined(A31_1A)
#if defined(ATMOS) && !defined(GLOBAL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine INBOUNDA
!LL
!LL Purpose : Takes as input,the code defining whether updates of
!LL  boundary data are required.  The physical files required are
!LL  identified, and the headers lookup tables are read into common
!LL  blocks.  Reads the update intervals from the boundary datasets.
!LL  Where the update interval is in months or years, the check will be
!LL  made daily.  This is a new routine at 5.3 which separates out the
!LL  atmosphere from INBOUND1
!LL
!LL Control routine for CRAY YMP
!LL
!LL CW, SI, ZG      <- programmer of some or all of previous code or
!LL                    changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  5.3     27/04/01 Separate routine for atmosphere.  Was previously
!LL                   part of INBOUND1.  Z Gardner
!LL  5.3     04/08/01 Allow sub-hourly input of atmosphere LBCs.
!LL                   Pete Clark.
!LL  5.3     31/10/01 Remove RimWeightsA_Orog. D.Robinson
!LL  5.4     17/09/02 Allow use of two boundary files, plus tidying.
!LL                   Adam Clayton/Dave Robinson
!    5.4     29/08/02 Include cruntimc for L_Fixed_lbcs logical
!                     and if true set BOUND_FIELDCODE(1) = 0 so
!                     that atmos lbc files are not read. R.M.Forbes
!   5.5     17/02/03 Move position of c_submodl.h  -  M.Hughes
!    5.5     03/02/03 Include additional microphysics lbcs in
!                     ITEM_BOUNDA array if active.   R.M.Forbes
!    5.5     30/01/03 CURRENT_TIME_DAYS and STEPS_TO_BDI_START
!                     calculations fixed. Adam Clayton
!    6.0     29/07/03 Include cloud fractions in lbcs. Damian Wilson
!    6.1     28/05/04 Allow start from any step. N. Christidis
!    6.1     02/08/04 Remove 5.4 mod which switched off reading of lbc
!                     file for fixed lbcs             Terry  Davies
!    6.2     24/03/05 Revise calculation of STEP_TO_BDI_START.
!                     Initialise BNDARY_OFFSETim and Current_LBC_Step.
!                     Dave Robinson.
!    6.2     01/07/05 Reinstate 5.4 mod which switched off reading of
!                     lbc file for fixed lbcs. R. Forbes.
!    6.2     29/01/06 Extra L_force_lbc logical.  Yongming Tang
!    6.2     01/10/04 Include murk aerosol lbcs in
!                     ITEM_BOUNDA array if active.   R.M.Forbes
!    6.2     21/03/06 Included nstypes.h J Ridley.
!LL
!LL Programming standard; Unified Model Documentation Paper No. 3
!LL version no. 1, dated 15/01/90
!LL
!LL Logical components covered : C720
!LL
!LL System task : C7
!LL
!LL Documentation : Unified Model Documentation Paper No C7
!LLEND
!
!*L  Arguments
! *****************************COPYRIGHT*******************************
      SUBROUTINE INBOUNDA(                                              &
#include "argduma.h"
#include "argsts.h"
#include "argptra.h"
#include "argbnd.h"
#include "argppx.h"
     &   A_LEN1_LEVDEPCDA,A_LEN2_LEVDEPCDA,                             &
     &   A_LEN1_ROWDEPCDA,A_LEN2_ROWDEPCDA,                             &
     &   A_LEN1_COLDEPCDA,A_LEN2_COLDEPCDA)

      IMPLICIT NONE

#include "cmaxsize.h"
#include "parvars.h"
#include "typsize.h"
#include "typduma.h"
#include "csubmodl.h"
#include "typsts.h"
#include "typptra.h"
#include "typbnd.h"
#include "nstypes.h"

      INTEGER :: A_LEN1_LEVDEPCDA   ! IN : copy of A_LEN1_LEVDEPC
      INTEGER :: A_LEN2_LEVDEPCDA   ! IN : copy of A_LEN2_LEVDEPC
      INTEGER :: A_LEN1_ROWDEPCDA   ! IN : copy of A_LEN1_ROWDEPC
      INTEGER :: A_LEN2_ROWDEPCDA   ! IN : copy of A_LEN2_ROWDEPC
      INTEGER :: A_LEN1_COLDEPCDA   ! IN : copy of A_LEN1_COLDEPC
      INTEGER :: A_LEN2_COLDEPCDA   ! IN : copy of A_LEN2_COLDEPC
      
#include "ppxlook.h"
#include "clookadd.h"
#include "chsunits.h"
#include "chistory.h"
#include "ccontrol.h"
#include "ctime.h"
#include "c_mdi.h"
#include "cenvir.h"
#include "ctracera.h"
#include "cruntimc.h"
#include "cprintst.h"

!    Local variables

      INTEGER                                                           &
     &        I,                                                        &
     &        J,                                                        &
     &        J1,                                                       &
     &        lbc_num,                                                  &
     &        START_BLOCK,                                              &
     &        NFTIN,                                                    &
     &        NFTENV,                                                   &
                                   ! Index in FT_ENV arrays
     &        im_index,                                                 &
                                   ! Internal model index
     &        ELAPSED_DAYS,                                             &
                                   ! Days since basis time
     &        ELAPSED_SECS,                                             &
                                   ! Secs since basis time
     &        CURRENT_TIME_DAYS,                                        &
                                   ! No. of days to current time
     &        CURRENT_TIME_SECS,                                        &
                                   ! No. of secs-in-day to current time
     &        DAYS_TO_DATA_START,                                       &
                                   ! Days  to start of boundary data
     &        SECS_TO_DATA_START,                                       &
                                   ! Secs  to start of boundary data
     &        DAYS_TO_DATA_END,                                         &
                                   ! Days  to end   of boundary data
     &        SECS_TO_DATA_END,                                         &
                                   ! Secs  to end   of boundary data
     &        STEPS_TO_DATA_START,                                      &
                                   ! Steps to start of boundary data
     &        STEPS_TO_DATA_END,                                        &
                                   ! Steps to end   of boundary data
     &        RIM_STEPSA_OLD,                                           &
                                   ! Data interval for last bndy file
     &        STEPS_TO_BDI_START,                                       &
                                   ! Steps to start/end of
     &        STEPS_TO_BDI_END,                                         &
                                   ! current boundary data interval
     &        BASIS_TO_DATA_START_STEPS,                                &
                                         ! Steps from basis time to
                                         ! start of boundary data
     &        ITEM_BOUNDA(RIM_LOOKUPSA)  ! Boundary updatable item list

      LOGICAL :: THIS_ALBC_FOR_BDI_END ! True if same boundary file to
                                       ! be used for end of current
                                       ! boundary data interval

      REAL                                                              &
     &      A_LEVDEPC_BO(A_LEN1_LEVDEPCDA,A_LEN2_LEVDEPCDA),            &
     &      A_ROWDEPC_BO(A_LEN1_ROWDEPCDA,A_LEN2_ROWDEPCDA),            &
     &      A_COLDEPC_BO(A_LEN1_COLDEPCDA,A_LEN2_COLDEPCDA)
      
      INTEGER FULL_LOOKUP_BOUNDA(LEN1_LOOKUP,BOUND_LOOKUPSA)

      INTEGER, Parameter :: DUMMY =1

      INTEGER             :: ErrorStatus      ! Return code
      CHARACTER (Len=256) :: CMESSAGE         ! Error message
      CHARACTER (Len=*), Parameter :: RoutineName = 'inbounda'

      NAMELIST/BOUNCNST/BOUND_FIELDCODE,RIMWEIGHTSA

      REAL P1,P2
      LOGICAL LNER
      LNER(P1,P2) = ((ABS(P1-P2))  >   (1.E-6*ABS(P1+P2)))

      LOGICAL, SAVE :: L_FirstCall = .TRUE.

!L      Internal Structure

      ErrorStatus=0
      CMESSAGE=' '

      IF (L_FirstCall) THEN

!  Rewind file and read namelist for boundary updating constants

        REWIND 5
        READ(5,BOUNCNST)

! If fixed lbcs is set then lateral boundaries will be fixed at
! the initial conditions. Turn off code to read in atmos lbc file
      If (L_Fixed_lbcs) BOUND_FIELDCODE(1) = 0

! If idealised forced lbcs option is on then lateral boundaries will
! be forced with a timeseries of single profiles from the idealised
! namelist. Turn off code to read in atmos lbc file
      If (L_Force_lbc) BOUND_FIELDCODE(1) = 0

!L  1.0 Initialise variables in COMMON/CBND/, lest undefined for some
!L      choices of boundary updating, eg. -DEF,GLOBAL but not DEF,FLOOR,
!L      as used in section 2 to set BOUNDARY_STEPSim for atmos
        RIM_STEPSA = 0
!  Initialise BNDARY_OFFSETim & BOUNDARY_STEPSim in CTIME comdeck
        DO I=1,INTERNAL_ID_MAX
          BNDARY_OFFSETim(I) = 0
          BOUNDARY_STEPSim(I) = 0
        END DO

      END IF ! (L_FirstCall)

!L  1.1 Update interval for lateral boundaries for atmosphere
!L      Read headers and test whether boundary updating required

      IF ( BOUND_FIELDCODE(1) <= 0) THEN
        RIM_STEPSA=0
      ELSE

!L      Open input boundary file and read headers

        NBOUND_LOOKUP(1)=1

        ! Use same unit number for first and second boundary files
        NFTIN = 125

        IF (ALBC_num == 2) THEN
          NFTENV = 126 ! Env Var for 2nd boundary file
        ELSE
          NFTENV = 125 ! Env Var for 1st boundary file
        END IF

! DEPENDS ON: file_open
        CALL FILE_OPEN(NFTIN,FT_ENVIRON(NFTENV),                        &
     &                 LEN_FT_ENVIR(NFTENV),0,0,ErrorStatus)
        IF(ErrorStatus /= 0) THEN
          Write (Cmessage,*)                                            &
     &    'Failure opening boundary file : No ', ALBC_num
! DEPENDS ON: ereport
          Call Ereport(RoutineName,ErrorStatus,Cmessage)
        ENDIF

!       Read in fixed header to get array dimensions
! DEPENDS ON: read_flh
        CALL READ_FLH(NFTIN,FIXHD_BOUNDA(1,1),                          &
     &                       LEN_FIXHD,ErrorStatus,CMESSAGE)
        IF (ErrorStatus >  0) THEN
          WRITE (Cmessage,*)                                            &
     &          'INBOUNDA : Error in READ_FLH for BOUNDA(1,1)'
! DEPENDS ON: ereport
           Call Ereport(RoutineName,ErrorStatus,Cmessage)
        ENDIF

!       Check for negative dimensions
        IF (FIXHD_BOUNDA(101,1) <= 0) FIXHD_BOUNDA(101,1)=1
        IF (FIXHD_BOUNDA(106,1) <= 0) FIXHD_BOUNDA(106,1)=1
        IF (FIXHD_BOUNDA(111,1) <= 0) FIXHD_BOUNDA(111,1)=1
        IF (FIXHD_BOUNDA(112,1) <= 0) FIXHD_BOUNDA(112,1)=1
        IF (FIXHD_BOUNDA(116,1) <= 0) FIXHD_BOUNDA(116,1)=1
        IF (FIXHD_BOUNDA(117,1) <= 0) FIXHD_BOUNDA(117,1)=1
        IF (FIXHD_BOUNDA(121,1) <= 0) FIXHD_BOUNDA(121,1)=1
        IF (FIXHD_BOUNDA(122,1) <= 0) FIXHD_BOUNDA(122,1)=1
        IF (FIXHD_BOUNDA(151,1) <= 0) FIXHD_BOUNDA(151,1)=1
        IF (FIXHD_BOUNDA(152,1) <= 0) FIXHD_BOUNDA(152,1)=1
        IF (FIXHD_BOUNDA(161,1) <= 0) FIXHD_BOUNDA(161,1)=1

!       Check if sufficient space allocated for LOOKUP table
        IF (FIXHD_BOUNDA(152,1) >  BOUND_LOOKUPSA) THEN
      write(6,*)' INBOUNDA; not enough space for LBC lookup headers.'
      write(6,*)'           try increasing value specified in umui'
      write(6,*)'           window atmos_Infile_Options_Headers'
          write(CMESSAGE,*)                                             &
     &        'INBOUNDA: Insufficient space for Lookup Table'
          Errorstatus = 2
! DEPENDS ON: ereport
          Call Ereport(RoutineName,ErrorStatus,Cmessage)
        ENDIF

! DEPENDS ON: setpos
        CALL SETPOS (NFTIN,0,ErrorStatus)
        IF (ErrorStatus >  0) THEN
          WRITE (6,*) 'INBOUNDA: Problem with SETPOS for BOUNDA(1,1)'
          WRITE (6,*) 'ErrorStatus ',ErrorStatus,' NFTIN ',NFTIN
          Write(Cmessage,*) 'Problem with Setpos for Bounda(1,1)'
! DEPENDS ON: ereport
          Call Ereport(RoutineName,Errorstatus,Cmessage)
        ENDIF

! DEPENDS ON: readhead
        CALL READHEAD(NFTIN,                                            &
     &                FIXHD_BOUNDA(1,1),LEN_FIXHD,                      &
     &                INTHD_BOUNDA(1,1),FIXHD_BOUNDA(101,1),            &
     &                REALHD_BOUNDA(1,1),FIXHD_BOUNDA(106,1),           &
     &                A_LEVDEPC_BO(1,1),                                &
     &                FIXHD_BOUNDA(111,1),FIXHD_BOUNDA(112,1),          &
     &                A_ROWDEPC_BO(1,1),                                &
     &                FIXHD_BOUNDA(116,1),FIXHD_BOUNDA(117,1),          &
     &                A_COLDEPC_BO(1,1),                                &
     &                FIXHD_BOUNDA(121,1),FIXHD_BOUNDA(122,1),          &
     &                DUMMY,DUMMY,DUMMY,                                &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                FULL_LOOKUP_BOUNDA,                               &
     &                FIXHD_BOUNDA(151,1),FIXHD_BOUNDA(152,1),          &
     &                FIXHD_BOUNDA(161,1),                              &
#include "argppx.h"
     &                START_BLOCK,ErrorStatus,CMESSAGE)

        IF (ErrorStatus >  0) THEN
          WRITE (6,*) 'INBOUNDA: Problem with READHEAD for BOUNDA(1,1)'
          WRITE (6,*) 'ErrorStatus ',ErrorStatus,' CMESSAGE ',CMESSAGE
          Write(cmessage,*) 'Problem with READHEAD for BOUNDA(1,1)'
! DEPENDS ON: ereport
          Call Ereport(RoutineName,errorStatus,cmessage)
        ENDIF

! Copy the first set of headers into LOOKUP_BOUNDA
        DO i=1,RIM_LOOKUPSA
          DO j=1,LEN1_LOOKUP
            LOOKUP_BOUNDA(j,i)=FULL_LOOKUP_BOUNDA(j,i)
          ENDDO ! j
        ENDDO ! i

! Copy all the varying items from the LOOKUP into COMP_LOOKUP_BOUNDA
        DO i=1,BOUND_LOOKUPSA

          LOOKUP_COMP_BOUNDA(LBCC_LBYR,i)=FULL_LOOKUP_BOUNDA(LBYR,i)
          LOOKUP_COMP_BOUNDA(LBCC_LBMON,i)=FULL_LOOKUP_BOUNDA(LBMON,i)
          LOOKUP_COMP_BOUNDA(LBCC_LBDAT,i)=FULL_LOOKUP_BOUNDA(LBDAT,i)
          LOOKUP_COMP_BOUNDA(LBCC_LBHR,i)=FULL_LOOKUP_BOUNDA(LBHR,i)
          LOOKUP_COMP_BOUNDA(LBCC_LBMIN,i)=FULL_LOOKUP_BOUNDA(LBMIN,i)
          LOOKUP_COMP_BOUNDA(LBCC_LBDAY,i)=FULL_LOOKUP_BOUNDA(LBDAY,i)
          LOOKUP_COMP_BOUNDA(LBCC_LBEGIN,i)=FULL_LOOKUP_BOUNDA(LBEGIN,i)
          LOOKUP_COMP_BOUNDA(LBCC_NADDR,i)=FULL_LOOKUP_BOUNDA(NADDR,i)

        ENDDO ! i

! Check validity of headers

! Integer headers
        IF (FIXHD_BOUNDA(100,1)  >   0) THEN

          IF (INTHD_BOUNDA(6,1) /= glsize(1,fld_type_p)) THEN
            WRITE(6,*) 'LBC Integer Header Mismatch:'
            WRITE(6,*) 'ROW_LENGTH from INTHD: ',INTHD_BOUNDA(6,1)
            WRITE(6,*) 'Model ROW_LENGTH: ',glsize(1,fld_type_p)

            ErrorStatus=3
            Write(CMESSAGE,*)'Integer header (ROW_LENGTH) mismatch'
! DEPENDS ON: ereport
            Call Ereport(RoutineName,ErrorStatus,Cmessage)
          ENDIF

          IF (INTHD_BOUNDA(7,1) /= glsize(2,fld_type_p)) THEN
            WRITE(6,*) 'LBC Integer Header Mismatch:'
            WRITE(6,*) 'Number of rows from INTHD: ',INTHD_BOUNDA(7,1)
            WRITE(6,*) 'Model number of rows: ',glsize(2,fld_type_p)

            ErrorStatus=4
            Write(CMESSAGE,*)'Integer header (N_ROWS) mismatch'
! DEPENDS ON: ereport
            Call Ereport(Routinename,Errorstatus,cmessage)
          ENDIF

          If ( INTHD_BOUNDA(17,1) /= A_INTHD(17) ) Then
            WRITE(6,*) 'LBC Integer Header Mismatch:'
            WRITE(6,*) 'LBC   : Height Generator Method : ',            &
     &                  INTHD_BOUNDA(17,1)
            WRITE(6,*) 'Model : Height Generator Method : ',            &
     &                  A_INTHD(17)

            ErrorStatus=5
            Write (CMESSAGE,*) 'INBOUNDA : Mis-match in height ',       &
     &                         'generator method.'
! DEPENDS ON: ereport
            Call Ereport(RoutineName,ErrorStatus,Cmessage)
          End If

          If ( INTHD_BOUNDA(24,1) /= A_INTHD(24) ) Then
            WRITE(6,*) 'LBC Integer Header Mismatch:'
            WRITE(6,*) 'LBC  : First rho level with constant height: ', &
     &                  INTHD_BOUNDA(24,1)
            WRITE(6,*) 'Model: First rho level with constant height: ', &
     &                  A_INTHD(24)

            ErrorStatus=6
            Write (CMESSAGE,*) 'INBOUNDA : Mis-match in height ',       &
     &                         'generator method.'
! DEPENDS ON: ereport
            Call Ereport(RoutineName,ErrorStatus,Cmessage)
          End if

        ENDIF ! IF (FIXHD_BOUNDA(100,1)  >   0)

! Real constants

        IF (FIXHD_BOUNDA(105,1)  >   0) THEN

          ! Check real headers in LBC file against those in dump to
          ! 32 bit accuracy
          DO j=1,6
            IF (LNER(REALHD_BOUNDA(j,1),A_REALHD(j))) THEN
              WRITE(6,*) 'LBC Real Header mismatch at position ',j,':'
              WRITE(6,*) 'Value from LBC file is ',REALHD_BOUNDA(j,1)
              WRITE(6,*) 'Value from model dump is ',A_REALHD(j)

              ErrorStatus=7
              Write(CMESSAGE,*) 'INBOUNDA : Real header mismatch'
! DEPENDS ON: ereport
              Call Ereport(RoutineName,ErrorStatus,Cmessage)
            ENDIF
          ENDDO ! j

          IF (LNER(REALHD_BOUNDA(16,1),A_REALHD(16))) THEN
            WRITE(6,*) 'LBC Real Header mismatch at position 16 :'
            WRITE(6,*) 'LBC file : Height at top of model ',            &
     &      REALHD_BOUNDA(16,1)
            WRITE(6,*) 'Model    : Height at top of model ',            &
     &      A_REALHD(16)

            ErrorStatus=8
            Write(CMESSAGE,*) 'INBOUNDA : Model Height mismatch'
! DEPENDS ON: ereport
            Call Ereport(RoutineName,ErrorStatus,Cmessage)
          ENDIF

        ENDIF ! IF (FIXHD_BOUNDA(105,1)  >   0)

! Level dependent constants

        If (FIXHD_BOUNDA(110,1) > 0) Then

! ---------------------------------
! Check eta values for theta levels
! ---------------------------------

          Do j=1, model_levels+1
            If ( LNER( A_LEVDEPC_BO(j,1),                               &
     &                 A_LEVDEPC(jetatheta+j-1) ) ) Then

              WRITE(6,*) 'LBC Level Dependent Constants Mismatch'
              WRITE(6,*) 'Eta values for theta levels mismatch for ',   &
     &                   'level ',J
              WRITE(6,*) 'Value from LBC file   : ',A_LEVDEPC_BO(j,1)
              WRITE(6,*) 'Value from model dump : ',                    &
     &                    A_LEVDEPC(jetatheta+j-1)
              ErrorStatus=9
              Write(CMESSAGE,*)                                         &
     &              'INBOUNDA : Level dependent constant mismatch'
! DEPENDS ON: ereport
              Call Ereport(RoutineName,ErrorStatus,Cmessage)

            End If
          End Do ! j

! -------------------------------
! Check eta values for rho levels
! -------------------------------

          Do j=1, model_levels
            If ( LNER( A_LEVDEPC_BO(j,2),                               &
     &                 A_LEVDEPC(jetarho+j-1) ) ) Then

              WRITE(6,*) 'LBC Level Dependent Constants Mismatch'
              WRITE(6,*) 'Eta values for rho levels mismatch for ',     &
     &                   'level ',J
              WRITE(6,*) 'Value from LBC file   : ',A_LEVDEPC_BO(j,2)
              WRITE(6,*) 'Value from model dump : ',                    &
     &                    A_LEVDEPC(jetarho+j-1)
              ErrorStatus=10
              Write(CMESSAGE,*)                                         &
     &            'INBOUNDA : Level dependent constant mismatch'
! DEPENDS ON: ereport
              Call Ereport(RoutineName,ErrorStatus,Cmessage)

            End If
          End Do ! j

        End If ! If (FIXHD_BOUNDA(110,1) > 0)

!L      Set update interval
!       If update interval includes months or years, a 360 day
!       calender assumed.

        ! Save previous value of RIM_STEPSA:
        IF (.NOT.L_FirstCall) RIM_STEPSA_OLD = RIM_STEPSA

        RIM_STEPSA=((FIXHD_BOUNDA(35,1)*8640+FIXHD_BOUNDA(36,1)*720     &
     &   +FIXHD_BOUNDA(37,1)*24+FIXHD_BOUNDA(38,1))*3600                &
     &   +FIXHD_BOUNDA(39,1)*60+FIXHD_BOUNDA(40,1))                     &
     &   *STEPS_PER_PERIODim(a_im)/SECS_PER_PERIODim(a_im)

        ! Check that RIM_STEPSA has not changed:
        IF (.NOT.L_FirstCall) THEN
          IF (RIM_STEPSA /= RIM_STEPSA_OLD) THEN
            ErrorStatus = 1
            WRITE (CMessage,*)                                          &
     &      'Boundary updating period (RIM_STEPSA) has changed from ',  &
     &      RIM_STEPSA_OLD, ' to ', RIM_STEPSA, '. Not allowed!'
! DEPENDS ON: ereport
            CALL EReport (RoutineName, ErrorStatus, CMessage)
          END IF
        END IF

        ! Initialise Current_LBC_Step

        If (L_FirstCall) Then
          Current_LBC_Step = 1 + STEPim(atmos_im)
          If (PrintStatus >= PrStatus_Normal) Then
            write (6,*) ' INBOUNDA : Timestep ',STEPim(a_im),           &
     &                  ' Current_LBC_Step ',Current_LBC_Step
          End If
        End If

        IF (Num_ALBCs == 2 .AND. L_FirstCall) THEN

          ! Calculate step on which to swap boundary files:
          ALBC_SwapStep = ALBC2_StartTime_steps - RIM_STEPSA

          ! Check that ALBC_SwapStep is non-negative:
          IF (ALBC_SwapStep < 0) THEN
            ErrorStatus = 1
            WRITE (CMessage,*)                                          &
     &        'Step on which to swap boundary files (ALBC_SwapStep = ', &
     &        ALBC_SwapStep, ') must be non-negative'
! DEPENDS ON: ereport
            CALL EReport (RoutineName, ErrorStatus, CMessage)
          END IF

          ! Check that ALBC_SwapStep is a multiple of RIM_STEPSA:
          IF (MOD(ALBC_SwapStep, RIM_STEPSA) /= 0) THEN
            ErrorStatus = 1
            WRITE (CMessage,*)                                          &
     &        'Step on which to swap boundary files (ALBC_SwapStep = ', &
     &        ALBC_SwapStep, ') must be a a multiple of the boundary '//&
     &        'updating period (RIM_STEPSA = ', RIM_STEPSA, ')'
            ErrorStatus=101
! DEPENDS ON: ereport
            CALL EReport (RoutineName, ErrorStatus, CMessage)
          END IF

        END IF

      END IF

!L  2   Set interval for setting any boundary field

      BOUNDARY_STEPSim(a_im) = RIM_STEPSA

!L  3    Check LOOKUP Table

      J1=0
      IF (BOUND_FIELDCODE(1) /= 0) THEN
        J1=FIXHD_BOUNDA(152,1)

        IF(FIXHD_BOUNDA(150,1) >  0) THEN

! Set up list of variables expected to be boundary updated.
          ITEM_BOUNDA(1)  = 31001 ! Orography
          ITEM_BOUNDA(2)  = 31002 ! U
          ITEM_BOUNDA(3)  = 31003 ! V
          ITEM_BOUNDA(4)  = 31004 ! W
          ITEM_BOUNDA(5)  = 31005 ! Density
          ITEM_BOUNDA(6)  = 31006 ! Potential temperature
          ITEM_BOUNDA(7)  = 31007 ! Specific humidity
          ITEM_BOUNDA(8)  = 31008 ! QCL
          ITEM_BOUNDA(9)  = 31009 ! QCF
          ITEM_BOUNDA(10) = 31010 ! Exner
          ITEM_BOUNDA(11) = 31011 ! U_Adv
          ITEM_BOUNDA(12) = 31012 ! V_Adv
          ITEM_BOUNDA(13) = 31013 ! W Adv

! Include additional microphysics lbcs if active

          lbc_num = 13
          IF (L_mcr_qcf2_lbc) THEN  ! qcf2 lbcs active
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31014
          ENDIF
          IF (L_mcr_qrain_lbc) THEN  ! qrain lbcs active
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31015
          ENDIF
          IF (L_mcr_qgraup_lbc) THEN  ! qgraup lbcs active
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31016
          ENDIF

          ! Setup for additional cloud fraction lbcs if active
          If (L_pc2_lbc) Then
           lbc_num = lbc_num+1
           ITEM_BOUNDA(lbc_num) = 31017 ! cf_bulk
           lbc_num = lbc_num+1
           ITEM_BOUNDA(lbc_num) = 31018 ! cf_liquid
           lbc_num = lbc_num+1
           ITEM_BOUNDA(lbc_num) = 31019 ! cf_frozen
          EndIf

          ! Include murk aerosol lbcs if in input lbc file
          If (L_murk_lbc) Then
            lbc_num = lbc_num + 1
            ITEM_BOUNDA(lbc_num) = 31020 ! murk aerosol
          EndIf

          IF (TR_VARS  >   0) THEN

            ! Find STASH item no. for each tracer in use
            i=0
            im_index=internal_model_index(A_IM)
            DO j = A_TRACER_FIRST,A_TRACER_LAST
              IF (SI(j,0,im_index) /= 1) THEN  ! tracer is in use
                i=i+1
                ITEM_BOUNDA(lbc_num+i)=31100+(j-A_TRACER_FIRST)
              ENDIF
            ENDDO ! j

          ENDIF ! IF (TR_VARS  >   0)

! DEPENDS ON: chk_look_bounda
          CALL CHK_LOOK_BOUNDA(                                         &
     &      ITEM_BOUNDA,FULL_LOOKUP_BOUNDA,                             &
#include "argbnd.h"
#include "argppx.h"
     &                         ErrorStatus,CMESSAGE)

          IF (ErrorStatus  /=  0) Then
            Write(CMessage,*) 'Error in CHK_LOOK_BOUNDA'
! DEPENDS ON: ereport
            Call Ereport(RoutineName,ErrorStatus,Cmessage)
          End If

!L Find start position in lookup tables

          ! Get days/seconds since basis time:
! DEPENDS ON: time2sec
          CALL TIME2SEC (I_YEAR,             I_MONTH,                   &
     &                   I_DAY,              I_HOUR,                    &
     &                   I_MINUTE,           I_SECOND,                  &
     &                   BASIS_TIME_DAYS,    BASIS_TIME_SECS,           &
     &                   ELAPSED_DAYS,       ELAPSED_SECS,              &
     &                   LCAL360)

          ! Get current model time in days/seconds-in-day:
          CURRENT_TIME_DAYS = BASIS_TIME_DAYS + ELAPSED_DAYS +          &
     &                       (BASIS_TIME_SECS + ELAPSED_SECS)/86400

          CURRENT_TIME_SECS = MOD(BASIS_TIME_SECS+ELAPSED_SECS, 86400)

          ! Get days/seconds to start of boundary data:
! DEPENDS ON: time2sec
          CALL TIME2SEC (FIXHD_BOUNDA(21,1), FIXHD_BOUNDA(22,1),        &
     &                   FIXHD_BOUNDA(23,1), FIXHD_BOUNDA(24,1),        &
     &                   FIXHD_BOUNDA(25,1), FIXHD_BOUNDA(26,1),        &
     &                   CURRENT_TIME_DAYS,  CURRENT_TIME_SECS,         &
     &                   DAYS_TO_DATA_START, SECS_TO_DATA_START,        &
     &                   LCAL360)

          ! Get days/seconds to end of boundary data:
! DEPENDS ON: time2sec
          CALL TIME2SEC (FIXHD_BOUNDA(28,1), FIXHD_BOUNDA(29,1),        &
     &                   FIXHD_BOUNDA(30,1), FIXHD_BOUNDA(31,1),        &
     &                   FIXHD_BOUNDA(32,1), FIXHD_BOUNDA(33,1),        &
     &                   CURRENT_TIME_DAYS,  CURRENT_TIME_SECS,         &
     &                   DAYS_TO_DATA_END,   SECS_TO_DATA_END,          &
     &                   LCAL360)

          ! Get steps to start of boundary data:
! DEPENDS ON: tim2step
          CALL TIM2STEP (DAYS_TO_DATA_START,                            &
     &                   SECS_TO_DATA_START,                            &
     &                   STEPS_PER_PERIODim(a_im),                      &
     &                   SECS_PER_PERIODim(a_im),                       &
     &                   STEPS_TO_DATA_START)

          ! Get steps to end of boundary data:
! DEPENDS ON: tim2step
          CALL TIM2STEP (DAYS_TO_DATA_END,                              &
     &                   SECS_TO_DATA_END,                              &
     &                   STEPS_PER_PERIODim(a_im),                      &
     &                   SECS_PER_PERIODim(a_im),                       &
     &                   STEPS_TO_DATA_END)

          ! Get steps from basis time to start of boundary data:
          BASIS_TO_DATA_START_STEPS = STEPim(a_im) + STEPS_TO_DATA_START

          ! Check that the above is a multiple of RIM_STEPSA only
          ! when two boundary files exist
          IF ((Num_ALBCs == 2) .AND.                                    &
     &        (MOD(BASIS_TO_DATA_START_STEPS, RIM_STEPSA) /= 0)) THEN
            ErrorStatus = 1
            WRITE (CMessage,*)                                          &
     &        'Steps from basis time to start of boundary data (',      &
     &        BASIS_TO_DATA_START_STEPS, ') must be a multiple of '//   &
     &        'the boundary updating period (', RIM_STEPSA, ')'
! DEPENDS ON: ereport
            CALL EReport (RoutineName, ErrorStatus, CMessage)
          END IF

          ! There are two situations we need to cater for here:
          !
          !   1. There is no boundary data in memory.
          !   2. The boundary data valid at the start of the current
          !      boundary data interval is already in memory, but
          !      data for the end of the interval needs to be read in.
          !
          ! The first case applies if and only if we are on the first
          ! call to this routine.
          !
          ! In the first case, we need to make sure that there is
          ! boundary data valid at the start of the current boundary
          ! data interval, and then point to it. We also need to make
          ! sure that data is available for the end of the interval. The
          ! only exception to the latter is the case where the data for
          ! the end of the interval is to come from a different boundary
          ! file, which may be the case if two boundary files are being
          ! used.
          !
          ! In the second case, we need to make sure that there is
          ! boundary data valid at the end of the current boundary data
          ! interval, and then point to it.

          ! Steps to start of current boundary data interval:
          STEPS_TO_BDI_START = MOD(BASIS_TO_DATA_START_STEPS,           &
     &                             RIM_STEPSA)

          ! Steps to end   of current boundary data interval:
          STEPS_TO_BDI_END   = STEPS_TO_BDI_START + RIM_STEPSA

          ! Initialise BNDARY_OFFSET for use in CTIME comdeck
          BNDARY_OFFSETim(a_im) = MOD(-BASIS_TO_DATA_START_STEPS,       &
     &                                 RIM_STEPSA)

          IF (L_FirstCall) THEN ! Case 1

            IF (STEPS_TO_DATA_START <= STEPS_TO_BDI_START .AND.         &
     &          STEPS_TO_DATA_END   >= STEPS_TO_BDI_START) THEN
              NBOUND_LOOKUP(1) = (-STEPS_TO_DATA_START / RIM_STEPSA)    &
     &                         * (RIM_LOOKUPSA-1) + 2
            ELSE IF (STEPS_TO_DATA_START > STEPS_TO_BDI_START) THEN
              ErrorStatus = 101
              WRITE (CMessage,*)                                        &
     &          'Boundary data starts after start of current '//        &
     &          'boundary data interval'
! DEPENDS ON: ereport
              CALL EReport (RoutineName, ErrorStatus, CMessage)
            ELSE
              ErrorStatus = 102
              WRITE (CMessage,*)                                        &
     &          'Boundary data ends before start of current '//         &
     &          'boundary data interval'
! DEPENDS ON: ereport
              CALL EReport (RoutineName, ErrorStatus, CMessage)
            END IF

            ! Is the data for the end of the current boundary data
            ! interval to come from the same boundary file?
            THIS_ALBC_FOR_BDI_END =.TRUE.
            IF (NUM_ALBCs == 2) THEN
              IF (STEPim(a_im) <= ALBC_SwapStep) THEN
                THIS_ALBC_FOR_BDI_END =.FALSE.
              END IF
            END IF

            ! If so, check that the data exists:
            IF (THIS_ALBC_FOR_BDI_END) THEN
              IF (STEPS_TO_DATA_START > STEPS_TO_BDI_END) THEN
                ErrorStatus = 101
                WRITE (CMessage,*)                                      &
     &            'Boundary data starts after end of current '//        &
     &            'boundary data interval'
! DEPENDS ON: ereport
                CALL EReport (RoutineName, ErrorStatus, CMessage)
              ELSE IF (STEPS_TO_DATA_END < STEPS_TO_BDI_END) THEN
                ErrorStatus = 102
                WRITE (CMessage,*)                                      &
     &            'Boundary data ends before end of current '//         &
     &            'boundary data interval'
! DEPENDS ON: ereport
                CALL EReport (RoutineName, ErrorStatus, CMessage)
              END IF
            END IF

          ELSE ! Not first call (Case 2)

            IF (STEPS_TO_DATA_START <= STEPS_TO_BDI_END .AND.           &
     &          STEPS_TO_DATA_END   >= STEPS_TO_BDI_END) THEN
              NBOUND_LOOKUP(1) = (1 - STEPS_TO_DATA_START / RIM_STEPSA) &
     &                         * (RIM_LOOKUPSA-1) + 2
            ELSE IF (STEPS_TO_DATA_START > STEPS_TO_BDI_END) THEN
              ErrorStatus = 101
              WRITE (CMessage,*)                                        &
     &          'Boundary data starts after end of current '//          &
     &          'boundary data interval'
! DEPENDS ON: ereport
              CALL EReport (RoutineName, ErrorStatus, CMessage)
            ELSE
              ErrorStatus = 102
              WRITE (CMessage,*)                                        &
     &          'Boundary data ends before end of current '//           &
     &          'boundary data interval'
! DEPENDS ON: ereport
              CALL EReport (RoutineName, ErrorStatus, CMessage)
            END IF

          END IF ! (L_FirstCall)

        END IF

      END IF ! lateral boundary

      L_FirstCall = .FALSE.

!L  4   End of routine

      RETURN
      END SUBROUTINE INBOUNDA
#endif
#endif
