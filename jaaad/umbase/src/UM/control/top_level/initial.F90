#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: INITIAL --------------------------------------------------
!LL
!LL  Purpose: Initialises the model ready for integration/assimilation.
!LL This involves reading the model control files and setting up STASH,
!LL reading the initial or restart dump,
!LL initialising the ancillary, boundary and interface
!LL field control routines and updating the ancillary fields on restart
!LL if time to do so, exchanging coupling fields and swapping dumps (if
!LL a coupled model), and initialising the assimilation package if
!LL necessary.  Subsidiary control routines are called to perform these
!LL functions.
!LL
!LL  Tested under compiler:  sxmpif90
!LL  Tested under OS version: 
!LL
!LL  Author:   T.C.Johns
!LL
!LL  Programming standard: UM Doc Paper 3, version 8 (01/06/2007)
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
      SUBROUTINE INITIAL(                                               &
#include "argppx.h"
#include "argszsp.h"
#include "argszspa.h"
#include "argszspc.h"
#include "argsp.h"
#include "argspa.h"
#include "argspc.h"
     &     IAU_lookup,                                                  &
     &     D1_IAU_k4,                                                   &
     &     D1_TFilt,                                                    &
     &     ngrgas,grgas_addr,                                           &
     &     internal_model,submodel,NGROUP,MEANLEV,                      &
     &     isSTASH)  ! FLUME-STASH UM vn7.1
!

! FLUME module
use atm_fields_mod

#if defined(ACCESS)
USE auscom_cpl_data_mod, ONLY : l_auscom
#endif

! Declarations:

      IMPLICIT NONE

! Common blocks:

#include "cmaxsize.h"
!CL Super array lengths
#include "typszsp.h"
#include "typszspa.h"
#include "typszspc.h"
!CL Super arrays
#include "typspd1.h"
#include "typspdua.h"
#include "typspst.h"
#include "typsppta.h"
#include "typspcoa.h"
#include "typspina.h"
#include "typspana.h"
#include "typspbo.h"
#include "typspboa.h"
#include "typspcpl.h"
!L       Model sizes
#include "parvars.h"
#include "typsize.h"
#include "nstypes.h"
!L       Addresses of arrays within super arrays
#include "spindex.h"
#include "typocdpt.h"
#include "typwvdpt.h"
#include "csubmodl.h"
! Rick's mods have *CALL CSUBMODL here
#include "cppxref.h"
#include "ppxlook.h"

#include "chsunits.h"
#include "chistory.h"
#include "ccontrol.h"
#include "cruntimc.h"
#include "ctime.h"
#include "cpprint.h"
#include "cenvir.h"
#include "c_global.h"
#include "c_writd.h"
#include "decomptp.h"
#include "decompdb.h"
#include "ctfilt.h"
#include "lbc_coup.h"
#include "cbound.h"
#include "cprintst.h"
#include "c_kinds.h"

! Subroutine arguments:
      LOGICAL       isSTASH  ! FLUME-STASH UM vn7.1
      INTEGER,      INTENT(INOUT) ::                                    &
                                     ! Lookup tables of IAU inc file
      
     &  IAU_lookup(IAU_Len1Lookup, IAU_Len2Lookup)

      REAL(KIND=real32), INTENT(INOUT) ::                               &
                                          ! Array for packed IAU fields
      
     &  D1_IAU_k4(D1_IAU_k4_len)

      REAL,         INTENT(INOUT) ::                                    &
                                     ! Array for unpacked IAU/TDF fields
      
     &  D1_TFilt(D1_TFilt_len)

      INTEGER internal_model ! OUT internal model identifier:
!                            !   1:Atmos; 2:Ocean; 3:Slab ; etc
      INTEGER submodel       ! OUT submodel partition (dump) identifier:
!                            !   1:Atmos; 2:Ocean; etc
      INTEGER NGROUP         ! OUT   - No of steps in "group"n
      INTEGER MEANLEV        ! OUT - Mean level indicator

! 3-D fields of species to be passed down to radiation
      INTEGER, INTENT(IN) :: ngrgas
      INTEGER, INTENT(OUT) :: grgas_addr(ngrgas)

! Local variables:

      INTEGER  IMEAN      ! Loop index over mean periods
      INTEGER  I          ! Loop index
      INTEGER  ISM        ! Loop index over submodels
      INTEGER  submodel_next      ! Submodel identifier
      INTEGER  NFTASWAP,NFTOSWAP  ! Fortran units for coupling swapfiles
      INTEGER  NFTSWAP    ! General Fortran unit for coupling swapfiles
      INTEGER  FTN_UNIT   ! Fortran unit for pp files
      INTEGER  TRANSALEN,TRANSOLEN !data length for TRANSOUT/IN
      INTEGER  TRANS_LEN  ! General data length for TRANSOUT/IN
      INTEGER NDEV      ! Unit no.
      INTEGER LL        ! Counter
      INTEGER LEN_FILENAME ! Length of FILENAME array
      INTEGER  G_THETA_FIELD_SIZE                                       &
                                    ! Sizes for dynamic allocation
     &        ,G_IMTJMT          ! in A-O coupling routines
      INTEGER CO2_DIMA,                                                 &
                                   ! CO2 array dimensions
     &        CO2_DIMO,                                                 &
     &        CO2_DIMO2
      INTEGER DMS_DIMA,                                                 &
                                   ! DMS array dimensions
     &        DMS_DIMO,                                                 &
     &        DMS_DIMO2
      INTEGER Dummyarg  ! Not used, needed to end arg list
      INTEGER SECS_PER_STEPA ! Atmos timestep length in seconds
      CHARACTER*14 PPNAME ! Dummy PP file name returned by PPCTL
      CHARACTER*80 FILENAME
      CHARACTER*2 ENVALUE
!
      integer len_wait_tot     ! Total wait time for boundary data
      character*8 ch_date2     ! Date from date_and_time
      character*10 ch_time2    ! Time from date_and_time
      integer*8 sleep          ! SLEEP Function to make UM wait
      integer*8 ISLEEP         ! SLEEP Function to make UM wait
#if defined(ATMOS) && !defined(GLOBAL)
      integer info             ! Return Code from GCOM routine.
#endif
      integer lbc_ntimes       ! No of BC's in communication file
      integer ms_ntimes        ! No of BC's required in mesoscale

      LOGICAL                                                           &
     &  not_enough_lbcs                                                 &
                          ! True if more LBCs are required
     &, read_more_values  ! True if we are to read more values in

! Error reporting
      INTEGER       ICODE       ! =0 normal exit; >0 error exit
      CHARACTER*256 Cmessage    ! Error message
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='INITIAL')
! Mixing ratio code
      Logical, Parameter ::                                             &
     & l_mr_tfiltctl = .false.  ! Use mixing ratio code (if available)

!- End of header ------------------------------------------------------

! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INITIAL ',3)
      ICODE=0
      Cmessage=''
!L
!L----------------------------------------------------------------------
!L
!L 1.2 Set FT units as inactive on first step of the integration
!L     and set last field written/read to zero for each unit
!L
      IF (H_STEPim(a_im) == 0) THEN
        DO I=20,NUNITS
          FT_ACTIVE(I)='N'
          FT_LASTFIELD(I)=0
        ENDDO
      ENDIF
!L
!L 1.3 Option to write RADINCS array to file cache2 on unit 16 removed
!L
      IF (PrintStatus  >=  PrStatus_Normal) THEN
      WRITE(6,*) 'Fast i/o of radincs directly to core memory'
      ENDIF
!L
!L Number of CPUs attached to program (NCPU) is hard-wired to 1.
!L
      NCPU = 1
!L
!L---------------------------------------------------------------------
!L 2. Initialise STASH control arrays from STASH control file.
!L---------------------------------------------------------------------

        IF(LTIMER) THEN
! DEPENDS ON: timer
          CALL TIMER('INITCTL ',3)
        END IF

! Note that NSECTS=NSECTP, NITEMS=NITEMP : set in WSTLST

! DEPENDS ON: initctl
      CALL INITCTL(                                                     &
#include "artsts.h"
#include "argppx.h"
#include "artd1.h"
     &                   ICODE,CMESSAGE )

        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to INITCTL'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF
        IF(LTIMER) THEN
! DEPENDS ON: timer
          CALL TIMER('INITCTL ',4)
        END IF

!L
!L----------------------------------------------------------------------
!L 3. Read appropriate submodel partition dump to memory.  If coupled,
!L    page out the D1 part of each dump to its 'swap' file and read the
!L    other dump(s) into memory.  Write temporary history file if dumps
!L    read successfully on timestep 0.
!L
!L
!L 3.1  Loop over submodel partition dumps
!L
      DO ISM=1,N_SUBMODEL_PARTITION

        submodel=SUBMODEL_PARTITION_LIST(ISM)
        IF (LTIMER) THEN
! DEPENDS ON: timer
          CALL TIMER('INITDUMP',5)
! DEPENDS ON: timer
          CALL TIMER('INITDUMP',3)
        ENDIF
! DEPENDS ON: initdump
        CALL INITDUMP(                                                  &
#include "artd1.h"
#include "artduma.h"
#include "artsts.h"
#include "artptra.h"
#include "artlndm.h"
#include "argppx.h"
     &               submodel,ICODE,CMESSAGE)
        IF (LTIMER) THEN
! DEPENDS ON: timer
          CALL TIMER('INITDUMP',4)
! DEPENDS ON: timer
          CALL TIMER('INITDUMP',6)
        ENDIF
        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to INITDUMP'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF

      ENDDO  ! ISM over submodel partitions

#if defined(ATMOS)
! SET_ATM_FIELDS points the fields in atm_fields_mod to the correct parts of D1
! It should be called after INITDUMP (which calls SET_ATM_POINTERS)
! and before INITDIAG (which calls ST_DIAG<n>).

! DEPENDS ON: set_atm_fields
      CALL Set_Atm_Fields(                                              &
#include "artptra.h"
#include "artsts.h"
     &     SPD1(IXD1(2)), SPD1(IXD1(2)), SPD1(IXD1(2)))
#endif

!L
#if defined(ATMOS)
!L
!L      Set RUN indicator in atmosphere dump header
!L
! DEPENDS ON: set_run_indic_op
      CALL SET_RUN_INDIC_OP(                                            &
#include "artduma.h"
     &              ICODE,CMESSAGE)
      IF (ICODE  /=  0) THEN
        WRITE(6,*) 'Failure in call to SET_RUN_INDIC_OP'
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ICODE,Cmessage)
      ENDIF
#endif
!L
!L 3.3  On NRUN initialise dump LOOKUP headers associated with
!L      diagnostic fields with the bare essentials needed to read and
!L      write dumps - the rest to be filled in by STASH during the run.
!L
      IF (H_STEPim(a_im)  == 0) THEN
! DEPENDS ON: inithdrs
        CALL INITHDRS(                                                  &
#include "artduma.h"
#include "artsts.h"
#include "argppx.h"
     &                ICODE,CMESSAGE)
        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to INITHDRS'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF

      ENDIF    ! End test for NRUN
!
#if defined(ATMOS)
      if(L_perturb_IC_theta) then
! DEPENDS ON: perturb_theta_ctl
       call perturb_theta_ctl(                                          &
#include "artd1.h"
#include "artptra.h"
     &                         row_length, rows, model_levels,          &
     &                         global_row_length, global_rows,          &
     &                         model_domain, at_extremity,              &
     &                         offx, offy, IntRand_Seed,   datastart    &
     &                        )
     endif
#endif

!L
!L 3.4 Write out temporary copy of D1 for current submodel
!L

      icode=0
      IF (L_WRIT_INIT) THEN

        DO ISM=1,N_SUBMODEL_PARTITION

          submodel=SUBMODEL_PARTITION_LIST(ISM)
        if (submodel  ==  atmos_sm) then
! DEPENDS ON: dumpctl
           CALL DUMPCTL (                                               &
#include "artd1.h"
#include "artduma.h"
#include "artlndm.h"
#include "artptra.h"
#include "artsts.h"
#include "argppx.h"
     &          atmos_sm,0,.TRUE.,'atminitial',0,                       &
     &          ICODE,CMESSAGE)
          IF (ICODE  /=  0) THEN
            WRITE(6,*) 'Failure in call to DUMPCTL (Atmos)'
! DEPENDS ON: ereport
            CALL Ereport(RoutineName,ICODE,Cmessage)
          ENDIF

        endif
        ENDDO ! ISM

      END IF

!L----------------------------------------------------------------------
!L 6.  Initialise means program control block
!L
      DO ISM=1,N_SUBMODEL_PARTITION

        submodel=SUBMODEL_PARTITION_LIST(ISM)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INITMEAN',3)
! DEPENDS ON: initmean
        CALL INITMEAN(                                                  &
#include "artduma.h"
     &                submodel,ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INITMEAN',4)
        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to INITMEAN for submodel ',       &
     &               ISM
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF

      ENDDO ! ISM over submodel partition dumps
!L----------------------------------------------------------------------
!L 4. Set up other control blocks and physical constants
!L
!L 4.1  Initialise the model time and check that history file data time
!L      matches dump(s); set derived time/step information
!L
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INITTIME',3)
! DEPENDS ON: inittime
      CALL INITTIME(                                                    &
#include "artduma.h"
     &              submodel,ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INITTIME',4)
        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to INITTIME'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF
!L
!L 4.2  Write up temporary history file after successfully reading
!L      initial dumps on timestep 0 and setting model_data_time if
!L      assimilation run, to allow CRUN from initial dumps.
!L
#if defined(MPP)
      IF (mype  ==  0) THEN
#endif
      IF (H_STEPim(a_im) == 0) THEN
! DEPENDS ON: temphist
        CALL TEMPHIST(THIST_UNIT,ICODE,CMESSAGE)
        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to TEMPHIST'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF
      END IF
#if defined(MPP)
      ENDIF
#endif
!L
!L 4.3  Set up control block for updating ancillary fields
!L
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INANCCTL',3)
! DEPENDS ON: inancctl
      CALL INANCCTL(                                                    &
#include "artduma.h"
#include "artsts.h"
#include "artptra.h"
#include "artanc.h"
#include "argppx.h"
     &           ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INANCCTL',4)
        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to INANCCTL'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF
!L
!L 4.4  Set up control block for updating boundary fields

#if defined(ATMOS) && !defined(GLOBAL)

      ALBC_num = 1

      ! Setup for use of two atmos boundary files:
      IF (Num_ALBCs == 2) THEN

        ! Check namelist variable ALBC2_StartTime_mins:
        SECS_PER_STEPA = NINT(SECS_PER_STEPim(a_im))

        IF (MOD(ALBC2_StartTime_mins*60, SECS_PER_STEPA) /= 0) THEN
          ICODE = 1
          WRITE (CMessage,*)                                            &
     &      'ALBC2_StartTime_mins (', ALBC2_StartTime_mins,             &
     &      ') does not coincide with the start of a timestep'
! DEPENDS ON: ereport
          CALL EReport (RoutineName, ICODE, CMessage)
        ELSE
          ! Convert into steps:
          ALBC2_StartTime_steps = ALBC2_StartTime_mins*60/SECS_PER_STEPA
        END IF

        ! If on a continuation run, we may be able to go straight to
        ! the second boundary file:
        IF (STEPim(a_im) >= ALBC2_StartTime_steps) ALBC_num = 2

      END IF

      ! NOTE: If using two boundary files, coupling is only activated
      !       for the second.
      IF (l_lbc_coup .AND.                                              &
     &    .NOT.(Num_ALBCs == 2 .AND. ALBC_num == 1)) THEN ! Global/Mes
                                                          ! coupling

        ms_ntimes=2  ! For now we assume the minimum requirement -
                     ! two sets of boundary conditions are required.
                     ! Once INBOUND has been called we may find that
                     ! we need more than this, in which case we will
                     ! jump back to line 100

      ENDIF

 100  CONTINUE  ! Return here if IN_BOUND has been called and there
                ! are insufficient boundary conditions to proceed.
                ! This condition may arise in CRUNS

      IF (l_lbc_coup .AND.                                              &
     &    .NOT.(Num_ALBCs == 2 .AND. ALBC_num == 1)) THEN ! Global/Mes
                                                          ! coupling

        CALL DATE_AND_TIME(ch_date2, ch_time2)

        IF (PrintStatus  >=  PrStatus_Normal) THEN
          WRITE(6,*) 'LBC_COUP: ',                                      &
     &      ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),' on ',   &
     &      ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4),          &
     &      ' Wait to call INBOUND in INITIAL.'
        ENDIF ! (PrintStatus  >=  PrStatus_Normal)

! DEPENDS ON: timer
        Call Timer ('LBC_WAIT',3)

        IF (mype  ==  0) THEN  ! Only do the check on the communication
                               ! file on PE 0
          WRITE(6,*) 'ms_ntimes in INITIAL ',ms_ntimes

          not_enough_lbcs=.TRUE.
          len_wait_tot=0

          DO WHILE (not_enough_lbcs) ! loop until the communication
                                     ! file indicates we have enough
                                     ! lbcs to start running with
            read_more_values=.TRUE.

            CLOSE (190) ! We need to close and reopen the communication
                        ! file to ensure we see its latest state
            OPEN (190,FILE=lbc_filename,ACTION="read",IOSTAT=ICODE)

            IF (ICODE  /=  0) THEN ! Something went wrong with the OPEN
              WRITE(6,*) 'Return code from OPEN of communication file ',&
     &                   'in INITIAL: ',ICODE
              ICODE=401
              WRITE(CMESSAGE,*) 'INITIAL : OPEN failed for Unit 190'

              read_more_values=.FALSE. ! Jump out of the loops
              not_enough_lbcs=.FALSE.  ! and call Ereport
            ENDIF ! IF (ICODE  /=  0)

            DO WHILE (read_more_values) ! loop while we need to read
                                        ! more values from the
                                        ! communication file

              READ (190,*,IOSTAT=ICODE) lbc_ntimes ! read the next value

              IF (ICODE  /=  0) THEN ! Problem with the READ...

                WRITE(6,*) 'ms : Return code from READ ',ICODE

                IF (len_wait_tot  >=  um_lbc_wait_max) THEN

                  ! Maximum wait time has been exceeded.
                  ! Insufficient Boundary Conditions to proceed.
                  ! Likely cause is delay in job generating the BC's.

                  WRITE(6,*) 'ms: Maximum wait time reached ',          &
     &                       'after ',um_lbc_wait_max,' seconds.'
                  ICODE=402
                  WRITE(CMESSAGE,*) 'INITIAL : Failed to find ',        &
     &              'required value in LBC_FILE.'

                  read_more_values=.FALSE. ! Jump out of the loops
                  not_enough_lbcs=.FALSE.  ! and call Ereport

                ELSE ! We've not exceeded the time limit

                  ! Insufficient BCs to proceed.
                  ! Wait for um_lbc_wait seconds before another
                  ! attempt to read the file to see if more BCs
                  ! have been written.

                  WRITE(6,*) 'ms: Wait for ',um_lbc_wait,               &
     &                       ' seconds and retry.'

                  isleep=SLEEP(um_lbc_wait)
                  len_wait_tot=len_wait_tot+um_lbc_wait

                  WRITE(6,*) 'ms : Total wait so far ',len_wait_tot,    &
     &                       ' seconds.'

                  read_more_values=.FALSE.
                  ! No more values in the file as it stands, so this
                  ! forces us to close and reopen the file and read
                  ! again to see if any new values have arrived.

                ENDIF ! IF (len_wait_tot  >=  um_lbc_wait_max)

              ELSE ! the READ(190) was successful - now we must
                   ! interpret the value pulled from the file

                IF (lbc_ntimes  >   1000) THEN
                  ! First value in the file is always >1000

                  read_more_values=.TRUE. ! Go round and read next value

                ELSEIF (lbc_ntimes  <   ms_ntimes) THEN
                  ! Value is not required. Proceed to read next value

                  WRITE(6,*) 'ms : lbc_ntimes = ',lbc_ntimes,           &
     &                       ' read in.'//                              &
     &                       ' lbc_ntimes >= ',ms_ntimes,               &
     &                       ' is required. Read next value.'
                  read_more_values=.TRUE.

                ELSEIF (lbc_ntimes  >=  ms_ntimes) THEN
                  ! Required value read in. Sufficient BCs to proceed

                  WRITE(6,*) 'ms : lbc_ntimes = ',lbc_ntimes,           &
     &                       ' read in.'//                              &
     &                       ' lbc_ntimes >= ',ms_ntimes,               &
     &                       ' is required. Proceed.'

                  CALL DATE_AND_TIME (ch_date2, ch_time2)

                  IF (PrintStatus  >=  PrStatus_Normal) THEN
                    WRITE(6,*) 'LBC_COUP: ',                            &
     &              ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),  &
     &              ' on ',                                             &
     &              ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4),  &
     &              ' Proceed to call INBOUND in INITIAL.'
                  ENDIF ! (PrintStatus  >=  PrStatus_Normal)

                  read_more_values=.FALSE. ! Don't read any more values
                  not_enough_lbcs=.FALSE.  ! Don't re-examine the file

                ENDIF ! Test for value of lbc_ntimes

              ENDIF ! IF (ICODE  /=  0)

            ENDDO ! DO WHILE (read_more_values)

          ENDDO ! DO WHILE (not_enough_lbcs)

        ENDIF ! IF (mype  ==  0)

! DEPENDS ON: timer
        Call Timer ('LBC_WAIT',4)

        CALL GC_IBCAST(100,1,0,nproc,info,ICODE) ! PE 0 broadcasts ICODE
                                                 ! to all processors
        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'INITIAL - Error detected in LBC coupling.'
          WRITE(6,*) 'ICODE: ',ICODE
          WRITE(6,*) 'CMESSAGE : ',CMESSAGE

! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,CMESSAGE)
        ENDIF ! IF (ICODE  /=  0)

        CALL GC_IBCAST(100,1,0,nproc,info,lbc_ntimes) ! PE 0 broadcasts
                                                      ! lbc_ntimes to
                                                      ! all processors
      ENDIF ! (l_lbc_coup .AND.
            !  .NOT.(Num_ALBCs == 2 .AND. ALBC_num == 1))
#endif

!L
      ! FLUME_STASH UM vn7.1
      ! isSTASH defaults to false in a non-flume run
      IF (.NOT.isSTASH) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('IN_BOUND',3)
! DEPENDS ON: in_bound
      CALL IN_BOUND(                                                    &
#include "artduma.h"
#include "artsts.h"
#include "artptra.h"
#include "artbnd.h"
     &   A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,                                 &
     &   A_LEN1_ROWDEPC,A_LEN2_ROWDEPC,                                 &
     &   A_LEN1_COLDEPC,A_LEN2_COLDEPC,                                 &
                                          ! for dynamic array
#include "argppx.h"
     &                   ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('IN_BOUND',4)
        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to IN_BOUND'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF

      END IF  ! .NOT.isSTASH

#if defined(ATMOS) && !defined(GLOBAL)

      IF (l_lbc_coup .AND.                                              &
     &    .NOT.(Num_ALBCs == 2 .AND. ALBC_num == 1)) THEN ! Global/Mes
                                                          ! coupling

!       Now that IN_BOUND has been called for the first time
!       double check that there are sufficient BCs to proceed.
!
!       Determine which boundary data is required to proceed
!       the next period.
        if (boundary_stepsim(a_im) >  0) then
          IF (ALBC_num == 2) THEN
            ms_ntimes = 1 + ( (stepim(a_im)-ALBC_SwapStep)              &
     &                      / boundary_stepsim(a_im) )
          ELSE
            ms_ntimes = 2 + (stepim(a_im)/boundary_stepsim(a_im))
          END IF
        endif

        if (lbc_ntimes <  ms_ntimes) then

!         There are insufficient BCs to proceed. Go back and wait
!         for sufficient BCs to proceed.
          WRITE(6,*) 'ms : lbc_ntimes = ',lbc_ntimes,                   &
     &               ' lbc_ntimes >= ',ms_ntimes,' is required. '//     &
     &               'Insufficient BCs to proceed. Wait.'

          GOTO 100

        endif

      endif  ! (l_lbc_coup .AND.
             ! .NOT.(Num_ALBCs == 2 .AND. ALBC_num == 1))

!
!    4.4.1  Update atmosphere boundary fields at step zero
!

      IF (BOUNDARY_STEPSim(a_im)  /=  0) THEN ! If there are BCs

        IF (l_lbc_coup .AND.                                            &
     &      .NOT.(Num_ALBCs == 2 .AND. ALBC_num == 1)) THEN ! Global/Mes
                                                            ! coupling

          IF (PrintStatus  >=  PrStatus_Normal) THEN
            WRITE(6,*) 'LBC_COUP: ',                                    &
     &        ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),        &
     &        ' on ',                                                   &
     &        ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4),        &
     &        ' Proceed to call UPBOUND in INITIAL.'
          ENDIF ! (PrintStatus  >=  PrStatus_Normal)

        ENDIF

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('UP_BOUND',3)

! DEPENDS ON: up_bound
        CALL UP_BOUND(submodel,                                         &
#include "artd1.h"
#include "artduma.h"
#include "artptra.h"
#include "artbnd.h"
#include "argppx.h"
     &              ICODE,CMESSAGE)

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('UP_BOUND',4)

        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'INITIAL: Failure in call to atmosphere UP_BOUND'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF

        IF (Num_ALBCs == 2) THEN

          ! If the data for the start and end of the current boundary
          ! data interval is to come from different boundary files, we
          ! need to make additional calls to INBOUNDA/UP_BOUND:
          IF (ALBC_num == 1 .AND. STEPim(a_im) >= ALBC_SwapStep) THEN

            IF (PrintStatus >= PrStatus_Normal) THEN
              WRITE (6,*) ''
              WRITE (6,*) 'INITIAL: Swapping to 2nd atmos boundary file'
              WRITE (6,*) ''
            END IF

            ALBC_num = 2

            IF (l_lbc_coup) THEN

              ms_ntimes = 1
              GOTO 100

            ELSE

! DEPENDS ON: inbounda
              CALL INBOUNDA(                                            &
#include "artduma.h"
#include "artsts.h"
#include "artptra.h"
#include "artbnd.h"
#include "argppx.h"
     &          A_LEN1_LEVDEPC,A_LEN2_LEVDEPC)

! DEPENDS ON: timer
              IF (LTIMER) CALL TIMER('UP_BOUND',3)

! DEPENDS ON: up_bound
              CALL UP_BOUND(submodel,                                   &
#include "artd1.h"
#include "artduma.h"
#include "artptra.h"
#include "artbnd.h"
#include "argppx.h"
     &          ICODE,CMESSAGE)

! DEPENDS ON: timer
              IF (LTIMER) CALL TIMER('UP_BOUND',4)

              IF (ICODE /= 0) THEN
                WRITE(6,*)                                              &
     &            'INITIAL: Failure in call to atmosphere UP_BOUND'
! DEPENDS ON: ereport
                CALL EReport (RoutineName, ICODE, CMESSAGE)
              END IF

            END IF ! (l_lbc_coup)

          END IF ! (ALBC_num == 1 .AND. STEPim(a_im) >= ALBC_SwapStep)

        END IF ! (Num_ALBCs == 2)

      ENDIF ! IF (BOUNDARY_STEPSim(a_im)  /=  0)

#endif

#if defined(ATMOS)
!
!    4.4.2  Call SETCONA
!

      IF (submodel /= atmos_sm) THEN  ! If we've not got the atmosphere
                                      ! dump in memory, then read it in
        TRANS_LEN=TRANSALEN
        NFTSWAP  =NFTASWAP
        submodel=atmos_sm    ! new submodel will be atmosphere

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('TRANSIN',3)
! DEPENDS ON: transin
          CALL TRANSIN(                                                 &
#include "artd1.h"
     &            TRANS_LEN,                                            &
     &            NFTSWAP,submodel,ICODE,CMESSAGE)
! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('TRANSIN',4)

        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to TRANSIN'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF
      ENDIF ! End check on submodel

! DEPENDS ON: setcona_ctl
      CALL SETCONA_CTL(                                                 &
#include "artd1.h"
#include "artduma.h"
#include "artptra.h"
#include "artlndm.h"
     &  ICODE,CMESSAGE,isSTASH)

      IF (ICODE  /=  0) THEN
        WRITE(6,*) 'INITIAL: Failure in call to SETCONA'
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ICODE,Cmessage)
      ENDIF
#endif


!L
!L  4.5  Set up control block for writing interface fields.
!L
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INTF_CTL',3)
! DEPENDS ON: intf_ctl
      CALL  INTF_CTL (                                                  &
#include "artinfa.h"
#include "artduma.h"
#include "artptra.h"
     &                ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INTF_CTL',4)
        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to INTF_CTL'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF
!L
#if defined(ATMOS)
!L
!L 4.6  Initialise physical constants used in main physics
!L      packages - includes radiation lookup tables
!L

      IF (L_Physics) THEN

! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INITPHYS',3)
! DEPENDS ON: initphys
      CALL INITPHYS(ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INITPHYS',4)
        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to INITPHYS'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF

      ENDIF ! on L_Physics

      IF (L_emcorr) THEN ! Energy correction enabled
!L--------------------------------------------------------------------
!L 4.7  Initialise total atmospheric energy & energy correction
!L--------------------------------------------------------------------
!L  Only done for a new run and then only if the header values in the
!L  dump are set to missing data. If the arrays were set at the
!L  beginning of all NRUNS (new runs) then bit reproducibility
!L  would be lost for short reruns.
!L  The energy correction applied in any day comes from the
!L  change calculated for the previous day. Initialisation for a NRUN
!L  sets the energy correction to zero for the first day.
!L
!L    A_REALHD(18) - total energy flux into the atmosphere
!L                   (used to accumulate change in energy throughout
!L                   a day from physics). Value at start of run zero
!L    A_REALHD(19) - total mass of the atmosphere (wet + dry)
!L                   Note this is not conserved by the dynamics
!L                   The model from UM 5.0 only conserves dry mass.
!L    A_REALHD(20) - total energy of the atmosphere calculated from
!L                   fields in dump.
!L    A_REALHD(21) - energy correction evaluated from previous day of
!L                   run (ie previous run or needs setting to zero).

      IF (STEPim(a_im) == 0 ) THEN

! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INEMCORR',3)

! DEPENDS ON: init_emcorr
         CALL INIT_EMCORR(                                              &
#include "artd1.h"
#include "artduma.h"
#include "artptra.h"
#include "artlndm.h"
     &                    ICODE,CMESSAGE)

! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INEMCORR',4)
        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to INIT_EMCORR'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF

        END IF      ! end test on T+0
!
      END IF    !    LEMCORR
#endif
!L----------------------------------------------------------------------
!L 5. Set timestep group control switches for initial step
!L
!L 5.1 Set timestep control switches for initial step
!L
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('SETTSCTL',3)
! DEPENDS ON: settsctl
      CALL SETTSCTL (                                                   &
#include "artduma.h"
#include "artsts.h"
#include "artinfa.h"
     &             internal_model,.TRUE.,MEANLEV,ICODE,CMESSAGE,        &
     &             isSTASH)   ! FLUME-STASH UM vn7.1
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('SETTSCTL',4)
        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to SETTSCTL'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF
!L
!L 5.2 Initialise PP files at step 0
!L
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('PPCTL_INIT',3)
! DEPENDS ON: ppctl_init
      CALL PPCTL_INIT(                            &
#include "artduma.h"
#include "artinfa.h"
     &           submodel,PPNAME,ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('PPCTL_INIT',4)
        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to PPCTL'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF
!L
!L 5.3  Initialise assimilation package (not if assimilation completed)
!L
#if defined(A18_2A) 
      IF ( (ASSIM_STEPSim(a_im)+ASSIM_EXTRASTEPSim(a_im) >  0  .AND.    &
     &   (MODEL_ASSIM_MODE == "Atmosphere"       .OR.                   &
     &    MODEL_ASSIM_MODE == "Coupled   ")      .AND.                  &
     &   (RUN_ASSIM_MODE   == "Atmosphere"       .OR.                   &
     &    RUN_ASSIM_MODE   == "Coupled   ")      .AND.                  &
     &    STEPim(a_im)  <   ASSIM_FIRSTSTEPim(a_im) +                   &
     &              ASSIM_STEPSim(a_im) + ASSIM_EXTRASTEPSim(a_im))     &
     & ) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('IN_ACCTL',3)
! DEPENDS ON: in_acctl
      CALL IN_ACCTL(                                                    &
#include "artduma.h"
#include "artptra.h"
#include "argppx.h"
     &                  ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('IN_ACCTL',4)
        IF (ICODE  /=  0) THEN
          WRITE(6,*) 'Failure in call to IN_ACCTL'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF
      ENDIF
#endif
!L
!L 5.4  Open unit for model increments diagnostics if requested
!L
#if defined(ATMOS)
      IF(LPRFLD)THEN
        LEN_FILENAME=LEN(FILENAME)
#if defined(MPP)
        CALL FORT_GET_ENV(FT_ENVIRON(NDEV_FLD),LEN_FT_ENVIR(NDEV_FLD),  &
     &                    FILENAME,LEN_FILENAME,ICODE)

        IF (ICODE  /=  0) THEN
          Cmessage='Failure in call to FORT_GET_ENV'
          WRITE(6,*) 'Attempting to get value of environment variable ',&
     &      FT_ENVIRON(NDEV_FLD)
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF
!         Search for end of filename
        LL=0
        DO I=1,LEN_FILENAME
          IF(FILENAME(I:I) /= ' ') THEN
             LL=LL+1
          ENDIF
        ENDDO    ! I over characters

!         Construct filename with PE no. appended
        FILENAME(LL+1:LL+1)='.'
        WRITE(FILENAME(LL+2:LL+5),'(i4.4)') mype
#endif

#if defined(MPP)
        LEN_FLD_FILENAME=LL+5
        FLD_FILENAME=FILENAME
        CALL OPEN_SINGLE(NDEV_FLD,FLD_FILENAME,                         &
     &                   LEN_FLD_FILENAME,1,1,ICODE)
#else
! DEPENDS ON: file_open
        CALL FILE_OPEN(NDEV_FLD,FT_ENVIRON(NDEV_FLD),                   &
     &                 LEN_FT_ENVIR(NDEV_FLD),1,0,ICODE)
#endif
        IF(ICODE /= 0) THEN
            Cmessage='Error opening cached file'
            WRITE(6,*) 'Attempting to open file on unit ',NDEV_FLD
! DEPENDS ON: ereport
            CALL Ereport(RoutineName,ICODE,Cmessage)
        ENDIF
      ENDIF
#endif

#if defined(ATMOS)
!L----------------------------------------------------------------------
!L
!L 6.   Call to temporal filtering control routine (section 18)
!L
      IF (L_IAU .OR. L_TDF) THEN
! DEPENDS ON: tfilt_cntl
        CALL TFilt_cntl (                                               &
       U, V, W, U_adv, V_adv, W_adv,                                    &
       Theta, Exner_rho_levels, Rho,                                    &
       Q, Qcl, Qcf, Murk, ozone_tracer,                                 &
       Deep_soil_temp,                                                  &
       P, P_theta_levels, Exner_theta_levels,                           &
       snodep,                                                          &
       cf_area, cf_bulk, cf_liquid, cf_frozen,                          &
       Pstar, Tstar, Tstar_tile,                                        &
#include "artduma.h"
#include "artlndm.h"
#include "argppx.h"
     &                    l_mr_tfiltctl,                                &
                                      ! use mixing ratio code ?
     &                    IAU_lookup,                                   &
                                      ! inout
     &                    D1_IAU_k4,                                    &
                                      ! inout
     &                    D1_TFilt )  ! inout

! 6.1   If required, write out TS0 dump, including any IAU increments
!       which may have just been added.

        IF (L_IAU .AND. STEPim(a_im) == 0) THEN
          IF (L_IAU_DumpTS0State) THEN

! DEPENDS ON: dumpctl
          CALL DUMPCTL (                                                &
#include "artd1.h"
#include "artduma.h"
#include "artlndm.h"
#include "artptra.h"
#include "artsts.h"
#include "argppx.h"
     &                   atmos_sm, MEANLEV, .FALSE., '           ', 0,  &
     &                   ICODE, CMESSAGE)

            IF (ICODE /= 0) THEN
              WRITE (6,*) 'Failure in call to DUMPCTL (Atmos)'
! DEPENDS ON: ereport
              CALL Ereport(RoutineName, ICODE, CMESSAGE)
            END IF

          END IF
        END IF

      END IF

#endif
!L----------------------------------------------------------------------

#if defined(ATMOS)
!L
!L 7.1  Get derived diagnostics from start fields (atmosphere)
!L
      IF (STEPim(a_im) == 0) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INITDIAG',3)
! DEPENDS ON: initdiag
        CALL INITDIAG(                                                  &
#include "artd1.h"
#include "artduma.h"
#include "artsts.h"
#include "arg_atm_fields.h"
#include "artlndm.h"
#include "argppx.h"
     & Dummyarg)

! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INITDIAG',4)
          IF (ICODE  /=  0) THEN
            WRITE(6,*) 'Failure in call to INITDIAG'
! DEPENDS ON: ereport
            CALL Ereport(RoutineName,ICODE,Cmessage)
          ENDIF
      ENDIF
!
! 7.2 Code to update the boundaries now moved forwards to 4.4.1
!
!L
!L 7.3 Update ancillary fields in dump if start time corresponds to
!L     an ancillary field update time. Also done at T+0 with values
!L     updated to half a period back from first standard update time
!L     to ensure reproducibility between long runs and new runs
!L     started from dump at any time.
!L
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('UP_ANCIL',3)
      IF (ANCILLARY_STEPSim(a_im) >  0) THEN
        IF (STEPim(a_im) == 0 .OR.                                      &
     &      MOD(STEPim(a_im),ANCILLARY_STEPSim(a_im)) == 0)             &
! DEPENDS ON: up_ancil
     &   CALL UP_ANCIL (                                                &
#include "artd1.h"
#include "artduma.h"
#include "artptra.h"
#include "artanc.h"
     &                  submodel,                                       &
#include "argppx.h"
     &                  ICODE,CMESSAGE)
      ENDIF
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('UP_ANCIL',4)
          IF (ICODE  /=  0) THEN
            WRITE(6,*) 'Failure in call to UP_ANCIL'
! DEPENDS ON: ereport
            CALL Ereport(RoutineName,ICODE,Cmessage)
          ENDIF
!L
!L
!L 7.3.1 Initialize tiled prognostics, gridbox mean vegetation
!L       parameters and TRIFFID accumulation prognostics.
!L
!      L_VEG_FRACS=.FALSE. ! Hardwire until it's put back in CNTLATM at
!                          ! vn5.3 or 5.4. Comes originally from umui.
!                         ! The same applies to L_TRIFFID (unused here).
      IF (L_VEG_FRACS .AND. STEPim(a_im) == 0 ) THEN
! DEPENDS ON: timer
                                 IF(LTIMER) CALL TIMER('INIT_VEG',3)
!  Skip INIT_VEG if LAND_FIELD=0 for this PE.
        IF (LAND_FIELD  >   0) THEN
! DEPENDS ON: init_veg
          CALL INIT_VEG(STEPim(a_im),                                   &
#include "artd1.h"
#include "artduma.h"
#include "artptra.h"
     &                  ICODE,CMESSAGE)
          IF (ICODE  /=  0) THEN
            WRITE(6,*) 'Failure in call to INIT_VEG'
! DEPENDS ON: ereport
            CALL Ereport(RoutineName,ICODE,Cmessage)
          ENDIF
        ELSE
          WRITE(6,*)'INITIAL; skip INIT_VEG, LAND_FIELD=0 for this PE'
        END IF
! DEPENDS ON: timer
                                 IF(LTIMER) CALL TIMER('INIT_VEG',4)
      END IF
!L 7.3.2 Ensure that canopy water does not exceed canopy
!L       capacity at step zero (this may be a problem when
!L       using interpolated fields
!L
      IF (STEPim(a_im) == 0) THEN
! DEPENDS ON: timer
                                 IF(LTIMER) CALL TIMER('INIT_HYD',3)
      IF (.NOT.L_VEG_FRACS) THEN
#if defined(MPP)
!  Skip INIT_HYD if LAND_FIELD=0 for this PE.
      IF (LAND_FIELD  >   0) THEN
#endif
! DEPENDS ON: init_hyd
         CALL INIT_HYD(                                                 &
#include "artd1.h"
#include "artptra.h"
     &                 ICODE,CMESSAGE)
          IF (ICODE  /=  0) THEN
            WRITE(6,*) 'Failure in call to INIT_HYD'
! DEPENDS ON: ereport
            CALL Ereport(RoutineName,ICODE,Cmessage)
          ENDIF
#if defined(MPP)
      ELSE
      write(6,*)' INITIAL; skip INIT_HYD, LAND_FIELD=0 for this PE'
      END IF
#endif
      ENDIF
! DEPENDS ON: timer
                                 IF(LTIMER) CALL TIMER('INIT_HYD',4)
      END IF
!
!L
!L 7.3.3 Ensure that convective cloud cover and liquid water path
!L       are consistent with convective cloud base & top. (Corrects
!L       for occasional problems caused by reconfiguration.)
!L
      IF (STEPim(a_im) == 0) THEN
! DEPENDS ON: timer
                                 IF(LTIMER) CALL TIMER('INIT_CNV',3)
! DEPENDS ON: init_cnv
         CALL INIT_CNV(                                                 &
#include "artd1.h"
#include "artptra.h"
     &                 ICODE,CMESSAGE)
          IF (ICODE  /=  0) THEN
            WRITE(6,*) 'Failure in call to INIT_CNV'
! DEPENDS ON: ereport
            CALL Ereport(RoutineName,ICODE,Cmessage)
          ENDIF
! DEPENDS ON: timer
                                 IF(LTIMER) CALL TIMER('INIT_CNV',4)
      END IF
!
! River routing
       IF(L_RIVERS)THEN
! Initialise the step counter for river routing.
        IF (STEPim(a_im) == 0) THEN
! DEPENDS ON: timer
                                 IF(LTIMER) CALL TIMER('INIT_RIV',3)
! DEPENDS ON: init_riv
         CALL INIT_RIV(                                                 &
#include "artd1.h"
#include "artduma.h"
#include "artptra.h"
     &                 ICODE,CMESSAGE)
          IF (ICODE  /=  0) THEN
! DEPENDS ON: ereport
            CALL Ereport(RoutineName,ICODE,Cmessage)
          ENDIF
! DEPENDS ON: timer
                                 IF(LTIMER) CALL TIMER('INIT_RIV',4)
        END IF
      ENDIF                   ! L_RIVERS
!
!L
!L 7.3.4 ! initialization of radiative feedback
!L  Initialise and check address array to feed chemical tracers from
!L  UKCA into radiation scheme.  This needs to be done even for a CRUN
!L  so no need to check for STEPim(a_im) == 0.
!L
      grgas_addr = -1
      IF (L_UKCA) THEN
! DEPENDS ON: timer
                                 IF(LTIMER) CALL TIMER('INIT_RADUKCA',3)
! DEPENDS ON: init_radukca
        CALL INIT_RADUKCA(                                             &
#include "artd1.h"
#include "artptra.h"
     &       ngrgas,grgas_addr)
! DEPENDS ON: timer
                                 IF(LTIMER) CALL TIMER('INIT_RADUKCA',4)
      END IF
!!
!L
!L 7.4 Generate interface fields at step zero if required
!L
      IF (LINTERFACE .AND. STEPim(a_im) == 0) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('GEN_INTF',3)
! DEPENDS ON: gen_intf
        CALL GEN_INTF (                                                 &
#include "artd1.h"
#include "artduma.h"
#include "artsts.h"
#include "artptra.h"
#include "artinfa.h"
#include "argppx.h"
     &          submodel,ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('GEN_INTF',4)
          IF (ICODE  /=  0) THEN
            WRITE(6,*) 'Failure in call to GEN_INTF'
! DEPENDS ON: ereport
            CALL Ereport(RoutineName,ICODE,Cmessage)
          ENDIF
      ENDIF

#endif
!
!L----------------------------------------------------------------------
!L 8. If coupled model, initialise addresses of coupling fields,
!L    and if model has restarted at the end of a coupling period
!L    exchange coupling fields and swap data (full ocean model)
!L    or both models are at step 0, exchange coupling fields and
!L    swap data (in sense O-A at step 0).
!L
#if defined(ATMOS)

#if defined(ACCESS)
      IF (L_OASIS .and. l_auscom) THEN
#else
      IF (L_OASIS) THEN
#endif

      ! If Running with OASIS3 or OASIS4 we need to set up 
      ! pointers to relevant coupling fields. 

! DEPENDS ON: timer
      If (ltimer) Call timer('oasis_inita2o',3)
! DEPENDS ON: oasis_inita2o
      Call oasis_inita2o(                                               &
#include "artd1.h"
#include "artsts.h"
#include "artduma.h"
#include "artptra.h"
     &                icode,                                            &
     &                cmessage)
      If (icode/=0) Then
        Write(6,*) 'Failure in call to inita2nemo'
! DEPENDS ON: ereport
        Call Ereport(RoutineName,icode,cmessage)
      Endif
! DEPENDS ON: timer
      If (ltimer) Call timer('oasis_inita2o',4)

      ENDIF



#if defined(A26_1A)
! Set up the ATMOS grid lat./long.coords for regridding to river
! routing grid
      IF(L_RIVERS)THEN
! DEPENDS ON: init_a2t
        CALL INIT_A2T(                                                  &
#include "artd1.h"
#include "artsts.h"
#include "artduma.h"
#include "artptra.h"
#include "artatcpl.h"
     &              ICODE,CMESSAGE)
          IF (ICODE  /=  0) THEN
! DEPENDS ON: ereport
            CALL Ereport(RoutineName,ICODE,Cmessage)
          ENDIF
      ENDIF
#endif
#endif
!L----------------------------------------------------------------------
!L 9. Print formatted diagnostics from initial dump
!L

! Printctl: printing of climate global and zonal diagnostics is no
!           longer supported

!L----------------------------------------------------------------------
!L 10. Initialisation complete - return to master routine
!L
! Check that operational model running on MPP has finished
! initialisation and write a message to the operator
#if defined(MPP)
      IF(mype == 0) THEN
         IF(MODEL_STATUS  ==  'Operational') THEN
! DEPENDS ON: operatormessage
            CALL OperatorMessage(nproc)
         ENDIF
      ENDIF
#endif
 999  CONTINUE
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INITIAL ',4)
      RETURN
!L----------------------------------------------------------------------
      END SUBROUTINE INITIAL
#endif
