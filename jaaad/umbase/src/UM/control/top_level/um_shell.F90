#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine: UM_SHELL -------------------------------------------------
!LL
!LL  Purpose: Control routine for the Atm Model.
!LL           Acquires size information needed for dynamic allocation of
!LL           configuration-dependent arrays and calls U_MODEL (the
!LL           master control routine) to allocate the arrays and perform
!LL           the top-level control functions and timestepping.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 6.1.5A
!LL
!LL  Programming standard: UM Doc Paper 3, version 8 (01/06/2007)
!LL
!LL  Logical components covered: C0
!LL
!LL  Project task: C0
!LL
!LL  External documentation: On-line UM document C1 - The top-level
!LL                          dynamic allocation
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!     FLUME-STASH  um_shell coverted from program to subroutine
!
      Subroutine UM_SHELL(isSTASH)

#if defined(FLUME)
      Use MSTASHFlumeModel, Only: MSTASHFlumeModel_initGatherData
      Use flume
#endif
      Use flumerun
!
#if defined(OASIS3)
      USE mod_prism
      USE oasis3_atmos_init
#endif
#if defined(ACCESS)
      USE auscom_cpl_data_mod, ONLY : l_auscom
#endif


!*----------------------------------------------------------------------
      IMPLICIT NONE
!
!  Subroutines called
!
#if defined(T3E)
      External                                                          &
     & PXFCONST,PXFSYSCONF
#endif
!
!  Local parameters
!
      CHARACTER*(*) RoutineName
      PARAMETER (RoutineName = 'UM_SHELL')
!
!  Local variables
!
#if defined(CRI_FFIO)
      real secondr, close_time
#endif
      INTEGER ICODE             ! Work - Internal return code
      INTEGER(kind=ip_intwp_p) ::  ICODE1
      INTEGER ISTATUS       ! RETURN STATUS FROM OPEN
      Integer loop_pe

      CHARACTER*80 FILENAME ! RETURN FILENAME FROM GET_FILE
      CHARACTER*256 CMESSAGE    ! Work - Internal error message
      INTEGER                                                           &
#if defined(ATMOS)
     &        atm_nprocx,                                               &
                              ! number of procs EW for atmosphere
     &        atm_nprocy,                                               &
                              ! number of procs NS for atmosphere
#endif
     &        err             ! error return from FORT_GET_ENV
      INTEGER :: nproc_um_npes       ! the verson of nproc from UM_NPES

      CHARACTER*8 c_nproc            ! to get nproc_x and nproc_y from
!                                    ! environment variables.
      CHARACTER*100 parexe_env       ! to hold the name of the parallel
!                                    ! executable script
      CHARACTER*200 stdout_filename  ! file to write stdout to
      CHARACTER*180 stdout_basename  ! base of filename
      CHARACTER*170 dataw_char       ! value of $DATAW
      CHARACTER*5  runid_char        ! value of $RUNID
#ifdef ACCESS
      CHARACTER*8  env_auscom       ! value of $RUNID
      CHARACTER*8  auscom_flag       ! value of $RUNID
#endif

#include "mppconf.h"
#include "decomptp.h"
#include "parvars.h"
#include "gccom.h"
!
      integer um_nam_max_seconds
!
      character*8 c_nam_max_seconds
!
!  Configuration-dependent sizes for dynamic arrays
!
#include "csubmodl.h"
#include "typsize.h"
! For STASH sizes
#include "typstsz.h"
!
!  Super array sizes for dynamic allocation in U_MODEL
!
#include "typszsp.h"
#include "typszspa.h"
#include "typszspc.h"
!
!   Declare ppxref look-up arrays, pointer array, and associated
!     sizes. The lengths (ppxRecs) of the ppx look-up arrays are
!     dynamically allocated.
#include "cppxref.h"
#include "version.h"
#include "cstash.h"
! Length of ppx look-up arrays - for dynamic allocation
      INTEGER ppxRecs
#include "chsunits.h"
#include "ccontrol.h"
#include "cntl_io.h"
!L
!L Control for temporal filtering.
!L
#include "ctfilt.h"
!L
!
!  Localized sizes for ocean decomposition:
      Integer                                                           &
     &     row_length_oce , rows_oce

      integer sect_err, rnl_err, um_rnl_skip
      integer io_timer_activate          ! Switch for IO timer.
!
      character*8 c_um_sector_size, c_um_rnl_skip
      character*8 ch_date2   !  Date returned from date_and_time
      character*10 ch_time2  !  Time returned from date_and_time

#if defined(T3E)
#include "t3eclktk.h"
      INTEGER iclktck,ierr
#endif
!   Fortran unit numbers
      INTEGER NFTPPXREF
      INTEGER NFTSTMSTU
      DATA NFTPPXREF/22/,NFTSTMSTU/2/
! Variables needed to close all the files

      INTEGER :: I
      INTEGER :: UNIT_STATUS
      CHARACTER*8 dummy_name

#if defined(FLUME)
! FLUME-STASH
! Communicator containing this process
      INTEGER :: commWorld
#endif
! Flag to indicate whether this is STASH process 
!   (subroutine argument - is set to FALSE in non-Flume run)
      LOGICAL :: isSTASH

#if defined(OASIS3)
      integer :: comm_in  ! Local communicator defined by OASIS      
      integer :: dummy_comm     ! Dummy communicator for OASIS
      integer  :: comp_id       ! Component ID for OASIS
      character*32 :: grid_name ! Grid name for OASIS
      integer  :: out_unit      

      ! gol124: use own mpi error handler to enforce coredumps
      integer(kind=4) :: l_comm, core_dumper, ierr32

      external mpi_errors_coredump
#endif

      cmessage = ' '
!
#if defined(FLUME)
      IF (Flume_run) THEN
        print *,"Start of UM_Shell: Current active model info is..."
        CALL printActiveModelDetails()
      END IF
#endif


! For IBM platform, turn on signal handling
#if defined(IBM)
      CALL signal_trap(0)
#endif
!
!L----------------------------------------------------------------------
!L 0. Start Timer running
!L
#if defined(T3E)
! Find out the number of clock ticks per second on this machine.
! This information is required to calculate the wallclock times
! in TIMER

      iclktck=0
      ticks_per_second=0

      CALL PXFCONST('CLK_TCK',iclktck,ierr)

      IF (ierr  /=  0) THEN
        WRITE(6,*) 'UMSHELL : Failure in PXFCONST, err= ',ierr
        ICODE=1
        GO TO 9999
      END IF

      CALL PXFSYSCONF(iclktck,ticks_per_second,ierr)

      IF (ierr  /=  0) THEN
        WRITE(6,*) 'UMSHELL : Failure in PXFSYSCONF, err= ',ierr
        ICODE=1
        GO TO 9999
      END IF

#endif
!   Open file for UNIT 5 before initialisation of model. All runtime
!   control variables subsequently read in from UNIT 5 by namelist.
      CALL GET_FILE(5,FILENAME,80,ICODE)
      OPEN(5,FILE=FILENAME,IOSTAT=ISTATUS)
      IF(ISTATUS /= 0) THEN
        ICODE=500
        WRITE(6,*) ' ERROR OPENING FILE ON UNIT 5'
        WRITE(6,*) ' FILENAME =',FILENAME
        WRITE(6,*) ' IOSTAT =',ISTATUS
        GO TO 9999
      END IF
!L------------------------------------------------------------------
!L 0.1 Get submodel/internal model components of model run.
!L
      ICODE=0
! DEPENDS ON: um_submodel_init
      CALL UM_Submodel_Init(ICODE)

!L----------------------------------------------------------------------
!----------------------------------------------------------------------
! 1.0 Initialise Message Passing Libraries
!

#if defined(ATMOS)
! Get the atmosphere decomposition

      CALL FORT_GET_ENV('UM_ATM_NPROCX',13,c_nproc,8,err)
      IF (err  /=  0) THEN
        WRITE(6,*) 'Warning: Environment variable UM_ATM_NPROCX has ',  &
     &             'not been set.'
        WRITE(6,*) 'Setting nproc_x to 1'
        atm_nprocx=1
      ELSE
        READ(c_nproc,'(I4)') atm_nprocx
      ENDIF
      CALL FORT_GET_ENV('UM_ATM_NPROCY',13,c_nproc,8,err)
      IF (err  /=  0) THEN
        WRITE(6,*) 'Warning: Environment variable UM_ATM_NPROCY has ',  &
     &             'not been set.'
        WRITE(6,*) 'Setting nproc_y to 1'
        atm_nprocy=1
      ELSE
        READ(c_nproc,'(I4)') atm_nprocy
      ENDIF
#endif

! Find out the maximum number of processors to be used in this
! run of the model

      CALL FORT_GET_ENV('UM_NPES',7,c_nproc,8,err)
      IF ( (err  /=  0) .OR. (c_nproc  ==  'UNSET') ) THEN
        WRITE(6,*) 'Error : Environment variable UM_NPES has ',         &
     &             'not been set.'
        WRITE(6,*) 'Exiting'
        GO TO 9999
      END IF

      READ(c_nproc,'(I4)') nproc_max

! Check MAXPROC is big enough for nproc_max

      IF (nproc_max  >   MAXPROC) THEN
        WRITE(6,*) 'Error : MAXPROC is not big enough.'
        WRITE(6,*) 'You will need to edit the parameter in comdeck ',   &
     &             'PARPARM.'
        WRITE(6,*) 'MAXPROC= ',MAXPROC,' nproc_max= ',nproc_max
        WRITE(6,*) 'Exiting'
        GO TO 9999
      END IF

#if defined(ATMOS)
! Check that there are enough processors to support the
! decompositions that have been requested.

#if defined(OASIS3) || defined(OASIS4)
      ! Ignore the cross check against the total number of PES
      ! since this is not a valid test when running with OASIS.
      nproc_max = atm_nprocx*atm_nprocy
#endif

      IF ((atm_nprocx*atm_nprocy)  /=  nproc_max ) THEN
        WRITE(6,*) 'Error : Atmosphere decomposition of ',atm_nprocx,   &
     &             ' x ',atm_nprocy,' processors cannot be supported ', &
     &             'using ',nproc_max,' processors.'
        WRITE(6,*) 'Exiting'
        GO TO 9999
      END IF
#endif

#if !defined(T3D) && !defined(T3E) && !defined(NEC)
      CALL FORT_GET_ENV('PAREXE',6,parexe_env,100,err)
      IF (err  /=  0) THEN
        WRITE(6,*) 'Failed to get the name of the parallel executable ',&
     &             'script from $PAREXE'
        WRITE(6,*) '*** Model Exiting.'
        GO TO 9999
      END IF
#else
      parexe_env=' '
#endif

#if defined(FLUME)
      IF (Flume_run) THEN
        ! FLUME-STASH MPI initialisation done in FLUME layer
        ! Call gc_init_final to set up values required by GCOM 
        commWorld=flumeGetCommunicator()
        print *,"flumegetcommunicator: commworld ",commWorld
        call gc_init_final(mype,nproc_max,commWorld)
      ELSE
        ! Non-Flume UM
        CALL GC_INIT(parexe_env,mype,nproc_max)        
      END IF
#elif defined(OASIS3) || defined(OASIS4)
!----------------------------------------------------------------------
!L
!L Call routine to initialise OASIS
!L
!----------------------------------------------------------------------
      comm_in=-999
#if defined(ACCESS)
      ! goal: same executable can be used for stand-alone and
      ! coupled runs
      ! L_oasis is defined in namelist files but not available yet
      ! use an environment variable instead
      CALL FORT_GET_ENV('AUSCOM_CPL',10,env_auscom,8,err)
      auscom_flag = trim(env_auscom)
      IF ( (auscom_flag ==  'TRUE') .or. (auscom_flag == 'true') ) THEN
        L_auscom=.TRUE.
! l_oasis will be set below by reading namelists below um_setup
! it is not needed before that, so we can leave it untouched
!        L_oasis=.true.
      else
        l_auscom=.false.
      END IF
!      IF ( (err  /=  0) .OR. (auscom_flag ==  'UNSET') ) THEN
!        L_OASIS=.FALSE.
!      ELSE
!        L_OASIS=.TRUE.
!      END IF
#endif
      L_OASIS=.TRUE.
  
      ! The key thing here is to get hold of the
      ! communicator defined for us by PRISM and then
      ! use that in GCOM rather than letting GCOM define
      ! its own MPI_COMM_WORLD.

#if defined(OASIS3)
#if defined(ACCESS)
      IF (L_auscom) THEN
! DEPENDS ON: OASIS3_ATMOS_INIT
        CALL OASIS3_UM_INIT(comm_in,ICODE,CMESSAGE)
      END IF
#else
! DEPENDS ON: OASIS3_ATMOS_INIT
      CALL OASIS3_UM_INIT(comm_in,ICODE,CMESSAGE)
#endif
#endif

#if defined(OASIS4)
! DEPENDS ON: OASIS4_ATMOS_INIT
      CALL OASIS4_UM_INIT(comm_in,ICODE,CMESSAGE)
#endif
      ! Check that MPI (or other) communication method
      ! is initialised. Discard the communicator returned
      ! from this call since we'll use an OASIS defined one.
#ifdef ACCESS
      ! do normal init in stand-alone mode
      IF (L_auscom) THEN
        CALL GC_INIT_INTRO(DUMMY_COMM)
      ELSE
        CALL GC_INIT_INTRO(comm_in)
      END IF
#else
      CALL GC_INIT_INTRO(DUMMY_COMM)
#endif

      ! Do all the initialisation with the correct communicator
      CALL GC_INIT_FINAL(mype,nproc_max,comm_in)

      WRITE(6,*) "UM_SHELL: GCOM for OASIS", mype,nproc_max,comm_in
#ifdef ACCESS
      IF (L_auscom) THEN
        WRITE(6,*) "UM_SHELL: running in coupled mode"
      ELSE
        WRITE(6,*) "UM_SHELL: running in stand alone mode"
      END IF

      ! gol124: use own mpi error handler to enforce coredumps
      l_comm = comm_in
      !depends on: mpi_errors_coredump
      call MPI_Comm_create_errhandler (mpi_errors_coredump,             &
                   core_dumper, ierr32)
      call MPI_Comm_set_errhandler (l_comm, core_dumper, ierr32)
#endif

#else
! Standard UM GCOM initialisation when OASIS is not used

      nproc_um_npes = nproc_max
      CALL GC_INIT(parexe_env,mype,nproc_max)

! Check the number of processors in UM_NPES matches that obtained
! by GCOM
      IF (nproc_max /= nproc_um_npes) THEN
       
!DEPENDS ON: ereport
        ICODE = 100
        WRITE(CMESSAGE,*) 'UM started on ', nproc_max,              &
                          ' PEs but ', nproc_um_npes,               &
                          ' asked for. Please adjust decomposition'
        CALL EREPORT(ROUTINENAME,ICODE,CMESSAGE)
      END IF
#endif


      IF (nproc_max  <   0) THEN
        WRITE(6,*) 'Parallel initialisation failed'
        GO TO 9999
      ELSE

! Send output to unique filename on every PE

        CALL FORT_GET_ENV('UM_STDOUT_FILE',14,stdout_basename,180,err)
        IF (err  /=  0) THEN
! Environment variable UM_STDOUT_FILE has not been set, so we will
! construct a default stdout_basename of $DATAW/$RUNID.fort6.pe
          CALL FORT_GET_ENV('DATAW',5,dataw_char,170,err)
          IF (err  /=  0) THEN
            WRITE(6,*) 'UMSHELL : Failed to get value of $DATAW'
            WRITE(6,*) '*** Model Exiting.'
            GO TO 9999
          END IF
          CALL FORT_GET_ENV('RUNID',5,runid_char,5,err)
          IF (err  /=  0) THEN
            WRITE(6,*) 'UMSHELL : Failed to get value of $RUNID'
            WRITE(6,*) '*** Model Exiting.'
            GO TO 9999
          END IF
          stdout_basename=trim(dataw_char)//'/'//                       &
     &                    trim(runid_char)//'.fort6.pe'
        ENDIF

! Now add PE number (mype) to stdout_basename to get the complete
! stdout_filename for this PE.

        IF (mype  <   10) THEN
          WRITE(stdout_filename,'(A,I1)') trim(stdout_basename),mype
        ELSEIF (mype  <   100) THEN
          WRITE(stdout_filename,'(A,I2)') trim(stdout_basename),mype
        ELSEIF (mype  <   1000) THEN
          WRITE(stdout_filename,'(A,I3)') trim(stdout_basename),mype
        ELSE
          WRITE(stdout_filename,'(A,I4)') trim(stdout_basename),mype
        ENDIF

! and close unit 6, then reopen to new filename

! FLUME-STASH Do not want open/close unit 6 in STASH proc
!   isSTASH defaults to FALSE in non-Flume run

        IF (.not.isSTASH) CLOSE(6)
        IF (.not.isSTASH) OPEN(6,FILE=stdout_filename)
! Force a close with a delete action - so if there is an existing
! unit6 output file it will be deleted, and the output from this
! run will go to a fresh file
        IF (.not.isSTASH) CLOSE(6,STATUS='DELETE')

        IF (.not.isSTASH) OPEN(6,FILE=stdout_filename)

        WRITE(6,*) nproc_max,' Processors initialised.'
        WRITE(6,*) 'I am PE ',mype
      ENDIF
!
! DEPENDS ON: timer
      IF (.not.isSTASH) CALL TIMER('UM_SHELL',1)

      if (mype == 0) then
        call date_and_time(ch_date2, ch_time2)
        write(6,*) 'Start of UM Job : ',                                &
     &  ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),' on ',       &
     &  ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4)
      endif


!L
!L Get the current sector size for disk I/O
!L

      CALL FORT_GET_ENV('UM_SECTOR_SIZE',14,c_um_sector_size,8,sect_err)
      IF (sect_err  /=  0) THEN
        WRITE(6,*) 'Warning: Environment variable UM_SECTOR_SIZE has ', &
     &             'not been set.'
        WRITE(6,*) 'Setting um_sector_size to 2048'
        um_sector_size=2048
      ELSE
        READ(c_um_sector_size,'(I4)') um_sector_size
      ENDIF
#if defined(CRAY)

!L
!L Get the current NAMELIST Skip value
!L

      CALL FORT_GET_ENV('UM_RNL_SKIP',11,c_um_rnl_skip,8,rnl_err)
      IF (rnl_err  /=  0) THEN
        WRITE(6,*) 'Warning: Environment variable UM_RNL_SKIP has ',    &
     &             'not been set.'
        WRITE(6,*) 'Setting um_rnl_skip to 0 - Omit Skipped Messages'
        um_rnl_skip=0
      ELSE
        READ(c_um_rnl_skip,'(I4)') um_rnl_skip
      ENDIF
      call rnlskip(um_rnl_skip)
#endif

!----------------------------------------------------------------------
!L
!L    Open unit 8 for server requests and send wakeup message
!L
      IF (mype == 0) THEN
        CALL GET_FILE(8,FILENAME,80,ICODE)
        OPEN(8,FILE=FILENAME)
        WRITE(8,10)
      ENDIF
   10 FORMAT('** WAKEUP **')
!----------------------------------------------------------------------
!L 1.1 Get configuration-dependent sizes needed for dynamic allocation.
!L
! DEPENDS ON: readsize
      CALL READSIZE()

!L   Read history and control files for NRUN; also interim control
!L   file for CRUN, and housekeeping file for operational run.
!L
! DEPENDS ON: um_setup
      CALL UM_SETUP(                                                    &
     &              ICODE,CMESSAGE)

      IF (ICODE >  0) GO TO 9999

#ifdef ACCESS
! supported modes:
! l_oasis       l_auscom        mode
! .false.       .false.         ATM only with ATM settings
! .true.        .false.         ATM only with ACCESS settings
! .true.        .true.          ACCESS mode with ACCESS settings
      IF (L_auscom.and.l_oasis) THEN
        WRITE(6,*) "UM_SHELL: ACCESS mode"
      ELSE IF (l_oasis.and.(.not.l_auscom)) THEN
        WRITE(6,*) "UM_SHELL: ATM mode with ACCESS settings"
      ELSE IF (.not.(l_oasis.or.l_auscom)) THEN
        WRITE(6,*) "UM_SHELL: ATM mode"
      ELSE
        WRITE(6,*) "UM_SHELL: l_oasis/l_auscom invalid combination"
        GOTO 9999
      END IF
#endif

!L----------------------------------------------------------------------
! Initialise print status for standard output
! DEPENDS ON: initprintstatus
      CALL InitPrintStatus

! Initialise IO Timer
      IF ( L_IO_Timer ) THEN
        io_timer_activate = 1
      ELSE
        io_timer_activate = 0
      END IF
      CALL io_timing_init( io_timer_activate )

! Initialise PP buffers
! DEPENDS ON: init_buffers_pp
         Call init_buffers_pp()
#if defined(ATMOS)

! Decompose atmosphere data and find new local data size


! DEPENDS ON: decompose_atmos
      CALL DECOMPOSE_ATMOS(global_ROW_LENGTH,global_ROWS,MODEL_LEVELS,  &
     &                     RIVER_ROWS, RIVER_ROW_LENGTH,                &
     &                     model_domain,                                &
     &                     atm_nprocx, atm_nprocy,                      &
     &                     extended_halo_size_EW,                       &
     &                     extended_halo_size_NS,                       &
     &                     RIMWIDTHA, NRIMA_MAX,ROW_LENGTH,ROWS )

#if defined(FLUME)
      ! FLUME-STASH 
      ! Note - isSTASH defaults to false in a non-Flume run
      IF (isSTASH) THEN
         ! Switch STASH to look like it's running on pe 0 on its own.
         nproc=1

         ! Must set mype to 0 before calling change_decomp because it
         ! populates the array based on this value. Just get rubbish
         ! if you pass actual PE of STASH process
         mype=0
         nproc_max=1

         CALL CHANGE_DECOMPOSITION(decomp_standard_atmos,ICODE)

         ! Now store the values for gathering decomposed data
         !  in the Flume layer.
         CALL MSTASHFlumeModel_initGatherData()

         ! Now switch STASH to a single processor decomposition 
         atm_nprocx=1
         atm_nprocy=1
         CALL DECOMPOSE_ATMOS(global_ROW_LENGTH,global_ROWS,         &
     &                        MODEL_LEVELS,                          &
     &                        RIVER_ROWS, RIVER_ROW_LENGTH,          &
     &                        model_domain,                          &
     &                        atm_nprocx, atm_nprocy,                &
     &                        extended_halo_size_EW,                 &
     &                        extended_halo_size_NS,                 &
     &                        RIMWIDTHA, NRIMA_MAX,ROW_LENGTH,ROWS)
         ! Set current decomposition to unset in order to force 
         ! CHANGE_DECOMPOSITION to update the arrays. Otherwise
         ! it just returns because it thinks it's already in that
         ! decomposition.
         current_decomp_type=DECOMP_UNSET
      END IF
#endif

! Set up the atmosphere decomposition in PARVARS
! DEPENDS ON: change_decomposition
      CALL CHANGE_DECOMPOSITION(decomp_standard_atmos,ICODE)

      IF (ICODE  /=  0) GO TO 9999
! Output range of gridpoints for each PE
      if (MYPE  ==  0) then
       do loop_pe = 0,ATM_NPROCX*ATM_NPROCY-1
       WRITE(6,*)'Range of gridpoints for PE',loop_pe,' :-'
       write(6,*)'EW direction: ',g_datastart_f(1,1,loop_pe),' to ',    &
     & g_datastart_f(1,1,loop_pe)+g_blsize(1,1,loop_pe)-1,              &
     & '  NS direction: ',g_datastart_f(2,1,loop_pe),' to ',            &
     & g_datastart_f(2,1,loop_pe)+g_blsize(2,1,loop_pe)-1
       end do
      else
       WRITE(6,*)'Range of gridpoints for PE',mype,' :-'
       write(6,*)'EW direction: ',g_datastart_f(1,1,mype),' to ',       &
     & g_datastart_f(1,1,mype)+g_blsize(1,1,mype)-1,                    &
     & '  NS direction: ',g_datastart_f(2,1,mype),' to ',               &
     & g_datastart_f(2,1,mype)+g_blsize(2,1,mype)-1
      end if

#endif

! Call DERVSIZE (the call in READSIZE has been deleted)

      ICODE=0
! DEPENDS ON: dervsize
      CALL DERVSIZE(ICODE,CMESSAGE)
      IF (ICODE  /=  0) GO TO 9999

#if defined(ATMOS)

!     Ensure that domain decomposition is set for Atmosphere
! DEPENDS ON: change_decomposition
      call change_decomposition (decomp_standard_atmos,icode)
      if (icode /= 0) then
        write (6,*) ' Error returned in CHANGE_DECOMPOSITION',          &
     &              ' before DERV_LAND_FIELD.'
        write (6,*) ' Error code ',icode
        write (cmessage,*) 'UM_SHELL : Error in CHANGE_DECOMPOSITION.'
        go to 9999   !  Exit
      endif

!     For MPP jobs, calculate the no of land-points on each PE.
! DEPENDS ON: derv_land_field
      CALL DERV_LAND_FIELD (21,icode,cmessage)
      if (icode >  0) then
        write (6,*) 'Error returned from DERV_LAND_FIELD.'
        write (6,*) 'Error code ',icode
        go to 9999   !  Exit
      endif

! Derive lengths involved with output boundary files - atmos.
! DEPENDS ON: derv_intf_a
      CALL DERV_INTF_A (TOT_LEN_INTFA_P,TOT_LEN_INTFA_U,                &
     &     MAX_INTF_MODEL_LEVELS,MAX_LBCROW_LENGTH,MAX_LBCROWS,         &
     &     N_INTF_A,U_FIELD_SIZE,U_FIELD_INTFA)

#endif
!
      IF (ICODE >  0) GO TO 9999
!-----------------------------------------------------------------------
! 1.2 Call STASH_PROC: top level control routine for processing of
!                      STASH requests and STASH addressing.

! Open STASHmaster file(s) and count number of records
!   This number is assigned to ppxRecs and used to dynamically
!   allocate the PPX_ arrays in which stash master records are held
      ppxRecs = 1
      ICODE   = 0
      IF (INTERNAL_MODEL_INDEX(A_IM) >  0)                              &
! DEPENDS ON: hdppxrf
     &   CALL HDPPXRF                                                   &
     &(NFTPPXREF,'STASHmaster_A',ppxRecs,ICODE,CMESSAGE)
      IF (ICODE /= 0) GO TO 9999
! Add number of user stash records
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(0,'             ',ppxRecs,ICODE,CMESSAGE)

      IF (icode  <   0) then
         IF (mype  ==  0) then
            write (0,*) 'WARNING : Problem in STASHmaster file(s)'
            WRITE (0,*) '        ',TRIM(cmessage)
         END IF
      ELSE IF (icode  >   0) then
         IF (mype  ==  0) then
            write (0,*) 'ERROR : Problem in STASHmaster files(s)'
            WRITE (0,*) '      ',TRIM(cmessage)
         END IF
         GO TO 9999  ! Always abort on fatal error.
      END IF

! DEPENDS ON: stash_proc
      CALL STASH_PROC(NFTPPXREF,NFTSTMSTU,.FALSE.,                      &
     &                ppxRecs,ICODE,CMESSAGE  )
      IF (ICODE >  0) GO TO 9999

! Total number of entries (N_PPXRECS) in STASH-addresses array IN_S has
!  obtained by WSTLST in STASH_PROC. Reset ppxRecs to equal this value.
! This is used to dynamically
!  allocate the ppx look-up arrays PPXI, PPXC in U_MODEL.

      ppxRecs = N_PPXRECS

!L----------------------------------------------------------------------
!L 1.3 Calculate addresses of super arrays passed down for dynamic
!L     allocation.
!L
      ICODE=0
! DEPENDS ON: um_index
      CALL UM_INDEX(                                                    &
#include "argszsp.h"
#include "argszspa.h"
#include "argszspc.h"
     &              ICODE,CMESSAGE)

      IF (ICODE >  0) GO TO 9999

!L----------------------------------------------------------------------
!L 1.4 Allow Override of namelist input in operational environment.
!L
! DEPENDS ON: oper_emergency
      CALL Oper_Emergency

!L----------------------------------------------------------------------
!L 1.5 Set up temporal filtering (or set defaults).
!L

! DEPENDS ON: setup_tfilt
      CALL Setup_TFilt

!L----------------------------------------------------------------------
!L 2. Call U_MODEL master routine to allocate the main data arrays
!L    and do the calculations.
!L
! DEPENDS ON: u_model
      CALL U_MODEL (                                                    &
     & NFTPPXREF,NFTSTMSTU,                                             &
#include "argszsp.h"
#include "argszspa.h"
#include "argszspc.h"
     &      ppxRecs,isSTASH)  ! FLUME-STASH
!
 9999 CONTINUE
      CLOSE(5)
!L----------------------------------------------------------------------
!L

! DEPENDS ON: timer
      if (.not.isSTASH) CALL TIMER('UM_SHELL',2)

      ! FLUME-STASH  
      ! For STASH procees do not execute the rest of the code otherwise
      !   files are closed before they've been written to.
      ! isSTASH defaults to FALSE in non-Flume run
      IF (isSTASH) then
         RETURN
      END IF      
!
! Make sure all the Files are properly Closed (and hence buffers
! are flushed correctly)
!
      DO I=20, NUNITS
        CALL IS_UNIT_OPEN(I,UNIT_STATUS)
        IF (UNIT_STATUS == 0) THEN   ! Unit is open
          ICODE=0
! DEPENDS ON: file_close
          CALL FILE_CLOSE(I, DUMMY_NAME, 5, 1, 0, ICODE)
          IF (ICODE /= 0) THEN
            WRITE (0,*)'UMSHELL1: Final Closing of Unit ',              &
     &       i, ' Returns Error Code ', icode, ' - Ignored'
            ICODE = 0
          END IF
          CALL GC_SSYNC(NPROC, ICODE)
        END IF
      END DO


      ! Write out PE 0 IO timings if required
      IF ( L_IO_Timer .AND. mype == 0) THEN
        Write (6,*) 'IO timings for PE 0'
        Call io_total_timings()
      END IF


      WRITE(6,*) 'Process ',mype,' has exited.'
#if defined(CRI_FFIO)
      call barrier()
      close_time=secondr()
      call close_all_files()
      call barrier()
      close_time=secondr()-close_time
      if (mype == 0) write(0,9976) close_time
      if (mype == 0) write(6,9976) close_time
9976  format(/'Time to Close All Files was ',f7.3,' Seconds'/)
#endif
!
      if (mype == 0) then
        call date_and_time(ch_date2, ch_time2)
        write(6,*) 'End of UM Job : ',                                  &
     &  ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),' on ',       &
     &  ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4)
      end if

      IF(ICODE /= 0) THEN
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ICODE,Cmessage)
      ENDIF


#if defined(OASIS3)
!----------------------------------------------------------------------
!
! Call prism routine to close OASIS3 instead of GC_EXIT
! 
!----------------------------------------------------------------------
#ifdef ACCESS
      IF (L_auscom) THEN
         CALL prism_terminate_proto(ICODE1)
         write(6,*)'UM_SHELL: Called prism_terminate_proto',ICODE1
      ELSE
         CALL GC_EXIT()
      END IF
#else
      CALL prism_terminate_proto(ICODE1)
      write(6,*)'UM_SHELL: Called prism_terminate_proto',ICODE1
#endif

#else
      ! Only bother calling GC_EXIT  when PRISM isnt going to shut
      ! things down for this component.
! Close down parallel process communication
!      FLUME-STASH  This is done in the Flume layer
       IF (.NOT.Flume_run) CALL GC_EXIT()
#endif

      END Subroutine UM_SHELL
#endif
