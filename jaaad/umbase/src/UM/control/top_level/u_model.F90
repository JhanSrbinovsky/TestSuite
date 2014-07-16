#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine: U_MODEL -----------------------------------------------
!LL
!LL  Purpose: High level control program for the Unified Model
!LL           (master routine).  Calls lower level control routines
!LL           according to top level switch settings. Called by
!LL           top level routine UMSHELL which provides dimension sizes
!LL           for dynamic allocation of data arrays.
!LL
!LL  Tested under compiler:   sxmpif90
!LL  Tested under OS version: 
!LL
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
      SUBROUTINE U_MODEL(                                               &
     &       NFT,NFTU,                                                  &
#include "argszsp.h"
#include "argszspa.h"
#include "argszspc.h"
     & ppxRecs, isSTASH) ! FLUME-STASH

! D1 replacement module
use atm_fields_mod

! OASIS Modules

#if defined(OASIS3)
      USE oasis3_atm_data_mod
#endif
#if defined(OASIS4)
      USE oasis4_atmos_init, ONLY : set_oasis_window32
#endif
#if defined(ACCESS)
      USE auscom_cpl_data_mod, ONLY : l_auscom, access_tfs, ocn_sss
#endif


! FLUME-STASH
#if defined(FLUME)
Use flume
Use SharedPPX  ! Declares allocatable versions of PPXI, PPXC, PPXPTR
Use MFlmModel
#endif   
Use flumerun       
      IMPLICIT NONE

!*L  Interface and arguments: ------------------------------------------
!L       Sizes of super arrays
#include "typszsp.h"
#include "typszspa.h"
#include "typszspc.h"
!L
!L       Model sizes
#include "parvars.h"
#include "typsize.h"
!L
!L       Addresses of component arrays within super arrays
#include "spindex.h"
!L
      LOGICAL isSTASH   ! true for FLUME parallel-STASH process
!L
!*----------------------------------------------------------------------
!
!  Common blocks
!
#include "cmaxsize.h"
#include "csubmodl.h"
#include "chsunits.h"
#include "chistory.h"
#include "ccontrol.h"
#include "ctime.h"
#include "cintfa.h"
#include "c_global.h"
!L
!  Dynamic allocation of ppxref look-up arrays and declaration of
!                                            ppxref pointer array.
#include "cppxref.h"
#if defined(FLUME)
! FLUME-STASH Include ppxlook explicitly so it can be altered to
!  avoid clash with module SharedPPX, which declares PPXI, PPXC, PPXPTR
! COMDECK PPXLOOK
#include "version.h"
! No. of STASH items per section
      INTEGER      PPXREF_ITEMS
        PARAMETER (PPXREF_ITEMS    =NITEMP)
! No. of STASH sections per internal model
      INTEGER      PPXREF_SECTIONS
        PARAMETER (PPXREF_SECTIONS =NSECTP-55)
! Max. number of non-null records in ppxref file (>1200)
      INTEGER      NUM_DIAG_MAX
        PARAMETER (NUM_DIAG_MAX    =NDIAGP)
! Max. number of user-defined ppxref records allowed
      INTEGER      NUM_USR_DIAG_MAX
        PARAMETER (NUM_USR_DIAG_MAX=450)
! No. of ppxref records read into PPXI,PPXC (for dyn. allocation)
      INTEGER      ppxRecs
! Global arrays:
! ppxref look-up arrays
      !INTEGER   PPXI(ppxRecs,PPXREF_CODELEN)
      !CHARACTER PPXC(NUM_DIAG_MAX,PPXREF_CHARLEN)
! Arrays for temporary storage of user-ppxref records -
!   used to transfer these records from STASH_PROC into U_MODEL
      INTEGER   PPXI_U(NUM_USR_DIAG_MAX,PPXREF_CODELEN)
      CHARACTER PPXC_U(NUM_USR_DIAG_MAX,PPXREF_CHARLEN)
! Array of flags to indicate origin of ppxref record
! 'P' for ppxref file; 'U' for user-stash master file
      CHARACTER OriginFlag(NUM_DIAG_MAX)
! Array of indices to identify which ppxref record corresponds to
!   any given row of PPXI, PPXC
      INTEGER   RowIndex(NUM_DIAG_MAX)
! Pointer array for PPXI, PPXC arrays
      !INTEGER PPXPTR                                                    &
      !& (N_INTERNAL_MODEL    ,0:PPXREF_SECTIONS ,PPXREF_ITEMS)

! Common block:
      COMMON/PPX_INT/ RowIndex,PPXI_U
      COMMON/PPX_CHA/ OriginFlag,PPXC_U
! - End ppxlook -----------------------
#else
#include "ppxlook.h"
#endif

#include "decomptp.h"
#include "decompdb.h"
! Data kind parameters
#include "c_kinds.h"
!
!L
!L  DYNAMIC ALLOCATION OF SUPER ARRAYS:
!L
!L       Main D1 data array
#include "typspd1.h"
#if defined(NECSX6)
!dir$ cache_align spd1
#endif
!L
!L       STASH related arrays
#include "typspst.h"
!L
!L       Dump headers and lookups
#include "typspdua.h"
!L
!L       Pointers (addresses) of model variables and constants
#include "typsppta.h"
!L Maximum sizes of fields limited by User Interface
!L CMAXSIZE now included earlier in routine
!L
!L       Model derived constants arrays
#include "typspcoa.h"
!L
!L       Generation of output interface fields
#include "typspina.h"
!L
!L       Updating of model from ancillary files
#include "typspana.h"
!L
!L       Boundary updating for Limited Area Models
#include "typspbo.h"
#include "typspboa.h"
!L
!L       Coupled model arrays (atmosphere-ocean)
#include "typspcpl.h"
!L
!L       Control for temporal filtering.
#include "ctfilt.h"
#if defined(ATMOS)
!  Sizes for allocation of AC scheme arrays
#include "csizeobs.h"
      INTEGER :: obs_flag_len,obs_len
      INTEGER, ALLOCATABLE, DIMENSION(:) :: OBS_FLAG
      REAL,    ALLOCATABLE, DIMENSION(:) :: OBS
#endif
!L       Control for boundary updating.
#include "cbound.h"
!L
!
!  Local variables
!
      INTEGER internal_model    ! Work - Internal model identifier
      INTEGER internal_model_prev!Work - Previous internal model ident
      INTEGER submodel          ! Work - Submodel id for dump partition
      INTEGER submodel_prev     ! Work - Previous submodel dump id
      INTEGER NGROUP            ! Work - Number of steps in "group"
      INTEGER MEANLEV           ! Work - Mean level indicator
      INTEGER IABORT            ! Work - Internal return code
      INTEGER I_STEP            ! Work - Loop counter over timesteps
      INTEGER G_THETA_FIELD_SIZE                                        &
                                   ! Sizes for MPP dynamic allocation
     &       ,G_IMTJMT          ! in A-O coupling routines
!
! River routing
      INTEGER G_RIVER_FIELD_SIZE   ! Sizes for MPP dynamic allocation
!
      LOGICAL LEXITNOW          ! Work - Immediate exit indicator
      CHARACTER*14 PPNAME       ! Work - Dummy PP filename
      INTEGER NFT           ! Unit no. for standard STASHmaster files
      INTEGER NFTU          ! Do. user STASH files (for GET_FILE)
      INTEGER RowNumber     ! Row no. counter for PPXI, PPXC arrays
      INTEGER I,J,K         ! Loop counters
      INTEGER CO2_DIMA,                                                 &
                                   ! CO2 array dimensions
     &        CO2_DIMO,                                                 &
     &        CO2_DIMO2
      INTEGER DMS_DIMA,                                                 &
                                   ! DMS array dimensions
     &        DMS_DIMO,                                                 &
     &        DMS_DIMO2
      Integer info   ! Return code from GCom routines

      integer len_runid       !  No of chars in RUNID

      character*4 runtype_char  !  Run Type (ie. NRUN, CRUN)

! 3-D fields of species to be passed down to radiation
      INTEGER, PARAMETER :: ngrgas = 8
      INTEGER, SAVE :: grgas_addr(ngrgas)

      ! Lookup tables of IAU increment file:
      INTEGER IAU_lookup(IAU_Len1Lookup, IAU_Len2Lookup)

      ! Arrays for storing IAU/TDF fields, one for packed fields, one
      ! for unpacked fields:
      REAL(KIND=real32) D1_IAU_k4(D1_IAU_k4_len)
      REAL         D1_TFilt (D1_TFilt_len )
      LOGICAL PUT_STEP ! True when we're going to do a put
                       ! for OASIS4 purposes, false if we use dummy data.
      LOGICAL GET_STEP ! True if we need to get data from OASIS4, false
                       ! if we perform a get but discard results.
      LOGICAL CPL_UPDATE_STEP ! True when we need to make sure
                              ! the D1 prognostic data is updated
                              ! in preparation for dump creation or
                              ! coupling actions.  
      INTEGER OASIS_COUPLE_TS ! OASIS_COUPLE_FREQ (hours) converted to
                              ! a number of timesteps.     
 
      INTEGER stepno
      INTEGER oasis_error  ! Just a flag
      LOGICAL L_OASIS3, L_OASIS4

#include "cprintst.h"
! Error reporting
      INTEGER       ICODE       ! =0 normal exit; >0 error exit
      CHARACTER*256 Cmessage    ! Error message
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='U_MODEL')

! Needed for non-coupled ACCESS initialisation
#include "c_0_dg_c.h"

#if defined(FLUME)
! FLUME-STASH  PPX arrays in SharedPPX module
      ALLOCATE(PPXI(ppxRecs,PPXREF_CODELEN))
      ALLOCATE(PPXC(NUM_DIAG_MAX,PPXREF_CHARLEN))
      ALLOCATE(PPXPTR                                                    &
     & (N_INTERNAL_MODEL    ,0:PPXREF_SECTIONS ,PPXREF_ITEMS))
#endif

      ICODE=0
      CMESSAGE=''

!L----------------------------------------------------------------------
!L 0. Start Timer call for U_MODEL (NB: not conditional on LTIMER)
!L
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('U_MODEL ',3)
      END IF

! DEPENDS ON: lbc_coup_setup
      Call lbc_coup_setup(LTIMER)

      ICODE=0
#if defined(T3E)
!
!--find the start address of spd1
      i=loc(spd1(ixd1(2)))
!--find the offset to the nearest SCACHE line boundary upwards
      j=((i+63)/64)*64-i
!--compute the offset to be added to the index values
      j=j/8
!
!--add this offset on to the current addresses
      do k=1, ixd1_len
        ixd1(k)=ixd1(k)+j
      end do
!
      if (mype == 0) then
      IF (PrintStatus  >=  PrStatus_Diag) THEN
        WRITE(0,*) 'Memory address of submodel starts:'
        write(0,'(4z17)') (loc(spd1(ixd1(i))), i=1,4)
      END IF
      end if
#endif

!  Routine GETPPX_PART reads those ppxref records which correspond to
!  entries in the stash list into the ppx look-up arrays PPXI, PPXC.
!  It also sets the ppx pointer array PPXPTR. The lengths of PPXI, PPXC
!  have been dynamically allocated to the value of ppxRecs.

! Initialise row number in PPXI, PPXC arrays
      RowNumber = 1

! Initialise lookup and pointer array
      DO I=1,ppxRecs
        DO J=1,PPXREF_CODELEN
          PPXI(I,J)=0
        END DO
        DO J=1,PPXREF_CHARLEN
          PPXC(I,J) = ' '
        END DO
      END DO
      DO I = 1,N_INTERNAL_MODEL
        DO J   = 0,PPXREF_SECTIONS
          DO K = 1,PPXREF_ITEMS
            PPXPTR(I,J,K)=0
          END DO
        END DO
      END DO

! Read in STASHmaster records
      IF (INTERNAL_MODEL_INDEX(A_IM) >  0) THEN
! DEPENDS ON: getppx_part
      CALL GETPPX_PART(NFT,NFTU,'STASHmaster_A',A_IM,RowNumber,         &
#include "argppx.h"
     &                        ICODE,CMESSAGE)
! DEPENDS ON: ereport
      IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
      END IF

!L----------------------------------------------------------------------
!L 1. General initialisation of control and physical data blocks
!L

      ICODE=0

! DEPENDS ON: initial
      CALL INITIAL(                                                     &
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
     &     isSTASH)  ! FLUME-STASH
! DEPENDS ON: ereport
      IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)

! DEPENDS ON: icopenoutput
      CALL ICOPENOUTPUT(runtype_char)

!  Allocate AC scheme arrays using sizes from AC_INIT
      IF (L_AC) THEN
        obs_flag_len = A_MAX_NO_OBS
        ALLOCATE (OBS_FLAG(obs_flag_len))
! 2048 gives enough space in WORK in RDOBS2 when no or very few obs.
        obs_len = A_MAX_OBS_SIZE+2048
        ALLOCATE (OBS(obs_len))
      write(6,*)'U_MODEL - OBS arrays allocated with sizes ',           &
     & A_MAX_NO_OBS,A_MAX_OBS_SIZE+2048
      ELSE
        A_MAX_NO_OBS = 1
        A_MAX_OBS_SIZE = 1
        obs_flag_len = A_MAX_NO_OBS
        ALLOCATE (OBS_FLAG(obs_flag_len))
        obs_len = A_MAX_OBS_SIZE
        ALLOCATE (OBS(obs_len))
      write(6,*)'U_MODEL - OBS arrays allocated with length 1'
      END IF

!L----------------------------------------------------------------------
!L 2. Check for nothing-to-do
!L
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EXITCHEK',3)
! DEPENDS ON: exitchek
      CALL EXITCHEK( internal_model, LEXITNOW)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EXITCHEK',4)
      IF (LEXITNOW) GO TO 9999

#if defined(FLUME)
      ! FLUME-STASH isSTASH defaults to FALSE in non-Flume run
      IF (isSTASH) then
      ! Call flumeInterface to take copies of all the data structures
      !  which will bw needed when calling stash
! DEPENDS ON: flumeInterface
         CALL flumeInterface(  &
#include "artd1.h"
#include "artduma.h"
#include "artsts.h"
#include "argppx.h"
     &     ICODE,CMESSAGE)
         RETURN  ! For STASH process do not enter timestep loop
      endif 
#endif

!L----------------------------------------------------------------------
!L 3. Start group of timesteps

       ! This is the best place to initialise the OASIS3 
       ! grid details and transient variables because by this
       ! stage in the code, we should have read the start data and we'll
       ! have the things we need such as land-sea masks available
#if defined(ACCESS)
       IF (L_OASIS .and. l_auscom) THEN
#else
       IF (L_OASIS) THEN
#endif

          PUT_STEP=.FALSE.
          GET_STEP=.FALSE.
          CPL_UPDATE_STEP=.FALSE. 

          ! Work out the frequency of coupling timesteps
          OASIS_COUPLE_TS = STEPS_PER_PERIODim(1) *                        &
     &                  (OASIS_COUPLE_FREQ*3600)/SECS_PER_PERIODim(1)

#if defined(ACCESS)
! gol124: auscom coupling
          WRITE(6,*)'U_MODEL: STEPS_PER_PERIODim=',STEPS_PER_PERIODim(1)
          WRITE(6,*)'U_MODEL: OASIS_COUPLE_FREQ=',OASIS_COUPLE_FREQ
          WRITE(6,*)'U_MODEL: SECS_PER_PERIODim=',SECS_PER_PERIODim(1)
          WRITE(6,*)'U_MODEL: OASIS_COUPLE_TS=',OASIS_COUPLE_TS
#endif

          ! Save our model start date for comparison with PRISM start 
          ! date. The easiest thing is to copy it to an array which we 
          ! use in a special module as follows:
          ! um_start_date(1:6) = MODEL_BASIS_TIME(1:6)
          ! But for now this is deactivated until we sort out similar
          ! checks in other components. 


#if defined(OASIS3)       
          L_OASIS3=.TRUE.
          ! Set timestep for use in OASIS4 calls to be
          ! real version of the standard atmos timestep
          PRISM_Timestep = 1.0*SECS_PER_STEPim(atmos_sm)
#if defined(ACCESS)
          IF (l_auscom) THEN
             WRITE(6,*)'U_MODEL: Calling GRID64_OASIS3'

! DEPENDS ON: OASIS3_GRID64
             CALL OASIS3_GRID64(OASIS_comp_id,MODEL_BASIS_TIME,ICODE)

             WRITE(6,*) 'U_MODEL: PRISM_TIMESTEP=',PRISM_Timestep

             WRITE(6,*) 'U_MODEL: Completed GRID64_OASIS3',ICODE
          else
             WRITE(6,*)'U_MODEL: Skipping GRID64_OASIS3'
          endif
#else
          WRITE(6,*)'U_MODEL: Calling GRID64_OASIS3'

! DEPENDS ON: OASIS3_GRID64
          CALL OASIS3_GRID64(OASIS_comp_id,MODEL_BASIS_TIME,ICODE)

          ! Set timestep for use in OASIS4 calls to be
          ! real version of the standard atmos timestep
          PRISM_Timestep = 1.0*SECS_PER_STEPim(atmos_sm)

          WRITE(6,*) 'U_MODEL: PRISM_TIMESTEP=',PRISM_Timestep

          WRITE(6,*) 'U_MODEL: Completed GRID64_OASIS3',ICODE
#endif

#endif

#if defined(OASIS4)
          L_OASIS4=.TRUE.
          WRITE(6,*)'U_MODEL: Calling GRID64_OASIS4'

! DEPENDS ON: OASIS4_GRID64
          CALL OASIS4_GRID64(OASIS_comp_id,MODEL_BASIS_TIME,ICODE)

          ! Work out the time window for the zeroth timestep
          CALL set_oasis_window32(icode)

          WRITE(6,*)'U_MODEL: Completed GRID64_OASIS4',ICODE
#endif
#if defined(ACCESS)
       ELSE
          ! Should come from c_0_dg_c.h 
          access_tfs = 271.35
          ocn_sss=.false.
          ! initilise um variables by "input_atm.nml" file
! DEPENDS ON: OASIS3_UMVARS_INIT
          CALL OASIS3_UMVARS_INIT()

#endif

       END IF

#if defined(ACCESS)
       ! Extra initialisation for running coupled model stand-alone
       if (.not. L_AUSCOM ) then
          ocn_sss = .false.
          access_tfs = tfs
       end if
#endif

!L----------------------------------------------------------------------
!L 3. Start group of timesteps
!L
   1  CONTINUE
!L----------------------------------------------------------------------
!L 3.1. Start main timestep loop
!L
!L 3.1.1 Increment model time ..
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INCRTIME',3)
! DEPENDS ON: incrtime
      CALL INCRTIME (                                                   &
#include "artduma.h"
     &       internal_model,ICODE,CMESSAGE)
! DEPENDS ON: ereport
      IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INCRTIME',4)

! Keep tabs on PRISM PUT/GET timesteps.
! At the moment we say a put and get timestep are one and the same.
! We're also hard coding to do operations every 48th timestep i.e.
! assuming a half hour atmos TS!. 
#if defined(ACCESS)
      IF (L_OASIS .and. l_auscom) THEN
#else
      IF (L_OASIS) THEN
#endif

         ! Is this a genuine exchange timestep.
#if defined(ACCESS)
! gol124: auscom coupling
         PUT_STEP = (MOD(STEPim(atmos_im),OASIS_COUPLE_TS).EQ.0)
#else
         PUT_STEP = (MOD(STEPim(atmos_im),OASIS_COUPLE_TS).EQ.1)
#endif
         GET_STEP = (MOD(STEPim(atmos_im),OASIS_COUPLE_TS).EQ.1)
         CPL_UPDATE_STEP=(MOD(STEPim(atmos_im),OASIS_COUPLE_TS).EQ.0)
         stepno = STEPim(atmos_im)
         if (PUT_STEP ) write(6,*) "time_step to couple PUT:",stepno
         if (GET_STEP ) write(6,*) "time_step to couple GET:",stepno
         if (CPL_UPDATE_STEP ) write(6,*) "time_step to couple UPDATE:",stepno


         ! Perform coupling exchanges relating to TS-1
#if defined(ACCESS)
         ! avoid unnecessary calls to prism and scatter
         IF (L_OASIS3 .and. get_step) THEN
#else
         IF (L_OASIS3) THEN
#endif
      
! DEPENDS ON: timer
            IF (LTIMER) CALL TIMER('oasis3get',3)
    
! DEPENDS ON: oasis3_geto2a
            CALL oasis3_geto2a(                                         &
#include "artd1.h"
#include "arg_atm_fields.h"
     &               GET_STEP,cmessage)

! DEPENDS ON: timer
            IF (LTIMER) CALL TIMER('oasis3get',4)

#if !defined(ACCESS)
! DEPENDS ON: timer
            IF (LTIMER) CALL TIMER('oasis3put',3)

! DEPENDS ON: oasis3_puta2o
            CALL oasis3_puta2o(                                         &
#include "artd1.h"
#include "arg_atm_fields.h"
     &              PUT_STEP,cmessage)

! DEPENDS ON: timer
            IF (LTIMER) CALL TIMER('oasis3put',4)

         ! Advance date ready for next timestep (if there is one)
! DEPENDS ON: OASIS3_ADVANCE_DATE64
            CALL OASIS3_ADVANCE_DATE64(oasis_error)
#endif

         END IF  ! L_OASIS3=true
      END IF

!L 3.1.2 .. set timestep control switches
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('SETTSCTL',3)
! DEPENDS ON: settsctl
      CALL SETTSCTL (                                                   &
#include "artduma.h"
#include "artsts.h"
#include "artinfa.h"
     &         internal_model,.FALSE.,MEANLEV,ICODE,CMESSAGE,&
     &         isSTASH)       ! FLUME-STASH
! DEPENDS ON: ereport
      IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('SETTSCTL',4)
!L 3.1.3 If PPfile initialisation time call PP control routine
!L          for instantaneous data (MEANLEV=0)
! FLUME-STASH - Flume_run is set to false in non-flume runs
      IF (LPP .and. .not.Flume_run) THEN
        write(6,*)  "PPCTL_REINIT timestep", STEPim(atmos_im)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('PPCTL_REINIT',3)
! DEPENDS ON: ppctl_reinit
        CALL PPCTL_REINIT(                                              &
#include "artduma.h"
#include "artinfa.h"
     &         internal_model,PPNAME,ICODE,CMESSAGE)
! DEPENDS ON: ereport
      IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('PPCTL_REINIT',4)
      END IF  ! FLUME-STASH

!L       Integrate atmosphere or ocean by 1 timestep
          IF (internal_model == atmos_im) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('ATM_STEP',3)
! Synchronize before the timestep starts
      CALL GC_GSYNC(nproc,info)

! River routing
      IF (L_RIVERS) THEN
! Get 'global' atmos horizontal domain sizes from database
! in DECOMPDB to set dynamic allocation in ATM_STEP for River routing
! on PE 0                                                   .
        G_THETA_FIELD_SIZE=                                             &
     &  decomp_db_glsize(1,fld_type_p,decomp_standard_atmos) *          &
     &  decomp_db_glsize(2,fld_type_p,decomp_standard_atmos)
        G_RIVER_FIELD_SIZE=                                             &
     &  decomp_db_glsize(1,fld_type_r,decomp_standard_atmos) *          &
     &  decomp_db_glsize(2,fld_type_r,decomp_standard_atmos)
      ELSE
        G_THETA_FIELD_SIZE=1
        G_RIVER_FIELD_SIZE=1
      END IF
!

! DEPENDS ON: atm_step
         CALL ATM_STEP (                                                &
#include "artd1.h"
#include "artduma.h"
#include "artlndm.h"
#include "arg_atm_fields.h"
#include "artbnd.h"
#include "artsts.h"
#include "argppx.h"
! River routing
#include "artatcpl.h"
     & G_THETA_FIELD_SIZE,                                              &
     & G_RIVER_FIELD_SIZE,                                              &
     & OBS_FLAG,OBS,obs_flag_len,obs_len,                               &
!
     &      ngrgas,grgas_addr,                                          &
     &                     IAU_lookup,                                  &
                                       ! inout
     &                     D1_IAU_k4,                                   &
                                       ! inout
     &                     D1_TFilt,   &! inout
          !jhan:pass thru end step to CABLE
     &    TARGET_END_STEPim(a_im) )

! DEPENDS ON: timer
      IF (LTIMER) CALL TIMER('ATM_STEP',4)

!L Generate Atmosphere lateral boundary values

      If (linterface) Then

! DEPENDS ON: timer
                            If (LTIMER) Call TIMER ('GEN_INTF',3)

! DEPENDS ON: gen_intf
        Call GEN_INTF (                                                 &
#include "artd1.h"
#include "artduma.h"
#include "artsts.h"
#include "artptra.h"
#include "artinfa.h"
#include "argppx.h"
     &       submodel, ICode, CMessage)

! DEPENDS ON: timer
                            If (LTIMER) Call TIMER ('GEN_INTF',4)

        If (ICode /= 0) Then
! DEPENDS ON: ereport
          CALL Ereport(RoutineName, ICode ,Cmessage)
        End If

      End If  ! LINTERFACE

      IF (L_ukca) THEN
        WRITE(6,*) ' Before call to UKCA_MAIN1'
! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('UKCA_MAIN1',5)
! DEPENDS ON: ukca_main1
        CALL UKCA_MAIN1(                                                &
#include "artd1.h"
#include "artduma.h"
#include "artsts.h"
#include "artlndm.h"
#include "argppx.h"
#include "artptra.h"
     &   I)                  ! dummy to terminate call
! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('UKCA_MAIN1',6)
      END IF

#if defined(A25_1A)
        WRITE(6,*) ' REDIST_STOCHEM called'
! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('REDIST_STOCHEM',5)
! DEPENDS ON: redist_stochem
        CALL REDIST_STOCHEM(                                            &
#include "artd1.h"
#include "artduma.h"
#include "artsts.h"
     &   L_USE_STOCHEM_CH4,L_USE_STOCHEM_O3,                            &
     &   I)                  ! dummy to terminate call
! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('REDIST_STOCHEM',6)
#endif
         END IF        ! internal_model = atmos_im

         IF (L_OASIS) THEN

            ! Applicable to OASIS3 or OASIS4

#if defined(ACCESS)
            IF (CPL_UPDATE_STEP .and. l_auscom) THEN
#else
            IF (CPL_UPDATE_STEP) THEN
#endif
!---------------------------------------------------------------------
! Ensure atmos coupling data in prognostic areas of D1 are up to date.
! The logic here is that CPL_UPDATE_STEP will be TRUE on the timestep
! BEFORE coupling is due to take place. 
!
! Newly generated coupling data is intercepted after ATM_STEP and 
! copied to  D1 prognostics.
!
! On the next timestep i.e. a coupling timestep, the prognostic
! contents will be sent to the other components in the coupling process.
!
! If there is no subsequent timestep (i.e. if this is the last model
! timestep then the D1 contents will be written to the dump
! ready for any future restart). There is no "end-of-model"
! coupling exchange. 
!---------------------------------------------------------------------

               If (ltimer) Call timer('updatecpl',3)

! DEPENDS ON: oasis_updatecpl
               CALL OASIS_UPDATECPL(                                           &
#include "artd1.h"
#include "arg_atm_fields.h"
     &               cmessage)

               If (ltimer) Call timer('updatecpl',4)

            END IF

#if defined(ACCESS)
! gol124: auscom coupling
          IF (l_oasis3) THEN
! DEPENDS ON: timer
            IF (LTIMER) CALL TIMER('oasis3put',3)

! DEPENDS ON: oasis3_puta2o
            CALL oasis3_puta2o(                                         &
#include "artd1.h"
#include "arg_atm_fields.h"
     &              PUT_STEP,cmessage)

! DEPENDS ON: timer
            IF (LTIMER) CALL TIMER('oasis3put',4)

         ! Advance date ready for next timestep (if there is one)
! DEPENDS ON: OASIS3_ADVANCE_DATE64
            CALL OASIS3_ADVANCE_DATE64(oasis_error)
           end if ! l_oasis3
#endif
         END IF  ! l_oasis


!L 3.1.4 If dump time, call dump control routine
          IF (LDUMP) THEN
            IF (LTIMER) THEN
! DEPENDS ON: timer
              CALL TIMER('DUMPCTL',5)
! DEPENDS ON: timer
              CALL TIMER('DUMPCTL ',3)
            END IF
! DEPENDS ON: dumpctl
            CALL DUMPCTL (                                              &
#include "artd1.h"
#include "artduma.h"
#include "artlndm.h"
#include "artptra.h"
#include "artsts.h"
#include "argppx.h"
     &          submodel,MEANLEV,.false.,'           ',0,               &
     &          ICODE,CMESSAGE)
! DEPENDS ON: ereport
            IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)

            IF (LTIMER) THEN
! DEPENDS ON: timer
              CALL TIMER('DUMPCTL',4)
! DEPENDS ON: timer
              CALL TIMER('DUMPCTL ',6)
            END IF
!L 3.1.4.1 Update interim history file unless means are to follow
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('TEMPHIST',3)
            IF (.NOT.LMEAN) THEN
              IF (.NOT. (N_SUBMODEL_PARTITION  >  1 .AND. submodel      &
     &          /=  SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION) .AND. &
     &          steps_per_periodim(submodel)  /=                        &
     &          dumpfreqim(submodel) )) THEN
! DEPENDS ON: set_history_values
                CALL SET_HISTORY_VALUES
                CALL GC_SSYNC(nproc,info)
                IF (MYPE  ==  0) THEN
! DEPENDS ON: temphist
                  CALL TEMPHIST(PHIST_UNIT,ICODE,CMESSAGE)
! DEPENDS ON: ereport
            IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
                END IF
              END IF 
            END IF
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('TEMPHIST',4)
!L 3.1.4.2 If atmosphere timestep recalculate prognostic data and
!L         wrap-around fields using rounded off values
            IF (submodel == atmos_sm) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('RESETATM',3)
! DEPENDS ON: resetatm
       CALL RESETATM(                                                   &
#include "artd1.h"
#include "artptra.h"
     &                                ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('RESETATM',4)
              IF (ICODE >  0) GO TO 9999
            END IF
          END IF ! L_DUMP

!L 3.1.5 If printed output time, call print control routine
          IF (LPRINT) THEN

          IF (PrintStatus >= PrStatus_Oper) THEN
         WRITE(6,*) RoutineName,':Warning, Printing of climate global ' &
     &             ,'and zonal diagnostics no longer supported'
          END IF  ! PrintStatus test

          END IF
!L 3.1.6 If interface generation time, generate interface fields
          IF ((internal_model == ocean_im) .or.                         &
     &               (internal_model == wave_im)) THEN

            IF (LINTERFACE) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('GEN_INTF',3)
! DEPENDS ON: gen_intf
              CALL GEN_INTF (                                           &
#include "artd1.h"
#include "artduma.h"
#include "artsts.h"
#include "artptra.h"
#include "artinfa.h"
#include "argppx.h"
     &              submodel,ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('GEN_INTF',4)
! DEPENDS ON: ereport
            IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
            END IF

          END IF
!L 3.1.6.1 Release job to process output created so far, if selected
          IF (LJOBRELEASE) THEN
! DEPENDS ON: flush_all_pp
            CALL FLUSH_ALL_PP()
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('JOBCTL  ',3)
! DEPENDS ON: jobctl
            CALL JOBCTL(internal_model,ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('JOBCTL  ',4)
! DEPENDS ON: ereport
          IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
          END IF
!L 3.1.7 If partial sum/mean creation time, call means control routine
!L       (calls mean PPfield and diagnostic print routines internally)
          IF (LMEAN) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('MEANCTL ',3)
! DEPENDS ON: meanctl
            CALL MEANCTL (                                              &
#include "artd1.h"
#include "artduma.h"
#include "artptra.h"
#include "artsts.h"
#include "artlndm.h"
#include "artinfa.h"
#include "argppx.h"
     &                  submodel,MEANLEV,ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('MEANCTL ',4)
            IF (ICODE >  0) THEN
! DEPENDS ON: del_hist
              CALL DEL_HIST(PHIST_UNIT)
      WRITE(6,*)'U_MODEL: interim history file deleted due to failu     &
     &re writing partial sum files'
! DEPENDS ON: ereport
              CALL Ereport(RoutineName,ICODE,Cmessage)
            END IF
!L 3.1.7.1 On successful completion, update interim history file
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('TEMPHIST',3)

            IF (.NOT. (N_SUBMODEL_PARTITION  >  1 .AND. submodel        &
     &          /=  SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION) .AND. &
     &          steps_per_periodim(submodel)  /=                        &
     &          dumpfreqim(submodel) )) THEN

! DEPENDS ON: set_history_values
              CALL SET_HISTORY_VALUES
              CALL GC_SSYNC(nproc,info)
              IF (MYPE  ==  0) THEN
! DEPENDS ON: temphist
                CALL TEMPHIST(PHIST_UNIT,ICODE,CMESSAGE)
              END IF
! DEPENDS ON: ereport
            IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
            END IF
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('TEMPHIST',4)
          END IF
!L 3.1.8 Update temporary history file if at a 'safe' restart point
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('TEMPHIST',3)
          IF (LHISTORY) THEN
!          In coupled model do not update history file until both
!          submodels have reached the safe restart point

            IF (.NOT. (N_SUBMODEL_PARTITION  >  1 .AND. submodel        &
     &          /=  SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION) .AND. &
     &          steps_per_periodim(submodel)  /=                        &
     &          dumpfreqim(submodel) )) THEN
! DEPENDS ON: set_history_values
              CALL SET_HISTORY_VALUES
            CALL GC_SSYNC(nproc,info)
            IF (MYPE  ==  0) THEN
! DEPENDS ON: temphist
            CALL TEMPHIST(THIST_UNIT,ICODE,CMESSAGE)
            END IF
          END IF
        END IF
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('TEMPHIST',4)
! DEPENDS ON: ereport
      IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
!L 3.1.9 If exit check time, check for immediate exit
          IF (LEXIT) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EXITCHEK',3)
! DEPENDS ON: exitchek
            CALL EXITCHEK(internal_model, LEXITNOW)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EXITCHEK',4)
            IF (LEXITNOW) THEN
              IF (.NOT.LDUMP) THEN

                IF (.NOT. (N_SUBMODEL_PARTITION  >  1 .AND. submodel    &
     &           /=  SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION).AND. &
     &          steps_per_periodim(submodel)  /=                        &
     &          dumpfreqim(submodel) )) THEN

! DEPENDS ON: set_history_values
                CALL SET_HISTORY_VALUES
                CALL GC_SSYNC(nproc,info)
                IF (MYPE  ==  0) THEN
! DEPENDS ON: temphist
                CALL TEMPHIST(PHIST_UNIT,ICODE,CMESSAGE)
                END IF
!               Exit model so no need to set ppflush
         END IF
              END IF
              GO TO 9999
      END IF
          END IF
!L 3.1.10 Update ancillary fields if necessary
          IF (LANCILLARY) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('UP_ANCIL',3)
! DEPENDS ON: up_ancil
         CALL UP_ANCIL (                                                &
#include "artd1.h"
#include "artduma.h"
#include "artptra.h"
#include "artanc.h"
     &                   submodel,                                      &
#include "argppx.h"
     &                   ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('UP_ANCIL',4)
! DEPENDS ON: ereport
      IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
          END IF
!L 3.1.11 Update boundary fields if necessary
          IF (LBOUNDARY) THEN
             ! DEPENDS ON: lbc_coup_update
             call lbc_coup_update( &
#include "artd1.h"
#include "artduma.h"
#include "artptra.h"
#include "artbnd.h"
#include "argppx.h"
                  submodel,ICODE,CMESSAGE)
          END IF
!L
!L      End main timestep loop
!L----------------------------------------------------------------------

#if defined(FLUME)
      IF (Flume_run) THEN
        ! FLUME-STASH  Mark the end of the timestep
        CALL flumeEOTS()  ! sends end of timestep tag
      END IF
#endif

      GO TO 1
!
 9999 CONTINUE

#if defined(FLUME)
      IF (Flume_run) THEN
        ! FLUME-STASH  Mark the end of processing
        CALL flumeEOM() ! sends end of processing tag
      END IF
#endif

!L----------------------------------------------------------------------
!L 4. Exit processing: Output error messages and perform tidy-up
!L
#if defined(ACCESS)
      IF (L_OASIS .and. l_auscom) THEN
#else
      IF (L_OASIS) THEN
#endif
! DEPENDS ON: oasis_tidy
         CALL oasis_tidy(                                        &
#include "artd1.h"
#include "artsts.h"
#include "artduma.h"
#include "artptra.h"
     &        icode,cmessage)
      ENDIF


! DEPENDS ON: iccloseoutput
      CALL ICCLOSEOUTPUT()

!L 4.1 Exit processing: If abnormal completion, output error message
      IABORT = ICODE
      IF (ICODE /= 0) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EREPORT ',3)
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ICODE,Cmessage)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EREPORT ',4)
      END IF
!L 4.2 Exit processing: Perform tidy-up
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EXITPROC',3)
! DEPENDS ON: exitproc
      CALL EXITPROC(ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EXITPROC',4)
!L 4.3 Exit processing: If error in exit processing, output error mess
      IF (ICODE /= 0) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EREPORT ',3)
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ICODE,Cmessage)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EREPORT ',4)
      END IF

#if defined(FLUME)
! FLUME-STASH  PPX arrays in SharedPPX module
      DEALLOCATE(PPXI)
      DEALLOCATE(PPXC)
      DEALLOCATE(PPXPTR)  
#endif

!L----------------------------------------------------------------------
!L 5. Complete Timer call and return
!L
      ICODE=IABORT
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('U_MODEL ',4)
      END IF
      RETURN
      END SUBROUTINE U_MODEL
#endif
