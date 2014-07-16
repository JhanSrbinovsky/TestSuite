
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Initialise model for submodel and internal model coupling
!
! Subroutine Interface:
      SUBROUTINE UM_Submodel_Init(ErrorStatus)

      IMPLICIT NONE
!
! Description:
!   UM_Submodel_Init initialises the model with information specifying
!   internal model and submodel partitions for the run, which is
!   required for control of coupling when more than one internal model
!   is present.
!
! Method:
!   The routine reads information from the user interface, providing
!   lists of internal models and their associated submodel data
!   partitions. This is required in both the reconfiguration and the
!   model as a prior step to calculating addressing in STASH_PROC.
!
! Current Code Owner: R. Rawlins
! History:
! Version   Date     Comment
! -------   ----     -------
! 3.5    07/04/95   Original code. R. Rawlins.
!LL 4.3-4.4   16/09/97 D1 addressing change and subsequent correction
!LL                    S.D.Mullerworth
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered:
! System Task:
! Declarations:
!
!
! Global variables (*CALLed COMDECKs etc...):
! CSUBMODL start
!
! Description:
!    Describes the number and identity of submodels available
!    within the system, and those included in the current
!    experiment.  Parameters set by the User Interface give
!    the relevant array sizes; other submodel configuration
!    information is either read from NAMELIST input, or
!    derived from dump header information.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! pre 3.0           Original code. T. Johns
! 3.5    07/04/95   Expansion for stage 1 of submodel project, allowing
!                   flexible specification of internal models within
!                   submodel partitions. R. Rawlins
!
! Declarations:
!
!  1. Internal model and submodel dump partition identifiers - fixed
!     for all experiments.
! CSMID start
!
! Description:
!    Hold parameters defining internal model identifiers and submodel
!    data partition (ie main D1 data array and consequent dump), both
!    short and long form.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! pre 3.0           Original code. T. Johns
! 3.3    26/10/93   M. Carter. Part of an extensive mod that:
!                    1.Removes the limit on primary STASH item numbers.
!                    2.Removes the assumption that (section,item)
!                      defines the sub-model.
!                    3.Thus allows for user-prognostics.
!                    Add index to submodel home dump.
! 3.5    13/03/95   Expansion for stage 1 of submodel project, allowing
!                   flexible specification of internal models within
!                   submodel partitions. R. Rawlins
! 6.0    02/07/03   Add X_IM and X_SM for small exec.      E.Leung
!
! Declarations:
!
!   Hold parameters defining internal model identifiers and submodel
!   data partition (ie main D1 data array and consequent dump), both
!   short and long form
      ! Internal models
      INTEGER,PARAMETER:: A_IM      = 1 ! Atmosphere internal model
      INTEGER,PARAMETER:: ATMOS_IM  = 1 ! Atmosphere internal model
      INTEGER,PARAMETER:: O_IM      = 2 ! Ocean internal model
      INTEGER,PARAMETER:: OCEAN_IM  = 2 ! Ocean internalmodel
      INTEGER,PARAMETER:: S_IM      = 3 ! Slab internal model
      INTEGER,PARAMETER:: SLAB_IM   = 3 ! Slab internal model
      INTEGER,PARAMETER:: W_IM      = 4 ! Wave internal model
      INTEGER,PARAMETER:: WAVE_IM   = 4 ! Wave internal model
      INTEGER,PARAMETER:: I_IM      = 5 ! Sea=ice internal model
      INTEGER,PARAMETER:: SEAICE_IM = 5 ! Sea=ice internal model
      ! New dynamics (Charney-Phillips grid)
      INTEGER,PARAMETER:: N_IM      = 6 ! ND internal model
      INTEGER,PARAMETER:: NATMOS_IM = 6 ! ND internal model
      ! Small Executables
      INTEGER,PARAMETER:: X_IM      = 7 ! SX indicator

      ! Submodels
      INTEGER,PARAMETER:: A_SM      = 1 ! Atmosphere submodel
      INTEGER,PARAMETER:: ATMOS_SM  = 1 ! Atmosphere submodel
      INTEGER,PARAMETER:: O_SM      = 2 ! Ocean submodel
      INTEGER,PARAMETER:: OCEAN_SM  = 2 ! Ocean submodel
      INTEGER,PARAMETER:: W_SM      = 4 ! Wave submodel
      INTEGER,PARAMETER:: WAVE_SM   = 4 ! Wave submodel
      ! New dynamics (Charney-Phillips grid)
      INTEGER,PARAMETER:: N_SM      = 6 ! ND submodel
      INTEGER,PARAMETER:: NATMOS_SM = 6 ! ND submodel
      ! Small Executables
      INTEGER,PARAMETER:: X_SM      = 7 ! SX indicator

! CSMID end

!
!  2. Maximum internal model/submodel array sizes for this version.
!
! CSUBMAX start
!
! Description:
!    Describes the number and identity of submodels available
!    within the system, and those included in the current
!    experiment.  Parameters set by the User Interface give
!    the relevant array sizes; other submodel configuration
!    information is either read from NAMELIST input, or
!    derived from dump header information.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 3.5    13/07/95   Original code. D.M. Goddard
! 4.0     3/11/95   Reduce max internal model, submodel from 10 to 4
!                   to save space in model. At 4.0 the max no of
!                   supported models is 3, 1 slot is reserved for
!                   expansion. Rick Rawlins.
!  4.1  21/02/96  Wave model introduced as 4th sub-model.  RTHBarnes
!
! Declarations:
!
!
!  1. Maximum internal model/submodel array sizes for this version.
!
      ! Max no. of internal models
      INTEGER,PARAMETER:: N_INTERNAL_MODEL_MAX=4

      ! Max no. of submodel dump partitions
      INTEGER,PARAMETER:: N_SUBMODEL_PARTITION_MAX=4

      ! Max value of internal model id
      INTEGER,PARAMETER:: INTERNAL_ID_MAX=N_INTERNAL_MODEL_MAX

      ! Max value of submodel dump id
      INTEGER,PARAMETER:: SUBMODEL_ID_MAX=N_SUBMODEL_PARTITION_MAX

! CSUBMAX end
!
!  3. Lists of internal models and their submodel dump partitions -
!     initialised by the user interface - experiment specific.
      INTEGER :: N_INTERNAL_MODEL          ! No. of internal models
      INTEGER :: N_SUBMODEL_PARTITION      ! No. of submodel partitions

      ! Internal models
      INTEGER :: INTERNAL_MODEL_LIST(N_INTERNAL_MODEL_MAX)

      ! Submodel identifier for each internal model in list
      INTEGER :: SUBMODEL_FOR_IM    (N_INTERNAL_MODEL_MAX)

      ! Submodel number for each submodel id
      INTEGER :: SUBMODEL_FOR_SM(N_INTERNAL_MODEL_MAX)

      ! Namelist for information in 3.
      NAMELIST/NSUBMODL/N_INTERNAL_MODEL,N_SUBMODEL_PARTITION,          &
     &  INTERNAL_MODEL_LIST,SUBMODEL_FOR_IM

      ! 4. Lists calculated in model from user interface supplied arrays
      ! experiment specific.

      ! No of internal models in each submodel partition indexed by sm
      !  identifier
      INTEGER :: N_INTERNAL_FOR_SM(SUBMODEL_ID_MAX)

      ! List of  submodel partition identifiers
      INTEGER :: SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION_MAX)

      ! Submodel partition identifier indexed by internal model identifie
      INTEGER :: SUBMODEL_PARTITION_INDEX(INTERNAL_ID_MAX)

      ! Sequence number of internal model indexed by internal model
      ! identifier: required to map from id to STASH internal model
      ! sequence
      INTEGER :: INTERNAL_MODEL_INDEX(INTERNAL_ID_MAX)


      ! Last internal model within a submodel partition if .TRUE.,
      ! indexed by internal model id.
      LOGICAL :: LAST_IM_IN_SM(INTERNAL_ID_MAX)

      ! Common block for information in 3. and 4.
      COMMON/SUBMODL/N_INTERNAL_MODEL,N_SUBMODEL_PARTITION,             &
     &  INTERNAL_MODEL_LIST,SUBMODEL_FOR_IM,SUBMODEL_FOR_SM,            &
     &  N_INTERNAL_FOR_SM,SUBMODEL_PARTITION_LIST,                      &
     &  SUBMODEL_PARTITION_INDEX,                                       &
     &  INTERNAL_MODEL_INDEX,                                           &
     &  LAST_IM_IN_SM

!
!  5. Time information specifying coupling frequencies between internal
!     models and submodels, and multipliers, indexed by sequence of
!     internal models and submodels (ie left to right along node tree).
!     {Not required at this release}.
!
! Namelists for information in 5. {Not required at this release}
!
!
!  6. Lists of coupling nodes defining coupling frequencies between
!     internal models and between submodel partitions. (Not defined
!     yet at this release).
!CALL CNODE
!
!  7. Variables dealing with general coupling switches at the control
!     level. {These will require revision at the next release when
!     coupling between internal models is dealt with more generally.
!     Logicals below are set in routine SETGRCTL.}

      ! new internal model next group of timesteps if .true.
      LOGICAL :: new_im

      ! new submodel dump  next group of timesteps if .true.
      LOGICAL :: new_sm

      COMMON/CSUBMGRP/new_im,new_sm

      INTEGER SUBMODEL_IDENT
      COMMON/SUBMODID/SUBMODEL_IDENT
! CSUBMODL end
! An alternative common block required by TYPD1
! CALTSUBM
! TYPD1 needs access to N_SUBMODEL_PARTITION/_MAX in CSUBMODL. However,
! they are not always called in the same decks and in the right order.
! Therefore, copy the values to another file and include it from TYPD1

      INTEGER ALT_N_SUBMODEL_PARTITION

      INTEGER, PARAMETER :: ALT_N_SUBMODEL_PARTITION_MAX=4

      COMMON/CALTSUBM/ALT_N_SUBMODEL_PARTITION
! CALTSUBM end
! Subroutine arguments
!   Scalar arguments with intent(in):

!   Array  arguments with intent(in):

!   Scalar arguments with intent(InOut):

!   Array  arguments with intent(InOut):

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):

!   ErrorStatus
      INTEGER      ErrorStatus          ! Error flag (0 = OK)

! Local parameters:

! Local scalars:
      INTEGER                                                           &
     & s                                                                &
                       ! submodel loop
     &,i                                                                &
                       ! internal model loop
     &,sm                                                               &
                       ! submodel identifier
     &,im                                                               &
                       ! internal model identifier
     &,sm_prev                                                          &
                       ! previous submodel identifier
     &,im_prev         ! previous internal model identifier

! Local dynamic arrays:

! Function & Subroutine calls: None

!- End of header
!
! 1. Initialise lists before obtaining values for this experiment.
!
      do i=1,N_INTERNAL_MODEL_MAX
         INTERNAL_MODEL_LIST(i)      = 0
         SUBMODEL_FOR_IM(i)          = 0
      enddo   ! i over internal model list

      do im=1,INTERNAL_ID_MAX
         SUBMODEL_PARTITION_INDEX(im) = 0
         INTERNAL_MODEL_INDEX(im) = 0
         LAST_IM_IN_SM(im)=.false.
      enddo   ! im over internal model ids

      do s=1,N_SUBMODEL_PARTITION_MAX
         SUBMODEL_PARTITION_LIST(s)= 0
         SUBMODEL_FOR_SM(s)=0
      enddo  ! s over submodel list

      do sm=1,SUBMODEL_ID_MAX
         N_INTERNAL_FOR_SM(sm)      = 0
      enddo  ! sm over submodel ids

!
! 2. Obtain internal model and submodel identifiers from umui
!    generated namelist.
!
      read(5,NSUBMODL)
!
!
! 3. Check umui supplied values.
!
!
! 3.1 Check for umui supplied dimensions against parameter maxima.

      if(N_INTERNAL_MODEL >  N_INTERNAL_MODEL_MAX) then
         write(6,*) 'UM_Submodel_In: FATAL ERROR. Too many internal ',  &
     &   'models =',N_INTERNAL_MODEL,                                   &
     &   ' :You need to increase N_INTERNAL_MODEL_MAX'
         ErrorStatus=1       ! Set error flag
      endif
!
! 3.2 Check umui suppiled values are valid
!
      do i=1,N_INTERNAL_MODEL ! loop over internal models

        im = INTERNAL_MODEL_LIST(i) ! internal model identifier
        if(im <= 0.or.im >  INTERNAL_ID_MAX) then
         write(6,*) 'UM_Submodel_In: FATAL ERROR. Illegal internal ',   &
     &   'model identifier=',im,                                        &
     &   ' :Check values in namelist NSUBMODL supplied by umui'
         ErrorStatus=1       ! Set error flag
        endif

        sm = SUBMODEL_FOR_IM(i)     ! submodel for this internal model
        if(sm <= 0.or.sm >  SUBMODEL_ID_MAX) then
         write(6,*) 'UM_Submodel_In: FATAL ERROR. Illegal submodel ',   &
     &   'dump identifier=',sm,                                         &
     &   ' :Check values in namelist NSUBMODL supplied by umui'
         ErrorStatus=1       ! Set error flag
        endif

      enddo ! i=1,N_INTERNAL_MODEL
!
! 4. Form internal model and submodel description arrays.
!
      sm_prev = 0             ! Null value of submodel identifier
      N_SUBMODEL_PARTITION=0  ! Count no. of submodel partitions

      do i=1,N_INTERNAL_MODEL ! loop over internal models

        im = INTERNAL_MODEL_LIST(i) ! internal model identifier
        sm = SUBMODEL_FOR_IM(i)     ! submodel for this internal model
        INTERNAL_MODEL_INDEX(im)=i  ! sequence no. for STASH arrays

        if(sm /= sm_prev) then  ! new submodel

           N_SUBMODEL_PARTITION = N_SUBMODEL_PARTITION+1
           SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION) = sm

!   Since this is a new submodel, the previous internal model must be
!   the last internal model in its submodel partition.
           IF(N_SUBMODEL_PARTITION >  1) THEN ! Not first dump
              LAST_IM_IN_SM(im_prev) = .true.
           ENDIF

        endif                   ! test on new submodel
        SUBMODEL_FOR_SM(IM) = N_SUBMODEL_PARTITION

        SUBMODEL_PARTITION_INDEX(im)=sm
        N_INTERNAL_FOR_SM(sm)=N_INTERNAL_FOR_SM(sm)+1

        im_prev=im
        sm_prev=sm

      enddo ! i=1,N_INTERNAL_MODEL

      LAST_IM_IN_SM(im) = .true.  ! last im in list is last im in sm

!
! 5. Check calculated dimensions against parameter maxima.

      if(N_SUBMODEL_PARTITION >  N_SUBMODEL_PARTITION_MAX) then
         write(6,*) 'UM_Submodel_In: FATAL ERROR. Too many submodels =',&
     &   N_SUBMODEL_PARTITION,                                          &
     &   ' You need to increase N_SUBMODEL_PARTITION_MAX'
         ErrorStatus=1       ! Set error flag
      endif
!
!     Need a copy of No of submodels for use by TYPD1.
      ALT_N_SUBMODEL_PARTITION=N_SUBMODEL_PARTITION

      if (ALT_N_SUBMODEL_PARTITION_MAX /= N_SUBMODEL_PARTITION_MAX)THEN
        write(6,*)'UM_Submodel_In: Mismatch in parameters '
        WRITE(6,*)'N_SUBMODEL_PARTITION_MAX and '
        WRITE(6,*)'ALT_N_SUBMODEL_PARTITION_MAX. '
        WRITE(6,*)'They should be identical '
        ErrorStatus=1
        endif
      return
      END SUBROUTINE UM_Submodel_Init
