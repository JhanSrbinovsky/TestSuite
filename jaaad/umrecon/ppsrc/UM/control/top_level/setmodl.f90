
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Defines submodel and section/version configuration
!
! Subroutine Interface:
      SUBROUTINE SETMODL(ErrorStatus,CMESSAGE)


      Use Rcf_Submodel_Mod, Only :                                      &
     &    Internal_Model_List,                                          &
     &    Atmos_IM,         A_IM,                                       &
     &    A_SM

      Use Rcf_Recon_Mod, Only :                                         &
     &    TR_VARS

      Use Rcf_Model_Mod  ! Currently uses most of this module...

      Use Rcf_CntlAtm_Mod, Only :                                       &
     &    A_SW_RadStep,                                                 &
     &    A_LW_RadStep,                                                 &
     &    A_SW_RadStep_diag,                                            &
     &    A_SW_Radstep_prog,                                            &
     &    A_LW_RadStep_diag,                                            &
     &    A_LW_Radstep_prog,                                            &
     &    H_SWBands,                                                    &
     &    H_LWBands

      Use Rcf_Grid_Type_Mod, Only :                                     &
     &    Output_Grid

      Use Rcf_Lsm_Mod, Only :                                           &
     &    Glob_Land_out


      IMPLICIT NONE
! Description:
!
! Method:
!
! Current code owner:  S.J.Swarbrick
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Global variables:

! Array arguments with intent(out):
      CHARACTER*(80) CMESSAGE    ! Error return message

! Error status:
      INTEGER        ErrorStatus ! +ve = fatal error


! Local scalars
      REAL    ASteps !Atmos timesteps per hour
      INTEGER I,J
      INTEGER Im,Is
      INTEGER Obs_type(MAX_AOBS)

! Local parameters
      INTEGER NOASFLDS
      PARAMETER(NOASFLDS=4)

! Function and subroutine calls:
      INTEGER  CTOITST
!- End of Header ---------------------------------------------------

!  Define submodel configuration

      H_ATMOS = 'N'
      H_OCEAN = 'N'
      H_SLAB  = 'N'
      H_WAVE  = 'N'

      DO I = 1,N_INTERNAL_MODEL_MAX
        IF (INTERNAL_MODEL_LIST(I) == ATMOS_IM) THEN
          H_ATMOS ='Y'
        END IF
      END DO

        DO I = 1,N_INTERNAL_MODEL_MAX
          H_VERS(I,0)= 1
        END DO

      IF (H_ATMOS  ==  'Y') THEN
        Im = ATMOS_IM
        DO Is = 1,NSECTP
! DEPENDS ON: ctoitst
          H_VERS(Im,Is)= CTOITST(ATMOS_SR(Is)(1:1))
        END DO
      END IF
! Submodel independent

      Im = ATMOS_IM
      IF (H_ATMOS  ==  'Y') THEN
! Atmos model included

       IF(OCALB == 1) THEN
          H_FLOOR=FLOOR
          H_STRAT='N'
        ELSE
          H_FLOOR='Y'
          H_STRAT='Y'
        END IF

        DO I=1,A_MAX_TRVARS   ! Up to A_MAX_TRVARS=150 free tracers
          IF (TCA(I) == 0) THEN
            TRACER_A(I)=.FALSE.
          ELSE
            TRACER_A(I)=.TRUE.
          END IF
        END DO

        DO I=1,A_MAX_UKCAVARS ! Up to A_MAX_UKCAVARS=150 UKCA tracers
          IF (TC_UKCA(I) == 0) THEN
            TR_UKCA_A(I)=.FALSE.
          ELSE
            TR_UKCA_A(I)=.TRUE.
          END IF
        END DO

        IF(TOTAE == 'Y') THEN
          IF(TOTEM == 'Y') THEN
            H_TOTEM='Y'
          ELSE
            H_TOTEM='N'
          END IF
        ELSE
          H_TOTEM='N'
        END IF



        H_OROG_GRAD= (H_VERS(Im,6) == 3 .OR. H_VERS(Im,6) == 4)

        IF(OCAAA == 1) THEN
! Atmos global model
          H_GLOBAL(A_IM)='Y'

          H_A_EWSPACE=360.0/ Output_Grid % GLOB_P_ROW_LENGTH
          H_A_NSSPACE=180.0/( Output_Grid % GLOB_P_ROWS-1)

          H_A_FIRSTLAT=-90.0       ! S to N
          H_A_FIRSTLONG=0.0
          H_A_POLELAT=90.0
          H_A_POLELONG=0.0
          LMESO=.FALSE.
        ELSE IF (OCAAA == 2) THEN
! Atmos LAM
          H_GLOBAL(A_IM)='N'
          LMESO=(MESO == 'Y')
          H_A_EWSPACE=EWSPACEA
          H_A_NSSPACE=NSSPACEA
          H_A_FIRSTLAT=FRSTLATA
          H_A_FIRSTLONG=FRSTLONA
          IF (H_A_FIRSTLONG <  0) H_A_FIRSTLONG=H_A_FIRSTLONG+360.0
          H_A_POLELAT=POLELATA
          H_A_POLELONG=POLELONA
        ELSE IF (OCAAA == 3) THEN
! Atmos single column
          H_GLOBAL(A_IM)='N'
          H_A_EWSPACE=360.
          H_A_NSSPACE=180.
          H_A_FIRSTLAT=LATS
          H_A_FIRSTLONG=LONS
          H_A_POLELAT=90.0
          H_A_POLELONG=0.0
          LMESO=.FALSE.
        ELSE IF (OCAAA /= 0) THEN
          write(6,*)                                                    &
     &   'Setmodl: UNEXPECTED ATMOSPHERIC AREA CODE OCAAA',OCAAA
        END IF

        LEXTRA(A_SM)=(1+ Output_Grid % MODEL_LEVELS + 2 *               &
     &                Output_Grid % WET_LEVELS) *                       &
     &                Output_Grid % GLOB_P_ROWS *                       &
     &                Output_Grid % GLOB_P_ROW_LENGTH

      ELSE   ! Atmosphere model not included

        ZonAvOzone    =.FALSE.
!       L_MICROPHY    =.FALSE. ! Now hard-wired parameter
        H_OROG_GRAD   =.FALSE.
        LMESO         =.FALSE.
        H_FLOOR       ='N'
        H_STRAT       ='N'
        TOTAE         ='N'
        H_TOTEM       ='N'
        H_GLOBAL(A_IM)='N'

        GLOB_LAND_OUT =0
        Output_Grid % WET_LEVELS         =0
        Output_Grid % CLOUD_LEVELS       =0
        Output_Grid % TR_LEVELS          =0
        TR_VARS                          =0
        DO I=1,A_MAX_TRVARS   ! Up to A_MAX_TRVARS=150 free tracers
          TRACER_A(I) =.FALSE.
        END DO
        DO I=1,A_MAX_UKCAVARS ! Up to A_MAX_UKCAVARS=150 UKCA tracers
          TR_UKCA_A(I) =.FALSE.
        END DO
        StLevGWdrag   =0
        BotVDiffLev   =0
        TopVDifflev   =0
        A_SW_RADSTEP  =0
        A_LW_RADSTEP  =0
        A_SW_RADSTEP_DIAG  = 0
        A_SW_RADSTEP_PROG  = 0
        A_LW_RADSTEP_DIAG  = 0
        A_LW_RADSTEP_PROG  = 0
        H_OROG_ROUGH  =0
!       A_CONV_STEP   =0       ! Now hard-wired parameter
        H_SWBANDS     =0
        H_LWBANDS     =0
        H_A_EWSPACE   =0.
        H_A_NSSPACE   =0.
        H_A_FIRSTLAT  =0.
        H_A_FIRSTLONG =0.
        H_A_POLELAT   =0.
        H_A_POLELONG  =0.
!       PHENOL_PERIOD =0.      ! Now hard-wired parameter
!       TRIFFID_PERIOD=0.      ! Now hard-wired parameter
      END IF

      H_SLAB_CAL   ='N'


      CLOSE(UNIT=10)
      RETURN
      END SUBROUTINE SETMODL

!- End of Subroutine code ----------------------------------------------

! Function Interface:

