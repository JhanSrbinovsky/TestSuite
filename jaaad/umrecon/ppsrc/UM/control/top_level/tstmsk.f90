
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Checks version mask and option code in a ppx record
!
! Subroutine Interface:

      SUBROUTINE TSTMSK(MODL,ISEC,LMASK,LADRES,ErrorStatus,CMESSAGE)


      Use Rcf_Submodel_Mod, Only :                                      &
     &    A_IM,                                                         &
     &    Internal_Model_List,                                          &
     &    N_Internal_Model_Max

      Use Rcf_Address_Vars_Mod, Only :                                  &
     &    IOPN,                                                         &
     &    ISPACE,                                                       &
     &    VMSK

      Use Rcf_CntlAtm_Mod, Only :                                       &
     &    LMOSES,              L_SNOW_ALBEDO,                           &
     &    L_DUST,              L_USE_GRAD_CORR,     L_USE_BIOGENIC,     &
     &    L_USE_ARCLBIOM,      L_USE_ARCLBLCK,      L_USE_ARCLSSLT,     &
     &    L_USE_ARCLSULP,      L_USE_ARCLDUST,      L_USE_ARCLOCFF,     &
     &    L_USE_ARCLDLTA,                                               &
     &    L_USE_STOCHEM_CH4,   L_USE_STOCHEM_O3,                        &
     &    L_SULPC_SO2,         L_SO2_SURFEM,                            &
     &    L_SO2_HILEM,         L_SO2_NATEM,                             &
     &    L_SULPC_DMS,         L_DMS_EM,            L_DMS_OINTER,       &
     &    L_SULPC_OZONE,       L_SULPC_NH3,                             &
     &    L_NH3_EM,            L_VEG_FRACS,                             &
     &    L_TRIFFID,           L_LSPICE,            L_USE_CARIOLLE

      Use Rcf_CntlAtm_Mod, Only :                                       &
     &    L_CCRAD,             L_3D_CCA,            L_3D_CCW,           &
     &    L_SOOT,              L_SOOT_SUREM,        L_SOOT_HILEM,       &
     &    L_BIOMASS,           L_BMASS_SUREM,                           &
     &    L_BMASS_HILEM,                                                &
     &    L_OCFF,              L_OCFF_SUREM,                            &
     &    L_OCFF_HILEM,                                                 &
     &    L_RHCPT,             L_CO2_INTERACTIVE,                       &
     &    L_MICROPHY,          L_USE_SULPC_INDIRECT_SW,                 &
     &    L_CLD_AREA,                                                   &
     &    L_USE_SEASALT_DIRECT,L_USE_SEASALT_INDIRECT                   &
     &   ,L_GWD,               L_USE_USSP                               &
     &   ,L_CTILE,             L_RIVERS,            L_INLAND            &
     &   ,L_OASIS

      Use Rcf_CntlAtm_Mod, Only :                                       &
     &    i_ozone_int,                                                  &
     &    L_USE_TPPS_OZONE,    L_PC2                                    &
     &   ,L_mcr_qcf2,          L_mcr_qrain                              &
     &   ,L_mcr_qgraup                                                  &
     &   ,CAN_MODEL,           L_TOP

      Use Rcf_Model_Mod, Only :                                         &
     &    MEAN_NUMBER,         H_GLOBAL,                                &
     &    AASSET,                                                       &
     &    A_Max_TrVars,        TRACER_A,                                &
     &    A_Max_UKCAVars,      TR_UKCA_A,                               &
     &    H_VERS,              ATMOS_SR,                                &
     &    TOTAE,               SSTANOM,                                 &
     &    OROGR,                                                        &
     &    H_SLAB_CAL,                                                   &
     &    H_OROG_GRAD,         H_STRAT,                                 &
     &    H_TOTEM,                                                      &
     &    H_FLOOR

      Use Rcf_Model_Mod, Only :                                         &
     &    H_SLAB,                                                       &
     &    H_OCEAN,                                                      &
     &    OCBOHANEY

      Use Rcf_o3intp_Mod


      Use Rcf_Recon_Mod, Only :                                         &
     &    Var_Recon,                                                    &
     &    RIMWIDTHA,                                                    &
     &    TR_VARS,                                                      &
     &    NICE


      IMPLICIT NONE

! Description:
!   Determines whether a diagnostic is available to a particular
!   (internal model,section) version. Also checks the option code IOPN
!   Called by INACTR, PRELIM.
!
! Method:
!
! The decimal value of the version mask, VMSK, was read from the
! ppxref file by GETPPX, each (model, section, item). The version
! number of the (model, section) - NMASK - is obtained from the
! H_VERS array.
!
! The procedure for checking diagnostic availability is as follows:
! (1) Check whether the relevant internal model is included in the
!     submodel configuration, by examining the INTERNAL_MODEL_LIST
!     array; if not, return.
! (2) Check whether NMASK=0 - this implies that the diag is
!      unavailable to any version. If so, return.
! (3) Check whether the diag is available to the specified version.
!     The 'version mask' binary code in the STASHmaster specifies
!     which versions the diag is available to. Eg., if it is available
!     to vns. 1 and 3, but not vn.2, the version mask is 101. VMSK
!     is the decimal equivalent of the version mask (5 in this example).
!     TSTMSK therefore calculates the quantity
!
!             IMOD = MOD(VMSK,2**NMASK)/2**(NMASK-1)
!
!       If IMOD=1, the diag is available; if IMOD=0, it isn't.
!
!   The option code is checked in accordance with the ppxref option
!   code definitions.
!
!   Note for code developers adding new tests: LMASK is initialised to
!   .TRUE. at top of deck. A series of tests is made and if any one
!   test fails LMASK is reset to .FALSE. Therefore do not reinitialise
!   LMASK to .TRUE. anywhere otherwise you may inadvertently overwrite
!   a preceding .FALSE. setting or you may mislead future code writers
!   who may do so.
!    For this reason, all unnecessary LMASK=.TRUE. were removed at 5.1
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Global variables:

! Subroutine arguments:
!   Scalar arguments with intent(in):
      INTEGER MODL    ! Internal model number
      INTEGER ISEC    ! Section number

!   Scalar arguments with intent(out):
      LOGICAL LMASK   ! T if diag is available to specified version
      LOGICAL LADRES  ! T if diag is available for addressing primary
      CHARACTER*80 CMESSAGE

! Local scalars
      INTEGER NMASK   ! Version number for (model,section)
      INTEGER IMOD    ! Determines whether diag is available
      INTEGER TWONM   ! Used in calculation of IMOD
      INTEGER TWONM1  ! Used in calculation of IMOD
      INTEGER COUNT
      INTEGER I
      INTEGER N0,N1,N2,N2N1,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13
      INTEGER N14,N15,N16,N17,N18
      INTEGER N19,N20,N21,N22,N23,N24,N25,N26,N27,N28,N29,N30
      INTEGER SUM_IOPN

! ErrorStatus
      INTEGER ErrorStatus

!- End of Header ---------------------------------------------------

      SUM_IOPN=IOPN(1)+IOPN(2)+IOPN(3)+IOPN(4)+IOPN(5)+IOPN(6)
!--------------------------------------------------------------------
!Check whether internal model is included in submodel configuration
!--------------------------------------------------------------------
      COUNT = 0
      DO I = 1,N_INTERNAL_MODEL_MAX
      IF (INTERNAL_MODEL_LIST(I) == MODL) THEN
        LMASK =.TRUE.
        LADRES=.TRUE.
        COUNT = COUNT + 1
      END IF
      END DO
      IF (COUNT == 0) THEN
        LMASK =.FALSE.
        LADRES=.FALSE.
        GO TO 9999
      END IF

!--------------------------------------------------------------------
!Check whether diagnostic is unavailable to any version
!--------------------------------------------------------------------
      NMASK=H_VERS(MODL,ISEC)
      IF (NMASK == 0) THEN
        LMASK =.FALSE.
        LADRES=.FALSE.
        GOTO 9999
      END IF

! Determine whether the diag is available to the specified version
      TWONM  = 2**NMASK
      TWONM1 = 2**(NMASK-1)
      IMOD   = MOD(VMSK,TWONM)/TWONM1
      IF(IMOD == 1) THEN
        LMASK =.TRUE.
      ELSE IF(IMOD == 0) THEN
        LMASK =.FALSE.
        LADRES=.FALSE.
        GOTO 9999
      ELSE
        WRITE(6,*)'S: TSTMSK INVALID DECODING OF VMSK',VMSK
        WRITE(6,*)'S: ... IMOD=',IMOD,'     NMASK=',NMASK
      END IF

!-------------------------------------------------------------------
!Check option codes
!-------------------------------------------------------------------
      LMASK=.TRUE.
      IF (MODL == A_IM) THEN
!-------------------------------------------------------------------
!Var reconfiguration now is not restricted to primary fields.  
!Therefore need to perform on all sections.
!-------------------------------------------------------------------
        IF ( VAR_RECON ) THEN                                           
          N17 = MOD((IOPN(4)/10),10)
          IF (N17 /= 1) THEN
            LMASK = .FALSE.
          END IF
        END IF

        IF((ISEC == 0) .OR. (ISEC == 31) .OR.                           &
     &     (ISEC == 32) .OR. (ISEC == 33) .OR.                          &
     &     (ISEC == 34) ) THEN
!Atmosphere primary field
            IF(SUM_IOPN /= 0) THEN
              N2N1=MOD (IOPN(1),100)
! Up to A_MAX_TRVARS=150 free tracers now allowed in section 33
              IF (ISEC == 33)   N2N1=MOD (IOPN(1),1000)
! Up to A_MAX_UKCAVARS=150 UKCA tracers now allowed in section 34
              IF (ISEC == 34)   N2N1=MOD (IOPN(1),1000)
              N3  =MOD((IOPN(1)/100),10)
              N4  =MOD((IOPN(1)/1000),10)
              N5  =MOD((IOPN(1)/10000),10)
              N6  =MOD( IOPN(2),10)
              N7  =MOD((IOPN(2)/10),10)
              N8  =MOD((IOPN(2)/100),10)
              N9  =MOD((IOPN(2)/1000),10)
              N10 =MOD((IOPN(2)/10000),10)
              N11 =MOD( IOPN(3),10)
              N12 =MOD((IOPN(3)/10),10)
              N13 =MOD((IOPN(3)/100),10)
              N14 =MOD((IOPN(3)/1000),10)
              N15 =MOD((IOPN(3)/10000),10)
              N16 =MOD( IOPN(4),10)
              N17 =MOD((IOPN(4)/10),10)
              N18 =MOD((IOPN(4)/100),10)
              N19 =MOD((IOPN(4)/1000),10)
              N20 =MOD((IOPN(4)/10000),10)
              N21 =MOD( IOPN(5),10)
              N22 =MOD((IOPN(5)/10),10)
              N23 =MOD((IOPN(5)/100),10)
              N24 =MOD((IOPN(5)/1000),10)
              N25 =MOD((IOPN(5)/10000),10)
              N26 =MOD( IOPN(6),10)
              N27 =MOD((IOPN(6)/10),10)
              N28 =MOD((IOPN(6)/100),10)
              N29 =MOD((IOPN(6)/1000),10)
              N30 =MOD((IOPN(6)/10000),10)
              IF( ISEC == 33 .AND.                                      &
     &            (N2N1 >  0) .AND. (N2N1 <= A_MAX_TRVARS) .AND.        &
     &            .NOT. TRACER_A(N2N1) ) THEN
                LMASK=.FALSE.
              ELSE IF( ISEC == 34 .AND.                                 & 
     &                 (N2N1 >  0) .AND. (N2N1 <= A_MAX_UKCAVARS) .AND. &
     &                 .NOT. TR_UKCA_A(N2N1) ) THEN
                LMASK=.FALSE.
              
!             N3=1 not used
!             N3=2 not used

              ELSE IF((N3 == 3).AND.(.NOT.LMOSES)) THEN
                LMASK=.FALSE.
              ELSE IF((N3 == 4).AND.(.NOT.LMOSES) )THEN
                LMASK=.FALSE.
              ELSE IF((N3 == 5).AND.(.NOT.L_TOP)) THEN
                LMASK=.FALSE.
              ELSE IF((N4 == 1).AND.(.NOT.L_USE_BIOGENIC)) THEN
                LMASK=.FALSE.
              ELSE IF((N4 == 2).AND.(.NOT.L_USE_ARCLBIOM)) THEN
                LMASK=.FALSE.
              ELSE IF((N4 == 3).AND.(.NOT.L_USE_ARCLBLCK)) THEN
                LMASK=.FALSE.
              ELSE IF((N4 == 4).AND.(.NOT.L_USE_ARCLSSLT)) THEN
                LMASK=.FALSE.
              ELSE IF((N4 == 5).AND.(.NOT.L_USE_ARCLSULP)) THEN
                LMASK=.FALSE.
              ELSE IF((N4 == 6).AND.(.NOT.L_USE_ARCLDUST)) THEN
                LMASK=.FALSE.
              ELSE IF((N4 == 7).AND.(.NOT.L_USE_ARCLOCFF)) THEN
                LMASK=.FALSE.
              ELSE IF((N4 == 8).AND.(.NOT.L_USE_ARCLDLTA)) THEN
                LMASK=.FALSE.
              ELSE IF((N4 > 8)) THEN
                LMASK=.FALSE.  ! these are not yet defined
              ELSE IF((N5 == 1).AND.(OROGR == 'N')) THEN
                LMASK=.FALSE.
              ELSE IF((N6 == 1) .AND.                                   &
     &               (RIMWIDTHA == 0 .OR. H_GLOBAL(A_IM) == 'Y')) THEN
                LMASK=.FALSE.
              ELSE IF((N6 == 2).AND.(H_FLOOR == 'N')) THEN
                LMASK=.FALSE.
              ELSE IF((N7 == 1).AND.(.NOT.L_OASIS)) THEN
                LMASK=.FALSE.
              ELSE IF((N7 == 2)) THEN
                LMASK=.FALSE.
              ELSE IF((N7 == 7).AND.(.NOT.(L_DMS_OINTER))) THEN
                LMASK=.FALSE.
              ELSE IF((N8 == 1).AND.(SSTAnom == 'N')) THEN
                LMASK=.FALSE.
              ELSE IF((N8 == 4).AND.(TOTAE /= 'Y')) THEN
                LMASK=.FALSE.
              ELSE IF((N8 == 5).AND.(H_TOTEM /= 'Y')) THEN
                LMASK=.FALSE.
              ELSE IF((N8 == 6).AND.(.NOT.L_SNOW_ALBEDO)) THEN
                LMASK=.FALSE.

!             N8=7 not used

              ELSEIF ( ( N8  ==  8) .AND. (ATMOS_SR(14) == '0A') ) THEN
                LMASK=.FALSE.
              ELSE IF((N9 == 1).AND.(.NOT.H_OROG_GRAD)) THEN
                LMASK=.FALSE.
              ELSE IF((N9 == 2).AND.(.NOT.L_USE_GRAD_CORR)) THEN
                LMASK=.FALSE.   ! Gradient correction for radiation off
              ELSE IF((N10 == 1).AND.(.NOT.L_SULPC_SO2)) THEN
                LMASK=.FALSE.
              ELSE IF((N10 == 2).AND.((.NOT.L_SULPC_SO2)                &
     &                          .OR.  (.NOT.L_SO2_SURFEM )) ) THEN
                LMASK=.FALSE.
              ELSE IF((N10 == 3).AND.((.NOT.L_SULPC_SO2)                &
     &                          .OR.  (.NOT.L_SO2_HILEM )) ) THEN
                LMASK=.FALSE.
              ELSE IF((N10 == 4).AND.((.NOT.L_SULPC_SO2)                &
     &                          .OR.  (.NOT.L_SO2_NATEM )) ) THEN
                LMASK=.FALSE.
              ELSE IF((N10 == 5).AND.((.NOT.L_SULPC_SO2)                &
     &                          .OR.  (.NOT.L_SULPC_DMS)) ) THEN
                LMASK=.FALSE.
              ELSE IF((N10 == 6).AND.((.NOT.L_SULPC_SO2)                &
     &                          .OR.  (.NOT.L_SULPC_DMS)                &
     &                          .OR.  (.NOT.L_DMS_EM    )) ) THEN
                LMASK=.FALSE.
              ELSE IF((N10 == 7).AND.((.NOT.L_SULPC_SO2)                &
     &                          .OR.  (.NOT.L_SULPC_OZONE )) ) THEN
                LMASK=.FALSE.
             ELSE IF((N10 == 8).AND.((.NOT.L_SULPC_SO2)                 &
     &                         .OR.  (.NOT.L_SULPC_OZONE)               &
     &                         .OR.  (.NOT.L_SULPC_NH3))  )  THEN
               LMASK=.FALSE.
             ELSE IF((N10 == 9).AND.((.NOT.L_SULPC_SO2)                 &
     &                         .OR.  (.NOT.L_SULPC_OZONE)               &
     &                         .OR.  (.NOT.L_SULPC_NH3)                 &
     &                         .OR.  (.NOT.L_NH3_EM))  )  THEN
               LMASK=.FALSE.
              ELSE IF((N11 == 1).AND.(L_VEG_FRACS.OR.L_TRIFFID) )THEN
                LMASK=.FALSE.
              ELSE IF((N11 == 2).AND.(.NOT.L_VEG_FRACS))THEN
                LMASK=.FALSE.
              ELSE IF((N11 == 3).AND.                                   &
     &                (.NOT.L_VEG_FRACS.OR..NOT.L_TRIFFID) )THEN
                LMASK=.FALSE.
              ELSE IF((N12 == 1) .AND. (L_LSPICE)) THEN
                LMASK=.FALSE.   ! QCF is secondary 403 (not prognostic)
              ELSE IF((N12 == 2) .AND. (.NOT.L_LSPICE)) THEN
                LMASK=.FALSE.   ! QCF is primary 12 (prognostic)
              ELSE IF((N12 == 3) .AND. (.NOT.L_mcr_qcf2)) THEN
                LMASK=.FALSE.   ! QCF2 is prognostic
              ELSE IF((N12 == 4) .AND. (.NOT.L_mcr_qrain)) THEN
                LMASK=.FALSE.   ! QRAIN is prognostic
              ELSE IF((N12 == 5) .AND. (.NOT.L_mcr_qgraup)) THEN
                LMASK=.FALSE.   ! QGRAUP is prognostic
              ELSE IF((N12 == 6) .AND. (.NOT.L_pc2)) THEN
                LMASK=.FALSE.   ! PC2 is active
              ELSE IF ((N13 == 1).AND.(L_3D_CCA)) THEN
                LMASK=.FALSE.    !  CCA is 2D
              ELSE IF ((N13 == 2).AND.(.NOT.L_3D_CCA)) THEN
                LMASK=.FALSE.    !  CCA is 3D
              ELSE IF((N14 == 1).AND.(.NOT.L_SOOT))  THEN
                LMASK=.FALSE.
              ELSE IF((N14 == 2).AND.((.NOT.L_SOOT)                     &
     &                          .OR.  (.NOT.L_SOOT_SUREM)) )  THEN
                LMASK=.FALSE.
              ELSE IF((N14 == 3).AND.((.NOT.L_SOOT)                     &
     &                         .OR.  (.NOT.L_SOOT_HILEM)) )  THEN
                LMASK=.FALSE.
              ELSE IF ((N16 == 1).AND.(.NOT.L_CO2_INTERACTIVE)) THEN
                LMASK=.FALSE.    !  Carbon cycle switched off
!N17 is left out here due to test is now on all sections not just primarys.
              ELSE IF ((N18 == 1).AND.(.NOT.L_CTILE)) THEN
                LMASK=.FALSE.    !  coastal tiling switched off
              ELSE IF ((N19 == 1).AND.(.NOT.L_USE_TPPS_OZONE)) THEN
                LMASK=.FALSE.    !  `tpps_ozone' switched off
              ELSE IF ((N20 == 1).AND.(CAN_MODEL /= 4)) THEN
                LMASK=.FALSE.    !  Snow canopy switched off
              ELSE IF ((N21 == 1).AND.(.NOT.L_BIOMASS)) THEN
                LMASK=.FALSE.
              ELSE IF ((N21 == 2).AND.((.NOT.L_BIOMASS)                 &
     &                           .OR. (.NOT.L_BMASS_SUREM)) ) THEN
                LMASK=.FALSE.
              ELSE IF ((N21 == 3).AND.((.NOT.L_BIOMASS)                 &
     &                           .OR. (.NOT.L_BMASS_HILEM)) ) THEN
                LMASK=.FALSE.
              ELSE IF ((N22 == 1).AND.(.NOT.L_RIVERS)) THEN
                LMASK=.FALSE.    !  River routing switched off
              ELSE IF ((N22 == 2).AND.(.NOT.L_INLAND)) THEN
                LMASK=.FALSE.    !Inland re-routing switched off
              ELSE IF ((N23 == 1).AND.(NICE == 1)) THEN
                LMASK=.FALSE.    ! Sea ice catagories switched off
              ELSE IF ((N24 == 1).AND.(.NOT.L_DUST)) THEN
                LMASK=.FALSE.    !  mineral dust switched off
              ELSE IF ((N25 == 1).AND. (.NOT.L_USE_STOCHEM_CH4)) THEN
                LMASK=.FALSE.   ! STOCHEM methane switched off
              ELSE IF ((N25 == 2).AND. (.NOT.L_USE_STOCHEM_O3)) THEN
                LMASK=.FALSE.   ! STOCHEM ozone switched off
              ELSE IF ((N26 == 1).AND.(ATMOS_SR(25) == '0A' .AND.       &
     &          ATMOS_SR(34) == '0A')) THEN
                LMASK=.FALSE.   ! Direct PAR prognostic switched off
              ELSE IF ((N27 == 1).AND.(.NOT.L_CCRAD)) THEN
                LMASK=.FALSE.   ! CCRad is off
              ELSE IF ((N27 == 2).AND.                                  &
                      (.NOT.(L_CCRAD.AND.L_3D_CCW))) THEN      
                LMASK=.FALSE.   ! Either L_CCRAD or L_3d_ccw is false
              ELSE IF ((N28 == 1).AND.(.NOT.L_USE_CARIOLLE)) THEN      
                LMASK=.FALSE.   ! L_USE_CARIOLLE is false
              ELSE IF ((N29 == 1).AND.(.NOT.L_OCFF)) THEN
                LMASK=.FALSE.
              ELSE IF ((N29 == 2).AND.((.NOT.L_OCFF)                    &
     &                           .OR. (.NOT.L_OCFF_SUREM)) ) THEN
                LMASK=.FALSE.
              ELSE IF ((N29 == 3).AND.((.NOT.L_OCFF)                    &
     &                           .OR. (.NOT.L_OCFF_HILEM)) ) THEN
                LMASK=.FALSE.
              END IF
            END IF
        END IF
!End of atmosphere primary block
      END IF

      IF (MODL == A_IM) THEN
!Atmos diagnostics
      IF(ISEC == 1) THEN
!Shortwave radiation
        IF (SUM_IOPN /= 0) THEN
          N1=MOD (IOPN(1)    ,10)
          N2=MOD((IOPN(1)/10),10)
          N3 = MOD((IOPN(1)/100),10)
          N4 = MOD((IOPN(1)/1000),10)
          N5 = MOD((IOPN(1)/10000),10)
          N6 = MOD( IOPN(2),10)
          IF((N1 == 1).AND.(H_GLOBAL(A_IM) /= 'Y')) THEN
            LMASK=.FALSE.
          ELSE IF((N2 == 1).AND.(.NOT.L_MICROPHY)) THEN
            LMASK=.FALSE.
          ELSEIF ((N4  ==  1).AND.(.NOT. L_USE_SULPC_INDIRECT_SW)) THEN
            LMASK = .FALSE.
          ELSEIF ((N5  ==  1).AND.                                      &
     &           (.NOT.(L_USE_SEASALT_DIRECT                            &
     &             .OR. L_USE_SEASALT_INDIRECT))) THEN
            LMASK = .FALSE.
          ELSEIF ( ( N6  ==  1 ) .AND. ( .NOT. L_PC2 ) ) THEN
            LMASK = .FALSE.
          END IF
        END IF
      ELSE IF(ISEC == 2) THEN
!Longwave  radiation
        IF (SUM_IOPN /= 0) THEN
          N1=MOD(IOPN(1)   ,10)
          IF ((N1 == 1).AND.                                            &
     &        (I_OZONE_INT /= IO3_TROP_MAP_MASSCON)) THEN
             LMASK=.FALSE.
          END IF
          N6 = MOD( IOPN(2),10)
          IF ( ( N6  ==  1 ) .AND. ( .NOT. L_PC2 ) ) THEN
            LMASK = .FALSE.
          END IF
        END IF
      ELSE IF(ISEC == 3) THEN
!Boundary layer
        IF (SUM_IOPN /= 0) THEN
          N1=MOD (IOPN(1)    ,10)
          N2=MOD((IOPN(1)/10),10)
          N3=MOD((IOPN(1)/100),10)
          N4=MOD((IOPN(1)/1000),10)
          N6 = MOD( IOPN(2),10)
          IF((N1 == 1).AND.(OROGR == 'N')) THEN
            LMASK=.FALSE.
          END IF
          IF((N2 == 1).AND.(.NOT.L_SULPC_SO2)) THEN
            LMASK=.FALSE.
          ELSE IF ((N2 == 2).AND.(.NOT.L_SULPC_NH3)) THEN
            LMASK=.FALSE.
          ELSE IF ((N2 == 3).AND.(.NOT.L_SOOT)) THEN
            LMASK=.FALSE.
          ELSE IF ((N2 == 4).AND.(.NOT.L_BIOMASS)) THEN
            LMASK=.FALSE.
          ELSE IF ((N2 == 5).AND.(.NOT.L_DUST)) THEN
            LMASK=.FALSE.
          ELSE IF ((N2 == 6).AND.(.NOT.L_OCFF)) THEN
            LMASK=.FALSE.
          END IF
          IF((N3 == 1).AND.(.NOT.L_CO2_INTERACTIVE)) THEN
            LMASK=.FALSE.
          ELSEIF ( ( N6  ==  1 ) .AND. ( .NOT. L_PC2 ) ) THEN
            LMASK = .FALSE.
          END IF
          IF((N4 == 1) .AND. (nice == 1))THEN
            LMASK = .FALSE.
          ENDIF
        END IF
      ELSE IF(ISEC == 4) THEN
!Large-scale precipitation
          N2=MOD((IOPN(1)/10),10)
          N6 = MOD( IOPN(2),10)
          IF((N2 == 1).AND.(.NOT.L_SULPC_SO2)) THEN
            LMASK=.FALSE.
          ELSE IF ((N2 == 2).AND.(.NOT.L_SULPC_NH3)) THEN
            LMASK=.FALSE.
          ELSE IF ((N2 == 3).AND.(.NOT.L_SOOT)) THEN
            LMASK=.FALSE.
          ELSE IF ((N2 == 4).AND.(.NOT.L_BIOMASS)) THEN
            LMASK=.FALSE.
          ELSE IF ((N2 == 5).AND.(.NOT.L_DUST)) THEN
            LMASK=.FALSE.
          ELSE IF ((N2 == 6).AND.(.NOT.L_OCFF)) THEN
            LMASK=.FALSE.
          END IF
          IF ( ( N6  ==  1 ) .AND. ( .NOT. L_PC2 ) ) THEN
            LMASK = .FALSE.
          END IF
      ELSE IF(ISEC == 5) THEN
!Convection
          N2=MOD((IOPN(1)/10),10)
          N3=MOD((IOPN(1)/100),10)
          N6=MOD( IOPN(2),10)

          IF((N2 == 1).AND.(.NOT.L_SULPC_SO2)) THEN
            LMASK=.FALSE.
          ELSE IF ( (N3 == 1).AND.(L_3D_CCA) ) THEN
            LMASK=.FALSE.
          ELSE IF ( (N3 == 2).AND.(.NOT.L_3D_CCA) ) THEN
            LMASK=.FALSE.
          ELSE IF ((N2 == 2).AND.(.NOT.L_SULPC_NH3)) THEN
            LMASK=.FALSE.
          ELSE IF ((N2 == 3).AND.(.NOT.L_SOOT)) THEN
            LMASK=.FALSE.
          ELSE IF ((N2 == 4).AND.(.NOT.L_BIOMASS)) THEN
            LMASK=.FALSE.
          ELSE IF ((N2 == 5).AND.(.NOT.L_DUST)) THEN
            LMASK=.FALSE.
          ELSE IF ((N2 == 6).AND.(.NOT.L_OCFF)) THEN
            LMASK=.FALSE.
          END IF
          IF ( ( N6  ==  1 ) .AND. ( .NOT. L_PC2 ) ) THEN
            LMASK = .FALSE.
          END IF

      ELSE IF(ISEC == 6) THEN
!Gravity Wave Drag parametrization
        IF (SUM_IOPN /= 0) THEN
          N2=MOD((IOPN(1)/10),10)
          N3=MOD((IOPN(1)/100),10)
          IF((N2 == 1).AND.(.NOT.L_GWD)) THEN
            LMASK=.FALSE.
          ELSE IF ((N3 == 1).AND.(.NOT.L_USE_USSP)) THEN
            LMASK=.FALSE.
          END IF
        END IF


      ELSE IF(ISEC == 8) THEN
!Hydrology
        IF (SUM_IOPN /= 0) THEN
          N1=MOD (IOPN(1)    ,10)
          N22=MOD((IOPN(5)/10),10)
          IF((N1 == 1).AND.(.NOT.L_TOP)) THEN
            LMASK=.FALSE.
          ELSE IF ((N22 == 1).AND.(.NOT.L_RIVERS)) THEN
                LMASK=.FALSE.    !River routing switched off
          ELSE IF ((N22 == 2).AND.(.NOT.L_INLAND)) THEN
                LMASK=.FALSE.    !Inland re-routing switched off
          END IF
        END IF

      ELSE IF(ISEC == 9) THEN
!Cloud parametrization
        IF (SUM_IOPN /= 0) THEN
          N2=MOD((IOPN(1)/10),10)
          N3=MOD((IOPN(1)/100),10)
          N6 = MOD( IOPN(2),10)
          IF((N2 == 1).AND.(.NOT.L_CLD_AREA)) THEN
            LMASK=.FALSE.
          ELSE IF((N3 == 1).AND.(.NOT.L_RHCPT)) THEN
            LMASK=.FALSE.
          ELSEIF ( ( N6  ==  1 ) .AND. ( .NOT. L_PC2 ) ) THEN
            LMASK = .FALSE.
          END IF
        END IF
      ELSE IF(ISEC == 12) THEN
!Dynamics Advection
        IF (SUM_IOPN /= 0) THEN
          N6 = MOD( IOPN(2),10)
          IF ( ( N6  ==  1 ) .AND. ( .NOT. L_PC2 ) ) THEN
            LMASK = .FALSE.
          END IF
          END IF
      ELSE IF(ISEC == 16) THEN
!Extra physics
        IF (SUM_IOPN /= 0) THEN
          N2N1=MOD(IOPN(1),100)
          N6 = MOD( IOPN(2),10)
          IF((N2N1 /= 0).AND.(.NOT.TRACER_A(N2N1))) THEN
            LMASK=.FALSE.
          END IF
          IF ( ( N6  ==  1 ) .AND. ( .NOT. L_PC2 ) ) THEN
            LMASK = .FALSE.
          END IF
        END IF
      ELSE IF (ISEC == 17) THEN
! Chemistry section
        N1=MOD(IOPN(1),     10)
        N2=MOD((IOPN(1)/10),10)
        IF ((N1 == 1).AND.(.NOT.L_SULPC_SO2)) THEN
            LMASK=.FALSE.
        ELSE IF ((N1 == 2).AND.(.NOT.L_SOOT)) THEN
            LMASK=.FALSE.
        ELSE IF ((N1 == 3).AND.(.NOT.L_BIOMASS)) THEN
            LMASK=.FALSE.
        END IF
        IF ((N2 == 1).AND.(.NOT.L_SULPC_DMS)) THEN
            LMASK=.FALSE.
        ELSE IF ((N2 == 2).AND.(.NOT.L_SULPC_OZONE)) THEN
            LMASK=.FALSE.
       ELSE IF (N2 == 3) THEN
           LMASK = .FALSE.
        END IF
      ELSE IF(ISEC == 18) THEN
!Data assimilation
        IF (SUM_IOPN /= 0) THEN
          N1=MOD (IOPN(1)    ,10)
          N2=MOD((IOPN(1)/10),10)
! This code appears to be obselete
!   Section 18 history:
!     Used for assimilation at vn 4.5
!     Not used for vn 5.x
!     Reintroduced at vn 6.1
!          IF(.NOT.AASSET(N1)) THEN
!            LMASK=.FALSE.
!          ELSE IF((N2 == 1).AND.(TOTAE == 'N')) THEN
!            LMASK=.FALSE.
!          END IF
        END IF
      ELSE IF((ISEC >= 1).AND.(ISEC <= 20)) THEN  ! BUT NOT 1,3,18
        IF (SUM_IOPN /= 0) THEN
          WRITE(6,*)                                                    &
     &   'MESSAGE FROM ROUTINE TSTMSK: UNEXPECTED OPTION CODE',         &
     &    IOPN(1)
          WRITE(6,*)'IN ATMOSPHERE SECTION ',ISEC
        END IF
      ELSE IF((ISEC >= 21).AND.(ISEC <= 24)) THEN
!Atmos climate mean diagnostics
        N0=ISEC-20
        IF(N0 >  MEAN_NUMBER(A_IM)) THEN
          LMASK=.FALSE.
        ELSE
          IF(SUM_IOPN /= 0) THEN
          N2N1=MOD( IOPN(1),1000)  ! Up to A_MAX_TRVARS=150 free tracers
! Not sure whether this will work for free tracers in section 33;
! Not set up for UKCA tracers in section 34.
            N3=MOD((IOPN(1)/100),10)
            N4=MOD((IOPN(1)/1000),10)
            N5=MOD((IOPN(1)/10000),10)
            N6=MOD( IOPN(2),10)
            N7=MOD((IOPN(2)/10),10)
            N8=MOD((IOPN(2)/100),10)
            N9=MOD((IOPN(2)/1000),10)
           N10=MOD((IOPN(2)/10000),10)
           N11=MOD( IOPN(3),10)
           N12=MOD((IOPN(3)/10),10)
           N13 =MOD((IOPN(3)/100),10)
           N14 =MOD((IOPN(3)/1000),10)
           N15 =MOD((IOPN(3)/10000),10)
           N16 =MOD( IOPN(4),10)
            IF((N2N1 /= 0).AND.(.NOT.TRACER_A(N2N1))) THEN
              LMASK=.FALSE.
            ELSE IF((N3 == 3).AND.(.NOT.LMOSES)) THEN
              LMASK=.FALSE.
            ELSE IF((N3 == 4).AND.(.NOT.LMOSES) ) THEN
              LMASK=.FALSE.
            ELSE IF((N4 == 1).AND.(H_STRAT == 'Y')) THEN
              LMASK=.FALSE.
            ELSE IF((N5 == 1).AND.(OROGR == 'N')) THEN
              LMASK=.FALSE.
            ELSE IF((N6 == 1).AND.(H_GLOBAL(A_IM) == 'Y')) THEN
              LMASK=.FALSE.
            ELSE IF((N6 == 2).AND.(H_FLOOR == 'N')) THEN
              LMASK=.FALSE.
            ELSE IF((N7 == 1).AND.(.NOT.L_OASIS)) THEN
              LMASK=.FALSE.
            ELSE IF((N8 == 1).AND.(SSTAnom == 'N')) THEN
              LMASK=.FALSE.
            ELSE IF((N8 == 4).AND.(TOTAE /= 'Y')) THEN
              LMASK=.FALSE.
            ELSE IF((N8 == 5).AND.(H_TOTEM /= 'Y')) THEN
              LMASK=.FALSE.
            ELSE IF((N8 == 6).AND.(.NOT.L_SNOW_ALBEDO)) THEN
              LMASK=.FALSE.
            ELSE IF((N9 == 1).AND.(.NOT.H_OROG_GRAD)) THEN
              LMASK=.FALSE.
            ELSE IF((N10 == 1).AND.(.NOT.L_SULPC_SO2)) THEN
              LMASK=.FALSE.
            ELSE IF((N10 == 2).AND.((.NOT.L_SULPC_SO2)                  &
     &                        .OR.  (.NOT.L_SO2_SURFEM )) ) THEN
              LMASK=.FALSE.
            ELSE IF((N10 == 3).AND.((.NOT.L_SULPC_SO2)                  &
     &                        .OR.  (.NOT.L_SO2_HILEM )) ) THEN
              LMASK=.FALSE.
            ELSE IF((N10 == 4).AND.((.NOT.L_SULPC_SO2)                  &
     &                        .OR.  (.NOT.L_SO2_NATEM )) ) THEN
              LMASK=.FALSE.
            ELSE IF((N10 == 5).AND.((.NOT.L_SULPC_SO2)                  &
     &                        .OR.  (.NOT.L_SULPC_DMS)) ) THEN
              LMASK=.FALSE.
            ELSE IF((N10 == 6).AND.((.NOT.L_SULPC_SO2)                  &
     &                        .OR.  (.NOT.L_SULPC_DMS)                  &
     &                        .OR.  (.NOT.L_DMS_EM    )) ) THEN
              LMASK=.FALSE.
            ELSE IF((N10 == 7).AND.((.NOT.L_SULPC_SO2)                  &
     &                        .OR.  (.NOT.L_SULPC_OZONE )) ) THEN
              LMASK=.FALSE.
            ELSE IF((N10 == 8).AND.((.NOT.L_SULPC_SO2)                  &
     &                        .OR.  (.NOT.L_SULPC_OZONE)                &
     &                        .OR.  (.NOT.L_SULPC_NH3))  )  THEN
              LMASK=.FALSE.
            ELSE IF((N10 == 9).AND.((.NOT.L_SULPC_SO2)                  &
     &                        .OR.  (.NOT.L_SULPC_OZONE)                &
     &                        .OR.  (.NOT.L_SULPC_NH3)                  &
     &                        .OR.  (.NOT.L_NH3_EM))  )  THEN
              LMASK=.FALSE.
            ELSE IF((N11 == 1).AND.(L_VEG_FRACS.OR.L_TRIFFID) )THEN
              LMASK=.FALSE.
            ELSE IF((N11 == 2).AND.(.NOT.L_VEG_FRACS))THEN
              LMASK=.FALSE.
            ELSE IF((N11 == 3).AND.                                     &
     &              (.NOT.L_VEG_FRACS.OR..NOT.L_TRIFFID) )THEN
              LMASK=.FALSE.
            ELSE IF((N12 == 2) .AND. (.NOT.L_LSPICE)) THEN
              LMASK=.FALSE.   ! QCF is primary 12 (prognostic)
            ELSE IF ((N13 == 1).AND.(L_3D_CCA)) THEN
              LMASK=.FALSE.    !  CCA is 2D
            ELSE IF ((N13 == 2).AND.(.NOT.L_3D_CCA)) THEN
              LMASK=.FALSE.    !  CCA is 3D
            ELSE IF((N14 == 1).AND.(.NOT.L_SOOT))  THEN
              LMASK=.FALSE.
            ELSE IF((N14 == 2).AND.((.NOT.L_SOOT)                       &
     &                        .OR.  (.NOT.L_SOOT_SUREM)) )  THEN
              LMASK=.FALSE.
            ELSE IF((N14 == 3).AND.((.NOT.L_SOOT)                       &
     &                        .OR.  (.NOT.L_SOOT_HILEM)) )  THEN
              LMASK=.FALSE.
            ELSE IF ((N16 == 1).AND.(.NOT.L_CO2_INTERACTIVE)) THEN
              LMASK=.FALSE.    !  Carbon cycle switched off
            ELSE IF((N21 == 1).AND.(.NOT.L_BIOMASS))  THEN
              LMASK=.FALSE.
            ELSE IF((N21 == 2).AND.((.NOT.L_BIOMASS)                    &
     &                        .OR.  (.NOT.L_BMASS_SUREM)) )  THEN
              LMASK=.FALSE.
            ELSE IF((N21 == 3).AND.((.NOT.L_BIOMASS)                    &
     &                        .OR.  (.NOT.L_BMASS_HILEM)) )  THEN
              LMASK=.FALSE.
            ELSE IF ((N24 == 1).AND.(.NOT.L_DUST)) THEN
              LMASK=.FALSE.    !  mineral dust switched off
            ELSE IF ((N25 == 1).AND. (.NOT.L_USE_STOCHEM_CH4)) THEN
              LMASK=.FALSE.   ! STOCHEM methane switched off
            ELSE IF ((N25 == 2).AND. (.NOT.L_USE_STOCHEM_O3)) THEN
              LMASK=.FALSE.   ! STOCHEM ozone switched off
            END IF
          END IF
        ENDIF
      ELSE IF(ISEC  ==  26) THEN  !  River Routing Diagnostics
        N22 = MOD((IOPN(5)/10),10)
        IF ((N22 == 1).AND.(.NOT.L_RIVERS)) THEN
          LMASK=.FALSE.    !  River routing switched off
          ELSE IF ((N22 == 2).AND.(.NOT.L_INLAND)) THEN
                LMASK=.FALSE.    !Inland re-routing switched off
      ENDIF
!*********************************************************************
      END IF
!End of atmos diagnostic block
      END IF

      IF(LMASK) THEN
        LADRES=.TRUE.
        IF((ISPACE == 3).OR.(ISPACE == 5))                              &
     &    LMASK=.FALSE.
      ELSE
        LADRES=.FALSE.
      END IF
!
 9999 RETURN
      END SUBROUTINE TSTMSK
