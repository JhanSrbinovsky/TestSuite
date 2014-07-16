#if defined(C84_1A) || defined(FLDOP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE PPHEAD------------------------------------------
!LL
!LL  Creates a 64 word PP header from the the following:-
!LL  1)  PP_XREF (PP cross-reference array record for this sect/item)
!LL  2)  FIXED length header
!LL  3)  INTEGER constants array
!LL  4)  REAL constants array
!LL  5)  Some input arguments
!LL
!LL  Tested under compiler CFT77
!LL  Tested under OS version 5.1
!LL
!LL T.Johns     <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL   3.2  27/05/93  Code for new real missing data indicator. (TCJ)
!LL   3.5  05/06/95  Remove PP_XREF from argument list and call
!LL                  EXPPXI instead, for submodels work. K Rogers
!LL   4.0  12/09/95  LBUSER(3) [PP_INT_HEAD(LBUSER3)] set to 0.
!LL                  LBCODE set to 31300 + 20 (for Gregorian calendar)
!LL                  or 31300 + 23 (for any other calendar type), if
!LL                  the field is a timeseries.
!LL                  BRLEV, BHRLEV, BULEV[BRSVD1] and BHULEV[BRSVD2]
!LL                  contain lower level boundary and upper level bndry
!LL                  information. Above changes agreed by the WGDUM in
!LL                  first half of 1994. Code for new LBEXP experiment
!LL                  name encoding. Also removed RUN_INDIC_OP from arg
!LL                  list as it is called from CHISTORY  (Andy Brady)
!LL  4.0  12/10/95  Set Lookup(model_code) to internal model ident. RTHB
!LL  4.1  18/04/96  RUN_ID now declared in CHISTORY.  RTHBarnes.
!LL  4.1    Apr. 96  Rationalise *CALLs  S.J.Swarbrick
!LL  4.3    14/02/97 Correct bug where ocean models can try to access
!LL                  uninitialised BKH array               P.Burton
!LL  4.5    14/05/98 Put the correct data type into PP header
!LL                                                  P.Burton
!LL  4.5  14/10/97   Set correct packing type for platform
!LL                  Author D.M. Goddard
!LL  4.5    02/09/98 Set Projection No for High Res Global. D. Robinson.
!LL  5.0    09/11/99 Take into account the reversed orientation of the
!LL                  atmos grid : realhd(2) has opposite sign from its
!LL                  value before 5.0, no more offset for fields not
!LL                  starting from the origin.  JC Thil.
!LL  5.0    30/06/99 Replace references to hybrid vertical coords
!LL                  (ak,bk) by eta values. R Rawlins
!LL  5.1    13/04/00 Hack verification time if running for less than
!LL                  1 hour. Stuart Bell
!LL  5.1    20/06/00 Set sensible eta values in pp header for new
!LL                  p level above top of model. R Rawlins
!LL  5.2  7/11/99  Enable run length encoding of ocean fieldsfiles to
!LL                compress the sequences of mdi values that represent
!LL                the land points. Ian Edmond
!LL  5.2    26/10/00 Set alternative lbvc value for hybrid height.
!LL                  R Rawlins
!LL  5.2    04/09/00 (1) Set LBLEV header to surface value for model
!LL                  level=0 case.
!LL                  (2) New definition of brlev/bhrlev/blev/bhlev/
!LL                   bulev/bhulev for definition of model levels on
!LL                   Charney-Philips vertical grid. R Rawlins
!LL  5.2    06/11/00 Correct BZX,BZY headers for sub-areas. R Rawlins
!    5.3    23/07/01 Improve error trap message for inconsistent
!                    lbvc/lvcode and replace magic nos. BLEV now
!                    set from STLEVELS routine for all special
!                    level fields, except surface. R Rawlins
!LL  5.3    03/01/02 Revise hack that allows us to obtain "T+0" fields
!LL  5.3             from the physics modules. Adam Clayton
!    5.3    05/06/01 Correct BRLEV,BHRLEV lookup items for level 1 of
!                    diagnostics on model theta levels. R Rawlins
!    5.3    29/01/02 Set LBSRCE to UM version ID xxxxyyyy,
!                    where xxxx is UM version number e.g. 0503
!                    and yyyy is the model identifier e.g. 1111 for
!                    unified model. D.M. Goddard
!    5.4    11/04/02 Set up details for variable horizontal
!                    grids.                         R. Hill
!    5.4    14/05/02 Reverse some of GRR1F503 as it adversely affects
!                    VER and Horace.  P.Selwood
!    5.4    19/08/02 Further revise "T+0" hack. Adam Clayton
!     5.5    17/01/03  Update LBCODE for variable grids
!                      to indicate to MASS that grid data
!                      is in extra data vectors. R. Hill
!     6.0   06/10/03 Cater for river grid 23. C.Bunton
!     6.0   17/11/03   Call expt_enc only once and store answer
!                          A. A. Dickinson
!     6.2   15/12/05   Correct BZY in PPHEAD for sub-area velocities.
!                      R Barnes.
!     6.2   19/11/04   Reverse ORH3F505 at request of users since the
!                      11110 LBCODE now seems to cause problems
!                      downstream. In fact LBCODE now appears
!                      virtually irrelevant to the UM directly.
!                      R. Hill
!LL
!LL  Programming standard: U M DOC  Paper NO. 4,
!LL
!LL  Logical components covered: D40
!LL
!LL  Project TASK: C4
!LL
!LL  External documentation  C4
!LL
!LLEND-------------------------------------------------------------

!
!*L  INTERFACE and ARGUMENTS:------------------------------------------
      SUBROUTINE PP_HEAD(                                               &
#include "argppx.h"
     &    im_ident,FIXHD,INTHD,REALHD,                                  &
     &    LEN_INTHD,LEN_REALHD,IE,IS,GR,                                &
     &    lfullfield,real_level,pseudo_level,                           &
     &    samples,start,start_or_verif_time,end_or_data_time,pp_len,    &
     &    extraw,PP_INT_HEAD,PP_REAL_HEAD,N_COLS_OUT,NUM_WORDS,         &
     &    LEN_BUF_WORDS,N_ROWS_OUT,NROW_IN,SROW_IN,WCOL_IN,ECOL_IN,     &
     &    lbproc_comp,                                                  &
     &    sample_prd,FCST_PRD,COMP_ACCRCY,PACKING_TYPE,                 &
     &    st_grid,IWA,zseak_rho,Ck_rho,zseak_theta,Ck_theta,            &
     &    model_levels,LevIndex,ROTATE,ELF,                             &
     &    OCEAN,OCN_DZ,OCN_KM,                                          &
     &    ICODE,CMESSAGE)
!*----------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, PARAMETER :: LEN_FIXHD = 256

      CHARACTER*(80) CMESSAGE !OUT OUT MESSAGE FROM ROUTINE
!
      LOGICAL                                                           &
     &  start                                                           &
                      ! IN flag to control update for verif/start time
     &, OCEAN                                                           &
                       !IN TRUE if processing an ocean diagnostic
     &, lfullfield     !IN TRUE if output field on full horiz domain
!
      INTEGER                                                           &
     &  start_or_verif_time(7)                                          &
                               ! IN verif time/start time for means etc
     &, end_or_data_time(7)                                             &
                               ! IN data time/end time for means etc
     &, samples                ! IN no of samples in period (timeseries)
!
      INTEGER                                                           &
     &  ICODE                                                           &
                          !IN    Return code from the routine
     &, im_ident                                                        &
                          !IN    Internal model identifier
     &, PP_LEN                                                          &
                          !IN    Length of the lookup table
     &, LEN_INTHD                                                       &
                          !IN    Length of the Integer Constants
     &, LEN_REALHD                                                      &
                          !IN    Length of the Real Constants
     &, FIXHD(LEN_FIXHD)                                                &
                          !IN    Array of Fixed Constants
     &, INTHD(LEN_INTHD)                                                &
                          !IN    Array of Integer Constants
     &, OCN_KM            !IN    number of ocean model levels
!
      INTEGER                                                           &
     &  st_grid                                                         &
                          !IN    STASH horizontal grid type
     &, model_levels                                                    &
                          !IN    No of model levels
     &, LevIndex                                                        &
                          !IN    level index
     &, N_ROWS_OUT                                                      &
                          !IN    PPHORIZ_OUT=N_ROWS_OUT*N_COLS_OUT+extra
     &, N_COLS_OUT                                                      &
                          !IN    PPHORIZ_OUT=N_COLS_OUT*N_ROWS_OUT+extra
     &, NROW_IN,SROW_IN                                                 &
                          !IN    The most nrthrly/southerly row.
     &, WCOL_IN,ECOL_IN                                                 &
                          !IN    The most westerly/easterly column
     &, pseudo_level                                                    &
                          !IN    Output PP pseudo-level
     &, COMP_ACCRCY                                                     &
                          !IN    PACKING ACCURACY IN POWER OF 2
     &, PACKING_TYPE      !IN   0 = No packing, 1 = WGDOS, 3 = GRIB
      INTEGER                                                           &
     &  NUM_WORDS                                                       &
                          !IN    Number of 64 Bit words to hold DATA
     &, extraw                                                          &
                          !IN    Number of extra-data words
     &, LEN_BUF_WORDS                                                   &
                          !IN    Number of 64 Bit words (rounded to 512)
     &, IWA                                                             &
                          !IN    Start word address.
     &, IE                                                              &
                          !IN    Item Number
     &, IS                                                              &
                          !IN    Section Number
     &, GR                                                              &
                          !IN    Grid point code
     &, LBPROC_COMP(14)                                                 &
                          !IN    Subcomponents(0/1) to make up LBPROC
     &, PP_INT_HEAD(PP_LEN)          !OUT  Integer Lookup table
!
! UM6.5 -  MODEL_ANALYSIS_HRS changed to REAL,  
!             requires FCST_PRD changed to REAL also
      REAL                                                              &
     &  FCST_PRD                                                        &
                            !IN    Forecast period
     &, PP_REAL_HEAD(PP_LEN)                                            &
                            !OUT Real Lookup table
     &, REALHD(LEN_REALHD)                                              &
                            !IN  Real header
     &, real_level                                                      &
                            !IN  Output PP level(REAL)
     &, sample_prd                                                      &
                            !IN  Sampling period in hours for time mean
     &, zseak_rho    (model_levels)                                     &
                                   !IN vert coeff zsea on rho levels
     &, Ck_rho       (model_levels)                                     &
                                   !IN vert coeff Ck   on rho levels
     &, zseak_theta(0:model_levels)                                     &
                                   !IN vert coeff zsea on theta levels
     &, Ck_theta   (0:model_levels)                                     &
                                   !IN vert coeff Ck   on theta levels
     &, OCN_DZ(OCN_KM)      !IN  ocean depths at KM levels
!
!*---------------------------------------------------------------------
#include "clookadd.h"
#include "stparam.h"
#include "csubmodl.h"
#include "cppxref.h"
! Contains *CALL VERSION
#include "ppxlook.h"
#include "c_mdi.h"
#include "chsunits.h"
#include "cntlall.h"
#include "chistory.h"
#include "parvars.h"
#include "comvgrid.h"
#include "cmaxsize.h"
#include "ctime.h"
#include "ctfilt.h"
#include "c_model_id.h"

        EXTERNAL EXPPXI
        EXTERNAL EXPT_ENC

!*L  WORKSPACE USAGE:-------------------------------------------------
!   DEFINE LOCAL WORKSPACE ARRAYS: None
!
!*---------------------------------------------------------------------
!    DEFINE LOCAL VARIABLES

      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='PP_head')
      REAL                                                              &
     &  ocn_depth                                                       &
                      !     depth of ocean at level
     &, ocn_depth_h   !     depth of ocean at half level
      INTEGER                                                           &
     &  PP_LBFC                                                         &
                      !     M08 Level code
     &, PP_LBTYP                                                        &
                      !     M08 Field type code
     &, PP_LBLEV                                                        &
                      !     M08 Field level code
     &, PP_IPROJ                                                        &
                      !     M08 Projection number
     &, PP_LBVC                                                         &
                      !     Vertical coord type
     &, II                                                              &
                      !     Local Counter
     &, int_level                                                       &
                      !     integer value of level
     &, K                                                               &
                      !     local counter
     &, IA,IB,IC                                                        &
                      !     Component codes to make up LBTIM
     &, mean_code                                                       &
                      !     spatial averaging code derived from GR
     &, lvcode                                                          &
                      !     lv code
     &, EXPPXI                                                          &
                      !     Function to extract ppxref info
     &, EXPTCODE      !     integer coded experiment name
      INTEGER :: get_um_version_id

! Local Parameters
      INTEGER, PARAMETER :: unset_exptenc = -5555

! Local scalars:
      INTEGER,SAVE  ::  expt_enc_answer=unset_exptenc

      LOGICAL                                                           &
     &  ELF,                                                            &
     &  ROTATE
!
!LL   Construct PP header
!
!  Timestamps ----------------------------------------------------------
!
!L
!L Set up time info dependent on start flag.
!L For all but time series start will be TRUE so all time information
!L will be set up from FIXHD in effect, but for time series start
!L will be set up by TEMPORAL and passed in, so that dump headers are
!L set correctly for such fields.
!L Note: end_or_data_time will be updated from current model time in
!L       FIXHD(28-34) for time means/accumulations etc.
!L
      IF (start) THEN    ! start timestep so update start time
        PP_INT_HEAD(LBYR)=start_or_verif_time(1)
        PP_INT_HEAD(LBMON)=start_or_verif_time(2)
        PP_INT_HEAD(LBDAT)=start_or_verif_time(3)
        PP_INT_HEAD(LBHR)=start_or_verif_time(4)
        PP_INT_HEAD(LBMIN)=start_or_verif_time(5)
        PP_INT_HEAD(LBDAY)=start_or_verif_time(7)
! If the IAU scheme was used on the previous timestep to add a complete
! increment at the nominal analysis time, reset verification minute to
! zero for fields not in sections 0, 15 or 16. Required so we can obtain
! "analysis" fields from the physics modules.
! UM6.5 - MODEL_ANALYSIS_HRS replaced by MODEL_ANALYSIS_MINS
        IF (L_IAU .AND. im_ident == A_IM) THEN
          IF ( IAU_StartMin == IAU_EndMin              .AND.            &
     &         IAU_StartMin == MODEL_ANALYSIS_MINS     .AND.            &
     &         STEPIM(A_IM) == IAU_DTResetStep + 1     .AND.            &
     &         IS /= 0                                 .AND.            &
     &         IS /= 15                                .AND.            &
     &         IS /= 16 ) THEN
            PP_INT_HEAD(LBMIN) = 0
          END IF
        END IF
      END IF
      PP_INT_HEAD(LBYRD)=end_or_data_time(1)
      PP_INT_HEAD(LBMOND)=end_or_data_time(2)
      PP_INT_HEAD(LBDATD)=end_or_data_time(3)
      PP_INT_HEAD(LBHRD)=end_or_data_time(4)
      PP_INT_HEAD(LBMIND)=end_or_data_time(5)
      PP_INT_HEAD(LBDAYD)=end_or_data_time(7)
!
!  Secondary time information ------------------------------------------
!
! LBTIM is 100*IA+10*IB+IC - this encodes the time processing type
!
      IA=INT(sample_prd)           ! Sampling period in whole hours
      IF(sample_prd == 0.0) THEN   ! NB: may be a fraction of an hour
        IB=1                       ! Forecast field
      ELSE
        IF (IA == 0) THEN
          IA=1                     ! 0 < sample_prd < 1 counts as 1 hour
        ENDIF
        IB=2                       ! Time mean or accumulation
      ENDIF
      IC=FIXHD(8)                  ! Calendar (1: Gregorian, 2: 360 day)
!
      PP_INT_HEAD(LBTIM)=100*IA+10*IB+IC
      PP_INT_HEAD(LBFT)=FCST_PRD
!
!  Data length ---------------------------------------------------------
!
      PP_INT_HEAD(LBLREC)=NUM_WORDS
!
!  Grid code (determined from dump fixed-length header) ----------------
!
      IF (samples == 0) THEN
!       Field is not a timeseries
        IF(FIXHD(4) <  100) THEN
          IF (VAR_GRID_TYPE == 0) THEN
             PP_INT_HEAD(LBCODE)=1   ! Regular lat/long grid
          ELSE
             ! Set variable  grid code to the same as regular
             ! grid which effectively renders this peice of
             ! information meaningless
             PP_INT_HEAD(LBCODE) = 1   ! Variable grid
          ENDIF
        ELSE
          IF (VAR_GRID_TYPE == 0) THEN
             PP_INT_HEAD(LBCODE)=101   ! Regular lat/long grid
          ELSE
             ! Set variable  grid code to the same as regular
             ! grid which effectively renders this peice of
             ! information meaningless
             PP_INT_HEAD(LBCODE) = 101 ! Variable grid
          ENDIF

        ENDIF
      ELSE
!       Field is a timeseries
        PP_INT_HEAD(LBCODE)=31300
        IF (FIXHD(8) == 1) THEN
!         Calendar --  1: Gregorian
          PP_INT_HEAD(LBCODE)=PP_INT_HEAD(LBCODE)+20
        ELSEIF (FIXHD(8) == 2) THEN
!         Calendar -- 360 day (Model Calendar)
          PP_INT_HEAD(LBCODE)=PP_INT_HEAD(LBCODE)+23
        ELSE
!         Unknown calendar. Fail.
          ICODE=2
      CMESSAGE='PPHEAD: unknown calender type in fixhd(8)'
        ENDIF
      ENDIF
!
!  Hemispheric subregion indicator -------------------------------------
!
      IF (samples >  0 .OR. .NOT.lfullfield) THEN
!  Field is a timeseries/trajectory or subdomain of the full model area
        PP_INT_HEAD(LBHEM)=3
      ELSEIF (FIXHD(4) <  100) THEN
!  Otherwise, use the value for the full model area encoded in the dump
        PP_INT_HEAD(LBHEM)=FIXHD(4)
      ELSE
        PP_INT_HEAD(LBHEM)=FIXHD(4)-100
      ENDIF
!
!  Field dimensions (rows x cols) --------------------------------------
!
      PP_INT_HEAD(LBROW)=N_ROWS_OUT
      PP_INT_HEAD(LBNPT)=N_COLS_OUT
!
!  'Extra data' length (now accomodates timeseries sampling data) ------
!
      PP_INT_HEAD(LBEXT)=extraw
!
!  Packing method indicator (new definition introduced at vn2.8)--------
#if defined(CRAY) && defined(T3E)
       IF(PACKING_TYPE == 1)THEN    ! WGDOS packing
         PP_INT_HEAD(LBPACK)=03001
       ELSEIF(PACKING_TYPE == 4)THEN ! Run length encoding
         PP_INT_HEAD(LBPACK)=03004
       ELSEIF(PACKING_TYPE == 3)THEN ! GRIB packing
         PP_INT_HEAD(LBPACK)=04003
       ELSEIF(PACKING_TYPE == 0)THEN ! No packing
         PP_INT_HEAD(LBPACK)=03000
       ELSE
         ICODE=1
         CMESSAGE='PPHEAD  Packing type undefined'
         PP_INT_HEAD(LBPACK)=03000
      ENDIF
#endif
#if defined(CRAY) && !defined(T3E)
       IF(PACKING_TYPE == 1)THEN    ! WGDOS packing
         PP_INT_HEAD(LBPACK)=02001
       ELSEIF(PACKING_TYPE == 4)THEN ! Run length encoding
         PP_INT_HEAD(LBPACK)=02004
       ELSEIF(PACKING_TYPE == 3)THEN ! GRIB packing
         PP_INT_HEAD(LBPACK)=04003
       ELSEIF(PACKING_TYPE == 0)THEN ! No packing
         PP_INT_HEAD(LBPACK)=02000
       ELSE
         ICODE=1
         CMESSAGE='PPHEAD  Packing type undefined'
         PP_INT_HEAD(LBPACK)=02000
      ENDIF
#endif
#if !defined(CRAY)
       IF(PACKING_TYPE == 1)THEN    ! WGDOS packing
         PP_INT_HEAD(LBPACK)=00001
       ELSEIF(PACKING_TYPE == 4)THEN ! Run length encoding
         PP_INT_HEAD(LBPACK)=00004
       ELSEIF(PACKING_TYPE == 3)THEN ! GRIB packing
         PP_INT_HEAD(LBPACK)=00003
       ELSEIF(PACKING_TYPE == 0)THEN ! No packing
         PP_INT_HEAD(LBPACK)=00000
       ELSE
         ICODE=1
         CMESSAGE='PPHEAD  Packing type undefined'
         PP_INT_HEAD(LBPACK)=00000
      ENDIF
#endif
!
!  PP header release no ------------------------------------------------
!
      PP_INT_HEAD(LBREL)=2
!
!  Primary fieldcode (some hardwiring for ELF winds) -------------------
!  Secondary fieldcode not used currently
!
! DEPENDS ON: exppxi
      PP_LBFC=EXPPXI(im_ident, is, ie, ppx_field_code,                  &
#include "argppx.h"
     &               icode, cmessage)
      IF(ELF.AND..NOT.ROTATE) THEN  ! ELF winds are in x,y direction
        IF(PP_LBFC == 56) PP_LBFC=48
        IF(PP_LBFC == 57) PP_LBFC=49
      ENDIF
      PP_INT_HEAD(LBFC)=PP_LBFC
      PP_INT_HEAD(LBCFC)=0
!
!  Processing code (encodes several things in one field) ---------------
!
      PP_INT_HEAD(LBPROC)=0
      DO II=14,1,-1
        PP_INT_HEAD(LBPROC)=PP_INT_HEAD(LBPROC)*2+LBPROC_COMP(II)
      ENDDO
!
!  Vertical coordinate type --------------------------------------------
!  Vertical coordinate type for reference level not coded
!
! DEPENDS ON: exppxi
      PP_LBVC=EXPPXI(im_ident, is, ie, ppx_lbvc_code,                   &
#include "argppx.h"
     &               icode, cmessage)
      PP_INT_HEAD(LBVC)=PP_LBVC
      PP_INT_HEAD(LBRVC)=0

! [Note that, although most diagnostics are defined over int_level=
! (1:model_levels), and hence int_level=LevIndex, some variables are
! defined over int_level=(0:model_levels). For this special case,
! LevIndex is (1:model_levels+1), since LevIndex is defined starting at
! 1, and int_level != LevIndex.]
      int_level = real_level+0.00001  ! ensure no rounding problems
!
!  Experiment number coded from EXPT_ID and JOB_ID for non
!  operational set to RUN_INDIC_OP for operational use.
!
      IF (MODEL_STATUS /= 'Operational') THEN
        RUN_ID(1:4)=EXPT_ID
        RUN_ID(5:5)=JOB_ID
        IF (expt_enc_answer == unset_exptenc) THEN
!  Function EXPT_ENC will encode the run_id into a unique integer
! DEPENDS ON: expt_enc
           CALL EXPT_ENC(RUN_ID,EXPTCODE,ICODE,CMESSAGE)
           expt_enc_answer=EXPTCODE

        ELSE
          EXPTCODE=expt_enc_answer
        ENDIF


        PP_INT_HEAD(LBEXP)=EXPTCODE          ! LBEXP
      ELSE
        PP_INT_HEAD(LBEXP)=RUN_INDIC_OP      ! LBEXP (ITAB)
      ENDIF
!
!  Direct access dataset start address and no of records ---------------
!
      PP_INT_HEAD(LBEGIN)=IWA
      PP_INT_HEAD(LBNREC)=LEN_BUF_WORDS
!
!  Operational fieldsfile projection no, fieldtype + level codes -------
!  These are hardwired according to model resolution
!
      IF(INTHD(6) == 192) THEN
        PP_IPROJ=802
      ELSE IF(INTHD(6) == 288) THEN
        PP_IPROJ=800
      ELSE IF(INTHD(6) == 96) THEN
        PP_IPROJ=870
       ELSE IF(INTHD(6) == 432) THEN
         PP_IPROJ=800
      ELSE
        PP_IPROJ=900
      ENDIF
! DEPENDS ON: exppxi
      PP_LBTYP=EXPPXI(im_ident, is, ie, ppx_meto8_fieldcode,            &
#include "argppx.h"
     &               icode, cmessage)
! DEPENDS ON: exppxi
      lvcode=EXPPXI(im_ident, is, ie, ppx_lv_code,                      &
#include "argppx.h"
     &             icode, cmessage)
      IF(real_level == -1.0) THEN
! DEPENDS ON: exppxi
        PP_LBLEV=EXPPXI(im_ident, is, ie, ppx_meto8_levelcode,          &
#include "argppx.h"
     &               icode, cmessage)   ! levelcode 9999 or 8888
      ELSE
        IF (im_ident  ==  atmos_im) THEN
          IF (lvcode == ppx_half_level .AND. int_level == 0) THEN
!             This is a surface level: reset lblev
            PP_LBLEV=ppx_meto8_surf
          ELSE
            PP_LBLEV=int_level
          ENDIF
        ELSE
          PP_LBLEV=int_level
        ENDIF
      ENDIF
      PP_INT_HEAD(LBPROJ)=PP_IPROJ
      PP_INT_HEAD(LBTYP)=PP_LBTYP
      PP_INT_HEAD(LBLEV)=PP_LBLEV
!
!  Reserved slots for future expansion ---------------------------------
!
      PP_INT_HEAD(LBRSVD1)=0
      PP_INT_HEAD(LBRSVD2)=0
      PP_INT_HEAD(LBRSVD3)=0
      PP_INT_HEAD(LBRSVD4)=0
!
! Generate model version_id
! DEPENDS ON: get_um_version_id
      PP_INT_HEAD(LBSRCE)=get_um_version_id(model_id)
!
! Data type - extract from PPXREF
! DEPENDS ON: exppxi
      PP_INT_HEAD(DATA_TYPE)=EXPPXI(im_ident, is, ie, ppx_data_type,    &
#include "argppx.h"
     &                              icode, cmessage)
!
!  Address within dump or PP file --------------------------------------
!
      PP_INT_HEAD(NADDR)=IWA
!
!  LBUSER3 is not currently used (ie set to 0).
!
      PP_INT_HEAD(LBUSER3)=0
!
!  STASH section/item code ---------------------------------------------
!
      PP_INT_HEAD(ITEM_CODE)=IS*1000+IE
!
!  STASH pseudo-level (for fields which have pseudo-levels defined) ----
!
      PP_INT_HEAD(LBPLEV)=pseudo_level
!
!  Spare for user's use ------------------------------------------------
!
      PP_INT_HEAD(LBUSER6)=0
      PP_INT_HEAD(MODEL_CODE) = im_ident
!
!  Reserved for future PP package use ----------------------------------
!
      PP_REAL_HEAD(BRSVD3)=0.0
      PP_REAL_HEAD(BRSVD4)=0.0
      PP_REAL_HEAD(BDATUM)=0.0
      PP_REAL_HEAD(BACC)=COMP_ACCRCY ! packing accuracy stored as real
!
!  Vertical grid description -------------------------------------------
!  Level and reference level
!
      IF(PP_LBVC >= 126.AND.PP_LBVC <= 139) THEN ! Special codes
!                                                  (surf botttom,
!                                                   top all zero)
        PP_REAL_HEAD(BLEV)=0.0
        PP_REAL_HEAD(BHLEV)=0.0
        PP_REAL_HEAD(BRLEV)=0.0
        PP_REAL_HEAD(BHRLEV)=0.0
        PP_REAL_HEAD(BULEV)=0.0
        PP_REAL_HEAD(BHULEV)=0.0
      ELSEIF(PP_LBVC == 9.OR.PP_LBVC == 65) THEN ! Hybrid/ETA levels
! DEPENDS ON: exppxi
        lvcode=EXPPXI(im_ident, is, ie, ppx_lv_code,                    &
#include "argppx.h"
     &    icode, cmessage)

! From vn5.2:
! height of model level k above mean sea level is
!       z(i,j,k) = zsea(k) + C(k)*zorog(i,j)
! bulev,bhulev      zsea,C of upper layer boundary
! blev ,bhlev       zsea,C of level
! brlev,bhrlev      zsea,C of lower level boundary
! The level here can refer to either a theta or rho level, with
! layer boundaries defined by surrounding rho or theta levels.
!
! [Note assumption that the top of the atmosphere model,
!  ie eta_theta_levels(model_levels) = 1.0, is equidistant from a
!  bounding p level above (not represented explicitly) and the p
!  level below (at eta_rho_levels(model_levels) ).]

        IF (lvcode == ppx_half_level) THEN ! theta level (& w)

          IF(int_level >= model_levels) THEN ! top level
            PP_REAL_HEAD(bulev) =  zseak_theta(model_levels) * 2.0      &
     &                                   - zseak_rho(model_levels)
            PP_REAL_HEAD(bhulev)=     Ck_theta(model_levels) * 2.0      &
     &                                   -    Ck_rho(model_levels)
          ELSE
            PP_REAL_HEAD(bulev) = zseak_rho(int_level+1)
            PP_REAL_HEAD(bhulev)=    Ck_rho(int_level+1)
          ENDIF                             ! top level

          PP_REAL_HEAD(blev) = zseak_theta(int_level)
          PP_REAL_HEAD(bhlev)=    Ck_theta(int_level)

! Note that the lowest theta level (=1) has a lower layer boundary at
! the surface, as set explicitly in the model interface to physics.
          IF(int_level <= 1) THEN            ! bottom level
             PP_REAL_HEAD(brlev) = 0.     ! zsea at/below surface
             PP_REAL_HEAD(bhrlev)= 1.     ! C    at/below surface
          ELSE
            PP_REAL_HEAD(brlev) = zseak_rho(int_level)
            PP_REAL_HEAD(bhrlev)=    Ck_rho(int_level)
          ENDIF                              ! bottom level

        ELSEIF(lvcode == ppx_full_level) THEN ! rho level (& u,v,p)

          IF(int_level >  model_levels) THEN ! p above top level
            PP_REAL_HEAD(bulev) = zseak_theta(model_levels) * 2.0       &
     &                                  - zseak_rho(model_levels)
            PP_REAL_HEAD(bhulev)=    Ck_theta(model_levels) * 2.0       &
     &                                   -   Ck_rho(model_levels)
            PP_REAL_HEAD(blev) = PP_REAL_HEAD(bulev)
           PP_REAL_HEAD(bhlev)= PP_REAL_HEAD(bhulev)
          ELSE
            PP_REAL_HEAD(bulev) = zseak_theta(int_level)
            PP_REAL_HEAD(bhulev)=    Ck_theta(int_level)
            PP_REAL_HEAD(blev)  = zseak_rho(int_level)
            PP_REAL_HEAD(bhlev) =    Ck_rho(int_level)
          ENDIF                              ! p above top level

          IF(int_level <= 0) THEN            ! bottom level
            PP_REAL_HEAD(brlev) = 0.    ! zsea at/below surface
            PP_REAL_HEAD(bhrlev)= 1.    ! C    at/below surface
          ELSE
            PP_REAL_HEAD(brlev) = zseak_theta(int_level-1)
            PP_REAL_HEAD(bhrlev)=    Ck_theta(int_level-1)
          ENDIF                              ! bottom level

        ELSE                ! Illegal lvcode
          ICODE=1
          CMESSAGE=' Inconsistent vertical coordinate codes in STASHmas&
     &ter for this output field: LevelT indicates that this variable is&
     & on model levels, but LBVC used for pp header label does not.'
          WRITE(6,*) RoutineName,CMESSAGE                               &
     &            ,' section,item,LevelT,pp_lbvc=',is,ie,lvcode,pp_lbvc

        ENDIF               ! Test on lvcode

      ELSEIF (PP_LBVC == 2.AND.OCEAN) THEN ! Depth levels
        PP_REAL_HEAD(BHRLEV)=0.0
        PP_REAL_HEAD(BHULEV)=0.0
! DEPENDS ON: exppxi
        lvcode=EXPPXI(im_ident, is, ie, ppx_lv_code,                    &
#include "argppx.h"
     &               icode, cmessage)
! ocn_depth defined for ocean full levels (e.g. temperature), and
! ocn_depth_h for ocean half-levels (e.g. vertical velocity)
        ocn_depth=0.5*OCN_DZ(1)
        ocn_depth_h=0.
        IF (int_level >  1) THEN
          DO K=2,int_level
!           Loop over levels calculating half levels as we go.
            ocn_depth=ocn_depth+0.5*(OCN_DZ(K-1)+OCN_DZ(K))

            ocn_depth_h=ocn_depth_h+OCN_DZ(K-1)
          END DO
        ENDIF
        IF (lvcode == ppx_half_level) THEN
          PP_REAL_HEAD(BLEV)=ocn_depth_h
          PP_REAL_HEAD(BRLEV)=ocn_depth
          IF (int_level == 1) THEN
            PP_REAL_HEAD(BULEV)=0.0 ! This level would be a
!                                     half level above the ocean.
!                                     Set to zero.
          ELSE
            PP_REAL_HEAD(BULEV)=ocn_depth_h-0.5*OCN_DZ(int_level-1)
          ENDIF
        ELSE
          PP_REAL_HEAD(BLEV)=ocn_depth
          PP_REAL_HEAD(BRLEV)=ocn_depth_h+OCN_DZ(int_level)
          PP_REAL_HEAD(BULEV)=ocn_depth_h
        ENDIF

        PP_REAL_HEAD(BHLEV)=0.0
      ELSE
        PP_REAL_HEAD(BLEV)=real_level
        PP_REAL_HEAD(BHLEV)=0.0
        PP_REAL_HEAD(BRLEV)=0.0  ! The boundary levels
        PP_REAL_HEAD(BHRLEV)=0.0 ! are not known
        PP_REAL_HEAD(BULEV)=0.0  ! for pressure
        PP_REAL_HEAD(BHULEV)=0.0 ! levels.
      ENDIF
!
!  Horizontal grid description -----------------------------------------
!  Position of pole (from dump fixed-length header)
!  Grid orientation (hardwired 0.0)
!  Origin and spacing of grid (depends on output grid type)
!
      PP_REAL_HEAD(BPLAT)=REALHD(5)
      PP_REAL_HEAD(BPLON)=REALHD(6)
      PP_REAL_HEAD(BGOR)=0.0
      IF (samples >  0) THEN   ! Indicates a timeseries/trajectory
        PP_REAL_HEAD(BZX)=0.0
        PP_REAL_HEAD(BDX)=0.0
        PP_REAL_HEAD(BZY)=0.0
        PP_REAL_HEAD(BDY)=0.0
      ELSE
        IF (OCEAN) THEN       !   set BZY,BZX,BDY,BDX for ocean
          IF (st_grid == st_uv_grid .OR. st_grid == st_zu_grid          &
     &        .OR. st_grid == st_mu_grid) THEN
            PP_REAL_HEAD(BZY)=REALHD(3)-REALHD(2)/2.0
            PP_REAL_HEAD(BZX)=REALHD(4)-REALHD(1)/2.0
          ELSEIF (st_grid == st_tp_grid .OR. st_grid == st_zt_grid      &
     &       .OR. st_grid == st_mt_grid .OR. st_grid == st_scalar) THEN
            PP_REAL_HEAD(BZY)=REALHD(3)-REALHD(2)
            PP_REAL_HEAD(BZX)=REALHD(4)-REALHD(1)
          ELSEIF (st_grid == st_cu_grid) THEN
            PP_REAL_HEAD(BZY)=REALHD(3)-REALHD(2)
            PP_REAL_HEAD(BZX)=REALHD(4)-REALHD(1)/2.0
          ELSEIF (st_grid == st_cv_grid) THEN
            PP_REAL_HEAD(BZY)=REALHD(3)-REALHD(2)/2.0
            PP_REAL_HEAD(BZX)=REALHD(4)-REALHD(1)
          ENDIF
          IF (REALHD(32) >  REALHD(29)) THEN !   greater than RMDI
            PP_REAL_HEAD(BDY)=0.0
            PP_REAL_HEAD(BDX)=REALHD(32)
          ELSE
            PP_REAL_HEAD(BDY)=REALHD(2)
            PP_REAL_HEAD(BDX)=REALHD(1)
          ENDIF
      ! If this is variable horizontal grid information
      ! then the horizontal spacings in the headers must be
      ! set to zero in the appropriate direction.
      IF (VAR_GRID_TYPE >  0) THEN
         IF (X_VAR_GRID) PP_REAL_HEAD(BDX) = 0.0
         IF (Y_VAR_GRID) PP_REAL_HEAD(BDY) = 0.0
      ENDIF
        ELSE                 !   set BZY,BZX,BDY,BDX for atmos
          IF(st_grid == st_riv_grid)THEN
            PP_REAL_HEAD(BDY) = -180.0/PP_INT_HEAD(LBROW)
            PP_REAL_HEAD(BZY) = -REALHD(3) - PP_REAL_HEAD(BDY)*0.5
            PP_REAL_HEAD(BDX) = 360.0/PP_INT_HEAD(LBNPT)
            PP_REAL_HEAD(BZX) = REALHD(4) - PP_REAL_HEAD(BDX)*0.5
          ELSE
           IF(st_grid == st_uv_grid.OR.st_grid == st_cv_grid .OR.       &
     &       st_grid == st_zu_grid.OR.st_grid == st_mu_grid) THEN
            PP_REAL_HEAD(BZY)=REALHD(3)-REALHD(2)/2.0 ! UV pts
           ELSE
            PP_REAL_HEAD(BZY)=REALHD(3)-REALHD(2) ! Zeroth Lat BZY
           ENDIF
!
           IF(st_grid == st_uv_grid.OR.st_grid == st_cu_grid .OR.       &
     &       st_grid == st_zu_grid.OR.st_grid == st_mu_grid) THEN
            PP_REAL_HEAD(BZX)=REALHD(4)-REALHD(1)/2.0 !UV points
           ELSE
            PP_REAL_HEAD(BZX)=REALHD(4)-REALHD(1) ! Zeroth Long BZX
           ENDIF
           PP_REAL_HEAD(BDX)=REALHD(1) ! Long intvl BDX
           PP_REAL_HEAD(BDY)=REALHD(2) ! Lat intvl BDY
          ENDIF
        ENDIF
!
! Add on offset for fields not starting from the origin (sub-areas)
!
      ! If this is variable horizontal grid information
      ! then the horizontal start points are NOT equally spaced.
      ! The existing code is meaningless in this case so
      ! we must get our actual start point some other way!
      IF (VAR_GRID_TYPE >  0) THEN
         IF (X_VAR_GRID) PP_REAL_HEAD(BZX) =                            &
     &                       X_BOUNDARY(WCOL_IN,VAR_GRID_TYPE)
         IF (Y_VAR_GRID) PP_REAL_HEAD(BZY) =                            &
     &                       Y_BOUNDARY(SROW_IN,VAR_GRID_TYPE)
      ELSE
          IF(st_grid == st_uv_grid.OR.st_grid == st_cv_grid .OR.        &
     &       st_grid == st_zu_grid.OR.st_grid == st_mu_grid) THEN
            IF(SROW_IN /= 1) THEN
! v area shifted up a row so diagnostic same as pre new dynamics
! when the rows were reversed. ie. N-S
              PP_REAL_HEAD(BZY)=PP_REAL_HEAD(BZY)                       &
     &                         +(SROW_IN-2)*PP_REAL_HEAD(BDY)
             END IF
          ELSE
            PP_REAL_HEAD(BZY)=PP_REAL_HEAD(BZY)                         &
     &                       +(SROW_IN-1)*PP_REAL_HEAD(BDY)
          END IF
        PP_REAL_HEAD(BZX)=PP_REAL_HEAD(BZX)                             &
     &                    +(WCOL_IN-1)*PP_REAL_HEAD(BDX)

      ENDIF
        IF(PP_REAL_HEAD(BZX) >= 360.0) THEN
           PP_REAL_HEAD(BZX)=PP_REAL_HEAD(BZX)-360.0
        ENDIF

!
! If horizontal averaging has been applied to the output field,
! set BDX and/or BDY to the full (sub)domain extent which was processed.
! If the input field was intrinsically non-2D (eg. zonal), assume that
! the collapsed dimension(s) covered the full model domain.
!
        mean_code=(GR/block_size)*block_size
        IF (st_grid == st_zt_grid .OR. st_grid == st_zu_grid            &
     &      .OR. st_grid == st_scalar) THEN
          PP_REAL_HEAD(BDX)=REAL(INTHD(6))*PP_REAL_HEAD(BDX)
        ELSEIF (mean_code == zonal_mean_base .OR.                       &
     &      mean_code == field_mean_base .OR.                           &
     &      mean_code == global_mean_base) THEN
          PP_REAL_HEAD(BDX)=ABS(REAL(ECOL_IN-WCOL_IN))*PP_REAL_HEAD(BDX)
        ENDIF
!
        IF (st_grid == st_mt_grid .OR. st_grid == st_mu_grid            &
     &      .OR. st_grid == st_scalar) THEN
          PP_REAL_HEAD(BDY)=REAL(INTHD(7))*PP_REAL_HEAD(BDY)
        ELSEIF (mean_code == merid_mean_base .OR.                       &
     &      mean_code == field_mean_base .OR.                           &
     &      mean_code == global_mean_base) THEN
          PP_REAL_HEAD(BDY)=ABS(REAL(NROW_IN-SROW_IN))*PP_REAL_HEAD(BDY)

        ENDIF
      ENDIF
      IF(.NOT. OCEAN .AND. (FIXHD(117) > 1 .AND. FIXHD(122) > 1)) THEN 
        PP_REAL_HEAD(BDX) = RMDI
        PP_REAL_HEAD(BDY) = RMDI
        PP_REAL_HEAD(BZX) = RMDI
        PP_REAL_HEAD(BZY) = RMDI
      ENDIF
!
! Missing data indicator (from PARAMETER) ------------------------------
! MKS scaling factor (unity as model uses SI units throughout)
!
      PP_REAL_HEAD(BMDI)=RMDI
      PP_REAL_HEAD(BMKS)=1.0
!

      RETURN
      END SUBROUTINE PP_HEAD
#endif
