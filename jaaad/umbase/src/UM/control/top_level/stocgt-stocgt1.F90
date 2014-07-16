#if defined(CONTROL) && defined(OCEAN)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine STOCGT   -----------------------------------------------
!LL
!LL Purpose: To extract 3D primary data and 2D fields (after removal of
!LL          redundant columns when cyclic boundary conditions are in
!LL          use) for use by STASH
!LL
!LL Author: N.K.Taylor      Date: 25/01/91
!LL
!LL Tested under compiler: cft77
!LL Tested under OS version: UNICOS 5.1.10
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.2  02/07/93  Dynamic allocation changes - R.T.H.Barnes.
!LL  3.3  22/11/93  Remove hard coded validation of sub-model no.
!LL                 Replace with check using parameters. - R.S.R.Hill
!LL  3.5  19/05/95  Remove ARGSIZE, pass in equivalent args. K Rogers
!LL  4.1  22/03/96  Increased no. of dimensions in array SI
!LL                 to include internal model id. G Henderson
!     4.1    Apr. 96  Rationalise *CALLs     S.J.Swarbrick
!LL  4.2  26/11/96  Allow uncompressed ocean dumps
!    4.3  18/04/97  Include calls to OD12SLAB for mpp code.
!    5.0  21/05/99  STASH sizes now held separately. R.Rawlins
!LL  5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!LL  5.3  13/09/01  Allow ocean tracer pointers to work for
!LL                 non-consceutive tracers                 S.Spall
!LL  5.4  22/10/02  Change CMESSAGE size to 80 for consistency
!LL                 with calling and called routines. R. Hill
!    5.5  20/02/03  Removed comments from #includes         P.Dando
!    6.2  31/05/06  Remove MPP and fix externals for FCM. P.Selwood
!LL
!LL Programming Standard: UM Doc Paper 3, version 2 (10/08/90)
!LL
!LL Logical components covered: CO
!LL
!LL Project Task:
!LL
!LL External documentation:
!LL
!LLEND -----------------------------------------------------------------
!*L Arguments

      SUBROUTINE STOCGT (                                               &
#include "argppx.h"
     &     ROWS_ARG, ROW_LEN, LEVELS, im_ident, submodel,               &
     &DATA,STNO,SECTION,LEVEL1,LEVEL2,base_level,                       &
     &     RMDI,VALUES,IX,IY,IZ,NT_DIM,                                 &
     &     SI,JOC_TRACER,JOC_U,JOC_V,O_CFI1,                            &
     &     O_CFI2,O_CFI3,JOC_NO_SEAPTS,JOC_NO_SEGS,ICODE,CMESSAGE)

      IMPLICIT NONE

!L  Common Blocks

#include "parparm.h"
#include "typsize.h"
! For STASH sizes
#include "typstsz.h"
#include "csubmodl.h"
#include "cppxref.h"
! Contains *CALL VERSION
#include "ppxlook.h"

      INTEGER                                                           &
     &       ROWS_ARG,                                                  &
                               ! (IN) No of rows
     &       ROW_LEN,                                                   &
                           ! (IN) No of points per row
     &       LEVELS,                                                    &
                           ! (IN) No of levels
     &       im_ident,                                                  &
                           ! (IN) Internal model id
     &       submodel,                                                  &
                           ! (IN) Submodel id
     &       STNO ,                                                     &
                           !(IN) stash identifier for extracted variable
     &       SECTION,                                                   &
                           ! (IN)
     &       LEVEL1,LEVEL2                                              &
                           ! (IN) levels between which data is extracted
     &      ,base_level                                                 &
                           ! (IN) first ocean level diagnosed
     &      ,NT_DIM                                                     &
                           ! (IN) number of tracers
     &      ,SI(NITEMS,0:NSECTS,N_INTERNAL_MODEL)                       &
                                                  !STASH IN ADDRESS
     &      ,JOC_TRACER(NT_DIM,2)                                       &
                                         ! (IN) ocean tracer pointers
     &      ,JOC_U(2),JOC_V(2)                                          &
                                         ! (IN) ocean pointers
     &      ,JOC_NO_SEAPTS,JOC_NO_SEGS                                  &
     &      ,O_CFI1(O_LEN_CFI1+1)                                       &
                                    ! (IN) ocean compressed field index
     &      ,O_CFI2(O_LEN_CFI2+1),O_CFI3(O_LEN_CFI3+1)
!
      INTEGER                                                           &
     &       IX,                                                        &
                           ! (OUT)
     &       IY,                                                        &
                           ! (OUT)
     &       IZ,                                                        &
                           ! (OUT)
     &       ICODE         ! (OUT)

      REAL                                                              &
     &    DATA(*),                                                      &
                          ! (IN) input1-d array of data, possibly compr
     &    RMDI,                                                         &
                          ! (IN) missing data indicator
     &    VALUES(*) ! (OUT) output array of data for diagnostic

      CHARACTER*80                                                      &
     &       CMESSAGE     ! (OUT) Error message if return code >0

!L  External subroutines called

!*--------------------------------------------------------------------
!
! Local variables
!
      INTEGER                                                           &
     &        ROW1,                                                     &
                        ! 1st row to be extracted from compressed data
     &        ROW2,                                                     &
                        ! Last  "   "  "      "      "       "       "
     &        PT,                                                       &
                        ! Local pointer to correct position in array D1
     &        INDEX,                                                    &
                        ! Local pointer to position in output array valu
     &        JJ,KK,KK1,                                                &
                        ! Temporary inner loop counters
     &        I,J,K,                                                    &
                        ! Do loop indices
     &       im_index,                                                  &
                           ! Internal model index number
     &        GR,                                                       &
                        ! Local value of PP cross reference grid code
     &        IXDATA                                                    &
                        ! No of x-points in input data array
     &         ,EXPPXI                                                  &
                          ! Function to extract ppxref info
     &         ,COX_TRACER  ! Index for JOC_TRACER
!
!L--------------------------------------------------------------------
! Get internal model index
      im_index = internal_model_index(im_ident)
!L   Check input values for validity
!
!    Check levels
!
       IF (LEVEL1 >  LEVEL2) THEN
          ICODE=3
          CMESSAGE='STOCGT  : Incorrect levels: level1 > level2'
          GOTO 999
       ENDIF
!
       IF (LEVEL1 <  1) THEN
          ICODE=3
          CMESSAGE='STOCGT  : Incorrect levels: level1 < 1'
          GOTO 999
       ENDIF
!
       IF (LEVEL1 >  LEVELS) THEN
          ICODE=3
          CMESSAGE='STOCGT  : Incorrect levels: level1 > LEVELS'
          GOTO 999
       ENDIF
!
       IF (LEVEL2 <  1) THEN
          ICODE=3
          CMESSAGE='STOCGT  : Incorrect levels: level2 < 1'
          GOTO 999
       ENDIF
!
       IF (LEVEL2 >  LEVELS) THEN
          ICODE=3
          CMESSAGE='STOCGT  : Incorrect levels: level2 > LEVELS'
          GOTO 999
       ENDIF
!
!    Check STASH identifier
!
       IF ( im_ident  /= OCEAN_IM) THEN
          ICODE=1
          CMESSAGE='STOCGT : Invalid SUB-MODEL'
          GOTO 999
       ENDIF
!L----------------------------------------------------------------------
!L   Calculate IX,IY,IZ - the output field dimensions, and
!L             IXDATA   - the input field x-dimension
!
!    Look up local PP cross reference
!
! DEPENDS ON: exppxi
        GR = EXPPXI(im_ident, section, stno, ppx_grid_type,             &
#include "argppx.h"
     &            icode, cmessage)
!
!     Extract all physical points - number depends on EW bdry condition
!                                   and the input grid type
!
      IF (CYCLIC_OCEAN) THEN     ! CYCLIC_OCEAN is set by user interface
        IF ((GR == ppx_ocn_tfield).OR.(GR == ppx_ocn_ufield)            &
     &       .OR.(GR == ppx_ocn_tmerid).OR.(GR == ppx_ocn_umerid)) THEN
          IXDATA=IMTM2
          IX=IMTM2
        ELSEIF ((GR == ppx_ocn_tzonal).OR.(GR == ppx_ocn_uzonal)        &
     &      .OR.(GR == ppx_ocn_scalar)) THEN
          IXDATA=1
          IX=1
        ELSE
          IXDATA=IMTM2 +2
          IX=IMTM2
        ENDIF
      ELSE
        IF ((GR == ppx_ocn_tzonal).OR.(GR == ppx_ocn_uzonal)            &
     &      .OR.(GR == ppx_ocn_scalar)) THEN
          IXDATA=1
          IX=1
        ELSE
          IXDATA=ROW_LEN
          IX=ROW_LEN
        ENDIF
      ENDIF
!
!    Remove row ROWS at U,V points in output field
!
      IF    ((GR == ppx_ocn_uall)                                       &
     &   .OR.(GR == ppx_ocn_ucomp)                                      &
     &   .OR.(GR == ppx_ocn_ufield)                                     &
     &   .OR.(GR == ppx_ocn_uzonal)) THEN
        IY = ROWS_ARG
      ELSEIF ((GR == ppx_ocn_tmerid).OR.(GR == ppx_ocn_umerid)          &
     &         .OR.(GR == ppx_ocn_scalar)) THEN
           IY=1
      ELSE
        IY = ROWS_ARG
      ENDIF
!
      IZ=LEVEL2-LEVEL1+1
!
!L----------------------------------------------------------------------
!L   Set pointer to stash address

      IF (STNO >= 101.AND.STNO <= 122) THEN

!    For all dual time level variables, except climate mean sections,
!    select 'update' level (joc_variable(2)).

        IF (SECTION >= 41.AND.SECTION <= 44) THEN
          PT = SI(STNO,SECTION,im_index)
        ELSE
          IF (STNO >= 101.AND.STNO <= 102) THEN   ! T+S tracers
             PT = JOC_TRACER(STNO-100,2)
          ELSE IF (STNO >= 103.AND.STNO <= 120) THEN   ! other tracers
          COX_TRACER=2
          DO I = 103,STNO
             IF (SI(I,0,im_index) >  1) THEN
                COX_TRACER = COX_TRACER + 1
             ENDIF
          ENDDO
          PT=JOC_TRACER(COX_TRACER,2)
          ELSE IF (STNO == 121) THEN         ! zonal velocity
            PT = JOC_U(2)
          ELSE IF (STNO == 122) THEN         ! meridional velocity
            PT = JOC_V(2)
          ENDIF
        ENDIF

      ELSE

        PT = SI(STNO,SECTION,im_index)

      ENDIF

!L----------------------------------------------------------------------
!L   Extract data
!
!
      IF (GR == ppx_ocn_tcomp) THEN
!    Extract tracers.
!    Need to access 3d primary data, using UNPACK to expand from
!    compressed form
!    Extract all rows
         ROW1=1
         ROW2=IY
! DEPENDS ON: od12slab
           CALL OD12SLAB(1,IMT,ROW1,ROW2,1,KM,IMT,1,KM,IMT,JMT,KM,      &
     &        DATA(PT),VALUES)
      ELSE IF (GR == ppx_ocn_ucomp) THEN
!    Extract velocities
!    Extract all rows
         ROW1=1
         ROW2=IY
!
! DEPENDS ON: od12slab
           CALL OD12SLAB(1,IMT,ROW1,ROW2,1,KM,IMT,1,KM,IMT,JMT,KM,      &
     &        DATA(PT),VALUES)
      ELSE IF((GR == ppx_ocn_tall).OR.(GR == ppx_ocn_uall).OR.          &
     &        (GR == ppx_ocn_tfield).OR.(GR == ppx_ocn_ufield) .OR.     &
     &        (GR == ppx_ocn_tzonal).OR.(GR == ppx_ocn_uzonal) .OR.     &
     &        (GR == ppx_ocn_tmerid).OR.(GR == ppx_ocn_umerid) .OR.     &
     &        (GR == ppx_ocn_scalar)) THEN
!     Extract uncompressed data, excluding cyclic points if present
!     in input field
         DO K=1,IZ
            KK=(K-1)*IX*IY
            KK1=(K-1)*IXDATA*IY+(LEVEL1-base_level)*IXDATA*IY
            DO J=1,IY
               JJ=(J-1)*IXDATA
               INDEX=KK+(J-1)*IX
               DO I=1,IX
                  INDEX=INDEX+1
                  VALUES(INDEX)=DATA(PT+(I-1)+JJ+KK1)
               ENDDO
            ENDDO
         ENDDO
      ELSE
         ICODE=1
         CMESSAGE='STOCGT  : Unrecognised grid type'
         GOTO 999
      ENDIF
!
      ICODE=0
      CMESSAGE='Normal Exit'
!
  999 CONTINUE
      RETURN
      END SUBROUTINE STOCGT
#endif
