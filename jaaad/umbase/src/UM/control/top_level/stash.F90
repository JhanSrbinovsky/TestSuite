#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine STASH --------------------------------------------------
!LL
!LL Purpose: Control routine for diagnostic processing step-by-step.
!LL          Called after each code section to process diagnostic fields
!LL          from workspace STASH_WORK to their final destination in D1
!LL          or PP file.  This routine loops over raw input fields and
!LL          calls a service routine STWORK to do the actual processing.
!LL
!LL TJ, SI      <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.1    2/02/93 : Add NUNITS to argument list of STWORK to increase
!LL                   no. of i/o units for 'C' portable code.
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL  3.2  13/04/93  Dynamic allocation of main arrays. R T H Barnes.
!LL  3.3  29/03/94  Correct serious error in computing LENOUT.  This can
!LL                 be longer than the input length for timeseries.  TCJ
!LL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
!LL                  1.Removes the limit on primary STASH item numbers.
!LL                  2.Removes the assumption that (section,item)
!LL                    defines the sub-model.
!LL                  3.Thus allows for user-prognostics.
!LL                  Re-Index PPXREF.
!LL  3.5  07/04/95   Recoded for first stage of submodels.  Added
!LL                  routine to pack internal model specific data into
!LL                  a superarray.  K Rogers
!LL  4.0  06/09/95   Pass in stash superarrays from above. K Rogers
!LL  4.1  03/05/96   Remove unnecessary comdecks ARGPTRA,ARGPTRO,
!LL                  ARGCONA,TYPPTRA and TYPPTRO. Add in wave
!LL                  comdecks and wave call to STWORK. K Rogers.
!    4.1  03/04/96   Pass DUMP_PACKim to STWORK. D. Robinson.
!LL  4.3  17/02/97   MPP code : Added code to change to correct
!LL                  decomposition before calling STASH.
!LL                  Add MOS_MASK_LEN to STWORK argument list
!LL                                                    P.Burton
!LL  4.3  16/05/97   Fix for MPP use of MOS diagnostics -
!LL                  calculate LENOUT correctly.    P.Burton
!LL  4.4   7/10/97   FT_FIRSTSTEP needed in STWORK to check whether
!LL                  re-initialised file stream is written to before
!LL                  being opened. R.Rawlins
!LL  4.4  13/10/97   Pass LEN_A/O/W_SPSTS to STWORK. D. Robinson.
!LL  4.5  13/01/98   Add global_LENOUT argument to STWORK   P.Burton
!LL  5.2  11/02/00 Pass Run Length Encoding logical to STWORK I. Edmond
!LL  5.0  19/05/99   Changed variable names
!LL                                                        P.Burton
!LL  5.0  11/05/99   Change STASH argument list to generic form.
!LL                  Replace im/jm/km references by im/jm/km _ui.
!LL                  R. Rawlins
!LL  5.1  25/01/00   Corrections for uv points on B grid: Amend
!LL                  field_len2 fed into atmos STWORK. R.Rawlins
!    5.3  05/10/01   Remove superfluous defined(mpp). R Rawlins
!    5.5  17/02/03   Fix for Wave model: pass NGX into row_len
!                    D.Holmes-Bell
!    5.5  17/02/03   Amendment for Wave model D.Holmes-Bell
!    5.5  03/02/03   River routing support. P.Selwood
!    6.1  18/08/04   Remove repeated declarations. Declare LEN2_LOOKUP
!                    before the arrays it dimensions.  P.Dando
!    6.2  23/11/05   Removed all references to the wavemodel.
!                    T.Edwards
!    6.2  10/08/05   Make global_lenout>=lenout. R Barnes.
!LL
!LL
!LL Programming standard : UM Doc Paper no 3
!LL
!LL Logical components covered : C3, C4, C8
!LL
!LL Project task : C4
!LL
!LL External documentation : UMDP no C4
!LL
!L* Interface and arguments --------------------------------------------
!
      SUBROUTINE STASH(sm_ident,im_ident,IS,STASH_WORK,                 &
#include "argd1.h"
#include "argdumg.h"
#include "argsts.h"
#include "argppx.h"
     &                 ICODE,CMESSAGE)
!
Use flumerun  ! FLUME-STASH

      IMPLICIT NONE

      integer LEN2_LOOKUP            ! length 2nd dim LOOKUP
#include "cmaxsize.h"
#include "csubmodl.h"
#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
#include "typdumg.h"
#include "typsts.h"
!!! #include "typcona.h" !!! not used, but cconsts.h is !!! 
#include "cconsts.h"
#include "typlndm.h"

      INTEGER                                                           &
     &    sm_ident                                                      &
                         !IN      Submodel identifier
     &   ,im_ident                                                      &
                         !IN      Internal model identifier
     &   ,IS                                                            &
                         !IN      SECTION NUMBER
     &   ,ICODE          !OUT     RETURN CODE
      REAL                                                              &
     &    STASH_WORK(*)  !IN     Area holding the data if not in D1
      CHARACTER*(80)                                                    &
     &    CMESSAGE       !OUT     ANY ERROR MESSAGE PASSED BACK

!*---------------------------------------------------------------------
!
! Common blocks and PARAMETERs
!
#include "ppxlook.h"
#include "chsunits.h"
#include "clookadd.h"
#include "stparam.h"
#include "ccontrol.h"
#include "ctime.h"
#include "chistory.h"
#include "cntl_io.h"
#if defined(MPP)
#include "decomptp.h"
#endif
!
! Subroutines called
!
      EXTERNAL STWORK
!
! Local variables.
!
      LOGICAL                                                           &
     &   LCYCLIC                       ! TRUE if submodel is cyclic

      INTEGER                                                           &
     &   IE,                                                            &
                                       ! Index over items in section
     &   ILEND,                                                         &
                                       ! End point in STASHlist
     &   ILSTART,                                                       &
                                       ! Start point in STASHlist
     &   IL,                                                            &
                                       ! STASHlist index
     &   IM,                                                            &
                                       ! Item number in section
     &   IPPX,                                                          &
                                       ! Index to record in PP_XREF
     &   LENOUT                                                         &
                                       ! Maximum output length
     &   ,STEP                                                          &
                                 ! Step number in integration
     &   ,STEPS_PER_PERIOD                                              &
                                 ! No of steps per period
     &   ,SECS_PER_PERIOD                                               &
                                 ! No of seconds per period
     &   ,im_index                                                      &
                                 ! Internal model index number
     &   ,num_rows1                                                     &
                                 ! Number of rows in field type 1
     &   ,num_rows2                                                     &
                                 ! Number of rows in field type 1
     &   ,row_len                                                       &
                                 ! Row length
     &   ,field_len1                                                    &
                                 ! Length of field type 1
     &   ,field_len2                                                    &
                                 ! Length of field type 2
     &   ,num_levels                                                    &
                                 ! Number of levels
#if defined(OCEAN)
     &   ,LEN_OCNWORK                                                   &
                                 ! Length of ocean work array
#endif
     &   ,orig_decomp                                                   &
                                 ! Decomposition on entry
     &,   global_LENOUT          ! Size of output field on disk
#if !defined(ATMOS)
#if defined(OCEAN)
      LOGICAL ELF                ! True if input grid rotated
                                 ! equatorial
#endif
#endif



!L---------------------------------------------------------------------
!L 0. Initialise.
!L    Set LCYCLIC to indicate EW boundary condition.
!L
! Find current decomposition, so we can return to this after STASH
      orig_decomp=current_decomp_type

      ICODE = 0
#if defined(ATMOS)
      IF (sm_ident == atmos_sm) LCYCLIC = .NOT. ELF
#endif
#if defined(OCEAN)
      IF (sm_ident == ocean_sm) LCYCLIC = CYCLIC_OCEAN
#endif

      im_index = internal_model_index(im_ident)

      ! FLUME-STASH
      ! In a Flume run, ATM processes do not call STASH - 
      ! data is sent from copydiag (or dyndiag, phydiag, etc.) 
      ! and received in Flume STASH wrapper which calls stwork 
      IF(Flume_run) return

!L---------------------------------------------------------------------
!L 1. Loop over items within this section and call STWORK with
!L    appropriate argument list depending on whether atmosphere/ocean
!L
      DO 100 IE=1,NITEMS    ! max no of items in this section
        IF(STINDEX(2,IE,IS,im_index) >  0) THEN ! Entries for this
                                                ! SECTION/ITEM
          ILSTART=STINDEX(1,IE,IS,im_index)
          ILEND=ILSTART+STINDEX(2,IE,IS,im_index)-1
          IM=STLIST(1,ILSTART)             ! item number

          IF(SF(IM,IS)) THEN       ! item/section reqd for this t/s
            IF(STLIST(st_proc_no_code,STINDEX(1,IE,IS,im_index))        &
     &                                                  /= 0) THEN
! required by STASH so continue.
! It should not be possible to have any multiple entries for a given
! ITEM SECTION number and one of those not be required by STASH

#if defined(PRINT84)
              WRITE(6,109) IM,IS
  109         FORMAT(' ITEM',I4,' SECTION',I4,'REQUIRED BY STASH')
#endif

! NB: max poss output length for the item/sect can be LONGER than
!     the input length in the case of timeseries.
              LENOUT=0
              global_LENOUT=0
              DO IL=ILSTART,ILEND
                IF (STLIST(st_output_code,IL)  ==  -89) THEN ! MOS
                  LENOUT=MAX(LENOUT,MOS_OUTPUT_LENGTH*MODEL_LEVELS)
                  global_LENOUT=MAX(global_LENOUT,                      &
     &                              MOS_OUTPUT_LENGTH*MODEL_LEVELS)
                ELSE
                  LENOUT=MAX(LENOUT,STLIST(st_output_length,IL))
                  global_LENOUT=                                        &
     &              MAX(global_LENOUT,                                  &
     &                  STLIST(st_dump_level_output_length,IL))
                ENDIF
              ENDDO
! Add on an extra UM_SECTOR_SIZE on the end, ensuring array is
! big enough for data+extra space to round up to the next
! UM_SECTOR_SIZE
              global_LENOUT=global_LENOUT+UM_SECTOR_SIZE
! Make sure global_lenout at least as large as lenout
! so that time series data are not corrupted
              global_LENOUT = MAX(global_LENOUT,LENOUT)

! Set general variables

              STEP = STEPim(im_ident)
              STEPS_PER_PERIOD = STEPS_PER_PERIODim(im_ident)
              SECS_PER_PERIOD = SECS_PER_PERIODim(im_ident)

!    DUMP_PACKim controls the packing of fields in a dump
!    1 : Use PPXREF file to control packing
!    2 : Do not pack prognostics, as 1 for diagnostics
!    3 : Do not pack prognostics or diagnostics

! Make superarrays to pass into STWORK
#if defined(ATMOS)
              if (im_ident  ==  atmos_sm) then

! Change to atmosphere decomposition
                IF (current_decomp_type  /=                             &
     &              decomp_standard_atmos) THEN
! DEPENDS ON: change_decomposition
                  CALL CHANGE_DECOMPOSITION(decomp_standard_atmos,      &
     &                                      ICODE)
                ENDIF
                IF (ICODE  /=  0) THEN
                  CMESSAGE='STASH : Unsupported MPP submodel : atmos'
                  GOTO 999
                ENDIF
                num_rows1  = ROWS
                num_rows2  = n_rows
                row_len    = row_length
                field_len1 = THETA_FIELD_SIZE
                field_len2 = V_FIELD_SIZE
                num_levels = MODEL_LEVELS

! DEPENDS ON: stwork
                CALL STWORK(                                            &
#include "argppx.h"
     &           D1,LEN_TOT,STASH_WORK,STASH_MAXLEN(IS,im_index),LENOUT,&
     &           global_LENOUT,                                         &
     &           1,IS,IM,ILSTART,ILEND,STEP,STEPS_PER_PERIOD,           &
     &           SECS_PER_PERIOD, PREVIOUS_TIME,                        &
     &           STLIST,LEN_STLIST,TOTITEMS,SI,NSECTS,NITEMS,           &
     &           STASH_LEVELS,NUM_STASH_LEVELS,NUM_LEVEL_LISTS,         &
     &           STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,NUM_PSEUDO_LISTS, &
     &           MAX_STASH_LEVS,STTABL,NSTTIMS,NSTTABL,                 &
     &           STASH_SERIES,nstash_series_records,time_series_rec_len,&
     &           stash_series_index,nstash_series_block,                &
     &           MOS_MASK,MOS_MASK_LEN,MOS_OUTPUT_LENGTH,               &
     &           PP_PACK_CODE,MODEL_FT_UNIT,FT_STEPS,FT_FIRSTSTEP,      &
     &           FIXHD,INTHD,                                           &
     &           REALHD,LEN_FIXHD,LEN_INTHD,LEN_REALHD,                 &
     &           LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC,                     &
     &           LOOKUP,LOOKUP,                                         &
                                 ! 2nd copy used as REAL in PP_HEAD
     &           LEN1_LOOKUP,LEN2_LOOKUP,PP_LEN2_LOOKUP,                &
     &           NUNITS,PP_LEN2_LOOK,                                   &
     &           lcyclic,lrle,num_rows1,num_rows2,                      &
     &           row_len,field_len1,field_len2,num_levels,              &
     &           river_rows, river_row_length,                          &
     &           FORECAST_HRS,RUN_INDIC_OP,ELF,FT_LASTFIELD,            &
     &           sm_ident,im_ident,DUMP_PACKim(sm_ident),               &
     &           len_spsts, spsts, spsts, ixsts, len_ixsts,             &
     &           ldump,                                                 &
     &           ICODE,CMESSAGE)
              endif
#endif
#if defined(OCEAN)
              if (im_ident  ==  ocean_sm) then

! Change to ocean (no wrap around points) decomposition
                IF (current_decomp_type  /=                             &
     &              decomp_nowrap_ocean) THEN
! DEPENDS ON: change_decomposition
                  CALL CHANGE_DECOMPOSITION(decomp_nowrap_ocean,        &
     &                                      ICODE)
                ENDIF
                IF (ICODE  /=  0) THEN
                  CMESSAGE='STASH : Unsupported MPP submodel : ocean'
                  GOTO 999
                ENDIF
                IF (CYCLIC_OCEAN) THEN
                  row_len = imtm2
                ELSE
                  row_len = imt_ui
                ENDIF
                len_ocnwork = spsts(ixsts(9))


                num_rows1  = jmt_ui
                num_rows2  = jmtm1
                field_len1 = jmt_ui * row_len
                field_len2 = jmtm1 * row_len
                num_levels = km_ui
                ELF = .FALSE.

! DEPENDS ON: stwork
                CALL STWORK(                                            &
#include "argppx.h"
     &           D1,LEN_TOT,STASH_WORK,STASH_MAXLEN(IS,im_index),LENOUT,&
     &           global_LENOUT,                                         &
     &           LEN_OCNWORK,                                           &
     &           IS,IM,ILSTART,ILEND,STEP,STEPS_PER_PERIOD,             &
     &           SECS_PER_PERIOD, PREVIOUS_TIME,                        &
     &           STLIST,LEN_STLIST,TOTITEMS,SI,NSECTS,NITEMS,           &
     &           STASH_LEVELS,NUM_STASH_LEVELS,NUM_LEVEL_LISTS,         &
     &           STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,NUM_PSEUDO_LISTS, &
     &           MAX_STASH_LEVS,STTABL,NSTTIMS,NSTTABL,                 &
     &           STASH_SERIES,nstash_series_records,time_series_rec_len,&
     &           stash_series_index,nstash_series_block,                &
     &           MOS_MASK,MOS_MASK_LEN,MOS_OUTPUT_LENGTH,               &
     &           PP_PACK_CODE,MODEL_FT_UNIT,FT_STEPS,FT_FIRSTSTEP,      &
     &           FIXHD,INTHD,                                           &
     &           REALHD,LEN_FIXHD,LEN_INTHD,LEN_REALHD,                 &
     &           LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC,                     &
     &           LOOKUP,LOOKUP,                                         &
                                 ! 2nd copy used as REAL in PP_HEAD
     &           LEN1_LOOKUP,LEN2_LOOKUP,PP_LEN2_LOOKUP,                &
     &           NUNITS,PP_LEN2_LOOK,                                   &
     &           lcyclic,lrle,num_rows1,num_rows2,                      &
     &           row_len,field_len1,field_len2,num_levels,              &
     &           river_rows, river_row_length,                          &
     &           FORECAST_HRS,RUN_INDIC_OP,ELF,FT_LASTFIELD,            &
     &           sm_ident,im_ident,DUMP_PACKim(sm_ident),               &
     &           len_spsts, spsts, spsts, ixsts, len_ixsts,             &
     &           ldump,                                                 &
     &           ICODE,CMESSAGE)

              endif
#endif
#if defined(SLAB)
              if (im_ident  ==  slab_im) then


                num_rows1  = ROWS
                num_rows2  = n_rows
                row_len    = row_length
                field_len1 = THETA_FIELD_SIZE
                field_len2 = U_FIELD_SIZE
                num_levels = MODEL_LEVELS

! DEPENDS ON: stwork
                CALL STWORK(                                            &
#include "argppx.h"
     &           D1,LEN_TOT,STASH_WORK,STASH_MAXLEN(IS,im_index),LENOUT,&
     &           global_LENOUT,                                         &
     &           1,IS,IM,ILSTART,ILEND,STEP,STEPS_PER_PERIOD,           &
     &           SECS_PER_PERIOD, PREVIOUS_TIME,                        &
     &           STLIST,LEN_STLIST,TOTITEMS,SI,NSECTS,NITEMS,           &
     &           STASH_LEVELS,NUM_STASH_LEVELS,NUM_LEVEL_LISTS,         &
     &           STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,NUM_PSEUDO_LISTS, &
     &           MAX_STASH_LEVS,STTABL,NSTTIMS,NSTTABL,                 &
     &           STASH_SERIES,nstash_series_records,time_series_rec_len,&
     &           stash_series_index,nstash_series_block,                &
     &           MOS_MASK,MOS_MASK_LEN,MOS_OUTPUT_LENGTH,               &
     &           PP_PACK_CODE,MODEL_FT_UNIT,FT_STEPS,FT_FIRSTSTEP,      &
     &           FIXHD,INTHD,                                           &
     &           REALHD,LEN_FIXHD,LEN_INTHD,LEN_REALHD,                 &
     &           LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC,                     &
     &           LOOKUP,LOOKUP,                                         &
                                 ! 2nd copy used as REAL in PP_HEAD
     &           LEN1_LOOKUP,LEN2_LOOKUP,PP_LEN2_LOOKUP,                &
     &           NUNITS,PP_LEN2_LOOK,                                   &
     &           lcyclic,lrle,num_rows1,num_rows2,                      &
     &           row_len,field_len1,field_len2,num_levels,              &
     &           river_rows, river_row_length,                          &
     &           FORECAST_HRS,RUN_INDIC_OP,ELF,FT_LASTFIELD,            &
     &           sm_ident,im_ident,DUMP_PACKim(sm_ident),               &
     &           len_spsts, spsts, spsts, ixsts, len_ixsts,             &
     &           ldump,                                                 &
     &           ICODE,CMESSAGE)

              endif
#endif


            ENDIF

          ENDIF

        ENDIF
! Handle error/warning conditions on return from STWORK
        IF (icode >  0) THEN
          WRITE(6,*)'STASH    : Error processing diagnostic section ',  &
     &            is,', item ',im,', code ',icode
          WRITE(6,*)'  ',cmessage
          goto 999
        ELSEIF (icode <  0) THEN
          WRITE(6,*)'STASH    : Warning processing diagnostic section ',&
     &            is,', item ',im,', code ',icode
          WRITE(6,*)'  ',cmessage
          icode=0
        ENDIF
  100 CONTINUE

 999    CONTINUE
      IF (current_decomp_type  /=  orig_decomp) THEN
! DEPENDS ON: change_decomposition
        CALL CHANGE_DECOMPOSITION(orig_decomp,                          &
     &                                      ICODE)
        IF (ICODE  /=  0) THEN
          CMESSAGE='STASH : Unsupported MPP submodel'
        ENDIF
      ENDIF
        RETURN
        END SUBROUTINE STASH
#endif
