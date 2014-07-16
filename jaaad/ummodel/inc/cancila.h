!*L--------------------COMDECK  CANCILA ---------------------------
!
! Purpose : Contains index blocks for control of update of
!           ancillary fields.
!
! System component F0171
!
! History :
! Version   Date   Comment
! -------   ----   -------
!  3.4   23/06/94  Update comments. D. Robinson
!  3.4   04/10/94  Increase NANCIL_FIELDS from 43 to 71. RTHBarnes
!  4.1   22/05/96  Move NANCIL_FIELDS to comdeck CANCMAXA. D. Robinson.
!  4.4   28/07/97  Add LAMIPII to common block. R A Stratton
!  6.2   22/08/05  Remove un-needed ampersand. P.Selwood
!
! -------------------------------------------------------------------
!
#include "cancmaxa.h"
! Type Declarations

      INTEGER                                                           &
     &  FILEANCIL,                                                      &
                         ! File number associated with ancillary fields
     &  NLOOKUP,                                                        &
                         ! Position of ancillary field in lookup tables.
     &  LOOKUP_STEP,                                                    &
                         ! Interval between PP Headers refering to
!                        ! to the same ancillary fields at diferent time
     &  LEVELS,                                                         &
                         ! Number of levels of data in each ancillary
!                        ! field.
     &  STASHANCIL,                                                     &
                         ! Stash codes for ancillary files
     &  D1_ANCILADD      ! Address of ancillary field in main data block


      COMMON/IXANCILA/ FILEANCIL(NANCIL_FIELDS),                        &
     &           NLOOKUP(NANCIL_FIELDS),                                &
     &           LOOKUP_STEP(NANCIL_FIELDS),                            &
     &           LEVELS(NANCIL_FIELDS),                                 &
     &           STASHANCIL(NANCIL_FIELDS),                             &
     &           D1_ANCILADD(NANCIL_FIELDS)

!*L---------- Control data calculated from NAMELIST-------------------
      LOGICAL                                                           &
     &         UPDATE                                                   &
     &,      L_SSTANOM                                                  &
                                ! Indicator if SST anom to be formed
                                ! (RECON=T) or used (-DEF,RECON)
     & ,     LAMIPII            ! True if special AMIP II updating

      INTEGER  FIELDCODE,                                               &
     &         STEPS

    !kdcorbin, 05/10
    INTEGER :: ANC_FILE_FINPUT
    CHARACTER*120 :: ANC_FILE_FNAME    
    
!*----------------------------------------------------------------------
      COMMON/CTANCILA/                                                  &
     &                L_SSTANOM,LAMIPII,                                &
#if !defined(RECON)
     &         FIELDCODE(2,NANCIL_FIELDS),                              &
#else
     &         FIELDCODE(NANCIL_FIELDS),                                &
#endif
     &         STEPS(NANCIL_FIELDS),UPDATE(NANCIL_FIELDS),              &
               ANC_FILE_FINPUT(NANCIL_FIELDS),ANC_FILE_FNAME(NANCIL_FIELDS)
!kdcorbin, 05/10 - added anc_file_finput and anc_file_fname
