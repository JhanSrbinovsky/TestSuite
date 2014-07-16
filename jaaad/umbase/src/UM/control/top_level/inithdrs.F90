#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE INITHDRS -----------------------------------------------
!LL
!LL  PURPOSE:   Initialises dump LOOKUP headers reserved for diagnostic
!LL             fields with size and other basic information to allow
!LL             dump IO routines to work correctly before STASH has
!LL             updated the addressed fields.
!LL
!LL  PROGRAMMING STANDARD: UNIFIED MODEL DP NO. 3, VERSION 3
!LL
!LL  SYSTEM TASK: C4
!LL
!LL  SYSTEM COMPONENTS: C401
!LL
!LL  EXTERNAL DOCUMENTATION: UMDP NO. C4 VERSION NO. 8
!LL
!LLEND-------------------------------------------------------------

      SUBROUTINE INITHDRS(                                              &
#include "argduma.h"
#include "argsts.h"
#include "argppx.h"
     &                  ICODE,CMESSAGE)

      IMPLICIT NONE

!*L Arguments
!L
#include "csubmodl.h"
#include "parparm.h"
#include "typsize.h"
#include "typduma.h"
! Contains *CALL CPPXREF
#include "typsts.h"
      INTEGER                                                           &
     &    ICODE                  ! OUT: Error return code
!
      CHARACTER*80  CMESSAGE     ! OUT: Error return message

! PARAMETERs

#include "chsunits.h"
#include "ccontrol.h"
#include "c_mdi.h"
#include "clookadd.h"
#include "stparam.h"

! PPXREF lookup arrays

#include "ppxlook.h"

! Local variables

      REAL                                                              &
     &        r_eqv_rmdi
      INTEGER                                                           &
     &        i_eqv_rmdi                                                &
     &,disk_address                                                     &
                                       ! Current rounded disk address
     &,number_of_data_words_on_disk                                     &
                                       ! Number of data words on disk
     &,number_of_data_words_in_memory  ! Number of Data Words in memory
      INTEGER                                                           &
     &        i,ii,is,im,iproc,ilookup,iheaders,ilength,imean,j         &
     &       ,im_ident                                                  &
                              ! internal model identifier
     &       ,sm_ident                                                  &
                              ! submodel partition (dump) identifier
     &       ,N1              ! Packing Indicator for Lookup(21)

      EQUIVALENCE(r_eqv_rmdi,i_eqv_rmdi)

! Function and subroutine calls:
      INTEGER  EXPPXI
      EXTERNAL EXPPXI

!L----------------------------------------------------------------------
!L  1. Set dump LOOKUP headers with basic information needed by
!L     READDUMP and WRITDUMP, by scanning STASHlist items for
!L     diagnostics destined for dump addresses.  NB: timeseries
!L     fields cannot be 32-bit packed in dumps as the extra data
!L     will contain integers.
!L
      r_eqv_rmdi=rmdi
!     gol124: ensure sm_ident always initialized to prevent
!     problems when running under debugger
      sm_ident = 0
!
      DO II=1,totitems
        IF (stlist(st_output_code,II) == st_dump) THEN
!       output is to addressed D1

          im_ident  =stlist(st_model_code,II)          ! internal model
          sm_ident  =SUBMODEL_PARTITION_INDEX(im_ident)! submodel

          is     =stlist(st_sect_no_code,II)
          im     =stlist(st_item_code,II)
          iproc  =stlist(st_proc_no_code,II)
          ilookup=stlist(st_lookup_ptr,II)

          IF (sm_ident /= atmos_sm) THEN
            icode=111
            cmessage='INITHDRS: Only Atmos diagnostic requests possible'
            GO TO 9999
          END IF

!         Calculate the total number of headers required by this
!         stashlist record.
          imean=(stlist(st_gridpoint_code,II)/block_size)*block_size
          IF (stlist(st_output_bottom,II) == 100) THEN
!           single level
            iheaders=1
          ELSE IF (stlist(st_proc_no_code,II)  ==                       &
     &              st_time_series_code.OR.                             &
     &      stlist(st_proc_no_code,II) == st_time_series_mean) THEN
!           time series
            iheaders=1
          ELSE IF (stlist(st_proc_no_code,II)  ==                       &
     &             st_append_traj_code) THEN
!           append trajectories
            iheaders=1
          ELSE IF (imean == vert_mean_base) THEN
!           vertical mean
            iheaders=1
          ELSE IF (imean == global_mean_base) THEN
!           total 3-D mean
            iheaders=1
          ELSE IF (stlist(st_output_bottom,II) <  0) THEN
!           level list, not vertical mean.
            iheaders=STASH_LEVELS(1, -stlist(st_output_bottom,II) )
          ELSE
!           level range, not vertical mean.
            iheaders=stlist(st_output_top,II)-                          &
     &               stlist(st_output_bottom,II)+1
          END IF
          IF(stlist(st_pseudo_out,II) >  0) THEN !Output pseudo levs
            iheaders=iheaders*                                          &
     &      STASH_PSEUDO_LEVELS(1,stlist(st_pseudo_out,II))
          END IF

#if !defined(MPP)
          ilength=stlist(st_output_length,II) / iheaders
#else
          ilength=stlist(st_dump_output_length,II) / iheaders
#endif
!         Loop down the headers.
          DO I=0,iheaders-1
#if defined(ATMOS)
            IF (sm_ident == atmos_sm) THEN
              DO j=1,len1_lookup
                a_lookup(j,ilookup+I)=imdi
              END DO
              DO j=46,len1_lookup
                a_lookup(j,ilookup+I)=i_eqv_rmdi
              END DO
              a_lookup(lbnrec   ,ilookup+I)=0
              a_lookup(item_code,ilookup+I)=is*1000+im
              a_lookup(model_code,ilookup+I)=im_ident
              a_lookup(data_type,ilookup+I)=                            &
! DEPENDS ON: exppxi
     &                          EXPPXI(im_ident,is,im,ppx_data_type,    &
#include "argppx.h"
     &                                             ICODE,CMESSAGE)
              a_lookup(lblrec,ilookup+I)=ilength
#if !defined(MPP)
              a_lookup(naddr ,ilookup+I)=stlist(st_output_addr  ,II)+   &
     &                                   ( ilength * I )
#else
              a_lookup(naddr ,ilookup+I)=                               &
     &          stlist(st_dump_output_addr  ,II)+( ilength * I )
#endif
              IF (iproc == st_time_series_code .OR.                     &
     &            iproc == st_time_series_mean .OR.                     &
     &        iproc == st_append_traj_code) THEN
                a_lookup(lbpack,ilookup+I)=2000
              ELSE
                a_lookup(lbpack,ilookup+I)=2000+                        &
! DEPENDS ON: exppxi
     &                          EXPPXI(im_ident,is,im,ppx_dump_packing, &
#include "argppx.h"
     &                                             ICODE,CMESSAGE)
                IF (DUMP_PACKim(sm_ident) == 3 ) THEN
!                 Do not pack data ; Override PPXREF packing indicator
                  N1 = 0   !   No packing
                  a_lookup(lbpack,ilookup+I) =                          &
     &           (a_lookup(lbpack,ilookup+I)/10)*10 + N1
                END IF
              END IF
            END IF
#endif
          END DO ! I, Loop over headers for this STASHlist entry

        END IF
      END DO  ! II LOOP OVER TOTITEMS

!
!--reset the disk addresses and lengths for well-formed I/O
#if defined(ATMOS)
      if (sm_ident == atmos_sm) then
! DEPENDS ON: set_dumpfile_address
        call set_dumpfile_address(a_fixhd, len_fixhd,                   &
     &                            a_lookup, len1_lookup,                &
     &                            a_len2_lookup,                        &
     &                            number_of_data_words_in_memory,       &
     &                            number_of_data_words_on_disk,         &
     &                            disk_address)
      endif
#endif
 9999 CONTINUE
      RETURN
      END SUBROUTINE INITHDRS
#endif
