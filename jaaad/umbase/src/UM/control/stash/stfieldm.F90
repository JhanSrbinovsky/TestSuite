#if defined(C84_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: STFIELDM -------------------------------------------------
!LL
!LL  Purpose: Calculate weighted field mean within a region specified
!LL           by a lower left hand and upper right hand corner.
!LL           Single level fields only.
!LL           (STASH service routine).
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Author:   T.Johns/S.Tett
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.3  16/09/93  Allow level-by-level mass-weighting if mass-weights
!LL                  are so defined, otherwise use P*.
!LL   4.3  28/01/97  Moved weighting and masking calculations up to
!LL                  SPATIAL.
!LL                  Significantly rewritten for MPP mode - data
!LL                  must be gathered to a processor for
!LL                  reproducible sums to be calculated.   P.Burton
!LL   4.4  13/06/97  MPP: Set fieldout to zero for processors in
!LL                  subdomain area which will not otherwise receive
!LL                  the result of the field mean.
!LL                  MPP: Correct bug in calculating SUMGBOT in non
!LL                       reproducible code                   P.Burton
!LL   4.5  12/01/98  Replaced usage of shmem common block by a
!LL                  dynamic array.                   P.Burton
!LL   4.5  09/01/98  Correct calculation of sum_pe      P.Burton
!LL   5.0  13/07/99  Changes for C-P C grid upgrade.
!LL                  R Rawlins.
!LL   5.0  17/11/99  Removed row_length,p_rows arguments
!LL                                             P.Burton
!LL   5.0  22/06/99  Added halo_type argument            P.Burton
!LL   5.0  15/9/99   Changes to South->North grid         P.Burton
!LL   5.1  15/05/00  S-N ordering consistency correction. R Rawlins
!     6.1  13/07/04  Add packing in 2 stage gather.
!                    P.Selwood/B. Carruthers (Cray)
!     6.2   18/10/05 Fix bugs introduced by gan1f601      Andy Malcolm
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: D715
!LL
!LL  Project task: D7
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage handling and diagnostic
!LL                                 system (STASH)
!LL
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE STFIELDM(fieldin,vx,vy,fld_type,gr,halo_type,          &
     &                    lwrap,lmasswt,                                &
     &                  xstart,ystart,xend,yend,                        &
     &                  global_xstart,global_ystart,                    &
     &                  global_xend,global_yend,                        &
     &                  fieldout,                                       &
     &                  pstar_weight,                                   &
     &                  area_weight,mask,                               &
     &                  level_code,mask_code,weight_code,rmdi,          &
     &                  icode,cmessage)
!
      IMPLICIT NONE
!
      INTEGER                                                           &
     &    vx,vy,                                                        &
                                                ! IN  input field size
     &    fld_type,                                                     &
                                                ! IN field type (u/v/p)
     &    gr,                                                           &
                                                ! IN input fld grid
     &    halo_type,                                                    &
                                                ! IN halo type
     &    xstart,ystart,                                                &
                                                ! IN  lower LH corner
     &    xend,yend,                                                    &
                                                ! IN  upper RH corner
     &    global_xstart,global_ystart,                                  &
                                                ! IN global versions of
     &    global_xend,  global_yend,                                    &
                                                ! IN xstart etc.
     &    level_code,                                                   &
                                                ! IN  input level code
     &    mask_code,                                                    &
                                                ! IN  masking code
     &    weight_code,                                                  &
                                                ! IN  weighting code
     &    icode                                 ! OUT error return code
      CHARACTER*(*)                                                     &
     &    cmessage                              ! OUT error return msg
      LOGICAL                                                           &
     &    lwrap,                                                        &
                                                ! IN  TRUE if wraparound
     &    lmasswt,                                                      &
                                                ! IN  TRUE if masswts OK
     &    mask(vx+1,vy)                         ! IN  mask array
      REAL                                                              &
     &    fieldin(vx,vy),                                               &
                                                ! IN  input field
     &    fieldout,                                                     &
                                                ! OUT output field
     &    pstar_weight(vx+1,vy),                                        &
                                                ! IN  mass weight factor
     &    area_weight(vx+1,vy),                                         &
                                                ! IN  area weight factor
! (already interpolated to the correct grid and
!  set to 1.0 where no area weighting is required)
     &    rmdi                                  ! IN  missing data indic
!*----------------------------------------------------------------------
!
! External subroutines called
!
!
#include "stparam.h"
#include "sterr.h"
#include "i_stgfld.h"
!
! Local variables
!
        INTEGER i,ii,j ! ARRAY INDICES FOR VARIABLE


#include "parvars.h"


#if defined(REPROD)

      INTEGER                                                           &
! Co-ords to PE at top left of subarea
     &  proc_top_left_x , proc_top_left_y                               &

! unused return values from GLOBAL_TO_LOCAL_RC
     &, dummy1 , dummy2                                                 &

! PE number of PE at top left of subarea
     &, sum_pe                                                          &

! size of local and global arrays
     &, local_size,global_size

! Weighted version of fieldin
      REAL local_sum_array_top(xstart:xend,ystart:yend)
! Weights applied to fieldin
      REAL local_sum_array_bot(xstart:xend,ystart:yend)

#else

      INTEGER                                                           &
! limits of local data to be summed
     &  local_sum_xstart,local_sum_xend                                 &
     &, local_sum_ystart,local_sum_yend                                 &

! return code from GCOM routines
     &, info

#endif


#if defined(REPROD)
      INTEGER                                                           &
! Sizes of the global_sum_arrays defined below
     &  global_sum_array_sizex,global_sum_array_sizey

      REAL                                                              &
! Collected versions of fieldin and the weights containing
! whole (subarea) columns of meridional data
     &  global_sum_array_top(global_xstart:global_xend,                 &
     &                       global_ystart:global_yend)                 &
     &, global_sum_array_bot(global_xstart:global_xend,                 &
     &                       global_ystart:global_yend)

#else

! sum(1) is equivalenced to SUMFBOT
! sum(2) is equivalenced to SUMFTOP
      REAL sum(2)

      EQUIVALENCE                                                       &
     &  (sum(1) , SUMFBOT ) , (sum(2) , SUMFTOP)

#endif

        REAL SUMFTOP
        REAL SUMFBOT

!L----------------------------------------------------------------------
!L 0. Initialise sums
!L
      SUMFTOP=0.0
      SUMFBOT=0.0
!L----------------------------------------------------------------------

#if defined(REPROD)
! Create arrays of weighted data suitable to be summed

! Only do the calculations if some of the subarea is contained
! within this processor
      IF ((xstart  /=  st_no_data) .AND. (xend  /=  st_no_data) .AND.   &
     &    (ystart  /=  st_no_data) .AND. (yend  /=  st_no_data)) THEN


        DO i=xstart,xend
            IF ( lwrap .AND.                                            &
     &          (i  >   (lasize(1,fld_type,halo_type)-                  &
     &                   halosize(1,halo_type)))) THEN
! miss halos on wrap around
              ii=i-blsize(1,fld_type)
          ELSE
            ii=i
          ENDIF
          DO j=ystart,yend
            IF (mask(ii,j)) THEN
                local_sum_array_bot(i,j)=                               &
     &            pstar_weight(ii,j)*area_weight(ii,j)
                local_sum_array_top(i,j)=                               &
     &            fieldin(ii,j)*pstar_weight(ii,j)*area_weight(ii,j)
            ELSE
              local_sum_array_bot(i,j)=0.0
              local_sum_array_top(i,j)=0.0
            ENDIF
          ENDDO
        ENDDO

      ENDIF  ! if this processor contains any of the subarea

! Initialise fieldout - so all PE's have valid data
! (Only PEs on top left of subdomain get the field mean)

      fieldout=0.0


! The local_sum_arrays must be distributed so that the complete
! sub-area exists on a single processor, so that a reproducible sum
! can be carried out.

! 0.0 : Initialise variables defining the size of the arrays
!       global_sum_arrays

      global_sum_array_sizex=global_xend-global_xstart+1
      global_sum_array_sizey=global_yend-global_ystart+1

! 1.0 Gather the fields to a single processor

! DEPENDS ON: global_to_local_rc
      CALL GLOBAL_TO_LOCAL_RC(gr,halo_type,                             &
     &  global_xstart , global_ystart,                                  &
     &  proc_top_left_x, proc_top_left_y,                               &
     &  dummy1,dummy2)

      sum_pe=proc_top_left_x + nproc_x*proc_top_left_y

      local_size=(xend-xstart+1)*(yend-ystart+1)
      global_size=global_sum_array_sizex*global_sum_array_sizey

! DEPENDS ON: stash_gather_field
      CALL STASH_GATHER_FIELD (                                         &
     &  local_sum_array_top , global_sum_array_top ,                    &
     &  local_size          , global_size,                              &
     &  1,                                                              &
            ! 1 level
     &  global_yend, global_xend, global_ystart, global_xstart,         &
     &  gr , halo_type, sum_pe,                                         &
     &  .TRUE.,                                                         &
                ! data has been extracted
     &  ICODE=ICODE, CMESSAGE=CMESSAGE)

      IF (ICODE  /=  0) THEN
        WRITE(6,*) 'STFIELDM : MPP Error in STASH_GATHER_FIELD'
        WRITE(6,*) CMESSAGE
        GOTO 999
      ENDIF

! DEPENDS ON: stash_gather_field
      CALL STASH_GATHER_FIELD (                                         &
     &  local_sum_array_bot , global_sum_array_bot ,                    &
     &  local_size          , global_size,                              &
     &  1,                                                              &
            ! 1 level
     &  global_yend, global_xend, global_ystart, global_xstart,         &
     &  gr , halo_type, sum_pe,                                         &
     &  .TRUE.,                                                         &
                ! data has been extracted
     &  ICODE=ICODE, CMESSAGE=CMESSAGE)

      IF (ICODE  /=  0) THEN
        WRITE(6,*) 'STFIELDM : MPP Error in STASH_GATHER_FIELD'
        WRITE(6,*) CMESSAGE
        GOTO 999
      ENDIF

! 2.0 Calculate the sums

      IF (mype  ==  sum_pe) THEN

        DO i=global_xstart,global_xend
          DO j=global_ystart,global_yend
            SUMFTOP=SUMFTOP+global_sum_array_top(i,j)
            SUMFBOT=SUMFBOT+global_sum_array_bot(i,j)
          ENDDO
        ENDDO

        IF (SUMFBOT  ==  0.0) THEN
          fieldout=rmdi
        ELSE
         fieldout=SUMFTOP/SUMFBOT
        ENDIF

      ENDIF

#else

! 1.0 Find the bounds of the actual data required in the summation
!    (ie. excluding the halos, contained within
!    xstart,xend,ystart,yend.

! DEPENDS ON: global_to_local_subdomain
      CALL GLOBAL_TO_LOCAL_SUBDOMAIN(.FALSE.,.FALSE.,                   &
     &  gr,halo_type,mype,                                              &
     &  global_ystart,global_xend,                                      &
     &  global_yend,global_xstart,                                      &
     &  local_sum_ystart,local_sum_xend,                                &
     &  local_sum_yend,local_sum_xstart)

      IF (local_sum_xstart  >   local_sum_xend)                         &
     &  local_sum_xend=local_sum_xend+vx-2*Offx

! 2.0 Calculate the partial sums

! Only do the calculations if some of the subdomain exists on this
! processor

      IF ( (local_sum_xstart  /=  st_no_data) .AND.                     &
     &     (local_sum_xend    /=  st_no_data) .AND.                     &
     &     (local_sum_ystart  /=  st_no_data) .AND.                     &
     &     (local_sum_yend    /=  st_no_data)) THEN

! 2.2 Do the actual sum

        DO i=local_sum_xstart,local_sum_xend

            IF ( lwrap .AND.                                            &
     &          (i  >   (lasize(1,fld_type,halo_type)-                  &
     &                   halosize(1,halo_type)))) THEN
! miss halos on wrap around
              ii=i-blsize(1,fld_type)
          ELSE
            ii=i
          ENDIF

          DO j=local_sum_ystart,local_sum_yend
            IF (mask(ii,j)) THEN

                SUMFBOT=SUMFBOT+                                        &
     &            pstar_weight(ii,j)*area_weight(ii,j)
                SUMFTOP=SUMFTOP+                                        &
     &            fieldin(ii,j)*pstar_weight(ii,j)*area_weight(ii,j)

            ENDIF ! if this point is to be processed
          ENDDO ! j : loop over rows
        ENDDO ! i : loop over columns
      ENDIF ! if subdomain covers this processor

! 3.0  add all the partial sums together, and store

! sum(1) is equivalenced to SUMFTOP
! sum(2) is equivalenced to SUMFBOT

      CALL GC_RSUM(2,nproc,info,sum)

      IF ( (local_sum_xstart  /=  st_no_data) .AND.                     &
     &     (local_sum_xend    /=  st_no_data) .AND.                     &
     &     (local_sum_ystart  /=  st_no_data) .AND.                     &
     &     (local_sum_yend    /=  st_no_data)) THEN

        IF (SUMFBOT  ==  0.0) THEN
          fieldout=rmdi
        ELSE
          fieldout=SUMFTOP/SUMFBOT
        ENDIF
      ENDIF

#endif
  999 CONTINUE
      RETURN
      END SUBROUTINE STFIELDM
#endif
