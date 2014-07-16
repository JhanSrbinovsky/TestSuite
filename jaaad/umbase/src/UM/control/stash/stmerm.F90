#if defined(C84_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: STMERM ---------------------------------------------------
!LL
!LL  Purpose: Calculate weighted meridional mean within a region
!LL           specified by lower left hand and upper right hand corner.
!LL           (STASH service routine).
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Author:   T.Johns/S.Tett
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: D714
!LL
!LL  Project task: D7
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage handling and diagnostic
!LL                                 system (STASH)
!LL
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE STMERM(fieldin,vx,vy,fld_type,gr,halo_type,            &
     &                  lwrap,lmasswt,                                  &
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
     &    fieldout(xstart:xend),                                        &
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
!
! Local variables
!
        INTEGER i,ii,j   ! ARRAY INDICES FOR VARIABLE


#include "parvars.h"

#if defined(REPROD)
      INTEGER                                                           &
! Processor co-ordinates of processors at the corners of the
! processed subdomain
     &  proc_top_left_x,proc_top_left_y                                 &
     &, proc_bot_right_x,proc_bot_right_y                               &

! size of the full subarea in both horizonal dimensions
     &, merid_sum_global_len_x                                          &
     &, merid_sum_global_len_y                                          &

! loop variables for loops over processors in subdomain
     &, proc_x,proc_y                                                   &

! real processor x co-ordinate - when proc_x > nproc_x is just
! proc_x-nproc_x
     &, eff_proc_x                                                      &

! processor id of processor (proc_x,proc_y)
     &, proc_id                                                         &

! definition of the extracted subarea array on processor proc_id
     &, local_array_top_left_x,local_array_top_left_y                   &
     &, local_array_bot_right_x,local_array_bot_right_y                 &

! definition of the real data contained within the extracted
! subarea on processor proc_id (ie. not including halos)
     &, local_data_top_left_x,local_data_top_left_y                     &
     &, local_data_bot_right_x,local_data_bot_right_y                   &

! size in the x dimension of the subarea array on proc_id
     &, local_array_size_x                                              &

! length of (partial) local meridional data to be sent
     &, local_send_size_y                                               &

! offset of data to be sent from start of local row
     &, local_send_off_x                                                &

! position in the full meridional column of (partial)
! data to be sent
     &, pos_in_merid_array                                              &

! number of (partial) meridional mean columns on
! processor proc_id
     &, local_n_cols_to_send                                            &

! the first point on the column to be sent to processor proc_id
     &, local_send_off_y                                                &

! the global meridional mean number of the first local
! (partial) meridional mean column to be sent from proc_id
     &, global_merid_col_number_start                                   &

! loop counter for loop over columns to send
     &, col                                                             &

! index of column in proc_id's array of local data
     &, local_col                                                       &

! global meridional column number of a column
     &, global_merid_col_id                                             &

! processor which this meridional mean column will be sent to
     &, dest_proc                                                       &

! column number on the destination processor
     &, work_dest_col_id                                                &

! number of items of meridional data to send and receive
     &, n_send_data , n_rec_data                                        &

! number of final meridional means to send and receive
     &, n_send_means , n_rec_means                                      &

! number of columns of (full) meridional data on this processor
     &, n_cols_full_merid_data                                          &

! size of local_sum_arrays and global_sum_arrays
     &, local_sum_array_len                                             &
     &, global_sum_array_len                                            &


! arguments for GCOM routines
     &, flag , info                                                     &

! dummy variables (unused return values from subroutine calls)
     &, dummy1,dummy2
#else
      INTEGER                                                           &
! size of subarea on this processor, not including halo areas
     &  local_sum_xstart,local_sum_xend                                 &
     &, local_sum_ystart,local_sum_yend                                 &

! number of columns of meridional data to sum on this processor
! (xend-xstart+1)
     &, partial_sum_data_sizex                                          &

! return code from GCOM routines
     &, info
#endif


#if defined(REPROD)
      LOGICAL                                                           &

! indicates if the subarea requested for meridional meaning wraps
! over zero longitude
     &  lwrap_merid_mean                                                &

! indicates if the subdomain contains processors which hold both
! the start and end of the subdomain, which wraps over zero
! longitude
     &, lwrap_proc                                                      &

! indicates that a full field is being zonal meaned
     &, fullfield

      REAL                                                              &
! temporary variables used in calculation of meridional means
     &  merid_sum_top,merid_sum_bot

      INTEGER                                                           &

! Send/receive maps for meridional data arrays to be summed
     &  send_data_map(7,(xend-xstart+1))                                &
     &, rec_data_map(7,(global_xend-global_xstart+1)*nproc)             &

! send/receive maps for meridional means
     &, send_means_map(7,(global_xend-global_xstart+1))                 &
     &, rec_means_map(7,(xend-xstart+1))

! Weighted version of fieldin
      REAL local_sum_array_top(xstart:xend,ystart:yend)
! Weights applied to fieldin
      REAL local_sum_array_bot(xstart:xend,ystart:yend)
#endif


#if defined(REPROD)
      INTEGER                                                           &
! Sizes of the global_sum_arrays defined below
     &  global_sum_array_sizex,global_sum_array_sizey

      REAL                                                              &
! Collected versions of fieldin and the weights containing
! whole (subarea) columns of meridional data
     &  global_sum_array_top(global_ystart:global_yend,                 &
     &                       global_xend-global_xstart+1)               &
     &, global_sum_array_bot(global_ystart:global_yend,                 &
     &                       global_xend-global_xstart+1)               &

! Calculated meridional means on the calculating processor
     &, merid_mean_array(global_xend-global_xstart+1)

#else
      INTEGER                                                           &
! Size of local partial sum arrays defined below
     &  partial_sum_array_sizex

      REAL                                                              &
! Partial meridional sums of subarea columns
     &  partial_SUMMTOP(vx*2)                                           &
     &, partial_SUMMBOT(vx*2)

#endif

#if defined(REPROD)

! Integer function used for obtaining field type
      INTEGER GET_FLD_TYPE

#endif

#include "gccom.h"

        REAL SUMMTOP(xstart:xend)
        REAL SUMMBOT(xstart:xend)

!L----------------------------------------------------------------------
!L 0. Initialise sums
!L
!FPP$ NOINNER R
      DO i=xstart,xend
        SUMMTOP(i)=0.0
        SUMMBOT(i)=0.0
      ENDDO
!L----------------------------------------------------------------------

! pstar_weight and area_weight arrays contain appropriate
! weighting factors, interpolated to the correct grid, for
! mass weighting and area weighting respectively. If either type
! of weighting is not required, the relevant array is set to 1.0
! The mask array contains appropriate masking

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

! Initialise fieldout array - so all PE's have valid data
! (Only PEs on top of subdomain get the meridional means)
      DO i=xstart,xend
        fieldout(i)=0.0
      ENDDO

! The local_sum_arrays must be distributed so that complete
! sub-area columns exist on processors, so that a reproducible sum
! can be carried out.
! The following code calculates where the local_sum_array data
! must be sent to, and where the final answers must be sent back to

! 0.0 : Initialise variables defining the size of the arrays
!       global_sum_arrays

      global_sum_array_sizex=global_xend-global_xstart+1
      global_sum_array_sizey=global_yend-global_ystart+1

      local_sum_array_len=((xend-xstart)+1)*((yend-ystart)+1)
      global_sum_array_len=global_sum_array_sizex*                      &
     &                     global_sum_array_sizey

! Set a logicial indicating if the area being meaned is the
! full field

      fullfield= ( (global_xstart  ==  1                 ) .AND.        &
     &             (global_xend    ==  glsize(1,fld_type)) .AND.        &
     &             (global_yend    ==  glsize(2,fld_type))      )

! Calculate the length of the full meridional subarea

      merid_sum_global_len_x=global_xend-global_xstart+1
      merid_sum_global_len_y=global_yend-global_ystart+1

! 1.0 Find the set of processors covering the requested sub-area

! DEPENDS ON: global_to_local_rc
      CALL GLOBAL_TO_LOCAL_RC(gr,halo_type,                             &
     &   global_xstart , global_ystart,                                 &
     &   proc_top_left_x, proc_top_left_y,                              &
     &   dummy1,dummy2)

! DEPENDS ON: global_to_local_rc
      CALL GLOBAL_TO_LOCAL_RC(gr,halo_type,                             &
     &   global_xend,global_yend,                                       &
     &   proc_bot_right_x, proc_bot_right_y,                            &
     &   dummy1,dummy2)

! Set a logical to indicate if the meridional mean area required
! wraps over zero longitude

      lwrap_merid_mean=                                                 &
     & ( (global_xend  >   glsize(1,fld_type)) .OR.                     &
     &   (global_xend  <   global_xstart))

! If there is a wrap around over 0 longitude, ensure that
! proc_bot_right_x > proc_top_left_x

      IF (lwrap_merid_mean)                                             &
     &  proc_bot_right_x=proc_bot_right_x+nproc_x

! Set up a logical to indicate if a processor in the subdomain
! contains both the start and end of a meridional mean which wraps
! over zero longitude. If TRUE, some extra work is required at
! this processor as it contains data for two non-contiguous parts
! of the meridional mean

      lwrap_proc=(proc_bot_right_x  ==  proc_top_left_x+nproc_x)

! 2.0 Loop over all the processors in the subdomain, and set
!     up the send/receive maps defining the redistribution
!     of data

      n_send_data=0            ! number of items of data to send
      n_rec_data=0             ! number of items of data to receive
      n_send_means=0           ! number of merid. means I will send
      n_rec_means=0            ! number of merid. means I will receive
      n_cols_full_merid_data=0 ! number of cols. of data I will mean

      DO proc_y=proc_top_left_y , proc_bot_right_y

        DO proc_x=proc_top_left_x , proc_bot_right_x

          eff_proc_x=MOD(proc_x,nproc_x)
          proc_id=eff_proc_x+proc_y*nproc_x

! 2.1  Find the size of the array containing the meridional
!      arrays on processor proc_id

! DEPENDS ON: global_to_local_subdomain
          CALL GLOBAL_TO_LOCAL_SUBDOMAIN(.TRUE.,.TRUE.,                 &
     &      gr,halo_type,proc_id,                                       &
     &      global_ystart,global_xend,                                  &
     &      global_yend,global_xstart,                                  &
     &      local_array_top_left_y,local_array_bot_right_x,             &
     &      local_array_bot_right_y,local_array_top_left_x)

! 2.2 Using this information, calculate the size of this array in
!     the x dimension. If the data is wrapped round, the calculation
!     is done differently:

          IF (local_array_top_left_x  <=  local_array_bot_right_x)      &
     &    THEN
            local_array_size_x=                                         &
     &        local_array_bot_right_x-local_array_top_left_x+1
          ELSE
            local_array_size_x=                                         &
     &        local_array_bot_right_x-local_array_top_left_x+1+         &
     &        g_blsize(1,fld_type,proc_id)
          ENDIF

! 2.3 Find out the size of the actual meridional mean data within the
!     subarea array on processor proc_id

! DEPENDS ON: global_to_local_subdomain
          CALL GLOBAL_TO_LOCAL_SUBDOMAIN(.FALSE.,.FALSE.,               &
     &      gr,halo_type,proc_id,                                       &
     &      global_ystart,global_xend,                                  &
     &      global_yend,global_xstart,                                  &
     &      local_data_top_left_y,local_data_bot_right_x,               &
     &      local_data_bot_right_y,local_data_top_left_x)

! 2.4 Calculate various quantities:
!     local_send_size_y  : the length of data to be sent
!     local_send_off_y   : the offset of this data from the
!                          start of column
!     pos_in_merid_array : position of this data in the full
!                          meridional array

          local_send_size_y=                                            &
     &      local_data_bot_right_y-local_data_top_left_y+1
          local_send_off_y=                                             &
     &      local_data_top_left_y-local_array_top_left_y
          pos_in_merid_array=                                           &
     &      g_datastart(2,proc_id)+local_data_top_left_y-               &
     &      halosize(2,halo_type)-                                      &
     &      global_ystart

! 2.5 Find the number of meridional mean columns to be sent
!     from this processor, the first column to be sent,
!     and which global meridional mean this is

          IF ((LWRAP_PROC) .AND. (proc_x  ==  proc_top_left_x)) THEN
! Processor containing start and end of sumdomain - but here
! we're interested only in the start segment

            local_n_cols_to_send=                                       &
     &        g_lasize(1,fld_type,halo_type,proc_id) -                  &
     &                               local_data_top_left_x -            &
     &                               halosize(1,halo_type)
            local_send_off_x=                                           &
     &        local_data_top_left_x-local_array_top_left_x
            global_merid_col_number_start=                              &
     &        g_datastart(1,proc_id)+local_data_top_left_x-             &
     &        halosize(1,halo_type)-                                    &
     &        global_xstart

          ELSEIF ((LWRAP_PROC) .AND.                                    &
     &            (proc_x  ==  proc_bot_right_x)) THEN
! Processor containing start and end of subdomain - but here
! we're interested only in the end segment

            local_n_cols_to_send=local_data_bot_right_x-                &
     &      halosize(1,halo_type)
            local_send_off_x=local_array_size_x-local_n_cols_to_send
            global_merid_col_number_start=                              &
     &        merid_sum_global_len_x-local_n_cols_to_send+1

          ELSE
! all other processors

            local_n_cols_to_send=                                       &
     &        local_data_bot_right_x-local_data_top_left_x+1
            local_send_off_x=                                           &
     &        local_data_top_left_x-local_array_top_left_x
            global_merid_col_number_start=                              &
     &        g_datastart(1,proc_id)+local_data_top_left_x-             &
     &        halosize(1,halo_type)-                                    &
     &        global_xstart
          ENDIF

          IF (global_merid_col_number_start  <   1) THEN
! This means the sub-area wraps over zero longitude - so to get
! the correct position in the array we add the global row length
            global_merid_col_number_start=                              &
     &        global_merid_col_number_start+glsize(1,fld_type)
          ENDIF


! 2.6 Loop over columns and construct send/receive maps

          DO col=1,local_n_cols_to_send

! 2.6.1 Find the local column index on proc_id, and the global
!       meridional column index of this column

            local_col=col+local_send_off_x
            global_merid_col_id=global_merid_col_number_start+col-1

! 2.6.2 and find the destination processor of this column, and
!       where on this processor it will be sent to

            dest_proc=MOD(global_merid_col_id-1,nproc)
            work_dest_col_id=((global_merid_col_id-1)/nproc)+1

! 2.6.3 If this processor is proc_id construct a send_data_map
!       entry for this column of data

            IF (mype  ==  proc_id) THEN

              n_send_data = n_send_data+1

              send_data_map(S_DESTINATION_PE,n_send_data)=              &
     &          dest_proc
              send_data_map(S_BASE_ADDRESS_IN_SEND_ARRAY,               &
     &                      n_send_data)=                               &
     &          local_send_off_y*local_array_size_x +                   &
     &          local_send_off_x+col
              send_data_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,               &
     &                      n_send_data)=                               &
     &          local_send_size_y
              send_data_map(S_STRIDE_IN_SEND_ARRAY,n_send_data)=        &
     &          local_array_size_x
              send_data_map(S_ELEMENT_LENGTH,n_send_data)=1
              send_data_map(S_BASE_ADDRESS_IN_RECV_ARRAY,               &
     &                      n_send_data)=                               &
     &          (work_dest_col_id-1)*global_sum_array_sizey +           &
     &          pos_in_merid_array
              send_data_map(S_STRIDE_IN_RECV_ARRAY,n_send_data)=1

! 2.6.3.1 If this processor is at the top of the subarea, then it is
!         responsible for holding the final meridional mean values.
!         So we must set up a rec_means_map entry to allow the
!         meridional mean value for this column to be returned.

              IF (proc_y  ==  proc_top_left_y) THEN

                n_rec_means = n_rec_means + 1

                rec_means_map(R_SOURCE_PE,n_rec_means)=dest_proc
                IF (fullfield) THEN ! We don't want halos
                  rec_means_map(R_BASE_ADDRESS_IN_RECV_ARRAY,           &
     &                          n_rec_means)=                           &
     &              local_col-halosize(1,halo_type)
                ELSE ! halos are automatically removed
                  rec_means_map(R_BASE_ADDRESS_IN_RECV_ARRAY,           &
     &                          n_rec_means)=                           &
     &              local_col
                ENDIF
                rec_means_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,             &
     &                        n_rec_means)= 1
                rec_means_map(R_STRIDE_IN_RECV_ARRAY,                   &
     &                        n_rec_means)= 1
                rec_means_map(R_ELEMENT_LENGTH,n_rec_means)=1
                rec_means_map(R_BASE_ADDRESS_IN_SEND_ARRAY,             &
     &                        n_rec_means)=                             &
     &            work_dest_col_id
                rec_means_map(R_STRIDE_IN_SEND_ARRAY,                   &
     &                        n_rec_means)=1
              ENDIF

            ENDIF

! 2.6.4 If this processor is dest_proc construct a rec_data_map
!       entry for this column of data

            IF (mype  ==  dest_proc) THEN

              IF (proc_y  ==  proc_top_left_y)                          &
! increment counter of full meridional columns on this processor
     &          n_cols_full_merid_data=n_cols_full_merid_data+1

              n_rec_data = n_rec_data+1

              rec_data_map(R_SOURCE_PE,n_rec_data)=                     &
     &          proc_id
              rec_data_map(R_BASE_ADDRESS_IN_RECV_ARRAY,n_rec_data)=    &
     &          (work_dest_col_id-1)*global_sum_array_sizey +           &
     &          pos_in_merid_array
              rec_data_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,n_rec_data)=    &
     &          local_send_size_y
              rec_data_map(R_STRIDE_IN_RECV_ARRAY,n_rec_data)=1
              rec_data_map(R_ELEMENT_LENGTH,n_rec_data)=1
              rec_data_map(R_BASE_ADDRESS_IN_SEND_ARRAY,n_rec_data)=    &
     &          local_send_off_y*local_array_size_x +                   &
     &          local_send_off_x+col
              rec_data_map(R_STRIDE_IN_SEND_ARRAY,n_rec_data)=          &
     &          local_array_size_x

! 2.6.4.1 Set up the send_means_map entry for sending the
!         resulting meridional mean of this column back to
!         the processor at the top of the subarea.
!         We only need to do this once per column (not for
!         each value of proc_y).

              IF (proc_y  ==  proc_top_left_y) THEN

                n_send_means = n_send_means+1

                send_means_map(S_DESTINATION_PE,n_send_means)=          &
     &            proc_id
                send_means_map(S_BASE_ADDRESS_IN_SEND_ARRAY,            &
     &                         n_send_means)=work_dest_col_id
                send_means_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,            &
     &                         n_send_means)=1
                send_means_map(S_STRIDE_IN_SEND_ARRAY,                  &
     &                         n_send_means)=1
                send_means_map(S_ELEMENT_LENGTH,n_send_means)=1
                IF (fullfield) THEN ! we don't want halos
                  send_means_map(S_BASE_ADDRESS_IN_RECV_ARRAY,          &
     &                           n_send_means)=local_col-               &
     &                           halosize(1,halo_type)
                ELSE ! halos are automatically removed
                  send_means_map(S_BASE_ADDRESS_IN_RECV_ARRAY,          &
     &                           n_send_means)=local_col
                ENDIF
                send_means_map(S_STRIDE_IN_RECV_ARRAY,                  &
     &                         n_send_means)=1
              ENDIF ! if at top of subarea

            ENDIF ! if mype  ==  dest_proc

          ENDDO ! col : loop over local columns on proc_id

        ENDDO ! proc_x : loop over processors in x dimension

      ENDDO ! proc_y : loop over processors in y dimension

! 3.0 Now the send and receive maps are set up, use
!     GCG_RALLTOALLE to redistribute the data
!     from the local_sum_arrays to the global_sum_arrays
      flag=GC_NONE ! flag argument is currently ignored by GCOM
      info=GC_NONE

      CALL GCG_RALLTOALLE(                                              &
     &  local_sum_array_top , send_data_map , n_send_data ,             &
     &  local_sum_array_len ,                                           &
     &  global_sum_array_top , rec_data_map , n_rec_data ,              &
     &  global_sum_array_len ,                                          &
     &  gc_all_proc_group , flag , info)

      info=GC_NONE

      CALL GCG_RALLTOALLE(                                              &
     &  local_sum_array_bot , send_data_map , n_send_data ,             &
     &  local_sum_array_len ,                                           &
     &  global_sum_array_bot , rec_data_map , n_rec_data ,              &
     &  global_sum_array_len ,                                          &
     &  gc_all_proc_group , flag , info)

! 4.0 Calculate mean of any meridional data on this processor

      DO i=1,n_cols_full_merid_data

        merid_sum_top=0.0
        merid_sum_bot=0.0

        DO j=global_ystart,global_yend

          merid_sum_top=merid_sum_top+                                  &
     &                     global_sum_array_top(j,i)
          merid_sum_bot=merid_sum_bot+                                  &
     &                     global_sum_array_bot(j,i)
        ENDDO

        IF (merid_sum_bot  ==  0.0) THEN
          merid_mean_array(i)=rmdi
        ELSE
          merid_mean_array(i)=merid_sum_top/merid_sum_bot
        ENDIF

      ENDDO

! 5.0 Send the calculated means back to the processors at the
!     top of the subarea, into the fieldout array

      flag=GC_NONE ! flag argument is currently ignored by GCOM
      info=GC_NONE

      CALL GCG_RALLTOALLE(                                              &
     &  merid_mean_array , send_means_map , n_send_means ,              &
     &  global_sum_array_sizex,                                         &
     &  fieldout , rec_means_map , n_rec_means,                         &
     &  (xend-xstart)+1,                                                &
     &  gc_all_proc_group , flag , info)

#else

! 0.0 : Initialise variables defining the size of the arrays
!       partial_sum_arrays
      partial_sum_array_sizex=vx*2

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
     &  local_sum_xend=local_sum_xend+vx-                               &
     &                 2*halosize(1,halo_type)

! 1.1 And the number of partial sums on this processor
!     (If there is a wrap around, with the subarea starting
!      and ending on this processor, the halo cols are included
!      in this number - although no sums will actually be
!      carried out there)

      partial_sum_data_sizex=                                           &
     &  local_sum_xend-local_sum_xstart+1


! 1.2 Initialise the sum arrays

      IF ((local_sum_xstart  /=  st_no_data) .AND.                      &
     &    (local_sum_xend  /=  st_no_data)) THEN

        DO i=local_sum_xstart,local_sum_xend
          partial_SUMMBOT(i)=0.0
          partial_SUMMTOP(i)=0.0
        ENDDO

      ENDIF

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

                partial_SUMMBOT(i)=partial_SUMMBOT(i)+                  &
     &            pstar_weight(ii,j)*area_weight(ii,j)
                partial_SUMMTOP(i)=partial_SUMMTOP(i)+                  &
     &            fieldin(ii,j)*pstar_weight(ii,j)*area_weight(ii,j)
            ENDIF ! if this point is to be processed
          ENDDO ! j : loop over rows
        ENDDO ! i : loop over columns
      ENDIF ! if subdomain covers this processor

! 3.0 Sums up the partial sums down each column to make a full sum

! So a sum down the processor column if the subdomain covers any
! processor(s) along the column

      IF ((local_sum_xstart  /=  st_no_data) .AND.                      &
     &    (local_sum_xend  /=  st_no_data)) THEN

        CALL GCG_RSUM(partial_sum_data_sizex,gc_proc_col_group,         &
     &                info,partial_SUMMBOT(local_sum_xstart))
        CALL GCG_RSUM(partial_sum_data_sizex,gc_proc_col_group,         &
     &                info,partial_SUMMTOP(local_sum_xstart))

! So now the partial_* arrays actually contain the full sums
! along the column

! 3.1 And put the mean meridional values into the fieldout array

! Only processors in the subdomain area need to record the
! results
        IF ((local_sum_ystart  /=  st_no_data) .AND.                    &
     &      (local_sum_yend  /=  st_no_data)) THEN

          DO i=local_sum_xstart,local_sum_xend
            IF (partial_SUMMBOT(i)  ==  0.0) THEN
              fieldout(i)=rmdi
            ELSE
              fieldout(i)=partial_SUMMTOP(i)/partial_SUMMBOT(i)
            ENDIF
          ENDDO

        ENDIF ! is this processor in the subdomain

      ENDIF ! does the subdomain intersect with this processor
!           ! column

#endif
!L
  999 CONTINUE
      RETURN
      END SUBROUTINE STMERM
#endif
