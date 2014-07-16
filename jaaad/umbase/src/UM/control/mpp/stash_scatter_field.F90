#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(UTILIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Scatters STASHed data from one processor to many processors
!
! Subroutine interface:
      SUBROUTINE STASH_SCATTER_FIELD (                                  &
     &  LOCAL_FIELD , GLOBAL_FIELD ,                                    &
     &  LOCAL_SIZE, GLOBAL_SIZE, LEVELS,                                &
     &  GLOBAL_NORTH , GLOBAL_EAST_IN , GLOBAL_SOUTH , GLOBAL_WEST,     &
     &  GRIDTYPE_CODE, HALO_TYPE,                                       &
     &  SCATTER_PE,                                                     &
     &  ICODE, CMESSAGE)

      IMPLICIT NONE


! Description:
! Takes a decomposed, STASH processed field and gathers
! it to a single processor, ready for I/O,
!
! Method:
! See in-line documentation
!
! Current code owner : P.Burton
!
!
! Subroutine arguments:

      INTEGER                                                           &
     &  LOCAL_SIZE                                                      &
                        ! IN: size of level of LOCAL_FIELD
     &, GLOBAL_SIZE                                                     &
                        ! IN: size of level of GLOBAL_FIELD
     &, LEVELS                                                          &
                        ! IN: number of levels
     &, GLOBAL_NORTH                                                    &
                        ! IN: specification of subdomain boundaries
     &, GLOBAL_EAST_IN                                                  &
                        ! IN: ""
     &, GLOBAL_SOUTH                                                    &
                        ! IN: ""
     &, GLOBAL_WEST                                                     &
                        ! IN: ""
     &, GRIDTYPE_CODE                                                   &
                        ! IN: indicates the type of grid output
     &, HALO_TYPE                                                       &
                        ! IN: type of halo on this field
     &, SCATTER_PE                                                      &
                        ! IN: the PE to scatter global field from
     &, ICODE           ! OUT: return code, 0=OK

      REAL                                                              &
     &  LOCAL_FIELD(LOCAL_SIZE,LEVELS)                                  &
!        ! OUT : local scattered data
     &, GLOBAL_FIELD(GLOBAL_SIZE,LEVELS)
!        ! IN : (PE SCATTER_PE only) - full field

      CHARACTER*80                                                      &
     &  CMESSAGE        ! OUT: Error message if ICODE  /=  0

! Parameters and common blocks
#include "stparam.h"
#include "cppxref.h"
#include "parvars.h"
#include "gccom.h"
#include "c_mdi.h"

! Local variables

      INTEGER                                                           &
     &  GLOBAL_EAST                                                     &
                      ! copy of GLOBAL_EAST_IN with wrap around s.t.
!                     ! GLOBAL_EAST > GLOBAL_ROW_LEN
     &, global_x                                                        &
                      ! size of global data EW
     &, global_y                                                        &
                      ! size of global data NS
     &, fld_type                                                        &
                      ! indicates if field is on P or U grid
     &, level                                                           &
                      ! loop index for loop over levels
     &, i                                                               &
                      ! loop counter
     &, proc_topleft_x,proc_topleft_y                                   &
                                        ! processors at corners of
     &, proc_botright_x,proc_botright_y                                 &
                                        ! the subarea
     &, dummy1,dummy2                                                   &
                      ! ignored return arguments
     &, procx,procy                                                     &
                      ! loop indexes for loops over processors
     &, eff_procx                                                       &
                      ! real x co-ord of processor column procx
     &, procid                                                          &
                      ! processor id of (procx,procy)
     &, local_xstart,local_xend                                         &
                                 ! boundaries of subdomain for
     &, local_ystart,local_yend                                         &
                                 ! processor procid
     &, local_start_row                                                 &
                          ! first row to receive on procid
     &, local_start_col                                                 &
                          ! first column to receive on procid
     &, sendsize_x                                                      &
                          ! number of points on each row to send
!                         ! to procid
     &, nrows_to_send                                                   &
                          ! number of rows to send to procid
     &, local_row_length                                                &
                          ! size of receiving array EW
     &, global_start_row                                                &
                          ! first row to send on SCATTER_PE
     &, global_start_col                                                &
                          ! first col. to send on SCATTER_PE
     &, global_row_length                                               &
                          ! size of sending array EW
     &, flag,info         ! GCOM arguments

      LOGICAL                                                           &
     &  l_vec             ! Indicates if a field is a vector quantity
                          ! (SWAPBOUNDS argument)
! Copies of arguments / variables used to decide if we can use the
! send/receive maps used in the last call

      INTEGER                                                           &
     &  old_LOCAL_SIZE , old_GLOBAL_SIZE                                &
     &, old_GLOBAL_NORTH , old_GLOBAL_EAST_IN                           &
     &, old_GLOBAL_SOUTH , old_GLOBAL_WEST                              &
     &, old_GRIDTYPE_CODE , old_SCATTER_PE                              &
     &, old_HALO_TYPE,old_current_decomp_type

      INTEGER                                                           &
! variables defining send and receive maps to be passed to
! GC_RALL_TO_ALL, defining the data transposition
     &  send_map(7,MAXPROC*2)                                           &
     &, receive_map(7,2)                                                &
     &, n_sends,n_recvs  ! number of sends and receives


      LOGICAL                                                           &
     &  wrap                                                            &
              ! if the subdomain wraps over 0 degree meridion
     &, wrap_special                                                    &
                     ! if there is a wrap around, which starts and
!                      ends on the same processor
     &, zonal_data                                                      &
                     ! if this is a zonal data grid
     &, fullfield    ! if this is a full field - NOT a subarea

! Save all the variables that may be used in the next call
      SAVE                                                              &
     &  old_LOCAL_SIZE , old_GLOBAL_SIZE                                &
     &, old_GLOBAL_NORTH , old_GLOBAL_EAST_IN                           &
     &, old_GLOBAL_SOUTH , old_GLOBAL_WEST                              &
     &, old_GRIDTYPE_CODE , old_SCATTER_PE                              &
     &, old_HALO_TYPE,old_current_decomp_type                           &
     &, send_map,receive_map,n_sends,n_recvs

! Set all the old_* variables to a number indicating they've
! not been used yet

      DATA                                                              &
     &  old_LOCAL_SIZE , old_GLOBAL_SIZE                                &
     &, old_GLOBAL_NORTH , old_GLOBAL_EAST_IN                           &
     &, old_GLOBAL_SOUTH , old_GLOBAL_WEST                              &
     &, old_GRIDTYPE_CODE , old_SCATTER_PE                              &
     &, old_HALO_TYPE,old_current_decomp_type                           &
     &  / -1,-1,-1,-1,-1,-1,-1,-1,-1,-1 /

! Functions

      INTEGER GET_FLD_TYPE
! ------------------------------------------------------------------

! DEPENDS ON: get_fld_type
      fld_type=GET_FLD_TYPE(GRIDTYPE_CODE)

      l_vec=((fld_type  ==  fld_type_u) .OR.                            &
     &       (fld_type  ==  fld_type_v))

! See if there is wrap around over meridion, and if so make
! sure that GLOBAL_EAST is > glsize(1)

      GLOBAL_EAST=GLOBAL_EAST_IN
      IF (fld_type  ==  fld_type_unknown) THEN
        WRITE(6,*)                                                      &
     &  'STASH_SCATTER_FIELD cannot process field with ppx gridtype ',  &
     &  GRIDTYPE_CODE
        ICODE=1
        CMESSAGE='STASH_SCATTER_FIELD : Incompatible GRIDTYPE code'
        GOTO 9999
      ENDIF

      IF (GLOBAL_EAST  >   glsize(1,fld_type)) THEN
        wrap=.TRUE.
      ELSEIF (GLOBAL_EAST  <   GLOBAL_WEST) THEN
        wrap=.TRUE.
        GLOBAL_EAST=GLOBAL_EAST_IN+glsize(1,fld_type)
      ELSE
        wrap=.FALSE.
      ENDIF


      IF ((GRIDTYPE_CODE  ==  ppx_atm_tzonal) .OR.                      &
     &    (GRIDTYPE_CODE  ==  ppx_ocn_uzonal) .OR.                      &
     &    (GRIDTYPE_CODE  ==  ppx_ocn_tzonal)) THEN

        zonal_data=.TRUE.
        global_x=1
      ELSE

! This is a normal field

        zonal_data=.FALSE.
        global_x=glsize(1,fld_type)

      ENDIF

      global_y=glsize(2,fld_type)

! Set up logical indicating if this is a full field, or just
! a subdomain

      IF (zonal_data) THEN

        fullfield= ( ( GLOBAL_NORTH  ==  global_y) .AND.                &
     &             ( GLOBAL_SOUTH  ==  1))

      ELSE

        fullfield = (( GLOBAL_WEST  ==  1) .AND.                        &
     &               ( GLOBAL_EAST  ==  global_x) .AND.                 &
     &               ( GLOBAL_NORTH  ==  global_y) .AND.                &
     &               ( GLOBAL_SOUTH  ==  1))

      ENDIF

! Dealing with fields not in model grid

      if((global_x == 0).or.(global_x == imdi)) then
        do level=1,levels
          do i=1,global_size
            local_field(i,level)=global_field(i,level)
          enddo
        enddo
      else

! If this is a fullfield, we can simply use the standard
! SCATTER_FIELD routine

      IF (fullfield) THEN

        IF (zonal_data) THEN

! DEPENDS ON: scatter_zonal_field
          CALL SCATTER_ZONAL_FIELD( LOCAL_FIELD,GLOBAL_FIELD,           &
     &                              lasize(2,fld_type,halo_type),       &
     &                              global_y,                           &
     &                              LEVELS,GRIDTYPE_CODE,fld_type,      &
     &                              halo_type,                          &
     &                              SCATTER_PE)

! Don't call swapbounds for ocean zonal fields which currently
! do not have halos

          IF ((GRIDTYPE_CODE  /=  ppx_ocn_tzonal).AND.                  &
     &        (GRIDTYPE_CODE  /=  ppx_ocn_uzonal)) THEN
! DEPENDS ON: swap_bounds
            CALL SWAP_BOUNDS(LOCAL_FIELD,                               &
     &                      1,blsize(2,fld_type),LEVELS,                &
     &                      0,halosize(2,halo_type),                    &
     &                      fld_type,l_vec)
          ENDIF

        ELSE

          DO level=1,LEVELS

! DEPENDS ON: scatter_field
            CALL SCATTER_FIELD( LOCAL_FIELD(1,level) ,                  &
     &                          GLOBAL_FIELD(1,level),                  &
     &                          lasize(1,fld_type,halo_type),           &
     &                          lasize(2,fld_type,halo_type),           &
     &                          global_x,global_y,                      &
     &                          fld_type,halo_type,                     &
     &                          SCATTER_PE,GC_ALL_PROC_GROUP,           &
     &                          ICODE,CMESSAGE)

            IF (ICODE  /=  0) THEN
              WRITE(6,*) 'STASH_SCATTER_FIELD : SCATTER_FIELD failed'
              WRITE(6,*) 'Return code was : ',ICODE
              WRITE(6,*) 'Error message was : ',CMESSAGE
              CMESSAGE='SCATTER_FIELD failed'
              GOTO 9999
            ENDIF


          ENDDO

! DEPENDS ON: swap_bounds
          CALL SWAP_BOUNDS(LOCAL_FIELD,                                 &
     &                    blsize(1,fld_type),                           &
     &                    blsize(2,fld_type),                           &
     &                    LEVELS,                                       &
     &                    halosize(1,halo_type),                        &
     &                    halosize(2,halo_type),                        &
     &                    fld_type,l_vec)

         ENDIF
       ELSE
!
! Check for WAM Wave fields - no code exists now
!
        if(fld_type  ==  fld_type_comp_wave .or.                        &
     &     fld_type  ==  fld_type_full_wave) then
          write(6,'(''STASH_SCATTER_FIELD:'',                           &
     &     '' No Code Exists to Generate Sub-domains for the'',         &
     &     '' WAM Wave Model'')')
          WRITE(6,*) 'Unable to process this field.'
          cmessage='STASH_SCATTER_FIELD:'//                             &
     &     ' Unable to process this field.'
          icode=100
          goto 9999
        endif
!
! for subdomains, life is not so easy - we must explicitly
! calculate our own send and receive maps, and use GCG_RALLTOALLE
! to shift the data around.

! If the same arguments are used as were used in the last call
! to this routine, we can just use the previously calculated
! send and receive maps, otherwise we need to calculate new maps

        IF (.NOT. (                                                     &
     &    (LOCAL_SIZE  ==  old_LOCAL_SIZE) .AND.                        &
     &    (GLOBAL_SIZE  ==  old_GLOBAL_SIZE) .AND.                      &
     &    (GLOBAL_NORTH  ==  old_GLOBAL_NORTH) .AND.                    &
     &    (GLOBAL_EAST_IN  ==  old_GLOBAL_EAST_IN) .AND.                &
     &    (GLOBAL_SOUTH  ==  old_GLOBAL_SOUTH) .AND.                    &
     &    (GLOBAL_WEST  ==  old_GLOBAL_WEST) .AND.                      &
     &    (GRIDTYPE_CODE  ==  old_GRIDTYPE_CODE) .AND.                  &
     &    (HALO_TYPE  ==  old_HALO_TYPE) .AND.                          &
     &    (SCATTER_PE  ==  old_SCATTER_PE) .AND.                        &
     &    (current_decomp_type  ==  old_current_decomp_type ))) THEN

          old_LOCAL_SIZE=LOCAL_SIZE
          old_GLOBAL_SIZE=GLOBAL_SIZE
          old_GLOBAL_NORTH=GLOBAL_NORTH
          old_GLOBAL_EAST_IN=GLOBAL_EAST_IN
          old_GLOBAL_SOUTH=GLOBAL_SOUTH
          old_GLOBAL_WEST=GLOBAL_WEST
          old_GRIDTYPE_CODE=GRIDTYPE_CODE
          old_HALO_TYPE=HALO_TYPE
          old_SCATTER_PE=SCATTER_PE
          old_current_decomp_type=current_decomp_type

! Find out what the boundaries of the subdomain area

! DEPENDS ON: global_to_local_rc
          CALL GLOBAL_TO_LOCAL_RC(GRIDTYPE_CODE,HALO_TYPE,              &
     &                            GLOBAL_WEST,GLOBAL_NORTH,             &
     &                            proc_topleft_x,proc_topleft_y,        &
     &                            dummy1,dummy2)
! DEPENDS ON: global_to_local_rc
          CALL GLOBAL_TO_LOCAL_RC(GRIDTYPE_CODE,HALO_TYPE,              &
     &                            GLOBAL_EAST,GLOBAL_SOUTH,             &
     &                            proc_botright_x,proc_botright_y,      &
     &                            dummy1,dummy2)

! Ensure that the processor x co-ords are such that the botright_x is
! always greater than (or equal to) top_left_x.
          IF (wrap) proc_botright_x=gridsize(1)+proc_botright_x

! wrap_special is set to true if there is a wrap around which starts
! and ends on the same processor. This case requires extra work as
! the processor in question
          IF (wrap .AND. (proc_topleft_x+gridsize(1)  ==                &
     &                    proc_botright_x)) THEN
            wrap_special=.TRUE.
          ELSE
            wrap_special=.FALSE.
          ENDIF

          n_sends=0
          n_recvs=0

          DO procy=proc_botright_y,proc_topleft_y
            DO procx=proc_topleft_x,proc_botright_x

              eff_procx=MOD(procx,gridsize(1))
              procid=eff_procx+procy*gridsize(1)

! DEPENDS ON: global_to_local_subdomain
              CALL GLOBAL_TO_LOCAL_SUBDOMAIN(                           &
     &          .TRUE.,.TRUE.,                                          &
     &          GRIDTYPE_CODE,HALO_TYPE,procid,                         &
     &          GLOBAL_SOUTH,GLOBAL_EAST,                               &
     &          GLOBAL_NORTH,GLOBAL_WEST,                               &
     &          local_ystart,local_xend,                                &
     &          local_yend  ,local_xstart)

! Calculate the shape of the arrays, and where to start sending/
! receiving data, and how many rows to send

              local_start_row=1
              nrows_to_send=local_yend-local_ystart+1

              global_start_row=g_datastart(2,procid)+local_ystart -     &
     &                         halosize(2,halo_type) -                  &
     &                         GLOBAL_SOUTH
              global_row_length=GLOBAL_EAST-GLOBAL_WEST+1

! Calculate the following variables:
! local_row_length : the X dimension size of the local array
! local_send_offx  : the offset into each row to start sending from
! sendsize_x       : the number of points on each row to send
! The calculation of these numbers is different for processors
! at the start and end of a wrap_special case

              IF (wrap_special .AND. procx  ==  proc_topleft_x) THEN
                local_row_length=g_blsize(1,fld_type,procid) +          &
     &                           local_xend - local_xstart + 1
                local_start_col=1
                sendsize_x=g_lasize(1,fld_type,halo_type,procid) -      &
     &                     local_xstart
                global_start_col=1

              ELSEIF (wrap_special .AND. procx  ==  proc_botright_x)    &
     &        THEN
                local_row_length=g_blsize(1,fld_type,procid) +          &
     &                           local_xend - local_xstart + 1
                local_start_col=local_row_length - local_xend +         &
     &                          halosize(1,halo_type) + 1
                sendsize_x=local_xend - halosize(1,halo_type)
                global_start_col=global_row_length-sendsize_x+1

              ELSE
                local_row_length=local_xend-local_xstart+1
                local_start_col=1
                sendsize_x=local_xend-local_xstart+1
                global_start_col=local_xstart -                         &
     &                           (halosize(1,halo_type) + 1 ) +         &
     &                           g_datastart(1,procid)-GLOBAL_WEST+1
              ENDIF

              IF (global_start_col  <   0) THEN
! Wrapped around field, but this processor is not start or end
! processor
                global_start_col=global_start_col+glsize(1,fld_type)
              ENDIF

! Now we can set up the send and receive map entries for the data on
! this processor

              IF (mype  ==  procid) THEN  ! I need to receive some data

                  n_recvs=n_recvs+1

                receive_map(R_SOURCE_PE,n_recvs) = SCATTER_PE
                receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,n_recvs) =     &
     &            (local_start_row-1)*local_row_length +                &
     &            local_start_col
                receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,n_recvs) =     &
     &            nrows_to_send
                receive_map(R_STRIDE_IN_RECV_ARRAY,n_recvs) =           &
     &            local_row_length
                receive_map(R_ELEMENT_LENGTH,n_recvs) = sendsize_x
                receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,n_recvs) =     &
     &            (global_start_row-1)*global_row_length +              &
     &            global_start_col
                receive_map(R_STRIDE_IN_SEND_ARRAY,n_recvs) =           &
     &            global_row_length

              ENDIF ! if I'm receiving data

              IF (mype  ==  SCATTER_PE) THEN ! I need to send data

                n_sends=n_sends+1

                send_map(S_DESTINATION_PE,n_sends) = procid
                send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,n_sends) =        &
     &            (global_start_row-1)*global_row_length +              &
     &            global_start_col
                send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,n_sends) =        &
     &            nrows_to_send
                send_map(S_STRIDE_IN_SEND_ARRAY,n_sends) =              &
     &            global_row_length
                send_map(S_ELEMENT_LENGTH,n_sends) = sendsize_x
                send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,n_sends) =        &
     &            (local_start_row-1)*local_row_length +                &
     &            local_start_col
                send_map(S_STRIDE_IN_RECV_ARRAY,n_sends) =              &
     &            local_row_length

              ENDIF ! if I'm sending data

            ENDDO ! procx : loop along processor row

          ENDDO ! procy : loop down processor column

        ENDIF ! if I need to recalculate my send/receive maps

! Send / receive the data using GCG_RALLTOALLE

        DO level=1,LEVELS

          flag=0  ! This is currently ignored at GCG v1.1
          info=GC_NONE

          CALL GCG_RALLTOALLE(                                          &
     &      GLOBAL_FIELD(1,level)  ,                                    &
     &      send_map    , n_sends  ,GLOBAL_SIZE  ,                      &
     &      LOCAL_FIELD(1,level) ,                                      &
     &      receive_map , n_recvs , LOCAL_SIZE ,                        &
     &      GC_ALL_PROC_GROUP , flag, info)

        ENDDO

      ENDIF ! if this is a full or extracted field

      endif

 9999 CONTINUE

      RETURN
      END SUBROUTINE STASH_SCATTER_FIELD

#endif
