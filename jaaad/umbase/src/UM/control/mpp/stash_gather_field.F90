#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(UTILIO) || defined(FLDIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Gathers STASHed data from many processors to one processor
!
! Subroutine interface:
      SUBROUTINE STASH_GATHER_FIELD (                                   &
     &  LOCAL_FIELD , GLOBAL_FIELD ,                                    &
     &  LOCAL_SIZE, GLOBAL_SIZE, LEVELS,                                &
     &  GLOBAL_NORTH , GLOBAL_EAST_IN , GLOBAL_SOUTH , GLOBAL_WEST,     &
     &  GRIDTYPE_CODE ,HALO_TYPE,                                       &
     &  GATHER_PE,                                                      &
     &  DATA_EXTRACTED,                                                 &
     &  PACKING, IM_IDENT, LRLE, PACKING_TYPE,                          &
     &  NUM_OUT,                                                        &
     &  COMP_ACCRCY, loc_RMDI,                                          &
     &  ICODE, CMESSAGE)

      IMPLICIT NONE

! Description:
! Takes a decomposed, STASH processed field and gathers
! it to a single processor, ready for I/O,
!
! Method:
! See in-line documentation
!
! Current code owner : P.Selwood
!
! Subroutine arguments:


      INTEGER, INTENT(IN) ::                                            &
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
     &, GATHER_PE       ! IN: the PE to gather the global field to

      INTEGER, INTENT(OUT) ::                                           &
     &  ICODE           ! OUT: return code, 0=OK
!
! Optional Arguments to handle the COEX packing if necessary
!
      LOGICAL, INTENT(IN), OPTIONAL ::                                  &
     &  PACKING                                                         &
                        ! IN: Set .true. if packing of the input
                        !     field is to be packed
     &, LRLE            ! IN: True if Run Length Encoding is required

      INTEGER, INTENT(IN), OPTIONAL ::                                  &
     &  IM_IDENT        ! IN: Internal model identifier

      INTEGER, INTENT(INOUT), OPTIONAL ::                               &
     &  PACKING_TYPE    ! IN/OUT: This flag is zero on input,
                        !         then stash packing is selected,
                        !         and the routine returns the
                        !         packing flag.
                        !
                        !         If the variable is set to 1 on input
                        !         then 32-bit packing for dumpfiles
                        !         is selected

      INTEGER, INTENT(OUT), OPTIONAL ::                                 &
     &  NUM_OUT         ! OUT: Number of 32-bit IBM words in the Packed
                        !      field for WDGOS packing

      INTEGER, INTENT(IN), OPTIONAL ::                                  &
     &  COMP_ACCRCY     ! IN: Packing Accuracy in Power of 2

      REAL, INTENT(IN), OPTIONAL ::                                     &
     &  loc_RMDI        ! IN: Missing data indicator
!
! Remaining Non-Optional Arguments
!
      LOGICAL, INTENT(IN) ::                                            &
     &  DATA_EXTRACTED  ! IN: TRUE if the data in LOCAL_FIELD has
                        !     already been extracted, or FALSE if
                        !     the extraction must be done here.

      REAL, INTENT(IN) ::                                               &
     &  LOCAL_FIELD(LOCAL_SIZE,LEVELS)
                        ! IN : local data

      REAL, INTENT(OUT) ::                                              &
     &  GLOBAL_FIELD(GLOBAL_SIZE,LEVELS)
                        ! OUT : (PE GATHER_PE only) - full gathered
                        !       field

      CHARACTER*(*), INTENT(OUT) ::                                     &
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
                          ! first row to send from procid
     &, local_start_col                                                 &
                          ! first column to send from procid
     &, sendsize_x                                                      &
                          ! number of points on each row to send
     &, nrows_to_send                                                   &
                          ! number of rows to send from procid
     &, local_row_length                                                &
                          ! size of sending array EW
     &, global_start_row                                                &
                          ! first row to receive at on GATHER_PE
     &, global_start_col                                                &
                          ! first col. to recieve on GATHER_PE
     &, global_row_length                                               &
                          ! size of receiving array EW
     &, flag,info         ! GCOM arguments

! Copies of arguments / variables used to decide if we can use the
! send/receive maps used in the last call

      INTEGER                                                           &
     &  old_LOCAL_SIZE , old_GLOBAL_SIZE                                &
     &, old_GLOBAL_NORTH , old_GLOBAL_EAST_IN                           &
     &, old_GLOBAL_SOUTH , old_GLOBAL_WEST                              &
     &, old_GRIDTYPE_CODE , old_GATHER_PE                               &
     &, old_current_decomp_type, old_HALO_TYPE

      INTEGER                                                           &
! variables defining send and receive maps to be passed to
! GC_RALL_TO_ALL, defining the data transposition
     &  send_map(7,2)                                                   &
     &, receive_map(7,2*MAXPROC)                                        &
     &, n_sends,n_recvs  ! number of sends and receives


      LOGICAL                                                           &
     &  wrap                                                            &
              ! if the subdomain wraps over 0 degree meridion
     &, wrap_special                                                    &
                     ! if there is a wrap around, which starts and
!                      ends on the same processor
     &, zonal_data                                                      &
                     ! if this is a zonal data grid
     &, fullfield                                                       &
                     ! if this is a full field - NOT a subarea
     &, l_packing    ! if packing of data is required

! Save all the variables that may be used in the next call
      SAVE                                                              &
     &  old_LOCAL_SIZE , old_GLOBAL_SIZE                                &
     &, old_GLOBAL_NORTH , old_GLOBAL_EAST_IN                           &
     &, old_GLOBAL_SOUTH , old_GLOBAL_WEST                              &
     &, old_GRIDTYPE_CODE , old_GATHER_PE                               &
     &, old_current_decomp_type                                         &
     &, send_map,receive_map,n_sends,n_recvs                            &
     &, old_HALO_TYPE

! Set all the old_* variables to a number indicating they've
! not been used yet

      DATA                                                              &
     &  old_LOCAL_SIZE , old_GLOBAL_SIZE                                &
     &, old_GLOBAL_NORTH , old_GLOBAL_EAST_IN                           &
     &, old_GLOBAL_SOUTH , old_GLOBAL_WEST                              &
     &, old_GRIDTYPE_CODE , old_GATHER_PE                               &
     &, old_current_decomp_type, old_HALO_TYPE                          &
     &  / -1,-1,-1,-1,-1,-1,-1,-1,-1,-1 /

! Functions

      INTEGER GET_FLD_TYPE
! ------------------------------------------------------------------

      ICODE=0
      IF (PRESENT(PACKING)) THEN
        L_PACKING = PACKING
      ELSE
        L_PACKING = .FALSE.
      END IF

! DEPENDS ON: get_fld_type
       fld_type=GET_FLD_TYPE(GRIDTYPE_CODE)
! See if there is wrap around over meridion, and if so make
! sure that GLOBAL_EAST is > glsize(1)

      GLOBAL_EAST=GLOBAL_EAST_IN
       IF (GLOBAL_EAST  >   glsize(1,fld_type)) THEN
        wrap=.TRUE.
      ELSEIF (GLOBAL_EAST  <   GLOBAL_WEST) THEN
        wrap=.TRUE.
        GLOBAL_EAST=GLOBAL_EAST_IN+glsize(1,fld_type)
      ELSE
        wrap=.FALSE.
      ENDIF

      IF ((GRIDTYPE_CODE  ==  ppx_atm_tzonal) .OR.                      &
                                                   ! Atmos T zonal
     &   ( GRIDTYPE_CODE  ==  ppx_atm_uzonal) .OR.                      &
                                                   ! Atmos U zonal
     &   ( GRIDTYPE_CODE  ==  ppx_ocn_tzonal) .OR.                      &
                                                   ! Ocean T zonal
     &   ( GRIDTYPE_CODE  ==  ppx_ocn_uzonal))                          &
                                                   ! Atmos U zonal
     &  THEN

! This is a zonal field

        zonal_data=.TRUE.
        global_x=1

        IF ((GRIDTYPE_CODE  ==  ppx_atm_tzonal) .OR.                    &
                                                     ! Atmos T zonal
     &      ( GRIDTYPE_CODE  ==  ppx_ocn_tzonal))                       &
                                                     ! Ocean T zonal
     &  THEN
          fld_type=fld_type_p
        ELSE
          fld_type=fld_type_u
        ENDIF
      ELSE

! This is a normal field

        zonal_data=.FALSE.
        global_x=glsize(1,fld_type)
      ENDIF

      global_y=glsize(2,fld_type)

! Set up logical indicating if this is a full field, or just
! a subdomain

      IF (zonal_data) THEN

        fullfield= ( ( GLOBAL_SOUTH  ==  1) .AND.                       &
     &             ( GLOBAL_NORTH  ==  global_y))

      ELSE

        fullfield = (( GLOBAL_WEST  ==  1) .AND.                        &
     &               ( GLOBAL_EAST  ==  global_x) .AND.                 &
     &               ( GLOBAL_SOUTH  ==  1) .AND.                       &
     &               ( GLOBAL_NORTH  ==  global_y))

      ENDIF

! Dealing with fields not in model grid

      if((global_x == 0).or.(global_x == imdi)) then
        write(6,*)'local_size=',local_size
        write(6,*)'global_size=',global_size
        do level=1,levels
          do i=1,global_size
            global_field(i,level)=local_field(i,level)
          enddo
        enddo
      else

! If this a fullfield, we can simply use the standard
! GATHER_FIELD routine

      IF (fullfield) THEN

        IF (zonal_data) THEN

! DEPENDS ON: gather_zonal_field
          CALL GATHER_ZONAL_FIELD( LOCAL_FIELD,GLOBAL_FIELD,            &
     &                          lasize(2,fld_type,halo_type),global_y,  &
     &                            LEVELS,GRIDTYPE_CODE,fld_type,        &
     &                            halo_type,GATHER_PE)

        ELSE

          DO level=1,LEVELS

            IF (L_PACKING) THEN
! DEPENDS ON: gather_pack_field
              CALL GATHER_PACK_FIELD(                                   &
     &                        LOCAL_FIELD(1,level),                     &
     &                        GLOBAL_FIELD(1,level),                    &
     &                        lasize(1,fld_type,halo_type),             &
     &                        lasize(2,fld_type,halo_type),             &
     &                        global_x,global_y,                        &
     &                        fld_type,halo_type,                       &
     &                        GATHER_PE,GC_ALL_PROC_GROUP,              &
     &                        PACKING, IM_IDENT, LRLE, PACKING_TYPE,    &
     &                        NUM_OUT,                                  &
     &                        COMP_ACCRCY, loc_RMDI)
            ELSE
! DEPENDS ON: gather_field
              CALL GATHER_FIELD( LOCAL_FIELD(1,level),                  &
     &                           GLOBAL_FIELD(1,level),                 &
     &                           lasize(1,fld_type,halo_type),          &
     &                           lasize(2,fld_type,halo_type),          &
     &                           global_x, global_y,                    &
     &                           fld_type, halo_type,                   &
     &                           GATHER_PE, GC_ALL_PROC_GROUP,          &
     &                           ICODE, CMESSAGE)
            END IF

            IF (ICODE  /=  0) THEN
              WRITE(6,*)                                                &
     &         'STASH_GATHER_FIELD : Failed during GATHER_FIELD'
              WRITE(6,*) 'Error code : ',ICODE
              WRITE(6,*) 'Message : ',CMESSAGE

              ICODE=2
              CMESSAGE='Failed to gather field'
              GOTO 9999
            ENDIF

           ENDDO

         ENDIF
       ELSE
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
     &    (GATHER_PE  ==  old_GATHER_PE) .AND.                          &
     &    (HALO_TYPE  ==  old_HALO_TYPE) .AND.                          &
     &    (current_decomp_type  ==  old_current_decomp_type ))) THEN

          old_LOCAL_SIZE=LOCAL_SIZE
          old_GLOBAL_SIZE=GLOBAL_SIZE
          old_GLOBAL_NORTH=GLOBAL_NORTH
          old_GLOBAL_EAST_IN=GLOBAL_EAST_IN
          old_GLOBAL_SOUTH=GLOBAL_SOUTH
          old_GLOBAL_WEST=GLOBAL_WEST
          old_GRIDTYPE_CODE=GRIDTYPE_CODE
          old_GATHER_PE=GATHER_PE
          old_current_decomp_type=current_decomp_type
          old_HALO_TYPE=HALO_TYPE

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

              IF (DATA_EXTRACTED) THEN
                local_start_row=1
              ELSE
                local_start_row=local_ystart
              ENDIF
              nrows_to_send=local_yend-local_ystart+1

             global_start_row=g_datastart(2,procid)+local_ystart-       &
     &                         halosize(2,halo_type) - GLOBAL_SOUTH
              global_row_length=GLOBAL_EAST-GLOBAL_WEST+1

! Calculate the following variables:
! local_row_length : the X dimension size of the local array
! local_send_offx  : the offset into each row to start sending from
! sendsize_x       : the number of points on each row to send
! The calculation of these numbers is different for processors
! at the start and end of a wrap_special case
! Note that when DATA_EXTRACTED is true, then local_field has no
! halos.

              IF (wrap_special .AND. procx  ==  proc_topleft_x) THEN
                IF (DATA_EXTRACTED) THEN
      local_row_length = g_lasize(1,fld_type,halo_type_no_halo,procid)  &
     &                   + local_xend - local_xstart + 1
      sendsize_x       = g_lasize(1,fld_type,halo_type_no_halo,procid)  &
     &                   - local_xstart + 1
      local_start_col  = 1
                ELSE
      local_row_length = g_lasize(1,fld_type,halo_type,procid)
      sendsize_x       = g_lasize(1,fld_type,halo_type,procid)          &
     &                   - local_xstart
      local_start_col  = local_xstart

                ENDIF
                global_start_col=1

              ELSEIF (wrap_special .AND. procx  ==  proc_botright_x)    &
     &        THEN
                IF (DATA_EXTRACTED) THEN
      local_row_length = g_lasize(1,fld_type,halo_type_no_halo,procid)  &
     &                   + local_xend - local_xstart + 1
      local_start_col  = local_row_length - local_xend + 1
      sendsize_x       = local_xend
                ELSE
      local_row_length = g_lasize(1,fld_type,halo_type,procid)
      local_start_col  = Offx + 1
      sendsize_x       = local_xend - Offx
                ENDIF
                global_start_col=global_row_length-sendsize_x+1

              ELSE
                IF (DATA_EXTRACTED) THEN
                  local_row_length=local_xend-local_xstart+1
                  local_start_col=1
                ELSE
          local_row_length=g_lasize(1,fld_type,halo_type,procid)
                  local_start_col=local_xstart
                ENDIF
                sendsize_x=local_xend-local_xstart+1
                global_start_col=local_xstart-halosize(1,halo_type)+    &
     &                           g_datastart(1,procid)-GLOBAL_WEST
              ENDIF

              IF (global_start_col  <   0) THEN
! Wrapped around field, but this processor is not start or end
! processor
        global_start_col=global_start_col+glsize(1,fld_type)
              ENDIF

! Now we can set up the send and receive map entries for the data on
! this processor

              IF (mype  ==  procid) THEN  ! I need to send some data


                n_sends=n_sends+1

                send_map(S_DESTINATION_PE,n_sends) = GATHER_PE
                send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,n_sends) =        &
     &            (local_start_row-1)*local_row_length +                &
     &            local_start_col
                send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,n_sends) =        &
     &            nrows_to_send
                send_map(S_STRIDE_IN_SEND_ARRAY,n_sends) =              &
     &            local_row_length
                send_map(S_ELEMENT_LENGTH,n_sends) = sendsize_x
                send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,n_sends) =        &
     &            (global_start_row-1)*global_row_length +              &
     &            global_start_col
                send_map(S_STRIDE_IN_RECV_ARRAY,n_sends) =              &
     &            global_row_length

              ENDIF ! if I'm sending data

              IF (mype  ==  GATHER_PE) THEN  ! I need to receive data

                n_recvs=n_recvs+1

                receive_map(R_SOURCE_PE,n_recvs) = procid
                receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,n_recvs) =     &
     &            (global_start_row-1)*global_row_length +              &
     &            global_start_col
                receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,n_recvs) =     &
     &            nrows_to_send
                receive_map(R_STRIDE_IN_RECV_ARRAY,n_recvs) =           &
     &            global_row_length
                receive_map(R_ELEMENT_LENGTH,n_recvs) = sendsize_x
                receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,n_recvs) =     &
     &            (local_start_row-1)*local_row_length +                &
     &            local_start_col
                receive_map(R_STRIDE_IN_SEND_ARRAY,n_recvs) =           &
     &            local_row_length

              ENDIF ! if I'm receiving data

            ENDDO ! procx : loop along processor row

          ENDDO ! procy : loop down processor column

        ENDIF ! if I need to recalculate my send/receive maps

! Send / receive the data using GCG_RALLTOALLE


        DO level=1,LEVELS

          flag=0  ! This is currently ignored at GCG v1.1
          info=GC_NONE

          CALL GCG_RALLTOALLE(                                          &
     &      LOCAL_FIELD(1,level)  ,                                     &
     &      send_map    , n_sends  ,LOCAL_SIZE  ,                       &
     &      GLOBAL_FIELD(1,level) ,                                     &
     &      receive_map , n_recvs , GLOBAL_SIZE ,                       &
     &      GC_ALL_PROC_GROUP , flag, info)

        ENDDO

        ENDIF ! if this is a full or extracted field

      endif


 9999 CONTINUE

      RETURN
      END SUBROUTINE STASH_GATHER_FIELD

#endif
