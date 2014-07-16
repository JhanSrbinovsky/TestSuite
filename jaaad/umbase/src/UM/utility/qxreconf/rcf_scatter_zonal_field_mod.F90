#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! + Scatters zonal field from one processor to many processors

Module Rcf_Scatter_Zonal_Field_Mod

!  Subroutine Rcf_Scatter_Zonal - scatters a zonal field
!
! Description:
! Takes a zonal field on a single processor, and decomposes it over
! many processors.
!
! Method:
!   Derived from UM 4.5 code.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Scatter_Zonal_Field (  LOCAL_FIELD , GLOBAL_FIELD , &
                                      LOCAL_SIZE  , GLOBAL_SIZE  , &
                                      LEVELS, GRID_TYPE , &
                                      SCATTER_PE)

Use Rcf_Parvars_Mod     !  Lots used.

IMPLICIT NONE

! Subroutine arguments:
Integer, Intent(In)   :: Local_Size    ! size of local_field level
Integer, Intent(In)   :: Global_Size   ! size of global_filed level
Integer, Intent(In)   :: Levels        ! number of levels in field
Integer, Intent(In)   :: Grid_Type     ! grid type of field
Integer, Intent(In)   :: Scatter_PE    ! PE on which Global_field lives

                                       ! Field to scatter
Real, Intent(In)      :: Global_Field( Global_Size, Levels )
                                       ! Local part of field
Real, Intent(Out)     :: Local_Field( Local_Size, Levels )

! Comdecks
#include "cppxref.h"
#include "gccom.h"

! Local variables
Integer               :: fld_type            ! P or U
Integer               :: info                ! return from GCOM
Integer               :: send_map(7,MAXPROC) ! send map
Integer               :: receive_map(7,1)    ! receive map
Integer               :: flag                ! dummy arg for GCOM
Integer               :: n_mess_to_send      ! number of messages
Integer               :: k                   ! loop counter
Integer               :: iproc               ! loop counter

Logical               :: mead_fld            ! is ocean field a
                                             ! Mead diagnostic
!====================================================================

IF ((grid_type .EQ. ppx_atm_tzonal) .OR. &
    (grid_type .EQ. ppx_ocn_tzonal)) THEN
  fld_type=fld_type_p
ELSE
  fld_type=fld_type_u
ENDIF

mead_fld=((grid_type .EQ. ppx_ocn_uzonal) .OR. &
          (grid_type .EQ. ppx_ocn_tzonal))

!--------------------------------------------------------------------

n_mess_to_send=0

IF (mype .EQ. SCATTER_PE) THEN
  DO iproc=0,nproc-1
    send_map(S_DESTINATION_PE,iproc+1) = iproc
    send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) = &
             g_datastart(2,iproc)
    send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = 1
    send_map(S_STRIDE_IN_SEND_ARRAY,iproc+1) = 0
    IF (fld_type .EQ. fld_type_p) THEN
      send_map(S_ELEMENT_LENGTH,iproc+1) = g_blsizep(2,iproc)
    ELSE
      send_map(S_ELEMENT_LENGTH,iproc+1) = g_blsizeu(2,iproc)
    ENDIF
    IF (mead_fld) THEN
      send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = 1
    ELSE
      send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = Offy+1
    ENDIF
    send_map(S_STRIDE_IN_RECV_ARRAY,iproc+1) = 0
  ENDDO
  n_mess_to_send=nproc
ENDIF

receive_map(R_SOURCE_PE,1) = SCATTER_PE
IF (mead_fld) THEN
  receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,1) = 1
ELSE
  receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,1) = Offy+1
ENDIF
receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,1) = 1
receive_map(R_STRIDE_IN_RECV_ARRAY,1) = 0
IF (fld_type .EQ. fld_type_p) THEN
  receive_map(R_ELEMENT_LENGTH,1) = blsizep(2)
ELSE
  receive_map(R_ELEMENT_LENGTH,1) = blsizeu(2)
ENDIF
receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,1) = datastart(2)
receive_map(R_STRIDE_IN_SEND_ARRAY,1) = 0

DO k=1,LEVELS

  info=GC_NONE
  flag=GC_NONE

  IF (fld_type .EQ. fld_type_p) THEN
    CALL GCG_RALLTOALLE( GLOBAL_FIELD(1,k),    send_map,    &
                         n_mess_to_send,       glsize(2),   &
                         LOCAL_FIELD(1,k),     receive_map, &
                         1,                    lasize(2),   &
                         GC_ALL_PROC_GROUP,    flag,    info)
  ELSE
    CALL GCG_RALLTOALLE( GLOBAL_FIELD(1,k),     send_map,    &
                         n_mess_to_send,        glsize(2)-1, &
                         LOCAL_FIELD(1,k),      receive_map, &
                         1,                     lasize(2),   &
                         GC_ALL_PROC_GROUP,     flag,    info)
  ENDIF

ENDDO

RETURN

END Subroutine Rcf_Scatter_Zonal_Field
END Module Rcf_Scatter_Zonal_Field_Mod


#endif
