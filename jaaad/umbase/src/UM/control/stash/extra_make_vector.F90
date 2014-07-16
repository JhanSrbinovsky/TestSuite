#if defined(C84_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL -----------------------------------------------------------
!LL Stash routine
!LL purpose: Generate extra data for the timeseries.
!LL This extra data provides information about what processing was done
!LL to produce the timeseries. This information will hopefully be of som
!LL use to users doing further processing of the timeseries data.
!LL This deck contains two subroutines
!LL (1) EXTRA_TS_INFO : which generates the codes and sets up the space
!LL                   : for the extra data.
!LL
!LL (2) EXTRA_MAKE_VECTOR: which computes the long/latt ht domain info
!LL                   : and puts that into the correct place in the
!LL                   : extra data
!LL Routines are  called by stmulspa1.
!LL
!LL To some extent this routine has much in common with the
!LL multi_spatial routine but as it has a different function
!LL viz generate info on timeseries rather than generating a single time
!LL for the timeseries it is coded separately.
!LL when modifying multi_spatial be sure also to modify this routine and
!LL vice versa
!LL
!LL 16/3/92 Written by Simon Tett
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL
!LL Programming Standard: UM DOC Paper3, Verion 4 (05/02/92)
!LL
!LL System Component Covered: D711
!LL
!LL System Task:C4
!LL

!*L Interface and arguments ------------------------------------
!


!*L Interface and arguments: -----------------------------
      SUBROUTINE EXTRA_MAKE_VECTOR(control,control_len,record_cnt,      &
     &  no_records,extra_data,extra_data_len,                           &
     &   bzx,bzy,bdx,bdy)
      implicit none
      integer control_len ! IN size of control record
      integer control(control_len) ! IN stash control record
      integer record_cnt ! IN record that is being processed
      integer no_records ! IN total number of records
      integer extra_data_len ! IN size of extra data
      real extra_data(extra_data_len) !IN/OUT extra data
      real bdx,bdy,bzx,bzy ! IN grid descriptors
!* ------------------------------------------------------
! Parameters
!
#include "stparam.h"
#include "sterr.h"
!*L
! Subroutines called: none
!
!*L Local variables
      integer addr ! what address in extra data are we at
      integer record_len ! how many words in a block ?
!*L
      record_len=no_records+1
      addr=1+record_cnt
!L put in the first latitude
      extra_data(addr)=control(st_south_code)*bdy+bzy
      addr=addr+record_len
!L put in the first long
      extra_data(addr)=control(st_west_code)*bdx+bzx
      addr=addr+record_len
!L put in the second lat
      extra_data(addr)=control(st_north_code)*bdy+bzy
      addr=addr+record_len
!L put in the second long
      extra_data(addr)=control(st_east_code)*bdx+bzx
      addr=addr+record_len
!L put in the lowest level
      extra_data(addr)=control(st_input_bottom)
      addr=addr+record_len
!L and now the highest  level
      extra_data(addr)=control(st_input_top)

      RETURN
      END SUBROUTINE EXTRA_MAKE_VECTOR


!LL    Subroutine: STUFF_INT-----------------------------------------
!LL
!LL    Purpose: To put the binary representation of an integer into a
!LL    real variable through hidden equivalencing via argument passing
!LL
!LL  Model              Modification history from model version 3.0:
!LL version  date
!LL   5.3    14/05/01   Changes data type of data_in to integer.
!LL                     E.Leung
!LL   5.4    10/04/02   Reverse the above change. S.D.Mullerworth

#endif
