#if defined(A32_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate a pressure field for the source model levels

      Subroutine Calc_P_Src (                                           &
#include "arginfa.h"
     &           ip_p, len_p, len_intf_p,                               &
     &           gather_PE, pe_for_levels, local_level,                 &
     &           ex_p, p_src )

! Description:
!   This routine calculates a pressure field for the source model levels
!   at the LBC points. Required for interpolation of LBC prognostics
!   from model levels to LBC levels.
!
! Method:
!   1. Set up Global fields of Exner_P
!   2. Convert Exner_P to P
!   3. Horizontal interpolation of P from global grid to the LBC grid.
!   4. Distribute P to the PEs.
!
! Owner: Dave Robinson
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

      Implicit None

#include "cmaxsize.h"
#include "parvars.h"
#include "typsize.h"
#include "cintfa.h"
#include "typinfa.h"

      Integer :: ip_p
      Integer :: len_p
      Integer :: len_intf_p
      Integer :: gather_PE
      Integer :: pe_for_levels(model_levels)
      Integer :: local_level(model_levels)

      Real :: ex_p  (len_p     , model_levels)
      Real :: p_src (len_intf_p, model_levels)

#include "gccom.h"
#include "c_r_cp.h"

      INTEGER                                                           &
     &  send_map(7,model_levels)                                        &
     &, recv_map(7,model_levels)

      INTEGER                                                           &
     &  map(model_levels)                                               &
                          ! mapping of processors holding which level
     &, level_on_gat_pe(model_levels)                                   &
                                      ! what is the level number on the
                                      ! gathering PE of this level
     &, level_count(0:nproc-1)                                          &
                                      ! Count of levels on each PE
     &, max_level_count               ! Max #levels on PE

      INTEGER                                                           &
     &  iproc                                                           &
                    ! loop counter over processors
     &, iseg                                                            &
                    ! loop counter over segments
     &, n_send                                                          &
                    ! number of items of data to send
     &, n_recv                                                          &
                    ! number of items of data to receive
     &, level                                                           &
                    ! loop counter over levels
     &, info                                                            &
                    ! return code
     &, flag                                                            &
                    ! GCOM input code
     &, ErrorStatus ! Return Code

      Integer, Parameter :: theta_halo_type = 1
      Integer :: i,k
      Integer :: off_x, off_y
      Logical :: l_include_halos

      Real, allocatable :: ex_p_global (:)
      Real, allocatable :: p_global (:)
      Real, allocatable :: p_lbc_pts (:,:)

      Character (Len=*), Parameter :: RoutineName = 'Calc_P_Src'
      Character (Len=80)           :: CMessage

! ----------------------------------------
! Find maximum number of levels on each PE
! ----------------------------------------

      map(:) = 0
      level_count(:) = 0
      Do k = 1, model_levels
        map(k) = mod ( k-1, nproc)
        level_count(map(k)) = level_count(map(k)) + 1
      End Do

      max_level_count = max( 1, level_count(mype) )

! ----------------------------
! Get global field for exner_p
! ----------------------------

!     Allocate space for :
!     global field of exner p
!     global field of p
!     p at LBC points

      allocate (ex_p_global(glsize(1,fld_type_p)*glsize(2,fld_type_p)))
      allocate (p_global (glsize(1,fld_type_p)*glsize(2,fld_type_p)))
      allocate (p_lbc_pts(len_intf_p,max_level_count))

      level_count(:) = 0

      Do level = 1, model_levels

! DEPENDS ON: gather_field
        Call gather_field ( ex_p(1,level),                              &
     &                      ex_p_global,                                &
     &                      lasize(1,fld_type_p,theta_halo_type),       &
     &                      lasize(2,fld_type_p,theta_halo_type),       &
     &                      glsize(1,fld_type_p),                       &
     &                      glsize(2,fld_type_p),                       &
     &                      fld_type_p,                                 &
     &                      theta_halo_type,                            &
     &                      pe_for_levels(level),                       &
     &                      gc_all_proc_group,                          &
     &                      errorstatus,                                &
     &                      cmessage )

        If ( mype == pe_for_levels(level) ) Then

          level_count(mype) = level_count(mype) + 1

!         --------------------
!         Convert exner_p to p
!         --------------------

          off_x = 0
          off_y = 0
          l_include_halos = .false.

! DEPENDS ON: calc_p_from_exner
          Call Calc_P_from_Exner (                                      &
     &         p_global, kappa, p_zero,                                 &
     &         glsize(1,fld_type_p), glsize(2,fld_type_p), 1,           &
     &         off_x, off_y,                                            &
     &         ex_p_global, L_include_halos)

! P is in p_global

!         ---------------------------
!         Interpolate p to lbc points
!         ---------------------------

! DEPENDS ON: h_int_bl
          Call H_INT_BL(glsize(2,fld_type_p),                           &
     &                  glsize(1,fld_type_p),LEN_INTF_P                 &
     &,                 AP_INDEX_B_L(IP_P),AP_INDEX_B_R(IP_P)           &
     &,                 p_global                                        &
     &,                 AP_WEIGHT_B_L(IP_P),AP_WEIGHT_B_R(IP_P)         &
     &,                 AP_WEIGHT_T_L(IP_P),AP_WEIGHT_T_R(IP_P)         &
     &,                 p_lbc_pts(1,level_count(mype)) )

! Interpolated P is in p_lbc_pts

        End If
      End Do

! Ensure that P has been set up for all levels

      Call GC_SSync (nproc, info)

! --------------------------------------
! Gather all levels of P on PE GATHER_PE
! --------------------------------------

!     Set up SEND and RECEIVE maps for GC_RALLTOALL

      n_send=0
      n_recv=0
      map (:) = 0
      send_map (:,:) = 0
      recv_map (:,:) = 0
      level_count(:) = 0
      level_on_gat_pe(:) = 0

      Do k=1,model_levels
        map(k)              = MOD(k-1,nproc)
        level_count(map(k)) = level_count(map(k))+1
        level_on_gat_pe(k)  = level_count(map(k))
      Enddo

      Do k=1,model_levels

        If (Map(k) == mype) Then ! I need to send data
          n_send=n_send+1
          send_map(S_DESTINATION_PE,n_send)=GATHER_PE
          send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,n_send)=                &
     &    1+((LEVEL_ON_GAT_PE(k)-1)*len_intf_p)
          send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,n_send)=1
          send_map(S_STRIDE_IN_SEND_ARRAY,n_send)=1
          send_map(S_ELEMENT_LENGTH,n_send)=len_intf_p
          send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,n_send)=                &
     &      1+(k-1)*len_intf_p
          send_map(S_STRIDE_IN_RECV_ARRAY,n_send)=1
        Endif

        If (mype == GATHER_PE) Then ! I need to receive data
          n_recv=n_recv+1
          recv_map(R_SOURCE_PE,n_recv)=MAP(k)
          recv_map(R_BASE_ADDRESS_IN_RECV_ARRAY,n_recv)=                &
     &      1+(k-1)*len_intf_p
          recv_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,n_recv)=1
          recv_map(R_STRIDE_IN_RECV_ARRAY,n_recv)=1
          recv_map(R_ELEMENT_LENGTH,n_recv)=len_intf_p
          recv_map(R_BASE_ADDRESS_IN_SEND_ARRAY,n_recv)=                &
     &      1+((LEVEL_ON_GAT_PE(k)-1)*len_intf_p)
          recv_map(R_STRIDE_IN_SEND_ARRAY,n_recv)=1
        Endif

      Enddo ! k

      flag=GC_NONE
      info=GC_NONE

      Call GCG_RALLTOALLE(p_lbc_pts, send_map, n_send,                  &
     &                    len_intf_p,                                   &
     &                    p_src, recv_map, n_recv,                      &
     &                    len_intf_p*model_levels,                      &
     &                    gc_all_proc_group,flag,info)


      If (info /= GC_NONE) THEN
        ErrorStatus = 30
        Write (cmessage,*) 'GCG_RALLTOALLE failed'
! DEPENDS ON: ereport
        Call Ereport ( RoutineName, ErrorStatus, cmessage)
      Endif

      deallocate (ex_p_global)
      deallocate (p_global)
      deallocate (p_lbc_pts)

      Return
      END SUBROUTINE Calc_P_Src
#endif
