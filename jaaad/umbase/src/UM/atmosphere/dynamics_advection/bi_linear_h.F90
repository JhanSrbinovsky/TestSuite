#if defined(A12_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Bi_Linear_H.

      Subroutine Bi_Linear_H(                                           &
     &                       Data_in, lambda_out, phi_out,              &
     &                       dim_i_in, dim_j_in, dim_k_in,              &
     &                       dim_i_out, dim_j_out, dim_k_out,           &
     &                       row_length, rows,                          &
     &                       i_out_in, j_out_in,                        &
     &                       weight_lambda, weight_phi,                 &
     &                       model_domain,                              &
     &                       me, n_procx,                               &
     &                       halo_i, halo_j, datastart,                 &
     &                       g_row_length, g_i_pe, at_extremity,        &
     &                       pole_handling, proc_row_group,             &
     &                       L_sl_halo_reprod, L_regular,               &

     &                       Data_out)

! Purpose:
!          Performs bi-linear horizontal interpolation of a field
!          defined on a w grid to another grid. Input data can be on a
!          sphere or be a rectangular box. Requested output points must
!          lie inside the region defined for the input Data.
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
! Original Programmer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
! Version     Date     Comment
! -------     ----     -------
!   5.1     22/03/00   Fix to stop out of memory address calls,
!                                                 Andy Malcolm
!LL   5.1   11/02/00  Use DOMTYP parameters                    P.Burton
!   5.3     15/09/01  add mt_bi_cyclic_LAM code            A. Malcolm
!  5.3     14/12/01  Make GCOM tags fit within MPI limits. P.Selwood.
!  6.2     21/10/05  Remove commented out code. P.Selwood
!  6.2     25/12/05  Variable resolution changes          Yongming Tang
!  6.4     27/11/06  Fix to integer arithmetic bug          A.Malcolm
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  dim_i_in                                                        &
                    ! Dimension of Data_in in i direction.
     &, dim_j_in                                                        &
                    ! Dimension of Data_in in j direction.
     &, dim_k_in                                                        &
                    ! Dimension of Data_in in k direction.
     &, dim_i_out                                                       &
                    ! Dimension of Data_out in i direction.
     &, dim_j_out                                                       &
                    ! Dimension of Data_out in j direction.
     &, dim_k_out                                                       &
                    ! Dimension of Data_out in j direction.
     &, me                                                              &
                    ! My processor number
     &, n_procx                                                         &
                    ! Number of processors in longitude
     &, halo_i                                                          &
                    ! Size of halo in i direction.
     &, halo_j                                                          &
                    ! Size of halo in j direction.
     &, datastart(3)                                                    &
                     ! First gridpoints held by this processor.
     &, g_row_length                                                    &
                     ! global number of points on a row
     &, g_i_pe(1-halo_i:g_row_length+halo_i)                            &
                                             ! processor on my
                   ! processor-row holding a given value in i direction
     &, pole_handling                                                   &
                       ! How to treat the poles:
                       !   0 - no calculations at the poles
                       !   1 - poles in one point
                       !   2 - poles in all points
     &, proc_row_group                                                  &
                   ! Group id for processors on the same row
     &, row_length                                                      &
!                  ! points in a row on this pe
     &, rows
!                  !  rows on this pe

      Logical                                                           &
     &  L_sl_halo_reprod                                                &
                         ! if true then sl code bit repoducible with
                         ! any sensible halo size
     &, L_regular  ! false if variable resolution

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Real                                                              &
     &  Data_in (1-halo_i:dim_i_in+halo_i, 1-halo_j:dim_j_in+halo_j,    &
     &           dim_k_in)                 ! data to be interpolated

      Real                                                              &
     &  lambda_out (dim_i_out, dim_j_out,dim_k_in)                      &
                                                    ! Lambda
                                           ! co-ordinate of
                                           ! output data on
                                           ! input.
     &, phi_out (dim_i_out, dim_j_out, dim_k_in)     ! Phi Co-ordinate
                                           ! of output data
                                           ! on input.
      Real                                                              &
     &  weight_lambda (dim_i_out, dim_j_out, dim_k_out)                 &
     &, weight_phi (dim_i_out, dim_j_out, dim_k_out)

      Integer                                                           &
     &  i_out_in (dim_i_out, dim_j_out, dim_k_out)                      &
     &, j_out_in (dim_i_out, dim_j_out, dim_k_out)

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid


! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  Data_out (dim_i_out, dim_j_out, dim_k_out) ! data interpolated
                                                   ! to desired locatns.

! Local Variables.

! scalars

      Integer                                                           &
     &  i, j, k, kk                                                     &
                             ! Loop indices
     &, j0, j1                                                          &
     &, index

      Integer :: ErrorStatus

! arrays

      Integer                                                           &
     &  i_out (dim_i_out, dim_j_out, dim_k_out)                         &
     &, j_out (dim_i_out, dim_j_out, dim_k_out)

! Varibles applied in the "compute-on-demand" strategy

      Integer                                                           &
     &  ime, ibase, irecv, my_imin, my_imax, dim_e_out                  &
     &, nsend, nrecv, info, len, itmp, sender                           &
     &, my_iminp, my_imaxp

      Integer                                                           &
     &  n_sendto(0:n_procx-1), n_recvfrom(0:n_procx-1)                  &
     &, i_store(dim_k_out*dim_i_out*8,0:n_procx-1)                      &
     &, j_store(dim_k_out*dim_i_out*8,0:n_procx-1)                      &
     &, k_store(dim_k_out*dim_i_out*8,0:n_procx-1)                      &
     &, i_out_e(dim_k_out*dim_i_out*8)                                  &
     &, j_out_e(dim_k_out*dim_i_out*8)                                  &
     &, isend_arr(2,dim_k_out*dim_i_out*8,0:n_procx-1)                  &
     &, irecv_arr(2,dim_k_out*dim_i_out*8,0:n_procx-1)                  &
     &, sp_send(0:n_procx-1), sp_levels(0:n_procx-1,dim_k_out)          &
     &, np_send(0:n_procx-1), np_levels(0:n_procx-1,dim_k_out)
      Real                                                              &
     &  rsend_arr(2,dim_k_out*dim_i_out*8,0:n_procx-1)                  &
     &, rrecv_arr(2,dim_k_out*dim_i_out*8,0:n_procx-1)                  &
     &, weight_lambda_e(dim_k_out*dim_i_out*8)                          &
     &, weight_phi_e(dim_k_out*dim_i_out*8)                             &
     &, Data_out_e(dim_k_out*dim_i_out*8)                               &
     &, recv_data(dim_k_out*dim_i_out*8,0:n_procx-1)                    &
     &, bcast_data(4*dim_k_out)

#include "parparm.h"
#include "domtyp.h"
! External Routines:
      External bi_linear

! Functions: None

! ----------------------------------------------------------------------
! Section 1.   Extend input data and r arrays to bigger area to
!              allow interpolation to be done without having to re-do
!              any end points.
! ----------------------------------------------------------------------

! No extension required in the parallel version - taken care of by
! swap_bounds.

! ----------------------------------------------------------------------
! Section 2.   For each output point find i,j so that the point on the
!              output grid lies between i and i+1, j and j+1
! ----------------------------------------------------------------------

! i_out and j_out must be copied because they change in this subroutine
      Do k = 1, dim_k_out
        Do j = 1, dim_j_out
          Do i = 1, dim_i_out
            i_out(i,j,k) = i_out_in(i,j,k)
            j_out(i,j,k) = j_out_in(i,j,k)
          End Do
        End Do
      End Do

      If ( n_procx > 1 ) then
      If ( model_domain == mt_Global ) then
! Send the points outside my region to the appropriate processor for
! interpolation. Only performed if the domain is decomposed in the
! i direction and not performed for LAM versions of the model.

! The first and last point I can interpolate in, based on available
! data on this processor

        my_imin = datastart(1) - halo_i + 2
        my_imax = datastart(1) + dim_i_out - 1 + halo_i - 2

! values for use in polar row to ensure pole is only calculated on one
! processor.
        my_iminp = datastart(1)
        my_imaxp = datastart(1)+dim_i_out-1
        If(at_extremity(PWest) ) my_iminp = my_imin
        If(at_extremity(PEast) ) my_imaxp = my_imax

! The base processor on this row, and my address relative to that
! processor

        ibase = (me/n_procx) * n_procx
        ime = me - ibase

        Do i = 0, n_procx-1
          n_sendto(i) = 0
        End Do

        j0 = 1
        j1 = dim_j_out

! If the pole values not are going to be used, we we can save a
! large amount of communication by evaluating local dummys at the poles

        If (pole_handling  ==  0) then

          If (at_extremity(PSouth)) then
            j0 = 2
            Do k = 1, dim_k_out
              Do i = 1, dim_i_out
                i_out(i,1,k) = i
              End Do
            End Do
          End If

          If (at_extremity(PNorth)) then
            j1 = dim_j_out-1
            Do k = 1, dim_k_out
              Do i = 1, dim_i_out
                i_out(i,dim_j_out,k) = i
              End Do
            End Do
          End If

        End If

! If all values along a polar row is evaluated in one point, we can
! save a large amount of communication by evaluating only one point
! correctly and then broadcast this value.

        If (pole_handling  ==  1) then

          Do i = 0, n_procx-1
            sp_send(i) = 0
          End Do
          If (at_extremity(PSouth)) then
            j0 = 2
            Do k = 1, dim_k_out
              If (i_out(1,1,k)  >=  my_iminp .and.                      &
     &              i_out(1,1,k)  <=  my_imaxp) then
                i_out(1,1,k) = i_out(1,1,k) - datastart(1) + 1
                sp_send(ime) = sp_send(ime) + 1
                sp_levels(ime,sp_send(ime)) = k
                Do i = 2, dim_i_out
                  i_out(i,1,k) = i ! i_out(i,1,k) - datastart(1) + 1
                End Do
              Else
                sender = g_i_pe(i_out(1,1,k))
                sp_send(sender) = sp_send(sender) + 1
                sp_levels(sender,sp_send(sender)) = k
                Do i = 1, dim_i_out
                  i_out(i,1,k) = i
                End Do
              End If
            End Do
          End If

          Do i = 0, n_procx-1
            np_send(i) = 0
          End Do
          If (at_extremity(PNorth)) then
            j1 = dim_j_out-1
            Do k = 1, dim_k_out
              If (i_out(1,dim_j_out,k)  >=  my_iminp .and.              &
     &            i_out(1,dim_j_out,k)  <=  my_imaxp)  then
                np_send(ime) = np_send(ime) + 1
                np_levels(ime,np_send(ime)) = k
                i_out(1,dim_j_out,k) =                                  &
     &                 i_out(1,dim_j_out,k) - datastart(1) + 1
                Do i = 2, dim_i_out
                  i_out(i,dim_j_out,k) = i
                End Do
              Else
                sender = g_i_pe(i_out(1,dim_j_out,k))
                np_send(sender) = np_send(sender) + 1
                np_levels(sender,np_send(sender)) = k
                Do i = 1, dim_i_out
                  i_out(i,dim_j_out,k) = i
                End Do
              End If
            End Do
          End If

        End If

        If( L_sl_halo_reprod) Then

! On the global boundaries, use i_out < 1 or i_out > g_row_length
! if that makes local computation possible. Not required when
! L_sl_halo_reprod is false is other logic ensures this is done.

! This code unsafe if applied at poles, where it isn't required.

          If (at_extremity(PWest)) then
            Do k = 1, dim_k_out
               Do j = j0, j1
                  Do i = 1, halo_i
                     If (i_out(i,j,k)  >   g_row_length-halo_i+2)       &
     &                    i_out(i,j,k) = i_out(i,j,k) - g_row_length
                  End Do
               End Do
            End Do
          End If
          If (at_extremity(PEast)) then
            Do k = 1, dim_k_out
               Do j = j0, j1
                  Do i = dim_i_out-halo_i+1, dim_i_out
                     If (i_out(i,j,k)  <   halo_i-1)                    &
     &                    i_out(i,j,k) = i_out(i,j,k) + g_row_length
                  End Do
               End Do
            End Do
          End If

        End If ! on L_sl_halo_reprod

        Do k = 1, dim_k_out
          Do j = j0, j1
            Do i = 1, dim_i_out
              If (i_out(i,j,k)  >=  my_imin .and.                       &
     &            i_out(i,j,k)  <=  my_imax) then
! Process locally, so find the local destination
                i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
              Else
! Send to a remote processor, given by the array g_i_pe
                irecv = g_i_pe(i_out(i,j,k))
             if(irecv > 1000) then
                print*,'stop here', irecv
             endif
                n_sendto(irecv) = n_sendto(irecv) + 1
                itmp = n_sendto(irecv)
                isend_arr(1,itmp,irecv) = i_out(i,j,k)
                isend_arr(2,itmp,irecv) = j_out(i,j,k)
                rsend_arr(1,itmp,irecv) = weight_lambda(i,j,k)
                rsend_arr(2,itmp,irecv) = weight_phi(i,j,k)
                i_store(itmp,irecv) = i
                j_store(itmp,irecv) = j
                k_store(itmp,irecv) = k
!!! fix to stop out of memory address calls
                i_out(i,j,k)=i
                j_out(i,j,k)=j
              End If
            End Do
          End Do
        End Do

        nsend = 0
        Do i = 0,n_procx-1
          call gc_isend(10*(me+1)+ibase+i, 1, ibase+i, info,            &
     &                  n_recvfrom(ime), n_sendto(i))
          nsend = nsend + n_sendto(i)
        End Do

        Call gcg_ssync(proc_row_group, info)

        nrecv = 0
        Do i = 0, n_procx-1
           call gc_irecv(10*(ibase+i+1)+me, 1, ibase+i, info,           &
     &           n_recvfrom(i), n_sendto(ime))
           nrecv = nrecv + n_recvfrom(i)
        End Do

!        Call gcg_ssync(proc_row_group, info)

        Do i = 0,n_procx-1
          len = 2*n_sendto(i)
          If (n_sendto(i)  >   0) then
             call gc_rsend(20*(me+1)+ibase+i, len, ibase+i, info,       &
     &                     rrecv_arr(1,1,ime), rsend_arr(1,1,i))
          End If
        End Do

        Call gcg_ssync(proc_row_group, info)

        Do i = 0,n_procx-1
          len = 2*n_recvfrom(i)
          If (n_recvfrom(i)  >   0) then
            call gc_rrecv(20*(ibase+i+1)+me, len, ibase+i, info,        &
     &              rrecv_arr(1,1,i), rsend_arr(1,1,ime))
          End If
        End Do

!        Call gcg_ssync(proc_row_group, info)

        Do i = 0,n_procx-1
          len = 2*n_sendto(i)
          If (n_sendto(i)  >   0) then
            call gc_isend(30*(me+1)+ibase+i, len, ibase+i, info,        &
     &          irecv_arr(1,1,ime), isend_arr(1,1,i))
          End If
        End Do

        Call gcg_ssync(proc_row_group, info)

        Do i = 0,n_procx-1
          len = 2*n_recvfrom(i)
          If (n_recvfrom(i)  >   0) then
            call gc_irecv(30*(ibase+i+1)+me, len, ibase+i, info,        &
     &              irecv_arr(1,1,i), isend_arr(1,1,ime))
          End If
        End Do

!      Else If ( L_regular ) then   ! for LAMs only
      Else      ! for LAMs only

        k = 1
        j = 1
        Do i = 1, dim_i_out * dim_j_out * dim_k_out
          i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
        End Do  ! i = 1, dim_i_out * dim_j_out * dim_k_out

      EndIf ! model_domain == mt_Global

      EndIf ! n_procx > 1

! ----------------------------------------------------------------------
! Section 3.   Perform required Interpolation.
! ----------------------------------------------------------------------

! DEPENDS ON: bi_linear
      call bi_linear (dim_i_out, dim_j_out, dim_k_out,                  &
     &                dim_i_in, dim_j_in, dim_k_in,                     &
     &                halo_i, halo_j, Data_in,                          &
     &                i_out, j_out, weight_lambda, weight_phi,          &
     &                Data_out)

      If ( n_procx > 1 .and. model_domain == mt_Global ) then

        dim_e_out = 0
        Do i = 0, n_procx-1
          If (n_recvfrom(i)  >   0) then
            Do j = 1, n_recvfrom(i)
              dim_e_out = dim_e_out + 1
              i_out_e(dim_e_out) = irecv_arr(1,j,i)-datastart(1)+1
              j_out_e(dim_e_out) = irecv_arr(2,j,i)
              weight_lambda_e(dim_e_out) = rrecv_arr(1,j,i)
              weight_phi_e(dim_e_out) = rrecv_arr(2,j,i)
            End Do
          End If
        End Do

        If (dim_e_out  >   dim_k_out*dim_i_out*8) Then
          ErrorStatus = 10
! DEPENDS ON: Ereport
          Call Ereport("Bi_linear_h", ErrorStatus,                      &
     &         "over-writing due to dim_e_out size" )
        End If 

        If (dim_e_out  >   0) then
! DEPENDS ON: bi_linear
          call bi_linear (dim_e_out, 1, 1,                              &
     &                    dim_i_in, dim_j_in, dim_k_in,                 &
     &                    halo_i, halo_j, Data_in,                      &
     &                    i_out_e, j_out_e,                             &
     &                    weight_lambda_e, weight_phi_e,                &
     &                    Data_out_e)
        End If


        nsend = 1
        Do i = 0, n_procx-1
          If (n_recvfrom(i)  >   0) then
            len = n_recvfrom(i)
            call gc_rsend(40*(me+1)+ibase+i, len, ibase+i, info,        &
     &             recv_data(1,ime), Data_out_e(nsend))
            nsend = nsend + n_recvfrom(i)
          End If
        End Do

        Call gcg_ssync(proc_row_group, info)

        Do i = 0, n_procx-1
          If (n_sendto(i)  >   0) then
            len = n_sendto(i)
            call gc_rrecv(40*(ibase+i+1)+me, len, ibase+i, info,        &
     &             recv_data(1,i), Data_out_e)
            Do j = 1, n_sendto(i)
              Data_out(i_store(j,i), j_store(j,i),                      &
     &                 k_store(j,i)) = recv_data(j,i)
            End Do
          End If
        End Do

! Now distribute the pole values.

        If (pole_handling  ==  1) then

          If (at_extremity(PSouth)) then
            Do j = 0, n_procx-1
              If (sp_send(j)  >   0) then
                Do kk = 1, sp_send(j)
                  k = sp_levels(j,kk)
                  bcast_data(sp_send(j)+kk) =                           &
     &                     Data_out(1,1,k)
                End Do
                len = sp_send(j)
                call gcg_rbcast(201, len, ibase+j,                      &
     &                 proc_row_group, info, bcast_data)
                Do kk = 1, sp_send(j)
                  k = sp_levels(j,kk)
                  Do i = 1, dim_i_out
                    Data_out(i,1,k) =                                   &
     &                       bcast_data(sp_send(j)+kk)
                  End Do
                End Do
              End If
            End Do
          End If

          If (at_extremity(PNorth)) then
            Do j = 0, n_procx-1
              If (np_send(j)  >   0) then
                Do kk = 1, np_send(j)
                  k = np_levels(j,kk)
                  bcast_data(np_send(j)+kk) =                           &
     &                     Data_out(1,dim_j_out,k)
                End Do
                len = np_send(j)
                call gcg_rbcast(201, len, ibase+j,                      &
     &                 proc_row_group, info, bcast_data)
                Do kk = 1, np_send(j)
                  k = np_levels(j,kk)
                  Do i = 1, dim_i_out
                    Data_out(i,dim_j_out,k) =                           &
     &                       bcast_data(np_send(j)+kk)
                  End Do
                End Do
              End If
            End Do
          End If

        End If

      End If !  n_procx > 1 .and. model_domain == mt_Global

! End of routine.
      return
      END SUBROUTINE Bi_Linear_H

#endif
