
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine GCR_abs_print

      subroutine GCR_abs_print(                                         &
     &                           Error, off_x, off_y,                   &
     &                            gc_proc_row_group,                    &
     &                            global_row_length, global_rows,       &
     &                            g_rows, g_row_length, g_datastart,    &
     &                            n_proc, n_procx, n_procy, me,         &
     &                            row_length, rows,                     &
     &                            model_levels, model_domain,           &
     &                            at_extremity)

! Purpose:
!          Diagnostic print of absolute error norms
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
! Original Progammer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
! Date     Version     Comment
! ----     -------     -------
! 22/03/00 5.1        Changed j_start/end for cyclic LAM  Andy Malcolm
!LL   5.1   10/02/00  Use DOMTYP parameters                    P.Burton
!     5.3   19/10/01  Use appropriate gcg routines.   S. Cusack
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
      Implicit None

! Arguments with Intent IN. ie: Input variables.
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, global_row_length                                               &
     &, global_rows                                                     &
     &, gc_proc_row_group                                               &
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, me                                                              &
     &, off_x                                                           &
     &, off_y                                                           &
     &, model_domain

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid


      Integer                                                           &
     &  g_rows(0:n_proc-1)                                              &
     &, g_row_length(0:n_proc-1)                                        &
     &, g_datastart(3,0:n_proc-1)

      Real                                                              &
     &  Error (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)

! Local variables
      Real                                                              &
     &  error_sum(row_length, rows)                                     &
     &, mean(global_rows)                                               &
     &, rbuf(global_rows)                                               &
     &, max_val(model_levels)

      Integer d, i, j, k, istat, info, gj                               &
     &,       i_start, i_end, j_start, j_end

!========================== COMDECK PARPARM ====================
!   Description:
!
!   This COMDECK contains PARAMETERs for the mpp-UM
!
!   Two sets of parameters are set up -
!     i)  for the mpp-UM itself.
!     ii) for the interface to the Message Passing Software.
!
      !=================================================================
      ! Parameters needed for the mpp-UM
      !=================================================================
      ! maximum number of spatial dimensions
      INTEGER,PARAMETER:: Ndim_max = 3 ! 3d data

      ! number of different halo types
      INTEGER,PARAMETER:: NHalo_max = 3 ! for N.D. atmos. model

      INTEGER,PARAMETER:: halo_type_single   = 1
      INTEGER,PARAMETER:: halo_type_extended = 2
      INTEGER,PARAMETER:: halo_type_no_halo  = 3

! FLDTYPE definitions for the different field types recognised on the
! decomposition
      INTEGER,PARAMETER:: Nfld_max=7 ! maximum number of field types
      INTEGER,PARAMETER:: fld_type_p=1       ! grid on P points
      INTEGER,PARAMETER:: fld_type_u=2       ! grid on U points
      INTEGER,PARAMETER:: fld_type_v=3       ! grid on V points
      INTEGER,PARAMETER:: fld_type_comp_wave  = 4
                              ! Compressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_full_wave  = 5
                              ! Uncompressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_rim_wave   = 6
                              ! Boundary data for WAM Wave Field
      INTEGER,PARAMETER:: fld_type_r=7       ! grid on river points
      INTEGER,PARAMETER:: fld_type_unknown=-1! non-standard grid
! FLDTYPE end

      ! Used in addressing to indicate if calculation is for a local or
      ! global (ie. disk dump) size

      INTEGER,PARAMETER:: local_data=1
      INTEGER,PARAMETER:: global_dump_data=2

      ! maximum permitted size of a halo
      INTEGER,PARAMETER:: Max_Halo_Size=10

      !=================================================================
      ! Parameters needed for the Message Passing Software
      !=================================================================
      INTEGER,PARAMETER:: Maxproc = 512 ! Max number of processors

      ! Processor addresses in the neighbour array
      INTEGER,PARAMETER:: PNorth   = 1
      INTEGER,PARAMETER:: PEast    = 2
      INTEGER,PARAMETER:: PSouth   = 3
      INTEGER,PARAMETER:: PWest    = 4

      ! Value in neighbour array if the domain has  no neighbour in this
      ! direction. Otherwise the value will be the tid of the neighbor
      INTEGER,PARAMETER:: NoDomain = -1

      INTEGER,PARAMETER:: BC_STATIC   = 1 ! Static boundary conditions
      INTEGER,PARAMETER:: BC_CYCLIC   = 2 ! Cyclic boundary conditions
      INTEGER,PARAMETER:: BC_OVERPOLE = 3 ! Transfer over pole
! PARPARM end
! DOMTYP contains different model domain types
!
! Author : P.Burton
! History:
! Version  Date      Comment.
! 5.0      15/04/99  New comdeck
! 5.2      15/11/00  add bi_cyclic_lam domain   A. Malcolm

      INTEGER,PARAMETER:: mt_global        = 1
      INTEGER,PARAMETER:: mt_lam           = 2
      INTEGER,PARAMETER:: mt_cyclic_lam    = 3
      INTEGER,PARAMETER:: mt_bi_cyclic_lam = 4
      INTEGER,PARAMETER:: mt_single_column = 5
! DOMTYP end
!    External Routines:
      External                                                          &
     &  gcg_rvecsumr, gc_imax, gc_rsend, gc_rrecv, gc_ssync

! ----------------------------------------------------------------------
! Section 1.   Calculate error norm over rows and levels.
! ----------------------------------------------------------------------

      i_start = 1
      i_end = row_length
      j_start = 1
      j_end = rows
      If (model_domain  ==  mt_lam) Then
        if(at_extremity(PSouth)) j_start = 2
        if(at_extremity(PNorth)) j_end = rows-1
        if(at_extremity(PEast)) i_end = row_length-1
        if(at_extremity(PWest)) i_start = 2
      End If
      If (model_domain  ==  mt_cyclic_LAM) Then
        if(at_extremity(PSouth)) j_start = 2
        if(at_extremity(PNorth)) j_end = rows-1
      End If

      Do k = 1, model_levels
        max_val(k) = 0.0
      End Do
      k = 1
      Do j = 1, rows
        Do i = 1, row_length
          error_sum(i,j) = 0.0
        End Do
      End Do
      Do k = 1, model_levels
        Do j = j_start, j_end
          Do i = i_start, i_end
            error_sum(i,j) = error_sum(i,j) + error(i,j,k)
            If (error(i,j,k)  >   max_val(k) )                          &
     &        max_val(k) = error(i,j,k)
          End Do
        End Do
      End Do

      call gcg_rvecsumr(row_length,row_length,1,                        &
     &   rows,error_sum,gc_proc_row_group,istat,mean)

      If (model_domain  /=  mt_lam) Then
        Do j = 1, rows
          mean(j) = mean(j) / (model_levels*global_row_length)
        End Do
      Else
        Do j = 1, rows
          mean(j) = mean(j) / (model_levels*(global_row_length-2))
        End Do
      End If

! Gather all data to processor 0 from its column and form
! global rows fields

      Call gc_imax(model_levels, n_proc, info, max_val)

! step 1 send data to processor zero
      If (me  /=  0) then
        Do j = 1, rows
          rbuf(j) = mean(j)
        End Do

        call gc_rsend(100*me, g_rows(me),                               &
     &                  0, info, rbuf, rbuf)
      End If

      call gc_ssync(n_proc,info)

! step 2 processor zero receives data and puts in correct location

      If (me  ==  0) then
        Do d = 1, n_procx*n_procy-1
          call gc_rrecv(100*d, g_rows(d),                               &
     &                  d, info, rbuf, rbuf)
          Do j = 1, g_rows(d)
            gj = g_datastart(2,d) + j - 1
            mean(gj) = rbuf(j)
          End Do
        End Do

        print*,' mean approximate absolute density change left '
        Do j = 1, global_rows
          write(*,100) j, mean(j)
        End Do
        print*,' max value left per level '
        Do k = 1, model_levels
            write(*,100) k, max_val(k)
        End Do

      End If

      call gc_ssync(n_proc,info)

 100  Format(1x,I3,1x,E10.2)

      return
      END SUBROUTINE GCR_abs_print

