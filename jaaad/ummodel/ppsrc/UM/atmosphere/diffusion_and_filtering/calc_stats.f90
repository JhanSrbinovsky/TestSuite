
! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT*****************************
!
! Subroutine calc_stats

      subroutine calc_stats(                                            &
     &                      field,                                      &
     &                      row_length, rows, levels,                   &
     &                      halo_x, halo_y, n_proc, size,               &
     &                      mean, variance, std_dev, skew, kurtosis)
      implicit none
!
! description:
! Standard statistics of a field
! see Numerical recipes p604-608
!
! method:
! Do global sums  on processor 0 only.
! Results from different processor configuration will not
! bit-reproduce but model fields are not affected (diagnostic only)
!
! current code owner : Terry Davies
!
! history
!  model    date      modification history from model version 6.2
!  version
!    6.2    25/10/05  new code               Terry Davies
!
! subroutine arguments:

      integer                                                           &
     &  row_length                                                      &
                           ! in columns
     &, rows                                                            &
                           ! in rows
     &, levels                                                          &
                           ! no. of levels
     &, size                                                            & 
                           ! array size for summing over pe's
     &, halo_x                                                          &
                           ! in halo in i direction ew
     &, halo_y                                                          &
                           ! in halo in j direction ns
     &, n_proc             ! Total number of processors

! inout
      real                                                              &
     & field(1-halo_x: row_length+halo_x, 1-halo_y: rows+halo_y         &
     &       , levels)            ! in array containing fields

! inout
      real                                                              &
     &  mean(levels)                                                    &
     &, std_dev(levels)                                                 &
     &, variance(levels)                                                &
     &, skew(levels)                                                    &
     &, kurtosis(levels)
!
! parameters and common

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

! local variables
      real                                                              &
     &  sum(size)                                                       &
                   ! various sums
     &, points                                                          &
                         ! field length (local and global)
     &, recip_points                                                    &
                         ! reciprocal field length
     &, diff                                                            &
                       ! difference
     &, diff_power     ! various powers of difference

      integer                                                           &
     &  i, j, k, k1                                                     &
                       ! loop variables
     &, info           ! return code from gc stuff

!-----------------------------------------------------------------------
! mpp code for summation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! 1. Do the sums
!-----------------------------------------------------------------------

        do k = 1, levels
          sum(k) = 0.0   ! to sum field
          do j = 1, rows
            do i = 1, row_length
              sum(k) = sum(k) + field(i,j,k)
            enddo
          enddo
        enddo  ! k = 1, levels

!  Put local number of points in sum(levels + 1)
        k1 = levels + 1
        sum(k1) = row_length * rows  ! to field size

        call GC_RSUMR(k1, n_proc, info, sum)

!  Total number of points is in sum(levels + 1)
        points = sum(levels+1)
        recip_points = 1.0 / points

        do k = 1, levels
          mean(k) = sum(k) * recip_points
        enddo

!   sum(1) for difference from mean, sum(2) for variance
!   sum(3) for skew, sum(4) for kurtosis
        do i = 1, size
          sum(i) = 0.0   ! 4 sums required
        enddo

        do k = 1, levels
          k1 = 4 * k - 3
          sum(k1) = 0.0
          sum(k1+1) = 0.0
          sum(k1+2) = 0.0
          sum(k1+3) = 0.0
          do j = 1, rows
            do i = 1, row_length
              diff  = field(i,j,k) - mean(k)
              sum(k1) = sum(k1) + diff
              diff_power  = diff * diff
              sum(k1+1) = sum(k1+1) + diff_power  !  variance
              diff_power  = diff_power * diff
              sum(k1+2) = sum(k1+2) + diff_power  !  skew
              diff_power  = diff_power * diff
              sum(k1+3) = sum(k1+3) + diff_power  !  kurtosis
            enddo
          enddo
        enddo

      k1 = 4 * levels
      call GC_RSUMR(k1, n_proc, info, sum)

      do k = 1, levels
        k1 = 4 * k - 3
        variance(k) = (sum(k1+1) - sum(k1) * sum(k1) * recip_points) /  &
     &                                                    (points - 1.0)
        if( variance(k) > 0 ) then
          std_dev(k) = sqrt(variance(k))
          skew(k) = sum(k1+2) * recip_points / std_dev(k)**3
          kurtosis(k) = sum(k1+3) * recip_points /                      &
     &                                   (variance(k) * variance(k)) - 3
        else
          std_dev(k) = 0.0
          skew(k) = 0.0
          kurtosis(k) = 0.0
        endif  !  variance(k) > 0
      enddo !  k = 1, levels

      return      ! End of routine calc_stats
      END SUBROUTINE calc_stats

