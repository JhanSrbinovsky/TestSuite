#if defined(A14_1B) || defined(A11_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ General purpose global sum routine for calculating the
!  sum of a horizontal field
!
! Subroutine Interface
      SUBROUTINE DO_SUMS(SUM_ARRAY,                                     &
     &                   row_length, rows, halo_il,halo_jl,             &
     &                   grid_type,halo_type,                           &
     &                   N_SUMS,SUM_RESULTS)
      IMPLICIT NONE
!
! Description:
! Primarily written for the energy correction suite of routines:
! Calculate N_SUMS global sums of the SUM_ARRAY field between
! the points on the actual grid ignoring any halos adding results
! onto SUM_RESULTS
!
! Method:
!#if defined(REPROD)
! Farm out N_SUMS global sums : 1 global sum/processor - do the
! sums and then return the results to all processors
!#else
! Every processor does its local part of the sum - then all these
! sub-sums are summed together.
!#endif
!
! Current code owner : Paul Burton
!
! History
!  Model    Date      Modification history from model version 4.1
!  version
!    4.1    9/11/95   New DECK created to make EMDIAG suitable for
!                     MPP use. P.Burton
!    4.2    18/11/96  *CALL to AMAXSIZE for MaxFieldSize  P.Burton
!    4.2    18/10/96  New name for group of processors in gather_field
!                     P.Burton
!    4.3    29/04/97  Correct call to GATHER_FIELD. D Robinson.
!    5.1    25/01/00  Version for New dynamics code. R A Stratton
!    5.3    25/10/01  Redimension global_sum_data to a more
!                     appropriate size to prevent memory problems.
!    6.4    14/11/06  Add if defined A11_2A              Andy Malcolm
!
! Subroutine Arguments:

      INTEGER                                                           &
     &  row_length                                                      &
                           ! IN columns
     &, rows                                                            &
                           ! IN rows
     &, halo_il                                                         &
                           ! IN halo in i direction EW
     &, halo_jl                                                         &
                           ! IN halo in j direction NS
     &, grid_type                                                       &
                           ! IN type of grid
     &, halo_type                                                       &
                           ! IN type of halo i=on input grid
     &, N_SUMS             ! IN number of sums to do

      REAL                                                              &
     & SUM_ARRAY((row_length+2*halo_il)*(rows+2*halo_jl),N_SUMS),       &
!                                 ! IN array containing fields
!                                 !    to be summed
     & SUM_RESULTS(N_SUMS)        ! INOUT sum of SUM_ARRAY added onto
!                                 !       initial value of SUM_RESULTS

! Parameters and COMMON

#include "parvars.h"
#include "gccom.h"

! Local variabels
      REAL sum_results_tmp(N_SUMS) ! actual sum which will eventually
!                                    ! be added to SUM_RESULTS

      INTEGER info  ! return code from GC stuff

#if defined(REPROD)
      INTEGER MAP(N_SUMS),                                              &
                            ! processor number for sum
     &        N_SUMS_ON_PROC(0:MAXPROC),                                &
                                          ! number of sums to do on pe
     &        RESULT_NUMBER,                                            &
                             ! result index in SUM_RESULTS
     &        ICODE          ! return code from gather

      REAL global_sum_data(glsize(1,grid_type)*glsize(2,grid_type),     &
     &                     (N_SUMS/nproc)+1)
                            ! area for doing global sums in

      CHARACTER*(80)                                                    &
     &  CMESSAGE        ! Error message
#endif

      INTEGER I,J,K,ij  ! loop variables

!------------------------------------------------------------------------
! Non MPP code for summation - removed
!------------------------------------------------------------------------
! MPP code for summation
!------------------------------------------------------------------------
#if defined(REPROD)
! Reproducible parallel global sums
! We assume all the fields are standard P_FIELDS mapping
! onto the full global grid

! 3. Calculate mapping - which sum is done on which processor

      DO K=1,N_SUMS
        SUM_RESULTS_TMP(K)=0.0
        MAP(K)=first_comp_pe+MOD((K-1),nproc)
        IF (mype  ==  MAP(K)) SUM_RESULTS_TMP(K)=SUM_RESULTS(K)
      ENDDO

! 4. Distribute the sums

      DO K=0,nproc-1
        N_SUMS_ON_PROC(K)=0
      ENDDO

      DO K=1,N_SUMS
        N_SUMS_ON_PROC(MAP(K))=N_SUMS_ON_PROC(MAP(K))+1

! DEPENDS ON: gather_field
        CALL GATHER_FIELD(SUM_ARRAY(1,K),                               &
     &                    global_sum_data(1,N_SUMS_ON_PROC(MAP(K))),    &
     &                    lasize(1,grid_type,halo_type),                &
     &                    lasize(2,grid_type,halo_type),                &
     &                    glsize(1,grid_type),glsize(2,grid_type),      &
     &                    grid_type,halo_type,                          &
     &                    MAP(K),GC_ALL_PROC_GROUP,                     &
     &                    icode,cmessage)

      ENDDO ! K : loop over N_SUMS

! 5. And do the sums

      DO K=1,N_SUMS_ON_PROC(mype)
        RESULT_NUMBER=(K-1)*nproc+mype+1

        Do i=1,glsize(2,grid_type)*glsize(1,grid_type) ! global rows
!          Do i=1,glsize(1,grid_type)    ! global row_length
           SUM_RESULTS_TMP(RESULT_NUMBER)=                              &
     &      SUM_RESULTS_TMP(RESULT_NUMBER)+ global_sum_data(I,K)
!          Enddo  ! i loop over columns
        Enddo  ! j loop over rows

      ENDDO ! K : loop over number of sums I must do

! 6.  Broadcast the results to everyone
!     Rather than do a bcast for each sum, we'll do a
!     parallel sum. Only the processor doing a particular
!     sum will contribute, ie.:
!     SUM  PE0  PE1  PE2  PE3           PE0  PE1  PE2  PE3
!     1    3.2  0.0  0.0  0.0  --SUM--> 3.2  3.2  3.2  3.2
!     2    0.0  9.2  0.0  0.0  --SUM--> 9.2  9.2  9.2  9.2
!     3    0.0  0.0  5.7  0.0  --SUM--> 5.7  5.7  5.7  5.7

       CALL GC_RSUM(N_SUMS,nproc,info,SUM_RESULTS_TMP)

       DO K=1,N_SUMS
         SUM_RESULTS(K)= SUM_RESULTS_TMP(K)
       ENDDO
#else
! This is the faster version of the global sum.
! Each processor works out its local sum across its part of the field
! and then these are summed up.
! This can give non-reproducible answers of two kinds:
! 1) If sums across processors are always done in the same order then
!    the same answer will always be obtained for the same processor
!    arrangement. However, the answer will be different on different
!    processor arrangements
! 2) If the sums across processor are not done in any particular order
!    then different results will be obtained even when the same
!    processor arrangement is used.
! Which case is true depends on how GC_RSUM has been implemented on
! this particular platform. It is probably safest to assume that (2)
! is the case, and this is the method that generally gives the
! fastest sum.

! 2 Do the local sum of my part of SUM_ARRAY

      Do k=1,n_sums
        sum_results_tmp(k)=0.0
        Do j=1,rows
          Do i=1,row_length
            ij = i + halo_il + (j-1)*(row_length + 2*halo_il)
            sum_results_tmp(k)=sum_results_tmp(k)+sum_array(ij,k)
          Enddo
        Enddo
      Enddo  ! K : loop over N_SUMS


! 3.  and now sum up all the local sums, and give everyone
!     the total sum


      CALL GC_RSUM(N_SUMS,nproc,info,sum_results_tmp)

      Do k=1,n_sums
        sum_results(k)=sum_results(k)+sum_results_tmp(k)
      Enddo

#endif

      RETURN
      END SUBROUTINE DO_SUMS
#endif
