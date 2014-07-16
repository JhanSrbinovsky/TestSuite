#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
! Fast-j: Copyright and license yet to be agreed
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Fast-j routine for calculating online photolysis rates
!
!     This subroutine computes a send/receive matrix on all PEs
!
!     Given a number of points to compute which vary on existing PEs,
!     a matrix is computed which contains the number of points to send
!     (negative number) or receive (positiv number) to other PEs
!     The matrix is computed on all PEs.

!     First dimension of matrix is local, second dimesion is target PE,s
!     the new number of points on existing PE is returned. (np_new)
!     a PE is either sender PE, ie. the original number of points is
!     above average or receiver PE.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!  Current Code Owner:       Colin Johsnon/Fiona O'Connor
!                            Oliver Wild
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
      SUBROUTINE FASTJ_CALC_NEW_NUMBER_OF_POINTS (                   &
     &        np,npes,npoints,sr_matrix,me,np_new,sender_PE)

      USE MPL, ONLY :     &
          MPL_COMM_WORLD, &
          MPL_SUM,        &
          MPL_INTEGER

      IMPLICIT   NONE

      INTEGER,INTENT(IN)                       :: np
      INTEGER,INTENT(IN)                       :: npes
      INTEGER,INTENT(OUT),DIMENSION(0:npes-1)  :: npoints
      INTEGER,INTENT(OUT),                                           &
     &    DIMENSION(0:npes-1,0:npes-1)         :: sr_matrix
      INTEGER,INTENT(OUT)                      :: me
      INTEGER,INTENT(OUT)                      :: np_new 

      LOGICAL,INTENT(OUT)                      :: sender_PE

!     Local variables

      INTEGER                                  :: istat
      INTEGER,DIMENSION(0:npes-1)              :: np_in

!     Get MPI Rank, might be different from GC rank

      CALL MPL_COMM_RANK( MPL_COMM_WORLD, me, istat )
      
!     Get old Number of Points from all PEs

      np_in     = 0
      np_in(me) = np

      CALL  MPL_ALLREDUCE (np_in(0),npoints(0),npes,MPL_INTEGER,     &
     &      MPL_SUM, MPL_COMM_WORLD, istat)

      CALL  FASTJ_LOADBAL_CALC_DIST                                  &
     &                      (npes,me,npoints,sr_matrix,sender_PE)

      np_new = npoints(me)

      RETURN
 
      CONTAINS

        SUBROUTINE FASTJ_LOADBAL_CALC_DIST(                          &
     &             n,me,npoints,sendrecv_matrix,sender_PE)
        IMPLICIT NONE

        INTEGER                         :: n
        INTEGER                         :: me
        INTEGER,DIMENSION(0:n-1)        :: npoints
        INTEGER,DIMENSION(0:n-1,0:n-1)  :: sendrecv_matrix

        LOGICAL                         :: sender_PE

!       Local variables
!       Min nunber of points is set to 10. Works for fastj

        INTEGER,PARAMETER               :: MIN_POINTS_2_MOVE=10
        INTEGER                         :: nt
        INTEGER                         :: max_send_PE,max_recv_PE
        INTEGER,DIMENSION(1)            :: xxx
        INTEGER,DIMENSION(0:n-1)        :: Ntrans

        LOGICAL                         :: Not_all_done
        LOGICAL,DIMENSION(0:n-1)        :: s_PE,r_PE

        REAL                            :: average

        s_PE = .false.
        r_PE = .false.

!       Calculate avarage Number of points

        average = sum(float(npoints(0:n-1)))/float(n)

!       Sender_PE or Receiver_PE ?

!CDIR Novector
        WHERE (npoints > average)
          s_PE = .true.
        ELSEWHERE
          r_PE = .true.
        END WHERE

!CDIR Novector
        Ntrans = average-npoints

!       Algrorihtm to compute new Number of Points 

!       Final distribution

        sendrecv_matrix = 0
        Not_all_done    = .true.

        DO WHILE (Not_all_done)

!         Find larges chunk

          xxx         = MINLOC (ntrans(0:n-1))-1
          max_recv_PE = xxx(1)
          xxx         = MAXLOC (ntrans(0:n-1))-1
          max_send_PE = xxx(1)

!         Smaller chunk will be selected and added into sendrecv_matrix
!         This PE then is balanced

          nt = MIN(abs(ntrans(max_recv_PE)),ntrans(max_send_PE))

 
          IF (nt <= MIN_POINTS_2_MOVE)   EXIT             !Done
  
          sendrecv_matrix(max_send_PE,max_recv_PE)                      &
     &              = sendrecv_matrix(max_send_PE,max_recv_PE) + nt 
          sendrecv_matrix(max_recv_PE,max_send_PE)                      &
     &              = sendrecv_matrix(max_recv_PE,max_send_PE) - nt 

          ntrans(max_recv_PE)  = ntrans(max_recv_PE) + nt
          ntrans(max_send_PE)  = ntrans(max_send_PE) - nt

          npoints(max_recv_PE) = npoints(max_recv_PE) - nt
          npoints(max_send_PE) = npoints(max_send_PE) + nt

        END DO  

        sender_PE = s_PE(me)

        RETURN
        END SUBROUTINE FASTJ_LOADBAL_CALC_DIST
      END SUBROUTINE FASTJ_CALC_NEW_NUMBER_OF_POINTS 
#endif
