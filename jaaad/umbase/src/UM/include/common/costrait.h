!LL COMDECK COSTRAIT ---------------------------------------------
!LL
!LL Contains dimensioning parameters and small arrays required for
!LL generalised ocean straits exchange flow parametrisation.
!LL
!LL History:
!LL
!LL  Version   Date     Modification history
!LL   5.4    29/08/02   New COMDECK created. D. Storkey
!LL   6.2    11/08/05   Fix to continuations. P.Selwood.
!LL---------------------------------------------------------------
!
! These parameters used here and also in the routines O_READ_STRAIT
! and O_SET_STRAIT.
      INTEGER,PARAMETER:: max_strait=10 ! maximum no of straits
      INTEGER,PARAMETER:: kmax=100      ! maximum no of levels
      INTEGER,PARAMETER:: ntmax=20      ! maximum no of tracers
      INTEGER,PARAMETER:: max_months=12 ! maximum no of months (=12!!)

      INTEGER :: i_strait_u(max_strait,2) ! i coord of strait points
                                          ! on velocity grid
      INTEGER :: j_strait_u(max_strait,2) ! global j coord of strait
                                          ! points on velocity grid
      INTEGER :: i_strait_t(max_strait,4) ! i coord of strait points
                                          ! on tracer grid
      INTEGER :: j_strait_t(max_strait,4) ! global j coord of strait
                                          ! points on tracer grid
      INTEGER :: face_strait(max_strait,2) ! which way does end of strait
                                           ! face? 1=S,2=W,3=N,4=E
      INTEGER :: indx_strait(max_strait,2) ! updating index
      INTEGER :: nt_strait(max_strait)     ! no of tracers advected
                                           ! through strait
      INTEGER :: kmax_strait(max_strait)   ! (effective) depth of strait
                                           ! in model levels
      INTEGER :: iproc_strait_u(max_strait,2)! iproc for velocity
                                             ! grid points
      INTEGER :: iproc_strait_t(max_strait,4)! iproc for tracer
                                             ! grid points

      COMMON/COSTRAIT/ i_strait_u, j_strait_u, i_strait_t, j_strait_t,  &
     &                 face_strait, indx_strait, nt_strait,             &
     &                 kmax_strait, iproc_strait_u, iproc_strait_t
!LL---------------------------------------------------------------
