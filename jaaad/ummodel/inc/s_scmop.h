#if defined(SCMA)
! Start of include file: s_scmop.h
! Description:
!  Declares and defines some parameters necessary for calling SCMoutput
! Current Code Owner: Luke Jones
!
! History:
! Version  Date      Comment
! =======  ====      =======
! 5.5      06/02/03  Original code (Luke Jones)
! 6.0      29/08/03  Removed everything apart from what is necessary
!                    for calling SCMoutput from any routine. Luke Jones.
! 6.2      29/11/05  Add d_cloud for cloud_levels profile.  Add
!                    L_SCMDiags logicals for diagnostics packages
!                                                  A Kerr-Munslow.
!
!
      ! Integers to represent the different time profiles. All must
      ! be non-negative and less than "only_radsteps".
      integer, parameter ::     &
     &     t_inst=1,            & ! Give the instantaneous value
     &     t_avg =2,            & ! Construct the average value
     &     t_max =3,            & ! " the maximum value
     &     t_min =4,            & ! " the minimum value
     &     t_acc =5,            & ! " the accumulated value
     &     t_div =7,            & ! " the average value divided by
                                  ! another diagnostic
     &     t_mult=8,            & ! " the average value multiplied by
                                  ! another diagnostic
     &     t_acc_div=9,         & ! " the accumulated value divided by
                                  ! another diagnostic
     &     t_acc_mult=10,       & ! " the accumulated value multiplied by
                                  ! another diagnostic
     &     t_const=11,          & ! The value is constant.
     &     only_radsteps=100      ! When added to one of the above
                                  ! parameters, flags that the diagnostic
                                  ! is only available on radiation
                                  ! timesteps

      ! Integers to represent the different domain profiles
      integer, parameter ::                                             &
     &     d_sl=1,                                                      &
     &     d_soilt=2,                                                   &
     &     d_bl=3,                                                      &
     &     d_wet=4,                                                     &
     &     d_all=5,                                                     &
     &     d_soilm=6,                                                   &
     &     d_tile=7,                                                    &
     &     d_vis=9,                                                     &
     &     d_point=13,                                                  &
     &     d_allxtra=14,                                                &
     &     d_land=15,                                                   &
     &     d_cloud=16

      ! Statement function to encode a stream number into an integer
      integer Stream,strm
      Stream(strm)=2**(strm-1)

      ! The default streams for diagnostics to go to will be 1,2,3,4,5
      ! and 6. The following should thus be equal to
      ! Stream(1) [2^0=1] + Stream(2) [2^1=2]  + Stream(3) [2^2=4]
      ! Stream(4) [2^3=8] + Stream(5) [2^4=16] + Stream(6) [2^5=32]
      ! Total = 63
      ! where Stream() is the statement function defined above.
      ! default is 63 (all)
      integer, parameter :: default_streams=63

      ! Integers to represent the different diagnostics packages
      Integer, Parameter ::      &
     &     SCMDiag_gen   = 1,    & ! General diagnostics
     &     SCMDiag_rad   = 2,    & ! Radiation
     &     SCMDiag_bl    = 3,    & ! Boundary layer
     &     SCMDiag_surf  = 4,    & ! Surface
     &     SCMDiag_land  = 5,    & ! Land points only
     &     SCMDiag_sea   = 6,    & ! Sea points only
     &     SCMDiag_lsp   = 7,    & ! Large scale precip
     &     SCMDiag_conv  = 8,    & ! Convection
     &     SCMDiag_lscld = 9,    & ! Large scale cloud
     &     SCMDiag_pc2   = 10,   & ! PC2
     &     SCMDiag_forc  = 11,   & ! Forcing
     &     SCMDiag_incs  = 12      ! Increments

! End of include file: s_scmop.h
#endif
