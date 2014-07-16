! Description: COMDECK containing different preconditioner options
!
! Author : Andy Malcolm
! History:
! Version  Date      Comment.
! 5.3      23/10/01  New comdeck
!
      INTEGER                                                           &
     &     no_precon                                                    &
     &,    vert_precon                                                  &
     &,    vert_plus_xyz_ADI_precon                                     &
     &,    xyz_ADI_precon                                               &
     &,    vert_plus_xz_ADI_precon                                      &
     &,    xz_ADI_precon                                                &
     &,    Dufort_Frankel_precon

      PARAMETER(                                                        &
     &     no_precon               = 0                                  &
     &,    vert_precon             = 1                                  &
     &,    vert_plus_xyz_ADI_precon = 2                                 &
     &,    xyz_ADI_precon          = 3                                  &
     &,    vert_plus_xz_ADI_precon = 4                                  &
     &,    xz_ADI_precon           = 5                                  &
     &,    Dufort_Frankel_precon   = 6 )
