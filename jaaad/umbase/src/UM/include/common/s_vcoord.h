!include file: s_vcoord.h
! Description:
!   Defines vertical coordinates.
!
! Current Code Owner: Z Gardner
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 5.3      27/04/01 Original code. Z Gardner
!
! Declarations:

      INTEGER                                                           &
     &  first_constant_r_rho_level

      Real                                                              &
     &  eta_rho(max_model_levels)                                       &
     &, eta_theta(max_model_levels+1)                                   &
     &, z_top_of_model

!- End of COMDECK declaration

      NAMELIST /VERTLEVS/ Z_TOP_OF_MODEL, FIRST_CONSTANT_R_RHO_LEVEL,   &
     &  ETA_RHO, ETA_THETA

      COMMON /VERTLEVS/ Z_TOP_OF_MODEL, FIRST_CONSTANT_R_RHO_LEVEL,     &
     &  ETA_RHO, ETA_THETA
!---------------------------------------------------------------------
