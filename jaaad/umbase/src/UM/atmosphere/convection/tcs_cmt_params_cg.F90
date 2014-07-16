#if defined(A05_4A) || defined(A05_5A)
MODULE tcs_cmt_params_cg

! Description:
!  Contains congestus CMT settings for turbulence CMT scheme
!
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------

! Parameters for calculation of cloud base stress

  Real, parameter ::                &
    gamma_cmt_cong = 1.0/2.3        & ! Value from fit to CRM results 
   ,delta_cmt_cong = (1.0-1.63/2.3)   ! value from fit to CRM results
    

END MODULE tcs_cmt_params_cg
#endif
