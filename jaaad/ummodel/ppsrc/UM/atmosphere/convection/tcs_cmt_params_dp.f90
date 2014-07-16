
MODULE tcs_cmt_params_dp

! Description:
!  Contains deep CMT settings for turbulence CMT scheme
!
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------

! Parameters for calculation of cloud base stress

  Real, parameter ::       &
    gamma_cmt_deep = 0.3   & ! Value from fit to CRM results 
                             ! (Alan Grant had 0.5)      
   ,delta_cmt_deep = 0.7     ! value from fit to CRM results
                             ! (Alan Grant had 0.185)
    
! Parameters for calculation of wup/wcld

  Real, parameter ::        &
    a_wup_deep = 0.95       & ! Value from fit to CRM results fit WG1 
   ,b_wup_deep = 0.67       & ! value from fit to CRM results      
   ,c_wup_deep = 3.0          ! value from fit to CRM results (guess)     

! Parameters for calculation of K

  Real, parameter ::        &
    a_kcmt_deep = 0.3         ! Value from fit to CRM results (anywhere 
                              ! between (0.2 - 0.5 Or possibly 1.)


END MODULE tcs_cmt_params_dp
