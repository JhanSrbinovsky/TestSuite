! --------------------- COMDECK LBC_COUP --------------------------
!    Description:
!       This COMDECK stores the variables connected with the
!       Parallel Running between the Global and Mesoscale models.
!       The Mesoscale has to be held back until there are sufficient
!       Boundary Conditions (BCs) to proceed.
!
!   History:
!
!   Model    Date     Modification history
!  version
!   4.5    17/08/98 New comdeck. D. Robinson
!
      logical l_lbc_coup        ! T : Global/Meso Coupling on
      integer um_lbc_stream     ! Output Stream generating BCs.
      
      ! UM 6.5 -  MODEL_ANALYSIS_HRS changed to REAL - 
      !                requires LBC_FC_HRS to REAL also
      real    lbc_fc_hrs        ! Forecast time w.r.t analysis time
      integer um_lbc_wait       ! Wait time between re-tries if BCs
                                ! not available.
      integer um_lbc_wait_max   ! Maximum wait time.
      character*80 lbc_filename ! Name of file that communicates between
                                ! Global and Meso.

      COMMON /LBC_COUP/ l_lbc_coup, um_lbc_stream, lbc_fc_hrs,          &
     &                  um_lbc_wait, um_lbc_wait_max, lbc_filename
! ----------------------------------------------------------------------
