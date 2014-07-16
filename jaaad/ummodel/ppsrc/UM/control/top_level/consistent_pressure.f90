
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Given the exner field at rho levels derive all other pressure fields
!
      SUBROUTINE Consistent_Pressure (                                  &
     &           exner_rho_levels,                                      &
     &           offx,offy,halo_i,halo_J,                               &
     &           row_length,rows,model_levels,                          &
     &           kappa, g, r_theta_levels, r_rho_levels, rho, p_zero,   &
     &           p, pstar, p_theta_levels,exner_theta_levels)          

      IMPLICIT NONE

      Integer, Intent(IN) :: row_length
      Integer, Intent(IN) :: rows
      Integer, Intent(IN) :: model_levels
      Integer, Intent(IN) :: offx        
      Integer, Intent(IN) :: offy        
      Integer, Intent(IN) :: halo_i      
      Integer, Intent(IN) :: halo_j      

      Real, Intent(IN) :: exner_rho_levels(1-offx:row_length+offx,      &
     &                    1-offy:rows+offy, model_levels+1)
      Real, Intent(IN) :: kappa
      Real, Intent(IN) :: g    
      Real, Intent(IN) :: p_zero
      Real, Intent(IN) :: rho(1-offx:row_length+offx,                   &
     &                        1-offy:rows+offy, model_levels)
      Real, Intent(IN) :: r_rho_levels(1-halo_i:row_length+halo_i,      &
     &                              1-halo_j:rows+halo_j, model_levels)
      Real, Intent(IN) :: r_theta_levels(1-halo_i:row_length+halo_i,    &
     &                            1-halo_j:rows+halo_j, 0:model_levels)

      Real, Intent(OUT) :: p(1-offx:row_length+offx,                    &
     &                       1-offy:rows+offy, model_levels+1)
      Real, Intent(OUT) :: pstar(row_length,rows)
      Real, Intent(OUT) :: p_theta_levels(1-offx:row_length+offx,       &
     &                                  1-offy:rows+offy, model_levels) 
      Real, Intent(OUT) :: exner_theta_levels(1-offx:row_length+offx,   &
     &                                  1-offy:rows+offy, model_levels) 

      Integer :: i,j,k
      Real :: recip_kappa

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

! Calculate pressure from Exner.
      recip_kappa = 1./ kappa
      Do k = 1, model_levels+1
        Do j = 1, rows
!CDIR NODEP
          Do i = 1, row_length
            p(i,j,k)= (exner_rho_levels(i,j,k)**recip_kappa) * p_zero
          End Do
        End Do
      End Do

! Halos updated 
! DEPENDS ON: swap_bounds
      call Swap_Bounds(P,                                               &
     &                 row_length, rows, model_levels+1,                &
     &                 offx, offy, fld_type_p, .false.)

! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(P,row_length, rows,                      &
     &                     model_levels+1,offx,offy)

! DEPENDS ON: calc_p_star
      Call Calc_P_star(                                                 &
     &                   r_theta_levels, r_rho_levels,                  &
     &                   P, RHO,                                        &
     &                   g, row_length, rows, model_levels,             &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   PSTAR )

! DEPENDS ON: calc_exner_at_theta
      Call Calc_Exner_at_theta(                                         &
     &                   r_theta_levels, r_rho_levels,                  &
     &                   EXNER_RHO_LEVELS,                              &
     &                   row_length, rows, model_levels,                &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   EXNER_THETA_LEVELS, .FALSE.)

! Calculate pressure from Exner at theta levels.
! DEPENDS ON: calc_p_from_exner
      call Calc_P_from_Exner(                                           &
     &                   P_THETA_LEVELS, kappa, p_zero,                 &
     &                   row_length, rows, model_levels,                &
     &                   offx, offy,                                    &
     &                   EXNER_THETA_LEVELS,.FALSE.)

      RETURN
      END SUBROUTINE Consistent_pressure

