
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initiate cloud and liquid water within the PC2 Cloud Scheme.

      Subroutine pc2_initiation_ctl (                                   &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y                                    &
!
! model dimensions.
     &, row_length, rows, rhc_row_length, rhc_rows                      &
     &, model_levels, wet_model_levels                                  &
!
! Model switches
     &, Ltimer, l_mixing_ratio                                          &
     &, L_ACF_CUSACK                                                    &
!
! in time stepping information.
     &, timestep                                                        &
!
! SCM diagnostics switches
     &, nSCMDpkgs,L_SCMDiags                                            &
!
! Primary fields passed IN/OUT
     &,  t,q,qcl,qcf,cf,cfl,cff,rhts,tlts,qtts,ptts,CF_area             &
!
! Primary fields passed IN
     &, p,pstar,p_theta_levels,ccb,cumulus,rhcrit                       &
!
! Output increments
     &, t_work,q_work,qcl_work,qcf_work,cf_work,cfl_work,cff_work       &
!
     & )

      Implicit None
!
! Description:
!   Initiate a small amount of cloud fraction and liquid
!   water content for the PC2 Cloud Scheme. Check that moisture
!   variables are consistent with each other.
!
! Method:
!   See the PC2 documentation.
!
! Current Code Owner:  Damian Wilson.
!
! History:
! Version   Date     Comment.
! -------   ----     --------
! 5.4       22/07/02 Original code.  Damian Wilson.
! 5.5       16/12/02 Tidy up unused code and bring up to date with
!                    latest PC2 formulation. Damian Wilson
! 6.1       27/07/04 Pass convective cloud information to pc2_initiate
!                    and correct pressure dimension error  Damian Wilson
! 6.2       23/01/06 Improvements to Single Column Model Diagnostics
!                    System                          A. Kerr-Munslow
! 6.2       16/02/06 Note model levels issue. A. Kerr-Munslow
! 6.2       31/03/05 Pass in RH(TL) to pc2_initiate   Damian Wilson
! 6.4       18/08/06 Use mixing ratio formulation.  Damian Wilson
! 7.3       25/11/09 Added call to new pc2_arcld routine   Ian Boutle
!
! Code Description:
!   Language: fortran 77 + common extensions
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
! Arguments with intent in. ie: input variables.
!
! Parallel setup variables
      Integer                                                           &
     &  halo_i                                                          &
                   ! Size of halo in i direction.
     &, halo_j                                                          &
                   ! Size of halo in j direction.
     &, off_x                                                           &
                   ! Size of small halo in i
     &, off_y      ! Size of small halo in j.
!
! Model dimensions
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, rhc_row_length                                                  &
                         ! Dimensions of RHcrit variable
     &, rhc_rows                                                        &
                         ! Dimensions of RHcrit variable
     &, model_levels                                                    &
     &, wet_model_levels
! Note that in the full UM pc2_initiation_ctl is called from
! pc2_pressure_forcing, where the only vertical dimension available
! is wet_levels, so that model_levels here is actually just a copy
! of wet_levels.  In the SCM pc2_initiation_ctl is called from
! scm_main, where model_levels exists, and so model_levels here is
! as expected.
!
      Logical                                                           &
     &  Ltimer                                                          &
                 ! true then output some timing information
     &, cumulus(row_length,rows)                                        &
                                  ! convection is occurring
     &, L_ACF_CUSACK                                                    &
     &, l_mixing_ratio            ! Use mixing ratio formulation        
!
! time information for current timestep
      Real                                                              &
     &  timestep
!
! Primary fields passed in/out
      Real                                                              &
     &  T(row_length, rows, model_levels)                               &
     &, q(row_length, rows, wet_model_levels)                           &
     &, qcl(row_length, rows, wet_model_levels)                         &
     &, qcf(row_length, rows, wet_model_levels)                         &
     &, cf(row_length, rows, wet_model_levels)                          &
     &, cfl(row_length, rows, wet_model_levels)                         &
     &, cff(row_length, rows, wet_model_levels)                         &
     &, CF_area(row_length, rows, wet_model_levels)
!
! Primary fields passed in
      Real                                                              &
     &  p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &    wet_model_levels+1)                                           &
     &, pstar(row_length, rows)                                         &
     &, p_theta_levels(row_length, rows, wet_model_levels)              &
     &, rhcrit(rhc_row_length,rhc_rows,wet_model_levels)                &
     &, rhts(row_length,rows,wet_model_levels)                          &
     &, tlts(row_length,rows,wet_model_levels)                          &
!       TL at start of timestep
     &, qtts(row_length,rows,wet_model_levels)                          &
!       qT at start of timestep
     &, ptts(row_length,rows,wet_model_levels)
!       pressure at theta levels at start of timestep
!
      Integer                                                           &
     &  ccb(row_length,rows)
!
      Integer                                                           &
     &  nSCMDpkgs             ! No of SCM diagnostics packages

      Logical                                                           &
     &  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages
!
! Local variables
!

! Output increment diagnostics
      Real                                                              &
     &  t_work(row_length,rows,model_levels)                            &
     &, q_work(row_length,rows,wet_model_levels)                        &
     &, qcl_work(row_length,rows,wet_model_levels)                      &
     &, qcf_work(row_length,rows,wet_model_levels)                      &
     &, cf_work(row_length,rows,wet_model_levels)                       &
     &, cfl_work(row_length,rows,wet_model_levels)                      &
     &, cff_work(row_length,rows,wet_model_levels)                      &
     &, p_layer_boundaries(row_length, rows, 0:model_levels)            &
              ! pressure at layer boundaries. Same as p except at
              ! bottom level = pstar, and at top = 0.
     &, p_layer_centres(row_length, rows, 0:model_levels)
              ! pressure at layer centres. Same as p_theta_levels
              !except bottom level = pstar, and at top = 0.

      Integer                                                           &
     & i,j,k                                                            &
     &, large_levels                                                    &
!       Total no. of sub-levels being processed by cloud scheme.
!       Currently ((wet_model_levels - 2)*levels_per_level) + 2
     &, levels_per_level
!       No. of sub-levels being processed by area cloud scheme.
!       Want an odd number of sublevels per level.
!       NB: levels_per_level = 3 is currently hardwired in the do loops
! Loop counters

      Character*(*), Parameter ::  RoutineName = 'pc2_initiation_ctl'

!
! External Functions:
!
!- End of header
!




!
! Call timer for PC2 initiation code
! DEPENDS ON: timer
      If (Ltimer) Call timer ('PC2 Initiation',3)
!
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
!
! Work fields are set to starting fields so diagnostics
! can be calculated
!
              q_work(i,j,k)   = q(i,j,k)
              qcl_work(i,j,k) = qcl(i,j,k)
              qcf_work(i,j,k) = qcf(i,j,k)
              cf_work(i,j,k)  = cf(i,j,k)
              cfl_work(i,j,k) = cfl(i,j,k)
              cff_work(i,j,k) = cff(i,j,k)
            End Do
          End Do
        End Do
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              T_work(i,j,k)   = T(i,j,k)
            End Do
          End Do
        End Do
!
! Call checking routine
!
! DEPENDS ON: pc2_checks
        CALL PC2_CHECKS(p_theta_levels,wet_model_levels,                &
     &      row_length, rows, t, cf, cfl, cff, q, qcl, qcf,             &
     &      l_mixing_ratio)
!
! Call initiation routine
!
! use Cusack interpolation option
        IF (L_ACF_CUSACK) THEN
!
! set p at layer boundaries.
          Do j = 1, rows
            Do i = 1, row_length
              p_layer_boundaries(i,j,0) = pstar(i,j)
              p_layer_centres(i,j,0) = pstar(i,j)
            End Do
          End Do
          Do k = 1, model_levels - 1
            Do j = 1, rows
              Do i = 1, row_length
                p_layer_boundaries(i,j,k) = p(i,j,k+1)
                p_layer_centres(i,j,k) = p_theta_levels(i,j,k)
              End Do
            End Do
          End Do
          Do j = 1, rows
            Do i = 1, row_length
              p_layer_boundaries(i,j,model_levels) = 0.0
              p_layer_centres(i,j,model_levels) =                       &
     &                           p_theta_levels(i,j,model_levels)
            End Do
          End Do
!
! Determine number of sublevels for vertical gradient area cloud
! Want an odd number of sublevels per level: 3 is hardwired in do loops
          levels_per_level = 3
          large_levels = ((model_levels - 2)*levels_per_level) + 2
!
! DEPENDS ON: pc2_arcld
          CALL PC2_ARCLD(p_layer_centres,p_layer_boundaries,            &
     &      ccb,cumulus,rhcrit,wet_model_levels,                        &
     &      row_length, rows, rhc_row_length,rhc_rows,                  &
     &      large_levels,levels_per_level,CF_area,                      &
     &      t,cf,cfl,cff,q,qcl,qcf,rhts,tlts,qtts,ptts,l_mixing_ratio)
!
        ELSE !l_acf_cusack
!
! DEPENDS ON: pc2_initiate
          CALL PC2_INITIATE(p_theta_levels,ccb,cumulus,rhcrit,          &
     &      wet_model_levels, row_length, rows, rhc_row_length,rhc_rows,&
     &      t,cf,cfl,cff,q,qcl,rhts,l_mixing_ratio)

        END IF !l_acf_cusack
!
! Call second checking routine
!
! DEPENDS ON: pc2_checks2
        CALL PC2_CHECKS2(p_theta_levels,rhcrit,                         &
     &      wet_model_levels, row_length, rows, rhc_row_length,rhc_rows,&
     &      t, cf, cfl, cff, q, qcl, l_mixing_ratio)
!
! Call first checking routine again
!
! DEPENDS ON: pc2_checks
        CALL PC2_CHECKS(p_theta_levels,wet_model_levels,                &
     &      row_length, rows, t, cf, cfl, cff, q, qcl, qcf,             &
     &      l_mixing_ratio)
!
! use Cusack interpolation option
        IF (L_ACF_CUSACK) THEN
! depends on: pc2_hom_arcld
          CALL pc2_hom_arcld(p_layer_centres,p_layer_boundaries,        &
     &       wet_model_levels,row_length,rows,                          &
     &       large_levels,levels_per_level,                             &
     &       CF_area,t,cf,cfl,cff,q,qcl,qcf,                            &
     &       l_mixing_ratio)
        End if
!
! Update work array to hold net increment from the above routines
!
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              q_work(i,j,k)   = q(i,j,k)   - q_work(i,j,k)
              qcl_work(i,j,k) = qcl(i,j,k) - qcl_work(i,j,k)
              qcf_work(i,j,k) = qcf(i,j,k) - qcf_work(i,j,k)
              cf_work(i,j,k)  = cf(i,j,k)  - cf_work(i,j,k)
              cfl_work(i,j,k) = cfl(i,j,k) - cfl_work(i,j,k)
              cff_work(i,j,k) = cff(i,j,k) - cff_work(i,j,k)
            End Do
          End Do
        End Do
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              t_work(i,j,k)   = t(i,j,k)   - t_work(i,j,k)
            End Do
          End Do
        End Do
!
!
! Call timer for PC2 initiation code
! DEPENDS ON: timer
        If (Ltimer) Call timer ('PC2 Initiation',4)
!
! End of routine initiation_ctl
      Return
      END SUBROUTINE pc2_initiation_ctl
